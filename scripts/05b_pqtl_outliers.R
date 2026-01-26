#!/usr/bin/env Rscript
# ==============================================================================
# 05b_pqtl_outliers.R - pQTL-Based Outlier Detection
# ==============================================================================
#
# Purpose:
#   Detects sample mismatches by comparing observed protein levels with
#   genotype-predicted protein levels using pQTL variants. Calculates Z-scores
#   per protein-pQTL pair and uses Mean Absolute Z-score (MeanAbsZ) as the
#   primary metric for outlier assignment. Uses PCA-cleaned matrix for robust
#   comparisons.
#
# Author: Reza Jabal, PhD (rjabal@broadinstitute.org)
# Date: December 2025
# ==============================================================================

suppressPackageStartupMessages({
    library(data.table)
    library(tidyverse)
    library(ggplot2)
    library(yaml)
    library(logger)
    library(stringr)
    library(digest)
    library(jsonlite)
})

# Ensure PLINK is in PATH (adding common locations and user's conda env)
old_path <- Sys.getenv("PATH")
Sys.setenv(PATH = paste(old_path, "/home/rjabal/miniconda3/envs/gwas_env/bin", "/mnt/longGWAS_disk_100GB/long_gwas/bin", sep = ":"))

# Source path utilities
script_dir <- tryCatch(
    {
        cmd_args <- commandArgs(trailingOnly = FALSE)
        file_arg <- grep("^--file=", cmd_args, value = TRUE)
        if (length(file_arg) > 0) {
            script_path <- sub("^--file=", "", file_arg)
            dirname(normalizePath(script_path))
        } else {
            getwd()
        }
    },
    error = function(e) getwd()
)
source(file.path(script_dir, "path_utils.R"))

# Get config path
config_file <- Sys.getenv("PIPELINE_CONFIG", "")
if (config_file == "" || !file.exists(config_file)) {
  stop("PIPELINE_CONFIG environment variable not set or config file not found. Please provide path to config file.")
}
config <- read_yaml(config_file)

# Get batch context
batch_id <- Sys.getenv("PIPELINE_BATCH_ID", config$batch$default_batch_id %||% "batch_01")
step_num <- get_step_number()

# Set up logging
log_path <- get_log_path(step_num, batch_id, config)
ensure_output_dir(log_path)
log_appender(appender_file(log_path))
log_info("Starting pQTL-based outlier detection for batch: {batch_id}")

# Helper to list GS files
list_gs_files <- function(path, pattern = NULL) {
    cmd <- paste("gsutil ls", path)
    files <- tryCatch(system(cmd, intern = TRUE), error = function(e) character(0))
    if (!is.null(pattern)) {
        files <- files[grepl(pattern, files)]
    }
    return(files)
}

# Helper to read GS file (copies to temp)
read_gs_file <- function(gs_path) {
    tmp <- tempfile()
    # Register cleanup
    on.exit(if (file.exists(tmp)) unlink(tmp), add = TRUE)

    system(paste("gsutil cp", gs_path, tmp))
    if (!file.exists(tmp)) {
        return(NULL)
    }

    dt <- fread(tmp)
    return(dt)
}

# Helper function for Inverse Rank Normalization (IRN) per protein
# Transforms each protein column to follow N(0,1) distribution
# This makes z-score calculation more robust to outliers and non-normal distributions
inverse_rank_normalize <- function(npx_matrix) {
    log_info("Applying Inverse Rank Normalization (IRN) per protein...")

    # Store original dimensions and names
    original_rownames <- rownames(npx_matrix)
    original_colnames <- colnames(npx_matrix)

    # Convert to matrix if not already
    if (!is.matrix(npx_matrix)) {
        npx_matrix <- as.matrix(npx_matrix)
    }

    # Apply IRN column-wise (per protein)
    npx_irn <- apply(npx_matrix, 2, function(x) {
        # Get non-NA values
        valid_idx <- !is.na(x)
        n_valid <- sum(valid_idx)

        # If too few valid values, return original (with NAs)
        if (n_valid < 3) {
            return(x)
        }

        # Extract valid values
        x_valid <- x[valid_idx]

        # Calculate ranks (average for ties)
        # rank() handles ties by default with average method
        ranks <- rank(x_valid, ties.method = "average", na.last = "keep")

        # Convert ranks to quantiles: (rank - 0.5) / n
        # This ensures quantiles are in (0, 1) range
        quantiles <- (ranks - 0.5) / n_valid

        # Apply inverse normal CDF (qnorm) to get N(0,1) distributed values
        # Handle edge cases: quantiles exactly 0 or 1 would give -Inf/Inf
        # Use small epsilon to avoid infinite values
        epsilon <- 1e-6
        quantiles <- pmax(epsilon, pmin(1 - epsilon, quantiles))
        x_irn_valid <- qnorm(quantiles)

        # Create output vector with same length as input
        x_irn <- rep(NA_real_, length(x))
        x_irn[valid_idx] <- x_irn_valid

        return(x_irn)
    })

    # Ensure matrix format and restore names
    npx_irn <- as.matrix(npx_irn)
    rownames(npx_irn) <- original_rownames
    colnames(npx_irn) <- original_colnames

    # Log summary statistics
    npx_irn_flat <- as.vector(npx_irn)
    npx_irn_flat <- npx_irn_flat[!is.na(npx_irn_flat)]
    if (length(npx_irn_flat) > 0) {
        log_info("IRN transformation complete:")
        log_info("  Mean: {round(mean(npx_irn_flat), 4)} (expected ~0)")
        log_info("  SD: {round(sd(npx_irn_flat), 4)} (expected ~1)")
        log_info("  Range: [{round(min(npx_irn_flat), 4)}, {round(max(npx_irn_flat), 4)}]")
    }

    return(npx_irn)
}

# Helper: Extract variant IDs from genotype data
extract_variant_ids <- function(dt_geno, local_genotype_path) {
    # Try .bim file first (most reliable)
    bim_file <- paste0(local_genotype_path, ".bim")
    if (file.exists(bim_file)) {
        dt_bim <- fread(bim_file, col.names = c("chr", "rsid", "cm", "pos", "a1", "a2"))
        return(unique(dt_bim$rsid))
    }

    # Fallback: extract from column names
    standard_cols <- c("FID", "IID", "PAT", "MAT", "SEX", "PHENOTYPE")
    variant_cols <- setdiff(names(dt_geno), standard_cols)
    if (length(variant_cols) > 0) {
        rsids <- unique(gsub("_[AGCT]$", "", variant_cols))
        return(rsids)
    }

    return(character(0))
}

# Helper: Apply MAF filter with logging
apply_maf_filter <- function(variants, maf_threshold, context = "selection") {
    if (!("maf" %in% names(variants))) {
        log_warn("MAF column not found. Skipping MAF filter for {context}.")
        return(variants)
    }

    n_before <- nrow(variants)
    variants_filtered <- variants[!is.na(maf) & maf > maf_threshold]
    n_after <- nrow(variants_filtered)
    n_filtered_out <- n_before - n_after

    log_info("Applied MAF > {maf_threshold*100}% filter for {context}:")
    log_info("  Variants before: {n_before}, after: {n_after}, filtered out: {n_filtered_out}")

    if (n_after > 0 && context == "z-score calculation") {
        log_info("  P-value range: {sprintf('%.2e', min(variants_filtered$p, na.rm=TRUE))} - {sprintf('%.2e', max(variants_filtered$p, na.rm=TRUE))}")
        log_info("  Beta range: {round(min(abs(variants_filtered$beta), na.rm=TRUE), 4)} - {round(max(abs(variants_filtered$beta), na.rm=TRUE), 4)}")
        log_info("  MAF range: {round(min(variants_filtered$maf, na.rm=TRUE), 4)} - {round(max(variants_filtered$maf, na.rm=TRUE), 4)}")
        if ("heterozygosity" %in% names(variants_filtered)) {
            log_info("  Heterozygosity range: {round(min(variants_filtered$heterozygosity, na.rm=TRUE), 4)} - {round(max(variants_filtered$heterozygosity, na.rm=TRUE), 4)}")
        }
    }

    return(variants_filtered)
}

# Helper: Load consensus/trained pQTLs or provided pQTL list
load_trained_pqtls <- function(pqtl_config, config, batch_id) {
    use_consensus_pqtls <- tryCatch(isTRUE(pqtl_config$use_consensus_pqtls), error = function(e) FALSE)
    training_config <- config$pqtl_training
    trained_pqtls <- NULL
    trained_thresholds <- NULL

    # First, check for provided pQTL list file (when training is disabled)
    provided_pqtl_file <- tryCatch(pqtl_config$provided_pqtl_list_file, error = function(e) NULL)
    if (!is.null(provided_pqtl_file) && file.exists(provided_pqtl_file)) {
        log_info("Loading provided pQTL list from: {provided_pqtl_file}")
        provided_pqtls_dt <- tryCatch(fread(provided_pqtl_file), error = function(e) NULL)
        if (!is.null(provided_pqtls_dt) && nrow(provided_pqtls_dt) > 0) {
            # Handle different column names (rsid, variant, variant_id, etc.)
            if ("rsid" %in% names(provided_pqtls_dt)) {
                trained_pqtls <- data.table(rsid = provided_pqtls_dt$rsid)
            } else if ("variant" %in% names(provided_pqtls_dt)) {
                trained_pqtls <- data.table(rsid = provided_pqtls_dt$variant)
            } else if ("variant_id" %in% names(provided_pqtls_dt)) {
                trained_pqtls <- data.table(rsid = provided_pqtls_dt$variant_id)
            } else if (ncol(provided_pqtls_dt) == 1) {
                # Single column - assume it's rsid
                trained_pqtls <- data.table(rsid = unlist(provided_pqtls_dt[, 1]))
            } else {
                log_warn("Provided pQTL file has unexpected format. Expected 'rsid', 'variant', or 'variant_id' column, or single column. Falling back to other sources.")
            }
            if (!is.null(trained_pqtls) && nrow(trained_pqtls) > 0) {
                log_info("Loaded {nrow(trained_pqtls)} pQTLs from provided list")
                return(list(pqtls = trained_pqtls, thresholds = trained_thresholds))
            }
        } else {
            log_warn("Provided pQTL file exists but could not be read or is empty: {provided_pqtl_file}")
        }
    }

    # Second, check for consensus pQTLs from step 05a training
    if (use_consensus_pqtls) {
        training_output_dir <- file.path(config$output$base_dir, training_config$output_dir %||% "pqtl-training", batch_id)
        consensus_pqtls_path <- file.path(training_output_dir, "05_05a_consensus_pqtls.tsv")
        consensus_config_path <- file.path(training_output_dir, "05_05a_consensus_config.yaml")

        # Fallback to legacy names
        if (!file.exists(consensus_pqtls_path)) {
            consensus_pqtls_path <- file.path(training_output_dir, training_config$outputs$selected_pqtls %||% "05_05a_selected_pqtls.tsv")
        }
        if (!file.exists(consensus_config_path)) {
            consensus_config_path <- file.path(training_output_dir, training_config$outputs$trained_thresholds %||% "05_05a_trained_thresholds.yaml")
        }

        if (file.exists(consensus_pqtls_path) && file.exists(consensus_config_path)) {
            log_info("Loading consensus pQTLs from step 05a")
            consensus_pqtls_dt <- fread(consensus_pqtls_path)
            if ("rsid" %in% names(consensus_pqtls_dt)) {
                trained_pqtls <- consensus_pqtls_dt
                log_info("Loaded {nrow(trained_pqtls)} consensus pQTLs")
            }

            consensus_config <- read_yaml(consensus_config_path)
            trained_thresholds <- consensus_config$thresholds

            if (!is.null(consensus_config$pqtls) && is.null(trained_pqtls)) {
                trained_pqtls <- data.table(rsid = unlist(consensus_config$pqtls))
                log_info("Loaded {nrow(trained_pqtls)} consensus pQTLs from YAML")
            }
        }
    } else if (!is.null(training_config) && isTRUE(training_config$enabled)) {
        # Legacy format
        training_output_dir <- file.path(config$output$base_dir, training_config$output_dir %||% "pqtl-training", batch_id)
        selected_pqtls_path <- file.path(training_output_dir, training_config$outputs$selected_pqtls %||% "05_05a_selected_pqtls.tsv")
        thresholds_path <- file.path(training_output_dir, training_config$outputs$trained_thresholds %||% "05_05a_trained_thresholds.yaml")

        if (file.exists(selected_pqtls_path) && file.exists(thresholds_path)) {
            log_info("Loading trained parameters (legacy format)")
            trained_pqtls <- fread(selected_pqtls_path)
            trained_thresholds <- read_yaml(thresholds_path)
            log_info("Loaded {nrow(trained_pqtls)} LASSO-selected pQTLs")
        }
    }

    return(list(pqtls = trained_pqtls, thresholds = trained_thresholds))
}

# Helper: Merge sex predictions with pQTL stats (simplified)
merge_sex_predictions <- function(sample_stats, sex_preds) {
    # Normalize column names
    if ("SAMPLE_ID" %in% names(sex_preds) && !("SampleID" %in% names(sex_preds))) {
        setnames(sex_preds, "SAMPLE_ID", "SampleID")
    }

    # Determine merge column
    if ("FINNGENID" %in% names(sample_stats) && "FINNGENID" %in% names(sex_preds)) {
        merge_col <- "FINNGENID"
    } else if ("SampleID" %in% names(sample_stats) && "SampleID" %in% names(sex_preds)) {
        merge_col <- "SampleID"
    } else {
        stop("Cannot merge: missing identifier columns")
    }

    # Select columns to merge
    sex_cols <- c(merge_col, "predicted_prob", "genetic_sex", "mismatch", "sex_outlier")
    available_cols <- intersect(sex_cols, names(sex_preds))

    if (length(available_cols) == 0 || !merge_col %in% available_cols) {
        stop("Cannot merge: required columns missing in sex_preds")
    }

    # Perform merge
    plot_data <- merge(sample_stats, sex_preds[, ..available_cols],
                      by = merge_col, all = TRUE, suffixes = c("", ".sex"))

    # Ensure SampleID exists
    if (!("SampleID" %in% names(plot_data))) {
        if ("FINNGENID" %in% names(plot_data)) {
            plot_data[, SampleID := FINNGENID]
        } else {
            stop("Cannot create SampleID: missing identifier columns")
        }
    }

    return(plot_data)
}

# Main function
main <- function() {
    # 1. Load Inputs
    # ----------------
    log_info("Loading inputs...")

    # Use step 01 cleaned matrix (PCA outliers removed) to match step 04, with fallback to step 00
    # This ensures consistency between sex predictions (step 04) and pQTL stats (step 05b)
        npx_path <- get_output_path("01", "npx_matrix_pca_cleaned", batch_id, "outliers", config = config)
    if (!file.exists(npx_path)) {
        log_warn("step 01 cleaned matrix not found: {npx_path}. Falling back to step 00 matrix.")
        npx_path <- get_output_path("00", "npx_matrix_analysis_ready", batch_id, "qc", config = config)
    }
    if (!file.exists(npx_path)) {
        stop("No input NPX matrix found (tried step 01 and step 00): {npx_path}")
    }
    npx_matrix <- readRDS(npx_path)
    log_info("Using NPX matrix: {npx_path}")
    log_info("Loaded NPX matrix: {nrow(npx_matrix)} samples x {ncol(npx_matrix)} proteins")

    # Load sample mapping to convert Olink IDs to FINNGENIDs
    # The genotype data uses FINNGENIDs, but NPX matrix uses Olink SampleIDs
    sample_mapping_path <- get_output_path("00", "sample_mapping", batch_id, "qc", "tsv", config = config)

    if (!file.exists(sample_mapping_path)) {
        stop("Sample mapping file not found: {sample_mapping_path}. Run step 00 first.")
    }

    dt_mapping <- fread(sample_mapping_path)
    log_info("Loaded sample mapping: {nrow(dt_mapping)} samples")

    # Convert NPX matrix row names from SampleID to FINNGENID
    # Match row names to SampleID column
    idx <- match(rownames(npx_matrix), dt_mapping$SampleID)

    # Filter to samples that have FINNGENIDs
    valid_idx <- which(!is.na(idx) & !is.na(dt_mapping$FINNGENID[idx]))

    if (length(valid_idx) == 0) {
        stop("No samples with FINNGENIDs found. Check sample mapping.")
    }

    npx_matrix <- npx_matrix[valid_idx, ]
    rownames(npx_matrix) <- dt_mapping$FINNGENID[idx[valid_idx]]

    log_info("Mapped to FINNGENIDs: {nrow(npx_matrix)} samples with valid IDs")

    # Get pQTL config
    pqtl_config <- config$parameters$pqtl_outliers
    if (is.null(pqtl_config)) stop("No pqtl_outliers configuration found")

    finemap_path <- pqtl_config$finemap_path
    if (is.null(finemap_path)) stop("finemap_path is not defined in pqtl_config")

    genotype_path <- pqtl_config$genotype_path
    if (is.null(genotype_path) || length(genotype_path) == 0) {
        stop("genotype_path is not defined in pqtl_config. Check config file.")
    }

    z_threshold <- pqtl_config$z_threshold %||% 4

    log_info("Genotype path: {genotype_path}")

    # Apply Inverse Rank Normalization (IRN) per protein if enabled
    # IRN transforms each protein to N(0,1) distribution, making z-score calculation more robust
    apply_irn <- tryCatch(
        pqtl_config$apply_irn %||% TRUE,  # Default: enabled
        error = function(e) TRUE
    )

    if (isTRUE(apply_irn)) {
        npx_matrix <- inverse_rank_normalize(npx_matrix)
        log_info("IRN applied: NPX matrix transformed to N(0,1) per protein")
    } else {
        log_info("IRN disabled: Using raw NPX values for z-score calculation")
    }

    log_info("Configuration:")
    log_info("  Fine-mapping path: {finemap_path}")
    log_info("  Genotype path: {genotype_path}")
    log_info("  Z-score threshold: {z_threshold}")

    # =========================================================================
    # Check for consensus pQTLs from step 05a (or trained parameters)
    # =========================================================================
    trained_result <- load_trained_pqtls(pqtl_config, config, batch_id)
    trained_pqtls <- trained_result$pqtls
    trained_thresholds <- trained_result$thresholds
    use_trained_params <- !is.null(trained_pqtls) && nrow(trained_pqtls) > 0

    if (use_trained_params && !is.null(trained_thresholds)) {
        log_info("Using trained/consensus pQTLs: {nrow(trained_pqtls)} variants")
        if (!is.null(trained_thresholds$mean_abs_z)) {
            if (!is.null(trained_thresholds$mean_abs_z$k_sd)) {
                log_info("  MeanAbsZ: k_sd = {round(trained_thresholds$mean_abs_z$k_sd, 4)}")
            } else if (!is.null(trained_thresholds$mean_abs_z$threshold)) {
                log_info("  MeanAbsZ: threshold = {round(trained_thresholds$mean_abs_z$threshold, 4)}")
                z_threshold <- trained_thresholds$mean_abs_z$threshold
            }
        }
    }

    # 2. Parse Fine-mapping Results
    # -----------------------------
    log_info("Parsing fine-mapping results...")

    # Create a collated file if not exists
    collated_path <- get_output_path(step_num, "05b_finemap_collated", batch_id, "pqtl", "tsv", config = config)
    ensure_output_dir(collated_path)

    if (file.exists(collated_path) && file.size(collated_path) > 0) {
        log_info("Loading existing collated fine-mapping results from: {collated_path}")
        # Explicitly specify column classes to prevent fread from misinterpreting values as dates
        # This is critical to prevent POSIXt errors in downstream calculations
        top_variants <- fread(collated_path, colClasses = list(
            character = c("trait", "rsid"),
            numeric = c("p", "beta", "sd", "cs_avg_r2", "maf", "heterozygosity")
        ))
        log_info("Loaded {nrow(top_variants)} variants from cache")

        # Ensure numeric columns are properly typed (prevent POSIXt issues)
        # Convert any non-numeric columns that should be numeric
        # CRITICAL: fread may auto-detect certain values as dates (POSIXt)
        numeric_cols <- c("p", "beta", "sd", "cs_avg_r2", "maf", "heterozygosity")
        for (col in numeric_cols) {
            if (col %in% names(top_variants)) {
                col_class <- class(top_variants[[col]])[1]
                col_is_numeric <- is.numeric(top_variants[[col]]) && !inherits(top_variants[[col]], "POSIXt")
                # Force conversion if not properly numeric (character, POSIXt, or any other non-numeric type)
                if (!col_is_numeric) {
                    log_warn("Column {col} has unexpected type ({col_class}). Converting to numeric.")
                    # For POSIXt, as.numeric gives seconds since epoch - we need to handle this specially
                    if (inherits(top_variants[[col]], "POSIXt")) {
                        log_warn("Column {col} was read as date/time. Setting to NA and will be recalculated.")
                        top_variants[, (col) := NA_real_]
                    } else {
                        top_variants[, (col) := as.numeric(get(col))]
                    }
                }
            }
        }
        log_info("Type checking complete for cached columns")
    } else {
        log_info("Collating fine-mapping results from GS...")

    files <- list_gs_files(finemap_path, "\\.SUSIE\\.cred\\.summary\\.tsv$")
    log_info("Found {length(files)} fine-mapping result files")

    if (length(files) == 0) stop("No fine-mapping files found")

    # Function to read and extract top variant
    get_top_variant <- function(f) {
        tryCatch(
            {
                dt <- read_gs_file(f)
                    if (is.null(dt) || nrow(dt) == 0) {
                    return(NULL)
                }

                # Filter for high confidence
                # Sort by p-value (asc) and cs_avg_r2 (desc)
                setorder(dt, p, -cs_avg_r2)

                # Take top 1
                top <- dt[1]

                    # Select informative columns

                    cols <- c("trait", "rsid", "p", "beta", "sd", "cs_avg_r2")

                    # Ensure required columns exist
                    missing <- setdiff(cols, names(top))
                    if (length(missing) > 0) {
                        for (m in missing) top[[m]] <- NA
                    }
                    return(top[, ..cols])
            },
            error = function(e) {
                return(NULL)
            }
        )
    }

        dts <- list()
        for (i in seq_along(files)) {
            if (i %% 50 == 0) log_info("Processed {i}/{length(files)} files")
            res <- get_top_variant(files[i])
            if (!is.null(res)) dts[[i]] <- res
        }
        top_variants <- rbindlist(dts)
        fwrite(top_variants, collated_path, sep = "\t")
        log_info("Saved collated results to: {collated_path}")
    }

    log_info("Identified {nrow(top_variants)} top variants")

    # Note: We DON'T initialize heterozygosity column here.
    # It will be calculated and merged later after genotype extraction.
    # If loaded from cache, the heterozygosity column may already exist with values.
    if ("heterozygosity" %in% names(top_variants)) {
        log_info("Heterozygosity column exists in cached file: {sum(!is.na(top_variants$heterozygosity))} variants with values")
    }

    # 2.5. Apply Top N Selection and Sex Chromosome Handling
    # --------------------------------------------------------
    sex_chr_enabled <- tryCatch(isTRUE(pqtl_config$sex_chromosome_pQTLs), error = function(e) FALSE)
    top_n <- tryCatch(pqtl_config$top_n_pQTLs, error = function(e) NULL)

    # Handle sex chromosome variants (chr23 → chrX/chrx)
    if (sex_chr_enabled) {
        log_info("Sex chromosome pQTL handling enabled: Converting chr23 to chrX/chrx for matching")
        # Count chr23 variants
        chr23_count <- sum(grepl("^chr23_", top_variants$rsid))
        log_info("Found {chr23_count} chr23 variants in fine-mapping results")

        # Create a mapping column for matching (convert chr23_* to chrX_* or chrx_*)
        # We'll try both chrX and chrx when matching
        top_variants[, rsid_original := rsid]
        top_variants[, rsid_for_matching := ifelse(grepl("^chr23_", rsid),
                                                    gsub("^chr23_", "chrX_", rsid),
                                                    rsid)]
    } else {
        top_variants[, rsid_original := rsid]
        top_variants[, rsid_for_matching := rsid]
        # Count chr23 variants that will be skipped
        chr23_count <- sum(grepl("^chr23_", top_variants$rsid))
        if (chr23_count > 0) {
            log_info("Found {chr23_count} chr23 variants (will be skipped unless sex_chromosome_pQTLs is enabled)")
        }
    }

    # 2.6. Unified Genotype Extraction, MAF, and Heterozygosity Calculation
    # ------------------------------------------------------------------------
    # REFACTORED: Extract genotypes ONCE for ALL variants and ALL NPX samples
    # Then calculate MAF and heterozygosity from exported data
    # This eliminates redundant sparse BED creation and PLINK --freq runs

    log_info("=== UNIFIED GENOTYPE EXTRACTION ===")
    log_info("Extracting genotypes for ALL {nrow(top_variants)} variants and ALL {nrow(npx_matrix)} NPX samples...")
    log_info("Will calculate MAF and heterozygosity from exported genotype data")

    # Get all NPX sample FINNGENIDs (for keep file)
    npx_sample_finngenids <- rownames(npx_matrix)
    log_info("NPX matrix samples: {length(npx_sample_finngenids)} samples")

    # Create keep file for PLINK (all NPX samples)
    keep_file <- tempfile(pattern = "npx_samples_", fileext = ".txt")
    keep_dt <- data.table(FID = npx_sample_finngenids, IID = npx_sample_finngenids)
    fwrite(keep_dt, keep_file, sep = "\t", col.names = FALSE)

    # Create variant list file (all variants)
    var_file <- tempfile(pattern = "all_variants_", fileext = ".snplist")
    writeLines(top_variants$rsid_for_matching, var_file)
    log_info("Variant list: {length(top_variants$rsid_for_matching)} variants")

    # Generate cache key based on variant list + sample list + collated file mtime
    variant_hash <- digest::digest(sort(top_variants$rsid_for_matching), algo = "md5")
    sample_hash <- digest::digest(sort(npx_sample_finngenids), algo = "md5")
    collated_mtime <- if (file.exists(collated_path)) file.mtime(collated_path) else Sys.time()
    cache_key <- paste0(variant_hash, "_", sample_hash, "_", as.integer(collated_mtime))

    # Check for cached genotype files and exported .raw file
    base_dir <- config$output$base_dir %||% Sys.getenv("PIPELINE_OUTPUT_DIR", "output")
    cache_base <- file.path(base_dir, "temp_work", "pqtl_cache", batch_id)
    if (!dir.exists(cache_base)) dir.create(cache_base, recursive = TRUE)

    cache_dir <- file.path(cache_base, cache_key)

    # Safety check: ensure genotype_path is defined before using it
    if (!exists("genotype_path") || is.null(genotype_path)) {
        stop("genotype_path is not defined. This should have been set earlier in the script.")
    }

    base_name <- basename(genotype_path)

    cached_bed <- file.path(cache_dir, paste0(base_name, ".bed"))
    cached_bim <- file.path(cache_dir, paste0(base_name, ".bim"))
    cached_fam <- file.path(cache_dir, paste0(base_name, ".fam"))
    cached_raw <- file.path(cache_dir, paste0(base_name, "_exported.raw"))

    cache_hit <- (file.exists(cached_bed) && file.exists(cached_bim) && file.exists(cached_fam) &&
                  file.size(cached_bed) > 0 && file.size(cached_bim) > 0 && file.size(cached_fam) > 0)
    exported_cache_hit <- (file.exists(cached_raw) && file.size(cached_raw) > 0)

    local_genotype_path <- NULL
    temp_geno_dir <- NULL
    dt_geno_exported <- NULL

    if (cache_hit && exported_cache_hit) {
        log_info("Cache HIT: Found cached sparse BED and exported genotypes")
        log_info("Cache directory: {cache_dir}")
        local_genotype_path <- file.path(cache_dir, base_name)

        # Load exported genotypes directly
        log_info("Loading cached exported genotypes from: {cached_raw}")
        dt_geno_exported <- fread(cached_raw)
        log_info("Loaded {nrow(dt_geno_exported)} samples from cached exported genotypes")
    } else {
        log_info("Cache MISS: Will extract genotypes and export to cache")
        log_info("Cache key: {cache_key}")

        # Extract sparse BED (reuse existing logic from main extraction)
        if (startsWith(genotype_path, "gs://")) {
            log_info("Detected GS genotype path. Using Smart Extraction (Sparse Download)...")

            # Create temp dir for genotypes
            large_temp_base <- file.path(config$output$base_dir, "temp_work")
            if (!dir.exists(large_temp_base)) dir.create(large_temp_base, recursive = TRUE)

            temp_geno_dir <- tempfile(pattern = "plink_input_", tmpdir = large_temp_base)
            dir.create(temp_geno_dir)

            # Download .bim and .fam
            log_info("Downloading metadata (.bim, .fam)...")
            bim_gs <- paste0(genotype_path, ".bim")
            fam_gs <- paste0(genotype_path, ".fam")
            bed_gs <- paste0(genotype_path, ".bed")

            bim_local <- file.path(temp_geno_dir, paste0(base_name, ".bim"))
            fam_local <- file.path(temp_geno_dir, paste0(base_name, ".fam"))

            system(paste("gsutil cp", bim_gs, bim_local))
            system(paste("gsutil cp", fam_gs, fam_local))

            if (!file.exists(bim_local) || !file.exists(fam_local)) {
                stop("Failed to download .bim or .fam files. Check paths.")
            }

            # Read .bim and match variants
            log_info("Matching {length(top_variants$rsid_for_matching)} variants in .bim file...")
            dt_bim <- fread(bim_local, header = FALSE)
            setnames(dt_bim, c("chr", "rsid", "cm", "pos", "a1", "a2"))
            dt_bim[, idx := .I]

            target_rsids <- top_variants$rsid_for_matching
            if (sex_chr_enabled) {
                dt_targets <- dt_bim[rsid %in% target_rsids]
                if (nrow(dt_targets) < length(target_rsids)) {
                    missing_targets <- setdiff(target_rsids, dt_targets$rsid)
                    missing_chrx <- gsub("^chrX_", "chrx_", missing_targets)
                    dt_targets_chrx <- dt_bim[rsid %in% missing_chrx]
                    if (nrow(dt_targets_chrx) > 0) {
                        dt_targets_chrx[, rsid := gsub("^chrx_", "chrX_", rsid)]
                        dt_targets <- rbind(dt_targets, dt_targets_chrx)
                    }
                }
            } else {
                dt_targets <- dt_bim[rsid %in% target_rsids]
            }

            if (nrow(dt_targets) == 0) {
                stop("No matching variants found in .bim file. Check RSIDs.")
            }

            log_info("Found {nrow(dt_targets)} matching variants in .bim file")

            # Read .fam to get sample count
            dt_fam <- fread(fam_local, header = FALSE)
            n_samples <- nrow(dt_fam)
            bytes_per_variant <- ceiling(n_samples / 4)

            # Create sparse BED file
            log_info("Downloading variant data (Sparse BED creation) for {nrow(dt_targets)} variants...")
            sparse_bed_local <- file.path(temp_geno_dir, paste0(base_name, ".bed"))

            f_out <- file(sparse_bed_local, "wb")
            writeBin(as.raw(c(0x6c, 0x1b, 0x01)), f_out)

            setorder(dt_targets, idx)
            for (i in seq_len(nrow(dt_targets))) {
                var_idx <- dt_targets[i]$idx
                start_byte <- 3 + (var_idx - 1) * bytes_per_variant
                end_byte <- start_byte + bytes_per_variant - 1
                range_str <- paste0(format(start_byte, scientific = FALSE), "-", format(end_byte, scientific = FALSE))

                chunk_file <- file.path(temp_geno_dir, "chunk.bin")
                cmd <- paste0("gsutil cat -r ", range_str, " ", bed_gs, " > ", chunk_file)
                ret <- system(cmd)

                if (ret == 0 && file.exists(chunk_file) && file.size(chunk_file) == bytes_per_variant) {
                    chunk_data <- readBin(chunk_file, "raw", n = bytes_per_variant + 100)
                    writeBin(chunk_data[1:bytes_per_variant], f_out)
                    unlink(chunk_file)
                    dt_targets[i, downloaded := TRUE]
                }
            }
            close(f_out)

            # Create filtered .bim
            dt_final_bim <- dt_targets[downloaded == TRUE, .(chr, rsid, cm, pos, a1, a2)]
            fwrite(dt_final_bim, file.path(temp_geno_dir, paste0(base_name, ".bim")), sep = "\t", col.names = FALSE)

            local_genotype_path <- file.path(temp_geno_dir, base_name)
            log_info("Sparse dataset created at: {local_genotype_path}")

            # Cache sparse BED files
            if (!dir.exists(cache_dir)) dir.create(cache_dir, recursive = TRUE)
            log_info("Caching sparse BED files to: {cache_dir}")
            file.copy(sparse_bed_local, cached_bed, overwrite = TRUE)
            file.copy(file.path(temp_geno_dir, paste0(base_name, ".bim")), cached_bim, overwrite = TRUE)
            file.copy(fam_local, cached_fam, overwrite = TRUE)
            local_genotype_path <- file.path(cache_dir, base_name)
            log_info("Sparse BED files cached successfully")
        } else {
            # Local genotype path
            local_genotype_path <- genotype_path
            log_info("Using local genotype path: {local_genotype_path}")
        }

        # Run PLINK --export A to export genotypes for all variants and all NPX samples
        log_info("Running PLINK --export A for all variants and all NPX samples...")
        out_prefix <- tempfile(pattern = "geno_export_")

        is_bed <- file.exists(paste0(local_genotype_path, ".bed"))
        flag <- if (is_bed) "--bfile" else "--pfile"

    plink_cmd <- sprintf(
            "plink2 %s %s --extract %s --keep %s --export A --out %s 2>&1",
            flag, local_genotype_path, var_file, keep_file, out_prefix
        )

        log_info("PLINK command: {plink_cmd}")
        plink_result <- system(plink_cmd, intern = TRUE)
        plink_exit_code <- attr(plink_result, "status") %||% 0

        if (plink_exit_code != 0) {
            log_warn("PLINK command failed with exit code: {plink_exit_code}")
            log_warn("PLINK output (last 20 lines): {paste(tail(plink_result, 20), collapse = '\\n')}")
        }

    raw_file <- paste0(out_prefix, ".raw")
        if (!file.exists(raw_file) || file.size(raw_file) == 0) {
            stop("PLINK export failed. Output file not found or empty: {raw_file}")
        }

        # Load exported genotypes
        log_info("Loading exported genotypes from: {raw_file}")
        dt_geno_exported <- fread(raw_file)
        log_info("Loaded {nrow(dt_geno_exported)} samples from exported genotypes")

        # Cache exported .raw file
        if (!dir.exists(cache_dir)) dir.create(cache_dir, recursive = TRUE)
        file.copy(raw_file, cached_raw, overwrite = TRUE)
        log_info("Cached exported genotypes to: {cached_raw}")

        # Cleanup temp files
        if (file.exists(var_file)) unlink(var_file)
        if (file.exists(keep_file)) unlink(keep_file)
        if (file.exists(out_prefix)) unlink(paste0(out_prefix, "*"))
        if (!is.null(temp_geno_dir) && dir.exists(temp_geno_dir)) {
            # Don't cleanup if we cached the files
            if (!cache_hit) {
                unlink(temp_geno_dir, recursive = TRUE)
            }
        }
    }

    # Calculate MAF from exported genotype data
    log_info("=== CALCULATING MAF FROM EXPORTED GENOTYPE DATA ===")
    maf_data <- NULL
    if (!is.null(dt_geno_exported) && nrow(dt_geno_exported) > 0) {
        # Get variant columns (exclude standard PLINK columns)
        standard_cols <- c("FID", "IID", "PAT", "MAT", "SEX", "PHENOTYPE")
        variant_cols <- setdiff(names(dt_geno_exported), standard_cols)

        # Get variant IDs from .bim file
        bim_file <- paste0(local_genotype_path, ".bim")
        variant_id_map <- data.table()
        if (file.exists(bim_file)) {
            dt_bim_maf <- fread(bim_file, col.names = c("chr", "rsid", "cm", "pos", "a1", "a2"))
            for (rsid_bim in dt_bim_maf$rsid) {
                matching_cols <- grep(paste0("^", rsid_bim, "_"), variant_cols, value = TRUE)
                if (length(matching_cols) > 0) {
                    variant_id_map <- rbind(variant_id_map, data.table(
                        column_name = matching_cols[1],
                        rsid = rsid_bim
                    ))
                }
            }
        }

        log_info("Calculating MAF for {nrow(variant_id_map)} variants from exported genotypes...")
        maf_results <- list()

        for (i in seq_len(nrow(variant_id_map))) {
            var_col <- variant_id_map[i]$column_name
            rsid <- variant_id_map[i]$rsid

            # Get genotype values (0, 1, 2 = hom ref, het, hom alt)
            geno_values <- dt_geno_exported[[var_col]]
            geno_values <- geno_values[!is.na(geno_values)]

            if (length(geno_values) == 0) next

            # Calculate allele frequency: mean(geno) / 2 gives alternate allele frequency
            alt_freq <- mean(geno_values) / 2
            maf <- pmin(alt_freq, 1 - alt_freq)

            # Handle edge cases
            if (is.na(maf) || maf <= 0) {
                maf <- 0.001  # Minimum MAF for rare variants
            }

            maf_results[[length(maf_results) + 1]] <- data.table(rsid = rsid, maf = maf)
        }

        if (length(maf_results) > 0) {
            maf_data <- rbindlist(maf_results)
            log_info("Calculated MAF for {nrow(maf_data)} variants")
            log_info("MAF range: {round(min(maf_data$maf, na.rm=TRUE), 4)} - {round(max(maf_data$maf, na.rm=TRUE), 4)}")
        } else {
            log_warn("No MAF values calculated. Proceeding without MAF weighting.")
        }
    } else {
        log_warn("Exported genotype data not available. Cannot calculate MAF.")
    }

    # Calculate heterozygosity from exported genotype data
    log_info("=== CALCULATING HETEROZYGOSITY FROM EXPORTED GENOTYPE DATA ===")
    heterozygosity_data <- NULL
    if (!is.null(dt_geno_exported) && nrow(dt_geno_exported) > 0) {
        # Use same variant_id_map from MAF calculation
        log_info("Calculating heterozygosity for {nrow(variant_id_map)} variants (samples in NPX matrix)...")

        # Filter to samples that are in NPX matrix
        geno_samples <- dt_geno_exported$IID
        npx_samples <- rownames(npx_matrix)
        common_samples_het <- intersect(geno_samples, npx_samples)

        if (length(common_samples_het) == 0) {
            log_warn("No common samples between genotype data and NPX matrix for heterozygosity calculation")
        } else {
            log_info("Calculating heterozygosity for {length(common_samples_het)} common samples")

            # Filter genotype data to common samples
            dt_geno_common <- dt_geno_exported[IID %in% common_samples_het]

            het_results <- list()
            for (i in seq_len(nrow(variant_id_map))) {
                var_col <- variant_id_map[i]$column_name
                rsid <- variant_id_map[i]$rsid

                geno_values <- dt_geno_common[[var_col]]
                geno_values <- geno_values[!is.na(geno_values)]

                if (length(geno_values) == 0) next

                # Count heterozygous (genotype == 1)
                n_het <- sum(geno_values == 1, na.rm = TRUE)
                n_total <- length(geno_values)
                het_prop <- if (n_total > 0) n_het / n_total else 0

                het_results[[length(het_results) + 1]] <- data.table(rsid = rsid, heterozygosity = het_prop)
            }

            if (length(het_results) > 0) {
                heterozygosity_data <- rbindlist(het_results)
                log_info("Calculated heterozygosity for {nrow(heterozygosity_data)} variants")
                log_info("Heterozygosity range: {round(min(heterozygosity_data$heterozygosity, na.rm=TRUE), 4)} - {round(max(heterozygosity_data$heterozygosity, na.rm=TRUE), 4)}")
                n_zero_het <- sum(heterozygosity_data$heterozygosity == 0, na.rm = TRUE)
                if (n_zero_het > 0) {
                    log_info("Variants with 0% heterozygosity (all homozygous): {n_zero_het}")
                }
            } else {
                log_warn("No heterozygosity values calculated.")
            }
        }
    } else {
        log_warn("Exported genotype data not available. Cannot calculate heterozygosity.")
    }

    # Merge MAF and heterozygosity into top_variants
    log_info("=== MERGING MAF AND HETEROZYGOSITY INTO COLLATED FILE ===")

    # Save a copy of full top_variants before merging (for updating collated file)
    top_variants_full <- copy(top_variants)

    # Merge MAF
    if (!is.null(maf_data) && nrow(maf_data) > 0) {
        if (!("rsid_for_matching" %in% names(top_variants))) {
            log_warn("rsid_for_matching column missing. Creating from rsid.")
            top_variants[, rsid_for_matching := rsid]
        }
        if ("maf" %in% names(top_variants)) {
            top_variants[, maf := NULL]
        }
        top_variants <- merge(top_variants, maf_data, by.x = "rsid_for_matching", by.y = "rsid", all.x = TRUE)
        n_with_maf <- sum(!is.na(top_variants$maf))
        n_without_maf <- sum(is.na(top_variants$maf))
        if (n_without_maf > 0) {
            log_warn("{n_without_maf} variants have missing MAF. Using default MAF (0.01) for ranking.")
            top_variants[is.na(maf), maf := 0.01]
        }
        log_info("MAF merged: {n_with_maf} variants with MAF, {n_without_maf} variants using default")
    } else {
        top_variants[, maf := 0.01]
        log_warn("MAF calculation failed or unavailable. Using default MAF (0.01) for all variants.")
    }

    # Merge heterozygosity
    if (!is.null(heterozygosity_data) && nrow(heterozygosity_data) > 0) {
        if (!("rsid_for_matching" %in% names(top_variants))) {
            top_variants[, rsid_for_matching := rsid]
        }
        if ("heterozygosity" %in% names(top_variants)) {
            top_variants[, heterozygosity := NULL]
        }
        top_variants <- merge(top_variants, heterozygosity_data, by.x = "rsid_for_matching", by.y = "rsid", all.x = TRUE)
        n_with_het <- sum(!is.na(top_variants$heterozygosity))
        n_without_het <- sum(is.na(top_variants$heterozygosity))
        if (n_without_het > 0) {
            log_warn("{n_without_het} variants have missing heterozygosity. Using default (0.5) for ranking.")
            top_variants[is.na(heterozygosity), heterozygosity := 0.5]
        }
        log_info("Heterozygosity merged: {n_with_het} variants with heterozygosity, {n_without_het} variants using default")
    } else {
        top_variants[, heterozygosity := 0.5]
        log_warn("Heterozygosity calculation failed or unavailable. Using default (0.5) for all variants.")
    }

    # Update top_variants_full with MAF and heterozygosity
    # IMPORTANT: De-duplicate by rsid first to avoid cartesian join
    # Multiple proteins can share the same variant (rsid), but MAF and heterozygosity are per-variant
    # Take the first occurrence of each rsid (they should all have the same MAF/heterozygosity)
    maf_het_by_rsid <- unique(top_variants[, .(rsid, maf, heterozygosity)], by = "rsid")
    log_info("Created unique rsid mapping: {nrow(maf_het_by_rsid)} unique rsids from {nrow(top_variants)} variant rows")

    # Merge with unique values to avoid cartesian product
    top_variants_full <- merge(top_variants_full, maf_het_by_rsid, by = "rsid", all.x = TRUE, suffixes = c("", ".new"))

    # Handle potential duplicate columns from merge
    if ("maf.new" %in% names(top_variants_full)) {
        # Use new values, preferring non-NA
        top_variants_full[, maf := fifelse(is.na(maf.new), maf, maf.new)]
        top_variants_full[, maf.new := NULL]
    }
    if ("heterozygosity.new" %in% names(top_variants_full)) {
        top_variants_full[, heterozygosity := fifelse(is.na(heterozygosity.new), heterozygosity, heterozygosity.new)]
        top_variants_full[, heterozygosity.new := NULL]
    }

    # Apply heterozygosity-weighted ranking
    log_info("=== APPLYING HETEROZYGOSITY-WEIGHTED RANKING ===")

    # Handle zero p-values before ranking
    min_p_value <- .Machine$double.xmin * 10
    zero_p_count <- sum(top_variants$p == 0, na.rm=TRUE)
    if (zero_p_count > 0) {
        log_warn("Found {zero_p_count} variants with p=0. Replacing with minimum p-value ({min_p_value}) to prevent Inf in composite score.")
        top_variants[p == 0, p := min_p_value]
    }
    tiny_p_count <- sum(top_variants$p > 0 & top_variants$p < 1e-300, na.rm=TRUE)
    if (tiny_p_count > 0) {
        log_warn("Found {tiny_p_count} variants with p < 1e-300. Replacing with minimum p-value to prevent numerical issues.")
        top_variants[p > 0 & p < 1e-300, p := min_p_value]
    }

    # Calculate composite score with heterozygosity weighting
    if ("maf" %in% names(top_variants) && "heterozygosity" %in% names(top_variants)) {
        log_info("Applying heterozygosity-weighted ranking: (-log10(p) * abs(beta) * maf) / (heterozygosity + 0.01)")
        top_variants[, heterozygosity_adj := heterozygosity + 0.01]
        top_variants[, composite_score := (-log10(p) * abs(beta) * maf) / heterozygosity_adj]
        log_info("Ranked variants using heterozygosity-weighted composite score")
    } else if ("maf" %in% names(top_variants)) {
        log_info("Applying MAF-weighted ranking: -log10(p) * abs(beta) * maf")
        top_variants[, composite_score := -log10(p) * abs(beta) * maf]
        log_info("Ranked variants using MAF-weighted composite score")
    } else {
        log_info("Applying standard ranking: -log10(p) * abs(beta)")
        top_variants[, composite_score := -log10(p) * abs(beta)]
        log_info("Ranked variants using standard composite score")
    }

    # Sort by composite score (descending)
    setorder(top_variants, -composite_score)

    # =========================================================================
    # Apply MAF filter for selection (if configured)
    # =========================================================================
    min_maf_for_selection <- tryCatch(pqtl_config$min_maf_for_selection, error = function(e) NULL)
    if (!is.null(min_maf_for_selection) && is.numeric(min_maf_for_selection) && min_maf_for_selection > 0) {
        top_variants_backup <- copy(top_variants)
        top_variants <- apply_maf_filter(top_variants, min_maf_for_selection, "pQTL selection")
        if (nrow(top_variants) == 0) {
            log_warn("No variants remain after MAF filter. Restoring all variants.")
            top_variants <- top_variants_backup
        }
    }

    # =========================================================================
    # Apply pQTL selection: Use trained pQTLs if available, else use top_n
    # =========================================================================
    selected_rsids <- character(0)

    if (use_trained_params && !is.null(trained_pqtls) && nrow(trained_pqtls) > 0) {
        # Use consensus/LASSO-selected pQTLs from training step
        log_info("")
        log_info("Using consensus pQTLs from training step (05a)")

        # Extract rsid column (handle both TSV and YAML list formats)
        if ("rsid" %in% names(trained_pqtls)) {
            trained_rsids <- trained_pqtls$rsid
        } else if (ncol(trained_pqtls) == 1) {
            # Single column data.table (from YAML list conversion)
            trained_rsids <- unlist(trained_pqtls)
        } else {
            log_warn("Unexpected format in trained pQTLs. Expected 'rsid' column or single column. Falling back to top_n selection.")
            use_trained_params <- FALSE
            trained_rsids <- character(0)
        }

        if (length(trained_rsids) > 0) {
            # Filter top_variants to only include consensus pQTLs
            # Match by rsid or rsid_for_matching (for sex chromosome variants)
            matching_rsids <- intersect(trained_rsids, top_variants$rsid)
            if (length(matching_rsids) == 0 && "rsid_for_matching" %in% names(top_variants)) {
                # Try matching with rsid_for_matching
                matching_rsids <- intersect(trained_rsids, top_variants$rsid_for_matching)
                if (length(matching_rsids) > 0) {
                    top_variants <- top_variants[rsid_for_matching %in% matching_rsids]
                    selected_rsids <- top_variants$rsid
                }
            } else if (length(matching_rsids) > 0) {
                top_variants <- top_variants[rsid %in% matching_rsids]
                selected_rsids <- matching_rsids
            }

            if (length(matching_rsids) > 0) {
                # IMPORTANT: Deduplicate by rsid to ensure we only use unique consensus variants
                # The same rsid can appear multiple times if it's a pQTL for multiple proteins
                # When using consensus pQTLs, we want to use each variant only once
                n_before_dedup <- nrow(top_variants)
                if ("composite_score" %in% names(top_variants)) {
                    # Keep the variant-protein combination with highest composite score
                    setorder(top_variants, -composite_score)
                    top_variants <- top_variants[, .SD[1], by = rsid]
                } else if ("p" %in% names(top_variants)) {
                    # Fallback: keep the one with lowest p-value
                    setorder(top_variants, p)
                    top_variants <- top_variants[, .SD[1], by = rsid]
                } else {
                    # Last resort: just take first occurrence
                    top_variants <- top_variants[, .SD[1], by = rsid]
                }
                n_after_dedup <- nrow(top_variants)

                # Update selected_rsids to match deduplicated variants
                selected_rsids <- top_variants$rsid

                log_info("Selected {length(unique(selected_rsids))} unique consensus pQTLs (out of {length(trained_rsids)} from training)")
                if (n_before_dedup > n_after_dedup) {
                    log_info("  Deduplicated: {n_before_dedup} variant-protein combinations -> {n_after_dedup} unique variants")
                }
                if (length(matching_rsids) < length(trained_rsids)) {
                    log_warn("  {length(trained_rsids) - length(matching_rsids)} consensus pQTLs not found in current dataset")
                    missing <- setdiff(trained_rsids, matching_rsids)
                    if (length(missing) <= 10) {
                        log_warn("  Missing pQTLs: {paste(missing, collapse=', ')}")
                    } else {
                        log_warn("  Missing pQTLs: {paste(head(missing, 10), collapse=', ')}, ... ({length(missing)-10} more)")
                    }
                }
            } else {
                log_warn("No consensus pQTLs found in current dataset. Falling back to top_n selection.")
                use_trained_params <- FALSE
            }
        } else {
            log_warn("No consensus pQTLs to use. Falling back to top_n selection.")
            use_trained_params <- FALSE
        }
    }

    # Fall back to top_n selection if trained params not available
    if (!use_trained_params || length(selected_rsids) == 0) {
        if (!is.null(top_n) && is.numeric(top_n) && top_n >= 100 && top_n < nrow(top_variants)) {
            log_info("Selecting top {top_n} pQTLs from {nrow(top_variants)} available...")
            top_variants <- top_variants[1:top_n]
            selected_rsids <- top_variants$rsid
            log_info("Selected top {nrow(top_variants)} variants based on composite score ranking")
            if ("heterozygosity" %in% names(top_variants)) {
                log_info("Selected variants heterozygosity range: {round(min(top_variants$heterozygosity, na.rm=TRUE), 4)} - {round(max(top_variants$heterozygosity, na.rm=TRUE), 4)}")
            }
        } else if (!is.null(top_n)) {
            if (top_n < 100) {
                log_warn("top_n_pQTLs ({top_n}) is less than 100. Using all variants instead.")
                selected_rsids <- top_variants$rsid
            } else if (top_n >= nrow(top_variants)) {
                log_info("top_n_pQTLs ({top_n}) is >= total variants ({nrow(top_variants)}). Using all variants.")
                selected_rsids <- top_variants$rsid
            } else {
                selected_rsids <- top_variants$rsid
            }
        } else {
            log_info("Using all {nrow(top_variants)} variants (top_n_pQTLs not specified)")
            selected_rsids <- top_variants$rsid
        }
    }

    # Update included_in_top_n flag
    if (length(selected_rsids) > 0) {
        top_variants_full[, included_in_top_n := rsid %in% selected_rsids]
    } else {
        top_variants_full[, included_in_top_n := TRUE]
    }

    # Update composite_score in top_variants_full (for variants that were selected)
    # De-duplicate by rsid to avoid cartesian join (multiple proteins can share same variant)
    composite_by_rsid <- unique(top_variants[, .(rsid, composite_score)], by = "rsid")
    top_variants_full <- merge(top_variants_full, composite_by_rsid, by = "rsid", all.x = TRUE, suffixes = c("", ".new"))

    # Handle potential duplicate columns from merge
    if ("composite_score.new" %in% names(top_variants_full)) {
        top_variants_full[, composite_score := fifelse(is.na(composite_score.new), composite_score, composite_score.new)]
        top_variants_full[, composite_score.new := NULL]
    }

    # Update collated file with MAF, heterozygosity, composite_score, and inclusion flag
    output_cols <- c("trait", "rsid", "p", "beta", "sd", "cs_avg_r2")
    if ("maf" %in% names(top_variants_full)) {
        output_cols <- c(output_cols, "maf")
    }
    if ("heterozygosity" %in% names(top_variants_full)) {
        output_cols <- c(output_cols, "heterozygosity")
    }
    if ("composite_score" %in% names(top_variants_full)) {
        output_cols <- c(output_cols, "composite_score")
    }
    if ("included_in_top_n" %in% names(top_variants_full)) {
        output_cols <- c(output_cols, "included_in_top_n")
    }

    # Ensure all output columns exist
    for (col in output_cols) {
        if (!(col %in% names(top_variants_full))) {
            top_variants_full[[col]] <- NA
        }
    }

    # Select only output columns and save
    top_variants_output <- top_variants_full[, ..output_cols]
    fwrite(top_variants_output, collated_path, sep = "\t")
    log_info("Updated collated fine-mapping file with MAF, heterozygosity, composite_score, and top_n inclusion flag: {collated_path}")
    log_info("  Variants with MAF: {sum(!is.na(top_variants_output$maf))}")
    if ("heterozygosity" %in% names(top_variants_output)) {
        log_info("  Variants with heterozygosity: {sum(!is.na(top_variants_output$heterozygosity))}")
        log_info("  Heterozygosity range: {round(min(top_variants_output$heterozygosity, na.rm=TRUE), 4)} - {round(max(top_variants_output$heterozygosity, na.rm=TRUE), 4)}")
    }
    if ("composite_score" %in% names(top_variants_output)) {
        log_info("  Composite score range: {round(min(top_variants_output$composite_score, na.rm=TRUE), 2)} - {round(max(top_variants_output$composite_score, na.rm=TRUE), 2)}")
    }
    log_info("  Variants included in top_n: {sum(top_variants_output$included_in_top_n, na.rm=TRUE)}")

    # Generate output file with variants used for z-score calculation (ordered by composite score)
    log_info("=== GENERATING Z-SCORE VARIANT LIST ===")
    if (length(selected_rsids) > 0) {
        zscore_variants <- top_variants[, .(trait, rsid, p, beta, sd, cs_avg_r2, maf, heterozygosity, composite_score)]
        # Sort by composite score (descending)
        setorder(zscore_variants, -composite_score)

        zscore_variants_path <- get_output_path(step_num, "05b_zscore_variants", batch_id, "pqtl", "tsv", config = config)
        ensure_output_dir(zscore_variants_path)
        fwrite(zscore_variants, zscore_variants_path, sep = "\t")
        log_info("Saved z-score variant list (ordered by composite score) to: {zscore_variants_path}")
        log_info("  Variants: {nrow(zscore_variants)}")
        log_info("  Composite score range: {round(min(zscore_variants$composite_score, na.rm=TRUE), 2)} - {round(max(zscore_variants$composite_score, na.rm=TRUE), 2)}")
    }

    # Store exported genotypes for z-score calculation (will be reused later)
    # dt_geno_exported is already loaded and cached
    log_info("=== UNIFIED GENOTYPE EXTRACTION COMPLETE ===")
    log_info("Exported genotypes cached and ready for z-score calculation")

    # 3. Prepare Genotype Data for Z-score Calculation
    # -------------------------------------------------
    # REFACTORED: Reuse cached exported genotypes (dt_geno_exported) from unified extraction
    log_info("=== PREPARING GENOTYPE DATA FOR Z-SCORE CALCULATION ===")

    # dt_geno_exported is already loaded from unified extraction above
    if (is.null(dt_geno_exported) || nrow(dt_geno_exported) == 0) {
        stop("Exported genotype data not available. Cannot proceed with z-score calculation.")
    }

    log_info("Using cached exported genotypes: {nrow(dt_geno_exported)} samples")

    # Filter to selected variants (top_n) if needed
    # Get variant columns from exported data
    standard_cols <- c("FID", "IID", "PAT", "MAT", "SEX", "PHENOTYPE")
    all_variant_cols <- setdiff(names(dt_geno_exported), standard_cols)

    # Get selected variant rsids (already filtered to top_n above)
    selected_rsids_for_matching <- top_variants$rsid_for_matching

    # Get variant columns for selected variants
    selected_variant_cols <- character(0)
    for (rsid_match in selected_rsids_for_matching) {
        matching_cols <- grep(paste0("^", rsid_match, "_"), all_variant_cols, value = TRUE)
        if (length(matching_cols) > 0) {
            selected_variant_cols <- c(selected_variant_cols, matching_cols[1])  # Use first match
        }
    }

    log_info("Filtering to {length(selected_variant_cols)} selected variant columns (from {length(all_variant_cols)} total)")

    # Keep standard columns + selected variant columns
    cols_to_keep <- c(standard_cols, selected_variant_cols)
    cols_to_keep <- intersect(cols_to_keep, names(dt_geno_exported))
    dt_geno <- dt_geno_exported[, ..cols_to_keep]

    log_info("Filtered genotype data to {length(selected_variant_cols)} variant columns")

    # Extract variant IDs and filter to successfully extracted variants
    geno_rsids <- extract_variant_ids(dt_geno, local_genotype_path)
    top_variants_extracted <- top_variants[rsid_for_matching %in% geno_rsids]
    log_info("Found {length(geno_rsids)} extracted variants, {nrow(top_variants_extracted)} match selected variants")

    # Apply MAF filter for z-score calculation (MAF > 5%)
    maf_threshold_zscore <- 0.05
    top_variants_extracted <- apply_maf_filter(top_variants_extracted, maf_threshold_zscore, "z-score calculation")

    if (nrow(top_variants_extracted) == 0) {
        log_warn("No variants remaining after filtering. Z-score calculation will be skipped.")
        return(NULL)
    }

    # 4. Calculate Z-scores
    # ---------------------
    log_info("Calculating Z-scores for {nrow(top_variants_extracted)} variants...")

    # Get available proteins in NPX matrix
    available_proteins <- colnames(npx_matrix)
    log_info("NPX matrix contains {length(available_proteins)} proteins")

    results <- list()
    failed_variants <- data.table()

    for (i in 1:nrow(top_variants_extracted)) {
        prot <- top_variants_extracted[i]$trait
        # Use rsid_for_matching for finding genotype column (may have chrX conversion)
        rsid_match <- top_variants_extracted[i]$rsid_for_matching
        # Use rsid_original for storing results
        rsid_orig <- top_variants_extracted[i]$rsid_original

        # Check if protein exists in NPX matrix (fixes "subscript out of bounds" error)
        if (!prot %in% available_proteins) {
            log_warn("Protein '{prot}' not found in NPX matrix for variant {rsid_orig}. Skipping.")
            failed_variants <- rbind(failed_variants, data.table(
                trait = prot,
                rsid = rsid_orig,
                reason = "Protein not found in NPX matrix"
            ))
            next
        }

        # Find column in geno (PLINK export A format: rsid_A)
        # Try matching with rsid_for_matching first, then try chrX/chrx variants
        geno_col <- grep(paste0("^", rsid_match), names(dt_geno), value = TRUE)
        if (length(geno_col) == 0 && sex_chr_enabled && grepl("^chrX_", rsid_match)) {
            # Try chrx (lowercase) as fallback
            rsid_chrx <- gsub("^chrX_", "chrx_", rsid_match)
            geno_col <- grep(paste0("^", rsid_chrx), names(dt_geno), value = TRUE)
        }
        if (length(geno_col) == 0) {
            log_warn("Variant {rsid_orig} (matching ID: {rsid_match}) not found in genotype columns. Skipping.")
            failed_variants <- rbind(failed_variants, data.table(
                trait = prot,
                rsid = rsid_orig,
                reason = "Variant not found in genotype columns"
            ))
            next
        }

        # Get data
        # Match samples
        common_samples <- intersect(dt_geno$IID, rownames(npx_matrix))
        if (length(common_samples) == 0) {
            log_warn("No common samples for variant {rsid_orig}. Skipping.")
            failed_variants <- rbind(failed_variants, data.table(
                trait = prot,
                rsid = rsid_orig,
                reason = "No common samples"
            ))
            next
        }

        # Extract protein values (now safe since we checked prot exists)
        protein_values <- npx_matrix[common_samples, prot, drop = FALSE]

        df <- data.table(
            SampleID = common_samples,
            Geno = dt_geno[match(common_samples, IID)][[geno_col[1]]],
            Protein = as.numeric(protein_values[, 1])
        )

        # Flip genotype to match reference implementation:
        # PLINK --export A: 0 = hom ref, 1 = het, 2 = hom alt
        # After flip: 2 = hom ref, 1 = het, 0 = hom alt
        # This aligns with beta interpretation (per-alternate-allele effect)
        df[, Geno := 2 - Geno]

        df <- df[!is.na(Protein) & !is.na(Geno)]

        # Calculate stats per genotype
        stats <- df[, .(Mean = mean(Protein), SD = sd(Protein)), by = Geno]

        # Check for valid SD (non-zero) before calculating Z
        stats <- stats[!is.na(SD) & SD > 0]
        if (nrow(stats) == 0) {
            log_warn("No valid genotype groups with non-zero SD for variant {rsid_orig}. Skipping.")
            failed_variants <- rbind(failed_variants, data.table(
                trait = prot,
                rsid = rsid_orig,
                reason = "No valid genotype groups"
            ))
            next
        }

        # Calc Z
        df <- merge(df, stats, by = "Geno", all.x = FALSE)
        if (nrow(df) == 0) {
            log_warn("No matching genotype groups after merge for variant {rsid_orig}. Skipping.")
            failed_variants <- rbind(failed_variants, data.table(
                trait = prot,
                rsid = rsid_orig,
                reason = "No matching genotype groups"
            ))
            next
        }

        df[, Z := (Protein - Mean) / SD]

        # Calculate beta-weighted residual for MAR metric
        beta_val <- top_variants_extracted[i]$beta
        if (!is.na(beta_val) && is.finite(beta_val) && nrow(df) > 3) {
            pop_mean_protein <- mean(df$Protein, na.rm = TRUE)
            mean_geno <- mean(df$Geno, na.rm = TRUE)

            # Predicted protein = population_mean + beta * (genotype - mean_genotype)
            df[, Predicted := pop_mean_protein + beta_val * (Geno - mean_geno)]
            df[, Residual := Protein - Predicted]

            # Standardize residual
            resid_sd <- sd(df$Residual, na.rm = TRUE)
            if (!is.na(resid_sd) && resid_sd > 0) {
                df[, StdResidual := Residual / resid_sd]
            } else {
                df[, StdResidual := NA_real_]
            }
        } else {
            df[, StdResidual := NA_real_]
        }

        results[[i]] <- df[, .(SampleID, Protein, rsid = rsid_orig, Z, StdResidual)]
    }

    # Save failed variants log
    if (nrow(failed_variants) > 0) {
        failed_variants_log <- get_output_path(step_num, "05b_failed_zscore_variants", batch_id, "pqtl", "tsv", config = config)
        ensure_output_dir(failed_variants_log)
        fwrite(failed_variants, failed_variants_log, sep = "	")
        log_warn("Failed to calculate Z-scores for {nrow(failed_variants)} variants. Log saved to: {failed_variants_log}")
    }

    # Remove NULL entries from results
    results <- results[!sapply(results, is.null)]

    if (length(results) == 0) {
        log_error("No Z-scores calculated. All variants failed.")
        return(NULL)
    }

    dt_z <- rbindlist(results)
    log_info("Successfully calculated Z-scores for {length(unique(dt_z$rsid))} variants")

    # 5. Identify Outliers
    # --------------------
    log_info("Identifying outliers...")

    if (nrow(dt_z) == 0) {
        log_warn("No Z-scores calculated. Check overlap between genotypes and phenotypes.")
        return(NULL)
    }

    sample_stats <- dt_z[, .(
        MeanAbsZ = mean(abs(Z), na.rm = TRUE),
        MedianAbsZ = median(abs(Z), na.rm = TRUE),
        MaxAbsZ = max(abs(Z), na.rm = TRUE),
        MedianAbsResidual = median(abs(StdResidual), na.rm = TRUE),
        N_Prots = .N,
        N_ValidResiduals = sum(!is.na(StdResidual))
    ), by = SampleID]

    # Define thresholds
    # Calculate population statistics first (needed for k_sd/k_mad conversion)
    pop_mean <- mean(sample_stats$MeanAbsZ, na.rm = TRUE)
    pop_sd <- sd(sample_stats$MeanAbsZ, na.rm = TRUE)
    pop_median_medabsz <- median(sample_stats$MedianAbsZ, na.rm = TRUE)
    pop_mad_medabsz <- median(abs(sample_stats$MedianAbsZ - pop_median_medabsz), na.rm = TRUE)
    mar_median <- median(sample_stats$MedianAbsResidual, na.rm = TRUE)
    mar_mad <- median(abs(sample_stats$MedianAbsResidual - mar_median), na.rm = TRUE)

    # MeanAbsZ threshold
    if (use_trained_params && !is.null(trained_thresholds)) {

        # Handle consensus config format (k_sd/k_mad) vs legacy format (absolute threshold)
        if (!is.null(trained_thresholds$mean_abs_z)) {
            if (!is.null(trained_thresholds$mean_abs_z$k_sd) && !is.na(trained_thresholds$mean_abs_z$k_sd)) {
                # Convert k_sd to absolute threshold
                k_sd <- trained_thresholds$mean_abs_z$k_sd
                cutoff_meanabsz <- pop_mean + k_sd * pop_sd
                log_info("Using consensus threshold (k_sd format):")
                log_info("  MeanAbsZ: {round(cutoff_meanabsz, 4)} (mean + {round(k_sd, 4)}*SD)")
            } else if (!is.null(trained_thresholds$mean_abs_z$threshold)) {
                # Use absolute threshold directly
                cutoff_meanabsz <- trained_thresholds$mean_abs_z$threshold
                log_info("Using trained threshold (absolute format):")
                log_info("  MeanAbsZ: {round(cutoff_meanabsz, 4)}")
            } else {
                # Fallback to default
                cutoff_meanabsz <- pop_mean + z_threshold * pop_sd
                log_warn("Trained MeanAbsZ threshold not found. Using default: {round(cutoff_meanabsz, 4)}")
            }
        } else {
            cutoff_meanabsz <- pop_mean + z_threshold * pop_sd
            log_warn("Trained MeanAbsZ threshold structure not found. Using default: {round(cutoff_meanabsz, 4)}")
        }

        # Handle MAR threshold
        if (!is.null(trained_thresholds$mar)) {
            if (!is.null(trained_thresholds$mar$k_mad) && !is.na(trained_thresholds$mar$k_mad)) {
                # Convert k_mad to absolute threshold
                k_mad <- trained_thresholds$mar$k_mad
                cutoff_mar <- mar_median + k_mad * mar_mad * 1.4826
                log_info("  MAR: {round(cutoff_mar, 4)} (median + {round(k_mad, 4)}*MAD*1.4826)")
            } else if (!is.null(trained_thresholds$mar$threshold)) {
                # Use absolute threshold directly
                cutoff_mar <- trained_thresholds$mar$threshold
                log_info("  MAR: {round(cutoff_mar, 4)}")
            } else {
                # Fallback to default
                cutoff_mar <- mar_median + z_threshold * mar_mad * 1.4826
                log_warn("Trained MAR threshold not found. Using default: {round(cutoff_mar, 4)}")
            }
        } else if (!is.null(trained_thresholds$median_abs_residual) && !is.null(trained_thresholds$median_abs_residual$threshold)) {
            # Legacy format
            cutoff_mar <- trained_thresholds$median_abs_residual$threshold
            log_info("  MAR: {round(cutoff_mar, 4)} (legacy format)")
        } else {
            cutoff_mar <- mar_median + z_threshold * mar_mad * 1.4826
            log_warn("Trained MAR threshold structure not found. Using default: {round(cutoff_mar, 4)}")
        }
    } else {
        # Fall back to population-based threshold
        cutoff_meanabsz <- pop_mean + z_threshold * pop_sd
        cutoff_mar <- mar_median + z_threshold * mar_mad * 1.4826

        log_info("Using population-based thresholds:")
        log_info("  MeanAbsZ: {round(cutoff_meanabsz, 4)} (mean + {z_threshold}*SD)")
        log_info("  MAR: {round(cutoff_mar, 4)} (median + {z_threshold}*MAD*1.4826)")
    }

    # Calculate threshold for MedianAbsZ (MAD-based, similar to MAR)
    cutoff_medianabsz <- pop_median_medabsz + z_threshold * pop_mad_medabsz * 1.4826

    # Round all thresholds to 1 decimal place
    cutoff_meanabsz <- round(cutoff_meanabsz, digits = 1)
    cutoff_medianabsz <- round(cutoff_medianabsz, digits = 1)
    cutoff_mar <- round(cutoff_mar, digits = 1)

    log_info("Thresholds (rounded to 1 decimal place):")
    log_info("  MeanAbsZ: {cutoff_meanabsz}")
    log_info("  MedianAbsZ: {cutoff_medianabsz}")
    log_info("  MAR: {cutoff_mar}")

    # Flag outliers by each metric (using rounded thresholds)
    sample_stats[, Outlier_MeanAbsZ := MeanAbsZ > cutoff_meanabsz]
    sample_stats[, Outlier_MedianAbsZ := MedianAbsZ > cutoff_medianabsz]
    sample_stats[, Outlier_MAR := MedianAbsResidual > cutoff_mar & N_ValidResiduals >= 10]

    # Count outliers for each metric
    n_outliers_meanabsz <- sum(sample_stats$Outlier_MeanAbsZ, na.rm = TRUE)
    n_outliers_medianabsz <- sum(sample_stats$Outlier_MedianAbsZ, na.rm = TRUE)
    n_outliers_mar <- sum(sample_stats$Outlier_MAR, na.rm = TRUE)

    log_info("Outliers detected:")
    log_info("  By MeanAbsZ > {cutoff_meanabsz}: {n_outliers_meanabsz}")
    log_info("  By MedianAbsZ > {cutoff_medianabsz}: {n_outliers_medianabsz}")
    log_info("  By MAR > {cutoff_mar}: {n_outliers_mar}")
    log_info("Using MeanAbsZ for final outlier assignment: {n_outliers_meanabsz} outliers")

    # 6. Save Outputs
    # ---------------
    # Add FINNGENID column to sample_stats
    # Note: SampleID column currently contains FINNGENIDs (from rownames conversion)
    # We need to preserve original SampleID and add FINNGENID
    log_info("Adding FINNGENID mapping to output tables...")
    sample_stats <- add_finngenid_column(sample_stats, batch_id = batch_id, config = config,
                                         sample_id_col = "SampleID", preserve_original = TRUE)

    out_file <- get_output_path(step_num, "05b_pqtl_outliers", batch_id, "outliers", "tsv", config = config)
    ensure_output_dir(out_file)
    # Only include samples flagged by MeanAbsZ (not MAR)
    fwrite(sample_stats[Outlier_MeanAbsZ == TRUE], out_file, sep = "\t")

    # Save full stats
    stats_file <- get_output_path(step_num, "05b_pqtl_stats", batch_id, "outliers", "tsv", config = config)
    fwrite(sample_stats, stats_file, sep = "\t")

    # Also save RDS version with FINNGENID
    stats_file_rds <- get_output_path(step_num, "05b_pqtl_stats", batch_id, "outliers", "rds", config = config)
    saveRDS(sample_stats, stats_file_rds)

    # 7. Visualization
    # ----------------
    log_info("Generating plots...")

    # Plot 1: Mean Absolute Z-score
    p1 <- ggplot(sample_stats, aes(x = MeanAbsZ, fill = !Outlier_MeanAbsZ)) +
        geom_histogram(bins = 50) +
        geom_vline(xintercept = cutoff_meanabsz, linetype = "dashed", color = "red", linewidth = 1) +
        scale_fill_manual(
            values = c("TRUE" = "steelblue", "FALSE" = "firebrick"),
            labels = c("TRUE" = "Non-outlier", "FALSE" = "Outlier"),
            name = "Sample Status"
        ) +
        labs(
            title = "Distribution of Mean Absolute Z-scores",
            subtitle = paste("Threshold:", cutoff_meanabsz, "| Outliers detected:", n_outliers_meanabsz),
            x = "Mean |Z|", y = "Count"
        ) +
        theme_bw() +
        theme(plot.title = element_text(size = 14, face = "bold"),
              plot.subtitle = element_text(size = 11))

    plot_file <- get_output_path(step_num, "05b_mean_abs_z_dist", batch_id, "outliers", "pdf", config = config)
    ensure_output_dir(plot_file)
    ggsave(plot_file, p1, width = 8, height = 6)
    log_info("Saved Mean Absolute Z-score plot to: {plot_file}")

    # Plot 2: Median Absolute Z-score
    p2 <- ggplot(sample_stats, aes(x = MedianAbsZ, fill = !Outlier_MedianAbsZ)) +
        geom_histogram(bins = 50) +
        geom_vline(xintercept = cutoff_medianabsz, linetype = "dashed", color = "red", linewidth = 1) +
        scale_fill_manual(
            values = c("TRUE" = "steelblue", "FALSE" = "firebrick"),
            labels = c("TRUE" = "Non-outlier", "FALSE" = "Outlier"),
            name = "Sample Status"
        ) +
        labs(
            title = "Distribution of Median Absolute Z-scores",
            subtitle = paste("Threshold:", cutoff_medianabsz, "| Samples above threshold:", n_outliers_medianabsz),
            x = "Median |Z|", y = "Count"
        ) +
        theme_bw() +
        theme(plot.title = element_text(size = 14, face = "bold"),
              plot.subtitle = element_text(size = 11))

    plot_file2 <- get_output_path(step_num, "05b_median_abs_z_dist", batch_id, "outliers", "pdf", config = config)
    ensure_output_dir(plot_file2)
    ggsave(plot_file2, p2, width = 8, height = 6)
    log_info("Saved Median Absolute Z-score plot to: {plot_file2}")

    # Plot 3: Median Absolute Residual
    p3 <- ggplot(sample_stats, aes(x = MedianAbsResidual, fill = !Outlier_MAR)) +
        geom_histogram(bins = 50) +
        geom_vline(xintercept = cutoff_mar, linetype = "dashed", color = "red", linewidth = 1) +
        scale_fill_manual(
            values = c("TRUE" = "steelblue", "FALSE" = "firebrick"),
            labels = c("TRUE" = "Non-outlier", "FALSE" = "Outlier"),
            name = "Sample Status"
        ) +
        labs(
            title = "Distribution of Median Absolute Residuals",
            subtitle = paste("Threshold:", cutoff_mar, "| Outliers detected:", n_outliers_mar),
            x = "Median |Standardized Residual|", y = "Count"
        ) +
        theme_bw() +
        theme(plot.title = element_text(size = 14, face = "bold"),
              plot.subtitle = element_text(size = 11))

    plot_file3 <- get_output_path(step_num, "05b_median_abs_residual_dist", batch_id, "outliers", "pdf", config = config)
    ensure_output_dir(plot_file3)
    ggsave(plot_file3, p3, width = 8, height = 6)
    log_info("Saved Median Absolute Residual plot to: {plot_file3}")

    # 7.5. Cross-Reference Plot: pQTL Z-scores vs Sex Predictions
    # ------------------------------------------------------------
    log_info("Creating cross-reference plot: pQTL Z-scores vs Sex Predictions...")

    # Load sex predictions from step 04
    sex_predictions_path <- get_output_path("04", "sex_predictions", batch_id, "outliers", "tsv", config = config)

    if (file.exists(sex_predictions_path)) {
        sex_preds <- fread(sex_predictions_path)
        log_info("Loaded sex predictions: {nrow(sex_preds)} samples")

        # Normalize column names: sex_preds uses SAMPLE_ID, but we need SampleID for consistency
        if ("SAMPLE_ID" %in% names(sex_preds) && !("SampleID" %in% names(sex_preds))) {
            setnames(sex_preds, "SAMPLE_ID", "SampleID")
        }

        # Extract Youden J threshold from step 04 log file
        # Look for format: "Threshold-based mismatches (Youden J = 0.674): 19" or "Mismatches (Youden J = 0.674): 19"
        youden_threshold <- 0.5  # Default fallback

        # Construct log file paths to try (failure-proof for both single and multi-batch modes)
        # Primary: Use get_log_path() helper (works when run through pipeline runner)
        log_file_04_primary <- get_log_path("04", batch_id, config)

        # Get log directory for fallback searches
        base_dir <- config$output$base_dir %||% Sys.getenv("PIPELINE_OUTPUT_DIR", "output")
        log_dir <- file.path(base_dir, config$output$logs_dir %||% "logs")
        multi_batch_mode <- tryCatch(isTRUE(config$parameters$normalization$multi_batch_mode), error = function(e) FALSE)
        if (multi_batch_mode && !is.null(batch_id)) {
            log_dir <- file.path(log_dir, batch_id)
        }

        # Fallback log file names (for direct script execution in single-batch mode)
        log_file_04_fallback1 <- file.path(log_dir, "04_sex_outliers.log")
        log_file_04_fallback2 <- file.path(log_dir, "04_script.log")
        log_file_04_fallback3 <- file.path(log_dir, "04_pipeline.log")

        # Try all possible log file paths
        log_file_04 <- NULL
        log_file_tried <- character()

        for (candidate_file in c(log_file_04_primary, log_file_04_fallback1, log_file_04_fallback2, log_file_04_fallback3)) {
            log_file_tried <- c(log_file_tried, candidate_file)
            if (file.exists(candidate_file)) {
                log_file_04 <- candidate_file
                log_info("Found step 04 log file: {log_file_04}")
                break
            }
        }

        # If still not found, search for any 04_*.log file in log directory
        if (is.null(log_file_04) && dir.exists(log_dir)) {
            log_files_04 <- list.files(log_dir, pattern = "^04_.*\\.log$", full.names = TRUE)
            if (length(log_files_04) > 0) {
                log_file_04 <- log_files_04[1]
                log_info("Found step 04 log file via pattern search: {log_file_04}")
            }
        }

        # Extract Youden J threshold from log file
        if (!is.null(log_file_04) && file.exists(log_file_04)) {
            log_lines <- readLines(log_file_04)
            # Look specifically for the "Threshold-based mismatches (Youden J = X.XXX)" or "Mismatches (Youden J = X.XXX)" pattern
            youden_line <- grep("(Threshold-based mismatches|Mismatches).*Youden J", log_lines, value = TRUE, ignore.case = TRUE)
            if (length(youden_line) == 0) {
                # Fallback: look for any line with "Youden J"
                youden_line <- grep("Youden J", log_lines, value = TRUE, ignore.case = TRUE)
            }
            if (length(youden_line) > 0) {
                # Extract threshold value - handle formats like:
                # "Mismatches (Youden J = 0.674): 19"
                # "Youden J = 0.674"
                # Try pattern: (Youden J = X.XXX) or Youden J = X.XXX
                youden_pattern <- "\\(Youden J\\s*=\\s*([0-9]+\\.[0-9]+)\\)|Youden J\\s*=\\s*([0-9]+\\.[0-9]+)"
                youden_match <- regmatches(youden_line[1], regexpr(youden_pattern, youden_line[1], ignore.case = TRUE, perl = TRUE))
                if (length(youden_match) > 0) {
                    # Extract the number part
                    num_match <- regmatches(youden_match[1], regexpr("[0-9]+\\.[0-9]+", youden_match[1]))
                    if (length(num_match) > 0) {
                        youden_val <- as.numeric(num_match[1])
                        if (!is.na(youden_val) && youden_val > 0 && youden_val <= 1) {
                            youden_threshold <- youden_val
                            log_info("Extracted Youden J threshold from step 04 log: {youden_threshold}")
                        } else {
                            log_warn("Extracted Youden J value ({youden_val}) is invalid. Using default: {youden_threshold}")
                        }
                    }
                } else {
                    log_warn("Could not extract Youden J threshold from log line: {youden_line[1]}. Using default: {youden_threshold}")
                }
            } else {
                log_warn("No Youden J threshold found in step 04 log file. Using default: {youden_threshold}")
            }
        } else {
            log_warn("Step 04 log file not found. Tried: {paste(log_file_tried, collapse=', ')}. Using default Youden J threshold: {youden_threshold}")
        }

        # Merge pQTL stats with sex predictions (using helper function)
        if (!("SampleID" %in% names(sample_stats)) && "FINNGENID" %in% names(sample_stats)) {
            sample_stats[, SampleID := FINNGENID]
        }

        plot_data <- merge_sex_predictions(sample_stats, sex_preds)

        log_info("Merged data for plotting: {nrow(plot_data)} total samples")
        log_info("  Samples with pQTL stats: {sum(!is.na(plot_data$MeanAbsZ))}")
        log_info("  Samples with sex predictions: {sum(!is.na(plot_data$predicted_prob))}")
        log_info("  Samples with both: {sum(!is.na(plot_data$predicted_prob) & !is.na(plot_data$MeanAbsZ))}")

        # Keep all samples, but handle missing values for plotting
        # Set default values for missing data
        plot_data[is.na(predicted_prob), predicted_prob := 0.5]  # Default to middle if missing
        plot_data[is.na(MeanAbsZ), MeanAbsZ := 0]  # Default to 0 if missing
        plot_data[is.na(genetic_sex), genetic_sex := "unknown"]
        plot_data[is.na(mismatch), mismatch := FALSE]
        plot_data[is.na(sex_outlier), sex_outlier := FALSE]

        if (nrow(plot_data) > 0) {
            # Normalize genetic_sex to lowercase for consistent matching
            plot_data[, genetic_sex := tolower(as.character(genetic_sex))]
            plot_data[!genetic_sex %in% c("male", "female"), genetic_sex := "unknown"]

            # Create color mapping based on genetic sex
            plot_data[, sex_color := ifelse(genetic_sex == "male", "male", ifelse(genetic_sex == "female", "female", "unknown"))]

            # Determine shape: triangle for samples above MeanAbsZ threshold, circle otherwise
            plot_data[, point_shape := ifelse(MeanAbsZ > cutoff_meanabsz, "Above threshold", "Below threshold")]

            # Add border indicator for sex outliers
            plot_data[, has_border := ifelse(mismatch == TRUE | sex_outlier == TRUE, TRUE, FALSE)]

            # NEW: Color gradient based on predicted_prob outlier distance (not MeanAbsZ)
            # Calculate distance from expected predicted_prob for each sex
            # Males should have low predicted_prob (< youden_threshold), females should have high (> youden_threshold)
            # Outliers are those far from expected values
            # For males: outliers have high predicted_prob (closer to 1 = more outlier)
            # For females: outliers have low predicted_prob (closer to 0 = more outlier)
            plot_data[, prob_outlier_dist := ifelse(
                genetic_sex == "male",
                # For males: distance from 0 (expected low prob) - higher prob = more outlier
                # Use predicted_prob directly (0 = normal, 1 = extreme outlier)
                predicted_prob,
                ifelse(
                    genetic_sex == "female",
                    # For females: distance from 1 (expected high prob) - lower prob = more outlier
                    # Use 1 - predicted_prob (0 = normal, 1 = extreme outlier)
                    1 - predicted_prob,
                    0  # Unknown sex: no outlier distance
                )
            )]

            # Normalize outlier distance for color mapping (0 = normal, 1 = extreme outlier)
            valid_dist <- plot_data[!is.na(prob_outlier_dist) & prob_outlier_dist > 0]$prob_outlier_dist
            if (length(valid_dist) > 0) {
                dist_min <- min(valid_dist, na.rm = TRUE)
                dist_max <- max(valid_dist, na.rm = TRUE)
                dist_range <- max(dist_max - dist_min, 0.001)
            } else {
                dist_min <- 0
                dist_max <- 1
                dist_range <- 1
            }

            plot_data[, prob_norm := ifelse(!is.na(prob_outlier_dist) & prob_outlier_dist > 0,
                                           (prob_outlier_dist - dist_min) / dist_range,
                                           0)]

            # Create color palette functions
            # Order: light (normal) → dark (outlier)
            # For males: blue to green gradient (darker shades for sex outliers)
            # Light blue (normal) → medium blue → teal → green → dark green (outliers)
            male_pal <- colorRampPalette(c("#87CEEB", "#5F9EA0", "#4682B4", "#1E90FF", "#00CED1", "#20B2AA", "#2E8B57", "#228B22", "#006400"))
            # For females: orange to pink gradient (darker shades for sex outliers)
            # Light orange (normal) → orange → coral → pink → dark pink/magenta (outliers)
            female_pal <- colorRampPalette(c("#FFE4B5", "#FFA500", "#FF8C00", "#FF6347", "#FF69B4", "#FF1493", "#DC143C", "#C71585", "#8B008B"))
            gray_pal <- colorRampPalette(c("#808080", "#A0A0A0", "#C0C0C0"))  # For unknown sex

            # Assign colors based on sex and normalized predicted_prob outlier distance
            # prob_norm = 0 (normal) → light color (index 1)
            # prob_norm = 1 (extreme outlier) → dark color (index 100)
            # Higher prob_norm (more outlier) = darker shade
            plot_data[, point_color := ifelse(
                genetic_sex == "male",
                male_pal(100)[pmin(100, pmax(1, round(prob_norm * 99 + 1)))],
                ifelse(
                    genetic_sex == "female",
                    female_pal(100)[pmin(100, pmax(1, round(prob_norm * 99 + 1)))],
                    gray_pal(10)[5]  # Gray for unknown
                )
            )]

            # Create the plot
            # Add border color column for mapping
            plot_data[, border_col := ifelse(has_border, "black", "transparent")]

            p_cross <- ggplot(plot_data, aes(x = predicted_prob, y = MeanAbsZ)) +
                # All points with color gradient (using fill for filled shapes)
                # Use shape 21 (filled circle) and 24 (filled triangle) which support fill aesthetic
                geom_point(aes(fill = point_color, shape = point_shape, color = border_col),
                          size = 2.5, alpha = 0.7, stroke = 0.5) +
                # Overlay thicker black borders for sex outliers
                geom_point(data = plot_data[has_border == TRUE],
                          aes(shape = point_shape),
                          fill = NA, color = "black",
                          size = 3, stroke = 1.5, alpha = 1) +
                # Use identity scale for fill (colors already assigned) - no legend
                scale_fill_identity(guide = "none") +
                # Use identity scale for border color (transparent for normal, black for outliers)
                scale_color_identity(guide = "none") +
                # Shape scale: triangle for above threshold, circle for below
                scale_shape_manual(
                    values = c("Above threshold" = 24, "Below threshold" = 21),
                    name = "MeanAbsZ Threshold",
                    labels = c("Above threshold" = paste("Above", cutoff_meanabsz),
                               "Below threshold" = paste("Below", cutoff_meanabsz))
                ) +
                # Vertical line for Youden J threshold (pale red)
                geom_vline(xintercept = youden_threshold, linetype = "dashed", color = "red", linewidth = 0.8, alpha = 0.5) +
                # Horizontal line for MeanAbsZ threshold (pale blue)
                geom_hline(yintercept = cutoff_meanabsz, linetype = "dashed", color = "blue", linewidth = 0.8, alpha = 0.5) +
                # Labels and theme
                labs(
                    title = "pQTL Z-scores vs Sex Predictions",
                    subtitle = paste("MeanAbsZ threshold:", cutoff_meanabsz, "| Youden J threshold:", round(youden_threshold, 3)),
                    x = "Predicted Female Probability",
                    y = "Mean Absolute Z-score",
                    shape = "Z-score Threshold"
                ) +
                theme_bw() +
                theme(
                    plot.title = element_text(size = 14, face = "bold"),
                    plot.subtitle = element_text(size = 11),
                    legend.position = "right"
                ) +
                # Annotate thresholds
                annotate("text", x = youden_threshold, y = max(plot_data$MeanAbsZ, na.rm = TRUE) * 0.95,
                        label = paste("Youden J =", round(youden_threshold, 3)),
                        hjust = -0.1, vjust = 1, color = "red", size = 3) +
                annotate("text", x = max(plot_data$predicted_prob, na.rm = TRUE) * 0.95, y = cutoff_meanabsz,
                        label = paste("MeanAbsZ threshold =", cutoff_meanabsz),
                        hjust = 1, vjust = -0.5, color = "blue", size = 3)

            # Save original plot
            cross_plot_file <- get_output_path(step_num, "05b_pqtl_sex_crossref", batch_id, "outliers", "pdf", config = config)
            ensure_output_dir(cross_plot_file)
            ggsave(cross_plot_file, p_cross, width = 10, height = 8)
            log_info("Saved cross-reference plot to: {cross_plot_file}")

            # Create faceted plot: separate panels for males and females
            # Filter to only samples with known genetic sex (male or female)
            plot_data_faceted <- plot_data[genetic_sex %in% c("male", "female")]

            # Ensure factor levels: female first (top), then male (bottom)
            if (nrow(plot_data_faceted) > 0) {
                plot_data_faceted[, genetic_sex := factor(genetic_sex, levels = c("female", "male"))]

                # Create faceted version with same color logic
                p_cross_faceted <- ggplot(plot_data_faceted, aes(x = predicted_prob, y = MeanAbsZ)) +
                    # All points with color gradient (using fill for filled shapes)
                    geom_point(aes(fill = point_color, shape = point_shape, color = border_col),
                              size = 2.5, alpha = 0.7, stroke = 0.5) +
                    # Overlay thicker black borders for sex outliers
                    geom_point(data = plot_data_faceted[has_border == TRUE],
                              aes(shape = point_shape),
                              fill = NA, color = "black",
                              size = 3, stroke = 1.5, alpha = 1) +
                    # Facet by genetic sex: vertical layout (nrow = 2), female on top, male on bottom
                    facet_wrap(~ genetic_sex, nrow = 2, scales = "free_x") +
                    # Use identity scale for fill (colors already assigned) - no legend
                    scale_fill_identity(guide = "none") +
                    # Use identity scale for border color (transparent for normal, black for outliers)
                    scale_color_identity(guide = "none") +
                    # Shape scale: triangle for above threshold, circle for below
                    scale_shape_manual(
                        values = c("Above threshold" = 24, "Below threshold" = 21),
                        name = "MeanAbsZ Threshold",
                        labels = c("Above threshold" = paste("Above", cutoff_meanabsz),
                                   "Below threshold" = paste("Below", cutoff_meanabsz))
                    ) +
                    # Vertical line for Youden J threshold (pale red)
                    geom_vline(xintercept = youden_threshold, linetype = "dashed", color = "red", linewidth = 0.8, alpha = 0.5) +
                    # Horizontal line for MeanAbsZ threshold (pale blue)
                    geom_hline(yintercept = cutoff_meanabsz, linetype = "dashed", color = "blue", linewidth = 0.8, alpha = 0.5) +
                    # Labels and theme
                    labs(
                        title = "pQTL Z-scores vs Sex Predictions (Faceted by Genetic Sex)",
                        subtitle = paste("MeanAbsZ threshold:", cutoff_meanabsz, "| Youden J threshold:", round(youden_threshold, 3)),
                        x = "Predicted Female Probability",
                        y = "Mean Absolute Z-score",
                        shape = "Z-score Threshold"
                    ) +
                    theme_bw() +
                    theme(
                        plot.title = element_text(size = 14, face = "bold"),
                        plot.subtitle = element_text(size = 11),
                        legend.position = "right",
                        strip.text = element_text(size = 12, face = "bold")
                    ) +
                    # Annotate thresholds (per facet)
                    annotate("text", x = youden_threshold, y = max(plot_data_faceted$MeanAbsZ, na.rm = TRUE) * 0.95,
                            label = paste("Youden J =", round(youden_threshold, 3)),
                            hjust = -0.1, vjust = 1, color = "red", size = 3) +
                    annotate("text", x = max(plot_data_faceted$predicted_prob, na.rm = TRUE) * 0.95, y = cutoff_meanabsz,
                            label = paste("MeanAbsZ threshold =", cutoff_meanabsz),
                            hjust = 1, vjust = -0.5, color = "blue", size = 3)

                # Save faceted plot
                cross_plot_faceted_file <- get_output_path(step_num, "05b_pqtl_sex_crossref_faceted", batch_id, "outliers", "pdf", config = config)
                ensure_output_dir(cross_plot_faceted_file)
                ggsave(cross_plot_faceted_file, p_cross_faceted, width = 14, height = 7)
                log_info("Saved faceted cross-reference plot to: {cross_plot_faceted_file}")
            } else {
                log_warn("No samples with known genetic sex for faceted plot. Skipping faceted plot generation.")
            }
        } else {
            log_warn("No overlapping samples between pQTL stats and sex predictions. Skipping cross-reference plot.")
        }
    } else {
        log_warn("Sex predictions file not found: {sex_predictions_path}. Skipping cross-reference plot.")
    }

    # 8. Generate Final Cleaned Matrix (All Outliers Removed)
    # --------------------------------------------------------
    log_info("Generating final cleaned matrix with all outliers removed...")

    # Load base matrix from step 00
    base_matrix_path <- get_output_path("00", "npx_matrix_analysis_ready", batch_id, "qc", config = config)
    if (!file.exists(base_matrix_path)) {
        log_warn("Base matrix not found: {base_matrix_path}. Skipping final cleaned matrix generation.")
    } else {
        base_matrix <- readRDS(base_matrix_path)
        log_info("Loaded base matrix: {nrow(base_matrix)} samples x {ncol(base_matrix)} proteins")

        # Collect all outlier SampleIDs from all QC steps
        all_outlier_ids <- character()

        # step 01: PCA outliers
        pca_file <- get_output_path("01", "pca_outliers_by_source", batch_id, "outliers", "tsv", config = config)
        if (file.exists(pca_file)) {
            pca_outliers <- try(fread(pca_file), silent = TRUE)
            if (!inherits(pca_outliers, "try-error") && nrow(pca_outliers) > 0 && "SampleID" %in% names(pca_outliers)) {
                pca_outlier_ids <- pca_outliers[Any == 1]$SampleID
                all_outlier_ids <- unique(c(all_outlier_ids, pca_outlier_ids))
                log_info("Added {length(pca_outlier_ids)} PCA outliers")
            }
        }

        # step 04: Sex outliers
        sex_file <- get_output_path("04", "sex_mismatches", batch_id, "outliers", "tsv", config = config)
        if (file.exists(sex_file)) {
            sex_outliers <- try(fread(sex_file), silent = TRUE)
            if (!inherits(sex_outliers, "try-error") && nrow(sex_outliers) > 0) {
                id_col <- if ("SAMPLE_ID" %in% names(sex_outliers)) "SAMPLE_ID" else "SampleID"
                if (id_col %in% names(sex_outliers)) {
                    sex_outlier_ids <- sex_outliers[[id_col]]
                    all_outlier_ids <- unique(c(all_outlier_ids, sex_outlier_ids))
                    log_info("Added {length(sex_outlier_ids)} sex outliers")
                }
            }
        }

        # step 02: Technical outliers
        tech_file <- get_output_path("02", "technical_outlier_summary", batch_id, "outliers", "tsv", config = config)
        if (file.exists(tech_file)) {
            tech_outliers <- try(fread(tech_file), silent = TRUE)
            if (!inherits(tech_outliers, "try-error") && nrow(tech_outliers) > 0 && "SampleID" %in% names(tech_outliers)) {
                # Build filter condition based on available columns
                tech_cols <- names(tech_outliers)
                is_outlier <- rep(FALSE, nrow(tech_outliers))

                if ("Tech_Any" %in% tech_cols) {
                    is_outlier <- is_outlier | (tech_outliers$Tech_Any == 1)
                }
                if ("is_plate_outlier" %in% tech_cols) {
                    is_outlier <- is_outlier | (tech_outliers$is_plate_outlier == TRUE)
                }
                if ("is_batch_outlier" %in% tech_cols) {
                    is_outlier <- is_outlier | (tech_outliers$is_batch_outlier == TRUE)
                }
                if ("is_processing_outlier" %in% tech_cols) {
                    is_outlier <- is_outlier | (tech_outliers$is_processing_outlier == TRUE)
                }
                if ("is_sample_outlier" %in% tech_cols) {
                    is_outlier <- is_outlier | (tech_outliers$is_sample_outlier == TRUE)
                }

                tech_outlier_ids <- tech_outliers[is_outlier]$SampleID
                if (length(tech_outlier_ids) > 0) {
                    all_outlier_ids <- unique(c(all_outlier_ids, tech_outlier_ids))
                    log_info("Added {length(tech_outlier_ids)} technical outliers")
                }
            }
        }

        # step 03: Z-score outliers
        zscore_file <- get_output_path("03", "zscore_outliers_list", batch_id, "outliers", "tsv", config = config)
        if (file.exists(zscore_file)) {
            zscore_outliers <- try(fread(zscore_file), silent = TRUE)
            if (!inherits(zscore_outliers, "try-error") && nrow(zscore_outliers) > 0 && "SampleID" %in% names(zscore_outliers)) {
                zscore_outlier_ids <- zscore_outliers$SampleID
                all_outlier_ids <- unique(c(all_outlier_ids, zscore_outlier_ids))
                log_info("Added {length(zscore_outlier_ids)} Z-score outliers")
            }
        }

        # step 05b: pQTL outliers (only MeanAbsZ-based)
        # Use SampleID if available, otherwise use FINNGENID
        id_col_pqtl <- if ("SampleID" %in% names(sample_stats)) "SampleID" else if ("FINNGENID" %in% names(sample_stats)) "FINNGENID" else NULL
        if (!is.null(id_col_pqtl)) {
            pqtl_outlier_ids <- sample_stats[Outlier_MeanAbsZ == TRUE][[id_col_pqtl]]
            pqtl_outlier_ids <- pqtl_outlier_ids[!is.na(pqtl_outlier_ids)]
            if (length(pqtl_outlier_ids) > 0) {
                all_outlier_ids <- unique(c(all_outlier_ids, pqtl_outlier_ids))
                log_info("Added {length(pqtl_outlier_ids)} pQTL outliers (MeanAbsZ-based)")
            }
        } else {
            log_warn("Cannot extract pQTL outlier IDs: neither SampleID nor FINNGENID found in sample_stats")
        }

        # Remove all outliers
        all_outlier_ids <- unique(all_outlier_ids)
        log_info("Total unique outliers across all QC steps: {length(all_outlier_ids)}")

        # Filter base matrix
        keep_samples <- setdiff(rownames(base_matrix), all_outlier_ids)
        final_cleaned_matrix <- base_matrix[keep_samples, ]

        log_info("Final cleaned matrix: {nrow(final_cleaned_matrix)} samples x {ncol(final_cleaned_matrix)} proteins")
        log_info("Removed {length(all_outlier_ids)} samples ({round(100*length(all_outlier_ids)/nrow(base_matrix), 2)}%)")

        # Save final cleaned matrix
        final_matrix_path <- get_output_path("05d", "05d_npx_matrix_all_qc_passed", batch_id, "phenotypes", "rds", config = config)
        ensure_output_dir(final_matrix_path)
        saveRDS(final_cleaned_matrix, final_matrix_path)
        log_info("Saved final cleaned matrix to: {final_matrix_path}")

        # Also save pcnorm version if available
        pcnorm_base_path <- get_output_path("00", "pcnorm_matrix_analysis_ready", batch_id, "qc", config = config)
        if (file.exists(pcnorm_base_path)) {
            pcnorm_base <- readRDS(pcnorm_base_path)
            pcnorm_cleaned <- pcnorm_base[keep_samples, ]
            pcnorm_final_path <- get_output_path("05b", "05b_pcnorm_matrix_all_outliers_removed", batch_id, "outliers", "rds", config = config)
            ensure_output_dir(pcnorm_final_path)
            saveRDS(pcnorm_cleaned, pcnorm_final_path)
            log_info("Saved final cleaned PCNorm matrix to: {pcnorm_final_path}")
        }
    }

    log_info("Done.")
}

if (!interactive()) {
    main()
}
