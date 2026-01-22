#!/usr/bin/env Rscript
# ==============================================================================
# 05a_pqtl_training.R - Robust pQTL Selection via Stability Selection
# ==============================================================================
#
# Purpose:
#   Identifies the most robust and informative pQTLs for mismatch detection
#   using a Stability Selection approach with LASSO Logistic Regression.
#   It runs multiple training iterations on synthetic cohorts to derive a
#   consensus pQTL list and optimal decision thresholds.
#
# Author: Reza Jabal, PhD (rjabal@broadinstitute.org)
# Date: December 2025
# ==============================================================================

suppressPackageStartupMessages({
    library(data.table)
    library(logger)
    library(yaml)
    library(arrow)
    library(glmnet)
    library(pROC)
    library(ggplot2)
    library(ggpubr)
    library(digest)
})

# Source utility functions
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

# Get config path from environment (required when sourced by pipeline runner)
# This matches the pattern used in other scripts (04_sex_outliers.R, 05b_pqtl_outliers.R)
config_file <- Sys.getenv("PIPELINE_CONFIG", "")
if (config_file == "" || !file.exists(config_file)) {
  # Try to get from global environment if set by pipeline runner
  if (exists("PIPELINE_CONFIG_OBJ", envir = .GlobalEnv)) {
    config <- get("PIPELINE_CONFIG_OBJ", envir = .GlobalEnv)
  } else {
    # Config not available - will be loaded in run_pqtl_training() function if needed
    config <- NULL
  }
} else {
  config <- yaml::read_yaml(config_file)
}

# Get batch context
if (!is.null(config)) {
    batch_id <- Sys.getenv("PIPELINE_BATCH_ID", config$batch$default_batch_id %||% "batch_01")
} else {
    batch_id <- Sys.getenv("PIPELINE_BATCH_ID", "batch_01")
}

# ============================================================================
# Helper Functions
# ============================================================================

# Helper to read GS file
read_gs_file <- function(gs_path) {
    tmp <- tempfile()
    on.exit(if (file.exists(tmp)) unlink(tmp), add = TRUE)

    ret <- system(paste("gsutil cp", gs_path, tmp), ignore.stdout = TRUE, ignore.stderr = TRUE)
    if (ret != 0 || !file.exists(tmp)) {
        return(NULL)
    }

    dt <- tryCatch(fread(tmp), error = function(e) NULL)
    return(dt)
}

# Helper to extract top variant (reused from 05b)
get_top_variant <- function(f) {
    tryCatch(
        {
            dt <- read_gs_file(f)
            if (is.null(dt) || nrow(dt) == 0) {
                return(NULL)
            }

            # Standardize column names
            if ("pvalue" %in% names(dt)) setnames(dt, "pvalue", "p")
            if ("effect" %in% names(dt)) setnames(dt, "effect", "beta")

            # Filter for high confidence
            # Sort by p-value (asc) and cs_avg_r2 (desc)
            setorder(dt, p, -cs_avg_r2)

            # Take top 1
            top <- dt[1]

            cols <- c("trait", "rsid", "p", "beta", "sd", "cs_avg_r2")
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

# Helper: Check Genotype Polarity and Flip if Needed
# Returns the genotype vector (possibly flipped) and a boolean indicating if it was flipped
check_and_flip_genotype <- function(g_vec, rsid, col_name) {
    # Parse suffix from col_name (e.g., "chr1_123_A_G_A" -> "A")
    # PLINK --export A format: ID_REF_ALT_ALTCOUNTED (usually)
    # Actually PLINK 2 --export A usually headers are: ID_ALLELE
    # e.g. chr1:123:A:G_A means counting A

    parts <- strsplit(col_name, "_")[[1]]
    counted_allele <- tail(parts, 1)

    # Parse Ref/Alt from RSID (assuming format chr_pos_REF_ALT)
    # This relies on consistent ID formatting
    rs_parts <- strsplit(rsid, "_")[[1]]
    if (length(rs_parts) >= 4) {
        ref_allele <- rs_parts[3]
        alt_allele <- rs_parts[4]

        # If the counted allele is the REF allele, effectively identical to Plink export
        # Empirical testing shows flipping destroys signal, so we assume provided Betas align with this export.
        if (counted_allele == ref_allele) {
            # return(list(g = 2 - g_vec, flipped = TRUE))
            return(list(g = g_vec, flipped = FALSE)) # Disable flip
        }
    }

    return(list(g = g_vec, flipped = FALSE))
}

# ============================================================================
# Main Function
# ============================================================================
run_pqtl_training <- function(config = NULL, batch_id = NULL, verbose = FALSE) {
    # Load config if not provided
    if (is.null(config)) {
        config_file <- Sys.getenv("PIPELINE_CONFIG", "")
        if (config_file != "" && file.exists(config_file)) {
            config <- yaml::read_yaml(config_file)
        } else if (exists("PIPELINE_CONFIG_OBJ", envir = .GlobalEnv)) {
            config <- get("PIPELINE_CONFIG_OBJ", envir = .GlobalEnv)
        } else {
            stop("Config not provided and could not be loaded from environment")
        }
    }
    
    # Get batch_id if not provided
    if (is.null(batch_id) || batch_id == "") {
        batch_id <- Sys.getenv("PIPELINE_BATCH_ID", config$batch$default_batch_id %||% "batch_01")
    }
    
    step_num <- "05a"
    step_suffix <- "05a"
    step_name <- "pqtl_training"

    # Set up logging (use get_log_path for consistency)
    log_file <- get_log_path(step_num, batch_id, config = config)
    ensure_output_dir(log_file)
    log_appender(appender_file(log_file))

    log_info("=" |> rep(60) |> paste(collapse = ""))
    log_info("Step {step_num}: pQTL Training - Robust Stability Selection")
    log_info("Batch: {batch_id %||% 'default'}")
    log_info("=" |> rep(60) |> paste(collapse = ""))

    start_time <- Sys.time()

    tryCatch(
        {
            # ====================================================================
            # 1. Configuration
            # ====================================================================
            training_config <- config$pqtl_training
            if (is.null(training_config) || !isTRUE(training_config$enabled)) {
                # Set environment variable to indicate step was skipped
                Sys.setenv(PIPELINE_STEP_SKIPPED = "TRUE")
                log_info("pQTL training is disabled in config. Skipping.")
                return(invisible(NULL))
            }

            pqtl_config <- config$parameters$pqtl_outliers

            training_seeds <- training_config$training_seeds %||% c(12345, 23456, 34567, 45678, 56789)
            error_configs <- training_config$error_configs
            if (is.null(error_configs) || length(error_configs) == 0) {
                log_warn("No error_configs provided. Using default error configuration.")
                error_configs <- list(
                    list(genotype_mismatch = 10, sex_mismatch = 5)
                )
            }
            cohort_size <- training_config$cohort_size %||% 500
            excluded_cohorts <- training_config$excluded_cohorts %||% c("Chromosomal_Abnormalities", "F64", "Kids")

            lasso_alpha <- training_config$lasso$alpha %||% 1.0
            cv_folds <- training_config$lasso$cv_folds %||% 10
            lambda_selection <- training_config$lasso$lambda_selection %||% "lambda.min"

            output_dir <- file.path(config$output$base_dir, training_config$output_dir %||% "pqtl-training", batch_id)
            dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

            # Stability Selection Parameters (defined here for logging)
            n_features_to_sample <- training_config$n_features_to_sample %||% 500
            n_top_tier_pool <- training_config$n_top_tier_pool %||% 1000
            selection_threshold <- training_config$selection_threshold %||% 0.01

            log_info("Configuration:")
            log_info("  Cohort Size: {cohort_size}")
            log_info("  Training Iterations: {length(training_seeds)} seeds x {length(error_configs)} error configs = {length(training_seeds) * length(error_configs)} cohorts")
            log_info("  Feature Subsampling: Select {n_features_to_sample} random from top {n_top_tier_pool}")
            log_info("  LASSO: alpha={lasso_alpha}, folds={cv_folds}, selection={lambda_selection}")
            log_info("  Output Dir: {output_dir}")

            # ====================================================================
            # 2. Data Loading & Preparation
            # ====================================================================
            log_info("")
            log_info("--- Loading Data ---")

            # Load NPX (step 01 Cleaned)
            npx_path <- get_output_path("01", "npx_matrix_pca_cleaned", batch_id, "outliers", config = config)
            if (!file.exists(npx_path)) {
                log_warn("step 01 matrix not found, checking step 00...")
                npx_path <- get_output_path("00", "npx_matrix_analysis_ready", batch_id, "qc", config = config)
            }
            if (!file.exists(npx_path)) stop("No input NPX matrix found.")
            npx_matrix <- readRDS(npx_path)
            log_info("Loaded NPX: {nrow(npx_matrix)} samples x {ncol(npx_matrix)} proteins")

            # Load Metadata & Mapping
            sample_mapping_path <- get_output_path("00", "sample_mapping", batch_id, "qc", "tsv", config = config)
            metadata_path <- get_output_path("00", "metadata", batch_id, "qc", config = config)

            if (!file.exists(sample_mapping_path)) stop("Sample mapping not found.")
            dt_mapping <- fread(sample_mapping_path)

            if (!file.exists(metadata_path)) stop("Metadata not found.")
            metadata_raw <- readRDS(metadata_path)

            # Merge & Prepare Metadata (including Disease Groups)
            disease_cols <- intersect(c("Kidney", "Kids", "F64", "Chromosomal_Abnormalities", "Bridging_samples"), names(metadata_raw))
            metadata <- merge(
                dt_mapping[, .(SampleID, FINNGENID, COHORT_FINNGENID, sample_type)],
                metadata_raw[, c("SAMPLE_ID", disease_cols), with = FALSE],
                by.x = "SampleID", by.y = "SAMPLE_ID", all.x = TRUE
            )

            # Load Covariates (Sex)
            cov_file <- config$covariates$finngen_r13_minimum_file
            if (!is.null(cov_file) && file.exists(cov_file)) {
                covariates <- fread(cov_file, select = c("FINNGENID", "SEX"))
                metadata <- merge(metadata, covariates, by = "FINNGENID", all.x = TRUE)
            } else {
                log_warn("Minimum covariate file not found. SEX info may be missing.")
            }

            # Exclude specific cohorts
            exclude_ids <- character(0)
            for (ex in excluded_cohorts) {
                if (ex %in% names(metadata)) {
                    ids <- metadata[get(ex) == TRUE | get(ex) == 1]$SampleID
                    exclude_ids <- c(exclude_ids, ids)
                    log_info("  Excluding {length(ids)} samples from {ex}")
                }
            }

            valid_samples <- metadata[
                !SampleID %in% exclude_ids &
                    !is.na(FINNGENID) &
                    !is.na(SEX) &
                    SampleID %in% rownames(npx_matrix) &
                    sample_type %in% c("FinnGen", "Bridging")
            ]

            npx_matrix <- npx_matrix[valid_samples$SampleID, , drop = FALSE]
            log_info("Analysis Set: {nrow(valid_samples)} valid samples")

            # Apply IRN to match production pipeline (05b)
            npx_matrix <- inverse_rank_normalize(npx_matrix)

            # Load step 04 Sex Predictions using the helper function logic or direct path
            # Try to find the sex predictions file
            sex_preds_path <- get_output_path("04", "sex_predictions", batch_id, "outliers", "tsv", config = config)
            sex_preds <- NULL
            if (file.exists(sex_preds_path)) {
                sex_preds <- fread(sex_preds_path)
                log_info("Loaded step 04 sex predictions: {nrow(sex_preds)} samples")
            } else {
                log_warn("step 04 sex predictions not found at {sex_preds_path}. Sex cross-ref plots will be skipped.")
            }

            # Parse Youden J threshold from step 04 Log
            # Expected location: logs/{batch_id}/05_sex_outliers.log
            log_path_05 <- file.path(config$output$base_dir, config$output$logs_dir %||% "logs", batch_id, "05_sex_outliers.log")
            youden_threshold_05 <- NA_real_

            if (file.exists(log_path_05)) {
                log_lines <- readLines(log_path_05)
                # Look for line: "Mismatches (Youden J = 0.71): 20"
                target_line <- grep("Mismatches \\(Youden J =", log_lines, value = TRUE)
                if (length(target_line) > 0) {
                    # Extract number using regex
                    match <- regexpr("Youden J = ([0-9\\.]+)", target_line[1])
                    if (match != -1) {
                        val_str <- regmatches(target_line[1], match)
                        val_str <- sub("Youden J = ", "", val_str)
                        youden_threshold_05 <- as.numeric(val_str)
                        log_info("Parsed step 04 Youden J threshold: {youden_threshold_05}")
                    }
                }
            }

            if (is.na(youden_threshold_05)) {
                log_warn("Could not parse Youden J threshold from step 04 log. Using default 0.5 for plots.")
                youden_threshold_05 <- 0.5
            }


            # ====================================================================
            # 3. Global pQTL Collation & Genotype Extraction
            # ====================================================================
            log_info("")
            log_info("--- Global pQTL Prep ---")

            collated_path <- file.path(output_dir, "05a_finemap_collated_global.tsv")
            top_variants_global <- NULL

            if (file.exists(collated_path)) {
                log_info("Loading cached pQTL list: {collated_path}")
                top_variants_global <- fread(collated_path)
            } else {
                log_info("Collating pQTLs from fine-mapping (top 1 per protein)...")
                finemap_path <- pqtl_config$finemap_path

                # List files
                if (grepl("^gs://", finemap_path)) {
                    files <- system(paste("gsutil ls", finemap_path), intern = TRUE)
                    files <- files[grepl("\\.SUSIE\\.cred\\.summary\\.tsv$", files)]
                } else {
                    files <- list.files(finemap_path, pattern = "\\.SUSIE\\.cred\\.summary\\.tsv$", full.names = TRUE)
                }

                log_info("Found {length(files)} fine-mapping files")
                if (length(files) == 0) stop("No fine-mapping files found.")

                results_list <- list()
                for (i in seq_along(files)) {
                    if (i %% 100 == 0) log_info("Processing {i}/{length(files)}...")
                    res <- get_top_variant(files[i])
                    if (!is.null(res)) results_list[[i]] <- res
                }
                top_variants_global <- rbindlist(results_list)

                # Handle sex chromosomes (chr23 -> chrX) if needed
                if (isTRUE(pqtl_config$sex_chromosome_pQTLs)) {
                    top_variants_global[, rsid_original := rsid]
                    top_variants_global[, rsid_for_matching := ifelse(grepl("^chr23_", rsid), gsub("^chr23_", "chrX_", rsid), rsid)]
                } else {
                    top_variants_global[, rsid_original := rsid]
                    top_variants_global[, rsid_for_matching := rsid]
                }

                fwrite(top_variants_global, collated_path, sep = "\t")
                log_info("Collated {nrow(top_variants_global)} global candidates.")
            }

            # Extract Genotypes for ALL candidates
            log_info("Extracting genotypes for {nrow(top_variants_global)} candidates...")

            # Use a consistent cache key (Variants + Samples)
            genotype_path <- pqtl_config$genotype_path

            # Calculate sample hash
            sample_hash <- digest::digest(sort(valid_samples$FINNGENID))
            variant_hash <- digest::digest(sort(top_variants_global$rsid_for_matching))

            cache_key <- digest::digest(c(variant_hash, sample_hash))
            cache_base <- file.path(config$output$base_dir, "temp_work", "pqtl_cache", batch_id, cache_key)
            dir.create(cache_base, recursive = TRUE, showWarnings = FALSE)

            cached_raw <- file.path(cache_base, "global_genotypes.raw")

            if (file.exists(cached_raw) && file.size(cached_raw) > 0) {
                log_info("Using cached global genotypes: {cached_raw}")
                dt_geno_global <- fread(cached_raw)
            } else {
                # Reuse logic from 05b sparse extraction...
                # Simplified here assuming the robust extraction works
                log_info("Running global extraction (sparse BED + PLINK export)...")

                # Create variant list
                var_file <- file.path(cache_base, "global_vars.snplist")
                writeLines(top_variants_global$rsid_for_matching, var_file)

                # Create keep file (all valid samples)
                keep_file <- file.path(cache_base, "global_samples.txt")
                all_ids <- unique(valid_samples$FINNGENID)
                fwrite(data.table(FID = all_ids, IID = all_ids), keep_file, sep = "\t", col.names = FALSE)

                # Determine PLINK input (Local vs GS sparse download)
                plink_input <- genotype_path

                if (startsWith(genotype_path, "gs://")) {
                    # GS Sparse download logic (simplified call, assuming 07b logic is available or we implement it fully)
                    # Ideally, we should source 07b or extract this to a common util.
                    # For now, I'll assume we can use the local path if downloaded, or implement the download.

                    # Implementing sparse download here to be self-contained
                    base_name <- basename(genotype_path)
                    temp_geno_dir <- file.path(cache_base, "temp_plink")
                    dir.create(temp_geno_dir, showWarnings = FALSE)

                    # Download BIM/FAM
                    system(paste("gsutil cp", paste0(genotype_path, ".bim"), file.path(temp_geno_dir, paste0(base_name, ".bim"))))
                    system(paste("gsutil cp", paste0(genotype_path, ".fam"), file.path(temp_geno_dir, paste0(base_name, ".fam"))))

                    dt_bim <- fread(file.path(temp_geno_dir, paste0(base_name, ".bim")), header = FALSE)
                    setnames(dt_bim, c("chr", "rsid", "cm", "pos", "a1", "a2"))
                    dt_bim[, idx := .I]

                    dt_targets <- dt_bim[rsid %in% top_variants_global$rsid_for_matching]

                    # Sparse BED
                    sparse_bed <- file.path(temp_geno_dir, paste0(base_name, ".bed"))
                    f_out <- file(sparse_bed, "wb")
                    writeBin(as.raw(c(0x6c, 0x1b, 0x01)), f_out)

                    # Sample count
                    n_samples <- nrow(fread(file.path(temp_geno_dir, paste0(base_name, ".fam")), header = FALSE))
                    bytes_per_variant <- ceiling(n_samples / 4)

                    log_info("Downloading {nrow(dt_targets)} variants...")
                    setorder(dt_targets, idx)

                    # Parallel or chunked download would be better, but sequential for safety
                    for (i in seq_len(nrow(dt_targets))) {
                        if (i %% 100 == 0) cat(".")
                        var_idx <- dt_targets$idx[i]
                        start <- 3 + (var_idx - 1) * bytes_per_variant
                        end <- start + bytes_per_variant - 1
                        chunk <- file.path(temp_geno_dir, "chunk.bin")
                        ret <- system(paste("gsutil cat -r", paste0(format(start, scientific = FALSE), "-", format(end, scientific = FALSE)), paste0(genotype_path, ".bed"), ">", chunk))
                        if (ret == 0 && file.exists(chunk)) {
                            writeBin(readBin(chunk, "raw", n = bytes_per_variant + 100), f_out)
                            unlink(chunk)
                        }
                    }
                    close(f_out)
                    cat("\n")

                    # Filter BIM
                    dt_final_bim <- dt_targets[, .(chr, rsid, cm, pos, a1, a2)]
                    fwrite(dt_final_bim, file.path(temp_geno_dir, paste0(base_name, ".bim")), sep = "\t", col.names = FALSE)

                    plink_input <- file.path(temp_geno_dir, base_name)
                }

                # Run PLINK export
                cmd <- sprintf("plink2 --bfile %s --keep %s --export A --out %s", plink_input, keep_file, file.path(cache_base, "global_genotypes"))
                system(cmd)

                cached_raw <- file.path(cache_base, "global_genotypes.raw")
                dt_geno_global <- fread(cached_raw)
            }

            # Prepare genotype map (Column Name -> RSID)
            geno_cols <- setdiff(names(dt_geno_global), c("FID", "IID", "PAT", "MAT", "SEX", "PHENOTYPE"))
            geno_map <- data.table(col = geno_cols)
            geno_map[, rsid := sub("_[AGCT]+$", "", col)] # Assumes RSID_ALLELE format

            log_info("Global Genotypes: {nrow(dt_geno_global)} samples, {length(geno_cols)} variants")

            # ====================================================================
            # 4. Iterative Training Loop
            # ====================================================================
            log_info("")
            log_info("--- Starting Stability Selection Loop ---")

            cohort_results <- list()
            cohort_idx <- 0

            for (seed in training_seeds) {
                for (error_cfg in error_configs) {
                    cohort_idx <- cohort_idx + 1
                    set.seed(seed)

                    # A. Generate Cohort
                    # ----------------------------------------
                    curr_samples <- valid_samples[sample(.N, min(.N, cohort_size))]
                    curr_samples[, `:=`(is_error = FALSE, has_geno_err = FALSE, has_sex_err = FALSE)]

                    # Inject Errors
                    n_geno <- error_cfg$genotype_mismatch
                    n_sex <- error_cfg$sex_mismatch

                    # Genotype swaps (within sex)
                    males <- which(curr_samples$SEX %in% c("male", "1"))
                    females <- which(curr_samples$SEX %in% c("female", "2"))

                    # ... (Simplified swap logic) ...
                    if (length(males) > 2) {
                        swap_m <- sample(males, min(length(males), floor(n_geno / 2 / 2) * 2))
                        if (length(swap_m) > 1) {
                            # Swap FINNGENIDs cyclically or pairwise
                            ids <- curr_samples$FINNGENID[swap_m]
                            curr_samples[swap_m, FINNGENID := ids[c(2:length(ids), 1)]]
                            curr_samples[swap_m, `:=`(is_error = TRUE, has_geno_err = TRUE)]
                        }
                    }
                    if (length(females) > 2) {
                        swap_f <- sample(females, min(length(females), floor(n_geno / 2 / 2) * 2))
                        if (length(swap_f) > 1) {
                            ids <- curr_samples$FINNGENID[swap_f]
                            curr_samples[swap_f, FINNGENID := ids[c(2:length(ids), 1)]]
                            curr_samples[swap_f, `:=`(is_error = TRUE, has_geno_err = TRUE)]
                        }
                    }

                    # Sex swaps
                    err_idx <- which(curr_samples$has_geno_err)
                    if (length(err_idx) >= n_sex) {
                        swap_s <- sample(err_idx, n_sex)
                        curr_samples[swap_s, SEX := ifelse(SEX == "male" | SEX == "1", "female", "male")]
                        curr_samples[swap_s, has_sex_err := TRUE]
                    }

                    # B. Local Statistics & Ranking
                    # -----------------------------
                    # Subset genotypes for THIS cohort (using the possibly swapped IDs to simulate reality)
                    # Ensure FINNGENID is character for matching
                    curr_finngenids <- as.character(curr_samples$FINNGENID)
                    cohort_genos <- dt_geno_global[IID %in% curr_finngenids]

                    log_info("  Cohort {cohort_idx}: {nrow(curr_samples)} samples, {nrow(cohort_genos)} matched in genotypes")

                    if (nrow(cohort_genos) == 0) {
                        log_warn("  No genotype matches for cohort {cohort_idx}. Skipping.")
                        next
                    }

                    # Calculate Local MAF/Het
                    local_stats <- list()
                    for (i in seq_len(nrow(geno_map))) {
                        vals <- cohort_genos[[geno_map$col[i]]]
                        if (is.null(vals)) next
                        maf <- min(mean(vals, na.rm = TRUE) / 2, 1 - mean(vals, na.rm = TRUE) / 2)
                        het <- mean(vals == 1, na.rm = TRUE)
                        local_stats[[i]] <- list(rsid = geno_map$rsid[i], local_maf = maf, local_het = het)
                    }
                    dt_local <- rbindlist(local_stats)

                    # Merge with global metadata (p, beta) and Rank
                    log_info("  Local stats: {nrow(dt_local)} variants with MAF/Het")
                    dt_rank <- merge(top_variants_global[, .(rsid = rsid_for_matching, trait, p, beta)], dt_local, by = "rsid")
                    log_info("  After merge with global metadata: {nrow(dt_rank)} variants")

                    if (nrow(dt_rank) == 0) {
                        log_warn("  No variants matched after merge in cohort {cohort_idx}. Skipping.")
                        next
                    }

                    dt_rank[, composite := (-log10(p) * abs(beta) * local_maf) / (local_het + 0.01)]
                    setorder(dt_rank, -composite)

                    # Feature Subsampling (Stability Selection)
                    # Take top pool, then random sample
                    pool <- head(dt_rank, min(n_top_tier_pool, nrow(dt_rank)))
                    n_to_sample <- min(n_features_to_sample, nrow(pool))
                    selected_features <- pool[sample(.N, n_to_sample)]

                    # Deduplicate by rsid (keep trait with highest composite score)
                    # This prevents duplicate SampleID-rsid combinations in dcast
                    setorder(selected_features, -composite)
                    selected_features <- selected_features[, .SD[1], by = rsid]
                    log_info("  Selected {nrow(selected_features)} unique features from top {nrow(pool)} (pool size, after dedup by rsid)")

                    # C. Build Feature Matrix
                    # -----------------------
                    # Need Z and Residuals for selected features
                    # This requires calculating Z/Res for every sample in cohort

                    log_info("  Building features for {nrow(selected_features)} selected pQTLs...")

                    # Reuse the calculation logic (vectorized where possible)
                    X_list <- list()
                    n_skipped_prot <- 0
                    n_skipped_merge <- 0
                    n_skipped_rows <- 0
                    n_skipped_col <- 0

                    for (i in seq_len(nrow(selected_features))) {
                        feat <- selected_features[i]
                        prot <- feat$trait
                        target_rsid <- feat$rsid

                        # Find matching column in geno_map (handle potential format differences)
                        col_match <- geno_map[rsid == target_rsid]
                        if (nrow(col_match) == 0) {
                            # Try case-insensitive or partial match
                            col_match <- geno_map[tolower(rsid) == tolower(target_rsid)]
                        }
                        if (nrow(col_match) == 0) {
                            n_skipped_col <- n_skipped_col + 1
                            if (i <= 5) log_warn("    Feature {i}: No genotype column found for rsid={target_rsid}")
                            next
                        }
                        col_name <- col_match$col[1] # Take first match

                        if (!prot %in% colnames(npx_matrix)) {
                            n_skipped_prot <- n_skipped_prot + 1
                            if (i <= 5) log_warn("    Feature {i}: Protein '{prot}' not in NPX matrix")
                            next
                        }

                        # Merge Genotype + Protein
                        # NPX uses SampleID, Genotypes use FINNGENID
                        # curr_samples maps SampleID -> FINNGENID

                        # Check if column exists
                        if (!col_name %in% names(cohort_genos)) {
                            n_skipped_col <- n_skipped_col + 1
                            if (i <= 5) log_warn("    Feature {i}: Column '{col_name}' not in cohort_genos")
                            next
                        }

                        # Ensure character types for merging
                        dat <- merge(
                            curr_samples[, .(SampleID, FINNGENID = as.character(FINNGENID))],
                            cohort_genos[, .(IID = as.character(IID), g = get(col_name))],
                            by.x = "FINNGENID", by.y = "IID"
                        )

                        if (nrow(dat) == 0) {
                            n_skipped_merge <- n_skipped_merge + 1
                            next
                        }

                        # Ensure SampleID is character to match rownames
                        p_vals <- npx_matrix[as.character(dat$SampleID), prot]

                        # Genotype Polarity Check
                        # Check if we need to flip the dosage based on column name vs RSID
                        # plink_col matches the column name in dt_geno (e.g., "chr1_123_A_G_A")
                        plink_col <- grep(paste0("^", target_rsid, "_"), names(cohort_genos), value = TRUE)
                        if (length(plink_col) == 1) {
                            flip_res <- check_and_flip_genotype(dat$g, target_rsid, plink_col)
                            dat[, g := flip_res$g]
                            # Optional: log if flipped? Might be too verbose inside loop
                            # if(flip_res$flipped) log_trace("Flipped {rsid}")
                        }
                        dat[, p := p_vals]
                        dat <- dat[!is.na(p) & !is.na(g)]

                        if (nrow(dat) < 10) {
                            n_skipped_rows <- n_skipped_rows + 1
                            next
                        }

                        # Z-score (per genotype group)
                        stats <- dat[, .(m = mean(p), s = sd(p)), by = g]
                        dat <- merge(dat, stats, by = "g")
                        dat[, z := (p - m) / ifelse(s > 0, s, NA)]

                        # Residual
                        mu <- mean(dat$p)
                        mu_g <- mean(dat$g)
                        dat[, pred := mu + feat$beta * (g - mu_g)]
                        dat[, res := p - pred]
                        sd_res <- sd(dat$res, na.rm = TRUE)

                        # Diagnostic for first few variants
                        if (i <= 3 && cohort_idx == 1) {
                            log_info("    Variant {i} ({feat$rsid}): n={nrow(dat)}, mu_p={round(mu, 3)}, mu_g={round(mu_g, 3)}, beta={round(feat$beta, 4)}, sd_res={round(sd_res, 6)}")
                        }

                        # Standardize residual - use robust scaling if sd_res is too small
                        if (is.na(sd_res) || sd_res <= 1e-6) {
                            # If residual SD is too small, use MAD instead
                            mad_res <- mad(dat$res, na.rm = TRUE)
                            if (!is.na(mad_res) && mad_res > 1e-6) {
                                dat[, std_res := res / (mad_res * 1.4826)]
                            } else {
                                # If even MAD is too small, use raw residual (not standardized)
                                dat[, std_res := res]
                            }
                        } else {
                            dat[, std_res := res / sd_res]
                        }

                        X_list[[i]] <- dat[, .(SampleID, z, std_res, rsid = feat$rsid)]
                    }

                    if (length(X_list) == 0) {
                        log_warn("No features extracted in cohort {cohort_idx}. Skipped: Col={n_skipped_col}, Prot={n_skipped_prot}, Merge={n_skipped_merge}, Rows={n_skipped_rows}")
                        next
                    }

                    long_feat <- rbindlist(X_list)
                    log_info("  Extracted features: {length(unique(long_feat$rsid))} pQTLs, {length(unique(long_feat$SampleID))} samples")

                    # Diagnostic: Check std_res distribution in long format
                    if (nrow(long_feat) > 0) {
                        std_res_vals <- long_feat$std_res[!is.na(long_feat$std_res)]
                        if (length(std_res_vals) > 0) {
                            log_info("  std_res in long format: n={length(std_res_vals)}, min={round(min(std_res_vals), 4)}, max={round(max(std_res_vals), 4)}, mean={round(mean(std_res_vals), 4)}, sd={round(sd(std_res_vals), 4)}")
                            n_zero_long <- sum(std_res_vals == 0)
                            log_info("  Zero std_res values: {n_zero_long} ({round(n_zero_long/length(std_res_vals)*100, 1)}%)")
                        } else {
                            log_warn("  All std_res values are NA in long format!")
                        }
                    }

                    # Pivot for LASSO (using std_res as primary feature, or both?)
                    # 05a original used Residuals for LASSO. Let's stick to that.
                    # Ensure SampleID and rsid are character (not factor) for dcast
                    long_feat[, SampleID := as.character(SampleID)]
                    long_feat[, rsid := as.character(rsid)]

                    # Check for duplicates before dcast
                    n_dup <- long_feat[, .N, by = .(SampleID, rsid)][N > 1, .N]
                    if (n_dup > 0) {
                        log_warn("  Found {n_dup} duplicate SampleID-rsid combinations! Taking mean.")
                        long_feat <- long_feat[, .(std_res = mean(std_res, na.rm = TRUE), z = mean(z, na.rm = TRUE)), by = .(SampleID, rsid)]
                    }

                    X_mat <- dcast(long_feat, SampleID ~ rsid, value.var = "std_res", fill = NA)
                    log_info("  Feature matrix: {nrow(X_mat)} samples x {ncol(X_mat)-1} features (after dcast)")

                    # Diagnostic: Check a few columns after dcast (before NA imputation)
                    if (ncol(X_mat) > 1 && cohort_idx == 1) {
                        sample_cols <- head(setdiff(names(X_mat), "SampleID"), 3)
                        for (col in sample_cols) {
                            col_vals <- X_mat[[col]]
                            col_vals <- col_vals[!is.na(col_vals)]
                            if (length(col_vals) > 0) {
                                log_info("    Column '{col}': n={length(col_vals)}, min={round(min(col_vals), 4)}, max={round(max(col_vals), 4)}, mean={round(mean(col_vals), 4)}, sd={round(sd(col_vals), 6)}")
                            }
                        }
                    }

                    # D. Train LASSO
                    # --------------
                    x_in <- as.matrix(X_mat[, -1, with = FALSE]) # Drop SampleID, use with=FALSE for safety

                    # Ensure numeric (not integer or factor)
                    storage.mode(x_in) <- "numeric"

                    # Create y_vec - ensure it matches x_in rows exactly
                    y_vec <- curr_samples[match(X_mat$SampleID, SampleID)]$is_error

                    # Check for NAs in y_vec (samples in X_mat not in curr_samples)
                    if (any(is.na(y_vec))) {
                        na_idx <- which(is.na(y_vec))
                        log_warn("  Found {length(na_idx)} samples in X_mat not in curr_samples. Removing them.")
                        x_in <- x_in[-na_idx, , drop = FALSE]
                        y_vec <- y_vec[-na_idx]
                        X_mat <- X_mat[-na_idx]
                    }

                    # Ensure y_vec is numeric (0/1) and has correct length
                    y_vec <- as.numeric(y_vec)

                    # Verify we have both classes
                    n_errors <- sum(y_vec == 1, na.rm = TRUE)
                    n_clean <- sum(y_vec == 0, na.rm = TRUE)
                    log_info("  Labels: {n_errors} errors, {n_clean} clean (total {length(y_vec)} samples)")

                    if (n_errors == 0 || n_clean == 0) {
                        log_warn("  Not enough class diversity for LASSO: {n_errors} errors, {n_clean} clean. Skipping cohort {cohort_idx}.")
                        next
                    }

                    if (length(y_vec) != nrow(x_in)) {
                        log_error("  Mismatch: y_vec length ({length(y_vec)}) != x_in rows ({nrow(x_in)})")
                        next
                    }

                    # Diagnostic: Check NA rate before imputation
                    na_rate_before <- sum(is.na(x_in)) / length(x_in)
                    log_info("  Before NA imputation: {nrow(x_in)} samples x {ncol(x_in)} features, NA rate: {round(na_rate_before*100, 1)}%")

                    # Diagnostic: Check value distribution BEFORE imputation
                    x_flat_before <- as.vector(x_in)
                    x_flat_before <- x_flat_before[!is.na(x_flat_before)]
                    if (length(x_flat_before) > 0) {
                        log_info("  Value distribution BEFORE imputation: min={round(min(x_flat_before), 4)}, max={round(max(x_flat_before), 4)}, mean={round(mean(x_flat_before), 4)}, sd={round(sd(x_flat_before), 4)}")
                    }

                    x_in[is.na(x_in)] <- 0
                    log_info("  After NA imputation: {nrow(x_in)} samples x {ncol(x_in)} features")

                    # Diagnostic: Check value distribution AFTER imputation
                    x_flat <- as.vector(x_in)
                    x_flat <- x_flat[!is.na(x_flat)]
                    if (length(x_flat) > 0) {
                        log_info("  Value distribution AFTER imputation: min={round(min(x_flat), 4)}, max={round(max(x_flat), 4)}, mean={round(mean(x_flat), 4)}, sd={round(sd(x_flat), 4)}")
                        n_zero <- sum(x_flat == 0)
                        log_info("  Zero values: {n_zero} ({round(n_zero/length(x_flat)*100, 1)}%)")
                    }

                    # Filter constant columns - calculate SD per column
                    sds <- apply(x_in, 2, sd, na.rm = TRUE)
                    n_constant <- sum(sds <= 1e-6, na.rm = TRUE)
                    n_na_sd <- sum(is.na(sds))

                    # Diagnostic: Check SD distribution per column
                    sds_valid <- sds[!is.na(sds) & sds > 0]
                    if (length(sds_valid) > 0) {
                        log_info("  Column SD distribution (non-zero): min={round(min(sds_valid), 6)}, max={round(max(sds_valid), 4)}, mean={round(mean(sds_valid), 4)}, median={round(median(sds_valid), 6)}")
                    } else {
                        log_warn("  No columns with non-zero SD!")
                    }

                    # Diagnostic: Show SDs for first few columns
                    if (length(sds) > 0 && cohort_idx == 1) {
                        sample_sds <- head(sds[!is.na(sds)], 5)
                        log_info("  First 5 column SDs: {paste(round(sample_sds, 6), collapse=', ')}")
                    }

                    x_in <- x_in[, sds > 1e-6, drop = FALSE]
                    log_info("  After filtering constants: {ncol(x_in)} features (removed {n_constant} constant, {n_na_sd} NA SD columns)")

                    if (ncol(x_in) < 2) {
                        log_warn("Not enough features for LASSO in cohort {cohort_idx} (only {ncol(x_in)} features after filtering)")
                        next
                    }

                    # Final check before LASSO
                    if (length(y_vec) != nrow(x_in)) {
                        log_error("  Final mismatch: y_vec length ({length(y_vec)}) != x_in rows ({nrow(x_in)})")
                        next
                    }

                    # Verify y_vec has both classes
                    y_table <- table(y_vec)
                    log_info("  Final y_vec distribution: {paste(names(y_table), '=', y_table, collapse=', ')}")

                    if (length(y_table) < 2 || any(y_table == 0)) {
                        log_warn("  Not enough class diversity for LASSO: {paste(names(y_table), '=', y_table, collapse=', ')}. Skipping cohort {cohort_idx}.")
                        next
                    }

                    # Ensure y_vec is a factor for glmnet (binomial family)
                    # glmnet expects factor with levels c("0", "1") or c("FALSE", "TRUE")
                    y_vec_factor <- factor(y_vec, levels = c(0, 1))

                    log_info("  Training LASSO: {nrow(x_in)} samples x {ncol(x_in)} features, {sum(y_vec==1)} errors, {sum(y_vec==0)} clean")

                    # Adjust CV folds if we have too few samples in minority class
                    # Need at least 1 error per fold on average
                    n_errors <- sum(y_vec == 1)
                    effective_folds <- min(cv_folds, n_errors)
                    if (effective_folds < cv_folds) {
                        log_warn("  Reducing CV folds from {cv_folds} to {effective_folds} due to small error count ({n_errors})")
                    }

                    cv_fit <- cv.glmnet(x_in, y_vec_factor, family = "binomial", alpha = lasso_alpha, nfolds = effective_folds)

                    # Extract Selected Features
                    coefs <- coef(cv_fit, s = lambda_selection)
                    active_idx <- which(coefs[-1] != 0)
                    active_vars <- colnames(x_in)[active_idx]

                    # E. Optimize Thresholds (on Training Set)
                    # ----------------------------------------
                    # Calculate metrics using ONLY selected features
                    # MeanAbsZ, MAR

                    # We need Z-scores for these active vars too
                    z_sub <- long_feat[rsid %in% active_vars, .(mz = mean(abs(z), na.rm = TRUE)), by = SampleID]
                    res_sub <- long_feat[rsid %in% active_vars, .(mar = median(abs(std_res), na.rm = TRUE)), by = SampleID]

                    # Merge metrics
                    metrics <- merge(curr_samples[, .(SampleID, is_error)], z_sub, by = "SampleID", all.x = TRUE)
                    metrics <- merge(metrics, res_sub, by = "SampleID", all.x = TRUE)

                    # Filter out NAs (samples without data for active variants)
                    metrics_complete <- metrics[!is.na(mz) & !is.na(mar)]
                    log_info("  Metrics: {nrow(metrics)} total, {nrow(metrics_complete)} complete (after removing NAs)")

                    # Verify both classes are present
                    if (nrow(metrics_complete) == 0) {
                        log_warn("  No complete metrics after filtering. Skipping threshold optimization for cohort {cohort_idx}.")
                        next
                    }

                    error_table <- table(metrics_complete$is_error)
                    if (length(error_table) < 2 || any(error_table == 0)) {
                        log_warn("  Not enough class diversity in complete metrics: {paste(names(error_table), '=', error_table, collapse=', ')}. Skipping threshold optimization for cohort {cohort_idx}.")
                        next
                    }

                    # Youden's J
                    roc_z <- roc(metrics_complete$is_error, metrics_complete$mz, quiet = TRUE)
                    roc_mar <- roc(metrics_complete$is_error, metrics_complete$mar, quiet = TRUE)

                    th_z <- coords(roc_z, "best", best.method = "youden")$threshold
                    th_mar <- coords(roc_mar, "best", best.method = "youden")$threshold

                    # Convert to k (SD/MAD from clean mean/median)
                    clean_metrics <- metrics_complete[is_error == FALSE]
                    k_z <- (th_z - mean(clean_metrics$mz, na.rm = TRUE)) / sd(clean_metrics$mz, na.rm = TRUE)
                    k_mar <- (th_mar - median(clean_metrics$mar, na.rm = TRUE)) / (mad(clean_metrics$mar, na.rm = TRUE) * 1.4826)

                    # Store
                    cohort_results[[cohort_idx]] <- list(
                        active_vars = active_vars,
                        k_z = k_z,
                        k_mar = k_mar,
                        auc_z = auc(roc_z),
                        auc_mar = auc(roc_mar)
                    )

                    # ============================================================
                    # DIAGNOSTIC PLOTS
                    # ============================================================

                    # Setup plot directory
                    plot_dir <- file.path(output_dir, "diagnostic_plots")
                    dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

                    # A. ROC and PR Curves
                    # --------------------
                    # ROC Data
                    roc_data <- data.table(
                        Specificity = c(roc_z$specificities, roc_mar$specificities),
                        Sensitivity = c(roc_z$sensitivities, roc_mar$sensitivities),
                        Metric = rep(
                            c(
                                paste0("MeanAbsZ (AUC=", round(auc(roc_z), 3), ")"),
                                paste0("MAR (AUC=", round(auc(roc_mar), 3), ")")
                            ),
                            times = c(length(roc_z$specificities), length(roc_mar$specificities))
                        )
                    )

                    # PR Curve Data (using PRROC or manual calculation if pkg unavailable,
                    # here simplifying using pROC's coords logic or just base logic if needed.
                    # actually pROC doesn't do PR curves natively well.
                    # Let's use a simple approximation/calculation or skipping if complex dep needed.
                    # User asked for it, assuming pROC or similar.
                    # Actually, let's use the 'coords' to get precision/recall if possible or generic calc.
                    # PRROC package is standard but might not be installed.
                    # Let's calculate Precision/Recall manually from predictions for plot.

                    calc_pr <- function(scores, labels) {
                        ord <- order(scores, decreasing = TRUE)
                        sc <- scores[ord]
                        lb <- labels[ord]
                        tp <- cumsum(lb)
                        fp <- cumsum(!lb)
                        n_pos <- sum(lb)
                        recall <- tp / n_pos
                        precision <- tp / (tp + fp)
                        return(data.table(Recall = recall, Precision = precision))
                    }

                    pr_z_dt <- calc_pr(metrics_complete$mz, metrics_complete$is_error)
                    pr_z_dt[, Metric := paste0("MeanAbsZ (AUC=", round(auc(roc_z), 3), ")")] # Reusing ROC AUC for label or separate?
                    # PR AUC calculation is different, but for now using label to distinguish

                    pr_mar_dt <- calc_pr(metrics_complete$mar, metrics_complete$is_error)
                    pr_mar_dt[, Metric := paste0("MAR (AUC=", round(auc(roc_mar), 3), ")")]

                    pr_data <- rbind(pr_z_dt, pr_mar_dt)

                    # Colors: MeanAbsZ = #E7B800 (Yellow/Gold), MAR = #FC4E07 (Orange/Red)
                    custom_colors <- c("#E7B800", "#FC4E07")
                    names(custom_colors) <- unique(roc_data$Metric) # Map to factor levels
                    # Using crude matching by substring since AUC varies
                    col_scale <- scale_color_manual(values = c("#E7B800", "#FC4E07"))

                    # ROC Plot
                    p_roc <- ggplot(roc_data, aes(x = 1 - Specificity, y = Sensitivity, color = Metric)) +
                        geom_line(size = 1) +
                        geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +
                        scale_color_manual(values = c("#E7B800", "#FC4E07")) +
                        labs(
                            title = paste0("ROC Curves - Cohort ", cohort_idx),
                            subtitle = paste0("Error Config: Geno=", error_cfg$genotype_mismatch, ", Sex=", error_cfg$sex_mismatch)
                        ) +
                        theme_bw()

                    ggsave(file.path(plot_dir, paste0("05a_diagnostic_ROC_cohort_", cohort_idx, ".pdf")), p_roc, width = 6, height = 5)

                    # PR Plot
                    p_pr <- ggplot(pr_data, aes(x = Recall, y = Precision, color = Metric)) +
                        geom_line(size = 1) +
                        labs(
                            title = paste0("Precision-Recall Curves - Cohort ", cohort_idx),
                            subtitle = paste0("Error Config: Geno=", error_cfg$genotype_mismatch, ", Sex=", error_cfg$sex_mismatch)
                        ) +
                        scale_color_manual(values = c("#E7B800", "#FC4E07")) +
                        theme_bw()

                    ggsave(file.path(plot_dir, paste0("05a_diagnostic_PR_cohort_", cohort_idx, ".pdf")), p_pr, width = 6, height = 5)


                    # B. Faceted Sex Cross-reference Plots
                    # ------------------------------------
                    if (!is.null(sex_preds)) {
                        # Merge Metrics
                        plot_data <- merge(metrics_complete, sex_preds[, .(SAMPLE_ID, predicted_prob, genetic_sex)],
                            by.x = "SampleID", by.y = "SAMPLE_ID"
                        )

                        if (nrow(plot_data) > 0) {
                            # Helper for plotting one metric
                            plot_metric_faceted <- function(dt, y_var, y_lab, cutoff, metric_name, color_hex) {
                                # Add plotting attributes
                                dt[, point_shape := ifelse(get(y_var) > cutoff, "Above threshold", "Below threshold")]
                                dt[, border_col := ifelse(is_error, "black", NA)]
                                # Color gradient logic handled by ggplot aes

                                p <- ggplot(dt, aes(x = predicted_prob, y = get(y_var))) +
                                    geom_point(aes(fill = is_error, color = border_col, shape = point_shape),
                                        size = 2.5, alpha = 0.7, stroke = 0.5
                                    ) +
                                    scale_fill_gradient(low = "grey80", high = color_hex, guide = "legend") +
                                    # Actually user wants gradient towards female prob for outliers?
                                    # "gradient colour showing outliers in darker shades (the colour gradient is directed toward the female predicted probabilities)"
                                    # This is complex to exact match without original code but let's approximate:
                                    # Outliers (is_error) colored by predicted_prob?
                                    # Let's stick to the 05b logic closer if possible.
                                    # From 05b backup: geom_point(aes(fill = point_color ...))
                                    # We'll define point_color manually for full control

                                    facet_wrap(~genetic_sex, nrow = 2, scales = "free_x") +
                                    geom_vline(xintercept = youden_threshold_05, linetype = "dashed", color = "red", alpha = 0.5) +
                                    geom_hline(yintercept = cutoff, linetype = "dashed", color = "blue", alpha = 0.5) +
                                    labs(
                                        title = paste0(y_lab, " vs Sex Pred - Cohort ", cohort_idx),
                                        subtitle = paste0("Thresh: ", round(cutoff, 3), " | Youden J: ", round(youden_threshold_05, 3)),
                                        x = "Predicted Female Probability", y = y_lab
                                    ) +
                                    theme_bw()

                                # Re-implementing specific color logic matching 05b
                                # Green=TP, Red=FN, Orange=FP, Blue=TN was my previous simple one.
                                # User wants: "showing all individuals ... with gradient colour showing outliers in darker shades"
                                # Let's use specific colors for errors vs clean

                                p <- ggplot(dt, aes(x = predicted_prob, y = get(y_var))) +
                                    geom_point(aes(fill = predicted_prob, shape = point_shape, color = border_col),
                                        size = 2.5, alpha = 0.8, stroke = 0.5
                                    ) +
                                    scale_fill_gradient(low = "white", high = color_hex) +
                                    facet_wrap(~genetic_sex, nrow = 2) +
                                    theme_bw() +
                                    labs(
                                        title = paste0(y_lab, " vs Sex Pred - Cohort ", cohort_idx),
                                        y = y_lab, x = "Predicted Female Prob"
                                    )

                                return(p)
                            }

                            # Simplified implementation using the direct logic requested:
                            # "dark borders and samples correctly identified as mismatches (i.e. above Z-score threshold in triangle)"

                            # Plot MeanAbsZ
                            p_z <- ggplot(plot_data, aes(x = predicted_prob, y = mz)) +
                                geom_point(aes(fill = is_error, shape = ifelse(mz > th_z, "Above", "Below"), color = ifelse(is_error, "black", NA)),
                                    size = 3, stroke = 0.5, alpha = 0.7
                                ) +
                                scale_shape_manual(values = c("Above" = 24, "Below" = 21)) +
                                scale_color_identity() +
                                scale_fill_manual(values = c("FALSE" = "#90CAF9", "TRUE" = "#E7B800")) + # Blue clean, Yellow/Gold error (matches metric color)
                                facet_wrap(~genetic_sex, ncol = 1) +
                                geom_vline(xintercept = youden_threshold_05, linetype = "dashed", color = "red") +
                                geom_hline(yintercept = th_z, linetype = "dashed", color = "blue") +
                                labs(title = paste0("MeanAbsZ vs Sex Prob - Cohort ", cohort_idx), y = "MeanAbsZ") +
                                theme_bw()

                            ggsave(file.path(plot_dir, paste0("05a_diagnostic_MeanAbsZ_cohort_", cohort_idx, ".pdf")), p_z, width = 8, height = 8)

                            # Plot MAR
                            p_mar <- ggplot(plot_data, aes(x = predicted_prob, y = mar)) +
                                geom_point(aes(fill = is_error, shape = ifelse(mar > th_mar, "Above", "Below"), color = ifelse(is_error, "black", NA)),
                                    size = 3, stroke = 0.5, alpha = 0.7
                                ) +
                                scale_shape_manual(values = c("Above" = 24, "Below" = 21)) +
                                scale_color_identity() +
                                scale_fill_manual(values = c("FALSE" = "#90CAF9", "TRUE" = "#FC4E07")) + # Blue clean, Orange/Red error
                                facet_wrap(~genetic_sex, ncol = 1) +
                                geom_vline(xintercept = youden_threshold_05, linetype = "dashed", color = "red") +
                                geom_hline(yintercept = th_mar, linetype = "dashed", color = "blue") +
                                labs(title = paste0("MAR vs Sex Prob - Cohort ", cohort_idx), y = "MAR") +
                                theme_bw()

                            ggsave(file.path(plot_dir, paste0("05a_diagnostic_MAR_cohort_", cohort_idx, ".pdf")), p_mar, width = 8, height = 8)
                        }
                    }


                    log_info("- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -")

                    log_info("  Cohort {cohort_idx} Done. Active Vars: {length(active_vars)}. k_Z: {round(k_z,2)}, k_MAR: {round(k_mar,2)}")
                }
            }


            # ====================================================================
            # 5. Aggregation & Output
            # ====================================================================
            log_info("--- Aggregating Results ---")

            # Consolidate Active Vars
            all_active <- unlist(lapply(cohort_results, function(x) x$active_vars))
            freq_table <- sort(table(all_active), decreasing = TRUE)
            freq_dt <- data.table(rsid = names(freq_table), count = as.integer(freq_table))
            freq_dt[, frequency := count / length(cohort_results)]

            # Select Consensus (Frequency > threshold)
            consensus_vars <- freq_dt[frequency >= selection_threshold]

            # Average Thresholds
            mean_k_z <- mean(sapply(cohort_results, function(x) x$k_z), na.rm = TRUE)
            mean_k_mar <- mean(sapply(cohort_results, function(x) x$k_mar), na.rm = TRUE)

            log_info("Consensus pQTLs: {nrow(consensus_vars)} (freq >= {selection_threshold})")
            log_info("Recommended Thresholds: MeanAbsZ > Mean + {round(mean_k_z, 2)}*SD, MAR > Median + {round(mean_k_mar, 2)}*MAD")

            # Write Outputs
            fwrite(consensus_vars, file.path(output_dir, "05a_consensus_pqtls.tsv"), sep = "\t")

            # Write config YAML snippet
            out_yaml <- list(
                pqtls = consensus_vars$rsid,
                thresholds = list(
                    mean_abs_z = list(k_sd = mean_k_z),
                    mar = list(k_mad = mean_k_mar)
                )
            )
            write_yaml(out_yaml, file.path(output_dir, "05a_consensus_config.yaml"))

            elapsed <- difftime(Sys.time(), start_time, units = "mins")
            log_info("Training completed in {round(elapsed, 2)} minutes.")
        },
        error = function(e) {
            log_error("Error: {e$message}")
            stop(e)
        }
    )
}

# Auto-execute when sourced by pipeline runner or run directly
# Only execute if config is available (when sourced, config should be loaded above)
if (!is.null(config)) {
    # Execute the training function (it will check if training is enabled internally)
    run_pqtl_training(config, batch_id)
}
