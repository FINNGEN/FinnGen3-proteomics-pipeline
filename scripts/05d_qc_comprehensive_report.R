#!/usr/bin/env Rscript

#################################################
# Script: 05d_qc_comprehensive_report.R
# Author: Reza Jabal, PhD (rjabal@broadinstitute.org)
# Description: Generate comprehensive QC reports after all outlier detection steps
#              Refactored version - integrates outliers from Steps 01-05b
#              Creates integrated outlier tracking, annotated metadata, and clean matrices
# Date: December 2025
#################################################

suppressPackageStartupMessages({
    library(data.table)
    library(tidyverse)
    library(yaml)
    library(logger)
    library(arrow)  # For parquet file support
})

# Source path utilities
script_dir <- tryCatch({
    cmd_args <- commandArgs(trailingOnly = FALSE)
    file_arg <- grep("^--file=", cmd_args, value = TRUE)
    if (length(file_arg) > 0) {
        script_path <- sub("^--file=", "", file_arg)
        dirname(normalizePath(script_path))
    } else {
        getwd()
    }
}, error = function(e) getwd())
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
log_info("Starting comprehensive QC report generation for batch: {batch_id}")

# Function to identify control probes
identify_control_probes <- function(protein_names) {
    control_patterns <- c(
        "Incubation control",
        "Extension control",
        "Amplification control"
    )

    control_probes <- character(0)
    for (pattern in control_patterns) {
        matches <- grep(pattern, protein_names, ignore.case = TRUE, value = TRUE)
        control_probes <- c(control_probes, matches)
    }

    return(control_probes)
}

# Helper function to convert matrix to data.table with SampleID
matrix_to_dt <- function(mat) {
    dt <- as.data.table(mat)
    dt[, SampleID := rownames(mat)]
    setcolorder(dt, c("SampleID", setdiff(names(dt), "SampleID")))
    return(dt)
}

# Main function
main <- function() {
    # ========================================================================
    # PHASE 1: DATA COLLECTION
    # ========================================================================

    log_info("=" |> rep(70) |> paste(collapse = ""))
    log_info("PHASE 1: LOADING BASE DATA AND SAMPLE MAPPINGS")
    log_info("=" |> rep(70) |> paste(collapse = ""))

    # 1.1: Load metadata (AG samples and blank samples already removed in step 00)
    metadata_path <- get_output_path("00", "metadata", batch_id, "qc", config = config)
    if (!file.exists(metadata_path)) {
        stop("Metadata file not found: {metadata_path}. Run step 00 first.")
    }
    metadata_full <- readRDS(metadata_path)
    log_info("Loaded metadata: {nrow(metadata_full)} samples")
    log_info("  Note: AG samples (53) and blank/control samples (20) were already excluded in step 00")

    # 1.2: Load sample mapping (for FINNGENID lookup)
    mapping_path <- get_output_path("00", "sample_mapping", batch_id, "qc", "tsv", config = config)
    if (!file.exists(mapping_path)) {
        stop("Sample mapping not found: {mapping_path}. Run step 00 first.")
    }
    sample_mapping <- fread(mapping_path)
    log_info("Loaded sample mapping: {nrow(sample_mapping)} samples")

    # 1.3: Load base NPX matrix (step 00 - analysis ready, includes bridging samples)
    base_npx_path <- get_output_path("00", "npx_matrix_analysis_ready", batch_id, "qc", config = config)
    if (!file.exists(base_npx_path)) {
        stop("Base NPX matrix not found: {base_npx_path}. Run step 00 first.")
    }
    base_npx_matrix <- readRDS(base_npx_path)
    log_info("Loaded base NPX matrix: {nrow(base_npx_matrix)} samples x {ncol(base_npx_matrix)} proteins")
    log_info("  Note: Blank/control samples (20) were removed in step 00; only biological samples included")

    # ========================================================================
    # PHASE 2: COLLECT OUTLIER LISTS FROM ALL QC STEPS
    # ========================================================================

    log_info("")
    log_info("=" |> rep(70) |> paste(collapse = ""))
    log_info("PHASE 2: COLLECTING OUTLIER LISTS FROM ALL QC STEPS")
    log_info("=" |> rep(70) |> paste(collapse = ""))

    # Load AG samples (excluded samples) to filter them out from step 00 failures
    excluded_samples_path <- get_output_path("00", "excluded_samples", batch_id, "qc", config = config)
    ag_samples <- character(0)
    if (file.exists(excluded_samples_path)) {
        excluded_dt <- readRDS(excluded_samples_path)
        # Handle different column names (SAMPLE_ID or SampleID)
        id_col_excluded <- if ("SAMPLE_ID" %in% names(excluded_dt)) "SAMPLE_ID" else "SampleID"
        ag_samples <- excluded_dt[[id_col_excluded]]
        log_info("Loaded {length(ag_samples)} AG samples (to be excluded from step 00 failures)")
    } else {
        log_warn("Excluded samples file not found: {excluded_samples_path}. AG samples may not be filtered correctly.")
    }

    # Note: AG samples are already removed in pre-filtered input, so no initial QC failures to track
    step00_failures <- character(0)
    log_info("step 00: No initial QC failures (AG samples already removed in pre-filtered input)")

    # step 01: PCA outliers
    pca_path <- get_output_path("01", "pca_outliers_by_source", batch_id, "outliers", "tsv", config = config)
    if (file.exists(pca_path)) {
        pca_dt <- fread(pca_path)
        step01_outliers <- pca_dt[Any == 1]$SampleID
        log_info("step 01 (PCA): {length(step01_outliers)} outliers")
    } else {
        step01_outliers <- character(0)
        log_warn("step 01 PCA outliers not found: {pca_path}")
    }

    # step 02: Technical outliers
    tech_path <- get_output_path("02", "technical_outlier_summary", batch_id, "outliers", "tsv", config = config)
    if (file.exists(tech_path)) {
        tech_dt <- fread(tech_path)
        # Multiple ways to identify technical outliers
        tech_cols <- names(tech_dt)
        is_outlier <- rep(FALSE, nrow(tech_dt))
        if ("Tech_Any" %in% tech_cols) is_outlier <- is_outlier | (tech_dt$Tech_Any == 1)
        if ("is_sample_outlier" %in% tech_cols) is_outlier <- is_outlier | (tech_dt$is_sample_outlier == TRUE)
        if ("is_plate_outlier" %in% tech_cols) is_outlier <- is_outlier | (tech_dt$is_plate_outlier == TRUE)
        if ("is_batch_outlier" %in% tech_cols) is_outlier <- is_outlier | (tech_dt$is_batch_outlier == TRUE)
        step02_outliers <- tech_dt[is_outlier]$SampleID
        log_info("step 02 (Technical): {length(step02_outliers)} outliers")
    } else {
        step02_outliers <- character(0)
        log_warn("step 02 technical outliers not found: {tech_path}")
    }

    # step 03: Z-score outliers
    zscore_path <- get_output_path("03", "zscore_outliers_list", batch_id, "outliers", "tsv", config = config)
    if (file.exists(zscore_path)) {
        zscore_dt <- fread(zscore_path)
        step03_outliers <- zscore_dt$SampleID
        log_info("step 03 (Z-score): {length(step03_outliers)} outliers")
    } else {
        step03_outliers <- character(0)
        log_warn("step 03 Z-score outliers not found: {zscore_path}")
    }

    # step 04: Sex mismatches (STRICT - predicted_sex != genetic_sex)
    sex_mismatch_path <- get_output_path("04", "sex_mismatches", batch_id, "outliers", "tsv", config = config)
    if (file.exists(sex_mismatch_path)) {
        sex_mismatch_dt <- fread(sex_mismatch_path)
        id_col <- if ("SAMPLE_ID" %in% names(sex_mismatch_dt)) "SAMPLE_ID" else "SampleID"
        step04_mismatches <- sex_mismatch_dt[[id_col]]
        log_info("step 04 (Sex Mismatch - Strict): {length(step04_mismatches)} samples")
    } else {
        step04_mismatches <- character(0)
        log_warn("step 04 sex mismatches not found: {sex_mismatch_path}")
    }

    # step 04: Sex outliers (threshold-based, not strict mismatch)
    sex_outlier_path <- get_output_path("04", "sex_outliers", batch_id, "outliers", "tsv", config = config)
    if (file.exists(sex_outlier_path)) {
        sex_outlier_dt <- fread(sex_outlier_path)
        id_col <- if ("SAMPLE_ID" %in% names(sex_outlier_dt)) "SAMPLE_ID" else "SampleID"
        step04_outliers <- sex_outlier_dt[[id_col]]
        log_info("step 04 (Sex Outlier - Threshold): {length(step04_outliers)} samples")
    } else {
        step04_outliers <- character(0)
        log_warn("step 04 sex outliers not found: {sex_outlier_path}")
    }

    # step 05b: pQTL outliers
    pqtl_path <- get_output_path("05b", "pqtl_outliers", batch_id, "outliers", "tsv", config = config)
    if (file.exists(pqtl_path)) {
        pqtl_dt <- fread(pqtl_path)
        # Can be either SampleID or FINNGENID
        id_col <- if ("SampleID" %in% names(pqtl_dt)) {
            "SampleID"
        } else if ("FINNGENID" %in% names(pqtl_dt)) {
            "FINNGENID"
        } else {
            NULL
        }

        if (!is.null(id_col)) {
            step05_outliers_raw <- pqtl_dt[[id_col]]
            # If using FINNGENID, convert to SampleID for consistency
            if (id_col == "FINNGENID") {
                step05_outliers <- sample_mapping[FINNGENID %in% step05_outliers_raw]$SampleID
                log_info("step 05b (pQTL): {length(step05_outliers_raw)} outliers (converted from FINNGENID to SampleID)")
            } else {
                step05_outliers <- step05_outliers_raw
                log_info("step 05b (pQTL): {length(step05_outliers)} outliers")
            }
        } else {
            step05_outliers <- character(0)
            log_warn("step 05b pQTL outliers missing identifier column (expected SampleID or FINNGENID)")
        }
    } else {
        step05_outliers <- character(0)
        log_warn("step 05b pQTL outliers not found: {pqtl_path}")
    }

    # ========================================================================
    # PHASE 3: BUILD COMPREHENSIVE TRACKING TABLE
    # ========================================================================

    log_info("")
    log_info("=" |> rep(70) |> paste(collapse = ""))
    log_info("PHASE 3: BUILDING COMPREHENSIVE QC TRACKING TABLE")
    log_info("=" |> rep(70) |> paste(collapse = ""))

    # 3.1: Create master sample list from base NPX matrix
    all_samples <- data.table(SampleID = rownames(base_npx_matrix))
    log_info("Starting with {nrow(all_samples)} samples from base NPX matrix")

    # Add FINNGENID and other columns from sample mapping (includes bridging samples)
    # Sample mapping from step 00 has: SampleID, FINNGENID, COHORT_FINNGENID, BIOBANK_PLASMA, sample_type
    mapping_cols_to_merge <- intersect(c("SampleID", "FINNGENID", "COHORT_FINNGENID", "BIOBANK_PLASMA", "sample_type"),
                                       names(sample_mapping))
    if (length(mapping_cols_to_merge) > 0) {
        all_samples <- merge(all_samples, sample_mapping[, mapping_cols_to_merge, with = FALSE],
                             by = "SampleID", all.x = TRUE)
    }
    n_with_finngenid <- sum(!is.na(all_samples$FINNGENID))
    n_bridging <- sum(all_samples$sample_type == "Bridging", na.rm = TRUE)
    log_info("Mapped FINNGENIDs: {n_with_finngenid}/{nrow(all_samples)} samples")
    log_info("Bridging samples identified: {n_bridging}")

    # Add metadata columns (BIOBANK_PLASMA, DISEASE_GROUP) from step 00 metadata
    # Note: Bridging samples may not be in step 00 metadata, so we use sample_mapping as primary source
    id_col_meta <- if ("SAMPLE_ID" %in% names(metadata_full)) "SAMPLE_ID" else "SampleID"
    metadata_subset <- if ("BIOBANK_PLASMA" %in% names(metadata_full) && "DISEASE_GROUP" %in% names(metadata_full)) {
        metadata_full[, c(id_col_meta, "BIOBANK_PLASMA", "DISEASE_GROUP"), with = FALSE]
    } else if ("BIOBANK_PLASMA" %in% names(metadata_full)) {
        metadata_full[, c(id_col_meta, "BIOBANK_PLASMA"), with = FALSE]
    } else {
        NULL
    }

    if (!is.null(metadata_subset)) {
        # Merge metadata, but preserve values from sample_mapping for bridging samples
        all_samples <- merge(all_samples, metadata_subset,
                             by.x = "SampleID", by.y = id_col_meta, all.x = TRUE, suffixes = c("", "_meta"))
        # Use metadata values only if sample_mapping values are missing
        if ("BIOBANK_PLASMA_meta" %in% names(all_samples)) {
            all_samples[is.na(BIOBANK_PLASMA) & !is.na(BIOBANK_PLASMA_meta),
                       BIOBANK_PLASMA := BIOBANK_PLASMA_meta]
            all_samples[, BIOBANK_PLASMA_meta := NULL]
        }
        if ("DISEASE_GROUP" %in% names(metadata_subset)) {
            if ("DISEASE_GROUP_meta" %in% names(all_samples)) {
                all_samples[is.na(DISEASE_GROUP) & !is.na(DISEASE_GROUP_meta),
                           DISEASE_GROUP := DISEASE_GROUP_meta]
                all_samples[, DISEASE_GROUP_meta := NULL]
            }
        }
        log_info("Added metadata columns from original metadata file")
    }

    # Ensure all samples from base NPX matrix are included (including bridging samples)
    # This is already handled by starting with base_npx_matrix row names, but verify
    n_samples_in_matrix <- nrow(base_npx_matrix)
    n_samples_in_tracking <- nrow(all_samples)
    if (n_samples_in_matrix != n_samples_in_tracking) {
        log_warn("Sample count mismatch: NPX matrix has {n_samples_in_matrix} samples, tracking table has {n_samples_in_tracking}")
    } else {
        log_info("Sample count verified: {n_samples_in_tracking} samples in tracking table matches NPX matrix")
    }

    # 3.2: Note: No step 00 failures (AG samples already removed in pre-filtered input)

    # 3.3: Create binary QC flag columns
    log_info("Creating QC flag columns for each step...")
    all_samples[, QC_pca := as.integer(SampleID %in% step01_outliers)]
    all_samples[, QC_technical := as.integer(SampleID %in% step02_outliers)]
    all_samples[, QC_zscore := as.integer(SampleID %in% step03_outliers)]
    all_samples[, QC_sex_mismatch := as.integer(SampleID %in% step04_mismatches)]
    all_samples[, QC_sex_outlier := as.integer(SampleID %in% step04_outliers)]
    all_samples[, QC_pqtl := as.integer(SampleID %in% step05_outliers)]

    # 3.4: Create summary columns
    log_info("Creating summary columns (QC_flag, N_Methods, Detection_Steps)...")

    # QC_flag: 1 if flagged by ANY step
    all_samples[, QC_flag := as.integer(
        QC_pca == 1 | QC_technical == 1 | QC_zscore == 1 |
        QC_sex_mismatch == 1 | QC_sex_outlier == 1 | QC_pqtl == 1
    )]

    # N_Methods: Count of flagging methods
    all_samples[, N_Methods := QC_pca + QC_technical + QC_zscore +
                               QC_sex_mismatch + QC_sex_outlier + QC_pqtl]

    # Detection_Steps: Comma-separated list
    all_samples[, Detection_Steps := {
        steps <- character(0)
        if (QC_pca == 1) steps <- c(steps, "PCA")
        if (QC_technical == 1) steps <- c(steps, "Technical")
        if (QC_zscore == 1) steps <- c(steps, "Zscore")
        if (QC_sex_mismatch == 1) steps <- c(steps, "SexMismatch")
        if (QC_sex_outlier == 1) steps <- c(steps, "SexOutlier")
        if (QC_pqtl == 1) steps <- c(steps, "pQTL")
        if (length(steps) == 0) NA_character_ else paste(steps, collapse = ",")
    }, by = SampleID]

    log_info("Comprehensive tracking table complete: {nrow(all_samples)} samples")

    # ========================================================================
    # PHASE 3.5: LOAD AND MERGE RAW QC METRICS
    # ========================================================================

    log_info("")
    log_info("=" |> rep(70) |> paste(collapse = ""))
    log_info("PHASE 3.5: LOADING RAW QC METRICS FROM ALL STEPS")
    log_info("=" |> rep(70) |> paste(collapse = ""))

    # step 00: Missing rate (if available)
    # Note: Step 00 may not output a QC summary file, so this is optional
    qc_summary_path <- get_output_path("00", "qc_summary", batch_id, "qc", "tsv", config = config)
    if (file.exists(qc_summary_path)) {
        log_info("Loading step 00 metrics (missing_rate)...")
        qc_summary <- fread(qc_summary_path)
        step00_metrics <- qc_summary[, .(SampleID, QC_initial_qc_missing_rate = missing_rate)]
        all_samples <- merge(all_samples, step00_metrics, by = "SampleID", all.x = TRUE)
        log_info("  Merged missing_rate for {sum(!is.na(all_samples$QC_initial_qc_missing_rate))} samples")
    } else {
        log_warn("step 00 QC summary not found - skipping missing_rate metrics")
        all_samples[, QC_initial_qc_missing_rate := NA_real_]
    }

    # step 01: PCA scores (PC1-PC4 for all samples)
    pca_result_path <- get_output_path("01", "pca_result", batch_id, "outliers", config = config)
    if (file.exists(pca_result_path)) {
        log_info("Loading step 01 metrics (PC1-PC4 scores)...")
        pca_result <- readRDS(pca_result_path)
        if (!is.null(pca_result$scores) && nrow(pca_result$scores) > 0) {
            pca_scores <- as.data.table(pca_result$scores)
            pca_scores[, SampleID := rownames(pca_result$scores)]
            # Extract PC1-PC4
            pca_cols <- c("PC1", "PC2", "PC3", "PC4")
            available_pcs <- intersect(pca_cols, names(pca_scores))
            if (length(available_pcs) > 0) {
                pca_metrics <- pca_scores[, c("SampleID", available_pcs), with = FALSE]
                setnames(pca_metrics, available_pcs, paste0("QC_pca_", tolower(available_pcs)))
                all_samples <- merge(all_samples, pca_metrics, by = "SampleID", all.x = TRUE)
                log_info("  Merged PCA scores for {sum(!is.na(all_samples$QC_pca_pc1))} samples")
            } else {
                log_warn("  No PC1-PC4 columns found in PCA scores")
                all_samples[, QC_pca_pc1 := NA_real_]
                all_samples[, QC_pca_pc2 := NA_real_]
                all_samples[, QC_pca_pc3 := NA_real_]
                all_samples[, QC_pca_pc4 := NA_real_]
            }
        } else {
            log_warn("  PCA result file exists but scores are missing")
            all_samples[, QC_pca_pc1 := NA_real_]
            all_samples[, QC_pca_pc2 := NA_real_]
            all_samples[, QC_pca_pc3 := NA_real_]
            all_samples[, QC_pca_pc4 := NA_real_]
        }
    } else {
        log_warn("step 01 PCA result not found - skipping PCA score metrics")
        all_samples[, QC_pca_pc1 := NA_real_]
        all_samples[, QC_pca_pc2 := NA_real_]
        all_samples[, QC_pca_pc3 := NA_real_]
        all_samples[, QC_pca_pc4 := NA_real_]
    }

    # step 04: Sex prediction metrics (shared for both mismatch and outlier)
    if (file.exists(sex_mismatch_path) || file.exists(sex_outlier_path)) {
        log_info("Loading step 04 metrics (predicted_prob, predicted_sex, genetic_sex)...")
        sex_predictions_path <- get_output_path("04", "sex_predictions", batch_id, "outliers", "tsv", config = config)
        if (file.exists(sex_predictions_path)) {
            sex_preds <- fread(sex_predictions_path)
            id_col_sex <- if ("SAMPLE_ID" %in% names(sex_preds)) "SAMPLE_ID" else "SampleID"
            sex_metrics <- sex_preds[, .(
                SampleID = get(id_col_sex),
                QC_sex_predicted_prob = predicted_prob,
                QC_sex_predicted_sex = predicted_sex,
                QC_sex_genetic_sex = genetic_sex
            )]
            all_samples <- merge(all_samples, sex_metrics, by = "SampleID", all.x = TRUE)
            log_info("  Merged sex prediction metrics for {sum(!is.na(all_samples$QC_sex_predicted_prob))} samples")
        } else {
            log_warn("  step 04 sex predictions file not found - skipping sex metrics")
            all_samples[, QC_sex_predicted_prob := NA_real_]
            all_samples[, QC_sex_predicted_sex := NA_character_]
            all_samples[, QC_sex_genetic_sex := NA_character_]
        }
    } else {
        log_warn("step 04 sex files not found - skipping sex metrics")
        all_samples[, QC_sex_predicted_prob := NA_real_]
        all_samples[, QC_sex_predicted_sex := NA_character_]
        all_samples[, QC_sex_genetic_sex := NA_character_]
    }

    # step 02: Technical outlier metrics (all four metrics for all samples)
    tech_stats_path <- get_output_path("02", "sample_technical_stats", batch_id, "outliers", "tsv", config = config)
    if (file.exists(tech_stats_path)) {
        log_info("Loading step 02 metrics (mean_npx, sd_npx, missing_rate, qc_fail_rate)...")
        tech_stats <- fread(tech_stats_path)
        tech_metrics <- tech_stats[, .(
            SampleID,
            QC_technical_mean_npx = mean_npx,
            QC_technical_sd_npx = sd_npx,
            QC_technical_missing_rate = missing_rate,
            QC_technical_qc_fail_rate = qc_fail_rate
        )]
        all_samples <- merge(all_samples, tech_metrics, by = "SampleID", all.x = TRUE)
        log_info("  Merged technical metrics for {sum(!is.na(all_samples$QC_technical_mean_npx))} samples")
    } else {
        log_warn("step 02 technical stats not found - skipping technical metrics")
        all_samples[, QC_technical_mean_npx := NA_real_]
        all_samples[, QC_technical_sd_npx := NA_real_]
        all_samples[, QC_technical_missing_rate := NA_real_]
        all_samples[, QC_technical_qc_fail_rate := NA_real_]
    }

    # step 03: Z-score outlier metrics
    zscore_summary_path <- get_output_path("03", "zscore_outlier_summary", batch_id, "outliers", "tsv", config = config)
    if (file.exists(zscore_summary_path)) {
        log_info("Loading step 03 metrics (n_outlier_proteins, pct_outlier_proteins, max_abs_zscore)...")
        zscore_summary <- fread(zscore_summary_path)
        zscore_metrics <- zscore_summary[, .(
            SampleID,
            QC_zscore_n_outlier_proteins = n_outlier_proteins,
            QC_zscore_pct_outlier_proteins = pct_outlier_proteins,
            QC_zscore_max_abs_zscore = max_abs_zscore
        )]
        all_samples <- merge(all_samples, zscore_metrics, by = "SampleID", all.x = TRUE)
        log_info("  Merged Z-score metrics for {sum(!is.na(all_samples$QC_zscore_n_outlier_proteins))} samples")
    } else {
        log_warn("step 03 Z-score summary not found - skipping Z-score metrics")
        all_samples[, QC_zscore_n_outlier_proteins := NA_integer_]
        all_samples[, QC_zscore_pct_outlier_proteins := NA_real_]
        all_samples[, QC_zscore_max_abs_zscore := NA_real_]
    }

    # step 05b: pQTL outlier metrics
    pqtl_stats_path <- get_output_path("05b", "pqtl_stats", batch_id, "outliers", "tsv", config = config)
    if (file.exists(pqtl_stats_path)) {
        log_info("Loading step 05b metrics (mean_abs_z, max_abs_z, median_abs_residual, n_prots)...")
        pqtl_stats <- fread(pqtl_stats_path)
        # Convert SampleID if needed (pQTL might use FINNGENID)
        id_col_pqtl <- if ("SampleID" %in% names(pqtl_stats)) {
            "SampleID"
        } else if ("FINNGENID" %in% names(pqtl_stats)) {
            # Convert FINNGENID to SampleID
            pqtl_stats <- merge(pqtl_stats, sample_mapping[, .(SampleID, FINNGENID)],
                               by = "FINNGENID", all.x = TRUE)
            "SampleID"
        } else {
            NULL
        }

        if (!is.null(id_col_pqtl) && "SampleID" %in% names(pqtl_stats)) {
            pqtl_metrics <- pqtl_stats[, .(
                SampleID,
                QC_pqtl_mean_abs_z = MeanAbsZ,
                QC_pqtl_max_abs_z = MaxAbsZ,
                QC_pqtl_median_abs_residual = MedianAbsResidual,
                QC_pqtl_n_prots = N_Prots
            )]
            all_samples <- merge(all_samples, pqtl_metrics, by = "SampleID", all.x = TRUE)
            log_info("  Merged pQTL metrics for {sum(!is.na(all_samples$QC_pqtl_mean_abs_z))} samples")
        } else {
            log_warn("  pQTL stats file exists but SampleID column not found")
            all_samples[, QC_pqtl_mean_abs_z := NA_real_]
            all_samples[, QC_pqtl_max_abs_z := NA_real_]
            all_samples[, QC_pqtl_median_abs_residual := NA_real_]
            all_samples[, QC_pqtl_n_prots := NA_integer_]
        }
    } else {
        log_warn("step 05b pQTL stats not found - adding columns with NA (will be populated when pQTL step runs)")
        all_samples[, QC_pqtl_mean_abs_z := NA_real_]
        all_samples[, QC_pqtl_max_abs_z := NA_real_]
        all_samples[, QC_pqtl_median_abs_residual := NA_real_]
        all_samples[, QC_pqtl_n_prots := NA_integer_]
    }

    log_info("Raw QC metrics loading complete")

    # ========================================================================
    # PHASE 4: GENERATE OUTPUTS
    # ========================================================================

    log_info("")
    log_info("=" |> rep(70) |> paste(collapse = ""))
    log_info("PHASE 4: GENERATING OUTPUT FILES")
    log_info("=" |> rep(70) |> paste(collapse = ""))

    # 4.1: Output 1 - Comprehensive Outlier List (flagged samples only)
    log_info("Generating Output 1: Comprehensive Outlier List...")

    outlier_list <- all_samples[QC_flag == 1]

    # Note: No step 00 failures (AG samples already removed in pre-filtered input)

    # Reorder columns: base columns, then each QC flag followed by its metrics
    col_order <- c("SampleID", "FINNGENID", "BIOBANK_PLASMA", "DISEASE_GROUP",
                   # step 01: PCA
                   "QC_pca", "QC_pca_pc1", "QC_pca_pc2", "QC_pca_pc3", "QC_pca_pc4",
                   # step 02: Technical
                   "QC_technical", "QC_technical_mean_npx", "QC_technical_sd_npx",
                   "QC_technical_missing_rate", "QC_technical_qc_fail_rate",
                   # step 03: Z-score
                   "QC_zscore", "QC_zscore_n_outlier_proteins", "QC_zscore_pct_outlier_proteins",
                   "QC_zscore_max_abs_zscore",
                   # step 04: Sex (shared metrics after mismatch)
                   "QC_sex_mismatch", "QC_sex_predicted_prob", "QC_sex_predicted_sex", "QC_sex_genetic_sex",
                   "QC_sex_outlier",
                   # step 05b: pQTL
                   "QC_pqtl", "QC_pqtl_mean_abs_z", "QC_pqtl_max_abs_z",
                   "QC_pqtl_median_abs_residual", "QC_pqtl_n_prots",
                   # Summary columns
                   "QC_flag", "N_Methods", "Detection_Steps")
    # Only keep columns that exist
    col_order <- intersect(col_order, names(outlier_list))
    outlier_list <- outlier_list[, ..col_order]

    # Save using path_utils naming convention
    output1_path <- get_output_path("05d", "comprehensive_outliers_list",
                                    batch_id, "phenotypes", "tsv", config = config)
    output1_parquet_path <- get_output_path("05d", "comprehensive_outliers_list",
                                           batch_id, "phenotypes", "parquet", config = config)
    ensure_output_dir(output1_path)
    ensure_output_dir(output1_parquet_path)
    fwrite(outlier_list, output1_path, sep = "\t")
    write_parquet(outlier_list, output1_parquet_path)
    log_info("Saved comprehensive outlier list:")
    log_info("  TSV: {output1_path}")
    log_info("  Parquet: {output1_parquet_path}")
    log_info("  Outliers: {nrow(outlier_list)} unique samples")

    # 4.2: Output 2 - Annotated Metadata (all samples with QC flags)
    log_info("")
    log_info("Generating Output 2: Annotated Metadata...")

    # Start with all_samples (includes bridging samples from base NPX matrix)
    # This ensures all 2,527 samples (2,477 FinnGen + 50 bridging) are included
    annotated_metadata <- copy(all_samples)

    # Rename SampleID to match metadata column name if needed
    if (id_col_meta != "SampleID" && "SampleID" %in% names(annotated_metadata)) {
        setnames(annotated_metadata, "SampleID", id_col_meta)
    }

    # Merge with full metadata to get all original metadata columns
    # Use all.x = TRUE to keep bridging samples that might not be in metadata_full
    if (nrow(metadata_full) > 0) {
        # Get all columns from metadata_full except those already in annotated_metadata
        meta_cols_to_add <- setdiff(names(metadata_full), names(annotated_metadata))
        if (length(meta_cols_to_add) > 0) {
            metadata_to_merge <- metadata_full[, c(id_col_meta, meta_cols_to_add), with = FALSE]
            annotated_metadata <- merge(annotated_metadata, metadata_to_merge,
                                        by = id_col_meta, all.x = TRUE)
            log_info("Merged additional metadata columns: {length(meta_cols_to_add)} columns")
        }
    }

    # Fill NAs with 0 for binary flags (samples not flagged)
    binary_qc_cols <- c("QC_initial_qc", "QC_pca", "QC_sex_mismatch", "QC_sex_outlier",
                        "QC_technical", "QC_zscore", "QC_pqtl", "QC_flag")
    binary_qc_cols <- intersect(binary_qc_cols, names(annotated_metadata))
    annotated_metadata[, (binary_qc_cols) := lapply(.SD, function(x) fifelse(is.na(x), 0L, x)),
                       .SDcols = binary_qc_cols]

    # Reorder columns: original metadata, then QC flags with metrics, then summary
    original_meta_cols <- names(metadata_full)
    qc_flag_cols <- c("QC_initial_qc", "QC_initial_qc_missing_rate",
                      "QC_pca", "QC_pca_pc1", "QC_pca_pc2", "QC_pca_pc3", "QC_pca_pc4",
                      "QC_sex_mismatch", "QC_sex_predicted_prob", "QC_sex_predicted_sex", "QC_sex_genetic_sex",
                      "QC_sex_outlier",
                      "QC_technical", "QC_technical_mean_npx", "QC_technical_sd_npx",
                      "QC_technical_missing_rate", "QC_technical_qc_fail_rate",
                      "QC_zscore", "QC_zscore_n_outlier_proteins", "QC_zscore_pct_outlier_proteins",
                      "QC_zscore_max_abs_zscore",
                      "QC_pqtl", "QC_pqtl_mean_abs_z", "QC_pqtl_max_abs_z",
                      "QC_pqtl_median_abs_residual", "QC_pqtl_n_prots",
                      "QC_flag", "N_Methods", "Detection_Steps")
    # Only keep columns that exist
    qc_flag_cols <- intersect(qc_flag_cols, names(annotated_metadata))
    final_col_order <- c(original_meta_cols, qc_flag_cols)
    # Add any remaining columns that weren't in the specified order
    remaining_cols <- setdiff(names(annotated_metadata), final_col_order)
    final_col_order <- c(final_col_order, remaining_cols)
    # Only keep columns that actually exist
    final_col_order <- intersect(final_col_order, names(annotated_metadata))
    annotated_metadata <- annotated_metadata[, ..final_col_order]

    # Save using path_utils naming convention
    output2_path <- get_output_path("05d", "qc_annotated_metadata",
                                    batch_id, "phenotypes", "tsv", config = config)
    output2_parquet_path <- get_output_path("05d", "qc_annotated_metadata",
                                           batch_id, "phenotypes", "parquet", config = config)
    ensure_output_dir(output2_path)
    ensure_output_dir(output2_parquet_path)
    fwrite(annotated_metadata, output2_path, sep = "\t")
    write_parquet(annotated_metadata, output2_parquet_path)
    log_info("Saved annotated metadata:")
    log_info("  TSV: {output2_path}")
    log_info("  Parquet: {output2_parquet_path}")
    log_info("  Total samples: {nrow(annotated_metadata)}")
    log_info("  Samples with QC flags: {sum(annotated_metadata$QC_flag == 1, na.rm = TRUE)}")

    # Verify sample count matches NPX matrix (should be 2,527: 2,477 FinnGen + 50 bridging)
    expected_total_samples <- 2527
    if (nrow(annotated_metadata) == expected_total_samples) {
        log_info("  ✓ Sample count verified: {nrow(annotated_metadata)} samples (matches NPX matrix)")
    } else {
        log_warn("  Sample count mismatch: {nrow(annotated_metadata)} samples, expected {expected_total_samples}")
    }

    # Count bridging samples
    n_bridging_in_meta <- sum(annotated_metadata$sample_type == "Bridging", na.rm = TRUE)
    if (n_bridging_in_meta == 50) {
        log_info("  ✓ Bridging samples verified: {n_bridging_in_meta} bridging samples included")
    } else {
        log_warn("  Bridging sample count mismatch: {n_bridging_in_meta} found, expected 50")
    }

    # 4.3: Output 3 - Clean NPX Matrix (all QC failures removed)
    log_info("")
    log_info("Generating Output 3: Clean NPX Matrix...")

    # Get list of all flagged samples
    flagged_samples <- all_samples[QC_flag == 1]$SampleID

    # Remove flagged samples from base NPX matrix
    clean_samples <- setdiff(rownames(base_npx_matrix), flagged_samples)
    clean_npx_matrix <- base_npx_matrix[clean_samples, ]

    # Save using path_utils naming convention
    # Output 3a: Clean NPX Matrix with all proteins (including control probes)
    output3_path <- get_output_path("05d", "npx_matrix_all_qc_passed",
                                    batch_id, "phenotypes", "rds", config = config)
    ensure_output_dir(output3_path)
    saveRDS(clean_npx_matrix, output3_path)
    log_info("Saved clean NPX matrix (all proteins): {output3_path}")
    log_info("  Clean samples: {nrow(clean_npx_matrix)} ({ncol(clean_npx_matrix)} proteins)")
    log_info("  Removed samples: {length(flagged_samples)}")
    pct_removed <- round(100 * length(flagged_samples) / nrow(base_npx_matrix), 2)
    log_info("  Removal rate: {pct_removed}%")

    # 4.4: Output 3b - Biological-only NPX Matrix (control probes removed)
    # Check if control probe removal is enabled (default: true)
    remove_control_probes <- tryCatch(
        isTRUE(config$parameters$qc$remove_control_probes),
        error = function(e) TRUE  # Default to true if not specified
    )

    # Initialize output paths (will be set if control probe removal is enabled)
    output3b_rds_path <- NULL
    output3b_parquet_path <- NULL
    output3b_tsv_path <- NULL

    if (remove_control_probes) {
        log_info("")
        log_info("Generating Output 3b: Biological-only NPX Matrix (control probes removed)...")

        # Identify control probes
        protein_names <- colnames(clean_npx_matrix)
        control_probes <- identify_control_probes(protein_names)

        if (length(control_probes) > 0) {
            log_info("  Control probes identified: {length(control_probes)}")
            incubation <- sum(grepl("Incubation control", control_probes, ignore.case = TRUE))
            extension <- sum(grepl("Extension control", control_probes, ignore.case = TRUE))
            amplification <- sum(grepl("Amplification control", control_probes, ignore.case = TRUE))
            log_info("    - Incubation controls: {incubation}")
            log_info("    - Extension controls: {extension}")
            log_info("    - Amplification controls: {amplification}")

            # Remove control probes
            biological_proteins <- setdiff(protein_names, control_probes)
            biological_npx_matrix <- clean_npx_matrix[, biological_proteins, drop = FALSE]
            log_info("  Biological proteins: {length(biological_proteins)}")
            log_info("  Control probes removed: {length(control_probes)}")
        } else {
            log_warn("  No control probes found in matrix - all proteins will be retained")
            biological_npx_matrix <- clean_npx_matrix
        }

        # Save biological-only RDS
        output3b_rds_path <- get_output_path("05d", "npx_matrix_all_qc_passed_biological_only",
                                            batch_id, "phenotypes", "rds", config = config)
        ensure_output_dir(output3b_rds_path)
        saveRDS(biological_npx_matrix, output3b_rds_path)
        log_info("  Saved biological-only RDS: {output3b_rds_path}")
        log_info("    Dimensions: {nrow(biological_npx_matrix)} samples × {ncol(biological_npx_matrix)} proteins")

        # Convert to data.table for Parquet/TSV export
        biological_dt <- matrix_to_dt(biological_npx_matrix)

        # Save biological-only Parquet
        output3b_parquet_path <- get_output_path("05d", "npx_matrix_all_qc_passed_biological_only",
                                                batch_id, "phenotypes", "parquet", config = config)
        ensure_output_dir(output3b_parquet_path)
        write_parquet(biological_dt, output3b_parquet_path)
        log_info("  Saved biological-only Parquet: {output3b_parquet_path}")

        # Save biological-only TSV
        output3b_tsv_path <- get_output_path("05d", "npx_matrix_all_qc_passed_biological_only",
                                            batch_id, "phenotypes", "tsv", config = config)
        ensure_output_dir(output3b_tsv_path)
        fwrite(biological_dt, output3b_tsv_path, sep = "\t")
        log_info("  Saved biological-only TSV: {output3b_tsv_path}")

        log_info("  ✓ Biological-only outputs created (RDS, Parquet, TSV)")
    } else {
        log_info("")
        log_info("Control probe removal is disabled (parameters.qc.remove_control_probes: false)")
        log_info("  Skipping biological-only output generation")
    }

    # ========================================================================
    # PHASE 5: ENHANCED LOGGING & SUMMARY
    # ========================================================================

    log_info("")
    log_info("=" |> rep(70) |> paste(collapse = ""))
    log_info("COMPREHENSIVE QC SUMMARY REPORT")
    log_info("=" |> rep(70) |> paste(collapse = ""))
    log_info("")
    log_info("Batch: {batch_id}")
    log_info("Total samples tracked: {nrow(all_samples)}")
    log_info("  Note: Blank/control samples (20) were removed in step 00 and are not tracked here")
    log_info("")

    # Unique samples flagged
    n_flagged_unique <- sum(all_samples$QC_flag == 1)
    n_passed <- sum(all_samples$QC_flag == 0)
    pct_flagged <- round(100 * n_flagged_unique / nrow(all_samples), 2)

    log_info("QC Results:")
    log_info("  Unique samples flagged: {n_flagged_unique} ({pct_flagged}%)")
    log_info("  Samples passing all QC: {n_passed} ({round(100 - pct_flagged, 2)}%)")
    log_info("")

    # Per-step breakdown
    log_info("Samples flagged by each QC step:")
    log_info("  step 00 (Initial QC):   {sum(all_samples$QC_initial_qc == 1)} samples")
    log_info("  step 01 (PCA):          {sum(all_samples$QC_pca == 1)} samples")
    log_info("  step 04 (Sex Mismatch): {sum(all_samples$QC_sex_mismatch == 1)} samples")
    log_info("  step 04 (Sex Outlier):  {sum(all_samples$QC_sex_outlier == 1)} samples")
    log_info("  step 02 (Technical):    {sum(all_samples$QC_technical == 1)} samples")
    log_info("  step 03 (Z-score):      {sum(all_samples$QC_zscore == 1)} samples")
    log_info("  step 05b (pQTL):       {sum(all_samples$QC_pqtl == 1)} samples")
    log_info("")

    # Breakdown by number of methods
    log_info("Samples flagged by multiple methods:")
    method_counts <- all_samples[QC_flag == 1, .N, by = N_Methods][order(N_Methods)]
    if (nrow(method_counts) > 0) {
        for (i in seq_len(nrow(method_counts))) {
            n_methods <- method_counts[i]$N_Methods
            n_samples <- method_counts[i]$N
            pct <- round(100 * n_samples / n_flagged_unique, 1)
            log_info("  {n_methods} method(s): {n_samples} samples ({pct}%)")
        }
    } else {
        log_info("  No samples flagged by multiple methods")
    }
    log_info("")

    log_info("Output files generated:")
    log_info("  1. Comprehensive outlier list:")
    log_info("     TSV:     {output1_path}")
    log_info("     Parquet: {output1_parquet_path}")
    log_info("  2. Annotated metadata:")
    log_info("     TSV:     {output2_path}")
    log_info("     Parquet: {output2_parquet_path}")
    log_info("  3. Clean NPX matrix (all proteins): {output3_path}")
    if (remove_control_probes && !is.null(output3b_rds_path)) {
        log_info("  3b. Biological-only NPX matrix (control probes removed):")
        log_info("       RDS:     {output3b_rds_path}")
        log_info("       Parquet: {output3b_parquet_path}")
        log_info("       TSV:     {output3b_tsv_path}")
    }
    log_info("")
    log_info("=" |> rep(70) |> paste(collapse = ""))
    log_info("QC comprehensive report generation completed successfully")
    log_info("=" |> rep(70) |> paste(collapse = ""))
}

# Run main function
if (!interactive()) {
    main()
}

