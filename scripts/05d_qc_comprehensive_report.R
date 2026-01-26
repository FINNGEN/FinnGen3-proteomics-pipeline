#!/usr/bin/env Rscript
# ==============================================================================
# 05d_qc_comprehensive_report.R - Comprehensive QC Report Generation
# ==============================================================================
#
# Purpose:
#   Generates comprehensive QC reports after all outlier detection steps. Integrates
#   outliers from all QC methods (Initial QC, PCA, Technical, Z-score, Sex, pQTL)
#   using union logic. Creates integrated outlier tracking visualisations, annotated
#   metadata with QC flags, and clean matrices with all QC-passed samples.
#
# Author: Reza Jabal, PhD (rjabal@broadinstitute.org)
# Date: December 2025 (Updated: January 2026 - Initial QC integration and visualisation)
# ==============================================================================

suppressPackageStartupMessages({
    library(data.table)
    library(tidyverse)
    library(yaml)
    library(logger)
    library(arrow)  # For parquet file support
    library(gridExtra)  # For combining plots
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

# Function to create comprehensive outlier tracking visualisation
# Similar to the plot generated in step 03, but includes all QC steps including initial QC
create_comprehensive_outlier_tracking_plot <- function(all_samples) {
    log_info("Creating comprehensive outlier tracking visualisation...")

    # Ensure QC flag columns exist and are integers (handle NAs)
    qc_cols <- c("QC_initial_qc", "QC_pca", "QC_sex_mismatch", "QC_sex_outlier",
                 "QC_technical", "QC_zscore", "QC_pqtl")
    qc_cols <- intersect(qc_cols, names(all_samples))

    # Convert to integer, treating NA as 0
    for (col in qc_cols) {
        if (col %in% names(all_samples)) {
            all_samples[, (col) := as.integer(fifelse(is.na(get(col)), 0L, get(col)))]
        }
    }

    # 1. Outlier method overlap barplot (separate sex mismatch and sex outlier, include initial QC and pQTL)
    # Safely get counts for each method (handle missing columns)
    get_count <- function(col_name) {
        if (col_name %in% names(all_samples)) {
            sum(all_samples[[col_name]] == 1, na.rm = TRUE)
        } else {
            0L
        }
    }

    method_counts <- data.table(
        Method = c("Initial QC", "PCA", "Sex Mismatch", "Sex Outlier", "Technical", "Z-score", "pQTL"),
        Count = c(
            get_count("QC_initial_qc"),
            get_count("QC_pca"),
            get_count("QC_sex_mismatch"),
            get_count("QC_sex_outlier"),
            get_count("QC_technical"),
            get_count("QC_zscore"),
            get_count("QC_pqtl")
        ),
        Type = c("QC", "Primary", "Sex", "Sex", "Primary", "Primary", "Primary")
    )

    # Filter out methods with 0 counts for cleaner visualisation
    method_counts <- method_counts[Count > 0]

    p1 <- ggplot(method_counts, aes(x = reorder(Method, Count), y = Count, fill = Type)) +
        geom_bar(stat = "identity") +
        geom_text(aes(label = Count), hjust = -0.2, size = 4) +
        scale_fill_manual(values = c("Primary" = "#3A5F8A", "Sex" = "#E74C3C", "QC" = "#F39C12")) +
        coord_flip() +
        labs(title = "Outliers Detected by Each Method",
             subtitle = "All QC Steps Integration (Steps 00-05b)",
             x = "Detection Method", y = "Number of Samples") +
        theme_bw() +
        theme(legend.position = "none",
              plot.title = element_text(size = 14, face = "bold"),
              plot.subtitle = element_text(size = 11))

    # 2. Multi-method detection distribution
    # CRITICAL: Ensure all QC columns exist before using them
    # Create missing columns with default value 0
    required_qc_cols <- c("QC_initial_qc", "QC_pca", "QC_sex_mismatch", "QC_sex_outlier",
                         "QC_technical", "QC_zscore", "QC_pqtl")
    for (col in required_qc_cols) {
        if (!col %in% names(all_samples)) {
            all_samples[, (col) := 0L]
            log_warn("Column {col} not found - initialized to 0")
        }
    }

    # Ensure N_Methods exists and handle NAs
    if (!"N_Methods" %in% names(all_samples)) {
        all_samples[, N_Methods := QC_initial_qc + QC_pca + QC_technical + QC_zscore +
                                   QC_sex_mismatch + QC_sex_outlier + QC_pqtl]
    }
    all_samples[, N_Methods := as.integer(fifelse(is.na(N_Methods), 0L, N_Methods))]

    # Ensure QC_flag exists
    if (!"QC_flag" %in% names(all_samples)) {
        all_samples[, QC_flag := as.integer(
            QC_initial_qc == 1 | QC_pca == 1 | QC_technical == 1 | QC_zscore == 1 |
            QC_sex_mismatch == 1 | QC_sex_outlier == 1 | QC_pqtl == 1
        )]
    }

    n_methods_dist <- all_samples[QC_flag == 1, .N, by = N_Methods][order(N_Methods)]
    if (nrow(n_methods_dist) == 0) {
        # If no outliers, create empty plot
        n_methods_dist <- data.table(N_Methods = 0L, N = 0L)
    }

    p2 <- ggplot(n_methods_dist, aes(x = factor(N_Methods), y = N, fill = factor(N_Methods))) +
        geom_bar(stat = "identity") +
        geom_text(aes(label = N), vjust = -0.5, size = 4) +
        scale_fill_brewer(palette = "YlOrRd") +
        labs(title = "Samples Flagged by Multiple Methods",
             subtitle = "Higher counts = higher confidence outliers",
             x = "Number of Methods Detecting Sample", y = "Number of Samples") +
        theme_bw() +
        theme(legend.position = "none",
              plot.title = element_text(size = 14, face = "bold"),
              plot.subtitle = element_text(size = 11))

    # 3. Overlap patterns (Venn diagram style) - including all methods
    # All QC columns should exist at this point (created/initialized above)
    # Use direct column access - columns are guaranteed to exist after initialization
    overlap_summary <- all_samples[, .(
        InitialQC_only = sum(QC_initial_qc == 1 & QC_pca == 0 & QC_sex_mismatch == 0 &
                            QC_sex_outlier == 0 & QC_technical == 0 & QC_zscore == 0 & QC_pqtl == 0, na.rm = TRUE),
        PCA_only = sum(QC_initial_qc == 0 & QC_pca == 1 & QC_sex_mismatch == 0 &
                      QC_sex_outlier == 0 & QC_technical == 0 & QC_zscore == 0 & QC_pqtl == 0, na.rm = TRUE),
        SexMismatch_only = sum(QC_initial_qc == 0 & QC_pca == 0 & QC_sex_mismatch == 1 &
                               QC_sex_outlier == 0 & QC_technical == 0 & QC_zscore == 0 & QC_pqtl == 0, na.rm = TRUE),
        SexOutlier_only = sum(QC_initial_qc == 0 & QC_pca == 0 & QC_sex_mismatch == 0 &
                             QC_sex_outlier == 1 & QC_technical == 0 & QC_zscore == 0 & QC_pqtl == 0, na.rm = TRUE),
        Technical_only = sum(QC_initial_qc == 0 & QC_pca == 0 & QC_sex_mismatch == 0 &
                           QC_sex_outlier == 0 & QC_technical == 1 & QC_zscore == 0 & QC_pqtl == 0, na.rm = TRUE),
        Zscore_only = sum(QC_initial_qc == 0 & QC_pca == 0 & QC_sex_mismatch == 0 &
                         QC_sex_outlier == 0 & QC_technical == 0 & QC_zscore == 1 & QC_pqtl == 0, na.rm = TRUE),
        pQTL_only = sum(QC_initial_qc == 0 & QC_pca == 0 & QC_sex_mismatch == 0 &
                       QC_sex_outlier == 0 & QC_technical == 0 & QC_zscore == 0 & QC_pqtl == 1, na.rm = TRUE),
        Two_Methods = sum(N_Methods == 2, na.rm = TRUE),
        Three_Methods = sum(N_Methods == 3, na.rm = TRUE),
        Four_Plus_Methods = sum(N_Methods >= 4, na.rm = TRUE)
    )]

    overlap_long <- data.table(
        Category = names(overlap_summary),
        Count = as.numeric(overlap_summary[1,])
    )
    overlap_long <- overlap_long[Count > 0]

    if (nrow(overlap_long) > 0) {
        p3 <- ggplot(overlap_long, aes(x = reorder(Category, Count), y = Count)) +
            geom_bar(stat = "identity", fill = "#5C9EAD") +
            geom_text(aes(label = Count), hjust = -0.2, size = 3.5) +
            coord_flip() +
            labs(title = "Outlier Detection Overlap Patterns",
                 x = "Overlap Category", y = "Number of Samples") +
            theme_bw() +
            theme(plot.title = element_text(size = 14, face = "bold"),
                  axis.text.y = element_text(size = 8))
    } else {
        # Create empty plot if no overlaps
        p3 <- ggplot(data.frame(x = 1, y = 1), aes(x = x, y = y)) +
            geom_blank() +
            labs(title = "Outlier Detection Overlap Patterns",
                 subtitle = "No overlap patterns detected",
                 x = "Overlap Category", y = "Number of Samples") +
            theme_bw() +
            theme(plot.title = element_text(size = 14, face = "bold"))
    }

    # Combine plots
    combined <- grid.arrange(p1, p2, p3, ncol = 2, nrow = 2)

    log_info("Comprehensive outlier tracking visualisation created")

    return(list(
        method_counts = p1,
        n_methods_dist = p2,
        overlap_patterns = p3,
        combined = combined
    ))
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

    # step 00: Initial QC failures (samples with >10% missing data)
    # Load from QC summary file created in Step 00
    qc_summary_path <- get_output_path("00", "qc_summary", batch_id, "qc", "tsv", config = config)
    if (file.exists(qc_summary_path)) {
        log_info("Loading step 00 initial QC failures from QC summary...")
        qc_summary <- fread(qc_summary_path)
        step00_failures <- qc_summary[QC_initial_qc == 1]$SampleID
        log_info("  Found {length(step00_failures)} samples failing initial QC (>10% missing data)")
        if(length(step00_failures) > 0) {
            log_info("  Failed samples: {paste(step00_failures, collapse=', ')}")
        }
    } else {
        step00_failures <- character(0)
        log_warn("step 00 QC summary not found - assuming no initial QC failures")
    }

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
    # Try multiple possible file naming patterns (step_num can be "05" or "05b")
    pqtl_path <- get_output_path("05b", "05b_pqtl_outliers", batch_id, "outliers", "tsv", config = config)
    if (!file.exists(pqtl_path)) {
        # Fallback 1: Try with step prefix using "05b" (backward compatibility)
        pqtl_path <- get_output_path("05b", "pqtl_outliers", batch_id, "outliers", "tsv", config = config)
    }
    if (!file.exists(pqtl_path)) {
        # Fallback 2: Try with step "05" and old format (backward compatibility)
        pqtl_path <- get_output_path("05", "05b_pqtl_outliers", batch_id, "outliers", "tsv", config = config)
    }
    if (!file.exists(pqtl_path)) {
        # Fallback 3: Try with step "05b" and old format (backward compatibility)
        pqtl_path <- get_output_path("05b", "pqtl_outliers", batch_id, "outliers", "tsv", config = config)
    }
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

    # 3.2: Add initial QC failures to tracking table (even if not in analysis-ready matrix)
    # CRITICAL FIX: Initial QC failures are removed from matrix before tracking table is created,
    # so we must add them explicitly to ensure all evaluated samples are tracked
    if (length(step00_failures) > 0) {
        log_info("Adding {length(step00_failures)} initial QC failures to tracking table (not in analysis-ready matrix)")

        # Check which initial QC failures are already in tracking table
        missing_from_tracking <- setdiff(step00_failures, all_samples$SampleID)

        if (length(missing_from_tracking) > 0) {
            log_info("  Adding {length(missing_from_tracking)} initial QC failures to tracking table")

            # Create data.table for missing initial QC failures
            # Load QC summary to get metadata for these samples (qc_summary from Phase 2 may not be in scope)
            qc_summary_path <- get_output_path("00", "qc_summary", batch_id, "qc", "tsv", config = config)
            if (file.exists(qc_summary_path)) {
                qc_summary_all <- fread(qc_summary_path)
                initial_qc_rows <- qc_summary_all[SampleID %in% missing_from_tracking]

                # Create minimal tracking table rows for initial QC failures
                initial_qc_tracking <- data.table(
                    SampleID = missing_from_tracking,
                    FINNGENID = NA_character_,
                    COHORT_FINNGENID = NA_character_,
                    BIOBANK_PLASMA = NA_character_,
                    sample_type = NA_character_,
                    DISEASE_GROUP = NA_character_
                )

                # Try to get FINNGENID from sample mapping if available
                if (nrow(sample_mapping) > 0 && "FINNGENID" %in% names(sample_mapping)) {
                    mapping_subset <- sample_mapping[SampleID %in% missing_from_tracking,
                                                     .(SampleID, FINNGENID, COHORT_FINNGENID, BIOBANK_PLASMA, sample_type)]
                    if (nrow(mapping_subset) > 0) {
                        initial_qc_tracking <- merge(initial_qc_tracking, mapping_subset,
                                                     by = "SampleID", all.x = TRUE, suffixes = c("", "_map"))
                        # Use mapping values where available
                        if ("FINNGENID_map" %in% names(initial_qc_tracking)) {
                            initial_qc_tracking[!is.na(FINNGENID_map), FINNGENID := FINNGENID_map]
                            initial_qc_tracking[, FINNGENID_map := NULL]
                        }
                        if ("COHORT_FINNGENID_map" %in% names(initial_qc_tracking)) {
                            initial_qc_tracking[!is.na(COHORT_FINNGENID_map), COHORT_FINNGENID := COHORT_FINNGENID_map]
                            initial_qc_tracking[, COHORT_FINNGENID_map := NULL]
                        }
                        if ("BIOBANK_PLASMA_map" %in% names(initial_qc_tracking)) {
                            initial_qc_tracking[!is.na(BIOBANK_PLASMA_map), BIOBANK_PLASMA := BIOBANK_PLASMA_map]
                            initial_qc_tracking[, BIOBANK_PLASMA_map := NULL]
                        }
                        if ("sample_type_map" %in% names(initial_qc_tracking)) {
                            initial_qc_tracking[!is.na(sample_type_map), sample_type := sample_type_map]
                            initial_qc_tracking[, sample_type_map := NULL]
                        }
                    }
                }

                # Add to tracking table
                all_samples <- rbindlist(list(all_samples, initial_qc_tracking), fill = TRUE)
                log_info("  Added {length(missing_from_tracking)} initial QC failures to tracking table")
                log_info("  Tracking table now has {nrow(all_samples)} samples (includes {length(missing_from_tracking)} initial QC failures)")
            } else {
                log_warn("  QC summary file not found - cannot add initial QC failures with metadata")
            }
        } else {
            log_info("  All initial QC failures already in tracking table")
        }
    }

    # 3.3: Create binary QC flag columns
    log_info("Creating QC flag columns for each step...")
    # CRITICAL FIX: Include initial QC flag (from step00_failures loaded earlier)
    all_samples[, QC_initial_qc := as.integer(SampleID %in% step00_failures)]
    all_samples[, QC_pca := as.integer(SampleID %in% step01_outliers)]
    all_samples[, QC_technical := as.integer(SampleID %in% step02_outliers)]
    all_samples[, QC_zscore := as.integer(SampleID %in% step03_outliers)]
    all_samples[, QC_sex_mismatch := as.integer(SampleID %in% step04_mismatches)]
    all_samples[, QC_sex_outlier := as.integer(SampleID %in% step04_outliers)]
    all_samples[, QC_pqtl := as.integer(SampleID %in% step05_outliers)]

    # 3.4: Create summary columns
    log_info("Creating summary columns (QC_flag, N_Methods, Detection_Steps)...")

    # QC_flag: 1 if flagged by ANY step (including initial QC)
    all_samples[, QC_flag := as.integer(
        QC_initial_qc == 1 | QC_pca == 1 | QC_technical == 1 | QC_zscore == 1 |
        QC_sex_mismatch == 1 | QC_sex_outlier == 1 | QC_pqtl == 1
    )]

    # N_Methods: Count of flagging methods (including initial QC)
    all_samples[, N_Methods := QC_initial_qc + QC_pca + QC_technical + QC_zscore +
                               QC_sex_mismatch + QC_sex_outlier + QC_pqtl]

    # Detection_Steps: Comma-separated list (including initial QC)
    all_samples[, Detection_Steps := {
        steps <- character(0)
        if (QC_initial_qc == 1) steps <- c(steps, "InitialQC")
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

    # step 00: Initial QC metrics (missing_rate and QC_initial_qc flag)
    # Load from QC summary file created in Step 00
    # CRITICAL: QC_initial_qc flag was already created in Phase 3.3 (line 546)
    # This merge only adds missing_rate metric - preserve existing QC_initial_qc flag
    qc_summary_path <- get_output_path("00", "qc_summary", batch_id, "qc", "tsv", config = config)
    if (file.exists(qc_summary_path)) {
        log_info("Loading step 00 metrics (missing_rate, QC_initial_qc)...")
        qc_summary <- fread(qc_summary_path)
        step00_metrics <- qc_summary[, .(SampleID,
                                         QC_initial_qc_missing_rate = missing_rate)]
        # Only merge missing_rate - QC_initial_qc flag already exists from Phase 3.3
        all_samples <- merge(all_samples, step00_metrics, by = "SampleID", all.x = TRUE)
        log_info("  Merged initial QC metrics for {sum(!is.na(all_samples$QC_initial_qc_missing_rate))} samples")
        log_info("  Samples flagged by initial QC: {sum(all_samples$QC_initial_qc == 1, na.rm = TRUE)}")
    } else {
        log_warn("step 00 QC summary not found - skipping initial QC metrics")
        all_samples[, QC_initial_qc_missing_rate := NA_real_]
        # QC_initial_qc flag already exists from Phase 3.3 - don't overwrite
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
    # PHASE 3.6: GENERATE OUTLIER TRACKING VISUALISATION
    # ========================================================================

    log_info("")
    log_info("=" |> rep(70) |> paste(collapse = ""))
    log_info("PHASE 3.6: GENERATING OUTLIER TRACKING VISUALISATION")
    log_info("=" |> rep(70) |> paste(collapse = ""))

    # Create comprehensive outlier tracking plot (similar to step 03)
    outlier_tracking_plots <- create_comprehensive_outlier_tracking_plot(all_samples)

    # Save individual plots
    method_counts_path <- get_output_path(step_num, "05d_outlier_method_counts", batch_id, "outliers", "pdf", config = config)
    n_methods_dist_path <- get_output_path(step_num, "05d_outlier_multi_method_dist", batch_id, "outliers", "pdf", config = config)
    overlap_patterns_path <- get_output_path(step_num, "05d_outlier_overlap_patterns", batch_id, "outliers", "pdf", config = config)
    tracking_combined_path <- get_output_path(step_num, "05d_outlier_tracking_combined", batch_id, "outliers", "pdf", config = config)

    ensure_output_dir(method_counts_path)
    ensure_output_dir(n_methods_dist_path)
    ensure_output_dir(overlap_patterns_path)
    ensure_output_dir(tracking_combined_path)

    ggsave(method_counts_path, outlier_tracking_plots$method_counts, width = 8, height = 6)
    ggsave(n_methods_dist_path, outlier_tracking_plots$n_methods_dist, width = 8, height = 6)
    ggsave(overlap_patterns_path, outlier_tracking_plots$overlap_patterns, width = 10, height = 6)
    ggsave(tracking_combined_path, outlier_tracking_plots$combined, width = 14, height = 10)

    log_info("Saved comprehensive outlier tracking visualisations:")
    log_info("  Method counts: {method_counts_path}")
    log_info("  Multi-method distribution: {n_methods_dist_path}")
    log_info("  Overlap patterns: {overlap_patterns_path}")
    log_info("  Combined plot: {tracking_combined_path}")

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

    # Note: Initial QC failures are now included in tracking table and comprehensive list

    # Reorder columns: base columns, then each QC flag followed by its metrics
    col_order <- c("SampleID", "FINNGENID", "BIOBANK_PLASMA", "DISEASE_GROUP",
                   # step 00: Initial QC
                   "QC_initial_qc", "QC_initial_qc_missing_rate",
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
    output1_path <- get_output_path("05d", "05d_comprehensive_outliers_list",
                                    batch_id, "phenotypes", "tsv", config = config)
    output1_parquet_path <- get_output_path("05d", "05d_comprehensive_outliers_list",
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
    # Note: QC_initial_qc is already set in tracking table creation, but ensure it's included
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
    output2_path <- get_output_path("05d", "05d_qc_annotated_metadata",
                                    batch_id, "phenotypes", "tsv", config = config)
    output2_parquet_path <- get_output_path("05d", "05d_qc_annotated_metadata",
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

    # Verify sample count matches expected (2,527 original - 5 initial QC = 2,522 in matrix + 5 initial QC in tracking = 2,527 total tracked)
    # OR: 2,522 in analysis-ready matrix + 5 initial QC failures = 2,527 total tracked
    expected_total_samples <- 2527  # Original input samples
    expected_in_matrix <- 2522  # After initial QC removal
    if (nrow(annotated_metadata) == expected_total_samples) {
        log_info("  ✓ Sample count verified: {nrow(annotated_metadata)} samples (includes {length(step00_failures)} initial QC failures)")
    } else if (nrow(annotated_metadata) == expected_in_matrix) {
        log_warn("  Sample count: {nrow(annotated_metadata)} samples (missing {length(step00_failures)} initial QC failures in tracking)")
    } else {
        log_warn("  Sample count mismatch: {nrow(annotated_metadata)} samples, expected {expected_total_samples} (or {expected_in_matrix} if initial QC not tracked)")
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

    # Get list of all flagged samples (includes initial QC failures even if not in base matrix)
    flagged_samples <- all_samples[QC_flag == 1]$SampleID
    log_info("  Total flagged samples: {length(flagged_samples)} (includes initial QC failures)")

    # Remove flagged samples from base NPX matrix
    # Note: Initial QC failures are already removed from base_npx_matrix, so they won't be in it
    # This is correct - we just need to track them for reporting purposes
    clean_samples <- setdiff(rownames(base_npx_matrix), flagged_samples)
    clean_npx_matrix <- base_npx_matrix[clean_samples, ]

    # Verify: Initial QC failures should not be in base matrix (they were removed in Step 00)
    initial_qc_in_base <- intersect(step00_failures, rownames(base_npx_matrix))
    if (length(initial_qc_in_base) > 0) {
        log_warn("  WARNING: {length(initial_qc_in_base)} initial QC failures found in base matrix (should be 0)")
    } else {
        log_info("  ✓ Verified: Initial QC failures already removed from base matrix (not present)")
    }

    # Save using path_utils naming convention
    # Output 3a: Clean NPX Matrix with all proteins (including control probes)
    output3_path <- get_output_path("05d", "05d_npx_matrix_all_qc_passed",
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
        output3b_rds_path <- get_output_path("05d", "05d_npx_matrix_all_qc_passed_biological_only",
                                            batch_id, "phenotypes", "rds", config = config)
        ensure_output_dir(output3b_rds_path)
        saveRDS(biological_npx_matrix, output3b_rds_path)
        log_info("  Saved biological-only RDS: {output3b_rds_path}")
        log_info("    Dimensions: {nrow(biological_npx_matrix)} samples × {ncol(biological_npx_matrix)} proteins")

        # Convert to data.table for Parquet/TSV export
        biological_dt <- matrix_to_dt(biological_npx_matrix)

        # Save biological-only Parquet
        output3b_parquet_path <- get_output_path("05d", "05d_npx_matrix_all_qc_passed_biological_only",
                                                batch_id, "phenotypes", "parquet", config = config)
        ensure_output_dir(output3b_parquet_path)
        write_parquet(biological_dt, output3b_parquet_path)
        log_info("  Saved biological-only Parquet: {output3b_parquet_path}")

        # Save biological-only TSV
        output3b_tsv_path <- get_output_path("05d", "05d_npx_matrix_all_qc_passed_biological_only",
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
    log_info("  Note: Initial QC failures ({length(step00_failures)} samples) are included in tracking table for comprehensive reporting")
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

