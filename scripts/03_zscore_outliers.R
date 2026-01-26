#!/usr/bin/env Rscript
# ==============================================================================
# 03_zscore_outliers.R - Z-Score Based Outlier Detection
# ==============================================================================
#
# Purpose:
#   Detects outliers using per-protein Z-score calculation with iterative refinement.
#   Flags samples with >10% of proteins having |Z| > 4. Uses iterative detection
#   to recalculate Z-scores after removing outliers. Operates on the base
#   analysis-ready matrix for parallel flagging architecture.
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
  library(jsonlite)
})

# Source path utilities
# Get script directory first (before sourcing path_utils)
script_dir <- tryCatch({
  env_script <- Sys.getenv("SCRIPT_NAME", "")
  if (env_script != "" && file.exists(env_script)) {
    dirname(normalizePath(env_script))
  } else {
    args <- commandArgs(trailingOnly = FALSE)
    file_arg <- grep("^--file=", args, value = TRUE)
    if (length(file_arg) > 0) {
      script_path <- sub("^--file=", "", file_arg)
      dirname(normalizePath(script_path))
    } else {
      getwd()
    }
  }
}, error = function(e) {
  getwd()
})
source(file.path(script_dir, "path_utils.R"))

# Get config path from environment (required, no default)
config_file <- Sys.getenv("PIPELINE_CONFIG", "")
if (config_file == "" || !file.exists(config_file)) {
  stop("PIPELINE_CONFIG environment variable not set or config file not found. Please provide path to config file.")
}
config <- read_yaml(config_file)

# Get batch context
batch_id <- Sys.getenv("PIPELINE_BATCH_ID", config$batch$default_batch_id %||% "batch_01")
step_num <- get_step_number()

# Set up logging with batch-aware path
log_path <- get_log_path(step_num, batch_id, config)
ensure_output_dir(log_path)
log_appender(appender_file(log_path))
log_info("Starting Z-score based outlier detection for batch: {batch_id}")

# Set theme for plots
theme_set(theme_bw())

# Function to calculate Z-scores for proteins
# Input: NPX matrix (samples x proteins) cleaned from technical outliers
# Method: For each protein (column), standardize values: Z = (X - mean) / SD
# Purpose: Identify samples with extreme protein expression values (|Z| > threshold)
calculate_protein_zscores <- function(npx_matrix) {
  log_info("Calculating Z-scores for each protein (per-protein standardization)")
  log_info("Input matrix: npx_matrix_technical_cleaned.rds from step 02")
  log_info("Formula: Z-score = (NPX_value - protein_mean) / protein_SD")

  # Calculate Z-scores for each protein (column)
  zscore_matrix <- apply(npx_matrix, 2, function(x) {
    # Remove NA values for calculation
    x_clean <- x[!is.na(x)]
    if(length(x_clean) < 3) return(rep(NA, length(x)))

    # Calculate mean and SD for this protein across all samples
    x_mean <- mean(x_clean)
    x_sd <- sd(x_clean)

    # Avoid division by zero
    if(x_sd == 0) return(rep(0, length(x)))

    # Calculate Z-scores: how many SDs each sample deviates from protein mean
    (x - x_mean) / x_sd
  })

  # Ensure matrix format is preserved
  zscore_matrix <- as.matrix(zscore_matrix)
  rownames(zscore_matrix) <- rownames(npx_matrix)

  log_info("Z-score matrix created: {nrow(zscore_matrix)} x {ncol(zscore_matrix)}")

  return(zscore_matrix)
}

# Function to detect outliers based on Z-scores
detect_zscore_outliers <- function(zscore_matrix, threshold = 4) {
  log_info("Detecting outliers with Z-score threshold: {threshold}")

  # Identify outliers (abs(Z) > threshold)
  outlier_matrix <- abs(zscore_matrix) > threshold

  # Count outliers per sample
  sample_outlier_counts <- rowSums(outlier_matrix, na.rm = TRUE)

  # Count outliers per protein
  protein_outlier_counts <- colSums(outlier_matrix, na.rm = TRUE)

  # Identify samples with excessive outliers
  outlier_threshold_pct <- 0.1  # Flag samples with >10% proteins as outliers
  n_proteins <- ncol(zscore_matrix)
  sample_outliers <- names(sample_outlier_counts)[
    sample_outlier_counts > n_proteins * outlier_threshold_pct
  ]

  # Create detailed outlier summary
  outlier_summary <- data.table(
    SampleID = rownames(zscore_matrix),
    n_outlier_proteins = sample_outlier_counts,
    pct_outlier_proteins = sample_outlier_counts / n_proteins * 100,
    max_abs_zscore = apply(abs(zscore_matrix), 1, max, na.rm = TRUE),
    is_outlier = rownames(zscore_matrix) %in% sample_outliers
  )

  log_info("Samples with >10% outlier proteins: {length(sample_outliers)}")

  return(list(
    outlier_matrix = outlier_matrix,
    sample_outlier_counts = sample_outlier_counts,
    protein_outlier_counts = protein_outlier_counts,
    sample_outliers = sample_outliers,
    outlier_summary = outlier_summary
  ))
}

# Function to perform iterative outlier detection
iterative_outlier_detection <- function(npx_matrix, max_iterations = 5, threshold = 4) {
  log_info("Starting iterative outlier detection")

  current_matrix <- npx_matrix
  all_outliers <- character()
  iteration_results <- list()

  for(i in 1:max_iterations) {
    log_info("Iteration {i}")

    # Calculate Z-scores
    zscore_matrix <- calculate_protein_zscores(current_matrix)

    # Detect outliers
    outlier_result <- detect_zscore_outliers(zscore_matrix, threshold)

    if(length(outlier_result$sample_outliers) == 0) {
      log_info("No more outliers detected. Stopping iterations.")
      break
    }

    # Store results
    iteration_results[[paste0("iteration_", i)]] <- outlier_result
    all_outliers <- c(all_outliers, outlier_result$sample_outliers)

    # Remove outliers for next iteration
    current_matrix <- current_matrix[!rownames(current_matrix) %in% outlier_result$sample_outliers, ]

    log_info("Removed {length(outlier_result$sample_outliers)} outliers. Matrix size: {nrow(current_matrix)} x {ncol(current_matrix)}")
  }

  log_info("Iterative detection complete. Total outliers: {length(unique(all_outliers))}")

  return(list(
    all_outliers = unique(all_outliers),
    iteration_results = iteration_results,
    final_matrix = current_matrix
  ))
}

# Function to analyze outlier patterns
analyze_outlier_patterns <- function(zscore_matrix, outlier_result, metadata = NULL) {
  log_info("Analyzing outlier patterns")

  # Get top outlier proteins
  top_outlier_proteins <- names(sort(outlier_result$protein_outlier_counts, decreasing = TRUE)[1:20])

  # Get top outlier samples
  top_outlier_samples <- outlier_result$outlier_summary[order(-n_outlier_proteins)][1:min(20, .N)]

  # If metadata available, check for cohort patterns
  if(!is.null(metadata)) {
    outlier_cohorts <- merge(
      outlier_result$outlier_summary[is_outlier == TRUE],
      metadata[, .(SAMPLE_ID, COHORT_FINNGENID, BIOBANK_PLASMA)],
      by.x = "SampleID",
      by.y = "SAMPLE_ID",
      all.x = TRUE
    )

    # Summarize by cohort
    cohort_summary <- outlier_cohorts[, .(n_outliers = .N), by = COHORT_FINNGENID]
    setorder(cohort_summary, -n_outliers)
  } else {
    cohort_summary <- NULL
  }

  return(list(
    top_outlier_proteins = top_outlier_proteins,
    top_outlier_samples = top_outlier_samples,
    cohort_summary = cohort_summary
  ))
}

# Function to create integrated outlier tracking across all steps (04-05b)
create_integrated_outlier_tracking <- function(zscore_outliers, metadata = NULL, batch_id = NULL, config = NULL) {
  log_info("Creating integrated outlier tracking across steps 04-05b")

  # Get batch_id and config if not provided
  if (is.null(batch_id)) {
    batch_id <- Sys.getenv("PIPELINE_BATCH_ID", "batch_01")
  }
  if (is.null(config)) {
    config_file <- Sys.getenv("PIPELINE_CONFIG",
                              "/mnt/longGWAS_disk_100GB/long_gwas/Github_clones/fg3_Olink_analysis/config/batch2_config.yaml")
    if (file.exists(config_file)) {
      config <- yaml::read_yaml(config_file)
    } else {
      stop("Config file not found: {config_file}")
    }
  }

  # File paths for outlier detection outputs (batch-aware)
  pca_file <- get_output_path("01", "pca_outliers_by_source", batch_id, "outliers", "tsv", config = config)
  sex_file <- get_output_path("04", "sex_mismatches", batch_id, "outliers", "tsv", config = config)
  tech_file <- get_output_path("02", "technical_outlier_summary", batch_id, "outliers", "tsv", config = config)
  pqtl_file <- get_output_path("05b", "05b_pqtl_outliers", batch_id, "outliers", "tsv", config = config)

  # Initialize empty integrated table
  all_samples <- character()

  # Load PCA outliers
  pca_outliers <- NULL
  if (file.exists(pca_file)) {
    pca_outliers <- try(fread(pca_file), silent = TRUE)
    if (!inherits(pca_outliers, "try-error") && nrow(pca_outliers) > 0) {
      all_samples <- unique(c(all_samples, pca_outliers$SampleID))
      log_info("Loaded {nrow(pca_outliers)} PCA outliers")
    }
  }

  # Load sex mismatches (strict - from sex_mismatches.tsv)
  sex_mismatches <- NULL
  if (file.exists(sex_file)) {
    sex_mismatches <- try(fread(sex_file), silent = TRUE)
    if (!inherits(sex_mismatches, "try-error") && nrow(sex_mismatches) > 0) {
      # Rename SAMPLE_ID to SampleID for consistency
      if ("SAMPLE_ID" %in% names(sex_mismatches)) {
        sex_mismatches[, SampleID := SAMPLE_ID]
      }
      all_samples <- unique(c(all_samples, sex_mismatches$SampleID))
      log_info("Loaded {nrow(sex_mismatches)} sex mismatches (strict)")
    }
  }

  # Load sex outliers (relaxed - from sex_predictions.tsv with sex_outlier flag)
  sex_pred_file <- get_output_path("04", "sex_predictions", batch_id, "outliers", "tsv", config = config)
  sex_outliers <- NULL
  if (file.exists(sex_pred_file)) {
    sex_predictions <- try(fread(sex_pred_file), silent = TRUE)
    if (!inherits(sex_predictions, "try-error") && nrow(sex_predictions) > 0) {
      # Filter for sex_outlier == TRUE
      if ("sex_outlier" %in% names(sex_predictions)) {
        # Handle both SAMPLE_ID and SampleID column names
        id_col <- if ("SAMPLE_ID" %in% names(sex_predictions)) "SAMPLE_ID" else "SampleID"
        sex_outliers <- sex_predictions[sex_outlier == TRUE, .(SampleID = get(id_col))]
        all_samples <- unique(c(all_samples, sex_outliers$SampleID))
        log_info("Loaded {nrow(sex_outliers)} sex outliers (relaxed threshold)")
      }
    }
  }

  # Load pQTL outliers
  pqtl_outliers <- NULL
  if (file.exists(pqtl_file)) {
    pqtl_data <- try(fread(pqtl_file), silent = TRUE)
    if (!inherits(pqtl_data, "try-error") && nrow(pqtl_data) > 0) {
      # pQTL outliers file should have SampleID or SAMPLE_ID column
      id_col <- if ("SAMPLE_ID" %in% names(pqtl_data)) "SAMPLE_ID" else "SampleID"
      if (id_col %in% names(pqtl_data)) {
        pqtl_outliers <- pqtl_data[, .(SampleID = get(id_col))]
        all_samples <- unique(c(all_samples, pqtl_outliers$SampleID))
        log_info("Loaded {nrow(pqtl_outliers)} pQTL outliers")
      }
    }
  }

  # Load technical outliers
  tech_outliers <- NULL
  if (file.exists(tech_file)) {
    tech_outliers <- try(fread(tech_file), silent = TRUE)
    if (!inherits(tech_outliers, "try-error") && nrow(tech_outliers) > 0) {
      all_samples <- unique(c(all_samples, tech_outliers$SampleID))
      log_info("Loaded {nrow(tech_outliers)} technical outliers")
    }
  }

  # Add z-score outliers
  all_samples <- unique(c(all_samples, zscore_outliers))

  # Create base table with all outlier samples
  integrated <- data.table(SampleID = all_samples)

  # Add FINNGENID, BIOBANK_PLASMA, DISEASE_GROUP from metadata or PCA file
  # #region agent log
  log_file_path <- "/mnt/longGWAS_disk_100GB/long_gwas/.cursor/debug.log"
  log_debug_data1 <- list(
    sessionId = "debug-session",
    runId = "pre-fix",
    hypothesisId = "B",
    location = "07_zscore_outliers.R:290",
    message = "Before FINNGENID merge - checking sources",
    data = list(
      total_samples = nrow(integrated),
      pca_outliers_exists = !is.null(pca_outliers) && nrow(pca_outliers) > 0,
      pca_outliers_count = if(!is.null(pca_outliers)) nrow(pca_outliers) else 0,
      metadata_exists = !is.null(metadata),
      metadata_cols = if(!is.null(metadata)) names(metadata) else character(0),
      all_samples_count = length(all_samples)
    ),
    timestamp = as.integer(as.numeric(Sys.time()) * 1000)
  )
  cat(jsonlite::toJSON(log_debug_data1, auto_unbox=TRUE), "\n", file=log_file_path, append=TRUE)
  # #endregion

  if (!is.null(pca_outliers) && nrow(pca_outliers) > 0) {
    # Make DISEASE_GROUP optional (may not exist for all batches)
    available_cols <- c("SampleID", "FINNGENID", "BIOBANK_PLASMA")
    if ("DISEASE_GROUP" %in% colnames(pca_outliers)) {
      available_cols <- c(available_cols, "DISEASE_GROUP")
    }
    meta_cols <- pca_outliers[, ..available_cols]
    integrated <- merge(integrated, meta_cols, by = "SampleID", all.x = TRUE)

    # #region agent log
    log_debug_data2 <- list(
      sessionId = "debug-session",
      runId = "pre-fix",
      hypothesisId = "B",
      location = "07_zscore_outliers.R:292",
      message = "After PCA merge - FINNGENID population",
      data = list(
        total_samples = nrow(integrated),
        finngenid_filled = sum(!is.na(integrated$FINNGENID)),
        finngenid_empty = sum(is.na(integrated$FINNGENID)),
        pca_samples_with_fgid = sum(integrated$SampleID %in% pca_outliers$SampleID & !is.na(integrated$FINNGENID)),
        non_pca_samples_with_fgid = sum(!(integrated$SampleID %in% pca_outliers$SampleID) & !is.na(integrated$FINNGENID))
      ),
      timestamp = as.integer(as.numeric(Sys.time()) * 1000)
    )
    cat(jsonlite::toJSON(log_debug_data2, auto_unbox=TRUE), "\n", file=log_file_path, append=TRUE)
    # #endregion
  } else if (!is.null(metadata)) {
    # Try to get from metadata
    if (all(c("SAMPLE_ID", "FINNGENID", "BIOBANK_PLASMA") %in% names(metadata))) {
      meta_cols <- metadata[SAMPLE_ID %in% all_samples, .(SampleID = SAMPLE_ID, FINNGENID, BIOBANK_PLASMA)]
      integrated <- merge(integrated, meta_cols, by = "SampleID", all.x = TRUE)
    }
  }

  # Add PCA outlier flags
  if (!is.null(pca_outliers) && nrow(pca_outliers) > 0) {
    pca_flags <- pca_outliers[, .(
      SampleID,
      PCA_PC1_PC2 = PC1_PC2,
      PCA_PC3_PC4 = PC3_PC4,
      PCA_Median = Median,
      PCA_IQR = IQR,
      PCA_Any = Any
    )]
    integrated <- merge(integrated, pca_flags, by = "SampleID", all.x = TRUE)
    # Fill NAs with 0
    for (col in c("PCA_PC1_PC2", "PCA_PC3_PC4", "PCA_Median", "PCA_IQR", "PCA_Any")) {
      integrated[is.na(get(col)), (col) := 0]
    }
  } else {
    integrated[, `:=` (PCA_PC1_PC2 = 0, PCA_PC3_PC4 = 0, PCA_Median = 0, PCA_IQR = 0, PCA_Any = 0)]
  }

  # Add sex mismatch flag (strict threshold)
  if (!is.null(sex_mismatches) && nrow(sex_mismatches) > 0) {
    integrated[, SexMismatch := as.integer(SampleID %in% sex_mismatches$SampleID)]
  } else {
    integrated[, SexMismatch := 0]
  }

  # Add sex outlier flag (relaxed threshold - more comprehensive)
  if (!is.null(sex_outliers) && nrow(sex_outliers) > 0) {
    integrated[, SexOutlier := as.integer(SampleID %in% sex_outliers$SampleID)]
  } else {
    integrated[, SexOutlier := 0]
  }

  # Sex_Any: flagged by either strict or relaxed
  integrated[, Sex_Any := as.integer(SexMismatch == 1 | SexOutlier == 1)]

  # Add technical outlier flags
  if (!is.null(tech_outliers) && nrow(tech_outliers) > 0) {
    tech_flags <- tech_outliers[, .(
      SampleID,
      Tech_Plate = as.integer(is_plate_outlier),
      Tech_Batch = as.integer(is_batch_outlier),
      Tech_Processing = as.integer(is_processing_outlier),
      Tech_Sample = as.integer(is_sample_outlier),
      Tech_Any = as.integer(is_plate_outlier | is_batch_outlier | is_processing_outlier | is_sample_outlier)
    )]
    integrated <- merge(integrated, tech_flags, by = "SampleID", all.x = TRUE)
    # Fill NAs with 0
    for (col in c("Tech_Plate", "Tech_Batch", "Tech_Processing", "Tech_Sample", "Tech_Any")) {
      integrated[is.na(get(col)), (col) := 0]
    }
  } else {
    integrated[, `:=` (Tech_Plate = 0, Tech_Batch = 0, Tech_Processing = 0, Tech_Sample = 0, Tech_Any = 0)]
  }

  # Add z-score outlier flag
  integrated[, Zscore := as.integer(SampleID %in% zscore_outliers)]

  # Add pQTL outlier flag
  if (!is.null(pqtl_outliers) && nrow(pqtl_outliers) > 0) {
    integrated[, pQTL := as.integer(SampleID %in% pqtl_outliers$SampleID)]
  } else {
    integrated[, pQTL := 0]
  }

  # Calculate cumulative flags (now including pQTL)
  integrated[, Any_Outlier := as.integer(PCA_Any == 1 | Sex_Any == 1 | Tech_Any == 1 | Zscore == 1 | pQTL == 1)]
  integrated[, N_Methods := PCA_Any + Sex_Any + Tech_Any + Zscore + pQTL]

  # Create detection step summary
  integrated[, Detection_Steps := ""]
  integrated[PCA_Any == 1, Detection_Steps := paste0(Detection_Steps, "PCA,")]
  integrated[Sex_Any == 1, Detection_Steps := paste0(Detection_Steps, "Sex,")]
  integrated[Tech_Any == 1, Detection_Steps := paste0(Detection_Steps, "Technical,")]
  integrated[Zscore == 1, Detection_Steps := paste0(Detection_Steps, "Zscore,")]
  integrated[pQTL == 1, Detection_Steps := paste0(Detection_Steps, "pQTL,")]
  # Remove trailing comma
  integrated[, Detection_Steps := gsub(",$", "", Detection_Steps)]

  log_info("Integrated outlier tracking complete: {nrow(integrated)} samples")

  return(integrated)
}

# Function to create integrated outlier tracking visualizations
create_integrated_plots <- function(integrated_tracking) {
  log_info("Creating integrated outlier tracking visualizations")

  library(pheatmap)
  library(gridExtra)

  # 1. Outlier method overlap barplot (separate sex mismatch and sex outlier, include pQTL)
  method_counts <- data.table(
    Method = c("PCA", "Sex Mismatch", "Sex Outlier", "Technical", "Z-score", "pQTL"),
    Count = c(
      sum(integrated_tracking$PCA_Any),
      sum(integrated_tracking$SexMismatch),
      sum(integrated_tracking$SexOutlier),
      sum(integrated_tracking$Tech_Any),
      sum(integrated_tracking$Zscore),
      sum(integrated_tracking$pQTL %||% 0)
    ),
    Type = c("Primary", "Sex", "Sex", "Primary", "Primary", "Primary")
  )

  p1 <- ggplot(method_counts, aes(x = reorder(Method, Count), y = Count, fill = Type)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = Count), hjust = -0.2, size = 4) +
    scale_fill_manual(values = c("Primary" = "#3A5F8A", "Sex" = "#E74C3C")) +
    coord_flip() +
    labs(title = "Outliers Detected by Each Method",
         subtitle = "Steps 01-05b Integration (Sex Mismatch & Sex Outlier shown separately)",
         x = "Detection Method", y = "Number of Samples") +
    theme_bw() +
    theme(legend.position = "none")

  # 2. Multi-method detection distribution
  n_methods_dist <- integrated_tracking[, .N, by = N_Methods][order(N_Methods)]

  p2 <- ggplot(n_methods_dist, aes(x = factor(N_Methods), y = N, fill = factor(N_Methods))) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = N), vjust = -0.5, size = 4) +
    scale_fill_brewer(palette = "YlOrRd") +
    labs(title = "Samples Flagged by Multiple Methods",
         subtitle = "Higher counts = higher confidence outliers",
         x = "Number of Methods Detecting Sample", y = "Number of Samples") +
    theme_bw() +
    theme(legend.position = "none")

  # 3. Overlap patterns (Venn diagram style) - now including pQTL
  # Handle pQTL column (may not exist if step 05b hasn't run yet)
  if (!"pQTL" %in% names(integrated_tracking)) {
    integrated_tracking[, pQTL := 0]
  }

  overlap_summary <- integrated_tracking[, .(
    PCA_only = sum(PCA_Any == 1 & Sex_Any == 0 & Tech_Any == 0 & Zscore == 0 & pQTL == 0),
    Sex_only = sum(PCA_Any == 0 & Sex_Any == 1 & Tech_Any == 0 & Zscore == 0 & pQTL == 0),
    Tech_only = sum(PCA_Any == 0 & Sex_Any == 0 & Tech_Any == 1 & Zscore == 0 & pQTL == 0),
    Zscore_only = sum(PCA_Any == 0 & Sex_Any == 0 & Tech_Any == 0 & Zscore == 1 & pQTL == 0),
    pQTL_only = sum(PCA_Any == 0 & Sex_Any == 0 & Tech_Any == 0 & Zscore == 0 & pQTL == 1),
    PCA_Sex = sum(PCA_Any == 1 & Sex_Any == 1 & Tech_Any == 0 & Zscore == 0 & pQTL == 0),
    PCA_Tech = sum(PCA_Any == 1 & Sex_Any == 0 & Tech_Any == 1 & Zscore == 0 & pQTL == 0),
    PCA_Zscore = sum(PCA_Any == 1 & Sex_Any == 0 & Tech_Any == 0 & Zscore == 1 & pQTL == 0),
    PCA_pQTL = sum(PCA_Any == 1 & Sex_Any == 0 & Tech_Any == 0 & Zscore == 0 & pQTL == 1),
    Sex_Tech = sum(PCA_Any == 0 & Sex_Any == 1 & Tech_Any == 1 & Zscore == 0 & pQTL == 0),
    Sex_Zscore = sum(PCA_Any == 0 & Sex_Any == 1 & Tech_Any == 0 & Zscore == 1 & pQTL == 0),
    Sex_pQTL = sum(PCA_Any == 0 & Sex_Any == 1 & Tech_Any == 0 & Zscore == 0 & pQTL == 1),
    Tech_Zscore = sum(PCA_Any == 0 & Sex_Any == 0 & Tech_Any == 1 & Zscore == 1 & pQTL == 0),
    Tech_pQTL = sum(PCA_Any == 0 & Sex_Any == 0 & Tech_Any == 1 & Zscore == 0 & pQTL == 1),
    Zscore_pQTL = sum(PCA_Any == 0 & Sex_Any == 0 & Tech_Any == 0 & Zscore == 1 & pQTL == 1),
    Three_plus = sum(N_Methods >= 3)
  )]

  overlap_long <- data.table(
    Category = names(overlap_summary),
    Count = as.numeric(overlap_summary[1,])
  )
  overlap_long <- overlap_long[Count > 0]

  p3 <- ggplot(overlap_long, aes(x = reorder(Category, Count), y = Count)) +
    geom_bar(stat = "identity", fill = "#5C9EAD") +
    geom_text(aes(label = Count), hjust = -0.2, size = 3.5) +
    coord_flip() +
    labs(title = "Outlier Detection Overlap Patterns",
         x = "Overlap Category", y = "Number of Samples") +
    theme_bw()

  # Combine some plots
  combined <- grid.arrange(p1, p2, p3, ncol = 2, nrow = 2)

  return(list(
    method_counts = p1,
    n_methods_dist = p2,
    overlap_patterns = p3,
    combined = combined
  ))
}

# Function to create Z-score plots
create_zscore_plots <- function(zscore_matrix, outlier_result, npx_matrix, metadata = NULL) {
  log_info("Creating Z-score visualization plots")

  # Z-score distribution
  zscore_values <- as.vector(zscore_matrix)
  zscore_values <- zscore_values[!is.na(zscore_values)]

  p1 <- ggplot(data.frame(z = zscore_values), aes(x = z)) +
    geom_histogram(bins = 100, fill = "#3A5F8A", alpha = 0.7) +
    geom_vline(xintercept = c(-4, 4), color = "red", linetype = "dashed", linewidth = 1) +
    annotate("text", x = 4, y = Inf, label = "+4 SD", vjust = 1.5, hjust = -0.1, color = "red", size = 3.5) +
    annotate("text", x = -4, y = Inf, label = "-4 SD", vjust = 1.5, hjust = 1.1, color = "red", size = 3.5) +
    labs(title = "Z-score Distribution",
         subtitle = paste("Outlier threshold: |Z-score| > ±4 SD (mean ± 4 standard deviations)"),
         x = "Z-score", y = "Count") +
    theme_bw()

  # Sample outlier counts
  p2 <- ggplot(outlier_result$outlier_summary, aes(x = pct_outlier_proteins)) +
    geom_histogram(bins = 50, fill = "#5C9EAD", alpha = 0.7) +
    geom_vline(xintercept = 10, color = "red", linetype = "dashed") +
    labs(title = "Distribution of Outlier Proteins per Sample",
         x = "Percentage of Proteins with |Z| > 4",
         y = "Number of Samples") +
    theme_bw()

  # Scatter plot: max |Z-score| vs number of outlier proteins with top 20 highlighted
  # Add FINNGENID to outlier_summary for labeling
  if (!is.null(metadata) && "SAMPLE_ID" %in% names(metadata) && "FINNGENID" %in% names(metadata)) {
    outlier_with_id <- merge(
      outlier_result$outlier_summary,
      metadata[, .(SampleID = SAMPLE_ID, FINNGENID)],
      by = "SampleID",
      all.x = TRUE
    )
  } else {
    outlier_with_id <- copy(outlier_result$outlier_summary)
    outlier_with_id[, FINNGENID := NA_character_]
  }

  # Identify top 20 samples
  outlier_with_id[, is_top20 := FALSE]
  top20_ids <- outlier_with_id[order(-n_outlier_proteins)][1:min(20, .N)]$SampleID
  outlier_with_id[SampleID %in% top20_ids, is_top20 := TRUE]

  library(ggrepel)
  # Annotated version with FINNGENIDs
  p3 <- ggplot(outlier_with_id, aes(x = n_outlier_proteins, y = max_abs_zscore)) +
    geom_point(data = outlier_with_id[is_top20 == FALSE], aes(color = "Other"), alpha = 0.5, size = 2) +
    geom_point(data = outlier_with_id[is_top20 == TRUE], aes(color = "Top 20"), size = 3, alpha = 0.8) +
    geom_text_repel(data = outlier_with_id[is_top20 == TRUE & !is.na(FINNGENID)],
                    aes(label = FINNGENID),
                    size = 2.5,
                    max.overlaps = 20,
                    box.padding = 0.5,
                    point.padding = 0.3) +
    scale_color_manual(values = c("Other" = "#95A5A6", "Top 20" = "#E74C3C"),
                       name = "Samples") +
    labs(title = "Sample Outlier Characteristics (Annotated)",
         subtitle = "Top 20 samples by outlier protein count highlighted with FINNGENIDs",
         x = "Number of Outlier Proteins (|Z| > 4)",
         y = "Maximum Absolute Z-score") +
    theme_bw() +
    theme(legend.position = "bottom")

  # Unannotated version with color gradient based on max absolute z-score
  library(viridis)
  # Calculate thresholds
  n_proteins <- ncol(zscore_matrix)
  outlier_threshold_pct <- 0.1
  protein_count_threshold <- n_proteins * outlier_threshold_pct
  zscore_threshold <- 4

  # Count how many samples would be flagged
  n_flagged <- sum(outlier_with_id$n_outlier_proteins > protein_count_threshold)

  p3_unannotated <- ggplot(outlier_with_id, aes(x = n_outlier_proteins, y = max_abs_zscore, color = max_abs_zscore)) +
    geom_point(alpha = 0.7, size = 2.5) +
    # Add threshold lines
    geom_hline(yintercept = zscore_threshold, linetype = "dashed", color = "red", linewidth = 1) +
    geom_vline(xintercept = protein_count_threshold, linetype = "dashed", color = "blue", linewidth = 1) +
    # Add threshold annotations
    annotate("text", x = Inf, y = zscore_threshold, label = paste0("Z-score threshold = ", zscore_threshold),
             hjust = 1.05, vjust = -0.5, color = "red", size = 3.5, fontface = "bold") +
    annotate("text", x = protein_count_threshold, y = Inf,
             label = paste0("Sample threshold\n(10% of ", n_proteins, " proteins = ", round(protein_count_threshold), ")"),
             hjust = -0.05, vjust = 1.2, color = "blue", size = 3.5, fontface = "bold") +
    # Add shaded region for flagged samples
    annotate("rect", xmin = protein_count_threshold, xmax = Inf, ymin = zscore_threshold, ymax = Inf,
             alpha = 0.1, fill = "red") +
    scale_color_viridis_c(option = "plasma", direction = 1) +
    labs(title = "Sample Outlier Characteristics (Color Gradient)",
         subtitle = paste0("Samples flagged as outliers: ", n_flagged, " (requires >", round(protein_count_threshold), " proteins with |Z|>4)"),
         x = "Number of Outlier Proteins (|Z| > 4)",
         y = "Maximum Absolute Z-score",
         color = "Max |Z|") +
    theme_bw() +
    theme(legend.position = "right")

  # Maximum absolute Z-score per sample
  p4 <- ggplot(outlier_result$outlier_summary, aes(x = max_abs_zscore, fill = is_outlier)) +
    geom_histogram(bins = 50, alpha = 0.7) +
    scale_fill_manual(values = c("FALSE" = "#3A5F8A", "TRUE" = "#FF6B6B")) +
    labs(title = "Maximum Absolute Z-score per Sample",
         x = "Max |Z-score|", y = "Count",
         fill = "Outlier") +
    theme_bw()

  # Per-protein outlier frequency plot (TOP 10 ONLY)
  log_info("Creating protein outlier frequency plot (top 10)")
  protein_stats <- data.table(
    Protein = colnames(zscore_matrix),
    n_outlier_samples = outlier_result$protein_outlier_counts,
    mean_zscore = apply(zscore_matrix, 2, mean, na.rm = TRUE)
  )

  # Identify top 10 proteins by outlier count
  protein_stats <- protein_stats[order(-n_outlier_samples)]
  protein_stats[, is_top10 := FALSE]
  protein_stats[1:min(10, .N), is_top10 := TRUE]

  library(paletteer)
  p5 <- ggplot(protein_stats, aes(x = mean_zscore, y = n_outlier_samples)) +
    geom_point(data = protein_stats[is_top10 == FALSE],
               aes(color = n_outlier_samples), size = 2.5, alpha = 0.7) +
    geom_point(data = protein_stats[is_top10 == TRUE],
               aes(color = n_outlier_samples), size = 3.5, alpha = 0.9, shape = 21,
               stroke = 1.5, fill = "white") +
    geom_text_repel(data = protein_stats[is_top10 == TRUE],
                    aes(label = Protein),
                    size = 3,
                    max.overlaps = 15,
                    box.padding = 0.5,
                    point.padding = 0.3,
                    fontface = "bold") +
    scale_color_paletteer_c("ggthemes::Orange", direction = 1) +
    labs(title = "Protein Outlier Frequency",
         subtitle = "Top 10 most frequently flagged proteins annotated with black borders",
         x = "Mean Z-score",
         y = "Number of Outlier Samples",
         color = "Outlier Count") +
    theme_bw() +
    theme(legend.position = "right")

  # NPX vs Z-score comparison plot with MARGINAL HISTOGRAMS
  log_info("Creating NPX vs Z-score comparison plot with marginal histograms")

  # Flatten matrices to vectors (sample randomly to avoid overplotting)
  npx_values <- as.vector(npx_matrix)
  zscore_values_all <- as.vector(zscore_matrix)

  # Remove NA pairs
  valid_idx <- !is.na(npx_values) & !is.na(zscore_values_all)
  npx_clean <- npx_values[valid_idx]
  zscore_clean <- zscore_values_all[valid_idx]

  # Sample for plotting (to avoid too many points)
  if (length(npx_clean) > 50000) {
    sample_idx <- sample(length(npx_clean), 50000)
    npx_plot <- npx_clean[sample_idx]
    zscore_plot <- zscore_clean[sample_idx]
  } else {
    npx_plot <- npx_clean
    zscore_plot <- zscore_clean
  }

  # Calculate regression statistics
  lm_fit <- lm(zscore_plot ~ npx_plot)
  lm_summary <- summary(lm_fit)
  r_squared <- lm_summary$r.squared
  p_value <- lm_summary$coefficients[2, 4]

  # Calculate summary statistics
  npx_stats <- data.frame(
    Mean = mean(npx_clean, na.rm = TRUE),
    Median = median(npx_clean, na.rm = TRUE),
    SD = sd(npx_clean, na.rm = TRUE),
    Min = min(npx_clean, na.rm = TRUE),
    Max = max(npx_clean, na.rm = TRUE)
  )

  zscore_stats <- data.frame(
    Mean = mean(zscore_clean, na.rm = TRUE),
    Median = median(zscore_clean, na.rm = TRUE),
    SD = sd(zscore_clean, na.rm = TRUE),
    Min = min(zscore_clean, na.rm = TRUE),
    Max = max(zscore_clean, na.rm = TRUE)
  )

  # Format statistics text
  npx_stats_text <- sprintf(
    "NPX: Mean=%.2f, Median=%.2f, SD=%.2f",
    npx_stats$Mean, npx_stats$Median, npx_stats$SD
  )

  zscore_stats_text <- sprintf(
    "Z-score: Mean=%.2f, Median=%.2f, SD=%.2f",
    zscore_stats$Mean, zscore_stats$Median, zscore_stats$SD
  )

  # Format p-value
  if (p_value < 2.2e-16) {
    p_text <- "p < 2.2e-16"
  } else if (p_value < 0.001) {
    p_text <- sprintf("p = %.2e", p_value)
  } else {
    p_text <- sprintf("p = %.4f", p_value)
  }

  regression_text <- sprintf("R² = %.4f, %s", r_squared, p_text)

  # Create comparison data frame
  comparison_df <- data.frame(
    NPX = npx_plot,
    Zscore = zscore_plot
  )

  # Main scatter plot with density contours
  library(ggExtra)
  p6_base <- ggplot(comparison_df, aes(x = NPX, y = Zscore)) +
    stat_density_2d(aes(fill = after_stat(level)), geom = "polygon", alpha = 0.3, bins = 20, color = "gray70", linewidth = 0.3) +
    geom_point(alpha = 0.03, size = 0.3, color = "gray30") +
    geom_smooth(method = "lm", color = "#1f4788", linetype = "dashed", linewidth = 0.8, se = TRUE, alpha = 0.15) +
    scale_fill_viridis_c(option = "plasma") +
    annotate("text", x = Inf, y = Inf, label = regression_text,
             hjust = 1.05, vjust = 1.5, size = 4, fontface = "bold",
             color = "#1f4788") +
    annotate("label", x = -Inf, y = Inf, label = npx_stats_text,
             hjust = -0.05, vjust = 1.1, size = 3, fontface = "plain",
             color = "black", fill = "white", alpha = 0.9) +
    annotate("label", x = Inf, y = -Inf, label = zscore_stats_text,
             hjust = 1.05, vjust = -0.1, size = 3, fontface = "plain",
             color = "black", fill = "white", alpha = 0.9) +
    labs(title = "NPX vs Z-score Comparison with Marginal Distributions",
         subtitle = "Linear regression (dashed dark blue) with 2D density contours and marginal histograms",
         x = "Original NPX Value",
         y = "Z-score",
         fill = "Density") +
    theme_bw() +
    theme(legend.position = "right",
          panel.grid.minor = element_blank())

  # Add marginal histograms
  p6 <- ggMarginal(p6_base, type = "histogram", fill = "#3A5F8A", alpha = 0.7, bins = 50)

  log_info("NPX vs Z-score plot created: R² = {round(r_squared, 4)}, p {ifelse(p_value < 0.001, '< 0.001', paste('=', round(p_value, 4)))}")

  return(list(
    zscore_dist = p1,
    outlier_pct = p2,
    sample_scatter = p3,
    sample_scatter_unannotated = p3_unannotated,
    max_zscore = p4,
    protein_frequency = p5,
    npx_zscore_comparison = p6
  ))
}

# Function to analyze disease group distribution for top outlier proteins
analyze_protein_outlier_disease_groups <- function(zscore_matrix, outlier_result, metadata, top_n = 40, threshold = 4) {
  log_info("Analyzing disease group distribution for top {top_n} outlier proteins")

  # Get top N proteins by outlier count
  protein_stats <- data.table(
    Protein = colnames(zscore_matrix),
    n_outlier_samples = outlier_result$protein_outlier_counts
  )
  protein_stats <- protein_stats[order(-n_outlier_samples)]
  top_proteins <- protein_stats[1:min(top_n, .N)]$Protein

  log_info("Top {length(top_proteins)} proteins selected: {paste(head(top_proteins, 5), collapse=', ')}...")

  # Identify disease group columns in metadata (binary columns)
  disease_cols <- c("Kidney", "Kids", "F64", "MFGE8", "Parkinsons", "Metabolic",
                    "AMD", "Rheuma", "Pulmo", "Chromosomal_Abnormalities",
                    "Blood_donors", "Bridging_samples")
  available_disease_cols <- disease_cols[disease_cols %in% names(metadata)]

  if (length(available_disease_cols) == 0) {
    log_warn("No disease group columns found in metadata")
    return(NULL)
  }

  log_info("Found {length(available_disease_cols)} disease group columns")

  # For each protein, identify outlier samples and their disease groups
  disease_group_summary <- rbindlist(lapply(top_proteins, function(prot) {
    # Get outlier samples for this protein
    prot_zscores <- zscore_matrix[, prot]
    outlier_samples <- names(prot_zscores)[!is.na(prot_zscores) & abs(prot_zscores) > threshold]

    # Match with metadata
    if (!is.null(metadata) && "SAMPLE_ID" %in% names(metadata)) {
      sample_metadata <- metadata[SAMPLE_ID %in% outlier_samples]

      # Assign each sample to its disease group(s)
      if (nrow(sample_metadata) > 0) {
        # For each sample, find which disease groups it belongs to
        sample_disease_groups <- lapply(1:nrow(sample_metadata), function(i) {
          row <- sample_metadata[i, ]
          groups <- character()
          for (col in available_disease_cols) {
            if (!is.na(row[[col]]) && row[[col]] == 1) {
              groups <- c(groups, col)
            }
          }
          if (length(groups) == 0) groups <- "No_Disease_Group"
          return(data.table(SAMPLE_ID = row$SAMPLE_ID, DISEASE_GROUP = groups))
        })

        sample_disease_dt <- rbindlist(sample_disease_groups)

        # Count by disease group
        disease_counts <- sample_disease_dt[, .N, by = DISEASE_GROUP]
        disease_counts[, Protein := prot]
        disease_counts[, Total_Outliers := length(outlier_samples)]
        disease_counts[, Percentage := (N / Total_Outliers) * 100]
        return(disease_counts)
      }
    }
    return(NULL)
  }), fill = TRUE)

  # Filter out proteins with no disease group info
  if (is.null(disease_group_summary) || nrow(disease_group_summary) == 0) {
    log_warn("No disease group information available for outlier proteins")
    return(NULL)
  }

  # Reorder proteins by total outliers
  disease_group_summary[, Protein := factor(Protein, levels = unique(protein_stats[Protein %in% top_proteins]$Protein))]

  # Order disease groups by total number of outlier samples (descending)
  disease_group_totals <- disease_group_summary[, .(Total = sum(N)), by = DISEASE_GROUP][order(-Total)]
  disease_group_summary[, DISEASE_GROUP := factor(DISEASE_GROUP, levels = disease_group_totals$DISEASE_GROUP)]

  # Create heatmap-style visualization
  p_heatmap <- ggplot(disease_group_summary, aes(x = DISEASE_GROUP, y = Protein, fill = N)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = N), size = 2.5, color = "white", fontface = "bold") +
    scale_fill_viridis_c(option = "rocket", direction = -1) +
    labs(title = "Disease Group Distribution for Top 40 Outlier Proteins",
         subtitle = "Number of samples per disease group (ordered by total outlier count, |Z| > 4)",
         x = "Disease Group",
         y = "Protein",
         fill = "Count") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
          axis.text.y = element_text(size = 7),
          legend.position = "right",
          panel.grid = element_blank())

  # Create stacked bar chart showing percentage composition with ordered bars
  # Order disease groups within each protein by frequency (descending)
  disease_group_summary <- disease_group_summary[order(Protein, -N)]
  disease_group_summary[, DISEASE_GROUP := factor(DISEASE_GROUP, levels = unique(DISEASE_GROUP))]

  # Custom color palette as requested
  custom_colors <- c("#DF2E28FF", "#FE801AFF", "#E9BF35FF", "#81BB42FF", "#32C7A9FF", "#4A9BDCFF")

  # Expand palette if more disease groups than colors
  n_disease_groups <- length(unique(disease_group_summary$DISEASE_GROUP))
  if (n_disease_groups > length(custom_colors)) {
    custom_colors <- colorRampPalette(custom_colors)(n_disease_groups)
  }

  p_stacked <- ggplot(disease_group_summary, aes(x = Protein, y = N, fill = DISEASE_GROUP)) +
    geom_bar(stat = "identity", position = "fill") +
    scale_y_continuous(labels = scales::percent) +
    scale_fill_manual(values = custom_colors) +
    labs(title = "Disease Group Composition for Top 40 Outlier Proteins",
         subtitle = "Percentage composition of disease groups (ordered by frequency within each protein)",
         x = "Protein",
         y = "Percentage of Outlier Samples",
         fill = "Disease Group") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
          legend.position = "bottom",
          legend.text = element_text(size = 8))

  log_info("Disease group analysis complete for {length(unique(disease_group_summary$Protein))} proteins")

  return(list(
    summary_table = disease_group_summary,
    heatmap = p_heatmap,
    stacked_bar = p_stacked
  ))
}

# Main execution
main <- function() {
  # Load data from previous steps
  log_info("Loading data from previous steps")
  # CRITICAL: Use base matrix from step 00 (analysis_ready) for parallel flagging
  # All outlier detection steps (01, 02, 03) now use the same base matrix for parallel flagging
  # This matches the original pipeline design where technical and Z-score detection operate
  # on the same base matrix independently, enabling parallel flagging with union logic in final QC
  prev_step00_num <- "00"
  npx_matrix_path <- get_output_path(prev_step00_num, "npx_matrix_analysis_ready", batch_id, "qc", config = config)
  metadata_path <- get_output_path(prev_step00_num, "metadata", batch_id, "qc", config = config)

  if (!file.exists(npx_matrix_path)) {
    stop("NPX matrix file not found: ", npx_matrix_path, ". Run Step 00 first.")
  }
  if (!file.exists(metadata_path)) {
    stop("Metadata file not found: ", metadata_path, ". Run Step 00 first.")
  }

  log_info("Using base NPX matrix (step 00) for parallel outlier flagging: {npx_matrix_path}")
  npx_matrix <- readRDS(npx_matrix_path)
  metadata <- readRDS(metadata_path)
  log_info("Loaded NPX matrix: {nrow(npx_matrix)} samples x {ncol(npx_matrix)} proteins")

  # Calculate Z-scores
  zscore_matrix <- calculate_protein_zscores(npx_matrix)

  # Detect outliers using configured threshold
  outlier_result <- detect_zscore_outliers(
    zscore_matrix,
    threshold = config$parameters$outliers$zscore_threshold
  )

  # Perform iterative outlier detection
  iterative_result <- iterative_outlier_detection(
    npx_matrix,
    max_iterations = 5,
    threshold = config$parameters$outliers$zscore_threshold
  )

  # Analyze outlier patterns
  pattern_analysis <- analyze_outlier_patterns(zscore_matrix, outlier_result, metadata)

  # Create integrated outlier tracking across all steps (01-05)
  integrated_tracking <- create_integrated_outlier_tracking(outlier_result$sample_outliers, metadata, batch_id = batch_id, config = config)

  # Create integrated tracking visualizations
  integrated_plots <- create_integrated_plots(integrated_tracking)

  # Create plots
  zscore_plots <- create_zscore_plots(zscore_matrix, outlier_result, npx_matrix, metadata)

  # Analyze disease group distribution for top 40 outlier proteins
  disease_group_analysis <- analyze_protein_outlier_disease_groups(zscore_matrix, outlier_result, metadata, top_n = 40, threshold = 4)

  # Remove Z-score outliers from matrix
  npx_clean <- npx_matrix[!rownames(npx_matrix) %in% outlier_result$sample_outliers, ]

  # Also create version with iterative outliers removed
  npx_iterative_clean <- npx_matrix[!rownames(npx_matrix) %in% iterative_result$all_outliers, ]

  log_info("Matrix after Z-score outlier removal: {nrow(npx_clean)} x {ncol(npx_clean)}")
  log_info("Matrix after iterative outlier removal: {nrow(npx_iterative_clean)} x {ncol(npx_iterative_clean)}")

  # Save outputs with batch-aware paths
  log_info("Saving Z-score outlier detection results")

  # Always save these
  zscore_matrix_path <- get_output_path(step_num, "zscore_matrix", batch_id, "outliers", config = config)
  pattern_analysis_path <- get_output_path(step_num, "zscore_pattern_analysis", batch_id, "outliers", config = config)
  ensure_output_dir(zscore_matrix_path)
  ensure_output_dir(pattern_analysis_path)

  saveRDS(zscore_matrix, zscore_matrix_path)
  saveRDS(pattern_analysis, pattern_analysis_path)

  # Conditional saves based on whether outliers were detected
  if (length(outlier_result$sample_outliers) > 0) {
    log_info("Saving outlier-specific files (outliers detected: {length(outlier_result$sample_outliers)})")
    zscore_outliers_path <- get_output_path(step_num, "zscore_outliers", batch_id, "outliers", config = config)
    npx_zscore_cleaned_path <- get_output_path(step_num, "npx_matrix_zscore_cleaned", batch_id, "outliers", config = config)
    ensure_output_dir(zscore_outliers_path)
    ensure_output_dir(npx_zscore_cleaned_path)

    saveRDS(outlier_result, zscore_outliers_path)
    saveRDS(npx_clean, npx_zscore_cleaned_path)
  } else {
    log_info("No z-score outliers detected in single-pass detection - skipping outlier-specific RDS files")
  }

  if (length(iterative_result$all_outliers) > 0) {
    log_info("Saving iterative outlier files (iterative outliers detected: {length(iterative_result$all_outliers)})")
    zscore_iterative_path <- get_output_path(step_num, "zscore_outliers_iterative", batch_id, "outliers", config = config)
    npx_iterative_cleaned_path <- get_output_path(step_num, "npx_matrix_iterative_cleaned", batch_id, "outliers", config = config)
    ensure_output_dir(zscore_iterative_path)
    ensure_output_dir(npx_iterative_cleaned_path)

    saveRDS(iterative_result, zscore_iterative_path)
    saveRDS(npx_iterative_clean, npx_iterative_cleaned_path)
  } else {
    log_info("No z-score outliers detected in iterative detection - skipping iterative-specific RDS files")
  }

  # Save tables with FINNGENID
  log_info("Adding FINNGENID mapping to output tables...")

  # Always save the summary (even if no outliers, it contains all sample stats)
  outlier_summary_with_fgid <- add_finngenid_column(outlier_result$outlier_summary, batch_id = batch_id, config = config)
  outlier_summary_path <- get_output_path(step_num, "zscore_outlier_summary", batch_id, "outliers", "tsv", config = config)
  ensure_output_dir(outlier_summary_path)
  fwrite(outlier_summary_with_fgid, outlier_summary_path, sep = "\t")

  # Update outlier_result for RDS
  outlier_result$outlier_summary <- outlier_summary_with_fgid

  # Conditional outlier list files
  if (length(outlier_result$sample_outliers) > 0) {
    outliers_list_with_fgid <- add_finngenid_column(
      data.table(SampleID = outlier_result$sample_outliers),
      batch_id = batch_id, config = config
    )
    outliers_list_path <- get_output_path(step_num, "zscore_outliers_list", batch_id, "outliers", "tsv", config = config)
    ensure_output_dir(outliers_list_path)
    fwrite(outliers_list_with_fgid, outliers_list_path, sep = "\t")
  } else {
    log_info("No z-score outliers detected - skipping zscore_outliers_list.tsv")
  }

  if (length(iterative_result$all_outliers) > 0) {
    iterative_list_with_fgid <- add_finngenid_column(
      data.table(SampleID = iterative_result$all_outliers),
      batch_id = batch_id, config = config
    )
    iterative_list_path <- get_output_path(step_num, "zscore_outliers_iterative_list", batch_id, "outliers", "tsv", config = config)
    ensure_output_dir(iterative_list_path)
    fwrite(iterative_list_with_fgid, iterative_list_path, sep = "\t")
  } else {
    log_info("No iterative z-score outliers detected - skipping zscore_outliers_iterative_list.tsv")
  }

  # Save integrated outlier tracking across all steps
  # Ensure FINNGENID is present (it may already be there from PCA file merge)

  if (!"FINNGENID" %in% names(integrated_tracking) || sum(!is.na(integrated_tracking$FINNGENID)) == 0) {
    integrated_tracking <- add_finngenid_column(integrated_tracking, batch_id = batch_id, config = config)
  }
  integrated_tracking_path <- get_output_path(step_num, "outliers_by_step_integrated", batch_id, "outliers", "tsv", config = config)
  ensure_output_dir(integrated_tracking_path)
  fwrite(integrated_tracking, integrated_tracking_path, sep = "\t")
  log_info("Saved integrated outlier tracking: {nrow(integrated_tracking)} samples tracked across steps 01-05b")

  # Only save cohort summary if there are z-score outliers
  if (length(outlier_result$sample_outliers) > 0 && !is.null(pattern_analysis$cohort_summary)) {
    # Check if cohort_summary has SampleID
    if ("SampleID" %in% names(pattern_analysis$cohort_summary)) {
      cohort_summary_with_fgid <- add_finngenid_column(pattern_analysis$cohort_summary, batch_id = batch_id, config = config)
    } else {
      cohort_summary_with_fgid <- pattern_analysis$cohort_summary
    }
    cohort_summary_path <- get_output_path(step_num, "zscore_outlier_cohorts", batch_id, "outliers", "tsv", config = config)
    ensure_output_dir(cohort_summary_path)
    fwrite(cohort_summary_with_fgid, cohort_summary_path, sep = "\t")
    log_info("Saved Z-score outlier cohort summary")
  } else {
    log_info("No Z-score outliers detected - skipping zscore_outlier_cohorts.tsv")
  }

  # Save plots with batch-aware paths
  zscore_dist_path <- get_output_path(step_num, "zscore_distribution", batch_id, "outliers", "pdf", config = config)
  zscore_outlier_pct_path <- get_output_path(step_num, "zscore_outlier_percentage", batch_id, "outliers", "pdf", config = config)
  sample_scatter_path <- get_output_path(step_num, "sample_outlier_scatter", batch_id, "outliers", "pdf", config = config)
  sample_scatter_gradient_path <- get_output_path(step_num, "sample_outlier_scatter_gradient", batch_id, "outliers", "pdf", config = config)
  zscore_max_path <- get_output_path(step_num, "zscore_max_per_sample", batch_id, "outliers", "pdf", config = config)
  protein_frequency_path <- get_output_path(step_num, "protein_outlier_frequency", batch_id, "outliers", "pdf", config = config)
  npx_zscore_comparison_path <- get_output_path(step_num, "npx_zscore_comparison", batch_id, "outliers", "pdf", config = config)

  ensure_output_dir(zscore_dist_path)
  ensure_output_dir(zscore_outlier_pct_path)
  ensure_output_dir(sample_scatter_path)
  ensure_output_dir(sample_scatter_gradient_path)
  ensure_output_dir(zscore_max_path)
  ensure_output_dir(protein_frequency_path)
  ensure_output_dir(npx_zscore_comparison_path)

  ggsave(zscore_dist_path, zscore_plots$zscore_dist, width = 10, height = 6)
  ggsave(zscore_outlier_pct_path, zscore_plots$outlier_pct, width = 8, height = 6)
  ggsave(sample_scatter_path, zscore_plots$sample_scatter, width = 12, height = 8)
  ggsave(sample_scatter_gradient_path, zscore_plots$sample_scatter_unannotated, width = 10, height = 8)
  ggsave(zscore_max_path, zscore_plots$max_zscore, width = 8, height = 6)
  ggsave(protein_frequency_path, zscore_plots$protein_frequency, width = 11, height = 8)
  ggsave(npx_zscore_comparison_path, zscore_plots$npx_zscore_comparison, width = 10, height = 8)

  # Save integrated tracking plots with batch-aware paths
  method_counts_path <- get_output_path(step_num, "outlier_method_counts", batch_id, "outliers", "pdf", config = config)
  n_methods_dist_path <- get_output_path(step_num, "outlier_multi_method_dist", batch_id, "outliers", "pdf", config = config)
  overlap_patterns_path <- get_output_path(step_num, "outlier_overlap_patterns", batch_id, "outliers", "pdf", config = config)
  tracking_combined_path <- get_output_path(step_num, "outlier_tracking_combined", batch_id, "outliers", "pdf", config = config)

  ensure_output_dir(method_counts_path)
  ensure_output_dir(n_methods_dist_path)
  ensure_output_dir(overlap_patterns_path)
  ensure_output_dir(tracking_combined_path)

  ggsave(method_counts_path, integrated_plots$method_counts, width = 8, height = 6)
  ggsave(n_methods_dist_path, integrated_plots$n_methods_dist, width = 8, height = 6)
  ggsave(overlap_patterns_path, integrated_plots$overlap_patterns, width = 10, height = 6)
  ggsave(tracking_combined_path, integrated_plots$combined, width = 14, height = 10)
  log_info("Saved integrated tracking visualizations")

  # Save disease group analysis plots and table
  if (!is.null(disease_group_analysis)) {
    disease_heatmap_path <- get_output_path(step_num, "protein_disease_group_heatmap", batch_id, "outliers", "pdf", config = config)
    disease_composition_path <- get_output_path(step_num, "protein_disease_group_composition", batch_id, "outliers", "pdf", config = config)
    disease_summary_path <- get_output_path(step_num, "protein_disease_group_summary", batch_id, "outliers", "tsv", config = config)

    ensure_output_dir(disease_heatmap_path)
    ensure_output_dir(disease_composition_path)
    ensure_output_dir(disease_summary_path)

    ggsave(disease_heatmap_path, disease_group_analysis$heatmap, width = 12, height = 16)
    ggsave(disease_composition_path, disease_group_analysis$stacked_bar, width = 16, height = 8)

    # Disease group summary doesn't have SampleID, so no FINNGENID needed
    fwrite(disease_group_analysis$summary_table, disease_summary_path, sep = "\t")
    log_info("Saved disease group analysis for top 40 outlier proteins")
  }

  # Print summary
  cat("\n=== Z-SCORE OUTLIER DETECTION SUMMARY ===\n")
  cat("Z-score threshold:", config$parameters$outliers$zscore_threshold, "\n")
  cat("Single-pass outliers:", length(outlier_result$sample_outliers), "\n")
  cat("Iterative outliers:", length(iterative_result$all_outliers), "\n")
  cat("\nTop outlier proteins:\n")
  cat(paste("  -", pattern_analysis$top_outlier_proteins[1:min(5, length(pattern_analysis$top_outlier_proteins))]), sep = "\n")
  cat("\nMatrix after single-pass removal:", nrow(npx_clean), "x", ncol(npx_clean), "\n")
  cat("Matrix after iterative removal:", nrow(npx_iterative_clean), "x", ncol(npx_iterative_clean), "\n")

  cat("\n=== INTEGRATED OUTLIER TRACKING (STEPS 04-05b) ===\n")
  cat("Total samples flagged by any method:", nrow(integrated_tracking), "\n")
  cat("Samples flagged by multiple methods:", sum(integrated_tracking$N_Methods > 1), "\n")
  cat("  - 2 methods:", sum(integrated_tracking$N_Methods == 2), "\n")
  cat("  - 3 methods:", sum(integrated_tracking$N_Methods == 3), "\n")
  cat("  - 4 methods:", sum(integrated_tracking$N_Methods == 4), "\n")
  cat("  - 5 methods:", sum(integrated_tracking$N_Methods == 5), "\n")
  cat("\nOutliers by step:\n")
  cat("  - PCA:", sum(integrated_tracking$PCA_Any), "\n")
  cat("  - Sex (any):", sum(integrated_tracking$Sex_Any), "\n")
  cat("    - Sex mismatch (strict):", sum(integrated_tracking$SexMismatch), "\n")
  cat("    - Sex outlier (relaxed):", sum(integrated_tracking$SexOutlier), "\n")
  cat("  - Technical:", sum(integrated_tracking$Tech_Any), "\n")
  cat("  - Z-score:", sum(integrated_tracking$Zscore), "\n")
  cat("  - pQTL:", sum(integrated_tracking$pQTL %||% 0), "\n")

  cat("\nResults saved to: ../output/outliers/\n")
  cat("  - Integrated tracking: 07_outliers_by_step_integrated.tsv\n")

  log_info("Z-score outlier detection completed")

  return(list(
    zscore_matrix = zscore_matrix,
    outliers = outlier_result,
    iterative_outliers = iterative_result,
    clean_matrix = npx_clean
  ))
}

# Run if executed directly
if (!interactive()) {
  result <- main()
}
