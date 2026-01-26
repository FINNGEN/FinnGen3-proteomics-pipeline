#!/usr/bin/env Rscript
# ==============================================================================
# 06_normalize_data.R - Proteomics Data Normalisation
# ==============================================================================
#
# Purpose:
#   Normalises proteomics data to remove technical variability and ensure sample
#   comparability. Mode-dependent behaviour:
#   - Single-batch mode: Uses median normalisation only (standard intra-batch step).
#     Bridge and ComBat normalisations are NOT applicable in single-batch mode.
#   - Multi-batch mode: Uses bridge normalisation (primary) with fallback to ComBat
#     or median if bridge samples are insufficient. Bridge and ComBat are ONLY
#     applicable in multi-batch (cross-batch) mode. Evaluates normalisation
#     effectiveness using SD, MAD, and IQR reduction metrics.
#
# Author: Reza Jabal, PhD (rjabal@broadinstitute.org)
# Date: December 2025
# ==============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(tidyverse)
  library(OlinkAnalyze)
  library(sva)  # For ComBat
  library(yaml)
  library(logger)
  library(ggplot2)
})

# Source path utilities for batch-aware paths
# Get script directory safely (handles both direct execution and sourcing)
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
}, error = function(e) getwd())
source(file.path(script_dir, "path_utils.R"))

# Get batch context
batch_id <- Sys.getenv("PIPELINE_BATCH_ID", "batch_01")
step_num <- get_step_number()

# Load configuration
config_file <- Sys.getenv("PIPELINE_CONFIG", "")
if (config_file == "" || !file.exists(config_file)) {
  stop("PIPELINE_CONFIG environment variable not set or config file not found. Please provide path to config file.")
}
config <- read_yaml(config_file)

# Set up logging with batch-aware path
log_path <- get_log_path(step_num, batch_id, config = config)
log_appender(appender_file(log_path))
log_info("Starting data normalisation for batch: {batch_id}")

# Function to check protein consistency between batches
check_protein_consistency <- function(batch2_proteins, batch1_proteins) {
  log_info("Checking protein consistency between batches")

  common_proteins <- intersect(batch2_proteins, batch1_proteins)
  batch2_only <- setdiff(batch2_proteins, batch1_proteins)
  batch1_only <- setdiff(batch1_proteins, batch2_proteins)

  log_info("Protein consistency:")
  log_info("  Common proteins: {length(common_proteins)}")
  log_info("  Batch 2 only: {length(batch2_only)}")
  log_info("  Batch 1 only: {length(batch1_only)}")

  if (length(batch2_only) > 0) {
    log_warn("Batch 2 has {length(batch2_only)} proteins not in batch 1")
    log_warn("These will be excluded from cross-batch normalization")
  }

  if (length(batch1_only) > 0) {
    log_warn("Batch 1 has {length(batch1_only)} proteins not in batch 2")
    log_warn("These will be excluded from cross-batch normalization")
  }

  return(list(
    common_proteins = common_proteins,
    batch2_only = batch2_only,
    batch1_only = batch1_only
  ))
}

# Function for bridge sample normalisation
# NOTE: Bridge normalization is ONLY applicable in multi-batch (cross-batch) mode
normalize_bridge <- function(npx_matrix, bridge_samples, reference_batch = NULL, bridge_mapping = NULL, multi_batch_mode = FALSE) {
  if (!multi_batch_mode) {
    log_error("Bridge normalization is only applicable in multi-batch (cross-batch) mode")
    stop("Bridge normalization cannot be used in single-batch mode")
  }
  log_info("Performing bridge sample normalisation")

  # Identify bridge samples in matrix
  bridge_in_matrix <- intersect(rownames(npx_matrix), bridge_samples)

  if(length(bridge_in_matrix) < 10) {
    log_warn("Too few bridge samples ({length(bridge_in_matrix)}) for normalization")
    return(NULL)
  }

  log_info("Using {length(bridge_in_matrix)} bridge samples for normalisation")

  # Calculate scaling factors based on bridge samples
  bridge_data <- npx_matrix[bridge_in_matrix, ]

  # If multi-batch mode and reference batch provided, use cross-batch normalisation
  if(multi_batch_mode && !is.null(reference_batch) && !is.null(bridge_mapping)) {
    log_info("Multi-batch mode: Performing cross-batch bridge normalisation")

    # Check protein consistency
    protein_check <- check_protein_consistency(colnames(npx_matrix), colnames(reference_batch))
    common_proteins <- protein_check$common_proteins

    if (length(common_proteins) < 100) {
      log_error("Too few common proteins ({length(common_proteins)}) for cross-batch normalization")
      return(NULL)
    }

    # Use only common proteins
    bridge_data <- bridge_data[, common_proteins, drop = FALSE]
    npx_matrix_subset <- npx_matrix[, common_proteins, drop = FALSE]

    # Map bridge samples from batch 2 to batch 1
    bridge_mapping_subset <- bridge_mapping[SAMPLE_ID_batch2 %in% bridge_in_matrix & !is.na(SAMPLE_ID_batch1)]
    batch1_bridge_ids <- bridge_mapping_subset$SAMPLE_ID_batch1
    batch1_bridge_in_matrix <- intersect(batch1_bridge_ids, rownames(reference_batch))

    if (length(batch1_bridge_in_matrix) < 10) {
      log_warn("Too few bridge samples in batch 1 ({length(batch1_bridge_in_matrix)}) for cross-batch normalization")
      log_warn("Falling back to single-batch normalization")
    } else {
      log_info("Found {length(batch1_bridge_in_matrix)} bridge samples in batch 1")

      # Extract batch 1 bridge samples (common proteins only)
      batch1_bridge_data <- reference_batch[batch1_bridge_in_matrix, common_proteins, drop = FALSE]

      # Combine bridge samples from both batches
      # Use median across both batches for each protein
      batch2_medians <- apply(bridge_data, 2, median, na.rm = TRUE)
      batch1_medians <- apply(batch1_bridge_data, 2, median, na.rm = TRUE)

      # Use combined median (average of batch medians) as reference
      combined_medians <- (batch2_medians + batch1_medians) / 2

      # Calculate global median across all proteins
      global_median <- median(combined_medians, na.rm = TRUE)

      # Calculate scaling factors for batch 2
      scaling_factors <- global_median / batch2_medians
      scaling_factors[is.na(scaling_factors) | is.infinite(scaling_factors)] <- 1

      # Cap extreme scaling factors
      scaling_factors[scaling_factors < 0.1] <- 0.1
      scaling_factors[scaling_factors > 10] <- 10

      log_info("Cross-batch scaling factors - Range: [{round(min(scaling_factors), 3)}, {round(max(scaling_factors), 3)}]")

      # Apply scaling to batch 2 matrix (common proteins only)
      normalized_matrix_subset <- sweep(npx_matrix_subset, 2, scaling_factors, "*")

      # Reconstruct full matrix (non-common proteins unchanged)
      normalized_matrix <- npx_matrix
      normalized_matrix[, common_proteins] <- normalized_matrix_subset

      return(list(
        normalized_matrix = normalized_matrix,
        scaling_factors = scaling_factors,
        bridge_samples_used = bridge_in_matrix,
        bridge_samples_batch1 = batch1_bridge_in_matrix,
        common_proteins = common_proteins,
        method = "bridge_cross_batch"
      ))
    }
  }

  # Single-batch normalization (fallback or default)
  log_info("Performing single-batch bridge normalisation")

  # Calculate median for each protein in bridge samples
  bridge_medians <- apply(bridge_data, 2, median, na.rm = TRUE)

  # Calculate global median across all proteins in bridge samples
  global_median <- median(bridge_medians, na.rm = TRUE)

  # Calculate scaling factors
  scaling_factors <- global_median / bridge_medians
  scaling_factors[is.na(scaling_factors) | is.infinite(scaling_factors)] <- 1

  # Cap extreme scaling factors to prevent outliers (keep within 0.1 to 10 range)
  scaling_factors[scaling_factors < 0.1] <- 0.1
  scaling_factors[scaling_factors > 10] <- 10

  log_info("Scaling factors - Range: [{round(min(scaling_factors), 3)}, {round(max(scaling_factors), 3)}]")

  # Apply scaling to entire matrix
  normalized_matrix <- sweep(npx_matrix, 2, scaling_factors, "*")

  # Note: NPX data is already log2-transformed by Olink, so we do NOT apply additional log transformation

  return(list(
    normalized_matrix = normalized_matrix,
    scaling_factors = scaling_factors,
    bridge_samples_used = bridge_in_matrix,
    method = "bridge"
  ))
}

# Function for median normalisation
normalize_median <- function(npx_matrix, by_plate = FALSE, plate_info = NULL) {
  log_info("Performing median normalisation")

  if(by_plate && !is.null(plate_info)) {
    log_info("Normalising by plate")

    # Group samples by plate
    normalized_list <- list()

    for(plate in unique(plate_info$PlateID)) {
      plate_samples <- plate_info[PlateID == plate]$SampleID
      plate_samples <- intersect(plate_samples, rownames(npx_matrix))

      if(length(plate_samples) > 0) {
        plate_matrix <- npx_matrix[plate_samples, , drop = FALSE]

        # Median normalise within plate
        col_medians <- apply(plate_matrix, 2, median, na.rm = TRUE)
        global_median <- median(col_medians, na.rm = TRUE)
        scaling_factors <- global_median / col_medians
        scaling_factors[is.na(scaling_factors) | is.infinite(scaling_factors)] <- 1

        normalized_list[[plate]] <- sweep(plate_matrix, 2, scaling_factors, "*")
      }
    }

    # Combine normalised plates
    normalized_matrix <- do.call(rbind, normalized_list)

  } else {
    log_info("Performing global median normalization")

    # Calculate column medians
    col_medians <- apply(npx_matrix, 2, median, na.rm = TRUE)

    # Calculate global median
    global_median <- median(col_medians, na.rm = TRUE)

    # Calculate scaling factors
    scaling_factors <- global_median / col_medians
    scaling_factors[is.na(scaling_factors) | is.infinite(scaling_factors)] <- 1

    # Cap extreme scaling factors to prevent outliers (keep within 0.1 to 10 range)
    scaling_factors[scaling_factors < 0.1] <- 0.1
    scaling_factors[scaling_factors > 10] <- 10

    # Apply scaling
    normalized_matrix <- sweep(npx_matrix, 2, scaling_factors, "*")
  }

  # Note: NPX data is already log2-transformed by Olink, so we do NOT apply additional log transformation

  return(list(
    normalized_matrix = normalized_matrix,
    method = ifelse(by_plate, "median_by_plate", "median_global")
  ))
}

# Function for ComBat normalisation
# NOTE: ComBat normalization is ONLY applicable in multi-batch (cross-batch) mode
# ComBat is designed to correct for batch effects between multiple batches
# When used as a fallback, it may normalize within a single batch (removing plate/technical effects)
normalize_combat <- function(npx_matrix, batch_info, covariates = NULL) {
  if (is.null(batch_info)) {
    log_error("Batch info is required for ComBat normalization")
    stop("Missing batch info for ComBat")
  }
  log_info("Performing ComBat batch correction")

  # Prepare batch vector
  sample_order <- rownames(npx_matrix)
  batch_vector <- batch_info[match(sample_order, batch_info$SampleID)]$batch

  if(any(is.na(batch_vector))) {
    log_warn("Missing batch information for some samples")
    batch_vector[is.na(batch_vector)] <- "Unknown"
  }

  # Prepare model matrix if covariates provided
  if(!is.null(covariates)) {
    log_info("Including covariates in ComBat")
    # Align covariates with samples
    cov_aligned <- covariates[match(sample_order, covariates$SampleID), ]

    # Create model matrix (example with age and sex)
    mod <- model.matrix(~ age + sex, data = cov_aligned)
  } else {
    mod <- NULL
  }

  # Transpose matrix for ComBat (needs proteins as rows)
  expr_matrix <- t(npx_matrix)

  # Remove rows with all NAs
  complete_proteins <- rowSums(is.na(expr_matrix)) < ncol(expr_matrix)
  expr_matrix <- expr_matrix[complete_proteins, ]

  # Impute remaining NAs with row means
  for(i in 1:nrow(expr_matrix)) {
    expr_matrix[i, is.na(expr_matrix[i, ])] <- mean(expr_matrix[i, ], na.rm = TRUE)
  }

  # Initialize normalized_matrix and method to handle error case
  normalized_matrix <- NULL
  method_used <- "combat"  # Default method

  # Apply ComBat
  tryCatch({
    combat_expr <- ComBat(dat = expr_matrix,
                         batch = batch_vector,
                         mod = mod,
                         par.prior = TRUE,
                         prior.plots = FALSE)

    # Transpose back
    normalized_matrix <- t(combat_expr)

    log_info("ComBat normalisation successful")
    method_used <- "combat"

  }, error = function(e) {
    log_error("ComBat failed: {e$message}")
    log_warn("ComBat requires multiple batches to correct for batch effects")
    log_warn("  Current batch_info has only {length(unique(batch_vector))} unique batch(es)")
    log_warn("  ComBat cannot work with a single batch - this is expected when used as fallback")
    log_info("Falling back to median normalisation")

    # Fallback to median normalisation
    result <- normalize_median(npx_matrix)
    # Use <<- to assign to parent scope
    normalized_matrix <<- result$normalized_matrix
    method_used <<- "median"  # Changed from "combat" since we fell back
  })

  # Verify normalized_matrix was created (either by ComBat or fallback)
  if(is.null(normalized_matrix)) {
    log_error("ComBat failed and fallback also failed - using original matrix")
    normalized_matrix <- npx_matrix
    method_used <- "none"  # Indicates failure
  }

  return(list(
    normalized_matrix = normalized_matrix,
    batch_vector = batch_vector,
    method = method_used
  ))
}

# Master normalization function
normalize_data <- function(npx_matrix, method = "bridge",
                          bridge_samples = NULL,
                          plate_info = NULL,
                          batch_info = NULL,
                          covariates = NULL,
                          reference_batch = NULL,
                          bridge_mapping = NULL,
                          multi_batch_mode = FALSE) {

  log_info("Normalising data using method: {method}")

  result <- switch(method,
    "bridge" = {
      if (!multi_batch_mode) {
        log_error("Bridge normalization is only applicable in multi-batch (cross-batch) mode")
        stop("Bridge normalization cannot be used in single-batch mode")
      }
      if(is.null(bridge_samples)) {
        log_error("Bridge samples required for bridge normalization")
        stop("Missing bridge samples")
      }
      normalize_bridge(npx_matrix, bridge_samples, reference_batch, bridge_mapping, multi_batch_mode)
    },
    "median" = {
      normalize_median(npx_matrix, by_plate = FALSE, plate_info = plate_info)
    },
    "median_plate" = {
      if(is.null(plate_info)) {
        log_warn("Plate info not provided, using global median normalization")
        normalize_median(npx_matrix, by_plate = FALSE)
      } else {
        normalize_median(npx_matrix, by_plate = TRUE, plate_info = plate_info)
      }
    },
    "combat" = {
      if (!multi_batch_mode) {
        log_error("ComBat normalization is only applicable in multi-batch (cross-batch) mode")
        stop("ComBat normalization cannot be used in single-batch mode")
      }
      if(is.null(batch_info)) {
        log_error("Batch info required for ComBat normalization")
        stop("Missing batch info")
      }
      normalize_combat(npx_matrix, batch_info, covariates)
    },
    {
      log_error("Unknown normalization method: {method}")
      stop("Invalid normalization method")
    }
  )

  return(result)
}

# Function to evaluate normalisation
evaluate_normalization <- function(raw_matrix, normalized_matrix, metadata) {
  log_info("Evaluating normalisation performance")

  # Validate inputs
  if (is.null(normalized_matrix)) {
    log_error("normalized_matrix is NULL - cannot evaluate normalization")
    stop("normalized_matrix is NULL in evaluate_normalization()")
  }

  if (!is.matrix(normalized_matrix) && !is.data.frame(normalized_matrix)) {
    log_error("normalized_matrix is not a matrix or data.frame")
    stop("normalized_matrix must be a matrix or data.frame")
  }

  if (nrow(normalized_matrix) == 0 || ncol(normalized_matrix) == 0) {
    log_error("normalized_matrix has zero dimensions: {nrow(normalized_matrix)} x {ncol(normalized_matrix)}")
    stop("normalized_matrix has invalid dimensions")
  }

  if (nrow(normalized_matrix) != nrow(raw_matrix) || ncol(normalized_matrix) != ncol(raw_matrix)) {
    log_error("Dimension mismatch: raw_matrix {nrow(raw_matrix)}x{ncol(raw_matrix)} vs normalized_matrix {nrow(normalized_matrix)}x{ncol(normalized_matrix)}")
    stop("raw_matrix and normalized_matrix must have the same dimensions")
  }

  # For log-transformed data like NPX, use SD and MAD instead of CV
  # (CV is not meaningful when data is centered around 0)

  # Calculate per-protein SD before and after
  sd_before <- apply(raw_matrix, 2, sd, na.rm = TRUE)
  sd_after <- apply(normalized_matrix, 2, sd, na.rm = TRUE)

  # Calculate per-protein MAD (Median Absolute Deviation) before and after
  mad_before <- apply(raw_matrix, 2, mad, na.rm = TRUE)
  mad_after <- apply(normalized_matrix, 2, mad, na.rm = TRUE)

  # Calculate per-protein IQR before and after
  iqr_before <- apply(raw_matrix, 2, IQR, na.rm = TRUE)
  iqr_after <- apply(normalized_matrix, 2, IQR, na.rm = TRUE)

  # Summary statistics
  eval_stats <- data.table(
    metric = c("Mean SD", "Median SD", "SD reduction",
               "Mean MAD", "Median MAD", "MAD reduction",
               "Mean IQR", "Median IQR", "IQR reduction"),
    before = c(mean(sd_before, na.rm = TRUE),
              median(sd_before, na.rm = TRUE),
              NA,
              mean(mad_before, na.rm = TRUE),
              median(mad_before, na.rm = TRUE),
              NA,
              mean(iqr_before, na.rm = TRUE),
              median(iqr_before, na.rm = TRUE),
              NA),
    after = c(mean(sd_after, na.rm = TRUE),
             median(sd_after, na.rm = TRUE),
             mean(sd_before - sd_after, na.rm = TRUE),
             mean(mad_after, na.rm = TRUE),
             median(mad_after, na.rm = TRUE),
             mean(mad_before - mad_after, na.rm = TRUE),
             mean(iqr_after, na.rm = TRUE),
             median(iqr_after, na.rm = TRUE),
             mean(iqr_before - iqr_after, na.rm = TRUE))
  )

  log_info("Mean SD - Before: {round(eval_stats$before[1], 3)}, After: {round(eval_stats$after[1], 3)}")
  log_info("SD reduction: {round(eval_stats$after[3], 3)}")

  return(eval_stats)
}

# Function to create normalization plots
create_normalization_plots <- function(raw_matrix, normalized_matrix, title_suffix = "", n_samples = 100) {
  log_info("Creating normalisation visualisation plots")

  # Sample distributions before and after
  set.seed(123)
  n_samples_actual <- min(n_samples, nrow(raw_matrix))
  sample_idx <- sample(nrow(raw_matrix), n_samples_actual)

  # Before normalisation - convert to data.table first
  df_before <- as.data.table(raw_matrix[sample_idx, ], keep.rownames = TRUE)
  setnames(df_before, "rn", "Sample")
  df_before <- melt(df_before, id.vars = "Sample", variable.name = "Protein", value.name = "Expression")
  df_before$Stage <- "Before"

  # After normalisation - convert to data.table first
  df_after <- as.data.table(normalized_matrix[sample_idx, ], keep.rownames = TRUE)
  setnames(df_after, "rn", "Sample")
  df_after <- melt(df_after, id.vars = "Sample", variable.name = "Protein", value.name = "Expression")
  df_after$Stage <- "After"

  # Combine
  df_combined <- rbind(df_before, df_after)
  df_combined$Stage <- factor(df_combined$Stage, levels = c("Before", "After"))

  # Calculate summary statistics for before and after
  stats_before <- df_before[, .(
    Mean = round(mean(Expression, na.rm = TRUE), 2),
    Median = round(median(Expression, na.rm = TRUE), 2),
    SD = round(sd(Expression, na.rm = TRUE), 2),
    IQR = round(IQR(Expression, na.rm = TRUE), 2)
  )]

  stats_after <- df_after[, .(
    Mean = round(mean(Expression, na.rm = TRUE), 2),
    Median = round(median(Expression, na.rm = TRUE), 2),
    SD = round(sd(Expression, na.rm = TRUE), 2),
    IQR = round(IQR(Expression, na.rm = TRUE), 2)
  )]

  # Create subtitle with sample info
  total_samples <- nrow(raw_matrix)
  total_proteins <- ncol(raw_matrix)
  subtitle <- sprintf(
    "Showing %d randomly selected samples (out of %d total) × %d proteins\nBefore: Mean=%.2f, SD=%.2f, IQR=%.2f | After: Mean=%.2f, SD=%.2f, IQR=%.2f",
    n_samples_actual, total_samples, total_proteins,
    stats_before$Mean, stats_before$SD, stats_before$IQR,
    stats_after$Mean, stats_after$SD, stats_after$IQR
  )

  # Create boxplot
  p <- ggplot(df_combined, aes(x = as.factor(Sample), y = Expression, fill = Stage)) +
    geom_boxplot(alpha = 0.6, outlier.size = 0.5) +
    facet_wrap(~ Stage, scales = "free_y", ncol = 1) +
    scale_fill_manual(values = c("Before" = "#FF6B6B", "After" = "#3A5F8A")) +
    labs(title = paste("Normalization Effect", title_suffix),
         subtitle = subtitle,
         x = "Sample (randomly selected, ordered by index)",
         y = "NPX Expression Level (log2 scale)") +
    theme_bw() +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      plot.subtitle = element_text(size = 9, color = "gray40")
    )

  return(p)
}

# Main execution
main <- function() {

  # Load data
  log_info("Loading data from previous steps")
  # Use final cleaned matrix (all outliers removed) from step 05
  # Fallback to step 00 QC matrix if final cleaned matrix doesn't exist yet
  final_cleaned_path <- get_output_path("05d", "05d_npx_matrix_all_qc_passed", batch_id, "phenotypes", "rds", config = config)
  if (file.exists(final_cleaned_path)) {
    log_info("Using final cleaned matrix (all outliers removed) from step 05d")
    npx_matrix <- readRDS(final_cleaned_path)
  } else {
    log_warn("Final cleaned matrix not found: {final_cleaned_path}")
    log_warn("Falling back to step 00 QC matrix (outliers may not be removed)")
    npx_matrix_path <- get_output_path("00", "npx_matrix_qc", batch_id, "qc", config = config)
    if (!file.exists(npx_matrix_path)) {
      stop("NPX matrix file not found: {npx_matrix_path}")
    }
    npx_matrix <- readRDS(npx_matrix_path)
  }
  # Load metadata and samples_data with batch-aware paths
  metadata_path <- get_output_path("00", "metadata", batch_id, "qc", config = config)
  samples_data_path <- get_output_path("00", "samples_data_raw", batch_id, "qc", config = config)
  bridge_result_path <- get_output_path("00", "bridge_samples_identified", batch_id, "normalized", config = config)

  if (!file.exists(metadata_path)) {
    stop("Metadata file not found: {metadata_path}")
  }
  if (!file.exists(samples_data_path)) {
    stop("Samples data file not found: {samples_data_path}")
  }
  if (!file.exists(bridge_result_path)) {
    stop("Bridge result file not found: {bridge_result_path}")
  }

  metadata <- readRDS(metadata_path)
  samples_data <- readRDS(samples_data_path)
  bridge_result <- readRDS(bridge_result_path)

  # CRITICAL: Filter samples_data to only include samples present in npx_matrix
  # step 00 removed samples with >10% missingness, but samples_data_raw still contains all original samples
  # We must filter samples_data to match the samples in npx_matrix to avoid using samples that were already removed
  samples_in_matrix <- rownames(npx_matrix)
  n_samples_before <- length(unique(samples_data$SampleID))
  samples_data <- samples_data[SampleID %in% samples_in_matrix]
  n_samples_after <- length(unique(samples_data$SampleID))

  if (n_samples_before != n_samples_after) {
    log_info("Filtered samples_data to match npx_matrix:")
    log_info("  Samples in original samples_data: {n_samples_before}")
    log_info("  Samples in filtered samples_data: {n_samples_after}")
    log_info("  Samples removed (already filtered in step 00): {n_samples_before - n_samples_after}")
    log_info("  Note: These removed samples had >10% missingness and were excluded in step 00")
  } else {
    log_info("All samples in samples_data are present in npx_matrix (no filtering needed)")
  }

  # Check if multi-batch mode is enabled
  multi_batch_mode <- tryCatch(
    isTRUE(config$parameters$normalization$multi_batch_mode),
    error = function(e) FALSE
  )

  # Load batch 1 data and bridge mapping if in multi-batch mode
  batch1_npx_matrix <- NULL
  bridge_mapping <- NULL

  if (multi_batch_mode) {
    log_info("Multi-batch mode enabled: Loading batch 1 data for normalisation")

    # Load other batch NPX matrix (use reference batch's final cleaned matrix if available)
    # Get other batch ID from config (not hardcoded)
    other_batch_id <- get_other_batch_id(batch_id, config)
    if (is.null(other_batch_id)) {
      log_warn("Could not determine other batch ID for multi-batch normalization. Skipping cross-batch normalization.")
      batch1_npx_matrix <- NULL
    } else {
      batch1_npx_path <- get_output_path("05d", "05d_npx_matrix_all_qc_passed", other_batch_id, "phenotypes", "rds", config = config)
      if (!file.exists(batch1_npx_path)) {
        # Fallback to step 00 QC matrix
        batch1_npx_path <- get_output_path("00", "npx_matrix_qc", other_batch_id, "qc", config = config)
      }
      if (file.exists(batch1_npx_path)) {
        batch1_npx_matrix <- readRDS(batch1_npx_path)
        log_info("Loaded {other_batch_id} NPX matrix: {nrow(batch1_npx_matrix)} samples x {ncol(batch1_npx_matrix)} proteins")
      } else {
        log_warn("{other_batch_id} NPX matrix not found: {batch1_npx_path}")
        batch1_npx_matrix <- NULL
      }
    }

    # Load bridge mapping (from step 00)
    bridge_mapping_path <- get_output_path("00", "bridge_sample_mapping", batch_id, "normalized", config = config)
    if (file.exists(bridge_mapping_path)) {
      bridge_mapping <- readRDS(bridge_mapping_path)
      log_info("Loaded bridge mapping: {nrow(bridge_mapping)} samples")
    } else {
      log_warn("Bridge mapping file not found: {bridge_mapping_path}")
    }
  }

  # Prepare plate information
  plate_info <- unique(samples_data[, .(SampleID, PlateID)])

  # Prepare batch information for ComBat
  # In multi-batch mode, we need to distinguish between actual batches (batch_01, batch_02)
  # For ComBat to work properly across batches, we should use batch_id, not PlateID
  if (multi_batch_mode) {
    # Create batch_info with actual batch IDs
    # For current batch, all samples get the current batch_id
    # If we have reference batch data, we could combine, but for now use current batch only
    batch_info <- data.table(
      SampleID = rownames(npx_matrix),
      batch = batch_id
    )
    log_info("Prepared batch_info for ComBat: {nrow(batch_info)} samples from {batch_id}")
    
    # If we have reference batch matrix, we could combine both batches for better ComBat normalization
    # But for now, ComBat will normalize within the current batch only (as fallback when bridge fails)
    if (!is.null(batch1_npx_matrix) && other_batch_id != batch_id) {
      log_info("Note: ComBat is being applied to current batch only ({batch_id})")
      log_info("  For full cross-batch ComBat, both batches would need to be combined")
      log_info("  This is acceptable as a fallback when bridge normalization fails")
    }
  } else {
    # Single-batch mode: Use PlateID as batch (original behavior)
    batch_info <- copy(plate_info)
    setnames(batch_info, "PlateID", "batch")
  }

  # Get bridge sample IDs from current batch
  bridge_samples <- bridge_result$bridge_ids

  # In multi-batch mode, if current batch has no bridge samples, try to use reference batch's bridge samples
  # CRITICAL: Bridge samples have the SAME FINNGENIDs but DIFFERENT SampleIDs between batches
  # Solution: Map bridge samples using FINNGENIDs, not SampleIDs
  if (multi_batch_mode && length(bridge_samples) == 0 && !is.null(other_batch_id)) {
    log_warn("Current batch ({batch_id}) has no bridge samples. Attempting to use reference batch's bridge samples via FINNGENID mapping.")
    
    # Load reference batch's bridge result
    ref_bridge_result_path <- get_output_path("00", "bridge_samples_identified", other_batch_id, "normalized", config = config)
    if (file.exists(ref_bridge_result_path)) {
      ref_bridge_result <- readRDS(ref_bridge_result_path)
      ref_bridge_sampleids <- ref_bridge_result$bridge_ids
      
      if (length(ref_bridge_sampleids) > 0) {
        log_info("Found {length(ref_bridge_sampleids)} bridge samples in reference batch ({other_batch_id})")
        
        # Load sample mappings from both batches to map via FINNGENIDs
        ref_sample_mapping_path <- get_output_path("00", "sample_mapping", other_batch_id, "qc", config = config)
        current_sample_mapping_path <- get_output_path("00", "sample_mapping", batch_id, "qc", config = config)
        
        if (file.exists(ref_sample_mapping_path) && file.exists(current_sample_mapping_path)) {
          ref_sample_mapping <- readRDS(ref_sample_mapping_path)
          current_sample_mapping <- readRDS(current_sample_mapping_path)
          
          # Get FINNGENIDs of bridge samples from reference batch
          ref_bridge_mapping <- ref_sample_mapping[SampleID %in% ref_bridge_sampleids & !is.na(FINNGENID)]
          ref_bridge_finngenids <- unique(ref_bridge_mapping$FINNGENID)
          
          if (length(ref_bridge_finngenids) > 0) {
            log_info("Found {length(ref_bridge_finngenids)} unique FINNGENIDs from reference batch bridge samples")
            
            # Map FINNGENIDs to current batch's SampleIDs
            current_bridge_mapping <- current_sample_mapping[FINNGENID %in% ref_bridge_finngenids & !is.na(FINNGENID)]
            bridge_samples <- unique(current_bridge_mapping$SampleID)
            
            # Filter to only include samples present in current batch's matrix
            bridge_samples <- intersect(bridge_samples, rownames(npx_matrix))
            
            if (length(bridge_samples) > 0) {
              log_info("Mapped {length(bridge_samples)} bridge samples from reference batch to current batch using FINNGENIDs")
              log_info("  Reference batch bridge FINNGENIDs: {length(ref_bridge_finngenids)}")
              log_info("  Current batch bridge SampleIDs (mapped): {length(bridge_samples)}")
              
              # Create bridge_mapping on-the-fly for cross-batch normalization
              # CRITICAL: normalize_bridge function expects:
              # - SAMPLE_ID_batch2 = current batch (being normalized)
              # - SAMPLE_ID_batch1 = reference batch
              # This is because the function was written assuming batch_02 is always the current batch
              # So we need to map: SAMPLE_ID_batch2 (current batch) -> SAMPLE_ID_batch1 (reference batch)
              if (is.null(bridge_mapping) || nrow(bridge_mapping) == 0) {
                # Get reference batch bridge SampleIDs for the mapped FINNGENIDs
                ref_bridge_for_mapping <- ref_sample_mapping[FINNGENID %in% ref_bridge_finngenids & !is.na(FINNGENID)]
                
                # Create mapping: SAMPLE_ID_batch2 (current batch) -> SAMPLE_ID_batch1 (reference batch)
                # Note: Column names are swapped because function expects batch_02 as current batch
                bridge_mapping <- merge(
                  current_bridge_mapping[, .(SAMPLE_ID_batch2 = SampleID, FINNGENID)],  # Current batch = batch2
                  ref_bridge_for_mapping[, .(SAMPLE_ID_batch1 = SampleID, FINNGENID)],    # Reference batch = batch1
                  by = "FINNGENID",
                  all.x = TRUE
                )
                bridge_mapping[, FINNGENID := NULL]  # Remove FINNGENID column (not needed in mapping)
                bridge_mapping <- bridge_mapping[!is.na(SAMPLE_ID_batch1)]  # Only keep samples that exist in both batches
                
                if (nrow(bridge_mapping) > 0) {
                  log_info("Created bridge mapping: {nrow(bridge_mapping)} bridge samples mapped between batches")
                  log_info("  Mapping: SAMPLE_ID_batch2 ({batch_id}) -> SAMPLE_ID_batch1 ({other_batch_id})")
                } else {
                  log_warn("Bridge mapping creation resulted in 0 rows. Cross-batch normalization may fall back to single-batch mode.")
                }
              }
            } else {
              log_warn("No bridge samples found in current batch matrix after FINNGENID mapping. Bridge normalization may fail.")
              log_warn("  This may indicate that bridge samples were removed during QC steps or are not present in the data.")
            }
          } else {
            log_warn("Reference batch bridge samples have no FINNGENIDs. Cannot map to current batch.")
          }
        } else {
          log_warn("Sample mapping files not found. Cannot map bridge samples via FINNGENIDs.")
          log_warn("  Reference mapping: {ref_sample_mapping_path} (exists: {file.exists(ref_sample_mapping_path)})")
          log_warn("  Current mapping: {current_sample_mapping_path} (exists: {file.exists(current_sample_mapping_path)})")
        }
      } else {
        log_warn("Reference batch ({other_batch_id}) also has no bridge samples")
      }
    } else {
      log_warn("Reference batch bridge result file not found: {ref_bridge_result_path}")
    }
  }

  # Apply normalisation based on mode
  # Initialize variables for both modes
  normalization_method_used <- NULL
  normalized_result <- NULL
  norm_bridge <- NULL
  norm_combat <- NULL
  norm_median <- NULL
  norm_median_plate <- NULL

  if (multi_batch_mode) {
    log_info("=== MULTI-BATCH MODE: Applying normalization with fallback chain ===")
    log_info("Fallback order: Bridge (preferred) -> ComBat -> Median")

    # Step 1: Try bridge normalization first (preferred for cross-batch harmonisation)
    log_info("Attempting bridge normalization (preferred method)...")
    norm_bridge <- normalize_data(
      npx_matrix,
      method = "bridge",
      bridge_samples = bridge_samples,
      plate_info = plate_info,
      batch_info = batch_info,
      reference_batch = batch1_npx_matrix,
      bridge_mapping = bridge_mapping,
      multi_batch_mode = multi_batch_mode
    )

    if (!is.null(norm_bridge) && !is.null(norm_bridge$normalized_matrix)) {
      normalized_result <- norm_bridge
      normalization_method_used <- "bridge"
      log_info("✓ Bridge normalization succeeded - using as primary method")
      log_info("  Bridge samples used: {length(norm_bridge$bridge_samples_used %||% bridge_samples)}")
    } else {
      log_warn("Bridge normalization failed or insufficient bridge samples ({length(bridge_samples)} available)")
      log_warn("  Minimum required: 10 bridge samples")
      log_warn("  Falling back to ComBat normalization...")

      # Step 2: Try ComBat normalization (fallback 1)
      # NOTE: ComBat normalization is ONLY applicable in multi-batch (cross-batch) mode
      log_info("Attempting ComBat normalization (fallback method 1 for cross-batch correction)...")
      norm_combat <- tryCatch({
        normalize_data(npx_matrix, method = "combat", batch_info = batch_info, multi_batch_mode = TRUE)
      }, error = function(e) {
        log_warn("ComBat normalization failed: {e$message}")
        NULL
      })

      if (!is.null(norm_combat) && !is.null(norm_combat$normalized_matrix)) {
        # Check if ComBat actually succeeded or fell back to median
        if (norm_combat$method == "combat") {
          normalized_result <- norm_combat
          normalization_method_used <- "combat"
          log_info("✓ ComBat normalization succeeded - using as primary method")
        } else {
          # ComBat fell back to median, so use median as primary
          log_warn("ComBat failed and fell back to median - using median as primary method")
          normalized_result <- norm_combat
          normalization_method_used <- "median"
          log_info("✓ Median normalization (from ComBat fallback) - using as primary method")
        }
      } else {
        log_warn("ComBat normalization failed - falling back to median normalization...")

        # Step 3: Try median normalization (fallback 2)
        log_info("Attempting median normalization (fallback method 2)...")
        norm_median <- normalize_data(npx_matrix, method = "median")

        if (!is.null(norm_median) && !is.null(norm_median$normalized_matrix)) {
          normalized_result <- norm_median
          normalization_method_used <- "median"
          log_info("✓ Median normalization succeeded - using as primary method")
        } else {
          log_error("All normalization methods failed (bridge, ComBat, median)")
          stop("Normalization failed - cannot proceed")
        }
      }
    }

    # Generate alternative normalisations for comparison (only if not already computed)
    log_info("Generating alternative normalisations for comparison (multi-batch mode)")

    # Only compute methods that weren't used as primary
    # NOTE: Bridge and ComBat are ONLY applicable in multi-batch mode
    if (normalization_method_used != "bridge" && is.null(norm_bridge)) {
      norm_bridge <- tryCatch({
        normalize_data(npx_matrix, method = "bridge", bridge_samples = bridge_samples,
                      plate_info = plate_info, batch_info = batch_info,
                      reference_batch = batch1_npx_matrix, bridge_mapping = bridge_mapping,
                      multi_batch_mode = TRUE)  # Explicitly set to TRUE (we're in multi-batch mode)
      }, error = function(e) {
        log_warn("Bridge normalization (for comparison) failed: {e$message}")
        NULL
      })
    }

    if (normalization_method_used != "combat" && is.null(norm_combat)) {
      norm_combat <- tryCatch({
        normalize_data(npx_matrix, method = "combat", batch_info = batch_info, multi_batch_mode = TRUE)
      }, error = function(e) {
        log_warn("ComBat normalization (for comparison) failed: {e$message}")
        NULL
      })
    }

    if (normalization_method_used != "median" && is.null(norm_median)) {
      norm_median <- normalize_data(npx_matrix, method = "median")
    }

    # Log summary of normalization methods
    log_info("Normalization method summary:")
    log_info("  Primary method used: {normalization_method_used}")
    log_info("  Bridge available: {!is.null(norm_bridge) && !is.null(norm_bridge$normalized_matrix)}")
    log_info("  ComBat available: {!is.null(norm_combat) && !is.null(norm_combat$normalized_matrix)}")
    log_info("  Median available: {!is.null(norm_median) && !is.null(norm_median$normalized_matrix)}")

  } else {
    log_info("=== SINGLE-BATCH MODE: Applying median normalisation (standard intra-batch step) ===")
    log_info("Median normalisation ensures samples are comparable within the batch")
    log_info("Bridge normalisation and ComBat are only applicable in multi-batch mode")

    # Primary normalisation: Median (standard for single-batch)
    normalized_result <- normalize_data(npx_matrix, method = "median")

    if (is.null(normalized_result)) {
      log_error("Median normalization failed")
      stop("Normalization failed - cannot proceed")
    }

    # Set method used for single-batch mode
    normalization_method_used <- "median"
    
    # Set norm_median for consistency (it's the primary result)
    norm_median <- normalized_result

    # In single-batch mode, only median normalisation is applied
    # Bridge and ComBat don't make sense for single batch
    norm_bridge <- NULL
    norm_combat <- NULL
    norm_median_plate <- NULL

    log_info("Single-batch mode: Only median normalisation applied (no alternatives needed)")
  }

  # Evaluate normalisations based on mode
  if (multi_batch_mode) {
    # Multi-batch mode: Evaluate primary method and alternatives
    log_info("Evaluating normalisations (multi-batch mode)")
    log_info("Primary method: {normalization_method_used}")

    # Evaluate primary normalization method
    eval_primary <- evaluate_normalization(npx_matrix, normalized_result$normalized_matrix, metadata)

    # Evaluate alternative normalisations (only if they exist and weren't used as primary)
    eval_bridge <- NULL
    eval_combat <- NULL
    eval_median <- NULL

    if (!is.null(norm_bridge) && !is.null(norm_bridge$normalized_matrix)) {
      if (normalization_method_used == "bridge") {
        eval_bridge <- eval_primary
      } else {
        eval_bridge <- evaluate_normalization(npx_matrix, norm_bridge$normalized_matrix, metadata)
      }
    }

    if (!is.null(norm_combat) && !is.null(norm_combat$normalized_matrix)) {
      if (normalization_method_used == "combat") {
        eval_combat <- eval_primary
      } else {
        eval_combat <- evaluate_normalization(npx_matrix, norm_combat$normalized_matrix, metadata)
      }
    }

    if (!is.null(norm_median) && !is.null(norm_median$normalized_matrix)) {
      if (normalization_method_used == "median") {
        eval_median <- eval_primary
      } else {
        eval_median <- evaluate_normalization(npx_matrix, norm_median$normalized_matrix, metadata)
      }
    }

    eval_median_plate <- NULL  # Not used in multi-batch mode

  } else {
    # Single-batch mode: Only evaluate median normalisation
    log_info("Evaluating normalisation (single-batch mode: median only)")

    eval_bridge <- NULL  # Not applicable in single-batch mode
    eval_median <- evaluate_normalization(npx_matrix, normalized_result$normalized_matrix, metadata)
    eval_median_plate <- NULL  # Not used
    eval_combat <- NULL  # Not applicable in single-batch mode
  }

  # Combine evaluations based on mode
  if (multi_batch_mode) {
    # Multi-batch mode: Include all available methods (primary + alternatives)
    evaluation_list <- list()

    # Add primary method evaluation
    if (normalization_method_used == "bridge" && !is.null(eval_bridge)) {
      evaluation_list$bridge <- cbind(method = "bridge", eval_bridge)
    } else if (normalization_method_used == "combat" && !is.null(eval_combat)) {
      evaluation_list$combat <- cbind(method = "combat", eval_combat)
    } else if (normalization_method_used == "median" && !is.null(eval_median)) {
      evaluation_list$median <- cbind(method = "median", eval_median)
    }

    # Add alternative methods (if available and not primary)
    if (!is.null(eval_bridge) && normalization_method_used != "bridge") {
      evaluation_list$bridge <- cbind(method = "bridge", eval_bridge)
    }

    if (!is.null(eval_combat) && normalization_method_used != "combat") {
      evaluation_list$combat <- cbind(method = "combat", eval_combat)
    }

    if (!is.null(eval_median) && normalization_method_used != "median") {
      evaluation_list$median <- cbind(method = "median", eval_median)
    }

    if (length(evaluation_list) > 0) {
      all_evaluations <- do.call(rbind, evaluation_list)
    } else {
      log_warn("No normalization evaluations available - creating empty evaluation table")
      all_evaluations <- data.table(method = character(), stringsAsFactors = FALSE)
    }

  } else {
    # Single-batch mode: Only median normalisation
    all_evaluations <- cbind(method = "median", eval_median)
  }

  # Create plots based on mode
  if (multi_batch_mode) {
    # Multi-batch mode: Create plots for primary method and alternatives
    primary_label <- paste0("- ", toupper(substring(normalization_method_used, 1, 1)), substring(normalization_method_used, 2), " (Primary)")
    plot_primary <- create_normalization_plots(npx_matrix, normalized_result$normalized_matrix, primary_label)

    # Create plots for alternative methods (if available and not primary)
    plot_bridge <- NULL
    plot_combat <- NULL
    plot_median <- NULL

    if (!is.null(norm_bridge) && !is.null(norm_bridge$normalized_matrix)) {
      if (normalization_method_used == "bridge") {
        plot_bridge <- plot_primary
      } else {
        plot_bridge <- create_normalization_plots(npx_matrix, norm_bridge$normalized_matrix, "- Bridge (Comparison)")
      }
    }

    if (!is.null(norm_combat) && !is.null(norm_combat$normalized_matrix)) {
      if (normalization_method_used == "combat") {
        plot_combat <- plot_primary
      } else {
        plot_combat <- create_normalization_plots(npx_matrix, norm_combat$normalized_matrix, "- ComBat (Comparison)")
      }
    }

    if (!is.null(norm_median) && !is.null(norm_median$normalized_matrix)) {
      if (normalization_method_used == "median") {
        plot_median <- plot_primary
      } else {
        plot_median <- create_normalization_plots(npx_matrix, norm_median$normalized_matrix, "- Median (Comparison)")
      }
    }

  } else {
    # Single-batch mode: Only plot median normalisation
    plot_bridge <- NULL  # Not applicable
    plot_median <- create_normalization_plots(npx_matrix, normalized_result$normalized_matrix, "- Median")
    plot_combat <- NULL  # Not applicable
  }

  # Save outputs with batch-aware paths
  log_info("Saving normalisation results")

  # Save primary normalisation with 06_ prefix
  normalized_matrix_path <- get_output_path(step_num, "npx_matrix_normalized", batch_id, "normalized", config = config)
  normalization_result_path <- get_output_path(step_num, "normalization_result", batch_id, "normalized", config = config)
  ensure_output_dir(normalized_matrix_path)
  ensure_output_dir(normalization_result_path)

  saveRDS(normalized_result$normalized_matrix, normalized_matrix_path)
  saveRDS(normalized_result, normalization_result_path)

  method_name <- if (!is.null(normalization_method_used)) {
    normalization_method_used
  } else {
    normalized_result$method
  }
  log_info("Saved primary normalisation ({method_name}): {normalized_matrix_path}")

  # Prepare paths for alternative normalizations (only used in multi-batch mode)
  median_path <- get_output_path(step_num, "npx_matrix_normalized_median", batch_id, "normalized", config = config)
  combat_path <- get_output_path(step_num, "npx_matrix_normalized_combat", batch_id, "normalized", config = config)
  ensure_output_dir(median_path)
  ensure_output_dir(combat_path)

  # Save alternative normalisations based on mode
  if (multi_batch_mode) {
    # Multi-batch mode: Save median and ComBat for comparison
    log_info("Saving alternative normalisations (multi-batch mode)")

    if (!is.null(norm_median)) {
      saveRDS(norm_median$normalized_matrix, median_path)
      log_info("Saved median normalisation (comparison): {median_path}")
    } else {
      log_warn("Skipping median normalization save (normalization failed)")
    }

    if (!is.null(norm_combat)) {
      saveRDS(norm_combat$normalized_matrix, combat_path)
      log_info("Saved ComBat normalisation (comparison): {combat_path}")
    } else {
      log_warn("Skipping ComBat normalization save (normalization failed)")
    }

    # Median-plate not typically used in multi-batch mode
    log_info("Median-plate normalisation not saved (not used in multi-batch mode)")

  } else {
    # Single-batch mode: No alternative normalisations to save
    log_info("Single-batch mode: Only median normalisation saved (no alternatives needed)")
    log_info("Bridge and ComBat normalisations are not applicable in single-batch mode")
  }

  # Save evaluations with 06_ prefix
  evaluations_path <- get_output_path(step_num, "normalization_evaluations", batch_id, "normalized", "tsv", config = config)
  ensure_output_dir(evaluations_path)
  fwrite(all_evaluations, evaluations_path, sep = "\t")
  log_info("Saved normalisation evaluations: {evaluations_path}")

  # Prepare plot paths (mode-dependent usage)
  bridge_plot_path <- get_output_path(step_num, "normalization_effect_bridge", batch_id, "normalized", "pdf", config = config)
  median_plot_path <- get_output_path(step_num, "normalization_effect_median", batch_id, "normalized", "pdf", config = config)
  combat_plot_path <- get_output_path(step_num, "normalization_effect_combat", batch_id, "normalized", "pdf", config = config)
  ensure_output_dir(bridge_plot_path)
  ensure_output_dir(median_plot_path)
  ensure_output_dir(combat_plot_path)

  # Save plots based on mode
  if (multi_batch_mode) {
    # Multi-batch mode: Save bridge, median, and ComBat plots
    if (!is.null(plot_bridge)) {
      ggsave(bridge_plot_path, plot_bridge, width = 14, height = 9)
    }

    if (!is.null(plot_median)) {
      ggsave(median_plot_path, plot_median, width = 14, height = 9)
    } else {
      log_warn("Skipping median normalization plot (normalization failed)")
    }

    if (!is.null(plot_combat)) {
      ggsave(combat_plot_path, plot_combat, width = 14, height = 9)
    } else {
      log_warn("Skipping ComBat normalization plot (normalization failed)")
    }

  } else {
    # Single-batch mode: Only save median plot
    if (!is.null(plot_median)) {
      ggsave(median_plot_path, plot_median, width = 14, height = 9)
    }
  }

  # Print summary
  cat("\n=== NORMALIZATION SUMMARY ===\n")
  if (multi_batch_mode) {
    cat("Mode: MULTI-BATCH\n")
    cat("Primary method used:", normalization_method_used, "\n")
    cat("Fallback chain: Bridge -> ComBat -> Median\n")
    
    # Show bridge sample info if bridge was used
    if(normalization_method_used == "bridge" && !is.null(normalized_result$bridge_samples_used)) {
      cat("Bridge samples used:", length(normalized_result$bridge_samples_used), "\n")
      if (!is.null(normalized_result$bridge_samples_batch1)) {
        cat("Bridge samples in reference batch:", length(normalized_result$bridge_samples_batch1), "\n")
      }
    }
    
    # Show available alternative methods
    alt_methods <- character()
    if (!is.null(norm_bridge) && !is.null(norm_bridge$normalized_matrix) && normalization_method_used != "bridge") {
      alt_methods <- c(alt_methods, "Bridge")
    }
    if (!is.null(norm_combat) && !is.null(norm_combat$normalized_matrix) && normalization_method_used != "combat") {
      alt_methods <- c(alt_methods, "ComBat")
    }
    if (!is.null(norm_median) && !is.null(norm_median$normalized_matrix) && normalization_method_used != "median") {
      alt_methods <- c(alt_methods, "Median")
    }
    if (length(alt_methods) > 0) {
      cat("Alternative methods available:", paste(alt_methods, collapse = ", "), "\n")
    }
  } else {
    cat("Mode: SINGLE-BATCH\n")
    cat("Method: Median normalization (standard intra-batch step)\n")
    cat("Note: Bridge normalization and ComBat are only applicable in multi-batch mode\n")
  }
  cat("Matrix dimensions:", nrow(normalized_result$normalized_matrix), "x",
      ncol(normalized_result$normalized_matrix), "\n")

  cat("\n=== NORMALIZATION EVALUATION ===\n")
  print(all_evaluations)

  cat("\nNormalized data saved to: ../output/normalized/\n")

  log_info("Normalisation completed")

  # Build alternative normalizations list based on mode
  alternative_normalizations <- list()

  if (multi_batch_mode) {
    # Multi-batch mode: Include median and ComBat for comparison
    if (!is.null(norm_median)) {
      alternative_normalizations$median <- norm_median
    }
    if (!is.null(norm_combat)) {
      alternative_normalizations$combat <- norm_combat
    }
  } else {
    # Single-batch mode: No alternatives (only median is primary)
    log_info("Single-batch mode: No alternative normalisations (median is primary)")
  }

  return(list(
    normalized_matrix = normalized_result$normalized_matrix,
    normalization_result = normalized_result,
    evaluations = all_evaluations,
    alternative_normalizations = alternative_normalizations,
    mode = if (multi_batch_mode) "multi_batch" else "single_batch"
  ))
}

# Run if executed directly
if (!interactive()) {
  result <- main()
}
