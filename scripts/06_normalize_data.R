#!/usr/bin/env Rscript
# ==============================================================================
# 06_normalize_data.R - Proteomics Data Normalisation
# ==============================================================================
#
# Purpose:
#   Normalises proteomics data to remove technical variability and ensure sample
#   comparability. Mode-dependent behaviour: single-batch mode uses median normalisation
#   (standard intra-batch step), while multi-batch mode uses bridge normalisation
#   (primary) with ComBat and median for comparison. Evaluates normalisation
#   effectiveness using SD, MAD, and IQR reduction metrics.
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
normalize_bridge <- function(npx_matrix, bridge_samples, reference_batch = NULL, bridge_mapping = NULL, multi_batch_mode = FALSE) {
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
normalize_combat <- function(npx_matrix, batch_info, covariates = NULL) {
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

  # Initialize normalized_matrix to handle error case
  normalized_matrix <- NULL

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

  }, error = function(e) {
    log_error("ComBat failed: {e$message}")
    log_info("Falling back to median normalisation")

    # Fallback to median normalisation
    result <- normalize_median(npx_matrix)
    # Use <<- to assign to parent scope
    normalized_matrix <<- result$normalized_matrix
  })

  # Verify normalized_matrix was created (either by ComBat or fallback)
  if(is.null(normalized_matrix)) {
    log_error("ComBat failed and fallback also failed - using original matrix")
    normalized_matrix <- npx_matrix
  }

  return(list(
    normalized_matrix = normalized_matrix,
    batch_vector = batch_vector,
    method = "combat"
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
  final_cleaned_path <- get_output_path("05d", "npx_matrix_all_qc_passed", batch_id, "phenotypes", "rds", config = config)
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

    # Load batch 1 NPX matrix (use reference batch's final cleaned matrix if available)
    # In multi-batch mode, batch_01 is the other batch
    other_batch_id <- if (batch_id == "batch_02") "batch_01" else "batch_02"
    batch1_npx_path <- get_output_path("05d", "npx_matrix_all_qc_passed", other_batch_id, "phenotypes", "rds", config = config)
    if (!file.exists(batch1_npx_path)) {
      # Fallback to step 00 QC matrix
      batch1_npx_path <- get_output_path("00", "npx_matrix_qc", other_batch_id, "qc", config = config)
    }
    if (file.exists(batch1_npx_path)) {
      batch1_npx_matrix <- readRDS(batch1_npx_path)
      log_info("Loaded {other_batch_id} NPX matrix: {nrow(batch1_npx_matrix)} samples x {ncol(batch1_npx_matrix)} proteins")
    } else {
      log_warn("{other_batch_id} NPX matrix not found: {batch1_npx_path}")
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

  # Prepare batch information (using PlateID as batch for now)
  batch_info <- copy(plate_info)
  setnames(batch_info, "PlateID", "batch")

  # Get bridge sample IDs
  bridge_samples <- bridge_result$bridge_ids

  # Apply normalisation based on mode
  if (multi_batch_mode) {
    log_info("=== MULTI-BATCH MODE: Applying bridge normalisation (primary) ===")
    log_info("Bridge normalisation is required for cross-batch harmonisation")

    # Primary normalisation: Bridge (for cross-batch harmonisation)
    normalized_result <- normalize_data(
      npx_matrix,
      method = "bridge",
      bridge_samples = bridge_samples,
      plate_info = plate_info,
      batch_info = batch_info,
      reference_batch = batch1_npx_matrix,
      bridge_mapping = bridge_mapping,
      multi_batch_mode = multi_batch_mode
    )

    # Check if bridge normalisation succeeded
    if (is.null(normalized_result)) {
      log_error("Bridge normalization failed in multi-batch mode")
      log_error("This is required for cross-batch harmonization - cannot proceed")
      stop("Bridge normalization failed - required for multi-batch mode")
    }

    # Generate alternative normalisations for comparison (multi-batch mode)
    log_info("Generating alternative normalisations for comparison (multi-batch mode)")

    norm_median <- normalize_data(npx_matrix, method = "median")
    norm_combat <- normalize_data(npx_matrix, method = "combat", batch_info = batch_info)
    norm_median_plate <- NULL  # Not typically used in multi-batch mode

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

    # In single-batch mode, only median normalisation is applied
    # Bridge and ComBat don't make sense for single batch
    norm_median <- NULL
    norm_combat <- NULL
    norm_median_plate <- NULL

    log_info("Single-batch mode: Only median normalisation applied (no alternatives needed)")
  }

  # Evaluate normalisations based on mode
  if (multi_batch_mode) {
    # Multi-batch mode: Evaluate bridge (primary), median, and ComBat
    log_info("Evaluating normalisations (multi-batch mode)")

    eval_bridge <- evaluate_normalization(npx_matrix, normalized_result$normalized_matrix, metadata)

    # Validate alternative normalisations before evaluation
    if (is.null(norm_median) || is.null(norm_median$normalized_matrix)) {
      log_warn("Median normalization failed - skipping evaluation")
      norm_median <- NULL
      eval_median <- NULL
    } else {
      eval_median <- evaluate_normalization(npx_matrix, norm_median$normalized_matrix, metadata)
    }

    if (is.null(norm_combat) || is.null(norm_combat$normalized_matrix)) {
      log_warn("ComBat normalization failed - skipping evaluation")
      norm_combat <- NULL
      eval_combat <- NULL
    } else {
      eval_combat <- evaluate_normalization(npx_matrix, norm_combat$normalized_matrix, metadata)
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
    # Multi-batch mode: Include bridge (primary), median, and ComBat
    evaluation_list <- list(
      bridge = cbind(method = "bridge", eval_bridge)
    )

    if (!is.null(eval_median)) {
      evaluation_list$median <- cbind(method = "median", eval_median)
    }

    if (!is.null(eval_combat)) {
      evaluation_list$combat <- cbind(method = "combat", eval_combat)
    }

    all_evaluations <- do.call(rbind, evaluation_list)

  } else {
    # Single-batch mode: Only median normalisation
    all_evaluations <- cbind(method = "median", eval_median)
  }

  # Create plots based on mode
  if (multi_batch_mode) {
    # Multi-batch mode: Create plots for bridge (primary), median, and ComBat
    plot_bridge <- create_normalization_plots(npx_matrix, normalized_result$normalized_matrix, "- Bridge (Primary)")

    plot_median <- if (!is.null(norm_median)) {
      create_normalization_plots(npx_matrix, norm_median$normalized_matrix, "- Median (Comparison)")
    } else {
      NULL
    }

    plot_combat <- if (!is.null(norm_combat)) {
      create_normalization_plots(npx_matrix, norm_combat$normalized_matrix, "- ComBat (Comparison)")
    } else {
      NULL
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

  log_info("Saved primary normalisation ({normalized_result$method}): {normalized_matrix_path}")

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
    cat("Primary method:", normalized_result$method, "\n")
    if(normalized_result$method == "bridge" || normalized_result$method == "bridge_cross_batch") {
      cat("Bridge samples used:", length(normalized_result$bridge_samples_used), "\n")
      if (!is.null(normalized_result$bridge_samples_batch1)) {
        cat("Bridge samples in reference batch:", length(normalized_result$bridge_samples_batch1), "\n")
      }
    }
    cat("Alternative methods:", if (!is.null(norm_median)) "median" else "",
        if (!is.null(norm_combat)) "ComBat" else "", "\n")
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
