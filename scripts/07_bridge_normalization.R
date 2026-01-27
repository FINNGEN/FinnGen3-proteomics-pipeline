#!/usr/bin/env Rscript
# ==============================================================================
# 07_bridge_normalization.R - Enhanced Bridge Sample Normalisation
# ==============================================================================
#
# Purpose:
#   Enhanced bridge sample normalisation for multi-batch harmonisation. Provides
#   median-based and quantile-based normalisation methods using bridge samples
#   from multiple batches. Optional step that is only used in multi-batch mode
#   for cross-batch data integration.
#
# Author: Reza Jabal, PhD (rjabal@broadinstitute.org)
# Date: December 2025
# ==============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(tidyverse)
  library(OlinkAnalyze)
  library(ggplot2)
  library(yaml)
  library(logger)
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
config_file <- Sys.getenv("PIPELINE_CONFIG", "")
if (config_file == "" || !file.exists(config_file)) {
  stop("PIPELINE_CONFIG environment variable not set or config file not found. Please provide path to config file.")
}
config <- read_yaml(config_file)

# Set up logging with batch-aware path
log_path <- get_log_path("09_bridge_normalization", batch_id, config = config)
log_appender(appender_file(log_path))
log_info("Starting enhanced bridge normalisation for batch: {batch_id}")

# Set theme for plots
theme_set(theme_bw())

# Function to validate bridge samples for normalisation
validate_bridge_samples <- function(npx_matrix, bridge_samples) {
  log_info("Validating bridge samples for normalisation")

  # Check bridge sample availability
  available_bridges <- intersect(rownames(npx_matrix), bridge_samples)
  missing_bridges <- setdiff(bridge_samples, rownames(npx_matrix))

  log_info("Available bridge samples: {length(available_bridges)}")
  if(length(missing_bridges) > 0) {
    log_warn("Missing bridge samples: {length(missing_bridges)}")
  }

  # Extract bridge sample data
  bridge_data <- npx_matrix[available_bridges, ]

  # Check quality of bridge samples
  bridge_quality <- data.table(
    SampleID = available_bridges,
    missing_rate = rowSums(is.na(bridge_data)) / ncol(bridge_data),
    mean_expr = rowMeans(bridge_data, na.rm = TRUE),
    sd_expr = apply(bridge_data, 1, sd, na.rm = TRUE)
  )

  # Flag low-quality bridge samples
  bridge_quality[, high_quality := missing_rate < 0.2]

  high_quality_bridges <- bridge_quality[high_quality == TRUE]$SampleID

  log_info("High-quality bridge samples: {length(high_quality_bridges)}")

  return(list(
    available_bridges = available_bridges,
    high_quality_bridges = high_quality_bridges,
    bridge_quality = bridge_quality
  ))
}

# Function for enhanced bridge normalisation
enhanced_bridge_normalization <- function(npx_matrix, bridge_samples, method = "median") {
  log_info("Performing enhanced bridge normalisation with method: {method}")

  # Validate bridge samples
  validation <- validate_bridge_samples(npx_matrix, bridge_samples)

  if(length(validation$high_quality_bridges) < 10) {
    log_error("Insufficient high-quality bridge samples for normalisation")
    return(NULL)
  }

  # Use high-quality bridge samples
  bridge_data <- npx_matrix[validation$high_quality_bridges, ]

  if(method == "median") {
    # Median-based normalisation
    bridge_medians <- apply(bridge_data, 2, median, na.rm = TRUE)
    global_median <- median(bridge_medians, na.rm = TRUE)
    scaling_factors <- global_median / bridge_medians

  } else if(method == "mean") {
    # Mean-based normalisation
    bridge_means <- apply(bridge_data, 2, mean, na.rm = TRUE)
    global_mean <- mean(bridge_means, na.rm = TRUE)
    scaling_factors <- global_mean / bridge_means

  } else if(method == "quantile") {
    # Quantile normalisation using bridge samples
    log_info("Performing quantile normalisation")

    # Get reference distribution from bridge samples
    bridge_quantiles <- apply(bridge_data, 2, function(x) {
      quantile(x, probs = seq(0, 1, 0.01), na.rm = TRUE)
    })
    reference_quantiles <- rowMeans(bridge_quantiles, na.rm = TRUE)

    # Apply quantile normalisation to full matrix
    normalized_matrix <- apply(npx_matrix, 2, function(x) {
      x_ranks <- rank(x, na.last = "keep", ties.method = "average")
      x_probs <- x_ranks / (sum(!is.na(x)) + 1)
      x_normalized <- approx(seq(0, 1, 0.01), reference_quantiles, xout = x_probs)$y
      return(x_normalized)
    })

    return(list(
      normalized_matrix = normalized_matrix,
      method = "quantile",
      bridge_samples_used = validation$high_quality_bridges,
      reference_quantiles = reference_quantiles
    ))
  }

  # Apply scaling factors for median/mean methods
  if(method %in% c("median", "mean")) {
    scaling_factors[is.na(scaling_factors) | is.infinite(scaling_factors)] <- 1

    # Cap extreme scaling factors to prevent outliers (keep within 0.1 to 10 range)
    scaling_factors[scaling_factors < 0.1] <- 0.1
    scaling_factors[scaling_factors > 10] <- 10

    normalized_matrix <- sweep(npx_matrix, 2, scaling_factors, "*")

    return(list(
      normalized_matrix = normalized_matrix,
      scaling_factors = scaling_factors,
      method = method,
      bridge_samples_used = validation$high_quality_bridges
    ))
  }
}

# Function to harmonise with reference batch
harmonize_with_reference <- function(current_batch, reference_batch = NULL, bridge_samples) {
  log_info("Harmonising with reference batch")

  if(is.null(reference_batch)) {
    log_info("No reference batch provided, using internal normalisation")
    return(enhanced_bridge_normalization(current_batch, bridge_samples, method = "median"))
  }

  # This function would be implemented when batch 1 data is available
  # For now, return placeholder
  log_info("Reference batch harmonization to be implemented")

  return(list(
    harmonized_matrix = current_batch,
    message = "Reference batch harmonization pending implementation"
  ))
}

# Function to evaluate normalisation quality
evaluate_bridge_normalization <- function(raw_matrix, normalized_matrix, bridge_samples) {
  log_info("Evaluating bridge normalisation quality")

  # Get bridge sample indices
  bridge_idx <- which(rownames(raw_matrix) %in% bridge_samples)
  non_bridge_idx <- which(!rownames(raw_matrix) %in% bridge_samples)

  # Calculate CV before and after for bridge samples
  cv_bridge_before <- apply(raw_matrix[bridge_idx, ], 2, function(x) {
    sd(x, na.rm = TRUE) / mean(x, na.rm = TRUE)
  })

  cv_bridge_after <- apply(normalized_matrix[bridge_idx, ], 2, function(x) {
    sd(x, na.rm = TRUE) / mean(x, na.rm = TRUE)
  })

  # Calculate CV for non-bridge samples
  cv_non_bridge_before <- apply(raw_matrix[non_bridge_idx, ], 2, function(x) {
    sd(x, na.rm = TRUE) / mean(x, na.rm = TRUE)
  })

  cv_non_bridge_after <- apply(normalized_matrix[non_bridge_idx, ], 2, function(x) {
    sd(x, na.rm = TRUE) / mean(x, na.rm = TRUE)
  })

  # Summary statistics
  eval_stats <- data.table(
    sample_type = c("Bridge", "Bridge", "Non-Bridge", "Non-Bridge"),
    stage = c("Before", "After", "Before", "After"),
    mean_cv = c(mean(cv_bridge_before, na.rm = TRUE),
                mean(cv_bridge_after, na.rm = TRUE),
                mean(cv_non_bridge_before, na.rm = TRUE),
                mean(cv_non_bridge_after, na.rm = TRUE)),
    median_cv = c(median(cv_bridge_before, na.rm = TRUE),
                  median(cv_bridge_after, na.rm = TRUE),
                  median(cv_non_bridge_before, na.rm = TRUE),
                  median(cv_non_bridge_after, na.rm = TRUE))
  )

  # Calculate improvement
  bridge_improvement <- mean(cv_bridge_before - cv_bridge_after, na.rm = TRUE)
  overall_improvement <- mean(c(cv_bridge_before, cv_non_bridge_before) -
                             c(cv_bridge_after, cv_non_bridge_after), na.rm = TRUE)

  log_info("Bridge sample CV improvement: {round(bridge_improvement, 4)}")
  log_info("Overall CV improvement: {round(overall_improvement, 4)}")

  return(list(
    eval_stats = eval_stats,
    bridge_improvement = bridge_improvement,
    overall_improvement = overall_improvement
  ))
}

# Function to create bridge normalisation plots
create_bridge_plots <- function(raw_matrix, normalized_matrix, bridge_samples) {
  log_info("Creating bridge normalisation visualisation")

  # Get bridge and non-bridge samples
  bridge_idx <- rownames(raw_matrix) %in% bridge_samples

  # Sample data for plotting
  set.seed(123)
  sample_proteins <- sample(ncol(raw_matrix), min(100, ncol(raw_matrix)))

  # Prepare data
  plot_data <- rbind(
    data.table(
      value = as.vector(raw_matrix[, sample_proteins]),
      stage = "Before",
      is_bridge = rep(bridge_idx, length(sample_proteins))
    ),
    data.table(
      value = as.vector(normalized_matrix[, sample_proteins]),
      stage = "After",
      is_bridge = rep(bridge_idx, length(sample_proteins))
    )
  )

  plot_data[, sample_type := ifelse(is_bridge, "Bridge", "Non-Bridge")]

  # Create boxplot
  p1 <- ggplot(plot_data[!is.na(value)], aes(x = sample_type, y = value, fill = stage)) +
    geom_boxplot(alpha = 0.7, outlier.size = 0.5) +
    scale_fill_manual(values = c("Before" = "#FF6B6B", "After" = "#3A5F8A")) +
    labs(title = "Bridge Normalization Effect",
         x = "Sample Type", y = "Expression Value",
         fill = "Stage") +
    theme_bw() +
    facet_wrap(~ stage, scales = "free_y")

  # Density plot
  p2 <- ggplot(plot_data[!is.na(value)], aes(x = value, fill = sample_type)) +
    geom_density(alpha = 0.5) +
    scale_fill_manual(values = c("Bridge" = "#5C9EAD", "Non-Bridge" = "#87CEEB")) +
    labs(title = "Expression Distribution",
         x = "Expression Value", y = "Density",
         fill = "Sample Type") +
    theme_bw() +
    facet_wrap(~ stage, scales = "free")

  return(list(
    boxplot = p1,
    density = p2
  ))
}

# Main execution
main <- function() {

  # Check if multi-batch mode is enabled
  multi_batch_mode <- tryCatch(
    isTRUE(config$parameters$normalization$multi_batch_mode),
    error = function(e) FALSE
  )

  # Check if enhanced bridge normalisation should be run
  run_enhanced_bridge <- tryCatch(
    isTRUE(config$parameters$normalization$run_enhanced_bridge),
    error = function(e) FALSE
  )

  # Auto-enable in multi-batch mode if not explicitly disabled
  if (multi_batch_mode && !run_enhanced_bridge) {
    log_info("Multi-batch mode detected: Auto-enabling enhanced bridge normalization")
    run_enhanced_bridge <- TRUE
  }

  if(!run_enhanced_bridge) {
    # Set environment variable to indicate step was skipped
    Sys.setenv(PIPELINE_STEP_SKIPPED = "TRUE")
    log_info("Enhanced bridge normalization is disabled in config (run_enhanced_bridge: false)")
    log_info("Skipping step 07 - Enhanced bridge normalization")
    log_info("This step is only needed for multi-batch integration (FG2+FG3 or Batch1+Batch2)")
    cat("\n=== ENHANCED BRIDGE NORMALIZATION SKIPPED ===\n")
    cat("Reason: Single-batch analysis (config: run_enhanced_bridge = false)\n")
    cat("To enable: Set run_enhanced_bridge: true in batch2_config.yaml\n")
    cat("Use case: Multi-batch data integration only\n\n")
    return(NULL)
  }

  # Load data from previous steps
  log_info("Loading data from previous steps")

  # Determine input based on configuration
  input_choice <- tryCatch(config$parameters$normalization$enhanced_bridge_input, error = function(e) "raw")

  if(input_choice == "raw") {
    log_info("Using RAW NPX data as input (npx_matrix_qc.rds)")
    npx_matrix_path <- get_output_path("00", "npx_matrix_qc", batch_id, "qc", config = config)
    if (!file.exists(npx_matrix_path)) {
      log_error("NPX matrix file not found: {npx_matrix_path}")
      stop(paste0("NPX matrix file not found: ", npx_matrix_path))
    }
    npx_matrix <- readRDS(npx_matrix_path)
  } else if(input_choice == "step08_median") {
    log_info("Using step 06 MEDIAN normalized data as input (08_npx_matrix_normalized_median.rds)")
    npx_matrix_path <- get_output_path("06", "npx_matrix_normalized_median", batch_id, "normalized", config = config)
    if (!file.exists(npx_matrix_path)) {
      log_error("step 06 normalized matrix file not found: {npx_matrix_path}")
      stop(paste0("step 06 normalized matrix file not found: ", npx_matrix_path))
    }
    npx_matrix <- readRDS(npx_matrix_path)
  } else {
    log_error("Invalid enhanced_bridge_input: {input_choice}. Must be 'raw' or 'step08_median'")
    stop(paste0("Invalid configuration: enhanced_bridge_input must be 'raw' or 'step08_median', got: ", input_choice))
  }

  bridge_result_path <- get_output_path("00", "bridge_samples_identified", batch_id, "normalized", config = config)
  if (!file.exists(bridge_result_path)) {
    stop("Bridge result file not found: {bridge_result_path}")
  }
  bridge_result <- readRDS(bridge_result_path)

  # Get bridge sample IDs
  bridge_samples <- bridge_result$bridge_ids

  # Perform enhanced bridge normalisation
  norm_result <- enhanced_bridge_normalization(
    npx_matrix,
    bridge_samples,
    method = "median"  # Can be "median", "mean", or "quantile"
  )

  # Also try quantile normalisation
  norm_quantile <- enhanced_bridge_normalization(
    npx_matrix,
    bridge_samples,
    method = "quantile"
  )

  # Evaluate normalisation
  eval_median <- evaluate_bridge_normalization(npx_matrix, norm_result$normalized_matrix, bridge_samples)
  eval_quantile <- evaluate_bridge_normalization(npx_matrix, norm_quantile$normalized_matrix, bridge_samples)

  # Compare with step 06 normalisation if available (for multi-batch mode)
  step08_cv_reduction <- NULL
  if (multi_batch_mode) {
    log_info("Comparing with step 06 normalization for CV reduction")

    # Try to load step 06 normalised matrix
    step08_file <- get_output_path("06", "npx_matrix_normalized", batch_id, "normalized", config = config)
    if (file.exists(step08_file)) {
      step08_matrix <- readRDS(step08_file)
      eval_step08 <- evaluate_bridge_normalization(npx_matrix, step08_matrix, bridge_samples)
      step08_cv_reduction <- eval_step08$overall_improvement
      log_info("step 06 CV reduction: {round(step08_cv_reduction, 4)}")
    } else {
      log_warn("step 06 normalized matrix not found, skipping comparison")
    }
  }

  # Select best method based on CV reduction
  median_cv_reduction <- eval_median$overall_improvement
  quantile_cv_reduction <- eval_quantile$overall_improvement

  best_method <- "median"
  best_cv_reduction <- median_cv_reduction
  best_result <- norm_result

  if (quantile_cv_reduction > median_cv_reduction) {
    best_method <- "quantile"
    best_cv_reduction <- quantile_cv_reduction
    best_result <- norm_quantile
  }

  if (!is.null(step08_cv_reduction) && step08_cv_reduction > best_cv_reduction) {
    log_info("step 06 normalization has better CV reduction than step 07 methods")
    log_info("Recommendation: Use step 06 normalized matrix for downstream analysis")
    best_method <- "step08"
    best_cv_reduction <- step08_cv_reduction
  } else {
    log_info("step 07 {best_method} normalization has best CV reduction: {round(best_cv_reduction, 4)}")
  }

  # Create plots
  plots_median <- create_bridge_plots(npx_matrix, norm_result$normalized_matrix, bridge_samples)
  plots_quantile <- create_bridge_plots(npx_matrix, norm_quantile$normalized_matrix, bridge_samples)

  # Save outputs with 09_ prefix using batch-aware paths
  log_info("Saving enhanced bridge normalization results")

  norm_result_path <- get_output_path(step_num, "bridge_norm_result", batch_id, "normalized", config = config)
  norm_quantile_path <- get_output_path(step_num, "bridge_norm_quantile", batch_id, "normalized", config = config)
  enhanced_matrix_path <- get_output_path(step_num, "npx_matrix_bridge_enhanced", batch_id, "normalized", config = config)
  quantile_matrix_path <- get_output_path(step_num, "npx_matrix_bridge_quantile", batch_id, "normalized", config = config)

  ensure_output_dir(norm_result_path)
  ensure_output_dir(norm_quantile_path)
  ensure_output_dir(enhanced_matrix_path)
  ensure_output_dir(quantile_matrix_path)

  saveRDS(norm_result, norm_result_path)
  saveRDS(norm_quantile, norm_quantile_path)
  saveRDS(norm_result$normalized_matrix, enhanced_matrix_path)
  saveRDS(norm_quantile$normalized_matrix, quantile_matrix_path)

  # Save evaluation results with 09_ prefix
  eval_median_path <- get_output_path(step_num, "bridge_norm_eval_median", batch_id, "normalized", "tsv", config = config)
  eval_quantile_path <- get_output_path(step_num, "bridge_norm_eval_quantile", batch_id, "normalized", "tsv", config = config)
  ensure_output_dir(eval_median_path)
  ensure_output_dir(eval_quantile_path)

  fwrite(eval_median$eval_stats, eval_median_path, sep = "\t")
  fwrite(eval_quantile$eval_stats, eval_quantile_path, sep = "\t")

  # Save plots with 09_ prefix
  boxplot_median_path <- get_output_path(step_num, "bridge_norm_boxplot_median", batch_id, "normalized", "pdf", config = config)
  density_median_path <- get_output_path(step_num, "bridge_norm_density_median", batch_id, "normalized", "pdf", config = config)
  boxplot_quantile_path <- get_output_path(step_num, "bridge_norm_boxplot_quantile", batch_id, "normalized", "pdf", config = config)
  density_quantile_path <- get_output_path(step_num, "bridge_norm_density_quantile", batch_id, "normalized", "pdf", config = config)
  ensure_output_dir(boxplot_median_path)
  ensure_output_dir(density_median_path)
  ensure_output_dir(boxplot_quantile_path)
  ensure_output_dir(density_quantile_path)

  ggsave(boxplot_median_path, plots_median$boxplot, width = 10, height = 6)
  ggsave(density_median_path, plots_median$density, width = 10, height = 6)
  ggsave(boxplot_quantile_path, plots_quantile$boxplot, width = 10, height = 6)
  ggsave(density_quantile_path, plots_quantile$density, width = 10, height = 6)

  # Print summary
  cat("\n=== ENHANCED BRIDGE NORMALIZATION SUMMARY ===\n")
  cat("Bridge samples used:", length(norm_result$bridge_samples_used), "\n")
  cat("\nMedian normalization:\n")
  cat("  - Bridge CV improvement:", round(eval_median$bridge_improvement, 4), "\n")
  cat("  - Overall CV improvement:", round(eval_median$overall_improvement, 4), "\n")
  cat("\nQuantile normalization:\n")
  cat("  - Bridge CV improvement:", round(eval_quantile$bridge_improvement, 4), "\n")
  cat("  - Overall CV improvement:", round(eval_quantile$overall_improvement, 4), "\n")
  if (!is.null(step08_cv_reduction)) {
    cat("\nStep 08 normalization:\n")
    cat("  - Overall CV improvement:", round(step08_cv_reduction, 4), "\n")
  }
  cat("\n=== BEST METHOD SELECTION ===\n")
  cat("Best method:", best_method, "\n")
  cat("Best CV reduction:", round(best_cv_reduction, 4), "\n")
  if (best_method == "step08") {
    cat("Recommendation: Use 08_npx_matrix_normalized.rds for downstream analysis\n")
  } else {
    cat("Recommendation: Use 09_npx_matrix_bridge_{best_method}.rds for downstream analysis\n")
  }
  cat("\nResults saved to: ../output/normalized/\n")

  log_info("Enhanced bridge normalization completed")

  return(list(
    norm_result = norm_result,
    norm_quantile = norm_quantile,
    evaluation = list(median = eval_median, quantile = eval_quantile)
  ))
}

# Run if executed directly
if (!interactive()) {
  result <- main()
}
