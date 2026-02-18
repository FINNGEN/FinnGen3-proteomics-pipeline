#!/usr/bin/env Rscript
# ==============================================================================
# 11_rank_normalize.R - Inverse Rank Normalisation
# ==============================================================================
#
# Purpose:
#   Applies inverse rank normalisation (IRN) to the cleaned NPX matrix column-wise
#   (per protein). Transforms data to a normal distribution suitable for downstream
#   statistical analysis. Generates PLINK format phenotype files and distribution
#   visualisations comparing before and after normalisation.
#
# Author: Reza Jabal, PhD (rjabal@broadinstitute.org)
# Date: December 2025 (Updated: January 2026 - Fixed file path issues)
# ==============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(tidyverse)
  library(ggplot2)
  library(yaml)
  library(logger)
})

# Suppress linting warnings
utils::globalVariables(c("theoretical"))

# Source path utilities for batch-aware paths
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

# Load configuration first
config_file <- Sys.getenv("PIPELINE_CONFIG", "")
if (config_file == "" || !file.exists(config_file)) {
  stop("PIPELINE_CONFIG environment variable not set or config file not found. Please provide path to config file.")
}
config <- read_yaml(config_file)

# Get batch context
batch_id <- Sys.getenv("PIPELINE_BATCH_ID", config$batch$default_batch_id %||% "batch_01")
step_num <- get_step_number()

# Set up logging with batch-aware path
log_path <- get_log_path(step_num, batch_id, config = config)
ensure_output_dir(log_path)
log_appender(appender_file(log_path))
log_info("Starting rank normalisation for batch: {batch_id}")

# Set theme for plots
theme_set(theme_bw())

# Function for inverse rank normalisation
inverse_rank_normalize <- function(x) {
  # Remove NA values
  x_clean <- x[!is.na(x)]
  n <- length(x_clean)

  if(n < 3) {
    # Too few values to normalise
    return(x)
  }

  # Rank the values
  ranks <- rank(x_clean, ties.method = "average")

  # Convert to quantiles
  quantiles <- (ranks - 0.5) / n

  # Convert to normal distribution
  normalized <- qnorm(quantiles)

  # Put back in original positions
  result <- x
  result[!is.na(x)] <- normalized

  return(result)
}

# Function to apply rank normalisation to matrix
rank_normalize_matrix <- function(phenotype_matrix, by = "column") {
  log_info("Applying inverse rank normalisation by {by}")

  if(by == "column") {
    # Normalise each protein separately (column-wise)
    # Use sapply with simplify=FALSE then combine to preserve matrix structure
    normalized_list <- lapply(seq_len(ncol(phenotype_matrix)), function(i) {
      inverse_rank_normalize(phenotype_matrix[, i])
    })
    normalized_matrix <- do.call(cbind, normalized_list)

  } else if(by == "row") {
    # Normalise each sample separately (row-wise)
    normalized_list <- lapply(seq_len(nrow(phenotype_matrix)), function(i) {
      inverse_rank_normalize(phenotype_matrix[i, ])
    })
    normalized_matrix <- do.call(rbind, normalized_list)

  } else {
    stop("Invalid normalization direction. Use 'column' or 'row'")
  }

  # Preserve row and column names
  rownames(normalized_matrix) <- rownames(phenotype_matrix)
  colnames(normalized_matrix) <- colnames(phenotype_matrix)

  return(normalized_matrix)
}

# Function to residualise after rank normalisation
residualize_covariates <- function(normalized_matrix, covariates) {
  log_info("Residualising for covariates after rank normalisation")

  # This would typically be done in the GWAS software
  # Here we provide the option for pre-residualization if needed

  residualized_matrix <- normalized_matrix

  # For each protein
  for(i in seq_len(ncol(normalized_matrix))) {
    y <- normalized_matrix[, i]

    # Check for complete cases
    complete_idx <- complete.cases(cbind(y, covariates))

    if(sum(complete_idx) < 50) {
      next
    }

    # Fit linear model
    fit_data <- data.frame(y = y[complete_idx], covariates[complete_idx, ])

    tryCatch({
      model <- lm(y ~ ., data = fit_data)
      residuals <- residuals(model)

      # Rank normalise the residuals
      residuals_normalized <- inverse_rank_normalize(residuals)

      residualized_matrix[complete_idx, i] <- residuals_normalized

    }, error = function(e) {
      log_debug("Failed to residualize protein {i}: {e$message}")
    })
  }

  return(residualized_matrix)
}

# Function to compare distributions
compare_distributions <- function(original_matrix, normalized_matrix) {
  log_info("Comparing distributions before and after normalisation")

  # Sample proteins for comparison
  set.seed(123)
  sample_proteins <- sample(ncol(original_matrix), min(10, ncol(original_matrix)))

  comparison_stats <- data.table()

  for(i in sample_proteins) {
    protein_name <- colnames(original_matrix)[i]

    # Original distribution stats
    orig_values <- original_matrix[, i]
    orig_values <- orig_values[!is.na(orig_values)]

    # Normalised distribution stats
    norm_values <- normalized_matrix[, i]
    norm_values <- norm_values[!is.na(norm_values)]

    # Shapiro-Wilk test for normality
    shapiro_orig <- if(length(orig_values) > 3 && length(orig_values) < 5000) {
      shapiro.test(orig_values)$p.value
    } else { NA }

    shapiro_norm <- if(length(norm_values) > 3 && length(norm_values) < 5000) {
      shapiro.test(norm_values)$p.value
    } else { NA }

    stats <- data.table(
      protein = protein_name,
      mean_orig = mean(orig_values),
      sd_orig = sd(orig_values),
      skew_orig = moments::skewness(orig_values),
      kurt_orig = moments::kurtosis(orig_values),
      shapiro_p_orig = shapiro_orig,
      mean_norm = mean(norm_values),
      sd_norm = sd(norm_values),
      skew_norm = moments::skewness(norm_values),
      kurt_norm = moments::kurtosis(norm_values),
      shapiro_p_norm = shapiro_norm
    )

    comparison_stats <- rbind(comparison_stats, stats)
  }

  return(comparison_stats)
}

# Function to create distribution plots
create_distribution_plots <- function(original_matrix, normalized_matrix) {
  log_info("Creating distribution visualisation plots")

  # Sample proteins for plotting
  set.seed(123)
  sample_proteins <- sample(ncol(original_matrix), min(4, ncol(original_matrix)))

  plot_list <- list()

  for(i in seq_along(sample_proteins)) {
    protein_idx <- sample_proteins[i]
    protein_name <- colnames(original_matrix)[protein_idx]

    # Prepare data
    plot_data <- data.table(
      original = original_matrix[, protein_idx],
      normalized = normalized_matrix[, protein_idx]
    )

    plot_data_long <- melt(plot_data, measure.vars = c("original", "normalized"),
                          variable.name = "stage", value.name = "expression")

    # Create plot
    p <- ggplot(plot_data_long[!is.na(expression)], aes(x = expression, fill = stage)) +
      geom_histogram(aes(y = after_stat(density)), bins = 30, alpha = 0.6, position = "identity") +
      geom_density(alpha = 0.3) +
      scale_fill_manual(values = c("original" = "#FF6B6B", "normalised" = "#3A5F8A"),
                       labels = c("Original", "Rank Normalized")) +
      labs(title = paste("Protein:", substr(protein_name, 1, 20)),
           x = "Expression", y = "Density",
           fill = "Stage") +
      theme_bw() +
      facet_wrap(~ stage, scales = "free")

    plot_list[[i]] <- p
  }

  # Combine plots
  if(length(plot_list) >= 4) {
    combined_plot <- gridExtra::grid.arrange(
      plot_list[[1]], plot_list[[2]],
      plot_list[[3]], plot_list[[4]],
      ncol = 2
    )
  } else {
    combined_plot <- plot_list[[1]]
  }

  # QQ plot for checking normality
  sample_protein <- sample(ncol(normalized_matrix), 1)
  norm_values <- normalized_matrix[!is.na(normalized_matrix[, sample_protein]), sample_protein]

  qq_data <- data.table(
    theoretical = qnorm(ppoints(length(norm_values))),
    sample = sort(norm_values)
  )

  qq_plot <- ggplot(qq_data, aes(x = theoretical, y = sample)) +
    geom_point(alpha = 0.5, color = "#3A5F8A") +
    geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
    labs(title = "Q-Q Plot After Rank Normalisation",
         x = "Theoretical Quantiles",
         y = "Sample Quantiles") +
    theme_bw()

  return(list(
    distribution_plots = plot_list,
    qq_plot = qq_plot
  ))
}

# Main execution
main <- function() {

  # Check if aggregation is enabled
  multi_batch_mode <- tryCatch(
    isTRUE(config$parameters$normalization$multi_batch_mode),
    error = function(e) FALSE
  )

  aggregate_output <- tryCatch(
    isTRUE(config$parameters$aggregation$aggregate_output),
    error = function(e) FALSE
  )

  # Load data from previous steps using batch-aware paths
  log_info("Loading data from previous steps")

  # Check if batch correction is enabled (need to do this BEFORE loading data)
  adjust_for_batch <- tryCatch(
    isTRUE(config$parameters$covariate_adjustment$adjust_for_batch),
    error = function(e) FALSE
  )

  # Try to load batch-corrected data first if batch correction is enabled
  phenotype_unrelated_path <- NULL
  finngenid_unrelated_path <- NULL

  if (adjust_for_batch && multi_batch_mode) {
    phenotype_unrelated_path_bc <- get_output_path("10", "phenotype_matrix_unrelated_batch_corrected", batch_id, "phenotypes", config = config)
    finngenid_unrelated_path_bc <- get_output_path("10", "phenotype_matrix_finngenid_unrelated_batch_corrected", batch_id, "phenotypes", config = config)

    if (file.exists(finngenid_unrelated_path_bc)) {
      phenotype_unrelated_path <- phenotype_unrelated_path_bc
      finngenid_unrelated_path <- finngenid_unrelated_path_bc
      log_info("Batch correction enabled - using BATCH-CORRECTED kinship-filtered data")
    } else {
      log_warn("Batch correction enabled but batch-corrected data not found: {finngenid_unrelated_path_bc}")
      log_warn("Falling back to non-batch-corrected data")
    }
  }

  # Fallback to non-batch-corrected data
  if (is.null(phenotype_unrelated_path)) {
    phenotype_unrelated_path <- get_output_path("10", "phenotype_matrix_unrelated", batch_id, "phenotypes", config = config)
    finngenid_unrelated_path <- get_output_path("10", "phenotype_matrix_finngenid_unrelated", batch_id, "phenotypes", config = config)
  }

  # When batch correction is enabled, prioritize FINNGENID matrix over SampleID matrix
  if (adjust_for_batch && multi_batch_mode && file.exists(finngenid_unrelated_path)) {
    # Load FINNGENID matrix first
    finngenid_unrelated <- readRDS(finngenid_unrelated_path)
    log_info("Loaded BATCH-CORRECTED FINNGENID-indexed matrix: {nrow(finngenid_unrelated)} samples x {ncol(finngenid_unrelated)} proteins")

    # Create a dummy SampleID matrix with same structure (won't be used for aggregation)
    phenotype_unrelated <- finngenid_unrelated
    rownames(phenotype_unrelated) <- rownames(finngenid_unrelated)
    log_info("Using FINNGENID matrix as primary data source (batch-corrected)")
  } else {
    # Load SampleID matrix (standard path)
    if (!file.exists(phenotype_unrelated_path)) {
      stop(paste0("Unrelated phenotype matrix not found: ", phenotype_unrelated_path, ". Please run step 10 first."))
    }

    phenotype_unrelated <- readRDS(phenotype_unrelated_path)
    log_info("Loaded unrelated phenotype matrix: {nrow(phenotype_unrelated)} samples x {ncol(phenotype_unrelated)} proteins")

    finngenid_unrelated <- NULL
    if(file.exists(finngenid_unrelated_path)) {
      finngenid_unrelated <- readRDS(finngenid_unrelated_path)
      log_info("Loaded unrelated FINNGENID-indexed matrix: {nrow(finngenid_unrelated)} samples")
    }
  }

  # Check if sample_id matrix is empty (can happen if filtering used FINNGENIDs)
  if(nrow(phenotype_unrelated) == 0 && !is.null(finngenid_unrelated) && nrow(finngenid_unrelated) > 0) {
    log_warn("Sample ID matrix is empty, using FINNGENID matrix for rank normalisation")
    phenotype_unrelated <- finngenid_unrelated
    use_finngenid_as_primary <- TRUE
  } else {
    use_finngenid_as_primary <- FALSE
  }

  if(nrow(phenotype_unrelated) == 0) {
    log_error("Both phenotype matrices are empty. Cannot proceed with rank normalisation.")
    stop("No samples available for rank normalisation")
  }

  # Apply rank normalisation
  phenotype_rint <- rank_normalize_matrix(phenotype_unrelated, by = "column")

  if(!is.null(finngenid_unrelated) && !use_finngenid_as_primary) {
    finngenid_rint <- rank_normalize_matrix(finngenid_unrelated, by = "column")
  } else if(use_finngenid_as_primary) {
    # If we used FINNGENID matrix as primary, set finngenid_rint to the same
    finngenid_rint <- phenotype_rint
  } else {
    finngenid_rint <- NULL
  }

  # Compare distributions
  comparison_stats <- compare_distributions(phenotype_unrelated, phenotype_rint)

  # Create plots
  distribution_plots <- create_distribution_plots(phenotype_unrelated, phenotype_rint)

  # Format for PLINK
  # Use FINNGENID matrix for FID/IID columns (required for PLINK)
  if(!is.null(finngenid_rint)) {
    # Use FINNGENID matrix row names (FINNGENIDs) for FID and IID
    plink_format <- data.table(
      FID = rownames(finngenid_rint),
      IID = rownames(finngenid_rint)
    )
    # Use FINNGENID matrix for phenotype values
    plink_format <- cbind(plink_format, as.data.table(finngenid_rint))
  } else {
    # Fallback: use phenotype_rint if FINNGENID matrix not available
    # This should not happen in normal workflow
    log_warn("FINNGENID matrix not available, using phenotype matrix row names for PLINK format")
    plink_format <- data.table(
      FID = rownames(phenotype_rint),
      IID = rownames(phenotype_rint)
    )
    plink_format <- cbind(plink_format, as.data.table(phenotype_rint))
  }

  # Replace NA with -9 for PLINK
  phenotype_cols <- if(!is.null(finngenid_rint)) colnames(finngenid_rint) else colnames(phenotype_rint)
  for(col in phenotype_cols) {
    plink_format[is.na(get(col)), (col) := -9]
  }

  # Create protein list for PLINK
  protein_list <- data.table(
    protein = colnames(phenotype_rint),
    index = seq_len(ncol(phenotype_rint))
  )

  # Save outputs using batch-aware paths
  log_info("Saving rank normalised phenotypes")

  # Save matrices
  phenotype_rint_path <- get_output_path(step_num, "phenotype_matrix_rint", batch_id, "phenotypes", config = config)
  ensure_output_dir(phenotype_rint_path)
  saveRDS(phenotype_rint, phenotype_rint_path)

  if(!is.null(finngenid_rint)) {
    finngenid_rint_path <- get_output_path(step_num, "phenotype_matrix_finngenid_rint", batch_id, "phenotypes", config = config)
    ensure_output_dir(finngenid_rint_path)
    saveRDS(finngenid_rint, finngenid_rint_path)
  }

  # Save PLINK format
  plink_filename <- if (aggregate_output) {
    "batch2_Olink5k_unrel_rint.pheno"
  } else {
    "13_Olink5k_batch2_unrel_rint.pheno"
  }
  plink_path <- get_output_path(step_num, plink_filename, batch_id, "phenotypes", "pheno", config = config)
  ensure_output_dir(plink_path)
  fwrite(plink_format, plink_path, sep = "\t", na = "-9", quote = FALSE)

  # Save protein list
  protein_list_path <- get_output_path(step_num, "proteins_all", batch_id, "phenotypes", "txt", config = config)
  ensure_output_dir(protein_list_path)
  fwrite(protein_list, protein_list_path, sep = "\t", col.names = FALSE, quote = FALSE)

  # Save comparison statistics
  comparison_stats_path <- get_output_path(step_num, "rank_norm_comparison", batch_id, "phenotypes", "tsv", config = config)
  ensure_output_dir(comparison_stats_path)
  fwrite(comparison_stats, comparison_stats_path, sep = "\t")

  # Save plots
  if(length(distribution_plots$distribution_plots) >= 4) {
    distributions_path <- get_output_path(step_num, "rank_norm_distributions", batch_id, "phenotypes", "pdf", config = config)
    ensure_output_dir(distributions_path)
    pdf(distributions_path, width = 12, height = 10)
    gridExtra::grid.arrange(
      distribution_plots$distribution_plots[[1]],
      distribution_plots$distribution_plots[[2]],
      distribution_plots$distribution_plots[[3]],
      distribution_plots$distribution_plots[[4]],
      ncol = 2
    )
    dev.off()
  }

  qq_plot_path <- get_output_path(step_num, "rank_norm_qq_plot", batch_id, "phenotypes", "pdf", config = config)
  ensure_output_dir(qq_plot_path)
  ggsave(qq_plot_path, distribution_plots$qq_plot, width = 8, height = 8)

  # Create summary
  summary_stats <- data.table(
    metric = c(
      "Samples",
      "Proteins",
      "Total measurements",
      "Missing values",
      "Missing rate",
      "Mean skewness before",
      "Mean skewness after",
      "Mean kurtosis before",
      "Mean kurtosis after"
    ),
    value = c(
      nrow(phenotype_rint),
      ncol(phenotype_rint),
      length(phenotype_rint),
      sum(is.na(phenotype_rint)),
      round(sum(is.na(phenotype_rint)) / length(phenotype_rint) * 100, 2),
      round(mean(abs(comparison_stats$skew_orig), na.rm = TRUE), 3),
      round(mean(abs(comparison_stats$skew_norm), na.rm = TRUE), 3),
      round(mean(comparison_stats$kurt_orig, na.rm = TRUE), 3),
      round(mean(comparison_stats$kurt_norm, na.rm = TRUE), 3)
    )
  )

  summary_stats_path <- get_output_path(step_num, "rank_norm_summary", batch_id, "phenotypes", "tsv", config = config)
  ensure_output_dir(summary_stats_path)
  fwrite(summary_stats, summary_stats_path, sep = "\t")

  # Handle aggregation if enabled
  if (aggregate_output && multi_batch_mode && !is.null(finngenid_rint)) {
    log_info("Aggregation enabled: Attempting to rank-normalise batch 1 and create aggregate output")

    # Check if other batch processed data exists
    # Get other batch ID from config (not hardcoded)
    other_batch_id <- get_other_batch_id(batch_id, config)
    if (is.null(other_batch_id)) {
      log_warn("Could not determine other batch ID for aggregation. Skipping aggregation.")
      aggregate_output <- FALSE
    }

    # Check if batch correction is enabled
    adjust_for_batch <- tryCatch(
      isTRUE(config$parameters$covariate_adjustment$adjust_for_batch),
      error = function(e) FALSE
    )

    # Try to load batch-corrected kinship-filtered data first (if batch correction enabled)
    batch1_file <- NULL
    if (adjust_for_batch) {
      batch1_file_batch_corrected <- get_output_path("10", "phenotype_matrix_finngenid_unrelated_batch_corrected",
                                                      other_batch_id, "phenotypes", config = config)
      if (file.exists(batch1_file_batch_corrected)) {
        batch1_file <- batch1_file_batch_corrected
        log_info("Using BATCH-CORRECTED kinship-filtered data from batch 1 for rank normalization")
      }
    }

    # Fallback to non-batch-corrected data
    if (is.null(batch1_file)) {
      batch1_file <- get_output_path("10", "phenotype_matrix_finngenid_unrelated", other_batch_id, "phenotypes", config = config)
      adjust_for_batch <- FALSE  # Reset flag if batch-corrected data not available
    }

    if (file.exists(batch1_file)) {
      log_info("Found batch 1 kinship-filtered data: Creating aggregate rank-normalised output")
      if (adjust_for_batch) {
        log_info("  Data is BATCH-CORRECTED (adjust_for_batch = true)")
      }

      # Load batch 1 data
      batch1_unrelated <- readRDS(batch1_file)

      # Rank normalise batch 1
      batch1_rint <- rank_normalize_matrix(batch1_unrelated, by = "column")

      # Combine batch 1 and batch 2 (on FINNGENID, common proteins)
      common_proteins <- intersect(colnames(batch1_rint), colnames(finngenid_rint))
      if (length(common_proteins) > 100) {
        batch1_subset <- batch1_rint[, common_proteins, drop = FALSE]
        batch2_subset <- finngenid_rint[, common_proteins, drop = FALSE]

        # Check for duplicate FINNGENIDs (use batch 2 as reference)
        common_finngenids <- intersect(rownames(batch1_subset), rownames(batch2_subset))
        if (length(common_finngenids) > 0) {
          log_warn("Found {length(common_finngenids)} FINNGENIDs in both batches, using batch 2 data")
          batch1_subset <- batch1_subset[!rownames(batch1_subset) %in% common_finngenids, , drop = FALSE]
        }

        # Combine
        aggregate_rint <- rbind(batch1_subset, batch2_subset)
        log_info("Aggregate rank-normalised matrix: {nrow(aggregate_rint)} samples x {ncol(aggregate_rint)} proteins")

        # Format for PLINK
        aggregate_plink <- data.table(
          FID = rownames(aggregate_rint),
          IID = rownames(aggregate_rint)
        )
        aggregate_plink <- cbind(aggregate_plink, as.data.table(aggregate_rint))
        for(col in colnames(aggregate_rint)) {
          aggregate_plink[is.na(get(col)), (col) := -9]
        }

        # Save aggregate outputs
        log_info("Saving aggregate rank-normalised outputs with 'aggregate_' prefix")
        aggregate_dir <- file.path(config$output$base_dir %||% Sys.getenv("PIPELINE_OUTPUT_DIR", "output"),
                                   "phenotypes", "aggregate")
        dir.create(aggregate_dir, recursive = TRUE, showWarnings = FALSE)

        # Use different naming based on whether batch correction was applied
        if (adjust_for_batch) {
          aggregate_rint_path <- file.path(aggregate_dir, "aggregate_phenotype_matrix_rank_normalized_batch_corrected.rds")
          aggregate_plink_path <- file.path(aggregate_dir, "aggregate_phenotypes_rank_normalized_batch_corrected_plink.txt")
          aggregate_protein_index_path <- file.path(aggregate_dir, "aggregate_proteins_rank_normalized_batch_corrected_index.txt")
          log_info("Data is BATCH-CORRECTED - using '_batch_corrected' suffix")
        } else {
          aggregate_rint_path <- file.path(aggregate_dir, "aggregate_phenotype_matrix_rank_normalized.rds")
          aggregate_plink_path <- file.path(aggregate_dir, "aggregate_phenotypes_rank_normalized_plink.txt")
          aggregate_protein_index_path <- file.path(aggregate_dir, "aggregate_proteins_rank_normalized_index.txt")
        }

        saveRDS(aggregate_rint, aggregate_rint_path)
        fwrite(aggregate_plink, aggregate_plink_path, sep = "\t", na = "-9", quote = FALSE)
        fwrite(data.table(protein = colnames(aggregate_rint), index = seq_len(ncol(aggregate_rint))),
               aggregate_protein_index_path, sep = "\t", col.names = FALSE, quote = FALSE)

        log_info("Saved aggregate rank-normalised outputs:")
        log_info("  - Matrix: {aggregate_rint_path}")
        log_info("  - PLINK format: {aggregate_plink_path}")
        log_info("  - Protein index: {aggregate_protein_index_path}")
      } else {
        log_warn("Too few common proteins for aggregation")
      }
    } else {
      log_warn("Batch 1 kinship-filtered data not found. Aggregate output skipped.")
      log_warn("Expected file: {batch1_file}")
    }
  }

  # Print summary
  cat("\n=== RANK NORMALISATION SUMMARY ===\n")
  cat("Phenotype matrix:", nrow(phenotype_rint), "samples x", ncol(phenotype_rint), "proteins\n")
  cat("Missing values:", sum(is.na(phenotype_rint)),
      "(", round(sum(is.na(phenotype_rint)) / length(phenotype_rint) * 100, 2), "%)\n")
  cat("\nDistribution improvements:\n")
  cat("  Mean |skewness|: ", round(mean(abs(comparison_stats$skew_orig), na.rm = TRUE), 3),
      " -> ", round(mean(abs(comparison_stats$skew_norm), na.rm = TRUE), 3), "\n")
  cat("  Mean kurtosis: ", round(mean(comparison_stats$kurt_orig, na.rm = TRUE), 3),
      " -> ", round(mean(comparison_stats$kurt_norm, na.rm = TRUE), 3), "\n")
  cat("\nOutput files:\n")
  cat("  - PLINK phenotype file: 13_Olink5k_batch2_unrel_rint.pheno\n")
  cat("  - Protein list: 13_proteins_all.txt\n")
  cat("  - R matrices: 13_phenotype_matrix_rint.rds\n")
  cat("\nResults saved to: ../output/phenotypes/\n")

  log_info("Rank normalisation completed")

  return(list(
    phenotype_rint = phenotype_rint,
    finngenid_rint = finngenid_rint,
    comparison_stats = comparison_stats,
    summary = summary_stats
  ))
}

# Run if executed directly
if (!interactive()) {
  # Load moments package for skewness/kurtosis
  if(!requireNamespace("moments", quietly = TRUE)) {
    install.packages("moments")
  }
  library(moments)

  result <- main()
}






