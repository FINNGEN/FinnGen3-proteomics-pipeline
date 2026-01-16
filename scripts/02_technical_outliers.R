#!/usr/bin/env Rscript

#################################################
# Script: 02_technical_outliers.R
# Author: Reza Jabal, PhD (rjabal@broadinstitute.org)
# Description: Detect technical outliers based on plate and batch effects
#              Refactored version - loads from Step 01 (PCA-cleaned matrix)
# Date: December 2025
#################################################

suppressPackageStartupMessages({
  library(data.table)
  library(tidyverse)
  library(ggplot2)
  library(gridExtra)
  library(yaml)
  library(logger)
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
log_info("Starting technical outlier detection for batch: {batch_id}")

# Set theme for plots
theme_set(theme_bw())

# Function to detect plate outliers
detect_plate_outliers <- function(samples_data, npx_matrix) {
  log_info("Detecting plate-based outliers")

  # Calculate plate-level statistics
  plate_stats <- samples_data[, .(
    n_samples = length(unique(SampleID)),
    mean_npx = mean(NPX, na.rm = TRUE),
    sd_npx = sd(NPX, na.rm = TRUE),
    median_npx = median(NPX, na.rm = TRUE),
    missing_rate = sum(is.na(NPX)) / .N,
    qc_fail_rate = sum(SampleQC == "FAIL") / .N
  ), by = PlateID]

  # Identify outlier plates using 5-MAD (robust, aligned with 4-SD in z-score step)
  plate_center <- median(plate_stats$mean_npx, na.rm = TRUE)
  plate_mad <- mad(plate_stats$mean_npx, constant = 1.4826, na.rm = TRUE)
  plate_threshold <- 5 * plate_mad
  plate_stats[, outlier := abs(mean_npx - plate_center) > plate_threshold]

  outlier_plates <- plate_stats[outlier == TRUE]$PlateID

  log_info("Plate-level outlier detection:")
  log_info("  Total plates: {nrow(plate_stats)}")
  log_info("  Center (median): {sprintf('%.3f', plate_center)}")
  log_info("  MAD: {sprintf('%.3f', plate_mad)}")
  log_info("  Threshold: ±{sprintf('%.3f', plate_threshold)} (5×MAD ≈ 4×SD)")
  log_info("  Outlier plates: {length(outlier_plates)}")

  if(length(outlier_plates) > 0) {
    log_info("  Plate IDs: {paste(outlier_plates, collapse = ', ')}")
  }

  # Get samples from outlier plates
  samples_from_outlier_plates <- samples_data[PlateID %in% outlier_plates]$SampleID
  samples_from_outlier_plates <- unique(samples_from_outlier_plates)
  samples_from_outlier_plates <- intersect(samples_from_outlier_plates, rownames(npx_matrix))

  log_info("Samples from outlier plates: {length(samples_from_outlier_plates)}")

  return(list(
    plate_stats = plate_stats,
    outlier_plates = outlier_plates,
    samples_from_outlier_plates = samples_from_outlier_plates
  ))
}

# Function to detect batch effects
detect_batch_effects <- function(samples_data, metadata) {
  log_info("Detecting batch effects")

  # Merge samples with metadata
  sample_batch <- merge(
    samples_data[, .(SampleID, PlateID, NPX)],
    metadata[, .(SAMPLE_ID, APPROX_TIMESTAMP_COLLECTION, BIOBANK_PLASMA)],
    by.x = "SampleID",
    by.y = "SAMPLE_ID",
    all.x = TRUE
  )

  # Convert timestamp to date
  sample_batch[, collection_date := as.Date(APPROX_TIMESTAMP_COLLECTION)]
  sample_batch[, collection_month := format(collection_date, "%Y-%m")]

  # Calculate batch statistics by collection month
  batch_stats <- sample_batch[!is.na(collection_month), .(
    n_samples = .N,
    mean_npx = mean(NPX, na.rm = TRUE),
    sd_npx = sd(NPX, na.rm = TRUE)
  ), by = collection_month]

  # Identify outlier batches using 5-MAD (robust, aligned with 4-SD in z-score step)
  batch_center <- median(batch_stats$mean_npx, na.rm = TRUE)
  batch_mad <- mad(batch_stats$mean_npx, constant = 1.4826, na.rm = TRUE)
  batch_threshold <- 5 * batch_mad
  batch_stats[, outlier := abs(mean_npx - batch_center) > batch_threshold]

  outlier_batches <- batch_stats[outlier == TRUE]$collection_month

  log_info("Batch-level outlier detection (by collection month):")
  log_info("  Total batches: {nrow(batch_stats)}")
  log_info("  Center (median): {sprintf('%.3f', batch_center)}")
  log_info("  MAD: {sprintf('%.3f', batch_mad)}")
  log_info("  Threshold: ±{sprintf('%.3f', batch_threshold)} (5×MAD ≈ 4×SD)")
  log_info("  Outlier batches: {length(outlier_batches)}")

  if(length(outlier_batches) > 0) {
    log_info("  Batch IDs: {paste(outlier_batches, collapse = ', ')}")
  }

  # Get samples from outlier batches
  samples_from_outlier_batches <- sample_batch[collection_month %in% outlier_batches]$SampleID

  return(list(
    batch_stats = batch_stats,
    outlier_batches = outlier_batches,
    samples_from_outlier_batches = samples_from_outlier_batches
  ))
}

# Function to detect processing time outliers
detect_processing_outliers <- function(metadata) {
  log_info("Detecting processing time outliers")

  # Calculate processing time (time from collection to freezing)
  metadata_time <- copy(metadata)
  metadata_time[, collection_time := as.POSIXct(APPROX_TIMESTAMP_COLLECTION)]
  metadata_time[, freezing_time := as.POSIXct(APPROX_TIMESTAMP_FREEZING)]
  metadata_time[, processing_hours := as.numeric(difftime(freezing_time, collection_time, units = "hours"))]

  # Remove invalid processing times
  metadata_time <- metadata_time[!is.na(processing_hours) & processing_hours > 0 & processing_hours < 48]

  # Identify outliers based on processing time (5-MAD, aligned with 4-SD in z-score step)
  proc_center <- median(metadata_time$processing_hours, na.rm = TRUE)
  proc_mad <- mad(metadata_time$processing_hours, constant = 1.4826, na.rm = TRUE)
  proc_threshold <- 5 * proc_mad
  metadata_time[, processing_outlier := abs(processing_hours - proc_center) > proc_threshold]

  processing_outliers <- metadata_time[processing_outlier == TRUE]$SAMPLE_ID

  log_info("Processing time outlier detection:")
  log_info("  Total samples with valid processing time: {nrow(metadata_time)}")
  log_info("  Center (median): {sprintf('%.2f', proc_center)} hours")
  log_info("  MAD: {sprintf('%.2f', proc_mad)} hours")
  log_info("  Threshold: ±{sprintf('%.2f', proc_threshold)} hours (5×MAD ≈ 4×SD)")
  log_info("  Processing time outliers: {length(processing_outliers)} samples")

  return(list(
    processing_stats = metadata_time[, .(SAMPLE_ID, processing_hours, processing_outlier)],
    processing_outliers = processing_outliers
  ))
}

# Function to detect sample-level technical outliers
detect_sample_technical_outliers <- function(samples_data) {
  log_info("Detecting sample-level technical outliers")
  log_info("  Note: samples_data has been filtered to match npx_matrix from Step 01")

  # Calculate sample-level statistics
  sample_stats <- samples_data[, .(
    n_proteins = .N,
    mean_npx = mean(NPX, na.rm = TRUE),
    sd_npx = sd(NPX, na.rm = TRUE),
    missing_rate = sum(is.na(NPX)) / .N,
    qc_fail_rate = sum(SampleQC == "FAIL" | AssayQC == "FAIL") / .N
  ), by = SampleID]

  # Identify outliers based on multiple criteria (robust cutoffs, aligned with z-score step)
  center_mean <- median(sample_stats$mean_npx, na.rm = TRUE)
  mad_mean <- mad(sample_stats$mean_npx, constant = 1.4826, na.rm = TRUE)
  mean_threshold <- 5 * mad_mean
  sample_stats[, mean_outlier := abs(mean_npx - center_mean) > mean_threshold]

  # High-variance samples (robust one-sided threshold: 4×MAD, aligned with z-score 4×SD)
  center_sd <- median(sample_stats$sd_npx, na.rm = TRUE)
  mad_sd <- mad(sample_stats$sd_npx, constant = 1.4826, na.rm = TRUE)
  sd_threshold <- center_sd + 4 * mad_sd
  sample_stats[, sd_outlier := sd_npx > sd_threshold]

  # Missing and QC thresholds
  sample_stats[, missing_outlier := missing_rate > 0.05]
  sample_stats[, qc_outlier := qc_fail_rate > 0.3]

  # Log detailed statistics
  log_info("Sample-level outlier detection thresholds:")
  log_info("  Mean NPX: center={sprintf('%.3f', center_mean)}, MAD={sprintf('%.3f', mad_mean)}, threshold=±{sprintf('%.3f', mean_threshold)} (5×MAD ≈ 4×SD)")
  log_info("  SD NPX: center={sprintf('%.3f', center_sd)}, MAD={sprintf('%.3f', mad_sd)}, threshold={sprintf('%.3f', sd_threshold)} (median + 4×MAD)")
  log_info("  Missing rate: threshold=5%")
  log_info("  QC fail rate: threshold=30%")

  # Combined outlier flag
  sample_stats[, is_outlier := mean_outlier | sd_outlier | missing_outlier | qc_outlier]

  technical_outliers <- sample_stats[is_outlier == TRUE]$SampleID

  # Log breakdown by criterion
  log_info("Sample-level outlier breakdown:")
  log_info("  Mean outliers: {sum(sample_stats$mean_outlier)} samples")
  log_info("  SD outliers: {sum(sample_stats$sd_outlier)} samples")
  log_info("  Missing outliers: {sum(sample_stats$missing_outlier)} samples")
  log_info("  QC outliers: {sum(sample_stats$qc_outlier)} samples")
  log_info("  Total sample-level outliers (OR logic): {length(technical_outliers)} samples")

  return(list(
    sample_stats = sample_stats,
    technical_outliers = technical_outliers
  ))
}

# Function to create technical outlier plots
create_technical_plots <- function(plate_result, batch_result, sample_result) {
  log_info("Creating technical outlier visualizations")

  # Plate statistics plot
  p1 <- ggplot(plate_result$plate_stats, aes(x = PlateID, y = mean_npx)) +
    geom_bar(stat = "identity", aes(fill = outlier)) +
    scale_fill_manual(values = c("FALSE" = "#3A5F8A", "TRUE" = "#FF6B6B")) +
    labs(title = "Mean NPX by Plate",
         x = "Plate ID", y = "Mean NPX",
         fill = "Outlier") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  # Batch statistics plot
  if(nrow(batch_result$batch_stats) > 0) {
    p2 <- ggplot(batch_result$batch_stats, aes(x = collection_month, y = mean_npx)) +
      geom_point(aes(color = outlier), size = 3) +
      geom_line(group = 1, alpha = 0.5) +
      scale_color_manual(values = c("FALSE" = "#3A5F8A", "TRUE" = "#FF6B6B")) +
      labs(title = "Mean NPX by Collection Month",
           x = "Collection Month", y = "Mean NPX",
           color = "Outlier") +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  } else {
    p2 <- NULL
  }

  # Sample technical outliers
  p3 <- ggplot(sample_result$sample_stats, aes(x = mean_npx, y = sd_npx)) +
    geom_point(aes(color = is_outlier), alpha = 0.6) +
    scale_color_manual(values = c("FALSE" = "#3A5F8A", "TRUE" = "#FF6B6B")) +
    labs(title = "Sample Technical Characteristics",
         x = "Mean NPX", y = "SD NPX",
         color = "Outlier") +
    theme_bw()

  # Missing rate distribution with categorical bar chart
  stats_dt <- sample_result$sample_stats
  n_zero <- sum(stats_dt$missing_rate == 0)
  n_low <- sum(stats_dt$missing_rate > 0 & stats_dt$missing_rate <= 0.05)
  n_high <- sum(stats_dt$missing_rate > 0.05)

  cat_data <- data.table(
    Category = factor(c("0% Missing", "0-5% Missing", ">5% Missing"),
                     levels = c("0% Missing", "0-5% Missing", ">5% Missing")),
    Count = c(n_zero, n_low, n_high),
    Percentage = c(100*n_zero/nrow(stats_dt), 100*n_low/nrow(stats_dt), 100*n_high/nrow(stats_dt))
  )

  p4_cat <- ggplot(cat_data, aes(x = Category, y = Count, fill = Category)) +
    geom_bar(stat = "identity", alpha = 0.8) +
    geom_text(aes(label = sprintf("%d\n(%.1f%%)", Count, Percentage)),
              vjust = -0.5, size = 3.5, fontface = "bold") +
    scale_fill_manual(values = c("0% Missing" = "#2E7D32",
                                  "0-5% Missing" = "#FFA726",
                                  ">5% Missing" = "#E53935")) +
    labs(title = "Missing Rate Categories",
         subtitle = sprintf("Total: %d samples | Excellent quality: %.1f%% have zero missing",
                           nrow(stats_dt), 100*n_zero/nrow(stats_dt)),
         x = "", y = "Sample Count") +
    theme_bw() +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 0, hjust = 0.5),
          plot.subtitle = element_text(size = 9, color = "gray30")) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15)))

  non_zero_missing <- stats_dt[missing_rate > 0]
  if(nrow(non_zero_missing) > 0) {
    p4_hist <- ggplot(non_zero_missing, aes(x = missing_rate)) +
      geom_histogram(bins = 20, fill = "#5C9EAD", alpha = 0.8, boundary = 0) +
      geom_vline(xintercept = 0.05, color = "red", linetype = "dashed", linewidth = 1) +
      annotate("text", x = 0.05, y = Inf, label = "5% threshold",
               vjust = 1.5, hjust = -0.1, color = "red", size = 3) +
      labs(title = "Distribution of Non-Zero Missing Rates",
           subtitle = sprintf("Only %d samples (%.2f%%) have any missing data",
                             nrow(non_zero_missing), 100*nrow(non_zero_missing)/nrow(stats_dt)),
         x = "Missing Rate", y = "Count") +
      theme_bw() +
      theme(plot.subtitle = element_text(size = 9, color = "gray30"))

    p4 <- gridExtra::grid.arrange(p4_cat, p4_hist, ncol = 2, widths = c(1, 1.2))
  } else {
    p4 <- p4_cat
  }

  return(list(
    plate_plot = p1,
    batch_plot = p2,
    sample_plot = p3,
    missing_plot = p4
  ))
}

# Function to combine technical outliers
combine_technical_outliers <- function(plate_outliers, batch_outliers,
                                      processing_outliers, sample_outliers) {
  log_info("Combining technical outliers")

  all_technical_outliers <- unique(c(
    plate_outliers,
    batch_outliers,
    processing_outliers,
    sample_outliers
  ))

  outlier_summary <- data.table(
    SampleID = all_technical_outliers,
    is_plate_outlier = all_technical_outliers %in% plate_outliers,
    is_batch_outlier = all_technical_outliers %in% batch_outliers,
    is_processing_outlier = all_technical_outliers %in% processing_outliers,
    is_sample_outlier = all_technical_outliers %in% sample_outliers
  )

  outlier_summary[, outlier_score := is_plate_outlier + is_batch_outlier +
                    is_processing_outlier + is_sample_outlier]

  log_info("Total technical outliers: {length(all_technical_outliers)}")

  return(list(
    all_outliers = all_technical_outliers,
    outlier_summary = outlier_summary
  ))
}

# Main execution
main <- function() {

  # Load data from previous step (Step 00 - analysis-ready base matrix)
  # CRITICAL: Use base matrix for parallel flagging (matches original pipeline design)
  # All outlier detection steps (01, 02, 03) now use the same base matrix for parallel flagging
  log_info("Loading data from previous step (Step 00)")
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

  # Load samples_data if available (for plate/batch analysis)
  # Note: samples_data may not be available in refactored pipeline
  # If not available, we'll extract from metadata and npx_matrix
  samples_data_path <- get_output_path(prev_step00_num, "samples_data_raw", batch_id, "qc", config = config)
  if (file.exists(samples_data_path)) {
    samples_data <- readRDS(samples_data_path)
    # Filter to match samples in matrix
    samples_in_matrix <- rownames(npx_matrix)
    samples_data <- samples_data[SampleID %in% samples_in_matrix]
    log_info("Loaded samples_data: {length(unique(samples_data$SampleID))} samples")
  } else {
    log_warn("samples_data_raw not found. Creating minimal samples_data from metadata and matrix.")
    # Create minimal samples_data from metadata
    samples_data <- metadata[, .(SampleID = SAMPLE_ID, PlateID, NPX = NA_real_)]
    # This is a fallback - plate detection may be limited
  }

  # Detect plate outliers
  plate_result <- detect_plate_outliers(samples_data, npx_matrix)

  # Detect batch effects
  batch_result <- detect_batch_effects(samples_data, metadata)

  # Detect processing time outliers
  processing_result <- detect_processing_outliers(metadata)

  # Detect sample-level technical outliers
  sample_result <- detect_sample_technical_outliers(samples_data)

  # Combine all technical outliers
  combined_result <- combine_technical_outliers(
    plate_result$samples_from_outlier_plates,
    batch_result$samples_from_outlier_batches,
    processing_result$processing_outliers,
    sample_result$technical_outliers
  )

  # Cross-reference technical outliers with sex mismatches (if available)
  sex_mismatch_file <- get_output_path("04", "sex_mismatches", batch_id, "outliers", "tsv", config = config)
  if (file.exists(sex_mismatch_file)) {
    log_info("Cross-referencing technical outliers with sex mismatches")
    sex_mm <- try(fread(sex_mismatch_file), silent = TRUE)
    if (!inherits(sex_mm, "try-error")) {
      id_col <- if ("SampleID" %in% names(sex_mm)) "SampleID" else if ("SAMPLE_ID" %in% names(sex_mm)) "SAMPLE_ID" else NULL
      if (!is.null(id_col)) {
        crossref <- merge(
          combined_result$outlier_summary,
          sex_mm[, .(SampleID = get(id_col), sex_mismatch = TRUE)],
          by = "SampleID",
          all.x = TRUE
        )
        crossref[is.na(sex_mismatch), sex_mismatch := FALSE]

        crossref <- add_finngenid_column(crossref, batch_id = batch_id, config = config,
                                         sample_id_col = "SampleID", preserve_original = TRUE)

        crossref_path <- get_output_path(step_num, "technical_sex_mismatch_crossref", batch_id, "outliers", "tsv", config = config)
        ensure_output_dir(crossref_path)
        fwrite(crossref, crossref_path, sep = "\t")
      }
    }
  }

  # Create plots
  technical_plots <- create_technical_plots(plate_result, batch_result, sample_result)

  # Remove technical outliers from matrix
  npx_clean <- npx_matrix[!rownames(npx_matrix) %in% combined_result$all_outliers, ]

  log_info("Matrix after technical outlier removal: {nrow(npx_clean)} x {ncol(npx_clean)}")

  # Save outputs
  log_info("Saving technical outlier detection results")

  plate_result_path <- get_output_path(step_num, "plate_outliers", batch_id, "outliers", config = config)
  batch_result_path <- get_output_path(step_num, "batch_outliers", batch_id, "outliers", config = config)
  processing_result_path <- get_output_path(step_num, "processing_outliers", batch_id, "outliers", config = config)
  sample_result_path <- get_output_path(step_num, "sample_technical_outliers", batch_id, "outliers", config = config)
  combined_result_path <- get_output_path(step_num, "technical_outliers_combined", batch_id, "outliers", config = config)
  npx_clean_path <- get_output_path(step_num, "npx_matrix_technical_cleaned", batch_id, "outliers", config = config)

  ensure_output_dir(plate_result_path)
  ensure_output_dir(batch_result_path)
  ensure_output_dir(processing_result_path)
  ensure_output_dir(sample_result_path)
  ensure_output_dir(combined_result_path)
  ensure_output_dir(npx_clean_path)

  saveRDS(plate_result, plate_result_path)
  saveRDS(batch_result, batch_result_path)
  saveRDS(processing_result, processing_result_path)

  # Add FINNGENID to output tables
  log_info("Adding FINNGENID mapping to output tables...")
  outlier_summary_with_fgid <- add_finngenid_column(combined_result$outlier_summary, batch_id = batch_id, config = config)
  outlier_summary_path <- get_output_path(step_num, "technical_outlier_summary", batch_id, "outliers", "tsv", config = config)
  ensure_output_dir(outlier_summary_path)
  fwrite(outlier_summary_with_fgid, outlier_summary_path, sep = "\t")

  plate_stats_path <- get_output_path(step_num, "plate_statistics", batch_id, "outliers", "tsv", config = config)
  ensure_output_dir(plate_stats_path)
  fwrite(plate_result$plate_stats, plate_stats_path, sep = "\t")

  sample_stats_with_fgid <- add_finngenid_column(sample_result$sample_stats, batch_id = batch_id, config = config)
  sample_stats_path <- get_output_path(step_num, "sample_technical_stats", batch_id, "outliers", "tsv", config = config)
  ensure_output_dir(sample_stats_path)
  fwrite(sample_stats_with_fgid, sample_stats_path, sep = "\t")

  sample_result$sample_stats <- sample_stats_with_fgid
  saveRDS(sample_result, sample_result_path)

  combined_result$outlier_summary <- outlier_summary_with_fgid
  saveRDS(combined_result, combined_result_path)
  saveRDS(npx_clean, npx_clean_path)

  # Save plots
  plate_plot_path <- get_output_path(step_num, "technical_plate_stats", batch_id, "outliers", "pdf", config = config)
  batch_plot_path <- get_output_path(step_num, "technical_batch_stats", batch_id, "outliers", "pdf", config = config)
  sample_plot_path <- get_output_path(step_num, "technical_sample_stats", batch_id, "outliers", "pdf", config = config)
  missing_plot_path <- get_output_path(step_num, "technical_missing_dist", batch_id, "outliers", "pdf", config = config)

  ensure_output_dir(plate_plot_path)
  ensure_output_dir(sample_plot_path)
  ensure_output_dir(missing_plot_path)

  ggsave(plate_plot_path, technical_plots$plate_plot, width = 12, height = 6)
  if(!is.null(technical_plots$batch_plot)) {
    ensure_output_dir(batch_plot_path)
    ggsave(batch_plot_path, technical_plots$batch_plot, width = 10, height = 6)
  }
  ggsave(sample_plot_path, technical_plots$sample_plot, width = 8, height = 8)

  pdf(missing_plot_path, width = 14, height = 6)
  grid::grid.draw(technical_plots$missing_plot)
  dev.off()

  # Print summary
  cat("\n=== TECHNICAL OUTLIER DETECTION SUMMARY ===\n")
  cat("Input Matrix:\n")
  cat("  Samples (after PCA outlier removal): ", nrow(npx_matrix), "\n")
  cat("  Proteins: ", ncol(npx_matrix), "\n\n")

  cat("Detection Results by Method:\n")
  cat("  Plate outliers:              ", length(plate_result$outlier_plates), " plates\n")
  cat("  Batch outliers:              ", length(batch_result$outlier_batches), " batches\n")
  cat("  Processing time outliers:    ", length(processing_result$processing_outliers), " samples\n")
  cat("  Sample-level outliers:       ", length(sample_result$technical_outliers), " samples\n\n")

  cat("Combined Results (OR logic):\n")
  cat("  Total technical outliers:       ", length(combined_result$all_outliers),
      sprintf(" (%.1f%% of samples)\n", 100*length(combined_result$all_outliers)/nrow(npx_matrix)))
  cat("  Samples retained:               ", nrow(npx_clean),
      sprintf(" (%.1f%% of samples)\n\n", 100*nrow(npx_clean)/nrow(npx_matrix)))

  cat("Output Matrix:\n")
  cat("  Samples: ", nrow(npx_clean), " × Proteins: ", ncol(npx_clean), "\n\n")

  log_info("Technical outlier detection completed")

  return(list(
    plate_outliers = plate_result,
    batch_outliers = batch_result,
    technical_outliers = combined_result,
    clean_matrix = npx_clean
  ))
}

# Run if executed directly
if (!interactive()) {
  result <- main()
}




