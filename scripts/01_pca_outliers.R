#!/usr/bin/env Rscript

#################################################
# Script: 01_pca_outliers.R
# Author: Reza Jabal, PhD (rjabal@broadinstitute.org)
# Description: PCA-based outlier detection for proteomics data
#              Refactored version - loads from Step 00 (pre-filtered NPX matrix)
# Date: December 2025
#################################################

suppressPackageStartupMessages({
  library(data.table)
  library(tidyverse)
  library(ggplot2)
  library(ggrepel)
  library(yaml)
  library(logger)
  library(scales)
  library(wesanderson)
})
utils::globalVariables(c(
  "Feature","Type","collection_year","APPROX_TIMESTAMP_COLLECTION","PC1","PC2","PC3","PC4",
  "SampleID","is_outlier","outlier_score","SAMPLE_ID","FINNGENID","COHORT_FINNGENID","BIOBANK_PLASMA",
  "sample_source","point_shape","annotation_label","flag_s","flag","flag_num","flag_bool","Group",
  "DISEASE_GROUP","groups","PC1_PC2","PC3_PC4","Median","Any","PlateID","..meta_avail","WellID",
  "gPC1","gPC2","gPC3","gPC4","gPC2_q","gPC4_q","IID","..cols"
))

# Source path utilities
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
log_info("Starting PCA-based outlier detection for batch: {batch_id}")

# Set theme for plots
theme_set(theme_bw())

# Helper: HUSL-like palette (uses hsluv if available, otherwise hue palette)
get_husl_palette <- function(n) {
  if (requireNamespace("hsluv", quietly = TRUE)) {
    # Evenly spaced hues in HSLuv space with fixed saturation/lightness
    hues <- seq(0, 360, length.out = n + 1)[1:n]
    cols <- vapply(hues, function(h) {
      # hsluv returns hex given h,s,l
      hsluv::hsluv(h = h, s = 100, l = 65)
    }, character(1))
    return(cols)
  }
  # Fallback to evenly spaced hues (close to sns.husl)
  return(scales::hue_pal()(n))
}

# Helper: continuous Wes Anderson "Zissou1" palette scaled to n groups
get_wes_cont_palette <- function(n) {
  base10 <- wes_palette("Zissou1", 10, type = "continuous")
  return(colorRampPalette(base10)(n))
}

# Convert protein loadings to long/wide tables
make_protein_loadings_dt <- function(pca_result) {
  load_mat <- pca_result$loadings
  if (is.null(load_mat)) return(NULL)
  dt <- as.data.table(load_mat)
  dt[, Feature := rownames(load_mat)]
  setcolorder(dt, c("Feature", setdiff(colnames(dt), "Feature")))
  dt[, Type := "protein"]
  setcolorder(dt, c("Type", "Feature", setdiff(colnames(dt), c("Type", "Feature"))))
  return(dt)
}

# Compute technical covariate associations with PCs (R^2 for categorical, Pearson r for numeric)
compute_technical_associations <- function(pca_scores, metadata) {
  if (is.null(metadata)) return(NULL)
  pc_dt <- as.data.table(pca_scores)
  pc_dt$SampleID <- rownames(pca_scores)

  # Select and derive potential technical covariates
  meta <- as.data.table(metadata)
  # Derive collection_year if timestamp present
  if ("APPROX_TIMESTAMP_COLLECTION" %in% colnames(meta)) {
    suppressWarnings({
      meta[, collection_year := as.numeric(format(as.POSIXct(APPROX_TIMESTAMP_COLLECTION), "%Y"))]
    })
  }

  tech_candidates <- intersect(c("PlateID", "BIOBANK_PLASMA", "COHORT_FINNGENID", "collection_year"), colnames(meta))
  if (length(tech_candidates) == 0) return(NULL)

  dt <- merge(pc_dt, meta[, c("SAMPLE_ID", tech_candidates), with = FALSE], by.x = "SampleID", by.y = "SAMPLE_ID", all.x = TRUE)

  pcs <- colnames(pc_dt)[startsWith(colnames(pc_dt), "PC")]
  res_list <- list()

  for (var in tech_candidates) {
    v <- dt[[var]]
    # Decide numeric vs categorical
    is_num <- suppressWarnings(!any(is.na(as.numeric(as.character(v)))))
    # But treat clearly categorical strings as categorical
    if (is.character(v) || is.factor(v)) is_num <- FALSE

    row_vals <- c(Type = "technical", Variable = var)
    for (pc in pcs) {
      pc_vals <- dt[[pc]]
      val <- NA_real_
      if (is_num) {
        suppressWarnings({ val <- suppressWarnings(stats::cor(as.numeric(v), pc_vals, use = "pairwise.complete.obs")) })
      } else {
        # R^2 from linear model PC ~ factor(var)
        suppressWarnings({
          try({ val <- summary(lm(pc_vals ~ as.factor(v)))$r.squared }, silent = TRUE)
        })
      }
      row_vals <- c(row_vals, setNames(val, pc))
    }
    res_list[[length(res_list) + 1]] <- as.data.table(as.list(row_vals))
  }

  assoc_dt <- rbindlist(res_list, fill = TRUE)
  # Ensure numeric on PC columns
  for (pc in pcs) assoc_dt[[pc]] <- as.numeric(assoc_dt[[pc]])
  return(assoc_dt)
}

# Function to perform PCA
perform_pca <- function(npx_matrix, n_components = 10) {
  log_info("Performing PCA with {n_components} components")

  # Handle missing values
  npx_complete <- npx_matrix

  # Impute missing values with column median
  for(i in seq_len(ncol(npx_complete))) {
    col_median <- median(npx_complete[, i], na.rm = TRUE)
    if(is.na(col_median)) {
      # If all values are NA, impute with 0
      npx_complete[, i] <- 0
    } else {
      npx_complete[is.na(npx_complete[, i]), i] <- col_median
    }
  }

  # Remove constant columns (zero variance)
  col_vars <- apply(npx_complete, 2, var, na.rm = TRUE)
  constant_cols <- which(col_vars == 0 | is.na(col_vars))

  if(length(constant_cols) > 0) {
    log_info("Removing {length(constant_cols)} constant/zero-variance columns")
    npx_complete <- npx_complete[, -constant_cols]
  }

  # Center and scale
  npx_scaled <- scale(npx_complete)

  # Perform PCA
  pca_result <- prcomp(npx_scaled, center = FALSE, scale. = FALSE)

  # Extract components
  n_actual_components <- min(n_components, ncol(pca_result$x))
  pca_scores <- pca_result$x[, seq_len(n_actual_components)]

  # Calculate variance explained
  # n_actual_components already computed
  var_explained <- summary(pca_result)$importance[2, 1:n_actual_components]
  cum_var_explained <- summary(pca_result)$importance[3, 1:n_actual_components]

  log_info("Variance explained by first {n_actual_components} PCs: {round(cum_var_explained[n_actual_components] * 100, 2)}%")

  return(list(
    scores = pca_scores,
    loadings = pca_result$rotation[, 1:n_actual_components],
    var_explained = var_explained,
    cum_var_explained = cum_var_explained,
    pca_object = pca_result
  ))
}

# Function to detect outliers using PCA
detect_pca_outliers <- function(pca_scores, sd_threshold = 3) {
  log_info("Detecting outliers with SD threshold: {sd_threshold}")

  outliers <- list()

  # For each PC, identify outliers
  for(i in seq_len(ncol(pca_scores))) {
    pc_name <- paste0("PC", i)
    pc_values <- pca_scores[, i]

    # Calculate mean and SD
    pc_mean <- mean(pc_values)
    pc_sd <- sd(pc_values)

    # Identify outliers
    outlier_idx <- abs(pc_values - pc_mean) > sd_threshold * pc_sd

    outliers[[pc_name]] <- rownames(pca_scores)[outlier_idx]

    if(sum(outlier_idx) > 0) {
      log_info("{pc_name}: {sum(outlier_idx)} outliers detected")
    }
  }

  # Combined outliers (samples that are outliers in any PC)
  all_outliers <- unique(unlist(outliers))

  # Calculate outlier score (number of PCs where sample is outlier)
  outlier_score <- sapply(rownames(pca_scores), function(x) {
    sum(sapply(outliers, function(y) x %in% y))
  })

  # Create outlier summary
  outlier_summary <- data.table(
    SampleID = rownames(pca_scores),
    outlier_score = outlier_score,
    is_outlier = outlier_score > 0,
    PC1 = pca_scores[, 1],
    PC2 = pca_scores[, 2]
  )

  log_info("Total outliers detected: {length(all_outliers)}")

  return(list(
    outliers_by_pc = outliers,
    all_outliers = all_outliers,
    outlier_summary = outlier_summary
  ))
}

# Function to detect outliers using sample-level IQR (match original batch1)
detect_iqr_outliers <- function(npx_matrix, sd_threshold = 5) {
  log_info("Detecting outliers using sample-level IQR with SD threshold: {sd_threshold}")

  # IQR across proteins per sample (row-wise)
  sample_iqr <- apply(npx_matrix, 1, IQR, na.rm = TRUE)

  mean_iqr <- mean(sample_iqr, na.rm = TRUE)
  sd_iqr <- sd(sample_iqr, na.rm = TRUE)

  lower_bound <- mean_iqr - sd_threshold * sd_iqr
  upper_bound <- mean_iqr + sd_threshold * sd_iqr

  outlier_idx <- which(sample_iqr < lower_bound | sample_iqr > upper_bound)
  outliers <- rownames(npx_matrix)[outlier_idx]

  log_info("Sample IQR outliers detected: {length(outliers)}")
  log_info("IQR range: {round(min(sample_iqr, na.rm=TRUE), 3)} - {round(max(sample_iqr, na.rm=TRUE), 3)}")
  log_info("Outlier boundaries: {round(lower_bound, 3)} - {round(upper_bound, 3)}")

  return(list(
    outliers = outliers,
    sample_iqr = sample_iqr,
    mean_iqr = mean_iqr,
    sd_iqr = sd_iqr,
    lower_bound = lower_bound,
    upper_bound = upper_bound
  ))
}

# Original-batch1 compatible outlier detection (PC1/PC2 5SD -> median 5SD -> IQR 5SD)
detect_outliers_original <- function(npx_matrix, pca_object) {
  log_info("Detecting outliers using original batch1 procedure (PC1/PC2 -> median -> IQR)")

  # Scale lambda as in batch1 implementation
  n_samples <- nrow(npx_matrix)
  scale_lambda <- pca_object$sdev * sqrt(n_samples)

  # Adjust PC1 and PC2 as in the original code
  pcx <- pca_object$x
  pc1_adj <- pcx[, "PC1"] / scale_lambda[1]
  pc2_adj <- pcx[, "PC2"] / scale_lambda[2]
  pc3_adj <- if (ncol(pcx) >= 3) pcx[, "PC3"] / scale_lambda[3] else rep(0, nrow(pcx))
  pc4_adj <- if (ncol(pcx) >= 4) pcx[, "PC4"] / scale_lambda[4] else rep(0, nrow(pcx))

  dt_pc <- data.table(SampleID = rownames(pcx), PC1 = pc1_adj, PC2 = pc2_adj, PC3 = pc3_adj, PC4 = pc4_adj)

  # Step 1: PC1 & PC2 within mean ± 5 SD
  mean1 <- mean(dt_pc$PC1); sd1 <- sd(dt_pc$PC1)
  mean2 <- mean(dt_pc$PC2); sd2 <- sd(dt_pc$PC2)
  upper1 <- mean1 + 5 * sd1; lower1 <- mean1 - 5 * sd1
  upper2 <- mean2 + 5 * sd2; lower2 <- mean2 - 5 * sd2
  keep_idx1 <- dt_pc[PC1 >= lower1 & PC1 <= upper1 & PC2 >= lower2 & PC2 <= upper2, which = TRUE]
  npx_val1 <- npx_matrix[keep_idx1, , drop = FALSE]
  dt_batch1 <- dt_pc[keep_idx1]
  log_info("After PC1/PC2 filtering (±5 SD): {nrow(dt_batch1)} samples remain")

  # Step 1b: PC3 & PC4 within mean ± 5 SD (same logic extended)
  if (ncol(pca_object$x) >= 4) {
    mean3 <- mean(dt_batch1$PC3); sd3 <- sd(dt_batch1$PC3)
    mean4 <- mean(dt_batch1$PC4); sd4 <- sd(dt_batch1$PC4)
    upper3 <- mean3 + 5 * sd3; lower3 <- mean3 - 5 * sd3
    upper4 <- mean4 + 5 * sd4; lower4 <- mean4 - 5 * sd4
    keep_idx_pc34 <- dt_batch1[PC3 >= lower3 & PC3 <= upper3 & PC4 >= lower4 & PC4 <= upper4, which = TRUE]
    npx_val1b <- npx_val1[keep_idx_pc34, , drop = FALSE]
    dt_batch1b <- dt_batch1[keep_idx_pc34]
    log_info("After PC3/PC4 filtering (±5 SD): {nrow(dt_batch1b)} samples remain")
  } else {
    dt_batch1b <- dt_batch1
    npx_val1b <- npx_val1
  }

  # Step 2: Row median within mean ± 5 SD
  row_median <- apply(npx_val1b, 1, median, na.rm = TRUE)
  mean_med <- mean(row_median, na.rm = TRUE); sd_med <- sd(row_median, na.rm = TRUE)
  upper_med <- mean_med + 5 * sd_med; lower_med <- mean_med - 5 * sd_med
  keep_idx2 <- which(row_median >= lower_med & row_median <= upper_med)
  npx_val2 <- npx_val1b[keep_idx2, , drop = FALSE]
  dt_batch2 <- dt_batch1b[keep_idx2]
  log_info("After median filtering (±5 SD): {nrow(dt_batch2)} samples remain")

  # Step 3: Row IQR within mean ± 5 SD
  row_iqr <- apply(npx_val2, 1, IQR, na.rm = TRUE)
  mean_iqr <- mean(row_iqr, na.rm = TRUE); sd_iqr <- sd(row_iqr, na.rm = TRUE)
  upper_iqr <- mean_iqr + 5 * sd_iqr; lower_iqr <- mean_iqr - 5 * sd_iqr
  keep_idx3 <- which(row_iqr >= lower_iqr & row_iqr <= upper_iqr)
  npx_val3 <- npx_val2[keep_idx3, , drop = FALSE]
  dt_batch3 <- dt_batch2[keep_idx3]
  log_info("After IQR filtering (±5 SD): {nrow(dt_batch3)} samples remain")

  final_kept <- dt_batch3$SampleID
  all_outliers <- setdiff(rownames(npx_matrix), final_kept)

  # Build outlier summary with source flags and additional metrics
  flag_pc12 <- setdiff(dt_pc$SampleID, dt_batch1$SampleID)
  flag_pc34 <- setdiff(dt_batch1$SampleID, dt_batch1b$SampleID)
  flag_med  <- setdiff(dt_batch1$SampleID, dt_batch2$SampleID)
  flag_iqr  <- setdiff(dt_batch2$SampleID, dt_batch3$SampleID)

  # Map median and IQR where available (computed on filtered subsets); others NA
  row_median_map <- setNames(row_median, dt_batch1b$SampleID)
  row_iqr_map    <- setNames(row_iqr, dt_batch2$SampleID)

  outlier_summary <- data.table(
    SampleID    = dt_pc$SampleID,
    PC1         = dt_pc$PC1,
    PC2         = dt_pc$PC2,
    PC3         = dt_pc$PC3,
    PC4         = dt_pc$PC4,
    row_median  = unname(row_median_map[dt_pc$SampleID]),
    row_iqr     = unname(row_iqr_map[dt_pc$SampleID])
  )

  outlier_summary[, outlier_PC1_PC2 := as.integer(SampleID %in% flag_pc12)]
  outlier_summary[, outlier_PC3_PC4 := as.integer(SampleID %in% flag_pc34)]
  outlier_summary[, outlier_Median  := as.integer(SampleID %in% flag_med)]
  outlier_summary[, outlier_IQR     := as.integer(SampleID %in% flag_iqr)]
  outlier_summary[, outlier_score   := outlier_PC1_PC2 + outlier_PC3_PC4 + outlier_Median + outlier_IQR]
  outlier_summary[, is_outlier      := outlier_score > 0]

  return(list(
    outliers_pc = setdiff(dt_pc$SampleID, dt_batch1$SampleID),
    outliers_pc34 = setdiff(dt_batch1$SampleID, dt_batch1b$SampleID),
    outliers_median = setdiff(dt_batch1$SampleID, dt_batch2$SampleID),
    outliers_iqr = setdiff(dt_batch2$SampleID, dt_batch3$SampleID),
    outliers_all = all_outliers,
    kept_samples = final_kept,
    outlier_summary = outlier_summary
  ))
}

# Function to create PCA plots
create_pca_plots <- function(pca_scores, outlier_summary, metadata = NULL) {
  log_info("Creating PCA visualization plots")

  # Prepare plot data
  plot_data <- as.data.table(pca_scores[, seq_len(min(4, ncol(pca_scores)))] )
  plot_data$SampleID <- rownames(pca_scores)
  plot_data <- merge(plot_data, outlier_summary[, .(SampleID, is_outlier, outlier_score)], by = "SampleID")

  # Add metadata if available
  if(!is.null(metadata)) {
    # Match by SAMPLE_ID first
    plot_data <- merge(plot_data, metadata[, .(SAMPLE_ID, FINNGENID, COHORT_FINNGENID, BIOBANK_PLASMA)],
                       by.x = "SampleID", by.y = "SAMPLE_ID", all.x = TRUE)
  }

  # Use biobank for coloring (matching Python implementation). Drop NAs.
  plot_data[, sample_source := BIOBANK_PLASMA]
  plot_data <- plot_data[!is.na(sample_source)]

  # Build dynamic palette with all observed biobanks
  unique_biobanks <- sort(unique(plot_data$sample_source))
  biobank_palette <- setNames(get_husl_palette(length(unique_biobanks)), unique_biobanks)

  # Create shape mapping (different shapes for outliers vs normal)
  plot_data[, point_shape := ifelse(is_outlier, "Outlier", "Normal")]

  # Create annotation labels (FINNGENID if available, otherwise SampleID)
  plot_data[, annotation_label := ifelse(!is.na(FINNGENID), FINNGENID, SampleID)]

  # PC1 vs PC2
  p1 <- ggplot(plot_data, aes(x = PC1, y = PC2)) +
    geom_point(aes(color = sample_source, shape = point_shape), alpha = 0.7, size = 2) +
    scale_color_manual(values = biobank_palette, drop = FALSE) +
    scale_shape_manual(values = c("Normal" = 16, "Outlier" = 17)) +
    geom_text_repel(data = plot_data[is_outlier == TRUE],
                    aes(label = annotation_label), size = 3, max.overlaps = 15,
                    box.padding = 0.5, point.padding = 0.3) +
    labs(title = "PCA: PC1 vs PC2 by Biobank",
         subtitle = "Outliers labeled with FINNGENIDs and highlighted with different shapes",
         x = "PC1", y = "PC2",
         color = "Biobank", shape = "Sample Type") +
    theme_bw() +
    theme(legend.position = "right")

  # PC3 vs PC4 if available
  if(ncol(pca_scores) >= 4) {
    p2 <- ggplot(plot_data, aes(x = PC3, y = PC4)) +
      geom_point(aes(color = sample_source, shape = point_shape), alpha = 0.7, size = 2) +
      scale_color_manual(values = biobank_palette, drop = FALSE) +
      scale_shape_manual(values = c("Normal" = 16, "Outlier" = 17)) +
      geom_text_repel(data = plot_data[is_outlier == TRUE],
                      aes(label = annotation_label), size = 3, max.overlaps = 15,
                      box.padding = 0.5, point.padding = 0.3) +
      labs(title = "PCA: PC3 vs PC4 by Biobank",
           subtitle = "Outliers labeled with FINNGENIDs and highlighted with different shapes",
           x = "PC3", y = "PC4",
           color = "Biobank", shape = "Sample Type") +
      theme_bw() +
      theme(legend.position = "right")
  } else {
    p2 <- NULL
  }

  # Outlier score distribution
  p3 <- ggplot(outlier_summary, aes(x = outlier_score)) +
    geom_histogram(bins = max(outlier_summary$outlier_score) + 1,
                   fill = "#5C9EAD", alpha = 0.7) +
    labs(title = "Outlier Score Distribution",
         x = "Number of PCs flagged as outlier",
         y = "Count") +
    theme_bw()

  # PC1 distribution with outliers highlighted
  p4 <- ggplot(plot_data, aes(x = PC1, fill = sample_source)) +
    geom_histogram(bins = 50, alpha = 0.7, position = "identity") +
    scale_fill_manual(values = biobank_palette, drop = FALSE) +
    labs(title = "PC1 Distribution by Biobank",
         subtitle = "Histogram showing distribution of PC1 values",
         x = "PC1", y = "Count",
         fill = "Biobank") +
    theme_bw() +
    theme(legend.position = "right")

  return(list(
    pc1_pc2 = p1,
    pc3_pc4 = p2,
    outlier_score = p3,
    pc1_dist = p4
  ))
}

# Build disease group column from the last 12 indicator columns (if present)
add_disease_group <- function(metadata) {
  if (is.null(metadata)) return(NULL)
  disease_cols_expected <- c(
    "Kidney", "Kids", "F64", "MFGE8", "Parkinsons", "Metabolic",
    "AMD", "Rheuma", "Pulmo", "Chromosomal_Abnormalities",
    "Blood_donors", "Bridging_samples"
  )
  ind_cols <- intersect(disease_cols_expected, colnames(metadata))
  if (length(ind_cols) == 0) return(metadata)

  dt <- as.data.table(metadata)
  long <- melt(dt[, c("SAMPLE_ID", ind_cols), with = FALSE],
               id.vars = "SAMPLE_ID", variable.name = "Group",
               value.name = "flag", variable.factor = FALSE)

  long[, flag_s := tolower(as.character(flag))]
  suppressWarnings(long[, flag_num := as.numeric(flag_s)])
  long[, flag_bool := fifelse(flag_s %in% c("1", "yes", "true", "y", "t"), TRUE,
                              fifelse(!is.na(flag_num) & flag_num > 0, TRUE, FALSE))]

  grp <- long[flag_bool == TRUE, .(groups = list(unique(Group))), by = SAMPLE_ID]
  dt <- merge(dt, grp, by = "SAMPLE_ID", all.x = TRUE)
  dt[, DISEASE_GROUP := fifelse(is.null(groups) | is.na(groups), "None",
                                ifelse(lengths(groups) == 1, vapply(groups, function(x) x[1], ""), "Multiple"))]
  dt[, groups := NULL]
  return(dt)
}

# Create pair of PC plots colored by an arbitrary grouping column (no labels)
create_pc_plots_by_group <- function(pca_scores, outlier_summary, metadata, group_col, palette_fn = NULL) {
  if (is.null(metadata) || !(group_col %in% colnames(metadata))) return(list(pc1_pc2 = NULL, pc3_pc4 = NULL))

  plot_data <- as.data.table(pca_scores[, seq_len(min(4, ncol(pca_scores)))] )
  plot_data$SampleID <- rownames(pca_scores)
  plot_data <- merge(plot_data, outlier_summary[, .(SampleID, is_outlier)], by = "SampleID")
  plot_data <- merge(plot_data, metadata[, .(SAMPLE_ID, g = get(group_col))], by.x = "SampleID", by.y = "SAMPLE_ID", all.x = TRUE)
  setnames(plot_data, "g", group_col)
  plot_data <- plot_data[!is.na(get(group_col))]

  # Dynamic palette
  groups <- sort(unique(plot_data[[group_col]]))
  pal_vec <- if (is.null(palette_fn)) get_husl_palette(length(groups)) else palette_fn(length(groups))
  pal <- setNames(pal_vec, groups)

  p1 <- ggplot(plot_data, aes(x = PC1, y = PC2)) +
    geom_point(aes(color = .data[[group_col]], shape = ifelse(is_outlier, "Outlier", "Normal")), alpha = 0.7, size = 2) +
    scale_color_manual(values = pal, drop = FALSE) +
    scale_shape_manual(values = c("Normal" = 16, "Outlier" = 17)) +
    labs(title = paste0("PCA: PC1 vs PC2 by ", group_col), x = "PC1", y = "PC2", color = group_col, shape = "Sample Type") +
    theme_bw()

  p2 <- NULL
  if (ncol(pca_scores) >= 4) {
    p2 <- ggplot(plot_data, aes(x = PC3, y = PC4)) +
      geom_point(aes(color = .data[[group_col]], shape = ifelse(is_outlier, "Outlier", "Normal")), alpha = 0.7, size = 2) +
      scale_color_manual(values = pal, drop = FALSE) +
      scale_shape_manual(values = c("Normal" = 16, "Outlier" = 17)) +
      labs(title = paste0("PCA: PC3 vs PC4 by ", group_col), x = "PC3", y = "PC4", color = group_col, shape = "Sample Type") +
      theme_bw()
  }

  return(list(pc1_pc2 = p1, pc3_pc4 = p2))
}

# Function to remove outliers
remove_outliers <- function(npx_matrix, outliers_to_remove) {
  log_info("Removing {length(outliers_to_remove)} outliers from matrix")

  # Filter matrix
  clean_matrix <- npx_matrix[!rownames(npx_matrix) %in% outliers_to_remove, ]

  log_info("Matrix dimensions after outlier removal: {nrow(clean_matrix)} x {ncol(clean_matrix)}")

  return(clean_matrix)
}

# --- Utilities for outlier source matrix and QC summary (defined before main) ---

write_outlier_sources <- function(outliers_orig, metadata) {
  dt <- data.table(SampleID = unique(c(outliers_orig$outliers_pc, outliers_orig$outliers_pc34, outliers_orig$outliers_median, outliers_orig$outliers_iqr)))
  if (nrow(dt) == 0) return(invisible(NULL))
  dt[, PC1_PC2 := as.integer(SampleID %in% outliers_orig$outliers_pc)]
  dt[, PC3_PC4 := as.integer(SampleID %in% outliers_orig$outliers_pc34)]
  dt[, Median := as.integer(SampleID %in% outliers_orig$outliers_median)]
  dt[, IQR := as.integer(SampleID %in% outliers_orig$outliers_iqr)]
  dt[, Any := as.integer(PC1_PC2 + PC3_PC4 + Median + IQR > 0)]
  if (!is.null(metadata)) {
    # Ensure disease group is present where possible
    metadata_with_disease <- add_disease_group(metadata)
    cols <- intersect(c("SAMPLE_ID", "FINNGENID", "BIOBANK_PLASMA", "PlateID", "DISEASE_GROUP"), colnames(metadata_with_disease))
    if (length(cols) > 0 && "SAMPLE_ID" %in% cols) {
      meta_sel <- metadata_with_disease[, ..cols]
      dt <- merge(dt, meta_sel, by.x = "SampleID", by.y = "SAMPLE_ID", all.x = TRUE)
    }
  }
  setcolorder(dt, c("FINNGENID", setdiff(names(dt), c("FINNGENID"))))
  path <- get_output_path(step_num, "pca_outliers_by_source", batch_id, "outliers", "tsv", config = config)
  ensure_output_dir(path)
  fwrite(dt, path, sep = "\t")
}

write_initial_qced_matrix <- function(metadata) {
  # Write analysis-ready, pre-outlier NPX matrix with metadata columns
  prev_step00_num <- "00"
  mat_path <- get_output_path(prev_step00_num, "npx_matrix_analysis_ready", batch_id, "qc", config = config)
  mat <- readRDS(mat_path)
  if (is.null(dim(mat))) return(invisible(NULL))
  df <- as.data.table(mat)
  df[, SampleID := rownames(mat)]
  setcolorder(df, c("SampleID", setdiff(names(df), "SampleID")))
  meta_cols <- c("SAMPLE_ID", "FINNGENID", "APPROX_TIMESTAMP_COLLECTION", "PlateID", "WellID", "COHORT_FINNGENID", "BIOBANK_PLASMA")
  meta_avail <- intersect(meta_cols, colnames(metadata))
  meta <- metadata[, ..meta_avail]
  setnames(meta, old = "SAMPLE_ID", new = "SampleID")
  merged <- merge(df, meta, by = "SampleID", all.x = TRUE)
  # Reorder to put metadata columns first
  front <- c("SampleID", "FINNGENID", "APPROX_TIMESTAMP_COLLECTION", "PlateID", "WellID", "COHORT_FINNGENID", "BIOBANK_PLASMA")
  front <- intersect(front, names(merged))
  merged <- merged[, c(front, setdiff(names(merged), front)), with = FALSE]
  # Use batch designation for filename
  batch_designation <- config$batches[[batch_id]]$batch_designation %||% paste0("batch_", gsub("batch_", "", batch_id))
  path <- get_output_path(step_num, paste0("olink5k_", batch_designation, "_QCed"), batch_id, "qc", "tsv", config = config)
  ensure_output_dir(path)
  fwrite(merged, path, sep = "\t")
}

write_qc_summary <- function() {
  # Robustly assemble requested summary
  out_path <- get_output_path(step_num, "pipeline_qc_summary", batch_id, "reports", "txt", config = config)
  ensure_output_dir(out_path)

  # Load what's available (batch-aware paths)
  prev_step00_num <- "00"
  prev_step04_num <- "04"

  sample_mapping_path <- get_output_path(prev_step00_num, "sample_mapping", batch_id, "qc", config = config)
  analysis_samples_path <- get_output_path(prev_step00_num, "analysis_samples", batch_id, "qc", config = config)
  mapping_validation_path <- get_output_path(prev_step00_num, "mapping_validation", batch_id, "qc", config = config)
  outliers_orig_path <- get_output_path(step_num, "pca_outliers_original", batch_id, "outliers", config = config)
  sex_outliers_path <- get_output_path(prev_step04_num, "sex_outliers", batch_id, "outliers", "tsv", config = config)

  sample_mapping <- try(readRDS(sample_mapping_path), silent = TRUE)
  analysis_samples <- try(readRDS(analysis_samples_path), silent = TRUE)
  mapping_validation <- try(readRDS(mapping_validation_path), silent = TRUE)
  outliers_orig <- try(readRDS(outliers_orig_path), silent = TRUE)

  n_received <- if (!inherits(sample_mapping, "try-error")) nrow(sample_mapping) else NA_integer_
  n_valid <- if (!inherits(analysis_samples, "try-error")) nrow(analysis_samples) else NA_integer_
  n_dups <- if (!inherits(mapping_validation, "try-error")) nrow(mapping_validation$duplicate_finngenids) else NA_integer_
  pc1pc2_rm <- if (!inherits(outliers_orig, "try-error")) length(outliers_orig$outliers_pc) else NA_integer_
  pc3pc4_rm <- if (!inherits(outliers_orig, "try-error")) length(outliers_orig$outliers_pc34) else NA_integer_
  median_rm <- if (!inherits(outliers_orig, "try-error")) length(outliers_orig$outliers_median) else NA_integer_
  iqr_rm <- if (!inherits(outliers_orig, "try-error")) length(outliers_orig$outliers_iqr) else NA_integer_
  n_normal <- if (!inherits(outliers_orig, "try-error")) length(outliers_orig$kept_samples) else NA_integer_
  sex_outliers_n <- if (file.exists(sex_outliers_path)) nrow(fread(sex_outliers_path)) else NA_integer_

  txt <- c(
    "# Sample QC:",
    sprintf("%s samples received in this batch.", ifelse(is.na(n_received), "NA", n_received)),
    sprintf("    %s valida samples", ifelse(is.na(n_valid), "NA", n_valid)),
    sprintf("        %s normal samples", ifelse(is.na(n_normal), "NA", n_normal)),
    "",
    sprintf("Start from %s valid samples (AG samples already removed)", ifelse(is.na(n_valid), "NA", n_valid)),
    sprintf("    PCA to remove PC1 and PC2 out of range mean +- 5SD: %s", ifelse(is.na(pc1pc2_rm), "NA", pc1pc2_rm)),
    sprintf("    PCA to remove PC3 and PC4 out of range mean +- 5SD: %s", ifelse(is.na(pc3pc4_rm), "NA", pc3pc4_rm)),
    sprintf("    Remove median protein expression level +- 5SD: %s", ifelse(is.na(median_rm), "NA", median_rm)),
    sprintf("    Remove protein expression level IQR +- 5SD: %s", ifelse(is.na(iqr_rm), "NA", iqr_rm)),
    sprintf("    Initial QC remain sample: %s (%s duplicates), saved to file (olink5k_batch2_QCed.tsv)",
            ifelse(is.na(n_valid), "NA", n_valid), ifelse(is.na(n_dups), "NA", n_dups)),
    "",
    "Further outlier detectin:",
    sprintf("    -Construct sex prediction model ... outlier: %s", ifelse(is.na(sex_outliers_n), "NA", sex_outliers_n)),
    "    - Detect sample swaps by utilizing strong cis pQTLs ...: NA",
    "",
    "pQTL running set:",
    "    Remove duplicate samples: N = NA",
    "    Keep only unrelated samples: N = NA (olink5k1_recover_unrel_rint.pheno, preadjust with age, sex, cohort, plate, row, col, PC1 - PC10, and then perform RINT)",
    "",
    "## File formats for protein expression:",
    "tsv:",
    "    SampleID: sample raw ID",
    "    FINNGENID: mapped FinnGen ID",
    "    APPROX_TIMESTAMP_COLLECTION: time when the plasma sample was collected",
    "    PlateID: plate ID of Olink HT",
    "    WellID: the well ID where sample was located in the plate",
    "    CORHORT_FINNGEN: the cohort information where the sample come from",
    "    BIOBANK_PLASMA: Biobank handle the plasma sample",
    "    Probe name (e.g. A1BG  A1CF...):  the protein name, and the following rows as the protein expression level for that probe",
    "",
    "pheno:",
    "    FID: sample ID",
    "    IID: sample ID",
    "    Probe name (e.g. A1BG  A1CF...):  the protein name, and the following rows as the protein expression level for that probe",
    "",
    "# Probe QC:",
    "xxx valid probes excluding control probes (fg3_batch3_probe_all.txt)"
  )
  writeLines(txt, out_path)
}

# Overlay proteomic PCs with genetic PCs extracted like in 05_sex_outliers.R
create_overlay_genetic_pc_plots <- function(pca_scores, outlier_summary, metadata, covariate_file) {
  # Prepare proteomic PCs
  df <- as.data.table(pca_scores[, seq_len(min(4, ncol(pca_scores)))])
  df$SampleID <- rownames(pca_scores)
  df <- merge(df, outlier_summary[, .(SampleID, is_outlier)], by = "SampleID", all.x = TRUE)
  df <- merge(df, metadata[, .(SAMPLE_ID, FINNGENID, BIOBANK_PLASMA)], by.x = "SampleID", by.y = "SAMPLE_ID", all.x = TRUE)

  # Load genetic PCs (as in 05_sex_outliers.R)
  covar <- try(fread(cmd = paste("zcat", covariate_file)), silent = TRUE)
  if (inherits(covar, "try-error")) return(list(pc12 = NULL, pc34 = NULL))
  gen <- covar[, .(FINNGENID = IID, gPC1 = PC1, gPC2 = PC2, gPC3 = PC3, gPC4 = PC4)]

  df <- merge(df, gen, by = "FINNGENID", all.x = TRUE)

  # Quintiles for shapes
  qlab <- c("Q1","Q2","Q3","Q4","Q5")
  if (!all(is.na(df$gPC2))) {
    br2 <- quantile(df$gPC2, probs = seq(0,1,0.2), na.rm = TRUE)
    df[, gPC2_q := cut(gPC2, breaks = unique(br2), include.lowest = TRUE, labels = qlab)]
  } else df[, gPC2_q := NA_character_]
  if (!all(is.na(df$gPC4))) {
    br4 <- quantile(df$gPC4, probs = seq(0,1,0.2), na.rm = TRUE)
    df[, gPC4_q := cut(gPC4, breaks = unique(br4), include.lowest = TRUE, labels = qlab)]
  } else df[, gPC4_q := NA_character_]

  # Color palette for continuous genetic PCs
  vir_pal <- if (requireNamespace("viridisLite", quietly = TRUE)) viridisLite::viridis(256) else scales::hue_pal()(256)

  # PC1 vs PC2 colored by gPC1, shaped by gPC2 quintile
  p12 <- ggplot(df, aes(x = PC1, y = PC2)) +
    geom_point(aes(color = gPC1, shape = gPC2_q), alpha = 0.7, size = 2, na.rm = TRUE) +
    scale_color_gradientn(colors = vir_pal, na.value = "#CCCCCC") +
    scale_shape_manual(values = c(16, 17, 15, 3, 4), na.translate = TRUE) +
    geom_text_repel(data = df[is_outlier == TRUE & !is.na(FINNGENID)],
                    aes(label = FINNGENID), size = 3, max.overlaps = 20,
                    box.padding = 0.5, point.padding = 0.3) +
    labs(title = "Proteomic PC1 vs PC2 overlay with Genetic PCs",
         subtitle = "Color = Genetic PC1, Shape = Genetic PC2 quintile",
         x = "Proteomic PC1", y = "Proteomic PC2",
         color = "Genetic PC1", shape = "Genetic PC2 quintile") +
    theme_bw()

  # PC3 vs PC4 colored by gPC3, shaped by gPC4 quintile
  p34 <- NULL
  if (ncol(pca_scores) >= 4) {
    p34 <- ggplot(df, aes(x = PC3, y = PC4)) +
      geom_point(aes(color = gPC3, shape = gPC4_q), alpha = 0.7, size = 2, na.rm = TRUE) +
      scale_color_gradientn(colors = vir_pal, na.value = "#CCCCCC") +
      scale_shape_manual(values = c(16, 17, 15, 3, 4), na.translate = TRUE) +
      geom_text_repel(data = df[is_outlier == TRUE & !is.na(FINNGENID)],
                      aes(label = FINNGENID), size = 3, max.overlaps = 20,
                      box.padding = 0.5, point.padding = 0.3) +
      labs(title = "Proteomic PC3 vs PC4 overlay with Genetic PCs",
           subtitle = "Color = Genetic PC3, Shape = Genetic PC4 quintile",
           x = "Proteomic PC3", y = "Proteomic PC4",
           color = "Genetic PC3", shape = "Genetic PC4 quintile") +
      theme_bw()
  }

  list(pc12 = p12, pc34 = p34)
}

# Genetic-only PC scatter plots for proteomics-linked samples, colored by biobank/disease
create_genetic_pc_plots <- function(pca_scores, outliers_orig, metadata, covariate_file) {
  # Restrict to samples with proteomics PCs
  proteo_samples <- rownames(pca_scores)
  meta_sub <- metadata[SAMPLE_ID %in% proteo_samples, .(SAMPLE_ID, FINNGENID, BIOBANK_PLASMA)]
  meta_sub <- unique(meta_sub, by = "FINNGENID")

  # Add disease group
  md_dis <- add_disease_group(metadata)
  dis_sub <- md_dis[SAMPLE_ID %in% proteo_samples, .(SAMPLE_ID, FINNGENID, DISEASE_GROUP)]
  dis_sub <- unique(dis_sub, by = "FINNGENID")

  # Genetic PCs
  covar <- try(fread(cmd = paste("zcat", covariate_file)), silent = TRUE)
  if (inherits(covar, "try-error")) return(list(biobank_pc12 = NULL, biobank_pc34 = NULL, disease_pc12 = NULL, disease_pc34 = NULL))
  gen <- covar[, .(FINNGENID = IID, gPC1 = PC1, gPC2 = PC2, gPC3 = PC3, gPC4 = PC4)]

  df <- merge(gen, meta_sub[, .(FINNGENID, BIOBANK_PLASMA)], by = "FINNGENID", all.x = TRUE)
  df <- merge(df, dis_sub[, .(FINNGENID, DISEASE_GROUP)], by = "FINNGENID", all.x = TRUE)

  # Outliers by FINNGENID
  outlier_finngen <- unique(metadata[SAMPLE_ID %in% outliers_orig$outliers_all, FINNGENID])
  df[, is_outlier := FINNGENID %in% outlier_finngen]

  # Palettes
  biobanks <- sort(unique(na.omit(df$BIOBANK_PLASMA)))
  pal_biobank <- setNames(get_husl_palette(length(biobanks)), biobanks)
  dis <- sort(unique(na.omit(df$DISEASE_GROUP)))
  pal_dis <- setNames(get_wes_cont_palette(length(dis)), dis)

  p_bio_12 <- ggplot(df, aes(x = gPC1, y = gPC2)) +
    geom_point(data = df[is.na(BIOBANK_PLASMA)], color = "#BDBDBD", alpha = 0.15, size = 1.5, na.rm = TRUE) +
    geom_point(data = df[!is.na(BIOBANK_PLASMA)], aes(color = BIOBANK_PLASMA), alpha = 0.7, size = 2, na.rm = TRUE) +
    scale_color_manual(values = pal_biobank, drop = TRUE, na.translate = FALSE) +
    geom_text_repel(data = df[is_outlier == TRUE & !is.na(FINNGENID) & !is.na(BIOBANK_PLASMA)], aes(label = FINNGENID),
                    size = 3, max.overlaps = 20, box.padding = 0.5, point.padding = 0.3) +
    labs(title = "Genetic PC1 vs PC2 (proteomics-linked)", x = "Genetic PC1", y = "Genetic PC2", color = "Biobank") +
    theme_bw()

  p_bio_34 <- ggplot(df, aes(x = gPC3, y = gPC4)) +
    geom_point(data = df[is.na(BIOBANK_PLASMA)], color = "#BDBDBD", alpha = 0.15, size = 1.5, na.rm = TRUE) +
    geom_point(data = df[!is.na(BIOBANK_PLASMA)], aes(color = BIOBANK_PLASMA), alpha = 0.7, size = 2, na.rm = TRUE) +
    scale_color_manual(values = pal_biobank, drop = TRUE, na.translate = FALSE) +
    geom_text_repel(data = df[is_outlier == TRUE & !is.na(FINNGENID) & !is.na(BIOBANK_PLASMA)], aes(label = FINNGENID),
                    size = 3, max.overlaps = 20, box.padding = 0.5, point.padding = 0.3) +
    labs(title = "Genetic PC3 vs PC4 (proteomics-linked)", x = "Genetic PC3", y = "Genetic PC4", color = "Biobank") +
    theme_bw()

  p_dis_12 <- ggplot(df, aes(x = gPC1, y = gPC2)) +
    geom_point(data = df[is.na(DISEASE_GROUP)], color = "#BDBDBD", alpha = 0.15, size = 1.5, na.rm = TRUE) +
    geom_point(data = df[!is.na(DISEASE_GROUP)], aes(color = DISEASE_GROUP), alpha = 0.7, size = 2, na.rm = TRUE) +
    scale_color_manual(values = pal_dis, drop = TRUE, na.translate = FALSE) +
    geom_text_repel(data = df[is_outlier == TRUE & !is.na(FINNGENID) & !is.na(DISEASE_GROUP)], aes(label = FINNGENID),
                    size = 3, max.overlaps = 20, box.padding = 0.5, point.padding = 0.3) +
    labs(title = "Genetic PC1 vs PC2 by disease (proteomics-linked)", x = "Genetic PC1", y = "Genetic PC2", color = "Disease") +
    theme_bw()

  p_dis_34 <- ggplot(df, aes(x = gPC3, y = gPC4)) +
    geom_point(data = df[is.na(DISEASE_GROUP)], color = "#BDBDBD", alpha = 0.15, size = 1.5, na.rm = TRUE) +
    geom_point(data = df[!is.na(DISEASE_GROUP)], aes(color = DISEASE_GROUP), alpha = 0.7, size = 2, na.rm = TRUE) +
    scale_color_manual(values = pal_dis, drop = TRUE, na.translate = FALSE) +
    geom_text_repel(data = df[is_outlier == TRUE & !is.na(FINNGENID) & !is.na(DISEASE_GROUP)], aes(label = FINNGENID),
                    size = 3, max.overlaps = 20, box.padding = 0.5, point.padding = 0.3) +
    labs(title = "Genetic PC3 vs PC4 by disease (proteomics-linked)", x = "Genetic PC3", y = "Genetic PC4", color = "Disease") +
    theme_bw()

  list(biobank_pc12 = p_bio_12, biobank_pc34 = p_bio_34, disease_pc12 = p_dis_12, disease_pc34 = p_dis_34)
}

# Main execution
main <- function() {

  # Load data from previous step (step 00, batch-aware paths)
  log_info("Loading data from previous step (step 00)")
  prev_step00_num <- "00"

  npx_matrix_path <- get_output_path(prev_step00_num, "npx_matrix_analysis_ready", batch_id, "qc", config = config)
  metadata_path <- get_output_path(prev_step00_num, "metadata", batch_id, "qc", config = config)

  if (!file.exists(npx_matrix_path)) {
    stop("NPX matrix file not found: {npx_matrix_path}. Run Step 00 first.")
  }
  if (!file.exists(metadata_path)) {
    stop("Metadata file not found: {metadata_path}. Run Step 00 first.")
  }

  npx_matrix <- readRDS(npx_matrix_path)
  metadata <- readRDS(metadata_path)

  # Note: AG samples are already removed in Step 00 (pre-filtered input)
  log_info("NPX matrix loaded: {nrow(npx_matrix)} samples x {ncol(npx_matrix)} proteins")
  log_info("Note: AG samples were already removed in pre-filtered input")

  # Perform PCA
  pca_result <- perform_pca(npx_matrix, n_components = config$parameters$outliers$pca_components)

  # Original (batch1) sequential outlier detection
  outliers_orig <- detect_outliers_original(npx_matrix, pca_result$pca_object)

  # Create plots
  pca_plots <- create_pca_plots(pca_result$scores, outliers_orig$outlier_summary, metadata)

  # Also create no-label plots by biobank and by disease group
  metadata_disease <- add_disease_group(metadata)
  biobank_nolabel <- create_pc_plots_by_group(pca_result$scores, outliers_orig$outlier_summary, metadata, "BIOBANK_PLASMA")
  disease_nolabel <- if (!is.null(metadata_disease) && "DISEASE_GROUP" %in% colnames(metadata_disease)) {
    create_pc_plots_by_group(pca_result$scores, outliers_orig$outlier_summary, metadata_disease, "DISEASE_GROUP", palette_fn = get_wes_cont_palette)
  } else { list(pc1_pc2 = NULL, pc3_pc4 = NULL) }

  # Also create labeled disease plots using Wes Anderson palette
  disease_labeled <- disease_nolabel

  # Remove outliers from matrix
  npx_clean <- remove_outliers(npx_matrix, outliers_orig$outliers_all)

  # Compute and save feature loadings
  protein_loadings_dt <- make_protein_loadings_dt(pca_result)
  tech_assoc_dt <- compute_technical_associations(pca_result$scores, metadata)

  # Save outputs with batch-aware paths
  log_info("Saving PCA outlier detection results")

  pca_result_path <- get_output_path(step_num, "pca_result", batch_id, "outliers", config = config)
  pca_outliers_original_path <- get_output_path(step_num, "pca_outliers_original", batch_id, "outliers", config = config)
  npx_matrix_pca_cleaned_path <- get_output_path(step_num, "npx_matrix_pca_cleaned", batch_id, "outliers", config = config)
  pc_protein_loadings_path <- get_output_path(step_num, "pc_protein_loadings", batch_id, "outliers", "tsv", config = config)
  pc_technical_loadings_path <- get_output_path(step_num, "pc_technical_loadings", batch_id, "outliers", "tsv", config = config)

  ensure_output_dir(pca_result_path)
  ensure_output_dir(pca_outliers_original_path)
  ensure_output_dir(npx_matrix_pca_cleaned_path)
  ensure_output_dir(pc_protein_loadings_path)
  ensure_output_dir(pc_technical_loadings_path)

  saveRDS(pca_result, pca_result_path)
  saveRDS(outliers_orig, pca_outliers_original_path)
  saveRDS(npx_clean, npx_matrix_pca_cleaned_path)
  if (!is.null(protein_loadings_dt)) fwrite(protein_loadings_dt, pc_protein_loadings_path, sep = "\t")
  if (!is.null(tech_assoc_dt)) fwrite(tech_assoc_dt, pc_technical_loadings_path, sep = "\t")

  # Save tables
  pca_outlier_summary_path <- get_output_path(step_num, "pca_outlier_summary", batch_id, "outliers", "tsv", config = config)
  pca_outliers_list_path <- get_output_path(step_num, "pca_outliers_list", batch_id, "outliers", "tsv", config = config)
  pca_outliers_by_step_path <- get_output_path(step_num, "pca_outliers_by_step", batch_id, "outliers", "tsv", config = config)

  ensure_output_dir(pca_outlier_summary_path)
  ensure_output_dir(pca_outliers_list_path)
  ensure_output_dir(pca_outliers_by_step_path)

  # Add FINNGENID to all outputs with SampleID
  log_info("Adding FINNGENID mapping to output tables...")
  outlier_summary_with_fgid <- add_finngenid_column(outliers_orig$outlier_summary, batch_id = batch_id, config = config)
  outliers_list_with_fgid <- add_finngenid_column(data.table(SampleID = outliers_orig$outliers_all), batch_id = batch_id, config = config)
  outliers_by_step_with_fgid <- add_finngenid_column(
    data.table(
    step = c(rep("PC1_PC2", length(outliers_orig$outliers_pc)),
             rep("PC3_PC4", length(outliers_orig$outliers_pc34)),
             rep("Median", length(outliers_orig$outliers_median)),
             rep("IQR", length(outliers_orig$outliers_iqr))),
    SampleID = c(outliers_orig$outliers_pc, outliers_orig$outliers_pc34, outliers_orig$outliers_median, outliers_orig$outliers_iqr)
    ),
    batch_id = batch_id, config = config
  )

  fwrite(outlier_summary_with_fgid, pca_outlier_summary_path, sep = "\t")
  fwrite(outliers_list_with_fgid, pca_outliers_list_path, sep = "\t")
  fwrite(outliers_by_step_with_fgid, pca_outliers_by_step_path, sep = "\t")

  # Also update the RDS object with FINNGENID
  outliers_orig$outlier_summary <- outlier_summary_with_fgid
  saveRDS(outliers_orig, pca_outliers_original_path)

  # Save plots
  pca_pc1_pc2_path <- get_output_path(step_num, "pca_pc1_pc2", batch_id, "outliers", "pdf", config = config)
  pca_pc3_pc4_path <- get_output_path(step_num, "pca_pc3_pc4", batch_id, "outliers", "pdf", config = config)
  pca_outlier_score_path <- get_output_path(step_num, "pca_outlier_score", batch_id, "outliers", "pdf", config = config)
  pca_pc1_distribution_path <- get_output_path(step_num, "pca_pc1_distribution", batch_id, "outliers", "pdf", config = config)

  ensure_output_dir(pca_pc1_pc2_path)
  ensure_output_dir(pca_pc3_pc4_path)
  ensure_output_dir(pca_outlier_score_path)
  ensure_output_dir(pca_pc1_distribution_path)

  ggsave(pca_pc1_pc2_path, pca_plots$pc1_pc2, width = 10, height = 8)
  if(!is.null(pca_plots$pc3_pc4)) {
    ggsave(pca_pc3_pc4_path, pca_plots$pc3_pc4, width = 10, height = 8)
  }
  ggsave(pca_outlier_score_path, pca_plots$outlier_score, width = 8, height = 6)
  ggsave(pca_pc1_distribution_path, pca_plots$pc1_dist, width = 8, height = 6)

  # Save no-label plots for both coloring schemes
  if (!is.null(biobank_nolabel$pc1_pc2)) {
    path <- get_output_path(step_num, "pca_pc1_pc2_biobank_nolabel", batch_id, "outliers", "pdf", config = config)
    ensure_output_dir(path)
    ggsave(path, biobank_nolabel$pc1_pc2, width = 10, height = 8)
  }
  if (!is.null(biobank_nolabel$pc3_pc4)) {
    path <- get_output_path(step_num, "pca_pc3_pc4_biobank_nolabel", batch_id, "outliers", "pdf", config = config)
    ensure_output_dir(path)
    ggsave(path, biobank_nolabel$pc3_pc4, width = 10, height = 8)
  }
  if (!is.null(disease_nolabel$pc1_pc2)) {
    path <- get_output_path(step_num, "pca_pc1_pc2_disease_nolabel", batch_id, "outliers", "pdf", config = config)
    ensure_output_dir(path)
    ggsave(path, disease_nolabel$pc1_pc2, width = 10, height = 8)
  }
  if (!is.null(disease_nolabel$pc3_pc4)) {
    path <- get_output_path(step_num, "pca_pc3_pc4_disease_nolabel", batch_id, "outliers", "pdf", config = config)
    ensure_output_dir(path)
    ggsave(path, disease_nolabel$pc3_pc4, width = 10, height = 8)
  }

  # Save labeled disease plots
  if (!is.null(disease_labeled$pc1_pc2)) {
    path <- get_output_path(step_num, "pca_pc1_pc2_disease", batch_id, "outliers", "pdf", config = config)
    ensure_output_dir(path)
    ggsave(path, disease_labeled$pc1_pc2, width = 10, height = 8)
  }
  if (!is.null(disease_labeled$pc3_pc4)) {
    path <- get_output_path(step_num, "pca_pc3_pc4_disease", batch_id, "outliers", "pdf", config = config)
    ensure_output_dir(path)
    ggsave(path, disease_labeled$pc3_pc4, width = 10, height = 8)
  }

  # [Removed] Proteomic vs genetic PC overlay plot exports (not needed)

  # Genetic-only PC plots colored by biobank and disease
  genplots <- create_genetic_pc_plots(pca_result$scores, outliers_orig, metadata, config$covariates$covariate_file)
  if (!is.null(genplots$biobank_pc12)) {
    path <- get_output_path(step_num, "genetic_pc1_pc2_biobank", batch_id, "outliers", "pdf", config = config)
    ensure_output_dir(path)
    ggsave(path, genplots$biobank_pc12, width = 10, height = 8)
  }
  if (!is.null(genplots$biobank_pc34)) {
    path <- get_output_path(step_num, "genetic_pc3_pc4_biobank", batch_id, "outliers", "pdf", config = config)
    ensure_output_dir(path)
    ggsave(path, genplots$biobank_pc34, width = 10, height = 8)
  }
  if (!is.null(genplots$disease_pc12)) {
    path <- get_output_path(step_num, "genetic_pc1_pc2_disease", batch_id, "outliers", "pdf", config = config)
    ensure_output_dir(path)
    ggsave(path, genplots$disease_pc12, width = 10, height = 8)
  }
  if (!is.null(genplots$disease_pc34)) {
    path <- get_output_path(step_num, "genetic_pc3_pc4_disease", batch_id, "outliers", "pdf", config = config)
    ensure_output_dir(path)
    ggsave(path, genplots$disease_pc34, width = 10, height = 8)
  }

  # Save outlier sources matrix and initial QCed matrix
  write_outlier_sources(outliers_orig, metadata)
  write_initial_qced_matrix(metadata)
  write_qc_summary()

  # Print summary
  cat("\n=== PCA OUTLIER DETECTION SUMMARY ===\n")
  cat("PCA components used:", config$parameters$outliers$pca_components, "\n")
  cat("Variance explained by PCs:", round(pca_result$cum_var_explained[config$parameters$outliers$pca_components] * 100, 2), "%\n")
  cat("\nOutlier detection (original procedure):\n")
  cat("  - PC1/PC2 ±5 SD removals:", length(outliers_orig$outliers_pc), "\n")
  cat("  - Median ±5 SD removals:", length(outliers_orig$outliers_median), "\n")
  cat("  - IQR ±5 SD removals:", length(outliers_orig$outliers_iqr), "\n")
  cat("  - Total removed:", length(outliers_orig$outliers_all), "\n")
  cat("\nMatrix after outlier removal:", nrow(npx_clean), "samples x", ncol(npx_clean), "proteins\n")
  cat("Results saved to: output/outliers/", if(batch_id != config$batch$default_batch_id ||
    tryCatch(isTRUE(config$parameters$normalization$multi_batch_mode), error = function(e) FALSE))
    paste0(batch_id, "/") else "", "\n", sep = "")

  log_info("PCA outlier detection completed")

  return(list(
    pca_result = pca_result,
    outliers = outliers_orig,
    clean_matrix = npx_clean
  ))
}

# Run if executed directly
if (!interactive()) {
  result <- main()
}

