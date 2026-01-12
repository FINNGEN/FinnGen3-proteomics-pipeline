#!/usr/bin/env Rscript

#################################################
# Script: 08_covariate_adjustment.R
# Author: Reza Jabal, PhD (rjabal@broadinstitute.org)
# Description: Adjust proteomics data for covariates
#              Refactored version
#              Refactored version
# Date: December 2025
#################################################

suppressPackageStartupMessages({
  library(data.table)
  library(tidyverse)
  library(ggplot2)
  library(ggrepel)
  library(gridExtra)
  library(paletteer)
  library(yaml)
  library(logger)
})

# Source helper functions for multicollinearity-aware analysis
source("file.path(script_dir, "08_helper_multicollinearity_adjustment.R")")

# Suppress "no visible binding" warnings for data.table operations
utils::globalVariables(
  c("IID", "SAMPLE_ID", "FINNGENID", "AGE_AT_DEATH_OR_END_OF_FOLLOWUP",
    "SEX_IMPUTED", "BMI", "harmonized_current_smoker", "SampleID",
    "correlation", "t_statistic", "before", "after", ".", "stage",
    "protein", "reduction", "label", "effect", "reduction_pct",
    "covariate", "PC1", "PC2", "..level..", "component", "variance_pct",
    "protein_idx", "r2_age", "r2_sex", "r2_full", "r2_residual",
    "APPROX_BIRTH_DATE", "APPROX_TIMESTAMP_COLLECTION", "BL_AGE",
    "birth_date", "bl_age", "sample_age", "neg_log10_p", "direction",
    "n_proteins", "mean_r2", "t_stat", "p_value", "abs_corr",
    "r2_ppc1", "r2_ppc2", "beta", "se", "bonferroni_sig", "abs_beta",
    "age", "young_sample_ids", "age_filtered_samples", "group_id",
    "group_size", "selected", "r_squared", "n_samples", "ci_lower",
    "ci_upper", "bonferroni_flag")
)

# Set up logging
log_appender(appender_file())
log_info("Starting covariate adjustment")

# Load configuration
config <- read_yaml("")

# Set theme for plots
theme_set(theme_bw())

# Load proteomic PCs from 04_pca_outliers (scores matrix)
load_proteomic_pcs <- function(pca_result_file, n_pcs = 10) {
  log_info("Loading proteomic PCs from: {pca_result_file}")
  pr <- readRDS(pca_result_file)
  scores <- as.data.table(pr$scores)
  scores$SampleID <- rownames(pr$scores)
  pcs <- names(scores)[grepl("^PC[0-9]+$", names(scores))]
  pcs <- pcs[seq_len(min(n_pcs, length(pcs)))]
  out <- scores[, c("SampleID", pcs), with = FALSE]
  setnames(out, pcs, paste0("pPC", seq_along(pcs)))
  out
}

# Function to load and prepare covariates
prepare_covariates <- function(covariate_file, sample_ids, metadata, proteomic_pcs) {
  log_info("Preparing covariates for adjustment (with proteomic PCs)")

  # Load covariates
  covariates <- fread(cmd = paste("zcat", covariate_file))

  # Load sex information with birth date
  sex_file <- "/mnt/longGWAS_disk_100GB/long_gwas/db_config/finngen_R13_minimum_1.0.txt.gz"
  sex_data <- fread(cmd = paste("zcat", sex_file))

  # Extract birth date and calculate age at sampling
  sex_info <- data.table(
    FINNGENID = sex_data$FINNGENID,
    birth_date = if("APPROX_BIRTH_DATE" %in% names(sex_data)) as.POSIXct(sex_data$APPROX_BIRTH_DATE, format="%Y-%m-%d") else NA,
    bl_age = if("BL_AGE" %in% names(sex_data)) sex_data$BL_AGE else NA_real_
  )

  # Calculate sample collection age if metadata has timestamps
  if ("APPROX_TIMESTAMP_COLLECTION" %in% names(metadata) && !all(is.na(sex_info$birth_date))) {
    id_map <- metadata[, .(SAMPLE_ID, FINNGENID, APPROX_TIMESTAMP_COLLECTION)]
    age_calc <- merge(id_map, sex_info, by = "FINNGENID", all.x = TRUE, allow.cartesian = TRUE)

    # Calculate age at sample collection for each SAMPLE_ID
    age_calc[!is.na(APPROX_TIMESTAMP_COLLECTION) & !is.na(birth_date),
             sample_age := as.numeric(difftime(APPROX_TIMESTAMP_COLLECTION, birth_date, units = "days") / 365.25)]

    # Use sample age as primary, fallback to BL_AGE
    age_calc[is.na(sample_age), sample_age := bl_age]

    age_data <- age_calc[, .(SAMPLE_ID, FINNGENID, age = sample_age)]
  } else {
    # Fallback to BL_AGE if no timestamps
    sample_mapping_temp <- metadata[SAMPLE_ID %in% sample_ids, .(SAMPLE_ID, FINNGENID)]
    age_data <- merge(sample_mapping_temp, sex_info[, .(FINNGENID, age = bl_age)], by = "FINNGENID", all.x = TRUE)
  }

  # Map samples to FINNGENIDs
  sample_mapping <- metadata[SAMPLE_ID %in% sample_ids, .(SAMPLE_ID, FINNGENID)]

  # Merge covariates with sample IDs (WITHOUT genomic PCs - we'll use proteomic PCs instead)
  covariate_data <- merge(
    sample_mapping,
    covariates[, .(
      FINNGENID = IID,
      sex = SEX_IMPUTED,
      bmi = BMI,
      smoking = harmonized_current_smoker
    )],
    by = "FINNGENID",
    all.x = TRUE
  )

  # Add calculated age at sampling
  covariate_data <- merge(covariate_data, age_data[, .(SAMPLE_ID, age)], by = "SAMPLE_ID", all.x = TRUE)

  # Filter out samples with age <= 20 years (high leverage pediatric/young adult samples)
  # These samples have fundamentally different proteomes and cause extreme effect sizes in age regressions
  n_young <- sum(covariate_data$age <= 20, na.rm = TRUE)
  if (n_young > 0) {
    log_info("Excluding {n_young} samples with age <= 20 years for age association analyses")
    young_sample_ids <- covariate_data[age <= 20 & !is.na(age), SAMPLE_ID]
  } else {
    young_sample_ids <- character(0)
  }

  # Align with sample order
  covariate_data <- covariate_data[match(sample_ids, covariate_data$SAMPLE_ID)]

  # Merge with proteomic PCs (align by SAMPLE_ID)
  pcs_aligned <- merge(data.table(SAMPLE_ID = sample_ids), proteomic_pcs,
                       by.x = "SAMPLE_ID", by.y = "SampleID", all.x = TRUE)
  # Reorder to exactly match sample_ids
  pcs_aligned <- pcs_aligned[match(sample_ids, pcs_aligned$SAMPLE_ID)]

  # Add proteomic PCs to covariate_data
  pc_cols <- grep("^pPC[0-9]+$", names(pcs_aligned), value = TRUE)
  for (col in pc_cols) {
    covariate_data[[col]] <- pcs_aligned[[col]]
  }

  # Check completeness (with proteomic PCs)
  n_check_pcs <- min(5, length(pc_cols))
  check_cols <- if (n_check_pcs > 0) {
    c("age", "sex", pc_cols[seq_len(n_check_pcs)])
  } else {
    c("age", "sex")
  }
  complete_samples <- complete.cases(covariate_data[, ..check_cols])

  log_info("Samples with complete covariates: {sum(complete_samples)} out of {length(complete_samples)}")
  log_info("Using {length(pc_cols)} proteomic PCs for adjustment")

  return(list(
    covariates = covariate_data,
    complete_samples = complete_samples,
    young_sample_ids = young_sample_ids
  ))
}

# Function to adjust for covariates using linear regression
adjust_covariates_lm <- function(npx_matrix, covariates, adjust_for = c("age", "sex", "pPC1-10")) {
  log_info("Adjusting for covariates using linear regression")

  # Prepare covariate matrix
  if("age" %in% adjust_for) {
    cov_matrix <- cbind(age = covariates$age)
  }

  if("sex" %in% adjust_for) {
    cov_matrix <- cbind(cov_matrix, sex = covariates$sex)
  }

  if("bmi" %in% adjust_for && !all(is.na(covariates$bmi))) {
    cov_matrix <- cbind(cov_matrix, bmi = covariates$bmi)
  }

  if("smoking" %in% adjust_for && !all(is.na(covariates$smoking))) {
    cov_matrix <- cbind(cov_matrix, smoking = covariates$smoking)
  }

  # Use proteomic PCs (pPC) instead of genomic PCs
  if("pPC1-10" %in% adjust_for) {
    pc_cols <- paste0("pPC", 1:10)
    pc_cols <- pc_cols[pc_cols %in% names(covariates)]
    if (length(pc_cols) > 0) {
    cov_matrix <- cbind(cov_matrix, as.matrix(covariates[, ..pc_cols]))
    }
  } else if("pPC1-5" %in% adjust_for) {
    pc_cols <- paste0("pPC", 1:5)
    pc_cols <- pc_cols[pc_cols %in% names(covariates)]
    if (length(pc_cols) > 0) {
    cov_matrix <- cbind(cov_matrix, as.matrix(covariates[, ..pc_cols]))
    }
  }

  # Initialize adjusted matrix
  adjusted_matrix <- npx_matrix

  # Adjust each protein
  n_proteins <- ncol(npx_matrix)

  for(i in 1:n_proteins) {
    if(i %% 100 == 0) {
      log_info("Processing protein {i}/{n_proteins}")
    }

    # Get protein expression
    y <- npx_matrix[, i]

    # Only adjust for samples with complete data
    complete_idx <- complete.cases(cbind(y, cov_matrix))

    if(sum(complete_idx) < 50) {
      log_debug("Skipping protein {i} due to insufficient complete cases")
      next
    }

    # Fit linear model
    fit_data <- data.frame(y = y[complete_idx], cov_matrix[complete_idx, , drop = FALSE])

    tryCatch({
      model <- lm(y ~ ., data = fit_data)

      # Get residuals
      residuals <- residuals(model)

      # Add back the mean
      adjusted_values <- residuals + mean(y[complete_idx], na.rm = TRUE)

      # Update adjusted matrix
      adjusted_matrix[complete_idx, i] <- adjusted_values

    }, error = function(e) {
      log_debug("Failed to adjust protein {i}: {e$message}")
    })
  }

  log_info("Covariate adjustment complete")

  return(adjusted_matrix)
}

# Function to evaluate covariate effects
evaluate_covariate_effects <- function(npx_matrix, covariates, use_all_proteins = TRUE) {
  log_info("Evaluating covariate effects on protein expression")

  # Use ALL proteins for comprehensive evaluation
  if (use_all_proteins) {
    sample_proteins <- seq_len(ncol(npx_matrix))
  } else {
  set.seed(123)
  sample_proteins <- sample(ncol(npx_matrix), min(100, ncol(npx_matrix)))
  }

  covariate_effects <- list()

  # Test age effect
  age_correlations <- sapply(sample_proteins, function(i) {
    cor(npx_matrix[, i], covariates$age, use = "complete.obs")
  })

  covariate_effects$age <- data.table(
    protein_idx = sample_proteins,
    correlation = age_correlations,
    abs_correlation = abs(age_correlations)
  )

  # Test sex effect
  sex_effects <- sapply(sample_proteins, function(i) {
    complete_idx <- !is.na(npx_matrix[, i]) & !is.na(covariates$sex)
    if(sum(complete_idx) < 10) return(NA)

    t.test(npx_matrix[complete_idx, i] ~ covariates$sex[complete_idx])$statistic
  })

  covariate_effects$sex <- data.table(
    protein_idx = sample_proteins,
    t_statistic = sex_effects,
    abs_t_statistic = abs(sex_effects)
  )

  # Test pPC1 effect (proteomic PC1)
  ppc1_correlations <- sapply(sample_proteins, function(i) {
    if ("pPC1" %in% names(covariates)) {
      cor(npx_matrix[, i], covariates$pPC1, use = "complete.obs")
    } else {
      NA_real_
    }
  })

  covariate_effects$pPC1 <- data.table(
    protein_idx = sample_proteins,
    correlation = ppc1_correlations,
    abs_correlation = abs(ppc1_correlations)
  )

  # Summary
  summary_effects <- data.table(
    covariate = c("Age", "Sex", "pPC1"),
    mean_effect = c(
      mean(abs(age_correlations), na.rm = TRUE),
      mean(abs(sex_effects), na.rm = TRUE) / 10,  # Normalize t-statistic
      mean(abs(ppc1_correlations), na.rm = TRUE)
    ),
    max_effect = c(
      max(abs(age_correlations), na.rm = TRUE),
      max(abs(sex_effects), na.rm = TRUE) / 10,
      max(abs(ppc1_correlations), na.rm = TRUE)
    )
  )

  return(list(
    effects = covariate_effects,
    summary = summary_effects
  ))
}

# Enhanced plotting functions
create_enhanced_covariate_plots <- function(npx_matrix, adjusted_matrix, effects_before, effects_after, covariates,
                                            npx_matrix_age_filtered = NULL, adjusted_matrix_age_filtered = NULL,
                                            covariates_age_filtered = NULL) {
  log_info("Creating enhanced covariate adjustment visualizations")

  # Define colors
  col_before <- "#FC4E07"  # Before (red/orange)
  col_after <- "#00AFBB"   # After (blue/cyan)

  plots <- list()

  # Use age-filtered datasets for age-specific volcano plots if provided
  if (!is.null(npx_matrix_age_filtered) && !is.null(covariates_age_filtered)) {
    log_info("Using age-filtered datasets (age > 20) for age volcano plots to avoid pediatric sample leverage")
    npx_for_age_volcano <- npx_matrix_age_filtered
    adjusted_for_age_volcano <- adjusted_matrix_age_filtered
    covariates_for_age_volcano <- covariates_age_filtered
  } else {
    log_info("Using full datasets for age volcano plots (no age filtering)")
    npx_for_age_volcano <- npx_matrix
    adjusted_for_age_volcano <- adjusted_matrix
    covariates_for_age_volcano <- covariates
  }

  # ============================================================================
  # HIGH PRIORITY 1: Before/After Comparison Panels
  # ============================================================================

  # Age Before/After
  age_data <- data.table(
    before = effects_before$effects$age$correlation,
    after = effects_after$effects$age$correlation
  )
  age_long <- melt(age_data, measure.vars = c("before", "after"),
                   variable.name = "stage", value.name = "correlation")

  mean_before_age <- mean(abs(age_data$before), na.rm = TRUE)
  mean_after_age <- mean(abs(age_data$after), na.rm = TRUE)
  reduction_age <- round(100 * (mean_before_age - mean_after_age) / mean_before_age, 1)

  plots$age_comparison <- ggplot(age_long, aes(x = correlation, fill = stage)) +
    geom_histogram(aes(y = ..density..), bins = 50, alpha = 0.6, position = "identity") +
    geom_density(aes(color = stage), alpha = 0, linewidth = 1) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray30") +
    scale_fill_manual(values = c("before" = col_before, "after" = col_after),
                      labels = c("Before", "After")) +
    scale_color_manual(values = c("before" = col_before, "after" = col_after),
                       labels = c("Before", "After")) +
    labs(title = "Age Effect on Protein Expression: Before vs After Adjustment",
         subtitle = sprintf("Mean |r| reduction: %.3f -> %.3f (%s%% decrease)",
                           mean_before_age, mean_after_age, reduction_age),
         x = expression("Correlation with Age ("*italic(r)*")"),
         y = "Density",
         fill = "Stage", color = "Stage") +
    theme_bw() +
    theme(legend.position = "top")

  # Sex Before/After
  sex_data <- data.table(
    before = effects_before$effects$sex$t_statistic,
    after = effects_after$effects$sex$t_statistic
  )
  sex_long <- melt(sex_data[complete.cases(sex_data)], measure.vars = c("before", "after"),
                   variable.name = "stage", value.name = "t_statistic")

  mean_before_sex <- mean(abs(sex_data$before), na.rm = TRUE)
  mean_after_sex <- mean(abs(sex_data$after), na.rm = TRUE)
  reduction_sex <- round(100 * (mean_before_sex - mean_after_sex) / mean_before_sex, 1)

  plots$sex_comparison <- ggplot(sex_long, aes(x = t_statistic, fill = stage)) +
    geom_histogram(aes(y = ..density..), bins = 50, alpha = 0.6, position = "identity") +
    geom_density(aes(color = stage), alpha = 0, linewidth = 1) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray30") +
    scale_fill_manual(values = c("before" = col_before, "after" = col_after),
                      labels = c("Before", "After")) +
    scale_color_manual(values = c("before" = col_before, "after" = col_after),
                       labels = c("Before", "After")) +
    labs(title = "Sex Effect on Protein Expression: Before vs After Adjustment",
         subtitle = sprintf("Mean |t| reduction: %.3f -> %.3f (%s%% decrease)",
                           mean_before_sex, mean_after_sex, reduction_sex),
         x = expression("T-statistic ("*italic(t)*")"),
         y = "Density",
         fill = "Stage", color = "Stage") +
    theme_bw() +
    theme(legend.position = "top")

  # pPC1 Before/After
  ppc1_data <- data.table(
    before = effects_before$effects$pPC1$correlation,
    after = effects_after$effects$pPC1$correlation
  )
  ppc1_long <- melt(ppc1_data, measure.vars = c("before", "after"),
                    variable.name = "stage", value.name = "correlation")

  mean_before_ppc1 <- mean(abs(ppc1_data$before), na.rm = TRUE)
  mean_after_ppc1 <- mean(abs(ppc1_data$after), na.rm = TRUE)
  reduction_ppc1 <- round(100 * (mean_before_ppc1 - mean_after_ppc1) / mean_before_ppc1, 1)

  plots$ppc1_comparison <- ggplot(ppc1_long, aes(x = correlation, fill = stage)) +
    geom_histogram(aes(y = ..density..), bins = 50, alpha = 0.6, position = "identity") +
    geom_density(aes(color = stage), alpha = 0, linewidth = 1) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray30") +
    scale_fill_manual(values = c("before" = col_before, "after" = col_after),
                      labels = c("Before", "After")) +
    scale_color_manual(values = c("before" = col_before, "after" = col_after),
                       labels = c("Before", "After")) +
    labs(title = "Proteomic PC1 Effect: Before vs After Adjustment",
         subtitle = sprintf("Mean |r| reduction: %.3f -> %.3f (%s%% decrease)",
                           mean_before_ppc1, mean_after_ppc1, reduction_ppc1),
         x = expression("Correlation with pPC1 ("*italic(r)*")"),
         y = "Density",
         fill = "Stage", color = "Stage") +
    theme_bw() +
    theme(legend.position = "top")

  # pPC2 Before/After (requested by user)
  if ("pPC2" %in% names(covariates)) {
    ppc2_cor_before <- sapply(seq_len(ncol(npx_matrix)), function(i) {
      cor(npx_matrix[, i], covariates$pPC2, use = "complete.obs")
    })
    ppc2_cor_after <- sapply(seq_len(ncol(adjusted_matrix)), function(i) {
      cor(adjusted_matrix[, i], covariates$pPC2, use = "complete.obs")
    })

    ppc2_data <- data.table(before = ppc2_cor_before, after = ppc2_cor_after)
    ppc2_long <- melt(ppc2_data, measure.vars = c("before", "after"),
                      variable.name = "stage", value.name = "correlation")

    mean_before_ppc2 <- mean(abs(ppc2_data$before), na.rm = TRUE)
    mean_after_ppc2 <- mean(abs(ppc2_data$after), na.rm = TRUE)
    reduction_ppc2 <- round(100 * (mean_before_ppc2 - mean_after_ppc2) / mean_before_ppc2, 1)

    plots$ppc2_comparison <- ggplot(ppc2_long, aes(x = correlation, fill = stage)) +
      geom_histogram(aes(y = ..density..), bins = 50, alpha = 0.6, position = "identity") +
      geom_density(aes(color = stage), alpha = 0, linewidth = 1) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "gray30") +
      scale_fill_manual(values = c("before" = col_before, "after" = col_after),
                        labels = c("Before", "After")) +
      scale_color_manual(values = c("before" = col_before, "after" = col_after),
                         labels = c("Before", "After")) +
      labs(title = "Proteomic PC2 Effect: Before vs After Adjustment",
           subtitle = sprintf("Mean |r| reduction: %.3f -> %.3f (%s%% decrease)",
                             mean_before_ppc2, mean_after_ppc2, reduction_ppc2),
           x = expression("Correlation with pPC2 ("*italic(r)*")"),
           y = "Density",
           fill = "Stage", color = "Stage") +
      theme_bw() +
      theme(legend.position = "top")
  }

  # ============================================================================
  # HIGH PRIORITY 2: Effect Size Reduction Scatter Plot
  # ============================================================================

  log_info("Creating effect size reduction scatter plots")

  # Age effect reduction scatter
  age_scatter_data <- data.table(
    protein = colnames(npx_matrix),
    before = abs(effects_before$effects$age$correlation),
    after = abs(effects_after$effects$age$correlation)
  )
  age_scatter_data[, reduction := before - after]

  # Annotate top 20 proteins
  age_scatter_data[, label := ""]
  top_age <- age_scatter_data[order(-before)][1:min(20, .N)]
  age_scatter_data[protein %in% top_age$protein, label := protein]

  plots$age_scatter <- ggplot(age_scatter_data, aes(x = before, y = after)) +
    geom_point(alpha = 0.4, size = 1.5, color = "gray50") +
    geom_point(data = age_scatter_data[label != ""], aes(color = reduction), size = 2.5) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black") +
    geom_text_repel(data = age_scatter_data[label != ""], aes(label = label),
                    size = 2.5, max.overlaps = 15) +
    scale_color_gradient(low = col_after, high = col_before) +
    labs(title = "Age Effect Reduction Per Protein",
         subtitle = sprintf("All %s proteins; diagonal = no change", nrow(age_scatter_data)),
         x = expression("|"*italic(r)*"| Before Adjustment"),
         y = expression("|"*italic(r)*"| After Adjustment"),
         color = "Reduction") +
    theme_bw() +
    theme(legend.position = "right")

  # Sex effect reduction scatter
  sex_scatter_data <- data.table(
    protein = colnames(npx_matrix),
    before = abs(effects_before$effects$sex$t_statistic),
    after = abs(effects_after$effects$sex$t_statistic)
  )
  sex_scatter_data[, reduction := before - after]
  sex_scatter_data[, label := ""]
  top_sex <- sex_scatter_data[order(-before)][1:min(20, .N)]
  sex_scatter_data[protein %in% top_sex$protein, label := protein]

  plots$sex_scatter <- ggplot(sex_scatter_data[complete.cases(sex_scatter_data)], aes(x = before, y = after)) +
    geom_point(alpha = 0.4, size = 1.5, color = "gray50") +
    geom_point(data = sex_scatter_data[label != ""], aes(color = reduction), size = 2.5) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black") +
    geom_text_repel(data = sex_scatter_data[label != ""], aes(label = label),
                    size = 2.5, max.overlaps = 15) +
    scale_color_gradient(low = col_after, high = col_before) +
    labs(title = "Sex Effect Reduction Per Protein",
         subtitle = sprintf("All %s proteins; diagonal = no change", nrow(sex_scatter_data)),
         x = expression("|"*italic(t)*"| Before Adjustment"),
         y = expression("|"*italic(t)*"| After Adjustment"),
         color = "Reduction") +
    theme_bw() +
    theme(legend.position = "right")

  # pPC1 effect reduction scatter
  ppc1_scatter_data <- data.table(
    protein = colnames(npx_matrix),
    before = abs(effects_before$effects$pPC1$correlation),
    after = abs(effects_after$effects$pPC1$correlation)
  )
  ppc1_scatter_data[, reduction := before - after]
  ppc1_scatter_data[, label := ""]
  top_ppc1 <- ppc1_scatter_data[order(-before)][1:min(20, .N)]
  ppc1_scatter_data[protein %in% top_ppc1$protein, label := protein]

  plots$ppc1_scatter <- ggplot(ppc1_scatter_data, aes(x = before, y = after)) +
    geom_point(alpha = 0.4, size = 1.5, color = "gray50") +
    geom_point(data = ppc1_scatter_data[label != ""], aes(color = reduction), size = 2.5) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black") +
    geom_text_repel(data = ppc1_scatter_data[label != ""], aes(label = label),
                    size = 2.5, max.overlaps = 15) +
    scale_color_gradient(low = col_after, high = col_before) +
    labs(title = "Proteomic PC1 Effect Reduction Per Protein",
         subtitle = sprintf("All %s proteins; diagonal = no change", nrow(ppc1_scatter_data)),
         x = expression("|"*italic(r)*"| Before Adjustment"),
         y = expression("|"*italic(r)*"| After Adjustment"),
         color = "Reduction") +
    theme_bw() +
    theme(legend.position = "right")

  # ============================================================================
  # HIGH PRIORITY 3: Covariate Importance Comparison
  # ============================================================================

  log_info("Creating covariate importance comparison")

  importance_data <- rbind(
    data.table(covariate = "Age", before = mean_before_age, after = mean_after_age),
    data.table(covariate = "Sex", before = mean_before_sex/10, after = mean_after_sex/10),
    data.table(covariate = "pPC1", before = mean_before_ppc1, after = mean_after_ppc1)
  )
  importance_data[, reduction_pct := round(100 * (before - after) / before, 1)]

  importance_long <- melt(importance_data, id.vars = c("covariate", "reduction_pct"),
                          measure.vars = c("before", "after"),
                          variable.name = "stage", value.name = "effect")

  plots$importance_comparison <- ggplot(importance_long, aes(x = covariate, y = effect, fill = stage)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7, alpha = 0.8) +
    geom_text(data = importance_data, aes(x = covariate, y = before, label = paste0(reduction_pct, "% decrease")),
              position = position_nudge(y = 0.02), size = 3.5, fontface = "bold", inherit.aes = FALSE) +
    scale_fill_manual(values = c("before" = col_before, "after" = col_after),
                      labels = c("Before", "After")) +
    labs(title = "Covariate Importance: Before vs After Adjustment",
         subtitle = "Mean absolute effect sizes across all proteins",
         x = "Covariate", y = "Mean |Effect|",
         fill = "Stage") +
    theme_bw() +
    theme(legend.position = "top")

  # ============================================================================
  # VOLCANO PLOT: Age-Associated Proteins
  # ============================================================================

  log_info("Creating age-associated volcano plot - using age-filtered data (age > 20)")

  # Calculate age correlations on age-filtered dataset
  age_corr_filtered <- sapply(seq_len(ncol(npx_for_age_volcano)), function(i) {
    cor.test(npx_for_age_volcano[, i], covariates_for_age_volcano$age, use = "complete.obs")$estimate
  })

  age_pval_filtered <- sapply(seq_len(ncol(npx_for_age_volcano)), function(i) {
    cor.test(npx_for_age_volcano[, i], covariates_for_age_volcano$age, use = "complete.obs")$p.value
  })

  # Prepare volcano data (use age-filtered correlations and p-values)
  volcano_age_data <- data.table(
    protein = colnames(npx_for_age_volcano),
    correlation = age_corr_filtered,
    p_value = age_pval_filtered
  )
  volcano_age_data[, neg_log10_p := -log10(p_value + 1e-300)]

  # Annotate top 10 positive and top 10 negative
  volcano_age_data[, label := ""]
  volcano_age_data[, direction := "non-significant"]
  top_positive <- volcano_age_data[order(-correlation)][1:10]
  top_negative <- volcano_age_data[order(correlation)][1:10]
  volcano_age_data[protein %in% top_positive$protein, `:=`(label = protein, direction = "positive")]
  volcano_age_data[protein %in% top_negative$protein, `:=`(label = protein, direction = "negative")]

  # Color by absolute effect size (like other plots)
  volcano_age_data[, abs_corr := abs(correlation)]

  plots$age_volcano <- ggplot(volcano_age_data, aes(x = correlation, y = neg_log10_p)) +
    geom_point(aes(color = abs_corr), alpha = 0.5, size = 2) +
    geom_point(data = volcano_age_data[label != ""], aes(color = abs_corr), size = 3.5, alpha = 0.8) +
    geom_text_repel(data = volcano_age_data[label != ""], aes(label = label),
                    size = 2.8, max.overlaps = 20, box.padding = 0.5) +
    scale_color_gradient(low = col_after, high = col_before,
                         name = expression("|"*italic(r)*"|")) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray30") +
    labs(title = "Age-Associated Proteins (Before Adjustment)",
         subtitle = "Samples with age <= 20 years excluded | Top 10 positive + Top 10 negative labeled",
         x = expression("Correlation with Age ("*italic(r)*")"),
         y = "-log10(p-value)") +
    theme_bw() +
    theme(legend.position = "right")

  # ============================================================================
  # VOLCANO PLOT: Age-Associated Proteins (Multicollinearity-Aware)
  # ============================================================================

  # Get configuration parameters (passed from main function)
  use_multicollinearity <- getOption("age_assoc_use_multicollinearity", TRUE)
  cor_threshold <- getOption("age_assoc_cor_threshold", 0.80)
  sd_filter_threshold <- getOption("age_assoc_sd_filter", 3)

  if (use_multicollinearity) {
    log_info("========== AGE-ASSOCIATED PROTEINS: MULTICOLLINEARITY-AWARE MODE ==========")
    log_info("Configuration:")
    log_info("  - Correlation threshold: {cor_threshold}")
    log_info("  - SD filter threshold: ±{sd_filter_threshold}SD")

    log_info("Step 1: Identifying correlated protein groups...")

    # Step 1: Identify correlated protein groups
    protein_groups <- identify_protein_groups(npx_for_age_volcano, cor_threshold = cor_threshold)

    log_info("Step 2: Selecting best protein from each group...")

    # Step 2: Select best protein from each group (multicollinearity resolution)
    # Convert covariates_for_age_volcano to data.frame for compatibility
    covariates_df <- as.data.frame(covariates_for_age_volcano)

    age_adj_results_all <- select_best_proteins_per_group(
      npx_for_age_volcano,
      covariates_df,
      protein_groups
    )

    log_info("Step 3: Creating volcano plot with confidence intervals...")

    # Step 3: Create volcano plot (selected proteins only, with CI, outlier filtering)
    plots$age_volcano_adjusted <- create_volcano_with_ci(
      age_adj_results_all,
      selected_only = TRUE,
      sd_threshold = sd_filter_threshold
    )

    # Step 4: Save comprehensive results tables
    log_info("Step 4: Saving comprehensive results tables...")
  } else {
    log_info("========== AGE-ASSOCIATED PROTEINS: STANDARD MODE (NO MULTICOLLINEARITY ADJUSTMENT) ==========")
    log_info("Configuration: use_multicollinearity_aware = FALSE")
    log_info("Will use standard regression for all proteins (no grouping/selection)")

    # Standard approach: fit all proteins individually
    # This is a fallback - not implementing now since user requested multicollinearity-aware approach
    stop("Standard mode not implemented. Please set use_multicollinearity_aware = TRUE in config.")
  }

  # Add Bonferroni significance flag to full results
  bonferroni_threshold_full <- 0.05 / sum(!is.na(age_adj_results_all$p_value))
  age_adj_results_all[, bonferroni_sig := p_value < bonferroni_threshold_full]

  # Full results table (all proteins)
  full_results <- age_adj_results_all[order(p_value)]
  setcolorder(full_results, c("protein", "beta", "se", "t_stat", "p_value",
                              "r_squared", "n_samples", "bonferroni_sig",
                              "group_id", "group_size", "selected"))

  # Filtered table (p < 0.05 only)
  sig_results <- full_results[p_value < 0.05]

  # Mark Bonferroni significant proteins with asterisk
  full_results[, bonferroni_flag := ifelse(bonferroni_sig, "***", "")]
  sig_results[, bonferroni_flag := ifelse(bonferroni_sig, "***", "")]

  # Save tables
  fwrite(full_results,
         ,
         sep = "\t", quote = FALSE)

  fwrite(sig_results,
         ,
         sep = "\t", quote = FALSE)

  # Selected proteins only (for volcano plot)
  selected_results <- full_results[selected == TRUE]
  fwrite(selected_results,
         ,
         sep = "\t", quote = FALSE)

  log_info("Saved results tables:")
  log_info("  - Full results (all proteins): {nrow(full_results)} proteins")
  log_info("  - Significant (p<0.05): {nrow(sig_results)} proteins")
  log_info("  - Selected (multicollinearity-resolved): {nrow(selected_results)} proteins")
  log_info("  - Bonferroni significant: {sum(full_results$bonferroni_sig, na.rm=TRUE)} proteins")

  # ============================================================================
  # MEDIUM PRIORITY: Multi-panel All Proteomic PCs (pPC1-10)
  # ============================================================================

  log_info("Creating multi-panel proteomic PC effects")

  # Collect all pPC correlations before/after
  ppc_cols <- grep("^pPC[0-9]+$", names(covariates), value = TRUE)
  if (length(ppc_cols) > 0) {
    ppc_effects_list <- list()

    for (ppc in ppc_cols) {
      # Before
      cor_before <- sapply(seq_len(ncol(npx_matrix)), function(i) {
        cor(npx_matrix[, i], covariates[[ppc]], use = "complete.obs")
      })
      # After
      cor_after <- sapply(seq_len(ncol(adjusted_matrix)), function(i) {
        cor(adjusted_matrix[, i], covariates[[ppc]], use = "complete.obs")
      })

      ppc_effects_list[[ppc]] <- data.table(
        PC = ppc,
        stage = rep(c("Before", "After"), each = length(cor_before)),
        correlation = c(cor_before, cor_after)
      )
    }

    ppc_all_data <- rbindlist(ppc_effects_list)
    ppc_all_data[, PC := factor(PC, levels = ppc_cols)]
    ppc_all_data[, stage := factor(stage, levels = c("Before", "After"))]

    plots$all_ppcs_facet <- ggplot(ppc_all_data, aes(x = correlation, fill = stage)) +
      geom_histogram(aes(y = ..density..), bins = 30, alpha = 0.6, position = "identity") +
      geom_density(aes(color = stage), alpha = 0, linewidth = 0.7) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "gray30", linewidth = 0.3) +
      scale_fill_manual(values = c("Before" = col_before, "After" = col_after)) +
      scale_color_manual(values = c("Before" = col_before, "After" = col_after)) +
      facet_wrap(~ PC, ncol = 5, scales = "free_y") +
      labs(title = "All Proteomic PCs (pPC1-10): Effect Before vs After Adjustment",
           subtitle = "Compression toward zero validates technical variation removal",
           x = "Correlation with Proteomic PC", y = "Density",
           fill = "Stage", color = "Stage") +
      theme_bw() +
      theme(legend.position = "top", strip.text = element_text(size = 9))
  }

  # ============================================================================
  # MEDIUM PRIORITY: Variance Decomposition - Paired Before/After Bars
  # ============================================================================

  log_info("Creating variance decomposition analysis with before/after comparison")

  # Get Acadia palette from paletteer (4 colors)
  acadia_colors <- paletteer::paletteer_d("nord::aurora", n = 4)
  names(acadia_colors) <- c("pPC1", "pPC2", "Sex", "Age")

  # Compute R² for 200, 500, and 1000 proteins (BEFORE and AFTER)
  protein_sample_sizes <- c(200, 500, 1000)
  var_decomp_list <- list()

  for (n_prots in protein_sample_sizes) {
    set.seed(456)
    sample_prots <- sample(seq_len(ncol(npx_matrix)), min(n_prots, ncol(npx_matrix)))

    # BEFORE adjustment
    var_decomp_before <- rbindlist(lapply(sample_prots, function(i) {
      y_before <- npx_matrix[, i]

      fit_age <- lm(y_before ~ age, data = data.frame(y_before = y_before, age = covariates$age))
      fit_sex <- lm(y_before ~ sex, data = data.frame(y_before = y_before, sex = covariates$sex))
      fit_ppc1 <- lm(y_before ~ ppc1, data = data.frame(y_before = y_before, ppc1 = covariates$pPC1))
      fit_ppc2 <- lm(y_before ~ ppc2, data = data.frame(y_before = y_before, ppc2 = covariates$pPC2))

      data.table(
        protein_idx = i,
        r2_ppc1 = tryCatch(summary(fit_ppc1)$r.squared, error = function(e) 0),
        r2_ppc2 = tryCatch(summary(fit_ppc2)$r.squared, error = function(e) 0),
        r2_sex = tryCatch(summary(fit_sex)$r.squared, error = function(e) 0),
        r2_age = tryCatch(summary(fit_age)$r.squared, error = function(e) 0)
      )
    }))

    # AFTER adjustment
    var_decomp_after <- rbindlist(lapply(sample_prots, function(i) {
      y_after <- adjusted_matrix[, i]

      fit_age <- lm(y_after ~ age, data = data.frame(y_after = y_after, age = covariates$age))
      fit_sex <- lm(y_after ~ sex, data = data.frame(y_after = y_after, sex = covariates$sex))
      fit_ppc1 <- lm(y_after ~ ppc1, data = data.frame(y_after = y_after, ppc1 = covariates$pPC1))
      fit_ppc2 <- lm(y_after ~ ppc2, data = data.frame(y_after = y_after, ppc2 = covariates$pPC2))

      data.table(
        protein_idx = i,
        r2_ppc1 = tryCatch(summary(fit_ppc1)$r.squared, error = function(e) 0),
        r2_ppc2 = tryCatch(summary(fit_ppc2)$r.squared, error = function(e) 0),
        r2_sex = tryCatch(summary(fit_sex)$r.squared, error = function(e) 0),
        r2_age = tryCatch(summary(fit_age)$r.squared, error = function(e) 0)
      )
    }))

    # Calculate mean R² for each covariate, before and after
    var_summary_before <- data.table(
      covariate = c("pPC1", "pPC2", "Sex", "Age"),
      mean_r2 = c(
        mean(var_decomp_before$r2_ppc1, na.rm = TRUE),
        mean(var_decomp_before$r2_ppc2, na.rm = TRUE),
        mean(var_decomp_before$r2_sex, na.rm = TRUE),
        mean(var_decomp_before$r2_age, na.rm = TRUE)
      ),
      stage = "Before",
      n_proteins = n_prots
    )

    var_summary_after <- data.table(
      covariate = c("pPC1", "pPC2", "Sex", "Age"),
      mean_r2 = c(
        mean(var_decomp_after$r2_ppc1, na.rm = TRUE),
        mean(var_decomp_after$r2_ppc2, na.rm = TRUE),
        mean(var_decomp_after$r2_sex, na.rm = TRUE),
        mean(var_decomp_after$r2_age, na.rm = TRUE)
      ),
      stage = "After",
      n_proteins = n_prots
    )

    var_decomp_list[[paste0(n_prots, "_before")]] <- var_summary_before
    var_decomp_list[[paste0(n_prots, "_after")]] <- var_summary_after
  }

  var_decomp_combined <- rbindlist(var_decomp_list)
  var_decomp_combined[, covariate := factor(covariate, levels = c("pPC1", "pPC2", "Sex", "Age"))]
  var_decomp_combined[, stage := factor(stage, levels = c("Before", "After"))]
  var_decomp_combined[, n_proteins := factor(n_proteins, levels = protein_sample_sizes)]

  plots$variance_decomposition <- ggplot(var_decomp_combined, aes(x = covariate, y = mean_r2, fill = covariate, alpha = stage)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.75) +
    geom_text(aes(label = sprintf("%.1f%%", mean_r2 * 100)),
              position = position_dodge(width = 0.8), vjust = -0.5, size = 2.5) +
    scale_fill_manual(values = acadia_colors) +
    scale_alpha_manual(values = c("Before" = 0.5, "After" = 1.0)) +
    facet_wrap(~ n_proteins, ncol = 3, labeller = labeller(n_proteins = function(x) paste0("n = ", x))) +
    labs(title = "Variance Explained by Each Covariate: Before vs After Adjustment",
         subtitle = expression(R^2*" from individual linear models - 200, 500, and 1,000 proteins"),
         x = "Covariate",
         y = expression("Mean "*R^2*" (Variance Explained)"),
         fill = "Covariate",
         alpha = "Stage") +
    theme_bw() +
    theme(legend.position = "right",
          strip.background = element_rect(fill = "gray95"),
          strip.text = element_text(face = "bold"))

  return(plots)
}

# Function to create proteomic PC contour density plots using filled.contour
create_pca_biplot_comparison <- function(npx_matrix, adjusted_matrix, covariates) {
  log_info("Creating PCA contour density plots with filled.contour (before/after)")

  # Remove constant columns and handle missing values
  sds_before <- apply(npx_matrix, 2, sd, na.rm = TRUE)
  keep_cols <- which(sds_before > 0 & !is.na(sds_before))
  npx_matrix_clean <- npx_matrix[, keep_cols]
  adjusted_matrix_clean <- adjusted_matrix[, keep_cols]

  # Remove rows with any missing values for PCA
  complete_rows <- complete.cases(npx_matrix_clean)
  npx_matrix_clean <- npx_matrix_clean[complete_rows, ]
  adjusted_matrix_clean <- adjusted_matrix_clean[complete_rows, ]

  log_info("PCA: using {nrow(npx_matrix_clean)} complete samples and {ncol(npx_matrix_clean)} proteins")

  # Perform PCA on protein matrices (before and after)
  pca_before <- prcomp(npx_matrix_clean, center = TRUE, scale. = TRUE)
  pc1_before <- pca_before$x[, 1]
  pc2_before <- pca_before$x[, 2]
  var_exp_before <- round(100 * summary(pca_before)$importance[2, 1:2], 1)

  pca_after <- prcomp(adjusted_matrix_clean, center = TRUE, scale. = TRUE)
  pc1_after <- pca_after$x[, 1]
  pc2_after <- pca_after$x[, 2]
  var_exp_after <- round(100 * summary(pca_after)$importance[2, 1:2], 1)

  # Compute 2D kernel density for filled.contour
  # BEFORE - use full data range
  kde_before <- MASS::kde2d(pc1_before, pc2_before, n = 100)

  # AFTER - use quantile-based limits to match earlier implementation
  # After adjustment, samples cluster more tightly, so we need tighter limits
  # Use 2.5th to 97.5th percentile to capture 95% of data without extreme outliers
  # This matches the earlier implementation's reasonable appearance (PC1: -3 to 2, PC2: -1 to 4)
  pc1_after_q025 <- quantile(pc1_after, 0.025, na.rm = TRUE)
  pc1_after_q975 <- quantile(pc1_after, 0.975, na.rm = TRUE)
  pc2_after_q025 <- quantile(pc2_after, 0.025, na.rm = TRUE)
  pc2_after_q975 <- quantile(pc2_after, 0.975, na.rm = TRUE)

  # Set limits to 2.5th-97.5th percentile (captures 95% of data, excludes extreme outliers)
  # Add small padding (5% of range) for visual clarity
  pc1_range <- pc1_after_q975 - pc1_after_q025
  pc2_range <- pc2_after_q975 - pc2_after_q025
  xlim_after <- c(pc1_after_q025 - 0.05*pc1_range, pc1_after_q975 + 0.05*pc1_range)
  ylim_after <- c(pc2_after_q025 - 0.05*pc2_range, pc2_after_q975 + 0.05*pc2_range)

  log_info("After adjustment plot limits (quantile-based):")
  log_info("  PC1 range: [{round(xlim_after[1], 2)}, {round(xlim_after[2], 2)}]")
  log_info("  PC2 range: [{round(ylim_after[1], 2)}, {round(ylim_after[2], 2)}]")

  # Compute KDE with explicit limits to ensure plot uses these ranges
  kde_after <- MASS::kde2d(pc1_after, pc2_after, n = 100,
                          lims = c(xlim_after, ylim_after))

  # Return data for plotting (we'll create the actual plots in a separate function)
  return(list(
    kde_before = kde_before,
    kde_after = kde_after,
    pc1_before = pc1_before,
    pc2_before = pc2_before,
    pc1_after = pc1_after,
    pc2_after = pc2_after,
    var_exp_before = var_exp_before,
    var_exp_after = var_exp_after,
    pca_before = pca_before,
    pca_after = pca_after,
    xlim_after = xlim_after,
    ylim_after = ylim_after
  ))
}

# Main execution
main <- function() {

  # Load data from previous steps
  log_info("Loading data from previous steps")

  # Try loading from step 07 (enhanced bridge) first, fallback to step 06 (median normalization)
  npx_file_candidates <- c(
    ,
    ,
    
  )

  npx_file <- NULL
  for (candidate in npx_file_candidates) {
    if (file.exists(candidate)) {
      npx_file <- candidate
      break
    }
  }

  if (is.null(npx_file)) {
    stop("No normalized NPX matrix found. Please run step 06 (and optionally step 07) first.")
  }

  log_info("Using normalized matrix: {npx_file}")
  npx_matrix <- readRDS(npx_file)
  metadata <- readRDS()

  # Get sample IDs
  sample_ids <- rownames(npx_matrix)

  # Load proteomic PCs from step 01
  log_info("Loading proteomic PCs from PCA analysis")
  proteomic_pcs <- load_proteomic_pcs(
    ,
    n_pcs = 10
  )

  # Prepare covariates (including proteomic PCs)
  covariate_result <- prepare_covariates(
    config$covariates$covariate_with_smoking,
    sample_ids,
    metadata,
    proteomic_pcs
  )

  # Create filtered datasets for age-specific analyses (exclude age <= 20)
  young_sample_ids <- covariate_result$young_sample_ids
  if (length(young_sample_ids) > 0) {
    log_info("Creating age-filtered datasets: excluding {length(young_sample_ids)} samples with age <= 20 years")
    age_filtered_samples <- !(sample_ids %in% young_sample_ids)
    npx_matrix_age_filtered <- npx_matrix[age_filtered_samples, ]
    covariates_age_filtered <- covariate_result$covariates[age_filtered_samples, ]
    adjusted_matrix_age_filtered <- NULL  # Will be set after adjustment
  } else {
    log_info("No samples with age <= 20 found")
    age_filtered_samples <- rep(TRUE, length(sample_ids))
    npx_matrix_age_filtered <- npx_matrix
    covariates_age_filtered <- covariate_result$covariates
    adjusted_matrix_age_filtered <- NULL
  }

  # Evaluate covariate effects before adjustment (ALL proteins)
  effects_before <- evaluate_covariate_effects(npx_matrix, covariate_result$covariates, use_all_proteins = TRUE)

  # Adjust for covariates (using proteomic PCs)
  adjusted_matrix <- adjust_covariates_lm(
    npx_matrix,
    covariate_result$covariates,
    adjust_for = c("age", "sex", "bmi", "smoking", "pPC1-10")
  )

  # Create age-filtered adjusted matrix (exclude samples with age <= 20)
  if (length(young_sample_ids) > 0) {
    adjusted_matrix_age_filtered <- adjusted_matrix[age_filtered_samples, ]
    log_info("Age-filtered datasets created: {nrow(npx_matrix_age_filtered)} samples (excluded {length(young_sample_ids)} with age <= 20)")
  } else {
    adjusted_matrix_age_filtered <- adjusted_matrix
  }

  # Evaluate covariate effects after adjustment (ALL proteins)
  effects_after <- evaluate_covariate_effects(adjusted_matrix, covariate_result$covariates, use_all_proteins = TRUE)

  # Set age association analysis options from config
  use_multicollinearity <- config$parameters$age_association$use_multicollinearity_aware
  cor_threshold <- config$parameters$age_association$protein_correlation_threshold
  sd_filter_threshold <- config$parameters$age_association$sd_filter_threshold

  log_info("Age association configuration:")
  log_info("  - Use multicollinearity-aware approach: {use_multicollinearity}")
  log_info("  - Protein correlation threshold: {cor_threshold}")
  log_info("  - SD filter threshold: ±{sd_filter_threshold}")

  # Set as global options for plotting functions
  options(age_assoc_use_multicollinearity = use_multicollinearity)
  options(age_assoc_cor_threshold = cor_threshold)
  options(age_assoc_sd_filter = sd_filter_threshold)

  # Create enhanced plots (use age-filtered datasets for age-specific analyses)
  log_info("Generating enhanced visualization plots")
  covariate_plots <- create_enhanced_covariate_plots(
    npx_matrix, adjusted_matrix,
    effects_before, effects_after,
    covariate_result$covariates,
    npx_matrix_age_filtered, adjusted_matrix_age_filtered,
    covariates_age_filtered
  )

  # Create PCA biplot comparison
  pca_comparison <- create_pca_biplot_comparison(npx_matrix, adjusted_matrix,
                                                  covariate_result$covariates)

  # Save outputs
  log_info("Saving covariate adjustment results")

  saveRDS(adjusted_matrix, )
  saveRDS(covariate_result, )
  saveRDS(effects_before, )
  saveRDS(effects_after, )

  # Save tables
  fwrite(effects_before$summary, , sep = "\t")
  fwrite(effects_after$summary, , sep = "\t")

  # Save HIGH PRIORITY plots (Before/After Comparisons)
  log_info("Saving HIGH PRIORITY plots: Before/After comparisons")
  ggsave(,
         covariate_plots$age_comparison, width = 10, height = 6)
  ggsave(,
         covariate_plots$sex_comparison, width = 10, height = 6)
  ggsave(,
         covariate_plots$ppc1_comparison, width = 10, height = 6)

  if (!is.null(covariate_plots$ppc2_comparison)) {
    ggsave(,
           covariate_plots$ppc2_comparison, width = 10, height = 6)
  }

  # Save HIGH PRIORITY plots (Effect Reduction Scatter)
  log_info("Saving HIGH PRIORITY plots: Effect reduction scatter plots")
  ggsave(,
         covariate_plots$age_scatter, width = 10, height = 8)
  ggsave(,
         covariate_plots$sex_scatter, width = 10, height = 8)
  ggsave(,
         covariate_plots$ppc1_scatter, width = 10, height = 8)

  # Save HIGH PRIORITY plots (Covariate Importance)
  log_info("Saving HIGH PRIORITY plots: Covariate importance comparison")
  ggsave(,
         covariate_plots$importance_comparison, width = 10, height = 6)

  # Save Volcano plots
  log_info("Saving age-associated volcano plots")
  ggsave(,
         covariate_plots$age_volcano, width = 10, height = 8)
  ggsave(,
         covariate_plots$age_volcano_adjusted, width = 10, height = 8)

  # Save PCA biplot comparison using filled.contour (base R graphics)
  log_info("Saving PCA contour density plots with filled.contour")
  pdf(,
      width = 14, height = 7)

  par(mfrow = c(1, 2), mar = c(4, 4, 3, 6))

  # BEFORE plot
  filled.contour(pca_comparison$kde_before,
                 color.palette = colorRampPalette(c("lightyellow", "orange", "darkorange", "red", "darkred")),
                 main = sprintf("Before Adjustment\nPC1: %.1f%%, PC2: %.1f%%",
                               pca_comparison$var_exp_before[1], pca_comparison$var_exp_before[2]),
                 xlab = paste0("PC1 (", pca_comparison$var_exp_before[1], "%)"),
                 ylab = paste0("PC2 (", pca_comparison$var_exp_before[2], "%)"),
                 plot.axes = {
                   axis(1)
                   axis(2)
                   contour(pca_comparison$kde_before, add = TRUE, lwd = 2, drawlabels = FALSE)
                   points(pca_comparison$pc1_before, pca_comparison$pc2_before,
                          pch = 16, cex = 0.3, col = rgb(0, 0, 0, 0.05))
                 })

  # AFTER plot - use quantile-based limits to match earlier implementation
  # The KDE was computed with explicit lims, so the grid matches the desired range
  # This ensures the plot shows a comprehensible view without the dense blob appearance
  filled.contour(pca_comparison$kde_after,
                 color.palette = colorRampPalette(c("lightyellow", "orange", "darkorange", "red", "darkred")),
                 main = sprintf("After Adjustment\nPC1: %.1f%%, PC2: %.1f%%",
                               pca_comparison$var_exp_after[1], pca_comparison$var_exp_after[2]),
                 xlab = paste0("PC1 (", pca_comparison$var_exp_after[1], "%)"),
                 ylab = paste0("PC2 (", pca_comparison$var_exp_after[2], "%)"),
                 plot.axes = {
                   axis(1, at = pretty(pca_comparison$xlim_after, n = 5))
                   axis(2, at = pretty(pca_comparison$ylim_after, n = 5))
                   contour(pca_comparison$kde_after, add = TRUE, lwd = 2, drawlabels = FALSE)
                   # Only plot points within the zoomed range
                   in_range <- pca_comparison$pc1_after >= pca_comparison$xlim_after[1] &
                              pca_comparison$pc1_after <= pca_comparison$xlim_after[2] &
                              pca_comparison$pc2_after >= pca_comparison$ylim_after[1] &
                              pca_comparison$pc2_after <= pca_comparison$ylim_after[2]
                   points(pca_comparison$pc1_after[in_range], pca_comparison$pc2_after[in_range],
                          pch = 16, cex = 0.3, col = rgb(0, 0, 0, 0.05))
                 })

  dev.off()

  # Save MEDIUM PRIORITY plots
  log_info("Saving MEDIUM PRIORITY plots")
  if (!is.null(covariate_plots$all_ppcs_facet)) {
    ggsave(,
           covariate_plots$all_ppcs_facet, width = 16, height = 10)
  }

  ggsave(,
         covariate_plots$variance_decomposition, width = 8, height = 6)

  # Print summary
  cat("\n=== COVARIATE ADJUSTMENT SUMMARY ===\n")
  cat("Samples with complete covariates:", sum(covariate_result$complete_samples), "/",
      length(covariate_result$complete_samples), "\n")
  cat("\nCovariate effects before adjustment:\n")
  print(effects_before$summary)
  cat("\nCovariate effects after adjustment:\n")
  print(effects_after$summary)
  cat("\nAdjusted matrix dimensions:", nrow(adjusted_matrix), "x", ncol(adjusted_matrix), "\n")
  cat("Results saved to: ../output/normalized/\n")

  log_info("Covariate adjustment completed")

  return(list(
    adjusted_matrix = adjusted_matrix,
    covariates = covariate_result,
    effects = list(before = effects_before, after = effects_after)
  ))
}

# Run if executed directly
if (!interactive()) {
  result <- main()
}
