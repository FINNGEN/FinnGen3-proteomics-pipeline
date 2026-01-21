#!/usr/bin/env Rscript
# ==============================================================================
# 08_helper_multicollinearity_adjustment.R - Multicollinearity-Aware Age Association
# ==============================================================================
#
# Purpose:
#   Helper functions for multicollinearity-aware age association analysis. Addresses
#   extreme effect sizes caused by highly correlated proteins by identifying correlated
#   protein groups and selecting representative proteins using regression-based
#   selection (lowest p-value, highest R²). Used in age-association analysis to
#   prevent inflated effect sizes from multicollinearity.
#
# Author: Reza Jabal, PhD (rjabal@broadinstitute.org)
# Date: December 2025
# ==============================================================================

#' Identify correlated protein groups
#'
#' @param npx_matrix Matrix of protein expression values (samples × proteins)
#' @param cor_threshold Correlation threshold (default 0.80)
#' @return List with correlation matrix and group assignments
identify_protein_groups <- function(npx_matrix, cor_threshold = 0.80) {

  log_info("========== MULTICOLLINEARITY ANALYSIS: START ==========")
  log_info("Computing protein-protein correlation matrix for {ncol(npx_matrix)} proteins...")
  log_info("This may take several minutes...")

  # Compute correlation matrix
  start_time <- Sys.time()
  cor_matrix <- cor(npx_matrix, use = "pairwise.complete.obs")
  elapsed <- round(difftime(Sys.time(), start_time, units = "secs"), 1)
  log_info("Correlation matrix computed in {elapsed} seconds")

  # Set diagonal to 0 (don't consider self-correlation)
  diag(cor_matrix) <- 0

  # Identify groups using graph-based clustering
  # Proteins are connected if |correlation| > threshold
  log_info("Identifying correlated protein groups using threshold |r| > {cor_threshold}...")

  # Create adjacency matrix
  log_info("Creating adjacency matrix (proteins connected if |r| > {cor_threshold})...")
  adj_matrix <- abs(cor_matrix) > cor_threshold

  # Count connections
  n_connections <- sum(adj_matrix) / 2  # Divide by 2 because matrix is symmetric
  log_info("Found {n_connections} protein pairs with |r| > {cor_threshold}")

  # Assign groups using connected components algorithm
  log_info("Assigning proteins to correlation groups (connected components)...")
  n_proteins <- ncol(npx_matrix)
  groups <- rep(0, n_proteins)
  current_group <- 1

  for (i in 1:n_proteins) {
    if (i %% 500 == 0) {
      log_info("Progress: {i}/{n_proteins} proteins processed for grouping...")
    }

    if (groups[i] == 0) {
      # Start new group
      to_visit <- i
      while (length(to_visit) > 0) {
        current <- to_visit[1]
        to_visit <- to_visit[-1]

        if (groups[current] == 0) {
          groups[current] <- current_group
          # Add connected proteins
          connected <- which(adj_matrix[current, ] & groups == 0)
          to_visit <- c(to_visit, connected)
        }
      }
      current_group <- current_group + 1
    }
  }

  log_info("Grouping complete: Assigned {n_proteins} proteins to {current_group - 1} groups")

  # Create group data table
  group_dt <- data.table(
    protein = colnames(npx_matrix),
    group_id = groups
  )

  # Summary
  n_groups <- max(groups)
  group_sizes <- table(groups)
  n_singletons <- sum(group_sizes == 1)
  n_clustered <- n_groups - n_singletons

  log_info("Identified {n_groups} groups:")
  log_info("  - {n_singletons} singleton proteins (uncorrelated)")
  log_info("  - {n_clustered} groups with multiple proteins")
  log_info("  - Largest group: {max(group_sizes)} proteins")

  return(list(
    cor_matrix = cor_matrix,
    groups = group_dt,
    n_groups = n_groups
  ))
}


#' Select best protein from each correlated group
#'
#' @param npx_matrix NPX matrix (age-filtered)
#' @param covariates_df Covariates data frame
#' @param protein_groups Result from identify_protein_groups()
#' @return data.table with regression results for selected proteins
select_best_proteins_per_group <- function(npx_matrix, covariates_df, protein_groups) {

  log_info("========== PROTEIN SELECTION WITHIN GROUPS: START ==========")
  log_info("Selecting best protein from each correlated group...")

  groups_dt <- protein_groups$groups
  unique_groups <- unique(groups_dt$group_id)

  # Count group sizes
  group_sizes <- groups_dt[, .N, by = group_id]
  n_singletons <- sum(group_sizes$N == 1)
  n_multi <- sum(group_sizes$N > 1)
  max_size <- max(group_sizes$N)

  log_info("Group statistics:")
  log_info("  - Total groups: {length(unique_groups)}")
  log_info("  - Singleton groups (1 protein): {n_singletons}")
  log_info("  - Multi-protein groups: {n_multi}")
  log_info("  - Largest group size: {max_size} proteins")
  log_info("Estimating {n_singletons + sum(group_sizes[N > 1, N])} regressions needed...")

  # Results storage
  all_results <- list()
  start_time <- Sys.time()
  proteins_processed <- 0
  total_proteins <- nrow(groups_dt)

  for (gid in unique_groups) {
    # Get proteins in this group
    group_proteins <- groups_dt[group_id == gid, protein]

    if (length(group_proteins) == 1) {
      # Singleton - just fit regression
      prot_idx <- which(colnames(npx_matrix) == group_proteins[1])
      result <- fit_single_protein(npx_matrix[, prot_idx], covariates_df, group_proteins[1])
      result[, group_id := gid]
      result[, group_size := 1]
      result[, selected := TRUE]
      all_results[[length(all_results) + 1]] <- result
      proteins_processed <- proteins_processed + 1

    } else {
      # Multiple proteins - evaluate all and select best
      log_info("Processing group {gid}: {length(group_proteins)} correlated proteins (|r| > threshold)")
      log_info("  Correlated proteins in group: {paste(group_proteins, collapse=', ')}")
      log_info("  Fitting {length(group_proteins)} regression models to select best representative...")

      group_results <- list()

      for (prot_name in group_proteins) {
        prot_idx <- which(colnames(npx_matrix) == prot_name)
        result <- fit_single_protein(npx_matrix[, prot_idx], covariates_df, prot_name)
        result[, group_id := gid]
        result[, group_size := length(group_proteins)]
        group_results[[length(group_results) + 1]] <- result
        proteins_processed <- proteins_processed + 1
      }

      group_dt <- rbindlist(group_results)

      # Select best: lowest p-value, then highest R²
      group_dt <- group_dt[order(p_value, -r_squared)]
      best_protein <- group_dt[1, protein]

      # Mark selection
      group_dt[, selected := protein == best_protein]

      # Calculate excluded proteins (after 'selected' column is created)
      excluded_proteins <- setdiff(group_proteins, best_protein)

      log_info("  ✓ SELECTED: {best_protein} (p={format(group_dt[1, p_value], digits=3, scientific=TRUE)}, R²={format(group_dt[1, r_squared], digits=4)})")
      log_info("  ✗ EXCLUDED: {paste(excluded_proteins, collapse=', ')}")

      all_results[[length(all_results) + 1]] <- group_dt
    }

    if (gid %% 50 == 0) {
      elapsed <- round(difftime(Sys.time(), start_time, units = "secs"), 1)
      pct_complete <- round(100 * proteins_processed / total_proteins, 1)
      log_info("Progress: {gid}/{length(unique_groups)} groups | {proteins_processed}/{total_proteins} proteins ({pct_complete}%) | {elapsed}s elapsed")
    }
  }

  # Combine all results
  log_info("Combining results from all groups...")
  combined_results <- rbindlist(all_results)

  total_elapsed <- round(difftime(Sys.time(), start_time, units = "mins"), 1)

  # Statistics
  n_selected <- sum(combined_results$selected)
  n_total <- nrow(combined_results)
  n_excluded <- n_total - n_selected

  log_info("========== PROTEIN SELECTION: COMPLETE ==========")
  log_info("Selection completed in {total_elapsed} minutes")
  log_info("Final statistics:")
  log_info("  - Total proteins evaluated: {n_total}")
  log_info("  - Proteins selected (representatives): {n_selected}")
  log_info("  - Proteins excluded (multicollinearity): {n_excluded}")
  log_info("  - Exclusion rate: {round(100 * n_excluded / n_total, 1)}%")

  return(combined_results)
}


#' Fit regression for single protein
#'
#' @param protein_values Vector of protein expression values
#' @param covariates_df Data frame with covariates
#' @param protein_name Protein name
#' @return data.table with regression results
fit_single_protein <- function(protein_values, covariates_df, protein_name) {

  # Create model data
  model_data <- data.frame(
    age = covariates_df$age,
    protein = protein_values,
    bmi = covariates_df$bmi,
    sex = covariates_df$sex,
    pPC1 = covariates_df$pPC1,
    pPC2 = covariates_df$pPC2,
    pPC3 = covariates_df$pPC3,
    pPC4 = covariates_df$pPC4,
    pPC5 = covariates_df$pPC5
  )

  # Remove missing
  model_data_complete <- model_data[complete.cases(model_data), ]

  if (nrow(model_data_complete) < 50) {
    return(data.table(
      protein = protein_name,
      beta = NA_real_,
      se = NA_real_,
      t_stat = NA_real_,
      p_value = NA_real_,
      r_squared = NA_real_,
      n_samples = nrow(model_data_complete)
    ))
  }

  # Fit model
  fit <- tryCatch(
    lm(age ~ protein + bmi + sex + pPC1 + pPC2 + pPC3 + pPC4 + pPC5,
       data = model_data_complete),
    error = function(e) NULL
  )

  if (is.null(fit)) {
    return(data.table(
      protein = protein_name,
      beta = NA_real_,
      se = NA_real_,
      t_stat = NA_real_,
      p_value = NA_real_,
      r_squared = NA_real_,
      n_samples = nrow(model_data_complete)
    ))
  }

  # Extract coefficients
  coef_summary <- summary(fit)$coefficients
  r_squared <- summary(fit)$r.squared

  if ("protein" %in% rownames(coef_summary)) {
    data.table(
      protein = protein_name,
      beta = coef_summary["protein", "Estimate"],
      se = coef_summary["protein", "Std. Error"],
      t_stat = coef_summary["protein", "t value"],
      p_value = coef_summary["protein", "Pr(>|t|)"],
      r_squared = r_squared,
      n_samples = nrow(model_data_complete)
    )
  } else {
    data.table(
      protein = protein_name,
      beta = NA_real_,
      se = NA_real_,
      t_stat = NA_real_,
      p_value = NA_real_,
      r_squared = r_squared,
      n_samples = nrow(model_data_complete)
    )
  }
}


#' Create volcano plot with confidence intervals and outlier filtering
#'
#' @param results data.table with regression results
#' @param selected_only Logical, plot only selected proteins?
#' @param sd_threshold SD threshold for outlier filtering (default 3)
#' @return ggplot object
create_volcano_with_ci <- function(results, selected_only = TRUE, sd_threshold = 3) {

  log_info("========== VOLCANO PLOT CREATION: START ==========")

  # Filter to selected proteins if requested
  if (selected_only) {
    plot_data <- results[selected == TRUE & !is.na(p_value)]
    subtitle_text <- "Selected proteins only (multicollinearity resolved)"
    log_info("Using SELECTED proteins only: {nrow(plot_data)} proteins")
  } else {
    plot_data <- results[!is.na(p_value)]
    subtitle_text <- "All proteins"
    log_info("Using ALL proteins: {nrow(plot_data)} proteins")
  }

  # Filter outliers by effect size (mean ± SD threshold)
  log_info("Filtering extreme effect sizes (mean ± {sd_threshold}SD)...")
  mean_beta <- mean(plot_data$beta, na.rm = TRUE)
  sd_beta <- sd(plot_data$beta, na.rm = TRUE)
  beta_lower <- mean_beta - sd_threshold * sd_beta
  beta_upper <- mean_beta + sd_threshold * sd_beta

  n_before <- nrow(plot_data)
  plot_data <- plot_data[beta >= beta_lower & beta <= beta_upper]
  n_after <- nrow(plot_data)
  n_filtered <- n_before - n_after

  log_info("Effect size statistics:")
  log_info("  - Mean β: {round(mean_beta, 2)}")
  log_info("  - SD β: {round(sd_beta, 2)}")
  log_info("  - Range: [{round(beta_lower, 2)}, {round(beta_upper, 2)}]")
  log_info("  - Proteins before filtering: {n_before}")
  log_info("  - Proteins after filtering: {n_after}")
  log_info("  - Extreme outliers removed: {n_filtered}")

  # Add derived columns
  plot_data[, neg_log10_p := -log10(p_value + 1e-300)]
  plot_data[, abs_beta := abs(beta)]

  # Calculate Bonferroni threshold
  bonferroni_threshold <- 0.05 / nrow(plot_data)
  plot_data[, bonferroni_sig := p_value < bonferroni_threshold]

  # Calculate 95% confidence intervals
  plot_data[, ci_lower := beta - 1.96 * se]
  plot_data[, ci_upper := beta + 1.96 * se]

  # Annotate top proteins
  log_info("Selecting proteins for annotation (3 criteria)...")
  plot_data[, label := ""]

  # Top 10 positive
  top_pos <- plot_data[beta > 0][order(-beta)][1:min(10, .N)]
  log_info("  - Top 10 POSITIVE by β: {nrow(top_pos)} proteins")

  # Top 10 negative
  top_neg <- plot_data[beta < 0][order(beta)][1:min(10, .N)]
  log_info("  - Top 10 NEGATIVE by β: {nrow(top_neg)} proteins")

  # Top 10 by significance
  top_sig <- plot_data[order(p_value)][1:min(10, .N)]
  log_info("  - Top 10 by p-value: {nrow(top_sig)} proteins")

  proteins_to_label <- unique(c(top_pos$protein, top_neg$protein, top_sig$protein))
  plot_data[protein %in% proteins_to_label, label := protein]

  log_info("  - Total unique proteins to label: {length(proteins_to_label)}")
  log_info("Creating ggplot object with confidence intervals...")

  # Create plot
  col_after <- "#00AFBB"
  col_before <- "#FC4E07"

  p <- ggplot(plot_data, aes(x = beta, y = neg_log10_p)) +
    # Non-significant points
    geom_point(data = plot_data[bonferroni_sig == FALSE],
               color = "gray70", alpha = 0.4, size = 1.5) +
    # Significant points colored by effect size
    geom_point(data = plot_data[bonferroni_sig == TRUE],
               aes(color = abs_beta), alpha = 0.6, size = 2) +
    # Confidence intervals for labeled proteins
    geom_errorbarh(data = plot_data[label != ""],
                   aes(xmin = ci_lower, xmax = ci_upper, color = abs_beta),
                   height = 0, alpha = 0.5, size = 0.8) +
    # Highlighted labeled points
    geom_point(data = plot_data[label != ""], aes(color = abs_beta),
               size = 3.5, alpha = 0.9) +
    geom_text_repel(data = plot_data[label != ""], aes(label = label),
                    size = 2.8, max.overlaps = 30, box.padding = 0.5,
                    segment.alpha = 0.5) +
    scale_color_gradient(low = col_after, high = col_before,
                         name = expression("|"*beta*"|")) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray30") +
    geom_hline(yintercept = -log10(bonferroni_threshold),
               linetype = "dashed", color = "red", alpha = 0.5) +
    labs(title = "Age-Associated Proteins (Multicollinearity-Adjusted)",
         subtitle = sprintf("%s | Filtered to ±%dSD | Bonf: p < %.2e | n=%d",
                           subtitle_text, sd_threshold, bonferroni_threshold, n_after),
         x = expression("Effect Size ("*beta*" = years per NPX unit) with 95% CI"),
         y = "-log10(p-value)") +
    theme_bw() +
    theme(legend.position = "right")

  log_info("========== VOLCANO PLOT CREATION: COMPLETE ==========")

  return(p)
}

