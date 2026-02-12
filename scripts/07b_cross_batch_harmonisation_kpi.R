#!/usr/bin/env Rscript
# ==============================================================================
# 07b_cross_batch_harmonisation_kpi.R
#   Post-hoc evaluation of cross-batch bridge normalisation efficacy
# ==============================================================================
#
# Purpose:
#   Runs immediately after 07_bridge_normalization.R. Computes a comprehensive
#   set of Key Performance Indicators (KPIs) that quantify whether the two
#   batches are properly aligned after bridge normalisation, while checking
#   that biological signal has not been collapsed by over-correction.
#
# Metrics computed (6 evaluation steps):
#   Step 0 — Fixed PCA basis (avoids basis leakage)
#   Step 1 — Bridge pair collapse (Euclidean, Mahalanobis, ICC, CCC,
#            paired log-ratio with bootstrap CI, noise-floor ratio)
#   Step 2 — Pair identity test (rank-1 k-NN identification rate)
#   Step 3 — Batch separability (Silhouette, kBET, LISI)
#   Step 4 — Biology preservation (sex prediction accuracy, variance decomposition)
#   Step 5 — Outlier detection (bridge pair flags, per-protein concordance)
#   Step 6 — KPI summary table & multi-panel PDF dashboard
#
# Author: Reza Jabal, PhD (rjabal@broadinstitute.org)
# Date: February 2026
# ==============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(tidyverse)
  library(ggplot2)
  library(ggpubr)
  library(RColorBrewer)
  library(cluster)
  library(MASS)
  library(FNN)
  library(irr)
  library(DescTools)
  library(kBET)
  library(lisi)
  library(pROC)  # For AUC calculation
  library(pheatmap)
  library(yaml)
  library(logger)
})

# ==============================================================================
# Boilerplate: script directory, config, batch context, logging
# ==============================================================================
script_dir <- tryCatch({
  env_script <- Sys.getenv("SCRIPT_NAME", "")
  if (env_script != "" && file.exists(env_script)) {
    dirname(normalizePath(env_script))
  } else {
    args <- commandArgs(trailingOnly = FALSE)
    file_arg <- grep("^--file=", args, value = TRUE)
    if (length(file_arg) > 0) {
      dirname(normalizePath(sub("^--file=", "", file_arg)))
    } else {
      getwd()
    }
  }
}, error = function(e) getwd())
source(file.path(script_dir, "path_utils.R"))

batch_id <- Sys.getenv("PIPELINE_BATCH_ID", "batch_01")
step_num <- get_step_number()

config_file <- Sys.getenv("PIPELINE_CONFIG", "")
if (config_file == "" || !file.exists(config_file)) {
  stop("PIPELINE_CONFIG environment variable not set or config file not found.")
}
config <- read_yaml(config_file)

log_path <- get_log_path(step_num, batch_id, config = config)
ensure_output_dir(log_path)
log_appender(appender_file(log_path))
log_info("=" |> rep(70) |> paste(collapse = ""))
log_info("Step {step_num}: Cross-Batch Harmonisation KPI Evaluation")
log_info("Batch context: {batch_id}")
log_info("=" |> rep(70) |> paste(collapse = ""))

theme_set(theme_bw(base_size = 11))

# ==============================================================================
# Guard: skip gracefully if disabled or single-batch
# ==============================================================================
run_kpi <- tryCatch(
  isTRUE(config$parameters$bridge_normalization$run_kpi_evaluation),
  error = function(e) FALSE
)
multi_batch_mode <- tryCatch(
  isTRUE(config$parameters$normalization$multi_batch_mode),
  error = function(e) FALSE
)

if (!multi_batch_mode) {
  log_info("Single-batch mode — skipping cross-batch KPI evaluation")
  Sys.setenv(PIPELINE_STEP_SKIPPED = "TRUE")
  cat("\n=== CROSS-BATCH KPI EVALUATION SKIPPED (single-batch mode) ===\n")
  quit(save = "no", status = 0)
}
if (!run_kpi) {
  log_info("KPI evaluation disabled in config (bridge_normalization.run_kpi_evaluation: false)")
  Sys.setenv(PIPELINE_STEP_SKIPPED = "TRUE")
  cat("\n=== CROSS-BATCH KPI EVALUATION SKIPPED (disabled in config) ===\n")
  quit(save = "no", status = 0)
}

# ==============================================================================
# 0. DATA LOADING
# ==============================================================================
log_info("Phase 0: Loading data")

all_batch_ids <- names(config$batches)
if (length(all_batch_ids) < 2) {
  log_error("Cross-batch KPI requires >= 2 batches, found: {length(all_batch_ids)}")
  stop("Insufficient batches")
}
batch1_id <- all_batch_ids[1]
batch2_id <- all_batch_ids[2]

# -- Pre-normalisation matrices (step 05d QC-passed, same input as step 07) ---
load_matrix <- function(bid) {
  p <- get_output_path("05d", "05d_npx_matrix_all_qc_passed", bid, "phenotypes", "rds", config = config)
  if (!file.exists(p)) {
    p <- get_output_path("00", "npx_matrix_analysis_ready", bid, "qc", config = config)
  }
  if (!file.exists(p)) stop("Pre-norm matrix not found for ", bid, ": ", p)
  readRDS(p)
}
pre_b1 <- load_matrix(batch1_id)
pre_b2 <- load_matrix(batch2_id)
log_info("  Pre-norm B1: {nrow(pre_b1)} x {ncol(pre_b1)}")
log_info("  Pre-norm B2: {nrow(pre_b2)} x {ncol(pre_b2)}")

# -- Post-normalisation matrices (step 07 bridge-normalised) -----------------
post_b1_path <- get_output_path("07", "npx_matrix_cross_batch_bridge", batch1_id, "normalized", config = config)
post_b2_path <- get_output_path("07", "npx_matrix_cross_batch_bridge", batch2_id, "normalized", config = config)
if (!file.exists(post_b1_path) || !file.exists(post_b2_path)) {
  stop("Post-normalisation matrices not found. Has step 07 completed?")
}
post_b1 <- readRDS(post_b1_path)
post_b2 <- readRDS(post_b2_path)
log_info("  Post-norm B1: {nrow(post_b1)} x {ncol(post_b1)}")
log_info("  Post-norm B2: {nrow(post_b2)} x {ncol(post_b2)}")

# -- Sample mappings ----------------------------------------------------------
map_b1 <- readRDS(get_output_path("00", "sample_mapping", batch1_id, "qc", config = config))
map_b2 <- readRDS(get_output_path("00", "sample_mapping", batch2_id, "qc", config = config))

# -- Bridge pairs (24 direct B1↔B2, excluding EA5) --------------------------
bridging_file <- tryCatch(
  get_batch_input_path("bridging_samples_file", batch1_id, config),
  error = function(e) NULL
)
if (is.null(bridging_file) || !file.exists(bridging_file)) {
  stop("Bridging metadata file not found.")
}
bridging_meta <- fread(bridging_file, sep = "\t")

# Identify 24 direct B1↔B2 FINNGENIDs (excluding EA5 / Batch_00)
b01_flag <- bridging_meta$in_Batch_01 %in% c(TRUE, "TRUE", "T")
b02_flag <- bridging_meta$in_Batch_02 %in% c(TRUE, "TRUE", "T")
b00_flag <- bridging_meta$in_Batch_00 %in% c(FALSE, "FALSE", "F")
bridge_finngenids <- bridging_meta[b01_flag & b02_flag & b00_flag]$FINNGENID
log_info("  Direct B1-B2 bridge FINNGENIDs: {length(bridge_finngenids)}")

# Map FINNGENID → SampleID in each batch
bridge_sid_b1 <- map_b1[FINNGENID %in% bridge_finngenids & SampleID %in% rownames(pre_b1)]
bridge_sid_b2 <- map_b2[FINNGENID %in% bridge_finngenids & SampleID %in% rownames(pre_b2)]
# If a FINNGENID has >1 SampleID in a batch, keep the first (deterministic)
bridge_sid_b1 <- bridge_sid_b1[!duplicated(FINNGENID)]
bridge_sid_b2 <- bridge_sid_b2[!duplicated(FINNGENID)]
# Align on common FINNGENIDs
common_bridge_fg <- intersect(bridge_sid_b1$FINNGENID, bridge_sid_b2$FINNGENID)
bridge_sid_b1 <- bridge_sid_b1[match(common_bridge_fg, FINNGENID)]
bridge_sid_b2 <- bridge_sid_b2[match(common_bridge_fg, FINNGENID)]
n_bridge <- length(common_bridge_fg)
log_info("  Aligned bridge pairs: {n_bridge}")

# -- Common proteins ----------------------------------------------------------
common_proteins <- Reduce(intersect, list(
  colnames(pre_b1), colnames(pre_b2),
  colnames(post_b1), colnames(post_b2)
))
log_info("  Common proteins across all matrices: {length(common_proteins)}")

# -- Within-batch replicates (noise floor) ------------------------------------
find_within_batch_replicates <- function(mapping, npx_matrix) {
  dup_fgs <- mapping[!is.na(FINNGENID), .N, by = FINNGENID][N >= 2]$FINNGENID
  pairs <- list()
  for (fg in dup_fgs) {
    sids <- mapping[FINNGENID == fg]$SampleID
    sids <- intersect(sids, rownames(npx_matrix))
    if (length(sids) >= 2) {
      pairs[[length(pairs) + 1]] <- list(
        finngenid = fg, sid1 = sids[1], sid2 = sids[2]
      )
    }
  }
  rbindlist(pairs)
}
wb_rep_b1 <- find_within_batch_replicates(map_b1, pre_b1)
wb_rep_b2 <- find_within_batch_replicates(map_b2, pre_b2)
log_info("  Within-batch replicate pairs: B1={nrow(wb_rep_b1)}, B2={nrow(wb_rep_b2)}")

# ==============================================================================
# STEP 0: FIXED PCA BASIS (avoid basis leakage)
# ==============================================================================
log_info("Step 0: Building fixed PCA basis on combined pre-harmonisation data")

build_fixed_pca <- function(mat_b1, mat_b2, proteins, max_pcs = 50, var_threshold = 0.95) {
  combined <- rbind(
    mat_b1[, proteins, drop = FALSE],
    mat_b2[, proteins, drop = FALSE]
  )
  # Impute remaining NAs with column medians
  na_count <- sum(is.na(combined))
  if (na_count > 0) {
    log_info("  Imputing {na_count} NAs with column medians")
    for (j in seq_len(ncol(combined))) {
      nas <- is.na(combined[, j])
      if (any(nas)) combined[nas, j] <- median(combined[, j], na.rm = TRUE)
    }
  }
  # Centre and scale, then PCA
  pca_fit <- prcomp(combined, center = TRUE, scale. = TRUE)
  cum_var <- cumsum(pca_fit$sdev^2 / sum(pca_fit$sdev^2))
  n_pcs <- min(max_pcs, which(cum_var >= var_threshold)[1])
  if (is.na(n_pcs)) n_pcs <- max_pcs
  log_info("  Retaining {n_pcs} PCs ({round(cum_var[n_pcs]*100,1)}% variance)")
  list(
    pca_model  = pca_fit,
    n_pcs      = n_pcs,
    var_explained = cum_var[n_pcs],
    center     = pca_fit$center,
    scale      = pca_fit$scale,
    rotation   = pca_fit$rotation[, 1:n_pcs, drop = FALSE]
  )
}

project_onto_basis <- function(mat_b1, mat_b2, pca_basis, proteins) {
  project_one <- function(mat) {
    m <- mat[, proteins, drop = FALSE]
    # Impute NAs
    for (j in seq_len(ncol(m))) {
      nas <- is.na(m[, j])
      if (any(nas)) m[nas, j] <- pca_basis$center[j]
    }
    # Centre, scale, multiply by rotation
    m_scaled <- scale(m, center = pca_basis$center, scale = pca_basis$scale)
    m_scaled %*% pca_basis$rotation
  }
  list(
    b1 = project_one(mat_b1),
    b2 = project_one(mat_b2)
  )
}

pca_basis <- build_fixed_pca(pre_b1, pre_b2, common_proteins)

pc_pre  <- project_onto_basis(pre_b1, pre_b2, pca_basis, common_proteins)
pc_post <- project_onto_basis(post_b1, post_b2, pca_basis, common_proteins)

log_info("  PC scores pre: B1={nrow(pc_pre$b1)}, B2={nrow(pc_pre$b2)}")
log_info("  PC scores post: B1={nrow(pc_post$b1)}, B2={nrow(pc_post$b2)}")

# ==============================================================================
# STEP 1: BRIDGE PAIR COLLAPSE (Primary KPI)
# ==============================================================================
log_info("Step 1: Computing bridge pair collapse metrics")

# Helper: Euclidean distance for matched pairs
paired_euclidean <- function(pc_b1, pc_b2, sids_b1, sids_b2) {
  sqrt(rowSums((pc_b1[sids_b1, , drop = FALSE] - pc_b2[sids_b2, , drop = FALSE])^2))
}

# Helper: Mahalanobis distance for matched pairs (robust covariance)
paired_mahalanobis <- function(pc_b1, pc_b2, sids_b1, sids_b2) {
  combined_pc <- rbind(pc_b1, pc_b2)
  # Use robust covariance (shrinkage/classical fallback)
  cov_mat <- tryCatch(
    cov.shrink <- (t(combined_pc) %*% combined_pc) / (nrow(combined_pc) - 1),
    error = function(e) cov(combined_pc)
  )
  cov_inv <- tryCatch(solve(cov_mat), error = function(e) {
    log_warn("Covariance matrix singular, using pseudoinverse")
    MASS::ginv(cov_mat)
  })
  diffs <- pc_b1[sids_b1, , drop = FALSE] - pc_b2[sids_b2, , drop = FALSE]
  sqrt(rowSums((diffs %*% cov_inv) * diffs))
}

sids_b1 <- bridge_sid_b1$SampleID
sids_b2 <- bridge_sid_b2$SampleID

# Compute MAD-whitened distances (robust per-PC scaled Euclidean distance)
# This is the primary distance metric for harmonisation assessment
compute_mad_whitened_distance <- function(pc_b1, pc_b2, sids_b1, sids_b2) {
  # Combine all PC scores to compute robust scaling factors
  combined_pc <- rbind(pc_b1, pc_b2)
  # Compute MAD per PC dimension
  mad_per_pc <- apply(combined_pc, 2, function(x) mad(x, na.rm = TRUE))
  # Avoid division by zero (set to 1 if MAD is 0)
  mad_per_pc[mad_per_pc == 0] <- 1

  # Compute MAD-whitened distance for each bridge pair
  dists <- numeric(length(sids_b1))
  for (i in seq_along(sids_b1)) {
    diff_vec <- pc_b1[sids_b1[i], , drop = FALSE] - pc_b2[sids_b2[i], , drop = FALSE]
    # Scale by MAD: d_ij / s_j
    scaled_diff <- diff_vec / mad_per_pc
    dists[i] <- sqrt(sum(scaled_diff^2, na.rm = TRUE))
  }
  dists
}

# 1a. Distance metrics
d_eucl_pre  <- paired_euclidean(pc_pre$b1,  pc_pre$b2,  sids_b1, sids_b2)
d_eucl_post <- paired_euclidean(pc_post$b1, pc_post$b2, sids_b1, sids_b2)
d_maha_pre  <- paired_mahalanobis(pc_pre$b1,  pc_pre$b2,  sids_b1, sids_b2)
d_maha_post <- paired_mahalanobis(pc_post$b1, pc_post$b2, sids_b1, sids_b2)

# MAD-whitened distances (primary metric for harmonisation assessment)
d_mad_pre  <- compute_mad_whitened_distance(pc_pre$b1,  pc_pre$b2,  sids_b1, sids_b2)
d_mad_post <- compute_mad_whitened_distance(pc_post$b1, pc_post$b2, sids_b1, sids_b2)

log_info("  Euclidean — pre: median={round(median(d_eucl_pre),3)}, post: median={round(median(d_eucl_post),3)}")
log_info("  Mahalanobis — pre: median={round(median(d_maha_pre),3)}, post: median={round(median(d_maha_post),3)}")

# 1a. ICC and CCC per bridge pair (on raw protein space, not PC space)
compute_pair_agreement <- function(mat_b1, mat_b2, sids_b1, sids_b2, proteins) {
  icc_vals <- numeric(length(sids_b1))
  ccc_vals <- numeric(length(sids_b1))
  for (i in seq_along(sids_b1)) {
    v1 <- mat_b1[sids_b1[i], proteins]
    v2 <- mat_b2[sids_b2[i], proteins]
    ok <- !is.na(v1) & !is.na(v2)
    if (sum(ok) < 10) {
      icc_vals[i] <- NA_real_
      ccc_vals[i] <- NA_real_
      next
    }
    # ICC: two-way mixed, single measures, consistency
    icc_res <- tryCatch(
      irr::icc(cbind(v1[ok], v2[ok]), model = "twoway", type = "consistency", unit = "single"),
      error = function(e) list(value = NA_real_)
    )
    icc_vals[i] <- icc_res$value
    # CCC
    ccc_res <- tryCatch(
      DescTools::CCC(v1[ok], v2[ok])$rho.c$est,
      error = function(e) NA_real_
    )
    ccc_vals[i] <- ccc_res
  }
  list(icc = icc_vals, ccc = ccc_vals)
}

agree_pre  <- compute_pair_agreement(pre_b1,  pre_b2,  sids_b1, sids_b2, common_proteins)
agree_post <- compute_pair_agreement(post_b1, post_b2, sids_b1, sids_b2, common_proteins)

# Verify lengths match n_bridge
if (length(agree_pre$icc) != n_bridge || length(agree_pre$ccc) != n_bridge) {
  log_error("agree_pre length mismatch: icc={length(agree_pre$icc)}, ccc={length(agree_pre$ccc)}, expected={n_bridge}")
  stop("agree_pre length mismatch")
}
if (length(agree_post$icc) != n_bridge || length(agree_post$ccc) != n_bridge) {
  log_error("agree_post length mismatch: icc={length(agree_post$icc)}, ccc={length(agree_post$ccc)}, expected={n_bridge}")
  stop("agree_post length mismatch")
}

log_info("  ICC — pre: median={round(median(agree_pre$icc, na.rm=T),3)}, post: median={round(median(agree_post$icc, na.rm=T),3)}")
log_info("  CCC — pre: median={round(median(agree_pre$ccc, na.rm=T),3)}, post: median={round(median(agree_post$ccc, na.rm=T),3)}")

# 1b. Noise floor ratio (within-batch replicate distances)
compute_noise_floor <- function(wb_pairs, pc_scores) {
  if (nrow(wb_pairs) == 0) return(NA_real_)
  dists <- numeric(nrow(wb_pairs))
  for (i in seq_len(nrow(wb_pairs))) {
    s1 <- wb_pairs$sid1[i]; s2 <- wb_pairs$sid2[i]
    if (s1 %in% rownames(pc_scores) && s2 %in% rownames(pc_scores)) {
      dists[i] <- sqrt(sum((pc_scores[s1, ] - pc_scores[s2, ])^2))
    } else {
      dists[i] <- NA_real_
    }
  }
  median(dists, na.rm = TRUE)
}

# Compute noise floor from pre-norm within-batch replicates
# (this is the measurement noise, should not change with harmonisation)
noise_floor_b1 <- compute_noise_floor(wb_rep_b1, pc_pre$b1)
noise_floor_b2 <- compute_noise_floor(wb_rep_b2, pc_pre$b2)
# Use the pooled median
all_wb_pairs <- rbind(wb_rep_b1, wb_rep_b2)

# Compute global MAD scaling factors (same as used for bridge pair MAD-whitened distances)
# so that noise floor and bridge distances are on the same scale
global_combined_pc <- rbind(pc_pre$b1, pc_pre$b2)
global_mad_per_pc <- apply(global_combined_pc, 2, function(x) mad(x, na.rm = TRUE))
global_mad_per_pc[global_mad_per_pc == 0] <- 1

noise_floor_dists <- numeric(nrow(all_wb_pairs))
noise_floor_dists_mad <- numeric(nrow(all_wb_pairs))
for (i in seq_len(nrow(all_wb_pairs))) {
  s1 <- all_wb_pairs$sid1[i]; s2 <- all_wb_pairs$sid2[i]
  # Determine which batch this pair belongs to
  if (s1 %in% rownames(pc_pre$b1) && s2 %in% rownames(pc_pre$b1)) {
    diff_vec <- pc_pre$b1[s1, , drop = FALSE] - pc_pre$b1[s2, , drop = FALSE]
    noise_floor_dists[i] <- sqrt(sum(diff_vec^2))
    noise_floor_dists_mad[i] <- sqrt(sum((diff_vec / global_mad_per_pc)^2, na.rm = TRUE))
  } else if (s1 %in% rownames(pc_pre$b2) && s2 %in% rownames(pc_pre$b2)) {
    diff_vec <- pc_pre$b2[s1, , drop = FALSE] - pc_pre$b2[s2, , drop = FALSE]
    noise_floor_dists[i] <- sqrt(sum(diff_vec^2))
    noise_floor_dists_mad[i] <- sqrt(sum((diff_vec / global_mad_per_pc)^2, na.rm = TRUE))
  } else {
    noise_floor_dists[i] <- NA_real_
    noise_floor_dists_mad[i] <- NA_real_
  }
}
noise_floor <- median(noise_floor_dists, na.rm = TRUE)
noise_floor_mad <- median(noise_floor_dists_mad, na.rm = TRUE)
noise_floor_ratio <- median(d_mad_post) / noise_floor_mad

log_info("  Noise floor (median within-batch replicate distance): {round(noise_floor, 3)}")
log_info("  Noise floor ratio (bridge_post / noise): {round(noise_floor_ratio, 3)}")

# 1c. Paired log-ratio effect size with bootstrap CI
# Primary metric: MAD-whitened distances (robust, per-PC scaled)
log_ratio_mad <- log(d_mad_post / d_mad_pre)
median_log_ratio_mad <- median(log_ratio_mad)
# Also compute for Euclidean (for backward compatibility in tables)
log_ratio <- log(d_eucl_post / d_eucl_pre)
median_log_ratio <- median(log_ratio)

# Bootstrap 95% CI (10,000 resamples) - for MAD-whitened (primary)
set.seed(42)
n_boot <- 10000
boot_medians_mad <- replicate(n_boot, {
  idx <- sample.int(n_bridge, replace = TRUE)
  median(log_ratio_mad[idx])
})
boot_ci_mad <- quantile(boot_medians_mad, probs = c(0.025, 0.975))

# Also for Euclidean (for tables)
boot_medians <- replicate(n_boot, {
  idx <- sample.int(n_bridge, replace = TRUE)
  median(log_ratio[idx])
})
boot_ci <- quantile(boot_medians, probs = c(0.025, 0.975))

log_info("  Paired log-ratio (MAD-whitened): median={round(median_log_ratio_mad, 4)}, 95% CI=[{round(boot_ci_mad[1], 4)}, {round(boot_ci_mad[2], 4)}]")
log_info("  CI excludes 0: {boot_ci_mad[2] < 0}")

# Wilcoxon signed-rank test (paired) - use MAD-whitened (primary)
wilcox_result <- wilcox.test(d_mad_post, d_mad_pre, paired = TRUE, alternative = "less")
log_info("  Wilcoxon signed-rank p-value (MAD-whitened): {format.pval(wilcox_result$p.value, digits = 4)}")

# Build per-pair distance table
# Verify all vectors have correct length before creating data.table
expected_length <- n_bridge
vectors_to_check <- list(
  pair_id = seq_len(n_bridge),
  finngenid = common_bridge_fg,
  sid_b1 = sids_b1,
  sid_b2 = sids_b2,
  eucl_pre = d_eucl_pre,
  eucl_post = d_eucl_post,
  maha_pre = d_maha_pre,
  maha_post = d_maha_post,
  mad_pre = d_mad_pre,
  mad_post = d_mad_post,
  icc_pre = agree_pre$icc,
  icc_post = agree_post$icc,
  ccc_pre = agree_pre$ccc,
  ccc_post = agree_post$ccc,
  log_ratio = log_ratio,
  log_ratio_mad = log(d_mad_post / d_mad_pre),
  noise_floor_ratio = d_mad_post / noise_floor_mad
)

# Check lengths
for (name in names(vectors_to_check)) {
  actual_length <- length(vectors_to_check[[name]])
  if (actual_length != expected_length) {
    log_error("Length mismatch in bridge_pair_dt: {name} has {actual_length} elements, expected {expected_length}")
    stop("Vector length mismatch in bridge_pair_dt")
  }
}

# Create bridge_pair_dt with explicit length validation
bridge_pair_dt <- tryCatch({
  dt <- data.table(
    pair_id        = vectors_to_check$pair_id,
    finngenid      = vectors_to_check$finngenid,
    sid_b1         = vectors_to_check$sid_b1,
    sid_b2         = vectors_to_check$sid_b2,
    eucl_pre       = vectors_to_check$eucl_pre,
    eucl_post      = vectors_to_check$eucl_post,
    maha_pre       = vectors_to_check$maha_pre,
    maha_post      = vectors_to_check$maha_post,
  mad_pre        = vectors_to_check$mad_pre,
  mad_post       = vectors_to_check$mad_post,
  icc_pre        = vectors_to_check$icc_pre,
  icc_post       = vectors_to_check$icc_post,
  ccc_pre        = vectors_to_check$ccc_pre,
  ccc_post       = vectors_to_check$ccc_post,
  log_ratio      = vectors_to_check$log_ratio,
  log_ratio_mad  = vectors_to_check$log_ratio_mad,
  noise_floor_ratio = vectors_to_check$noise_floor_ratio
  )
  # Validate all columns have correct length
  col_lengths <- sapply(dt, length)
  if (any(col_lengths != n_bridge)) {
    log_error("bridge_pair_dt column length mismatch:")
    for (col in names(col_lengths)) {
      if (col_lengths[col] != n_bridge) {
        log_error("  {col}: {col_lengths[col]} (expected {n_bridge})")
      }
    }
    stop("bridge_pair_dt column length mismatch")
  }
  dt
}, error = function(e) {
  log_error("Failed to create bridge_pair_dt: {e$message}")
  stop(e)
})

# ==============================================================================
# STEP 2: PAIR IDENTITY TEST (Rank-1 identification rate)
# ==============================================================================
log_info("Step 2: Pair identity test (k-NN rank-1 identification rate)")

compute_rank1_rate <- function(pc_b1, pc_b2, sids_b1, sids_b2) {
  # For each bridge sample in B1, find nearest neighbour in B2
  query <- pc_b1[sids_b1, , drop = FALSE]
  ref   <- pc_b2  # all B2 samples
  # get.knnx returns indices into ref sorted by distance
  knn_res <- FNN::get.knnx(ref, query, k = nrow(ref))
  # For each bridge sample, what rank does its true pair achieve?
  true_b2_idx <- match(sids_b2, rownames(ref))
  ranks <- numeric(length(sids_b1))
  for (i in seq_along(sids_b1)) {
    ranks[i] <- which(knn_res$nn.index[i, ] == true_b2_idx[i])
  }
  list(
    rank1_hits = sum(ranks == 1),
    rank1_rate = mean(ranks == 1),
    median_rank = median(ranks),
    ranks = ranks
  )
}

id_pre  <- compute_rank1_rate(pc_pre$b1,  pc_pre$b2,  sids_b1, sids_b2)
id_post <- compute_rank1_rate(pc_post$b1, pc_post$b2, sids_b1, sids_b2)

# Verify ranks have correct length
if (length(id_pre$ranks) != n_bridge) {
  log_error("id_pre$ranks length mismatch: {length(id_pre$ranks)} vs expected {n_bridge}")
  stop("id_pre$ranks length mismatch")
}
if (length(id_post$ranks) != n_bridge) {
  log_error("id_post$ranks length mismatch: {length(id_post$ranks)} vs expected {n_bridge}")
  stop("id_post$ranks length mismatch")
}

log_info("  Rank-1 ID rate — pre: {id_pre$rank1_hits}/{n_bridge} ({round(id_pre$rank1_rate*100,1)}%), post: {id_post$rank1_hits}/{n_bridge} ({round(id_post$rank1_rate*100,1)}%)")
log_info("  Median rank of true pair — pre: {id_pre$median_rank}, post: {id_post$median_rank}")

# Validate ranks have correct length before assigning
if (length(id_pre$ranks) != n_bridge) {
  log_error("id_pre$ranks length mismatch: {length(id_pre$ranks)} vs expected {n_bridge}")
  stop("id_pre$ranks length mismatch")
}
if (length(id_post$ranks) != n_bridge) {
  log_error("id_post$ranks length mismatch: {length(id_post$ranks)} vs expected {n_bridge}")
  stop("id_post$ranks length mismatch")
}

# Assign ranks using column assignment (safer than := when there might be recycling issues)
bridge_pair_dt$rank_pre <- id_pre$ranks
bridge_pair_dt$rank_post <- id_post$ranks

# Verify assignment succeeded
if (length(bridge_pair_dt$rank_pre) != n_bridge || length(bridge_pair_dt$rank_post) != n_bridge) {
  log_error("Rank assignment failed: rank_pre length={length(bridge_pair_dt$rank_pre)}, rank_post length={length(bridge_pair_dt$rank_post)}, expected {n_bridge}")
  stop("Rank assignment length mismatch")
}

# ==============================================================================
# STEP 3: BATCH SEPARABILITY (full cohort, ~4K samples)
# ==============================================================================
log_info("Step 3: Batch separability metrics (full cohort)")

# Combine PC scores with batch labels
build_cohort_pc <- function(pc_b1, pc_b2, bid1, bid2) {
  scores <- rbind(pc_b1, pc_b2)
  labels <- c(rep(bid1, nrow(pc_b1)), rep(bid2, nrow(pc_b2)))
  list(scores = scores, labels = labels)
}

cohort_pre  <- build_cohort_pc(pc_pre$b1,  pc_pre$b2,  batch1_id, batch2_id)
cohort_post <- build_cohort_pc(pc_post$b1, pc_post$b2, batch1_id, batch2_id)

# 3a. Silhouette score (on first 2 PCs for interpretability)
compute_silhouette <- function(scores, labels) {
  d <- dist(scores[, 1:2])
  si <- cluster::silhouette(as.integer(factor(labels)), d)
  mean(si[, "sil_width"])
}

sil_pre  <- compute_silhouette(cohort_pre$scores,  cohort_pre$labels)
sil_post <- compute_silhouette(cohort_post$scores, cohort_post$labels)
log_info("  Silhouette (PC1-PC2) — pre: {round(sil_pre, 4)}, post: {round(sil_post, 4)}")

# 3b. kBET (k-Nearest Neighbour Batch Effect Test) acceptance rate
compute_kbet <- function(scores, labels, k0 = 25) {
  k0 <- min(k0, floor(nrow(scores) * 0.1))
  res <- tryCatch(
    kBET::kBET(scores, labels, k0 = k0, plot = FALSE, verbose = FALSE),
    error = function(e) {
      log_warn("kBET (k-Nearest Neighbour Batch Effect Test) failed: {e$message}")
      list(summary = data.frame(kBET.observed = NA_real_))
    }
  )
  1 - res$summary$kBET.observed[1]  # acceptance rate = 1 - rejection rate
}

kbet_pre  <- compute_kbet(cohort_pre$scores,  cohort_pre$labels)
kbet_post <- compute_kbet(cohort_post$scores, cohort_post$labels)
log_info("  kBET (k-Nearest Neighbour Batch Effect Test) acceptance — pre: {round(kbet_pre, 4)}, post: {round(kbet_post, 4)}")

# 3c. LISI (Local Inverse Simpson's Index) batch mixing score
compute_lisi_score <- function(scores, labels) {
  meta_df <- data.frame(batch = labels, row.names = rownames(scores))
  res <- tryCatch(
    lisi::compute_lisi(scores, meta_df, "batch"),
    error = function(e) {
      log_warn("LISI (Local Inverse Simpson's Index) failed: {e$message}")
      data.frame(batch = NA_real_)
    }
  )
  median(res$batch, na.rm = TRUE)
}

lisi_pre  <- compute_lisi_score(cohort_pre$scores,  cohort_pre$labels)
lisi_post <- compute_lisi_score(cohort_post$scores, cohort_post$labels)
log_info("  LISI (Local Inverse Simpson's Index) — pre: {round(lisi_pre, 4)}, post: {round(lisi_post, 4)}")

# ==============================================================================
# STEP 4: BIOLOGY PRESERVATION (guard against over-correction)
# ==============================================================================
log_info("Step 4: Biology preservation checks")

# 4a. Sex prediction accuracy (re-use top sex proteins from step 04)
sex_pred_b1_path <- get_output_path("04", "sex_predictions", batch1_id, "outliers", "tsv", config = config)
sex_pred_b2_path <- get_output_path("04", "sex_predictions", batch2_id, "outliers", "tsv", config = config)

sex_accuracy_pre  <- NA_real_
sex_accuracy_post <- NA_real_

if (file.exists(sex_pred_b1_path) && file.exists(sex_pred_b2_path)) {
  log_info("  Loading sex predictions from step 04")
  sp_b1 <- fread(sex_pred_b1_path)
  sp_b2 <- fread(sex_pred_b2_path)

  # Get top sex-associated proteins
  sex_prot_path <- get_output_path("04", "sex_associated_proteins", batch1_id, "outliers", "tsv", config = config)
  top_sex_proteins <- character(0)
  if (file.exists(sex_prot_path)) {
    sex_prots <- fread(sex_prot_path)
    top_sex_proteins <- head(sex_prots$protein, 20)
    top_sex_proteins <- intersect(top_sex_proteins, common_proteins)
  }

  if (length(top_sex_proteins) >= 5) {
    # Simple logistic regression check: predict sex from top proteins
    compute_sex_accuracy <- function(mat_b1, mat_b2, map_b1, map_b2, sex_dt_b1, sex_dt_b2, proteins) {
      # Build combined dataset with sex labels
      build_one <- function(mat, mapping, sex_dt) {
        sids_in_mat <- intersect(rownames(mat), mapping$SampleID)
        fg_map <- mapping[SampleID %in% sids_in_mat, .(SampleID, FINNGENID)]
        fg_map <- fg_map[!duplicated(SampleID)]
        # Match sex data
        id_col <- if ("SAMPLE_ID" %in% names(sex_dt)) "SAMPLE_ID" else "SampleID"
        sex_map <- merge(fg_map, sex_dt[, .(get(id_col), genetic_sex)], by.x = "SampleID", by.y = "V1", all.x = TRUE)
        sex_map <- sex_map[!is.na(genetic_sex) & genetic_sex %in% c("male", "female")]
        if (nrow(sex_map) == 0) return(NULL)
        X <- mat[sex_map$SampleID, proteins, drop = FALSE]
        y <- as.integer(sex_map$genetic_sex == "female")
        list(X = X, y = y)
      }
      d1 <- build_one(mat_b1, map_b1, sex_dt_b1)
      d2 <- build_one(mat_b2, map_b2, sex_dt_b2)
      if (is.null(d1) || is.null(d2)) return(NA_real_)
      X_all <- rbind(d1$X, d2$X)
      y_all <- c(d1$y, d2$y)

      # Handle NAs
      complete <- complete.cases(X_all)
      X_all <- X_all[complete, , drop = FALSE]
      y_all <- y_all[complete]

      if (length(unique(y_all)) < 2 || length(y_all) < 50) return(NA_real_)
      fit <- tryCatch(
        glm(y_all ~ ., data = as.data.frame(X_all), family = binomial()),
        error = function(e) NULL
      )
      if (is.null(fit)) return(NA_real_)
      pred <- ifelse(predict(fit, type = "response") >= 0.5, 1, 0)
      mean(pred == y_all)
    }

    sex_accuracy_pre <- compute_sex_accuracy(
      pre_b1, pre_b2, map_b1, map_b2, sp_b1, sp_b2, top_sex_proteins
    )
    sex_accuracy_post <- compute_sex_accuracy(
      post_b1, post_b2, map_b1, map_b2, sp_b1, sp_b2, top_sex_proteins
    )
    log_info("  Sex accuracy — pre: {round(sex_accuracy_pre, 4)}, post: {round(sex_accuracy_post, 4)}")
  } else {
    log_warn("  Insufficient sex-associated proteins ({length(top_sex_proteins)}) — skipping sex accuracy check")
  }
} else {
  log_warn("  Sex prediction files not found — skipping sex accuracy check")
}

# 4b. Variance decomposition (batch vs sex for a subset of proteins)
log_info("  Computing variance decomposition (batch vs sex)")

compute_variance_decomposition <- function(mat_b1, mat_b2, map_b1, map_b2,
                                           sex_dt_b1, sex_dt_b2, proteins,
                                           bid1, bid2) {
  # Build combined matrix with batch and sex labels
  build_labeled <- function(mat, mapping, sex_dt, bid) {
    sids <- intersect(rownames(mat), mapping$SampleID)
    fg_map <- mapping[SampleID %in% sids, .(SampleID, FINNGENID)]
    fg_map <- fg_map[!duplicated(SampleID)]
    id_col <- if ("SAMPLE_ID" %in% names(sex_dt)) "SAMPLE_ID" else "SampleID"
    sex_map <- merge(fg_map, sex_dt[, .(get(id_col), genetic_sex)], by.x = "SampleID", by.y = "V1", all.x = TRUE)
    sex_map <- sex_map[!is.na(genetic_sex) & genetic_sex %in% c("male", "female")]
    if (nrow(sex_map) == 0) return(NULL)
    X <- mat[sex_map$SampleID, proteins, drop = FALSE]
    data.table(
      SampleID = sex_map$SampleID,
      batch = bid,
      sex = sex_map$genetic_sex,
      as.data.table(X)
    )
  }

  d1 <- build_labeled(mat_b1, map_b1, sex_dt_b1, bid1)
  d2 <- build_labeled(mat_b2, map_b2, sex_dt_b2, bid2)
  if (is.null(d1) || is.null(d2)) return(NULL)
  combined <- rbind(d1, d2, fill = TRUE)

  # For each protein, fit NPX ~ batch + sex and extract R^2 components
  result <- rbindlist(lapply(proteins, function(p) {
    y <- combined[[p]]
    if (all(is.na(y))) return(NULL)
    ok <- !is.na(y)
    fit_full <- tryCatch(
      anova(lm(y[ok] ~ factor(combined$batch[ok]) + factor(combined$sex[ok]))),
      error = function(e) NULL
    )
    if (is.null(fit_full)) return(NULL)
    ss <- fit_full[["Sum Sq"]]
    total_ss <- sum(ss)
    data.table(
      protein = p,
      batch_var_frac = ss[1] / total_ss,
      sex_var_frac   = ss[2] / total_ss,
      residual_frac  = ss[3] / total_ss
    )
  }))
  result
}

# Use a representative subset of proteins for variance decomposition
var_decomp_proteins <- if (length(common_proteins) > 200) {
  sample(common_proteins, 200)
} else {
  common_proteins
}

var_decomp_pre  <- NULL
var_decomp_post <- NULL
if (file.exists(sex_pred_b1_path) && file.exists(sex_pred_b2_path)) {
  var_decomp_pre <- compute_variance_decomposition(
    pre_b1, pre_b2, map_b1, map_b2, sp_b1, sp_b2,
    var_decomp_proteins, batch1_id, batch2_id
  )
  var_decomp_post <- compute_variance_decomposition(
    post_b1, post_b2, map_b1, map_b2, sp_b1, sp_b2,
    var_decomp_proteins, batch1_id, batch2_id
  )
  if (!is.null(var_decomp_pre) && !is.null(var_decomp_post)) {
    log_info("  Batch variance fraction — pre: {round(median(var_decomp_pre$batch_var_frac, na.rm=T), 4)}, post: {round(median(var_decomp_post$batch_var_frac, na.rm=T), 4)}")
    log_info("  Sex variance fraction — pre: {round(median(var_decomp_pre$sex_var_frac, na.rm=T), 4)}, post: {round(median(var_decomp_post$sex_var_frac, na.rm=T), 4)}")

    # Over-correction check: flag if sex variance drops by >20%
    sex_var_drop <- (median(var_decomp_pre$sex_var_frac, na.rm = TRUE) -
                     median(var_decomp_post$sex_var_frac, na.rm = TRUE)) /
                    median(var_decomp_pre$sex_var_frac, na.rm = TRUE)
    if (!is.na(sex_var_drop) && sex_var_drop > 0.2) {
      log_warn("  OVER-CORRECTION ALARM: sex variance dropped by {round(sex_var_drop*100,1)}% after harmonisation")
    }
  }
}

# ==============================================================================
# STEP 5: OUTLIER DETECTION
# ==============================================================================
log_info("Step 5: Outlier detection")

# 5a. Bridge pair outlier flagging (post-harmonisation Mahalanobis)
q1 <- quantile(d_maha_post, 0.25)
q3 <- quantile(d_maha_post, 0.75)
iqr_val <- q3 - q1
outlier_threshold <- q3 + 1.5 * iqr_val
bridge_pair_dt[, outlier_maha := d_maha_post > outlier_threshold]
bridge_pair_dt[, outlier_noise := noise_floor_ratio > 3.0]
bridge_pair_dt[, outlier_flag := outlier_maha | outlier_noise]

n_outliers <- sum(bridge_pair_dt$outlier_flag)
log_info("  Bridge pair outliers (Mahalanobis + noise ratio): {n_outliers}/{n_bridge}")
if (n_outliers > 0) {
  log_warn("  Flagged bridge pairs: {paste(bridge_pair_dt[outlier_flag == TRUE]$finngenid, collapse=', ')}")
}

# 5b. Per-protein concordance profiling
# NOTE: Correlation is preserved under additive transformations (Olink standard method),
# so we use Mean Squared Error (MSE) instead to measure harmonisation success
log_info("  Computing per-protein harmonisation (bridge pairs, MSE-based)")

protein_concordance <- rbindlist(lapply(common_proteins, function(p) {
  v1_pre  <- pre_b1[sids_b1,  p]
  v2_pre  <- pre_b2[sids_b2,  p]
  v1_post <- post_b1[sids_b1, p]
  v2_post <- post_b2[sids_b2, p]

  ok_pre  <- !is.na(v1_pre) & !is.na(v2_pre)
  ok_post <- !is.na(v1_post) & !is.na(v2_post)

  # Correlation (for reference, but not used for flagging)
  cor_pre  <- if (sum(ok_pre) >= 5) cor(v1_pre[ok_pre], v2_pre[ok_pre]) else NA_real_
  cor_post <- if (sum(ok_post) >= 5) cor(v1_post[ok_post], v2_post[ok_post]) else NA_real_

  # Mean Squared Error (MSE) - should decrease after harmonisation
  mse_pre  <- if (sum(ok_pre) >= 5) mean((v1_pre[ok_pre] - v2_pre[ok_pre])^2, na.rm = TRUE) else NA_real_
  mse_post <- if (sum(ok_post) >= 5) mean((v1_post[ok_post] - v2_post[ok_post])^2, na.rm = TRUE) else NA_real_

  # Mean Absolute Error (MAE) - alternative metric
  mae_pre  <- if (sum(ok_pre) >= 5) mean(abs(v1_pre[ok_pre] - v2_pre[ok_pre]), na.rm = TRUE) else NA_real_
  mae_post <- if (sum(ok_post) >= 5) mean(abs(v1_post[ok_post] - v2_post[ok_post]), na.rm = TRUE) else NA_real_

  # R² (coefficient of determination) - should increase toward 1.0
  # R² = 1 - (SS_res / SS_tot), where SS_res = sum((y - y_pred)^2), SS_tot = sum((y - mean(y))^2)
  r2_pre  <- if (sum(ok_pre) >= 5) {
    ss_res <- sum((v2_pre[ok_pre] - v1_pre[ok_pre])^2, na.rm = TRUE)
    ss_tot <- sum((v2_pre[ok_pre] - mean(v2_pre[ok_pre], na.rm = TRUE))^2, na.rm = TRUE)
    if (ss_tot > 0) 1 - (ss_res / ss_tot) else NA_real_
  } else NA_real_
  r2_post <- if (sum(ok_post) >= 5) {
    ss_res <- sum((v2_post[ok_post] - v1_post[ok_post])^2, na.rm = TRUE)
    ss_tot <- sum((v2_post[ok_post] - mean(v2_post[ok_post], na.rm = TRUE))^2, na.rm = TRUE)
    if (ss_tot > 0) 1 - (ss_res / ss_tot) else NA_real_
  } else NA_real_

  data.table(
    protein = p,
    cor_pre = cor_pre, cor_post = cor_post,
    mse_pre = mse_pre, mse_post = mse_post,
    mae_pre = mae_pre, mae_post = mae_post,
    r2_pre = r2_pre, r2_post = r2_post,
    mse_reduction = ifelse(!is.na(mse_pre) && !is.na(mse_post) && mse_pre > 0,
                          (mse_pre - mse_post) / mse_pre, NA_real_)
  )
}))
# Flag proteins where MSE did not decrease (or increased) or R² remains low
protein_concordance[, poorly_harmonised := (
  (!is.na(mse_reduction) & mse_reduction < 0) |  # MSE increased
  (!is.na(r2_post) & r2_post < 0.5)  # R² remains low
)]

n_poor <- sum(protein_concordance$poorly_harmonised, na.rm = TRUE)
log_info("  Poorly harmonised proteins (MSE increased or R² < 0.5): {n_poor}/{length(common_proteins)}")
if (n_poor > 0 && n_poor <= 20) {
  poor_prots <- protein_concordance[poorly_harmonised == TRUE][order(r2_post, na.last = TRUE)]
  log_info("  Worst proteins (by R²): {paste(head(poor_prots$protein, 10), collapse=', ')}")
} else if (n_poor > 20) {
  poor_prots <- protein_concordance[poorly_harmonised == TRUE][order(r2_post, na.last = TRUE)]
  log_info("  Worst 10 proteins (by R²): {paste(head(poor_prots$protein, 10), collapse=', ')}")
  log_info("  Median R² (poor proteins): {round(median(poor_prots$r2_post, na.rm=TRUE), 3)}")
}

# ==============================================================================
# STEP 6: KPI SUMMARY TABLE & MULTI-PANEL PDF DASHBOARD
# ==============================================================================
log_info("Step 6: Building KPI summary and dashboard")

# 6a. KPI summary table
kpi_dt <- data.table(
  metric = c(
    "Median bridge Euclidean distance",
    "Median bridge Mahalanobis distance",
    "Median bridge MAD-whitened distance",
    "Paired log-ratio median (MAD-whitened)",
    "Paired log-ratio 95% CI lower (MAD-whitened)",
    "Paired log-ratio 95% CI upper (MAD-whitened)",
    "CI excludes zero (MAD-whitened)",
    "Noise floor (within-batch replicate dist)",
    "Noise floor ratio (bridge_post / noise)",
    "ICC (median across pairs)",
    "CCC (median across pairs)",
    "Rank-1 ID rate",
    "Median rank of true pair",
    "Silhouette score (PC1-PC2)",
    "kBET (k-Nearest Neighbour Batch Effect Test) acceptance rate",
    "LISI (Local Inverse Simpson's Index) batch mixing score",
    "Sex prediction accuracy",
    "Batch variance fraction (median)",
    "Sex variance fraction (median)",
    "Flagged bridge outliers",
    "Poorly harmonised proteins (MSE↑ or R²<0.5)",
    "Wilcoxon signed-rank p-value"
  ),
  pre_harmonisation = c(
    round(median(d_eucl_pre), 4),
    round(median(d_maha_pre), 4),
    NA_real_,
    NA_real_,
    NA_real_,
    NA_real_,
    NA_real_,
    NA_real_,
    NA_real_,
    round(median(agree_pre$icc, na.rm = TRUE), 4),
    round(median(agree_pre$ccc, na.rm = TRUE), 4),
    round(id_pre$rank1_rate, 4),
    id_pre$median_rank,
    round(sil_pre, 4),
    round(kbet_pre, 4),
    round(lisi_pre, 4),
    round(sex_accuracy_pre, 4),
    ifelse(!is.null(var_decomp_pre), round(median(var_decomp_pre$batch_var_frac, na.rm = TRUE), 4), NA_real_),
    ifelse(!is.null(var_decomp_pre), round(median(var_decomp_pre$sex_var_frac, na.rm = TRUE), 4), NA_real_),
    NA_real_,
    NA_real_,
    NA_real_
  ),
  post_harmonisation = c(
    round(median(d_eucl_post), 4),
    round(median(d_maha_post), 4),
    round(median(d_mad_post), 4),
    round(median_log_ratio_mad, 4),
    round(boot_ci_mad[1], 4),
    round(boot_ci_mad[2], 4),
    as.numeric(boot_ci_mad[2] < 0),
    round(noise_floor_mad, 4),
    round(noise_floor_ratio, 4),
    round(median(agree_post$icc, na.rm = TRUE), 4),
    round(median(agree_post$ccc, na.rm = TRUE), 4),
    round(id_post$rank1_rate, 4),
    id_post$median_rank,
    round(sil_post, 4),
    round(kbet_post, 4),
    round(lisi_post, 4),
    round(sex_accuracy_post, 4),
    ifelse(!is.null(var_decomp_post), round(median(var_decomp_post$batch_var_frac, na.rm = TRUE), 4), NA_real_),
    ifelse(!is.null(var_decomp_post), round(median(var_decomp_post$sex_var_frac, na.rm = TRUE), 4), NA_real_),
    n_outliers,
    n_poor,
    round(wilcox_result$p.value, 6)
  ),
  target = c(
    "decrease",
    "decrease",
    "decrease",
    "negative",
    "negative",
    "< 0",
    "1 (TRUE)",
    "reference",
    "close to 1.0",
    "> 0.9",
    "> 0.9",
    "> 0.9",
    "1",
    "decrease toward 0",
    "> 0.8",
    "close to 2.0",
    "stable",
    "decrease",
    "stable (not drop >20%)",
    "0 ideal",
    "0 ideal",
    "< 0.05"
  )
)

# Determine pass/fail
kpi_dt[, pass := NA_character_]
kpi_dt[metric == "Median bridge Euclidean distance",
       pass := ifelse(post_harmonisation < pre_harmonisation, "PASS", "FAIL")]
kpi_dt[metric == "Median bridge Mahalanobis distance",
       pass := ifelse(post_harmonisation < pre_harmonisation, "PASS", "FAIL")]
kpi_dt[metric == "Median bridge MAD-whitened distance",
       pass := ifelse(post_harmonisation < pre_harmonisation, "PASS", "FAIL")]
kpi_dt[metric == "CI excludes zero (MAD-whitened)",
       pass := ifelse(post_harmonisation == 1, "PASS", "FAIL")]
kpi_dt[metric == "Noise floor ratio (bridge_post / noise)",
       pass := ifelse(!is.na(post_harmonisation) & post_harmonisation < 2.0, "PASS",
                      ifelse(!is.na(post_harmonisation) & post_harmonisation < 3.0, "WARN", "FAIL"))]
kpi_dt[metric == "ICC (median across pairs)",
       pass := ifelse(!is.na(post_harmonisation) & post_harmonisation > 0.9, "PASS",
                      ifelse(!is.na(post_harmonisation) & post_harmonisation > 0.7, "WARN", "FAIL"))]
kpi_dt[metric == "CCC (median across pairs)",
       pass := ifelse(!is.na(post_harmonisation) & post_harmonisation > 0.9, "PASS",
                      ifelse(!is.na(post_harmonisation) & post_harmonisation > 0.7, "WARN", "FAIL"))]
kpi_dt[metric == "Rank-1 ID rate",
       pass := ifelse(!is.na(post_harmonisation) & post_harmonisation > 0.9, "PASS",
                      ifelse(!is.na(post_harmonisation) & post_harmonisation > 0.7, "WARN", "FAIL"))]
kpi_dt[metric == "Silhouette score (PC1-PC2)",
       pass := ifelse(post_harmonisation < pre_harmonisation, "PASS", "FAIL")]
kpi_dt[metric == "kBET (k-Nearest Neighbour Batch Effect Test) acceptance rate",
       pass := ifelse(!is.na(post_harmonisation) & post_harmonisation > 0.8, "PASS",
                      ifelse(!is.na(post_harmonisation) & post_harmonisation > 0.5, "WARN", "FAIL"))]
kpi_dt[metric == "LISI (Local Inverse Simpson's Index) batch mixing score",
       pass := ifelse(!is.na(post_harmonisation) & post_harmonisation > 1.8, "PASS",
                      ifelse(!is.na(post_harmonisation) & post_harmonisation > 1.5, "WARN", "FAIL"))]
kpi_dt[metric == "Wilcoxon signed-rank p-value",
       pass := ifelse(!is.na(post_harmonisation) & post_harmonisation < 0.05, "PASS", "FAIL")]

# ==============================================================================
# SAVE OUTPUTS
# ==============================================================================
log_info("Saving output files")

output_batch_id <- batch1_id  # Primary batch for output naming

# KPI summary table
kpi_path <- get_output_path("07b", "kpi_summary", output_batch_id, "normalized", "tsv", config = config)
ensure_output_dir(kpi_path)
fwrite(kpi_dt, kpi_path, sep = "\t")
log_info("  Saved KPI summary: {kpi_path}")

# Bridge pair distances
pair_dist_path <- get_output_path("07b", "bridge_pair_distances", output_batch_id, "normalized", "tsv", config = config)
ensure_output_dir(pair_dist_path)
fwrite(bridge_pair_dt, pair_dist_path, sep = "\t")
log_info("  Saved bridge pair distances: {pair_dist_path}")

# Noise floor analysis
# Verify lengths match before creating data.table
source_vec <- c(
  rep("B1_within_batch", nrow(wb_rep_b1)),
  rep("B2_within_batch", nrow(wb_rep_b2))
)
expected_noise_length <- nrow(wb_rep_b1) + nrow(wb_rep_b2)
if (length(noise_floor_dists) != expected_noise_length) {
  log_warn("noise_floor_dists length mismatch: {length(noise_floor_dists)} vs expected {expected_noise_length}")
  # Use only the valid (non-NA) distances if there's a mismatch
  if (length(noise_floor_dists) < expected_noise_length) {
    # Pad with NAs if shorter
    noise_floor_dists <- c(noise_floor_dists, rep(NA_real_, expected_noise_length - length(noise_floor_dists)))
  } else {
    # Truncate if longer
    noise_floor_dists <- noise_floor_dists[seq_len(expected_noise_length)]
  }
}
if (length(source_vec) != length(noise_floor_dists)) {
  log_error("Source and distance vector length mismatch: source={length(source_vec)}, distance={length(noise_floor_dists)}")
  stop("Noise floor data.table length mismatch")
}

noise_dt <- data.table(
  source = source_vec,
  distance = noise_floor_dists
)
noise_path <- get_output_path("07b", "noise_floor_analysis", output_batch_id, "normalized", "tsv", config = config)
ensure_output_dir(noise_path)
fwrite(noise_dt, noise_path, sep = "\t")
log_info("  Saved noise floor analysis: {noise_path}")

# Protein concordance
conc_path <- get_output_path("07b", "protein_concordance", output_batch_id, "normalized", "tsv", config = config)
ensure_output_dir(conc_path)
fwrite(protein_concordance, conc_path, sep = "\t")
log_info("  Saved protein concordance: {conc_path}")

# Poorly harmonised proteins (subset with only flagged proteins)
poorly_harmonised_proteins <- protein_concordance[poorly_harmonised == TRUE]
if (nrow(poorly_harmonised_proteins) > 0) {
  # Sort by R² (worst first) for easy identification
  poorly_harmonised_proteins <- poorly_harmonised_proteins[order(r2_post, na.last = TRUE)]
  # Select relevant columns for the output
  poor_prot_output <- poorly_harmonised_proteins[, .(
    protein,
    cor_pre, cor_post,
    mse_pre, mse_post, mse_reduction,
    mae_pre, mae_post,
    r2_pre, r2_post
  )]
  poor_path <- get_output_path("07b", "poorly_harmonised_proteins", output_batch_id, "normalized", "tsv", config = config)
  ensure_output_dir(poor_path)
  fwrite(poor_prot_output, poor_path, sep = "\t")
  log_info("  Saved poorly harmonised proteins list: {poor_path} ({nrow(poor_prot_output)} proteins)")
} else {
  log_info("  No poorly harmonised proteins to save (all proteins passed harmonisation criteria)")
}

# Outlier flags
outlier_path <- get_output_path("07b", "outlier_flags", output_batch_id, "normalized", "tsv", config = config)
ensure_output_dir(outlier_path)
fwrite(bridge_pair_dt[, .(finngenid, sid_b1, sid_b2, eucl_post, maha_post, noise_floor_ratio,
                          outlier_maha, outlier_noise, outlier_flag, rank_post)],
       outlier_path, sep = "\t")
log_info("  Saved outlier flags: {outlier_path}")

# Variance decomposition
if (!is.null(var_decomp_pre) && !is.null(var_decomp_post)) {
  vd_combined <- merge(var_decomp_pre, var_decomp_post,
                       by = "protein", suffixes = c("_pre", "_post"))
  vd_path <- get_output_path("07b", "variance_decomposition", output_batch_id, "normalized", "tsv", config = config)
  ensure_output_dir(vd_path)
  fwrite(vd_combined, vd_path, sep = "\t")
  log_info("  Saved variance decomposition: {vd_path}")
}

# ==============================================================================
# MULTI-PANEL PDF DASHBOARD
# ==============================================================================
log_info("Generating multi-panel PDF dashboard")

dashboard_path <- get_output_path("07b", "kpi_dashboard", output_batch_id, "normalized", "pdf", config = config)
ensure_output_dir(dashboard_path)

# Custom colour palette for KPI dashboard
kpi_pal <- c("#2A363BFF", "#019875FF", "#99B898FF", "#FECEA8FF",
             "#FF847CFF", "#E84A5FFF", "#C0392BFF", "#96281BFF")
col_pre   <- kpi_pal[6]  # "#E84A5FFF" warm red-pink for pre
col_post  <- kpi_pal[2]  # "#019875FF" teal green for post
col_accent <- kpi_pal[8] # "#96281BFF" dark red for outliers/flagged

# ==============================================================================
# VALIDATION: Quick sanity check before plotting
# ==============================================================================
log_info("Generating plots ({n_bridge} bridge pairs, {nrow(bridge_pair_dt)} rows in bridge_pair_dt)")

# --- Plot 0: Bridge Pair Euclidean Distance Boxplot with Wilcoxon Test -------
# Paired boxplot comparing pre vs post bridge distances (from old implementation)
boxplot_df <- data.frame(
  distance = c(d_mad_pre, d_mad_post),
  stage = factor(rep(c("Pre-harmonisation", "Post-harmonisation"), each = n_bridge),
                 levels = c("Pre-harmonisation", "Post-harmonisation")),
  pair = rep(seq_len(n_bridge), 2)
)

p0_boxplot <- ggpaired(boxplot_df, x = "stage", y = "distance", id = "pair",
                       fill = "stage", palette = c(col_pre, col_post),
                       line.color = "grey70", line.size = 0.4,
                       point.size = 2, width = 0.5) +
  stat_compare_means(method = "wilcox.test", paired = TRUE,
                     label = "p.format", label.x = 1.5, label.y = max(c(d_mad_pre, d_mad_post)) * 1.05,
                     size = 4) +
  labs(title = "Bridge Pair Distance (MAD-Whitened, Fixed PC Basis)",
       subtitle = sprintf("Wilcoxon p = %s | n = %d pairs",
                          format.pval(wilcox_result$p.value, digits = 3), n_bridge),
       y = "MAD-Whitened Distance", x = NULL) +
  theme(legend.position = "none")

# --- Plot 1: Bridge Paired Distance "Dumbbell" Plot (Pre → Post) ------------
# Enhanced version: sorted by improvement ratio, highlight outliers
dumbbell_dt <- data.table(
  pair_id = seq_len(n_bridge),
  finngenid = common_bridge_fg,
  d_pre = d_mad_pre,
  d_post = d_mad_post,
  improvement_ratio = d_mad_post / d_mad_pre
)
# Sort by improvement ratio (best improvements first)
dumbbell_dt <- dumbbell_dt[order(improvement_ratio)]
dumbbell_dt[, pair_order := seq_len(.N)]
# Flag outliers (D_post > median + 3*MAD) — use column reference, not global vector
outlier_thresh_dumbbell <- median(dumbbell_dt$d_post) + 3 * mad(dumbbell_dt$d_post)
dumbbell_dt[, is_outlier := d_post > outlier_thresh_dumbbell]

# Convert to data.frame for ggplot (avoids data.table [.data.table auto-indexing clashes)
dumbbell_df <- as.data.frame(dumbbell_dt)

# Add truncated sample IDs for labelling
dumbbell_df$finngenid_short <- substr(dumbbell_df$finngenid, 1, 10)

p1_dumbbell <- ggplot(dumbbell_df, aes(x = pair_order)) +
  geom_segment(aes(x = pair_order, xend = pair_order, y = d_pre, yend = d_post,
                   colour = is_outlier), linewidth = 1.2, alpha = 0.7) +
  geom_point(aes(y = d_pre, fill = "Pre"), shape = 21, size = 3, alpha = 0.8) +
  geom_point(aes(y = d_post, fill = "Post"), shape = 21, size = 3, alpha = 0.8) +
  geom_text(aes(y = d_post, label = finngenid_short),
            angle = 45, hjust = 1, vjust = 1.5, size = 2, colour = "grey30") +
  scale_colour_manual(values = c("FALSE" = kpi_pal[1], "TRUE" = col_accent),
                      guide = "none") +
  scale_fill_manual(values = c("Pre" = col_pre, "Post" = col_post),
                    name = "Stage") +
  labs(title = "Bridge Paired Distance 'Dumbbell' Plot (MAD-Whitened)",
       subtitle = sprintf("Sorted by improvement ratio | Outliers (red): %d pairs with D_post > median+3*MAD",
                          sum(dumbbell_df$is_outlier)),
       x = "Bridge Pair (sorted by improvement)", y = "MAD-Whitened Distance") +
  theme(legend.position = "bottom", axis.text.x = element_blank(), axis.ticks.x = element_blank())

# --- Panel 2: Paired distance reduction (spaghetti + log-ratio histogram) ----
spaghetti_dt <- data.table(
  pair = rep(seq_len(n_bridge), 2),
  stage = rep(c("Pre", "Post"), each = n_bridge),
  distance = c(d_mad_pre, d_mad_post)  # Use MAD-whitened for consistency
)
spaghetti_dt[, stage := factor(stage, levels = c("Pre", "Post"))]
spaghetti_df <- as.data.frame(spaghetti_dt)

p2a <- ggplot(spaghetti_df, aes(x = stage, y = distance, group = pair)) +
  geom_line(alpha = 0.4, colour = "grey40") +
  geom_point(aes(colour = stage), size = 2) +
  scale_colour_manual(values = c(col_pre, col_post)) +
  labs(title = "Paired Distance Reduction (MAD-Whitened)",
       subtitle = "Each line = one bridge pair",
       y = "MAD-Whitened Distance", x = NULL) +
  theme(legend.position = "none")

# --- Plot 3: Log-Ratio Strip/Dot Plot (Enhanced) ------------------------------
# Use MAD-whitened distances for log-ratio (already computed above)

# Create histogram plot only (scatter plot removed)
p2b_hist_df <- data.frame(lr = log_ratio_mad, pair_id = seq_len(n_bridge))

p2b <- ggplot(p2b_hist_df, aes(x = lr)) +
  geom_histogram(bins = 15, fill = col_post, alpha = 0.7, colour = "black", linewidth = 0.3) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "red", linewidth = 0.7) +
  geom_vline(xintercept = median_log_ratio_mad, colour = "blue", linewidth = 0.9) +
  annotate("rect", xmin = boot_ci_mad[1], xmax = boot_ci_mad[2], ymin = -Inf, ymax = Inf,
           fill = "blue", alpha = 0.15) +
  labs(title = "Paired Log-Ratio: log(D_post / D_pre) [MAD-Whitened]",
       subtitle = sprintf("Median = %.3f, 95%% CI [%.3f, %.3f] | n = %d pairs",
                          median_log_ratio_mad, boot_ci_mad[1], boot_ci_mad[2], n_bridge),
       x = "log(D_post / D_pre)", y = "Frequency") +
  theme(legend.position = "none")

# --- Plot 2: Distribution Overlay of Bridge Distances (Pre vs Post) ----------
# ECDF overlay (recommended for n=34)
ecdf_dt <- rbind(
  data.table(stage = "Pre-harmonisation", distance = d_mad_pre),
  data.table(stage = "Post-harmonisation", distance = d_mad_post)
)
ecdf_dt[, stage := factor(stage, levels = c("Pre-harmonisation", "Post-harmonisation"))]
ecdf_df <- as.data.frame(ecdf_dt)

# Compute percentiles
p90_pre <- quantile(d_mad_pre, 0.90, na.rm = TRUE)
p95_pre <- quantile(d_mad_pre, 0.95, na.rm = TRUE)
p90_post <- quantile(d_mad_post, 0.90, na.rm = TRUE)
p95_post <- quantile(d_mad_post, 0.95, na.rm = TRUE)

p2_dist <- ggplot(ecdf_df, aes(x = distance, colour = stage)) +
  stat_ecdf(linewidth = 1.2, alpha = 0.8) +
  geom_vline(xintercept = median(d_mad_pre), linetype = "dashed", colour = col_pre, alpha = 0.6) +
  geom_vline(xintercept = median(d_mad_post), linetype = "dashed", colour = col_post, alpha = 0.6) +
  geom_vline(xintercept = p90_pre, linetype = "dotted", colour = col_pre, alpha = 0.5) +
  geom_vline(xintercept = p95_pre, linetype = "dotted", colour = col_pre, alpha = 0.5) +
  geom_vline(xintercept = p90_post, linetype = "dotted", colour = col_post, alpha = 0.5) +
  geom_vline(xintercept = p95_post, linetype = "dotted", colour = col_post, alpha = 0.5) +
  annotate("text", x = median(d_mad_pre), y = 0.5, label = "Pre median",
           hjust = -0.1, size = 3, colour = col_pre) +
  annotate("text", x = median(d_mad_post), y = 0.5, label = "Post median",
           hjust = -0.1, size = 3, colour = col_post) +
  scale_colour_manual(values = c(col_pre, col_post)) +
  labs(title = "Distribution Overlay: Bridge Distances (MAD-Whitened)",
       subtitle = sprintf("Pre: median=%.2f, p90=%.2f, p95=%.2f | Post: median=%.2f, p90=%.2f, p95=%.2f",
                          median(d_mad_pre), p90_pre, p95_pre, median(d_mad_post), p90_post, p95_post),
       x = "MAD-Whitened Distance", y = "Cumulative Probability") +
  theme(legend.position = "bottom", legend.title = element_blank())

# --- Panel 3: Noise floor comparison -----------------------------------------
# Use MAD-whitened distances for BOTH bridge pairs and noise floor (consistent units)
nf_plot_dt <- rbind(
  data.table(group = "Bridge pairs (post)", distance = d_mad_post),
  data.table(group = "Within-batch replicates", distance = noise_floor_dists_mad[!is.na(noise_floor_dists_mad)])
)

nf_plot_df <- as.data.frame(nf_plot_dt)

p3 <- ggplot(nf_plot_df, aes(x = distance, fill = group)) +
  geom_density(alpha = 0.5) +
  geom_vline(xintercept = noise_floor_mad, linetype = "dashed", colour = "grey30") +
  annotate("text", x = noise_floor_mad, y = Inf, label = sprintf("Noise floor = %.2f", noise_floor_mad),
           vjust = 2, hjust = -0.1, size = 3) +
  labs(title = "Bridge Distance vs Within-Batch Noise Floor",
       subtitle = sprintf("Ratio = %.2f (target: ~1.0)", noise_floor_ratio),
       x = "MAD-Whitened Distance (Fixed PC Basis)", y = "Density") +
  scale_fill_manual(values = c(col_post, "grey60")) +
  theme(legend.position = "bottom", legend.title = element_blank())

# --- Panel 4: Rank-1 identification barplot ----------------------------------
rank_dt <- data.table(
  stage = c("Pre-harmonisation", "Post-harmonisation"),
  rate  = c(id_pre$rank1_rate, id_post$rank1_rate),
  hits  = c(id_pre$rank1_hits, id_post$rank1_hits)
)
rank_dt[, stage := factor(stage, levels = c("Pre-harmonisation", "Post-harmonisation"))]
rank_df <- as.data.frame(rank_dt)

p4 <- ggplot(rank_df, aes(x = stage, y = rate, fill = stage)) +
  geom_col(alpha = 0.8) +
  geom_text(aes(label = sprintf("%d/%d", hits, n_bridge)), vjust = -0.3, size = 3.5) +
  scale_fill_manual(values = c(col_pre, col_post)) +
  scale_y_continuous(limits = c(0, 1.1), labels = scales::percent_format()) +
  labs(title = "Rank-1 Identification Rate",
       subtitle = "Fraction of bridge B1 samples whose NN in B2 is their true pair",
       y = "Rank-1 Hit Rate", x = NULL) +
  theme(legend.position = "none")

# --- Plot 4: PC Scatter with Bridge Connectors (Pre and Post Panels) ----------
# Enhanced: Add connector lines for bridge pairs
pca_plot_data <- function(pc_b1, pc_b2, stage, bid1, bid2, bridge_sids_b1, bridge_sids_b2) {
  # Use data.frame to avoid data.table [.data.table issues in ggplot
  df <- rbind(
    data.frame(PC1 = pc_b1[, 1], PC2 = pc_b1[, 2], batch = bid1, sample = rownames(pc_b1)),
    data.frame(PC1 = pc_b2[, 1], PC2 = pc_b2[, 2], batch = bid2, sample = rownames(pc_b2))
  )
  df$bridge <- df$sample %in% c(bridge_sids_b1, bridge_sids_b2)
  df$stage <- stage
  df
}

pca_dt <- rbind(
  pca_plot_data(pc_pre$b1,  pc_pre$b2,  "Pre-harmonisation",  batch1_id, batch2_id, sids_b1, sids_b2),
  pca_plot_data(pc_post$b1, pc_post$b2, "Post-harmonisation", batch1_id, batch2_id, sids_b1, sids_b2)
)
pca_dt$stage <- factor(pca_dt$stage, levels = c("Pre-harmonisation", "Post-harmonisation"))

# Create connector lines data for bridge pairs
create_bridge_connectors <- function(pc_b1, pc_b2, sids_b1, sids_b2, stage_name) {
  # Use data.frame (not data.table) to avoid [.data.table issues in ggplot
  rows <- lapply(seq_along(sids_b1), function(i) {
    if (sids_b1[i] %in% rownames(pc_b1) && sids_b2[i] %in% rownames(pc_b2)) {
      data.frame(
        x = pc_b1[sids_b1[i], 1],
        y = pc_b1[sids_b1[i], 2],
        xend = pc_b2[sids_b2[i], 1],
        yend = pc_b2[sids_b2[i], 2],
        pair_id = i,
        stage = stage_name,
        stringsAsFactors = FALSE
      )
    } else {
      NULL
    }
  })
  do.call(rbind, Filter(Negate(is.null), rows))
}

bridge_connectors_pre  <- create_bridge_connectors(pc_pre$b1,  pc_pre$b2,  sids_b1, sids_b2, "Pre-harmonisation")
bridge_connectors_post <- create_bridge_connectors(pc_post$b1, pc_post$b2, sids_b1, sids_b2, "Post-harmonisation")
bridge_connectors <- rbind(bridge_connectors_pre, bridge_connectors_post)
bridge_connectors$stage <- factor(bridge_connectors$stage, levels = c("Pre-harmonisation", "Post-harmonisation"))

# pca_dt and bridge_connectors are already data.frames (not data.tables)
pca_df <- pca_dt
bridge_connectors_df <- bridge_connectors

# Downsample non-bridge samples if too many (for performance)
n_non_bridge <- sum(pca_df$bridge == FALSE)
if (n_non_bridge > 2000) {
  set.seed(42)
  non_bridge_ids <- pca_df$sample[pca_df$bridge == FALSE]
  sample_idx <- sample.int(length(non_bridge_ids), 2000)
  pca_df_plot <- rbind(
    pca_df[pca_df$bridge == TRUE, ],
    pca_df[pca_df$sample %in% non_bridge_ids[sample_idx], ]
  )
} else {
  pca_df_plot <- pca_df
}

p5 <- ggplot(pca_df_plot[pca_df_plot$bridge == FALSE, ], aes(x = PC1, y = PC2, colour = batch)) +
  geom_point(alpha = 0.2, size = 0.5) +
  geom_segment(data = bridge_connectors_df, aes(x = x, y = y, xend = xend, yend = yend),
               inherit.aes = FALSE,
               colour = "grey50", alpha = 0.4, linewidth = 0.5) +
  geom_point(data = pca_df_plot[pca_df_plot$bridge == TRUE, ], aes(shape = bridge), size = 2.5, alpha = 0.9) +
  facet_wrap(~stage) +
  scale_colour_manual(values = c("#2c7bb6", "#d7191c"), name = "Batch", labels = c(batch1_id, batch2_id)) +
  scale_shape_manual(values = c("TRUE" = 17)) +
  labs(title = "PCA (Fixed Basis): Before vs After Harmonisation",
       subtitle = "Triangles = bridge samples | Lines = bridge pair connectors") +
  guides(shape = "none") +
  theme(legend.position = "bottom")

# --- Panel 6: Batch separability metrics comparison (will be replaced by p6_enhanced) ---
# This is kept for backward compatibility but p6_enhanced (with AUC) is used in PDF

# --- Panel 7: Variance decomposition (ggpubr boxplot) -------------------------
# Use full per-protein distributions instead of medians for richer visualization
p7 <- NULL
if (!is.null(var_decomp_pre) && !is.null(var_decomp_post)) {
  # Build long format with protein names
  vd_long <- rbind(
    data.frame(
      stage = "Pre",
      component = "Batch",
      fraction = var_decomp_pre$batch_var_frac,
      protein = var_decomp_pre$protein
    ),
    data.frame(
      stage = "Pre",
      component = "Sex",
      fraction = var_decomp_pre$sex_var_frac,
      protein = var_decomp_pre$protein
    ),
    data.frame(
      stage = "Post",
      component = "Batch",
      fraction = var_decomp_post$batch_var_frac,
      protein = var_decomp_post$protein
    ),
    data.frame(
      stage = "Post",
      component = "Sex",
      fraction = var_decomp_post$sex_var_frac,
      protein = var_decomp_post$protein
    )
  )
  vd_long$stage <- factor(vd_long$stage, levels = c("Pre", "Post"))
  vd_long$component <- factor(vd_long$component, levels = c("Batch", "Sex"))

  # Identify proteins with >30% variance explained (changed from >50%)
  vd_long$high_var <- vd_long$fraction > 0.3
  high_var_proteins <- unique(vd_long$protein[vd_long$high_var == TRUE])

  # Split into low-variance (for jittered points) and high-variance (for triangles only)
  vd_low_var <- subset(vd_long, !high_var)  # fraction <= 0.3
  vd_high_var <- subset(vd_long, high_var)   # fraction > 0.3

  # Compute median percentages for labels on top of boxes
  vd_summary <- aggregate(fraction ~ stage + component, data = vd_long,
                          FUN = function(x) median(x, na.rm = TRUE))
  names(vd_summary)[names(vd_summary) == "fraction"] <- "median_pct"

  # Create high-variance subset with aligned positions (no jitter)
  # Since facet_wrap(~stage) separates stages, triangles only need the component x-position
  if (nrow(vd_high_var) > 0) {
    vd_high_var$aligned_x <- as.numeric(vd_high_var$component)
    vd_high_var$aligned_y <- vd_high_var$fraction
  }

  # Compute subtitle text before passing to labs() (avoid nested if() in function arguments)
  vd_subtitle_base <- "Batch variance should decrease; sex variance should remain stable | Red triangles = proteins with >30% variance explained\n"
  if(length(high_var_proteins) > 0) {
    high_var_text <- paste(head(high_var_proteins, 10), collapse = ", ")
    if(length(high_var_proteins) > 10) {
      high_var_text <- paste0(high_var_text, " ...")
    }
    vd_subtitle <- paste0(vd_subtitle_base,
                          "High-variance proteins (", length(high_var_proteins), "): ", high_var_text)
  } else {
    vd_subtitle <- paste0(vd_subtitle_base, "No high-variance proteins")
  }

  # Use only low-variance data for the main boxplot with jitter
  # This ensures high-variance proteins are NOT shown as jittered points
  p7 <- ggboxplot(vd_low_var, x = "component", y = "fraction",
                   fill = "stage", palette = c(col_pre, col_post),
                   add = "jitter", add.params = list(size = 0.8, alpha = 0.3),
                   outlier.shape = NA) +
    # Add percentage labels on top of boxes
    geom_text(data = vd_summary,
              aes(x = component, y = 1.05, label = sprintf("%.1f%%", median_pct * 100)),
              position = position_dodge(width = 0.75), size = 3.5, fontface = "bold") +
    # Highlight high-variance proteins (>30%) with triangles (shape 17)
    # Triangles are aligned with the boxplot groups (no jitter) for clarity
    geom_point(data = vd_high_var,
               aes(x = aligned_x, y = aligned_y),
               shape = 17, size = 3.5, colour = "#FF0000", alpha = 0.9) +
    # Add protein name labels for high-variance proteins using ggrepel for clean positioning
    {if(nrow(vd_high_var) > 0) {
      if (requireNamespace("ggrepel", quietly = TRUE)) {
        ggrepel::geom_text_repel(data = vd_high_var,
                                  aes(x = aligned_x, y = aligned_y, label = protein),
                                  size = 2.5, colour = "grey30", fontface = "bold",
                                  min.segment.length = 0, box.padding = 0.3, point.padding = 0.2,
                                  max.overlaps = Inf, force = 2)
      } else {
        geom_text(data = vd_high_var,
                  aes(x = aligned_x, y = aligned_y + 0.05, label = protein),
                  size = 2.5, hjust = 0, vjust = 0, angle = 45, colour = "grey30",
                  fontface = "bold", nudge_y = 0.02)
      }
    } else {
      geom_blank()
    }} +
    scale_y_continuous(labels = scales::percent_format(), limits = c(0, 1.2)) +
    facet_wrap(~stage, ncol = 2) +
    labs(title = "Variance Decomposition (per-protein distributions)",
         subtitle = vd_subtitle,
         y = "Fraction of Variance", x = NULL, fill = "Stage") +
    theme(legend.position = "bottom",
          plot.subtitle = element_text(size = 9, lineheight = 1.2))
}

# --- Plot 5: Nearest-Neighbor Pairing Heatmap (Swap Detector) ----------------
# Compute distance matrix: for each bridge A_i, distance to all B_j
compute_nn_heatmap <- function(pc_b1, pc_b2, sids_b1, sids_b2) {
  # Get MAD scaling factors
  combined_pc <- rbind(pc_b1, pc_b2)
  mad_per_pc <- apply(combined_pc, 2, function(x) mad(x, na.rm = TRUE))
  mad_per_pc[mad_per_pc == 0] <- 1

  # Compute distance matrix
  dist_mat <- matrix(NA_real_, nrow = length(sids_b1), ncol = length(sids_b2))
  for (i in seq_along(sids_b1)) {
    for (j in seq_along(sids_b2)) {
      if (sids_b1[i] %in% rownames(pc_b1) && sids_b2[j] %in% rownames(pc_b2)) {
        diff_vec <- pc_b1[sids_b1[i], , drop = FALSE] - pc_b2[sids_b2[j], , drop = FALSE]
        scaled_diff <- diff_vec / mad_per_pc
        dist_mat[i, j] <- sqrt(sum(scaled_diff^2, na.rm = TRUE))
      }
    }
  }
  dist_mat
}

nn_dist_mat_post <- compute_nn_heatmap(pc_post$b1, pc_post$b2, sids_b1, sids_b2)

# Verify matrix dimensions match
if (nrow(nn_dist_mat_post) != length(sids_b1) || ncol(nn_dist_mat_post) != length(sids_b2)) {
  log_error("Distance matrix dimension mismatch: matrix={nrow(nn_dist_mat_post)}x{ncol(nn_dist_mat_post)}, expected={length(sids_b1)}x{length(sids_b2)}")
  stop("Distance matrix dimension mismatch")
}

# Label rows/cols with truncated FINNGENID for readability
fg_short <- substr(common_bridge_fg, 1, 10)
colnames(nn_dist_mat_post) <- fg_short
rownames(nn_dist_mat_post) <- fg_short

# Convert to long format for ggplot
nn_dist_df <- as.data.frame(nn_dist_mat_post)
nn_dist_df$B1_idx <- seq_len(nrow(nn_dist_df))
nn_heatmap_dt <- melt(as.data.table(nn_dist_df),
                      id.vars = "B1_idx", variable.name = "B2_label", value.name = "distance")
nn_heatmap_dt[, B2_idx := match(B2_label, fg_short)]
nn_heatmap_dt[, is_diagonal := B1_idx == B2_idx]
nn_heatmap_dt[, B1_label := fg_short[B1_idx]]

# Identify outlier bridges for red border highlighting
outlier_fg_short <- substr(bridge_pair_dt$finngenid[bridge_pair_dt$outlier_flag == TRUE], 1, 10)
nn_heatmap_dt[, is_outlier_row := B1_label %in% outlier_fg_short]
nn_heatmap_dt[, is_outlier_col := B2_label %in% outlier_fg_short]
nn_heatmap_dt[, is_outlier_cell := is_outlier_row | is_outlier_col]

# Create ggplot heatmap with spectral color brewer, sample IDs, and red borders for outliers
midpoint_val <- median(nn_dist_mat_post, na.rm = TRUE)
if (is.na(midpoint_val)) midpoint_val <- 0

nn_heatmap_df <- as.data.frame(nn_heatmap_dt)

# Add outlier status column for legend
nn_heatmap_df$outlier_status <- "Normal"
nn_heatmap_df$outlier_status[nn_heatmap_df$is_diagonal == TRUE] <- "True Pair"
nn_heatmap_df$outlier_status[nn_heatmap_df$is_outlier_cell == TRUE & nn_heatmap_df$is_diagonal != TRUE] <- "Outlier Bridge"

# Add flagged status for legend
nn_heatmap_df$flagged_status <- "Normal"
nn_heatmap_df$flagged_status[nn_heatmap_df$is_diagonal == TRUE] <- "True Pair (Diagonal)"
nn_heatmap_df$flagged_status[nn_heatmap_df$is_outlier_cell == TRUE & nn_heatmap_df$is_diagonal != TRUE] <- "Flagged Outlier"

p5_heatmap <- ggplot(nn_heatmap_df, aes(x = B2_idx, y = B1_idx, fill = distance)) +
  geom_tile(colour = "black", linewidth = 0.1) +  # Thin black borders for all cells
  geom_tile(data = nn_heatmap_df[nn_heatmap_df$is_diagonal == TRUE, ],
            colour = "yellow", linewidth = 0.3) +  # Yellow border for diagonal (true pairs)
  geom_tile(data = nn_heatmap_df[nn_heatmap_df$is_outlier_cell == TRUE & nn_heatmap_df$is_diagonal != TRUE, ],
            colour = col_accent, linewidth = 0.5) +  # Red border for outlier bridges
  geom_text(aes(label = round(distance, 2)), size = 2.2, colour = "black", fontface = "bold") +
  scale_fill_distiller(palette = "RdBu", name = "Distance\n(MAD-Whitened)", direction = -1,
                       guide = guide_colorbar(title.position = "top", title.hjust = 0.5)) +
  scale_x_continuous(breaks = seq_along(fg_short), labels = fg_short, expand = c(0, 0)) +
  scale_y_continuous(breaks = seq_along(fg_short), labels = fg_short, expand = c(0, 0)) +
  labs(title = "Nearest-Neighbor Pairing Heatmap (Swap Detector)",
       subtitle = paste0("Pairing calculation: For each bridge sample in Batch 1 (rows), compute MAD-whitened Euclidean distance ",
                        "to all bridge samples in Batch 2 (columns) in post-harmonisation PC space. ",
                        "True pairs form the diagonal (i==j). Yellow border = true pairs | Red border = flagged outliers (",
                        length(outlier_fg_short), " flagged) | Dark = close, Light = far"),
       x = "Batch 2 Bridge Sample (FINNGENID)", y = "Batch 1 Bridge Sample (FINNGENID)") +
  guides(fill = guide_colorbar(title.position = "top", title.hjust = 0.5)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
        axis.text.y = element_text(size = 6),
        legend.position = "right",
        plot.subtitle = element_text(size = 9, lineheight = 1.2))

# --- Plot 5b: Protein Correlation Heatmap for Bridge Samples (Pre & Post) -----
# For each bridge pair, compute Pearson correlation of protein profiles (B1 vs B2)
log_info("  Computing protein correlation for bridge samples (pre & post)")

compute_bridge_protein_corr_matrix <- function(mat_b1, mat_b2, sids_b1_all, sids_b2_all, proteins) {
  n <- length(sids_b1_all)
  corr_mat <- matrix(NA_real_, nrow = n, ncol = n)
  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      if (sids_b1_all[i] %in% rownames(mat_b1) && sids_b2_all[j] %in% rownames(mat_b2)) {
        x <- as.numeric(mat_b1[sids_b1_all[i], proteins])
        y <- as.numeric(mat_b2[sids_b2_all[j], proteins])
        ok <- !is.na(x) & !is.na(y)
        if (sum(ok) > 10) {
          corr_mat[i, j] <- cor(x[ok], y[ok], method = "pearson")
        }
      }
    }
  }
  rownames(corr_mat) <- substr(common_bridge_fg, 1, 10)
  colnames(corr_mat) <- substr(common_bridge_fg, 1, 10)
  corr_mat
}

bridge_corr_pre  <- compute_bridge_protein_corr_matrix(pre_b1, pre_b2, sids_b1, sids_b2, common_proteins)
bridge_corr_post <- compute_bridge_protein_corr_matrix(post_b1, post_b2, sids_b1, sids_b2, common_proteins)

# Identify bridges with low diagonal (true pair) correlation (for reference only)
diag_corr_pre  <- diag(bridge_corr_pre)
diag_corr_post <- diag(bridge_corr_post)
low_corr_thresh <- 0.8
low_corr_fg <- substr(common_bridge_fg[diag_corr_post < low_corr_thresh], 1, 10)

# Get flagged samples based on MAD-whitened outlier detection (primary flagging method)
outlier_fg_short_mad <- substr(bridge_pair_dt$finngenid[bridge_pair_dt$outlier_flag == TRUE], 1, 10)

# Convert protein correlation matrices to long format for ggplot (separate plots)
convert_corr_to_long <- function(corr_mat, fg_labels, stage_label, flagged_fg) {
  corr_df <- as.data.frame(corr_mat)
  corr_df$B1_idx <- seq_len(nrow(corr_df))
  corr_long <- melt(as.data.table(corr_df),
                    id.vars = "B1_idx", variable.name = "B2_label", value.name = "correlation")
  corr_long[, B2_idx := match(B2_label, fg_labels)]
  corr_long[, is_diagonal := B1_idx == B2_idx]
  corr_long[, B1_label := fg_labels[B1_idx]]
  corr_long[, stage := stage_label]
  # Add flagged status based on MAD-whitened outlier detection (not low correlation)
  corr_long[, is_flagged := B1_label %in% flagged_fg | B2_label %in% flagged_fg]
  as.data.frame(corr_long)
}

corr_pre_long <- convert_corr_to_long(bridge_corr_pre, fg_short, "Pre-harmonisation", outlier_fg_short_mad)
corr_post_long <- convert_corr_to_long(bridge_corr_post, fg_short, "Post-harmonisation", outlier_fg_short_mad)

# Create separate ggplot heatmaps with BrBG palette, theme_bw, and side color labels for flagged samples
# Add side color indicator for flagged samples
corr_pre_long$flagged_colour <- ifelse(corr_pre_long$is_flagged, col_accent, "transparent")
corr_post_long$flagged_colour <- ifelse(corr_post_long$is_flagged, col_accent, "transparent")

# Compute caption text for both plots (avoid ifelse in labs() which causes "In index: 1." error)
caption_text <- if(length(outlier_fg_short_mad) > 0) {
  sprintf("Flagged bridges: %s", paste(head(outlier_fg_short_mad, 10), collapse = ", "))
} else {
  "No flagged bridges"
}

p5b_corr_pre <- ggplot(corr_pre_long, aes(x = B2_idx, y = B1_idx, fill = correlation)) +
  geom_tile(colour = "black", linewidth = 0.1) +  # Thin black borders for all cells
  geom_tile(data = corr_pre_long[corr_pre_long$is_diagonal == TRUE, ],
            colour = "yellow", linewidth = 0.3) +  # Yellow border for diagonal
  geom_tile(data = corr_pre_long[corr_pre_long$is_flagged == TRUE & corr_pre_long$is_diagonal != TRUE, ],
            colour = col_accent, linewidth = 0.4) +  # Red border for flagged bridges
  geom_text(aes(label = round(correlation, 2)),
            size = 2.8, colour = "black", fontface = "bold") +
  # Add side color bar for flagged samples using geom_rect
  {flagged_rows_pre <- unique(corr_pre_long$B1_idx[corr_pre_long$is_flagged == TRUE & !is.na(corr_pre_long$is_flagged)])
  if(length(flagged_rows_pre) > 0) {
    side_bar_data_pre <- data.frame(
      xmin = rep(-0.6, length(flagged_rows_pre)),
      xmax = rep(-0.3, length(flagged_rows_pre)),
      ymin = flagged_rows_pre - 0.5,
      ymax = flagged_rows_pre + 0.5
    )
    geom_rect(data = side_bar_data_pre,
              aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
              fill = col_accent, colour = NA, inherit.aes = FALSE)
  } else {
    geom_blank()
  }} +
  scale_fill_distiller(palette = "BrBG", name = "Correlation\n(Pearson r)", direction = 1, limits = c(-1, 1),
                       guide = guide_colorbar(title.position = "top", title.hjust = 0.5),
                       na.value = "transparent") +
  scale_x_continuous(breaks = seq_along(fg_short), labels = fg_short, expand = c(0.1, 0)) +
  scale_y_continuous(breaks = seq_along(fg_short), labels = fg_short, expand = c(0, 0)) +
  labs(title = "Pre-Harmonisation: Bridge Sample Protein Correlation (B1 vs B2)",
       subtitle = sprintf("Yellow border = true pairs (diagonal) | Red side bar = flagged bridges (MAD-whitened outliers, n=%d)",
                          length(outlier_fg_short_mad)),
       x = "Batch 2 Bridge Sample (FINNGENID)", y = "Batch 1 Bridge Sample (FINNGENID)",
       caption = caption_text) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
        axis.text.y = element_text(size = 6),
        legend.position = "right",
        plot.caption = element_text(size = 7, hjust = 0),
        panel.grid = element_blank())

p5b_corr_post <- ggplot(corr_post_long, aes(x = B2_idx, y = B1_idx, fill = correlation)) +
  geom_tile(colour = "black", linewidth = 0.1) +  # Thin black borders for all cells
  geom_tile(data = corr_post_long[corr_post_long$is_diagonal == TRUE, ],
            colour = "yellow", linewidth = 0.3) +  # Yellow border for diagonal
  geom_tile(data = corr_post_long[corr_post_long$is_flagged == TRUE & corr_post_long$is_diagonal != TRUE, ],
            colour = col_accent, linewidth = 0.4) +  # Red border for flagged bridges
  geom_text(aes(label = round(correlation, 2)),
            size = 2.8, colour = "black", fontface = "bold") +
  # Add side color bar for flagged samples using geom_rect
  {flagged_rows_post <- unique(corr_post_long$B1_idx[corr_post_long$is_flagged == TRUE & !is.na(corr_post_long$is_flagged)])
  if(length(flagged_rows_post) > 0) {
    side_bar_data_post <- data.frame(
      xmin = rep(-0.6, length(flagged_rows_post)),
      xmax = rep(-0.3, length(flagged_rows_post)),
      ymin = flagged_rows_post - 0.5,
      ymax = flagged_rows_post + 0.5
    )
    geom_rect(data = side_bar_data_post,
              aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
              fill = col_accent, colour = NA, inherit.aes = FALSE)
  } else {
    geom_blank()
  }} +
  scale_fill_distiller(palette = "BrBG", name = "Correlation\n(Pearson r)", direction = 1, limits = c(-1, 1),
                       guide = guide_colorbar(title.position = "top", title.hjust = 0.5),
                       na.value = "transparent") +
  scale_x_continuous(breaks = seq_along(fg_short), labels = fg_short, expand = c(0.1, 0)) +
  scale_y_continuous(breaks = seq_along(fg_short), labels = fg_short, expand = c(0, 0)) +
  labs(title = "Post-Harmonisation: Bridge Sample Protein Correlation (B1 vs B2)",
       subtitle = sprintf("Yellow border = true pairs (diagonal) | Red side bar = flagged bridges (MAD-whitened outliers, n=%d)",
                          length(outlier_fg_short_mad)),
       x = "Batch 2 Bridge Sample (FINNGENID)", y = "Batch 1 Bridge Sample (FINNGENID)",
       caption = caption_text) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
        axis.text.y = element_text(size = 6),
        legend.position = "right",
        plot.caption = element_text(size = 7, hjust = 0),
        panel.grid = element_blank())

log_info("  Bridge protein correlation: pre diagonal median={round(median(diag_corr_pre, na.rm = TRUE), 3)}, post diagonal median={round(median(diag_corr_post, na.rm = TRUE), 3)}")
log_info("  Bridges flagged by MAD-whitened outlier detection: {length(outlier_fg_short_mad)}/{n_bridge}")

# --- Plot 6: Pairing Margin Plot (Confidence of Correct Match) ---------------
# For each bridge i: margin = d_true - d_best_wrong
compute_pairing_margins <- function(dist_mat) {
  margins <- numeric(nrow(dist_mat))
  for (i in seq_len(nrow(dist_mat))) {
    d_true <- dist_mat[i, i]
    d_others <- dist_mat[i, -i, drop = FALSE]
    d_best_wrong <- min(d_others, na.rm = TRUE)
    margins[i] <- d_true - d_best_wrong
  }
  margins
}

pairing_margins <- compute_pairing_margins(nn_dist_mat_post)

margin_df <- data.frame(
  pair_id = seq_len(n_bridge),
  finngenid = common_bridge_fg,
  margin = pairing_margins,
  is_negative = pairing_margins < 0,
  plot_type = "Confidence of Correct Match"
)

# Create margin plot data
margin_plot_df <- margin_df
margin_plot_df$plot_type <- "Confidence of Correct Match"

# Create margin plot (kept separate for backward compatibility if needed)
p6_margin_plot <- ggplot(margin_plot_df, aes(x = reorder(finngenid, margin), y = margin, fill = is_negative)) +
  geom_col(alpha = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "red", linewidth = 0.7) +
  scale_fill_manual(values = c("FALSE" = col_post, "TRUE" = col_accent),
                    labels = c("Positive margin", "Negative margin (swap candidate)"),
                    name = "Margin") +
  labs(title = "Pairing Margin Plot (Confidence of Correct Match)",
       subtitle = sprintf("Margin = d_true - d_best_wrong | Negative margins: %d/%d (swap candidates)",
                          sum(margin_df$is_negative), n_bridge),
       x = "Bridge Pair (FINNGENID)", y = "Margin (Distance)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
        legend.position = "bottom")

# --- Plot 6: Enhanced Rank-1 Identification Rate with Distance Scatters (3-panel faceted) ---
# Define new color palette for consistent 3-color scheme
col_normal <- "#00AFBB"  # Normal/Pre
col_intermediate <- "#E7B800"  # Post/Intermediate
col_flagged <- "#FC4E07"  # Flagged

# Create rank-1 identification data
rank1_df <- data.frame(
  pair_id = seq_len(n_bridge),
  finngenid = common_bridge_fg,
  rank_pre = id_pre$ranks,
  rank_post = id_post$ranks,
  is_rank1_pre = id_pre$ranks == 1,
  is_rank1_post = id_post$ranks == 1,
  outlier_flag = bridge_pair_dt$outlier_flag
)
rank1_long <- rbind(
  data.frame(pair_id = rank1_df$pair_id, finngenid = rank1_df$finngenid,
             value = rank1_df$rank_pre, is_rank1 = rank1_df$is_rank1_pre,
             stage = "Pre-harmonisation", plot_type = "Rank-1 Identification Rate",
             outlier_flag = rank1_df$outlier_flag),
  data.frame(pair_id = rank1_df$pair_id, finngenid = rank1_df$finngenid,
             value = rank1_df$rank_post, is_rank1 = rank1_df$is_rank1_post,
             stage = "Post-harmonisation", plot_type = "Rank-1 Identification Rate",
             outlier_flag = rank1_df$outlier_flag)
)
rank1_long$stage <- factor(rank1_long$stage, levels = c("Pre-harmonisation", "Post-harmonisation"))

# Create Mahalanobis distance scatter data
maha_scatter_df <- data.frame(
  pair_id = seq_len(n_bridge),
  finngenid = common_bridge_fg,
  x = d_maha_pre,
  y = d_maha_post,
  outlier_flag = bridge_pair_dt$outlier_flag,
  plot_type = "Mahalanobis Distance"
)

# Create MAD-whitened distance scatter data
mad_scatter_df <- data.frame(
  pair_id = seq_len(n_bridge),
  finngenid = common_bridge_fg,
  x = d_mad_pre,
  y = d_mad_post,
  outlier_flag = bridge_pair_dt$outlier_flag,
  plot_type = "MAD-Whitened Distance"
)

# Create three separate plots, then combine with ggarrange for faceted-like appearance
# Panel 1: Rank-1 bar plot
p6_rank1_plot <- ggplot(rank1_long, aes(x = reorder(finngenid, value), y = value, fill = stage)) +
  geom_col(alpha = 0.8, position = "dodge") +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "green", linewidth = 0.7) +
  geom_point(data = rank1_long[rank1_long$is_rank1 == TRUE, ],
             aes(shape = "Rank-1"), size = 3, colour = "green", position = position_dodge(width = 0.9)) +
  scale_fill_manual(values = c(col_normal, col_intermediate), name = "Stage") +
  scale_shape_manual(name = "", values = c("Rank-1" = 17)) +
  labs(title = "Rank-1 Identification Rate (per Bridge Pair)",
       subtitle = sprintf("Pre: %d/%d rank-1 hits (%.1f%%) | Post: %d/%d rank-1 hits (%.1f%%)",
                          id_pre$rank1_hits, n_bridge, id_pre$rank1_rate*100,
                          id_post$rank1_hits, n_bridge, id_post$rank1_rate*100),
       x = "Bridge Pair (FINNGENID)", y = "Rank of True Pair") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
        legend.position = "bottom",
        plot.title = element_text(size = 12, face = "bold"),
        plot.subtitle = element_text(size = 9))

# Panel 2: Mahalanobis distance scatter
p6_maha_scatter <- ggplot(maha_scatter_df, aes(x = x, y = y, colour = outlier_flag)) +
  geom_point(size = 2.5, alpha = 0.7) +
  geom_smooth(aes(x = x, y = y), method = "lm", se = TRUE, colour = "black", linewidth = 0.8, inherit.aes = FALSE) +
  scale_colour_manual(values = c("FALSE" = col_normal, "TRUE" = col_flagged),
                       labels = c("Normal", "Flagged"),
                       name = "Status") +
  labs(title = "Mahalanobis Distance: Pre vs Post",
       subtitle = sprintf("Flagged samples: %d/%d", n_outliers, n_bridge),
       x = "Pre-harmonisation Mahalanobis Distance",
       y = "Post-harmonisation Mahalanobis Distance") +
  theme(legend.position = "bottom",
        plot.title = element_text(size = 12, face = "bold"),
        plot.subtitle = element_text(size = 9))

# Panel 3: MAD-whitened distance scatter
p6_mad_scatter <- ggplot(mad_scatter_df, aes(x = x, y = y, colour = outlier_flag)) +
  geom_point(size = 2.5, alpha = 0.7) +
  geom_smooth(aes(x = x, y = y), method = "lm", se = TRUE, colour = "black", linewidth = 0.8, inherit.aes = FALSE) +
  scale_colour_manual(values = c("FALSE" = col_normal, "TRUE" = col_flagged),
                       labels = c("Normal", "Flagged"),
                       name = "Status") +
  labs(title = "MAD-Whitened Distance: Pre vs Post",
       subtitle = sprintf("Flagged samples: %d/%d", n_outliers, n_bridge),
       x = "Pre-harmonisation MAD-Whitened Distance",
       y = "Post-harmonisation MAD-Whitened Distance") +
  theme(legend.position = "bottom",
        plot.title = element_text(size = 12, face = "bold"),
        plot.subtitle = element_text(size = 9))

# Combine in faceted style (1 column, 3 rows) using ggarrange
# This provides the same visual effect as facet_wrap but allows proper axis labels for each panel type
p6_margin <- ggarrange(p6_rank1_plot, p6_maha_scatter, p6_mad_scatter,
                       ncol = 1, nrow = 3, heights = c(1, 1, 1),
                       common.legend = FALSE, align = "v")

# --- Plot 7: Batch Separability AUC (Enhanced p6) ---------------------------
# Train logistic regression to predict batch from PCs, compute AUC
compute_batch_auc <- function(pc_scores, labels) {
  if (length(unique(labels)) < 2) return(NA_real_)
  # Use data.frame to avoid data.table [.data.table indexing issues
  pc_df <- as.data.frame(pc_scores[, seq_len(min(10, ncol(pc_scores))), drop = FALSE])
  pc_df$batch <- as.integer(factor(labels)) - 1L  # 0/1 encoding for binomial GLM
  complete <- complete.cases(pc_df)
  if (sum(complete) < 50) return(NA_real_)
  pc_df <- pc_df[complete, , drop = FALSE]

  # Simple logistic regression
  fit <- tryCatch(
    glm(batch ~ ., data = pc_df, family = binomial()),
    error = function(e) NULL
  )
  if (is.null(fit)) return(NA_real_)

  # Compute AUC
  pred_prob <- predict(fit, type = "response")
  auc <- tryCatch(
    pROC::auc(pc_df$batch, pred_prob, quiet = TRUE),
    error = function(e) NA_real_
  )
  as.numeric(auc)
}

batch_auc_pre  <- compute_batch_auc(cohort_pre$scores,  cohort_pre$labels)
batch_auc_post <- compute_batch_auc(cohort_post$scores, cohort_post$labels)

# Enhanced p6 with AUC
sep_dt_enhanced <- data.table(
  metric = rep(c("Silhouette", "kBET (k-NN Batch Effect Test) Acceptance",
                 "LISI (Local Inverse Simpson's Index) Score", "Batch Separability AUC"), each = 2),
  stage  = rep(c("Pre", "Post"), 4),
  value  = c(sil_pre, sil_post, kbet_pre, kbet_post, lisi_pre, lisi_post,
             ifelse(is.na(batch_auc_pre), 0, batch_auc_pre),
             ifelse(is.na(batch_auc_post), 0, batch_auc_post))
)
sep_dt_enhanced[, stage := factor(stage, levels = c("Pre", "Post"))]
sep_df_enhanced <- as.data.frame(sep_dt_enhanced)

p6_enhanced <- ggplot(sep_df_enhanced, aes(x = stage, y = value, fill = stage)) +
  geom_col(alpha = 0.8) +
  facet_wrap(~metric, scales = "free_y", ncol = 2) +
  scale_fill_manual(values = c(col_pre, col_post)) +
  labs(title = "Batch Separability Metrics (Enhanced)",
       subtitle = "Silhouette: lower = better | kBET: higher = better | LISI: 2.0 = perfect | AUC: 0.5 = no separation",
       y = "Value", x = NULL) +
  theme(legend.position = "none")

# --- Plot 8: PC-by-Batch Effect Size Reduction (R² per PC) ------------------
compute_pc_batch_effect <- function(pc_scores, labels) {
  n_pcs <- min(10, ncol(pc_scores))
  r2_per_pc <- numeric(n_pcs)
  for (j in seq_len(n_pcs)) {
    fit <- tryCatch(
      lm(pc_scores[, j] ~ factor(labels)),
      error = function(e) NULL
    )
    if (!is.null(fit)) {
      r2_per_pc[j] <- summary(fit)$r.squared
    } else {
      r2_per_pc[j] <- NA_real_
    }
  }
  r2_per_pc
}

pc_r2_pre  <- compute_pc_batch_effect(cohort_pre$scores,  cohort_pre$labels)
pc_r2_post <- compute_pc_batch_effect(cohort_post$scores, cohort_post$labels)

pc_effect_dt <- data.table(
  PC = rep(paste0("PC", seq_len(length(pc_r2_pre))), 2),
  stage = rep(c("Pre-harmonisation", "Post-harmonisation"), each = length(pc_r2_pre)),
  R2 = c(pc_r2_pre, pc_r2_post)
)
pc_effect_dt[, stage := factor(stage, levels = c("Pre-harmonisation", "Post-harmonisation"))]
pc_effect_df <- as.data.frame(pc_effect_dt)

p8_pc_effect <- ggplot(pc_effect_df, aes(x = PC, y = R2, colour = stage, group = stage)) +
  geom_line(linewidth = 1.2, alpha = 0.8) +
  geom_point(size = 3, alpha = 0.8) +
  scale_colour_manual(values = c(col_pre, col_post)) +
  labs(title = "PC-by-Batch Effect Size Reduction (R² per PC)",
       subtitle = "R² from PC ~ batch regression | Lower R² = less batch effect",
       x = "Principal Component", y = "R² (Batch Effect)") +
  theme(legend.position = "bottom", legend.title = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))

# --- Plot 9: Outlier Drill-Down Dashboard (For Flagged Bridges Only) ---------
# Create panels for flagged bridges
if (n_outliers > 0) {
  flagged_pairs <- as.data.frame(bridge_pair_dt[outlier_flag == TRUE])
  outlier_panels <- list()

  for (idx in seq_len(min(6, nrow(flagged_pairs)))) {  # Limit to 6 for space
    pair_idx <- which(bridge_pair_dt$finngenid == flagged_pairs$finngenid[idx])
    if (length(pair_idx) == 0 || length(pair_idx) > 1) {
      log_warn("  Skipping outlier pair {flagged_pairs$finngenid[idx]}: pair_idx length = {length(pair_idx)}")
      next
    }

    # Get PC coordinates for this pair (pair_idx should be length 1)
    sid_b1_pair <- sids_b1[pair_idx]
    sid_b2_pair <- sids_b2[pair_idx]

    # Verify we have valid sample IDs
    if (length(sid_b1_pair) != 1 || length(sid_b2_pair) != 1) {
      log_warn("  Skipping outlier pair {flagged_pairs$finngenid[idx]}: invalid sample IDs")
      next
    }

    if (sid_b1_pair %in% rownames(pc_pre$b1) && sid_b2_pair %in% rownames(pc_pre$b2) &&
        sid_b1_pair %in% rownames(pc_post$b1) && sid_b2_pair %in% rownames(pc_post$b2)) {
      # Extract PC coordinates (ensure single values, not vectors)
      pc1_pre_b1 <- as.numeric(pc_pre$b1[sid_b1_pair, 1])
      pc1_pre_b2 <- as.numeric(pc_pre$b2[sid_b2_pair, 1])
      pc1_post_b1 <- as.numeric(pc_post$b1[sid_b1_pair, 1])
      pc1_post_b2 <- as.numeric(pc_post$b2[sid_b2_pair, 1])
      pc2_pre_b1 <- as.numeric(pc_pre$b1[sid_b1_pair, 2])
      pc2_pre_b2 <- as.numeric(pc_pre$b2[sid_b2_pair, 2])
      pc2_post_b1 <- as.numeric(pc_post$b1[sid_b1_pair, 2])
      pc2_post_b2 <- as.numeric(pc_post$b2[sid_b2_pair, 2])

      # Ensure all are single values (not vectors) - take first element if vector
      if (length(pc1_pre_b1) > 1) pc1_pre_b1 <- pc1_pre_b1[1]
      if (length(pc1_pre_b2) > 1) pc1_pre_b2 <- pc1_pre_b2[1]
      if (length(pc1_post_b1) > 1) pc1_post_b1 <- pc1_post_b1[1]
      if (length(pc1_post_b2) > 1) pc1_post_b2 <- pc1_post_b2[1]
      if (length(pc2_pre_b1) > 1) pc2_pre_b1 <- pc2_pre_b1[1]
      if (length(pc2_pre_b2) > 1) pc2_pre_b2 <- pc2_pre_b2[1]
      if (length(pc2_post_b1) > 1) pc2_post_b1 <- pc2_post_b1[1]
      if (length(pc2_post_b2) > 1) pc2_post_b2 <- pc2_post_b2[1]

      pair_pc_df <- data.frame(
        stage = rep(c("Pre", "Post"), each = 2),
        batch = rep(c("B1", "B2"), 2),
        PC1 = c(pc1_pre_b1, pc1_pre_b2, pc1_post_b1, pc1_post_b2),
        PC2 = c(pc2_pre_b1, pc2_pre_b2, pc2_post_b1, pc2_post_b2)
      )

      # Connector line (use the extracted scalar values)
      connectors <- data.frame(
        x = c(pc1_pre_b1, pc1_post_b1),
        y = c(pc2_pre_b1, pc2_post_b1),
        xend = c(pc1_pre_b2, pc1_post_b2),
        yend = c(pc2_pre_b2, pc2_post_b2),
        stage = c("Pre", "Post")
      )

      # Safely extract distance metrics for subtitle (BEFORE ggplot chain)
      d_pre_val <- if (pair_idx >= 1 && pair_idx <= length(d_mad_pre)) d_mad_pre[pair_idx] else NA_real_
      d_post_val <- if (pair_idx >= 1 && pair_idx <= length(d_mad_post)) d_mad_post[pair_idx] else NA_real_
      l_ratio_val <- if (pair_idx >= 1 && pair_idx <= length(log_ratio_mad)) log_ratio_mad[pair_idx] else NA_real_
      margin_val <- if (pair_idx >= 1 && pair_idx <= length(pairing_margins)) pairing_margins[pair_idx] else NA_real_

      p_outlier <- ggplot(pair_pc_df, aes(x = PC1, y = PC2, colour = batch, shape = batch)) +
        geom_segment(data = connectors, aes(x = x, y = y, xend = xend, yend = yend),
                     inherit.aes = FALSE,
                     colour = "grey50", alpha = 0.6, linewidth = 0.8) +
        geom_point(size = 4, alpha = 0.9) +
        facet_wrap(~stage) +
        scale_colour_manual(values = c("B1" = "#2c7bb6", "B2" = "#d7191c")) +
        scale_shape_manual(values = c("B1" = 16, "B2" = 17)) +
        labs(title = sprintf("Outlier: %s", substr(flagged_pairs$finngenid[idx], 1, 10)),
             subtitle = sprintf("D_pre=%.2f, D_post=%.2f, L=%.3f, Margin=%.2f",
                               d_pre_val, d_post_val, l_ratio_val, margin_val),
             x = "PC1", y = "PC2") +
        theme(legend.position = "none", plot.title = element_text(size = 9),
              plot.subtitle = element_text(size = 7))

      outlier_panels[[length(outlier_panels) + 1]] <- p_outlier
    }
  }

  if (length(outlier_panels) > 0) {
    p9_outlier_dashboard <- ggarrange(plotlist = outlier_panels, ncol = 3, nrow = 2)
  } else {
    p9_outlier_dashboard <- NULL
  }
} else {
  p9_outlier_dashboard <- NULL
}

# --- Panel 8: Outlier flag barplot (keep existing) ---------------------------
# Create outlier plot data as data.frame to avoid data.table indexing issues
outlier_plot_df <- data.frame(
  finngenid = bridge_pair_dt$finngenid,
  eucl_post = bridge_pair_dt$eucl_post,
  maha_post = bridge_pair_dt$maha_post,
  noise_floor_ratio = bridge_pair_dt$noise_floor_ratio,
  outlier_flag = bridge_pair_dt$outlier_flag,
  finngenid_short = substr(bridge_pair_dt$finngenid, 1, 10)
)

p8 <- ggplot(outlier_plot_df, aes(x = reorder(finngenid_short, -eucl_post), y = eucl_post)) +
  geom_col(aes(fill = outlier_flag), alpha = 0.8) +
  scale_fill_manual(values = c("FALSE" = col_post, "TRUE" = col_accent),
                    labels = c("Normal", "Flagged")) +
  labs(title = "Bridge Pair Post-Harmonisation Distance",
       subtitle = sprintf("Flagged outliers: %d/%d (threshold: Mahal > Q3+1.5*IQR or noise ratio > 3)", n_outliers, n_bridge),
       x = "Bridge Pair (FINNGENID)", y = "Euclidean Distance") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
        legend.position = "bottom", legend.title = element_blank())

# --- Assemble multi-panel PDF ------------------------------------------------
# Save as multi-page PDF
# Wrap PDF generation in tryCatch to catch and log exact error location
tryCatch({
  pdf(dashboard_path, width = 16, height = 20)

  # Page 1: Top 3 Quick-Inspection Plots
  # Dumbbell + PCA scatter (pheatmap rendered on separate page)
  page1_top <- ggarrange(p1_dumbbell, p5, ncol = 1, nrow = 2,
                         labels = c("A", "B"), heights = c(1, 1.2))
  print(annotate_figure(page1_top,
    top = text_grob("Cross-Batch Harmonisation KPI Dashboard — Top Quick-Inspection Plots",
                    face = "bold", size = 14)))

  # Page 2: Distance Metrics (4 panels: Boxplot-Wilcoxon, ECDF, Spaghetti, Log-Ratio)
  page2_dist <- ggarrange(p0_boxplot, p2_dist, p2a, p2b,
                           ncol = 2, nrow = 2, labels = c("A", "B", "C", "D"))
  print(annotate_figure(page2_dist,
    top = text_grob("Cross-Batch Harmonisation KPI Dashboard — Distance Metrics",
                    face = "bold", size = 14)))

  # Page 3: NN Pairing Heatmap (ggplot, full-page)
  print(annotate_figure(p5_heatmap,
    top = text_grob("Cross-Batch Harmonisation KPI Dashboard — Nearest-Neighbor Pairing Heatmap",
                    face = "bold", size = 14)))

  # Page 4: Pairing & Identity Tests (Faceted: Margin + Rank-1)
  print(annotate_figure(p6_margin,
    top = text_grob("Cross-Batch Harmonisation KPI Dashboard — Pairing & Identity Tests",
                    face = "bold", size = 14)))

  # Page 5: PCA with Connectors & Batch Separability
  page5_sep <- ggarrange(p5, p6_enhanced, ncol = 1, nrow = 2,
                          labels = c("A", "B"), heights = c(1.2, 1))
  print(annotate_figure(page5_sep,
    top = text_grob("Cross-Batch Harmonisation KPI Dashboard — PCA & Batch Separability",
                    face = "bold", size = 14)))

  # Page 6: PC Effects & Additional Metrics
  page6_list <- list(p8_pc_effect, p3, p8)
  if (!is.null(p7)) {
    page6_list <- c(page6_list, list(p7))
    page6 <- ggarrange(plotlist = page6_list, ncol = 2, nrow = 2,
                        labels = c("A", "B", "C", "D"))
  } else {
    page6 <- ggarrange(plotlist = page6_list, ncol = 2, nrow = 2,
                        labels = c("A", "B", "C"))
  }
  print(annotate_figure(page6,
    top = text_grob("Cross-Batch Harmonisation KPI Dashboard — PC Effects & Additional Metrics",
                    face = "bold", size = 14)))

  # Page 7: Protein Correlation Heatmap (Pre-harmonisation)
  print(annotate_figure(p5b_corr_pre,
    top = text_grob("Cross-Batch Harmonisation KPI Dashboard — Bridge Sample Protein Correlation (Pre-Harmonisation)",
                    face = "bold", size = 14)))

  # Page 8: Protein Correlation Heatmap (Post-harmonisation)
  print(annotate_figure(p5b_corr_post,
    top = text_grob("Cross-Batch Harmonisation KPI Dashboard — Bridge Sample Protein Correlation (Post-Harmonisation)",
                    face = "bold", size = 14)))

  # Page 9: Outlier Drill-Down Dashboard (if any outliers)
  if (!is.null(p9_outlier_dashboard)) {
    tryCatch({
      print(annotate_figure(p9_outlier_dashboard,
        top = text_grob("Cross-Batch Harmonisation KPI Dashboard — Outlier Drill-Down",
                        face = "bold", size = 14)))
    }, error = function(e) {
      log_warn("Failed to generate outlier dashboard page: {e$message}")
    })
  }

  dev.off()
  log_info("  Saved dashboard PDF: {dashboard_path}")
}, error = function(e) {
  # Try to close PDF device if it's open
  tryCatch(dev.off(), error = function(e2) NULL)
  log_error("Failed to generate PDF dashboard: {e$message}")
  log_error("Error occurred at: {deparse(e$call)}")
  # Re-throw to stop execution
  stop(e)
})

# ==============================================================================
# FINAL SUMMARY
# ==============================================================================
log_info("=" |> rep(70) |> paste(collapse = ""))
log_info("HARMONISATION KPI SUMMARY")
log_info("=" |> rep(70) |> paste(collapse = ""))
log_info("  Bridge pairs analysed: {n_bridge}")
log_info("  MAD-whitened distance reduction: {round(median(d_mad_pre),2)} -> {round(median(d_mad_post),2)}")
log_info("  Paired log-ratio (MAD-whitened): {round(median_log_ratio_mad, 4)} [{round(boot_ci_mad[1],4)}, {round(boot_ci_mad[2],4)}]")
log_info("  Noise floor ratio: {round(noise_floor_ratio, 2)}")
log_info("")
log_info("  NOISE FLOOR CALCULATION DETAILS:")
log_info("    The noise floor represents the measurement noise inherent in the assay system.")
log_info("    It is computed as the median MAD-whitened Euclidean distance between within-batch")
log_info("    replicate pairs (same biological sample measured multiple times in the same batch).")
log_info("    - Within-batch replicate pairs: B1={nrow(wb_rep_b1)}, B2={nrow(wb_rep_b2)}")
log_info("    - Noise floor (MAD-whitened): {round(noise_floor_mad, 3)}")
log_info("    - Bridge pair post-harmonisation distance (median): {round(median(d_mad_post), 3)}")
log_info("    - Ratio (bridge_post / noise): {round(noise_floor_ratio, 3)}")
log_info("    Interpretation: Ratio < 2.0 indicates bridge pairs are approaching measurement noise floor.")
if (!is.null(var_decomp_pre) && !is.null(var_decomp_post)) {
  high_var_proteins <- unique(c(
    var_decomp_pre$protein[var_decomp_pre$batch_var_frac > 0.3 | var_decomp_pre$sex_var_frac > 0.3],
    var_decomp_post$protein[var_decomp_post$batch_var_frac > 0.3 | var_decomp_post$sex_var_frac > 0.3]
  ))
  if (length(high_var_proteins) > 0) {
    log_info("")
    log_info("  HIGH-VARIANCE PROTEINS (>30% variance explained by Batch or Sex):")
    for (prot in head(high_var_proteins, 20)) {
      pre_batch <- var_decomp_pre$batch_var_frac[var_decomp_pre$protein == prot]
      pre_sex <- var_decomp_pre$sex_var_frac[var_decomp_pre$protein == prot]
      post_batch <- var_decomp_post$batch_var_frac[var_decomp_post$protein == prot]
      post_sex <- var_decomp_post$sex_var_frac[var_decomp_post$protein == prot]
      if (length(pre_batch) > 0 && length(post_batch) > 0) {
        log_info("    {prot}: Pre (Batch={round(pre_batch[1]*100,1)}%, Sex={round(pre_sex[1]*100,1)}%) | Post (Batch={round(post_batch[1]*100,1)}%, Sex={round(post_sex[1]*100,1)}%)")
      }
    }
    if (length(high_var_proteins) > 20) {
      log_info("    ... and {length(high_var_proteins) - 20} more proteins")
    }
  }
}
log_info("  Rank-1 ID rate: {id_pre$rank1_hits}/{n_bridge} -> {id_post$rank1_hits}/{n_bridge}")
log_info("  Silhouette: {round(sil_pre,4)} -> {round(sil_post,4)}")
log_info("  kBET (k-Nearest Neighbour Batch Effect Test) acceptance: {round(kbet_pre,4)} -> {round(kbet_post,4)}")
log_info("  LISI (Local Inverse Simpson's Index) score: {round(lisi_pre,4)} -> {round(lisi_post,4)}")
log_info("  Bridge outliers: {n_outliers}")
log_info("  Poorly harmonised proteins: {n_poor}")
log_info("  Pass/Fail summary:")
for (i in seq_len(nrow(kpi_dt))) {
  if (!is.na(kpi_dt$pass[i])) {
    log_info("    [{kpi_dt$pass[i]}] {kpi_dt$metric[i]}")
  }
}
log_info("=" |> rep(70) |> paste(collapse = ""))
log_info("Step 07b complete")

cat("\n=== CROSS-BATCH HARMONISATION KPI EVALUATION COMPLETE ===\n")
cat(sprintf("  KPI summary: %s\n", kpi_path))
cat(sprintf("  Dashboard: %s\n", dashboard_path))
cat(sprintf("  Bridge pairs: %d | Outliers: %d | Poor proteins: %d\n", n_bridge, n_outliers, n_poor))
cat(sprintf("  Log-ratio: %.3f [%.3f, %.3f]\n", median_log_ratio, boot_ci[1], boot_ci[2]))
