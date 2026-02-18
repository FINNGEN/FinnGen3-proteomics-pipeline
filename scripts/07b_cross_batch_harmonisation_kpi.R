#!/usr/bin/env Rscript
# ==============================================================================
# 07b_cross_batch_harmonisation_kpi.R
# Post-hoc evaluation of cross-batch bridge normalisation efficacy
# ==============================================================================
#
# Purpose:
#   Runs immediately after 07_bridge_normalization.R. Computes KPIs that quantify
#   whether the two batches are properly aligned after bridge normalisation,
#   while checking that biological signal has not been collapsed by over-correction.
#
# Metrics computed (6 evaluation steps):
#   Step 0 — Fixed PCA basis (avoids basis leakage; trained on combined pre-
#            harmonisation data so post cannot influence basis)
#   Step 1 — Bridge pair collapse (MAD-whitened Euclidean [PRIMARY], Euclidean,
#            Mahalanobis, ICC, CCC, paired log-ratio with bootstrap CI,
#            noise-floor ratio)
#   Step 2 — Pair identity test (rank-1 identification rate)
#   Step 3 — Batch separability (Silhouette, kBET, LISI, Batch AUC)
#   Step 4 — Biology preservation (sex prediction accuracy, variance decomposition)
#   Step 5 — Outlier detection (bridge pair flags, per-protein concordance)
#   Step 6 — KPI summary table & multi-panel PDF dashboard
#
# Design decisions:
#   1) Pre-harmonisation matrices MUST come from QC-passed outputs (05d). A
#      fallback is allowed only if 05d is missing; this is logged loudly.
#   2) Primary bridge alignment distance is MAD-whitened Euclidean distance in
#      fixed PC space (robust to heavy tails; scale-invariant per PC).
#   3) Noise floor is computed from within-batch technical replicates (same
#      FINNGENID within batch) using the same MAD-whitening. Ratio ~1 ideal.
#   4) Batch separability uses multiple views: silhouette / kBET / LISI / AUC.
#      Some separation post-harmonisation can be "true biology".
#
# Outputs:
#   - 07b/kpi_summary.tsv
#   - 07b/bridge_pair_distances.tsv
#   - 07b/noise_floor_analysis.tsv
#   - 07b/protein_concordance.tsv
#   - 07b/poorly_harmonised_proteins.tsv (if any)
#   - 07b/outlier_flags.tsv
#   - 07b/variance_decomposition.tsv (if available)
#   - 07b/kpi_dashboard.pdf
#
# Author: Reza Jabal, PhD (rjabal@broadinstitute.org)
# Date: February 2026
# ==============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggpubr)
  library(cluster)
  library(MASS)
  library(FNN)
  library(irr)
  library(DescTools)
  library(kBET)
  library(lisi)
  library(pROC)
  library(yaml)
  library(logger)
})

# ------------------------------------------------------------------------------
# Boilerplate: script directory, config, batch context, logging
# ------------------------------------------------------------------------------
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

log_info(strrep("=", 80))
log_info("Step {step_num}: Cross-Batch Harmonisation KPI Evaluation")
log_info("Batch context: {batch_id}")
log_info(strrep("=", 80))

# ------------------------------------------------------------------------------
# Aesthetics (clean + publication-ish)
# ------------------------------------------------------------------------------
theme_set(theme_bw(base_size = 12))
theme_update(
  plot.title = element_text(face = "bold", size = 13),
  plot.subtitle = element_text(size = 10, color = "grey30", lineheight = 1.3),
  axis.title = element_text(face = "bold", size = 11),
  legend.title = element_text(face = "bold"),
  panel.grid.minor = element_blank()
)

# Palette (consistent across dashboard)
PAL <- list(
  pre = "#E76F51",      # warm (pre)
  post = "#2A9D8F",     # teal (post)
  accent = "#8B0000",   # dark red (flags)
  neutral = "#264653",  # dark slate
  grey = "grey70"
)

# ------------------------------------------------------------------------------
# Guard: skip gracefully if disabled or single-batch
# ------------------------------------------------------------------------------
run_kpi <- tryCatch(isTRUE(config$parameters$bridge_normalization$run_kpi_evaluation),
                    error = function(e) FALSE)
multi_batch_mode <- tryCatch(isTRUE(config$parameters$normalization$multi_batch_mode),
                             error = function(e) FALSE)

if (!multi_batch_mode) {
  log_info("Single-batch mode — skipping cross-batch KPI evaluation")
  Sys.setenv(PIPELINE_STEP_SKIPPED = "TRUE")
  cat("\n=== CROSS-BATCH KPI EVALUATION SKIPPED (single-batch mode) ===\n")
  stop("STEP_SKIPPED")
}
if (!run_kpi) {
  log_info("KPI evaluation disabled in config (bridge_normalization.run_kpi_evaluation: false)")
  Sys.setenv(PIPELINE_STEP_SKIPPED = "TRUE")
  cat("\n=== CROSS-BATCH KPI EVALUATION SKIPPED (disabled in config) ===\n")
  stop("STEP_SKIPPED")
}

# ==============================================================================
# Helpers: assertions, safe IO, alignment, robust stats
# ==============================================================================
stopifnot_msg <- function(cond, msg) {
  if (!isTRUE(cond)) stop(msg, call. = FALSE)
}

`%||%` <- function(a, b) if (!is.null(a)) a else b

safe_median <- function(x) median(x, na.rm = TRUE)

matrix_missingness <- function(M) {
  if (length(M) == 0) return(list(na = 0L, total = 0L, frac = 0))
  na_count <- sum(is.na(M))
  total <- length(M)
  list(na = na_count, total = total, frac = ifelse(total > 0, na_count / total, 0))
}

ensure_named_rownames <- function(mat, label) {
  stopifnot_msg(is.matrix(mat) || is.data.frame(mat), paste0(label, " must be matrix/data.frame"))
  mat <- as.matrix(mat)
  stopifnot_msg(!is.null(rownames(mat)) && all(rownames(mat) != ""), paste0(label, " missing rownames"))
  stopifnot_msg(!is.null(colnames(mat)) && all(colnames(mat) != ""), paste0(label, " missing colnames"))
  mat
}

read_rds_or_stop <- function(path, what) {
  if (!file.exists(path)) stop(paste0(what, " not found: ", path), call. = FALSE)
  readRDS(path)
}

# Strict: prefer 05d QC-passed; allow fallback only with loud logging.
load_pre_matrix_qc_passed <- function(bid, config) {
  p_qc <- get_output_path("05d", "05d_npx_matrix_all_qc_passed", bid, "phenotypes", "rds", config = config)
  if (file.exists(p_qc)) {
    log_info("Pre matrix for {bid}: using QC-passed (05d) at {p_qc}")
    return(ensure_named_rownames(readRDS(p_qc), paste0("pre_", bid)))
  }
  p_fallback <- get_output_path("00", "npx_matrix_analysis_ready", bid, "qc", config = config)
  if (file.exists(p_fallback)) {
    log_warn("Pre matrix for {bid}: 05d QC-passed missing; FALLBACK to 00 at {p_fallback}")
    log_warn("This fallback should be rare; verify upstream pipeline outputs.")
    return(ensure_named_rownames(readRDS(p_fallback), paste0("pre_", bid)))
  }
  stop(paste0("Pre-harmonisation matrix not found for ", bid,
              ". Tried: ", p_qc, " and ", p_fallback), call. = FALSE)
}

load_post_matrix <- function(bid, config) {
  p <- get_output_path("07", "npx_matrix_cross_batch_bridge", bid, "normalized", config = config)
  if (!file.exists(p)) stop("Post-harmonisation matrix not found for ", bid, ": ", p, call. = FALSE)
  log_info("Post matrix for {bid}: {p}")
  ensure_named_rownames(readRDS(p), paste0("post_", bid))
}

# Deterministic “first sample per FINNGENID” (stable across runs)
first_sample_per_finngenid <- function(map_dt, sids_available) {
  map_dt <- as.data.table(map_dt)
  stopifnot_msg(all(c("FINNGENID", "SampleID") %in% names(map_dt)), "sample_mapping must include FINNGENID and SampleID")
  map_dt <- map_dt[!is.na(FINNGENID) & !is.na(SampleID)]
  map_dt <- map_dt[SampleID %in% sids_available]
  setorder(map_dt, FINNGENID, SampleID)
  map_dt[!duplicated(FINNGENID)]
}

# Find within-batch replicate pairs (noise floor)
find_within_batch_replicates <- function(mapping, npx_matrix) {
  mapping <- as.data.table(mapping)
  mapping <- mapping[!is.na(FINNGENID) & !is.na(SampleID)]
  dup_fgs <- mapping[, .N, by = FINNGENID][N >= 2]$FINNGENID
  pairs <- vector("list", length(dup_fgs))
  k <- 0L
  for (fg in dup_fgs) {
    sids <- mapping[FINNGENID == fg, SampleID]
    sids <- intersect(sids, rownames(npx_matrix))
    if (length(sids) >= 2) {
      k <- k + 1L
      pairs[[k]] <- data.table(finngenid = fg, sid1 = sids[1], sid2 = sids[2])
    }
  }
  if (k == 0) return(data.table())
  rbindlist(pairs[seq_len(k)])
}

# Robust covariance inverse (stable; avoids singularity)
robust_cov_inv <- function(X, ridge = 1e-6) {
  X <- as.matrix(X)
  cov_mat <- tryCatch({
    # MASS::cov.rob can fail if n too small or singular; fallback below
    MASS::cov.rob(X, method = "mcd")$cov
  }, error = function(e) {
    cov(X)
  })
  cov_mat <- cov_mat + diag(ridge, nrow(cov_mat))
  cov_inv <- tryCatch(solve(cov_mat), error = function(e) MASS::ginv(cov_mat))
  list(cov = cov_mat, inv = cov_inv)
}

# Impute NAs column-wise using precomputed reference medians (named vector).
# Uses named matching for safety — column names must match ref_medians names.
impute_with_reference_medians <- function(M, ref_medians) {
  M <- as.matrix(M)
  stopifnot_msg(!is.null(colnames(M)), "Matrix has no colnames; cannot impute safely.")
  stopifnot_msg(!is.null(names(ref_medians)), "ref_medians must be a named vector.")
  idx <- match(colnames(M), names(ref_medians))
  stopifnot_msg(!anyNA(idx), "Some matrix columns missing from ref_medians names.")
  meds <- ref_medians[idx]
  for (j in seq_len(ncol(M))) {
    nas <- is.na(M[, j])
    if (any(nas)) M[nas, j] <- meds[j]
  }
  M
}

# ------------------------------------------------------------------------------
# KPI NA-policy and PCA parameters (configurable with sensible defaults)
# ------------------------------------------------------------------------------
pca_impute_max_frac <- tryCatch(config$parameters$kpi$na_impute_overall_max_frac, error = function(e) NULL) %||% 0.001
pca_min_complete_rows <- tryCatch(config$parameters$kpi$complete_case_pca_min_rows, error = function(e) NULL) %||% 200
kpi_max_pcs <- tryCatch(config$parameters$kpi$pca_max_pcs, error = function(e) NULL) %||% 50
kpi_var_threshold <- tryCatch(config$parameters$kpi$pca_var_threshold, error = function(e) NULL) %||% 0.95

log_info("PCA NA policy: complete-case preferred; impute only if missingness <= {pca_impute_max_frac*100}%")
log_info("PCA settings: max_pcs={kpi_max_pcs}, var_threshold={kpi_var_threshold}, min_complete_rows={pca_min_complete_rows}")

# ==============================================================================
# 0) DATA LOADING
# ==============================================================================
log_info("Phase 0: Loading data")
all_batch_ids <- names(config$batches)
stopifnot_msg(length(all_batch_ids) >= 2, paste0("Cross-batch KPI requires >=2 batches; found ", length(all_batch_ids)))

batch1_id <- all_batch_ids[1]
batch2_id <- all_batch_ids[2]
log_info("Batch1={batch1_id} | Batch2={batch2_id}")

# Pre (QC-passed)
pre_b1 <- load_pre_matrix_qc_passed(batch1_id, config)
pre_b2 <- load_pre_matrix_qc_passed(batch2_id, config)

# Post (bridge normalised)
post_b1 <- load_post_matrix(batch1_id, config)
post_b2 <- load_post_matrix(batch2_id, config)

log_info("Pre B1: {nrow(pre_b1)} x {ncol(pre_b1)} | Pre B2: {nrow(pre_b2)} x {ncol(pre_b2)}")
log_info("Post B1: {nrow(post_b1)} x {ncol(post_b1)} | Post B2: {nrow(post_b2)} x {ncol(post_b2)}")

# Sample mappings
map_b1 <- read_rds_or_stop(get_output_path("00", "sample_mapping", batch1_id, "qc", config = config), "sample_mapping B1")
map_b2 <- read_rds_or_stop(get_output_path("00", "sample_mapping", batch2_id, "qc", config = config), "sample_mapping B2")
map_b1 <- as.data.table(map_b1)
map_b2 <- as.data.table(map_b2)

# Bridge metadata
bridging_file <- tryCatch(get_batch_input_path("bridging_samples_file", batch1_id, config), error = function(e) NULL)
stopifnot_msg(!is.null(bridging_file) && file.exists(bridging_file), "Bridging metadata file not found.")
bridging_meta <- fread(bridging_file, sep = "\t")

# Identify direct B1<->B2 (excluding Batch_00 / EA5)
b01_flag <- bridging_meta$in_Batch_01 %in% c(TRUE, "TRUE", "T")
b02_flag <- bridging_meta$in_Batch_02 %in% c(TRUE, "TRUE", "T")
b00_flag <- bridging_meta$in_Batch_00 %in% c(FALSE, "FALSE", "F")
bridge_finngenids <- bridging_meta[b01_flag & b02_flag & b00_flag]$FINNGENID
bridge_finngenids <- unique(na.omit(bridge_finngenids))
log_info("Direct B1-B2 bridge FINNGENIDs: {length(bridge_finngenids)}")

# Bridge SampleIDs per batch (deterministic)
bridge_sid_b1 <- first_sample_per_finngenid(map_b1[FINNGENID %in% bridge_finngenids], rownames(pre_b1))
bridge_sid_b2 <- first_sample_per_finngenid(map_b2[FINNGENID %in% bridge_finngenids], rownames(pre_b2))

common_bridge_fg <- intersect(bridge_sid_b1$FINNGENID, bridge_sid_b2$FINNGENID)
stopifnot_msg(length(common_bridge_fg) >= 5, "Too few aligned bridge pairs; check bridging metadata and mappings.")

bridge_sid_b1 <- bridge_sid_b1[match(common_bridge_fg, FINNGENID)]
bridge_sid_b2 <- bridge_sid_b2[match(common_bridge_fg, FINNGENID)]
n_bridge <- length(common_bridge_fg)

sids_b1 <- bridge_sid_b1$SampleID
sids_b2 <- bridge_sid_b2$SampleID

log_info("Aligned bridge pairs: {n_bridge}")

# --- Configurable bridge sample exclusion ---
# The user may supply either FINNGENIDs or SampleIDs in the config list.
# Match against common_bridge_fg (FINNGENIDs) AND sids_b1/sids_b2 (SampleIDs)
# so that exclusion works regardless of identifier type.
exclude_ids <- tryCatch(
  config$parameters$bridge_normalization$exclude_bridge_sample_ids,
  error = function(e) NULL) %||% character(0)
if (length(exclude_ids) > 0) {
  log_info("Excluding {length(exclude_ids)} bridge IDs from config: {paste(exclude_ids, collapse=', ')}")
  keep <- !(common_bridge_fg %in% exclude_ids) &
          !(sids_b1 %in% exclude_ids) &
          !(sids_b2 %in% exclude_ids)
  n_before <- n_bridge
  common_bridge_fg <- common_bridge_fg[keep]
  sids_b1 <- sids_b1[keep]
  sids_b2 <- sids_b2[keep]
  n_bridge <- length(common_bridge_fg)
  log_info("Bridge pairs after exclusion: {n_bridge} (removed {n_before - n_bridge})")
  stopifnot_msg(n_bridge >= 3, "Too few bridge pairs remaining after exclusion.")
}

# Common proteins across all matrices
common_proteins <- Reduce(intersect, list(colnames(pre_b1), colnames(pre_b2), colnames(post_b1), colnames(post_b2)))
stopifnot_msg(length(common_proteins) >= 200, paste0("Too few common proteins (", length(common_proteins), ")."))
log_info("Common proteins across all matrices: {length(common_proteins)}")

# Within-batch replicates
wb_rep_b1 <- find_within_batch_replicates(map_b1, pre_b1)
wb_rep_b2 <- find_within_batch_replicates(map_b2, pre_b2)
log_info("Within-batch replicate pairs: B1={nrow(wb_rep_b1)}, B2={nrow(wb_rep_b2)}")

# ==============================================================================
# STEP 0: FIXED PCA BASIS (avoid basis leakage)
#   NA policy: complete-case first; impute only if missingness <= threshold.
# ==============================================================================
log_info("Step 0: Fixed PCA basis on combined PRE data (no leakage from POST)")

build_fixed_pca <- function(mat_b1, mat_b2, proteins,
                            max_pcs = 50, var_threshold = 0.95,
                            impute_max_frac = 0.001, min_complete_rows = 200) {
  combined <- rbind(mat_b1[, proteins, drop = FALSE],
                    mat_b2[, proteins, drop = FALSE])

  # Reference medians (always computed, for allowed tiny-imputation or downstream use)
  ref_medians <- apply(combined, 2, safe_median)
  names(ref_medians) <- colnames(combined)

  miss <- matrix_missingness(combined)
  log_info("PCA training missingness: {miss$na}/{miss$total} ({round(miss$frac*100,4)}%)")

  method <- NA_character_
  train <- combined

  if (miss$frac <= impute_max_frac) {
    # Allowed tiny-imputation (<=threshold)
    method <- paste0("median_impute (<=", impute_max_frac * 100, "%)")
    if (miss$na > 0) {
      log_info("Missingness <= threshold; imputing NAs with reference medians for PCA training")
      train <- impute_with_reference_medians(train, ref_medians)
    } else {
      log_info("No missingness; PCA training without imputation")
    }
  } else {
    # COMPLETE-CASE FIRST: drop rows if enough remain; else drop NA columns
    cc_rows <- complete.cases(train)
    n_cc <- sum(cc_rows)
    if (n_cc >= min_complete_rows) {
      method <- "complete_case_drop_rows"
      log_info("Using complete-case PCA by dropping rows: {nrow(train)} -> {n_cc}")
      train <- train[cc_rows, , drop = FALSE]
    } else {
      # Fallback: drop NA columns (avoid imputation entirely)
      na_cols <- which(colSums(is.na(train)) > 0)
      keep_cols <- setdiff(seq_len(ncol(train)), na_cols)
      stopifnot_msg(length(keep_cols) >= 50,
                    paste0("Too many proteins have NAs; dropping NA columns leaves only ",
                           length(keep_cols), " proteins."))
      method <- "complete_case_drop_cols"
      log_warn("Too few complete rows ({n_cc}); switching to complete-case PCA by dropping NA columns: {length(na_cols)} dropped")
      train <- train[, keep_cols, drop = FALSE]
      proteins <- colnames(train)
    }
  }

  # PCA
  pca_fit <- prcomp(train, center = TRUE, scale. = TRUE)
  cum_var <- cumsum(pca_fit$sdev^2 / sum(pca_fit$sdev^2))
  idx <- which(cum_var >= var_threshold)[1]
  n_pcs <- if (is.na(idx)) min(max_pcs, ncol(train)) else idx
  n_pcs <- min(n_pcs, max_pcs, ncol(pca_fit$rotation))

  log_info("PCA basis built using method: {method}")
  log_info("Retaining {n_pcs} PCs (variance explained={round(cum_var[n_pcs]*100,1)}%)")

  list(
    method = method,
    proteins = proteins,
    n_pcs = n_pcs,
    ref_medians = ref_medians[proteins],
    center = pca_fit$center,
    scale = pca_fit$scale,
    rotation = pca_fit$rotation[, seq_len(n_pcs), drop = FALSE],
    var_explained = cum_var[n_pcs],
    training_missing_frac = miss$frac
  )
}

project_onto_basis <- function(mat, basis, impute_max_frac = 0.001) {
  M <- mat[, basis$proteins, drop = FALSE]
  miss <- matrix_missingness(M)

  if (miss$frac == 0) {
    M2 <- M
  } else if (miss$frac <= impute_max_frac) {
    log_info("  Projection: missingness {round(miss$frac*100,4)}% <= threshold; imputing with reference medians")
    M2 <- impute_with_reference_medians(M, basis$ref_medians)
  } else {
    # Complete-case projection: drop rows with any NA
    cc <- complete.cases(M)
    dropped <- sum(!cc)
    log_warn("  Projection: missingness {round(miss$frac*100,4)}% > threshold; dropping {dropped} rows for complete-case projection")
    M2 <- M[cc, , drop = FALSE]
  }

  M_scaled <- scale(M2, center = basis$center[basis$proteins], scale = basis$scale[basis$proteins])
  scores <- M_scaled %*% basis$rotation
  rownames(scores) <- rownames(M2)

  proj_method <- if (miss$frac == 0) "no_missing" else if (miss$frac <= impute_max_frac) "median_impute" else "complete_case_drop_rows"
  list(scores = scores, method = proj_method, missing_frac = miss$frac, n = nrow(scores))
}

pca_basis <- build_fixed_pca(pre_b1, pre_b2, common_proteins,
                             max_pcs = kpi_max_pcs, var_threshold = kpi_var_threshold,
                             impute_max_frac = pca_impute_max_frac,
                             min_complete_rows = pca_min_complete_rows)

proj_pre_b1  <- project_onto_basis(pre_b1,  pca_basis, impute_max_frac = pca_impute_max_frac)
proj_pre_b2  <- project_onto_basis(pre_b2,  pca_basis, impute_max_frac = pca_impute_max_frac)
proj_post_b1 <- project_onto_basis(post_b1, pca_basis, impute_max_frac = pca_impute_max_frac)
proj_post_b2 <- project_onto_basis(post_b2, pca_basis, impute_max_frac = pca_impute_max_frac)

log_info("PC scores pre: B1={proj_pre_b1$n} ({proj_pre_b1$method}), B2={proj_pre_b2$n} ({proj_pre_b2$method})")
log_info("PC scores post: B1={proj_post_b1$n} ({proj_post_b1$method}), B2={proj_post_b2$n} ({proj_post_b2$method})")

# Extract score matrices (preserving variable names used throughout downstream code)
pc_pre_b1  <- proj_pre_b1$scores
pc_pre_b2  <- proj_pre_b2$scores
pc_post_b1 <- proj_post_b1$scores
pc_post_b2 <- proj_post_b2$scores

# Reconcile bridge pairs: complete-case projection may have dropped some samples
keep_pre  <- (sids_b1 %in% rownames(pc_pre_b1))  & (sids_b2 %in% rownames(pc_pre_b2))
keep_post <- (sids_b1 %in% rownames(pc_post_b1)) & (sids_b2 %in% rownames(pc_post_b2))
keep_pair <- keep_pre & keep_post

dropped_pairs <- sum(!keep_pair)
if (dropped_pairs > 0) {
  log_warn("Dropping {dropped_pairs}/{n_bridge} bridge pairs due to incomplete PCA projection (NA policy).")
  common_bridge_fg <- common_bridge_fg[keep_pair]
  sids_b1 <- sids_b1[keep_pair]
  sids_b2 <- sids_b2[keep_pair]
  n_bridge <- length(common_bridge_fg)
}
stopifnot_msg(n_bridge >= 5, paste0("Too few bridge pairs after PCA NA-policy filtering: ", n_bridge))
log_info("Aligned bridge pairs (after PCA availability filter): {n_bridge}")

# ==============================================================================
# STEP 1: BRIDGE PAIR COLLAPSE (Primary KPI)
# ==============================================================================
log_info("Step 1: Bridge pair collapse metrics")

paired_euclidean <- function(pc_b1, pc_b2, sids_b1, sids_b2) {
  A <- pc_b1[sids_b1, , drop = FALSE]
  B <- pc_b2[sids_b2, , drop = FALSE]
  sqrt(rowSums((A - B)^2))
}

paired_mahalanobis <- function(pc_b1, pc_b2, sids_b1, sids_b2) {
  combined <- rbind(pc_b1, pc_b2)
  rc <- robust_cov_inv(combined)
  diffs <- pc_b1[sids_b1, , drop = FALSE] - pc_b2[sids_b2, , drop = FALSE]
  sqrt(rowSums((diffs %*% rc$inv) * diffs))
}

mad_whitened_distance <- function(pc_b1, pc_b2, sids_b1, sids_b2, mad_ref = NULL) {
  if (is.null(mad_ref)) {
    combined <- rbind(pc_b1, pc_b2)
    mad_ref <- apply(combined, 2, function(x) mad(x, na.rm = TRUE))
  }
  mad_ref[mad_ref == 0 | is.na(mad_ref)] <- 1

  A <- pc_b1[sids_b1, , drop = FALSE]
  B <- pc_b2[sids_b2, , drop = FALSE]
  diffs <- sweep(A - B, 2, mad_ref, "/")
  sqrt(rowSums(diffs^2, na.rm = TRUE))
}

# Distances (PC-space)
d_eucl_pre  <- paired_euclidean(pc_pre_b1,  pc_pre_b2,  sids_b1, sids_b2)
d_eucl_post <- paired_euclidean(pc_post_b1, pc_post_b2, sids_b1, sids_b2)
d_maha_pre  <- paired_mahalanobis(pc_pre_b1,  pc_pre_b2,  sids_b1, sids_b2)
d_maha_post <- paired_mahalanobis(pc_post_b1, pc_post_b2, sids_b1, sids_b2)

# Fixed MAD reference from combined pre-harmonisation PCs (stable; not influenced by post)
global_combined_pre_pc <- rbind(pc_pre_b1, pc_pre_b2)
global_mad_per_pc <- apply(global_combined_pre_pc, 2, function(x) mad(x, na.rm = TRUE))
global_mad_per_pc[global_mad_per_pc == 0 | is.na(global_mad_per_pc)] <- 1

d_mad_pre   <- mad_whitened_distance(pc_pre_b1,  pc_pre_b2,  sids_b1, sids_b2, mad_ref = global_mad_per_pc)
d_mad_post  <- mad_whitened_distance(pc_post_b1, pc_post_b2, sids_b1, sids_b2, mad_ref = global_mad_per_pc)

log_info("Euclidean median pre={round(median(d_eucl_pre),3)} post={round(median(d_eucl_post),3)}")
log_info("Mahalanobis median pre={round(median(d_maha_pre),3)} post={round(median(d_maha_post),3)}")
log_info("MAD-whitened median pre={round(median(d_mad_pre),3)} post={round(median(d_mad_post),3)}")

# Agreement on raw protein space (per-pair ICC/CCC)
compute_pair_agreement <- function(mat_b1, mat_b2, sids_b1, sids_b2, proteins) {
  icc_vals <- rep(NA_real_, length(sids_b1))
  ccc_vals <- rep(NA_real_, length(sids_b1))
  for (i in seq_along(sids_b1)) {
    v1 <- as.numeric(mat_b1[sids_b1[i], proteins])
    v2 <- as.numeric(mat_b2[sids_b2[i], proteins])
    ok <- !is.na(v1) & !is.na(v2)
    if (sum(ok) < 25) next

    icc_vals[i] <- tryCatch(
      irr::icc(cbind(v1[ok], v2[ok]), model = "twoway", type = "consistency", unit = "single")$value,
      error = function(e) NA_real_
    )
    ccc_vals[i] <- tryCatch(
      DescTools::CCC(v1[ok], v2[ok])$rho.c$est,
      error = function(e) NA_real_
    )
  }
  list(icc = icc_vals, ccc = ccc_vals)
}

agree_pre  <- compute_pair_agreement(pre_b1,  pre_b2,  sids_b1, sids_b2, common_proteins)
agree_post <- compute_pair_agreement(post_b1, post_b2, sids_b1, sids_b2, common_proteins)

stopifnot_msg(length(agree_pre$icc) == n_bridge && length(agree_post$icc) == n_bridge, "ICC vector length mismatch")
stopifnot_msg(length(agree_pre$ccc) == n_bridge && length(agree_post$ccc) == n_bridge, "CCC vector length mismatch")

log_info("ICC median pre={round(median(agree_pre$icc, na.rm=TRUE),3)} post={round(median(agree_post$icc, na.rm=TRUE),3)}")
log_info("CCC median pre={round(median(agree_pre$ccc, na.rm=TRUE),3)} post={round(median(agree_post$ccc, na.rm=TRUE),3)}")

# Noise floor (within-batch replicates) computed in *same* fixed PC space and MAD-whitened
compute_noise_floor_mad <- function(wb_pairs, pc_scores, global_mad_per_pc) {
  if (nrow(wb_pairs) == 0) return(numeric(0))
  d <- rep(NA_real_, nrow(wb_pairs))
  for (i in seq_len(nrow(wb_pairs))) {
    s1 <- wb_pairs$sid1[i]; s2 <- wb_pairs$sid2[i]
    if (s1 %in% rownames(pc_scores) && s2 %in% rownames(pc_scores)) {
      diff_vec <- pc_scores[s1, , drop = FALSE] - pc_scores[s2, , drop = FALSE]
      diff_vec <- sweep(diff_vec, 2, global_mad_per_pc, "/")
      d[i] <- sqrt(sum(diff_vec^2, na.rm = TRUE))
    }
  }
  d
}

# global_mad_per_pc already computed above (fixed MAD reference from pre-combined PCs)

noise_dists_b1 <- compute_noise_floor_mad(wb_rep_b1, pc_pre_b1, global_mad_per_pc)
noise_dists_b2 <- compute_noise_floor_mad(wb_rep_b2, pc_pre_b2, global_mad_per_pc)
noise_floor_dists_mad <- c(noise_dists_b1, noise_dists_b2)
noise_floor_mad <- median(noise_floor_dists_mad, na.rm = TRUE)

noise_floor_ratio <- median(d_mad_post, na.rm = TRUE) / noise_floor_mad
log_info("Noise floor MAD median={round(noise_floor_mad,3)} | Ratio bridge_post/noise={round(noise_floor_ratio,3)}")

# Paired log-ratio (primary: MAD) — epsilon guards against log(0)
log_ratio_mad <- log((d_mad_post + 1e-12) / (d_mad_pre + 1e-12))
median_log_ratio_mad <- median(log_ratio_mad, na.rm = TRUE)

# Bootstrap CI (median log-ratio)
set.seed(42)
n_boot <- 10000
boot_medians_mad <- replicate(n_boot, {
  idx <- sample.int(n_bridge, replace = TRUE)
  median(log_ratio_mad[idx], na.rm = TRUE)
})
boot_ci_mad <- quantile(boot_medians_mad, probs = c(0.025, 0.975), na.rm = TRUE)

wilcox_result <- wilcox.test(d_mad_post, d_mad_pre, paired = TRUE, alternative = "less")
log_info("Paired log-ratio (MAD): median={round(median_log_ratio_mad,4)} CI=[{round(boot_ci_mad[1],4)},{round(boot_ci_mad[2],4)}]")
log_info("Wilcoxon paired p (MAD): {format.pval(wilcox_result$p.value, digits=4)}")

# Build per-pair table (strict length checks)
bridge_pair_dt <- data.table(
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
  log_ratio_mad = log_ratio_mad,
  noise_floor_ratio = d_mad_post / noise_floor_mad
)
stopifnot_msg(nrow(bridge_pair_dt) == n_bridge, "bridge_pair_dt row mismatch")

# ==============================================================================
# STEP 2: PAIR IDENTITY TEST (Rank-1 identification rate)
# ==============================================================================
log_info("Step 2: Pair identity test (rank-1 kNN rate)")

# Vectorised rank-of-true-pair using FNN::get.knnx (O(n log n) vs O(n²) manual loop).
compute_rank_of_true_pair <- function(pc_b1, pc_b2, sids_b1, sids_b2) {
  ref <- pc_b2
  ref_ids <- rownames(ref)
  true_idx <- match(sids_b2, ref_ids)
  stopifnot_msg(!any(is.na(true_idx)), "Some true B2 pair SampleIDs not found in B2 PC matrix rownames.")

  query <- pc_b1[sids_b1, , drop = FALSE]

  # FNN returns nn.index: for each query row, the indices of k nearest neighbours in ref
  knn_res <- FNN::get.knnx(ref, query, k = nrow(ref))

  # For each query, find the rank position of its true pair
  ranks <- vapply(seq_len(nrow(query)), function(i) {
    which(knn_res$nn.index[i, ] == true_idx[i])
  }, integer(1))

  list(
    ranks = ranks,
    rank1_hits = sum(ranks == 1),
    rank1_rate = mean(ranks == 1),
    median_rank = median(ranks)
  )
}

id_pre  <- compute_rank_of_true_pair(pc_pre_b1,  pc_pre_b2,  sids_b1, sids_b2)
id_post <- compute_rank_of_true_pair(pc_post_b1, pc_post_b2, sids_b1, sids_b2)

bridge_pair_dt[, rank_pre := id_pre$ranks]
bridge_pair_dt[, rank_post := id_post$ranks]

log_info("Rank-1 ID rate pre: {id_pre$rank1_hits}/{n_bridge} ({round(id_pre$rank1_rate*100,1)}%)")
log_info("Rank-1 ID rate post: {id_post$rank1_hits}/{n_bridge} ({round(id_post$rank1_rate*100,1)}%)")
log_info("Median rank pre={id_pre$median_rank} post={id_post$median_rank}")

# ==============================================================================
# STEP 3: BATCH SEPARABILITY (full cohort)
# ==============================================================================
log_info("Step 3: Batch separability metrics (full cohort, fixed PC basis)")

build_cohort <- function(pc_b1, pc_b2, bid1, bid2) {
  scores <- rbind(pc_b1, pc_b2)
  labels <- c(rep(bid1, nrow(pc_b1)), rep(bid2, nrow(pc_b2)))
  list(scores = scores, labels = labels)
}

cohort_pre  <- build_cohort(pc_pre_b1,  pc_pre_b2,  batch1_id, batch2_id)
cohort_post <- build_cohort(pc_post_b1, pc_post_b2, batch1_id, batch2_id)

# Post-harmonisation PCA basis (for Page 5 right panel and variance-explained plot)
pca_basis_post <- build_fixed_pca(post_b1, post_b2, common_proteins,
                                  max_pcs = kpi_max_pcs, var_threshold = kpi_var_threshold,
                                  impute_max_frac = pca_impute_max_frac,
                                  min_complete_rows = pca_min_complete_rows)
proj_post_space_b1  <- project_onto_basis(post_b1, pca_basis_post, impute_max_frac = pca_impute_max_frac)
proj_post_space_b2  <- project_onto_basis(post_b2, pca_basis_post, impute_max_frac = pca_impute_max_frac)
pc_post_space_b1    <- proj_post_space_b1$scores
pc_post_space_b2    <- proj_post_space_b2$scores
cohort_post_space   <- build_cohort(pc_post_space_b1, pc_post_space_b2, batch1_id, batch2_id)

compute_silhouette <- function(scores, labels) {
  X <- scores[, 1:2, drop = FALSE]
  d <- dist(X)
  si <- cluster::silhouette(as.integer(factor(labels)), d)
  mean(si[, "sil_width"])
}

compute_kbet_acceptance <- function(scores, labels, k0 = 25) {
  k0 <- min(k0, floor(nrow(scores) * 0.1))
  res <- tryCatch(
    kBET::kBET(scores, labels, k0 = k0, plot = FALSE, verbose = FALSE),
    error = function(e) {
      log_warn("kBET failed: {e$message}")
      list(summary = data.frame(kBET.observed = NA_real_))
    }
  )
  1 - res$summary$kBET.observed[1]
}

compute_lisi_median <- function(scores, labels) {
  meta <- data.frame(batch = labels, row.names = rownames(scores))
  res <- tryCatch(
    lisi::compute_lisi(scores, meta, "batch"),
    error = function(e) {
      log_warn("LISI failed: {e$message}")
      data.frame(batch = NA_real_)
    }
  )
  median(res$batch, na.rm = TRUE)
}

compute_batch_auc <- function(scores, labels, max_pcs = 10) {
  if (length(unique(labels)) < 2) return(NA_real_)
  X <- as.data.frame(scores[, seq_len(min(max_pcs, ncol(scores))), drop = FALSE])
  y <- as.integer(factor(labels)) - 1L
  df <- cbind(X, batch = y)
  df <- df[complete.cases(df), , drop = FALSE]
  if (nrow(df) < 100) return(NA_real_)
  fit <- tryCatch(glm(batch ~ ., data = df, family = binomial()), error = function(e) NULL)
  if (is.null(fit)) return(NA_real_)
  p <- predict(fit, type = "response")
  as.numeric(tryCatch(pROC::auc(df$batch, p, quiet = TRUE), error = function(e) NA_real_))
}

sil_pre  <- compute_silhouette(cohort_pre$scores,  cohort_pre$labels)
sil_post <- compute_silhouette(cohort_post$scores, cohort_post$labels)

kbet_pre  <- compute_kbet_acceptance(cohort_pre$scores,  cohort_pre$labels)
kbet_post <- compute_kbet_acceptance(cohort_post$scores, cohort_post$labels)

lisi_pre  <- compute_lisi_median(cohort_pre$scores,  cohort_pre$labels)
lisi_post <- compute_lisi_median(cohort_post$scores, cohort_post$labels)

auc_pre  <- compute_batch_auc(cohort_pre$scores,  cohort_pre$labels)
auc_post <- compute_batch_auc(cohort_post$scores, cohort_post$labels)

log_info("Silhouette (PC1-2): pre={round(sil_pre,4)} post={round(sil_post,4)}")
log_info("kBET acceptance: pre={round(kbet_pre,4)} post={round(kbet_post,4)}")
log_info("LISI median: pre={round(lisi_pre,4)} post={round(lisi_post,4)}")
log_info("Batch AUC (PC1-10): pre={round(auc_pre,4)} post={round(auc_post,4)}")

# ==============================================================================
# STEP 4: BIOLOGY PRESERVATION (guard against over-correction)
# ==============================================================================
log_info("Step 4: Biology preservation checks")

sex_pred_b1_path <- get_output_path("04", "sex_predictions", batch1_id, "outliers", "tsv", config = config)
sex_pred_b2_path <- get_output_path("04", "sex_predictions", batch2_id, "outliers", "tsv", config = config)
sex_accuracy_pre <- NA_real_
sex_accuracy_post <- NA_real_
var_decomp_pre <- NULL
var_decomp_post <- NULL

# Helper: standardize sex prediction tables
standardize_sex_dt <- function(dt) {
  dt <- as.data.table(dt)
  # Accept common ID variants
  id_col <- intersect(names(dt), c("SAMPLE_ID", "SampleID", "sample_id", "sampleid"))
  stopifnot_msg(length(id_col) >= 1, "Sex prediction file must contain SAMPLE_ID (or SampleID) column.")
  id_col <- id_col[1]
  setnames(dt, id_col, "SampleID")
  stopifnot_msg("genetic_sex" %in% names(dt), "Sex prediction file must contain genetic_sex column.")
  dt[, genetic_sex := tolower(genetic_sex)]
  dt[genetic_sex %in% c("male", "female")][, .(SampleID, genetic_sex)]
}

# Sex accuracy: simple logistic model on top sex-associated proteins (if available)
compute_sex_accuracy <- function(mat_b1, mat_b2, map_b1, map_b2, sp_b1, sp_b2, proteins) {
  build_one <- function(mat, mapping, sp_dt) {
    sids <- intersect(rownames(mat), mapping$SampleID)
    m <- mapping[SampleID %in% sids, .(SampleID)]
    m <- unique(m)
    d <- merge(m, sp_dt, by = "SampleID", all.x = TRUE)
    d <- d[!is.na(genetic_sex)]
    if (nrow(d) < 50) return(NULL)
    X <- mat[d$SampleID, proteins, drop = FALSE]
    # complete cases
    ok <- complete.cases(X)
    X <- X[ok, , drop = FALSE]
    y <- as.integer(d$genetic_sex[ok] == "female")
    if (length(unique(y)) < 2) return(NULL)
    list(X = X, y = y)
  }

  d1 <- build_one(mat_b1, map_b1, sp_b1)
  d2 <- build_one(mat_b2, map_b2, sp_b2)
  if (is.null(d1) || is.null(d2)) return(NA_real_)

  X <- rbind(d1$X, d2$X)
  y <- c(d1$y, d2$y)
  if (nrow(X) < 100) return(NA_real_)
  fit <- tryCatch(glm(y ~ ., data = as.data.frame(X), family = binomial()), error = function(e) NULL)
  if (is.null(fit)) return(NA_real_)
  pred <- as.integer(predict(fit, type = "response") >= 0.5)
  mean(pred == y)
}

compute_variance_decomposition <- function(mat_b1, mat_b2, map_b1, map_b2, sp_b1, sp_b2, proteins, bid1, bid2) {
  build_labeled <- function(mat, mapping, sp_dt, bid) {
    sids <- intersect(rownames(mat), mapping$SampleID)
    m <- mapping[SampleID %in% sids, .(SampleID)]
    m <- unique(m)
    d <- merge(m, sp_dt, by = "SampleID", all.x = TRUE)
    d <- d[!is.na(genetic_sex)]
    if (nrow(d) < 50) return(NULL)
    X <- mat[d$SampleID, proteins, drop = FALSE]
    dt <- data.table(
      SampleID = d$SampleID,
      batch = bid,
      sex = d$genetic_sex
    )
    cbind(dt, as.data.table(X))
  }

  d1 <- build_labeled(mat_b1, map_b1, sp_b1, bid1)
  d2 <- build_labeled(mat_b2, map_b2, sp_b2, bid2)
  if (is.null(d1) || is.null(d2)) return(NULL)
  combined <- rbind(d1, d2, fill = TRUE)

  rbindlist(lapply(proteins, function(p) {
    y <- combined[[p]]
    ok <- !is.na(y)
    if (sum(ok) < 50) return(NULL)
    fit <- tryCatch(anova(lm(y[ok] ~ factor(combined$batch[ok]) + factor(combined$sex[ok]))),
                    error = function(e) NULL)
    if (is.null(fit)) return(NULL)
    ss <- fit[["Sum Sq"]]
    total <- sum(ss)
    data.table(
      protein = p,
      batch_var_frac = ss[1] / total,
      sex_var_frac = ss[2] / total,
      residual_frac = ss[3] / total
    )
  }))
}

if (file.exists(sex_pred_b1_path) && file.exists(sex_pred_b2_path)) {
  sp_b1 <- standardize_sex_dt(fread(sex_pred_b1_path))
  sp_b2 <- standardize_sex_dt(fread(sex_pred_b2_path))

  sex_prot_path <- get_output_path("04", "sex_associated_proteins", batch1_id, "outliers", "tsv", config = config)
  top_sex_proteins <- character(0)
  if (file.exists(sex_prot_path)) {
    sex_prots <- fread(sex_prot_path)
    if ("protein" %in% names(sex_prots)) {
      top_sex_proteins <- intersect(head(sex_prots$protein, 20), common_proteins)
    }
  }

  if (length(top_sex_proteins) >= 5) {
    sex_accuracy_pre <- compute_sex_accuracy(pre_b1, pre_b2, map_b1, map_b2, sp_b1, sp_b2, top_sex_proteins)
    sex_accuracy_post <- compute_sex_accuracy(post_b1, post_b2, map_b1, map_b2, sp_b1, sp_b2, top_sex_proteins)
    log_info("Sex accuracy pre={round(sex_accuracy_pre,4)} post={round(sex_accuracy_post,4)}")
  } else {
    log_warn("Insufficient sex-associated proteins; skipping sex accuracy check.")
  }

  # Variance decomposition on subset
  set.seed(42)
  var_decomp_proteins <- if (length(common_proteins) > 200) sample(common_proteins, 200) else common_proteins

  var_decomp_pre <- compute_variance_decomposition(pre_b1, pre_b2, map_b1, map_b2, sp_b1, sp_b2,
                                                  var_decomp_proteins, batch1_id, batch2_id)
  var_decomp_post <- compute_variance_decomposition(post_b1, post_b2, map_b1, map_b2, sp_b1, sp_b2,
                                                   var_decomp_proteins, batch1_id, batch2_id)
  if (!is.null(var_decomp_pre) && !is.null(var_decomp_post)) {
    log_info("Variance decomposition medians: batch pre={round(median(var_decomp_pre$batch_var_frac, na.rm=TRUE),4)} post={round(median(var_decomp_post$batch_var_frac, na.rm=TRUE),4)}")
    log_info("Variance decomposition medians: sex   pre={round(median(var_decomp_pre$sex_var_frac, na.rm=TRUE),4)} post={round(median(var_decomp_post$sex_var_frac, na.rm=TRUE),4)}")

    sex_drop <- (median(var_decomp_pre$sex_var_frac, na.rm=TRUE) - median(var_decomp_post$sex_var_frac, na.rm=TRUE)) /
      median(var_decomp_pre$sex_var_frac, na.rm=TRUE)
    if (!is.na(sex_drop) && sex_drop > 0.20) {
      log_warn("OVER-CORRECTION alarm: sex variance dropped by {round(sex_drop*100,1)}% post-harmonisation.")
    }
  }
} else {
  log_warn("Sex prediction files missing; skipping biology-preservation checks that require sex labels.")
}

# ==============================================================================
# STEP 5: OUTLIER DETECTION
# ==============================================================================
log_info("Step 5: Outlier detection")

# 5a) Bridge pair outliers using Mahalanobis (post) + noise ratio
q1 <- quantile(d_maha_post, 0.25, na.rm = TRUE)
q3 <- quantile(d_maha_post, 0.75, na.rm = TRUE)
iqr_val <- q3 - q1
outlier_threshold <- q3 + 1.5 * iqr_val

bridge_pair_dt[, outlier_maha := maha_post > outlier_threshold]
bridge_pair_dt[, outlier_noise := noise_floor_ratio > 3.0]
bridge_pair_dt[, outlier_flag := outlier_maha | outlier_noise]

n_outliers <- sum(bridge_pair_dt$outlier_flag, na.rm = TRUE)
log_info("Flagged bridge outliers: {n_outliers}/{n_bridge}")

if (n_outliers > 0) {
  log_warn("Flagged FINNGENIDs: {paste(bridge_pair_dt[outlier_flag==TRUE]$finngenid, collapse=', ')}")
}

# 5b) Protein concordance profiling (MSE/MAE/R2 across bridge pairs)
log_info("Computing per-protein bridge concordance (MSE/MAE/R²)")

protein_concordance <- rbindlist(lapply(common_proteins, function(p) {
  v1_pre  <- pre_b1[sids_b1,  p]
  v2_pre  <- pre_b2[sids_b2,  p]
  v1_post <- post_b1[sids_b1, p]
  v2_post <- post_b2[sids_b2, p]

  ok_pre  <- !is.na(v1_pre) & !is.na(v2_pre)
  ok_post <- !is.na(v1_post) & !is.na(v2_post)

  cor_pre  <- if (sum(ok_pre)  >= 5) cor(v1_pre[ok_pre],  v2_pre[ok_pre])  else NA_real_
  cor_post <- if (sum(ok_post) >= 5) cor(v1_post[ok_post], v2_post[ok_post]) else NA_real_

  mse_pre  <- if (sum(ok_pre)  >= 5) mean((v1_pre[ok_pre]  - v2_pre[ok_pre])^2) else NA_real_
  mse_post <- if (sum(ok_post) >= 5) mean((v1_post[ok_post] - v2_post[ok_post])^2) else NA_real_

  mae_pre  <- if (sum(ok_pre)  >= 5) mean(abs(v1_pre[ok_pre]  - v2_pre[ok_pre])) else NA_real_
  mae_post <- if (sum(ok_post) >= 5) mean(abs(v1_post[ok_post] - v2_post[ok_post])) else NA_real_

  r2_pre <- if (sum(ok_pre) >= 5) {
    ss_res <- sum((v2_pre[ok_pre] - v1_pre[ok_pre])^2)
    ss_tot <- sum((v2_pre[ok_pre] - mean(v2_pre[ok_pre]))^2)
    if (ss_tot > 0) 1 - ss_res/ss_tot else NA_real_
  } else NA_real_

  r2_post <- if (sum(ok_post) >= 5) {
    ss_res <- sum((v2_post[ok_post] - v1_post[ok_post])^2)
    ss_tot <- sum((v2_post[ok_post] - mean(v2_post[ok_post]))^2)
    if (ss_tot > 0) 1 - ss_res/ss_tot else NA_real_
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

protein_concordance[, poorly_harmonised := (
  (!is.na(mse_reduction) & mse_reduction < 0) |
  (!is.na(r2_post) & r2_post < 0.5)
)]
n_poor <- sum(protein_concordance$poorly_harmonised, na.rm = TRUE)
log_info("Poorly harmonised proteins: {n_poor}/{length(common_proteins)}")

# ==============================================================================
# STEP 6: KPI SUMMARY TABLE + PDF DASHBOARD
# ==============================================================================
log_info("Step 6: KPI summary + dashboard")

kpi_dt <- data.table(
  metric = c(
    "Median bridge Euclidean distance",
    "Median bridge Mahalanobis distance",
    "Median bridge MAD-whitened distance (PRIMARY)",
    "Paired log-ratio median (MAD-whitened)",
    "Paired log-ratio 95% CI lower (MAD-whitened)",
    "Paired log-ratio 95% CI upper (MAD-whitened)",
    "CI excludes zero (MAD-whitened)",
    "Noise floor (within-batch replicate dist, MAD)",
    "Noise floor ratio (bridge_post / noise)",
    "ICC (median across pairs)",
    "CCC (median across pairs)",
    "Rank-1 ID rate",
    "Median rank of true pair",
    "Silhouette score (PC1-PC2)",
    "kBET acceptance rate",
    "LISI batch mixing score",
    "Batch separability AUC (PC1-10)",
    "Sex prediction accuracy",
    "Batch variance fraction (median)",
    "Sex variance fraction (median)",
    "Flagged bridge outliers",
    "Poorly harmonised proteins (MSE↑ or R²<0.5)",
    "Wilcoxon signed-rank p-value (MAD, post<pre)"
  ),
  pre_harmonisation = c(
    round(median(d_eucl_pre), 4),
    round(median(d_maha_pre), 4),
    round(median(d_mad_pre), 4),
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
    round(auc_pre, 4),
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
    round(auc_post, 4),
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
    "<0",
    "1 (TRUE)",
    "reference",
    "~1.0",
    ">0.9",
    ">0.9",
    ">0.9",
    "1",
    "decrease toward 0",
    ">0.8",
    "close to 2.0",
    "toward 0.5",
    "stable",
    "decrease",
    "stable (not drop >20%)",
    "0 ideal",
    "0 ideal",
    "<0.05"
  )
)

# Pass/Warn/Fail rules (same spirit, slightly clarified)
kpi_dt[, pass := NA_character_]
kpi_dt[metric %in% c("Median bridge Euclidean distance",
                     "Median bridge Mahalanobis distance",
                     "Median bridge MAD-whitened distance (PRIMARY)"),
       pass := ifelse(post_harmonisation < pre_harmonisation, "PASS", "FAIL")]
kpi_dt[metric == "CI excludes zero (MAD-whitened)",
       pass := ifelse(post_harmonisation == 1, "PASS", "FAIL")]
kpi_dt[metric == "Noise floor ratio (bridge_post / noise)",
       pass := ifelse(!is.na(post_harmonisation) & post_harmonisation < 2.0, "PASS",
                      ifelse(!is.na(post_harmonisation) & post_harmonisation < 3.0, "WARN", "FAIL"))]
kpi_dt[metric %in% c("ICC (median across pairs)", "CCC (median across pairs)", "Rank-1 ID rate"),
       pass := ifelse(!is.na(post_harmonisation) & post_harmonisation > 0.9, "PASS",
                      ifelse(!is.na(post_harmonisation) & post_harmonisation > 0.7, "WARN", "FAIL"))]
kpi_dt[metric == "kBET acceptance rate",
       pass := ifelse(!is.na(post_harmonisation) & post_harmonisation > 0.8, "PASS",
                      ifelse(!is.na(post_harmonisation) & post_harmonisation > 0.5, "WARN", "FAIL"))]
kpi_dt[metric == "LISI batch mixing score",
       pass := ifelse(!is.na(post_harmonisation) & post_harmonisation > 1.8, "PASS",
                      ifelse(!is.na(post_harmonisation) & post_harmonisation > 1.5, "WARN", "FAIL"))]
kpi_dt[metric == "Wilcoxon signed-rank p-value (MAD, post<pre)",
       pass := ifelse(!is.na(post_harmonisation) & post_harmonisation < 0.05, "PASS", "FAIL")]

# ==============================================================================
# SAVE OUTPUTS
# ==============================================================================
log_info("Saving outputs")

output_batch_id <- batch1_id

kpi_path <- get_output_path("07b", "kpi_summary", output_batch_id, "normalized", "tsv", config = config)
ensure_output_dir(kpi_path)
fwrite(kpi_dt, kpi_path, sep = "\t")
log_info("Saved KPI summary: {kpi_path}")

pair_dist_path <- get_output_path("07b", "bridge_pair_distances", output_batch_id, "normalized", "tsv", config = config)
ensure_output_dir(pair_dist_path)
fwrite(bridge_pair_dt, pair_dist_path, sep = "\t")
log_info("Saved bridge pair distances: {pair_dist_path}")

noise_path <- get_output_path("07b", "noise_floor_analysis", output_batch_id, "normalized", "tsv", config = config)
ensure_output_dir(noise_path)
noise_dt <- data.table(
  source = c(rep("B1_within_batch", length(noise_dists_b1)),
             rep("B2_within_batch", length(noise_dists_b2))),
  distance_mad = c(noise_dists_b1, noise_dists_b2)
)
fwrite(noise_dt, noise_path, sep = "\t")
log_info("Saved noise floor analysis: {noise_path}")

conc_path <- get_output_path("07b", "protein_concordance", output_batch_id, "normalized", "tsv", config = config)
ensure_output_dir(conc_path)
fwrite(protein_concordance, conc_path, sep = "\t")
log_info("Saved protein concordance: {conc_path}")

poorly_harmonised <- protein_concordance[poorly_harmonised == TRUE]
if (nrow(poorly_harmonised) > 0) {
  poor_path <- get_output_path("07b", "poorly_harmonised_proteins", output_batch_id, "normalized", "tsv", config = config)
  ensure_output_dir(poor_path)
  poorly_harmonised <- poorly_harmonised[order(r2_post, na.last = TRUE)]
  fwrite(poorly_harmonised[, .(protein, cor_pre, cor_post, mse_pre, mse_post, mse_reduction, mae_pre, mae_post, r2_pre, r2_post)],
         poor_path, sep = "\t")
  log_info("Saved poorly harmonised proteins: {poor_path} (n={nrow(poorly_harmonised)})")
} else {
  log_info("No poorly harmonised proteins to save.")
}

outlier_path <- get_output_path("07b", "outlier_flags", output_batch_id, "normalized", "tsv", config = config)
ensure_output_dir(outlier_path)
fwrite(bridge_pair_dt[, .(finngenid, sid_b1, sid_b2, mad_post, maha_post, noise_floor_ratio, outlier_maha, outlier_noise, outlier_flag, rank_post)],
       outlier_path, sep = "\t")
log_info("Saved outlier flags: {outlier_path}")

if (!is.null(var_decomp_pre) && !is.null(var_decomp_post)) {
  vd_combined <- merge(var_decomp_pre, var_decomp_post, by = "protein", suffixes = c("_pre", "_post"))
  vd_path <- get_output_path("07b", "variance_decomposition", output_batch_id, "normalized", "tsv", config = config)
  ensure_output_dir(vd_path)
  fwrite(vd_combined, vd_path, sep = "\t")
  log_info("Saved variance decomposition: {vd_path}")
}

# ==============================================================================
# ADDITIONAL COMPUTATIONS (for enhanced dashboard)
# ==============================================================================
log_info("Computing additional metrics for enhanced dashboard")

fg_short <- substr(common_bridge_fg, 1, 10)
outlier_fg_short <- fg_short[bridge_pair_dt$outlier_flag == TRUE]

# Full NN distance matrix (MAD-whitened, bridge B1 x bridge B2)
compute_nn_heatmap <- function(pc_b1, pc_b2, s_b1, s_b2, mad_scales) {
  dist_mat <- matrix(NA_real_, nrow = length(s_b1), ncol = length(s_b2))
  for (i in seq_along(s_b1)) {
    for (j in seq_along(s_b2)) {
      if (s_b1[i] %in% rownames(pc_b1) && s_b2[j] %in% rownames(pc_b2)) {
        diff_vec <- pc_b1[s_b1[i], , drop = FALSE] - pc_b2[s_b2[j], , drop = FALSE]
        dist_mat[i, j] <- sqrt(sum((diff_vec / mad_scales)^2, na.rm = TRUE))
      }
    }
  }
  dist_mat
}

nn_dist_mat_post <- compute_nn_heatmap(pc_post_b1, pc_post_b2, sids_b1, sids_b2, global_mad_per_pc)
rownames(nn_dist_mat_post) <- fg_short
colnames(nn_dist_mat_post) <- fg_short
log_info("NN distance matrix: {nrow(nn_dist_mat_post)}x{ncol(nn_dist_mat_post)}")

# Pairing margins (confidence of correct match)
compute_pairing_margins <- function(dist_mat) {
  margins <- numeric(nrow(dist_mat))
  for (i in seq_len(nrow(dist_mat))) {
    margins[i] <- dist_mat[i, i] - min(dist_mat[i, -i], na.rm = TRUE)
  }
  margins
}

pairing_margins <- compute_pairing_margins(nn_dist_mat_post)
log_info("Pairing margins: {sum(pairing_margins < 0)} negative (swap candidates) out of {length(pairing_margins)}")

# Bridge protein correlation matrices (Pearson r)
compute_bridge_protein_corr <- function(mat_b1, mat_b2, s_b1, s_b2, proteins, labels) {
  n <- length(s_b1)
  corr_mat <- matrix(NA_real_, nrow = n, ncol = n)
  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      if (s_b1[i] %in% rownames(mat_b1) && s_b2[j] %in% rownames(mat_b2)) {
        x <- as.numeric(mat_b1[s_b1[i], proteins])
        y <- as.numeric(mat_b2[s_b2[j], proteins])
        ok <- !is.na(x) & !is.na(y)
        if (sum(ok) > 10) corr_mat[i, j] <- cor(x[ok], y[ok], method = "pearson")
      }
    }
  }
  rownames(corr_mat) <- labels
  colnames(corr_mat) <- labels
  corr_mat
}

bridge_corr_pre  <- compute_bridge_protein_corr(pre_b1, pre_b2, sids_b1, sids_b2, common_proteins, fg_short)
bridge_corr_post <- compute_bridge_protein_corr(post_b1, post_b2, sids_b1, sids_b2, common_proteins, fg_short)
diag_corr_pre  <- diag(bridge_corr_pre)
diag_corr_post <- diag(bridge_corr_post)
log_info("Bridge protein correlation: pre diagonal median={round(median(diag_corr_pre, na.rm=TRUE),3)}, post diagonal median={round(median(diag_corr_post, na.rm=TRUE),3)}")

# PC-by-batch effect size (R-squared per PC)
compute_pc_batch_effect <- function(pc_scores, labels, max_pcs = 10) {
  n_pcs <- min(max_pcs, ncol(pc_scores))
  r2_vals <- numeric(n_pcs)
  for (j in seq_len(n_pcs)) {
    fit <- tryCatch(lm(pc_scores[, j] ~ factor(labels)), error = function(e) NULL)
    r2_vals[j] <- if (!is.null(fit)) summary(fit)$r.squared else NA_real_
  }
  r2_vals
}

pc_r2_pre  <- compute_pc_batch_effect(cohort_pre$scores, cohort_pre$labels)
pc_r2_post <- compute_pc_batch_effect(cohort_post$scores, cohort_post$labels)
log_info("PC batch effect R-squared: PC1 pre={round(pc_r2_pre[1],4)} post={round(pc_r2_post[1],4)}")

# ==============================================================================
# DASHBOARD PLOTS
# ==============================================================================
log_info("Generating dashboard PDF")

dashboard_path <- get_output_path("07b", "kpi_dashboard", output_batch_id, "normalized", "pdf", config = config)
ensure_output_dir(dashboard_path)

# Plot 1: Paired distance box/paired lines (MAD)
boxplot_df <- data.frame(
  distance = c(d_mad_pre, d_mad_post),
  stage = factor(rep(c("Pre-harmonisation", "Post-harmonisation"), each = n_bridge),
                 levels = c("Pre-harmonisation", "Post-harmonisation")),
  pair = rep(seq_len(n_bridge), 2)
)
p_box <- ggpaired(
  boxplot_df, x = "stage", y = "distance", id = "pair",
  fill = "stage", palette = c(PAL$pre, PAL$post),
  line.color = "grey70", line.size = 0.4,
  point.size = 2.2, width = 0.5
) +
  stat_compare_means(method = "wilcox.test", paired = TRUE, label = "p.format", size = 4) +
  labs(
    title = "Bridge Pair Distance (PRIMARY: MAD-whitened, Fixed PCA Basis)",
    subtitle = paste0(
      "Each point is a bridge sample measured in both batches. Lines connect the same individual.\n",
      "A downward shift indicates improved cross-batch agreement after harmonisation.\n",
      sprintf("Wilcoxon paired test (post < pre): p = %s | n = %d bridge pairs",
              format.pval(wilcox_result$p.value, digits = 3), n_bridge)),
    y = "MAD-whitened distance",
    x = NULL
  ) +
  theme(legend.position = "none")

# Plot 2: Dumbbell plot sorted by improvement
dumbbell_dt <- data.table(
  finngenid = common_bridge_fg,
  d_pre = d_mad_pre,
  d_post = d_mad_post
)
dumbbell_dt[, improvement := d_post / d_pre]
setorder(dumbbell_dt, improvement)
dumbbell_dt[, idx := seq_len(.N)]
dumbbell_dt[, is_outlier := finngenid %in% bridge_pair_dt[outlier_flag==TRUE]$finngenid]
db <- as.data.frame(dumbbell_dt)
db$fin_short <- substr(db$finngenid, 1, 10)

p_dumbbell <- ggplot(db, aes(x = idx)) +
  geom_segment(aes(y = d_pre, yend = d_post, xend = idx, colour = is_outlier), linewidth = 1.1, alpha = 0.7) +
  geom_point(aes(y = d_pre), shape = 21, size = 3, fill = PAL$pre, color = "black", alpha = 0.85) +
  geom_point(aes(y = d_post), shape = 21, size = 3, fill = PAL$post, color = "black", alpha = 0.85) +
  geom_text(aes(y = d_post, label = fin_short), angle = 45, hjust = 1, vjust = 1.4, size = 2.4, color = "grey25") +
  scale_colour_manual(values = c(`FALSE` = PAL$neutral, `TRUE` = PAL$accent), guide = "none") +
  labs(
    title = "Bridge Pair Improvement (MAD distance)",
    subtitle = paste0(
      "Each segment connects pre- (orange) and post-harmonisation (teal) distance for one bridge pair.\n",
      "Pairs sorted by improvement ratio D_post/D_pre (left = most improved). Red = flagged outliers.\n",
      sprintf("n = %d bridge pairs | Flagged outliers: %d", n_bridge, sum(dumbbell_dt$is_outlier))),
    x = "Bridge pairs (sorted)",
    y = "MAD-whitened distance"
  ) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

# Plot 3: ECDF overlay pre vs post (enhanced with percentiles)
ecdf_df <- rbind(
  data.frame(stage = "Pre", distance = d_mad_pre),
  data.frame(stage = "Post", distance = d_mad_post)
)
ecdf_df$stage <- factor(ecdf_df$stage, levels = c("Pre", "Post"))
p90_pre  <- quantile(d_mad_pre, 0.90, na.rm = TRUE)
p95_pre  <- quantile(d_mad_pre, 0.95, na.rm = TRUE)
p90_post <- quantile(d_mad_post, 0.90, na.rm = TRUE)
p95_post <- quantile(d_mad_post, 0.95, na.rm = TRUE)
p_ecdf <- ggplot(ecdf_df, aes(x = distance, colour = stage)) +
  stat_ecdf(linewidth = 1.2) +
  geom_vline(xintercept = median(d_mad_pre),  linetype = "dashed", colour = PAL$pre,  alpha = 0.6) +
  geom_vline(xintercept = median(d_mad_post), linetype = "dashed", colour = PAL$post, alpha = 0.6) +
  geom_vline(xintercept = p90_pre,  linetype = "dotted", colour = PAL$pre,  alpha = 0.5) +
  geom_vline(xintercept = p95_pre,  linetype = "dotted", colour = PAL$pre,  alpha = 0.5) +
  geom_vline(xintercept = p90_post, linetype = "dotted", colour = PAL$post, alpha = 0.5) +
  geom_vline(xintercept = p95_post, linetype = "dotted", colour = PAL$post, alpha = 0.5) +
  annotate("text", x = median(d_mad_pre),  y = 0.5, label = "Pre median",
           hjust = -0.1, size = 3, colour = PAL$pre) +
  annotate("text", x = median(d_mad_post), y = 0.5, label = "Post median",
           hjust = -0.1, size = 3, colour = PAL$post) +
  scale_colour_manual(values = c(Pre = PAL$pre, Post = PAL$post)) +
  labs(
    title = "Distribution Overlay (ECDF): Bridge Distances",
    subtitle = paste0(
      "Empirical cumulative distribution of MAD-whitened bridge pair distances.\n",
      "A leftward shift of the post curve indicates distance reduction across the full distribution.\n",
      sprintf("Pre: median=%.2f, p90=%.2f, p95=%.2f | Post: median=%.2f, p90=%.2f, p95=%.2f",
              median(d_mad_pre), p90_pre, p95_pre, median(d_mad_post), p90_post, p95_post)),
    x = "MAD-whitened distance",
    y = "Cumulative probability"
  ) +
  theme(legend.position = "bottom", legend.title = element_blank())

# Plot 4: Log-ratio histogram
lr_df <- data.frame(lr = log_ratio_mad)
p_lr <- ggplot(lr_df, aes(x = lr)) +
  annotate("rect", xmin = boot_ci_mad[1], xmax = boot_ci_mad[2],
           ymin = -Inf, ymax = Inf, fill = "#1119B0", alpha = 0.45) +
  geom_histogram(bins = 15, fill = PAL$post, alpha = 0.7, colour = "black", linewidth = 0.25) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "#d7191c", linewidth = 0.8) +
  geom_vline(xintercept = median_log_ratio_mad, colour = "#1A237E", linewidth = 1.0) +
  labs(
    title = "Paired Log-Ratio: log(D_post / D_pre)",
    subtitle = paste0(
      "Distribution of per-pair log-ratios: negative values mean post-distance < pre-distance.\n",
      "Blue band = bootstrap 95% CI of the median. Dashed red = zero (no change).\n",
      sprintf("Median = %.3f | 95%% CI [%.3f, %.3f] | CI entirely < 0 confirms improvement",
              median_log_ratio_mad, boot_ci_mad[1], boot_ci_mad[2])),
    x = "log(D_post / D_pre)",
    y = "Count"
  )

# Plot 5: Noise floor density
nf_df <- rbind(
  data.frame(group = "Bridge pairs (post)", distance = d_mad_post),
  data.frame(group = "Within-batch replicates (noise)", distance = noise_floor_dists_mad)
)
nf_df <- nf_df[is.finite(nf_df$distance), , drop = FALSE]
p_noise <- ggplot(nf_df, aes(x = distance, fill = group)) +
  geom_density(alpha = 0.45) +
  geom_vline(xintercept = noise_floor_mad, linetype = "dashed", linewidth = 0.8) +
  labs(
    title = "Bridge Distances vs Noise Floor (MAD units, consistent)",
    subtitle = paste0(
      "Noise floor = median MAD-whitened distance between within-batch technical replicates.\n",
      "Ratio near 1.0 means bridge pairs approach measurement noise — ideal harmonisation.\n",
      sprintf("Noise floor median = %.2f | Bridge post median = %.2f | Ratio = %.2f",
              noise_floor_mad, median(d_mad_post), noise_floor_ratio)),
    x = "MAD-whitened distance",
    y = "Density"
  ) +
  scale_fill_manual(values = c("Bridge pairs (post)" = PAL$post, "Within-batch replicates (noise)" = "grey60")) +
  theme(legend.position = "bottom", legend.title = element_blank())

# Plot 6: Rank-1 rate bars
rank_df <- data.frame(
  stage = factor(c("Pre", "Post"), levels = c("Pre", "Post")),
  rate = c(id_pre$rank1_rate, id_post$rank1_rate),
  hits = c(id_pre$rank1_hits, id_post$rank1_hits)
)
p_rank1 <- ggplot(rank_df, aes(x = stage, y = rate, fill = stage)) +
  geom_col(alpha = 0.85) +
  geom_text(aes(label = sprintf("%d/%d", hits, n_bridge)), vjust = -0.35, size = 3.8) +
  scale_fill_manual(values = c(Pre = PAL$pre, Post = PAL$post)) +
  scale_y_continuous(limits = c(0, 1.1), labels = scales::percent_format()) +
  labs(
    title = "Rank-1 Identification Rate",
    subtitle = paste0(
      "Fraction of bridge samples whose true cross-batch pair is the nearest neighbour.\n",
      "Higher is better; a rate > 90% indicates strong identity preservation.\n",
      sprintf("Pre: %d/%d (%.0f%%) | Post: %d/%d (%.0f%%)",
              id_pre$rank1_hits, n_bridge, id_pre$rank1_rate * 100,
              id_post$rank1_hits, n_bridge, id_post$rank1_rate * 100)),
    x = NULL,
    y = "Rank-1 hit rate"
  ) +
  theme(legend.position = "none")

# Plot 7: Separability metrics (incl AUC) — facet titles show actual pre/post values
sep_df <- data.frame(
  metric = rep(c("Silhouette (↓)", "kBET acceptance (↑)", "LISI (→2.0)", "Batch AUC (→0.5)"), each = 2),
  stage = rep(c("Pre", "Post"), times = 4),
  value = c(sil_pre, sil_post, kbet_pre, kbet_post, lisi_pre, lisi_post, auc_pre, auc_post)
)
sep_df$stage <- factor(sep_df$stage, levels = c("Pre", "Post"))
sep_facet_labels <- paste0(
  c("Silhouette", "kBET acceptance", "LISI", "Batch AUC"), " (",
  round(c(sil_pre, kbet_pre, lisi_pre, auc_pre), 3), " / ",
  round(c(sil_post, kbet_post, lisi_post, auc_post), 3), ")"
)
sep_df$facet_label <- factor(
  sep_df$metric,
  levels = c("Silhouette (↓)", "kBET acceptance (↑)", "LISI (→2.0)", "Batch AUC (→0.5)"),
  labels = sep_facet_labels
)
p_sep <- ggplot(sep_df, aes(x = stage, y = value, fill = stage)) +
  geom_col(alpha = 0.85) +
  facet_wrap(~facet_label, scales = "free_y", ncol = 2) +
  scale_fill_manual(values = c(Pre = PAL$pre, Post = PAL$post)) +
  labs(
    title = "Batch Separability Metrics",
    subtitle = paste0(
      "Four complementary views of batch separability in fixed PC space.\n",
      "Silhouette should decrease; kBET acceptance should increase; LISI should approach 2.0; AUC toward 0.5.\n",
      "Some residual separation post-harmonisation may reflect true biology rather than technical artefact."),
    x = NULL,
    y = "Value"
  ) +
  theme(legend.position = "none")

# Plot 8: Outlier barplot (post MAD)
out_df <- as.data.frame(bridge_pair_dt)
out_df$fin_short <- substr(out_df$finngenid, 1, 10)
p_out <- ggplot(out_df, aes(x = reorder(fin_short, -mad_post), y = mad_post, fill = outlier_flag)) +
  geom_col(alpha = 0.85) +
  scale_fill_manual(values = c(`FALSE` = PAL$post, `TRUE` = PAL$accent), labels = c("Normal", "Flagged")) +
  labs(
    title = "Flagged Bridge Pairs (Post-harmonisation)",
    subtitle = paste0(
      "Post-harmonisation MAD-whitened distance for each bridge pair, sorted by magnitude.\n",
      "Red bars indicate flagged outliers (Mahalanobis > Q3 + 1.5 x IQR or noise ratio > 3.0).\n",
      sprintf("Flagged: %d/%d bridge pairs | Outlier threshold (Mahalanobis): %.2f",
              n_outliers, n_bridge, outlier_threshold)),
    x = "Bridge pair (FINNGENID, truncated)",
    y = "MAD-whitened distance (post)"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        legend.position = "bottom",
        legend.title = element_blank())

# Variance decomposition plot (enhanced: high-variance protein annotations)
p_var <- NULL
if (!is.null(var_decomp_pre) && !is.null(var_decomp_post)) {
  vd_long <- rbind(
    data.frame(stage="Pre", component="Batch", fraction=var_decomp_pre$batch_var_frac, protein=var_decomp_pre$protein),
    data.frame(stage="Pre", component="Sex",  fraction=var_decomp_pre$sex_var_frac,   protein=var_decomp_pre$protein),
    data.frame(stage="Post", component="Batch", fraction=var_decomp_post$batch_var_frac, protein=var_decomp_post$protein),
    data.frame(stage="Post", component="Sex",  fraction=var_decomp_post$sex_var_frac,   protein=var_decomp_post$protein)
  )
  vd_long$stage <- factor(vd_long$stage, levels=c("Pre","Post"))
  vd_long$component <- factor(vd_long$component, levels=c("Batch","Sex"))

  # Identify high-variance proteins (>30%)
  vd_long$high_var <- vd_long$fraction > 0.3
  high_var_proteins_vd <- unique(vd_long$protein[vd_long$high_var])
  vd_low_var  <- subset(vd_long, !high_var)
  vd_high_var <- subset(vd_long, high_var)

  # Compute median percentages for box annotation
  vd_summary <- aggregate(fraction ~ stage + component, data = vd_long,
                          FUN = function(x) median(x, na.rm = TRUE))
  names(vd_summary)[names(vd_summary) == "fraction"] <- "median_pct"
  # Position median labels so they stay readable when near 85%: avoid top edge and overlap
  vd_summary$label_y <- ifelse(vd_summary$median_pct > 0.78,
                              pmin(vd_summary$median_pct - 0.02, 0.88),
                              pmin(vd_summary$median_pct + 0.04, 0.92))

  # Jitter high-variance points so they and their annotations don't overlap.
  # NOTE: ggboxplot creates a discrete x-axis, so jitter must stay compatible
  # with the factor scale. We jitter y only and use position_jitter for x nudge
  # inside the geom calls, with a fixed seed for reproducibility.
  if (nrow(vd_high_var) > 0) {
    set.seed(40721)
    jitter_y_amount <- 0.025
    vd_high_var$jitter_y <- vd_high_var$fraction + runif(nrow(vd_high_var), -jitter_y_amount, jitter_y_amount)
    vd_high_var$jitter_y <- pmax(0.01, pmin(0.99, vd_high_var$jitter_y))
  }

  # Build subtitle
  vd_subtitle_base <- paste0(
    "ANOVA decomposition of NPX ~ batch + sex per protein (200 proteins sampled randomly). Batch variance should decrease post-harmonisation;\n",
    "sex variance should remain stable (drop > 20% signals over-correction). Red triangles = > 30% variance.\n")
  if (length(high_var_proteins_vd) > 0) {
    hv_text <- paste(head(high_var_proteins_vd, 10), collapse = ", ")
    if (length(high_var_proteins_vd) > 10) hv_text <- paste0(hv_text, " ...")
    vd_subtitle <- paste0(vd_subtitle_base, sprintf("High-variance proteins (%d): %s",
                                                     length(high_var_proteins_vd), hv_text))
  } else {
    vd_subtitle <- paste0(vd_subtitle_base, "No high-variance proteins detected.")
  }

  p_var <- ggboxplot(vd_low_var, x = "component", y = "fraction",
                     fill = "stage", palette = c(PAL$pre, PAL$post),
                     add = "jitter", add.params = list(size = 0.8, alpha = 0.3),
                     outlier.shape = NA) +
    geom_text(data = vd_summary,
              aes(x = component, y = label_y, label = sprintf("%.1f%%", median_pct * 100)),
              position = position_dodge(width = 0.75), size = 3.5, fontface = "bold",
              vjust = -0.3) +
    {if (nrow(vd_high_var) > 0) {
      geom_point(data = vd_high_var, aes(x = component, y = jitter_y),
                 position = position_jitter(width = 0.18, height = 0, seed = 40721),
                 shape = 17, size = 3.5, colour = "#FF0000", alpha = 0.9)
    } else geom_blank()} +
    {if (nrow(vd_high_var) > 0 && requireNamespace("ggrepel", quietly = TRUE)) {
      ggrepel::geom_text_repel(data = vd_high_var,
                               aes(x = component, y = jitter_y, label = protein),
                               position = position_jitter(width = 0.18, height = 0, seed = 40721),
                               size = 2.5, colour = "grey30", fontface = "bold",
                               min.segment.length = 0, box.padding = 0.3,
                               point.padding = 0.2, max.overlaps = Inf, force = 2)
    } else if (nrow(vd_high_var) > 0) {
      geom_text(data = vd_high_var,
                aes(x = component, y = jitter_y + 0.05, label = protein),
                position = position_jitter(width = 0.18, height = 0, seed = 40721),
                size = 2.5, hjust = 0, vjust = 0, angle = 45, colour = "grey30",
                fontface = "bold")
    } else geom_blank()} +
    scale_y_continuous(labels = scales::percent_format()) +
    coord_cartesian(ylim = c(0, 1.0)) +
    facet_wrap(~stage, ncol = 2) +
    labs(title = "Variance Decomposition (NPX ~ batch + sex)",
         subtitle = vd_subtitle,
         y = "Fraction of variance", x = NULL, fill = "Stage") +
    theme(legend.position = "bottom",
          plot.subtitle = element_text(size = 9, lineheight = 1.2))
}

# ------------------------------------------------------------------------------
# NEW PLOTS: A1-A8, B1 chart
# ------------------------------------------------------------------------------

# A5: Spaghetti plot (paired distance reduction)
spaghetti_df <- data.frame(
  pair = rep(seq_len(n_bridge), 2),
  stage = factor(rep(c("Pre", "Post"), each = n_bridge), levels = c("Pre", "Post")),
  distance = c(d_mad_pre, d_mad_post)
)
p_spaghetti <- ggplot(spaghetti_df, aes(x = stage, y = distance, group = pair)) +
  geom_line(alpha = 0.4, colour = "grey40") +
  geom_point(aes(colour = stage), size = 2) +
  scale_colour_manual(values = c(Pre = PAL$pre, Post = PAL$post)) +
  labs(title = "Paired Distance Reduction (MAD-Whitened)",
       subtitle = paste0(
         "Each line traces one bridge pair from pre- to post-harmonisation distance.\n",
         "Downward-sloping lines indicate the pair moved closer after normalisation.\n",
         sprintf("n = %d pairs | Median reduction: %.2f -> %.2f",
                 n_bridge, median(d_mad_pre), median(d_mad_post))),
       y = "MAD-whitened distance", x = NULL) +
  theme(legend.position = "none")

# A1: PCA Scatter with Bridge Connectors
pca_build_df <- function(pc_b1, pc_b2, stage_label, bid1, bid2) {
  data.frame(
    PC1 = c(pc_b1[, 1], pc_b2[, 1]),
    PC2 = c(pc_b1[, 2], pc_b2[, 2]),
    batch = c(rep(bid1, nrow(pc_b1)), rep(bid2, nrow(pc_b2))),
    sample = c(rownames(pc_b1), rownames(pc_b2)),
    bridge = c(rownames(pc_b1), rownames(pc_b2)) %in% c(sids_b1, sids_b2),
    stage = stage_label
  )
}
pca_df <- rbind(
  pca_build_df(pc_pre_b1, pc_pre_b2, "Pre-harmonisation", batch1_id, batch2_id),
  pca_build_df(pc_post_b1, pc_post_b2, "Post-harmonisation", batch1_id, batch2_id)
)
pca_df$stage <- factor(pca_df$stage, levels = c("Pre-harmonisation", "Post-harmonisation"))

build_connectors <- function(pcb1, pcb2, sb1, sb2, stage_label) {
  rows <- lapply(seq_along(sb1), function(i) {
    if (sb1[i] %in% rownames(pcb1) && sb2[i] %in% rownames(pcb2)) {
      data.frame(x = pcb1[sb1[i], 1], y = pcb1[sb1[i], 2],
                 xend = pcb2[sb2[i], 1], yend = pcb2[sb2[i], 2],
                 pair_id = i, stage = stage_label)
    }
  })
  do.call(rbind, Filter(Negate(is.null), rows))
}
bridge_conn <- rbind(
  build_connectors(pc_pre_b1, pc_pre_b2, sids_b1, sids_b2, "Pre-harmonisation"),
  build_connectors(pc_post_b1, pc_post_b2, sids_b1, sids_b2, "Post-harmonisation")
)
bridge_conn$stage <- factor(bridge_conn$stage, levels = c("Pre-harmonisation", "Post-harmonisation"))

pca_plot_df <- pca_df
n_non_bridge <- sum(!pca_df$bridge)
if (n_non_bridge > 2000) {
  set.seed(42)
  keep_idx <- c(which(pca_df$bridge), sample(which(!pca_df$bridge), 2000))
  pca_plot_df <- pca_df[keep_idx, ]
}

p_pca <- ggplot(pca_plot_df[!pca_plot_df$bridge, ], aes(x = PC1, y = PC2, colour = batch)) +
  geom_point(alpha = 0.2, size = 0.5) +
  geom_segment(data = bridge_conn, aes(x = x, y = y, xend = xend, yend = yend),
               inherit.aes = FALSE, colour = "grey50", alpha = 0.4, linewidth = 0.5) +
  geom_point(data = pca_plot_df[pca_plot_df$bridge, ], size = 2.5, alpha = 0.9, shape = 17) +
  facet_wrap(~stage) +
  scale_colour_manual(values = c("#2c7bb6", "#d7191c"), name = "Batch") +
  labs(title = "PCA (Fixed Basis): Before vs After Harmonisation",
       subtitle = paste0(
         "All samples projected onto a fixed PCA basis trained on combined pre-harmonisation data.\n",
         "Triangles mark bridge samples; lines connect the same individual across batches.\n",
         sprintf("Shorter connectors post-harmonisation indicate improved alignment. PCs retained: %d",
                 pca_basis$n_pcs))) +
  theme(legend.position = "bottom")

# Page 5: Left = pre, Right = post (same convention as plot B: triangles, connectors;
# flagged bridge samples larger + border; PCA space derived separately before/after).
# Refactored 2026-02; convention aligned with "PCA (Fixed Basis): Before vs After Harmonisation".
batch_cols_page5 <- setNames(c("#2c7bb6", "#d7191c"), c(batch1_id, batch2_id))
flagged_finngenid <- bridge_pair_dt[outlier_flag == TRUE]$finngenid

# Left: Pre-harmonisation — PCA from aggregated pre data; all samples, bridge = triangles, connectors
pca_pre_all <- pca_df[pca_df$stage == "Pre-harmonisation", ]
pca_pre_all$finngenid <- NA_character_
for (i in seq_along(common_bridge_fg)) {
  pca_pre_all$finngenid[pca_pre_all$sample == sids_b1[i]] <- common_bridge_fg[i]
  pca_pre_all$finngenid[pca_pre_all$sample == sids_b2[i]] <- common_bridge_fg[i]
}
pca_pre_all$is_flagged <- !is.na(pca_pre_all$finngenid) & (pca_pre_all$finngenid %in% flagged_finngenid)
n_non_bridge_pre <- sum(!pca_pre_all$bridge)
if (n_non_bridge_pre > 2000) {
  set.seed(42)
  keep_pre <- c(which(pca_pre_all$bridge), sample(which(!pca_pre_all$bridge), 2000))
  pca_pre_all <- pca_pre_all[keep_pre, ]
}
bridge_conn_pre <- bridge_conn[bridge_conn$stage == "Pre-harmonisation", ]
p_pca_left <- ggplot(pca_pre_all[!pca_pre_all$bridge, ], aes(x = PC1, y = PC2, colour = batch)) +
  geom_point(alpha = 0.2, size = 0.5) +
  geom_segment(data = bridge_conn_pre, aes(x = x, y = y, xend = xend, yend = yend),
               inherit.aes = FALSE, colour = "grey50", alpha = 0.4, linewidth = 0.5) +
  geom_point(data = pca_pre_all[pca_pre_all$bridge & !pca_pre_all$is_flagged, ],
             size = 2.5, alpha = 0.9, shape = 17)
if (any(pca_pre_all$is_flagged, na.rm = TRUE)) {
  p_pca_left <- p_pca_left +
    geom_point(data = pca_pre_all[pca_pre_all$bridge & pca_pre_all$is_flagged, ],
               aes(x = PC1, y = PC2, fill = batch), size = 4, shape = 24, colour = PAL$accent, stroke = 1.5,
               inherit.aes = FALSE) +
    scale_fill_manual(values = batch_cols_page5, guide = "none")
}
p_pca_left <- p_pca_left +
  scale_colour_manual(values = batch_cols_page5, name = "Batch") +
  labs(
    title = "Pre-harmonisation: aggregated batches",
    subtitle = paste0(
      "PCA space derived from aggregated samples before harmonisation (same convention as page 1).\n",
      "All samples: small points by batch. Bridge samples: triangles; lines connect same individual across batches.\n",
      "Flagged bridge pairs shown larger with border for tracking before/after.")) +
  theme(legend.position = "bottom")

# Right: Post-harmonisation — PCA from aggregated post data; all samples, bridge = triangles, connectors
post_space_all <- data.frame(
  PC1 = c(pc_post_space_b1[, 1], pc_post_space_b2[, 1]),
  PC2 = c(pc_post_space_b1[, 2], pc_post_space_b2[, 2]),
  batch = c(rep(batch1_id, nrow(pc_post_space_b1)), rep(batch2_id, nrow(pc_post_space_b2))),
  sample = c(rownames(pc_post_space_b1), rownames(pc_post_space_b2)),
  stringsAsFactors = FALSE
)
post_space_all$batch <- factor(post_space_all$batch, levels = c(batch1_id, batch2_id))
post_space_all$bridge <- post_space_all$sample %in% c(sids_b1, sids_b2)
post_space_all$finngenid <- NA_character_
for (i in seq_along(common_bridge_fg)) {
  post_space_all$finngenid[post_space_all$sample == sids_b1[i]] <- common_bridge_fg[i]
  post_space_all$finngenid[post_space_all$sample == sids_b2[i]] <- common_bridge_fg[i]
}
post_space_all$is_flagged <- !is.na(post_space_all$finngenid) & (post_space_all$finngenid %in% flagged_finngenid)
n_non_bridge_post <- sum(!post_space_all$bridge)
if (n_non_bridge_post > 2000) {
  set.seed(42)
  keep_post <- c(which(post_space_all$bridge), sample(which(!post_space_all$bridge), 2000))
  post_space_all <- post_space_all[keep_post, ]
}
bridge_conn_post <- build_connectors(pc_post_space_b1, pc_post_space_b2, sids_b1, sids_b2, "Post")
if (nrow(post_space_all) > 0) {
  p_pca_right <- ggplot(post_space_all[!post_space_all$bridge, ], aes(x = PC1, y = PC2, colour = batch)) +
    geom_point(alpha = 0.2, size = 0.5) +
    geom_segment(data = bridge_conn_post, aes(x = x, y = y, xend = xend, yend = yend),
                 inherit.aes = FALSE, colour = "grey50", alpha = 0.4, linewidth = 0.5) +
    geom_point(data = post_space_all[post_space_all$bridge & !post_space_all$is_flagged, ],
               size = 2.5, alpha = 0.9, shape = 17) +
    scale_colour_manual(values = batch_cols_page5, name = "Batch", drop = FALSE)
  if (any(post_space_all$is_flagged, na.rm = TRUE)) {
    p_pca_right <- p_pca_right +
      geom_point(data = post_space_all[post_space_all$bridge & post_space_all$is_flagged, ],
                 aes(x = PC1, y = PC2, fill = batch), size = 4, shape = 24,
                 colour = PAL$accent, stroke = 1.5, inherit.aes = FALSE) +
      scale_fill_manual(values = batch_cols_page5, guide = "none", drop = FALSE)
  }
  p_pca_right <- p_pca_right +
    labs(
      title = "Post-harmonisation: aggregated batches",
      subtitle = paste0(
        "PCA space derived from aggregated samples after harmonisation (same convention as page 1).\n",
        "All samples: small points by batch. Bridge samples: triangles; lines connect same individual across batches.\n",
        "Flagged bridge pairs shown larger with border for visual inspection.")) +
    theme(legend.position = "bottom")
} else {
  p_pca_right <- ggplot(data.frame(PC1 = 0, PC2 = 0, batch = factor(batch1_id, levels = c(batch1_id, batch2_id))),
                        aes(x = PC1, y = PC2, colour = batch)) +
    geom_point(alpha = 0) +
    scale_colour_manual(values = batch_cols_page5, name = "Batch", drop = FALSE) +
    labs(
      title = "Post-harmonisation: aggregated batches",
      subtitle = paste0(
        "No samples in post-derived PC space (projection dropped all).\n",
        "PCA space is derived separately from aggregated data after harmonisation.")) +
    theme(legend.position = "bottom")
}
p_pca_page5 <- ggarrange(p_pca_left, p_pca_right, ncol = 2, widths = c(1, 1),
                         labels = c("A1", "A2"), common.legend = TRUE, legend = "bottom")

# Variance explained by top 10 PCs (pre vs post, same colour scheme)
n_pcs_plot <- min(10, ncol(cohort_pre$scores), ncol(cohort_post_space$scores))
tot_var_pre  <- sum(apply(cohort_pre$scores, 2, function(x) var(x, na.rm = TRUE)))
tot_var_post <- sum(apply(cohort_post_space$scores, 2, function(x) var(x, na.rm = TRUE)))
pre_pct  <- if (tot_var_pre > 0)  100 * apply(cohort_pre$scores[, seq_len(n_pcs_plot), drop = FALSE], 2, function(x) var(x, na.rm = TRUE)) / tot_var_pre else rep(NA_real_, n_pcs_plot)
post_pct <- if (tot_var_post > 0) 100 * apply(cohort_post_space$scores[, seq_len(n_pcs_plot), drop = FALSE], 2, function(x) var(x, na.rm = TRUE)) / tot_var_post else rep(NA_real_, n_pcs_plot)
pc_labels_10 <- paste0("PC", seq_len(n_pcs_plot))
var_df <- data.frame(
  PC = factor(rep(pc_labels_10, 2), levels = pc_labels_10),
  stage = factor(rep(c("Pre", "Post"), each = n_pcs_plot), levels = c("Pre", "Post")),
  var_explained = c(pre_pct, post_pct),
  stringsAsFactors = FALSE
)
var_df$var_explained[!is.finite(var_df$var_explained)] <- 0
p_var_explained <- ggplot(var_df, aes(x = PC, y = var_explained, colour = stage, group = stage)) +
  geom_line(linewidth = 1, alpha = 0.9, na.rm = TRUE) +
  geom_point(size = 3, alpha = 0.9, na.rm = TRUE) +
  scale_colour_manual(values = c(Pre = PAL$pre, Post = PAL$post), name = "Stage") +
  labs(title = "Variance explained by top 10 PCs (pre vs post)",
       subtitle = "Pre = PC space from aggregated batches before harmonisation; Post = PC space from aggregated batches after harmonisation. Same colour scheme as PCA panels.",
       x = "Principal component", y = "Variance explained (%)") +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 45, hjust = 1))

# A2: NN Pairing Heatmap (swap detector) — with row annotation strip (grey = normal, accent = flagged)
nn_heatmap_long <- as.data.frame(as.table(nn_dist_mat_post))
names(nn_heatmap_long) <- c("B1_label", "B2_label", "distance")
nn_heatmap_long$B1_idx <- match(nn_heatmap_long$B1_label, fg_short)
nn_heatmap_long$B2_idx <- match(nn_heatmap_long$B2_label, fg_short)
nn_heatmap_long$is_diagonal <- nn_heatmap_long$B1_idx == nn_heatmap_long$B2_idx
nn_heatmap_long$is_outlier_cell <- nn_heatmap_long$B1_label %in% outlier_fg_short |
                                    nn_heatmap_long$B2_label %in% outlier_fg_short

nn_row_annot_df <- data.frame(
  B1_idx = seq_along(fg_short),
  B1_label = fg_short,
  is_flagged = fg_short %in% outlier_fg_short
)
nn_row_annot_df$is_flagged <- factor(nn_row_annot_df$is_flagged, levels = c(FALSE, TRUE),
                                     labels = c("Normal", "Flagged"))
p_nn_row_annot <- ggplot(nn_row_annot_df, aes(x = 1, y = B1_idx, fill = is_flagged)) +
  geom_tile(colour = "white", linewidth = 0.5) +
  scale_fill_manual(values = c(Normal = "grey90", Flagged = PAL$accent), name = "Row") +
  scale_y_continuous(breaks = seq_along(fg_short), expand = c(0, 0), limits = c(0.5, length(fg_short) + 0.5)) +
  theme_bw() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.title = element_blank(), panel.grid = element_blank(),
        legend.position = "left", plot.margin = margin(5, 0, 5, 5))

p_nn_heatmap_main <- ggplot(nn_heatmap_long, aes(x = B2_idx, y = B1_idx, fill = distance)) +
  geom_tile(colour = "white", linewidth = 0.5) +
  geom_tile(data = nn_heatmap_long[nn_heatmap_long$is_outlier_cell & !nn_heatmap_long$is_diagonal, ],
            colour = PAL$accent, linewidth = 0.5) +
  geom_text(aes(label = round(distance, 2)), size = 3, colour = "black", fontface = "bold") +
  scale_fill_distiller(palette = "RdBu", name = "Distance\n(MAD-Whitened)", direction = -1) +
  scale_x_continuous(breaks = seq_along(fg_short), labels = fg_short, expand = c(0, 0)) +
  scale_y_continuous(breaks = seq_along(fg_short), labels = fg_short, expand = c(0, 0),
                     limits = c(0.5, length(fg_short) + 0.5)) +
  labs(title = "Nearest-Neighbour Pairing Heatmap (Swap Detector)",
       subtitle = paste0(
         "MAD-whitened distance between every B1-B2 bridge pair after harmonisation.\n",
         "Diagonal = true pairs; off-diagonal minima suggest potential sample swaps. Left strip = row outlier flag.\n",
         sprintf("Flagged outliers: %d | n = %d bridge pairs",
                 length(outlier_fg_short), n_bridge)),
       x = "Batch 2 Bridge Sample", y = "Batch 1 Bridge Sample") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        axis.text.y = element_text(size = 8), legend.position = "right",
        plot.subtitle = element_text(size = 9, lineheight = 1.2), plot.margin = margin(5, 5, 5, 0))

p_nn_heatmap <- ggarrange(p_nn_row_annot, p_nn_heatmap_main, widths = c(1, 12), align = "h", common.legend = FALSE)

# A3/A4: Protein Correlation Heatmaps (Pre and Post)
build_corr_heatmap <- function(corr_mat, labels, stage_label, flagged_fg, pal_accent) {
  corr_long <- as.data.frame(as.table(corr_mat))
  names(corr_long) <- c("B1_label", "B2_label", "correlation")
  corr_long$B1_idx <- match(corr_long$B1_label, labels)
  corr_long$B2_idx <- match(corr_long$B2_label, labels)
  corr_long$is_diagonal <- corr_long$B1_idx == corr_long$B2_idx
  corr_long$is_flagged <- corr_long$B1_label %in% flagged_fg |
                           corr_long$B2_label %in% flagged_fg
  caption_text <- if (length(flagged_fg) > 0) {
    sprintf("Flagged bridges: %s", paste(head(flagged_fg, 10), collapse = ", "))
  } else "No flagged bridges"

  diag_median <- median(diag(corr_mat), na.rm = TRUE)

  # Row annotation: which B1 rows (sample rows) are flagged as outliers
  row_annot_df <- data.frame(
    B1_idx = seq_along(labels),
    B1_label = labels,
    is_flagged = labels %in% flagged_fg
  )
  row_annot_df$is_flagged <- factor(row_annot_df$is_flagged, levels = c(FALSE, TRUE),
                                   labels = c("Normal", "Flagged"))
  p_row_annot <- ggplot(row_annot_df, aes(x = 1, y = B1_idx, fill = is_flagged)) +
    geom_tile(colour = "white", linewidth = 0.5) +
    scale_fill_manual(values = c(Normal = "grey90", Flagged = pal_accent), name = "Row") +
    scale_y_continuous(breaks = seq_along(labels), expand = c(0, 0), limits = c(0.5, length(labels) + 0.5)) +
    theme_bw() +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          axis.text.y = element_blank(), axis.ticks.y = element_blank(),
          axis.title = element_blank(), panel.grid = element_blank(),
          legend.position = "left", plot.margin = margin(5, 0, 5, 5))

  # Main heatmap: RdBu, intensity by |r| (midpoint = 0); thick white borders; no yellow diagonal
  p_main <- ggplot(corr_long, aes(x = B2_idx, y = B1_idx, fill = correlation)) +
    geom_tile(colour = "white", linewidth = 0.5) +
    geom_text(aes(label = round(correlation, 2)), size = 3.2, colour = "black", fontface = "bold") +
    scale_fill_gradient2(low = "#B2182B", mid = "white", high = "#2166AC",
                         midpoint = 0, limits = c(-1, 1), name = "Pearson r",
                         na.value = "transparent") +
    scale_x_continuous(breaks = seq_along(labels), labels = labels, expand = c(0.02, 0)) +
    scale_y_continuous(breaks = seq_along(labels), labels = labels, expand = c(0, 0),
                       limits = c(0.5, length(labels) + 0.5)) +
    labs(title = sprintf("%s: Bridge Sample Protein Correlation (B1 vs B2)", stage_label),
         subtitle = paste0(
           "Pearson correlation of full protein profiles between each B1-B2 bridge pair.\n",
           "Color intensity = |r|; red = negative, blue = positive. Left strip = row outlier flag.\n",
           sprintf("Diagonal median r = %.3f | Flagged rows: %d",
                   diag_median, length(flagged_fg))),
         x = "Batch 2 Bridge Sample", y = "Batch 1 Bridge Sample",
         caption = caption_text) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
          axis.text.y = element_text(size = 8), legend.position = "right",
          plot.caption = element_text(size = 8, hjust = 0), panel.grid = element_blank(),
          plot.margin = margin(5, 5, 5, 0))

  ggarrange(p_row_annot, p_main, widths = c(1, 12), align = "h", common.legend = FALSE)
}

p_corr_pre  <- build_corr_heatmap(bridge_corr_pre,  fg_short, "Pre-Harmonisation",  outlier_fg_short, PAL$accent)
p_corr_post <- build_corr_heatmap(bridge_corr_post, fg_short, "Post-Harmonisation", outlier_fg_short, PAL$accent)

# A6: Pairing Margin Plot
margin_df <- data.frame(
  finngenid = common_bridge_fg, margin = pairing_margins,
  is_negative = pairing_margins < 0
)
p_margin <- ggplot(margin_df, aes(x = reorder(finngenid, margin), y = margin, fill = is_negative)) +
  geom_col(alpha = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "red", linewidth = 0.7) +
  scale_fill_manual(values = c("FALSE" = PAL$post, "TRUE" = PAL$accent),
                    labels = c("Positive margin", "Negative (swap candidate)")) +
  labs(title = "Pairing Margin (Confidence of Correct Match)",
       subtitle = paste0(
         "Margin = distance to true pair minus distance to best wrong pair (MAD-whitened).\n",
         "Negative margins (red) indicate the true pair is NOT the nearest neighbour — potential swap.\n",
         sprintf("Negative margins: %d/%d (swap candidates)",
                 sum(margin_df$is_negative), n_bridge)),
       x = "Bridge Pair (FINNGENID)", y = "Margin (MAD-whitened distance)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
        legend.position = "bottom", legend.title = element_blank())

# A7: Enhanced 3-Panel Identity Test
# Panel A: Per-pair rank bar
rank1_long_df <- rbind(
  data.frame(finngenid = common_bridge_fg, value = id_pre$ranks,
             is_rank1 = id_pre$ranks == 1, stage = "Pre"),
  data.frame(finngenid = common_bridge_fg, value = id_post$ranks,
             is_rank1 = id_post$ranks == 1, stage = "Post")
)
rank1_long_df$stage <- factor(rank1_long_df$stage, levels = c("Pre", "Post"))
p_rank1_detail <- ggplot(rank1_long_df, aes(x = reorder(finngenid, value), y = value, fill = stage)) +
  geom_col(alpha = 0.8, position = "dodge") +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "green", linewidth = 0.7) +
  scale_fill_manual(values = c(Pre = PAL$pre, Post = PAL$post)) +
  labs(title = "Rank-1 Identification (per Bridge Pair)",
       subtitle = paste0(
         "For each bridge sample in B1, rank of its true B2 pair among all B2 candidates (kNN).\n",
         "Rank = 1 (green dashed line) means the true pair is the nearest neighbour.\n",
         sprintf("Pre: %d/%d (%.1f%%) rank-1 | Post: %d/%d (%.1f%%) rank-1",
                 id_pre$rank1_hits, n_bridge, id_pre$rank1_rate * 100,
                 id_post$rank1_hits, n_bridge, id_post$rank1_rate * 100)),
       x = "Bridge Pair", y = "Rank of True Pair") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7), legend.position = "bottom")

# Panel B: Mahalanobis scatter
maha_scatter_df <- data.frame(x = d_maha_pre, y = d_maha_post,
                              outlier_flag = bridge_pair_dt$outlier_flag)
p_maha_scatter <- ggplot(maha_scatter_df, aes(x = x, y = y, colour = outlier_flag)) +
  geom_point(size = 2.5, alpha = 0.7) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
  geom_smooth(aes(x = x, y = y), method = "lm", se = TRUE, colour = "black",
              linewidth = 0.8, inherit.aes = FALSE) +
  scale_colour_manual(values = c("FALSE" = PAL$post, "TRUE" = PAL$accent),
                      labels = c("Normal", "Flagged")) +
  labs(title = "Mahalanobis Distance: Pre vs Post",
       subtitle = paste0(
         "Each point is one bridge pair; axes show robust Mahalanobis distance in PC space.\n",
         "Points below the diagonal improved after harmonisation; trend line shows overall shift.\n",
         sprintf("Flagged outliers (red): %d/%d", n_outliers, n_bridge)),
       x = "Pre Mahalanobis", y = "Post Mahalanobis") +
  theme(legend.position = "bottom")

# Panel C: MAD scatter
mad_scatter_df <- data.frame(x = d_mad_pre, y = d_mad_post,
                             outlier_flag = bridge_pair_dt$outlier_flag)
p_mad_scatter <- ggplot(mad_scatter_df, aes(x = x, y = y, colour = outlier_flag)) +
  geom_point(size = 2.5, alpha = 0.7) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
  geom_smooth(aes(x = x, y = y), method = "lm", se = TRUE, colour = "black",
              linewidth = 0.8, inherit.aes = FALSE) +
  scale_colour_manual(values = c("FALSE" = PAL$post, "TRUE" = PAL$accent),
                      labels = c("Normal", "Flagged")) +
  labs(title = "MAD-Whitened Distance: Pre vs Post",
       subtitle = paste0(
         "Each point is one bridge pair; axes show MAD-whitened Euclidean distance (primary KPI).\n",
         "Points below the diagonal improved after harmonisation; trend line shows overall shift.\n",
         sprintf("Flagged outliers (red): %d/%d", n_outliers, n_bridge)),
       x = "Pre MAD-Whitened", y = "Post MAD-Whitened") +
  theme(legend.position = "bottom")

# B1 chart: PC-by-Batch Effect Size Reduction
pc_labels <- paste0("PC", seq_along(pc_r2_pre))
pc_effect_df <- data.frame(
  PC = factor(rep(pc_labels, 2), levels = pc_labels),
  stage = factor(rep(c("Pre", "Post"), each = length(pc_r2_pre)), levels = c("Pre", "Post")),
  R2 = c(pc_r2_pre, pc_r2_post)
)
p_pc_effect <- ggplot(pc_effect_df, aes(x = PC, y = R2, colour = stage, group = stage)) +
  geom_line(linewidth = 1.2, alpha = 0.8) +
  geom_point(size = 3, alpha = 0.8) +
  scale_colour_manual(values = c(Pre = PAL$pre, Post = PAL$post)) +
  labs(title = "PC-by-Batch Effect Size (R-squared per PC)",
       subtitle = paste0(
         "R-squared from regressing each PC on batch label — quantifies batch confounding per component.\n",
         "Lower values indicate less batch effect. A large drop on PC1 is the strongest sign of correction.\n",
         sprintf("PC1 R-squared: pre = %.4f, post = %.4f", pc_r2_pre[1], pc_r2_post[1])),
       x = "Principal Component", y = "R-squared (Batch Effect)") +
  theme(legend.position = "bottom", legend.title = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))

# A8: Outlier Drill-Down Dashboard
p_outlier_drill <- NULL
if (n_outliers > 0) {
  flagged_pairs <- as.data.frame(bridge_pair_dt[outlier_flag == TRUE])
  outlier_panels <- list()
  for (idx in seq_len(min(6, nrow(flagged_pairs)))) {
    pair_idx <- which(bridge_pair_dt$finngenid == flagged_pairs$finngenid[idx])
    if (length(pair_idx) != 1) next
    sid_b1_pair <- sids_b1[pair_idx]
    sid_b2_pair <- sids_b2[pair_idx]
    if (!(sid_b1_pair %in% rownames(pc_pre_b1) && sid_b2_pair %in% rownames(pc_pre_b2) &&
          sid_b1_pair %in% rownames(pc_post_b1) && sid_b2_pair %in% rownames(pc_post_b2))) next
    pair_pc_df <- data.frame(
      stage = rep(c("Pre", "Post"), each = 2),
      batch = rep(c("B1", "B2"), 2),
      PC1 = c(pc_pre_b1[sid_b1_pair, 1], pc_pre_b2[sid_b2_pair, 1],
              pc_post_b1[sid_b1_pair, 1], pc_post_b2[sid_b2_pair, 1]),
      PC2 = c(pc_pre_b1[sid_b1_pair, 2], pc_pre_b2[sid_b2_pair, 2],
              pc_post_b1[sid_b1_pair, 2], pc_post_b2[sid_b2_pair, 2])
    )
    conn_df <- data.frame(
      x    = c(pc_pre_b1[sid_b1_pair, 1], pc_post_b1[sid_b1_pair, 1]),
      y    = c(pc_pre_b1[sid_b1_pair, 2], pc_post_b1[sid_b1_pair, 2]),
      xend = c(pc_pre_b2[sid_b2_pair, 1], pc_post_b2[sid_b2_pair, 1]),
      yend = c(pc_pre_b2[sid_b2_pair, 2], pc_post_b2[sid_b2_pair, 2]),
      stage = c("Pre", "Post")
    )
    d_pre_val  <- if (pair_idx <= length(d_mad_pre))  d_mad_pre[pair_idx]  else NA_real_
    d_post_val <- if (pair_idx <= length(d_mad_post)) d_mad_post[pair_idx] else NA_real_
    lr_val     <- if (pair_idx <= length(log_ratio_mad)) log_ratio_mad[pair_idx] else NA_real_
    margin_val <- if (pair_idx <= length(pairing_margins)) pairing_margins[pair_idx] else NA_real_

    p_out_i <- ggplot(pair_pc_df, aes(x = PC1, y = PC2, colour = batch, shape = batch)) +
      geom_segment(data = conn_df, aes(x = x, y = y, xend = xend, yend = yend),
                   inherit.aes = FALSE, colour = "grey50", alpha = 0.6, linewidth = 0.8) +
      geom_point(size = 4, alpha = 0.9) +
      facet_wrap(~stage) +
      scale_colour_manual(values = c(B1 = "#2c7bb6", B2 = "#d7191c")) +
      scale_shape_manual(values = c(B1 = 16, B2 = 17)) +
      labs(title = sprintf("Outlier: %s", substr(flagged_pairs$finngenid[idx], 1, 10)),
           subtitle = paste0(
             "PC1-PC2 positions of this bridge pair before and after harmonisation.\n",
             sprintf("D_pre = %.2f | D_post = %.2f | LR = %.3f | Margin = %.2f",
                     d_pre_val, d_post_val, lr_val, margin_val))) +
      theme(legend.position = "none", plot.title = element_text(size = 9),
            plot.subtitle = element_text(size = 7))
    outlier_panels[[length(outlier_panels) + 1]] <- p_out_i
  }
  if (length(outlier_panels) > 0) {
    p_outlier_drill <- ggarrange(plotlist = outlier_panels, ncol = 3, nrow = 2)
  }
}

# ------------------------------------------------------------------------------
# Write PDF (multi-page, 9 pages)
# ------------------------------------------------------------------------------
tryCatch({
  pdf(dashboard_path, width = 16, height = 20)

  # Page 1: Dumbbell + PCA Scatter with Connectors
  page1 <- ggarrange(p_dumbbell, p_pca, ncol = 1, nrow = 2,
                     labels = c("A", "B"), heights = c(1, 1.2))
  print(annotate_figure(page1, top = text_grob(
    "Cross-Batch Harmonisation KPI Dashboard — Summary", face = "bold", size = 15)))

  # Page 2: Boxplot + ECDF + Spaghetti + Log-ratio
  page2 <- ggarrange(p_box, p_ecdf, p_spaghetti, p_lr,
                     ncol = 2, nrow = 2, labels = c("A", "B", "C", "D"))
  print(annotate_figure(page2, top = text_grob(
    "Cross-Batch Harmonisation KPI Dashboard — Distance Metrics", face = "bold", size = 15)))

  # Page 3: NN Pairing Heatmap (full page)
  print(annotate_figure(p_nn_heatmap, top = text_grob(
    "Cross-Batch Harmonisation KPI Dashboard — Nearest-Neighbour Pairing Heatmap",
    face = "bold", size = 15)))

  # Page 4: Identity Tests + Pairing Margins
  page4 <- ggarrange(p_rank1_detail, p_margin, p_maha_scatter, p_mad_scatter,
                     ncol = 2, nrow = 2, labels = c("A", "B", "C", "D"))
  print(annotate_figure(page4, top = text_grob(
    "Cross-Batch Harmonisation KPI Dashboard — Pairing & Identity Tests",
    face = "bold", size = 15)))

  # Page 5: PCA pair (pre all / post bridge flagged) + Variance explained + Batch Separability
  page5 <- ggarrange(p_pca_page5, p_var_explained, p_sep, ncol = 1, nrow = 3,
                     labels = c("A", "B", "C"), heights = c(1.2, 0.7, 1))
  print(annotate_figure(page5, top = text_grob(
    "Cross-Batch Harmonisation KPI Dashboard — PCA & Batch Separability",
    face = "bold", size = 15)))

  # Page 6: PC Effect + Noise + Outlier barplot + Var Decomp
  page6_list <- list(p_pc_effect, p_noise, p_out)
  page6_labels <- c("A", "B", "C")
  if (!is.null(p_var)) {
    page6_list <- c(page6_list, list(p_var))
    page6_labels <- c(page6_labels, "D")
  }
  page6 <- ggarrange(plotlist = page6_list, ncol = 2, nrow = 2, labels = page6_labels)
  print(annotate_figure(page6, top = text_grob(
    "Cross-Batch Harmonisation KPI Dashboard — PC Effects & Additional Metrics",
    face = "bold", size = 15)))

  # Page 7: Protein Correlation Heatmap (Pre)
  print(annotate_figure(p_corr_pre, top = text_grob(
    "Cross-Batch Harmonisation KPI Dashboard — Bridge Protein Correlation (Pre)",
    face = "bold", size = 15)))

  # Page 8: Protein Correlation Heatmap (Post)
  print(annotate_figure(p_corr_post, top = text_grob(
    "Cross-Batch Harmonisation KPI Dashboard — Bridge Protein Correlation (Post)",
    face = "bold", size = 15)))

  # Page 9: Outlier Drill-Down (if any)
  if (!is.null(p_outlier_drill)) {
    tryCatch({
      print(annotate_figure(p_outlier_drill, top = text_grob(
        "Cross-Batch Harmonisation KPI Dashboard — Outlier Drill-Down",
        face = "bold", size = 15)))
    }, error = function(e) log_warn("Failed to generate outlier drill-down page: {e$message}"))
  }

  dev.off()
  log_info("Saved dashboard PDF: {dashboard_path}")
}, error = function(e) {
  tryCatch(dev.off(), error = function(e2) NULL)
  log_error("Failed to generate dashboard PDF: {e$message}")
  stop(e)
})

# ==============================================================================
# FINAL SUMMARY
# ==============================================================================
log_info(strrep("=", 80))
log_info("HARMONISATION KPI SUMMARY")
log_info(strrep("=", 80))
log_info("Bridge pairs analysed: {n_bridge}")
log_info("MAD-whitened distance reduction: {round(median(d_mad_pre),3)} -> {round(median(d_mad_post),3)}")
log_info("Paired log-ratio median: {round(median_log_ratio_mad,4)} | CI [{round(boot_ci_mad[1],4)},{round(boot_ci_mad[2],4)}]")
log_info("Noise floor MAD: {round(noise_floor_mad,3)} | ratio={round(noise_floor_ratio,3)}")
log_info("")
log_info("NOISE FLOOR CALCULATION DETAILS:")
log_info("  The noise floor represents measurement noise inherent in the assay.")
log_info("  Computed as median MAD-whitened Euclidean distance between within-batch")
log_info("  replicate pairs (same sample measured multiple times in the same batch).")
log_info("  Within-batch replicate pairs: B1={nrow(wb_rep_b1)}, B2={nrow(wb_rep_b2)}")
log_info("  Noise floor (MAD-whitened): {round(noise_floor_mad, 3)}")
log_info("  Bridge post-harmonisation distance (median): {round(median(d_mad_post), 3)}")
log_info("  Ratio (bridge_post / noise): {round(noise_floor_ratio, 3)}")
log_info("  Interpretation: Ratio < 2.0 = bridge pairs approaching measurement noise.")
log_info("")
if (!is.null(var_decomp_pre) && !is.null(var_decomp_post)) {
  hv_proteins <- unique(c(
    var_decomp_pre$protein[var_decomp_pre$batch_var_frac > 0.3 | var_decomp_pre$sex_var_frac > 0.3],
    var_decomp_post$protein[var_decomp_post$batch_var_frac > 0.3 | var_decomp_post$sex_var_frac > 0.3]
  ))
  if (length(hv_proteins) > 0) {
    log_info("HIGH-VARIANCE PROTEINS (>30% variance explained by Batch or Sex):")
    for (prot in head(hv_proteins, 20)) {
      pre_b <- var_decomp_pre$batch_var_frac[var_decomp_pre$protein == prot]
      pre_s <- var_decomp_pre$sex_var_frac[var_decomp_pre$protein == prot]
      post_b <- var_decomp_post$batch_var_frac[var_decomp_post$protein == prot]
      post_s <- var_decomp_post$sex_var_frac[var_decomp_post$protein == prot]
      if (length(pre_b) > 0 && length(post_b) > 0) {
        log_info("  {prot}: Pre (Batch={round(pre_b[1]*100,1)}%, Sex={round(pre_s[1]*100,1)}%) | Post (Batch={round(post_b[1]*100,1)}%, Sex={round(post_s[1]*100,1)}%)")
      }
    }
    if (length(hv_proteins) > 20) log_info("  ... and {length(hv_proteins) - 20} more")
  }
}
log_info("Rank-1 ID rate: {id_pre$rank1_hits}/{n_bridge} -> {id_post$rank1_hits}/{n_bridge}")
log_info("Separability: Silhouette {round(sil_pre,4)}->{round(sil_post,4)} | kBET {round(kbet_pre,4)}->{round(kbet_post,4)} | LISI {round(lisi_pre,4)}->{round(lisi_post,4)} | AUC {round(auc_pre,4)}->{round(auc_post,4)}")
log_info("PC batch R-squared (PC1): {round(pc_r2_pre[1],4)} -> {round(pc_r2_post[1],4)}")
log_info("Bridge outliers: {n_outliers}")
log_info("Pairing margins: {sum(pairing_margins < 0)} swap candidates")
log_info("Poorly harmonised proteins: {n_poor}")
log_info("")
log_info("Pass/Fail summary:")
for (i in seq_len(nrow(kpi_dt))) {
  if (!is.na(kpi_dt$pass[i])) {
    log_info("[{kpi_dt$pass[i]}] {kpi_dt$metric[i]}")
  }
}
log_info(strrep("=", 80))
log_info("Step 07b complete")

cat("\n=== CROSS-BATCH HARMONISATION KPI EVALUATION COMPLETE ===\n")
cat(sprintf("  KPI summary: %s\n", kpi_path))
cat(sprintf("  Dashboard:   %s\n", dashboard_path))
cat(sprintf("  Bridge pairs: %d | Outliers: %d | Poor proteins: %d\n", n_bridge, n_outliers, n_poor))