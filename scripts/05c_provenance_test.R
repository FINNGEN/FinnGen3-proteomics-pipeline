#!/usr/bin/env Rscript

#################################################
# Script: 05c_provenance_test.R
# Author: Reza Jabal, PhD (rjabal@broadinstitute.org)
# Description: Unit test for governance components validation
#              Refactored version
#              Refactored version
#              Tests sex outlier detection (04) and pQTL outlier detection (05b)
# Date: December 2025
#################################################

suppressPackageStartupMessages({
    library(data.table)
    library(yaml)
    library(logger)
    library(digest)
    library(ggplot2)
    library(ggrepel)
    library(glmnet)
    library(pROC)
    library(PRROC)
    library(MASS)
    library(ggpubr)
})

# Source path utilities
script_dir <- tryCatch({
    # Try to get script path from environment variable (set by pipeline runner)
    if (exists("SCRIPT_NAME", envir = .GlobalEnv) || "SCRIPT_NAME" %in% names(Sys.getenv())) {
        script_path <- Sys.getenv("SCRIPT_NAME", unset = get("SCRIPT_NAME", envir = .GlobalEnv, ifnotfound = NULL))
        if (!is.null(script_path) && file.exists(script_path)) {
            return(dirname(normalizePath(script_path)))
        }
    }
    # Try commandArgs (works when script is executed directly)
    cmd_args <- commandArgs(trailingOnly = FALSE)
    file_arg <- grep("^--file=", cmd_args, value = TRUE)
    if (length(file_arg) > 0) {
        script_path <- sub("^--file=", "", file_arg)
        return(dirname(normalizePath(script_path)))
    }
    # Fallback: try to find script location from sys.frame
    frame_files <- sapply(sys.frames(), function(x) {
        if (exists("ofile", x, inherits = FALSE)) {
            return(x$ofile)
        }
        return(NULL)
    })
    frame_files <- frame_files[!sapply(frame_files, is.null)]
    if (length(frame_files) > 0) {
        return(dirname(normalizePath(frame_files[[1]])))
    }
    # Last resort: assume we're in the scripts directory or use getwd
    wd <- getwd()
    if (file.exists(file.path(wd, "path_utils.R"))) {
        return(wd)
    }
    # Try scripts subdirectory
    scripts_dir <- file.path(wd, "scripts")
    if (file.exists(file.path(scripts_dir, "path_utils.R"))) {
        return(scripts_dir)
    }
    # Default fallback
    return(wd)
}, error = function(e) {
    # Final fallback
    wd <- getwd()
    if (file.exists(file.path(wd, "path_utils.R"))) {
        return(wd)
    }
    return(file.path(wd, "scripts"))
})

source(file.path(script_dir, "path_utils.R"))

# Get config path
config_file <- Sys.getenv(
    "PIPELINE_CONFIG",
    ""
)
config_file <- Sys.getenv("PIPELINE_CONFIG", "")
if (config_file == "" || !file.exists(config_file)) {
  stop("PIPELINE_CONFIG environment variable not set or config file not found. Please provide path to config file.")
}
config <- read_yaml(config_file)

# Check if test mode is enabled
test_enabled <- tryCatch(
    isTRUE(config$test_case$enabled),
    error = function(e) FALSE
)

if (!test_enabled) {
    log_info("Test case mode is disabled. Set test_case.enabled: true in config to enable.")
    quit(status = 0)
}

# Get batch context
batch_id <- Sys.getenv("PIPELINE_BATCH_ID", config$batch$default_batch_id %||% "batch_01")
step_num <- "07c"

# Set up logging - write to logs directory (primary)
# Override log filename to use 07_07c format (matching other 07* scripts)
base_dir <- config$output$base_dir %||% Sys.getenv("PIPELINE_OUTPUT_DIR", "output")
log_dir <- file.path(base_dir, config$output$logs_dir %||% "logs")
if (tryCatch(isTRUE(config$parameters$normalization$multi_batch_mode), error = function(e) FALSE) && !is.null(batch_id)) {
    log_dir <- file.path(log_dir, batch_id)
}
log_path <- file.path(log_dir, "07_07c_governance_test.log")
ensure_output_dir(log_path)
log_appender(appender_file(log_path))

# Also create a copy in test-case directory for convenience
test_case_log_path <- get_output_path(step_num, "governance_test", batch_id, "test-case", "log", config = config)
ensure_output_dir(test_case_log_path)

# Create a custom appender that writes to both files
dual_appender <- function(msg) {
    # Write to primary log (logs directory)
    cat(msg, file = log_path, append = TRUE)
    # Also write to test-case log
    cat(msg, file = test_case_log_path, append = TRUE)
}

# Use file appender for primary, and we'll copy the log at the end
log_info("Starting governance component validation test for batch: {batch_id}")

# Get test configuration
test_config <- tryCatch(config$test_case, error = function(e) NULL)
test_seed <- tryCatch(test_config$random_seed %||% 12345, error = function(e) 12345)
set.seed(test_seed)
log_info("Random seed set to {test_seed} for reproducibility")

# Get cohort size
cohort_size <- tryCatch(
    test_config$cohort$default_size %||% 500,
    error = function(e) 500
)
cohort_size <- max(100, min(500, cohort_size))  # Clamp to 100-500

# Get error sample sizes (NEW: 3 categories with exact counts)
n_category1 <- tryCatch(
    test_config$error_samples$category1_within_sex_id_shuffle %||% 20,
    error = function(e) 20
)
n_category2 <- tryCatch(
    test_config$error_samples$category2_sex_label_misannotation %||% 10,
    error = function(e) 10
)
n_category3 <- tryCatch(
    test_config$error_samples$category3_cross_sex_id_shuffle %||% 5,
    error = function(e) 5
)
total_errors <- n_category1 + n_category2 + n_category3

# Get genetic distance configuration
distance_config <- tryCatch(test_config$genetic_distance, error = function(e) NULL)
if (is.null(distance_config)) {
    distance_config <- list(
        enabled = FALSE,
        percentile_75_fraction = 0.6,
        percentile_90_fraction = 0.4,
        auto_lower_threshold = TRUE,
        log_warnings = TRUE
    )
}

log_info("Test configuration:")
log_info("  Cohort size: {cohort_size} samples")
log_info("  Category 1 (Within-Sex ID Shuffle): {n_category1} mismatches")
log_info("  Category 2 (Sex Label Misannotation): {n_category2} mismatches")
log_info("  Category 3 (Cross-Sex ID Shuffle): {n_category3} mismatches")
log_info("  Total errors: {total_errors} mismatches")
log_info("  Genetic distance stratification: {if(distance_config$enabled %||% FALSE) 'enabled' else 'disabled'}")

# Helper function to add disease group column
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
                                  ifelse(lengths(groups) == 1,
                                         vapply(groups, function(x) if(length(x) > 0) x[1] else "None", character(1)),
                                         "Multiple"))]
    dt[, groups := NULL]
    return(dt)
}

# ============================================================================
# Genetic Distance Calculation Functions
# ============================================================================

# Helper function to read gzipped files robustly (handles /tmp full scenarios)
# Tries zcat first (faster), falls back to gzfile() if zcat fails
read_gzipped_file <- function(file_path) {
    # Check if file is gzipped
    is_gzipped <- grepl("\\.gz$", file_path, ignore.case = TRUE)

    if (!is_gzipped) {
        # Not gzipped, read directly
        return(fread(file_path))
    }

    # Try zcat first (faster, but requires /tmp space)
    result <- tryCatch({
        fread(cmd = paste("zcat", file_path))
    }, error = function(e) {
        # zcat failed (likely /tmp full), try gzfile() instead
        log_warn("zcat failed (likely /tmp full): {e$message}. Falling back to gzfile()...")
        return(NULL)
    })

    if (!is.null(result)) {
        return(result)
    }

    # Fallback: use gzfile() which decompresses in memory (no /tmp needed)
    log_info("Using gzfile() for decompression (no /tmp space required)")
    return(fread(gzfile(file_path)))
}

# Load genetic PCs from covariate file (same approach as 04_pca_outliers.R)
load_genetic_pcs <- function(covariate_file, finngenids) {
    # Load genetic PCs from covariate file
    # Returns: data.table with FINNGENID, gPC1, gPC2, gPC3, gPC4

    log_info("Loading genetic PCs from covariate file: {covariate_file}")

    covariates <- tryCatch({
        read_gzipped_file(covariate_file)
    }, error = function(e) {
        log_error("Failed to load covariate file: {e$message}")
        stop("Cannot load genetic PCs from covariate file")
    })

    # Extract genetic PCs (PC1-PC4 from covariate file)
    if (!all(c("IID", "PC1", "PC2", "PC3", "PC4") %in% names(covariates))) {
        stop("Covariate file missing required columns: IID, PC1, PC2, PC3, PC4")
    }

    gen_pcs <- covariates[, .(
        FINNGENID = IID,
        gPC1 = PC1,
        gPC2 = PC2,
        gPC3 = PC3,
        gPC4 = PC4
    )]

    # Filter to requested FINNGENIDs
    gen_pcs <- gen_pcs[FINNGENID %in% finngenids]

    log_info("Loaded genetic PCs for {nrow(gen_pcs)} samples (out of {length(finngenids)} requested)")

    # Check for missing values
    n_missing <- sum(is.na(gen_pcs[, .(gPC1, gPC2, gPC3, gPC4)]))
    if (n_missing > 0) {
        log_warn("Found {n_missing} missing genetic PC values. Will impute with median.")
    }

    return(gen_pcs)
}

# Calculate pairwise genetic distances using Euclidean distance in PC space
calculate_genetic_distances <- function(gen_pcs, distance_metric = "euclidean") {
    # Calculate pairwise Euclidean distances in genetic PC space
    # Returns: distance matrix (samples x samples) with FINNGENID as row/col names

    log_info("Calculating genetic distances using {distance_metric} metric...")

    # Extract PC matrix
    pc_matrix <- as.matrix(gen_pcs[, .(gPC1, gPC2, gPC3, gPC4)])
    rownames(pc_matrix) <- gen_pcs$FINNGENID

    # Handle missing values (impute with column median)
    for (i in seq_len(ncol(pc_matrix))) {
        missing_idx <- is.na(pc_matrix[, i])
        if (any(missing_idx)) {
            median_val <- median(pc_matrix[, i], na.rm = TRUE)
            pc_matrix[missing_idx, i] <- median_val
            log_warn("Imputed {sum(missing_idx)} missing values in gPC{i} with median: {round(median_val, 4)}")
        }
    }

    # Calculate pairwise distances
    if (distance_metric == "euclidean") {
        dist_matrix <- as.matrix(dist(pc_matrix, method = "euclidean"))
    } else {
        stop("Only euclidean distance metric is currently supported. Got: {distance_metric}")
    }

    # Remove diagonal (self-distances)
    diag(dist_matrix) <- NA

    log_info("Calculated genetic distances for {nrow(dist_matrix)} samples")
    log_info("Distance range: {round(min(dist_matrix, na.rm=TRUE), 3)} - {round(max(dist_matrix, na.rm=TRUE), 3)}")
    log_info("Distance median: {round(median(dist_matrix, na.rm=TRUE), 3)}")

    return(dist_matrix)
}

# Filter swap pool by genetic distance percentile threshold
filter_swap_pool_by_distance <- function(
    source_finngenid,        # FINNGENID of sample to swap
    candidate_pool,         # Vector of candidate FINNGENIDs
    distance_matrix,         # Distance matrix
    min_percentile,         # Minimum distance percentile (75 or 90)
    auto_lower = TRUE,      # Auto-lower threshold if insufficient candidates
    log_warnings = TRUE      # Log warnings
) {
    # Filter candidate pool to only include swaps above min_percentile distance
    # Returns: filtered candidate pool

    if (length(candidate_pool) == 0) {
        return(character(0))
    }

    if (!source_finngenid %in% rownames(distance_matrix)) {
        if (log_warnings) {
            log_warn("Source FINNGENID {source_finngenid} not found in distance matrix")
        }
        return(character(0))
    }

    # Get distances from source to all candidates
    available_candidates <- intersect(candidate_pool, colnames(distance_matrix))
    if (length(available_candidates) == 0) {
        return(character(0))
    }

    distances <- distance_matrix[source_finngenid, available_candidates]
    distances <- distances[!is.na(distances) & is.finite(distances)]

    if (length(distances) == 0) {
        return(character(0))
    }

    # Calculate percentile threshold
    threshold <- quantile(distances, probs = min_percentile / 100, na.rm = TRUE)

    # Filter candidates above threshold
    filtered_pool <- names(distances)[distances >= threshold]

    # Auto-lower threshold if insufficient candidates
    if (length(filtered_pool) < 3 && auto_lower) {
        if (log_warnings) {
            log_warn("Insufficient candidates ({length(filtered_pool)}) at {min_percentile}th percentile for {source_finngenid}. Lowering threshold...")
        }

        # Lower threshold by 10 percentile points, but not below 50
        new_percentile <- max(50, min_percentile - 10)
        new_threshold <- quantile(distances, probs = new_percentile / 100, na.rm = TRUE)
        filtered_pool <- names(distances)[distances >= new_threshold]

        if (log_warnings) {
            log_warn("Lowered threshold to {new_percentile}th percentile. Found {length(filtered_pool)} candidates.")
        }
    }

    return(filtered_pool)
}

# Store genetic distances for plotting (distance from population center)
store_genetic_distances <- function(metadata, gen_pcs, error_tracking) {
    # Store genetic distance for each sample (distance from population center)
    # This is used for the scatter plot

    # Calculate population center (median of all PCs)
    pc_center <- gen_pcs[, .(
        gPC1_center = median(gPC1, na.rm = TRUE),
        gPC2_center = median(gPC2, na.rm = TRUE),
        gPC3_center = median(gPC3, na.rm = TRUE),
        gPC4_center = median(gPC4, na.rm = TRUE)
    )]

    # Calculate distance from center for each sample
    gen_pcs[, genetic_distance := sqrt(
        (gPC1 - pc_center$gPC1_center)^2 +
        (gPC2 - pc_center$gPC2_center)^2 +
        (gPC3 - pc_center$gPC3_center)^2 +
        (gPC4 - pc_center$gPC4_center)^2
    )]

    # Merge into error_tracking
    error_tracking <- merge(
        error_tracking,
        gen_pcs[, .(FINNGENID, genetic_distance)],
        by = "FINNGENID",
        all.x = TRUE
    )

    log_info("Stored genetic distances for {sum(!is.na(error_tracking$genetic_distance))} samples")

    return(error_tracking)
}

# Helper function for Inverse Rank Normalization (IRN) per protein
# Transforms each protein column to follow N(0,1) distribution
# This makes z-score calculation more robust to outliers and non-normal distributions
inverse_rank_normalize <- function(npx_matrix) {
    log_info("Applying Inverse Rank Normalization (IRN) per protein...")

    # Store original dimensions and names
    original_rownames <- rownames(npx_matrix)
    original_colnames <- colnames(npx_matrix)

    # Convert to matrix if not already
    if (!is.matrix(npx_matrix)) {
        npx_matrix <- as.matrix(npx_matrix)
    }

    # Apply IRN column-wise (per protein)
    npx_irn <- apply(npx_matrix, 2, function(x) {
        # Get non-NA values
        valid_idx <- !is.na(x)
        n_valid <- sum(valid_idx)

        # If too few valid values, return original (with NAs)
        if (n_valid < 3) {
            return(x)
        }

        # Extract valid values
        x_valid <- x[valid_idx]

        # Calculate ranks (average for ties)
        # rank() handles ties by default with average method
        ranks <- rank(x_valid, ties.method = "average", na.last = "keep")

        # Convert ranks to quantiles: (rank - 0.5) / n
        # This ensures quantiles are in (0, 1) range
        quantiles <- (ranks - 0.5) / n_valid

        # Apply inverse normal CDF (qnorm) to get N(0,1) distributed values
        # Handle edge cases: quantiles exactly 0 or 1 would give -Inf/Inf
        # Use small epsilon to avoid infinite values
        epsilon <- 1e-6
        quantiles <- pmax(epsilon, pmin(1 - epsilon, quantiles))
        x_irn_valid <- qnorm(quantiles)

        # Create output vector with same length as input
        x_irn <- rep(NA_real_, length(x))
        x_irn[valid_idx] <- x_irn_valid

        return(x_irn)
    })

    # Ensure matrix format and restore names
    npx_irn <- as.matrix(npx_irn)
    rownames(npx_irn) <- original_rownames
    colnames(npx_irn) <- original_colnames

    # Log summary statistics
    npx_irn_flat <- as.vector(npx_irn)
    npx_irn_flat <- npx_irn_flat[!is.na(npx_irn_flat)]
    if (length(npx_irn_flat) > 0) {
        log_info("IRN transformation complete:")
        log_info("  Mean: {round(mean(npx_irn_flat), 4)} (expected ~0)")
        log_info("  SD: {round(sd(npx_irn_flat), 4)} (expected ~1)")
        log_info("  Range: [{round(min(npx_irn_flat), 4)}, {round(max(npx_irn_flat), 4)}]")
    }

    return(npx_irn)
}

# Helper functions for sex detection (extracted from 05_sex_outliers.R)
# --------------------------------------------------------------------
residualize_by_confounders <- function(Y, X) {
    XtX <- crossprod(X)
    XtX_inv <- tryCatch(solve(XtX), error = function(e) MASS::ginv(XtX))
    B <- XtX_inv %*% crossprod(X, Y)
    R <- Y - X %*% B
    list(residuals = R, coef = B)
}

apply_residualization <- function(Y_new, X_new, coef_B) {
    Y_new - X_new %*% coef_B
}

compute_impute_means <- function(mat) {
    means <- colMeans(mat, na.rm = TRUE)
    means[is.na(means)] <- 0
    means
}

impute_matrix_with_means <- function(mat, means) {
    if (length(means) == 0) return(mat)
    out <- mat
    idx <- which(is.na(out), arr.ind = TRUE)
    if (nrow(idx) > 0) {
        out[idx] <- means[idx[, 2]]
    }
    out
}

make_folds <- function(y, k = 5, seed = 1) {
    set.seed(seed)
    y <- as.factor(y)
    idx <- split(seq_along(y), y)
    folds <- vector("list", k)
    for (lev in names(idx)) {
        ids <- sample(idx[[lev]])
        parts <- split(ids, rep(1:k, length.out = length(ids)))
        for (i in 1:k) folds[[i]] <- c(folds[[i]], parts[[i]])
    }
    lapply(folds, sort)
}

pick_threshold <- function(probs, labels) {
    if (length(probs) == 0 || all(is.na(probs))) {
        return(list(th = 0.5, J = 0, F1 = 0, sens = 0, spec = 0, acc = 0))
    }
    roc_obj <- try(pROC::roc(response = labels, predictor = probs, quiet = TRUE), silent = TRUE)
    if (inherits(roc_obj, "try-error")) {
        return(list(th = 0.5, J = 0, F1 = 0, sens = 0, spec = 0, acc = 0))
    }
    coords <- pROC::coords(roc_obj, "best", ret = c("threshold", "sensitivity", "specificity"), transpose = FALSE)
    th <- coords$threshold[1]
    sens <- coords$sensitivity[1]
    spec <- coords$specificity[1]
    J <- sens + spec - 1
    pred_binary <- ifelse(probs >= th, 1, 0)
    tp <- sum(pred_binary == 1 & labels == 1, na.rm = TRUE)
    fp <- sum(pred_binary == 1 & labels == 0, na.rm = TRUE)
    fn <- sum(pred_binary == 0 & labels == 1, na.rm = TRUE)
    prec <- if (tp + fp > 0) tp / (tp + fp) else 0
    rec <- if (tp + fn > 0) tp / (tp + fn) else 0
    F1 <- if (prec + rec > 0) 2 * prec * rec / (prec + rec) else 0
    acc <- mean(pred_binary == labels, na.rm = TRUE)
    list(th = th, J = J, F1 = F1, sens = sens, spec = spec, acc = acc)
}

# Two-panel distribution plot function (from 05_sex_outliers.R)
create_prediction_distribution_panels <- function(df_pred, title, file_out,
                                                  youden_thresh,
                                                  annotate_ids = FALSE, max_labels = 80,
                                                  subtitle = NULL) {
    # df_pred: columns SAMPLE_ID, FINNGENID, genetic_sex ("male"/"female"), predicted_prob (0..1)
    dt <- as.data.table(df_pred)
    dt <- dt[!is.na(predicted_prob) & !is.na(genetic_sex)]
    dt[, genetic_sex := factor(genetic_sex, levels = c("male", "female"))]
    # Guard against empty data (prevents faceting error)
    if (nrow(dt) == 0) {
        pg <- ggplot() + theme_void() + ggtitle("No data available for prediction distribution")
        dir.create(dirname(file_out), recursive = TRUE, showWarnings = FALSE)
        try(ggsave(file_out, pg, width = 12, height = 6), silent = TRUE)
        return(pg)
    }
    # numeric y for jitter
    dt[, y_num := as.numeric(genetic_sex)]
    # deterministic vertical jitter (consistent overlay) - larger spread
    set.seed(123)
    dt[, y_jit := y_num + runif(.N, min = -0.22, max = 0.22)]

    # Outlier flags for plotting on the probability strip
    dt[, is_mismatch := FALSE]
    dt[, is_sex_outlier := FALSE]
    # Mismatch: crosses Youden's J threshold
    dt[genetic_sex == "male" & predicted_prob >= youden_thresh, is_mismatch := TRUE]
    dt[genetic_sex == "female" & predicted_prob < youden_thresh, is_mismatch := TRUE]
    # Sex outlier: crosses 0.5 but not Youden's J
    dt[genetic_sex == "male" & predicted_prob >= 0.5 & predicted_prob < youden_thresh, is_sex_outlier := TRUE]

    # Classification: normal, outlier (pale red, between 0.5 and Youden), mismatch (red, beyond Youden)
    dt[, outlier_class := "normal"]
    dt[is_sex_outlier == TRUE, outlier_class := "outlier"]
    dt[is_mismatch == TRUE, outlier_class := "mismatch"]

    # Left: density + histogram faceted by sex
    p_left <- ggplot(dt, aes(x = predicted_prob, fill = genetic_sex)) +
        geom_histogram(aes(y = after_stat(density)), bins = 30, alpha = 0.6, position = "identity") +
        geom_density(alpha = 0.25, color = "black", linewidth = 0.3) +
        geom_vline(xintercept = 0.5, linetype = "dashed", color = "gray50", linewidth = 0.4) +
        geom_vline(xintercept = youden_thresh, linetype = "dotted", color = "#D7301F", linewidth = 0.8) +
        scale_fill_manual(values = c("male" = "#3A5F8A", "female" = "#FF6B6B")) +
        facet_wrap(~ genetic_sex, ncol = 1, scales = "free_y") +
        labs(title = title,
             subtitle = if (is.null(subtitle)) sprintf("Gray dashed = Outlier threshold (0.5); Red dotted = Mismatch threshold (Youden J = %.3f)", youden_thresh) else subtitle,
             x = "Predicted Female Probability", y = "Density", fill = "Genetic Sex") +
        theme_bw()

    # Right: probability strip by sex with jitter
    p_right <- ggplot(dt, aes(x = predicted_prob, y = y_jit)) +
        geom_vline(xintercept = 0.5, linetype = "dashed", color = "gray50", linewidth = 0.4) +
        geom_vline(xintercept = youden_thresh, linetype = "dotted", color = "#D7301F", linewidth = 0.8) +
        geom_point(data = dt[outlier_class == "normal"],
                   size = 0.9, alpha = 0.22, color = "gray20", show.legend = FALSE) +
        geom_point(data = dt[outlier_class == "outlier"],
                   shape = 21, size = 2.0, stroke = 0.3,
                   fill = "#D7301F", alpha = 0.35, color = "black", show.legend = FALSE) +
        geom_point(data = dt[outlier_class == "mismatch"],
                   shape = 21, size = 2.4, stroke = 0.5,
                   fill = "#D7301F", color = "black", show.legend = FALSE) +
        scale_y_continuous(breaks = c(1, 2), labels = c("male", "female"), limits = c(0.6, 2.4), expand = expansion(mult = c(0.02, 0.02))) +
        coord_cartesian(xlim = c(0, 1)) +
        theme_bw() +
        labs(x = "Predicted Female Probability", y = NULL) +
        theme(legend.position = "none") +
        annotate("rect", xmin = 0.02, xmax = 0.30, ymin = 2.25, ymax = 2.38, fill = "white", color = "black", linewidth = 0.3) +
        annotate("text", x = 0.04, y = 2.35, label = "Sex QC Flags", hjust = 0, size = 2.8, fontface = "bold") +
        annotate("point", x = 0.05, y = 2.31, shape = 21, size = 2.4, stroke = 0.5, fill = "#D7301F", color = "black") +
        annotate("text", x = 0.07, y = 2.31, label = "Mismatches", hjust = 0, size = 2.4) +
        annotate("point", x = 0.05, y = 2.27, shape = 21, size = 2.0, stroke = 0.3, fill = "#D7301F", alpha = 0.35, color = "black") +
        annotate("text", x = 0.07, y = 2.27, label = "Outliers", hjust = 0, size = 2.4)

    if (annotate_ids) {
        lab_dt <- dt[is_mismatch == TRUE & !is.na(FINNGENID)]
        if (nrow(lab_dt) > max_labels) {
            lab_m <- lab_dt[genetic_sex == "male"][order(-predicted_prob)][1:ceiling(max_labels/2)]
            lab_f <- lab_dt[genetic_sex == "female"][order(predicted_prob)][1:floor(max_labels/2)]
            lab_dt <- rbind(lab_m, lab_f, fill = TRUE)
        }
        if (nrow(lab_dt) > 0) {
            p_right <- p_right +
                ggrepel::geom_text_repel(data = lab_dt,
                                         aes(y = y_jit, label = FINNGENID),
                                         size = 1.8, color = "gray20",
                                         max.overlaps = 60, segment.size = 0.2,
                                         seed = 7)
        }
    }

    # Arrange as two panels
    pg <- ggpubr::ggarrange(p_left, p_right, ncol = 2, widths = c(2, 1))
    dir.create(dirname(file_out), recursive = TRUE, showWarnings = FALSE)
    try(ggsave(file_out, pg, width = 12, height = 6), silent = TRUE)
    pg
}

# ============================================================================
# Plotting Functions
# ============================================================================

# Calculate p-values for distribution comparisons
calculate_distribution_pvalues <- function(plot_data) {
    # Calculate Mann-Whitney U test p-values for each error category vs Clean
    clean_data <- plot_data[error_category == "Clean"]$MeanAbsZ
    clean_data <- clean_data[!is.na(clean_data) & is.finite(clean_data)]

    p_values <- data.table()

    for (cat in c("Within-Sex ID Shuffle", "Sex Label Misannotation", "Cross-Sex ID Shuffle")) {
        cat_data <- plot_data[error_category == cat]$MeanAbsZ
        cat_data <- cat_data[!is.na(cat_data) & is.finite(cat_data)]

        if (length(cat_data) > 0 && length(clean_data) > 0) {
            test_result <- tryCatch({
                wilcox.test(cat_data, clean_data)
            }, error = function(e) {
                list(p.value = NA_real_)
            })

            p_values <- rbind(p_values, data.table(
                category = cat,
                p_value = test_result$p.value
            ))
        }
    }

    return(p_values)
}

# Create faceted distribution plot (density curves + violin plots with p-values)
create_meanabsz_distribution_plot <- function(
    sample_stats,          # From pQTL detection (has MeanAbsZ column)
    ground_truth,          # Error tracking with categories
    output_path            # File path for saving plot
) {
    log_info("Creating MeanAbsZ distribution plot...")

    # Merge data
    plot_data <- merge(
        sample_stats[, .(SampleID, MeanAbsZ)],
        ground_truth[, .(SAMPLE_ID, category, swap_type)],
        by.x = "SampleID", by.y = "SAMPLE_ID",
        all.x = TRUE
    )

    # Create category labels
    plot_data[, error_category := fifelse(
        category == "Category1_WithinSex", "Within-Sex ID Shuffle",
        fifelse(category == "Category2_SexLabel", "Sex Label Misannotation",
        fifelse(category == "Category3_CrossSex", "Cross-Sex ID Shuffle",
        "Clean"))
    )]

    plot_data <- plot_data[!is.na(MeanAbsZ) & is.finite(MeanAbsZ)]

    # Calculate p-values
    p_values <- calculate_distribution_pvalues(plot_data)

    # Create faceted plot: violin plots with density overlays
    # For each category, create a combined violin + density plot
    p <- ggplot(plot_data, aes(x = error_category, y = MeanAbsZ)) +
        # Violin plots
        geom_violin(aes(fill = error_category), alpha = 0.4, trim = FALSE, width = 0.7) +
        # Box plots (for quartiles)
        geom_boxplot(width = 0.15, alpha = 0.7, outlier.shape = NA) +
        # Add p-values as text
        geom_text(data = p_values,
                 aes(x = category, y = max(plot_data$MeanAbsZ, na.rm = TRUE) * 1.15,
                     label = ifelse(is.na(p_value), "N/A",
                                   ifelse(p_value < 0.001, "p < 0.001",
                                          paste0("p = ", format(p_value, scientific = FALSE, digits = 3))))),
                 inherit.aes = FALSE, size = 3.5, fontface = "bold") +
        facet_wrap(~ error_category, scales = "free_y", ncol = 2) +
        scale_fill_manual(values = c(
            "Within-Sex ID Shuffle" = "#2E7D32",      # Green
            "Sex Label Misannotation" = "#FF9800",    # Orange
            "Cross-Sex ID Shuffle" = "#9C27B0",        # Purple
            "Clean" = "#90CAF9"                        # Light blue
        )) +
        labs(
            title = "Mean Absolute Z-Score Distribution by Error Category",
            subtitle = "Violin plots with statistical tests (Mann-Whitney U vs Clean)",
            x = "Error Category",
            y = "Mean Absolute Z-Score",
            fill = "Category"
        ) +
        theme_bw() +
        theme(
            plot.title = element_text(size = 14, face = "bold"),
            plot.subtitle = element_text(size = 11),
            legend.position = "none",
            strip.text = element_text(size = 11, face = "bold"),
            axis.text.x = element_text(angle = 45, hjust = 1)
        )

    # Create separate density plot for overlay
    p_density <- ggplot(plot_data, aes(x = MeanAbsZ, fill = error_category)) +
        geom_density(alpha = 0.5, adjust = 1.5) +
        facet_wrap(~ error_category, ncol = 2, scales = "free_y") +
        scale_fill_manual(values = c(
            "Within-Sex ID Shuffle" = "#2E7D32",
            "Sex Label Misannotation" = "#FF9800",
            "Cross-Sex ID Shuffle" = "#9C27B0",
            "Clean" = "#90CAF9"
        )) +
        labs(
            title = "Mean Absolute Z-Score Density Curves by Error Category",
            x = "Mean Absolute Z-Score",
            y = "Density",
            fill = "Category"
        ) +
        theme_bw() +
        theme(
            plot.title = element_text(size = 14, face = "bold"),
            legend.position = "right"
        )

    # Create overlay density plot: Correct samples vs All errors and Categories separately
    plot_data_overlay <- copy(plot_data)
    plot_data_overlay[, overlay_group := fifelse(
        error_category == "Clean", "Correct Samples",
        error_category
    )]

    # Create overlay plot
    p_overlay <- ggplot(plot_data_overlay, aes(x = MeanAbsZ, fill = overlay_group, color = overlay_group)) +
        # Density curves
        geom_density(alpha = 0.4, adjust = 1.5, linewidth = 1.2) +
        scale_fill_manual(values = c(
            "Correct Samples" = "#90CAF9",              # Light blue
            "Within-Sex ID Shuffle" = "#2E7D32",        # Green
            "Sex Label Misannotation" = "#FF9800",      # Orange
            "Cross-Sex ID Shuffle" = "#9C27B0"          # Purple
        ), name = "Sample Group") +
        scale_color_manual(values = c(
            "Correct Samples" = "#1976D2",              # Darker blue
            "Within-Sex ID Shuffle" = "#1B5E20",        # Darker green
            "Sex Label Misannotation" = "#E65100",     # Darker orange
            "Cross-Sex ID Shuffle" = "#4A148C"         # Darker purple
        ), guide = "none") +
        labs(
            title = "Mean Absolute Z-Score Distribution: Correct Samples vs Error Categories",
            subtitle = "Overlay density curves comparing correct samples to each error category",
            x = "Mean Absolute Z-Score",
            y = "Density",
            fill = "Sample Group"
        ) +
        theme_bw() +
        theme(
            plot.title = element_text(size = 14, face = "bold"),
            plot.subtitle = element_text(size = 11),
            legend.position = "right"
        )

    # Save violin plots on separate page (full page)
    violin_output_path <- gsub("\\.pdf$", "_violin_only.pdf", output_path)
    ggsave(violin_output_path, p, width = 14, height = 8)
    log_info("Saved violin plots (separate page) to: {violin_output_path}")

    # Combine density plots only (violin plots are now separate)
    p_combined <- ggpubr::ggarrange(p_density, p_overlay, nrow = 2, heights = c(1, 1), common.legend = FALSE)

    ggsave(output_path, p_combined, width = 14, height = 10)
    log_info("Saved distribution plot (density plots only) to: {output_path}")
    return(p_combined)
}

# Create scatter plot: gPC1 vs MeanAbsZ
create_gpc1_vs_meanabsz_plot <- function(
    sample_stats,          # From pQTL detection
    ground_truth,          # Error tracking with categories
    gen_pcs,               # Genetic PCs data.table with FINNGENID, gPC1
    output_path            # File path for saving plot
) {
    log_info("Creating gPC1 vs MeanAbsZ scatter plot...")

    # Merge data
    plot_data <- merge(
        sample_stats[, .(SampleID, MeanAbsZ, Outlier)],
        ground_truth[, .(SAMPLE_ID, FINNGENID, category, swap_type)],
        by.x = "SampleID", by.y = "SAMPLE_ID",
        all.x = TRUE
    )

    # Merge with genetic PCs
    if (!is.null(gen_pcs) && nrow(gen_pcs) > 0) {
        plot_data <- merge(
            plot_data,
            gen_pcs[, .(FINNGENID, gPC1)],
            by = "FINNGENID",
            all.x = TRUE
        )
    } else {
        log_warn("Genetic PCs not available. Skipping gPC1 scatter plot.")
        return(NULL)
    }

    # Create category labels
    plot_data[, error_category := fifelse(
        category == "Category1_WithinSex", "Within-Sex ID Shuffle",
        fifelse(category == "Category2_SexLabel", "Sex Label Misannotation",
        fifelse(category == "Category3_CrossSex", "Cross-Sex ID Shuffle",
        "Clean"))
    )]

    plot_data[, is_cross_sex := category == "Category3_CrossSex"]
    plot_data[, is_error := category != "Clean"]

    plot_data <- plot_data[!is.na(MeanAbsZ) & !is.na(gPC1) &
                          is.finite(MeanAbsZ) & is.finite(gPC1)]

    # Create scatter plot
    p <- ggplot(plot_data, aes(x = gPC1, y = MeanAbsZ)) +
        # Clean samples (background)
        geom_point(data = plot_data[error_category == "Clean"],
                  color = "gray70", alpha = 0.5, size = 2, shape = 21) +
        # Error samples by category
        geom_point(data = plot_data[error_category != "Clean" & !is_cross_sex],
                  aes(color = error_category),
                  size = 3, alpha = 0.8, shape = 16) +
        # Cross-sex mismatches with special shape and color
        geom_point(data = plot_data[is_cross_sex == TRUE],
                  color = "purple", size = 4, alpha = 0.9, shape = 23, stroke = 2) +
        # Highlight outliers
        geom_point(data = plot_data[Outlier == TRUE],
                  fill = NA, color = "red", size = 3.5, shape = 24, stroke = 1.5, alpha = 0.8) +
        # Annotate cross-sex mismatches
        geom_text_repel(data = plot_data[is_cross_sex == TRUE & !is.na(FINNGENID)],
                       aes(label = FINNGENID),
                       color = "purple", size = 2.5, fontface = "bold",
                       box.padding = 0.5, point.padding = 0.3, max.overlaps = Inf) +
        scale_color_manual(values = c(
            "Within-Sex ID Shuffle" = "#2E7D32",      # Green
            "Sex Label Misannotation" = "#FF9800"    # Orange
        ), name = "Error Category") +
        labs(
            title = "gPC1 vs Mean Absolute Z-Score",
            subtitle = paste("Outliers highlighted in red. Cross-sex mismatches (purple diamonds) annotated.",
                           "| Total outliers:", sum(plot_data$Outlier, na.rm = TRUE)),
            x = "gPC1 (First Genetic Principal Component)",
            y = "Mean Absolute Z-Score",
            color = "Error Category"
        ) +
        theme_bw() +
        theme(
            plot.title = element_text(size = 14, face = "bold"),
            plot.subtitle = element_text(size = 11),
            legend.position = "right"
        )

    ggsave(output_path, p, width = 12, height = 8)
    log_info("Saved gPC1 vs MeanAbsZ scatter plot to: {output_path}")
    return(p)
}

# Create scatter plot: MeanAbsZ vs Genetic Distance
create_meanabsz_vs_distance_plot <- function(
    sample_stats,          # From pQTL detection
    ground_truth,          # Error tracking with genetic_distance
    output_path            # File path for saving plot
) {
    log_info("Creating MeanAbsZ vs Genetic Distance scatter plot...")

    # Merge data
    plot_data <- merge(
        sample_stats[, .(SampleID, MeanAbsZ)],
        ground_truth[, .(SAMPLE_ID, category, swap_type, genetic_distance)],
        by.x = "SampleID", by.y = "SAMPLE_ID",
        all.x = TRUE
    )

    # Create category labels and shapes
    plot_data[, error_category := fifelse(
        category == "Category1_WithinSex", "Within-Sex ID Shuffle",
        fifelse(category == "Category2_SexLabel", "Sex Label Misannotation",
        fifelse(category == "Category3_CrossSex", "Cross-Sex ID Shuffle",
        "Clean"))
    )]

    plot_data[, point_shape := fifelse(
        category == "Category3_CrossSex", 17,  # Triangle for cross-sex
        fifelse(category != "Clean", 16, 21)   # Circle for errors, open circle for clean
    )]

    plot_data <- plot_data[!is.na(MeanAbsZ) & !is.na(genetic_distance) &
                          is.finite(MeanAbsZ) & is.finite(genetic_distance)]

    # Create scatter plot
    p <- ggplot(plot_data, aes(x = genetic_distance, y = MeanAbsZ)) +
        # Clean samples (background)
        geom_point(data = plot_data[error_category == "Clean"],
                  aes(shape = factor(point_shape)), color = "gray70", alpha = 0.5, size = 2) +
        # Error samples (foreground)
        geom_point(data = plot_data[error_category != "Clean"],
                  aes(color = error_category, shape = factor(point_shape)),
                  size = 3, alpha = 0.8, stroke = 1.2) +
        # Annotate cross-sex mismatches
        geom_text_repel(data = plot_data[category == "Category3_CrossSex"],
                       aes(label = FINNGENID),
                       color = "purple", size = 2.5, fontface = "bold",
                       box.padding = 0.5, point.padding = 0.3, max.overlaps = Inf) +
        scale_color_manual(values = c(
            "Within-Sex ID Shuffle" = "#2E7D32",      # Green
            "Sex Label Misannotation" = "#FF9800",    # Orange
            "Cross-Sex ID Shuffle" = "#9C27B0"        # Purple
        )) +
        scale_shape_manual(values = c(
            "16" = 16,  # Circle
            "17" = 17,  # Triangle
            "21" = 21   # Open circle
        ), guide = "none") +
        labs(
            title = "Mean Absolute Z-Score vs Genetic Distance",
            subtitle = "Error samples highlighted by category. Cross-sex mismatches annotated.",
            x = "Genetic Distance (Euclidean in PC space)",
            y = "Mean Absolute Z-Score",
            color = "Error Category"
        ) +
        theme_bw() +
        theme(
            plot.title = element_text(size = 14, face = "bold"),
            plot.subtitle = element_text(size = 11),
            legend.position = "right"
        )

    ggsave(output_path, p, width = 12, height = 8)
    log_info("Saved scatter plot to: {output_path}")
    return(p)
}

# Function to run validation tests on synthetic dataset (DIRECT IMPLEMENTATION)
# This implements sex and pQTL detection directly without sourcing the production scripts
run_validation_tests <- function(synthetic_npx_path, synthetic_metadata_path, ground_truth,
                                  batch_id, config, test_output_dir) {
    log_info("=== GOVERNANCE TEST: Running Validation Tests (Direct Implementation) ===")
    log_info("Running validation tests on synthetic dataset...")

    # Direct implementation - no file replacement needed
    # We work directly with synthetic data

    # Run validation tests
    validation_result <- list(
        n_sex_detected = 0,
        n_sex_expected = sum(ground_truth$expect_sex_flag),
        n_pqtl_detected = 0,
        n_pqtl_expected = sum(ground_truth$expect_pqtl_flag),
        sex_detection_rate = 0,
        pqtl_detection_rate = 0,
        sex_detected_samples = character(0),
        pqtl_detected_samples = character(0)
    )

    tryCatch({
        # ====================================================================
        # PART 1: SEX DETECTION ON SYNTHETIC DATA (Direct Implementation)
        # ====================================================================
        log_info("--- PART 1: Sex Detection on Synthetic Data (Direct Implementation) ---")

        # Load synthetic data
        log_info("Loading synthetic NPX matrix and metadata...")
        npx_synthetic <- readRDS(synthetic_npx_path)
        metadata_synthetic <- readRDS(synthetic_metadata_path)

        log_info("Synthetic NPX matrix: {nrow(npx_synthetic)} samples x {ncol(npx_synthetic)} proteins")
        log_info("Synthetic metadata: {nrow(metadata_synthetic)} samples")

        # Apply Inverse Rank Normalization (IRN) per protein if enabled
        # IRN transforms each protein to N(0,1) distribution, making z-score calculation more robust
        pqtl_config <- config$parameters$pqtl_outliers
        apply_irn <- tryCatch(
            pqtl_config$apply_irn %||% TRUE,  # Default: enabled
            error = function(e) TRUE
        )

        if (isTRUE(apply_irn)) {
            npx_synthetic <- inverse_rank_normalize(npx_synthetic)
            log_info("IRN applied: Synthetic NPX matrix transformed to N(0,1) per protein")
        } else {
            log_info("IRN disabled: Using raw NPX values for z-score calculation")
        }

        # Identify samples with sex errors (Category 1 & 2) - exclude from training
        samples_with_sex_errors <- ground_truth[expect_sex_flag == TRUE]$SAMPLE_ID
        log_info("Samples with sex errors (excluded from training): {length(samples_with_sex_errors)}")

        # Load covariates to get original sex information
        log_info("Loading covariates for sex information...")
        covariates <- read_gzipped_file(config$covariates$covariate_file)

        # Extract sex information
        if ("SEX_IMPUTED" %in% names(covariates)) {
            sex_raw <- covariates$SEX_IMPUTED
            sex_col_used <- "SEX_IMPUTED"
        } else if ("SEX" %in% names(covariates)) {
            sex_raw <- covariates$SEX
            sex_col_used <- "SEX"
        } else {
            stop("No sex column found in covariates")
        }

        # Create sex_info mapping FINNGENID to genetic sex
        sex_info_cov <- data.table(
            FINNGENID = covariates$IID,
            genetic_sex_raw = sex_raw
        )

        # Convert to PLINK coding (1=male, 2=female)
        if (sex_col_used == "SEX_IMPUTED") {
            sex_info_cov[, genetic_sex := ifelse(genetic_sex_raw == 1, 2, ifelse(genetic_sex_raw == 0, 1, NA))]
        } else {
            sex_info_cov[, genetic_sex := ifelse(tolower(genetic_sex_raw) == "female", 2,
                                                ifelse(tolower(genetic_sex_raw) == "male", 1, NA))]
        }

        # ============================================================================
        # CRITICAL: Use genetic_sex from SYNTHETIC METADATA, not from covariates!
        # ============================================================================
        # The synthetic metadata has the MODIFIED sex labels (flipped for sex error samples).
        # If we look up sex from covariates based on SWAPPED FINNGENID, we get the wrong sex
        # because within-sex shuffling means the swap partner has the same sex.
        #
        # Example:
        # - Sample A (male proteomics) gets swapped FINNGENID from Sample B (also male)
        # - We flip Sample A's sex label to "female"
        # - If we look up sex from covariates for B's FINNGENID, we get "male"
        # - Model predicts "male" (correct for proteomics), compares to "male" → NO MISMATCH!
        # - But we SHOULD compare to "female" (the flipped label) → MISMATCH!
        # ============================================================================

        # Use genetic_sex from synthetic metadata (which has flipped sex labels)
        if ("genetic_sex" %in% names(metadata_synthetic)) {
            log_info("Using genetic_sex from synthetic metadata (with flipped sex labels for error samples)")
            metadata_with_sex <- copy(metadata_synthetic[, .(SAMPLE_ID, FINNGENID, genetic_sex)])
            # Convert to PLINK coding (1=male, 2=female) if needed
            if (is.character(metadata_with_sex$genetic_sex)) {
                metadata_with_sex[, genetic_sex := ifelse(tolower(genetic_sex) == "female", 2,
                                                         ifelse(tolower(genetic_sex) == "male", 1, NA))]
            }
        } else {
            # Fallback: merge with covariates (this is the OLD behavior, may not detect errors correctly)
            log_warn("genetic_sex not in synthetic metadata. Falling back to covariates lookup (may not detect sex errors correctly).")
            metadata_with_sex <- merge(
                metadata_synthetic[, .(SAMPLE_ID, FINNGENID)],
                sex_info_cov[, .(FINNGENID, genetic_sex)],
                by = "FINNGENID",
                all.x = TRUE
            )
        }

        # Add age, BMI, smoking from covariates
        age_col <- if("BL_AGE" %in% names(covariates)) "BL_AGE" else
                   if("AGE_AT_DEATH_OR_END_OF_FOLLOWUP" %in% names(covariates)) "AGE_AT_DEATH_OR_END_OF_FOLLOWUP" else NA_character_
        bmi_col <- if("BMI_IRN" %in% names(covariates)) "BMI_IRN" else if("BMI" %in% names(covariates)) "BMI" else NA_character_
        smoke_col <- if("harmonized_current_smoker" %in% names(covariates)) "harmonized_current_smoker" else
                     if("CURRENT_SMOKER" %in% names(covariates)) "CURRENT_SMOKER" else NA_character_

        if (!is.na(age_col)) {
            age_map <- setNames(covariates[[age_col]], covariates$IID)
            metadata_with_sex[, age := age_map[FINNGENID]]
        } else {
            metadata_with_sex[, age := NA_real_]
        }

        if (!is.na(bmi_col)) {
            bmi_map <- setNames(covariates[[bmi_col]], covariates$IID)
            metadata_with_sex[, bmi := bmi_map[FINNGENID]]
        } else {
            metadata_with_sex[, bmi := NA_real_]
        }

        if (!is.na(smoke_col)) {
            smoke_map <- setNames(covariates[[smoke_col]], covariates$IID)
            metadata_with_sex[, smoking := smoke_map[FINNGENID]]
        } else {
            metadata_with_sex[, smoking := NA_real_]
        }

        # Filter to samples in NPX matrix
        sample_ids <- rownames(npx_synthetic)
        metadata_with_sex <- metadata_with_sex[SAMPLE_ID %in% sample_ids]

        # Create pheno_dt for modeling
        pheno_dt <- merge(
            data.table(SAMPLE_ID = sample_ids),
            metadata_with_sex,
            by = "SAMPLE_ID",
            all.x = TRUE
        )

        # Exclude samples with sex errors from training (Category 1 & 2)
        training_mask <- !pheno_dt$SAMPLE_ID %in% samples_with_sex_errors

        log_info("Training samples (excluding sex errors): {sum(training_mask)}")
        log_info("Prediction samples (including sex errors): {nrow(pheno_dt)}")

        # Exclude samples with missing sex, F64, chromosomal abnormalities
        excl_missing_sex <- is.na(pheno_dt$genetic_sex)
        excl_f64 <- rep(FALSE, nrow(pheno_dt))
        if ("F64" %in% names(metadata_synthetic)) {
            f64_map <- setNames(metadata_synthetic$F64, metadata_synthetic$SAMPLE_ID)
            excl_f64 <- !is.na(f64_map[pheno_dt$SAMPLE_ID]) & f64_map[pheno_dt$SAMPLE_ID]
        }
        excl_chrom <- rep(FALSE, nrow(pheno_dt))
        chrom_col <- names(metadata_synthetic)[tolower(names(metadata_synthetic)) == "chromosomal_abnormalities"]
        if (length(chrom_col) > 0) {
            chrom_map <- setNames(metadata_synthetic[[chrom_col[1]]], metadata_synthetic$SAMPLE_ID)
            excl_chrom <- !is.na(chrom_map[pheno_dt$SAMPLE_ID]) & chrom_map[pheno_dt$SAMPLE_ID]
        }
        excl_any <- excl_missing_sex | excl_f64 | excl_chrom
        keep_mask <- !excl_any & training_mask

        y_raw <- pheno_dt$genetic_sex
        y_bin <- ifelse(y_raw == 2, 1, ifelse(y_raw == 1, 0, NA))
        keep <- which(!is.na(y_bin) & keep_mask)

        if (length(keep) < 10) {
            stop("Too few samples for training ({length(keep)}). Need at least 10.")
        }

        Y <- npx_synthetic[keep, , drop = FALSE]
        sample_ids_keep <- sample_ids[keep]
        cov_dt <- pheno_dt[keep, .(SAMPLE_ID, FINNGENID, age, bmi, smoking)]
        Xconf <- cbind(Intercept = 1, as.matrix(cov_dt[, .(age, bmi, smoking)]))
        Xconf[is.na(Xconf)] <- 0

        # Nested CV elastic-net (full implementation matching production)
        log_info("Running nested CV elastic-net for sex prediction...")
        folds <- make_folds(y_bin[keep], k = 5, seed = 7)
        prot_names <- colnames(Y)
        metrics <- list()
        preds_all <- rep(NA_real_, length(keep))
        coef_list <- vector("list", length(folds))

        for (fi in seq_along(folds)) {
            log_info("Outer CV fold {fi}/{length(folds)}")
            te_idx <- folds[[fi]]
            tr_idx <- setdiff(seq_along(keep), te_idx)

            Y_tr_raw <- Y[tr_idx, , drop = FALSE]
            Y_te_raw <- Y[te_idx, , drop = FALSE]

            # Impute and residualize
            col_means <- compute_impute_means(Y_tr_raw)
            Y_tr_imp <- impute_matrix_with_means(Y_tr_raw, col_means)
            Y_te_imp <- impute_matrix_with_means(Y_te_raw, col_means)
            res_tr <- residualize_by_confounders(Y_tr_imp, Xconf[tr_idx, , drop = FALSE])
            R_tr <- res_tr$residuals
            coef_B <- res_tr$coef
            R_te <- apply_residualization(Y_te_imp, Xconf[te_idx, , drop = FALSE], coef_B)

            # Scale
            mu <- matrix(colMeans(R_tr, na.rm = TRUE), nrow = 1)
            sdv <- matrix(apply(R_tr, 2, sd, na.rm = TRUE), nrow = 1)
            sdv[sdv == 0 | is.na(sdv)] <- 1
            Z_tr <- sweep(sweep(R_tr, 2, mu, "-"), 2, sdv, "/")
            Z_te <- sweep(sweep(R_te, 2, mu, "-"), 2, sdv, "/")
            conf_tr <- Xconf[tr_idx, -1, drop = FALSE]
            conf_te <- Xconf[te_idx, -1, drop = FALSE]
            conf_mu <- matrix(colMeans(conf_tr, na.rm = TRUE), nrow = 1)
            conf_sd <- matrix(apply(conf_tr, 2, sd, na.rm = TRUE), nrow = 1)
            conf_sd[conf_sd == 0 | is.na(conf_sd)] <- 1
            conf_tr_sc <- sweep(sweep(conf_tr, 2, conf_mu, "-"), 2, conf_sd, "/")
            conf_te_sc <- sweep(sweep(conf_te, 2, conf_mu, "-"), 2, conf_sd, "/")
            X_tr <- cbind(Z_tr, conf_tr_sc)
            X_te <- cbind(Z_te, conf_te_sc)
            y_tr <- y_bin[keep][tr_idx]
            y_te <- y_bin[keep][te_idx]

            # Tune alpha
            pen <- c(rep(1, ncol(Z_tr)), rep(1, ncol(conf_tr)))
            w_tr <- rep(1, length(y_tr))
            p1 <- mean(y_tr == 1); p0 <- mean(y_tr == 0)
            w_tr[y_tr == 1] <- ifelse(p1 > 0, 0.5/p1, 1)
            w_tr[y_tr == 0] <- ifelse(p0 > 0, 0.5/p0, 1)
            alphas <- seq(0, 1, by = 0.1)
            best_auc <- -Inf; best_alpha <- NA; best_lambda <- NA
            set.seed(13+fi)
            inner_foldid <- sample(rep(1:5, length.out = length(y_tr)))

            for (a in alphas) {
                fit_cv <- try(cv.glmnet(X_tr, y_tr, family = "binomial", alpha = a, foldid = inner_foldid,
                                      weights = w_tr, type.measure = "auc", penalty.factor = pen, standardize = FALSE), silent = TRUE)
                if (inherits(fit_cv, "try-error")) next
                auc <- max(fit_cv$cvm)
                lam <- fit_cv$lambda[which.max(fit_cv$cvm)]
                if (auc > best_auc) { best_auc <- auc; best_alpha <- a; best_lambda <- lam }
            }

            if (is.na(best_alpha)) { best_alpha <- 0.5; best_lambda <- NULL }
            fit <- glmnet(X_tr, y_tr, family = "binomial", alpha = best_alpha, lambda = best_lambda,
                         weights = w_tr, penalty.factor = pen, standardize = FALSE)
            pr_te <- as.numeric(predict(fit, newx = X_te, type = "response"))
            preds_all[te_idx] <- pr_te

            # Store coefficients
            co <- as.matrix(coef(fit))
            rn <- rownames(co)
            prot_vec <- setNames(rep(0, length(prot_names)), prot_names)
            common <- intersect(prot_names, rn)
            if (length(common) > 0) {
                prot_vec[common] <- as.numeric(co[common, 1])
            }
            coef_list[[fi]] <- prot_vec

            # Metrics
            roc_obj <- try(pROC::roc(response = y_te, predictor = pr_te, quiet = TRUE), silent = TRUE)
            auc <- if (inherits(roc_obj, "try-error")) NA_real_ else as.numeric(pROC::auc(roc_obj))
            pr_obj <- try(PRROC::pr.curve(scores.class0 = pr_te[y_te == 1], scores.class1 = pr_te[y_te == 0], curve = FALSE), silent = TRUE)
            pr_auc <- if (inherits(pr_obj, "try-error")) NA_real_ else pr_obj$auc.integral
            thr <- pick_threshold(pr_te, y_te)
            metrics[[fi]] <- data.table(fold = fi, auc = auc, pr_auc = pr_auc, th = thr$th, J = thr$J, F1 = thr$F1,
                                       sens = thr$sens, spec = thr$spec, acc = thr$acc,
                                       alpha = best_alpha, lambda = ifelse(is.null(best_lambda), NA_real_, best_lambda))
        }

        metrics_dt <- rbindlist(metrics)
        alpha_pick <- metrics_dt[, .(m_auc = mean(auc, na.rm = TRUE)), by = alpha][order(-m_auc)][1]$alpha

        # Fit final model on all training data
        log_info("Fitting final model on all training data...")
        col_means_all <- compute_impute_means(Y)
        Y_all_imp <- impute_matrix_with_means(Y, col_means_all)
        res_all <- residualize_by_confounders(Y_all_imp, Xconf)
        R_all <- res_all$residuals
        mu_all <- matrix(colMeans(R_all, na.rm = TRUE), nrow = 1)
        sd_all <- matrix(apply(R_all, 2, sd, na.rm = TRUE), nrow = 1)
        sd_all[sd_all == 0 | is.na(sd_all)] <- 1
        Z_all <- sweep(sweep(R_all, 2, mu_all, "-"), 2, sd_all, "/")
        conf_all <- Xconf[, -1, drop = FALSE]
        conf_mu_all <- matrix(colMeans(conf_all, na.rm = TRUE), nrow = 1)
        conf_sd_all <- matrix(apply(conf_all, 2, sd, na.rm = TRUE), nrow = 1)
        conf_sd_all[conf_sd_all == 0 | is.na(conf_sd_all)] <- 1
        conf_all_sc <- sweep(sweep(conf_all, 2, conf_mu_all, "-"), 2, conf_sd_all, "/")
        X_all <- cbind(Z_all, conf_all_sc)
        y_all <- y_bin[keep]
        w_all <- rep(1, length(y_all)); p1 <- mean(y_all == 1); p0 <- mean(y_all == 0)
        w_all[y_all == 1] <- ifelse(p1 > 0, 0.5/p1, 1); w_all[y_all == 0] <- ifelse(p0 > 0, 0.5/p0, 1)
        pen_all <- c(rep(1, ncol(Z_all)), rep(1, ncol(conf_all_sc)))
        fit_all <- cv.glmnet(X_all, y_all, family = "binomial", alpha = alpha_pick,
                            weights = w_all, penalty.factor = pen_all, standardize = FALSE,
                            type.measure = "auc", nfolds = 5)
        pr_all <- as.numeric(predict(fit_all, newx = X_all, s = "lambda.1se", type = "response"))
        th_global <- pick_threshold(preds_all, y_bin[keep])$th

        # Predict on ALL samples (including those with sex errors)
        log_info("Predicting on all samples (including sex error samples)...")
        all_sample_ids <- rownames(npx_synthetic)
        all_pheno_dt <- merge(
            data.table(SAMPLE_ID = all_sample_ids),
            metadata_with_sex,
            by = "SAMPLE_ID",
            all.x = TRUE
        )

        # Prepare data for all samples
        Y_all_samples <- npx_synthetic
        all_cov_dt <- all_pheno_dt[, .(SAMPLE_ID, FINNGENID, age, bmi, smoking)]
        all_Xconf <- cbind(Intercept = 1, as.matrix(all_cov_dt[, .(age, bmi, smoking)]))
        all_Xconf[is.na(all_Xconf)] <- 0

        # Impute and residualize all samples using training statistics
        Y_all_imp_full <- impute_matrix_with_means(Y_all_samples, col_means_all)
        R_all_full <- apply_residualization(Y_all_imp_full, all_Xconf, res_all$coef)
        Z_all_full <- sweep(sweep(R_all_full, 2, mu_all, "-"), 2, sd_all, "/")
        conf_all_full <- all_Xconf[, -1, drop = FALSE]
        conf_all_full_sc <- sweep(sweep(conf_all_full, 2, conf_mu_all, "-"), 2, conf_sd_all, "/")
        X_all_full <- cbind(Z_all_full, conf_all_full_sc)

        # Predict
        pr_all_full <- as.numeric(predict(fit_all, newx = X_all_full, s = "lambda.1se", type = "response"))

        # Create predictions data.table
        all_y_raw <- all_pheno_dt$genetic_sex
        all_y_bin <- ifelse(all_y_raw == 2, 1, ifelse(all_y_raw == 1, 0, NA))
        pred_label <- ifelse(pr_all_full >= th_global, 1, 0)
        predicted_sex_str <- ifelse(pred_label == 1, "female", "male")
        genetic_sex_str_all <- ifelse(all_y_raw == 2, "female",
                                     ifelse(all_y_raw == 1, "male", NA_character_))

        preds_dt <- data.table(
            SAMPLE_ID = all_sample_ids,
            FINNGENID = all_cov_dt$FINNGENID,
            predicted_prob = pr_all_full,
            predicted_sex = predicted_sex_str,
            genetic_sex = genetic_sex_str_all,
            mismatch = FALSE,
            sex_outlier = FALSE
        )

        # Identify mismatches using Youden J threshold
        preds_dt[!is.na(genetic_sex) & !is.na(predicted_sex),
                mismatch := (genetic_sex == "male" & predicted_prob >= th_global) |
                           (genetic_sex == "female" & predicted_prob < th_global)]
        preds_dt[, sex_outlier := (genetic_sex == "male" & predicted_prob >= 0.5 & predicted_prob < th_global) |
                            (genetic_sex == "female" & predicted_prob < 0.5 & predicted_prob >= th_global)]

        log_info("Sex predictions complete:")
        log_info("  Mismatches (Youden J threshold = {round(th_global, 3)}): {sum(preds_dt$mismatch, na.rm=TRUE)}")
        log_info("  Sex outliers (0.5 threshold): {sum(preds_dt$sex_outlier, na.rm=TRUE)}")

        # Save sex predictions
        sex_pred_path <- file.path(test_output_dir, "07_07c_sex_predictions_test_case.tsv")
        ensure_output_dir(sex_pred_path)
        fwrite(preds_dt, sex_pred_path, sep = "\t")
        log_info("Saved sex predictions to: {sex_pred_path}")

        # Identify detected sex outliers
        sex_detected_samples <- preds_dt[mismatch == TRUE | sex_outlier == TRUE]$SAMPLE_ID
        validation_result$sex_detected_samples <- unique(sex_detected_samples)
        validation_result$n_sex_detected <- length(validation_result$sex_detected_samples)

        # Generate sex prediction distribution plot (using same function as 05_sex_outliers.R)
        log_info("Generating sex prediction distribution plot...")
        sex_plot_path <- file.path(test_output_dir, "07_07c_sex_predict_distribution_test_case_fg3_batch_02.pdf")
        ensure_output_dir(sex_plot_path)

        create_prediction_distribution_panels(
            df_pred = preds_dt[, .(SAMPLE_ID, FINNGENID, genetic_sex, predicted_prob)],
            title = "Sex Prediction Distribution – Test Case (Elastic Net Logistic Regression)",
            file_out = sex_plot_path,
            youden_thresh = th_global,
            annotate_ids = FALSE
        )
        log_info("Saved sex prediction plot to: {sex_plot_path}")

        # ====================================================================
        # PART 2: pQTL OUTLIER DETECTION ON SYNTHETIC DATA (Direct Implementation)
        # ====================================================================
        log_info("--- PART 2: pQTL Outlier Detection on Synthetic Data (Direct Implementation) ---")

        # Load fine-mapping results
        log_info("Loading fine-mapping results...")
        pqtl_config <- config$parameters$pqtl_outliers
        if (is.null(pqtl_config)) {
            log_warn("No pqtl_outliers configuration found. Skipping pQTL detection.")
        } else {
            finemap_path <- pqtl_config$finemap_path
            genotype_path <- pqtl_config$genotype_path
            z_threshold <- pqtl_config$z_threshold %||% 4
            per_variant_z_threshold <- pqtl_config$per_variant_z_threshold %||% 3  # Threshold for per-variant outlier detection
            top_n <- tryCatch(pqtl_config$top_n_pQTLs, error = function(e) NULL)

            # Load collated fine-mapping results
            # Note: The actual file uses "07" as step_num, not "07b", and filename is "07.b_finemap_collated"
            # Try both possible paths
            collated_path <- get_output_path("05", "07.b_finemap_collated", batch_id, "pqtl", "tsv", config = config)
            if (!file.exists(collated_path) || file.size(collated_path) == 0) {
                # Try with step_num "07" instead
                collated_path_alt <- get_output_path("03", "07.b_finemap_collated", batch_id, "pqtl", "tsv", config = config)
                if (file.exists(collated_path_alt) && file.size(collated_path_alt) > 0) {
                    collated_path <- collated_path_alt
                    log_info("Found collated file at alternative path: {collated_path}")
                }
            }

            if (file.exists(collated_path) && file.size(collated_path) > 0) {
                top_variants <- fread(collated_path)
                log_info("Loaded {nrow(top_variants)} variants from collated fine-mapping results")

                # Apply MAF filter if configured (for governance test)
                pqtl_selection_config <- tryCatch(test_config$pqtl_selection, error = function(e) NULL)
                min_maf_threshold <- tryCatch(
                    pqtl_selection_config$min_maf %||% NULL,
                    error = function(e) NULL
                )

                if (!is.null(min_maf_threshold) && is.numeric(min_maf_threshold) && min_maf_threshold > 0) {
                    if ("maf" %in% names(top_variants)) {
                        n_before_maf <- nrow(top_variants)
                        top_variants <- top_variants[!is.na(maf) & maf > min_maf_threshold]
                        n_after_maf <- nrow(top_variants)
                        log_info("Applied MAF > {min_maf_threshold*100}% filter: {n_before_maf} -> {n_after_maf} variants")

                        if (n_after_maf == 0) {
                            log_warn("No variants remaining after MAF filter. Using all variants.")
                            top_variants <- fread(collated_path)  # Reload original
                        }
                    } else {
                        log_warn("MAF column not found in collated file. Cannot apply MAF filter.")
                    }
                } else {
                    log_info("MAF filter not configured or disabled. Using all variants.")
                }

                # Filter to top 100 ONLY if MAF filter is NOT configured
                # If MAF filter is configured, use all variants that pass the MAF threshold
                if (is.null(min_maf_threshold) || !is.numeric(min_maf_threshold) || min_maf_threshold <= 0) {
                    # Only apply top-100 filter when MAF filter is disabled
                    if (nrow(top_variants) > 100) {
                        # If included_in_top_n column exists, use it; otherwise take first 100
                        if ("included_in_top_n" %in% names(top_variants)) {
                            top_variants_filtered <- top_variants[included_in_top_n == TRUE]
                            if (nrow(top_variants_filtered) > 100) {
                                top_variants_filtered <- top_variants_filtered[1:100]
                            }
                            top_variants <- top_variants_filtered
                        } else {
                            top_variants <- top_variants[1:100]
                        }
                        log_info("Filtered to top {nrow(top_variants)} variants for pQTL detection")
                    }
                } else {
                    log_info("MAF filter is active. Using all {nrow(top_variants)} variants that pass MAF > {min_maf_threshold*100}% threshold")
                }

                # Apply max_pqtls limit if configured
                max_pqtls <- tryCatch(
                    pqtl_selection_config$max_pqtls %||% NULL,
                    error = function(e) NULL
                )
                if (!is.null(max_pqtls) && is.numeric(max_pqtls) && max_pqtls > 0 && nrow(top_variants) > max_pqtls) {
                    log_info("Applying max_pqtls limit: {nrow(top_variants)} -> {max_pqtls} variants")
                    top_variants <- top_variants[1:max_pqtls]
                }
            } else {
                log_warn("Collated fine-mapping results not found at: {collated_path}")
                log_warn("Skipping pQTL detection.")
                top_variants <- NULL
            }

            if (!is.null(top_variants) && nrow(top_variants) > 0) {
                # Apply proper ranking (matching 07b_pqtl_outliers.R)
                # If heterozygosity column exists, use heterozygosity-weighted ranking
                # Otherwise, use MAF-weighted ranking or standard ranking
                if ("heterozygosity" %in% names(top_variants) && "maf" %in% names(top_variants)) {
                    log_info("Applying heterozygosity-weighted ranking: (-log10(p) * abs(beta) * maf) / (heterozygosity + 0.01)")
                    # Add epsilon to heterozygosity to avoid division by zero
                    top_variants[, heterozygosity_adj := heterozygosity + 0.01]
                    # Calculate composite score with heterozygosity weighting
                    top_variants[, composite_score := (-log10(p) * abs(beta) * maf) / heterozygosity_adj]
                    log_info("Ranked variants using heterozygosity-weighted composite score")
                } else if ("maf" %in% names(top_variants)) {
                    log_info("Applying MAF-weighted ranking: -log10(p) * abs(beta) * maf")
                    top_variants[, composite_score := -log10(p) * abs(beta) * maf]
                    log_info("Ranked variants using MAF-weighted composite score")
                } else {
                    log_info("Applying standard ranking: -log10(p) * abs(beta)")
                    top_variants[, composite_score := -log10(p) * abs(beta)]
                    log_info("Ranked variants using standard composite score")
                }

                # Sort by composite score (descending)
                setorder(top_variants, -composite_score)

                # Apply top_n selection if specified (further filtering beyond top 100)
                if (!is.null(top_n) && is.numeric(top_n) && top_n >= 100 && top_n < nrow(top_variants)) {
                    log_info("Selecting top {top_n} pQTLs from {nrow(top_variants)} available...")
                    top_variants <- top_variants[1:top_n]
                    log_info("Selected top {nrow(top_variants)} variants based on composite score ranking")
                    if ("heterozygosity" %in% names(top_variants)) {
                        log_info("Selected variants heterozygosity range: {round(min(top_variants$heterozygosity, na.rm=TRUE), 4)} - {round(max(top_variants$heterozygosity, na.rm=TRUE), 4)}")
                    }
                }

                # Generate PLINK genotype files for test samples
                log_info("Generating PLINK genotype files for test samples...")

                # Get test sample FINNGENIDs from synthetic metadata (with swaps)
                test_sample_finngenids <- unique(metadata_synthetic$FINNGENID[!is.na(metadata_synthetic$FINNGENID)])
                log_info("Test samples: {length(test_sample_finngenids)} unique FINNGENIDs")

                # Create keep file for PLINK (test samples only)
                test_geno_dir <- file.path(test_output_dir, "genotype_files")
                if (!dir.exists(test_geno_dir)) {
                    dir.create(test_geno_dir, recursive = TRUE)
                }
                keep_file <- file.path(test_geno_dir, "test_samples_keep.txt")
                keep_dt <- data.table(FID = test_sample_finngenids, IID = test_sample_finngenids)
                fwrite(keep_dt, keep_file, sep = "\t", col.names = FALSE)

                # Create variant list file
                var_file <- file.path(test_geno_dir, "variants.snplist")
                # Use rsid column (should exist in collated file)
                if (!"rsid" %in% names(top_variants)) {
                    log_error("rsid column not found in top_variants. Available columns: {paste(names(top_variants), collapse=', ')}")
                    stop("Cannot create variant list: rsid column missing")
                }
                variant_list <- unique(top_variants$rsid[!is.na(top_variants$rsid)])
                if (length(variant_list) == 0) {
                    log_error("No valid variants found in top_variants")
                    stop("Cannot create variant list: no valid rsids")
                }
                writeLines(variant_list, var_file)
                log_info("Created variant list file with {length(variant_list)} variants: {var_file}")

                # Check for cached genotype files (same logic as main function)
                base_dir <- config$output$base_dir %||% Sys.getenv("PIPELINE_OUTPUT_DIR", "output")
                temp_work_dir <- file.path(base_dir, "temp_work")
                cached_geno_dir <- NULL
                genotype_path <- pqtl_config$genotype_path
                base_name <- basename(genotype_path)

                # Look for cached genotype files (check pqtl_cache FIRST - more reliable)
                # Note: freq_sparse_* directories may have corrupted files (full .bim but sparse .bed)
                if (dir.exists(temp_work_dir)) {
                    # Helper function to validate cache consistency
                    validate_cache <- function(cache_dir, base_nm) {
                        bed_f <- file.path(cache_dir, paste0(base_nm, ".bed"))
                        bim_f <- file.path(cache_dir, paste0(base_nm, ".bim"))
                        fam_f <- file.path(cache_dir, paste0(base_nm, ".fam"))

                        if (!file.exists(bed_f) || !file.exists(bim_f) || !file.exists(fam_f)) {
                            return(FALSE)
                        }

                        # Check file sizes are non-zero
                        if (file.size(bed_f) == 0 || file.size(bim_f) == 0 || file.size(fam_f) == 0) {
                            return(FALSE)
                        }

                        # Validate BED file size matches expected size based on .bim and .fam
                        # BED format: 3 magic bytes + M variants * ceil(N samples / 4) bytes
                        tryCatch({
                            n_samples <- length(readLines(fam_f))
                            n_variants <- length(readLines(bim_f))
                            bytes_per_variant <- ceiling(n_samples / 4)
                            expected_bed_size <- 3 + n_variants * bytes_per_variant
                            actual_bed_size <- file.size(bed_f)

                            # Allow some tolerance for rounding
                            is_valid <- abs(actual_bed_size - expected_bed_size) <= n_variants
                            if (!is_valid) {
                                log_warn("Cache validation failed: BED size {actual_bed_size} != expected {expected_bed_size} (variants: {n_variants}, samples: {n_samples})")
                            }
                            return(is_valid)
                        }, error = function(e) {
                            log_warn("Cache validation error: {e$message}")
                            return(FALSE)
                        })
                    }

                    # PRIORITY 1: Check pqtl_cache (has correctly matched sparse BED/BIM/FAM)
                    cache_base <- file.path(temp_work_dir, "pqtl_cache", batch_id)
                    if (dir.exists(cache_base)) {
                        cache_dirs <- list.dirs(cache_base, full.names = TRUE, recursive = FALSE)
                        if (length(cache_dirs) > 0) {
                            # Sort by modification time (most recent first)
                            cache_dirs <- cache_dirs[order(file.mtime(cache_dirs), decreasing = TRUE)]
                            for (cache_dir_candidate in cache_dirs) {
                                if (validate_cache(cache_dir_candidate, base_name)) {
                                    cached_geno_dir <- cache_dir_candidate
                                    log_info("Found valid cached genotype files in pqtl_cache: {cached_geno_dir}")
                                    break
                                }
                            }
                        }
                    }

                    # PRIORITY 2: Check plink_input_* directories (created during sparse extraction)
                    if (is.null(cached_geno_dir)) {
                        plink_dirs <- list.dirs(temp_work_dir, full.names = TRUE, recursive = FALSE)
                        plink_dirs <- plink_dirs[grepl("plink_input_", basename(plink_dirs))]
                        if (length(plink_dirs) > 0) {
                            plink_dirs <- plink_dirs[order(file.mtime(plink_dirs), decreasing = TRUE)]
                            for (plink_dir_candidate in plink_dirs) {
                                if (validate_cache(plink_dir_candidate, base_name)) {
                                    cached_geno_dir <- plink_dir_candidate
                                    log_info("Found valid cached genotype files in plink_input: {cached_geno_dir}")
                                    break
                                }
                            }
                        }
                    }

                    # DO NOT use freq_sparse_* directories - they often have corrupted files
                    # (full .bim from original data but sparse .bed from extraction)
                }

                # Determine genotype source
                if (startsWith(genotype_path, "gs://")) {
                    if (!is.null(cached_geno_dir)) {
                        local_genotype_path <- file.path(cached_geno_dir, base_name)
                        log_info("Using cached genotype files from: {local_genotype_path}")
                    } else {
                        log_warn("No cached genotype files found. pQTL detection may fail.")
                        log_warn("Expected cached files in: {temp_work_dir}/freq_sparse_* or {temp_work_dir}/pqtl_cache/{batch_id}/*")
                        local_genotype_path <- NULL
                    }
                } else {
                    local_genotype_path <- genotype_path
                }

                if (!is.null(local_genotype_path) && file.exists(paste0(local_genotype_path, ".bed"))) {
                    # Run PLINK to extract genotypes for test samples and variants
                    out_prefix <- file.path(test_geno_dir, "test_genotypes")
                    is_bed <- file.exists(paste0(local_genotype_path, ".bed"))
                    flag <- if (is_bed) "--bfile" else "--pfile"

                    plink_cmd <- sprintf(
                        "plink2 %s %s --keep %s --extract %s --export A --out %s",
                        flag, local_genotype_path, keep_file, var_file, out_prefix
                    )

                    log_info("Running PLINK: {plink_cmd}")
                    plink_result <- system(plink_cmd, intern = TRUE)
                    plink_exit_code <- attr(plink_result, "status") %||% 0

                    if (plink_exit_code != 0) {
                        log_warn("PLINK command failed with exit code: {plink_exit_code}")
                        log_warn("PLINK output (last 20 lines): {paste(tail(plink_result, 20), collapse = '\\n')}")
                    } else {
                        log_info("PLINK command completed successfully")
                    }

                    raw_file <- paste0(out_prefix, ".raw")
                    if (file.exists(raw_file) && file.size(raw_file) > 0) {
                        dt_geno <- fread(raw_file)
                        log_info("Loaded genotype data: {nrow(dt_geno)} samples")

                        # ================================================================
                        # ENHANCED pQTL Z-SCORE CALCULATION (matching production logic)
                        # ================================================================
                        log_info("Calculating pQTL z-scores (enhanced implementation)...")

                        # Add rsid_for_matching if not present (handles chr23→chrX conversion)
                        if (!"rsid_for_matching" %in% names(top_variants)) {
                            top_variants[, rsid_for_matching := rsid]
                        }

                        # Filter to variants with MAF > 5% (matching production)
                        maf_threshold <- 0.05
                        if ("maf" %in% names(top_variants)) {
                            n_before <- nrow(top_variants)
                            top_variants_filtered <- top_variants[!is.na(maf) & maf > maf_threshold]
                            n_after <- nrow(top_variants_filtered)
                            log_info("Applied MAF > {maf_threshold*100}% filter: {n_before} -> {n_after} variants")
                            if (n_after == 0) {
                                log_warn("No variants after MAF filter. Using all variants.")
                                top_variants_filtered <- top_variants
                            }
                        } else {
                            log_warn("MAF column not found. Using all variants.")
                            top_variants_filtered <- top_variants
                        }

                        results <- list()
                        n_variants_processed <- 0

                        for (i in seq_len(nrow(top_variants_filtered))) {
                            prot <- top_variants_filtered[i]$trait
                            rsid <- top_variants_filtered[i]$rsid
                            rsid_match <- if ("rsid_for_matching" %in% names(top_variants_filtered)) {
                                top_variants_filtered[i]$rsid_for_matching
                            } else {
                                rsid
                            }

                            if (!prot %in% colnames(npx_synthetic)) {
                                next
                            }

                            # Find genotype column (use rsid_for_matching for chr23→chrX handling)
                            geno_col <- grep(paste0("^", rsid_match, "_"), names(dt_geno), value = TRUE)
                            if (length(geno_col) == 0) {
                                # Try without underscore suffix
                                geno_col <- grep(paste0("^", rsid_match), names(dt_geno), value = TRUE)
                            }
                            if (length(geno_col) == 0) next

                            # Match samples (using swapped FINNGENIDs from synthetic metadata)
                            common_samples <- intersect(dt_geno$IID, metadata_synthetic$FINNGENID)
                            if (length(common_samples) == 0) next

                            # Get sample IDs (not FINNGENIDs) for matching with NPX matrix
                            sample_map <- metadata_synthetic[FINNGENID %in% common_samples, .(FINNGENID, SAMPLE_ID)]
                            sample_map <- sample_map[!is.na(SAMPLE_ID)]

                            if (nrow(sample_map) == 0) next

                            # Get protein values and genotypes
                            df <- data.table(
                                SAMPLE_ID = sample_map$SAMPLE_ID,
                                FINNGENID = sample_map$FINNGENID,
                                Geno = dt_geno[match(sample_map$FINNGENID, IID)][[geno_col[1]]],
                                Protein = as.numeric(npx_synthetic[sample_map$SAMPLE_ID, prot])
                            )

                            # FIX: Flip genotype to match production implementation
                            # PLINK --export A: 0 = hom ref, 1 = het, 2 = hom alt
                            # After flip: 2 = hom ref, 1 = het, 0 = hom alt
                            # This aligns with beta interpretation (per-alternate-allele effect)
                            df[, Geno := 2 - Geno]

                            df <- df[!is.na(Protein) & !is.na(Geno)]
                            if (nrow(df) == 0) next

                            # Calculate stats per genotype
                            stats <- df[, .(Mean = mean(Protein), SD = sd(Protein)), by = Geno]
                            stats <- stats[!is.na(SD) & SD > 0]
                            if (nrow(stats) == 0) next

                            # Calculate Z-scores
                            df <- merge(df, stats, by = "Geno", all.x = FALSE)
                            if (nrow(df) == 0) next
                            df[, Z := (Protein - Mean) / SD]

                            # Get beta for beta-weighted residual calculation
                            beta_val <- if ("beta" %in% names(top_variants_filtered)) {
                                top_variants_filtered[i]$beta
                            } else {
                                NA_real_
                            }

                            # ================================================================
                            # BETA-WEIGHTED RESIDUAL CALCULATION
                            # ================================================================
                            # For correctly matched samples: residual should be small
                            # For mismatched samples: residual will be large
                            # Predicted protein = population_mean + beta * (genotype - mean_genotype)
                            if (!is.na(beta_val) && is.finite(beta_val) && nrow(df) > 3) {
                                pop_mean_protein <- mean(df$Protein, na.rm = TRUE)
                                mean_geno <- mean(df$Geno, na.rm = TRUE)

                                # Calculate predicted protein based on genotype and beta
                                df[, Predicted := pop_mean_protein + beta_val * (Geno - mean_geno)]

                                # Calculate raw residual
                                df[, Residual := Protein - Predicted]

                                # Standardize residual using population SD of residuals
                                residual_sd <- sd(df$Residual, na.rm = TRUE)
                                if (!is.na(residual_sd) && residual_sd > 0) {
                                    df[, StdResidual := Residual / residual_sd]
                                } else {
                                    df[, StdResidual := NA_real_]
                                }
                            } else {
                                df[, Predicted := NA_real_]
                                df[, Residual := NA_real_]
                                df[, StdResidual := NA_real_]
                            }

                            results[[length(results) + 1]] <- df[, .(
                                SampleID = SAMPLE_ID,
                                Protein = Protein,
                                Geno = Geno,
                                rsid = rsid,
                                Z = Z,
                                beta = beta_val,
                                Predicted = Predicted,
                                Residual = Residual,
                                StdResidual = StdResidual
                            )]
                            n_variants_processed <- n_variants_processed + 1
                        }

                        log_info("Processed {n_variants_processed} variants for z-score calculation")

                        if (length(results) > 0) {
                            dt_z <- rbindlist(results)
                            log_info("Successfully calculated Z-scores for {length(unique(dt_z$rsid))} variants")

                            # Count variants with valid residuals
                            n_with_residuals <- sum(!is.na(dt_z$StdResidual))
                            log_info("Variants with valid beta-weighted residuals: {n_with_residuals}")

                            # ================================================================
                            # ROBUST MULTI-METRIC OUTLIER DETECTION
                            # ================================================================
                            log_info("=== ROBUST MULTI-METRIC OUTLIER DETECTION ===")

                            # Get residual detection config
                            residual_config <- tryCatch(pqtl_config$residual_detection, error = function(e) NULL)
                            if (is.null(residual_config)) {
                                residual_config <- list(
                                    enabled = TRUE,
                                    threshold_method = "mad",
                                    threshold_k = 3.0,
                                    min_variants = 10,
                                    aggregation = "msr"
                                )
                            }

                            # Per-variant Z-score threshold (configurable, default: 3)
                            z_outlier_threshold <- per_variant_z_threshold
                            log_info("Per-variant Z-score threshold: |Z| > {z_outlier_threshold}")

                            # Calculate standard metrics
                            sample_stats <- dt_z[, .(
                                MeanAbsZ = mean(abs(Z), na.rm = TRUE),
                                MaxAbsZ = max(abs(Z), na.rm = TRUE),
                                MedianAbsZ = median(abs(Z), na.rm = TRUE),
                                N_Variants = .N,
                                OutlierVariantCount = sum(abs(Z) > z_outlier_threshold, na.rm = TRUE),
                                # Beta-weighted residual metrics
                                N_ValidResiduals = sum(!is.na(StdResidual)),
                                MeanSquaredResidual = mean(StdResidual^2, na.rm = TRUE),
                                MedianAbsResidual = median(abs(StdResidual), na.rm = TRUE),
                                MaxAbsResidual = max(abs(StdResidual), na.rm = TRUE)
                            ), by = SampleID]

                            # Calculate OutlierVariantFraction
                            sample_stats[, OutlierVariantFraction := OutlierVariantCount / N_Variants]

                            # ================================================================
                            # ROBUST THRESHOLD CALCULATION (MAD-based)
                            # ================================================================
                            log_info("Calculating robust thresholds using {residual_config$threshold_method} method...")

                            threshold_k <- residual_config$threshold_k %||% 3.0
                            min_variants <- residual_config$min_variants %||% 10

                            # Helper function for MAD-based threshold
                            calc_mad_threshold <- function(x, k = 3) {
                                x <- x[!is.na(x) & is.finite(x)]
                                if (length(x) < 5) return(Inf)
                                med <- median(x)
                                mad_val <- median(abs(x - med))
                                # 1.4826 is the scaling factor to make MAD consistent with SD for normal distributions
                                threshold <- med + k * mad_val * 1.4826
                                return(threshold)
                            }

                            # Helper function for SD-based threshold
                            calc_sd_threshold <- function(x, k = 3) {
                                x <- x[!is.na(x) & is.finite(x)]
                                if (length(x) < 5) return(Inf)
                                threshold <- mean(x) + k * sd(x)
                                return(threshold)
                            }

                            # Select threshold function - use MAD by default (more robust)
                            threshold_method <- residual_config$threshold_method %||% "mad"
                            calc_threshold <- if (threshold_method == "mad") {
                                calc_mad_threshold
                            } else {
                                calc_sd_threshold
                            }

                            # ================================================================
                            # THRESHOLD CALCULATION FOR EACH METRIC
                            # ================================================================

                            # Threshold 1: MeanAbsZ
                            cutoff_mean_z <- calc_threshold(sample_stats$MeanAbsZ, threshold_k)

                            # Threshold 2: MaxAbsZ
                            cutoff_max_z <- calc_threshold(sample_stats$MaxAbsZ, threshold_k)

                            # Threshold 3: OutlierVariantFraction
                            cutoff_outlier_frac <- calc_threshold(sample_stats$OutlierVariantFraction, threshold_k)
                            # Ensure minimum of 10% for outlier fraction
                            cutoff_outlier_frac <- max(cutoff_outlier_frac, 0.10)

                            # Threshold 4: Mean Squared Residual (ROBUST - key metric)
                            # Only calculate for samples with enough valid residuals
                            valid_msr <- sample_stats[N_ValidResiduals >= min_variants]$MeanSquaredResidual
                            cutoff_msr <- calc_threshold(valid_msr, threshold_k)

                            # Threshold 5: Median Absolute Residual
                            valid_mar <- sample_stats[N_ValidResiduals >= min_variants]$MedianAbsResidual
                            cutoff_mar <- calc_threshold(valid_mar, threshold_k)

                            log_info("Robust Thresholds (method: {threshold_method}, k: {threshold_k}):")
                            log_info("  MeanAbsZ: > {round(cutoff_mean_z, 3)}")
                            log_info("  MaxAbsZ: > {round(cutoff_max_z, 3)}")
                            log_info("  OutlierVariantFraction: > {round(cutoff_outlier_frac*100, 1)}%")
                            log_info("  MeanSquaredResidual: > {round(cutoff_msr, 3)}")
                            log_info("  MedianAbsResidual: > {round(cutoff_mar, 3)}")

                            # ================================================================
                            # APPLY MULTI-METRIC DETECTION (ANY metric triggers outlier)
                            # ================================================================
                            sample_stats[, Outlier_MeanAbsZ := MeanAbsZ > cutoff_mean_z]
                            sample_stats[, Outlier_MaxAbsZ := MaxAbsZ > cutoff_max_z]
                            sample_stats[, Outlier_VariantFrac := OutlierVariantFraction > cutoff_outlier_frac]
                            sample_stats[, Outlier_MSR := N_ValidResiduals >= min_variants & MeanSquaredResidual > cutoff_msr]
                            sample_stats[, Outlier_MAR := N_ValidResiduals >= min_variants & MedianAbsResidual > cutoff_mar]

                            # ================================================================
                            # COMPOSITE ANOMALY SCORE (Robust multi-metric combination)
                            # ================================================================
                            # Normalize each metric to [0,1] using rank percentile
                            # This is robust to outliers and allows fair comparison across metrics

                            log_info("Calculating composite anomaly score...")

                            # Helper function to calculate rank percentile (0-1)
                            rank_percentile <- function(x) {
                                x[is.na(x) | !is.finite(x)] <- 0
                                ranks <- rank(x, ties.method = "average", na.last = "keep")
                                percentiles <- ranks / (length(x) + 1)
                                return(percentiles)
                            }

                            # Normalize each metric to rank percentile
                            sample_stats[, MeanAbsZ_Pctl := rank_percentile(MeanAbsZ)]
                            sample_stats[, MaxAbsZ_Pctl := rank_percentile(MaxAbsZ)]
                            sample_stats[, OutlierFrac_Pctl := rank_percentile(OutlierVariantFraction)]
                            sample_stats[, MSR_Pctl := rank_percentile(MeanSquaredResidual)]
                            sample_stats[, MAR_Pctl := rank_percentile(MedianAbsResidual)]

                            # Composite score = mean of all normalized metrics
                            # Higher score = more likely to be an outlier
                            sample_stats[, CompositeScore := (MeanAbsZ_Pctl + MaxAbsZ_Pctl + OutlierFrac_Pctl +
                                                              MSR_Pctl + MAR_Pctl) / 5]

                            # Apply threshold based on composite score
                            # Use MAD-based threshold on composite score
                            composite_threshold <- calc_threshold(sample_stats$CompositeScore, threshold_k)
                            # Also provide a percentile-based threshold (e.g., top 10%)
                            percentile_threshold <- quantile(sample_stats$CompositeScore, 0.90, na.rm = TRUE)

                            log_info("Composite Score Statistics:")
                            log_info("  Min: {round(min(sample_stats$CompositeScore, na.rm=TRUE), 3)}")
                            log_info("  Median: {round(median(sample_stats$CompositeScore, na.rm=TRUE), 3)}")
                            log_info("  Mean: {round(mean(sample_stats$CompositeScore, na.rm=TRUE), 3)}")
                            log_info("  Max: {round(max(sample_stats$CompositeScore, na.rm=TRUE), 3)}")
                            log_info("  MAD-based threshold (k={threshold_k}): {round(composite_threshold, 3)}")
                            log_info("  90th percentile threshold: {round(percentile_threshold, 3)}")

                            # ================================================================
                            # TIERED DETECTION (Addresses fundamental limitation)
                            # ================================================================
                            # Within-sex shuffling creates subtle mismatches that may not be
                            # fully detectable through pQTL analysis. We use a tiered approach:
                            # - Tier 1 (High confidence): Flag for removal
                            # - Tier 2 (Moderate confidence): Flag for review
                            # - Tier 3 (Low confidence): Note but don't auto-flag

                            tier1_pctl <- residual_config$tier1_percentile %||% 0.92
                            tier2_pctl <- residual_config$tier2_percentile %||% 0.80
                            tier3_pctl <- residual_config$tier3_percentile %||% 0.65

                            tier1_cutoff <- quantile(sample_stats$CompositeScore, tier1_pctl, na.rm = TRUE)
                            tier2_cutoff <- quantile(sample_stats$CompositeScore, tier2_pctl, na.rm = TRUE)
                            tier3_cutoff <- quantile(sample_stats$CompositeScore, tier3_pctl, na.rm = TRUE)

                            log_info("Tiered Detection Thresholds:")
                            log_info("  Tier 1 (High confidence, {(1-tier1_pctl)*100}% flagged): > {round(tier1_cutoff, 3)}")
                            log_info("  Tier 2 (Moderate confidence, {(1-tier2_pctl)*100}% flagged): > {round(tier2_cutoff, 3)}")
                            log_info("  Tier 3 (Low confidence, {(1-tier3_pctl)*100}% flagged): > {round(tier3_cutoff, 3)}")

                            # Assign tiers
                            sample_stats[, OutlierTier := fifelse(
                                CompositeScore > tier1_cutoff, 1L,
                                fifelse(CompositeScore > tier2_cutoff, 2L,
                                    fifelse(CompositeScore > tier3_cutoff, 3L, 0L))
                            )]

                            # For validation: flag Tier 1 and Tier 2 as outliers
                            sample_stats[, Outlier := OutlierTier %in% c(1L, 2L)]
                            sample_stats[, Outlier_Composite := CompositeScore > composite_threshold]

                            # Also flag by individual metrics for comparison
                            sample_stats[, Outlier_Individual := Outlier_MeanAbsZ | Outlier_MaxAbsZ | Outlier_VariantFrac |
                                                                  Outlier_MSR | Outlier_MAR]

                            log_info("Samples by Tier:")
                            log_info("  Tier 1 (High confidence): {sum(sample_stats$OutlierTier == 1)}")
                            log_info("  Tier 2 (Moderate confidence): {sum(sample_stats$OutlierTier == 2)}")
                            log_info("  Tier 3 (Low confidence): {sum(sample_stats$OutlierTier == 3)}")
                            log_info("  Normal (no flag): {sum(sample_stats$OutlierTier == 0)}")
                            log_info("Using Tier 1 + Tier 2 as outliers for validation")

                            # Count detections by each metric
                            n_by_mean <- sum(sample_stats$Outlier_MeanAbsZ, na.rm = TRUE)
                            n_by_max <- sum(sample_stats$Outlier_MaxAbsZ, na.rm = TRUE)
                            n_by_frac <- sum(sample_stats$Outlier_VariantFrac, na.rm = TRUE)
                            n_by_msr <- sum(sample_stats$Outlier_MSR, na.rm = TRUE)
                            n_by_mar <- sum(sample_stats$Outlier_MAR, na.rm = TRUE)
                            n_by_composite <- sum(sample_stats$Outlier_Composite, na.rm = TRUE)
                            n_by_percentile <- sum(sample_stats$Outlier_Percentile, na.rm = TRUE)
                            n_combined <- sum(sample_stats$Outlier, na.rm = TRUE)

                            log_info("Outliers detected by each metric:")
                            log_info("  MeanAbsZ > {round(cutoff_mean_z, 3)}: {n_by_mean} samples")
                            log_info("  MaxAbsZ > {round(cutoff_max_z, 3)}: {n_by_max} samples")
                            log_info("  OutlierVariantFraction > {round(cutoff_outlier_frac*100, 1)}%: {n_by_frac} samples")
                            log_info("  MeanSquaredResidual > {round(cutoff_msr, 3)}: {n_by_msr} samples")
                            log_info("  MedianAbsResidual > {round(cutoff_mar, 3)}: {n_by_mar} samples")
                            log_info("  CompositeScore > {round(composite_threshold, 3)} (MAD-based): {n_by_composite} samples")
                            log_info("  CompositeScore > {round(percentile_threshold, 3)} (90th percentile): {n_by_percentile} samples")
                            log_info("  Combined (ANY metric OR Composite): {n_combined} samples")

                            # Save pQTL results with all metrics
                            pqtl_outliers_path <- file.path(test_output_dir, "07_07c_pqtl_outliers_test_case.tsv")
                            ensure_output_dir(pqtl_outliers_path)
                            fwrite(sample_stats[Outlier == TRUE], pqtl_outliers_path, sep = "\t")

                            pqtl_stats_path <- file.path(test_output_dir, "07_07c_pqtl_stats_multi_metric_test_case.tsv")
                            fwrite(sample_stats, pqtl_stats_path, sep = "\t")
                            log_info("Saved multi-metric stats to: {pqtl_stats_path}")

                            # Identify detected pQTL outliers
                            pqtl_detected_samples <- sample_stats[Outlier == TRUE]$SampleID
                            validation_result$pqtl_detected_samples <- unique(pqtl_detected_samples)
                            validation_result$n_pqtl_detected <- length(validation_result$pqtl_detected_samples)

                            # ================================================================
                            # GENERATE FACETED PLOTS FOR EACH METRIC
                            # ================================================================
                            log_info("Generating faceted plots for each detection metric...")

                            # Merge with sex predictions for plotting
                            plot_data <- merge(
                                sample_stats,
                                preds_dt[, .(SAMPLE_ID, predicted_prob, genetic_sex, mismatch)],
                                by.x = "SampleID", by.y = "SAMPLE_ID",
                                all.x = TRUE
                            )

                            # Also merge with ground truth to identify error types (including swap_type and FINNGENID)
                            plot_data <- merge(
                                plot_data,
                                ground_truth[, .(SAMPLE_ID, FINNGENID, has_sex_error, has_genotype_error, category, swap_type)],
                                by.x = "SampleID", by.y = "SAMPLE_ID",
                                all.x = TRUE
                            )

                            plot_data <- plot_data[!is.na(genetic_sex) & !is.na(predicted_prob)]

                            if (nrow(plot_data) > 0) {
                                plot_data[, genetic_sex := factor(tolower(as.character(genetic_sex)), levels = c("female", "male"))]

                                # Mark ground truth error samples
                                plot_data[is.na(has_genotype_error), has_genotype_error := FALSE]
                                plot_data[is.na(has_sex_error), has_sex_error := FALSE]
                                plot_data[, is_error_sample := has_genotype_error | has_sex_error]

                                # Helper function to create faceted plot for a metric
                                create_metric_plot <- function(data, y_col, y_label, threshold, title_suffix, file_suffix) {
                                    # Create threshold column name for this specific metric
                                    # Map metric names to their corresponding outlier flag columns
                                    outlier_col <- switch(y_col,
                                        "MeanAbsZ" = "Outlier_MeanAbsZ",
                                        "MaxAbsZ" = "Outlier_MaxAbsZ",
                                        "OutlierVariantFraction" = "Outlier_VariantFrac",
                                        "MeanSquaredResidual" = "Outlier_MSR",
                                        "MedianAbsResidual" = "Outlier_MAR",
                                        "CompositeScore" = "Outlier_Composite",
                                        paste0("Outlier_", gsub("AbsZ", "AbsZ", y_col))  # Fallback: try to construct column name
                                    )

                                    # Verify the column exists, fallback to composite Outlier if not
                                    if (!(outlier_col %in% names(data))) {
                                        log_warn("Outlier column '{outlier_col}' not found for metric '{y_col}'. Using composite 'Outlier' flag.")
                                        outlier_col <- "Outlier"
                                    }

                                    # Calculate plot-specific point_type using the metric-specific outlier flag
                                    # This ensures False Positive/False Negative designations are test-specific
                                    data[, point_type_metric := fifelse(
                                        is_error_sample & get(outlier_col), "True Positive (Error Detected)",
                                        fifelse(is_error_sample & !get(outlier_col), "False Negative (Error Missed)",
                                        fifelse(!is_error_sample & get(outlier_col), "False Positive",
                                        "True Negative (Clean)"))
                                    )]

                                    # Mark cross-sex mismatches
                                    data[, is_cross_sex := category == "Category3_CrossSex"]
                                    data[, point_shape := fifelse(get(y_col) > threshold, "Above threshold", "Below threshold")]
                                    data[, has_border := is_error_sample]

                                    # Color by detection status with clearer labels
                                    color_map <- c(
                                        "True Positive (Error Detected)" = "#2E7D32",  # Green
                                        "False Negative (Error Missed)" = "#C62828",   # Red
                                        "False Positive" = "#FF9800",                   # Orange
                                        "True Negative (Clean)" = "#90CAF9"             # Light blue
                                    )

                                    p <- ggplot(data, aes(x = predicted_prob, y = get(y_col))) +
                                        # Base points with detection status colors (using metric-specific point_type)
                                        geom_point(aes(fill = point_type_metric, shape = point_shape),
                                                  size = 2.5, alpha = 0.7, stroke = 0.5, color = "black") +
                                        # Highlight cross-sex mismatches with purple border and diamond shape
                                        geom_point(data = data[is_cross_sex == TRUE],
                                                  aes(shape = "Cross-Sex Mismatch"),
                                                  fill = NA, color = "purple",
                                                  size = 4, stroke = 2, alpha = 1) +
                                        # Highlight other error samples with black border
                                        geom_point(data = data[is_error_sample == TRUE & is_cross_sex == FALSE],
                                                  aes(shape = point_shape),
                                                  fill = NA, color = "black",
                                                  size = 3.5, stroke = 1.5, alpha = 1) +
                                        # Annotate cross-sex mismatches with FINNGENID
                                        geom_text_repel(data = data[is_cross_sex == TRUE & !is.na(FINNGENID)],
                                                       aes(label = FINNGENID),
                                                       color = "purple", size = 2.5, fontface = "bold",
                                                       box.padding = 0.5, point.padding = 0.3, max.overlaps = Inf) +
                                        facet_wrap(~ genetic_sex, nrow = 2, scales = "free_x") +
                                        scale_fill_manual(
                                            values = color_map,
                                            name = paste0("Detection Status (", y_col, ")"),
                                            labels = c(
                                                "True Positive (Error Detected)" = "✓ Error Detected (Green)",
                                                "False Negative (Error Missed)" = "✗ Error Missed (Red)",
                                                "False Positive" = "⚠ False Positive (Orange)",
                                                "True Negative (Clean)" = "○ Clean Sample (Blue)"
                                            ),
                                            guide = guide_legend(
                                                override.aes = list(shape = 21, size = 4),
                                                title.position = "top",
                                                title.hjust = 0.5
                                            )
                                        ) +
                                        scale_shape_manual(
                                            values = c(
                                                "Above threshold" = 24,              # Triangle up
                                                "Below threshold" = 21,              # Circle
                                                "Cross-Sex Mismatch" = 23            # Diamond
                                            ),
                                            name = paste0("Threshold: ", round(threshold, 3)),
                                            labels = c(
                                                "Above threshold" = paste("Above threshold (", round(threshold, 3), ")"),
                                                "Below threshold" = paste("Below threshold (", round(threshold, 3), ")"),
                                                "Cross-Sex Mismatch" = "Cross-Sex Mismatch (Purple diamond)"
                                            ),
                                            guide = guide_legend(
                                                override.aes = list(fill = "gray70", color = "black", size = 3),
                                                title.position = "top",
                                                title.hjust = 0.5
                                            )
                                        ) +
                                        geom_vline(xintercept = th_global, linetype = "dashed",
                                                  color = "red", linewidth = 0.8, alpha = 0.5) +
                                        geom_hline(yintercept = threshold, linetype = "dashed",
                                                  color = "blue", linewidth = 0.8, alpha = 0.5) +
                                        # Add threshold annotations
                                        annotate("text", x = max(data$predicted_prob, na.rm = TRUE) * 0.95,
                                                y = threshold,
                                                label = paste0("Threshold = ", round(threshold, 3)),
                                                hjust = 1, vjust = -0.5, color = "blue", size = 3, fontface = "bold") +
                                        annotate("text", x = th_global,
                                                y = max(data[[y_col]], na.rm = TRUE) * 0.95,
                                                label = paste0("Youden J = ", round(th_global, 3)),
                                                hjust = -0.1, vjust = 1, color = "red", size = 3, fontface = "bold") +
                                        labs(
                                            title = paste("pQTL Detection:", title_suffix, "(Test Case)"),
                                            subtitle = paste("Threshold:", round(threshold, 3),
                                                           "| Detected outliers (", y_col, "):", sum(data[[outlier_col]], na.rm=TRUE),
                                                           "| Expected errors:", sum(data$is_error_sample),
                                                           "| Cross-Sex Mismatches:", sum(data$is_cross_sex, na.rm = TRUE)),
                                            x = "Predicted Female Probability",
                                            y = y_label
                                        ) +
                                        theme_bw() +
                                        theme(
                                            plot.title = element_text(size = 14, face = "bold"),
                                            plot.subtitle = element_text(size = 11),
                                            legend.position = "right",
                                            strip.text = element_text(size = 12, face = "bold")
                                        )

                                    plot_path <- file.path(test_output_dir,
                                        paste0("07_07c_pqtl_", file_suffix, "_test_case_fg3_batch_02.pdf"))
                                    ggsave(plot_path, p, width = 12, height = 8)
                                    log_info("Saved {file_suffix} plot to: {plot_path}")

                                    return(p)
                                }

                                # Plot 1: MeanAbsZ
                                # Note: Each plot uses its own metric-specific outlier flag for False Positive/Negative designation
                                p1 <- create_metric_plot(
                                    copy(plot_data), "MeanAbsZ", "Mean Absolute Z-score",
                                    cutoff_mean_z, "Mean Absolute Z-score", "mean_absz"
                                )

                                # Plot 2: MaxAbsZ
                                p2 <- create_metric_plot(
                                    copy(plot_data), "MaxAbsZ", "Maximum Absolute Z-score",
                                    cutoff_max_z, "Maximum Absolute Z-score", "max_absz"
                                )

                                # Plot 3: OutlierVariantFraction
                                p3 <- create_metric_plot(
                                    copy(plot_data), "OutlierVariantFraction",
                                    paste0("Fraction of Variants with |Z| > ", z_outlier_threshold),
                                    cutoff_outlier_frac, "Outlier Variant Fraction", "outlier_frac"
                                )

                                # Plot 4: MeanSquaredResidual (beta-weighted)
                                if ("MeanSquaredResidual" %in% names(plot_data)) {
                                    p4 <- create_metric_plot(
                                        copy(plot_data), "MeanSquaredResidual",
                                        "Mean Squared Residual (Beta-Weighted)",
                                        cutoff_msr, "Mean Squared Residual", "msr"
                                    )
                                }

                                # Plot 5: MedianAbsResidual (beta-weighted)
                                if ("MedianAbsResidual" %in% names(plot_data)) {
                                    p5 <- create_metric_plot(
                                        copy(plot_data), "MedianAbsResidual",
                                        "Median Absolute Residual (Beta-Weighted)",
                                        cutoff_mar, "Median Absolute Residual", "mar"
                                    )
                                }

                                # Plot 6: CompositeScore (rank-normalized combination)
                                if ("CompositeScore" %in% names(plot_data)) {
                                    p6 <- create_metric_plot(
                                        copy(plot_data), "CompositeScore",
                                        "Composite Anomaly Score (Rank-Normalized)",
                                        composite_threshold, "Composite Score", "composite"
                                    )
                                }

                                # Plot 7: Combined view (original format for backward compatibility)
                                # Enhanced to annotate cross-sex mismatches
                                # This plot shows MeanAbsZ, so use MeanAbsZ-specific outlier flag for consistency
                                plot_data[, point_shape := fifelse(Outlier_MeanAbsZ, "Outlier", "Normal")]

                                # Add cross-sex mismatch indicator
                                plot_data[, is_cross_sex := category == "Category3_CrossSex"]
                                plot_data[, point_shape_enhanced := fifelse(
                                    is_cross_sex, "Cross-Sex Mismatch",
                                    fifelse(Outlier_MeanAbsZ, "Outlier", "Normal")
                                )]

                                # Calculate MeanAbsZ-specific point_type for combined plot
                                plot_data[, point_type_meanabsz := fifelse(
                                    is_error_sample & Outlier_MeanAbsZ, "True Positive (Error Detected)",
                                    fifelse(is_error_sample & !Outlier_MeanAbsZ, "False Negative (Error Missed)",
                                    fifelse(!is_error_sample & Outlier_MeanAbsZ, "False Positive",
                                    "True Negative (Clean)"))
                                )]

                                p_combined <- ggplot(plot_data, aes(x = predicted_prob, y = MeanAbsZ)) +
                                    geom_point(aes(fill = point_type_meanabsz, shape = point_shape),
                                              size = 2.5, alpha = 0.7, stroke = 0.5, color = "black") +
                                    # Highlight cross-sex mismatches with thicker purple border
                                    geom_point(data = plot_data[is_cross_sex == TRUE],
                                              aes(shape = point_shape_enhanced),
                                              fill = NA, color = "purple",
                                              size = 4, stroke = 2, alpha = 1) +
                                    # Highlight other error samples with black border
                                    geom_point(data = plot_data[is_error_sample == TRUE & is_cross_sex == FALSE],
                                              aes(shape = point_shape),
                                              fill = NA, color = "black",
                                              size = 3.5, stroke = 1.5, alpha = 1) +
                                    # Annotate cross-sex mismatches
                                    geom_text_repel(data = plot_data[is_cross_sex == TRUE & !is.na(FINNGENID)],
                                                   aes(label = FINNGENID),
                                                   color = "purple", size = 2.5, fontface = "bold",
                                                   box.padding = 0.5, point.padding = 0.3, max.overlaps = Inf) +
                                    facet_wrap(~ genetic_sex, nrow = 2, scales = "free_x") +
                                    scale_fill_manual(
                                        values = c(
                                            "True Positive (Error Detected)" = "#2E7D32",
                                            "False Negative (Error Missed)" = "#C62828",
                                            "False Positive" = "#FF9800",
                                            "True Negative (Clean)" = "#90CAF9"
                                        ),
                                        name = "Detection Status (MeanAbsZ)",
                                        labels = c(
                                            "True Positive (Error Detected)" = "✓ Error Detected (Green)",
                                            "False Negative (Error Missed)" = "✗ Error Missed (Red)",
                                            "False Positive" = "⚠ False Positive (Orange)",
                                            "True Negative (Clean)" = "○ Clean Sample (Blue)"
                                        ),
                                        guide = guide_legend(
                                            override.aes = list(shape = 21, size = 4),
                                            title.position = "top",
                                            title.hjust = 0.5
                                        )
                                    ) +
                                    scale_shape_manual(
                                        values = c(
                                            "Outlier" = 24,              # Triangle up
                                            "Normal" = 21,               # Circle
                                            "Cross-Sex Mismatch" = 23    # Diamond
                                        ),
                                        name = "Sample Type",
                                        labels = c(
                                            "Outlier" = "Outlier (Triangle)",
                                            "Normal" = "Normal (Circle)",
                                            "Cross-Sex Mismatch" = "Cross-Sex Mismatch (Purple diamond)"
                                        ),
                                        guide = guide_legend(
                                            override.aes = list(fill = "gray70", color = "black", size = 3),
                                            title.position = "top",
                                            title.hjust = 0.5
                                        )
                                    ) +
                                    geom_vline(xintercept = th_global, linetype = "dashed",
                                              color = "red", linewidth = 0.8, alpha = 0.5) +
                                    geom_hline(yintercept = cutoff_mean_z, linetype = "dashed",
                                              color = "blue", linewidth = 0.8, alpha = 0.5) +
                                    labs(
                                        title = "pQTL Multi-Metric Detection (Test Case - Combined)",
                                        subtitle = paste(
                                            "MeanAbsZ Threshold:", round(cutoff_mean_z, 3),
                                            "| True Positives:", sum(plot_data$is_error_sample & plot_data$Outlier_MeanAbsZ),
                                            "| False Negatives:", sum(plot_data$is_error_sample & !plot_data$Outlier_MeanAbsZ),
                                            "| False Positives:", sum(!plot_data$is_error_sample & plot_data$Outlier_MeanAbsZ),
                                            "| Cross-Sex Mismatches:", sum(plot_data$is_cross_sex, na.rm = TRUE)
                                        ),
                                        x = "Predicted Female Probability",
                                        y = "Mean Absolute Z-score"
                                    ) +
                                    theme_bw() +
                                    theme(
                                        plot.title = element_text(size = 14, face = "bold"),
                                        plot.subtitle = element_text(size = 11),
                                        legend.position = "right",
                                        strip.text = element_text(size = 12, face = "bold")
                                    )

                                pqtl_plot_path <- file.path(test_output_dir,
                                    "07_07c_pqtl_sex_crossref_faceted_test_case_fg3_batch_02.pdf")
                                ggsave(pqtl_plot_path, p_combined, width = 12, height = 8)
                                log_info("Saved combined multi-metric plot to: {pqtl_plot_path}")

                                # Generate distribution plot
                                dist_plot_path <- file.path(test_output_dir,
                                    "07_07c_meanabsz_distribution_test_case_fg3_batch_02.pdf")
                                create_meanabsz_distribution_plot(sample_stats, ground_truth, dist_plot_path)

                                # Generate scatter plot (MeanAbsZ vs Genetic Distance)
                                if ("genetic_distance" %in% names(ground_truth) &&
                                    sum(!is.na(ground_truth$genetic_distance)) > 0) {
                                    scatter_plot_path <- file.path(test_output_dir,
                                        "07_07c_meanabsz_vs_distance_test_case_fg3_batch_02.pdf")
                                    create_meanabsz_vs_distance_plot(sample_stats, ground_truth, scatter_plot_path)
                                } else {
                                    log_warn("Genetic distances not available. Skipping scatter plot.")
                                }

                                # Generate scatter plot (gPC1 vs MeanAbsZ)
                                covariate_file <- tryCatch(config$covariates$covariate_file, error = function(e) NULL)
                                if (!is.null(covariate_file) && file.exists(covariate_file)) {
                                    # Load genetic PCs for scatter plot
                                    all_finngenids <- unique(ground_truth$FINNGENID)
                                    gen_pcs_plot <- load_genetic_pcs(covariate_file, all_finngenids)
                                    if (!is.null(gen_pcs_plot) && nrow(gen_pcs_plot) > 0) {
                                        gpc1_plot_path <- file.path(test_output_dir,
                                            "07_07c_gpc1_vs_meanabsz_test_case_fg3_batch_02.pdf")
                                        create_gpc1_vs_meanabsz_plot(sample_stats, ground_truth, gen_pcs_plot, gpc1_plot_path)
                                    } else {
                                        log_warn("Genetic PCs not available. Skipping gPC1 scatter plot.")
                                    }
                                } else {
                                    log_warn("Covariate file not available. Skipping gPC1 scatter plot.")
                                }

                                # ================================================================
                                # DETECTION SUMMARY
                                # ================================================================
                                log_info("=== DETECTION SUMMARY ===")

                                # Calculate detection rates by error type
                                error_samples <- plot_data[is_error_sample == TRUE]
                                genotype_only <- plot_data[has_genotype_error == TRUE & has_sex_error == FALSE]
                                both_errors <- plot_data[has_genotype_error == TRUE & has_sex_error == TRUE]

                                if (nrow(error_samples) > 0) {
                                    tp <- sum(error_samples$Outlier)
                                    fn <- sum(!error_samples$Outlier)
                                    log_info("Error sample detection: {tp}/{nrow(error_samples)} ({round(100*tp/nrow(error_samples), 1)}%)")

                                    if (nrow(genotype_only) > 0) {
                                        tp_geno <- sum(genotype_only$Outlier)
                                        log_info("  Genotype-only errors: {tp_geno}/{nrow(genotype_only)} ({round(100*tp_geno/nrow(genotype_only), 1)}%)")
                                    }

                                    if (nrow(both_errors) > 0) {
                                        tp_both <- sum(both_errors$Outlier)
                                        log_info("  Genotype+Sex errors: {tp_both}/{nrow(both_errors)} ({round(100*tp_both/nrow(both_errors), 1)}%)")
                                    }

                                    # List missed samples
                                    missed <- error_samples[Outlier == FALSE]$SampleID
                                    if (length(missed) > 0) {
                                        log_warn("Missed error samples ({length(missed)}): {paste(head(missed, 10), collapse=', ')}")
                                    }
                                }

                                clean_samples <- plot_data[is_error_sample == FALSE]
                                if (nrow(clean_samples) > 0) {
                                    fp <- sum(clean_samples$Outlier)
                                    log_info("False positives: {fp}/{nrow(clean_samples)} ({round(100*fp/nrow(clean_samples), 1)}%)")
                                }
                            }
                        } else {
                            log_warn("No z-scores calculated. Check genotype file generation.")
                        }
                    } else {
                        log_warn("PLINK output not found. Skipping pQTL detection.")
                    }
                } else {
                    log_warn("Genotype files not available. Skipping pQTL detection.")
                }
            }
        }

        # Calculate detection rates
        validation_result$sex_detection_rate <- if (validation_result$n_sex_expected > 0) {
            validation_result$n_sex_detected / validation_result$n_sex_expected
        } else { 1 }

        validation_result$pqtl_detection_rate <- if (validation_result$n_pqtl_expected > 0) {
            validation_result$n_pqtl_detected / validation_result$n_pqtl_expected
        } else { 1 }

        log_info("Validation results:")
        log_info("  Sex detection: {validation_result$n_sex_detected}/{validation_result$n_sex_expected} ({round(validation_result$sex_detection_rate * 100, 1)}%)")
        log_info("  pQTL detection: {validation_result$n_pqtl_detected}/{validation_result$n_pqtl_expected} ({round(validation_result$pqtl_detection_rate * 100, 1)}%)")

    }, error = function(e) {
        log_error("Error during validation tests: {e$message}")
        log_error("Traceback: {paste(capture.output(traceback()), collapse = '\\n')}")
    })

    return(validation_result)
}

# Function to extract actual genotypes for a cohort using sparse method
# Similar to 07b_pqtl_outliers.R but for a specific set of FINNGENIDs
extract_cohort_genotypes <- function(cohort_finngenids, top_variants, genotype_path, config, batch_id, cached_genotype_info = NULL) {
    log_info("Extracting actual genotypes for {length(cohort_finngenids)} samples using sparse method...")

    # Check for cached genotypes first
    local_genotype_path <- NULL
    temp_geno_dir <- NULL

    if (!is.null(cached_genotype_info) &&
        file.exists(cached_genotype_info$bed) &&
        file.exists(cached_genotype_info$bim) &&
        file.exists(cached_genotype_info$fam)) {
        log_info("Using cached genotype files from: {cached_genotype_info$dir}")
        local_genotype_path <- cached_genotype_info$prefix
    } else {
        # Need to extract genotypes using sparse method
        # This is a simplified version - in production, we'd use the full sparse extraction
        # For now, we'll use PLINK directly if we have a local path, or extract sparse if GS
        if (startsWith(genotype_path, "gs://")) {
            log_info("GS genotype path detected. Will use sparse extraction method.")
            log_warn("Sparse extraction for cohort not fully implemented. Using cached genotypes or skipping.")
            return(NULL)
        } else {
            # Local path - use PLINK directly
            local_genotype_path <- genotype_path
        }
    }

    # Create temp files for PLINK
    var_file <- tempfile(pattern = "variants_", fileext = ".snplist")
    out_prefix <- tempfile(pattern = "cohort_geno_")

    # Write variant list
    if (!is.null(top_variants) && nrow(top_variants) > 0) {
        # Use rsid_for_matching if available, otherwise rsid
        variant_ids <- if ("rsid_for_matching" %in% names(top_variants)) {
            top_variants$rsid_for_matching
        } else {
            top_variants$rsid
        }
        writeLines(variant_ids, var_file)
    } else {
        log_warn("No variants provided for genotype extraction")
        return(NULL)
    }

    # Create sample list file
    sample_file <- tempfile(pattern = "samples_", fileext = ".txt")
    writeLines(cohort_finngenids, sample_file)

    # Run PLINK to extract genotypes
    is_bed <- file.exists(paste0(local_genotype_path, ".bed"))
    flag <- if (is_bed) "--bfile" else "--pfile"

    plink_cmd <- sprintf(
        "plink2 %s %s --extract %s --keep %s --export A --out %s 2>&1",
        flag, local_genotype_path, var_file, sample_file, out_prefix
    )

    log_info("Running PLINK for cohort genotype extraction: {plink_cmd}")
    plink_result <- system(plink_cmd, intern = TRUE)
    plink_exit_code <- attr(plink_result, "status") %||% 0

    if (plink_exit_code != 0) {
        log_warn("PLINK extraction failed. Exit code: {plink_exit_code}")
        log_warn("PLINK output: {paste(head(plink_result, 10), collapse = '\\n')}")
        return(NULL)
    }

    # Read PLINK output
    raw_file <- paste0(out_prefix, ".raw")
    if (!file.exists(raw_file)) {
        log_warn("PLINK output file not found: {raw_file}")
        return(NULL)
    }

    dt_geno <- fread(raw_file)
    log_info("Extracted genotypes for {nrow(dt_geno)} samples and {length(variant_ids)} variants")

    # Cleanup temp files
    if (file.exists(var_file)) unlink(var_file)
    if (file.exists(sample_file)) unlink(sample_file)
    if (file.exists(out_prefix)) unlink(paste0(out_prefix, "*"))

    return(dt_geno)
}

# ============================================================================
# Category-Specific Error Introduction Functions
# ============================================================================

# Category 1: Within-Sex ID Shuffles (swap FINNGENID within same sex)
apply_category1_swaps <- function(
    error_samples,          # data.table with SAMPLE_ID, FINNGENID, genetic_sex
    metadata_synthetic,     # Full cohort metadata
    error_tracking,         # Error tracking data.table
    distance_matrix,        # Genetic distance matrix
    distance_config         # Distance configuration
) {
    log_info("--- Category 1: Within-Sex ID Shuffles ---")

    n_errors <- nrow(error_samples)
    if (n_errors == 0) return(list(metadata = metadata_synthetic, error_tracking = error_tracking))

    # Separate by sex
    male_samples <- error_samples[genetic_sex == "male"]
    female_samples <- error_samples[genetic_sex == "female"]

    # Determine percentile thresholds (mixed strategy: 60% use 75th, 40% use 90th)
    n_75 <- round(n_errors * (distance_config$percentile_75_fraction %||% 0.6))
    n_90 <- n_errors - n_75

    # Assign thresholds randomly
    thresholds <- c(rep(75, n_75), rep(90, n_90))
    thresholds <- sample(thresholds)  # Randomize order

    # Process male samples
    if (nrow(male_samples) > 0) {
        male_pool <- metadata_synthetic[genetic_sex == "male" & !FINNGENID %in% error_samples$FINNGENID]$FINNGENID

        for (i in seq_len(nrow(male_samples))) {
            sample_id <- male_samples$SAMPLE_ID[i]
            source_finngenid <- male_samples$FINNGENID[i]
            threshold_percentile <- thresholds[i]

            # Filter pool by distance
            filtered_pool <- filter_swap_pool_by_distance(
                source_finngenid, male_pool, distance_matrix,
                threshold_percentile,
                auto_lower = distance_config$auto_lower_threshold %||% TRUE,
                log_warnings = distance_config$log_warnings %||% TRUE
            )

            if (length(filtered_pool) > 0) {
                new_finngenid <- sample(filtered_pool, 1)
                idx <- which(metadata_synthetic$SAMPLE_ID == sample_id)
                if (length(idx) > 0) {
                    metadata_synthetic[idx, FINNGENID := new_finngenid]
                    error_tracking[SAMPLE_ID == sample_id, swapped_finngenid := new_finngenid]
                    error_tracking[SAMPLE_ID == sample_id, swap_type := "within_sex"]
                    error_tracking[SAMPLE_ID == sample_id, swap_distance := distance_matrix[source_finngenid, new_finngenid]]
                    error_tracking[SAMPLE_ID == sample_id, category := "Category1_WithinSex"]
                    error_tracking[SAMPLE_ID == sample_id, has_genotype_error := TRUE]
                }
            }
        }
    }

    # Process female samples (similar logic)
    if (nrow(female_samples) > 0) {
        female_pool <- metadata_synthetic[genetic_sex == "female" & !FINNGENID %in% error_samples$FINNGENID]$FINNGENID

        offset <- nrow(male_samples)
        for (i in seq_len(nrow(female_samples))) {
            sample_id <- female_samples$SAMPLE_ID[i]
            source_finngenid <- female_samples$FINNGENID[i]
            threshold_percentile <- thresholds[offset + i]

            filtered_pool <- filter_swap_pool_by_distance(
                source_finngenid, female_pool, distance_matrix,
                threshold_percentile,
                auto_lower = distance_config$auto_lower_threshold %||% TRUE,
                log_warnings = distance_config$log_warnings %||% TRUE
            )

            if (length(filtered_pool) > 0) {
                new_finngenid <- sample(filtered_pool, 1)
                idx <- which(metadata_synthetic$SAMPLE_ID == sample_id)
                if (length(idx) > 0) {
                    metadata_synthetic[idx, FINNGENID := new_finngenid]
                    error_tracking[SAMPLE_ID == sample_id, swapped_finngenid := new_finngenid]
                    error_tracking[SAMPLE_ID == sample_id, swap_type := "within_sex"]
                    error_tracking[SAMPLE_ID == sample_id, swap_distance := distance_matrix[source_finngenid, new_finngenid]]
                    error_tracking[SAMPLE_ID == sample_id, category := "Category1_WithinSex"]
                    error_tracking[SAMPLE_ID == sample_id, has_genotype_error := TRUE]
                }
            }
        }
    }

    log_info("Applied Category 1 swaps to {n_errors} samples (within-sex)")
    return(list(metadata = metadata_synthetic, error_tracking = error_tracking))
}

# Category 2: Sex Label Misannotation Only (no ID swap)
apply_category2_sex_labels <- function(
    error_samples,          # data.table with SAMPLE_ID, FINNGENID, genetic_sex
    metadata_synthetic,     # Full cohort metadata
    error_tracking         # Error tracking data.table
) {
    log_info("--- Category 2: Sex Label Misannotation Only ---")

    n_errors <- nrow(error_samples)
    if (n_errors == 0) return(list(metadata = metadata_synthetic, error_tracking = error_tracking))

    # Flip sex labels only (no FINNGENID swap)
    for (i in seq_len(n_errors)) {
        sample_id <- error_samples$SAMPLE_ID[i]
        idx <- which(metadata_synthetic$SAMPLE_ID == sample_id)

        if (length(idx) > 0) {
            # Flip sex label
            metadata_synthetic[idx, genetic_sex := ifelse(genetic_sex == "male", "female", "male")]
            error_tracking[SAMPLE_ID == sample_id, modified_sex := metadata_synthetic[idx]$genetic_sex]
            error_tracking[SAMPLE_ID == sample_id, category := "Category2_SexLabel"]
            error_tracking[SAMPLE_ID == sample_id, has_sex_error := TRUE]
            error_tracking[SAMPLE_ID == sample_id, swap_type := NA_character_]  # No swap
        }
    }

    log_info("Applied Category 2 sex label flips to {n_errors} samples")
    return(list(metadata = metadata_synthetic, error_tracking = error_tracking))
}

# Category 3: Cross-Sex ID Shuffles (swap FINNGENID across sexes)
apply_category3_swaps <- function(
    error_samples,          # data.table with SAMPLE_ID, FINNGENID, genetic_sex
    metadata_synthetic,     # Full cohort metadata
    error_tracking,         # Error tracking data.table
    distance_matrix,        # Genetic distance matrix
    distance_config         # Distance configuration
) {
    log_info("--- Category 3: Cross-Sex ID Shuffles ---")

    n_errors <- nrow(error_samples)
    if (n_errors == 0) return(list(metadata = metadata_synthetic, error_tracking = error_tracking))

    # Determine percentile thresholds (mixed strategy)
    n_75 <- round(n_errors * (distance_config$percentile_75_fraction %||% 0.6))
    n_90 <- n_errors - n_75
    thresholds <- c(rep(75, n_75), rep(90, n_90))
    thresholds <- sample(thresholds)

    for (i in seq_len(n_errors)) {
        sample_id <- error_samples$SAMPLE_ID[i]
        source_finngenid <- error_samples$FINNGENID[i]
        source_sex <- error_samples$genetic_sex[i]
        threshold_percentile <- thresholds[i]

        # Get opposite-sex pool
        opposite_sex_pool <- metadata_synthetic[genetic_sex != source_sex & !FINNGENID %in% error_samples$FINNGENID]$FINNGENID

        # Filter pool by distance
        filtered_pool <- filter_swap_pool_by_distance(
            source_finngenid, opposite_sex_pool, distance_matrix,
            threshold_percentile,
            auto_lower = distance_config$auto_lower_threshold %||% TRUE,
            log_warnings = distance_config$log_warnings %||% TRUE
        )

        if (length(filtered_pool) > 0) {
            new_finngenid <- sample(filtered_pool, 1)
            idx <- which(metadata_synthetic$SAMPLE_ID == sample_id)
            if (length(idx) > 0) {
                metadata_synthetic[idx, FINNGENID := new_finngenid]
                error_tracking[SAMPLE_ID == sample_id, swapped_finngenid := new_finngenid]
                error_tracking[SAMPLE_ID == sample_id, swap_type := "cross_sex"]
                error_tracking[SAMPLE_ID == sample_id, swap_distance := distance_matrix[source_finngenid, new_finngenid]]
                error_tracking[SAMPLE_ID == sample_id, category := "Category3_CrossSex"]
                error_tracking[SAMPLE_ID == sample_id, has_genotype_error := TRUE]
                error_tracking[SAMPLE_ID == sample_id, has_sex_error := TRUE]  # Cross-sex creates sex mismatch
            }
        }
    }

    log_info("Applied Category 3 swaps to {n_errors} samples (cross-sex)")
    return(list(metadata = metadata_synthetic, error_tracking = error_tracking))
}

main <- function() {
    log_info("=== GOVERNANCE TEST: Data Preparation ===")

    # 1. Load required data
    # --------------------
    log_info("Loading data from previous steps...")

    # Load NPX matrix (step 01 cleaned or step 00)
    npx_path_04 <- get_output_path("01", "npx_matrix_pca_cleaned", batch_id, "outliers", config = config)
    npx_path_03 <- get_output_path("00", "npx_matrix_analysis_ready", batch_id, "qc", config = config)

    if (file.exists(npx_path_04)) {
        npx_matrix <- readRDS(npx_path_04)
        log_info("Loaded NPX matrix from step 01: {nrow(npx_matrix)} samples")
    } else if (file.exists(npx_path_03)) {
        npx_matrix <- readRDS(npx_path_03)
        log_info("Loaded NPX matrix from step 00: {nrow(npx_matrix)} samples")
    } else {
        stop("No NPX matrix found. Run Steps 03-04 first.")
    }

    # Load metadata
    metadata_path <- get_output_path("00", "metadata", batch_id, "qc", config = config)
    if (!file.exists(metadata_path)) {
        stop("Metadata file not found: {metadata_path}")
    }
    metadata <- readRDS(metadata_path)
    metadata <- add_disease_group(metadata)
    log_info("Loaded metadata: {nrow(metadata)} samples")

    # Load sample mapping
    sample_mapping_path <- get_output_path("00", "sample_mapping", batch_id, "qc", "tsv", config = config)
    if (!file.exists(sample_mapping_path)) {
        stop("Sample mapping file not found: {sample_mapping_path}")
    }
    dt_mapping <- fread(sample_mapping_path)
    log_info("Loaded sample mapping: {nrow(dt_mapping)} samples")

    # Load pQTL collated file to get top variants
    pqtl_config <- config$parameters$pqtl_outliers
    if (is.null(pqtl_config)) {
        stop("No pqtl_outliers configuration found")
    }

    # Try both possible paths (07b and 07)
    collated_path <- get_output_path("05", "07.b_finemap_collated", batch_id, "pqtl", "tsv", config = config)
    if (!file.exists(collated_path) || file.size(collated_path) == 0) {
        # Try with step_num "07" instead
        collated_path_alt <- get_output_path("03", "07.b_finemap_collated", batch_id, "pqtl", "tsv", config = config)
        if (file.exists(collated_path_alt) && file.size(collated_path_alt) > 0) {
            collated_path <- collated_path_alt
            log_info("Found collated file at alternative path: {collated_path}")
        }
    }

    if (!file.exists(collated_path) || file.size(collated_path) == 0) {
        log_warn("pQTL collated file not found. Will use all available variants.")
        top_variants <- NULL
    } else {
        top_variants <- fread(collated_path)
        log_info("Loaded {nrow(top_variants)} variants from collated file")

        # Apply MAF filter if configured (for governance test)
        pqtl_selection_config <- tryCatch(test_config$pqtl_selection, error = function(e) NULL)
        min_maf_threshold <- tryCatch(
            pqtl_selection_config$min_maf %||% NULL,
            error = function(e) NULL
        )

        if (!is.null(min_maf_threshold) && is.numeric(min_maf_threshold) && min_maf_threshold > 0) {
            if ("maf" %in% names(top_variants)) {
                n_before_maf <- nrow(top_variants)
                top_variants <- top_variants[!is.na(maf) & maf > min_maf_threshold]
                n_after_maf <- nrow(top_variants)
                log_info("Applied MAF > {min_maf_threshold*100}% filter: {n_before_maf} -> {n_after_maf} variants")

                if (n_after_maf == 0) {
                    log_warn("No variants remaining after MAF filter. Using all variants.")
                    top_variants <- fread(collated_path)  # Reload original
                }
            } else {
                log_warn("MAF column not found in collated file. Cannot apply MAF filter.")
            }
        } else {
            log_info("MAF filter not configured or disabled. Using all variants.")
        }

        # Apply proper ranking (matching 07b_pqtl_outliers.R)
        # If heterozygosity column exists, use heterozygosity-weighted ranking
        # Otherwise, use MAF-weighted ranking or standard ranking
        if ("heterozygosity" %in% names(top_variants) && "maf" %in% names(top_variants)) {
            log_info("Applying heterozygosity-weighted ranking: (-log10(p) * abs(beta) * maf) / (heterozygosity + 0.01)")
            # Add epsilon to heterozygosity to avoid division by zero
            top_variants[, heterozygosity_adj := heterozygosity + 0.01]
            # Calculate composite score with heterozygosity weighting
            top_variants[, composite_score := (-log10(p) * abs(beta) * maf) / heterozygosity_adj]
            # Sort by composite score (descending)
            setorder(top_variants, -composite_score)
            log_info("Ranked variants using heterozygosity-weighted composite score")
        } else if ("maf" %in% names(top_variants)) {
            log_info("Applying MAF-weighted ranking: -log10(p) * abs(beta) * maf")
            top_variants[, composite_score := -log10(p) * abs(beta) * maf]
            setorder(top_variants, -composite_score)
            log_info("Ranked variants using MAF-weighted composite score")
        } else {
            log_info("Applying standard ranking: -log10(p) * abs(beta)")
            top_variants[, composite_score := -log10(p) * abs(beta)]
            setorder(top_variants, -composite_score)
            log_info("Ranked variants using standard composite score")
        }

        # Filter to top 100 ONLY if MAF filter is NOT configured
        # If MAF filter is configured, use all variants that pass the MAF threshold
        if (is.null(min_maf_threshold) || !is.numeric(min_maf_threshold) || min_maf_threshold <= 0) {
            # Only apply top-100 filter when MAF filter is disabled
            if (nrow(top_variants) > 100) {
                # If included_in_top_n column exists, prefer those variants but still limit to top 100 by ranking
                if ("included_in_top_n" %in% names(top_variants)) {
                    included_variants <- top_variants[included_in_top_n == TRUE]
                    if (nrow(included_variants) > 0 && nrow(included_variants) <= 100) {
                        top_variants <- included_variants
                        log_info("Using {nrow(top_variants)} variants marked as included_in_top_n")
                    } else {
                        # Use top 100 by ranking
                        top_variants <- top_variants[1:100]
                        log_info("Selected top 100 variants by composite score ranking")
                    }
                } else {
                    # Use top 100 by ranking
                    top_variants <- top_variants[1:100]
                    log_info("Selected top 100 variants by composite score ranking")
                }
            }
        } else {
            log_info("MAF filter is active. Using all {nrow(top_variants)} variants that pass MAF > {min_maf_threshold*100}% threshold")
        }

        # Apply max_pqtls limit if configured
        max_pqtls <- tryCatch(
            pqtl_selection_config$max_pqtls %||% NULL,
            error = function(e) NULL
        )
        if (!is.null(max_pqtls) && is.numeric(max_pqtls) && max_pqtls > 0 && nrow(top_variants) > max_pqtls) {
            log_info("Applying max_pqtls limit: {nrow(top_variants)} -> {max_pqtls} variants")
            top_variants <- top_variants[1:max_pqtls]
        }

        log_info("Using top {nrow(top_variants)} pQTL variants for genotype filtering")
        if ("composite_score" %in% names(top_variants)) {
            log_info("Composite score range: {round(min(top_variants$composite_score, na.rm=TRUE), 2)} - {round(max(top_variants$composite_score, na.rm=TRUE), 2)}")
        }
    }

    # Check for cached genotype files in temp_work directories
    # Look for freq_sparse_* directories and pqtl_cache directories
    base_dir <- config$output$base_dir %||% Sys.getenv("PIPELINE_OUTPUT_DIR", "output")
    temp_work_dir <- file.path(base_dir, "temp_work")

    cached_genotype_info <- NULL
    if (dir.exists(temp_work_dir)) {
        log_info("Checking for cached genotype files in temp_work directory...")

        # Get genotype path and base name (needed for both freq_sparse and pqtl_cache checks)
        genotype_path <- pqtl_config$genotype_path
        base_name <- basename(genotype_path)

        # Check for freq_sparse_* directories (from MAF extraction)
        freq_sparse_dirs <- list.dirs(temp_work_dir, full.names = TRUE, recursive = FALSE)
        freq_sparse_dirs <- freq_sparse_dirs[grepl("freq_sparse_", basename(freq_sparse_dirs))]

        if (length(freq_sparse_dirs) > 0) {
            # Use the most recent freq_sparse directory
            freq_sparse_dirs <- freq_sparse_dirs[order(file.info(freq_sparse_dirs)$mtime, decreasing = TRUE)]
            latest_freq_dir <- freq_sparse_dirs[1]

            # Check if it contains genotype files
            cached_bed <- file.path(latest_freq_dir, paste0(base_name, ".bed"))
            cached_bim <- file.path(latest_freq_dir, paste0(base_name, ".bim"))
            cached_fam <- file.path(latest_freq_dir, paste0(base_name, ".fam"))

            if (file.exists(cached_bed) && file.exists(cached_bim) && file.exists(cached_fam)) {
                cached_genotype_info <- list(
                    type = "freq_sparse",
                    dir = latest_freq_dir,
                    bed = cached_bed,
                    bim = cached_bim,
                    fam = cached_fam,
                    prefix = file.path(latest_freq_dir, base_name)
                )
                log_info("Found cached genotype files in freq_sparse directory: {latest_freq_dir}")
            }
        }

        # If no freq_sparse found, check pqtl_cache
        if (is.null(cached_genotype_info)) {
            cache_base <- file.path(temp_work_dir, "pqtl_cache", batch_id)
            if (dir.exists(cache_base)) {
                cache_dirs <- list.dirs(cache_base, full.names = TRUE, recursive = FALSE)
                if (length(cache_dirs) > 0) {
                    # Use the most recent cache directory
                    cache_dirs <- cache_dirs[order(file.info(cache_dirs)$mtime, decreasing = TRUE)]
                    latest_cache_dir <- cache_dirs[1]

                    cached_bed <- file.path(latest_cache_dir, paste0(base_name, ".bed"))
                    cached_bim <- file.path(latest_cache_dir, paste0(base_name, ".bim"))
                    cached_fam <- file.path(latest_cache_dir, paste0(base_name, ".fam"))

                    if (file.exists(cached_bed) && file.exists(cached_bim) && file.exists(cached_fam)) {
                        cached_genotype_info <- list(
                            type = "pqtl_cache",
                            dir = latest_cache_dir,
                            bed = cached_bed,
                            bim = cached_bim,
                            fam = cached_fam,
                            prefix = file.path(latest_cache_dir, base_name)
                        )
                        log_info("Found cached genotype files in pqtl_cache: {latest_cache_dir}")
                    }
                }
            }
        }

        if (is.null(cached_genotype_info)) {
            log_info("No cached genotype files found in temp_work. Will use GS paths directly.")
        }
    } else {
        log_info("temp_work directory does not exist. Will use GS paths directly.")
    }

    # 2. Select 110 samples meeting criteria with genotype stratification
    # --------------------------------------------------------------------
    log_info("=== GOVERNANCE TEST: Sample Selection and Stratification ===")

    # Filter metadata to samples in NPX matrix
    sample_ids <- rownames(npx_matrix)
    metadata_filtered <- metadata[metadata$SAMPLE_ID %in% sample_ids]

    # Exclude F64 and Chromosomal_Abnormalities
    excl_f64 <- rep(FALSE, nrow(metadata_filtered))
    if ("F64" %in% names(metadata_filtered)) {
        v <- metadata_filtered$F64
        if (is.logical(v)) {
            excl_f64 <- !is.na(v) & v
        } else if (is.numeric(v) || is.integer(v)) {
            excl_f64 <- !is.na(v) & (v != 0)
        }
    }
    if ("DISEASE_GROUP" %in% names(metadata_filtered)) {
        dg <- toupper(as.character(metadata_filtered$DISEASE_GROUP))
        excl_f64 <- excl_f64 | (!is.na(dg) & dg == "F64")
    }

    excl_chrom <- rep(FALSE, nrow(metadata_filtered))
    chrom_col <- names(metadata_filtered)[tolower(names(metadata_filtered)) == "chromosomal_abnormalities"]
    if (length(chrom_col) > 0) {
        v <- metadata_filtered[[chrom_col[1]]]
        if (is.logical(v)) {
            excl_chrom <- !is.na(v) & v
        } else if (is.numeric(v) || is.integer(v)) {
            excl_chrom <- !is.na(v) & (v != 0)
        }
    }
    if ("DISEASE_GROUP" %in% names(metadata_filtered)) {
        dg <- toupper(as.character(metadata_filtered$DISEASE_GROUP))
        excl_chrom <- excl_chrom | (!is.na(dg) & dg == "CHROMOSOMAL_ABNORMALITIES")
    }

    # Filter eligible samples
    eligible <- metadata_filtered[!excl_f64 & !excl_chrom]
    log_info("Eligible samples (excluding F64 and Chromosomal_Abnormalities): {nrow(eligible)}")

    # Filter to samples with FINNGENIDs (required for genotype matching)
    eligible <- eligible[eligible$SAMPLE_ID %in% dt_mapping$SampleID]

    # Remove FINNGENID from eligible if it exists (to avoid merge conflicts)
    if ("FINNGENID" %in% names(eligible)) {
        eligible[, FINNGENID := NULL]
    }

    # Merge with mapping to get FINNGENID
    eligible <- merge(eligible, dt_mapping[, .(SAMPLE_ID = SampleID, FINNGENID)], by = "SAMPLE_ID", all.x = FALSE)
    eligible <- eligible[!is.na(FINNGENID)]
    log_info("Eligible samples with FINNGENIDs: {nrow(eligible)}")

    # Check if we have enough samples for cohort
    if (nrow(eligible) < cohort_size) {
        log_warn("Only {nrow(eligible)} samples meet all criteria (requested cohort size: {cohort_size})")
        cohort_size <- nrow(eligible)
        log_info("Adjusting cohort size to {cohort_size}")
    }

    # Downsample cohort first
    log_info("Downsampling cohort: selecting {cohort_size} samples from {nrow(eligible)} eligible")
    cohort_sample_ids <- sample(eligible$SAMPLE_ID, cohort_size)
    cohort_eligible <- eligible[SAMPLE_ID %in% cohort_sample_ids]
    log_info("Cohort selected: {nrow(cohort_eligible)} samples")

    # Ensure we have enough samples for error introduction
    if (nrow(cohort_eligible) < total_errors) {
        log_warn("Cohort size ({nrow(cohort_eligible)}) is smaller than requested total errors ({total_errors})")
        # Proportionally reduce each category
        reduction_factor <- nrow(cohort_eligible) / total_errors
        n_category1 <- max(1, round(n_category1 * reduction_factor))
        n_category2 <- max(1, round(n_category2 * reduction_factor))
        n_category3 <- max(0, round(n_category3 * reduction_factor))
        total_errors <- n_category1 + n_category2 + n_category3
        log_info("Adjusting error sizes: Category1={n_category1}, Category2={n_category2}, Category3={n_category3}, Total={total_errors}")
    }

    # Select error samples from within cohort (exactly as specified)
    log_info("Selecting {total_errors} error samples from cohort:")
    log_info("  Category 1 (Within-Sex ID Shuffle): {n_category1} samples")
    log_info("  Category 2 (Sex Label Misannotation): {n_category2} samples")
    log_info("  Category 3 (Cross-Sex ID Shuffle): {n_category3} samples")

    all_error_sample_ids <- sample(cohort_eligible$SAMPLE_ID, total_errors)
    cat1_sample_ids <- all_error_sample_ids[1:n_category1]
    cat2_sample_ids <- all_error_sample_ids[(n_category1 + 1):(n_category1 + n_category2)]
    cat3_sample_ids <- all_error_sample_ids[(n_category1 + n_category2 + 1):total_errors]

    log_info("Error sample selection complete:")
    log_info("  Category 1 samples: {length(cat1_sample_ids)}")
    log_info("  Category 2 samples: {length(cat2_sample_ids)}")
    log_info("  Category 3 samples: {length(cat3_sample_ids)}")
    log_info("  Clean samples: {nrow(cohort_eligible) - total_errors}")

    # 3. Create synthetic dataset with three-category error structure
    # ----------------------------------------------------------------
    log_info("=== GOVERNANCE TEST: Introducing Controlled Errors ===")

    # Create synthetic dataset directory
    # Use dirname to get directory path (get_output_path returns a file path)
    test_output_file <- get_output_path(step_num, "test_data", batch_id, "test-case", config = config)
    test_output_dir <- dirname(test_output_file)
    ensure_output_dir(test_output_dir)

    # Subset NPX matrix to FULL COHORT (not just error samples)
    npx_synthetic <- npx_matrix[cohort_eligible$SAMPLE_ID, , drop = FALSE]

    # Create synthetic metadata for FULL COHORT
    metadata_synthetic <- copy(cohort_eligible)

    # Load sex information from covariates
    covariate_file <- tryCatch(config$covariates$covariate_file, error = function(e) NULL)
    sex_info <- NULL
    if (!is.null(covariate_file) && file.exists(covariate_file)) {
        log_info("Loading sex information from covariates...")
        covariates <- read_gzipped_file(covariate_file)
        if ("SEX_IMPUTED" %in% names(covariates)) {
            sex_info <- data.table(
                FINNGENID = covariates$IID,
                genetic_sex = ifelse(covariates$SEX_IMPUTED == 1, "female", "male")
            )
        } else if ("SEX" %in% names(covariates)) {
            sex_info <- data.table(
                FINNGENID = covariates$IID,
                genetic_sex = tolower(covariates$SEX)
            )
        }
    }

    # Merge sex information
    if (!is.null(sex_info)) {
        metadata_synthetic <- merge(metadata_synthetic, sex_info, by = "FINNGENID", all.x = TRUE)
    } else {
        # Fallback: use random assignment
        log_warn("Sex information not available. Using random assignment.")
        metadata_synthetic[, genetic_sex := sample(c("male", "female"), .N, replace = TRUE)]
    }

    # Initialize error tracking for ALL cohort samples
    error_tracking <- data.table(
        SAMPLE_ID = metadata_synthetic$SAMPLE_ID,
        FINNGENID = metadata_synthetic$FINNGENID,
        category = "Clean",  # Will be: "Category1_WithinSex", "Category2_SexLabel", "Category3_CrossSex", "Clean"
        has_sex_error = FALSE,
        has_genotype_error = FALSE,
        original_sex = metadata_synthetic$genetic_sex,
        modified_sex = metadata_synthetic$genetic_sex,
        original_finngenid = metadata_synthetic$FINNGENID,
        swapped_finngenid = metadata_synthetic$FINNGENID,
        # NEW COLUMNS for enhanced error tracking:
        swap_type = NA_character_,           # "within_sex", "cross_sex", or NA (for sex label only)
        swap_distance = NA_real_,             # Genetic distance of swap pair
        swap_distance_percentile = NA_real_,  # Percentile rank of swap distance
        genetic_distance = NA_real_          # Genetic distance from population center (for plotting)
    )

    # ============================================================================
    # NEW: 3-CATEGORY ERROR INTRODUCTION WITH GENETIC DISTANCE STRATIFICATION
    # ============================================================================

    # Load genetic PCs and calculate distances if enabled
    distance_matrix <- NULL
    gen_pcs <- NULL

    if (isTRUE(distance_config$enabled)) {
        log_info("Loading genetic PCs for distance calculation...")
        all_finngenids <- unique(metadata_synthetic$FINNGENID)
        gen_pcs <- load_genetic_pcs(covariate_file, all_finngenids)
        distance_matrix <- calculate_genetic_distances(gen_pcs, distance_config$distance_metric %||% "euclidean")
    } else {
        log_info("Genetic distance stratification disabled. Using random swaps.")
    }

    # Store original FINNGENIDs before any swaps
    error_tracking[SAMPLE_ID %in% c(cat1_sample_ids, cat2_sample_ids, cat3_sample_ids),
                   original_finngenid := metadata_synthetic[SAMPLE_ID %in% c(cat1_sample_ids, cat2_sample_ids, cat3_sample_ids)]$FINNGENID]

    # Category 1: Within-Sex ID Shuffles (20 mismatches)
    if (n_category1 > 0) {
        cat1_samples <- metadata_synthetic[SAMPLE_ID %in% cat1_sample_ids, .(SAMPLE_ID, FINNGENID, genetic_sex)]
        result <- apply_category1_swaps(cat1_samples, metadata_synthetic, error_tracking,
                                        distance_matrix, distance_config)
        metadata_synthetic <- result$metadata
        error_tracking <- result$error_tracking
    }

    # Category 2: Sex Label Misannotation Only (10 mismatches)
    if (n_category2 > 0) {
        cat2_samples <- metadata_synthetic[SAMPLE_ID %in% cat2_sample_ids, .(SAMPLE_ID, FINNGENID, genetic_sex)]
        result <- apply_category2_sex_labels(cat2_samples, metadata_synthetic, error_tracking)
        metadata_synthetic <- result$metadata
        error_tracking <- result$error_tracking
    }

    # Category 3: Cross-Sex ID Shuffles (5 mismatches)
    if (n_category3 > 0) {
        cat3_samples <- metadata_synthetic[SAMPLE_ID %in% cat3_sample_ids, .(SAMPLE_ID, FINNGENID, genetic_sex)]
        result <- apply_category3_swaps(cat3_samples, metadata_synthetic, error_tracking,
                                        distance_matrix, distance_config)
        metadata_synthetic <- result$metadata
        error_tracking <- result$error_tracking
    }

    # Store genetic distances for all samples (for plotting)
    if (!is.null(gen_pcs)) {
        error_tracking <- store_genetic_distances(metadata_synthetic, gen_pcs, error_tracking)
    }

    # Summary
    log_info("Error introduction complete:")
    log_info("  Category 1 (Within-Sex ID Shuffle): {sum(error_tracking$category == 'Category1_WithinSex')} samples")
    log_info("  Category 2 (Sex Label Misannotation): {sum(error_tracking$category == 'Category2_SexLabel')} samples")
    log_info("  Category 3 (Cross-Sex ID Shuffle): {sum(error_tracking$category == 'Category3_CrossSex')} samples")
    log_info("  Clean samples: {sum(error_tracking$category == 'Clean')} samples")
    log_info("  Total with sex errors: {sum(error_tracking$has_sex_error)}")
    log_info("  Total with genotype errors: {sum(error_tracking$has_genotype_error)}")
    log_info("  Total with both errors: {sum(error_tracking$has_sex_error & error_tracking$has_genotype_error)}")

    # Log swap distance statistics
    if (!is.null(distance_matrix)) {
        cat1_distances <- error_tracking[category == "Category1_WithinSex" & !is.na(swap_distance)]$swap_distance
        cat3_distances <- error_tracking[category == "Category3_CrossSex" & !is.na(swap_distance)]$swap_distance
        if (length(cat1_distances) > 0) {
            log_info("Category 1 swap distances: mean={round(mean(cat1_distances), 3)}, median={round(median(cat1_distances), 3)}")
        }
        if (length(cat3_distances) > 0) {
            log_info("Category 3 swap distances: mean={round(mean(cat3_distances), 3)}, median={round(median(cat3_distances), 3)}")
        }
    }

    # 4. Save synthetic dataset with explicit error indicators in filenames
    # --------------------------------------------------------------------
    log_info("=== GOVERNANCE TEST: Saving Synthetic Dataset ===")

    # Create synthetic dataset directory (if not already created)
    # Use dirname to get directory path (get_output_path returns a file path)
    test_output_file <- get_output_path(step_num, "test_data", batch_id, "test-case", config = config)
    test_output_dir <- dirname(test_output_file)
    ensure_output_dir(test_output_dir)

    # Save cached genotype info to test output for reference (if available)
    if (!is.null(cached_genotype_info)) {
        cached_info_path <- file.path(test_output_dir, "cached_genotype_info.txt")
        cached_info_text <- paste0(
            "Cached Genotype Files Information\n",
            "==================================\n\n",
            "Type: ", cached_genotype_info$type, "\n",
            "Directory: ", cached_genotype_info$dir, "\n",
            "BED file: ", cached_genotype_info$bed, "\n",
            "BIM file: ", cached_genotype_info$bim, "\n",
            "FAM file: ", cached_genotype_info$fam, "\n",
            "PLINK prefix: ", cached_genotype_info$prefix, "\n\n",
            "To use these cached files, update genotype_path in config to:\n",
            cached_genotype_info$prefix, "\n"
        )
        writeLines(cached_info_text, cached_info_path)
        log_info("Saved cached genotype info: {cached_info_path}")
    }

    # Determine error types for filename based on categories
    has_genotype_only <- sum(error_tracking$category == "GenotypeMismatch", na.rm = TRUE) > 0
    has_both_errors <- sum(error_tracking$category == "GenotypeAndSexMismatch", na.rm = TRUE) > 0

    # Create comprehensive error suffix
    error_suffix <- "_with_controlled_errors"
    if (has_genotype_only && has_both_errors) {
        error_suffix <- "_with_genotype_and_sex_mismatch_errors"
    } else if (has_both_errors) {
        error_suffix <- "_with_genotype_and_sex_mismatch_errors"
    } else if (has_genotype_only) {
        error_suffix <- "_with_genotype_mismatch_errors"
    }

    # Save synthetic NPX matrix (all samples, including errors)
    npx_synthetic_path <- file.path(test_output_dir, paste0("synthetic_npx_matrix", error_suffix, ".rds"))
    saveRDS(npx_synthetic, npx_synthetic_path)
    log_info("Saved synthetic NPX matrix: {npx_synthetic_path}")

    # Save synthetic metadata (all samples, including errors)
    metadata_synthetic_path <- file.path(test_output_dir, paste0("synthetic_metadata", error_suffix, ".rds"))
    saveRDS(metadata_synthetic, metadata_synthetic_path)
    log_info("Saved synthetic metadata: {metadata_synthetic_path}")

    # Save category-specific files
    if (has_genotype_only) {
        genotype_only_samples <- error_tracking[category == "GenotypeMismatch"]$SAMPLE_ID
        npx_genotype_only <- npx_synthetic[genotype_only_samples, , drop = FALSE]
        metadata_genotype_only <- metadata_synthetic[SAMPLE_ID %in% genotype_only_samples]

        npx_genotype_only_path <- file.path(test_output_dir, "synthetic_npx_matrix_genotype_mismatch_only.rds")
        metadata_genotype_only_path <- file.path(test_output_dir, "synthetic_metadata_genotype_mismatch_only.rds")

        saveRDS(npx_genotype_only, npx_genotype_only_path)
        saveRDS(metadata_genotype_only, metadata_genotype_only_path)
        log_info("Saved genotype mismatch only samples - NPX: {npx_genotype_only_path}, Metadata: {metadata_genotype_only_path}")
    }

    if (has_both_errors) {
        both_errors_samples <- error_tracking[category == "GenotypeAndSexMismatch"]$SAMPLE_ID
        npx_both_errors <- npx_synthetic[both_errors_samples, , drop = FALSE]
        metadata_both_errors <- metadata_synthetic[SAMPLE_ID %in% both_errors_samples]

        npx_both_errors_path <- file.path(test_output_dir, "synthetic_npx_matrix_genotype_and_sex_mismatch.rds")
        metadata_both_errors_path <- file.path(test_output_dir, "synthetic_metadata_genotype_and_sex_mismatch.rds")

        saveRDS(npx_both_errors, npx_both_errors_path)
        saveRDS(metadata_both_errors, metadata_both_errors_path)
        log_info("Saved genotype + sex mismatch samples - NPX: {npx_both_errors_path}, Metadata: {metadata_both_errors_path}")
    }

    # Save separate files for samples with specific error types
    has_sex_errors <- has_both_errors
    has_genotype_errors <- has_genotype_only || has_both_errors

    if (has_sex_errors) {
        sex_error_samples <- error_tracking[has_sex_error == TRUE]$SAMPLE_ID
        npx_sex_errors <- npx_synthetic[sex_error_samples, , drop = FALSE]
        metadata_sex_errors <- metadata_synthetic[SAMPLE_ID %in% sex_error_samples]

        npx_sex_path <- file.path(test_output_dir, "synthetic_npx_matrix_with_sex_errors_only.rds")
        metadata_sex_path <- file.path(test_output_dir, "synthetic_metadata_with_sex_errors_only.rds")

        saveRDS(npx_sex_errors, npx_sex_path)
        saveRDS(metadata_sex_errors, metadata_sex_path)
        log_info("Saved sex error samples - NPX: {npx_sex_path}, Metadata: {metadata_sex_path}")
    }

    if (has_genotype_errors) {
        genotype_error_samples <- error_tracking[has_genotype_error == TRUE]$SAMPLE_ID
        npx_genotype_errors <- npx_synthetic[genotype_error_samples, , drop = FALSE]
        metadata_genotype_errors <- metadata_synthetic[SAMPLE_ID %in% genotype_error_samples]

        npx_genotype_path <- file.path(test_output_dir, "synthetic_npx_matrix_with_genotype_errors_only.rds")
        metadata_genotype_path <- file.path(test_output_dir, "synthetic_metadata_with_genotype_errors_only.rds")

        saveRDS(npx_genotype_errors, npx_genotype_path)
        saveRDS(metadata_genotype_errors, metadata_genotype_path)
        log_info("Saved genotype error samples - NPX: {npx_genotype_path}, Metadata: {metadata_genotype_path}")
    }

    # Save clean samples (no errors) for reference
    clean_samples <- error_tracking[has_sex_error == FALSE & has_genotype_error == FALSE]$SAMPLE_ID
    if (length(clean_samples) > 0) {
        npx_clean <- npx_synthetic[clean_samples, , drop = FALSE]
        metadata_clean <- metadata_synthetic[SAMPLE_ID %in% clean_samples]

        npx_clean_path <- file.path(test_output_dir, "synthetic_npx_matrix_no_errors.rds")
        metadata_clean_path <- file.path(test_output_dir, "synthetic_metadata_no_errors.rds")

        saveRDS(npx_clean, npx_clean_path)
        saveRDS(metadata_clean, metadata_clean_path)
        log_info("Saved clean samples (no errors) - NPX: {npx_clean_path}, Metadata: {metadata_clean_path}")
    }

    # 5. Create Ground Truth Reference File
    # ---------------------------------------
    log_info("=== GOVERNANCE TEST: Creating Ground Truth Reference ===")

    # Create ground truth with expected detection flags
    ground_truth <- copy(error_tracking)
    ground_truth[, expect_sex_flag := has_sex_error]
    ground_truth[, expect_pqtl_flag := has_genotype_error]

    # Add expected detection summary (updated for 3 categories)
    ground_truth_summary <- data.table(
        category = c("Category1_WithinSex", "Category2_SexLabel", "Category3_CrossSex", "Clean", "Total"),
        n_samples = c(
            sum(ground_truth$category == "Category1_WithinSex", na.rm = TRUE),
            sum(ground_truth$category == "Category2_SexLabel", na.rm = TRUE),
            sum(ground_truth$category == "Category3_CrossSex", na.rm = TRUE),
            sum(ground_truth$category == "Clean", na.rm = TRUE),
            nrow(ground_truth)
        ),
        expect_sex_flag = c(FALSE, TRUE, TRUE, FALSE, NA),  # Category 2 & 3 have sex errors
        expect_pqtl_flag = c(TRUE, FALSE, TRUE, FALSE, NA),  # Category 1 & 3 have genotype errors
        description = c(
            "Within-sex FINNGENID shuffle (genotype mismatch only)",
            "Sex label misannotation only (no ID swap)",
            "Cross-sex FINNGENID shuffle (genotype + sex mismatch)",
            "Clean samples (no errors)",
            "All categories"
        )
    )

    # Save ground truth
    ground_truth_path <- file.path(test_output_dir, "ground_truth.tsv")
    fwrite(ground_truth, ground_truth_path, sep = "\t")
    log_info("Saved ground truth reference: {ground_truth_path}")

    # Save ground truth summary
    ground_truth_summary_path <- file.path(test_output_dir, "ground_truth_summary.tsv")
    fwrite(ground_truth_summary, ground_truth_summary_path, sep = "\t")
    log_info("Saved ground truth summary: {ground_truth_summary_path}")

    # Save error tracking
    error_tracking_path <- file.path(test_output_dir, "error_tracking.tsv")
    fwrite(error_tracking, error_tracking_path, sep = "\t")
    log_info("Saved error tracking: {error_tracking_path}")

    # Save comprehensive summary
    summary_path <- file.path(test_output_dir, "test_case_summary.txt")
    n_category1_actual <- sum(error_tracking$category == "Category1_WithinSex", na.rm = TRUE)
    n_category2_actual <- sum(error_tracking$category == "Category2_SexLabel", na.rm = TRUE)
    n_category3_actual <- sum(error_tracking$category == "Category3_CrossSex", na.rm = TRUE)
    n_clean_actual <- sum(error_tracking$category == "Clean", na.rm = TRUE)

    # Calculate totals for compatibility
    n_genotype_errors <- n_category1_actual + n_category3_actual  # Category 1 & 3 have genotype errors
    n_sex_errors <- n_category2_actual + n_category3_actual  # Category 2 & 3 have sex errors
    n_genotype_only_actual <- n_category1_actual  # Category 1: genotype only
    n_both_errors_actual <- n_category3_actual  # Category 3: both genotype and sex

    summary_text <- paste0(
        "Governance Test Case Summary\n",
        "===========================\n\n",
        "Test Configuration:\n",
        "  Random seed: ", test_seed, "\n",
        "  Cohort size: ", cohort_size, " samples\n",
        "  Category 1 (Within-Sex ID Shuffle): ", n_category1, " requested, ", n_category1_actual, " created\n",
        "  Category 2 (Sex Label Misannotation): ", n_category2, " requested, ", n_category2_actual, " created\n",
        "  Category 3 (Cross-Sex ID Shuffle): ", n_category3, " requested, ", n_category3_actual, " created\n",
        "  Clean samples: ", n_clean_actual, "\n",
        "  Total samples: ", nrow(cohort_eligible), "\n\n",
        "Error Distribution:\n",
        "  Category 1 (genotype mismatch only): ", n_category1_actual, "\n",
        "  Category 2 (sex label error only): ", n_category2_actual, "\n",
        "  Category 3 (genotype + sex mismatch): ", n_category3_actual, "\n",
        "  Samples with sex errors: ", sum(error_tracking$has_sex_error), "\n",
        "  Samples with genotype errors: ", sum(error_tracking$has_genotype_error), "\n",
        "  Samples with both errors: ", sum(error_tracking$has_sex_error & error_tracking$has_genotype_error), "\n",
        "  Samples with any error: ", sum(error_tracking$has_sex_error | error_tracking$has_genotype_error), "\n\n",
        "Expected Detection (Ground Truth):\n",
        "  - Sex mismatches: ", sum(error_tracking$has_sex_error), " samples (should be detected by step 04)\n",
        "    * Category 2: ", n_category2_actual, " samples (sex label only)\n",
        "    * Category 3: ", n_category3_actual, " samples (cross-sex swaps)\n",
        "  - pQTL outliers: ", sum(error_tracking$has_genotype_error), " samples (should be detected by step 05b)\n",
        "    * Category 1: ", n_category1_actual, " samples (within-sex swaps)\n",
        "    * Category 3: ", n_category3_actual, " samples (cross-sex swaps)\n",
        "  - Combined errors: ", sum(error_tracking$has_sex_error & error_tracking$has_genotype_error), " samples\n",
        "    * Should be flagged by BOTH step 04 and step 05b\n\n",
        "Success Criteria:\n",
        "  ✓ ", n_category1_actual, " samples flagged for genotype mismatch only (Category 1)\n",
        "  ✓ ", n_category2_actual, " samples flagged for sex label error only (Category 2)\n",
        "  ✓ ", n_category3_actual, " samples flagged for genotype mismatch AND sex mismatch (Category 3)\n",
        "  ✓ ", sum(error_tracking$has_sex_error), " total sex outliers detected\n",
        "  ✓ ", sum(error_tracking$has_genotype_error), " total pQTL outliers detected\n"
    )
    writeLines(summary_text, summary_path)
    log_info("Saved test case summary: {summary_path}")

    # Save sample selection log
    selection_log_path <- file.path(test_output_dir, "sample_selection_log.txt")
    selection_log_text <- paste0(
        "Sample Selection Log\n",
        "===================\n\n",
        "Selection Criteria:\n",
        "  - Excluded F64 (Gender identity disorders)\n",
        "  - Excluded Chromosomal_Abnormalities\n",
        "  - Required FINNGENIDs for genotype matching\n",
        "  - Downsampled cohort: ", cohort_size, " samples\n\n",
        "Selected Samples:\n",
        "  Cohort size: ", cohort_size, " samples\n",
        "  Category 1 (Within-Sex ID Shuffle): ", length(cat1_sample_ids), " samples\n",
        "    ", paste(cat1_sample_ids, collapse = ", "), "\n\n",
        "  Category 2 (Sex Label Misannotation): ", length(cat2_sample_ids), " samples\n",
        "    ", paste(cat2_sample_ids, collapse = ", "), "\n\n",
        "  Category 3 (Cross-Sex ID Shuffle): ", length(cat3_sample_ids), " samples\n",
        "    ", paste(cat3_sample_ids, collapse = ", "), "\n\n",
        "  Clean samples: ", nrow(cohort_eligible) - total_errors, " samples\n"
    )
    writeLines(selection_log_text, selection_log_path)
    log_info("Saved sample selection log: {selection_log_path}")

    # Save error introduction log
    error_log_path <- file.path(test_output_dir, "error_introduction_log.txt")
    error_log_lines <- c(
        "Error Introduction Log",
        "=====================",
        "",
        "Category 1: Within-Sex ID Shuffles (genotype mismatch only)",
        paste(rep("-", 50), collapse = ""),
        ""
    )

    cat1_log <- error_tracking[category == "Category1_WithinSex"]
    for (i in seq_len(nrow(cat1_log))) {
        error_log_lines <- c(error_log_lines,
            paste0("Sample: ", cat1_log[i]$SAMPLE_ID),
            paste0("  Original FINNGENID: ", cat1_log[i]$original_finngenid, " → Swapped FINNGENID: ", cat1_log[i]$swapped_finngenid),
            paste0("  Sex: ", cat1_log[i]$original_sex, " (unchanged, within-sex swap)"),
            if (!is.na(cat1_log[i]$swap_distance)) paste0("  Genetic distance: ", round(cat1_log[i]$swap_distance, 4)),
            ""
        )
    }

    error_log_lines <- c(error_log_lines,
        "Category 2: Sex Label Misannotation Only (no ID swap)",
        paste(rep("-", 50), collapse = ""),
        ""
    )

    cat2_log <- error_tracking[category == "Category2_SexLabel"]
    for (i in seq_len(nrow(cat2_log))) {
        error_log_lines <- c(error_log_lines,
            paste0("Sample: ", cat2_log[i]$SAMPLE_ID),
            paste0("  FINNGENID: ", cat2_log[i]$FINNGENID, " (unchanged)"),
            paste0("  Original sex: ", cat2_log[i]$original_sex, " → Modified sex: ", cat2_log[i]$modified_sex),
            ""
        )
    }

    error_log_lines <- c(error_log_lines,
        "Category 3: Cross-Sex ID Shuffles (genotype + sex mismatch)",
        paste(rep("-", 50), collapse = ""),
        ""
    )

    cat3_log <- error_tracking[category == "Category3_CrossSex"]
    for (i in seq_len(nrow(cat3_log))) {
        error_log_lines <- c(error_log_lines,
            paste0("Sample: ", cat3_log[i]$SAMPLE_ID),
            paste0("  Original FINNGENID: ", cat3_log[i]$original_finngenid, " → Swapped FINNGENID: ", cat3_log[i]$swapped_finngenid),
            paste0("  Original sex: ", cat3_log[i]$original_sex, " → Modified sex: ", cat3_log[i]$modified_sex),
            if (!is.na(cat3_log[i]$swap_distance)) paste0("  Genetic distance: ", round(cat3_log[i]$swap_distance, 4)),
            ""
        )
    }

    writeLines(error_log_lines, error_log_path)
    log_info("Saved error introduction log: {error_log_path}")

    log_info("=== GOVERNANCE TEST: Synthetic Dataset Created ===")
    log_info("Test dataset saved to: {test_output_dir}")
    log_info("")
    log_info("Generated Files:")
    log_info("  - Synthetic NPX matrix: synthetic_npx_matrix{error_suffix}.rds")
    log_info("  - Synthetic metadata: synthetic_metadata{error_suffix}.rds")
    log_info("  - Ground truth: ground_truth.tsv")
    log_info("  - Error tracking: error_tracking.tsv")
    log_info("  - Category-specific files: synthetic_*_category*.rds")
    log_info("")

    # Run validation tests on synthetic dataset
    log_info("=== GOVERNANCE TEST: Running Validation Tests ===")
    validation_result <- run_validation_tests(
        synthetic_npx_path = npx_synthetic_path,
        synthetic_metadata_path = metadata_synthetic_path,
        ground_truth = ground_truth,
        batch_id = batch_id,
        config = config,
        test_output_dir = test_output_dir
    )

    log_info("=== GOVERNANCE TEST: Validation Complete ===")
    log_info("Validation results:")
    log_info("  - Sex outliers detected: {validation_result$n_sex_detected} / {validation_result$n_sex_expected}")
    log_info("  - pQTL outliers detected: {validation_result$n_pqtl_detected} / {validation_result$n_pqtl_expected}")
    log_info("  - Sex detection rate: {round(validation_result$sex_detection_rate * 100, 1)}%")
    log_info("  - pQTL detection rate: {round(validation_result$pqtl_detection_rate * 100, 1)}%")

    if (validation_result$sex_detection_rate >= 0.95 && validation_result$pqtl_detection_rate >= 0.95) {
        log_info("✓ Validation PASSED: Both detection rates >= 95%")
    } else {
        log_warn("⚠ Validation WARNING: Some detection rates < 95%")
    }

    return(list(
        test_output_dir = test_output_dir,
        n_samples = nrow(cohort_eligible),
        n_genotype_only = n_genotype_only_actual,
        n_both_errors = n_both_errors_actual,
        n_clean = n_clean_actual,
        n_sex_errors = sum(error_tracking$has_sex_error),
        n_genotype_errors = sum(error_tracking$has_genotype_error),
        error_tracking = error_tracking,
        ground_truth = ground_truth,
        validation_result = validation_result
    ))
}

# Run if executed directly
if (!interactive()) {
    result <- tryCatch({
        main()
    }, error = function(e) {
        log_error("Test case generation failed: {e$message}")
        traceback()
        quit(status = 1)
    })

    log_info("Governance test case generation completed successfully")
    log_info("Test case directory: {result$test_output_dir}")
}
