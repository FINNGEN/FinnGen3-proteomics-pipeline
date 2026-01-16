#!/usr/bin/env Rscript

#################################################
# Script: 04_sex_outliers.R
# Author: Reza Jabal, PhD (rjabal@broadinstitute.org)
# Description: Sex mismatch detection using protein expression and genetic data
#              Refactored version
# Date: December 2025
#################################################

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggrepel)
  library(yaml)
  library(logger)
  library(paletteer)
  library(ggthemes)
  library(ggpubr)
  library(dichromat)
  library(glmnet)
  library(pROC)
  library(PRROC)
  library(pheatmap)
})
## Optional model libraries (loaded lazily if available)
has_rf <- requireNamespace("randomForest", quietly = TRUE)
has_xgb <- requireNamespace("xgboost", quietly = TRUE)
has_caret <- requireNamespace("caret", quietly = TRUE)
has_keras <- requireNamespace("keras3", quietly = TRUE) || requireNamespace("keras", quietly = TRUE)
# This block suppresses "no visible binding for global variable" warnings for the specified variable names,
# which can occur in R CMD check or in non-standard evaluation contexts, such as within data.table or dplyr code.
# It tells R that these variable names are intentionally used globally within this script or package.
utils::globalVariables(
  c(
    "IID", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10",
    "SAMPLE_ID", "FINNGENID",
    paste0("pPC", 1:10),
    "SEX_IMPUTED", "AGE_AT_DEATH_OR_END_OF_FOLLOWUP", "BMI_IRN", "harmonized_current_smoker",
    ".", "p_adj", "p_value", "neg_log10_p", "significant", "label", "protein", "genetic_sex",
    "predicted_sex", "mismatch", "genetic_sex_label", "mismatch_label", "predicted_prob",
    "genetic_sex_check", "confident_mismatch", "confidence", "high_conf_mismatch",
    "bonf_sig", "OR", "CI_low", "CI_high", "dir", "idx", "color", "mag", "pt_size",
    "beta_low", "beta_high", "enet_beta_mean", "enet_beta_low", "enet_beta_high", "value"
  )
)

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
log_info("Starting sex outlier detection (nested-CV elastic-net with residualization) for batch: {batch_id}")

# Set theme for plots
theme_set(theme_bw())

# Choose priority NPX cleaned matrix from previous steps
# Prefer PCA-cleaned NPX for modeling
pick_clean_matrix <- function() {
  # Try PCA-cleaned first (Step 01), then technical-cleaned (Step 02), then zscore-cleaned (Step 03)
  candidates <- c(
    get_output_path("01", "npx_matrix_pca_cleaned", batch_id, "outliers", config = config),
    get_output_path("02", "npx_matrix_technical_cleaned", batch_id, "outliers", config = config),
    get_output_path("03", "npx_matrix_zscore_cleaned", batch_id, "outliers", config = config),
    get_output_path("00", "npx_matrix_analysis_ready", batch_id, "qc", config = config)
  )
  found <- candidates[file.exists(candidates)][1]
  if (is.na(found)) stop("No cleaned NPX matrix found. Run previous QC steps first.")
  log_info("Using matrix: {found}")
  return(found)
}

# Matrix residualization by confounders (train/test safe)
residualize_by_confounders <- function(Y, X) {
  # Y: n x p matrix (e.g., protein abundances), X: n x k design matrix (covariates including intercept) #linting: ignore
  # Returns: list with components:
  #   residuals: n x p matrix, residualized Y values
  #   coef:      k x p matrix, estimated coefficients for each covariate/protein

  # Compute (X'X)^(-1)
  XtX <- crossprod(X)
  XtX_inv <- tryCatch(solve(XtX), error = function(e) MASS::ginv(XtX)) # Fallback to generalized inverse for rank-deficiency

  # Compute coefficients B = (X'X)^(-1) X'Y
  B <- XtX_inv %*% crossprod(X, Y)  # each column of B gives betas for one protein

  # Compute residuals: Y - XB
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

# Prepare stratified K-folds
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

# Threshold sweep to maximize Youden's J and F1
pick_threshold <- function(probs, y_true) {
  grid <- unique(sort(c(seq(0.05, 0.95, by = 0.01), probs)))
  best <- data.table(th = NA_real_, J = -Inf, F1 = -Inf, sens = NA_real_, spec = NA_real_, acc = NA_real_)
  for (t in grid) {
    pred <- as.integer(probs >= t)
    TP <- sum(pred == 1 & y_true == 1, na.rm = TRUE)
    FP <- sum(pred == 1 & y_true == 0, na.rm = TRUE)
    TN <- sum(pred == 0 & y_true == 0, na.rm = TRUE)
    FN <- sum(pred == 0 & y_true == 1, na.rm = TRUE)
    sens <- ifelse((TP+FN) > 0, TP/(TP+FN), NA_real_)
    spec <- ifelse((TN+FP) > 0, TN/(TN+FP), NA_real_)
    prec <- ifelse((TP+FP) > 0, TP/(TP+FP), NA_real_)
    rec  <- sens
    F1 <- ifelse((prec+rec) > 0, 2*prec*rec/(prec+rec), NA_real_)
    J <- sens + spec - 1
    acc <- (TP+TN)/max(1, TP+TN+FP+FN)
    if (!is.na(J) && (J > best$J || (J == best$J && F1 > best$F1))) {
      best <- data.table(th = t, J = J, F1 = F1, sens = sens, spec = spec, acc = acc)
    }
  }
  best
}

# Function to load covariates with sex information and calculate sample age at collection
load_covariates <- function(covariate_file, metadata = NULL, finngen_r13_minimum_file = NULL) {
  log_info("Loading sex, age, BMI, and smoking from: {covariate_file}")

  # Load covariates (contains SEX, SEX_IMPUTED, APPROX_BIRTH_DATE, BL_AGE, BMI, smoking, etc.)
  covariates <- fread(cmd = paste("zcat", covariate_file))

  # Extract sex information: prefer SEX_IMPUTED (female=1, male=0), fallback to SEX (string)
  sex_col_used <- NA_character_
  if ("SEX_IMPUTED" %in% names(covariates)) {
    # SEX_IMPUTED: female=1, male=0
    sex_raw <- covariates$SEX_IMPUTED
    sex_col_used <- "SEX_IMPUTED"
    log_info("Using SEX_IMPUTED column from covariate file (female=1, male=0)")
  } else if ("SEX" %in% names(covariates)) {
    # SEX: "female" or "male" (string)
    sex_raw <- covariates$SEX
    sex_col_used <- "SEX"
    log_info("Using SEX column from covariate file (string format)")
  } else {
    log_warn("Neither SEX_IMPUTED nor SEX found in covariate file; sex will be NA")
    sex_raw <- NA_character_
  }

  # Extract birth date from covariate file
  birth_date_col <- if("APPROX_BIRTH_DATE" %in% names(covariates)) "APPROX_BIRTH_DATE" else NA_character_

  # Extract age from covariate file (BL_AGE or AGE_AT_DEATH_OR_END_OF_FOLLOWUP)
  age_col <- if("BL_AGE" %in% names(covariates)) "BL_AGE" else if("AGE_AT_DEATH_OR_END_OF_FOLLOWUP" %in% names(covariates)) "AGE_AT_DEATH_OR_END_OF_FOLLOWUP" else NA_character_

  # Use specific covariates as requested
  bmi_col <- if ("BMI_IRN" %in% names(covariates)) "BMI_IRN" else if ("BMI" %in% names(covariates)) "BMI" else NA_character_
  smoke_col <- if ("harmonized_current_smoker" %in% names(covariates)) "harmonized_current_smoker" else if ("CURRENT_SMOKER" %in% names(covariates)) "CURRENT_SMOKER" else NA_character_

  # Create base sex_info data.table with FINNGENID from covariate file
  sex_info <- data.table(
    FINNGENID = covariates$IID,  # IID column contains FINNGENID
    genetic_sex_raw = sex_raw,
    birth_date = if(!is.na(birth_date_col)) {
      as.POSIXct(covariates[[birth_date_col]], format="%Y-%m-%d")
    } else {
      NA
    }
  )

  # Calculate sample collection age if metadata is provided
  if (!is.null(metadata) && "APPROX_TIMESTAMP_COLLECTION" %in% names(metadata) && !all(is.na(sex_info$birth_date))) {
    # Merge with metadata to get collection timestamps (keep all SAMPLE_IDs)
    id_map <- metadata[, .(SAMPLE_ID, FINNGENID, APPROX_TIMESTAMP_COLLECTION)]

    # Log how many samples in metadata vs covariate file
    n_metadata_samples <- nrow(id_map)
    n_unique_finngenids_metadata <- length(unique(id_map$FINNGENID[!is.na(id_map$FINNGENID)]))
    n_covariate_finngenids <- length(unique(sex_info$FINNGENID))
    n_metadata_finngenids_in_covariate <- sum(unique(id_map$FINNGENID[!is.na(id_map$FINNGENID)]) %in% sex_info$FINNGENID)

    log_info("Merging sex information: {n_metadata_samples} samples in metadata, {n_unique_finngenids_metadata} unique FINNGENIDs")
    log_info("Covariate file contains {n_covariate_finngenids} FINNGENIDs; {n_metadata_finngenids_in_covariate} metadata FINNGENIDs found in covariate file")

    # Identify missing FINNGENIDs (not in covariate file)
    missing_finngenids <- setdiff(unique(id_map$FINNGENID[!is.na(id_map$FINNGENID)]), sex_info$FINNGENID)

    # Try to supplement missing FINNGENIDs from finngen_R13_minimum_1.0.txt.gz if provided
    if (length(missing_finngenids) > 0 && !is.null(finngen_r13_minimum_file) && file.exists(finngen_r13_minimum_file)) {
      log_info("Attempting to supplement {length(missing_finngenids)} missing FINNGENIDs from finngen_R13_minimum file")

      tryCatch({
        # Load finngen_R13_minimum file
        r13_minimum <- fread(cmd = paste("zcat", finngen_r13_minimum_file))

        # Filter to only "FG" prefix FINNGENIDs
        r13_minimum <- r13_minimum[grepl("^FG", FINNGENID)]

        # Extract sex information for missing FINNGENIDs
        missing_finngenids_filtered <- missing_finngenids[grepl("^FG", missing_finngenids)]
        r13_subset <- r13_minimum[FINNGENID %in% missing_finngenids_filtered]

        if (nrow(r13_subset) > 0 && "SEX" %in% names(r13_subset)) {
          # Create supplemental sex_info for missing FINNGENIDs
          # Store raw sex for later conversion (consistent with main sex_info structure)
          supplemental_sex_info <- data.table(
            FINNGENID = r13_subset$FINNGENID,
            genetic_sex_raw = r13_subset$SEX,
            birth_date = if("APPROX_BIRTH_DATE" %in% names(r13_subset)) {
              as.POSIXct(r13_subset$APPROX_BIRTH_DATE, format="%Y-%m-%d")
            } else {
              NA
            }
          )

          # Add age if available
          if ("BL_AGE" %in% names(r13_subset)) {
            supplemental_sex_info[, age := r13_subset$BL_AGE]
          } else {
            supplemental_sex_info[, age := NA_real_]
          }

          # Add BMI and smoking as NA (not available in R13 minimum file)
          supplemental_sex_info[, bmi := NA_real_]
          supplemental_sex_info[, smoking := NA_real_]

          # Append to main sex_info (before conversion, so genetic_sex_raw is preserved)
          sex_info <- rbind(sex_info, supplemental_sex_info, fill = TRUE)

          n_supplemented <- nrow(supplemental_sex_info)
          log_info("Supplemented {n_supplemented} FINNGENIDs with sex information from finngen_R13_minimum file")

          # Update missing list
          missing_finngenids <- setdiff(missing_finngenids, supplemental_sex_info$FINNGENID)

          # Update sex_col_used to indicate we have both SEX_IMPUTED and SEX (from R13 minimum)
          if (sex_col_used == "SEX_IMPUTED") {
            # Main data uses SEX_IMPUTED, supplemental uses SEX - we'll handle both in conversion
            sex_col_used <- "MIXED"
          }
        }
      }, error = function(e) {
        log_warn("Failed to load finngen_R13_minimum file: {e$message}")
      })
    }

    if (length(missing_finngenids) > 0) {
      log_warn("{length(missing_finngenids)} FINNGENIDs from metadata not found in covariate file or finngen_R13_minimum file")
      log_warn("These samples will have NA genetic_sex and will be excluded from sex-associated protein analysis")

      # Print the missing FINNGENIDs and their corresponding SAMPLE_IDs
      missing_samples <- id_map[FINNGENID %in% missing_finngenids, .(SAMPLE_ID, FINNGENID)]
      if (nrow(missing_samples) > 0) {
        log_warn("Missing FINNGENIDs and their SAMPLE_IDs:")
        for (i in seq_len(nrow(missing_samples))) {
          log_warn("  SAMPLE_ID: {missing_samples$SAMPLE_ID[i]}, FINNGENID: {missing_samples$FINNGENID[i]}")
        }
      } else {
        log_warn("Missing FINNGENIDs: {paste(missing_finngenids, collapse=', ')}")
      }
    }

    sex_info <- merge(id_map, sex_info, by = "FINNGENID", all.x = TRUE, allow.cartesian = TRUE)

    # Calculate age at sample collection for each SAMPLE_ID
    sex_info[!is.na(APPROX_TIMESTAMP_COLLECTION) & !is.na(birth_date),
             sample_age := as.numeric(difftime(APPROX_TIMESTAMP_COLLECTION, birth_date, units = "days") / 365.25)]

    # Use sample age as primary, fallback to age from covariate file
    if (!is.na(age_col) && age_col %in% names(covariates)) {
      sex_info[is.na(sample_age), sample_age := covariates[[age_col]][match(FINNGENID, covariates$IID)]]
    }

    # Keep all SAMPLE_IDs with their calculated ages
    sex_info <- sex_info[, .(SAMPLE_ID, FINNGENID, genetic_sex_raw, age = sample_age)]
  } else {
    # Fallback to age from covariate file if metadata not provided
    if (!is.na(age_col) && age_col %in% names(covariates)) {
      sex_info[, age := covariates[[age_col]][match(FINNGENID, covariates$IID)]]
    } else {
      sex_info[, age := NA_real_]
    }
  }

  # Add BMI from covariate file
  if (!is.na(bmi_col) && bmi_col %in% names(covariates)) {
    sex_info[, bmi := covariates[[bmi_col]][match(FINNGENID, covariates$IID)]]
  } else {
    sex_info[, bmi := NA_real_]
  }

  # Add smoking from covariate file
  if (!is.na(smoke_col) && smoke_col %in% names(covariates)) {
    sex_info[, smoking := covariates[[smoke_col]][match(FINNGENID, covariates$IID)]]
  } else {
    sex_info[, smoking := NA_real_]
  }

  # Convert sex to PLINK coding: 1=male, 2=female
  # SEX_IMPUTED: female=1, male=0 -> convert to PLINK: 1=male, 2=female
  # SEX: "male"/"female" strings -> convert to PLINK: 1=male, 2=female
  # Handle mixed case (SEX_IMPUTED from covariate file, SEX from R13 minimum file)
  if (sex_col_used == "SEX_IMPUTED" || sex_col_used == "MIXED") {
    # For SEX_IMPUTED: 0=male, 1=female -> PLINK: 1=male, 2=female
    # For SEX strings: "male"/"female" -> PLINK: 1=male, 2=female
    # Use vectorized & instead of && for element-wise operations
    sex_info[, genetic_sex := ifelse(
      # Check if numeric (SEX_IMPUTED: 0 or 1) OR character matching "^[01]$"
      is.numeric(genetic_sex_raw) | (is.character(genetic_sex_raw) & grepl("^[01]$", genetic_sex_raw)),
      # SEX_IMPUTED conversion
      ifelse(genetic_sex_raw == 0 | genetic_sex_raw == "0", 1,
             ifelse(genetic_sex_raw == 1 | genetic_sex_raw == "1", 2, NA_real_)),
      # SEX string conversion
      ifelse(genetic_sex_raw == "male", 1,
             ifelse(genetic_sex_raw == "female", 2, NA_real_))
    )]
  } else if (sex_col_used == "SEX") {
    # SEX: "male"/"female" strings -> PLINK: 1=male, 2=female
    sex_info[, genetic_sex := ifelse(genetic_sex_raw == "male", 1,
                                     ifelse(genetic_sex_raw == "female", 2, NA_real_))]
  } else {
    sex_info[, genetic_sex := NA_real_]
  }

  # Remove temporary genetic_sex_raw column
  if ("genetic_sex_raw" %in% names(sex_info)) {
    sex_info[, genetic_sex_raw := NULL]
  }

  log_info("Loaded sex information for {nrow(sex_info)} individuals")
  log_info("Sex distribution - Male (1): {sum(sex_info$genetic_sex == 1, na.rm=TRUE)}, Female (2): {sum(sex_info$genetic_sex == 2, na.rm=TRUE)}")

  return(sex_info)
}

# Load top proteomic PCs from 04_pca_outliers (scores matrix)
load_proteomic_pcs <- function(pca_result_file, n_pcs = 10) {
  if (is.null(pca_result_file) || !file.exists(pca_result_file)) {
    log_warn("PCA result file not found: {if(is.null(pca_result_file)) 'NULL' else pca_result_file}")
    log_warn("Proteomic PCs will not be included as covariates (step 01 may have failed)")
    # Return empty data.table with expected structure (using SampleID to match merge expectations)
    return(data.table(SampleID = character(), pPC1 = numeric(), pPC2 = numeric(), pPC3 = numeric(),
                      pPC4 = numeric(), pPC5 = numeric(), pPC6 = numeric(), pPC7 = numeric(),
                      pPC8 = numeric(), pPC9 = numeric(), pPC10 = numeric()))
  }
  log_info("Loading proteomic PCs from: {pca_result_file}")
  pr <- tryCatch(readRDS(pca_result_file), error = function(e) {
    log_error("Failed to read PCA result file: {e$message}")
    return(NULL)
  })
  if (is.null(pr)) {
    log_warn("Proteomic PCs will not be included as covariates")
    return(data.table(SampleID = character(), pPC1 = numeric(), pPC2 = numeric(), pPC3 = numeric(),
                      pPC4 = numeric(), pPC5 = numeric(), pPC6 = numeric(), pPC7 = numeric(),
                      pPC8 = numeric(), pPC9 = numeric(), pPC10 = numeric()))
  }
  scores <- as.data.table(pr$scores)
  scores$SampleID <- rownames(pr$scores)
  pcs <- names(scores)[grepl("^PC[0-9]+$", names(scores))]
  pcs <- pcs[seq_len(min(n_pcs, length(pcs)))]
  out <- scores[, c("SampleID", pcs), with = FALSE]
  setnames(out, pcs, paste0("pPC", seq_along(pcs)))
  out
}

# Function to identify sex-associated proteins
identify_sex_proteins <- function(npx_matrix, sex_info, id_map, proteomic_pcs) {
  log_info("Identifying sex-associated proteins using logistic regression")

  # Merge sex information with samples
  sample_ids <- rownames(npx_matrix)
  base_dt <- data.table(SAMPLE_ID = sample_ids)

  # Check if sex_info already has SAMPLE_ID (from metadata merge)
  if ("SAMPLE_ID" %in% names(sex_info)) {
    # sex_info already has SAMPLE_ID, merge directly
    sample_sex <- merge(base_dt, sex_info, by = "SAMPLE_ID", all.x = TRUE)
  } else {
    # Map SAMPLE_ID -> FINNGENID (provided mapping)
    map_dt <- merge(base_dt, id_map[, .(SAMPLE_ID, FINNGENID)], by = "SAMPLE_ID", all.x = TRUE)
    # Merge with sex info (keep base order)
    sample_sex <- merge(map_dt, sex_info, by = "FINNGENID", all.x = TRUE)
  }

  # Reorder to base_dt order
  sample_sex <- sample_sex[match(sample_ids, sample_sex$SAMPLE_ID)]

  # Align proteomic PCs (by SAMPLE_ID)
  pcs_aligned <- merge(base_dt, proteomic_pcs, by.x = "SAMPLE_ID", by.y = "SampleID", all.x = TRUE)
  # Reorder to exactly match sample_ids length and order
  pcs_aligned <- pcs_aligned[match(sample_ids, pcs_aligned$SAMPLE_ID)]

  # Remove samples without sex information
  has_sex <- !is.na(sample_sex$genetic_sex)
  npx_matrix_sex <- npx_matrix[has_sex, ]
  sample_sex_clean <- sample_sex[has_sex]
  pcs_aligned <- pcs_aligned[has_sex]

  log_info("Samples with sex information: {sum(has_sex)} out of {length(has_sex)}")

  # Run logistic regression for each protein
  protein_associations <- data.table()

  for(i in seq_len(ncol(npx_matrix_sex))) {
    protein_name <- colnames(npx_matrix_sex)[i]

    # Prepare data for regression
    reg_data <- data.table(
      sex = sample_sex_clean$genetic_sex - 1,  # Convert from 1/2 to 0/1 (0=male, 1=female)
      protein = npx_matrix_sex[, i],
      age = sample_sex_clean$age,
      bmi = sample_sex_clean$bmi,
      smoking = sample_sex_clean$smoking
    )
    # Bind proteomic PCs
    if (ncol(pcs_aligned) > 1) {
      reg_data <- cbind(reg_data, pcs_aligned[, -1])
    }

    # Check minimal non-NA count on key variables; allow NAs in BMI/smoking/pPCs
    n_non_na <- sum(complete.cases(reg_data[, c("sex", "protein", "age")]))
    if(n_non_na < 30) next  # Skip if too few samples with key variables

    # Fit logistic regression with covariates
    tryCatch({
      # Include proteomic PCs as covariates if available
      model <- glm(sex ~ protein + age + bmi + smoking + .,  # remaining columns are pPCs; NAs omitted
                   data = as.data.frame(reg_data),
                   family = binomial())

      # Extract protein coefficient - check if "protein" exists in model coefficients
      # (it may be dropped due to zero variance, perfect separation, or numerical issues)
      coef_names <- rownames(summary(model)$coefficients)
      if (!"protein" %in% coef_names) {
        log_debug("Protein '{protein_name}' coefficient not found in model (likely dropped due to zero variance or convergence issues)")
        next
      }

      coef_summary <- summary(model)$coefficients["protein", ]

      protein_associations <- rbind(protein_associations, data.table(
        protein = protein_name,
        beta = coef_summary[1],
        se = coef_summary[2],
        z_value = coef_summary[3],
        p_value = coef_summary[4],
        n_samples = nrow(reg_data)
      ))
    }, error = function(e) {
      log_debug("Failed to fit model for protein {protein_name}: {e$message}")
    })
    if (i %% 500 == 0) {
      log_info("Logistic regression progress: {i}/{ncol(npx_matrix_sex)} proteins")
    }
  }

  # Calculate adjusted p-values; if empty, return a well-typed empty table (avoids fwrite warning)
  if (nrow(protein_associations) == 0 || !("p_value" %in% names(protein_associations))) {
    log_warn("No valid protein association fits; returning empty associations table")
    empty_assocs <- data.table(
      protein = character(),
      beta = numeric(),
      se = numeric(),
      z_value = numeric(),
      p_value = numeric(),
      n_samples = integer(),
      p_adj = numeric(),
      neg_log10_p = numeric(),
      significant = logical()
    )
    return(list(
      associations = empty_assocs,
      npx_matrix_sex = npx_matrix_sex,
      sample_sex = sample_sex_clean
    ))
  }
  protein_associations[, p_adj := p.adjust(p_value, method = "BH")]
  protein_associations[, neg_log10_p := -log10(p_value)]
  protein_associations[, bonf_sig := p_value < (0.05/5440)]
  protein_associations[, OR := exp(beta)]
  protein_associations[, CI_low := exp(beta - 1.96*se)]
  protein_associations[, CI_high := exp(beta + 1.96*se)]

  # Identify significant sex-associated proteins (FDR)
  sig_threshold <- 0.05
  protein_associations[, significant := p_adj < sig_threshold]

  # Verbose summary of results
  log_info("Logistic regression complete: fit {nrow(protein_associations)} proteins")
  log_info("Proteins with nominal p < 0.05: {sum(protein_associations$p_value < 0.05, na.rm=TRUE)}")
  log_info("Sex-associated proteins (FDR < {sig_threshold}): {sum(protein_associations$significant, na.rm=TRUE)}; Bonf. sig: {sum(protein_associations$bonf_sig, na.rm=TRUE)}")

  return(list(
    associations = protein_associations,
    npx_matrix_sex = npx_matrix_sex,
    sample_sex = sample_sex_clean
  ))
}

# Function to identify smoking-associated proteins
identify_smoking_proteins <- function(npx_matrix, sex_info, id_map, proteomic_pcs) {
  log_info("Identifying smoking-associated proteins using logistic regression")

  # Merge smoking information with samples
  sample_ids <- rownames(npx_matrix)
  base_dt <- data.table(SAMPLE_ID = sample_ids)

  # Check if sex_info already has SAMPLE_ID (from metadata merge)
  if ("SAMPLE_ID" %in% names(sex_info)) {
    # sex_info already has SAMPLE_ID, merge directly
    sample_smoking <- merge(base_dt, sex_info, by = "SAMPLE_ID", all.x = TRUE)
  } else {
    # Map SAMPLE_ID -> FINNGENID (provided mapping)
    map_dt <- merge(base_dt, id_map[, .(SAMPLE_ID, FINNGENID)], by = "SAMPLE_ID", all.x = TRUE)
    # Merge with sex info (keep base order)
    sample_smoking <- merge(map_dt, sex_info, by = "FINNGENID", all.x = TRUE)
  }

  # Reorder to base_dt order
  sample_smoking <- sample_smoking[match(sample_ids, sample_smoking$SAMPLE_ID)]

  # Align proteomic PCs (by SAMPLE_ID)
  pcs_aligned <- merge(base_dt, proteomic_pcs, by.x = "SAMPLE_ID", by.y = "SampleID", all.x = TRUE)
  # Reorder to exactly match sample_ids length and order
  pcs_aligned <- pcs_aligned[match(sample_ids, pcs_aligned$SAMPLE_ID)]

  # Remove samples without smoking information
  has_smoking <- !is.na(sample_smoking$smoking)
  npx_matrix_smoking <- npx_matrix[has_smoking, ]
  sample_smoking_clean <- sample_smoking[has_smoking]
  pcs_aligned <- pcs_aligned[has_smoking]

  log_info("Samples with smoking information: {sum(has_smoking)} out of {length(has_smoking)}")

  # Run logistic regression for each protein
  protein_associations <- data.table()

  for(i in seq_len(ncol(npx_matrix_smoking))) {
    protein_name <- colnames(npx_matrix_smoking)[i]

    # Prepare data for regression
    reg_data <- data.table(
      smoking = sample_smoking_clean$smoking,  # Binary: 0=non-smoker, 1=smoker
      protein = npx_matrix_smoking[, i],
      age = sample_smoking_clean$age,
      sex = sample_smoking_clean$genetic_sex,
      bmi = sample_smoking_clean$bmi
    )
    # Bind proteomic PCs
    if (ncol(pcs_aligned) > 1) {
      reg_data <- cbind(reg_data, pcs_aligned[, -1])
    }

    # Check minimal non-NA count on key variables; allow NAs in BMI/pPCs
    n_non_na <- sum(complete.cases(reg_data[, c("smoking", "protein", "age", "sex")]))
    if(n_non_na < 30) next  # Skip if too few samples with key variables

    # Fit logistic regression with covariates
    tryCatch({
      # Include proteomic PCs as covariates if available
      model <- glm(smoking ~ protein + age + sex + bmi + .,  # remaining columns are pPCs; NAs omitted
                   data = as.data.frame(reg_data),
                   family = binomial())

      # Extract protein coefficient - check if "protein" exists in model coefficients
      # (it may be dropped due to zero variance, perfect separation, or numerical issues)
      coef_names <- rownames(summary(model)$coefficients)
      if (!"protein" %in% coef_names) {
        log_debug("Protein '{protein_name}' coefficient not found in model (likely dropped due to zero variance or convergence issues)")
        next
      }

      coef_summary <- summary(model)$coefficients["protein", ]

      protein_associations <- rbind(protein_associations, data.table(
        protein = protein_name,
        beta = coef_summary[1],
        se = coef_summary[2],
        z_value = coef_summary[3],
        p_value = coef_summary[4],
        n_samples = nrow(reg_data)
      ))
    }, error = function(e) {
      log_debug("Failed to fit model for protein {protein_name}: {e$message}")
    })
    if (i %% 500 == 0) {
      log_info("Logistic regression progress: {i}/{ncol(npx_matrix_smoking)} proteins")
    }
  }

  # Calculate adjusted p-values; if empty, return a well-typed empty table
  if (nrow(protein_associations) == 0 || !("p_value" %in% names(protein_associations))) {
    log_warn("No valid protein association fits; returning empty associations table")
    empty_assocs <- data.table(
      protein = character(),
      beta = numeric(),
      se = numeric(),
      z_value = numeric(),
      p_value = numeric(),
      n_samples = integer(),
      p_adj = numeric(),
      neg_log10_p = numeric(),
      significant = logical(),
      bonf_sig = logical()
    )
    return(list(
      associations = empty_assocs,
      npx_matrix_smoking = npx_matrix_smoking,
      sample_smoking = sample_smoking_clean
    ))
  }
  protein_associations[, p_adj := p.adjust(p_value, method = "BH")]
  protein_associations[, neg_log10_p := -log10(p_value)]
  n_tests <- nrow(protein_associations)
  bonf_threshold <- 0.05 / max(1, n_tests)
  protein_associations[, bonf_sig := p_value < bonf_threshold]
  protein_associations[, OR := exp(beta)]
  protein_associations[, CI_low := exp(beta - 1.96*se)]
  protein_associations[, CI_high := exp(beta + 1.96*se)]

  # Identify significant smoking-associated proteins (FDR)
  sig_threshold <- 0.05
  protein_associations[, significant := p_adj < sig_threshold]

  # Verbose summary of results
  log_info("Logistic regression complete: fit {nrow(protein_associations)} proteins")
  log_info("Proteins with nominal p < 0.05: {sum(protein_associations$p_value < 0.05, na.rm=TRUE)}")
  log_info("Smoking-associated proteins (FDR < {sig_threshold}): {sum(protein_associations$significant, na.rm=TRUE)}; Bonf. sig: {sum(protein_associations$bonf_sig, na.rm=TRUE)}")

  return(list(
    associations = protein_associations,
    npx_matrix_smoking = npx_matrix_smoking,
    sample_smoking = sample_smoking_clean
  ))
}

# Function to predict sex using protein expression
predict_sex_protein <- function(npx_matrix_sex, sample_sex, sex_proteins) {
  log_info("Predicting sex using protein expression")

  # Select top sex-associated proteins
  top_proteins <- sex_proteins$associations[significant == TRUE][order(p_value)][1:min(20, sum(sex_proteins$associations$significant))]$protein

  if(length(top_proteins) < 5) {
    log_warn("Too few significant sex-associated proteins for prediction")
    return(NULL)
  }

  # Create training data
  X <- npx_matrix_sex[, top_proteins]
  y <- sample_sex$genetic_sex - 1  # Convert from 1/2 to 0/1 (0=male, 1=female)

  # Add covariates
  X_with_cov <- cbind(X,
                      age = sample_sex$age,
                      bmi = sample_sex$bmi,
                      smoking = sample_sex$smoking)

  # Remove missing values
  complete_idx <- complete.cases(X_with_cov) & !is.na(y)
  X_train <- X_with_cov[complete_idx, ]
  y_train <- y[complete_idx]

  # Fit logistic regression model
  model <- glm(y_train ~ ., data = as.data.frame(X_train), family = binomial())

  # Predict on all samples
  predictions <- predict(model, newdata = as.data.frame(X_with_cov), type = "response")

  # Create prediction results
  pred_label <- ifelse(predictions >= 0.5, 1, 0)  # 1 indicates genetically female under 0/1 coding
  predicted_sex_str <- ifelse(pred_label == 1, "female", "male")
  genetic_sex_str <- ifelse(sample_sex$genetic_sex == 2, "female",
                            ifelse(sample_sex$genetic_sex == 1, "male", NA_character_))

  prediction_results <- data.table(
    SAMPLE_ID = sample_sex$SAMPLE_ID,
    FINNGENID = sample_sex$FINNGENID,
    genetic_sex = genetic_sex_str,
    predicted_prob = predictions,
    predicted_sex = predicted_sex_str,
    mismatch = FALSE
  )

  # Identify mismatches
  prediction_results[!is.na(genetic_sex) & !is.na(predicted_sex),
                    mismatch := genetic_sex != predicted_sex]

  log_info("Sex prediction complete - Mismatches: {sum(prediction_results$mismatch, na.rm=TRUE)}")

  return(list(
    predictions = prediction_results,
    model = model,
    proteins_used = top_proteins
  ))
}

# Function to create volcano plot
create_volcano_plot <- function(associations) {
  log_info("Creating volcano plot for sex-associated proteins")

  # Prepare data for plot
  plot_data <- associations[!is.na(beta) & !is.na(p_value)]
  plot_data[, neg_log10_p := -log10(p_value)]
  # Label strategy: Bonferroni + top10 by p for each sex-direction
  plot_data[, label := ""]
  # Top 10 male (beta < 0) and top 10 female (beta > 0) by smallest p
  top_male <- plot_data[beta < 0][order(p_value)][1:min(.N, 10)]$protein
  top_fem  <- plot_data[beta > 0][order(p_value)][1:min(.N, 10)]$protein
  plot_data[protein %in% top_male | protein %in% top_fem | bonf_sig == TRUE, label := protein]
  plot_data[, dir := ifelse(beta < 0, "male", "female")]  # beta<0 => higher in males; beta>0 => higher in females

  # Palettes
  pal_male <- as.character(paletteer::paletteer_c("ggthemes::Green-Gold", 30))
  pal_fem  <- as.character(paletteer::paletteer_c("ggthemes::Orange-Gold", 30))

  # Map significance to palette index 1..30
  minp <- min(plot_data$neg_log10_p, na.rm = TRUE)
  maxp <- max(plot_data$neg_log10_p, na.rm = TRUE)
  rng <- ifelse(is.finite(maxp - minp) && (maxp - minp) > 0, maxp - minp, 1)
  plot_data[, idx := pmin(30L, pmax(1L, as.integer(1 + floor(29 * (neg_log10_p - minp) / rng))))]
  plot_data[, color := ifelse(dir == "male", pal_male[idx], pal_fem[idx])]

  # Point sizes by association magnitude (exp(|beta|)) rescaled
  plot_data[, mag := exp(abs(beta))]
  q90 <- as.numeric(quantile(plot_data$mag, 0.9, na.rm = TRUE))
  if (!is.finite(q90) || q90 <= 0) q90 <- max(plot_data$mag, na.rm = TRUE)
  if (!is.finite(q90) || q90 <= 0) q90 <- 1
  plot_data[, pt_size := pmin(4.5, 1 + 3.5 * (mag / q90))]

  p <- ggplot(plot_data, aes(x = beta, y = neg_log10_p)) +
    geom_point(aes(color = color, size = pt_size), alpha = 0.8) +
    scale_color_identity() +
    scale_size_identity() +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray50") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray70") +
    geom_text_repel(aes(label = label), size = 3, max.overlaps = 20) +
    labs(title = "Sex-Associated Proteins",
         subtitle = paste0("Dual palettes: male=Green-Gold, female=Orange-Gold; size=exp(|", "\u03B2", "|)"),
         x = paste0("Effect size (", "\u03B2", ")"),
         y = "-log10(p-value)") +
    theme_bw() +
    theme(legend.position = "none")

  return(p)
}

# Volcano colored by OR (exp(beta)) with Bonferroni annotations
create_volcano_plot_or <- function(associations) {
  log_info("Creating OR-colored volcano plot with Bonferroni labels")
  plot_data <- associations[!is.na(beta) & !is.na(p_value)]
  plot_data[, OR := exp(beta)]
  plot_data[, neg_log10_p := -log10(p_value)]
  n_tests <- nrow(plot_data)
  bonf <- 0.05 / max(1, n_tests)
  plot_data[, bonf_sig := p_value < bonf]
  # Palette
  pal <- as.character(paletteer::paletteer_d("dichromat::BluetoOrange_10"))
  p <- ggplot(plot_data, aes(x = beta, y = neg_log10_p)) +
    geom_point(aes(color = OR), alpha = 0.7, size = 2) +
    scale_color_gradientn(colors = pal) +
    geom_hline(yintercept = -log10(bonf), linetype = "dashed", color = "gray30") +
    geom_text_repel(data = plot_data[bonf_sig == TRUE], aes(label = protein),
                    size = 3, max.overlaps = 25, box.padding = 0.5, point.padding = 0.3) +
    labs(title = "Sex-Associated Proteins colored by Odds Ratio",
         subtitle = "Color = Odds Ratio (exp(beta)); dashed line = Bonferroni threshold",
         x = paste0("Effect size (", "\u03B2", "; >0 higher in genetically female)"),
         y = "-log10(p-value)", color = "Odds Ratio") +
    theme_bw() +
    theme(legend.position = "right")
  p
}

# Volcano plot for smoking-associated proteins aligned with sex-associated styling
create_smoking_volcano_plot <- function(associations) {
  log_info("Creating smoking-associated volcano plot (sex-plot styling)")
  plot_data <- associations[!is.na(beta) & !is.na(p_value)]
  if (nrow(plot_data) == 0) return(ggplot() + theme_void() + ggtitle("No smoking associations"))

  plot_data[, neg_log10_p := -log10(p_value)]

  # Bonferroni threshold and labels
  n_tests <- nrow(plot_data)
  bonf_threshold <- 0.05 / max(1, n_tests)
  if (!("bonf_sig" %in% names(plot_data))) {
    plot_data[, bonf_sig := p_value < bonf_threshold]
  }
  plot_data[, label := ""]
  plot_data[bonf_sig == TRUE, label := protein]

  # Directional palettes (mirroring sex-associated volcano)
  pal_neg <- as.character(paletteer::paletteer_c("ggthemes::Green-Gold", 30))
  pal_pos <- as.character(paletteer::paletteer_c("ggthemes::Orange-Gold", 30))

  plot_data[, dir := ifelse(beta < 0, "lower", "higher")]

  # Map -log10(p) to palette index for visual consistency
  minp <- min(plot_data$neg_log10_p, na.rm = TRUE)
  maxp <- max(plot_data$neg_log10_p, na.rm = TRUE)
  rng <- ifelse(is.finite(maxp - minp) && (maxp - minp) > 0, maxp - minp, 1)
  plot_data[, idx := pmin(30L, pmax(1L, as.integer(1 + floor(29 * (neg_log10_p - minp) / rng))))]
  plot_data[, color := ifelse(dir == "lower", pal_neg[idx], pal_pos[idx])]

  # Point sizes using |β| exp scaling
  plot_data[, mag := exp(abs(beta))]
  q90 <- as.numeric(quantile(plot_data$mag, 0.9, na.rm = TRUE))
  if (!is.finite(q90) || q90 <= 0) q90 <- max(plot_data$mag, na.rm = TRUE)
  if (!is.finite(q90) || q90 <= 0) q90 <- 1
  plot_data[, pt_size := pmin(4.5, 1 + 3.5 * (mag / q90))]

  p <- ggplot(plot_data, aes(x = beta, y = neg_log10_p)) +
    geom_point(aes(color = color, size = pt_size), alpha = 0.8) +
    scale_color_identity() +
    scale_size_identity() +
    geom_hline(yintercept = -log10(bonf_threshold),
               linetype = "dashed", color = "gray50") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray70") +
    geom_text_repel(data = plot_data[bonf_sig == TRUE],
                    aes(label = label),
                    size = 3, max.overlaps = 25, box.padding = 0.5, point.padding = 0.3) +
    labs(title = "Smoking-Associated Proteins",
         subtitle = sprintf("Color palette matches sex volcano; Bonferroni: p < %.2e (n = %d)", bonf_threshold, n_tests),
         x = expression("Effect size (" * beta * "; >0 higher in smokers)"),
         y = expression("-log"[10] * "(p-value)")) +
    theme_bw() +
    theme(legend.position = "none")

  return(p)
}

# Forest plot for significant proteins (Bonferroni)
create_forest_plot <- function(associations) {
  # Use beta scale (log-OR) so male (<0) and female (>0) effects are comparable to volcano
  sig <- associations[bonf_sig == TRUE]
  if (nrow(sig) == 0) return(NULL)
  # Select balanced top effects by |beta|
  pos <- sig[beta > 0][order(-abs(beta))][1:min(.N, 10)]
  neg <- sig[beta < 0][order(-abs(beta))][1:min(.N, 10)]
  top20 <- rbind(pos, neg, fill = TRUE)
  if (nrow(top20) == 0) return(NULL)
  top20[, sex_label := ifelse(beta > 0, "Higher in Female", "Higher in Male")]
  # 95% CI on beta scale
  top20[, beta_low := beta - 1.96 * se]
  top20[, beta_high := beta + 1.96 * se]
  # Order by beta (female to the right)
  top20 <- top20[order(beta)]
  top20[, protein := factor(protein, levels = protein)]
  p <- ggplot(top20, aes(x = beta, y = protein, color = sex_label)) +
    geom_point(size = 2) +
    geom_errorbarh(aes(xmin = beta_low, xmax = beta_high), height = 0.2) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    scale_color_manual(values = c("Higher in Male" = "#3A5F8A", "Higher in Female" = "#FF6B6B")) +
    labs(title = "Top sex-associated proteins (Bonferroni)",
         subtitle = "Top 10 higher in females (beta > 0) and top 10 higher in males (beta < 0)",
         x = paste0("Effect size (", "\u03B2", "; >0 higher in genetically female)"), y = "Protein", color = "Association") +
    theme_bw()
  p
}

# Forest plot for smoking-associated proteins
create_smoking_forest_plot <- function(associations) {
  sig <- associations[bonf_sig == TRUE]
  if (nrow(sig) == 0) return(NULL)

  pos <- sig[beta > 0][order(-abs(beta))][1:min(.N, 10)]
  neg <- sig[beta < 0][order(-abs(beta))][1:min(.N, 10)]
  top20 <- rbind(pos, neg, fill = TRUE)
  if (nrow(top20) == 0) return(NULL)

  top20[, smoke_label := ifelse(beta > 0, "Higher in Smokers", "Lower in Smokers")]
  top20[, beta_low := beta - 1.96 * se]
  top20[, beta_high := beta + 1.96 * se]
  top20 <- top20[order(beta)]
  top20[, protein := factor(protein, levels = protein)]

  ggplot(top20, aes(x = beta, y = protein, color = smoke_label)) +
    geom_point(size = 2) +
    geom_errorbarh(aes(xmin = beta_low, xmax = beta_high), height = 0.2) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    scale_color_manual(values = c("Lower in Smokers" = "#3A5F8A", "Higher in Smokers" = "#FF6B6B")) +
    labs(title = "Top smoking-associated proteins (Bonferroni)",
         subtitle = "Top 10 higher and top 10 lower in smokers",
         x = expression("Effect size (" * beta * "; >0 higher in smokers)"),
         y = "Protein", color = "Association") +
    theme_bw()
}

# Faceted histograms for the top 10 Bonferroni-significant proteins by |beta|
create_top_hist_facets <- function(npx_matrix, associations, sample_sex, feature_list = NULL, top_n = 10) {
  # If feature_list is provided, use it directly; otherwise fall back to Bonferroni-significant by |beta|
  if (is.null(feature_list)) {
  sig <- associations[bonf_sig == TRUE]
  if (nrow(sig) == 0) return(NULL)
    top10 <- sig[order(-abs(beta))][1:min(top_n, .N)]$protein
  } else {
    top10 <- feature_list[1:min(top_n, length(feature_list))]
  }
  if (length(top10) == 0) return(NULL)

  # Align SAMPLE_IDs between matrix and sample_sex; restrict to samples with genetic sex available
  all_ids <- rownames(npx_matrix)
  ss <- as.data.table(sample_sex[, c("SAMPLE_ID", "genetic_sex")])
  ss[, genetic_sex := ifelse(genetic_sex == 2, "female",
                             ifelse(genetic_sex == 1, "male", NA_character_))]
  ss <- ss[!is.na(genetic_sex)]
  ids <- intersect(all_ids, ss$SAMPLE_ID)
  if (length(ids) == 0) return(NULL)

  # Build long data only for aligned IDs
  ss_aligned <- ss[match(ids, ss$SAMPLE_ID)]
  long_dt <- rbindlist(lapply(top10, function(p) {
    vals <- as.numeric(npx_matrix[ids, p])
    data.table(SAMPLE_ID = ids,
               protein = p,
               value = vals,
               genetic_sex = ss_aligned$genetic_sex)
  }), use.names = TRUE)
  long_dt <- long_dt[!is.na(genetic_sex)]
  long_dt[, genetic_sex := factor(genetic_sex, levels = c("male", "female"))]

  p <- gghistogram(
    long_dt, x = "value", y = "..density..",
    add = "mean", rug = TRUE,
    fill = "genetic_sex", palette = c("#00AFBB", "#E7B800"),
    add_density = TRUE, facet.by = "protein"
  ) +
    labs(title = sprintf("Distribution of top %d sex-differentiating proteins", length(top10)),
         x = "NPX", y = "Density", fill = "Imputed sex")

  p
}

# Scatter comparing logistic vs elastic-net effect sizes for top-50 proteins
create_effect_scatter_top50 <- function(associations, enet_coef_summary) {
  if (nrow(associations) == 0 || nrow(enet_coef_summary) == 0) return(NULL)
  top50 <- associations[order(p_value)][1:min(50, .N)]
  top50[, beta_low := beta - 1.96 * se]
  top50[, beta_high := beta + 1.96 * se]
  dt <- merge(top50[, .(protein, beta, beta_low, beta_high)],
              enet_coef_summary, by = "protein", all.x = TRUE)
  if (nrow(dt) == 0) return(NULL)

  p <- ggplot(dt, aes(x = beta, y = enet_beta_mean)) +
    geom_point(alpha = 0.8) +
    geom_errorbarh(aes(xmin = beta_low, xmax = beta_high), alpha = 0.4, height = 0) +
    geom_errorbar(aes(ymin = enet_beta_low, ymax = enet_beta_high), alpha = 0.4, width = 0) +
    geom_smooth(method = "lm", se = TRUE, color = "#3A5F8A") +
    labs(title = "Logistic Regression vs Elastic Net effect sizes (Top 50 by p-value)",
         x = paste0("Logistic ", "\u03B2", " (95% CI)"),
         y = paste0("Elastic Net ", "\u03B2", " (outer cross-validation 95% range)")) +
    theme_bw()

  # Annotate top 10 male-associated (beta<0) and top 10 female-associated (beta>0)
  lab_male <- dt[beta < 0][order(-abs(beta))][1:min(.N, 10)]
  lab_fem  <- dt[beta > 0][order(-abs(beta))][1:min(.N, 10)]
  labs_dt <- rbind(lab_male, lab_fem, fill = TRUE)
  if (nrow(labs_dt) > 0) {
    p <- p + ggrepel::geom_text_repel(
      data = labs_dt,
      aes(label = protein),
      size = 2.6, max.overlaps = 30
    )
  }

  p
}

# Correlation heatmap of top-N associated proteins with significance stars
create_corr_heatmap <- function(npx_matrix, associations, id_map, sex_info, top_n = 50) {
  sig <- associations[bonf_sig == TRUE]
  if (nrow(sig) == 0) return(NULL)
  top <- sig[order(p_value)][1:min(top_n, .N)]$protein
  sub <- npx_matrix[, top, drop = FALSE]
  # Spearman correlation
  R <- suppressWarnings(cor(sub, use = "pairwise.complete.obs", method = "spearman"))
  # P-values via cor.test pairwise
  pmat <- matrix(NA_real_, nrow = ncol(sub), ncol = ncol(sub), dimnames = list(colnames(sub), colnames(sub)))
  for (i in seq_len(ncol(sub))) for (j in seq_len(ncol(sub))) {
    if (i == j) { pmat[i,j] <- 0; next }
    ct <- try(suppressWarnings(cor.test(sub[, i], sub[, j], method = "spearman")), silent = TRUE)
    pmat[i,j] <- if (inherits(ct, "try-error")) NA_real_ else ct$p.value
  }
  stars <- matrix("", nrow = nrow(pmat), ncol = ncol(pmat), dimnames = dimnames(pmat))
  stars[pmat < 0.05 & !is.na(pmat)] <- "*"
  stars[pmat < 0.01 & !is.na(pmat)] <- "**"
  stars[pmat < 0.001 & !is.na(pmat)] <- "***"
  cols <- colorRampPalette(c("#4393C3", "#F7F7F7", "#D6604D"))(50)
  pheatmap(R, color = cols, display_numbers = stars, number_color = "black",
           border_color = NA, main = "Correlation of sex-associated proteins (Spearman)")
}

# Function to create sex prediction plot
create_prediction_plot <- function(predictions) {
  log_info("Creating sex prediction visualization")

  # Prepare data
  plot_data <- predictions[!is.na(genetic_sex) & !is.na(predicted_prob)]
  plot_data[, genetic_sex_label := ifelse(genetic_sex == "male", "Male", "Female")]
  plot_data[, mismatch_label := ifelse(mismatch, "Mismatch", "Match")]

  # Create density plot
  p <- ggplot(plot_data, aes(x = predicted_prob, fill = genetic_sex_label)) +
    geom_histogram(aes(y = ..density..), bins = 30, alpha = 0.6, position = "identity") +
    geom_density(alpha = 0.3) +
    geom_vline(xintercept = 0.5, linetype = "dashed", color = "black") +
    facet_wrap(~ genetic_sex_label, ncol = 1) +
    scale_fill_manual(values = c("Male" = "#3A5F8A", "Female" = "#FF6B6B")) +
    labs(title = "Sex Prediction from Protein Expression",
         x = "Predicted Female Probability",
         y = "Density",
         fill = "Genetic Sex") +
    theme_bw()

  # Highlight mismatches
  mismatch_data <- plot_data[mismatch == TRUE]
  if(nrow(mismatch_data) > 0) {
    p <- p + geom_point(data = mismatch_data,
                       aes(x = predicted_prob, y = 0),
                       color = "red", size = 3, shape = 17)
  }

  return(p)
}

# ============================================================================
# DEPRECATED: 3-sigma threshold approach (COMMENTED OUT - DO NOT USE)
# ============================================================================
# This function computed sex-specific 3σ outlier thresholds from CV predictions.
# ISSUE: Using 3σ from mean probability distribution is statistically incorrect
# for binary classification tasks. It conflates prediction uncertainty with
# true misclassification.
#
# REPLACEMENT: Use Youden's J threshold (optimized decision boundary from ROC curve)
# which properly balances sensitivity and specificity.
# ============================================================================
# learn_sex_specific_thresholds <- function(y_true_01, y_pred_prob) {
#   # y_true_01: 0=male, 1=female; y_pred_prob: predicted female probability from CV (held-out)
#   dt <- data.table(y = as.integer(y_true_01), p = as.numeric(y_pred_prob))
#   dt <- dt[!is.na(y) & !is.na(p)]
#   male <- dt[y == 0]$p
#   fem  <- dt[y == 1]$p
#   mu_m <- mean(male, na.rm = TRUE); sd_m <- sd(male, na.rm = TRUE)
#   mu_f <- mean(fem,  na.rm = TRUE); sd_f <- sd(fem,  na.rm = TRUE)
#   # 3-sigma rule on CV predictions (model-informed)
#   male_upper <- mu_m + 3 * sd_m
#   female_lower <- mu_f - 3 * sd_f
#   # Clamp to valid [0,1]
#   male_upper <- min(1, max(0, male_upper))
#   female_lower <- min(1, max(0, female_lower))
#   separation <- mu_f - mu_m
#   list(male_upper = male_upper, female_lower = female_lower, separation = separation,
#        mu_male = mu_m, mu_female = mu_f, sd_male = sd_m, sd_female = sd_f, n = nrow(dt))
# }

# Two-panel distribution: left = density/hist by sex; right = probability strip with outliers and optional FINNGENID labels
# Thresholds: 0.5 (mismatch), Youden's J (sex_outlier)
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
  # ============================================================================
  # UPDATED LOGIC (per user request):
  # - mismatch: Deviates from Youden's J threshold (FULL RED - more stringent)
  # - sex_outlier: Deviates from 0.5 but not from Youden's J (PALE RED - mild)
  # ============================================================================
  dt[, is_mismatch := FALSE]
  dt[, is_sex_outlier := FALSE]
  # Mismatch: crosses Youden's J threshold
  dt[genetic_sex == "male" & predicted_prob >= youden_thresh, is_mismatch := TRUE]
  dt[genetic_sex == "female" & predicted_prob < youden_thresh, is_mismatch := TRUE]
  # Sex outlier: crosses 0.5 but not Youden's J
  dt[genetic_sex == "male" & predicted_prob >= 0.5 & predicted_prob < youden_thresh, is_sex_outlier := TRUE]
  # Note: For females, when youden_thresh > 0.5, there's no "between" range
  # Females with prob < youden_thresh are already mismatches

  # Classification: normal, outlier (pale red, between 0.5 and Youden), mismatch (red, beyond Youden)
  dt[, outlier_class := "normal"]
  dt[is_sex_outlier == TRUE, outlier_class := "outlier"]
  dt[is_mismatch == TRUE, outlier_class := "mismatch"]

  # Left: density + histogram faceted by sex
  p_left <- ggplot(dt, aes(x = predicted_prob, fill = genetic_sex)) +
    geom_histogram(aes(y = ..density..), bins = 30, alpha = 0.6, position = "identity") +
    geom_density(alpha = 0.25, color = "black", linewidth = 0.3) +
    geom_vline(xintercept = 0.5, linetype = "dashed", color = "gray50", linewidth = 0.4) +
    geom_vline(xintercept = youden_thresh, linetype = "dotted", color = "#D7301F", linewidth = 0.8) +
    scale_fill_manual(values = c("male" = "#3A5F8A", "female" = "#FF6B6B")) +
    facet_wrap(~ genetic_sex, ncol = 1, scales = "free_y") +
    labs(title = title,
         subtitle = if (is.null(subtitle)) sprintf("Gray dashed = Outlier threshold (0.5); Red dotted = Mismatch threshold (Youden J = %.3f)", youden_thresh) else subtitle,
         x = "Predicted Female Probability", y = "Density", fill = "Genetic Sex") +
    theme_bw()

  # Right: probability strip by sex with jitter; no boxplots
  # use numeric y with jitter to reduce stacking; format y axis as categorical ticks
  p_right <- ggplot(dt, aes(x = predicted_prob, y = y_jit)) +
    # decision thresholds
    geom_vline(xintercept = 0.5, linetype = "dashed", color = "gray50", linewidth = 0.4) +
    geom_vline(xintercept = youden_thresh, linetype = "dotted", color = "#D7301F", linewidth = 0.8) +
    # jittered points - normal samples
    geom_point(data = dt[outlier_class == "normal"],
               size = 0.9, alpha = 0.22, color = "gray20", show.legend = FALSE) +
    # outliers (between 0.5 and Youden J) - pale red
    geom_point(data = dt[outlier_class == "outlier"],
               shape = 21, size = 2.0, stroke = 0.3,
               fill = "#D7301F", alpha = 0.35, color = "black", show.legend = FALSE) +
    # mismatches (beyond Youden J) - full red
    geom_point(data = dt[outlier_class == "mismatch"],
               shape = 21, size = 2.4, stroke = 0.5,
               fill = "#D7301F", color = "black", show.legend = FALSE) +
    scale_y_continuous(breaks = c(1, 2), labels = c("male", "female"), limits = c(0.6, 2.4), expand = expansion(mult = c(0.02, 0.02))) +
    coord_cartesian(xlim = c(0, 1)) +
    theme_bw() +
    labs(x = "Predicted Female Probability", y = NULL) +
    theme(legend.position = "none") +
    # Add legend annotation in top-left corner
    annotate("rect", xmin = 0.02, xmax = 0.30, ymin = 2.25, ymax = 2.38, fill = "white", color = "black", linewidth = 0.3) +
    annotate("text", x = 0.04, y = 2.35, label = "Sex QC Flags", hjust = 0, size = 2.8, fontface = "bold") +
    annotate("point", x = 0.05, y = 2.31, shape = 21, size = 2.4, stroke = 0.5, fill = "#D7301F", color = "black") +
    annotate("text", x = 0.07, y = 2.31, label = "Mismatches", hjust = 0, size = 2.4) +
    annotate("point", x = 0.05, y = 2.27, shape = 21, size = 2.0, stroke = 0.3, fill = "#D7301F", alpha = 0.35, color = "black") +
    annotate("text", x = 0.07, y = 2.27, label = "Outliers", hjust = 0, size = 2.4)

  if (annotate_ids) {
    # Label mismatches (most severe outliers)
    lab_dt <- dt[is_mismatch == TRUE & !is.na(FINNGENID)]
    if (nrow(lab_dt) > max_labels) {
      # choose far-from-threshold labels for readability
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
  # Titles can be long; allow multi-line by using labs on the combined plot via annotate
  pg <- ggpubr::ggarrange(p_left, p_right, ncol = 2, widths = c(2, 1))
  dir.create(dirname(file_out), recursive = TRUE, showWarnings = FALSE)
  try(ggsave(file_out, pg, width = 12, height = 6), silent = TRUE)
  pg
}

# Map model code to descriptive display name
model_display_name <- function(code, ensemble_members = NULL) {
  cl <- tolower(code)
  if (cl == "enet") return("Elastic Net Logistic Regression (L1 + L2 regularization)")
  if (cl == "ridge") return("Ridge Logistic Regression (L2 regularization)")
  if (cl == "lasso") return("Lasso Logistic Regression (L1 regularization)")
  if (cl == "logl2") return("Logistic Regression (L2 regularization)")
  if (cl == "rf") return("Random Forest Classifier")
  if (cl == "xgb") return("Extreme Gradient Boosting (XGBoost) Classifier")
  if (cl == "mlp") return("Multilayer Perceptron Neural Network")
  if (cl == "ae") return("Autoencoder-based Model")
  if (cl == "enet_main") return("Main Pipeline (Elastic Net Logistic Regression)")
  if (cl == "ensemble") {
    if (!is.null(ensemble_members) && length(ensemble_members) > 0) {
      return(paste0("Ensemble of ", paste(ensemble_members, collapse = ", ")))
    }
    return("Ensemble of available base models")
  }
  toupper(code)
}
# Create a distribution plot like 05.1_sex_outliers.R (histogram with sex-specific thresholds)
create_prediction_distribution_legacy <- function(predicted_probs, genetic_sex_01, outdir) {
  dt_pred <- data.table(sex = as.integer(genetic_sex_01), predicted = as.numeric(predicted_probs))
  dt_thresh <- dt_pred[, .(mp = mean(predicted, na.rm = TRUE), sdp = sd(predicted, na.rm = TRUE)), by = sex]
  dt_thresh[, low_p := mp - 3 * sdp]
  dt_thresh[, up_p := mp + 3 * sdp]
  dt_thresh[, sex := factor(sex, levels = c(0, 1), labels = c("male", "female"))]
  dt_plot <- melt(dt_thresh, id.vars = "sex", measure.vars = c("low_p", "up_p"), variable.name = "threshold", value.name = "thresh")

  dt_pred_plot <- copy(dt_pred)
  dt_pred_plot[, sex := factor(sex, levels = c(0, 1), labels = c("male", "female"))]
  p <- ggplot(dt_pred_plot, aes(x = predicted, fill = sex)) +
    geom_histogram(binwidth = 0.05, alpha = 0.6, position = "identity") +
    labs(x = "Predicted value", fill = "Sex", y = "Count") +
    geom_vline(data = dt_plot, aes(xintercept = thresh), linetype = "dashed", alpha = 0.6) +
    theme_bw()
  try(ggsave(file.path(outdir, "sex_predict_distribution.pdf"), p, width = 5, height = 3), silent = TRUE)
}

# Lightweight utility used by enhanced models
standardize_features <- function(X_train, X_test = NULL) {
  means <- colMeans(X_train, na.rm = TRUE)
  sds <- apply(X_train, 2, sd, na.rm = TRUE)
  sds[sds == 0 | is.na(sds)] <- 1
  X_train_sc <- sweep(sweep(X_train, 2, means, "-"), 2, sds, "/")
  if (is.null(X_test)) return(list(train = X_train_sc, means = means, sds = sds))
  X_test_sc <- sweep(sweep(X_test, 2, means, "-"), 2, sds, "/")
  list(train = X_train_sc, test = X_test_sc, means = means, sds = sds)
}

# Simple within-fold feature selection (ensemble of ranks)
select_features_ensemble <- function(X, y, top_n = 150) {
  n_features <- ncol(X)
  if (n_features <= top_n) {
    return(list(selected = colnames(X)))
  }
  fs <- data.table(feature = colnames(X))
  # p-values
  p_values <- sapply(seq_len(ncol(X)), function(i) {
    xi <- X[, i]
    if (sd(xi, na.rm = TRUE) == 0) return(1)
    p <- tryCatch(t.test(xi[y == 0], xi[y == 1])$p.value, error = function(e) 1)
    ifelse(is.finite(p), p, 1)
  })
  fs[, p_value := p_values][, p_rank := rank(p_value, ties.method = "average")]
  # effect size
  eff <- sapply(seq_len(ncol(X)), function(i) {
    xi <- X[, i]
    if (sd(xi, na.rm = TRUE) == 0) return(0)
    m0 <- mean(xi[y == 0], na.rm = TRUE); m1 <- mean(xi[y == 1], na.rm = TRUE)
    s0 <- sd(xi[y == 0], na.rm = TRUE); s1 <- sd(xi[y == 1], na.rm = TRUE)
    n0 <- sum(y == 0); n1 <- sum(y == 1)
    pooled <- sqrt(((n0 - 1) * s0^2 + (n1 - 1) * s1^2) / max(1, n0 + n1 - 2))
    if (!is.finite(pooled) || pooled == 0) return(0)
    abs(m1 - m0) / pooled
  })
  fs[, effect_size := eff][, effect_rank := rank(-effect_size, ties.method = "average")]
  # correlation
  cs <- sapply(seq_len(ncol(X)), function(i) {
    xi <- X[, i]
    if (sd(xi, na.rm = TRUE) == 0) return(0)
    v <- tryCatch(abs(cor(xi, y, use = "complete.obs")), error = function(e) 0)
    ifelse(is.finite(v), v, 0)
  })
  fs[, cor_score := cs][, cor_rank := rank(-cor_score, ties.method = "average")]
  # variance
  vs <- apply(X, 2, var, na.rm = TRUE)
  fs[, variance := ifelse(is.finite(vs), vs, 0)][, var_rank := rank(-variance, ties.method = "average")]
  # ensemble
  fs[, ensemble_rank := (p_rank + effect_rank + cor_rank + var_rank) / 4]
  sel <- fs[order(ensemble_rank)][1:min(top_n, n_features)]$feature
  list(selected = sel, scores = fs)
}

# Evaluate predictions quickly
evaluate_predictions_quick <- function(y_true, y_pred) {
  roc_obj <- try(pROC::roc(response = y_true, predictor = y_pred, quiet = TRUE), silent = TRUE)
  auc_val <- if (inherits(roc_obj, "try-error")) NA_real_ else as.numeric(pROC::auc(roc_obj))
  pr_obj <- try(PRROC::pr.curve(scores.class0 = y_pred[y_true == 1], scores.class1 = y_pred[y_true == 0], curve = FALSE), silent = TRUE)
  pr_auc <- if (inherits(pr_obj, "try-error")) NA_real_ else pr_obj$auc.integral
  list(auc = auc_val, pr_auc = pr_auc)
}

# Function to integrate protein and genetic evidence
integrate_sex_evidence <- function(protein_predictions, genetic_sex) {
  log_info("Integrating protein-based and genetic sex evidence")

  # Combine evidence
  integrated <- merge(
    protein_predictions,
    genetic_sex[, .(FINNGENID, genetic_sex_check = genetic_sex)],
    by = "FINNGENID",
    all.x = TRUE
  )

  # Flag samples with consistent mismatches
  integrated[, confident_mismatch := mismatch & !is.na(genetic_sex_check)]

  # Calculate confidence score
  integrated[, confidence := abs(predicted_prob - 0.5) * 2]  # Scale to 0-1

  # Identify high-confidence mismatches
  high_conf_threshold <- 0.95
  integrated[, high_conf_mismatch := confident_mismatch & confidence > high_conf_threshold]

  log_info("High-confidence sex mismatches: {sum(integrated$high_conf_mismatch, na.rm=TRUE)}")

  return(integrated)
}

# Build disease group column from indicator columns (copied from 04_pca_outliers.R)
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

# Helper function to add BIOBANK_PLASMA and DISEASE_GROUP to predictions data.table
add_metadata_to_predictions <- function(preds_dt, metadata) {
  if (is.null(metadata) || nrow(preds_dt) == 0) return(preds_dt)

  # Add DISEASE_GROUP to metadata if not already present
  metadata_with_disease <- add_disease_group(metadata)

  # Select columns to merge
  meta_cols <- intersect(c("SAMPLE_ID", "BIOBANK_PLASMA", "DISEASE_GROUP"), colnames(metadata_with_disease))
  if (length(meta_cols) == 0 || !("SAMPLE_ID" %in% meta_cols)) {
    log_warn("Cannot add BIOBANK_PLASMA/DISEASE_GROUP: SAMPLE_ID not found in metadata")
    return(preds_dt)
  }

  # Merge with predictions
  meta_sel <- metadata_with_disease[, ..meta_cols]
  preds_with_meta <- merge(preds_dt, meta_sel, by.x = "SAMPLE_ID", by.y = "SAMPLE_ID", all.x = TRUE)

  # Rename SAMPLE_ID to SampleID to match output format (if needed)
  if ("SAMPLE_ID" %in% names(preds_with_meta)) {
    setnames(preds_with_meta, "SAMPLE_ID", "SampleID")
  }

  # Ensure column order: SampleID, FINNGENID, then rest
  if ("SampleID" %in% names(preds_with_meta) && "FINNGENID" %in% names(preds_with_meta)) {
    other_cols <- setdiff(names(preds_with_meta), c("SampleID", "FINNGENID"))
    setcolorder(preds_with_meta, c("SampleID", "FINNGENID", other_cols))
  }

  return(preds_with_meta)
}

# Main execution
main <- function() {

  # Load data from previous steps (batch-aware paths)
  log_info("Loading data from previous steps")
  prev_step01_num <- "01"
  prev_step00_num <- "00"
  prev_step03_num <- "03"  # Z-score outliers step
  prev_step04_num <- "01"   # PCA outliers step (for PCA results)

  # CRITICAL: Use PCA-cleaned matrix to match original pipeline implementation
  # Original pipeline uses PCA-cleaned matrix (2,498 samples) for sex detection
  # This ensures consistency with original pipeline and proper model training
  npx_matrix_path <- get_output_path(prev_step04_num, "npx_matrix_pca_cleaned", batch_id, "outliers", config = config)
  if (!file.exists(npx_matrix_path)) {
    # Fallback to analysis-ready if PCA-cleaned not available
    npx_matrix_path <- get_output_path(prev_step00_num, "npx_matrix_analysis_ready", batch_id, "qc", config = config)
  }
  if (!file.exists(npx_matrix_path)) {
    stop("No NPX matrix file found for sex outlier step")
  }
  log_info("Using NPX matrix: {npx_matrix_path}")
  npx_matrix <- readRDS(npx_matrix_path)

  metadata_path <- get_output_path(prev_step00_num, "metadata", batch_id, "qc", config = config)
  metadata <- readRDS(metadata_path)

  # Also load bridge_metadata_enhanced from step 00 to ensure all bridging samples are included
  prev_step00_num <- "00"
  bridge_metadata_path <- get_output_path(prev_step01_num, "bridge_metadata_enhanced", batch_id, "normalized", config = config)
  if (file.exists(bridge_metadata_path)) {
    log_info("Loading enhanced bridge metadata to ensure all bridging samples are included")
    bridge_metadata_enhanced <- readRDS(bridge_metadata_path)
    # Merge bridge metadata with main metadata, prioritizing main metadata but filling gaps
    if (nrow(bridge_metadata_enhanced) > 0) {
      # Find bridging samples not in main metadata
      bridge_not_in_metadata <- bridge_metadata_enhanced[!SAMPLE_ID %in% metadata$SAMPLE_ID]
      if (nrow(bridge_not_in_metadata) > 0) {
        log_info("Adding {nrow(bridge_not_in_metadata)} bridging samples from bridge_metadata_enhanced to metadata")
        # Ensure bridge_not_in_metadata has same structure as metadata
        common_cols <- intersect(names(metadata), names(bridge_not_in_metadata))
        if (length(common_cols) < length(names(metadata))) {
          # Add missing columns with NA
          missing_cols <- setdiff(names(metadata), names(bridge_not_in_metadata))
          for (col in missing_cols) {
            bridge_not_in_metadata[, (col) := NA]
          }
          # Reorder to match metadata column order
          setcolorder(bridge_not_in_metadata, names(metadata))
        }
        metadata <- rbind(metadata, bridge_not_in_metadata, fill = TRUE)
      }
    }
  }

  id_map <- metadata[, .(SAMPLE_ID, FINNGENID)]
  # Ensure unique SAMPLE_ID mapping
  if (any(duplicated(id_map$SAMPLE_ID))) {
    id_map <- id_map[!duplicated(SAMPLE_ID)]
  }
  # Prefer comprehensive sample mapping if available
  map_file <- get_output_path(prev_step00_num, "sample_mapping", batch_id, "qc", config = config)
  if (file.exists(map_file)) {
    sm <- readRDS(map_file)
    if (all(c("SampleID","FINNGENID") %in% names(sm))) {
      id_map <- data.table(SAMPLE_ID = sm$SampleID, FINNGENID = sm$FINNGENID)
      if (any(duplicated(id_map$SAMPLE_ID))) {
        id_map <- id_map[!duplicated(SAMPLE_ID)]
      }
    }
  }

  # Load covariates with sex information and sample collection age
  # Sex information is now loaded from the covariate file (SEX_IMPUTED preferred, SEX as fallback)
  # Fallback to finngen_R13_minimum_file for missing FINNGENIDs
  finngen_r13_minimum_file <- tryCatch(config$covariates$finngen_r13_minimum_file, error = function(e) NULL)
  sex_info <- load_covariates(config$covariates$covariate_file, metadata, finngen_r13_minimum_file)

  # Load proteomic PCs from 04_pca_outliers (batch-aware path)
  pca_result_path <- get_output_path(prev_step04_num, "pca_result", batch_id, "outliers", config = config)
  proteomic_pcs <- load_proteomic_pcs(pca_result_path, n_pcs = 10)

  # Check if proteomic PCs were loaded successfully
  if (nrow(proteomic_pcs) == 0) {
    log_warn("No proteomic PCs available - proceeding without PC covariates")
    log_warn("This may occur if step 01 (PCA outliers) failed or was not run")
  }

  # Identify sex-associated proteins (univariate + confounders)
  sex_proteins <- identify_sex_proteins(npx_matrix, sex_info, id_map, proteomic_pcs)

  # Identify smoking-associated proteins (univariate + confounders)
  smoking_proteins <- identify_smoking_proteins(npx_matrix, sex_info, id_map, proteomic_pcs)

  # Predictive modeling: nested CV elastic-net with residualization + unpenalized confounders
  sample_ids <- rownames(npx_matrix)

  # Merge sex_info with sample_ids - handle both SAMPLE_ID and FINNGENID cases
  if ("SAMPLE_ID" %in% names(sex_info)) {
    # sex_info has SAMPLE_ID, merge directly to preserve sample-specific ages
    pheno_dt <- merge(data.table(SAMPLE_ID = sample_ids),
                      sex_info[, .(SAMPLE_ID, FINNGENID, genetic_sex, age, bmi, smoking)],
                      by = "SAMPLE_ID", all.x = TRUE)
  } else {
    # sex_info only has FINNGENID, merge through id_map
  dt_map <- merge(data.table(SAMPLE_ID = sample_ids), id_map, by = "SAMPLE_ID", all.x = TRUE)
    pheno_dt <- merge(dt_map, sex_info[, .(FINNGENID, genetic_sex, age, bmi, smoking)],
                      by = "FINNGENID", all.x = TRUE)
  }

  # Ensure order matches sample_ids
  pheno_dt <- pheno_dt[match(sample_ids, pheno_dt$SAMPLE_ID)]

  # Exclude samples with missing genetic sex, F64, and chromosomal abnormality cohorts from training
  log_info("Applying exclusions: missing genetic sex, F64, Chromosomal_Abnormalities cohorts (if available)")
  # Attach cohort flags from metadata when present
  chrom_meta_col <- names(metadata)[tolower(names(metadata)) == "chromosomal_abnormalities"]
  meta_cols <- intersect(names(metadata), unique(c("F64", "DISEASE_GROUP", "BIOBANK_PLASMA", chrom_meta_col)))
  if (length(meta_cols) > 0) {
    add_cols <- c("SAMPLE_ID", meta_cols)
    add_cols <- add_cols[add_cols %in% names(metadata)]
    pheno_dt <- merge(pheno_dt, unique(metadata[, ..add_cols]), by = "SAMPLE_ID", all.x = TRUE)
  }
  # Build exclusion flags
  excl_missing_sex <- is.na(pheno_dt$genetic_sex)
  # F64 can appear as a numeric/logical indicator column or via DISEASE_GROUP == "F64"
  excl_f64 <- rep(FALSE, nrow(pheno_dt))
  if ("F64" %in% names(pheno_dt)) {
    v <- pheno_dt$F64
    if (is.logical(v)) {
      excl_f64 <- !is.na(v) & v
    } else if (is.numeric(v) || is.integer(v)) {
      excl_f64 <- !is.na(v) & (v != 0)
    } else if (is.character(v)) {
      lv <- tolower(v)
      excl_f64 <- !is.na(lv) & (lv %in% c("1","true","t","yes","y"))
    }
  }
  if ("DISEASE_GROUP" %in% names(pheno_dt)) {
    dg <- toupper(as.character(pheno_dt$DISEASE_GROUP))
    excl_f64 <- excl_f64 | (!is.na(dg) & dg == "F64")
  }
  # Chromosomal abnormalities as a dedicated boolean column or via BIOBANK_PLASMA text
  excl_chrom <- rep(FALSE, nrow(pheno_dt))
  chrom_pheno_col <- names(pheno_dt)[tolower(names(pheno_dt)) == "chromosomal_abnormalities"]
  if (length(chrom_pheno_col) > 0) {
    v <- pheno_dt[[chrom_pheno_col[1]]]
    if (is.logical(v)) {
      excl_chrom <- !is.na(v) & v
    } else if (is.numeric(v) || is.integer(v)) {
      excl_chrom <- !is.na(v) & (v != 0)
    } else if (is.character(v)) {
      lv <- tolower(v)
      excl_chrom <- !is.na(lv) & (lv %in% c("1","true","t","yes","y"))
    }
  } else if ("BIOBANK_PLASMA" %in% names(pheno_dt)) {
    bp <- as.character(pheno_dt$BIOBANK_PLASMA)
    lbp <- tolower(bp)
    excl_chrom <- !is.na(lbp) & (lbp %in% c("chromosomal_abnormalities","chromosomal abnormalities","chromosomal-abnormalities"))
  }
  excl_any <- (excl_missing_sex | excl_f64 | excl_chrom)
  log_info("Exclusion counts - missing_sex: {sum(excl_missing_sex, na.rm=TRUE)}, F64: {sum(excl_f64, na.rm=TRUE)}, Chromosomal_Abnormalities: {sum(excl_chrom, na.rm=TRUE)}")
  keep_mask <- !excl_any

  y_raw <- pheno_dt$genetic_sex
  # Convert from 1/2 to 0/1 for modeling (0=male, 1=female)
  y_bin <- ifelse(y_raw == 2, 1, ifelse(y_raw == 1, 0, NA))
  keep <- which(!is.na(y_bin) & keep_mask)
  Y <- npx_matrix[keep, , drop = FALSE]
  sample_ids_keep <- sample_ids[keep]
  cov_dt <- pheno_dt[keep, .(SAMPLE_ID, FINNGENID, age, bmi, smoking)]
  Xconf <- cbind(Intercept = 1, as.matrix(cov_dt[, .(age, bmi, smoking)]))
  Xconf[is.na(Xconf)] <- 0
  folds <- make_folds(y_bin[keep], k = 5, seed = 7)
  prot_names <- colnames(Y)
  metrics <- list()
  preds_all <- rep(NA_real_, length(keep))
  coef_list <- vector("list", length(folds))
  for (fi in seq_along(folds)) {
    log_info("Outer CV fold {fi}/{length(folds)}: preparing train/test sets")
    te_idx <- folds[[fi]]
    tr_idx <- setdiff(seq_along(keep), te_idx)
    # Log class balance for train/test
    y_tr_dbg <- y_bin[keep][tr_idx]; y_te_dbg <- y_bin[keep][te_idx]
    log_info("Fold {fi}: train N={length(y_tr_dbg)} (male={sum(y_tr_dbg==0,na.rm=TRUE)}, female={sum(y_tr_dbg==1,na.rm=TRUE)}); test N={length(y_te_dbg)} (male={sum(y_te_dbg==0,na.rm=TRUE)}, female={sum(y_te_dbg==1,na.rm=TRUE)})")
    Y_tr_raw <- Y[tr_idx, , drop = FALSE]
    Y_te_raw <- Y[te_idx, , drop = FALSE]
    log_info("Outer CV fold {fi}: imputing and residualizing ({ncol(Y_tr_raw)} proteins)")
    col_means <- compute_impute_means(Y_tr_raw)
    Y_tr_imp <- impute_matrix_with_means(Y_tr_raw, col_means)
    Y_te_imp <- impute_matrix_with_means(Y_te_raw, col_means)
    res_tr <- residualize_by_confounders(Y_tr_imp, Xconf[tr_idx, , drop = FALSE])
    R_tr <- res_tr$residuals
    coef_B <- res_tr$coef
    R_te <- apply_residualization(Y_te_imp, Xconf[te_idx, , drop = FALSE], coef_B)
    log_info("Outer CV fold {fi}: scaling residualized matrices and confounders")
    mu <- matrix(colMeans(R_tr, na.rm = TRUE), nrow = 1)
    sdv <- matrix(apply(R_tr, 2, sd, na.rm = TRUE), nrow = 1)
    sdv[sdv == 0 | is.na(sdv)] <- 1
    Z_tr <- sweep(sweep(R_tr, 2, mu, "-"), 2, sdv, "/")
    Z_te <- sweep(sweep(R_te, 2, mu, "-"), 2, sdv, "/")
    conf_tr <- Xconf[tr_idx, -1, drop = FALSE]
    conf_te <- Xconf[te_idx, -1, drop = FALSE]
    # Scale confounders within training fold to mean 0, sd 1
    conf_mu <- matrix(colMeans(conf_tr, na.rm = TRUE), nrow = 1)
    conf_sd <- matrix(apply(conf_tr, 2, sd, na.rm = TRUE), nrow = 1)
    conf_sd[conf_sd == 0 | is.na(conf_sd)] <- 1
    conf_tr_sc <- sweep(sweep(conf_tr, 2, conf_mu, "-"), 2, conf_sd, "/")
    conf_te_sc <- sweep(sweep(conf_te, 2, conf_mu, "-"), 2, conf_sd, "/")
    X_tr <- cbind(Z_tr, conf_tr_sc)
    X_te <- cbind(Z_te, conf_te_sc)
    y_tr <- y_bin[keep][tr_idx]
    y_te <- y_bin[keep][te_idx]
    # Penalize confounders as well (extend shrinkage to confounders)
    pen <- c(rep(1, ncol(Z_tr)), rep(1, ncol(conf_tr)))
    w_tr <- rep(1, length(y_tr))
    p1 <- mean(y_tr == 1); p0 <- mean(y_tr == 0)
    w_tr[y_tr == 1] <- ifelse(p1 > 0, 0.5/p1, 1)
    w_tr[y_tr == 0] <- ifelse(p0 > 0, 0.5/p0, 1)
    # Stronger tuning grid for alpha
    alphas <- seq(0, 1, by = 0.1)
    best_auc <- -Inf; best_alpha <- NA; best_lambda <- NA
    set.seed(13+fi)
    inner_foldid <- sample(rep(1:5, length.out = length(y_tr)))
    log_info("Outer CV fold {fi}: inner CV over alpha grid ({length(alphas)} values)")
    for (a in alphas) {
      fit_cv <- try(cv.glmnet(X_tr, y_tr, family = "binomial", alpha = a, foldid = inner_foldid,
                              weights = w_tr, type.measure = "auc", penalty.factor = pen, standardize = FALSE), silent = TRUE)
      if (inherits(fit_cv, "try-error")) next
      auc <- max(fit_cv$cvm)
      lam <- fit_cv$lambda[which.max(fit_cv$cvm)]
      if (auc > best_auc) { best_auc <- auc; best_alpha <- a; best_lambda <- lam }
    }
    log_info("Outer CV fold {fi}: best inner AUC={round(best_auc,4)} at alpha={best_alpha}")
    if (is.na(best_alpha)) { best_alpha <- 0.5; best_lambda <- NULL }
    fit <- glmnet(X_tr, y_tr, family = "binomial", alpha = best_alpha, lambda = best_lambda,
                  weights = w_tr, penalty.factor = pen, standardize = FALSE)
    pr_te <- as.numeric(predict(fit, newx = X_te, type = "response"))
    preds_all[te_idx] <- pr_te
    # Store per-fold protein coefficients for later comparison with logistic regression
    co <- as.matrix(coef(fit))
    rn <- rownames(co)
    # Initialize vector of protein coefs for all proteins
    prot_vec <- setNames(rep(0, length(prot_names)), prot_names)
    # Intersect available coefficient names with protein names
    common <- intersect(prot_names, rn)
    if (length(common) > 0) {
      prot_vec[common] <- as.numeric(co[common, 1])
    }
    coef_list[[fi]] <- prot_vec
    roc_obj <- try(pROC::roc(response = y_te, predictor = pr_te, quiet = TRUE), silent = TRUE)
    auc <- if (inherits(roc_obj, "try-error")) NA_real_ else as.numeric(pROC::auc(roc_obj))
    pr_obj <- try(PRROC::pr.curve(scores.class0 = pr_te[y_te == 1], scores.class1 = pr_te[y_te == 0], curve = FALSE), silent = TRUE)
    pr_auc <- if (inherits(pr_obj, "try-error")) NA_real_ else pr_obj$auc.integral
    thr <- pick_threshold(pr_te, y_te)
    metrics[[fi]] <- data.table(fold = fi, auc = auc, pr_auc = pr_auc, th = thr$th, J = thr$J, F1 = thr$F1, sens = thr$sens, spec = thr$spec, acc = thr$acc,
                                alpha = best_alpha, lambda = ifelse(is.null(best_lambda), NA_real_, best_lambda))
  }
  metrics_dt <- rbindlist(metrics)
  alpha_pick <- metrics_dt[, .(m_auc = mean(auc, na.rm = TRUE)), by = alpha][order(-m_auc)][1]$alpha

  # Summarize elastic-net coefficients across outer folds (for proteins only)
  if (length(coef_list) > 0) {
    coefs_mat <- do.call(cbind, lapply(coef_list, function(v) v[prot_names]))
    enet_coef_summary <- data.table(
      protein = prot_names,
      enet_beta_mean = rowMeans(coefs_mat, na.rm = TRUE),
      enet_beta_low = apply(coefs_mat, 1, function(x) stats::quantile(x, 0.025, na.rm = TRUE)),
      enet_beta_high = apply(coefs_mat, 1, function(x) stats::quantile(x, 0.975, na.rm = TRUE))
    )
  } else {
    enet_coef_summary <- data.table(protein = character(), enet_beta_mean = numeric(), enet_beta_low = numeric(), enet_beta_high = numeric())
  }
  col_means_all <- compute_impute_means(Y)
  Y_all_imp <- impute_matrix_with_means(Y, col_means_all)
  res_all <- residualize_by_confounders(Y_all_imp, Xconf)
  R_all <- res_all$residuals
  mu_all <- matrix(colMeans(R_all, na.rm = TRUE), nrow = 1)
  sd_all <- matrix(apply(R_all, 2, sd, na.rm = TRUE), nrow = 1)
  sd_all[sd_all == 0 | is.na(sd_all)] <- 1
  Z_all <- sweep(sweep(R_all, 2, mu_all, "-"), 2, sd_all, "/")
  conf_all <- Xconf[, -1, drop = FALSE]
  # Scale confounders on all-data using training-like transform
  conf_mu_all <- matrix(colMeans(conf_all, na.rm = TRUE), nrow = 1)
  conf_sd_all <- matrix(apply(conf_all, 2, sd, na.rm = TRUE), nrow = 1)
  conf_sd_all[conf_sd_all == 0 | is.na(conf_sd_all)] <- 1
  conf_all_sc <- sweep(sweep(conf_all, 2, conf_mu_all, "-"), 2, conf_sd_all, "/")
  log_info("Fitting final glmnet on full data with alpha={alpha_pick}")
  X_all <- cbind(Z_all, conf_all_sc)
  y_all <- y_bin[keep]
  w_all <- rep(1, length(y_all)); p1 <- mean(y_all == 1); p0 <- mean(y_all == 0)
  w_all[y_all == 1] <- ifelse(p1 > 0, 0.5/p1, 1); w_all[y_all == 0] <- ifelse(p0 > 0, 0.5/p0, 1)
  pen_all <- c(rep(1, ncol(Z_all)), rep(1, ncol(conf_all_sc)))
  # Choose lambda via CV on all-data (avoids predicting across the whole lambda path)
  fit_all <- cv.glmnet(X_all, y_all, family = "binomial", alpha = alpha_pick,
                       weights = w_all, penalty.factor = pen_all, standardize = FALSE,
                       type.measure = "auc", nfolds = 5)
  pr_all <- as.numeric(predict(fit_all, newx = X_all, s = "lambda.1se", type = "response"))
  th_global <- pick_threshold(preds_all, y_bin[keep])$th
  pred_label <- ifelse(pr_all >= th_global, 1, 0)
  predicted_sex_str <- ifelse(pred_label == 1, "female", "male")
  genetic_sex_str_all <- ifelse(y_raw[keep] == 2, "female",
                                ifelse(y_raw[keep] == 1, "male", NA_character_))

  preds_dt <- data.table(
    SAMPLE_ID = sample_ids_keep,
    FINNGENID = cov_dt$FINNGENID,
    predicted_prob = pr_all,
    predicted_sex = predicted_sex_str,
    genetic_sex = genetic_sex_str_all
  )
  # Ensure one row per SAMPLE_ID
  if (any(duplicated(preds_dt$SAMPLE_ID))) {
    preds_dt <- preds_dt[!duplicated(SAMPLE_ID)]
  }
  preds_dt[, mismatch := !is.na(genetic_sex) & (genetic_sex != predicted_sex)]

  # Also predict for ALL samples (including excluded/missing-sex) using training transforms
  log_info("Generating predictions for ALL samples (including excluded)")
  # Build full confounder design
  Xconf_full <- cbind(Intercept = 1, as.matrix(pheno_dt[, .(age, bmi, smoking)]))
  Xconf_full[is.na(Xconf_full)] <- 0
  # Impute using training (kept) protein means
  Y_full_imp <- impute_matrix_with_means(npx_matrix, col_means_all)
  # Residualize all samples using coefficients from kept
  R_full <- apply_residualization(Y_full_imp, Xconf_full, res_all$coef)
  # Scale proteins using training stats
  Z_full <- sweep(sweep(R_full, 2, mu_all, "-"), 2, sd_all, "/")
  # Scale confounders using training stats
  conf_full <- Xconf_full[, -1, drop = FALSE]
  conf_full_sc <- sweep(sweep(conf_full, 2, conf_mu_all, "-"), 2, conf_sd_all, "/")
  X_full_pred <- cbind(Z_full, conf_full_sc)
  pr_full_all <- as.numeric(predict(fit_all, newx = X_full_pred, s = "lambda.1se", type = "response"))
  y_raw_full <- pheno_dt$genetic_sex
  genetic_sex_str_full <- ifelse(y_raw_full == 2, "female", ifelse(y_raw_full == 1, "male", NA_character_))

  # ============================================================================
  # THRESHOLD DEFINITIONS (UPDATED LOGIC)
  # ============================================================================
  # Youden's J threshold: optimized decision boundary from ROC curve (th_global)
  # Already computed at line 1142: th_global <- pick_threshold(preds_all, y_bin[keep])$th
  #
  # Two distinct flags (SWAPPED per user request):
  # 1. mismatch: More stringent - deviates from Youden's J threshold (FULL RED)
  # 2. sex_outlier: Mild - deviates from 0.5 but not from Youden's J (PALE RED)
  # ============================================================================

  preds_dt_all <- data.table(
    SAMPLE_ID = pheno_dt$SAMPLE_ID,
    FINNGENID = pheno_dt$FINNGENID,
    predicted_prob = pr_full_all,
    # Predicted sex based on 0.5 threshold (for reference only)
    predicted_sex = ifelse(pr_full_all >= 0.5, "female", "male"),
    genetic_sex = genetic_sex_str_full,
    # Threshold-based flags (for reference and outlier detection):
    # mismatch: using Youden's J threshold (more stringent threshold-based flag)
    mismatch = ( (!is.na(genetic_sex_str_full) & genetic_sex_str_full == "male" & pr_full_all >= th_global) |
                 (!is.na(genetic_sex_str_full) & genetic_sex_str_full == "female" & pr_full_all < th_global) ),
    # Sex outlier: deviates from 0.5 but not from Youden's J (mild threshold-based outlier)
    sex_outlier = ( (!is.na(genetic_sex_str_full) & genetic_sex_str_full == "male" & pr_full_all >= 0.5 & pr_full_all < th_global) )
    # Note: For females when th_global > 0.5, there's no "between" range (prob < th_global are already mismatches)
  )
  # STRICT MISMATCH: predicted_sex != genetic_sex (actual label mismatch)
  # Must be calculated AFTER predicted_sex is created (data.table limitation)
  preds_dt_all[, strict_mismatch := (!is.na(genetic_sex) & !is.na(predicted_sex) &
                                     genetic_sex != predicted_sex)]
  if (any(duplicated(preds_dt_all$SAMPLE_ID))) preds_dt_all <- preds_dt_all[!duplicated(SAMPLE_ID)]
  log_info("Predictions generated for all samples: {nrow(preds_dt_all)}")
  log_info("Strict mismatches (predicted_sex != genetic_sex): {sum(preds_dt_all$strict_mismatch, na.rm=TRUE)}")
  log_info("Threshold-based mismatches (Youden J = {round(th_global, 3)}): {sum(preds_dt_all$mismatch, na.rm=TRUE)}")
  log_info("Sex outliers (0.5 threshold, mild): {sum(preds_dt_all$sex_outlier, na.rm=TRUE)}")
  preds_dt_out <- preds_dt_all

  # Save a distribution plot of predicted probabilities
  plot_with_ids <- tryCatch(isTRUE(config$parameters$outliers$plot_with_ids), error = function(e) FALSE)
  sex_predict_distribution_path <- get_output_path(step_num, "sex_predict_distribution", batch_id, "outliers", "pdf", config = config)
  sex_predict_distribution_with_ids_path <- get_output_path(step_num, "sex_predict_distribution_with_ids", batch_id, "outliers", "pdf", config = config)
  ensure_output_dir(sex_predict_distribution_path)
  ensure_output_dir(sex_predict_distribution_with_ids_path)

  create_prediction_distribution_panels(
    df_pred = preds_dt_out[, .(SAMPLE_ID, FINNGENID, genetic_sex, predicted_prob)],
    title = "Sex Prediction Distribution – Main Pipeline (Elastic Net Logistic Regression)",
    file_out = sex_predict_distribution_path,
    youden_thresh = th_global, annotate_ids = FALSE
  )
  # With FINNGENID annotations on the right panel (only if configured)
  if (isTRUE(plot_with_ids)) {
    create_prediction_distribution_panels(
      df_pred = preds_dt_out[, .(SAMPLE_ID, FINNGENID, genetic_sex, predicted_prob)],
      title = "Sex Prediction Distribution\n(with FINNGENID labels) – Main Pipeline (Elastic Net Logistic Regression)",
      file_out = sex_predict_distribution_with_ids_path,
      youden_thresh = th_global, annotate_ids = TRUE
    )
  }

  # Create visualizations
  if (nrow(sex_proteins$associations) > 0) {
    volcano_plot <- create_volcano_plot(sex_proteins$associations)
    volcano_plot_or <- create_volcano_plot_or(sex_proteins$associations)
    forest_plot <- create_forest_plot(sex_proteins$associations)
    hist_facets <- create_top_hist_facets(npx_matrix, sex_proteins$associations, sex_proteins$sample_sex)
    scatter_plot <- create_effect_scatter_top50(sex_proteins$associations, enet_coef_summary)
    heatmap_obj <- create_corr_heatmap(npx_matrix, sex_proteins$associations, id_map, sex_info, top_n = 50)
  } else {
    volcano_plot <- ggplot() + theme_void() + ggtitle("No associations")
    volcano_plot_or <- volcano_plot
    forest_plot <- NULL
    hist_facets <- NULL
    scatter_plot <- NULL
  }

  # Create smoking-associated proteins volcano/forest plots
  if (nrow(smoking_proteins$associations) > 0) {
    smoking_volcano_plot <- create_smoking_volcano_plot(smoking_proteins$associations)
    smoking_forest_plot <- create_smoking_forest_plot(smoking_proteins$associations)
  } else {
    smoking_volcano_plot <- ggplot() + theme_void() + ggtitle("No smoking associations")
    smoking_forest_plot <- NULL
  }

  # Save outputs with batch-aware paths
  log_info("Saving sex outlier detection results")

  sex_proteins_path <- get_output_path(step_num, "sex_proteins", batch_id, "outliers", config = config)
  smoking_proteins_path <- get_output_path(step_num, "smoking_proteins", batch_id, "outliers", config = config)
  sex_model_cv_metrics_path <- get_output_path(step_num, "sex_model_cv_metrics", batch_id, "outliers", "tsv", config = config)
  sex_predictions_path <- get_output_path(step_num, "sex_predictions", batch_id, "outliers", "tsv", config = config)
  sex_mismatches_path <- get_output_path(step_num, "sex_mismatches", batch_id, "outliers", "tsv", config = config)
  sex_outliers_path <- get_output_path(step_num, "sex_outliers", batch_id, "outliers", "tsv", config = config)

  ensure_output_dir(sex_proteins_path)
  ensure_output_dir(smoking_proteins_path)
  ensure_output_dir(sex_model_cv_metrics_path)
  ensure_output_dir(sex_predictions_path)
  ensure_output_dir(sex_mismatches_path)
  ensure_output_dir(sex_outliers_path)

  saveRDS(sex_proteins, sex_proteins_path)
  saveRDS(smoking_proteins, smoking_proteins_path)
  fwrite(metrics_dt, sex_model_cv_metrics_path, sep = "\t")
  fwrite(preds_dt_out, sex_predictions_path, sep = "\t")

  # Add BIOBANK_PLASMA and DISEASE_GROUP before writing mismatches/outliers
  preds_dt_out_with_meta <- add_metadata_to_predictions(preds_dt_out, metadata)

  # ============================================================================
  # CREATE TWO SEPARATE LISTS:
  # 1. Strict mismatches: predicted_sex != genetic_sex (actual label mismatch)
  # 2. Outliers: Threshold-based outliers that are NOT strict mismatches
  # ============================================================================

  # List 1: Strict mismatches (predicted_sex != genetic_sex)
  strict_mismatches <- preds_dt_out_with_meta[strict_mismatch == TRUE]
  fwrite(strict_mismatches, sex_mismatches_path, sep = "\t")
  log_info("Saved {nrow(strict_mismatches)} strict mismatches (predicted_sex != genetic_sex) to: {sex_mismatches_path}")

  # List 2: Threshold-based outliers that are NOT strict mismatches
  # These are samples that cross thresholds but predicted_sex == genetic_sex
  threshold_outliers <- preds_dt_out_with_meta[
    (mismatch == TRUE | isTRUE(sex_outlier)) &
    strict_mismatch == FALSE
  ]
  fwrite(threshold_outliers, sex_outliers_path, sep = "\t")
  log_info("Saved {nrow(threshold_outliers)} threshold-based outliers (not strict mismatches) to: {sex_outliers_path}")

  # Summary statistics
  log_info("=== SEX DETECTION SUMMARY ===")
  log_info("Strict mismatches (predicted_sex != genetic_sex): {nrow(strict_mismatches)}")
  log_info("Threshold-based outliers (not strict mismatches): {nrow(threshold_outliers)}")
  log_info("Total flagged samples: {nrow(strict_mismatches) + nrow(threshold_outliers)}")
  # Also write a unified fold-level metrics file for the main model (schema-compatible with enhanced)
  # Construct enhanced_dir path using batch-aware path construction
  base_dir <- config$output$base_dir %||% "output"
  batch_subdir <- if(batch_id != config$batch$default_batch_id ||
    tryCatch(isTRUE(config$parameters$normalization$multi_batch_mode), error = function(e) FALSE))
    batch_id else ""
  enhanced_dir <- file.path(base_dir, "outliers", batch_subdir, "enhanced_models")
  dir.create(enhanced_dir, recursive = TRUE, showWarnings = FALSE)
  # Ensure all columns exist
  required_cols <- c("fold","auc","pr_auc","th","J","F1","sens","spec","acc","alpha","lambda")
  out_metrics <- copy(metrics_dt)
  for (cc in required_cols) if (!(cc %in% names(out_metrics))) out_metrics[, (cc) := NA_real_]
  out_metrics <- out_metrics[, ..required_cols]
  fwrite(out_metrics, file.path(enhanced_dir, "all_models_auc.tsv"), sep = "\t")

  # Distribution of sex outliers by plate (if plate information is available)
  # We try to load samples_data_raw.rds to get PlateID per SampleID (batch-aware path)
  plate_map_file <- get_output_path(prev_step00_num, "samples_data_raw", batch_id, "qc", config = config)
  if (file.exists(plate_map_file)) {
    log_info("Computing sex outlier distribution by PlateID")
    samples_data_raw <- try(readRDS(plate_map_file), silent = TRUE)
    if (!inherits(samples_data_raw, "try-error")) {
      # Expecting columns SampleID, PlateID
      plate_dt <- unique(samples_data_raw[, .(SampleID, PlateID)])
      # Use sex_outlier flag to match red points on right panel (consistent with threshold-based outliers)
      use_dt <- copy(preds_dt_out)
      if (!("sex_outlier" %in% names(use_dt))) use_dt[, sex_outlier := NA]
      plate_merge <- merge(use_dt[, .(SampleID = SAMPLE_ID, sex_outlier)], plate_dt, by = "SampleID", all.x = TRUE)
      plate_merge[, sex_outlier := as.logical(sex_outlier)]
      plate_summary <- plate_merge[, .(
        sexOutRate = mean(sex_outlier, na.rm = TRUE),
        Noutlier = sum(sex_outlier, na.rm = TRUE),
        N = .N
      ), by = PlateID][order(-sexOutRate)]
      fwrite(plate_summary, , sep = "\t")
    } else {
      log_warn("Failed to read samples_data_raw.rds; skipping sex outlier by plate summary")
    }
  } else {
    log_warn("samples_data_raw.rds not found; skipping sex outlier by plate summary")
  }

  # Save tables
  sex_associated_proteins_path <- get_output_path(step_num, "sex_associated_proteins", batch_id, "outliers", "tsv", config = config)
  smoking_associated_proteins_path <- get_output_path(step_num, "smoking_associated_proteins", batch_id, "outliers", "tsv", config = config)
  ensure_output_dir(sex_associated_proteins_path)
  ensure_output_dir(smoking_associated_proteins_path)
  fwrite(sex_proteins$associations, sex_associated_proteins_path, sep = "\t")
  fwrite(smoking_proteins$associations, smoking_associated_proteins_path, sep = "\t")



  # Save plots with batch-aware paths
  sex_proteins_volcano_path <- get_output_path(step_num, "sex_proteins_volcano", batch_id, "outliers", "pdf", config = config)
  sex_proteins_volcano_or_path <- get_output_path(step_num, "sex_proteins_volcano_or", batch_id, "outliers", "pdf", config = config)
  sex_forest_significant_path <- get_output_path(step_num, "sex_forest_significant", batch_id, "outliers", "pdf", config = config)
  sex_top10_histograms_path <- get_output_path(step_num, "sex_top10_histograms", batch_id, "outliers", "pdf", config = config)
  sex_effects_scatter_top50_path <- get_output_path(step_num, "sex_effects_scatter_top50", batch_id, "outliers", "pdf", config = config)
  sex_corr_heatmap_path <- get_output_path(step_num, "sex_corr_heatmap", batch_id, "outliers", "pdf", config = config)
  smoking_associated_proteins_volcano_path <- get_output_path(step_num, "smoking_associated_proteins_volcano", batch_id, "outliers", "pdf", config = config)
  smoking_associated_proteins_forest_path <- get_output_path(step_num, "smoking_associated_proteins_forest", batch_id, "outliers", "pdf", config = config)

  ensure_output_dir(sex_proteins_volcano_path)
  ensure_output_dir(sex_proteins_volcano_or_path)
  ensure_output_dir(sex_forest_significant_path)
  ensure_output_dir(sex_top10_histograms_path)
  ensure_output_dir(sex_effects_scatter_top50_path)
  ensure_output_dir(sex_corr_heatmap_path)
  ensure_output_dir(smoking_associated_proteins_volcano_path)
  ensure_output_dir(smoking_associated_proteins_forest_path)

  ggsave(sex_proteins_volcano_path, volcano_plot, width = 10, height = 8)
  ggsave(sex_proteins_volcano_or_path, volcano_plot_or, width = 10, height = 8)

  if (!is.null(forest_plot)) ggsave(sex_forest_significant_path, forest_plot, width = 9, height = 12)
  if (!is.null(hist_facets)) ggsave(sex_top10_histograms_path, hist_facets, width = 12, height = 10)
  if (!is.null(scatter_plot)) ggsave(sex_effects_scatter_top50_path, scatter_plot, width = 10, height = 8)
  if (exists("heatmap_obj")) {
    pdf(sex_corr_heatmap_path, width = 10, height = 10)
    print(heatmap_obj)
    dev.off()
  }
  # Save smoking volcano plot
  ggsave(smoking_associated_proteins_volcano_path, smoking_volcano_plot, width = 10, height = 8)
  if (!is.null(smoking_forest_plot)) {
    ggsave(smoking_associated_proteins_forest_path, smoking_forest_plot, width = 9, height = 12)
  }

  # Optional: Enhanced modeling (classical + tree-based + NN if available) under config toggles
  run_nn <- tryCatch(isTRUE(config$parameters$outliers$run_nn), error = function(e) FALSE)
  run_ae <- tryCatch(isTRUE(config$parameters$outliers$run_autoencoder), error = function(e) FALSE)
  run_enhanced <- tryCatch(isTRUE(config$parameters$outliers$run_enhanced), error = function(e) FALSE) || run_nn || run_ae
  top_n_nn <- tryCatch(as.integer(config$parameters$outliers$top_n_proteins_nn), error = function(e) 150)
  if (!is.finite(top_n_nn) || top_n_nn <= 0) top_n_nn <- 150

  if (run_enhanced) {
    log_info("Running enhanced models block (classical + tree-based; NN if available)")
    # Prepare CV folds on the same keep-set
    y_cv <- y_bin[keep]
    set.seed(42)
    folds_cv <- if (has_caret) caret::createFolds(as.factor(y_cv), k = 5, list = TRUE, returnTrain = FALSE) else split(seq_along(y_cv), rep(1:5, length.out = length(y_cv)))
    all_preds <- list()
    metrics_by_model <- list()  # per-fold metrics, same columns as main
    models <- c("enet","ridge","lasso","logl2", if (has_rf) "rf" else NULL, if (has_xgb) "xgb" else NULL, if (has_keras && run_nn) "mlp" else NULL, if (has_keras && run_ae) "ae" else NULL)
    models <- unlist(models)
    for (mi in models) all_preds[[mi]] <- rep(NA_real_, length(y_cv))
    for (mi in models) metrics_by_model[[mi]] <- data.table()

    for (fi in seq_along(folds_cv)) {
      te <- folds_cv[[fi]]; tr <- setdiff(seq_along(y_cv), te)
      X_tr_raw <- Y[tr, , drop = FALSE]
      X_te_raw <- Y[te, , drop = FALSE]
      # Within-fold feature selection on raw (pre-residual) protein matrix already imputed above in pipeline
      sel <- select_features_ensemble(X_tr_raw, y_cv[tr], top_n = top_n_nn)$selected
      X_tr <- X_tr_raw[, sel, drop = FALSE]
      X_te <- X_te_raw[, sel, drop = FALSE]

      # Elastic-net (ridge/lasso variants)
      sc <- standardize_features(X_tr, X_te)
      # ENET alpha grid
      best_auc <- -Inf; best_alpha <- 0.5; best_model <- NULL
      for (a in seq(0, 1, by = 0.2)) {
        cvm <- try(cv.glmnet(sc$train, y_cv[tr], family = "binomial", alpha = a, type.measure = "auc", nfolds = 3), silent = TRUE)
        if (inherits(cvm, "try-error")) next
        auc_cand <- max(cvm$cvm, na.rm = TRUE)
        if (is.finite(auc_cand) && auc_cand > best_auc) { best_auc <- auc_cand; best_alpha <- a; best_model <- cvm }
      }
      if (!is.null(best_model)) {
        pr <- as.numeric(predict(best_model, sc$test, s = "lambda.1se", type = "response"))
        all_preds[["enet"]][te] <- pr
        thr <- pick_threshold(pr, y_cv[te])
        metrics_by_model[["enet"]] <- rbind(metrics_by_model[["enet"]],
          data.table(fold = fi,
                     auc = tryCatch(as.numeric(pROC::auc(pROC::roc(y_cv[te], pr, quiet = TRUE))), error = function(e) NA_real_),
                     pr_auc = tryCatch(PRROC::pr.curve(scores.class0 = pr[y_cv[te] == 1], scores.class1 = pr[y_cv[te] == 0], curve = FALSE)$auc.integral, error = function(e) NA_real_),
                     th = thr$th, J = thr$J, F1 = thr$F1, sens = thr$sens, spec = thr$spec, acc = thr$acc,
                     alpha = best_alpha,
                     lambda = tryCatch(best_model$lambda.1se, error = function(e) NA_real_)))
      }
      # Ridge
      rid <- try(cv.glmnet(sc$train, y_cv[tr], family = "binomial", alpha = 0, type.measure = "auc", nfolds = 5), silent = TRUE)
      if (!inherits(rid, "try-error")) {
        pr_r <- as.numeric(predict(rid, sc$test, s = "lambda.1se", type = "response"))
        all_preds[["ridge"]][te] <- pr_r
        thr <- pick_threshold(pr_r, y_cv[te])
        metrics_by_model[["ridge"]] <- rbind(metrics_by_model[["ridge"]],
          data.table(fold = fi,
                     auc = tryCatch(as.numeric(pROC::auc(pROC::roc(y_cv[te], pr_r, quiet = TRUE))), error = function(e) NA_real_),
                     pr_auc = tryCatch(PRROC::pr.curve(scores.class0 = pr_r[y_cv[te] == 1], scores.class1 = pr_r[y_cv[te] == 0], curve = FALSE)$auc.integral, error = function(e) NA_real_),
                     th = thr$th, J = thr$J, F1 = thr$F1, sens = thr$sens, spec = thr$spec, acc = thr$acc,
                     alpha = 0, lambda = tryCatch(rid$lambda.1se, error = function(e) NA_real_)))
      }
      # LASSO
      las <- try(cv.glmnet(sc$train, y_cv[tr], family = "binomial", alpha = 1, type.measure = "auc", nfolds = 5), silent = TRUE)
      if (!inherits(las, "try-error")) {
        pr_l <- as.numeric(predict(las, sc$test, s = "lambda.1se", type = "response"))
        all_preds[["lasso"]][te] <- pr_l
        thr <- pick_threshold(pr_l, y_cv[te])
        metrics_by_model[["lasso"]] <- rbind(metrics_by_model[["lasso"]],
          data.table(fold = fi,
                     auc = tryCatch(as.numeric(pROC::auc(pROC::roc(y_cv[te], pr_l, quiet = TRUE))), error = function(e) NA_real_),
                     pr_auc = tryCatch(PRROC::pr.curve(scores.class0 = pr_l[y_cv[te] == 1], scores.class1 = pr_l[y_cv[te] == 0], curve = FALSE)$auc.integral, error = function(e) NA_real_),
                     th = thr$th, J = thr$J, F1 = thr$F1, sens = thr$sens, spec = thr$spec, acc = thr$acc,
                     alpha = 1, lambda = tryCatch(las$lambda.1se, error = function(e) NA_real_)))
      }
      # Logistic L2 (same as ridge but explicit)
      if (!inherits(rid, "try-error")) {
        pr_lr <- as.numeric(predict(rid, sc$test, s = "lambda.1se", type = "response"))
        all_preds[["logl2"]][te] <- pr_lr
        thr <- pick_threshold(pr_lr, y_cv[te])
        metrics_by_model[["logl2"]] <- rbind(metrics_by_model[["logl2"]],
          data.table(fold = fi,
                     auc = tryCatch(as.numeric(pROC::auc(pROC::roc(y_cv[te], pr_lr, quiet = TRUE))), error = function(e) NA_real_),
                     pr_auc = tryCatch(PRROC::pr.curve(scores.class0 = pr_lr[y_cv[te] == 1], scores.class1 = pr_lr[y_cv[te] == 0], curve = FALSE)$auc.integral, error = function(e) NA_real_),
                     th = thr$th, J = thr$J, F1 = thr$F1, sens = thr$sens, spec = thr$spec, acc = thr$acc,
                     alpha = 0, lambda = tryCatch(rid$lambda.1se, error = function(e) NA_real_)))
      }
      # RF
      if (has_rf) {
        rf_fit <- try(randomForest::randomForest(X_tr, as.factor(y_cv[tr]), ntree = 300, mtry = max(1, floor(sqrt(ncol(X_tr))))), silent = TRUE)
        if (!inherits(rf_fit, "try-error")) {
          pr_rf <- tryCatch(as.numeric(predict(rf_fit, X_te, type = "prob")[, 2]), error = function(e) NA_real_)
          all_preds[["rf"]][te] <- pr_rf
          thr <- pick_threshold(pr_rf, y_cv[te])
          metrics_by_model[["rf"]] <- rbind(metrics_by_model[["rf"]],
            data.table(fold = fi,
                       auc = tryCatch(as.numeric(pROC::auc(pROC::roc(y_cv[te], pr_rf, quiet = TRUE))), error = function(e) NA_real_),
                       pr_auc = tryCatch(PRROC::pr.curve(scores.class0 = pr_rf[y_cv[te] == 1], scores.class1 = pr_rf[y_cv[te] == 0], curve = FALSE)$auc.integral, error = function(e) NA_real_),
                       th = thr$th, J = thr$J, F1 = thr$F1, sens = thr$sens, spec = thr$spec, acc = thr$acc,
                       alpha = NA_real_, lambda = NA_real_))
        }
      }
      # XGB
      if (has_xgb) {
        dtr <- xgboost::xgb.DMatrix(data = as.matrix(X_tr), label = y_cv[tr])
        dte <- xgboost::xgb.DMatrix(data = as.matrix(X_te), label = y_cv[te])
        spw <- sum(y_cv[tr] == 0) / max(1, sum(y_cv[tr] == 1))
        params <- list(objective = "binary:logistic", eval_metric = "auc", max_depth = 4, eta = 0.05, subsample = 0.7, colsample_bytree = 0.7, min_child_weight = 5, gamma = 1, alpha = 1, lambda = 1, scale_pos_weight = spw)
        xgb <- try(xgboost::xgb.train(params = params, data = dtr, nrounds = 100, watchlist = list(train = dtr), verbose = 0), silent = TRUE)
        if (!inherits(xgb, "try-error")) {
          pr_xgb <- tryCatch(as.numeric(predict(xgb, dte)), error = function(e) NA_real_)
          all_preds[["xgb"]][te] <- pr_xgb
          thr <- pick_threshold(pr_xgb, y_cv[te])
          metrics_by_model[["xgb"]] <- rbind(metrics_by_model[["xgb"]],
            data.table(fold = fi,
                       auc = tryCatch(as.numeric(pROC::auc(pROC::roc(y_cv[te], pr_xgb, quiet = TRUE))), error = function(e) NA_real_),
                       pr_auc = tryCatch(PRROC::pr.curve(scores.class0 = pr_xgb[y_cv[te] == 1], scores.class1 = pr_xgb[y_cv[te] == 0], curve = FALSE)$auc.integral, error = function(e) NA_real_),
                       th = thr$th, J = thr$J, F1 = thr$F1, sens = thr$sens, spec = thr$spec, acc = thr$acc,
                       alpha = NA_real_, lambda = NA_real_))
        }
      }
      # MLP / AE (only if keras available and toggled)
      if (has_keras && run_nn) {
        # Use whichever namespace is available
        K <- if (requireNamespace("keras3", quietly = TRUE)) getNamespace("keras3") else getNamespace("keras")
        sc_nn <- standardize_features(X_tr, X_te)
        model <- try({
          m <- K$keras_model_sequential()
          m$add(K$layer_dense(units = 128L, input_shape = as.integer(ncol(sc_nn$train)), activation = "relu"))
          m$add(K$layer_dropout(rate = 0.3))
          m$add(K$layer_dense(units = 64L, activation = "relu"))
          m$add(K$layer_dropout(rate = 0.3))
          m$add(K$layer_dense(units = 1L, activation = "sigmoid"))
          m$compile(optimizer = K$optimizer_adam(learning_rate = 0.001), loss = "binary_crossentropy")
          m$fit(as.matrix(sc_nn$train), y_cv[tr], epochs = 30L, batch_size = 64L, validation_split = 0.1, verbose = 0L)
          m
        }, silent = TRUE)
        if (!inherits(model, "try-error")) {
          pr_mlp <- tryCatch(as.numeric(model$predict(as.matrix(sc_nn$test))), error = function(e) NA_real_)
          all_preds[["mlp"]][te] <- pr_mlp
          thr <- pick_threshold(pr_mlp, y_cv[te])
          metrics_by_model[["mlp"]] <- rbind(metrics_by_model[["mlp"]],
            data.table(fold = fi,
                       auc = tryCatch(as.numeric(pROC::auc(pROC::roc(y_cv[te], pr_mlp, quiet = TRUE))), error = function(e) NA_real_),
                       pr_auc = tryCatch(PRROC::pr.curve(scores.class0 = pr_mlp[y_cv[te] == 1], scores.class1 = pr_mlp[y_cv[te] == 0], curve = FALSE)$auc.integral, error = function(e) NA_real_),
                       th = thr$th, J = thr$J, F1 = thr$F1, sens = thr$sens, spec = thr$spec, acc = thr$acc,
                       alpha = NA_real_, lambda = NA_real_))
        }
      }
      if (has_keras && run_ae) {
        # Skip heavy AE here to keep runtime moderate; MLP suffices as NN baseline
      }
    }
    # Evaluate and save comparison
    # Build ENSEMBLE CV prediction (average of available model CV preds)
    avail_cv <- all_preds[!sapply(all_preds, function(v) all(is.na(v)))]
    ensemble_cv <- if (length(avail_cv) > 0) rowMeans(do.call(cbind, avail_cv), na.rm = TRUE) else rep(NA_real_, length(y_cv))
    # Metrics per model + ENSEMBLE (aggregate for winner selection)
    # Only evaluate models that have valid predictions (not all NA)
    valid_models <- names(all_preds)[!sapply(all_preds, function(v) all(is.na(v)))]
    if (length(valid_models) == 0) {
      log_warn("No models produced valid predictions - all models failed to train")
      metrics_models <- data.table(model = character(), auc = numeric(), pr_auc = numeric())
    } else {
      metrics_models <- rbindlist(lapply(valid_models, function(nm) {
        pr <- all_preds[[nm]]
        # Additional check: ensure predictions are not constant (would cause ROC to fail)
        if (length(unique(pr[!is.na(pr)])) < 2) {
          log_warn("Model {nm} has constant predictions - skipping AUC evaluation")
          return(data.table(model = toupper(nm), auc = NA_real_, pr_auc = NA_real_))
        }
        ev <- evaluate_predictions_quick(y_cv, pr)
        if (is.na(ev$auc)) {
          log_warn("Model {nm} AUC evaluation failed - predictions may be invalid")
        }
        data.table(model = toupper(nm), auc = ev$auc, pr_auc = ev$pr_auc)
      }), use.names = TRUE, fill = TRUE)
    }
    # Only compute ensemble metrics if ensemble has valid predictions
    if (all(is.na(ensemble_cv)) || length(unique(ensemble_cv[!is.na(ensemble_cv)])) < 2) {
      log_warn("Ensemble predictions are invalid (all NA or constant) - skipping ensemble AUC evaluation")
      metrics_ensemble <- data.table(model = "ENSEMBLE", auc = NA_real_, pr_auc = NA_real_)
    } else {
      ev_ensemble <- evaluate_predictions_quick(y_cv, ensemble_cv)
      if (is.na(ev_ensemble$auc)) {
        log_warn("Ensemble AUC evaluation failed - predictions may be invalid")
      }
      metrics_ensemble <- data.table(model = "ENSEMBLE", auc = ev_ensemble$auc, pr_auc = ev_ensemble$pr_auc)
    }
    metrics <- rbind(metrics_models, metrics_ensemble, fill = TRUE)
    final_metrics <- metrics
    # Ensure enhanced_models directory exists before any writes (batch-aware path)
    base_dir <- config$output$base_dir %||% "/mnt/longGWAS_disk_100GB/long_gwas/Github_clones/fg3_Olink_analysis/output"
    batch_subdir <- if(batch_id != config$batch$default_batch_id ||
      tryCatch(isTRUE(config$parameters$normalization$multi_batch_mode), error = function(e) FALSE))
      batch_id else ""
    enhanced_dir <- file.path(base_dir, "outliers", batch_subdir, "enhanced_models")
    dir.create(enhanced_dir, recursive = TRUE, showWarnings = FALSE)
    # Build per-fold metrics table for ALL models with the same schema as the main metrics file
    required_cols <- c("fold","auc","pr_auc","th","J","F1","sens","spec","acc","alpha","lambda")
    # Add ensemble per-fold metrics using fold-wise thresholds on test preds
    metrics_by_model[["ensemble"]] <- data.table()
    for (fi in seq_along(folds_cv)) {
      te <- folds_cv[[fi]]
      pr_e <- ensemble_cv[te]
      thr <- pick_threshold(pr_e, y_cv[te])
      metrics_by_model[["ensemble"]] <- rbind(metrics_by_model[["ensemble"]],
        data.table(fold = fi,
                   auc = tryCatch(as.numeric(pROC::auc(pROC::roc(y_cv[te], pr_e, quiet = TRUE))), error = function(e) NA_real_),
                   pr_auc = tryCatch(PRROC::pr.curve(scores.class0 = pr_e[y_cv[te] == 1], scores.class1 = pr_e[y_cv[te] == 0], curve = FALSE)$auc.integral, error = function(e) NA_real_),
                   th = thr$th, J = thr$J, F1 = thr$F1, sens = thr$sens, spec = thr$spec, acc = thr$acc,
                   alpha = NA_real_, lambda = NA_real_))
    }
    metrics_all_models <- rbindlist(lapply(names(metrics_by_model), function(nm) {
      dtm <- copy(metrics_by_model[[nm]])
      for (cc in required_cols) if (!(cc %in% names(dtm))) dtm[, (cc) := NA_real_]
      cbind(data.table(model = toupper(nm)), dtm[, ..required_cols])
    }), use.names = TRUE, fill = TRUE)
    fwrite(metrics_all_models, file.path(enhanced_dir, "all_models_auc.tsv"), sep = "\t")
    # Save per-model distribution plots
    sex_lbl <- ifelse(y_cv == 1, "female", "male")
    for (nm in names(all_preds)) {
      pd <- data.table(SAMPLE_ID = sample_ids_keep, genetic_sex = sex_lbl, predicted_prob = all_preds[[nm]])
      # Compute Youden's J threshold for this model from CV predictions
      thr_nm <- pick_threshold(all_preds[[nm]], y_cv)$th
      create_prediction_distribution_panels(
        df_pred = pd,
        title = paste0(model_display_name(nm), " – Sex Prediction Distribution"),
        file_out = file.path(enhanced_dir, paste0("sex_predict_dist_", nm, ".pdf")),
        youden_thresh = thr_nm, annotate_ids = FALSE
      )
      if (isTRUE(plot_with_ids)) {
        create_prediction_distribution_panels(
          df_pred = pd[, .(SAMPLE_ID, FINNGENID = NA_character_, genetic_sex, predicted_prob)],
          title = paste0(model_display_name(nm), " – Sex Prediction Distribution\n(with FINNGENID labels)"),
          file_out = file.path(enhanced_dir, paste0("sex_predict_dist_", nm, "_with_ids.pdf")),
          youden_thresh = thr_nm, annotate_ids = TRUE
        )
      }
    }
    # Also produce ENSEMBLE plot from ensemble_cv
    # Compose descriptive ensemble member list
    ensemble_member_names <- character(0)
    if (length(avail_cv) > 0) {
      member_codes <- names(avail_cv)
      ensemble_member_names <- vapply(member_codes, function(c) model_display_name(c), character(1))
    }
    # Compute Youden's J threshold for ensemble from CV predictions
    thr_en <- pick_threshold(ensemble_cv, y_cv)$th
    create_prediction_distribution_panels(
      df_pred = data.table(SAMPLE_ID = sample_ids_keep, genetic_sex = sex_lbl, predicted_prob = ensemble_cv),
      title = paste0(model_display_name("ensemble", ensemble_member_names), " – Sex Prediction Distribution"),
      file_out = file.path(enhanced_dir, "sex_predict_dist_ensemble.pdf"),
      youden_thresh = thr_en, annotate_ids = FALSE
    )
    create_prediction_distribution_panels(
      df_pred = data.table(SAMPLE_ID = sample_ids_keep, FINNGENID = NA_character_, genetic_sex = sex_lbl, predicted_prob = ensemble_cv),
      title = paste0(model_display_name("ensemble", ensemble_member_names), " – Sex Prediction Distribution (with FINNGENID labels)"),
      file_out = file.path(enhanced_dir, "sex_predict_dist_ensemble_with_ids.pdf"),
      youden_thresh = thr_en, annotate_ids = TRUE
    )

    # =============================
    # Select best model and align outputs
    # =============================
    # Compare ENET main pipeline AUC vs enhanced models
    # ENET_MAIN should always have valid AUC (basic glmnet models should work)
    enet_main_auc <- try(mean(metrics_dt$auc, na.rm = TRUE), silent = TRUE)
    enet_main_auc <- if (inherits(enet_main_auc, "try-error") || !is.finite(enet_main_auc)) {
      log_error("ENET_MAIN AUC calculation failed - this should never happen! Check metrics_dt")
      log_error("  metrics_dt dimensions: {nrow(metrics_dt)} rows")
      log_error("  AUC values: {paste(metrics_dt$auc, collapse=', ')}")
      NA_real_
    } else {
      enet_main_auc
    }

    # Filter out models with NA or invalid AUCs from final_metrics before comparison
    final_metrics_valid <- final_metrics[!is.na(auc) & is.finite(auc)]

    # Handle case where final_metrics might be empty or all AUCs are NA
    if (nrow(final_metrics_valid) == 0) {
      log_warn("No enhanced models with valid AUCs available")
      log_warn("  Total models attempted: {nrow(final_metrics)}")
      log_warn("  Models with valid AUCs: {nrow(final_metrics_valid)}")
      if (nrow(final_metrics) > 0) {
        log_warn("  Model names: {paste(final_metrics$model, collapse=', ')}")
        log_warn("  AUC values: {paste(final_metrics$auc, collapse=', ')}")
      }
      best_enh_auc <- NA_real_
      best_model_name <- "ENET_MAIN"
      winner <- "ENET_MAIN"
      log_info("Defaulting to ENET_MAIN as no enhanced models produced valid AUCs")
    } else {
      best_enh <- final_metrics_valid[which.max(auc)]
      best_model_name <- best_enh$model
      best_enh_auc <- best_enh$auc
      # Decide winner - both AUCs should be finite at this point
      winner <- best_model_name
      if (is.finite(enet_main_auc) && is.finite(best_enh_auc) && enet_main_auc >= best_enh_auc) {
        winner <- "ENET_MAIN"
        log_info("Best model by CV AUC: ENET_MAIN (AUC={round(enet_main_auc,4)})")
      } else if (is.finite(best_enh_auc)) {
        log_info("Best model by CV AUC: {winner} (AUC={round(best_enh_auc,4)})")
      } else {
        log_warn("Enhanced model AUC is invalid despite filtering - defaulting to ENET_MAIN")
        winner <- "ENET_MAIN"
      }
    }

    # Helper to compute Youden J threshold from CV predictions
    compute_youden_threshold <- function(y_true, y_pred) {
      thr <- pick_threshold(y_pred, y_true)
      if (is.null(thr) || is.na(thr$th)) return(0.5)
      as.numeric(thr$th)
    }

    # Build full-data predictions for winner
    # For ENET_MAIN we already have pr_all and th_global
    preds_dt_final <- copy(preds_dt)
    if (winner != "ENET_MAIN") {
      # Prepare full-data feature selection on raw proteins
      y_full <- y_bin[keep]
      sel_all <- select_features_ensemble(Y, y_full, top_n = top_n_nn)$selected
      X_full <- Y[, sel_all, drop = FALSE]
      # Standardize on all-data (for final fit only)
      sc_all <- standardize_features(X_full)
      X_full_sc <- as.matrix(sc_all$train)

      # Fit per model on full-data
      pred_full <- rep(NA_real_, length(y_full))
      if (winner == "RIDGE" || winner == "LOGL2") {
        cvm <- try(cv.glmnet(X_full_sc, y_full, family = "binomial", alpha = 0, type.measure = "auc", nfolds = 5), silent = TRUE)
        if (!inherits(cvm, "try-error")) pred_full <- as.numeric(predict(cvm, X_full_sc, s = "lambda.1se", type = "response"))
      } else if (winner == "LASSO") {
        cvm <- try(cv.glmnet(X_full_sc, y_full, family = "binomial", alpha = 1, type.measure = "auc", nfolds = 5), silent = TRUE)
        if (!inherits(cvm, "try-error")) pred_full <- as.numeric(predict(cvm, X_full_sc, s = "lambda.1se", type = "response"))
      } else if (winner == "ENET") {
        # Use alpha grid briefly to approximate best alpha on full data
        best_auc2 <- -Inf; best_alpha2 <- 0.5; best_mod2 <- NULL
        for (a in seq(0, 1, by = 0.2)) {
          m <- try(cv.glmnet(X_full_sc, y_full, family = "binomial", alpha = a, type.measure = "auc", nfolds = 5), silent = TRUE)
          if (inherits(m, "try-error")) next
          auc_c <- max(m$cvm, na.rm = TRUE)
          if (is.finite(auc_c) && auc_c > best_auc2) { best_auc2 <- auc_c; best_alpha2 <- a; best_mod2 <- m }
        }
        if (!is.null(best_mod2)) pred_full <- as.numeric(predict(best_mod2, X_full_sc, s = "lambda.1se", type = "response"))
      } else if (winner == "RF") {
        if (has_rf) {
          rf_fit <- try(randomForest::randomForest(X_full, as.factor(y_full), ntree = 300, mtry = max(1, floor(sqrt(ncol(X_full))))), silent = TRUE)
          if (!inherits(rf_fit, "try-error")) pred_full <- tryCatch(as.numeric(predict(rf_fit, X_full, type = "prob")[, 2]), error = function(e) pred_full)
        }
      } else if (winner == "XGB") {
        if (has_xgb) {
          dtr <- xgboost::xgb.DMatrix(data = as.matrix(X_full), label = y_full)
          spw <- sum(y_full == 0) / max(1, sum(y_full == 1))
          params <- list(objective = "binary:logistic", eval_metric = "auc", max_depth = 4, eta = 0.05, subsample = 0.7, colsample_bytree = 0.7, min_child_weight = 5, gamma = 1, alpha = 1, lambda = 1, scale_pos_weight = spw)
          xgb <- try(xgboost::xgb.train(params = params, data = dtr, nrounds = 100, watchlist = list(train = dtr), verbose = 0), silent = TRUE)
          if (!inherits(xgb, "try-error")) pred_full <- tryCatch(as.numeric(predict(xgb, dtr)), error = function(e) pred_full)
        }
      } else if (winner == "MLP") {
        if (has_keras && run_nn) {
          K <- if (requireNamespace("keras3", quietly = TRUE)) getNamespace("keras3") else getNamespace("keras")
          model <- try({
            m <- K$keras_model_sequential()
            m$add(K$layer_dense(units = 128L, input_shape = as.integer(ncol(X_full_sc)), activation = "relu"))
            m$add(K$layer_dropout(rate = 0.3))
            m$add(K$layer_dense(units = 64L, activation = "relu"))
            m$add(K$layer_dropout(rate = 0.3))
            m$add(K$layer_dense(units = 1L, activation = "sigmoid"))
            m$compile(optimizer = K$optimizer_adam(learning_rate = 0.001), loss = "binary_crossentropy")
            m$fit(as.matrix(X_full_sc), y_full, epochs = 30L, batch_size = 64L, validation_split = 0.1, verbose = 0L)
            m
          }, silent = TRUE)
          if (!inherits(model, "try-error")) pred_full <- tryCatch(as.numeric(model$predict(as.matrix(X_full_sc))), error = function(e) pred_full)
        }
      } else if (winner == "ENSEMBLE") {
        # Build ensemble from available base full-data fits
        pred_components <- list()
        # Ridge
        cvm_r <- try(cv.glmnet(X_full_sc, y_full, family = "binomial", alpha = 0, type.measure = "auc", nfolds = 5), silent = TRUE)
        if (!inherits(cvm_r, "try-error")) pred_components[["ridge"]] <- as.numeric(predict(cvm_r, X_full_sc, s = "lambda.1se", type = "response"))
        # Lasso
        cvm_l <- try(cv.glmnet(X_full_sc, y_full, family = "binomial", alpha = 1, type.measure = "auc", nfolds = 5), silent = TRUE)
        if (!inherits(cvm_l, "try-error")) pred_components[["lasso"]] <- as.numeric(predict(cvm_l, X_full_sc, s = "lambda.1se", type = "response"))
        # Logistic L2 (same as ridge)
        if (!inherits(cvm_r, "try-error")) pred_components[["logl2"]] <- as.numeric(predict(cvm_r, X_full_sc, s = "lambda.1se", type = "response"))
        # RF
        if (has_rf) {
          rf_fit <- try(randomForest::randomForest(X_full, as.factor(y_full), ntree = 300, mtry = max(1, floor(sqrt(ncol(X_full))))), silent = TRUE)
          if (!inherits(rf_fit, "try-error")) pred_components[["rf"]] <- tryCatch(as.numeric(predict(rf_fit, X_full, type = "prob")[, 2]), error = function(e) NULL)
        }
        # XGB
        if (has_xgb) {
          dtr <- xgboost::xgb.DMatrix(data = as.matrix(X_full), label = y_full)
          spw <- sum(y_full == 0) / max(1, sum(y_full == 1))
          params <- list(objective = "binary:logistic", eval_metric = "auc", max_depth = 4, eta = 0.05, subsample = 0.7, colsample_bytree = 0.7, min_child_weight = 5, gamma = 1, alpha = 1, lambda = 1, scale_pos_weight = spw)
          xgb <- try(xgboost::xgb.train(params = params, data = dtr, nrounds = 100, watchlist = list(train = dtr), verbose = 0), silent = TRUE)
          if (!inherits(xgb, "try-error")) pred_components[["xgb"]] <- tryCatch(as.numeric(predict(xgb, dtr)), error = function(e) NULL)
        }
        # MLP optional
        if (has_keras && run_nn) {
          K <- if (requireNamespace("keras3", quietly = TRUE)) getNamespace("keras3") else getNamespace("keras")
          model <- try({
            m <- K$keras_model_sequential()
            m$add(K$layer_dense(units = 128L, input_shape = as.integer(ncol(X_full_sc)), activation = "relu"))
            m$add(K$layer_dropout(rate = 0.3))
            m$add(K$layer_dense(units = 64L, activation = "relu"))
            m$add(K$layer_dropout(rate = 0.3))
            m$add(K$layer_dense(units = 1L, activation = "sigmoid"))
            m$compile(optimizer = K$optimizer_adam(learning_rate = 0.001), loss = "binary_crossentropy")
            m$fit(as.matrix(X_full_sc), y_full, epochs = 30L, batch_size = 64L, validation_split = 0.1, verbose = 0L)
            m
          }, silent = TRUE)
          if (!inherits(model, "try-error")) pred_components[["mlp"]] <- tryCatch(as.numeric(model$predict(as.matrix(X_full_sc))), error = function(e) NULL)
        }
        avail <- pred_components[!sapply(pred_components, is.null)]
        if (length(avail) > 0) pred_full <- Reduce(`+`, avail) / length(avail)
      }

      # Compute threshold from CV predictions of the selected model
      get_cv_pred <- function(nm) {
        pred_col <- tolower(nm)
        if (nm == "LOGL2") pred_col <- "logl2"
        if (nm == "ENET") pred_col <- "enet"
        if (nm == "ENET_MAIN") return(NULL)
        return(all_preds[[pred_col]])
      }
      cv_pred <- if (winner == "ENSEMBLE") ensemble_cv else get_cv_pred(winner)
      cv_thr <- if (!is.null(cv_pred)) compute_youden_threshold(y_cv, cv_pred) else 0.5

      # Build final preds_dt from winner
      preds_dt_final <- data.table(
        SAMPLE_ID = sample_ids_keep,
        FINNGENID = cov_dt$FINNGENID,
        predicted_prob = pred_full,
        predicted_sex = ifelse(pred_full >= cv_thr, "female", "male"),
        genetic_sex = genetic_sex_str_all
      )
      preds_dt_final[, mismatch := !is.na(genetic_sex) & (genetic_sex != predicted_sex)]
      log_info("Selected model = {winner}; using threshold = {round(cv_thr, 3)}")
      # Also predict ALL samples under winner transforms/threshold
      log_info("Generating ALL-sample predictions for winner")
      # Prepare features for ALL samples using same selected features
      X_full_all <- npx_matrix[, sel_all, drop = FALSE]
      # Standardize using training means/sds
      X_full_all_sc <- sweep(sweep(X_full_all, 2, sc_all$means, "-"), 2, sc_all$sds, "/")
      pred_full_all <- rep(NA_real_, nrow(X_full_all_sc))
      if (winner == "RIDGE" || winner == "LOGL2") {
        pred_full_all <- as.numeric(predict(cvm, as.matrix(X_full_all_sc), s = "lambda.1se", type = "response"))
      } else if (winner == "LASSO") {
        pred_full_all <- as.numeric(predict(cvm, as.matrix(X_full_all_sc), s = "lambda.1se", type = "response"))
      } else if (winner == "ENET") {
        pred_full_all <- as.numeric(predict(best_mod2, as.matrix(X_full_all_sc), s = "lambda.1se", type = "response"))
      } else if (winner == "RF") {
        pred_full_all <- tryCatch(as.numeric(predict(rf_fit, X_full_all, type = "prob")[, 2]), error = function(e) pred_full_all)
      } else if (winner == "XGB") {
        dfull <- xgboost::xgb.DMatrix(data = as.matrix(X_full_all))
        pred_full_all <- tryCatch(as.numeric(predict(xgb, dfull)), error = function(e) pred_full_all)
      } else if (winner == "MLP") {
        pred_full_all <- tryCatch(as.numeric(model$predict(as.matrix(X_full_all_sc))), error = function(e) pred_full_all)
      } else if (winner == "ENSEMBLE") {
        comp_all <- list()
        comp_all[["ridge"]] <- tryCatch(as.numeric(predict(cvm_r, as.matrix(X_full_all_sc), s = "lambda.1se", type = "response")), error = function(e) NULL)
        comp_all[["lasso"]] <- tryCatch(as.numeric(predict(cvm_l, as.matrix(X_full_all_sc), s = "lambda.1se", type = "response")), error = function(e) NULL)
        comp_all[["logl2"]] <- tryCatch(as.numeric(predict(cvm_r, as.matrix(X_full_all_sc), s = "lambda.1se", type = "response")), error = function(e) NULL)
        comp_all[["rf"]] <- tryCatch(as.numeric(predict(rf_fit, X_full_all, type = "prob")[, 2]), error = function(e) NULL)
        comp_all[["xgb"]] <- tryCatch({
          dfull <- xgboost::xgb.DMatrix(data = as.matrix(X_full_all))
          as.numeric(predict(xgb, dfull))
        }, error = function(e) NULL)
        comp_all[["mlp"]] <- tryCatch(as.numeric(model$predict(as.matrix(X_full_all_sc))), error = function(e) NULL)
        comp_av <- comp_all[!sapply(comp_all, is.null)]
        if (length(comp_av) > 0) pred_full_all <- Reduce(`+`, comp_av) / length(comp_av)
      }
      # Compute Youden's J threshold from winner's CV predictions (BEFORE using it)
      thr_winner_youden <- if (winner == "ENSEMBLE") {
        pick_threshold(ensemble_cv, y_cv)$th
      } else {
        pred_col <- tolower(winner); if (winner == "LOGL2") pred_col <- "logl2"
        pick_threshold(all_preds[[pred_col]], y_cv)$th
      }

      preds_dt_final_all <- data.table(
        SAMPLE_ID = pheno_dt$SAMPLE_ID,
        FINNGENID = pheno_dt$FINNGENID,
        predicted_prob = pred_full_all,
        # Predicted sex based on 0.5 threshold
        predicted_sex = ifelse(pred_full_all >= 0.5, "female", "male"),
        genetic_sex = genetic_sex_str_full
      )
      # STRICT MISMATCH: predicted_sex != genetic_sex (actual label mismatch)
      # Must be calculated AFTER predicted_sex is created (data.table limitation)
      preds_dt_final_all[, strict_mismatch := (!is.na(genetic_sex) & !is.na(predicted_sex) &
                                               genetic_sex != predicted_sex)]
      if (any(duplicated(preds_dt_final_all$SAMPLE_ID))) preds_dt_final_all <- preds_dt_final_all[!duplicated(SAMPLE_ID)]
      # Threshold-based flags (for reference and outlier detection):
      # Mismatch: using Youden's J threshold (more stringent threshold-based flag)
      preds_dt_final_all[, mismatch := (!is.na(genetic_sex) & genetic_sex == "male" & predicted_prob >= thr_winner_youden) |
                                       (!is.na(genetic_sex) & genetic_sex == "female" & predicted_prob < thr_winner_youden)]
      # Sex outlier: deviates from 0.5 but not from Youden's J (mild threshold-based outlier)
      preds_dt_final_all[, sex_outlier := (!is.na(genetic_sex) & genetic_sex == "male" & predicted_prob >= 0.5 & predicted_prob < thr_winner_youden)]

      # Overwrite outputs with selected model results (ALL samples) - batch-aware paths
      fwrite(preds_dt_final_all, sex_predictions_path, sep = "\t")

      # Add BIOBANK_PLASMA and DISEASE_GROUP before writing mismatches/outliers
      preds_dt_final_all_with_meta <- add_metadata_to_predictions(preds_dt_final_all, metadata)

      # ============================================================================
      # CREATE TWO SEPARATE LISTS (Enhanced Models):
      # 1. Strict mismatches: predicted_sex != genetic_sex (actual label mismatch)
      # 2. Outliers: Threshold-based outliers that are NOT strict mismatches
      # ============================================================================

      # List 1: Strict mismatches (predicted_sex != genetic_sex)
      strict_mismatches_enhanced <- preds_dt_final_all_with_meta[strict_mismatch == TRUE]
      fwrite(strict_mismatches_enhanced, sex_mismatches_path, sep = "\t")
      log_info("Saved {nrow(strict_mismatches_enhanced)} strict mismatches (predicted_sex != genetic_sex) to: {sex_mismatches_path}")

      # List 2: Threshold-based outliers that are NOT strict mismatches
      threshold_outliers_enhanced <- preds_dt_final_all_with_meta[
        (mismatch == TRUE | isTRUE(sex_outlier)) &
        strict_mismatch == FALSE
      ]
      fwrite(threshold_outliers_enhanced, sex_outliers_path, sep = "\t")
      log_info("Saved {nrow(threshold_outliers_enhanced)} threshold-based outliers (not strict mismatches) to: {sex_outliers_path}")

      # Summary statistics
      log_info("=== SEX DETECTION SUMMARY (Enhanced Models) ===")
      log_info("Strict mismatches (predicted_sex != genetic_sex): {nrow(strict_mismatches_enhanced)}")
      log_info("Threshold-based outliers (not strict mismatches): {nrow(threshold_outliers_enhanced)}")
      log_info("Total flagged samples: {nrow(strict_mismatches_enhanced) + nrow(threshold_outliers_enhanced)}")
      preds_dt_out <- preds_dt_final_all

      # Winner annotated distribution panels (with and without FINNGENID labels)
      plot_with_ids <- tryCatch(isTRUE(config$parameters$outliers$plot_with_ids), error = function(e) FALSE)
      plot_title <- if (winner == "ENET_MAIN") {
        "Sex Prediction Distribution – Main Pipeline (Elastic Net Logistic Regression)"
      } else if (winner == "ENSEMBLE") {
        paste0(model_display_name("ensemble", ensemble_member_names), " – Sex Prediction Distribution")
      } else {
        paste0(model_display_name(winner), " – Sex Prediction Distribution")
      }
      create_prediction_distribution_panels(
        df_pred = preds_dt_out[, .(SAMPLE_ID, FINNGENID, genetic_sex, predicted_prob)],
        title = plot_title,
        file_out = sex_predict_distribution_path,
        youden_thresh = thr_winner_youden, annotate_ids = FALSE
      )
      if (isTRUE(plot_with_ids)) {
        create_prediction_distribution_panels(
          df_pred = preds_dt_out[, .(SAMPLE_ID, FINNGENID, genetic_sex, predicted_prob)],
          title = paste0(plot_title, "\n(with FINNGENID labels)"),
          file_out = sex_predict_distribution_with_ids_path,
          youden_thresh = thr_winner_youden, annotate_ids = TRUE
        )
      }
      # Update by-plate summary if possible
      if (file.exists(plate_map_file)) {
        samples_data_raw <- try(readRDS(plate_map_file), silent = TRUE)
        if (!inherits(samples_data_raw, "try-error")) {
          plate_dt <- unique(samples_data_raw[, .(SampleID, PlateID)])
          plate_merge <- merge(preds_dt_out[, .(SampleID = SAMPLE_ID, strict_mismatch)], plate_dt, by = "SampleID", all.x = TRUE)
          plate_merge[, strict_mismatch := as.logical(strict_mismatch)]
          plate_summary <- plate_merge[, .(
            sexOutRate = mean(strict_mismatch, na.rm = TRUE),
            Noutlier = sum(strict_mismatch, na.rm = TRUE),
            N = .N
          ), by = PlateID][order(-sexOutRate)]
          fwrite(plate_summary, sex_outlier_plate_path, sep = "\t")
        }
      }
    } else {
      # Winner is ENET_MAIN; ensure current files reflect ENET main results (already saved)
      log_info("Selected model = ENET_MAIN; keeping existing ENET-based outputs")
    }
    # Write winner per-fold metrics using the same schema as main (replace prior all_models file)
    winner_key <- if (winner == "ENET_MAIN") "enet" else tolower(winner)
    if (winner == "LOGL2") winner_key <- "logl2"
    winner_metrics <- if (winner == "ENET_MAIN") metrics_dt else {
      m <- metrics_by_model[[winner_key]]
      # Ensure all required columns exist
      need <- c("fold","auc","pr_auc","th","J","F1","sens","spec","acc","alpha","lambda")
      for (nm in need) if (!(nm %in% names(m))) m[, (nm) := NA_real_]
      m[, ..need]
    }
    fwrite(winner_metrics, file.path(enhanced_dir, "model_metrics.tsv"), sep = "\t")

    # Create top-20 feature histograms for the winner
    top_features <- character(0)
    if (winner == "ENET_MAIN") {
      if (exists("enet_coef_summary") && nrow(enet_coef_summary) > 0) {
        top_features <- enet_coef_summary[order(-abs(enet_beta_mean))][1:min(20, .N)]$protein
      }
    } else if (winner %in% c("RIDGE","LOGL2")) {
      co <- try(as.matrix(coef(cvm, s = "lambda.1se")), silent = TRUE)
      if (!inherits(co, "try-error")) {
        cf <- data.table(feature = rownames(co), beta = as.numeric(co[, 1]))
        cf <- cf[feature %in% colnames(Y) & feature != "(Intercept)"]
        top_features <- cf[order(-abs(beta))][1:min(20, .N)]$feature
      }
    } else if (winner == "LASSO") {
      co <- try(as.matrix(coef(cvm, s = "lambda.1se")), silent = TRUE)
      if (!inherits(co, "try-error")) {
        cf <- data.table(feature = rownames(co), beta = as.numeric(co[, 1]))
        cf <- cf[feature %in% colnames(Y) & feature != "(Intercept)"]
        top_features <- cf[order(-abs(beta))][1:min(20, .N)]$feature
      }
    } else if (winner == "ENET") {
      co <- try(as.matrix(coef(best_mod2, s = "lambda.1se")), silent = TRUE)
      if (!inherits(co, "try-error")) {
        cf <- data.table(feature = rownames(co), beta = as.numeric(co[, 1]))
        cf <- cf[feature %in% colnames(Y) & feature != "(Intercept)"]
        top_features <- cf[order(-abs(beta))][1:min(20, .N)]$feature
      }
    } else if (winner == "RF") {
      imp <- try(randomForest::importance(rf_fit), silent = TRUE)
      if (!inherits(imp, "try-error")) {
        vi <- data.table(feature = rownames(imp), score = as.numeric(imp[, 1]))
        vi <- vi[feature %in% colnames(Y)]
        top_features <- vi[order(-score)][1:min(20, .N)]$feature
      }
    } else if (winner == "XGB") {
      if (exists("sel_all")) {
        vi <- sapply(sel_all, function(f) {
          xi <- Y[, f]; v <- try(abs(cor(xi, y_bin[keep], use = "complete.obs")), silent = TRUE)
          if (inherits(v, "try-error") || !is.finite(v)) 0 else v
        })
        ord <- order(-vi)
        top_features <- sel_all[ord][1:min(20, length(sel_all))]
      }
    } else if (winner == "MLP") {
      if (exists("sel_all")) {
        vi <- sapply(sel_all, function(f) {
          xi <- Y[, f]; v <- try(abs(cor(xi, y_bin[keep], use = "complete.obs")), silent = TRUE)
          if (inherits(v, "try-error") || !is.finite(v)) 0 else v
        })
        ord <- order(-vi)
        top_features <- sel_all[ord][1:min(20, length(sel_all))]
      }
    } else if (winner == "ENSEMBLE") {
      if (exists("cvm_r")) {
        co <- try(as.matrix(coef(cvm_r, s = "lambda.1se")), silent = TRUE)
        if (!inherits(co, "try-error")) {
          cf <- data.table(feature = rownames(co), beta = as.numeric(co[, 1]))
          cf <- cf[feature %in% colnames(Y) & feature != "(Intercept)"]
          top_features <- cf[order(-abs(beta))][1:min(20, .N)]$feature
        }
      }
      if (length(top_features) == 0 && exists("sel_all")) {
        vi <- sapply(sel_all, function(f) {
          xi <- Y[, f]; v <- try(abs(cor(xi, y_bin[keep], use = "complete.obs")), silent = TRUE)
          if (inherits(v, "try-error") || !is.finite(v)) 0 else v
        })
        ord <- order(-vi)
        top_features <- sel_all[ord][1:min(20, length(sel_all))]
      }
    }
    if (length(top_features) > 0) {
      hist_top20 <- create_top_hist_facets(npx_matrix, sex_proteins$associations, sex_proteins$sample_sex,
                                           feature_list = top_features, top_n = 20)
      if (!is.null(hist_top20)) {
        sex_top20_histograms_path <- get_output_path(step_num, "sex_top20_histograms", batch_id, "outliers", "pdf", config = config)
        ensure_output_dir(sex_top20_histograms_path)
        try(ggsave(sex_top20_histograms_path, hist_top20, width = 14, height = 12), silent = TRUE)
      }
    }
  }

  # Print summary
  cat("\n=== SEX OUTLIER DETECTION SUMMARY ===\n")
  cat("Sex-associated proteins identified:", sum(sex_proteins$associations$bonf_sig, na.rm = TRUE), "\n")
  cat("Nested-CV mean AUC:", round(mean(metrics_dt$auc, na.rm = TRUE), 3), "; PR-AUC:", round(mean(metrics_dt$pr_auc, na.rm = TRUE), 3), "\n")
  cat("Optimized threshold (outer-CV, Youden J):", round(th_global, 3), "\n")
  cat("Total strict sex mismatches (predicted_sex != genetic_sex):", sum(preds_dt_out$strict_mismatch, na.rm = TRUE), "\n")
  cat("Total threshold-based outliers:", sum((preds_dt_out$mismatch == TRUE | isTRUE(preds_dt_out$sex_outlier)) & preds_dt_out$strict_mismatch == FALSE, na.rm = TRUE), "\n")
  cat("\n=== SMOKING-ASSOCIATED PROTEINS SUMMARY ===\n")
  cat("Smoking-associated proteins (Bonferroni):", sum(smoking_proteins$associations$bonf_sig, na.rm = TRUE), "\n")
  cat("Smoking-associated proteins (FDR < 0.05):", sum(smoking_proteins$associations$significant, na.rm = TRUE), "\n")

  cat("\nResults saved to: output/outliers/", if(batch_id != config$batch$default_batch_id ||
    tryCatch(isTRUE(config$parameters$normalization$multi_batch_mode), error = function(e) FALSE))
    paste0(batch_id, "/") else "", "\n", sep = "")

  log_info("Sex outlier detection completed")

  # Integrate into prior-step summaries (batch-aware paths)
  pos_file <- get_output_path(prev_step04_num, "pca_outlier_summary", batch_id, "outliers", "tsv", config = config)
  if (file.exists(pos_file)) {
    pos <- fread(pos_file)
    # Remove any prior prediction columns to avoid duplicate names across reruns
    drop_cols <- grep("^(predicted_prob(\\.|$)|predicted_sex(\\.|$)|sex_mismatch(\\.|$))", names(pos), value = TRUE)
    if (length(drop_cols) > 0) pos[, (drop_cols) := NULL]
    pos <- merge(pos,
                 preds_dt_out[, .(SampleID = SAMPLE_ID,
                                  predicted_prob,
                                  predicted_sex,
                                  sex_mismatch = strict_mismatch)],
                 by = "SampleID", all.x = TRUE)
    fwrite(pos, pos_file, sep = "\t")
  }
  step_file <- get_output_path(prev_step04_num, "pca_outliers_by_step", batch_id, "outliers", "tsv", config = config)
  if (file.exists(step_file)) {
    step_dt <- fread(step_file)
    sex_rows <- preds_dt_out[strict_mismatch == TRUE, .(step = "SexMismatch", SampleID = SAMPLE_ID)]
    step_dt <- rbind(step_dt, sex_rows, fill = TRUE)
    fwrite(step_dt, step_file, sep = "\t")
  }
  src_file <- get_output_path(prev_step04_num, "pca_outliers_by_source", batch_id, "outliers", "tsv", config = config)
  if (file.exists(src_file)) {
    src <- fread(src_file)
    # Drop any old SexMismatch columns to avoid .x/.y accumulation
    if ("SexMismatch" %in% names(src)) src[, SexMismatch := NULL]
    src <- merge(src,
                 preds_dt_out[, .(FINNGENID, SampleID = SAMPLE_ID, SexMismatch = as.integer(strict_mismatch))],
                 by = c("FINNGENID","SampleID"), all.x = TRUE)
    fwrite(src, src_file, sep = "\t")
  }

  return(list(
    sex_proteins = sex_proteins,
    smoking_proteins = smoking_proteins,
    metrics = metrics_dt,
    predictions = preds_dt_out
  ))
}

# Run if executed directly
if (!interactive()) {
  result <- main()
}
