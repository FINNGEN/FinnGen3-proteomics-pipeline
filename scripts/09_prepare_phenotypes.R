#!/usr/bin/env Rscript
# ==============================================================================
# 09_prepare_phenotypes.R - Phenotype Matrix Preparation
# ==============================================================================
#
# Purpose:
#   Prepares phenotype matrices for downstream GWAS analysis. Uses comprehensive
#   outlier list from Step 05d (includes all QC steps: Initial QC, PCA, Technical,
#   Z-score, Sex, pQTL) and converts outlier SampleIDs to match matrix format.
#   Creates both SampleID-indexed and FINNGENID-indexed matrices with outliers removed.
#
# Author: Reza Jabal, PhD (rjabal@broadinstitute.org)
# Date: December 2025 (Updated: January 2026 - Comprehensive outlier integration and SampleID conversion)
# ==============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(tidyverse)
  library(yaml)
  library(logger)
})

# Suppress "no visible binding" warnings for data.table operations
utils::globalVariables(
  c("is_outlier", "SampleID", "SAMPLE_ID", "FINNGENID", ".")
)

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
if (config_file == "" || !file.exists(config_file)) {
  stop("PIPELINE_CONFIG environment variable not set or config file not found. Please provide path to config file.")
}
config <- read_yaml(config_file)

# Set up logging with batch-aware path
log_path <- get_log_path(step_num, batch_id, config = config)
log_appender(appender_file(log_path))
log_info("Starting phenotype preparation for batch: {batch_id}")

# Function to load comprehensive outlier list from Step 05d
# CRITICAL: Use comprehensive outlier list from Step 05d which includes ALL outliers:
# - Initial QC failures (step 00)
# - PCA outliers (step 01)
# - Technical outliers (step 02)
# - Z-score outliers (step 03)
# - Sex mismatches and outliers (step 04)
# - pQTL outliers (step 05b)
# NOTE: Comprehensive list uses "P..." format SampleIDs, but matrices use "EA5_OLI_..." format
# We need to convert using sample mapping to ensure proper matching
load_comprehensive_outlier_list <- function(batch_id, config) {
  log_info("Loading comprehensive outlier list from Step 05d")

  # Use comprehensive outlier list from Step 05d (includes ALL QC steps)
  comprehensive_outliers_path <- get_output_path("05d", "comprehensive_outliers_list", batch_id, "phenotypes", "tsv", config = config)

  if (!file.exists(comprehensive_outliers_path)) {
    log_warn("Comprehensive outlier list not found: {comprehensive_outliers_path}")
    log_warn("Falling back to combining individual outlier lists (may miss pQTL and initial QC outliers)")
    return(combine_outlier_lists_fallback(batch_id, config))
  }

  # Read comprehensive outlier list
  comprehensive_outliers <- fread(comprehensive_outliers_path)

  # Extract SampleID column (should exist in comprehensive list)
  if (!"SampleID" %in% names(comprehensive_outliers)) {
    log_error("Comprehensive outlier list missing SampleID column")
    stop("Invalid comprehensive outlier list format")
  }

  # CRITICAL: Comprehensive list uses "P..." format SampleIDs, but matrices use "EA5_OLI_..." format
  # Load sample mapping to convert outlier SampleIDs to match matrix format
  sample_mapping_path <- get_output_path("00", "sample_mapping", batch_id, "qc", "tsv", config = config)
  if (!file.exists(sample_mapping_path)) {
    log_warn("Sample mapping not found: {sample_mapping_path}")
    log_warn("Cannot convert outlier SampleIDs - using as-is (may not match matrix)")
    all_outliers <- comprehensive_outliers$SampleID
  } else {
    sample_mapping <- fread(sample_mapping_path)
    log_info("Loaded sample mapping: {nrow(sample_mapping)} samples")

    # The comprehensive list SampleIDs are in "P..." format
    # Check if they match the mapping directly
    outliers_in_mapping <- comprehensive_outliers$SampleID %in% sample_mapping$SampleID
    log_info("Outliers matching sample mapping: {sum(outliers_in_mapping)}/{nrow(comprehensive_outliers)}")

    # Use SampleIDs directly from comprehensive list (they should match matrix row names)
    # If comprehensive list uses "P..." format but matrix uses "EA5_OLI_...", we need to check
    # Actually, the comprehensive list SampleIDs should match what's in the base matrix
    # Let's use them as-is first, and verify matching later
    all_outliers <- comprehensive_outliers$SampleID
    log_info("Using outlier SampleIDs from comprehensive list (format: {head(all_outliers, 1)})")
  }

  all_outliers <- unique(all_outliers[!is.na(all_outliers)])

  log_info("Loaded {length(all_outliers)} outliers from comprehensive list (Step 05d)")
  log_info("  This includes outliers from: Initial QC, PCA, Technical, Z-score, Sex, and pQTL")

  # Create outlier sources breakdown for reporting
  outlier_sources <- list()
  if ("QC_initial_qc" %in% names(comprehensive_outliers)) {
    outlier_sources$initial_qc <- comprehensive_outliers[QC_initial_qc == 1]$SampleID
  }
  if ("QC_pca" %in% names(comprehensive_outliers)) {
    outlier_sources$pca <- comprehensive_outliers[QC_pca == 1]$SampleID
  }
  if ("QC_technical" %in% names(comprehensive_outliers)) {
    outlier_sources$technical <- comprehensive_outliers[QC_technical == 1]$SampleID
  }
  if ("QC_zscore" %in% names(comprehensive_outliers)) {
    outlier_sources$zscore <- comprehensive_outliers[QC_zscore == 1]$SampleID
  }
  if ("QC_sex_mismatch" %in% names(comprehensive_outliers)) {
    outlier_sources$sex_mismatch <- comprehensive_outliers[QC_sex_mismatch == 1]$SampleID
  }
  if ("QC_sex_outlier" %in% names(comprehensive_outliers)) {
    outlier_sources$sex_outlier <- comprehensive_outliers[QC_sex_outlier == 1]$SampleID
  }
  if ("QC_pqtl" %in% names(comprehensive_outliers)) {
    outlier_sources$pqtl <- comprehensive_outliers[QC_pqtl == 1]$SampleID
  }

  return(list(
    all_outliers = all_outliers,
    outlier_sources = outlier_sources,
    comprehensive_outliers_dt = comprehensive_outliers  # Return full DT for ID conversion
  ))
}

# Function to convert outlier SampleIDs to match matrix format
# Comprehensive list may use "P..." format, but matrices use "EA5_OLI_..." format
convert_outlier_sampleids <- function(outlier_sampleids, sample_mapping, npx_matrix, comprehensive_outliers_dt = NULL) {
  log_info("Converting {length(outlier_sampleids)} outlier SampleIDs to match matrix format")

  # Get matrix row names (target format)
  matrix_sampleids <- rownames(npx_matrix)

  # Check if outliers already match matrix format
  direct_matches <- intersect(outlier_sampleids, matrix_sampleids)
  log_info("  Direct matches (no conversion needed): {length(direct_matches)}")

  if (length(direct_matches) == length(outlier_sampleids)) {
    log_info("  All outliers already match matrix format - no conversion needed")
    return(outlier_sampleids)
  }

  # For non-matching outliers, try to convert via sample mapping using FINNGENID
  if (!is.null(comprehensive_outliers_dt) && "FINNGENID" %in% names(comprehensive_outliers_dt)) {
    # Match outliers to comprehensive list to get FINNGENIDs
    outlier_dt <- data.table(SampleID_outlier = outlier_sampleids)
    outlier_dt <- merge(outlier_dt, comprehensive_outliers_dt[, .(SampleID, FINNGENID)],
                       by.x = "SampleID_outlier", by.y = "SampleID", all.x = TRUE)

    # Match FINNGENIDs to sample mapping to get matrix SampleIDs
    if ("FINNGENID" %in% names(outlier_dt) && "FINNGENID" %in% names(sample_mapping)) {
      outlier_dt <- merge(outlier_dt, sample_mapping[, .(SampleID, FINNGENID)],
                         by = "FINNGENID", all.x = TRUE, suffixes = c("", "_matrix"))

      # Use matrix SampleID where available, otherwise keep original
      converted <- ifelse(!is.na(outlier_dt$SampleID), outlier_dt$SampleID, outlier_dt$SampleID_outlier)

      n_converted <- sum(!is.na(outlier_dt$SampleID) & outlier_dt$SampleID != outlier_dt$SampleID_outlier)
      log_info("  Converted via FINNGENID: {n_converted} outliers")

      # Verify converted IDs are in matrix
      final_matches <- intersect(converted, matrix_sampleids)
      log_info("  Final matches in matrix: {length(final_matches)}/{length(converted)}")

      return(converted)
    }
  }

  # Fallback: return original (will be checked in prepare_phenotype_matrix)
  log_warn("  Could not convert outlier SampleIDs - using original format")
  return(outlier_sampleids)
}

# Fallback function to combine outlier lists from individual steps (if comprehensive list not available)
combine_outlier_lists_fallback <- function(batch_id, config) {
  log_info("Combining outlier lists from individual detection methods (fallback mode)")

  # Use batch-aware paths for outlier files
  outlier_files <- list(
    pca = get_output_path("01", "pca_outliers_list", batch_id, "outliers", "tsv", config = config),
    sex = get_output_path("04", "sex_mismatches", batch_id, "outliers", "tsv", config = config),
    technical = get_output_path("02", "technical_outlier_summary", batch_id, "outliers", "tsv", config = config),
    zscore = get_output_path("03", "zscore_outlier_summary", batch_id, "outliers", "tsv", config = config)
  )

  all_outliers <- character()
  outlier_sources <- list()

  # Read PCA outliers
  if(file.exists(outlier_files$pca)) {
    pca_outliers <- fread(outlier_files$pca)$SampleID
    all_outliers <- c(all_outliers, pca_outliers)
    outlier_sources$pca <- pca_outliers
    log_info("PCA outliers: {length(pca_outliers)}")
  }

  # Read sex mismatches
  if(file.exists(outlier_files$sex)) {
    sex_outliers <- fread(outlier_files$sex)$SAMPLE_ID
    all_outliers <- c(all_outliers, sex_outliers)
    outlier_sources$sex <- sex_outliers
    log_info("Sex mismatches: {length(sex_outliers)}")
  }

  # Read technical outliers
  if(file.exists(outlier_files$technical)) {
    tech_outliers <- fread(outlier_files$technical)$SampleID
    all_outliers <- c(all_outliers, tech_outliers)
    outlier_sources$technical <- tech_outliers
    log_info("Technical outliers: {length(tech_outliers)}")
  }

  # Read Z-score outliers
  if(file.exists(outlier_files$zscore)) {
    zscore_data <- fread(outlier_files$zscore)
    zscore_outliers <- zscore_data[is_outlier == TRUE, SampleID]
    all_outliers <- c(all_outliers, zscore_outliers)
    outlier_sources$zscore <- zscore_outliers
    log_info("Z-score outliers: {length(zscore_outliers)}")
  }

  # Get unique outliers
  all_outliers <- unique(all_outliers)

  log_info("Total unique outliers (fallback): {length(all_outliers)}")
  log_warn("  NOTE: Fallback mode may miss pQTL outliers and initial QC failures")

  return(list(
    all_outliers = all_outliers,
    outlier_sources = outlier_sources
  ))
}

# Function to remove non-FinnGen samples
remove_excluded_samples <- function(npx_matrix, excluded_samples) {
  log_info("Removing excluded samples (non-FinnGen samples)")

  # Get samples to exclude
  samples_to_remove <- intersect(rownames(npx_matrix), excluded_samples$SAMPLE_ID)

  if(length(samples_to_remove) > 0) {
    log_info("Removing {length(samples_to_remove)} excluded samples")
    npx_matrix <- npx_matrix[!rownames(npx_matrix) %in% samples_to_remove, ]
  }

  return(npx_matrix)
}

# Function to prepare phenotype matrix
prepare_phenotype_matrix <- function(npx_matrix, outliers_to_remove, sample_info) {
  log_info("Preparing phenotype matrix")

  # Remove outliers
  clean_samples <- setdiff(rownames(npx_matrix), outliers_to_remove)
  phenotype_matrix <- npx_matrix[clean_samples, ]

  log_info("Matrix after outlier removal: {nrow(phenotype_matrix)} x {ncol(phenotype_matrix)}")

  # Add FINNGENID as rownames if available
  if(!is.null(sample_info)) {
    # sample_mapping uses SampleID (not SAMPLE_ID)
    sample_mapping <- sample_info[SampleID %in% rownames(phenotype_matrix), .(SampleID, FINNGENID)]

    # Check for samples with FINNGENID
    samples_with_finngen <- sample_mapping[!is.na(FINNGENID)]

    if(nrow(samples_with_finngen) > 0) {
      # Create FINNGENID-indexed matrix
      finngen_matrix <- phenotype_matrix[samples_with_finngen$SampleID, ]
      rownames(finngen_matrix) <- samples_with_finngen$FINNGENID

      log_info("Created FINNGENID-indexed matrix: {nrow(finngen_matrix)} samples")

      return(list(
        sample_id_matrix = phenotype_matrix,
        finngenid_matrix = finngen_matrix,
        sample_mapping = samples_with_finngen
      ))
    }
  }

  return(list(
    sample_id_matrix = phenotype_matrix,
    finngenid_matrix = NULL,
    sample_mapping = NULL
  ))
}

# Function to create phenotype info file
create_phenotype_info <- function(phenotype_matrix, protein_info = NULL) {
  log_info("Creating phenotype information file")

  # Basic info
  pheno_info <- data.table(
    protein = colnames(phenotype_matrix),
    n_samples = colSums(!is.na(phenotype_matrix)),
    missing_rate = colSums(is.na(phenotype_matrix)) / nrow(phenotype_matrix),
    mean = colMeans(phenotype_matrix, na.rm = TRUE),
    sd = apply(phenotype_matrix, 2, sd, na.rm = TRUE),
    median = apply(phenotype_matrix, 2, median, na.rm = TRUE),
    min = apply(phenotype_matrix, 2, min, na.rm = TRUE),
    max = apply(phenotype_matrix, 2, max, na.rm = TRUE)
  )

  # Add protein annotations if available
  if(!is.null(protein_info)) {
    pheno_info <- merge(pheno_info, protein_info, by.x = "protein", by.y = "Assay", all.x = TRUE)
  }

  return(pheno_info)
}

# Function to format for different analysis tools
format_phenotypes <- function(phenotype_matrix, format = "plink") {
  log_info("Formatting phenotypes for {format}")

  if(format == "plink") {
    # PLINK format: FID IID phenotype1 phenotype2 ...
    plink_pheno <- data.table(
      FID = rownames(phenotype_matrix),
      IID = rownames(phenotype_matrix)
    )

    # Add phenotype columns
    plink_pheno <- cbind(plink_pheno, as.data.table(phenotype_matrix))

    # Replace NA with -9 for PLINK
    for(col in colnames(phenotype_matrix)) {
      plink_pheno[is.na(get(col)), (col) := -9]
    }

    return(plink_pheno)

  } else if(format == "matrix") {
    # Standard matrix format
    return(phenotype_matrix)

  } else if(format == "long") {
    # Long format for mixed models
    long_pheno <- melt(as.data.table(phenotype_matrix, keep.rownames = "sample_id"),
                      id.vars = "sample_id",
                      variable.name = "protein",
                      value.name = "expression")

    return(long_pheno)
  }
}

# Function to create QCed sample lists
create_sample_lists <- function(phenotype_result, outlier_result) {
  log_info("Creating QCed sample lists")

  # Samples passing all QC
  qc_pass_samples <- rownames(phenotype_result$sample_id_matrix)

  # Samples with FINNGENIDs
  finngen_samples <- if(!is.null(phenotype_result$finngenid_matrix)) {
    rownames(phenotype_result$finngenid_matrix)
  } else {
    character()
  }

  # Create summary
  sample_lists <- list(
    qc_pass_all = qc_pass_samples,
    qc_pass_finngen = finngen_samples,
    outliers_removed = outlier_result$all_outliers,
    n_original = NA,  # Would be set from input
    n_after_qc = length(qc_pass_samples),
    n_with_finngen = length(finngen_samples)
  )

  return(sample_lists)
}

# Function to merge batch matrices on FINNGENID (common proteins only)
merge_batch_matrices <- function(batch1_matrix, batch2_matrix, batch1_mapping, batch2_mapping) {
  log_info("Merging batch matrices on FINNGENID")

  # Check protein consistency
  common_proteins <- intersect(colnames(batch1_matrix), colnames(batch2_matrix))
  log_info("Common proteins for aggregation: {length(common_proteins)}")

  if (length(common_proteins) < 100) {
    log_error("Too few common proteins ({length(common_proteins)}) for aggregation")
    return(NULL)
  }

  # Use only common proteins
  batch1_subset <- batch1_matrix[, common_proteins, drop = FALSE]
  batch2_subset <- batch2_matrix[, common_proteins, drop = FALSE]

  # Get FINNGENIDs for each batch
  batch1_finngenids <- batch1_mapping[SampleID %in% rownames(batch1_subset) & !is.na(FINNGENID)]
  batch2_finngenids <- batch2_mapping[SampleID %in% rownames(batch2_subset) & !is.na(FINNGENID)]

  # Set FINNGENID as rownames
  batch1_finngen <- batch1_subset[batch1_finngenids$SampleID, , drop = FALSE]
  rownames(batch1_finngen) <- batch1_finngenids$FINNGENID

  batch2_finngen <- batch2_subset[batch2_finngenids$SampleID, , drop = FALSE]
  rownames(batch2_finngen) <- batch2_finngenids$FINNGENID

  # Find common FINNGENIDs (samples in both batches)
  common_finngenids <- intersect(rownames(batch1_finngen), rownames(batch2_finngen))
  log_info("Common FINNGENIDs (samples in both batches): {length(common_finngenids)}")

  if (length(common_finngenids) > 0) {
    log_warn("Found {length(common_finngenids)} samples present in both batches")
    log_warn("Using batch 2 data for common samples (batch 2 is reference)")
    # Remove common samples from batch 1 to avoid duplicates
    batch1_finngen <- batch1_finngen[!rownames(batch1_finngen) %in% common_finngenids, , drop = FALSE]
  }

  # Combine matrices
  aggregate_matrix <- rbind(batch1_finngen, batch2_finngen)
  log_info("Aggregate matrix: {nrow(aggregate_matrix)} samples x {ncol(aggregate_matrix)} proteins")

  return(aggregate_matrix)
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

  # Load data from previous steps with batch-aware paths
  log_info("Loading data from previous steps")
  # Step 08 saves as "npx_matrix_covariate_adjusted" (covariate-adjusted matrix)
  # NOTE: This matrix is adjusted for age, sex, BMI, and smoking ONLY
  # Proteomic PCs are NOT adjusted for (preserved to maintain biological signal)
  npx_adjusted_path <- get_output_path("08", "npx_matrix_covariate_adjusted", batch_id, "normalized", config = config)
  sample_mapping_path <- get_output_path("00", "sample_mapping", batch_id, "qc", config = config)
  excluded_samples_path <- get_output_path("00", "excluded_samples", batch_id, "qc", config = config)

  if (!file.exists(npx_adjusted_path)) {
    stop("Adjusted NPX matrix not found: {npx_adjusted_path}. Please run step 08 (covariate adjustment) first.")
  }
  log_info("Using covariate-adjusted NPX matrix: {npx_adjusted_path}")
  log_info("  Note: Matrix adjusted for age, sex, BMI, smoking only (proteomic PCs preserved)")
  if (!file.exists(sample_mapping_path)) {
    stop("Sample mapping not found: {sample_mapping_path}")
  }

  npx_adjusted <- readRDS(npx_adjusted_path)
  sample_mapping <- readRDS(sample_mapping_path)

  # excluded_samples is optional
  excluded_samples <- if (file.exists(excluded_samples_path)) {
    readRDS(excluded_samples_path)
  } else {
    log_warn("Excluded samples file not found: {excluded_samples_path}. Proceeding without exclusions.")
    data.table(SAMPLE_ID = character(0))
  }

  # Store original dimensions
  n_original <- nrow(npx_adjusted)

  # Load comprehensive outlier list from Step 05d (includes ALL outliers from all QC steps)
  outlier_result <- load_comprehensive_outlier_list(batch_id, config)

  # CRITICAL FIX: Convert outlier SampleIDs to match matrix format
  # Comprehensive list may use "P..." format, but matrix uses "EA5_OLI_..." format
  # Use sample mapping to convert via FINNGENID if needed
  log_info("Converting outlier SampleIDs to match matrix format...")
  outliers_original <- outlier_result$all_outliers
  outliers_converted <- convert_outlier_sampleids(
    outlier_result$all_outliers,
    sample_mapping,
    npx_adjusted,
    comprehensive_outliers_dt = outlier_result$comprehensive_outliers_dt
  )
  outlier_result$all_outliers <- outliers_converted

  n_matched <- sum(outliers_converted %in% rownames(npx_adjusted))
  log_info("Outlier ID conversion: {length(outliers_original)} original -> {length(outliers_converted)} converted")
  log_info("  Outliers matching matrix row names: {n_matched}/{length(outliers_converted)}")

  if (n_matched < length(outliers_converted)) {
    log_warn("  {length(outliers_converted) - n_matched} outliers not found in matrix (may have been removed earlier)")
  }

  # Remove non-FinnGen samples first
  npx_clean <- remove_excluded_samples(npx_adjusted, excluded_samples)
  log_info("After removing excluded samples: {nrow(npx_clean)} samples")

  # Prepare phenotype matrices
  phenotype_result <- prepare_phenotype_matrix(
    npx_clean,
    outlier_result$all_outliers,
    sample_mapping
  )

  # Create phenotype info
  pheno_info <- create_phenotype_info(phenotype_result$sample_id_matrix)

  # Format for different tools
  plink_format <- format_phenotypes(phenotype_result$sample_id_matrix, format = "plink")
  long_format <- format_phenotypes(phenotype_result$sample_id_matrix, format = "long")

  # Create sample lists
  sample_lists <- create_sample_lists(phenotype_result, outlier_result)
  sample_lists$n_original <- n_original

  # Save outputs
  log_info("Saving phenotype matrices and information")

  # Determine output prefix based on aggregation mode
  if (aggregate_output) {
    output_prefix <- "batch2_"
    log_info("Aggregation mode: Saving batch 2 outputs with 'batch2_' prefix")
  } else {
    output_prefix <- "11_"
    log_info("Single-batch mode: Saving outputs with '11_' prefix")
  }

  # Save matrices with batch-aware paths
  phenotype_matrix_path <- get_output_path(step_num, "phenotype_matrix", batch_id, "phenotypes", config = config)
  ensure_output_dir(phenotype_matrix_path)
  saveRDS(phenotype_result$sample_id_matrix, phenotype_matrix_path)

  if(!is.null(phenotype_result$finngenid_matrix)) {
    finngenid_matrix_path <- get_output_path(step_num, "phenotype_matrix_finngenid", batch_id, "phenotypes", config = config)
    ensure_output_dir(finngenid_matrix_path)
    saveRDS(phenotype_result$finngenid_matrix, finngenid_matrix_path)
  }

  # Save formatted versions
  plink_path <- get_output_path(step_num, "phenotypes_plink", batch_id, "phenotypes", "txt", config = config)
  long_path <- get_output_path(step_num, "phenotypes_long", batch_id, "phenotypes", "txt", config = config)
  ensure_output_dir(plink_path)
  ensure_output_dir(long_path)

  fwrite(plink_format, plink_path, sep = "\t", na = "-9")
  fwrite(long_format, long_path, sep = "\t")

  # Save info files
  pheno_info_path <- get_output_path(step_num, "phenotype_info", batch_id, "phenotypes", "tsv", config = config)
  ensure_output_dir(pheno_info_path)
  fwrite(pheno_info, pheno_info_path, sep = "\t")

  # Add FINNGENID to sample lists
  qc_pass_dt <- add_finngenid_column(
    data.table(SampleID = sample_lists$qc_pass_all),
    batch_id = batch_id, config = config
  )
  qc_pass_path <- get_output_path(step_num, "samples_qc_pass", batch_id, "phenotypes", "tsv", config = config)
  ensure_output_dir(qc_pass_path)
  fwrite(qc_pass_dt, qc_pass_path, sep = "\t")

  finngenids_path <- get_output_path(step_num, "finngenids_qc_pass", batch_id, "phenotypes", "txt", config = config)
  ensure_output_dir(finngenids_path)
  fwrite(data.table(FINNGENID = sample_lists$qc_pass_finngen), finngenids_path, col.names = FALSE)

  # Add FINNGENID to outlier list
  outliers_dt <- add_finngenid_column(
    data.table(SampleID = outlier_result$all_outliers),
    batch_id = batch_id, config = config
  )
  outliers_path <- get_output_path(step_num, "all_outliers_removed", batch_id, "phenotypes", "tsv", config = config)
  ensure_output_dir(outliers_path)
  fwrite(outliers_dt, outliers_path, sep = "\t")

  # Save summary
  summary_stats <- data.table(
    stage = c("Original", "After excluding non-FinnGen samples", "After removing outliers", "With FINNGENID"),
    n_samples = c(sample_lists$n_original, nrow(npx_clean),
                 sample_lists$n_after_qc, sample_lists$n_with_finngen)
  )
  summary_path <- get_output_path(step_num, "sample_summary", batch_id, "phenotypes", "tsv", config = config)
  ensure_output_dir(summary_path)
  fwrite(summary_stats, summary_path, sep = "\t")

  # Handle aggregation if enabled
  if (aggregate_output && multi_batch_mode && !is.null(phenotype_result$finngenid_matrix)) {
    log_info("Aggregation enabled: Attempting to merge batch 1 and batch 2 data")

    # Check if batch 1 processed data exists
    # Note: Batch 1 would need to be processed through steps 00-10 separately
    other_batch_id <- if (batch_id == "batch_02") "batch_01" else "batch_02"
    batch1_file <- get_output_path("09", "phenotype_matrix_finngenid", other_batch_id, "phenotypes", config = config)
    batch1_mapping_file <- get_output_path("00", "sample_mapping", other_batch_id, "qc", config = config)

    if (file.exists(batch1_file) && file.exists(batch1_mapping_file)) {
      log_info("Found batch 1 processed data: Creating aggregate outputs")

      # Load batch 1 data
      batch1_matrix <- readRDS(batch1_file)
      batch1_mapping <- readRDS(batch1_mapping_file)

      # Merge matrices
      aggregate_matrix <- merge_batch_matrices(
        batch1_matrix,
        phenotype_result$finngenid_matrix,
        batch1_mapping,
        sample_mapping
      )

      if (!is.null(aggregate_matrix)) {
        # Create aggregate phenotype info
        aggregate_pheno_info <- create_phenotype_info(aggregate_matrix)

        # Format for different tools
        aggregate_plink <- format_phenotypes(aggregate_matrix, format = "plink")
        aggregate_long <- format_phenotypes(aggregate_matrix, format = "long")

        # Save aggregate outputs with batch-aware paths
        log_info("Saving aggregate outputs with 'aggregate_' prefix")
        base_dir <- config$output$base_dir %||% Sys.getenv("PIPELINE_OUTPUT_DIR", "output")
        aggregate_dir <- file.path(base_dir, "phenotypes")
        dir.create(aggregate_dir, recursive = TRUE, showWarnings = FALSE)

        saveRDS(aggregate_matrix, file.path(aggregate_dir, "aggregate_phenotype_matrix.rds"))
        fwrite(aggregate_plink, file.path(aggregate_dir, "aggregate_phenotypes_plink.txt"), sep = "\t", na = "-9")
        fwrite(aggregate_long, file.path(aggregate_dir, "aggregate_phenotypes_long.txt"), sep = "\t")
        fwrite(aggregate_pheno_info, file.path(aggregate_dir, "aggregate_phenotype_info.tsv"), sep = "\t")
        fwrite(data.table(FINNGENID = rownames(aggregate_matrix)),
               file.path(aggregate_dir, "aggregate_finngenids_qc_pass.txt"), col.names = FALSE)

        log_info("Aggregate outputs saved: {nrow(aggregate_matrix)} samples x {ncol(aggregate_matrix)} proteins")
      } else {
        log_warn("Failed to create aggregate matrix")
      }
    } else {
      log_warn("Batch 1 processed data not found. Aggregation skipped.")
      log_warn("Expected files:")
      log_warn("  - {batch1_file}")
      log_warn("  - {batch1_mapping_file}")
      log_warn("Batch 1 must be processed through steps 00-10 separately before aggregation")
    }
  }

  # Print summary
  cat("\n=== PHENOTYPE PREPARATION SUMMARY ===\n")
  cat("Original samples:", sample_lists$n_original, "\n")
  cat("After excluding non-FinnGen samples:", nrow(npx_clean), "\n")
  cat("Outliers removed:", length(outlier_result$all_outliers), "\n")
  # Report breakdown by source (if available)
  if (!is.null(outlier_result$outlier_sources)) {
    cat("  - Initial QC failures:", length(outlier_result$outlier_sources$initial_qc %||% character(0)), "\n")
    cat("  - PCA outliers:", length(outlier_result$outlier_sources$pca %||% character(0)), "\n")
    cat("  - Sex mismatches:", length(outlier_result$outlier_sources$sex_mismatch %||% character(0)), "\n")
    cat("  - Sex outliers:", length(outlier_result$outlier_sources$sex_outlier %||% character(0)), "\n")
    cat("  - Technical outliers:", length(outlier_result$outlier_sources$technical %||% character(0)), "\n")
    cat("  - Z-score outliers:", length(outlier_result$outlier_sources$zscore %||% character(0)), "\n")
    cat("  - pQTL outliers:", length(outlier_result$outlier_sources$pqtl %||% character(0)), "\n")
  } else {
    cat("  (Outlier source breakdown not available)\n")
  }
  cat("\nFinal QC-passed samples:", sample_lists$n_after_qc, "\n")
  cat("Samples with FINNGENID:", sample_lists$n_with_finngen, "\n")
  cat("\nPhenotype matrix:", nrow(phenotype_result$sample_id_matrix), "x",
      ncol(phenotype_result$sample_id_matrix), "\n")
  cat("\nOutputs saved to: ../output/phenotypes/\n")

  log_info("Phenotype preparation completed")

  return(list(
    phenotype_matrix = phenotype_result$sample_id_matrix,
    finngenid_matrix = phenotype_result$finngenid_matrix,
    sample_lists = sample_lists,
    pheno_info = pheno_info
  ))
}

# Run if executed directly
if (!interactive()) {
  result <- main()
}






