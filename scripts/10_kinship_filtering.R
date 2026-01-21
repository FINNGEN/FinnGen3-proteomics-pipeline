#!/usr/bin/env Rscript
# ==============================================================================
# 10_kinship_filtering.R - Kinship Filtering and Quality-Based Sample Selection
# ==============================================================================
#
# Purpose:
#   Removes related individuals using kinship matrix (3rd degree threshold:
#   0.0884 KING kinship coefficient). For FINNGENIDs with multiple SampleIDs
#   (duplicate samples), selects the best quality sample using hierarchical
#   criteria: missing rate (lower) → n_detected_proteins (higher) → SD NPX
#   (lower) → SampleID (lexicographic). Ensures consistent matrix dimensions
#   (one sample per person).
#
# Author: Reza Jabal, PhD (rjabal@broadinstitute.org)
# Date: December 2025 (Updated: January 2026 - Quality-based sample selection implementation)
# Date: December 2025 (Updated: January 2026 - Quality-based sample selection implementation)
# ==============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(tidyverse)
  library(yaml)
  library(logger)
})

# Suppress linting warnings for data.table columns
utils::globalVariables(c(
  "KINSHIP", "relationship", "IID1", "IID2", "IID", "n_relatives",
  "missing_rate", "individuals", "FINNGENID1", "FINNGENID2", "MR1", "MR2",
  "ID1_kept", "ID2_kept", "correct", "."
))

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
log_info("Starting kinship filtering for batch: {batch_id}")

# Function to load kinship matrix
load_kinship <- function(kinship_file) {
  log_info("Loading kinship matrix from: {kinship_file}")

  # Read kinship file (KING format expected)
  kinship <- fread(kinship_file)

  # Check for required columns (KING format can have ID1/ID2 or IID1/IID2)
  # Kinship column can be "Kinship" or "KINSHIP"
  has_iid_format <- all(c("IID1", "IID2") %in% colnames(kinship))
  has_id_format <- all(c("ID1", "ID2") %in% colnames(kinship))

  if(!has_iid_format && !has_id_format) {
    log_error("Kinship file missing required ID columns (need IID1/IID2 or ID1/ID2)")
    stop("Invalid kinship file format")
  }

  # Standardize column names to IID1/IID2
  if(has_id_format && !has_iid_format) {
    setnames(kinship, c("ID1", "ID2"), c("IID1", "IID2"))
    log_info("Renamed ID1/ID2 to IID1/IID2")
  }

  # Check for kinship column
  if("Kinship" %in% colnames(kinship) && !"KINSHIP" %in% colnames(kinship)) {
    setnames(kinship, "Kinship", "KINSHIP")
    log_info("Renamed Kinship to KINSHIP")
  }

  if(!"KINSHIP" %in% colnames(kinship)) {
    log_error("Kinship file missing KINSHIP column")
    stop("Invalid kinship file format")
  }

  log_info("Kinship pairs loaded: {nrow(kinship)}")

  return(kinship)
}

# Function to identify related pairs
identify_related_pairs <- function(kinship, threshold = 0.0884) {
  log_info("Identifying related pairs with threshold: {threshold}")

  # Filter to related pairs (3rd degree or closer)
  # Kinship coefficients:
  # - Identical twins/duplicates: 0.5
  # - 1st degree (parent-child, full siblings): 0.25
  # - 2nd degree (grandparent, half-siblings): 0.125
  # - 3rd degree (cousins): 0.0625
  # Default threshold 0.0884 is between 2nd and 3rd degree

  related_pairs <- kinship[KINSHIP >= threshold]

  # Categorize relationships
  related_pairs[, relationship := case_when(
    KINSHIP >= 0.354 ~ "Duplicate/MZ twin",
    KINSHIP >= 0.177 ~ "1st degree",
    KINSHIP >= 0.0884 ~ "2nd degree",
    KINSHIP >= 0.0442 ~ "3rd degree",
    TRUE ~ "Unrelated"
  )]

  log_info("Related pairs found: {nrow(related_pairs)}")
  log_info("  - Duplicates/MZ twins: {sum(related_pairs$relationship == 'Duplicate/MZ twin')}")
  log_info("  - 1st degree: {sum(related_pairs$relationship == '1st degree')}")
  log_info("  - 2nd degree: {sum(related_pairs$relationship == '2nd degree')}")
  log_info("  - 3rd degree: {sum(related_pairs$relationship == '3rd degree')}")

  return(related_pairs)
}

# Function to select unrelated individuals
select_unrelated <- function(related_pairs, sample_ids, phenotype_matrix = NULL) {
  log_info("Selecting unrelated individuals")

  if(nrow(related_pairs) == 0) {
    log_info("No related pairs found, all samples are unrelated")
    return(sample_ids)
  }

  # Get all individuals involved in relationships
  related_individuals <- unique(c(related_pairs$IID1, related_pairs$IID2))
  related_in_data <- intersect(related_individuals, sample_ids)

  log_info("Individuals with relatives in dataset: {length(related_in_data)}")
  # related_in_data used in log message above

  # Strategy: Keep individual with less missing data if phenotype matrix provided
  if(!is.null(phenotype_matrix)) {
    # Calculate missing rate for each sample
    missing_rates <- data.table(
      IID = rownames(phenotype_matrix),
      missing_rate = rowSums(is.na(phenotype_matrix)) / ncol(phenotype_matrix)
    )
  } else {
    missing_rates <- data.table(
      IID = sample_ids,
      missing_rate = 0  # No preference if no phenotype data
    )
  }

  # Iteratively remove related individuals
  samples_to_remove <- character()
  remaining_pairs <- copy(related_pairs)

  while(nrow(remaining_pairs) > 0) {
    # Count how many relationships each individual has
    relationship_counts <- rbind(
      remaining_pairs[, .(IID = IID1, n_relatives = .N), by = IID1][, .(IID, n_relatives)],
      remaining_pairs[, .(IID = IID2, n_relatives = .N), by = IID2][, .(IID, n_relatives)]
    )[, .(n_relatives = sum(n_relatives)), by = IID]

    # Prioritize removal of individuals with most relatives
    setorder(relationship_counts, -n_relatives)

    # Among those with same number of relatives, remove one with more missing data
    top_related <- relationship_counts[n_relatives == max(n_relatives)]$IID

    if(length(top_related) > 1 && !is.null(phenotype_matrix)) {
      # Choose based on missing rate
      missing_info <- missing_rates[IID %in% top_related]
      setorder(missing_info, -missing_rate)
      to_remove <- missing_info[1]$IID
    } else {
      to_remove <- top_related[1]
    }

    # Remove this individual
    samples_to_remove <- c(samples_to_remove, to_remove)

    # Update remaining pairs
    remaining_pairs <- remaining_pairs[IID1 != to_remove & IID2 != to_remove]
  }

  # Get unrelated samples
  unrelated_samples <- setdiff(sample_ids, samples_to_remove)

  log_info("Samples removed due to relatedness: {length(samples_to_remove)}")
  log_info("Unrelated samples remaining: {length(unrelated_samples)}")

  return(list(
    unrelated_samples = unrelated_samples,
    removed_samples = samples_to_remove,
    n_removed = length(samples_to_remove)
  ))
}

# Function to create relationship summary
create_relationship_summary <- function(related_pairs, sample_ids) {
  log_info("Creating relationship summary")

  if(nrow(related_pairs) == 0) {
    return(data.table(
      relationship = character(),
      n_pairs = integer(),
      n_individuals = integer()
    ))
  }

  # Count pairs by relationship type
  pair_summary <- related_pairs[, .(n_pairs = .N), by = relationship]

  # Count unique individuals by relationship type
  individual_summary <- related_pairs[, .(
    individuals = unique(c(IID1, IID2))
  ), by = relationship][, .(n_individuals = length(individuals)), by = relationship]

  # Combine summaries
  summary <- merge(pair_summary, individual_summary, by = "relationship")
  setorder(summary, relationship)

  return(summary)
}

# Function to select best SampleID per FINNGENID based on quality metrics
# Uses hierarchical criteria: missing rate -> n_detected_proteins -> SD NPX -> lexicographic order
select_best_sample_per_finngenid <- function(unrelated_finngenids, sample_mapping, phenotype_matrix,
                                             comprehensive_qc_data = NULL, batch_id, config) {
  log_info("Selecting best SampleID per FINNGENID based on quality metrics")

  # Get all SampleIDs for unrelated FINNGENIDs
  all_candidates <- sample_mapping[FINNGENID %in% unrelated_finngenids & !is.na(FINNGENID)]

  # Filter to samples that exist in phenotype matrix
  all_candidates <- all_candidates[SampleID %in% rownames(phenotype_matrix)]

  # Load comprehensive QC data if available (contains pre-calculated metrics)
  if (is.null(comprehensive_qc_data)) {
    comprehensive_qc_path <- get_output_path("05d", "comprehensive_outliers_list", batch_id, "phenotypes", "tsv", config = config)
    if (file.exists(comprehensive_qc_path)) {
      comprehensive_qc_data <- fread(comprehensive_qc_path)
      log_info("Loaded comprehensive QC data: {nrow(comprehensive_qc_data)} samples")
    } else {
      log_warn("Comprehensive QC data not found - will calculate metrics from phenotype matrix")
      comprehensive_qc_data <- NULL
    }
  }

  # Calculate quality metrics for each SampleID
  quality_metrics <- data.table(
    SampleID = all_candidates$SampleID,
    FINNGENID = all_candidates$FINNGENID
  )

  # Calculate missing rate and number of detected proteins from phenotype matrix
  # Filter phenotype matrix to only candidates for efficiency
  candidate_sampleids <- quality_metrics$SampleID
  samples_in_matrix <- intersect(candidate_sampleids, rownames(phenotype_matrix))

  # Initialize metrics with worst-case values
  quality_metrics[, missing_rate := 1.0]
  quality_metrics[, n_detected_proteins := 0]
  quality_metrics[, sd_npx := Inf]

  # Calculate metrics for samples that exist in matrix
  if (length(samples_in_matrix) > 0) {
    candidate_matrix <- phenotype_matrix[samples_in_matrix, , drop = FALSE]
    n_proteins <- ncol(candidate_matrix)

    # Calculate all metrics efficiently
    missing_counts <- rowSums(is.na(candidate_matrix))
    detected_counts <- n_proteins - missing_counts
    missing_rates <- missing_counts / n_proteins

    # Calculate SD for each sample
    sd_values <- apply(candidate_matrix, 1, function(x) {
      result <- sd(x, na.rm = TRUE)
      if (is.na(result)) Inf else result
    })

    # Update quality metrics for samples in matrix
    quality_metrics[SampleID %in% samples_in_matrix, missing_rate := missing_rates[SampleID]]
    quality_metrics[SampleID %in% samples_in_matrix, n_detected_proteins := detected_counts[SampleID]]
    quality_metrics[SampleID %in% samples_in_matrix, sd_npx := sd_values[SampleID]]
  }

  # If comprehensive QC data available, use pre-calculated SD NPX (more accurate)
  if (!is.null(comprehensive_qc_data) && "QC_technical_sd_npx" %in% names(comprehensive_qc_data)) {
    qc_sd <- comprehensive_qc_data[, .(SampleID, QC_technical_sd_npx)]
    quality_metrics <- merge(quality_metrics, qc_sd, by = "SampleID", all.x = TRUE)
    # Use QC SD if available, otherwise use calculated SD
    quality_metrics[!is.na(QC_technical_sd_npx), sd_npx := QC_technical_sd_npx]
    quality_metrics[, QC_technical_sd_npx := NULL]
  }

  # Handle NA values in SD (set to Inf so they rank last)
  quality_metrics[is.na(sd_npx), sd_npx := Inf]

  # Select best SampleID per FINNGENID using hierarchical criteria:
  # 1. Lower missing rate (better)
  # 2. Higher n_detected_proteins (better)
  # 3. Lower SD NPX (more consistent, better)
  # 4. Lexicographic SampleID order (deterministic tie-breaker)
  setorder(quality_metrics, FINNGENID, missing_rate, -n_detected_proteins, sd_npx, SampleID)

  # Keep first SampleID for each FINNGENID (best quality)
  selected_samples <- quality_metrics[, .SD[1], by = FINNGENID]

  n_duplicates_handled <- sum(quality_metrics[, .N, by = FINNGENID]$N > 1)
  if (n_duplicates_handled > 0) {
    log_info("Selected best SampleID for {n_duplicates_handled} FINNGENIDs with multiple SampleIDs")
    log_info("  Selection criteria: missing_rate (lower) -> n_detected_proteins (higher) -> sd_npx (lower) -> SampleID (lexicographic)")

    # Log example selections
    dup_counts <- quality_metrics[, .N, by = FINNGENID][N > 1]
    if (nrow(dup_counts) > 0) {
      n_examples <- min(3, nrow(dup_counts))
      example_finngenids <- dup_counts[seq_len(n_examples)]$FINNGENID
      for (fgid in example_finngenids) {
        candidates <- quality_metrics[FINNGENID == fgid]
        selected <- selected_samples[FINNGENID == fgid]$SampleID
        log_info("  Example FINNGENID {fgid}:")
        for (j in seq_len(nrow(candidates))) {
          marker <- if (candidates[j]$SampleID == selected) "✓ SELECTED" else "  "
          log_info("    {marker} {candidates[j]$SampleID}: MR={round(candidates[j]$missing_rate, 4)}, N={candidates[j]$n_detected_proteins}, SD={round(candidates[j]$sd_npx, 3)}")
        }
      }
    }
  } else {
    log_info("No duplicate FINNGENIDs found - all samples unique")
  }

  return(selected_samples$SampleID)
}

# Function to filter phenotype matrix
filter_phenotype_matrix <- function(phenotype_matrix, unrelated_samples) {
  log_info("Filtering phenotype matrix to unrelated samples")

  # Filter to unrelated samples
  unrelated_in_matrix <- intersect(rownames(phenotype_matrix), unrelated_samples)
  filtered_matrix <- phenotype_matrix[unrelated_in_matrix, ]

  log_info("Phenotype matrix filtered: {nrow(phenotype_matrix)} -> {nrow(filtered_matrix)} samples")

  return(filtered_matrix)
}

# Function to verify duplicate handling
verify_duplicate_handling <- function(duplicate_finngenids_file, related_pairs, unrelated_samples,
                                      removed_samples, sample_ids, phenotype_matrix = NULL,
                                      kinship_file_path = NULL) {
  log_info(paste(rep("=", 50), collapse = ""))
  log_info("VERIFYING DUPLICATE SAMPLE HANDLING")
  log_info(paste(rep("=", 50), collapse = ""))

  # Load duplicate FINNGENIDs list (optional)
  if(is.null(duplicate_finngenids_file) || duplicate_finngenids_file == "" || !file.exists(duplicate_finngenids_file)) {
    if(is.null(duplicate_finngenids_file) || duplicate_finngenids_file == "") {
      log_info("Duplicate FINNGENIDs file not specified - skipping duplicate verification")
    } else {
      log_warn("Duplicate FINNGENIDs file not found: {duplicate_finngenids_file}")
      log_warn("Skipping duplicate verification")
    }
    return(NULL)
  }

  duplicate_list <- fread(duplicate_finngenids_file)
  expected_duplicates <- duplicate_list$FINNGENID
  n_expected <- length(expected_duplicates)

  log_info("Expected duplicate FINNGENIDs from QC file: {n_expected}")

  # Initialize variables for tracking
  expected_in_kinship_all <- NA
  expected_in_dataset <- NA

  # Check which expected duplicates are in current dataset
  expected_in_dataset <- intersect(expected_duplicates, sample_ids)
  log_info("Expected duplicate FINNGENIDs present in current dataset: {length(expected_in_dataset)}/{n_expected}")

  # Filter to duplicate pairs in kinship data (KINSHIP >= 0.354)
  duplicate_pairs <- related_pairs[relationship == "Duplicate/MZ twin"]
  log_info("Duplicate pairs detected in kinship file (where both samples are in current dataset): {nrow(duplicate_pairs)}")

  # Also check all duplicate pairs in full kinship file (not just those in our data)
  # This helps identify if duplicates were already removed in earlier steps
  if(!is.null(phenotype_matrix) && !is.null(kinship_file_path)) {
    # Load full kinship file to check for all duplicate pairs
    if(file.exists(kinship_file_path)) {
      kinship_full <- fread(kinship_file_path)
      # Standardize column names
      if("ID1" %in% colnames(kinship_full) && !"IID1" %in% colnames(kinship_full)) {
        setnames(kinship_full, c("ID1", "ID2"), c("IID1", "IID2"))
      }
      if("Kinship" %in% colnames(kinship_full) && !"KINSHIP" %in% colnames(kinship_full)) {
        setnames(kinship_full, "Kinship", "KINSHIP")
      }

      all_duplicate_pairs <- kinship_full[KINSHIP >= 0.354]
      duplicate_ids_in_kinship_all <- unique(c(all_duplicate_pairs$IID1, all_duplicate_pairs$IID2))
      expected_in_kinship_all <- intersect(expected_duplicates, duplicate_ids_in_kinship_all)

      log_info("Checking full kinship file for expected duplicates...")
      log_info("Expected duplicate FINNGENIDs found in full kinship file: {length(expected_in_kinship_all)}/{n_expected}")

      if(length(expected_in_kinship_all) > 0 && length(expected_in_dataset) > 0) {
        # Check if any duplicate pairs exist where at least one ID is in our dataset
        pairs_with_one_in_data <- all_duplicate_pairs[IID1 %in% sample_ids | IID2 %in% sample_ids]
        log_info("Duplicate pairs in kinship file where at least one ID is in our dataset: {nrow(pairs_with_one_in_data)}")

        if(nrow(pairs_with_one_in_data) > 0 && nrow(duplicate_pairs) == 0) {
          log_info("Note: Duplicate pairs exist in kinship file but both samples are not in current dataset.")
          log_info("This suggests duplicate samples may have been removed in earlier QC steps.")
        }
      }
    }
  }

  if(nrow(duplicate_pairs) == 0) {
    log_warn("No duplicate pairs found where BOTH samples are in current dataset.")
    log_warn("This may indicate:")
    log_warn("  1. Duplicate samples were already removed in earlier QC steps")
    log_warn("  2. Duplicate pairs don't meet kinship threshold (>= 0.354)")
    log_warn("  3. Duplicate FINNGENIDs don't form duplicate pairs in kinship file")

    # If no pairs to verify, return early but with informative summary
    return(list(
      expected_duplicates = n_expected,
      found_in_kinship = if(!is.na(expected_in_kinship_all[1])) length(expected_in_kinship_all) else NA,
      in_current_dataset = length(expected_in_dataset),
      duplicate_pairs_checked = 0,
      correctly_handled = 0,
      issues_found = 0,
      note = if(!is.na(expected_in_kinship_all[1])) "No duplicate pairs in current dataset to verify (may have been removed in earlier steps)" else "No duplicate pairs in current dataset to verify"
    ))
  }

  # Get all FINNGENIDs involved in duplicate pairs
  duplicate_ids_in_kinship <- unique(c(duplicate_pairs$IID1, duplicate_pairs$IID2))
  log_info("Unique FINNGENIDs in duplicate pairs: {length(duplicate_ids_in_kinship)}")

  # Check overlap with expected duplicates
  # Note: Each expected duplicate should appear in 2 samples, so we expect 2*n_expected IDs
  # But some may have been filtered out already, so we check what's in our data
  duplicate_ids_in_data <- intersect(duplicate_ids_in_kinship, sample_ids)
  log_info("Duplicate FINNGENIDs present in current dataset: {length(duplicate_ids_in_data)}")
  # duplicate_ids_in_data used in log message above

  # Check which expected duplicates are found
  expected_in_kinship <- intersect(expected_duplicates, duplicate_ids_in_kinship)
  expected_not_in_kinship <- setdiff(expected_duplicates, duplicate_ids_in_kinship)

  log_info("Expected duplicates found in kinship file: {length(expected_in_kinship)}/{n_expected}")
  if(length(expected_not_in_kinship) > 0) {
    log_warn("Expected duplicates NOT found in kinship file: {length(expected_not_in_kinship)}")
    log_warn("Missing FINNGENIDs: {paste(head(expected_not_in_kinship, 10), collapse = ', ')}")
    if(length(expected_not_in_kinship) > 10) {
      log_warn("  ... and {length(expected_not_in_kinship) - 10} more")
    }
  }

  # Verify that for each duplicate pair, only one is kept (the one with lower missing rate)
  if(is.null(phenotype_matrix)) {
    log_warn("Phenotype matrix not provided, cannot verify missing rate selection")
    return(list(
      expected_duplicates = n_expected,
      found_in_kinship = length(expected_in_kinship),
      missing_from_kinship = length(expected_not_in_kinship),
      duplicate_pairs_checked = nrow(duplicate_pairs)
    ))
  }

  # Calculate missing rates
  missing_rates <- data.table(
    IID = rownames(phenotype_matrix),
    missing_rate = rowSums(is.na(phenotype_matrix)) / ncol(phenotype_matrix)
  )

  # Check each duplicate pair
  verification_results <- list()
  issues_found <- 0

  for(i in seq_len(nrow(duplicate_pairs))) {
    pair <- duplicate_pairs[i]
    id1 <- pair$IID1
    id2 <- pair$IID2

    # Get missing rates
    mr1 <- missing_rates[IID == id1, missing_rate]
    mr2 <- missing_rates[IID == id2, missing_rate]

    # Check which one is kept
    id1_kept <- id1 %in% unrelated_samples
    id2_kept <- id2 %in% unrelated_samples

    # Determine which should be kept (lower missing rate)
    if(length(mr1) == 0 || length(mr2) == 0) {
      log_warn("Missing rate not found for pair: {id1} - {id2}")
      next
    }

    should_keep_id1 <- mr1 < mr2 || (mr1 == mr2 && id1 < id2)  # Tie-break by ID
    should_keep_id2 <- mr2 < mr1 || (mr1 == mr2 && id2 < id1)

    # Verify correctness
    if(id1_kept && id2_kept) {
      issues_found <- issues_found + 1
      log_error("ISSUE: Both samples kept for duplicate pair: {id1} (MR={round(mr1, 6)}) - {id2} (MR={round(mr2, 6)})")
    } else if(!id1_kept && !id2_kept) {
      issues_found <- issues_found + 1
      log_error("ISSUE: Both samples removed for duplicate pair: {id1} (MR={round(mr1, 6)}) - {id2} (MR={round(mr2, 6)})")
    } else if(id1_kept && !should_keep_id1) {
      issues_found <- issues_found + 1
      log_error("ISSUE: Wrong sample kept for pair: {id1} (MR={round(mr1, 6)}) kept, but {id2} (MR={round(mr2, 6)}) has lower missing rate")
    } else if(id2_kept && !should_keep_id2) {
      issues_found <- issues_found + 1
      log_error("ISSUE: Wrong sample kept for pair: {id2} (MR={round(mr2, 6)}) kept, but {id1} (MR={round(mr1, 6)}) has lower missing rate")
    } else {
      # Correct handling
      kept_id <- if(id1_kept) id1 else id2
      removed_id <- if(id1_kept) id2 else id1
      kept_mr <- if(id1_kept) mr1 else mr2
      removed_mr <- if(id1_kept) mr2 else mr1
      log_debug("✓ Correct: {kept_id} kept (MR={round(kept_mr, 6)}), {removed_id} removed (MR={round(removed_mr, 6)})")
      # kept_id, removed_id, kept_mr, removed_mr used in log_debug above
    }

    verification_results[[i]] <- data.table(
      FINNGENID1 = id1,
      FINNGENID2 = id2,
      MR1 = mr1,
      MR2 = mr2,
      ID1_kept = id1_kept,
      ID2_kept = id2_kept,
      correct = (id1_kept && should_keep_id1) || (id2_kept && should_keep_id2)
    )
  }

  verification_dt <- rbindlist(verification_results)

  # Summary statistics
  n_pairs_checked <- nrow(verification_dt)
  n_correct <- sum(verification_dt$correct, na.rm = TRUE)
  n_issues <- n_pairs_checked - n_correct

  log_info(paste(rep("=", 50), collapse = ""))
  log_info("DUPLICATE VERIFICATION SUMMARY")
  log_info(paste(rep("=", 50), collapse = ""))
  log_info("Expected duplicate FINNGENIDs (from QC): {n_expected}")
  log_info("Duplicate pairs found in kinship file: {nrow(duplicate_pairs)}")
  log_info("Expected duplicates found in kinship: {length(expected_in_kinship)}/{n_expected}")
  log_info("Duplicate pairs verified: {n_pairs_checked}")
  log_info("Correctly handled: {n_correct}/{n_pairs_checked}")
  if(n_issues > 0) {
    log_error("ISSUES FOUND: {n_issues} duplicate pairs not handled correctly!")
  } else {
    log_info("✓ All duplicate pairs correctly handled (lowest missing rate kept)")
  }
  log_info(paste(rep("=", 50), collapse = ""))

  return(list(
    expected_duplicates = n_expected,
    found_in_kinship = length(expected_in_kinship),
    missing_from_kinship = length(expected_not_in_kinship),
    duplicate_pairs_checked = n_pairs_checked,
    correctly_handled = n_correct,
    issues_found = n_issues,
    verification_table = verification_dt
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
  phenotype_matrix_path <- get_output_path("09", "phenotype_matrix", batch_id, "phenotypes", config = config)
  finngenid_matrix_path <- get_output_path("09", "phenotype_matrix_finngenid", batch_id, "phenotypes", config = config)

  if (!file.exists(phenotype_matrix_path)) {
    stop("Phenotype matrix not found: {phenotype_matrix_path}. Please run step 09 first.")
  }

  phenotype_matrix <- readRDS(phenotype_matrix_path)
  finngenid_matrix <- NULL
  if(file.exists(finngenid_matrix_path)) {
    finngenid_matrix <- readRDS(finngenid_matrix_path)
    log_info("Loaded FINNGENID-indexed matrix: {nrow(finngenid_matrix)} samples")
  }
  log_info("Loaded phenotype matrix: {nrow(phenotype_matrix)} samples x {ncol(phenotype_matrix)} proteins")

  # Load kinship matrix
  kinship <- load_kinship(config$covariates$kinship_file)

  # Get sample IDs (use FINNGENIDs if available)
  if(!is.null(finngenid_matrix)) {
    sample_ids <- rownames(finngenid_matrix)
    matrix_for_filtering <- finngenid_matrix
    log_info("Using FINNGENID-indexed matrix for kinship filtering")
  } else {
    sample_ids <- rownames(phenotype_matrix)
    matrix_for_filtering <- phenotype_matrix
    log_info("Using sample ID matrix for kinship filtering")
  }

  # Identify related pairs
  related_pairs <- identify_related_pairs(
    kinship,
    threshold = config$parameters$pqtl$kinship_threshold
  )

  # Filter to pairs where both individuals are in our data
  related_pairs_in_data <- related_pairs[IID1 %in% sample_ids & IID2 %in% sample_ids]

  log_info("Related pairs in our data: {nrow(related_pairs_in_data)}")

  # Select unrelated individuals
  unrelated_result <- select_unrelated(
    related_pairs_in_data,
    sample_ids,
    matrix_for_filtering
  )

  # Verify duplicate handling (optional - requires duplicate FINNGENIDs file)
  # Check if duplicate file is specified in config, otherwise set to NULL
  duplicate_file <- tryCatch(
    config$covariates$duplicate_finngenids_file,
    error = function(e) NULL
  )

  # If not in config, try to find it in a standard location (optional)
  if (is.null(duplicate_file) || duplicate_file == "") {
    # Try to find duplicate file in output directory (if it exists from QC steps)
    duplicate_file_candidate <- get_output_path("00", "duplicate_finngenids", batch_id, "qc", "tsv", config = config)
    if (file.exists(duplicate_file_candidate)) {
      duplicate_file <- duplicate_file_candidate
      log_info("Found duplicate FINNGENIDs file: {duplicate_file}")
    } else {
      log_info("Duplicate FINNGENIDs file not specified or found - skipping duplicate verification")
      duplicate_file <- NULL
    }
  }

  verification_result <- verify_duplicate_handling(
    duplicate_file,
    related_pairs_in_data,
    unrelated_result$unrelated_samples,
    unrelated_result$removed_samples,
    sample_ids,
    matrix_for_filtering,
    config$covariates$kinship_file
  )
  # verification_result contains verification statistics (may be used for reporting)

  # Create relationship summary
  relationship_summary <- create_relationship_summary(related_pairs_in_data, sample_ids)

  # Filter phenotype matrices
  # CRITICAL: If using FINNGENID-indexed matrix, unrelated_samples are FINNGENIDs
  # We need to filter the FINNGENID matrix, and optionally convert back to SampleID matrix
  if(!is.null(finngenid_matrix)) {
    # Filter FINNGENID-indexed matrix (unrelated_samples are FINNGENIDs)
    finngenid_unrelated <- filter_phenotype_matrix(
      finngenid_matrix,
      unrelated_result$unrelated_samples
    )

    # Convert FINNGENID-unrelated samples back to SampleID format for phenotype_matrix
    # CRITICAL: Select ONE SampleID per FINNGENID based on quality metrics
    # This ensures SampleID matrix matches FINNGENID matrix dimensions (one sample per person)
    sample_mapping_path <- get_output_path("00", "sample_mapping", batch_id, "qc", "tsv", config = config)
    if (file.exists(sample_mapping_path)) {
      sample_mapping <- fread(sample_mapping_path)

      # Select best SampleID per FINNGENID using quality-based selection
      # Load comprehensive QC data for pre-calculated metrics
      comprehensive_qc_path <- get_output_path("05d", "comprehensive_outliers_list", batch_id, "phenotypes", "tsv", config = config)
      comprehensive_qc_data <- NULL
      if (file.exists(comprehensive_qc_path)) {
        comprehensive_qc_data <- fread(comprehensive_qc_path)
      }

      unrelated_sampleids <- select_best_sample_per_finngenid(
        unrelated_result$unrelated_samples,
        sample_mapping,
        phenotype_matrix,
        comprehensive_qc_data,
        batch_id,
        config
      )

      # Filter SampleID-indexed matrix using selected SampleIDs (one per FINNGENID)
      phenotype_unrelated <- filter_phenotype_matrix(
        phenotype_matrix,
        unrelated_sampleids
      )
      log_info("Selected {length(unrelated_sampleids)} SampleIDs (one per FINNGENID) from {length(unrelated_result$unrelated_samples)} unrelated FINNGENIDs")
    } else {
      log_warn("Sample mapping not found - cannot convert FINNGENIDs to SampleIDs")
      log_warn("Using empty SampleID matrix (will use FINNGENID matrix only)")
      phenotype_unrelated <- phenotype_matrix[character(0), ]  # Empty matrix with same structure
    }
  } else {
    # No FINNGENID matrix - use SampleIDs directly
    phenotype_unrelated <- filter_phenotype_matrix(
      phenotype_matrix,
      unrelated_result$unrelated_samples
    )
    finngenid_unrelated <- NULL
  }

  # Save outputs
  log_info("Saving kinship filtering results")

  # Save filtered matrices using batch-aware paths
  phenotype_unrelated_path <- get_output_path(step_num, "phenotype_matrix_unrelated", batch_id, "phenotypes", config = config)
  ensure_output_dir(phenotype_unrelated_path)
  saveRDS(phenotype_unrelated, phenotype_unrelated_path)

  if(!is.null(finngenid_unrelated)) {
    finngenid_unrelated_path <- get_output_path(step_num, "phenotype_matrix_finngenid_unrelated", batch_id, "phenotypes", config = config)
    ensure_output_dir(finngenid_unrelated_path)
    saveRDS(finngenid_unrelated, finngenid_unrelated_path)
  }

  # Save sample lists
  samples_unrelated_path <- get_output_path(step_num, "samples_unrelated", batch_id, "phenotypes", "txt", config = config)
  samples_removed_path <- get_output_path(step_num, "samples_related_removed", batch_id, "phenotypes", "txt", config = config)
  ensure_output_dir(samples_unrelated_path)
  ensure_output_dir(samples_removed_path)
  fwrite(data.table(IID = unrelated_result$unrelated_samples), samples_unrelated_path, col.names = FALSE)
  fwrite(data.table(IID = unrelated_result$removed_samples), samples_removed_path, col.names = FALSE)

  # Save relationship information with FINNGENID mapping
  # Source path utilities for FINNGENID mapping
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
  batch_id <- Sys.getenv("PIPELINE_BATCH_ID", config$batch$default_batch_id %||% "batch_01")

  # Add FINNGENID1 and FINNGENID2 to related_pairs
  # Load sample mapping
  sample_mapping_path <- get_output_path("00", "sample_mapping", batch_id, "qc", "tsv", config = config)
  if (file.exists(sample_mapping_path)) {
    sample_mapping <- fread(sample_mapping_path)
    # Map IID1 to FINNGENID1
    related_pairs_with_fgid <- merge(related_pairs_in_data,
                                     sample_mapping[, .(IID1 = SampleID, FINNGENID1 = FINNGENID)],
                                     by = "IID1", all.x = TRUE)
    # Map IID2 to FINNGENID2
    related_pairs_with_fgid <- merge(related_pairs_with_fgid,
                                     sample_mapping[, .(IID2 = SampleID, FINNGENID2 = FINNGENID)],
                                     by = "IID2", all.x = TRUE)
    # Reorder columns: IID1, FINNGENID1, IID2, FINNGENID2, ...
    col_order <- c("IID1", "FINNGENID1", "IID2", "FINNGENID2",
                   setdiff(names(related_pairs_with_fgid), c("IID1", "FINNGENID1", "IID2", "FINNGENID2")))
    setcolorder(related_pairs_with_fgid, col_order)
  } else {
    log_warn("Sample mapping file not found, saving related_pairs without FINNGENID")
    related_pairs_with_fgid <- related_pairs_in_data
  }

  related_pairs_path <- get_output_path(step_num, "related_pairs", batch_id, "phenotypes", "tsv", config = config)
  relationship_summary_path <- get_output_path(step_num, "relationship_summary", batch_id, "phenotypes", "tsv", config = config)
  ensure_output_dir(related_pairs_path)
  ensure_output_dir(relationship_summary_path)
  fwrite(related_pairs_with_fgid, related_pairs_path, sep = "\t")
  fwrite(relationship_summary, relationship_summary_path, sep = "\t")

  # Handle aggregation if enabled
  if (aggregate_output && multi_batch_mode && !is.null(finngenid_unrelated)) {
    log_info("Aggregation enabled: Attempting to filter batch 1 and create aggregate output")

    # Check if batch 1 processed data exists
    other_batch_id <- if (batch_id == "batch_02") "batch_01" else "batch_02"
    batch1_file <- get_output_path(step_num, "phenotype_matrix_finngenid_unrelated", other_batch_id, "phenotypes", config = config)

    if (file.exists(batch1_file)) {
      log_info("Found batch 1 kinship-filtered data: Creating aggregate output")

      # Load batch 1 data
      batch1_unrelated <- readRDS(batch1_file)

      # Combine batch 1 and batch 2 (on FINNGENID, common proteins)
      common_proteins <- intersect(colnames(batch1_unrelated), colnames(finngenid_unrelated))
      if (length(common_proteins) > 100) {
        batch1_subset <- batch1_unrelated[, common_proteins, drop = FALSE]
        batch2_subset <- finngenid_unrelated[, common_proteins, drop = FALSE]

        # Check for duplicate FINNGENIDs (use batch 2 as reference)
        common_finngenids <- intersect(rownames(batch1_subset), rownames(batch2_subset))
        if (length(common_finngenids) > 0) {
          log_warn("Found {length(common_finngenids)} FINNGENIDs in both batches, using batch 2 data")
          batch1_subset <- batch1_subset[!rownames(batch1_subset) %in% common_finngenids, , drop = FALSE]
        }

        # Combine
        aggregate_unrelated <- rbind(batch1_subset, batch2_subset)
        log_info("Aggregate unrelated matrix: {nrow(aggregate_unrelated)} samples x {ncol(aggregate_unrelated)} proteins")

        # Save aggregate output using batch-aware paths
        aggregate_dir <- file.path(config$output$base_dir, "phenotypes", "aggregate")
        dir.create(aggregate_dir, recursive = TRUE, showWarnings = FALSE)

        aggregate_matrix_path <- file.path(aggregate_dir, "aggregate_phenotype_matrix_finngenid_unrelated.rds")
        aggregate_samples_path <- file.path(aggregate_dir, "aggregate_samples_unrelated.txt")

        saveRDS(aggregate_unrelated, aggregate_matrix_path)
        fwrite(data.table(FINNGENID = rownames(aggregate_unrelated)), aggregate_samples_path, col.names = FALSE)
      } else {
        log_warn("Too few common proteins for aggregation")
      }
    } else {
      log_warn("Batch 1 kinship-filtered data not found. Aggregate output skipped.")
      log_warn("Expected file: {batch1_file}")
    }
  }

  # Create summary report
  summary_report <- data.table(
    metric = c(
      "Total samples",
      "Related individuals",
      "Samples removed",
      "Unrelated samples",
      "Duplicate/MZ pairs",
      "1st degree pairs",
      "2nd degree pairs",
      "3rd degree pairs"
    ),
    value = c(
      length(sample_ids),
      length(unique(c(related_pairs_in_data$IID1, related_pairs_in_data$IID2))),
      unrelated_result$n_removed,
      length(unrelated_result$unrelated_samples),
      sum(related_pairs_in_data$relationship == "Duplicate/MZ twin"),
      sum(related_pairs_in_data$relationship == "1st degree"),
      sum(related_pairs_in_data$relationship == "2nd degree"),
      sum(related_pairs_in_data$relationship == "3rd degree")
    )
  )

  fwrite(summary_report, , sep = "\t")

  # Print summary
  cat("\n=== KINSHIP FILTERING SUMMARY ===\n")
  cat("Kinship threshold:", config$parameters$pqtl$kinship_threshold, "\n")
  cat("Total samples:", length(sample_ids), "\n")
  cat("Related pairs in data:", nrow(related_pairs_in_data), "\n")
  print(relationship_summary)
  cat("\nSamples removed:", unrelated_result$n_removed, "\n")
  cat("Unrelated samples:", length(unrelated_result$unrelated_samples), "\n")
  cat("\nFiltered phenotype matrix:", nrow(phenotype_unrelated), "x", ncol(phenotype_unrelated), "\n")
  cat("Results saved to: ../output/phenotypes/\n")

  log_info("Kinship filtering completed")

  return(list(
    unrelated_samples = unrelated_result$unrelated_samples,
    removed_samples = unrelated_result$removed_samples,
    phenotype_unrelated = phenotype_unrelated,
    finngenid_unrelated = finngenid_unrelated,
    relationship_summary = relationship_summary
  ))
}

# Run if executed directly
if (!interactive()) {
  result <- main()
}






