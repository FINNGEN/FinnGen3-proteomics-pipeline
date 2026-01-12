#!/usr/bin/env Rscript

#################################################
# Script: 10_kinship_filtering.R
# Author: Reza Jabal, PhD (rjabal@broadinstitute.org)
# Description: Filter related individuals using kinship matrix
#              Refactored version
#              Refactored version
# Date: December 2025
#################################################

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

# Set up logging
log_appender(appender_file())
log_info("Starting kinship filtering")

# Load configuration
config <- read_yaml("")

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

  # Load duplicate FINNGENIDs list
  if(!file.exists(duplicate_finngenids_file)) {
    log_warn("Duplicate FINNGENIDs file not found: {duplicate_finngenids_file}")
    log_warn("Skipping duplicate verification")
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

  # Determine input prefix based on aggregation mode
  if (aggregate_output) {
    input_prefix <- "batch2_"
    output_prefix <- "batch2_"
    log_info("Aggregation mode: Using batch 2 inputs/outputs with 'batch2_' prefix")
  } else {
    input_prefix <- "11_"
    output_prefix <- "12_"
    log_info("Single-batch mode: Using standard '11_'/'12_' prefixes")
  }

  # Load data from previous steps
  log_info("Loading data from previous steps")
  phenotype_matrix <- readRDS(paste0(, input_prefix, "phenotype_matrix.rds"))
  finngenid_matrix <- NULL
  finngenid_file <- paste0(, input_prefix, "phenotype_matrix_finngenid.rds")
  if(file.exists(finngenid_file)) {
    finngenid_matrix <- readRDS(finngenid_file)
  }

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

  # Verify duplicate handling
  duplicate_file <-
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
  phenotype_unrelated <- filter_phenotype_matrix(
    phenotype_matrix,
    unrelated_result$unrelated_samples
  )

  if(!is.null(finngenid_matrix)) {
    finngenid_unrelated <- filter_phenotype_matrix(
      finngenid_matrix,
      unrelated_result$unrelated_samples
    )
  } else {
    finngenid_unrelated <- NULL
  }

  # Save outputs
  log_info("Saving kinship filtering results")

  # Save filtered matrices
  saveRDS(phenotype_unrelated,
          paste0(, output_prefix, "phenotype_matrix_unrelated.rds"))
  if(!is.null(finngenid_unrelated)) {
    saveRDS(finngenid_unrelated,
            paste0(, output_prefix, "phenotype_matrix_finngenid_unrelated.rds"))
  }

  # Save sample lists
  fwrite(data.table(IID = unrelated_result$unrelated_samples),
         paste0(, output_prefix, "samples_unrelated.txt"),
         col.names = FALSE)
  fwrite(data.table(IID = unrelated_result$removed_samples),
         paste0(, output_prefix, "samples_related_removed.txt"),
         col.names = FALSE)

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

  fwrite(related_pairs_with_fgid,
         paste0(, output_prefix, "related_pairs.tsv"),
         sep = "\t")
  fwrite(relationship_summary,
         paste0(, output_prefix, "relationship_summary.tsv"),
         sep = "\t")

  # Handle aggregation if enabled
  if (aggregate_output && multi_batch_mode && !is.null(finngenid_unrelated)) {
    log_info("Aggregation enabled: Attempting to filter batch 1 and create aggregate output")

    # Check if batch 1 processed data exists
    batch1_file <-

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

        # Save aggregate output
        saveRDS(aggregate_unrelated,
                )
        fwrite(data.table(FINNGENID = rownames(aggregate_unrelated)),
               ,
               col.names = FALSE)
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






