#!/usr/bin/env Rscript

#################################################
# Script: 00_data_loader.R
# Author: Reza Jabal, PhD (rjabal@broadinstitute.org)
# Description: Load pre-filtered NPX matrix and prepare analysis-ready data
#              Starting point: NPX matrix with AG samples already removed
#              Supports both single-batch and multi-batch modes
# Date: December 2025
#################################################

# Load required libraries
suppressPackageStartupMessages({
  library(arrow)
  library(data.table)
  library(tidyverse)
  library(yaml)
  library(logger)
})

# Source path utilities
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
}, error = function(e) getwd())
source(file.path(script_dir, "path_utils.R"))

# Get config path from environment or use default
config_file <- Sys.getenv("PIPELINE_CONFIG", "")
if (config_file == "") {
  stop("PIPELINE_CONFIG environment variable not set. Please provide path to config file.")
}
config <- read_yaml(config_file)

# Get batch context
batch_id <- Sys.getenv("PIPELINE_BATCH_ID", config$batch$default_batch_id %||% "batch_01")
step_num <- get_step_number()

# Set up logging with batch-aware path
log_path <- get_log_path(step_num, batch_id, config)
ensure_output_dir(log_path)
log_appender(appender_file(log_path))
log_info("Starting data loader for batch: {batch_id}")

# Function to load NPX matrix from parquet
load_npx_matrix_from_parquet <- function(parquet_file) {
  log_info("Loading NPX matrix from parquet: {parquet_file}")

  if (!file.exists(parquet_file)) {
    stop("NPX matrix parquet file not found: ", parquet_file)
  }

  # Read parquet file
  dt <- read_parquet(parquet_file)
  setDT(dt)

  # Check if SampleID column exists
  if (!"SampleID" %in% names(dt)) {
    # Try first column as SampleID
    if (ncol(dt) > 0) {
      setnames(dt, names(dt)[1], "SampleID")
      log_info("Using first column as SampleID")
    } else {
      stop("Could not identify SampleID column in parquet file")
    }
  }

  # Extract SampleID and convert to matrix
  sample_ids <- dt$SampleID
  protein_cols <- setdiff(names(dt), "SampleID")

  # Convert to matrix
  npx_matrix <- as.matrix(dt[, ..protein_cols])
  rownames(npx_matrix) <- sample_ids

  log_info("Loaded NPX matrix: {nrow(npx_matrix)} samples x {ncol(npx_matrix)} proteins")
  log_info("Missing values: {sum(is.na(npx_matrix))} ({round(sum(is.na(npx_matrix))/(nrow(npx_matrix)*ncol(npx_matrix))*100, 2)}%)")

  return(npx_matrix)
}

# Function to load metadata
load_metadata <- function(metadata_file) {
  if (is.null(metadata_file) || length(metadata_file) == 0 || !is.character(metadata_file)) {
    stop("metadata_file must be a non-empty character string")
  }

  if (!file.exists(metadata_file)) {
    stop("metadata_file does not exist: ", metadata_file)
  }

  log_info("Loading metadata from: {metadata_file}")

  # Explicitly use tab separator for metadata files
  metadata <- fread(metadata_file, sep = "\t")
  log_info("Metadata dimensions: {nrow(metadata)} rows x {ncol(metadata)} columns")

  return(metadata)
}

# Function to load bridging samples info
load_bridging_samples <- function(bridging_file) {
  if (is.null(bridging_file) || length(bridging_file) == 0 || !is.character(bridging_file)) {
    log_info("No bridging samples file provided, returning empty data.table")
    return(data.table())
  }

  if (!file.exists(bridging_file)) {
    log_warn("Bridging samples file does not exist: {bridging_file}, returning empty data.table")
    return(data.table())
  }

  log_info("Loading bridging samples from: {bridging_file}")

  bridging <- fread(bridging_file, sep = ";")
  log_info("Bridging samples: {nrow(bridging)}")

  return(bridging)
}

# Function to load EA5 to FG3 FINNGENID mapping
load_bridging_finngenid_map <- function(mapping_file) {
  if (is.null(mapping_file) || !file.exists(mapping_file)) {
    return(NULL)
  }

  log_info("Loading EA5 to FG3 FINNGENID mapping from: {mapping_file}")

  # Read the mapping file (TSV format with quoted headers)
  mapping <- fread(mapping_file, sep = "\t")

  # Clean column names (remove quotes if present)
  setnames(mapping, names(mapping), gsub('"', '', names(mapping)))

  # Expected columns: Tube_ID, FGID, PSEUDO_ID
  if (!all(c("PSEUDO_ID", "FGID") %in% names(mapping))) {
    log_warn("Expected columns PSEUDO_ID and FGID not found in mapping file")
    log_warn("Found columns: {paste(names(mapping), collapse=', ')}")
    return(NULL)
  }

  # Create clean mapping: PSEUDO_ID -> FINNGENID
  finngenid_map <- mapping[, .(SAMPLE_ID = PSEUDO_ID, FINNGENID = FGID)]
  finngenid_map <- finngenid_map[!is.na(SAMPLE_ID) & !is.na(FINNGENID)]

  # Filter to only include FINNGENIDs with "FG" prefix (universal rule)
  n_before <- nrow(finngenid_map)
  finngenid_map <- finngenid_map[grepl("^FG", FINNGENID)]
  n_after <- nrow(finngenid_map)

  if (n_before > n_after) {
    log_warn("Filtered out {n_before - n_after} FINNGENIDs that don't start with 'FG' prefix")
  }

  log_info("Loaded FINNGENID mapping for {nrow(finngenid_map)} bridging samples (all with 'FG' prefix)")

  return(finngenid_map)
}

# Function to create sample mapping
create_sample_mapping <- function(sample_ids, metadata, bridging_finngenid_map = NULL) {
  log_info("Creating sample mapping")

  # Create base mapping table
  mapping <- data.table(
    SampleID = sample_ids,
    has_metadata = sample_ids %in% metadata$SAMPLE_ID
  )

  # Add FINNGENID from metadata
  mapping <- merge(
    mapping,
    metadata[, .(SAMPLE_ID, FINNGENID, COHORT_FINNGENID, BIOBANK_PLASMA)],
    by.x = "SampleID",
    by.y = "SAMPLE_ID",
    all.x = TRUE
  )

  # Supplement with FINNGENIDs from EA5→FG3 mapping for bridging samples
  if (!is.null(bridging_finngenid_map)) {
    log_info("Supplementing sample mapping with FINNGENIDs from EA5→FG3 mapping")

    # Merge the EA5→FG3 mapping for samples without FINNGENID
    mapping <- merge(
      mapping,
      bridging_finngenid_map,
      by.x = "SampleID",
      by.y = "SAMPLE_ID",
      all.x = TRUE,
      suffixes = c("", "_ea5")
    )

    # Use EA5 mapping FINNGENID where regular metadata FINNGENID is missing
    n_before <- sum(!is.na(mapping$FINNGENID))
    mapping[is.na(FINNGENID) & !is.na(FINNGENID_ea5), FINNGENID := FINNGENID_ea5]
    mapping[, FINNGENID_ea5 := NULL]  # Remove temporary column
    n_after <- sum(!is.na(mapping$FINNGENID))

    log_info("Supplemented {n_after - n_before} samples with FINNGENIDs from EA5→FG3 mapping")
  }

  # Flag sample types (AG samples already removed, so only FinnGen and Bridging)
  mapping[, sample_type := case_when(
    grepl("CONTROL", SampleID) ~ "Control",
    grepl("EA5_OLI", SampleID) ~ "Bridging",
    !is.na(FINNGENID) ~ "FinnGen",
    TRUE ~ "Unknown"
  )]

  # Summary statistics
  log_info("Sample mapping summary:")
  log_info("  - Total samples: {nrow(mapping)}")
  log_info("  - FinnGen samples: {sum(mapping$sample_type == 'FinnGen')}")
  log_info("  - Bridging samples: {sum(mapping$sample_type == 'Bridging')}")
  log_info("  - Unknown samples: {sum(mapping$sample_type == 'Unknown')}")

  return(mapping)
}

# Function to validate sample mapping
validate_mapping <- function(mapping, npx_matrix) {
  log_info("Validating sample mapping")

  # Check if all matrix samples are mapped
  matrix_samples <- rownames(npx_matrix)
  unmapped <- setdiff(matrix_samples, mapping$SampleID)

  if(length(unmapped) > 0) {
    log_warn("{length(unmapped)} samples in matrix not found in mapping")
  }

  # Check for duplicate FINNGENIDs
  dup_finngen <- mapping[!is.na(FINNGENID), .(n = .N), by = FINNGENID][n > 1]

  if(nrow(dup_finngen) > 0) {
    log_warn("{nrow(dup_finngen)} duplicate FINNGENIDs found")
  }

  # Validation summary
  validation <- list(
    unmapped_samples = unmapped,
    duplicate_finngenids = dup_finngen,
    mapping_complete = length(unmapped) == 0 & nrow(dup_finngen) == 0
  )

  return(validation)
}

# Main execution
main <- function() {

  # Get input paths from config
  npx_matrix_file <- get_batch_input_path("npx_matrix_file", batch_id, config)
  metadata_file <- get_batch_input_path("metadata_file", batch_id, config)
  bridging_samples_file <- get_batch_input_path("bridging_samples_file", batch_id, config)
  bridging_samples_finngenid_map_file <- get_batch_input_path("bridging_samples_finngenid_map", batch_id, config)

  # Validate required inputs
  if (is.null(npx_matrix_file)) {
    stop("npx_matrix_file is required but not found in config for batch: ", batch_id)
  }
  if (is.null(metadata_file)) {
    stop("metadata_file is required but not found in config for batch: ", batch_id)
  }

  log_info("Input configuration:")
  log_info("  NPX matrix file: {npx_matrix_file}")
  log_info("  Metadata file: {metadata_file}")
  if (!is.null(bridging_samples_file)) {
    log_info("  Bridging samples file: {bridging_samples_file}")
  }
  if (!is.null(bridging_samples_finngenid_map_file)) {
    log_info("  Bridging FINNGENID mapping: {bridging_samples_finngenid_map_file}")
  }

  # Load data
  log_info("Loading data...")
  npx_matrix <- load_npx_matrix_from_parquet(npx_matrix_file)
  metadata <- load_metadata(metadata_file)

  # Load optional bridging sample information
  bridging_samples <- load_bridging_samples(bridging_samples_file)
  bridging_finngenid_map <- NULL
  if (!is.null(bridging_samples_finngenid_map_file)) {
    bridging_finngenid_map <- load_bridging_finngenid_map(bridging_samples_finngenid_map_file)
  }

  # Create sample mapping
  log_info("Creating sample mapping...")
  sample_ids <- rownames(npx_matrix)
  sample_mapping <- create_sample_mapping(sample_ids, metadata, bridging_finngenid_map)

  # Validate mapping
  validation <- validate_mapping(sample_mapping, npx_matrix)

  # Since AG samples are already removed, the analysis-ready matrix is the same as input
  # But we filter to ensure only FinnGen and Bridging samples are included
  analysis_samples <- sample_mapping[
    sample_type %in% c("FinnGen", "Bridging") &
    SampleID %in% rownames(npx_matrix)
  ]

  log_info("Analysis-ready samples: {nrow(analysis_samples)}")
  analysis_matrix <- npx_matrix[analysis_samples$SampleID, ]

  # Save outputs
  log_info("Saving outputs...")

  # Save metadata
  metadata_path <- get_output_path(step_num, "metadata", batch_id, "qc", config = config)
  ensure_output_dir(metadata_path)
  saveRDS(metadata, metadata_path)

  # Save sample mapping
  sample_mapping_path <- get_output_path(step_num, "sample_mapping", batch_id, "qc", config = config)
  sample_mapping_tsv_path <- get_output_path(step_num, "sample_mapping", batch_id, "qc", "tsv", config = config)
  ensure_output_dir(sample_mapping_path)
  ensure_output_dir(sample_mapping_tsv_path)
  saveRDS(sample_mapping, sample_mapping_path)
  fwrite(sample_mapping, sample_mapping_tsv_path, sep = "\t")

  # Save validation
  mapping_validation_path <- get_output_path(step_num, "mapping_validation", batch_id, "qc", config = config)
  ensure_output_dir(mapping_validation_path)
  saveRDS(validation, mapping_validation_path)

  # Save analysis-ready samples
  analysis_samples_path <- get_output_path(step_num, "analysis_samples", batch_id, "qc", config = config)
  analysis_samples_tsv_path <- get_output_path(step_num, "analysis_samples", batch_id, "qc", "tsv", config = config)
  ensure_output_dir(analysis_samples_path)
  ensure_output_dir(analysis_samples_tsv_path)
  saveRDS(analysis_samples, analysis_samples_path)
  fwrite(analysis_samples, analysis_samples_tsv_path, sep = "\t")

  # Save analysis-ready matrices
  npx_matrix_analysis_ready_path <- get_output_path(step_num, "npx_matrix_analysis_ready", batch_id, "qc", config = config)
  ensure_output_dir(npx_matrix_analysis_ready_path)
  saveRDS(analysis_matrix, npx_matrix_analysis_ready_path)

  # Create minimal samples_data_raw from metadata for downstream steps (needed for plate info)
  # Since we start from pre-filtered matrix, we create a minimal long-format data table
  log_info("Creating minimal samples_data_raw from metadata...")
  samples_data_raw <- data.table(
    SampleID = rownames(analysis_matrix)
  )
  # Add PlateID from metadata where available
  if("PlateID" %in% names(metadata)) {
    plate_map <- metadata[, .(SAMPLE_ID, PlateID)]
    samples_data_raw <- merge(
      samples_data_raw,
      plate_map,
      by.x = "SampleID",
      by.y = "SAMPLE_ID",
      all.x = TRUE
    )
  } else {
    samples_data_raw[, PlateID := NA_character_]
  }
  # Add NPX column (mean NPX per sample for compatibility)
  sample_means <- rowMeans(analysis_matrix, na.rm = TRUE)
  samples_data_raw[, NPX := sample_means[SampleID]]
  samples_data_raw[, SampleQC := "PASS"]  # Default QC pass
  samples_data_raw[, AssayQC := "PASS"]   # Default QC pass

  samples_data_raw_path <- get_output_path(step_num, "samples_data_raw", batch_id, "qc", config = config)
  ensure_output_dir(samples_data_raw_path)
  saveRDS(samples_data_raw, samples_data_raw_path)
  log_info("Saved minimal samples_data_raw: {samples_data_raw_path} ({nrow(samples_data_raw)} samples)")

  # Save duplicate FINNGENIDs if any
  if(nrow(validation$duplicate_finngenids) > 0) {
    duplicate_finngenids_path <- get_output_path(step_num, "duplicate_finngenids", batch_id, "qc", "tsv", config = config)
    ensure_output_dir(duplicate_finngenids_path)
    fwrite(validation$duplicate_finngenids, duplicate_finngenids_path, sep = "\t")
  }

  # Identify and save bridge samples (needed for normalization step)
  log_info("Identifying bridge samples...")
  bridge_samples <- sample_mapping[sample_type == "Bridging"]$SampleID
  if(length(bridge_samples) > 0) {
    log_info("Found {length(bridge_samples)} bridge samples")
    bridge_result <- list(
      bridge_ids = bridge_samples,
      bridge_summary = sample_mapping[sample_type == "Bridging"],
      n_bridge = length(bridge_samples)
    )
    bridge_result_path <- get_output_path(step_num, "bridge_samples_identified", batch_id, "normalized", config = config)
    ensure_output_dir(bridge_result_path)
    saveRDS(bridge_result, bridge_result_path)
    log_info("Saved bridge sample identification: {bridge_result_path}")
  } else {
    log_warn("No bridge samples found in the data")
    # Create empty bridge result for downstream steps
    bridge_result <- list(
      bridge_ids = character(0),
      bridge_summary = data.table(),
      n_bridge = 0
    )
    bridge_result_path <- get_output_path(step_num, "bridge_samples_identified", batch_id, "normalized", config = config)
    ensure_output_dir(bridge_result_path)
    saveRDS(bridge_result, bridge_result_path)
  }

  # Print summary
  cat("\n=== DATA LOADER SUMMARY ===\n")
  cat("NPX matrix loaded: ", nrow(npx_matrix), " samples x ", ncol(npx_matrix), " proteins\n", sep = "")
  cat("Analysis-ready samples: ", nrow(analysis_samples), "\n", sep = "")
  cat("Sample types:\n")
  print(table(sample_mapping$sample_type))
  if(nrow(validation$duplicate_finngenids) > 0) {
    cat("Duplicate FINNGENIDs: ", nrow(validation$duplicate_finngenids), "\n", sep = "")
  }
  cat("\nResults saved to: output/qc/", if(batch_id != config$batch$default_batch_id ||
    tryCatch(isTRUE(config$parameters$normalization$multi_batch_mode), error = function(e) FALSE))
    paste0(batch_id, "/") else "", "\n", sep = "")

  log_info("Data loader completed successfully")

  return(list(
    npx_matrix = npx_matrix,
    metadata = metadata,
    sample_mapping = sample_mapping,
    validation = validation,
    analysis_samples = analysis_samples,
    analysis_matrix = analysis_matrix
  ))
}

# Run if executed directly
if (!interactive()) {
  result <- main()
}


