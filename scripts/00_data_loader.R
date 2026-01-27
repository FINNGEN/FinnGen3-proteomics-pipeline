#!/usr/bin/env Rscript
# ==============================================================================
# 00_data_loader.R - Load NPX Matrix and Perform Initial QC
# ==============================================================================
#
# Purpose:
#   Loads NPX matrix from multiple formats (Parquet, RDS, or TSV - auto-detected),
#   performs initial QC by removing samples with >10% missing data, and prepares
#   analysis-ready matrices. Creates sample mapping with FINNGENIDs and long-format
#   sample data for downstream QC steps. Supports both single-batch and multi-batch modes.
#
# Author: Reza Jabal, PhD (rjabal@broadinstitute.org)
# Date: December 2025 (Updated: January 2026 - Added initial QC filtering)
# ==============================================================================

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

# Function to identify control probes (24 control aptamers: 8 incubation, 8 extension, 8 amplification)
# These are technical controls and should be removed early in the pipeline for robustness
identify_control_probes <- function(protein_names) {
  control_patterns <- c(
    "Incubation control",
    "Extension control",
    "Amplification control"
  )

  control_probes <- character(0)
  for (pattern in control_patterns) {
    matches <- grep(pattern, protein_names, ignore.case = TRUE, value = TRUE)
    control_probes <- c(control_probes, matches)
  }

  return(control_probes)
}

# Function to detect file format from extension
detect_file_format <- function(file_path) {
  ext <- tolower(tools::file_ext(file_path))
  if (ext %in% c("parquet")) {
    return("parquet")
  } else if (ext %in% c("rds")) {
    return("rds")
  } else if (ext %in% c("tsv", "txt", "csv")) {
    return("tsv")
  } else {
    # Try to infer from content
    if (grepl("\\.parquet$", file_path, ignore.case = TRUE)) {
      return("parquet")
    } else if (grepl("\\.rds$", file_path, ignore.case = TRUE)) {
      return("rds")
    } else {
      return("tsv")  # Default to TSV
    }
  }
}

# Function to load NPX matrix from various formats (parquet, RDS, TSV)
load_npx_matrix <- function(matrix_file) {
  log_info("Loading NPX matrix from: {matrix_file}")

  if (!file.exists(matrix_file)) {
    stop("NPX matrix file not found: ", matrix_file)
  }

  file_format <- detect_file_format(matrix_file)
  log_info("Detected file format: {file_format}")

  # Load based on format
  if (file_format == "parquet") {
    dt <- read_parquet(matrix_file)
    setDT(dt)
  } else if (file_format == "rds") {
    data <- readRDS(matrix_file)
    # Handle both matrix and data.table/data.frame
    if (is.matrix(data)) {
      dt <- as.data.table(data)
      dt[, SampleID := rownames(data)]
      setcolorder(dt, c("SampleID", setdiff(names(dt), "SampleID")))
    } else if (is.data.frame(data) || is.data.table(data)) {
      dt <- as.data.table(data)
    } else {
      stop("RDS file does not contain a matrix, data.frame, or data.table")
    }
  } else {  # TSV or CSV
    # Try to detect separator
    first_line <- readLines(matrix_file, n = 1)
    if (grepl("\t", first_line)) {
      sep <- "\t"
    } else if (grepl(",", first_line)) {
      sep <- ","
    } else {
      sep <- "\t"  # Default to tab
    }
    dt <- fread(matrix_file, sep = sep)
  }

  # Check if SampleID column exists
  if (!"SampleID" %in% names(dt)) {
    # Try common alternatives
    id_cols <- c("SAMPLE_ID", "sample_id", "Sample_ID", "ID", "id")
    found_id_col <- NULL
    for (col in id_cols) {
      if (col %in% names(dt)) {
        found_id_col <- col
        break
      }
    }

    if (!is.null(found_id_col)) {
      setnames(dt, found_id_col, "SampleID")
      log_info("Renamed column '{found_id_col}' to 'SampleID'")
    } else {
      # Try first column as SampleID
      if (ncol(dt) > 0) {
        setnames(dt, names(dt)[1], "SampleID")
        log_info("Using first column as SampleID")
      } else {
        stop("Could not identify SampleID column in matrix file")
      }
    }
  }

  # Detect if data is in long format (one row per sample-protein pair)
  # Long format indicators: Has SampleID, Assay/OlinkID, and NPX/PCNormalizedNPX columns
  has_assay_col <- "Assay" %in% names(dt) || "OlinkID" %in% names(dt)
  has_npx_col <- "NPX" %in% names(dt) || "PCNormalizedNPX" %in% names(dt)
  is_long_format <- has_assay_col && has_npx_col

  if (is_long_format) {
    log_info("Detected LONG format NPX file - converting to wide format matrix")
    
    # Identify protein identifier column (Assay or OlinkID)
    protein_id_col <- if ("Assay" %in% names(dt)) "Assay" else "OlinkID"
    log_info("Using '{protein_id_col}' as protein identifier")
    
    # Identify NPX value column (prefer NPX over PCNormalizedNPX for raw data)
    npx_value_col <- if ("NPX" %in% names(dt)) "NPX" else "PCNormalizedNPX"
    log_info("Using '{npx_value_col}' as NPX value column")
    
    # Check for duplicate sample-protein pairs
    key_cols <- c("SampleID", protein_id_col)
    n_duplicates <- sum(duplicated(dt[, ..key_cols]))
    if (n_duplicates > 0) {
      log_warn("Found {n_duplicates} duplicate sample-protein pairs - using mean value")
      # Aggregate duplicates by taking mean
      dt <- dt[, .(NPX = mean(get(npx_value_col), na.rm = TRUE)), 
                  by = c("SampleID", protein_id_col)]
      setnames(dt, "NPX", npx_value_col)
    }
    
    # Convert long to wide using dcast (data.table's pivot function)
    log_info("Converting long format to wide format matrix...")
    # Use formula interface for dcast - need to construct formula dynamically
    formula_str <- paste("SampleID ~", protein_id_col)
    dt_wide <- dcast(
      dt,
      as.formula(formula_str),
      value.var = npx_value_col,
      fun.aggregate = mean,  # Handle any remaining duplicates
      na.rm = TRUE
    )
    
    # Extract SampleID and protein columns
    sample_ids <- dt_wide$SampleID
    protein_cols <- setdiff(names(dt_wide), "SampleID")
    
    # Convert to matrix
    npx_matrix <- as.matrix(dt_wide[, ..protein_cols])
    rownames(npx_matrix) <- sample_ids
    
    log_info("Converted long format to wide: {nrow(npx_matrix)} samples x {ncol(npx_matrix)} proteins")
  } else {
    # WIDE format: Extract SampleID and convert to matrix
    sample_ids <- dt$SampleID
    protein_cols <- setdiff(names(dt), "SampleID")
    
    # Convert to matrix
    npx_matrix <- as.matrix(dt[, ..protein_cols])
    rownames(npx_matrix) <- sample_ids
    
    log_info("Loaded wide format NPX matrix: {nrow(npx_matrix)} samples x {ncol(npx_matrix)} proteins")
  }

  log_info("Final matrix dimensions: {nrow(npx_matrix)} samples x {ncol(npx_matrix)} proteins")
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

  # Try to detect file format (TSV or semicolon-separated)
  first_line <- readLines(bridging_file, n = 1)
  if (grepl("\t", first_line)) {
    # TSV format (new unified format)
    bridging <- fread(bridging_file, sep = "\t")
    log_info("Detected TSV format (unified bridging metadata)")
    
    # Validate required columns
    required_cols <- c("FINNGENID")
    if (!all(required_cols %in% names(bridging))) {
      log_warn("Unified bridging file missing required columns. Expected: {paste(required_cols, collapse=', ')}")
      log_warn("Found columns: {paste(names(bridging), collapse=', ')}")
      return(data.table())
    }
    
    # Filter to only include FINNGENIDs with "FG" prefix
    n_before <- nrow(bridging)
    bridging <- bridging[grepl("^FG", FINNGENID)]
    n_after <- nrow(bridging)
    
    if (n_before > n_after) {
      log_warn("Filtered out {n_before - n_after} entries that don't have FINNGENID starting with 'FG' prefix")
    }
    
    log_info("Loaded {nrow(bridging)} bridging samples from unified metadata (all with 'FG' prefix FINNGENIDs)")
  } else {
    # Legacy semicolon-separated format (for backward compatibility)
    bridging <- fread(bridging_file, sep = ";")
    log_info("Detected semicolon-separated format (legacy)")
    log_info("Bridging samples: {nrow(bridging)}")
  }

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
create_sample_mapping <- function(sample_ids, metadata, bridging_samples = NULL, bridging_finngenid_map = NULL, batch_id = NULL) {
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

  # Identify bridging samples using unified bridging metadata (new method)
  bridging_finngenids <- character(0)
  if (!is.null(bridging_samples) && nrow(bridging_samples) > 0) {
    # Check if this is the new unified format (has FINNGENID column)
    if ("FINNGENID" %in% names(bridging_samples)) {
      log_info("Using unified bridging metadata format (FINNGENID-based)")
      
      # Determine which batch flags to check based on batch_id
      # Handle both logical (TRUE/FALSE) and character ("TRUE"/"FALSE") boolean values
      if (!is.null(batch_id)) {
        if (batch_id == "batch_01") {
          # For batch_01, bridging samples are those in both batch_01 and batch_02
          # CRITICAL: Exclude EA5 samples (in_Batch_00=TRUE) - these are for Batch_00↔Batch_02 bridging, not Batch_01↔Batch_02
          # We want ONLY the 24 direct Batch_01↔Batch_02 samples (not the 6 EA5 samples in all three batches)
          if ("in_Batch_01" %in% names(bridging_samples) && "in_Batch_02" %in% names(bridging_samples)) {
            # Robust boolean filtering: handles both logical and character representations
            batch_01_true <- bridging_samples$in_Batch_01 == TRUE | bridging_samples$in_Batch_01 == "TRUE" | bridging_samples$in_Batch_01 == "T"
            batch_02_true <- bridging_samples$in_Batch_02 == TRUE | bridging_samples$in_Batch_02 == "TRUE" | bridging_samples$in_Batch_02 == "T"
            # Exclude EA5 samples (in_Batch_00=TRUE) - these are primarily for Batch_00↔Batch_02 bridging
            if ("in_Batch_00" %in% names(bridging_samples)) {
              batch_00_false <- bridging_samples$in_Batch_00 == FALSE | bridging_samples$in_Batch_00 == "FALSE" | bridging_samples$in_Batch_00 == "F"
              bridging_finngenids <- bridging_samples[batch_01_true & batch_02_true & batch_00_false]$FINNGENID
              log_info("Identified {length(bridging_finngenids)} bridging samples for batch_01 (direct Batch_01↔Batch_02, excluding EA5 samples)")
            } else {
              # Fallback if in_Batch_00 column doesn't exist
              bridging_finngenids <- bridging_samples[batch_01_true & batch_02_true]$FINNGENID
              log_info("Identified {length(bridging_finngenids)} bridging samples for batch_01 (present in both batch_01 and batch_02)")
            }
          }
        } else if (batch_id == "batch_02") {
          # For batch_02, bridging samples are those in batch_02 and either batch_00 or batch_01
          # This includes: 50 EA5 samples (Batch_00 ↔ Batch_02) + 24 direct Batch_01↔Batch_02 samples = 74 unique FINNGENIDs
          # Note: If duplicate FINNGENIDs exist in the sample mapping, the count of SampleIDs flagged as bridging
          # may exceed 74 (e.g., 75 if one bridging FINNGENID appears twice in the mapping)
          if ("in_Batch_02" %in% names(bridging_samples)) {
            # Robust boolean filtering: handles both logical and character representations
            batch_02_true <- bridging_samples$in_Batch_02 == TRUE | bridging_samples$in_Batch_02 == "TRUE" | bridging_samples$in_Batch_02 == "T"
            bridging_finngenids <- bridging_samples[batch_02_true]$FINNGENID
            log_info("Identified {length(bridging_finngenids)} unique bridging FINNGENIDs for batch_02")
            log_info("  Note: Actual number of SampleIDs flagged as bridging may be higher if duplicate FINNGENIDs exist in sample mapping")
          }
        }
      } else {
        # If batch_id not provided, use all bridging samples
        bridging_finngenids <- unique(bridging_samples$FINNGENID)
        log_info("Identified {length(bridging_finngenids)} total bridging samples (batch_id not specified)")
      }
    } else {
      # Legacy format - use PSEUDO_ID matching (for backward compatibility)
      log_info("Using legacy bridging samples format (PSEUDO_ID-based)")
      # This would require the old bridging_finngenid_map logic
    }
  }

  # Supplement with FINNGENIDs from EA5→FG3 mapping for bridging samples (legacy method, for backward compatibility)
  if (!is.null(bridging_finngenid_map)) {
    log_info("Supplementing sample mapping with FINNGENIDs from EA5→FG3 mapping (legacy method)")

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
    
    # Also add EA5 FINNGENIDs to bridging list if they match
    if (nrow(bridging_finngenid_map) > 0) {
      ea5_bridging_fgids <- unique(bridging_finngenid_map$FINNGENID)
      bridging_finngenids <- unique(c(bridging_finngenids, ea5_bridging_fgids))
    }
  }

  # Flag sample types (AG samples already removed, so only FinnGen and Bridging)
  # Bridging samples are identified by FINNGENID matching, not by SampleID pattern
  mapping[, sample_type := case_when(
    grepl("CONTROL", SampleID) ~ "Control",
    !is.na(FINNGENID) & FINNGENID %in% bridging_finngenids ~ "Bridging",
    grepl("EA5_OLI", SampleID) ~ "Bridging",  # Legacy: EA5_OLI pattern
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

# Function to create long-format samples_data_raw from wide matrix
# This reconstructs the long-format structure needed for technical outlier detection
# Format: Multiple rows per sample (one per protein), matching original pipeline structure
create_long_format_samples_data <- function(npx_matrix, metadata) {
  log_info("Reconstructing long-format samples_data_raw from wide matrix")
  log_info("  This is required for technical outlier detection (sd_npx calculation)")

  # Get sample IDs and protein names
  sample_ids <- rownames(npx_matrix)
  protein_names <- colnames(npx_matrix)
  n_samples <- length(sample_ids)
  n_proteins <- length(protein_names)

  log_info("  Matrix dimensions: {n_samples} samples × {n_proteins} proteins")
  log_info("  Expected long-format rows: {n_samples * n_proteins}")

  # Create PlateID mapping from metadata
  plate_map <- NULL
  if (!is.null(metadata) && "PlateID" %in% names(metadata) && "SAMPLE_ID" %in% names(metadata)) {
    plate_map <- metadata[, .(SAMPLE_ID, PlateID)]
    log_info("  PlateID mapping available for {nrow(plate_map)} samples")
  } else {
    log_warn("  PlateID not found in metadata - will be set to NA")
  }

  # Melt matrix to long format
  # Convert matrix to data.table with SampleID
  dt_matrix <- as.data.table(npx_matrix, keep.rownames = "SampleID")

  # Melt to long format: SampleID, Assay (protein), NPX
  samples_data_long <- melt(
    dt_matrix,
    id.vars = "SampleID",
    variable.name = "Assay",
    value.name = "NPX",
    variable.factor = FALSE
  )

  # Add PlateID from metadata
  if (!is.null(plate_map)) {
    samples_data_long <- merge(
      samples_data_long,
      plate_map,
      by.x = "SampleID",
      by.y = "SAMPLE_ID",
      all.x = TRUE
    )
  } else {
    samples_data_long[, PlateID := NA_character_]
  }

  # Add QC flags (default to PASS since we're working with pre-filtered data)
  samples_data_long[, SampleQC := "PASS"]
  samples_data_long[, AssayQC := "PASS"]

  # Add SampleType (all should be SAMPLE since controls are already filtered)
  samples_data_long[, SampleType := "SAMPLE"]

  # Set column order to match original pipeline structure
  setcolorder(samples_data_long, c("SampleID", "Assay", "NPX", "PlateID", "SampleQC", "AssayQC", "SampleType"))

  # Validate structure
  n_rows_expected <- n_samples * n_proteins
  n_rows_actual <- nrow(samples_data_long)

  if (n_rows_actual != n_rows_expected) {
    log_warn("  Row count mismatch: expected {n_rows_expected}, got {n_rows_actual}")
  } else {
    log_info("  Long-format structure validated: {n_rows_actual} rows")
  }

  # Validate that each sample has the correct number of proteins
  sample_counts <- samples_data_long[, .N, by = SampleID]
  expected_per_sample <- n_proteins
  mismatched_samples <- sample_counts[N != expected_per_sample]

  if (nrow(mismatched_samples) > 0) {
    log_warn("  {nrow(mismatched_samples)} samples have incorrect protein counts")
  } else {
    log_info("  All samples have {expected_per_sample} protein measurements")
  }

  # Verify that sd_npx can be calculated (critical for technical outlier detection)
  test_sample <- sample_ids[1]
  test_data <- samples_data_long[SampleID == test_sample]
  test_sd <- sd(test_data$NPX, na.rm = TRUE)

  if (is.na(test_sd) || length(test_data$NPX) < 2) {
    log_error("  CRITICAL: Cannot calculate sd_npx - long-format structure is incorrect!")
    stop("Long-format samples_data_raw structure validation failed")
  } else {
    log_info("  Validation passed: sd_npx calculation works (test sample SD: {sprintf('%.3f', test_sd)})")
  }

  return(samples_data_long)
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
  npx_matrix <- load_npx_matrix(npx_matrix_file)
  metadata <- load_metadata(metadata_file)

  # CRITICAL: Remove control aptamers (24 technical controls) early in the pipeline
  # This ensures both batches have the same number of biological proteins (5416)
  # and prevents dimension mismatches in cross-batch normalization
  log_info("Removing control aptamers (technical controls) from NPX matrix...")
  protein_names <- colnames(npx_matrix)
  control_probes <- identify_control_probes(protein_names)
  
  if (length(control_probes) > 0) {
    log_info("  Control aptamers identified: {length(control_probes)}")
    incubation <- sum(grepl("Incubation control", control_probes, ignore.case = TRUE))
    extension <- sum(grepl("Extension control", control_probes, ignore.case = TRUE))
    amplification <- sum(grepl("Amplification control", control_probes, ignore.case = TRUE))
    log_info("    - Incubation controls: {incubation}")
    log_info("    - Extension controls: {extension}")
    log_info("    - Amplification controls: {amplification}")
    
    # Remove control probes
    biological_proteins <- setdiff(protein_names, control_probes)
    npx_matrix <- npx_matrix[, biological_proteins, drop = FALSE]
    log_info("  Removed {length(control_probes)} control aptamers")
    log_info("  Biological proteins retained: {ncol(npx_matrix)}")
    log_info("  Matrix dimensions after control removal: {nrow(npx_matrix)} samples × {ncol(npx_matrix)} proteins")
  } else {
    log_info("  No control aptamers found in matrix (may have been removed previously)")
    log_info("  Biological proteins: {ncol(npx_matrix)}")
  }

  # Load optional bridging sample information
  bridging_samples <- load_bridging_samples(bridging_samples_file)
  bridging_finngenid_map <- NULL
  # Only load legacy EA5 mapping if bridging_samples_file doesn't contain FINNGENIDs (backward compatibility)
  if (!is.null(bridging_samples_finngenid_map_file) && 
      (is.null(bridging_samples) || nrow(bridging_samples) == 0 || !"FINNGENID" %in% names(bridging_samples))) {
    bridging_finngenid_map <- load_bridging_finngenid_map(bridging_samples_finngenid_map_file)
  }

  # Perform initial QC: Remove samples with >10% missing data
  # This matches the original pipeline's Step 02 (02_initial_qc.R) functionality
  log_info("Performing initial QC: filtering samples with >10% missing data")
  initial_qc_threshold <- config$parameters$qc$max_missing_per_sample %||% 0.10

  # Calculate missing rate per sample
  sample_missing_rates <- rowSums(is.na(npx_matrix)) / ncol(npx_matrix)

  # Identify samples failing initial QC
  initial_qc_failed <- names(sample_missing_rates)[sample_missing_rates > initial_qc_threshold]

  # Create QC summary for tracking (needed for comprehensive report)
  qc_summary <- data.table(
    SampleID = names(sample_missing_rates),
    missing_rate = sample_missing_rates,
    n_missing = rowSums(is.na(npx_matrix)),
    n_measured = rowSums(!is.na(npx_matrix)),
    QC_initial_qc = as.integer(sample_missing_rates > initial_qc_threshold),
    qc_pass = sample_missing_rates <= initial_qc_threshold
  )

  if(length(initial_qc_failed) > 0) {
    log_info("Removing {length(initial_qc_failed)} samples with >{initial_qc_threshold*100}% missing data")
    log_info("  Failed samples: {paste(initial_qc_failed, collapse=', ')}")

    # Remove failing samples from matrix
    npx_matrix <- npx_matrix[!rownames(npx_matrix) %in% initial_qc_failed, ]

    log_info("Matrix after initial QC: {nrow(npx_matrix)} samples (removed {length(initial_qc_failed)} samples)")
  } else {
    log_info("No samples failed initial QC (all samples have ≤{initial_qc_threshold*100}% missing data)")
  }

  # Save QC summary for comprehensive report integration
  qc_summary_path <- get_output_path(step_num, "qc_summary", batch_id, "qc", "tsv", config = config)
  ensure_output_dir(qc_summary_path)
  fwrite(qc_summary, qc_summary_path, sep = "\t")
  log_info("Saved initial QC summary: {qc_summary_path}")
  log_info("  - Total samples evaluated: {nrow(qc_summary)}")
  log_info("  - Samples passing QC: {sum(qc_summary$qc_pass)}")
  log_info("  - Samples failing QC: {sum(!qc_summary$qc_pass)}")

  # Create sample mapping
  log_info("Creating sample mapping...")
  sample_ids <- rownames(npx_matrix)
  sample_mapping <- create_sample_mapping(sample_ids, metadata, bridging_samples, bridging_finngenid_map, batch_id)

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

  # Create long-format samples_data_raw from wide matrix for downstream steps
  # This is CRITICAL for technical outlier detection which needs to calculate sd_npx
  # across all proteins per sample (requires multiple rows per sample)
  log_info("Creating long-format samples_data_raw from matrix...")
  samples_data_raw <- create_long_format_samples_data(analysis_matrix, metadata)

  samples_data_raw_path <- get_output_path(step_num, "samples_data_raw", batch_id, "qc", config = config)
  ensure_output_dir(samples_data_raw_path)
  saveRDS(samples_data_raw, samples_data_raw_path)
  log_info("Saved long-format samples_data_raw: {samples_data_raw_path}")
  log_info("  - Total rows: {nrow(samples_data_raw)} (long format: sample × protein pairs)")
  log_info("  - Unique samples: {length(unique(samples_data_raw$SampleID))}")
  log_info("  - Unique proteins: {length(unique(samples_data_raw$Assay))}")
  log_info("  - Average proteins per sample: {round(nrow(samples_data_raw) / length(unique(samples_data_raw$SampleID)), 1)}")

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

  # Validate long-format samples_data_raw structure
  log_info("Validating long-format samples_data_raw structure...")
  sample_stats_test <- samples_data_raw[, .(
    n_proteins = .N,
    mean_npx = mean(NPX, na.rm = TRUE),
    sd_npx = sd(NPX, na.rm = TRUE)
  ), by = SampleID]

  n_valid_sd <- sum(!is.na(sample_stats_test$sd_npx))
  n_samples_total <- nrow(sample_stats_test)

  if (n_valid_sd == n_samples_total) {
    log_info("  ✓ All {n_samples_total} samples have valid sd_npx (technical outlier detection will work)")
  } else {
    log_error("  ✗ Only {n_valid_sd}/{n_samples_total} samples have valid sd_npx")
    stop("Long-format samples_data_raw validation failed - technical outlier detection will not work")
  }

  # Print summary
  cat("\n=== DATA LOADER SUMMARY ===\n")
  cat("Input file format: ", detect_file_format(npx_matrix_file), "\n", sep = "")
  n_samples_loaded <- nrow(qc_summary)
  n_samples_after_qc <- nrow(npx_matrix)
  n_initial_qc_failed <- sum(qc_summary$QC_initial_qc == 1)
  cat("NPX matrix loaded: ", n_samples_loaded, " samples x ", ncol(npx_matrix), " proteins\n", sep = "")
  if(n_initial_qc_failed > 0) {
    cat("Initial QC: Removed ", n_initial_qc_failed, " samples with >10% missing data\n", sep = "")
    cat("  Matrix after initial QC: ", n_samples_after_qc, " samples\n", sep = "")
  } else {
    cat("Initial QC: All samples passed (≤10% missing data)\n", sep = "")
  }
  cat("Analysis-ready samples: ", nrow(analysis_samples), "\n", sep = "")
  cat("Long-format samples_data_raw: ", nrow(samples_data_raw), " rows (",
      length(unique(samples_data_raw$SampleID)), " samples × ",
      length(unique(samples_data_raw$Assay)), " proteins)\n", sep = "")
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


