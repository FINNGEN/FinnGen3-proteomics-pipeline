#!/usr/bin/env Rscript
# ==============================================================================
# path_utils.R - Batch-Aware Path Construction Utilities
# ==============================================================================
#
# Purpose:
#   Provides shared utility functions for constructing batch-aware output paths
#   across all pipeline scripts. Ensures consistent file naming and directory
#   structure for both single-batch and multi-batch modes. All paths are derived
#   from configuration files or environment variables with no hardcoded paths.
#
# Author: Reza Jabal, PhD (rjabal@broadinstitute.org)
# Date: December 2025
# ==============================================================================

# This file provides utility functions for constructing batch-aware output paths
# across all pipeline scripts. It ensures consistent file naming and directory
# structure for both single-batch and multi-batch modes.
# All paths are derived from configuration files or environment variables.

#' Get batch-aware output path
#'
#' Constructs output file paths with step prefixes and batch suffixes.
#' In single-batch mode, returns standard paths. In multi-batch mode,
#' includes batch-specific subdirectories and suffixes.
#'
#' @param step_num Character string of step number (e.g., "00", "05", "11")
#' @param filename Character string of filename (without extension)
#' @param batch_id Character string of batch identifier (e.g., "batch_01", "batch_02")
#'                 If NULL, uses default batch from config or single-batch mode
#' @param subdir Character string of output subdirectory (e.g., "qc", "outliers", "normalized")
#' @param extension Character string of file extension (default: "rds")
#' @param config Configuration object loaded from YAML
#' @param aggregate Logical indicating if this is an aggregate output (default: FALSE)
#'
#' @return Character string of full output path
#'
#' @examples
#' # Single-batch mode
#' get_output_path("05", "sex_proteins", subdir = "outliers", config = config)
#' # Returns: "output/outliers/05_sex_proteins.rds"
#'
#' # Multi-batch mode
#' get_output_path("05", "sex_proteins", batch_id = "batch_01", subdir = "outliers", config = config)
#' # Returns: "output/outliers/batch_01/05_sex_proteins_fg3_batch_01.rds"
#'
#' # Aggregate output
#' get_output_path("11", "phenotype_matrix", aggregate = TRUE, subdir = "phenotypes", config = config)
#' # Returns: "output/phenotypes/aggregate/11_aggregate_phenotype_matrix.rds"
get_output_path <- function(step_num, filename, batch_id = NULL, subdir = NULL,
                            extension = "rds", config = NULL, aggregate = FALSE) {
  # Load config if not provided
  if (is.null(config)) {
    config_file <- Sys.getenv("PIPELINE_CONFIG", "")
    if (config_file == "" || !file.exists(config_file)) {
      stop("Config file not found. Set PIPELINE_CONFIG environment variable or provide config object.")
    }
    config <- yaml::read_yaml(config_file)
  }

  # Get base output directory from config or environment
  base_dir <- config$output$base_dir
  if (is.null(base_dir)) {
    # Try environment variable as fallback
    base_dir <- Sys.getenv("PIPELINE_OUTPUT_DIR", "output")
  }

  # Check if multi-batch mode
  multi_batch_mode <- tryCatch(
    isTRUE(config$parameters$normalization$multi_batch_mode),
    error = function(e) FALSE
  )

  # Determine batch_id if not provided
  if (is.null(batch_id) && !aggregate) {
    # Try to get current batch from config
    batch_id <- tryCatch(
      config$batch$current_batch_id,
      error = function(e) NULL
    )

    # If still NULL and multi-batch mode, use default
    if (is.null(batch_id) && multi_batch_mode) {
      batch_id <- config$batch$default_batch_id %||% "batch_01"
    }
  }

  # Construct path components
  path_parts <- c(base_dir)

  # Add subdirectory
  if (!is.null(subdir)) {
    path_parts <- c(path_parts, subdir)
  }

  # Add batch subdirectory for multi-batch mode (not for aggregate)
  if (multi_batch_mode && !is.null(batch_id) && !aggregate) {
    path_parts <- c(path_parts, batch_id)
  }

  # Add aggregate subdirectory
  if (aggregate) {
    path_parts <- c(path_parts, "aggregate")
  }

  # Construct filename
  if (aggregate) {
    # Aggregate outputs: step_num_aggregate_filename.extension
    file_basename <- paste0(step_num, "_aggregate_", filename)
  } else if (multi_batch_mode && !is.null(batch_id)) {
    # Multi-batch outputs: step_num_filename_fg3_batch_id.extension
    # Convert batch_01 to fg3_batch_01
    batch_suffix <- gsub("^batch_", "fg3_batch_", batch_id)
    file_basename <- paste0(step_num, "_", filename, "_", batch_suffix)
  } else {
    # Single-batch outputs: step_num_filename.extension
    file_basename <- paste0(step_num, "_", filename)
  }

  # Add extension
  full_filename <- paste0(file_basename, ".", extension)

  # Combine path
  full_path <- do.call(file.path, as.list(c(path_parts, full_filename)))

  return(full_path)
}

#' Get batch-aware log path
#'
#' Constructs log file paths with batch-specific subdirectories.
#'
#' @param step_num Character string of step number
#' @param batch_id Character string of batch identifier (optional)
#' @param config Configuration object
#'
#' @return Character string of full log path
get_log_path <- function(step_num, batch_id = NULL, config = NULL) {
  # Load config if not provided
  if (is.null(config)) {
    config_file <- Sys.getenv("PIPELINE_CONFIG", "")
    if (config_file == "" || !file.exists(config_file)) {
      stop("Config file not found. Set PIPELINE_CONFIG environment variable or provide config object.")
    }
    config <- yaml::read_yaml(config_file)
  }

  # Get base output directory from config or environment
  base_dir <- config$output$base_dir
  if (is.null(base_dir)) {
    base_dir <- Sys.getenv("PIPELINE_OUTPUT_DIR", "output")
  }

  # Check if multi-batch mode
  multi_batch_mode <- tryCatch(
    isTRUE(config$parameters$normalization$multi_batch_mode),
    error = function(e) FALSE
  )

  # Determine batch_id if not provided
  if (is.null(batch_id)) {
    batch_id <- tryCatch(
      config$batch$current_batch_id,
      error = function(e) NULL
    )
  }

  # Construct log path
  log_dir <- file.path(base_dir, config$output$logs_dir %||% "logs")

  # Always add batch_id to log directory if provided (for batch-aware logging)
  # This ensures logs are organized by batch even in single-batch mode
  if (!is.null(batch_id) && batch_id != "") {
    log_dir <- file.path(log_dir, batch_id)
  }

  log_filename <- paste0(
    step_num, "_", gsub(
      "^[0-9]+[a-z]*_", "",  # Match step numbers like "05", "05a", "05b", "05c"
      gsub(
        "\\.R$", "",
        basename(Sys.getenv("SCRIPT_NAME", "script"))
      )
    ),
    ".log"
  )

  # If we can't determine script name, use step_num
  if (log_filename == paste0(step_num, "_.log")) {
    log_filename <- paste0(step_num, "_pipeline.log")
  }

  full_log_path <- file.path(log_dir, log_filename)
  return(full_log_path)
}

#' Get input path for a batch
#'
#' Retrieves input file path for a specific batch from config.
#'
#' @param file_type Character string of file type (e.g., "npx_matrix_file", "metadata_file")
#' @param batch_id Character string of batch identifier
#' @param config Configuration object
#'
#' @return Character string of input file path, or NULL if not found
get_batch_input_path <- function(file_type, batch_id, config) {
  # Check if multi-batch mode
  multi_batch_mode <- tryCatch(
    isTRUE(config$parameters$normalization$multi_batch_mode),
    error = function(e) FALSE
  )

  if (!multi_batch_mode) {
    # Single-batch mode: use data section
    return(tryCatch(config$data[[file_type]], error = function(e) NULL))
  }

  # Multi-batch mode: check batch-specific sections
  batch_config <- tryCatch(
    config$batches[[batch_id]],
    error = function(e) NULL
  )

  if (!is.null(batch_config)) {
    # First check inside data section
    path <- tryCatch(batch_config$data[[file_type]], error = function(e) NULL)
    if (!is.null(path)) {
      return(path)
    }
    # If not found in data section, check batch level
    path <- tryCatch(batch_config[[file_type]], error = function(e) NULL)
    if (!is.null(path)) {
      return(path)
    }
  }

  # Fallback to default data section
  return(tryCatch(config$data[[file_type]], error = function(e) NULL))
}

#' Ensure output directory exists
#'
#' Creates output directory structure if it doesn't exist.
#'
#' @param file_path Character string of full file path
ensure_output_dir <- function(file_path) {
  dir_path <- dirname(file_path)
  if (!dir.exists(dir_path)) {
    dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
  }
  return(invisible(NULL))
}

#' Get step number from script name
#'
#' Extracts step number from script filename (e.g., "05_sex_outliers.R" -> "05")
#'
#' @param script_name Character string of script name or path
#'
#' @return Character string of step number
get_step_number <- function(script_name = NULL) {
  # First, check environment variable (set by run_pipeline.R)
  step_num_env <- Sys.getenv("PIPELINE_STEP_NUM", "")
  if (step_num_env != "") {
    return(step_num_env)
  }

  if (is.null(script_name)) {
    # Try to get from environment variable
    script_name <- Sys.getenv("SCRIPT_NAME", "")

    # If still empty, try to get from command line (when executed directly)
    if (script_name == "") {
      args <- commandArgs(trailingOnly = FALSE)
      file_arg <- grep("^--file=", args, value = TRUE)
      if (length(file_arg) > 0) {
        script_name <- sub("^--file=", "", file_arg)
      }
    }

    # If still empty, try to infer from current script file (when sourced)
    if (script_name == "") {
      # Try to get from sys.frame() when script is sourced
      frames <- sys.frames()
      for (frame in frames) {
        if (exists("ofile", frame)) {
          script_name <- frame$ofile
          break
        }
      }
      # Alternative: check if we can get script path from source() call
      if (script_name == "") {
        # Look for script path in parent frames
        for (i in seq_len(sys.nframe())) {
          frame <- sys.frame(i)
          if (exists("script_path", frame, inherits = FALSE)) {
            script_name <- frame$script_path
            break
          }
        }
      }
    }
  }

  # Extract step number (first two digits) from script name
  if (script_name != "") {
    step_match <- regmatches(basename(script_name), regexpr("^[0-9]{2}", basename(script_name)))
    if (length(step_match) > 0) {
      return(step_match[1])
    }
  }

  # Default fallback - but this should rarely happen if environment is set correctly
  warning("Could not determine step number, defaulting to '00'. Set PIPELINE_STEP_NUM environment variable.")
  return("00")
}

# Helper function for NULL coalescing
`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

#' Get script directory safely (handles both direct execution and sourcing)
#'
#' @return Character string of script directory path
get_script_dir <- function() {
  # Try to get from environment (set by run_pipeline.R)
  env_script <- Sys.getenv("SCRIPT_NAME", "")
  if (env_script != "" && file.exists(env_script)) {
    return(dirname(normalizePath(env_script)))
  }

  # Try to get from command line (when executed directly)
  tryCatch(
    {
      cmd_args <- commandArgs(trailingOnly = FALSE)
      file_arg <- grep("^--file=", cmd_args, value = TRUE)
      if (length(file_arg) > 0) {
        script_path <- sub("^--file=", "", file_arg)
        return(dirname(normalizePath(script_path)))
      }
    },
    error = function(e) NULL
  )

  # Fallback: use current working directory
  return(getwd())
}

#' Add FINNGENID column to data.table/data.frame with SampleID
#'
#' This utility function adds a FINNGENID column to any data structure
#' containing SampleID. It loads the sample mapping from step 00 and merges
#' FINNGENID information. Handles cases where SampleID might already be FINNGENID.
#'
#' @param dt data.table or data.frame with SampleID column
#' @param batch_id Character string of batch identifier (e.g., "batch_01", "batch_02")
#' @param config Configuration object loaded from YAML (optional)
#' @param sample_id_col Character string of column name containing sample IDs (default: "SampleID")
#' @param preserve_original Logical indicating if original SampleID should be preserved (default: TRUE)
#'                           If FALSE and SampleID is already FINNGENID, it will be renamed
#'
#' @return data.table with FINNGENID column added (as second column after SampleID)
#'
#' @examples
#' # Add FINNGENID to outlier summary
#' dt_with_fgid <- add_finngenid_column(outlier_summary, batch_id = "batch_02", config = config)
add_finngenid_column <- function(dt, batch_id = NULL, config = NULL,
                                  sample_id_col = "SampleID", preserve_original = TRUE) {
  # Load config if not provided
  if (is.null(config)) {
    config_file <- Sys.getenv("PIPELINE_CONFIG", "")
    if (config_file == "" || !file.exists(config_file)) {
      warning("Config file not found. Cannot add FINNGENID column.")
      return(dt)
    }
    config <- yaml::read_yaml(config_file)
  }

  # Get batch_id if not provided
  if (is.null(batch_id)) {
    batch_id <- tryCatch(
      config$batch$current_batch_id,
      error = function(e) NULL
    )
    if (is.null(batch_id)) {
      batch_id <- config$batch$default_batch_id %||% "batch_01"
    }
  }

  # Convert to data.table if needed
  if (!inherits(dt, "data.table")) {
    dt <- as.data.table(dt)
  } else {
    dt <- copy(dt)  # Avoid modifying original
  }

  # Check if SampleID column exists (try both SampleID and SAMPLE_ID)
  if (!sample_id_col %in% names(dt)) {
    # Try alternative column name
    if (sample_id_col == "SampleID" && "SAMPLE_ID" %in% names(dt)) {
      sample_id_col <- "SAMPLE_ID"
    } else if (sample_id_col == "SAMPLE_ID" && "SampleID" %in% names(dt)) {
      sample_id_col <- "SampleID"
    } else {
      warning("Column '", sample_id_col, "' not found in data. Cannot add FINNGENID.")
      return(dt)
    }
  }

  # Load sample mapping from step 00
  sample_mapping_path <- get_output_path("00", "sample_mapping", batch_id, "qc", "tsv", config = config)

  if (!file.exists(sample_mapping_path)) {
    # Try RDS format
    sample_mapping_path <- get_output_path("00", "sample_mapping", batch_id, "qc", "rds", config = config)
    if (!file.exists(sample_mapping_path)) {
      warning("Sample mapping file not found. Cannot add FINNGENID column.")
      return(dt)
    }
    sample_mapping <- readRDS(sample_mapping_path)
  } else {
    sample_mapping <- fread(sample_mapping_path)
  }

  # Ensure sample_mapping has required columns
  if (!all(c("SampleID", "FINNGENID") %in% names(sample_mapping))) {
    warning("Sample mapping missing required columns (SampleID, FINNGENID). Cannot add FINNGENID.")
    return(dt)
  }

  # Check if SampleID in dt is already FINNGENID (starts with "FG")
  sample_ids <- dt[[sample_id_col]]
  is_finngenid <- grepl("^FG[0-9A-Z]+$", sample_ids)

  if (all(is_finngenid, na.rm = TRUE) && sum(!is.na(is_finngenid)) > 0) {
    # SampleID column already contains FINNGENIDs - need to map back to original SampleID
    # Create reverse mapping: FINNGENID -> SampleID (original)
    # Handle duplicates by taking first match
    reverse_map <- sample_mapping[!is.na(FINNGENID) & grepl("^FG", FINNGENID),
                                  .(Original_SampleID = SampleID[1]),
                                  by = FINNGENID]

    # Merge to get original SampleID
    dt <- merge(dt, reverse_map,
                by.x = sample_id_col, by.y = "FINNGENID",
                all.x = TRUE)

    # Rename columns appropriately
    if (preserve_original) {
      # Keep both: original SampleID (from mapping) and FINNGENID (current SampleID column)
      # Rename current SampleID column to FINNGENID
      setnames(dt, sample_id_col, "FINNGENID")
      # Rename Original_SampleID to SampleID
      setnames(dt, "Original_SampleID", sample_id_col)

      # Reorder columns: SampleID, FINNGENID, ...
      col_order <- c(sample_id_col, "FINNGENID", setdiff(names(dt), c(sample_id_col, "FINNGENID")))
      setcolorder(dt, col_order)
    } else {
      # Just rename: SampleID becomes FINNGENID, original SampleID is lost
      setnames(dt, sample_id_col, "FINNGENID")
      setnames(dt, "Original_SampleID", sample_id_col)

      # Reorder columns: SampleID, FINNGENID, ...
      col_order <- c(sample_id_col, "FINNGENID", setdiff(names(dt), c(sample_id_col, "FINNGENID")))
      setcolorder(dt, col_order)
    }
  } else {
    # SampleID is original Olink ID - add FINNGENID from mapping
    dt <- merge(dt,
                sample_mapping[, .(SampleID, FINNGENID)],
                by.x = sample_id_col, by.y = "SampleID",
                all.x = TRUE)

    # Reorder columns: SampleID, FINNGENID, ...
    col_order <- c(sample_id_col, "FINNGENID", setdiff(names(dt), c(sample_id_col, "FINNGENID")))
    setcolorder(dt, col_order)
  }

  return(dt)
}




