#!/usr/bin/env Rscript
# ==============================================================================
# run_pipeline.R - Master Pipeline Runner
# ==============================================================================
#
# Purpose:
#   Master pipeline runner for FinnGen 3 Olink analysis. Supports both single-batch
#   and multi-batch processing modes. Orchestrates execution of all pipeline steps
#   (00-11) with proper error handling, logging, and step dependencies. Provides
#   command-line interface for flexible execution control.
#
# Author: Reza Jabal, PhD (rjabal@broadinstitute.org)
# Date: December 2025
# ==============================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(yaml)
  library(logger)
})

# Source path utilities
script_dir <- tryCatch(
  {
    args <- commandArgs(trailingOnly = FALSE)
    file_arg <- grep("^--file=", args, value = TRUE)
    if (length(file_arg) > 0) {
      script_path <- sub("^--file=", "", file_arg)
      dirname(normalizePath(script_path))
    } else {
      getwd()
    }
  },
  error = function(e) getwd()
)
source(file.path(script_dir, "path_utils.R"))

# Command line options
option_list <- list(
  make_option(c("-c", "--config"),
    type = "character",
    default = NULL,
    help = "Configuration file path (required)"
  ),
  make_option(c("-s", "--step"),
    type = "character",
    default = "all",
    help = "Which step to run (all, qc, outlier, normalize, phenotype, pqtl, report) [default %default]"
  ),
  make_option(c("-f", "--from"),
    type = "character",
    default = NULL,
    help = "Start from specific step"
  ),
  make_option(c("-t", "--to"),
    type = "character",
    default = NULL,
    help = "Run until specific step"
  ),
  make_option(c("-d", "--dry_run"),
    action = "store_true",
    default = FALSE,
    help = "Dry run - show what would be executed without running"
  ),
  make_option(c("-v", "--verbose"),
    action = "store_true",
    default = FALSE,
    help = "Verbose output"
  ),
  make_option(c("-b", "--batch"),
    type = "character",
    default = NULL,
    help = "Process specific batch only (e.g., batch_01, batch_02). Overrides multi-batch mode."
  )
)

# Parse arguments
parser <- OptionParser(
  option_list = option_list,
  description = "FinnGen 3 Olink Analysis Pipeline (Single/Multi-Batch) - Refactored Version"
)
args <- parse_args(parser)

# Validate config file
if (is.null(args$config)) {
  # Try environment variable
  args$config <- Sys.getenv("PIPELINE_CONFIG", "")
  if (args$config == "" || !file.exists(args$config)) {
    stop("Configuration file required. Use --config PATH or set PIPELINE_CONFIG environment variable.")
  }
}

# Set up logging
if (args$verbose) {
  log_threshold(DEBUG)
} else {
  log_threshold(INFO)
}

# Load config early to determine mode
if (!file.exists(args$config)) {
  stop("Configuration file not found: ", args$config)
}
config <- read_yaml(args$config)

# Get base output directory from config or environment
base_dir <- config$output$base_dir
if (is.null(base_dir)) {
  base_dir <- Sys.getenv("PIPELINE_OUTPUT_DIR", "output")
}
log_dir <- file.path(base_dir, config$output$logs_dir %||% "logs")
try(dir.create(log_dir, recursive = TRUE, showWarnings = FALSE), silent = TRUE)

# Set up logging: file appender (always) + console appender (if verbose)
log_file_path <- file.path(log_dir, "pipeline_master.log")
if (args$verbose) {
  # In verbose mode: log to both file and console
  # Use appender_tee to write to both file and console
  log_appender(appender_tee(log_file_path, append = TRUE))
} else {
  # In normal mode: log only to file
  log_appender(appender_file(log_file_path, append = TRUE))
}

log_info("Starting FinnGen 3 Olink Analysis Pipeline (Refactored Version)")
log_info("Configuration loaded from: {args$config}")

# Define pipeline steps with new numbering
pipeline_steps <- list(
  "00_data_loader" = "Load pre-filtered NPX matrix and prepare analysis-ready data",
  "01_pca_outliers" = "PCA-based outlier detection",
  "02_technical_outliers" = "Technical outlier detection",
  "03_zscore_outliers" = "Z-score outlier detection",
  "04_sex_outliers" = "Sex mismatch detection",
  "05a_pqtl_training" = "pQTL training - LASSO feature selection (optional)",
  "05b_pqtl_outliers" = "pQTL-based outlier detection",
  "05c_provenance_test" = "Provenance component validation test (optional)",
  "05d_qc_comprehensive_report" = "Comprehensive QC report generation",
  "06_normalize_data" = "Data normalization",
  "07_bridge_normalization" = "Bridge sample normalization",
  "08_covariate_adjustment" = "Covariate adjustment",
  "09_prepare_phenotypes" = "Phenotype preparation",
  "10_kinship_filtering" = "Kinship filtering",
  "11_rank_normalize" = "Rank normalization"
)

# Steps that require cross-batch operations (run after all batches complete steps 00-05d)
# NOTE: 06_normalize_data is NOT a true cross-batch step - it runs per-batch but can use cross-batch data
cross_batch_steps <- c("07_bridge_normalization")

# Steps that require normalized data (run after normalization steps 06-07)
# NOTE: 06_normalize_data runs per-batch in Phase 3.5, before step 08
post_normalization_steps <- c("06_normalize_data", "08_covariate_adjustment")

# Steps that handle aggregation (run after all batches complete steps 00-08)
aggregation_steps <- c("09_prepare_phenotypes", "10_kinship_filtering", "11_rank_normalize")

# Function to check if script exists
check_script <- function(script_name) {
  script_path <- file.path(script_dir, paste0(script_name, ".R"))
  if (!file.exists(script_path)) {
    log_warn("Script not found: {script_path}")
    return(FALSE)
  }
  return(TRUE)
}

# Function to get batches to process
get_batches_to_process <- function(config, batch_arg = NULL) {
  multi_batch_mode <- tryCatch(
    isTRUE(config$parameters$normalization$multi_batch_mode),
    error = function(e) FALSE
  )

  if (!is.null(batch_arg)) {
    # Process specific batch if requested
    # Return as character vector (not list) for consistency
    return(c(batch_arg))
  }

  if (!multi_batch_mode) {
    # Single-batch mode: use default batch
    default_batch <- config$batch$default_batch_id %||% "batch_01"
    # Return as character vector (not list) for consistency
    return(c(default_batch))
  }

  # Multi-batch mode: get all batches from config
  batches <- names(config$batches)
  if (length(batches) == 0) {
    log_warn("Multi-batch mode enabled but no batches defined. Using default batch_01")
    # Return as character vector (not list) for consistency
    return(c("batch_01"))
  }

  # Sort batches to ensure consistent processing order
  batches <- sort(batches)
  log_info("Multi-batch mode: Processing {length(batches)} batches: {paste(batches, collapse=', ')}")
  return(batches)
}

# Function to create batch-specific directories
create_batch_directories <- function(config, batch_id = NULL) {
  base_dir <- config$output$base_dir
  if (is.null(base_dir)) {
    base_dir <- Sys.getenv("PIPELINE_OUTPUT_DIR", "output")
  }
  multi_batch_mode <- tryCatch(
    isTRUE(config$parameters$normalization$multi_batch_mode),
    error = function(e) FALSE
  )

  # Note: Aggregate outputs are saved to phenotypes/aggregate/ by scripts 09-11
  subdirs <- c("qc", "outliers", "normalized", "phenotypes", "pqtl", "reports", "logs", "pqtl-training")

  for (subdir in subdirs) {
    # Base directory
    dir_path <- file.path(base_dir, subdir)
    dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)

    # Batch-specific directory (if multi-batch mode)
    if (multi_batch_mode && !is.null(batch_id) && subdir != "aggregate") {
      batch_dir <- file.path(dir_path, batch_id)
      dir.create(batch_dir, recursive = TRUE, showWarnings = FALSE)
    }
  }
}

# Function to set batch context in environment
set_batch_context <- function(batch_id, config) {
  # Update config with current batch
  config$batch$current_batch_id <- batch_id

  # Set environment variables for scripts to access
  Sys.setenv(PIPELINE_BATCH_ID = batch_id)
  Sys.setenv(PIPELINE_CONFIG = args$config)

  # Store config in global environment (scripts will reload it)
  assign("PIPELINE_CONFIG_OBJ", config, envir = .GlobalEnv)

  return(config)
}

# Function to run a pipeline step for a specific batch
run_step_for_batch <- function(step_name, batch_id, config, dry_run = FALSE) {
  script_path <- file.path(script_dir, paste0(step_name, ".R"))

  if (!check_script(step_name)) {
    log_error("Cannot run step {step_name}: script not found")
    return(FALSE)
  }

  # Print step header (both to log and console in verbose mode)
  step_header <- paste0("\n", strrep("=", 60), "\n",
                       "Running: ", step_name, " for batch: ", batch_id, "\n",
                       "Description: ", pipeline_steps[[step_name]], "\n",
                       strrep("=", 60), "\n")
  cat(step_header)
  log_info(strrep("=", 60))
  log_info("Running: {step_name} for batch: {batch_id}")
  log_info("Description: {pipeline_steps[[step_name]]}")
  log_info(strrep("=", 60))

  if (dry_run) {
    log_info("[DRY RUN] Would execute: Rscript {script_path} (batch: {batch_id})")
    cat("[DRY RUN] Would execute:", script_path, "(batch:", batch_id, ")\n")
    return(TRUE)
  }

  # Set batch context
  config <- set_batch_context(batch_id, config)

  # Extract step number from step_name (e.g., "02_technical_outliers" -> "02")
  step_num_match <- regmatches(step_name, regexpr("^[0-9]{2}[a-z]?", step_name))
  step_num <- if (length(step_num_match) > 0) step_num_match[1] else "00"

  # Set step number in environment so scripts can access it
  Sys.setenv(PIPELINE_STEP_NUM = step_num)
  Sys.setenv(SCRIPT_NAME = script_path)

  # Create batch-specific directories
  create_batch_directories(config, batch_id)

  # Execute the script
  start_time <- Sys.time()

  tryCatch(
    {
      # Clear skip flag before sourcing (scripts can set this if they skip)
      Sys.unsetenv("PIPELINE_STEP_SKIPPED")

      # Source the script (it will read config and use batch context)
      source(script_path, local = FALSE)

      end_time <- Sys.time()
      duration <- difftime(end_time, start_time, units = "secs")

      # Check if step was skipped (scripts set PIPELINE_STEP_SKIPPED environment variable)
      was_skipped <- Sys.getenv("PIPELINE_STEP_SKIPPED", "FALSE") == "TRUE"

      if (was_skipped) {
        success_msg <- paste0("⊘ Step ", step_name, " skipped\n")
        cat(success_msg)
        log_info("Step {step_name} (batch {batch_id}) was skipped")
      } else {
        success_msg <- paste0("✓ Step ", step_name, " completed\n")
        cat(success_msg)
        log_info("Step {step_name} (batch {batch_id}) completed successfully in {round(duration, 2)} seconds")
      }
      return(TRUE)
    },
    error = function(e) {
      # Check if this is a skip signal (not a real error)
      if (e$message == "STEP_SKIPPED" || Sys.getenv("PIPELINE_STEP_SKIPPED", "FALSE") == "TRUE") {
        success_msg <- paste0("⊘ Step ", step_name, " skipped\n")
        cat(success_msg)
        log_info("Step {step_name} (batch {batch_id}) was skipped")
        return(TRUE)
      }
      # Real error - report it
      error_msg <- paste0("✗ Step ", step_name, " failed: ", e$message, "\n")
      cat(error_msg, file = stderr())
      log_error("Step {step_name} (batch {batch_id}) failed: {e$message}")
      log_error("Error traceback: {paste(capture.output(traceback()), collapse='\\n')}")
      return(FALSE)
    }
  )
}

# Function to run cross-batch step (after all batches complete prerequisite steps)
run_cross_batch_step <- function(step_name, batches, config, dry_run = FALSE) {
  step_header <- paste0("\n", strrep("=", 60), "\n",
                       "Running cross-batch step: ", step_name, "\n",
                       "Description: ", pipeline_steps[[step_name]], "\n",
                       "Batches involved: ", paste(batches, collapse = ", "), "\n",
                       strrep("=", 60), "\n")
  cat(step_header)
  log_info(strrep("=", 60))
  log_info("Running cross-batch step: {step_name}")
  log_info("Description: {pipeline_steps[[step_name]]}")
  log_info("Batches involved: {paste(batches, collapse=', ')}")
  log_info(strrep("=", 60))

  if (dry_run) {
    log_info("[DRY RUN] Would execute cross-batch step: {step_name}")
    return(TRUE)
  }

  # Create directories
  create_batch_directories(config)

  # Set batch context to first batch (scripts will handle cross-batch logic)
  config <- set_batch_context(batches[1], config)

  # Extract step number from step_name
  step_num_match <- regmatches(step_name, regexpr("^[0-9]{2}", step_name))
  step_num <- if (length(step_num_match) > 0) step_num_match[1] else "00"

  # Set step number in environment
  Sys.setenv(PIPELINE_STEP_NUM = step_num)
  script_path <- file.path(script_dir, paste0(step_name, ".R"))
  Sys.setenv(SCRIPT_NAME = script_path)

  start_time <- Sys.time()

  tryCatch(
    {
      # Clear skip flag before sourcing (scripts can set this if they skip)
      Sys.unsetenv("PIPELINE_STEP_SKIPPED")

      source(script_path, local = FALSE)

      end_time <- Sys.time()
      duration <- difftime(end_time, start_time, units = "secs")

      # Check if step was skipped
      was_skipped <- Sys.getenv("PIPELINE_STEP_SKIPPED", "FALSE") == "TRUE"

      if (was_skipped) {
        success_msg <- paste0("⊘ Cross-batch step ", step_name, " skipped\n")
        cat(success_msg)
        log_info("Cross-batch step {step_name} was skipped")
      } else {
        success_msg <- paste0("✓ Cross-batch step ", step_name, " completed\n")
        cat(success_msg)
        log_info("Cross-batch step {step_name} completed successfully in {round(duration, 2)} seconds")
      }
      return(TRUE)
    },
    error = function(e) {
      error_msg <- paste0("✗ Cross-batch step ", step_name, " failed: ", e$message, "\n")
      cat(error_msg, file = stderr())
      log_error("Cross-batch step {step_name} failed: {e$message}")
      log_error("Error traceback: {paste(capture.output(traceback()), collapse='\\n')}")
      return(FALSE)
    }
  )
}

# Function to get step range
get_step_range <- function(from = NULL, to = NULL) {
  all_steps <- names(pipeline_steps)

  if (is.null(from)) {
    from_idx <- 1
  } else {
    from_idx <- which(grepl(from, all_steps))
    if (length(from_idx) == 0) {
      stop(paste("Step not found:", from))
    }
    from_idx <- from_idx[1]
  }

  if (is.null(to)) {
    to_idx <- length(all_steps)
  } else {
    to_idx <- which(grepl(to, all_steps))
    if (length(to_idx) == 0) {
      stop(paste("Step not found:", to))
    }
    to_idx <- to_idx[1]
  }

  return(all_steps[from_idx:to_idx])
}

# Main execution
main <- function() {
  # Reload config (may have been modified)
  config <- read_yaml(args$config)

  # Determine mode
  multi_batch_mode <- tryCatch(
    isTRUE(config$parameters$normalization$multi_batch_mode),
    error = function(e) FALSE
  )

  aggregate_output <- tryCatch(
    isTRUE(config$parameters$aggregation$aggregate_output),
    error = function(e) FALSE
  )

  # Get batches to process
  batches <- get_batches_to_process(config, args$batch)

  # Determine which steps to run
  if (args$step == "all") {
    steps_to_run <- get_step_range(args$from, args$to)
  } else {
    matching_steps <- names(pipeline_steps)[grepl(args$step, names(pipeline_steps))]
    if (length(matching_steps) == 0) {
      stop(paste("No matching step found for:", args$step))
    }
    steps_to_run <- matching_steps
  }

  # Print initial pipeline information to console
  cat("\n", strrep("=", 60), "\n", sep = "")
  cat("FinnGen 3 Olink Analysis Pipeline (Refactored Version)\n")
  cat(strrep("=", 60), "\n", sep = "")
  cat("Mode:", if(multi_batch_mode) "Multi-batch" else "Single-batch", "\n")
  cat("Batches to process:", paste(batches, collapse = ", "), "\n")
  cat("Steps to run:", paste(steps_to_run, collapse = ", "), "\n")
  if (aggregate_output && multi_batch_mode) {
    cat("Aggregation enabled: Steps 09-11 will create aggregate outputs\n")
  }
  if (isTRUE(args$dry_run)) {
    cat("DRY RUN MODE - No actual execution\n")
  }
  cat(strrep("=", 60), "\n\n", sep = "")

  log_info("Mode: {if(multi_batch_mode) 'Multi-batch' else 'Single-batch'}")
  log_info("Batches to process: {paste(batches, collapse=', ')}")
  log_info("Steps to run: {paste(steps_to_run, collapse=', ')}")
  if (aggregate_output && multi_batch_mode) {
    log_info("Aggregation enabled: Steps 09-11 will create aggregate outputs")
  }
  if (isTRUE(args$dry_run)) {
    log_info("DRY RUN MODE - No actual execution")
  }

  # Track results
  batch_results <- list()
  success_count <- 0
  failed_steps <- character()

  if (multi_batch_mode && length(batches) > 1) {
    # MULTI-BATCH MODE: Process reference batch first, then other batches using QCed reference data

    # Get reference batch ID
    reference_batch_id <- tryCatch(
      config$parameters$normalization$reference_batch,
      error = function(e) NULL
    )
    if (is.null(reference_batch_id) || !reference_batch_id %in% batches) {
      # Default to first batch if not specified or invalid
      reference_batch_id <- batches[1]
      log_warn("Reference batch not specified or invalid, using: {reference_batch_id}")
    }
    log_info("Reference batch: {reference_batch_id}")
    cat("Reference batch:", reference_batch_id, "\n")

    # Separate reference batch from other batches
    other_batches <- setdiff(batches, reference_batch_id)
    other_batches_str <- if(length(other_batches) > 0) paste(other_batches, collapse = ", ") else "none"
    log_info("Other batches: {if(length(other_batches) > 0) paste(other_batches, collapse=', ') else 'none'}")
    cat("Other batches:", other_batches_str, "\n")

    # Phase 1: Process reference batch first through all independent steps (00-05d)
    # This ensures reference batch is fully QCed before being used by other batches
    # Note: Steps 06-07 are cross-batch, step 08 requires normalized data, steps 09-11 are aggregation
    independent_steps <- setdiff(steps_to_run, c(cross_batch_steps, post_normalization_steps, aggregation_steps))

    if (length(independent_steps) > 0) {
      phase1_msg <- paste0("\nPhase 1: Processing reference batch (", reference_batch_id, ") through all independent steps\n")
      cat(phase1_msg)
      log_info("Phase 1: Processing reference batch ({reference_batch_id}) through all independent steps")
      batch_results[[reference_batch_id]] <- list(success = 0, failed = character())

      for (step in independent_steps) {
        success <- run_step_for_batch(step, reference_batch_id, config, isTRUE(args$dry_run))

        if (success) {
          batch_results[[reference_batch_id]]$success <- batch_results[[reference_batch_id]]$success + 1
          success_count <- success_count + 1
        } else {
          failed_step <- paste0(step, " (", reference_batch_id, ")")
          batch_results[[reference_batch_id]]$failed <- c(batch_results[[reference_batch_id]]$failed, step)
          failed_steps <- c(failed_steps, failed_step)

          # Stop processing reference batch if it fails (other batches depend on it)
          log_error("Reference batch {reference_batch_id} failed at step {step}")
          if (!args$dry_run) {
            log_error("Cannot proceed with other batches - reference batch must complete successfully")
            break
          }
        }
      }
    }

    # Phase 2: Process other batches through independent steps (00-05d)
    # These batches can now use QCed data from reference batch
    if (length(other_batches) > 0 && length(independent_steps) > 0) {
      phase2_msg <- paste0("\nPhase 2: Processing other batches through independent steps (using QCed reference batch data)\n")
      cat(phase2_msg)
      log_info("Phase 2: Processing other batches through independent steps (using QCed reference batch data)")
      for (batch_id in other_batches) {
        cat("Processing batch:", batch_id, "\n")
        log_info("Processing batch: {batch_id}")
        batch_results[[batch_id]] <- list(success = 0, failed = character())

        for (step in independent_steps) {
          success <- run_step_for_batch(step, batch_id, config, isTRUE(args$dry_run))

          if (success) {
            batch_results[[batch_id]]$success <- batch_results[[batch_id]]$success + 1
            success_count <- success_count + 1
          } else {
            failed_step <- paste0(step, " (", batch_id, ")")
            batch_results[[batch_id]]$failed <- c(batch_results[[batch_id]]$failed, step)
            failed_steps <- c(failed_steps, failed_step)

            # Continue processing other batches even if one fails
            log_error("Batch {batch_id} failed at step {step}, but continuing with other batches")
          }
        }
      }
    }

    # Phase 3: Cross-batch steps (06-07) - run after all batches complete prerequisites
    cross_batch_to_run <- intersect(steps_to_run, cross_batch_steps)

    if (length(cross_batch_to_run) > 0) {
      phase3_msg <- paste0("\nPhase 3: Running cross-batch normalization steps (using QCed data from all batches)\n")
      cat(phase3_msg)
      log_info("Phase 3: Running cross-batch normalization steps (using QCed data from all batches)")
      for (step in cross_batch_to_run) {
        success <- run_cross_batch_step(step, batches, config, isTRUE(args$dry_run))
        if (success) {
          success_count <- success_count + 1
        } else {
          failed_steps <- c(failed_steps, step)
          log_error("Cross-batch step {step} failed")
        }
      }
    }

    # Phase 3.5: Normalization and post-normalization steps (06, 08) - run per-batch
    # Step 06 runs per-batch (can use cross-batch data for normalization)
    # Step 08 requires normalized data from step 06
    post_norm_to_run <- intersect(steps_to_run, post_normalization_steps)

    if (length(post_norm_to_run) > 0) {
      phase35_msg <- paste0("\nPhase 3.5: Running normalization and post-normalization steps (per-batch)\n")
      cat(phase35_msg)
      log_info("Phase 3.5: Running normalization and post-normalization steps (per-batch)")
      log_info("  Step 06 runs per-batch (can use cross-batch data for normalization)")
      log_info("  Step 08 requires normalized data from step 06")

      for (batch_id in batches) {
        for (step in post_norm_to_run) {
          success <- run_step_for_batch(step, batch_id, config, isTRUE(args$dry_run))
          if (success) {
            if (!batch_id %in% names(batch_results)) {
              batch_results[[batch_id]] <- list(success = 0, failed = character())
            }
            batch_results[[batch_id]]$success <- batch_results[[batch_id]]$success + 1
            success_count <- success_count + 1
          } else {
            failed_step <- paste0(step, " (", batch_id, ")")
            if (!batch_id %in% names(batch_results)) {
              batch_results[[batch_id]] <- list(success = 0, failed = character())
            }
            batch_results[[batch_id]]$failed <- c(batch_results[[batch_id]]$failed, step)
            failed_steps <- c(failed_steps, failed_step)
            log_error("Batch {batch_id} failed at step {step}, but continuing with other batches")
          }
        }
      }
    }

    # Phase 4: Steps 09-11 (phenotype preparation, kinship filtering, rank normalization)
    # CRITICAL FIX: Always run per-batch first, then scripts handle aggregation internally
    # When aggregate_output=true, scripts will automatically aggregate after processing their batch
    aggregation_steps_to_run <- intersect(steps_to_run, aggregation_steps)
    if (length(aggregation_steps_to_run) > 0) {
    if (aggregate_output) {
        phase4_msg <- paste0("\nPhase 4: Running steps 09-11 per-batch (scripts will aggregate outputs)\n")
        cat(phase4_msg)
        log_info("Phase 4: Running steps 09-11 per-batch (scripts will aggregate outputs)")
        log_info("  Scripts detect aggregate_output=true and merge batch outputs internally")
    } else {
        phase4_msg <- paste0("\nPhase 4: Processing batches through steps 09-11 (no aggregation)\n")
        cat(phase4_msg)
        log_info("Phase 4: Processing batches through steps 09-11 (no aggregation)")
      }

        for (batch_id in batches) {
          for (step in aggregation_steps_to_run) {
            success <- run_step_for_batch(step, batch_id, config, isTRUE(args$dry_run))
            if (success) {
            if (!batch_id %in% names(batch_results)) {
              batch_results[[batch_id]] <- list(success = 0, failed = character())
            }
              batch_results[[batch_id]]$success <- batch_results[[batch_id]]$success + 1
              success_count <- success_count + 1
            } else {
              failed_step <- paste0(step, " (", batch_id, ")")
            if (!batch_id %in% names(batch_results)) {
              batch_results[[batch_id]] <- list(success = 0, failed = character())
            }
            batch_results[[batch_id]]$failed <- c(batch_results[[batch_id]]$failed, step)
            failed_steps <- c(failed_steps, failed_step)
          }
        }
      }
    }
  } else {
    # SINGLE-BATCH MODE: Process one batch through all steps
    batch_id <- batches[1]
    cat("\nSingle-batch mode: Processing batch", batch_id, "\n\n")
    log_info("Single-batch mode: Processing batch {batch_id}")

    for (step in steps_to_run) {
      success <- run_step_for_batch(step, batch_id, config, isTRUE(args$dry_run))
      if (success) {
        success_count <- success_count + 1
      } else {
        failed_steps <- c(failed_steps, step)
        if (!args$dry_run) {
          log_error("Pipeline stopped due to failure in step: {step}")
          break
        }
      }
    }
  }

  # Print summary
  cat("\n", rep("=", 60), "\n", sep = "")
  cat("PIPELINE SUMMARY\n")
  cat(rep("=", 60), "\n", sep = "")
  cat("Mode:", if (multi_batch_mode) "Multi-batch" else "Single-batch", "\n")
  cat("Batches processed:", paste(batches, collapse = ", "), "\n")
  cat("Total steps attempted:", length(steps_to_run), "\n")
  cat("Successful steps:", success_count, "\n")
  cat("Failed steps:", length(failed_steps), "\n")

  if (multi_batch_mode && length(batch_results) > 0) {
    cat("\nPer-batch results:\n")
    for (batch_id in names(batch_results)) {
      cat("  Batch", batch_id, ":\n")
      cat("    Successful steps:", batch_results[[batch_id]]$success, "\n")
      if (length(batch_results[[batch_id]]$failed) > 0) {
        cat("    Failed steps:", paste(batch_results[[batch_id]]$failed, collapse = ", "), "\n")
      }
    }
  }

  if (length(failed_steps) > 0) {
    cat("\nFailed steps:\n")
    for (step in failed_steps) {
      cat("  -", step, "\n")
    }
  }

  if (success_count == length(steps_to_run)) {
    cat("\nPIPELINE COMPLETED SUCCESSFULLY!\n")
    log_info("Pipeline completed successfully")
  } else {
    cat("\nPIPELINE INCOMPLETE\n")
    log_warn("Pipeline incomplete - {length(failed_steps)} steps failed")
  }

  cat(rep("=", 60), "\n", sep = "")
}

# Helper function for NULL coalescing
`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

# Run pipeline
if (!interactive()) {
  main()
}




