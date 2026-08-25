#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------------------------------
# Step 00b: Build the Batch 03 sample metadata table
#
# Purpose
#   Batch 03 shipped no provider metadata TSV. Batches 01 and 02 each shipped one, and Step 00
#   (create_sample_mapping) merges four columns from it unconditionally:
#       SAMPLE_ID, FINNGENID, COHORT_FINNGENID, BIOBANK_PLASMA
#   The THL plate-2 volumes file supplies none of the last two and covers only 5,331 of 6,191 samples,
#   so Step 00 aborts with "object 'COHORT_FINNGENID' not found".
#
#   This step assembles an equivalent table covering EVERY delivered biological sample, from the
#   sources that do exist. It resolves three defects at once:
#     - the missing COHORT_FINNGENID / BIOBANK_PLASMA columns (Step 00 abort)
#     - the 860 P1516_* samples having no THL metadata row, which would otherwise be labelled Unknown
#       and silently dropped (13.9% of the batch)
#     - the two "-rep" technical replicates, whose base SampleID resolves but whose suffixed form does
#       not, reproducing the failure mode that lost 50 EA5 bridge samples in Batch 02
#
# Sources, in precedence order per field
#   SampleID / PlateID / WellID / OSI   step 00a sample index (authoritative: the delivered matrix)
#   FINNGENID                           THL combined crosswalk, with "-rep" resolved to its base ID
#   COHORT_FINNGENID                    R14 minimum_extended COHORT (covers all 4,635 individuals)
#   BIOBANK_PLASMA                      THL volumes BIOBANK; else plasma-selection BIOBANK_PLASMA;
#                                       else derived from COHORT for the Arctic/NFBC1966 stratum
#   pre-analytical fields               THL volumes file, then the Jan-2026 per-biobank availability
#                                       lists (which add CONTAINER_NAME and fill some VOLUME gaps)
#   SEX / BL_AGE / BMI / SMOKE3         R14 minimum_extended
#   disease-group eligibility           the derived registry-eligibility table (NOT provider labels)
#
# Outputs
#   00b_metadata_fg3_batch_03.tsv        pipeline-facing metadata (Step 00 consumes this)
#   00b_bridging_metadata_fg3_batch_03.tsv  bridge table carrying an in_Batch_03 column
#   00b_metadata_build_report.txt        provenance and coverage counts for the release note
#
# Note on identifiers: every ID column is read and written as character. Integer coercion of Olink
# SampleIDs strips significant leading zeros and produces silent zero-hit joins.
# ---------------------------------------------------------------------------------------------------

suppressPackageStartupMessages({ library(arrow); library(data.table); library(optparse) })

opt <- parse_args(OptionParser(option_list = list(
  make_option("--sample_index", type = "character"),
  make_option("--crosswalk", type = "character"),
  make_option("--volumes", type = "character"),
  make_option("--availability_dir", type = "character", default = NULL),
  make_option("--selection_dir", type = "character", default = NULL),
  make_option("--r14_minimum", type = "character"),
  make_option("--eligibility", type = "character", default = NULL),
  make_option("--bridging_in", type = "character", default = NULL),
  make_option("--outdir", type = "character"),
  make_option("--batch", type = "character", default = "fg3_batch_03")
)))

say <- function(...) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), paste0(...)))
op  <- function(n, e) file.path(opt$outdir, sprintf("00b_%s_%s.%s", n, opt$batch, e))
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)
CH <- "character"

# ---- delivered samples (authoritative row set) ----------------------------------------------------
si <- fread(opt$sample_index, colClasses = list(character = "SampleID"))
say("delivered biological samples: ", nrow(si))

M <- si[, .(SAMPLE_ID = SampleID, PlateID, WellID,
            OSICategory, OSISummary, OSITimeToCentrifugation, OSIPreparationTemperature,
            vendor_sample_verdict, vendor_subproject)]

# ---- FINNGENID, with -rep resolution --------------------------------------------------------------
xw <- fread(opt$crosswalk, colClasses = CH)
M[, base_sample_id := sub("-rep$", "", SAMPLE_ID)]
M[, is_technical_replicate := grepl("-rep$", SAMPLE_ID)]
M[, replicate_group_id := base_sample_id]
M <- merge(M, unique(xw[, .(base_sample_id = SAMPLE_ID, FINNGENID)]),
           by = "base_sample_id", all.x = TRUE)
say("FINNGENID resolved: ", sum(!is.na(M$FINNGENID)), " / ", nrow(M),
    "  (of which via -rep suffix stripping: ", M[is_technical_replicate == TRUE, .N], ")")
unres <- M[is.na(FINNGENID)]$SAMPLE_ID
if (length(unres)) say("UNRESOLVED (expected 0): ", paste(unres, collapse = ", "))

# ---- stratum -------------------------------------------------------------------------------------
M[, stratum := fifelse(grepl("^P1516_", SAMPLE_ID), "NFBC1966", "THL_routed")]

# ---- R14 cohort and covariates -------------------------------------------------------------------
r14 <- as.data.table(read_parquet(opt$r14_minimum,
        col_select = c("FINNGENID", "COHORT", "SEX", "BL_AGE", "BL_YEAR", "BMI", "SMOKE3")))
setnames(r14, "COHORT", "COHORT_FINNGENID")
M <- merge(M, r14, by = "FINNGENID", all.x = TRUE)
say("COHORT_FINNGENID populated: ", sum(!is.na(M$COHORT_FINNGENID)), " / ", nrow(M))

# ---- THL pre-analytical fields -------------------------------------------------------------------
v <- fread(opt$volumes, colClasses = CH, na.strings = c("NA", "", "empty"))
v <- v[!is.na(SAMPLE_ID)]
# two SAMPLE_IDs appear twice (genuine replicate aliquots); keep the first deterministically
v <- v[order(SAMPLE_ID, CONTAINER, POSITION)][!duplicated(SAMPLE_ID)]
vcols <- c("BIOBANK", "APPROX_TIMESTAMP_COLLECTION", "APPROX_TIMESTAMP_PROCESSING",
           "APPROX_TIMESTAMP_FREEZING", "BARCODE", "FREEZE_THAW_CYCLES", "HEMOLYSIS",
           "VOLUME", "CONTAINER", "CONTAINER_BARCODE", "POSITION")
M <- merge(M, v[, c("SAMPLE_ID", vcols), with = FALSE], by = "SAMPLE_ID", all.x = TRUE)

# ---- availability lists: CONTAINER_NAME and gap-filling -------------------------------------------
if (!is.null(opt$availability_dir) && dir.exists(opt$availability_dir)) {
  af <- list.files(opt$availability_dir, pattern = "^PROTEOMICS3_", full.names = TRUE)
  if (length(af)) {
    A <- unique(rbindlist(lapply(af, fread, colClasses = CH, na.strings = c("NA", "")), fill = TRUE),
                by = "SAMPLE_ID")
    M <- merge(M, A[, .(SAMPLE_ID, CONTAINER_NAME,
                        VOLUME_avail = VOLUME, FTC_avail = FREEZE_THAW_CYCLES)],
               by = "SAMPLE_ID", all.x = TRUE)
    M[is.na(VOLUME), VOLUME := VOLUME_avail]
    M[is.na(FREEZE_THAW_CYCLES), FREEZE_THAW_CYCLES := FTC_avail]
    M[, c("VOLUME_avail", "FTC_avail") := NULL]
    say("availability lists merged: CONTAINER_NAME on ", sum(!is.na(M$CONTAINER_NAME)), " samples")
  }
}

# ---- BIOBANK_PLASMA ------------------------------------------------------------------------------
M[, BIOBANK_PLASMA := BIOBANK]
if (!is.null(opt$selection_dir) && dir.exists(opt$selection_dir)) {
  sf <- list.files(opt$selection_dir, pattern = "plasma_selection|plasma_replacements",
                   full.names = TRUE)
  if (length(sf)) {
    S <- rbindlist(lapply(sf, fread, colClasses = CH), fill = TRUE)
    S <- unique(S[!is.na(BIOBANK_PLASMA), .(FINNGENID, BP_sel = BIOBANK_PLASMA)], by = "FINNGENID")
    M <- merge(M, S, by = "FINNGENID", all.x = TRUE)
    M[is.na(BIOBANK_PLASMA), BIOBANK_PLASMA := BP_sel]
    M[, BP_sel := NULL]
  }
}
# NFBC1966 is held by Arctic Biobank and never routed through THL, so it has no BIOBANK row
M[is.na(BIOBANK_PLASMA) & grepl("^ARCTIC BIOBANK", COHORT_FINNGENID), BIOBANK_PLASMA := "ARCTIC BIOBANK"]
M[is.na(BIOBANK_PLASMA), BIOBANK_PLASMA := "UNKNOWN"]
say("BIOBANK_PLASMA populated: ", sum(M$BIOBANK_PLASMA != "UNKNOWN"), " / ", nrow(M))

# ---- separate-project tag (Arctic Biobank, THL-routed) -------------------------------------------
M[, separate_project := fifelse(!is.na(BIOBANK) & BIOBANK == "Arctic Biobank", "ArcticBiobank", NA_character_)]

# ---- designated bridges --------------------------------------------------------------------------
bridge_fgids <- unique(v[BIOBANK == "BRIDGING_SAMPLE"]$FINNGENID)
M[, is_designated_bridge := FINNGENID %in% bridge_fgids]
say("designated bridge samples: ", M[is_designated_bridge == TRUE, .N],
    " rows / ", uniqueN(M[is_designated_bridge == TRUE]$FINNGENID), " individuals")

# ---- longitudinal structure ----------------------------------------------------------------------
M[, collection_date := as.Date(substr(APPROX_TIMESTAMP_COLLECTION, 1, 10))]
setorder(M, FINNGENID, collection_date, SAMPLE_ID)
M[, n_samples_this_individual := .N, by = FINNGENID]
M[, sample_index_within_individual := seq_len(.N), by = FINNGENID]
M[, days_since_previous_sample := as.integer(collection_date - shift(collection_date)), by = FINNGENID]
M[, collection_span_days_this_individual :=
    as.integer(max(collection_date, na.rm = TRUE) - min(collection_date, na.rm = TRUE)), by = FINNGENID]
M[is.infinite(collection_span_days_this_individual), collection_span_days_this_individual := NA_integer_]

# ---- metadata completeness -----------------------------------------------------------------------
M[, metadata_completeness := fifelse(is.na(BIOBANK), "thl_absent", "thl_complete")]
M[, QC_technical_evaluable := metadata_completeness == "thl_complete" &
                              !is.na(APPROX_TIMESTAMP_COLLECTION)]

# ---- registry-derived disease-group eligibility (NOT provider labels) ----------------------------
if (!is.null(opt$eligibility) && file.exists(opt$eligibility)) {
  E <- fread(opt$eligibility, colClasses = CH)
  ecols <- grep("^elig_", names(E), value = TRUE)
  E2 <- unique(E[, c("SAMPLE_ID", ecols, "n_eligible_groups", "eligible_groups",
                     "assignment_confidence"), with = FALSE], by = "SAMPLE_ID")
  M <- merge(M, E2, by = "SAMPLE_ID", all.x = TRUE)
  # Disease_Group is populated ONLY where registry eligibility is unambiguous. Provider selection
  # lists were not available for batch 03; see METADATA_GAPS_CHECKLIST.md.
  M[, Disease_Group := fifelse(assignment_confidence == "high_unique", eligible_groups, NA_character_)]
  M[, Disease_Group_source := fifelse(assignment_confidence == "high_unique",
                                      "derived_R14_registry_unique", "unavailable")]
  say("Disease_Group populated (high-confidence only): ", sum(!is.na(M$Disease_Group)), " / ", nrow(M))
}

# ---- write ---------------------------------------------------------------------------------------
M[, base_sample_id := NULL]
setcolorder(M, c("SAMPLE_ID", "FINNGENID", "COHORT_FINNGENID", "BIOBANK_PLASMA"))
setorder(M, SAMPLE_ID)
fwrite(M, op("metadata", "tsv"), sep = "\t")
say("wrote ", op("metadata", "tsv"), " : ", nrow(M), " rows x ", ncol(M), " cols")

# bridging table with an in_Batch_03 column, for the generalised Step 00 bridge branch
if (!is.null(opt$bridging_in) && file.exists(opt$bridging_in)) {
  B <- fread(opt$bridging_in, colClasses = CH)
  B[, in_Batch_03 := FINNGENID %in% bridge_fgids]
  extra <- setdiff(bridge_fgids, B$FINNGENID)
  if (length(extra)) B <- rbind(B, data.table(FINNGENID = extra, in_Batch_03 = TRUE), fill = TRUE)
  fwrite(B, op("bridging_metadata", "tsv"), sep = "\t")
  say("wrote bridging metadata: ", nrow(B), " rows, in_Batch_03 TRUE for ", sum(B$in_Batch_03))
}

rep_lines <- c(
  "Batch 03 metadata build report",
  paste0("generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  "",
  paste0("delivered biological samples        : ", nrow(M)),
  paste0("unique individuals                  : ", uniqueN(M$FINNGENID)),
  paste0("FINNGENID resolved                  : ", sum(!is.na(M$FINNGENID))),
  paste0("  of which via -rep suffix          : ", M[is_technical_replicate == TRUE, .N]),
  paste0("COHORT_FINNGENID populated          : ", sum(!is.na(M$COHORT_FINNGENID))),
  paste0("BIOBANK_PLASMA populated            : ", sum(M$BIOBANK_PLASMA != "UNKNOWN")),
  paste0("THL pre-analytical metadata present : ", M[metadata_completeness == "thl_complete", .N]),
  paste0("THL pre-analytical metadata absent  : ", M[metadata_completeness == "thl_absent", .N]),
  paste0("QC_technical_evaluable              : ", sum(M$QC_technical_evaluable)),
  paste0("stratum NFBC1966                    : ", M[stratum == "NFBC1966", .N]),
  paste0("stratum THL_routed                  : ", M[stratum == "THL_routed", .N]),
  paste0("designated bridge rows              : ", M[is_designated_bridge == TRUE, .N]),
  paste0("separate_project = ArcticBiobank    : ", M[!is.na(separate_project), .N]),
  paste0("technical replicates                : ", M[is_technical_replicate == TRUE, .N]),
  paste0("individuals with >1 sample          : ", uniqueN(M[n_samples_this_individual > 1]$FINNGENID)),
  paste0("Disease_Group populated             : ", sum(!is.na(M$Disease_Group))),
  "",
  "COHORT_FINNGENID composition:",
  paste0("  ", capture.output(print(M[, .N, by = COHORT_FINNGENID][order(-N)][1:15]))),
  "",
  "BIOBANK_PLASMA composition:",
  paste0("  ", capture.output(print(M[, .N, by = BIOBANK_PLASMA][order(-N)]))))
writeLines(rep_lines, file.path(opt$outdir, "00b_metadata_build_report.txt"))
say("done")
