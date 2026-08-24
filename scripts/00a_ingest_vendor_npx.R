#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------------------------------
# Step 00a: Vendor Extended NPX ingestion
#
# Purpose
#   Convert an Olink Explore HT vendor "Extended NPX" long-format delivery into the wide NPX matrix
#   that Step 00 (00_data_loader.R) consumes, while PRESERVING the sample- and assay-level metadata
#   that a naive long->wide pivot discards.
#
# Why this step exists
#   load_npx_matrix() in 00_data_loader.R can dcast a long-format file, but the pipeline then rebuilds
#   its long representation from the wide matrix (create_long_format_samples_data), fabricating
#   SampleQC = "PASS", AssayQC = "PASS" and SampleType = "SAMPLE". Feeding the raw vendor file straight
#   into Step 00 therefore silently loses: SampleType (so control wells are treated as samples), the
#   real vendor SampleQC/AssayQC verdicts, PlateID/WellID, per-sample QC warn/fail counters, assay
#   Count, and - new in NPX Map 2.0.0 - the four Olink Sample Index (OSI) columns.
#
#   This step emits those explicitly so no downstream stage has to reconstruct them.
#
# Inputs
#   --npx        Vendor Extended NPX parquet (long format)
#   --outdir     Output directory
#   --batch      Batch designation used in output filenames (e.g. fg3_batch_03)
#   --withdraw   Optional file of SampleIDs to hard-exclude (consent withdrawals), one per line
#
# Outputs (all prefixed 00a_, suffixed with the batch designation)
#   npx_matrix_raw          Wide biological NPX matrix, samples x biological assays (control probes removed)
#   pcnorm_matrix_raw       Wide biological PCNormalizedNPX matrix, same dimensions
#   controls_npx_matrix     Wide NPX matrix for control wells, retained for QC
#   sample_index            One row per biological sample: SampleType, PlateID, WellID, OSI, QC counters
#   assay_inventory         One row per delivered assay with control-probe and vendor-exclusion flags
#   lod_from_neg_controls   Per-assay LOD = mean(NC NPX) + 3 * SD(NC NPX)
#   withdrawal_register     SampleIDs excluded on consent grounds, with the reason
#   ingestion_summary       Counts for the release note and the validation gates
#
# Notes
#   - Identifier columns are read as character throughout. Integer coercion of Olink SampleIDs strips
#     significant leading zeros and produces silent zero-hit joins downstream.
#   - Assay columns are named by the Assay symbol. If the delivery contains duplicated Assay symbols
#     they are disambiguated as <Assay>_<Block>; the batch 03 delivery has none.
#   - Vendor-excluded assays (Count == 0 and NPX all NaN) are flagged, not dropped here. Step 00 applies
#     the named purge so the removal is explicit and logged in one place.
# ---------------------------------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(arrow)
  library(data.table)
  library(dplyr)
  library(optparse)
})

opt <- parse_args(OptionParser(option_list = list(
  make_option(c("-n", "--npx"), type = "character", help = "Vendor Extended NPX parquet (long format)"),
  make_option(c("-o", "--outdir"), type = "character", help = "Output directory"),
  make_option(c("-b", "--batch"), type = "character", default = "fg3_batch_03", help = "Batch designation"),
  make_option(c("-w", "--withdraw"), type = "character", default = NULL, help = "File of SampleIDs to hard-exclude")
)))

stopifnot(!is.null(opt$npx), !is.null(opt$outdir))
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

say <- function(...) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), paste0(...)))
op  <- function(name, ext) file.path(opt$outdir, sprintf("00a_%s_%s.%s", name, opt$batch, ext))

CONTROL_ASSAY_TYPES <- c("ext_ctrl", "inc_ctrl", "amp_ctrl")

say("opening ", opt$npx)
ds <- open_dataset(opt$npx)
say("total long-format rows: ", format(ds$num_rows, big.mark = ","))

# --- assay inventory -------------------------------------------------------------------------------
say("building assay inventory")
assays <- ds %>%
  select(OlinkID, UniProt, Assay, AssayType, Panel, Block) %>%
  distinct() %>% collect() %>% as.data.table()
setorder(assays, OlinkID)

assay_stats <- ds %>%
  filter(SampleType == "SAMPLE") %>%
  group_by(OlinkID) %>%
  summarise(n_obs        = sum(!is.na(NPX)),
            n_missing    = sum(is.na(NPX)),
            mean_count   = mean(Count, na.rm = TRUE),
            n_assayqc_warn = sum(AssayQC == "WARN", na.rm = TRUE)) %>%
  collect() %>% as.data.table()

assays <- merge(assays, assay_stats, by = "OlinkID", all.x = TRUE)
assays[, is_control_probe := AssayType %in% CONTROL_ASSAY_TYPES]
assays[, vendor_excluded  := !is_control_probe & n_obs == 0L]
assays[, missing_rate     := ifelse((n_obs + n_missing) > 0, n_missing / (n_obs + n_missing), NA_real_)]

# Assay symbols are unique in this delivery; disambiguate on Block only if that ever changes.
dups <- assays[is_control_probe == FALSE, .N, by = Assay][N > 1]$Assay
assays[, assay_colname := ifelse(Assay %in% dups, paste0(Assay, "_", Block), Assay)]

say("  assays delivered: ", nrow(assays),
    " | control probes: ", sum(assays$is_control_probe),
    " | vendor-excluded: ", sum(assays$vendor_excluded),
    " | duplicated symbols: ", length(dups))
fwrite(assays, op("assay_inventory", "tsv"), sep = "\t")

bio_assays <- assays[is_control_probe == FALSE]

# --- sample index ----------------------------------------------------------------------------------
say("building sample index")
# One row per sample. DataAnalysisRefID and the *Version fields vary by Block, so including them in the
# distinct() would emit one row per (sample, block) - eight rows per sample in an Explore HT delivery.
# Platform versions are a property of the delivery, not of a sample, so they are captured separately.
sidx <- ds %>%
  select(SampleID, SampleType, PlateID, WellID,
         OSICategory, OSISummary, OSITimeToCentrifugation, OSIPreparationTemperature) %>%
  distinct() %>% collect() %>% as.data.table()
sidx[, SampleID := as.character(SampleID)]

platform <- ds %>%
  select(SoftwareName, SoftwareVersion, PanelDataArchiveVersion, PreProcessingVersion,
         PreProcessingSoftware, InstrumentType) %>%
  distinct() %>% collect() %>% as.data.table()
fwrite(platform, op("platform_versions", "tsv"), sep = "\t")
say("  platform version combinations in delivery: ", nrow(platform))

qc_counters <- ds %>%
  group_by(SampleID) %>%
  summarise(n_sampleqc_fail    = sum(SampleQC == "FAIL", na.rm = TRUE),
            n_sampleqc_warn    = sum(SampleQC == "WARN", na.rm = TRUE),
            n_block_qc_fail    = sum(BlockQCFail, na.rm = TRUE),
            n_sampleblock_warn = sum(SampleBlockQCWarn, na.rm = TRUE),
            n_sampleblock_fail = sum(SampleBlockQCFail, na.rm = TRUE)) %>%
  collect() %>% as.data.table()
qc_counters[, SampleID := as.character(SampleID)]
sidx <- merge(sidx, qc_counters, by = "SampleID", all.x = TRUE)

# Vendor sample-level verdict: FAIL dominates WARN dominates PASS, across all blocks.
sidx[, vendor_sample_verdict := fifelse(n_sampleqc_fail > 0, "FAIL",
                                fifelse(n_sampleqc_warn > 0, "WARN", "PASS"))]

# derive the plate-name family, which is inconsistent across the two vendor sub-projects
sidx[, vendor_subproject := fifelse(grepl("90-1324676138", PlateID), "90-1324676138",
                             fifelse(grepl("90-1265725744", PlateID), "90-1265725744", "legacy_plate_label"))]
sidx[, plate_name_family := fifelse(grepl("^90-[0-9]+_B[0-9]+_P[AB]", PlateID), "subproject_underscore",
                             fifelse(grepl("^90-[0-9]+-B[0-9]+-P[AB]", PlateID), "subproject_hyphen",
                             fifelse(grepl("^Plate Layout", PlateID), "plate_layout_prefix", "other")))]

say("  distinct (SampleID, PlateID) rows: ", nrow(sidx))
say("  SampleType breakdown:")
print(sidx[, .(n_sampleIDs = uniqueN(SampleID)), by = SampleType][order(-n_sampleIDs)])

# --- withdrawals -----------------------------------------------------------------------------------
withdrawn <- character(0)
if (!is.null(opt$withdraw) && file.exists(opt$withdraw)) {
  withdrawn <- trimws(readLines(opt$withdraw))
  withdrawn <- withdrawn[nzchar(withdrawn) & !startsWith(withdrawn, "#")]
  say("withdrawals to hard-exclude: ", length(withdrawn))
}
reg <- data.table(SampleID = withdrawn,
                  reason = "consent_withdrawal_biobank_denial",
                  present_in_delivery = withdrawn %in% sidx$SampleID)
fwrite(reg, op("withdrawal_register", "tsv"), sep = "\t")
if (nrow(reg) && any(!reg$present_in_delivery))
  say("  WARNING: withdrawal IDs not found in delivery: ",
      paste(reg[present_in_delivery == FALSE]$SampleID, collapse = ", "))

bio_ids  <- setdiff(unique(sidx[SampleType == "SAMPLE"]$SampleID), withdrawn)
ctrl_ids <- unique(sidx[SampleType != "SAMPLE"]$SampleID)
say("biological samples after withdrawal removal: ", length(bio_ids), " | control wells: ", length(ctrl_ids))

fwrite(sidx[SampleID %in% bio_ids], op("sample_index", "tsv"), sep = "\t")
fwrite(sidx[SampleID %in% ctrl_ids], op("control_index", "tsv"), sep = "\t")

# --- LOD from negative controls, computed before controls are dropped ------------------------------
say("computing LOD from negative controls")
nc <- ds %>%
  filter(SampleType == "NEGATIVE_CONTROL") %>%
  select(SampleID, OlinkID, NPX) %>% collect() %>% as.data.table()
lod <- nc[!is.na(NPX), .(n_nc_obs = .N, nc_mean = mean(NPX), nc_sd = sd(NPX)), by = OlinkID]
lod[, LOD := nc_mean + 3 * nc_sd]
lod <- merge(lod, assays[, .(OlinkID, Assay, assay_colname, is_control_probe, vendor_excluded)],
             by = "OlinkID", all.y = TRUE)
fwrite(lod, op("lod_from_neg_controls", "tsv"), sep = "\t")
say("  LOD computed for ", lod[!is.na(LOD), .N], " of ", nrow(lod), " assays")

# --- wide matrices ---------------------------------------------------------------------------------
build_wide <- function(value_col, ids, label) {
  say("pivoting ", value_col, " for ", label, " (", length(ids), " samples)")
  d <- ds %>%
    select(SampleID, OlinkID, !!rlang::sym(value_col)) %>%
    collect() %>% as.data.table()
  setnames(d, value_col, "V")
  d[, SampleID := as.character(SampleID)]
  d <- d[SampleID %in% ids & OlinkID %in% bio_assays$OlinkID]
  # a (SampleID, OlinkID) pair is unique for biological samples; controls recur across plates
  w <- dcast(d, SampleID ~ OlinkID, value.var = "V", fun.aggregate = function(x) mean(x, na.rm = TRUE))
  ids_out <- w$SampleID
  w[, SampleID := NULL]
  cn <- bio_assays$assay_colname[match(names(w), bio_assays$OlinkID)]
  m <- as.matrix(w); rownames(m) <- ids_out; colnames(m) <- cn
  m[!is.finite(m)] <- NA_real_
  say("  -> ", nrow(m), " x ", ncol(m))
  m
}

npx_bio    <- build_wide("NPX", bio_ids, "biological")
pcnorm_bio <- build_wide("PCNormalizedNPX", bio_ids, "biological")
npx_ctrl   <- build_wide("NPX", ctrl_ids, "controls")

saveRDS(npx_bio,    op("npx_matrix_raw", "rds"))
saveRDS(pcnorm_bio, op("pcnorm_matrix_raw", "rds"))
saveRDS(npx_ctrl,   op("controls_npx_matrix", "rds"))

# --- summary ---------------------------------------------------------------------------------------
summ <- data.table(
  metric = c("long_format_rows", "distinct_sampleid_all", "biological_samples_delivered",
             "withdrawals_removed", "biological_samples_retained", "control_wells",
             "plates", "assays_delivered", "control_probes_removed", "biological_assays",
             "vendor_excluded_assays", "biological_assays_with_data",
             "vendor_sample_verdict_PASS", "vendor_sample_verdict_WARN", "vendor_sample_verdict_FAIL",
             "matrix_cells", "matrix_missing_cells", "matrix_missing_pct"),
  value = c(ds$num_rows,
            uniqueN(sidx$SampleID),
            uniqueN(sidx[SampleType == "SAMPLE"]$SampleID),
            length(intersect(withdrawn, sidx$SampleID)),
            nrow(npx_bio),
            length(ctrl_ids),
            uniqueN(sidx$PlateID),
            nrow(assays),
            sum(assays$is_control_probe),
            nrow(bio_assays),
            sum(assays$vendor_excluded),
            nrow(bio_assays) - sum(assays$vendor_excluded),
            sidx[SampleID %in% bio_ids & vendor_sample_verdict == "PASS", uniqueN(SampleID)],
            sidx[SampleID %in% bio_ids & vendor_sample_verdict == "WARN", uniqueN(SampleID)],
            sidx[SampleID %in% bio_ids & vendor_sample_verdict == "FAIL", uniqueN(SampleID)],
            length(npx_bio),
            sum(is.na(npx_bio)),
            round(100 * sum(is.na(npx_bio)) / length(npx_bio), 3)))
fwrite(summ, op("ingestion_summary", "tsv"), sep = "\t")
say("ingestion summary")
print(summ)
say("done; outputs in ", opt$outdir)
