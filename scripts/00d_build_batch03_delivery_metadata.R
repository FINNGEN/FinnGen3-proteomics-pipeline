#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------------------------------
# Step 00d: Build the single Batch 03 delivery metadata table
#
# Supersedes the earlier 00b metadata, which had to infer disease groups from R14 registry endpoints
# because the provider selection table had not been located. It has now been provided:
#
#   gs://fg-3/proteomics/genewiz_batch3/FG3_Batch_3_Plasma_Selections/Batch_3_final_plasma_selection_n_5308.tsv
#
# That file carries `Sample_set` — the authoritative task-force designation — for all 5,308 selected
# samples, plus a much richer pre-analytical block than the THL plate-2 volumes file (collection /
# processing / freezing times, hours-to-freezing, hemolysis, serial-sample flags, plasma index).
#
# The earlier 4-column `..._for_biobanks_set1.tsv` in the same directory is the stripped export that was
# sent TO the biobanks so they could pull tubes; it is not the selection record.
#
# Task-force assignment, in strict precedence order (recorded per row in Sample_set_source):
#   1. selection_exact      selection table matched on (FINNGENID, APPROX_DATE_COLLECTION)
#   2. selection_finngenid  matched on FINNGENID alone, used only where that individual has exactly one
#                           Sample_set in the selection table (true for 3,657 of 3,664 individuals)
#   3. supplement_107       the 107 Parkinson's priority-3 re-picks that replaced unavailable THL samples
#                           (TASK_FORCE column of the n_107 file)
#   4. replacement_2        the 2 Auria replacements, which the selection Notes record as drawn from the
#                           Atopic Dermatitis group
#   5. bridging             designated cross-batch bridging samples (THL BIOBANK = BRIDGING_SAMPLE)
#   6. nfbc1966_biobank     the 860 NFBC1966 samples, which Arctic Biobank selected and funded itself on
#                           the basis of incident disease. They are not part of the nine task forces and
#                           are retained in the release as their own stratum.
#   7. unassigned           no source resolves the row; Sample_set left NA
#
# Everything is kept in one table. Suitability for any given analysis is left to the analyst, per the
# project decision to deliver the batch whole rather than pre-filtered.
# ---------------------------------------------------------------------------------------------------

suppressPackageStartupMessages({ library(arrow); library(data.table); library(optparse) })

opt <- parse_args(OptionParser(option_list = list(
  make_option("--sample_index", type = "character"),
  make_option("--crosswalk", type = "character"),
  make_option("--volumes", type = "character"),
  make_option("--selection", type = "character", help = "Batch_3_final_plasma_selection_n_5308.tsv"),
  make_option("--supplement_107", type = "character"),
  make_option("--replacements_2", type = "character"),
  make_option("--r14_minimum", type = "character"),
  make_option("--qc_metadata", type = "character", default = NULL, help = "05d annotated metadata"),
  make_option("--outdir", type = "character"),
  make_option("--batch", type = "character", default = "fg3_batch_03")
)))
say <- function(...) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), paste0(...)))
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)
CH <- "character"; NAS <- c("NA", "", "empty")

# ---- 1. authoritative row set: every delivered biological sample -----------------------------------
si <- fread(opt$sample_index, colClasses = list(character = "SampleID"))
M <- si[, .(SAMPLE_ID = SampleID, PlateID, WellID, vendor_subproject, vendor_sample_verdict,
            OSICategory, OSISummary, OSITimeToCentrifugation, OSIPreparationTemperature)]
M[, base_sample_id := sub("-rep$", "", SAMPLE_ID)]
M[, is_technical_replicate := grepl("-rep$", SAMPLE_ID)]
M[, replicate_group_id := base_sample_id]
say("delivered biological samples: ", nrow(M))

# ---- 2. FINNGENID, resolving the -rep suffix ------------------------------------------------------
xw <- fread(opt$crosswalk, colClasses = CH)
M <- merge(M, unique(xw[, .(base_sample_id = SAMPLE_ID, FINNGENID)]), by = "base_sample_id", all.x = TRUE)
M[, stratum := fifelse(grepl("^P1516_", SAMPLE_ID), "NFBC1966", "THL_routed")]
say("FINNGENID resolved: ", sum(!is.na(M$FINNGENID)), " / ", nrow(M))

# ---- 3. THL pre-analytical (collection date is the join key to the selection table) ---------------
v <- fread(opt$volumes, colClasses = CH, na.strings = NAS)
v <- v[!is.na(SAMPLE_ID)][order(SAMPLE_ID, CONTAINER, POSITION)][!duplicated(SAMPLE_ID)]
M <- merge(M, v[, .(base_sample_id = SAMPLE_ID, BIOBANK, APPROX_TIMESTAMP_COLLECTION,
                    BARCODE, FREEZE_THAW_CYCLES, VOLUME, CONTAINER, CONTAINER_BARCODE, POSITION)],
           by = "base_sample_id", all.x = TRUE)
M[, collection_date := substr(APPROX_TIMESTAMP_COLLECTION, 1, 10)]

# ---- 4. the provider selection table — Sample_set and the rich pre-analytical block ---------------
S <- fread(opt$selection, colClasses = CH, na.strings = NAS)
sel_cols <- c("Sample_set", "PLASMA_INDEX.x", "duplicated_flag", "SERIAL_SAMPLES",
              "TIME_COLLECTION", "APPROX_DATE_PROCESSING", "TIME_PROCESSING",
              "APPROX_DATE_FREEZING", "TIME_FREEZING",
              "HOURS_FROM_COLLECTION_TO_FREEZING", "HOURS_FROM_COLLECTION_TO_FREEZING_BIN",
              "HEMOLYSIS", "FINNGENID_N_PLASMA_SAMPLES",
              "APPROX_BIRTH_DATE", "BL_YEAR", "BL_AGE", "SEX", "HEIGHT", "WEIGHT", "BMI",
              "SMOKE2", "SMOKE3", "SMOKE5", "CURRENT_SMOKER", "EVER_SMOKER",
              "regionofbirthname", "NUMBER_OF_OFFSPRING", "COHORT",
              "DEATH", "DEATH_FU_AGE", "AGE_AT_DEATH_OR_END_OF_FOLLOWUP")
sel_cols <- intersect(sel_cols, names(S))

# 4a. exact match on (FINNGENID, collection date)
Sx <- unique(S[, c("FINNGENID", "APPROX_DATE_COLLECTION", sel_cols), with = FALSE],
             by = c("FINNGENID", "APPROX_DATE_COLLECTION"))
M <- merge(M, Sx, by.x = c("FINNGENID", "collection_date"),
           by.y = c("FINNGENID", "APPROX_DATE_COLLECTION"), all.x = TRUE)
M[, Sample_set_source := fifelse(!is.na(Sample_set), "selection_exact", NA_character_)]
say("task force by exact (FINNGENID, date) match: ", sum(!is.na(M$Sample_set)))

# 4b. fall back to FINNGENID, only where the individual has a single Sample_set
uni <- S[, .(n = uniqueN(Sample_set), ss = Sample_set[1]), by = FINNGENID][n == 1, .(FINNGENID, ss)]
M <- merge(M, uni, by = "FINNGENID", all.x = TRUE)
n0 <- sum(!is.na(M$Sample_set))
M[is.na(Sample_set) & !is.na(ss), `:=`(Sample_set = ss, Sample_set_source = "selection_finngenid")]
say("  + by FINNGENID (single-task-force individuals): ", sum(!is.na(M$Sample_set)) - n0)
M[, ss := NULL]

# ---- 5. the two documented substitution rounds ----------------------------------------------------
A <- fread(opt$supplement_107, colClasses = CH, na.strings = NAS)
n0 <- sum(!is.na(M$Sample_set))
M[is.na(Sample_set) & FINNGENID %in% A$FINNGENID,
  `:=`(Sample_set = "Parkinsons", Sample_set_source = "supplement_107")]
say("  + 107 Parkinson's priority-3 supplement: ", sum(!is.na(M$Sample_set)) - n0)

R <- fread(opt$replacements_2, colClasses = CH, na.strings = NAS)
repl <- R[REMOVED_REPLACEMENT == "REPLACEMENT"]$FINNGENID
n0 <- sum(!is.na(M$Sample_set))
M[is.na(Sample_set) & FINNGENID %in% repl,
  `:=`(Sample_set = "Atopic_dermatitis", Sample_set_source = "replacement_2")]
say("  + 2 Auria replacements (Atopic Dermatitis, per Notes): ", sum(!is.na(M$Sample_set)) - n0)
M[, removed_replacement := R$REMOVED_REPLACEMENT[match(FINNGENID, R$FINNGENID)]]

# ---- 6. bridging and the NFBC1966 stratum ---------------------------------------------------------
bridge_fg <- unique(v[BIOBANK == "BRIDGING_SAMPLE"]$FINNGENID)
M[, is_designated_bridge := FINNGENID %in% bridge_fg]
n0 <- sum(!is.na(M$Sample_set))
M[is.na(Sample_set) & is_designated_bridge == TRUE,
  `:=`(Sample_set = "Bridging", Sample_set_source = "bridging")]
say("  + designated bridging samples: ", sum(!is.na(M$Sample_set)) - n0)

n0 <- sum(!is.na(M$Sample_set))
M[is.na(Sample_set) & stratum == "NFBC1966",
  `:=`(Sample_set = "NFBC1966", Sample_set_source = "nfbc1966_biobank")]
say("  + NFBC1966 biobank-selected stratum: ", sum(!is.na(M$Sample_set)) - n0)

M[is.na(Sample_set), Sample_set_source := "unassigned"]
M[, is_task_force_sample := Sample_set %in%
    c("Pulmo","Metabolic","Rheuma","IBD","Eye","Heart failure","Atopic_dermatitis",
      "Kidney","Parkinsons","Psychosis")]

# ---- 7. R14 demographics, filling the NFBC1966 arm which is absent from the selection table -------
r14 <- as.data.table(read_parquet(opt$r14_minimum,
        col_select = c("FINNGENID","COHORT","SEX","BL_AGE","BL_YEAR","BMI","SMOKE3","APPROX_BIRTH_DATE")))
setnames(r14, setdiff(names(r14), "FINNGENID"), paste0("r14_", setdiff(names(r14), "FINNGENID")))
M <- merge(M, r14, by = "FINNGENID", all.x = TRUE)
for (c in c("COHORT","SEX","BL_AGE","BL_YEAR","BMI","SMOKE3","APPROX_BIRTH_DATE")) {
  if (c %in% names(M)) M[is.na(get(c)), (c) := get(paste0("r14_", c))]
}
M[, COHORT_FINNGENID := fifelse(!is.na(COHORT), COHORT, r14_COHORT)]
M[, paste0("r14_", c("COHORT","SEX","BL_AGE","BL_YEAR","BMI","SMOKE3","APPROX_BIRTH_DATE")) := NULL]

# age at collection, guarding the one corrupt source collection year
M[, cy := suppressWarnings(as.integer(substr(collection_date, 1, 4)))]
M[, age_at_collection := fifelse(!is.na(cy) & cy >= 1990 & cy <= 2026,
                                 as.numeric(BL_AGE) + (cy - as.numeric(BL_YEAR)), NA_real_)]
M[, collection_year_implausible := !is.na(cy) & (cy < 1990 | cy > 2026)]

# ---- 8. longitudinal structure --------------------------------------------------------------------
setorder(M, FINNGENID, collection_date, SAMPLE_ID)
M[, n_samples_this_individual := .N, by = FINNGENID]
M[, sample_index_within_individual := seq_len(.N), by = FINNGENID]
M[, days_since_previous_sample :=
    as.integer(as.Date(collection_date) - shift(as.Date(collection_date))), by = FINNGENID]

# ---- 9. QC flags from step 05d --------------------------------------------------------------------
if (!is.null(opt$qc_metadata) && file.exists(opt$qc_metadata)) {
  q <- fread(opt$qc_metadata, na.strings = NAS)
  idc <- intersect(c("SampleID","SAMPLE_ID"), names(q))[1]
  keep <- c(idc, grep("^QC_", names(q), value = TRUE))
  M <- merge(M, unique(q[, ..keep], by = idc), by.x = "SAMPLE_ID", by.y = idc, all.x = TRUE)
  say("QC flags merged: ", length(keep) - 1, " columns")
}

# ---- 10. tidy and write ---------------------------------------------------------------------------
M[, c("base_sample_id","cy") := NULL]
setcolorder(M, c("SAMPLE_ID","FINNGENID","Sample_set","Sample_set_source","is_task_force_sample",
                 "stratum","BIOBANK_PLASMA_THL" <- "BIOBANK","COHORT_FINNGENID","PlateID","WellID"))
setorder(M, SAMPLE_ID)
out <- file.path(opt$outdir, sprintf("FG3_batch03_delivery_metadata_%s.tsv", opt$batch))
fwrite(M, out, sep = "\t")
write_parquet(M, sub("\\.tsv$", ".parquet", out))
say("wrote ", out, " : ", nrow(M), " rows x ", ncol(M), " cols")

say("=== task force assignment ===")
print(M[, .N, by = .(Sample_set, Sample_set_source)][order(-N)])
say("=== assignment source summary ===")
print(M[, .(samples = .N, individuals = uniqueN(FINNGENID)), by = Sample_set_source][order(-samples)])
