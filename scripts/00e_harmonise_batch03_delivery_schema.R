#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# 00e — Harmonise the Batch 03 delivery package with the Batch 02 schema, and
#       repair the defects found by the post-release adversarial audit.
#
# Batch 02 precedent (the schema this must be compatible with):
#   11.fg3_Proteomics/02.genewiz_batch2/02.QCed_Batch02_Release_Dec2025/
#     qc_annotated_metadata_all_2527_samples_fg3_batch_02.tsv   (56 cols)
#     86_outliers_list_fg3_batch_02.tsv                         (32 cols)
#   11.fg3_Proteomics/02.genewiz_batch2/00.metadata/
#     FG3_Batch_2_Olink_Proteomics_Metadata.tsv                 (26 cols)
#
# What this script fixes:
#   A. DISEASE_GROUP in the comprehensive outliers list is NA for all 263 rows
#      -> populate from the authoritative Sample_set.
#   B. qc_annotated_metadata carries 54 "_b0"-suffixed duplicate columns. 51 are
#      value-identical dead weight; 3 (FINNGENID, COHORT_FINNGENID,
#      BIOBANK_PLASMA) have a BLANK base value for the 5 initial-QC failures
#      while the _b0 copy holds the true value -> coalesce base <- _b0, then drop
#      every _b0 column. 140 -> 86 information-bearing columns.
#   C. The four-way technical decomposition (plate / temporal batch / processing
#      time / sample-level 5xMAD) exists only in un-delivered pipeline RDS files,
#      so an analyst cannot reproduce the per-method removal table -> emit the
#      four sub-flags into every metadata artefact.
#   D. Batch 02 shipped ONE-HOT task-force indicator columns; Batch 03 shipped a
#      single Sample_set factor plus registry-derived elig_* flags, which are a
#      different variable -> emit Batch-02-style one-hot indicators alongside
#      Sample_set, using the Batch 02 vocabulary wherever the concept matches.
#   E. Per-level EVALUABLE denominators are not recoverable from the delivered
#      files, so flagged/delivered rates are misleading for filters that could
#      not assess a given sample -> emit QC_*_evaluable columns.
#
# Read-only with respect to the QC results themselves: no sample's QC verdict
# changes, no matrix is rewritten. Only metadata artefacts are re-emitted.
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(arrow)
})

REL  <- "/mnt/longGWAS_disk_100GB/long_gwas/11.fg3_Proteomics/08.genewiz_batch3/01.QCed_Batch03_Release_Aug2026"
OUT  <- "/mnt/longGWAS_disk_100GB/long_gwas/Github_clones/fg3_olink_pipeline/output_batch03_phase1_qc"
B02  <- "/mnt/longGWAS_disk_100GB/long_gwas/11.fg3_Proteomics/02.genewiz_batch2"
SUF  <- "fg3_batch_03"
R <- function(...) file.path(REL, ...)
O <- function(...) file.path(OUT, ...)
# Pristine pre-00e copies of the three artefacts this script re-emits, so the
# script is re-runnable: it always READS the untouched originals and WRITES the
# harmonised versions. Created by hand before the first run; recreated here if
# absent so a fresh checkout still behaves.
PRIS <- file.path(dirname(REL), "05.audit_evidence_pre_00e")
Rin <- function(fn) { p <- file.path(PRIS, fn); if (file.exists(p)) p else file.path(REL, fn) }

say <- function(...) cat(sprintf(...), "\n", sep = "")
nz  <- function(x) !is.na(x) & trimws(as.character(x)) != ""

# --------------------------------------------------------------------------
# 1. Authoritative inputs
# --------------------------------------------------------------------------
say("== loading authoritative inputs ==")
dm <- fread(Rin(sprintf("FG3_batch03_delivery_metadata_%s.tsv", SUF)),
            colClasses = list(character = c("SAMPLE_ID", "FINNGENID")))
stopifnot(nrow(dm) == 6191, !anyDuplicated(dm$SAMPLE_ID))

passed   <- rownames(readRDS(R(sprintf("npx_matrix_5928_qc_passed_%s.rds", SUF))))
rm_truth <- setdiff(dm$SAMPLE_ID, passed)
stopifnot(length(passed) == 5928, length(rm_truth) == 263)

# four-way technical decomposition, from the step-02 artefacts
plt  <- readRDS(O("outliers/02_plate_outliers.rds"))
bat  <- readRDS(O("outliers/02_batch_outliers.rds"))
prc  <- readRDS(O("outliers/02_processing_outliers.rds"))
stec <- readRDS(O("outliers/02_sample_technical_outliers.rds"))
tcmb <- readRDS(O("outliers/02_technical_outliers_combined.rds"))
pca  <- readRDS(O("outliers/01_pca_outliers_original.rds"))
sxp  <- fread(O("outliers/04_sex_predictions.tsv"), colClasses = list(character = "SAMPLE_ID"))
pqs  <- fread(O("outliers/05b_05b_pqtl_stats.tsv"), colClasses = list(character = "FINNGENID"))

# The four QC_pqtl_* metric columns shipped 0-populated in all three metadata
# artefacts although the release note documents them (audit finding H1). The
# statistics exist per FINNGENID in 05b_pqtl_stats.tsv; join them back.
pq_cols <- intersect(c("MeanAbsZ", "MedianAbsZ", "MaxAbsZ", "MedianAbsResidual", "N_Prots", "N_ValidResiduals"),
                     names(pqs))
say("   pQTL stat columns available: %s", paste(pq_cols, collapse = ", "))
pq_map <- unique(pqs, by = "FINNGENID")

keep <- function(v) intersect(unique(as.character(v[nz(v)])), dm$SAMPLE_ID)
tech <- list(plate            = keep(plt$samples_from_outlier_plates),
             temporal_batch   = keep(bat$samples_from_outlier_batches),
             processing       = keep(prc$processing_outliers),
             sample_level     = keep(stec$technical_outliers))
say("   technical sub-flags: plate %d | temporal %d | processing %d | sample-level %d",
    length(tech$plate), length(tech$temporal_batch), length(tech$processing), length(tech$sample_level))
stopifnot(setequal(keep(tcmb$all_outliers), Reduce(union, tech)))
say("   union == step-02 combined (%d) OK", length(keep(tcmb$all_outliers)))

# --------------------------------------------------------------------------
# 2. Task-force one-hot indicators, Batch 02 style
# --------------------------------------------------------------------------
# Batch 02 vocabulary, in the order Batch 02 shipped it.
B02_GROUPS <- c("Kidney", "Kids", "F64", "MFGE8", "Parkinsons", "Metabolic", "AMD",
                "Rheuma", "Pulmo", "Chromosomal_Abnormalities", "Blood_donors",
                "Bridging_samples")
# Batch 03 Sample_set value -> indicator column name.
#   Five names are shared with Batch 02 verbatim. "Bridging" takes the Batch 02
#   name "Bridging_samples". The rest are Batch-03-only task forces.
#   "Eye" is deliberately NOT aliased to Batch 02's "AMD": AMD is age-related
#   macular degeneration specifically, a narrower concept than the Batch 03
#   "Eye diseases" task force. Conflating them would silently widen AMD.
SET2COL <- c("Kidney" = "Kidney", "Parkinsons" = "Parkinsons", "Metabolic" = "Metabolic",
             "Rheuma" = "Rheuma", "Pulmo" = "Pulmo", "Bridging" = "Bridging_samples",
             "IBD" = "IBD", "Eye" = "Eye", "Heart failure" = "Heart_failure",
             "Atopic_dermatitis" = "Atopic_dermatitis", "Psychosis" = "Psychosis",
             "NFBC1966" = "NFBC1966")
# Batch 02 groups with no Batch 03 counterpart:
#   0  = positive documented evidence of absence from Batch 03
#   NA = genuinely unknown for Batch 03; must not be asserted as 0
B02_ABSENT_ZERO <- c("F64", "Chromosomal_Abnormalities")
B02_ABSENT_NA   <- c("Kids", "MFGE8", "AMD", "Blood_donors")

ss <- dm[, .(SAMPLE_ID, Sample_set = fifelse(nz(Sample_set), Sample_set, NA_character_),
             Sample_set_source)]
stopifnot(all(ss[nz(Sample_set), Sample_set] %in% names(SET2COL)))

onehot <- data.table(SAMPLE_ID = ss$SAMPLE_ID)
for (s in names(SET2COL)) onehot[[SET2COL[[s]]]] <- as.integer(ss$Sample_set == s & !is.na(ss$Sample_set))
for (g in B02_ABSENT_ZERO) onehot[[g]] <- 0L
for (g in B02_ABSENT_NA)   onehot[[g]] <- NA_integer_
# Batch 02 column order first, then the Batch-03-only groups
oh_order <- c(B02_GROUPS, setdiff(unname(SET2COL), B02_GROUPS))
setcolorder(onehot, c("SAMPLE_ID", oh_order))
say("== one-hot indicators: %d columns (%d B02 vocabulary + %d B03-only) ==",
    ncol(onehot) - 1L, length(B02_GROUPS), length(setdiff(unname(SET2COL), B02_GROUPS)))
chk <- rowSums(as.matrix(onehot[, ..oh_order]), na.rm = TRUE)
say("   rows with exactly one indicator set: %d ; zero: %d (the %d unassigned)",
    sum(chk == 1), sum(chk == 0), ss[is.na(Sample_set), .N])
stopifnot(sum(chk > 1) == 0, sum(chk == 0) == ss[is.na(Sample_set), .N])

# --------------------------------------------------------------------------
# 3. Per-level evaluability
# --------------------------------------------------------------------------
iq_fail <- dm[QC_initial_qc %in% c(TRUE, 1), SAMPLE_ID]
ar      <- setdiff(dm$SAMPLE_ID, iq_fail)                    # 6,186 analysis-ready
ev <- data.table(SAMPLE_ID = dm$SAMPLE_ID)
ev[, QC_pca_evaluable        := as.integer(SAMPLE_ID %in% ar)]
ev[, QC_technical_evaluable_plate      := as.integer(SAMPLE_ID %in% ar & nz(dm$PlateID))]
ev[, QC_technical_evaluable_temporal   := as.integer(SAMPLE_ID %in% ar & nz(dm$collection_date))]
# The processing-time filter was SKIPPED, not passed. The metadata the run was
# given (output/qc/batch_03/00b_metadata) had APPROX_TIMESTAMP_PROCESSING for
# 149 of 6,191 samples (2.4%) and APPROX_TIMESTAMP_FREEZING for 0, and
# 02_technical_outliers.log records "No valid processing times found. Skipping
# processing time outlier detection." So nothing was evaluable AS RUN.
ev[, QC_technical_evaluable_processing := 0L]
# Forward-looking: the provider selection record, recovered after the run,
# does carry processing dates. This column says what a RE-RUN could evaluate.
ev[, QC_processing_date_available := as.integer(SAMPLE_ID %in% ar & nz(dm$APPROX_DATE_PROCESSING))]
ev[, QC_technical_evaluable_sample     := as.integer(SAMPLE_ID %in% ar)]
ev[, QC_zscore_evaluable     := as.integer(SAMPLE_ID %in% ar)]
# step 04 reads the PCA-CLEANED matrix, so sex is evaluated only on PCA survivors
ev[, QC_sex_evaluable        := as.integer(SAMPLE_ID %in% intersect(ar, sxp[nz(genetic_sex), SAMPLE_ID]))]
# step 05b keys on FINNGENID: one aliquot per individual is assessed, the repeats are not
ev[, QC_pqtl_evaluable       := as.integer(SAMPLE_ID %in% ar & dm$FINNGENID %in% pqs$FINNGENID)]
# NOTE: step 05b produces one statistic per FINNGENID and does not record which
# aliquot it was computed from, so which specific tube was assessed is unknowable
# from the artefacts. QC_pqtl_evaluable therefore means "this sample's individual
# has a pQTL statistic", NOT "this tube was checked". For a repeat-sampled
# individual only one tube's worth of evidence exists; see defect D-28.
ev[, QC_pqtl_n_aliquots_this_individual :=
     dm$n_samples_this_individual[match(SAMPLE_ID, dm$SAMPLE_ID)]]
say("== evaluability ==")
for (cn in setdiff(names(ev), "SAMPLE_ID")) {
  v <- ev[[cn]]
  if (all(v %in% c(0L, 1L, NA_integer_)))                       # 0/1 evaluability flag
    say("   %-38s %5d of 6191 (%5.1f%%)", cn, sum(v, na.rm = TRUE), 100 * sum(v, na.rm = TRUE) / 6191)
  else                                                          # a count column, not a flag
    say("   %-38s median %g, range %g-%g", cn, median(v, na.rm = TRUE),
        min(v, na.rm = TRUE), max(v, na.rm = TRUE))
}

fill_pqtl <- function(dt, fgid_col = "FINNGENID") {
  m <- match(dt[[fgid_col]], pq_map$FINNGENID)
  pick <- function(cands) { h <- intersect(cands, names(pq_map)); if (length(h)) pq_map[[h[1]]][m] else NA }
  dt[, QC_pqtl_mean_abs_z          := pick(c("MeanAbsZ", "mean_abs_z"))]
  dt[, QC_pqtl_max_abs_z           := pick(c("MaxAbsZ", "max_abs_z"))]
  dt[, QC_pqtl_median_abs_residual := pick(c("MedianAbsResidual", "median_abs_residual"))]
  dt[, QC_pqtl_n_prots             := pick(c("N_Prots", "n_prots"))]
  dt
}

# --------------------------------------------------------------------------
# 4. Fix A + C: comprehensive outliers list
# --------------------------------------------------------------------------
say("== fixing the comprehensive outliers list ==")
ol <- fread(Rin(sprintf("comprehensive_outliers_list_%s.tsv", SUF)),
            colClasses = list(character = c("SampleID", "FINNGENID")))
stopifnot(nrow(ol) == 263)
say("   before: DISEASE_GROUP non-empty = %d of %d", sum(nz(ol$DISEASE_GROUP)), nrow(ol))

# (H8) The same blanked-identity defect as D-34: FINNGENID and BIOBANK_PLASMA are
# empty for the 5 initial-QC failures in the outliers list too. Repair from the
# authoritative delivery metadata.
for (cn in c("FINNGENID", "BIOBANK_PLASMA")) {
  src <- if (cn == "BIOBANK_PLASMA") dm$BIOBANK else dm[[cn]]
  bad <- !nz(ol[[cn]])
  if (any(bad)) {
    ol[[cn]][bad] <- as.character(src)[match(ol$SampleID, dm$SAMPLE_ID)][bad]
    say("   repaired blank %s on %d outlier rows", cn, sum(bad))
  }
}
ol <- fill_pqtl(ol)
ol[, DISEASE_GROUP := ss$Sample_set[match(SampleID, ss$SAMPLE_ID)]]
ol[, Sample_set := DISEASE_GROUP]
ol[, Sample_set_source := ss$Sample_set_source[match(SampleID, ss$SAMPLE_ID)]]
for (k in names(tech)) ol[[paste0("QC_technical_", k)]] <- as.integer(ol$SampleID %in% tech[[k]])
ol <- merge(ol, ev, by.x = "SampleID", by.y = "SAMPLE_ID", all.x = TRUE, sort = FALSE)
say("   after:  DISEASE_GROUP non-empty = %d of %d ; %d distinct task forces",
    sum(nz(ol$DISEASE_GROUP)), nrow(ol), uniqueN(ol$DISEASE_GROUP[nz(ol$DISEASE_GROUP)]))
say("   technical sub-flags in list: plate %d temporal %d processing %d sample %d",
    sum(ol$QC_technical_plate), sum(ol$QC_technical_temporal_batch),
    sum(ol$QC_technical_processing), sum(ol$QC_technical_sample_level))

# --------------------------------------------------------------------------
# 5. Fix B + C + D: qc_annotated_metadata
# --------------------------------------------------------------------------
say("== fixing qc_annotated_metadata ==")
qa <- fread(Rin(sprintf("qc_annotated_metadata_all_6191_samples_%s.tsv", SUF)),
            colClasses = list(character = c("SAMPLE_ID", "FINNGENID", "FINNGENID_b0")))
stopifnot(nrow(qa) == 6191)
b0 <- grep("_b0$", names(qa), value = TRUE)
say("   before: %d columns, %d of them _b0 duplicates", ncol(qa), length(b0))

coalesced <- character(0)
for (col_b0 in b0) {
  base <- sub("_b0$", "", col_b0)
  if (!base %in% names(qa)) next
  bad <- !nz(qa[[base]]) & nz(qa[[col_b0]])
  if (any(bad)) {
    qa[[base]][bad] <- as.character(qa[[col_b0]])[bad]
    coalesced <- c(coalesced, sprintf("%s (%d rows)", base, sum(bad)))
  }
}
say("   coalesced base <- _b0 where base was blank: %s",
    if (length(coalesced)) paste(coalesced, collapse = ", ") else "none")
qa[, (b0) := NULL]

# The registry-derived pair is NOT the task force. Leaving it named `Disease_Group`
# next to a task-force column named `DISEASE_GROUP` — differing only in case — is a
# trap. Rename it explicitly. It stays in the file; only the name changes.
if ("Disease_Group" %in% names(qa))        setnames(qa, "Disease_Group", "Disease_Group_registry_derived")
if ("Disease_Group_source" %in% names(qa)) setnames(qa, "Disease_Group_source", "Disease_Group_registry_source")
qa <- fill_pqtl(qa)
qa[, Sample_set        := ss$Sample_set[match(SAMPLE_ID, ss$SAMPLE_ID)]]
qa[, Sample_set_source := ss$Sample_set_source[match(SAMPLE_ID, ss$SAMPLE_ID)]]
qa[, DISEASE_GROUP     := Sample_set]
for (k in names(tech)) qa[[paste0("QC_technical_", k)]] <- as.integer(qa$SAMPLE_ID %in% tech[[k]])
qa <- merge(qa, onehot, by = "SAMPLE_ID", all.x = TRUE, sort = FALSE)
qa <- merge(qa, ev,     by = "SAMPLE_ID", all.x = TRUE, sort = FALSE)
say("   after: %d columns", ncol(qa))
stopifnot(nrow(qa) == 6191, sum(!nz(qa$FINNGENID)) == 0)
say("   blank FINNGENID rows now: %d (was 5)", sum(!nz(qa$FINNGENID)))

# --------------------------------------------------------------------------
# 6. Fix C + D + E: delivery metadata
# --------------------------------------------------------------------------
say("== extending the delivery metadata ==")
dmo <- copy(dm)
# Batch 02 column-name aliases the Batch 03 delivery table was missing
dmo[, CONTAINER_NAME := CONTAINER]
dmo[, BIOBANK_PLASMA := BIOBANK]
ts <- function(d, t) {
  d <- trimws(as.character(d)); t <- trimws(as.character(t))
  out <- ifelse(nz(d), ifelse(nz(t), paste(d, t), d), NA_character_)
  out
}
dmo[, APPROX_TIMESTAMP_PROCESSING := ts(APPROX_DATE_PROCESSING, TIME_PROCESSING)]
dmo[, APPROX_TIMESTAMP_FREEZING   := ts(APPROX_DATE_FREEZING,   TIME_FREEZING)]
for (k in names(tech)) dmo[[paste0("QC_technical_", k)]] <- as.integer(dmo$SAMPLE_ID %in% tech[[k]])
dmo[, N_Methods := ol$N_Methods[match(SAMPLE_ID, ol$SampleID)]]
dmo[is.na(N_Methods), N_Methods := 0L]
dmo[, Detection_Steps := ol$Detection_Steps[match(SAMPLE_ID, ol$SampleID)]]
dmo <- merge(dmo, onehot, by = "SAMPLE_ID", all.x = TRUE, sort = FALSE)
dmo <- merge(dmo, ev,     by = "SAMPLE_ID", all.x = TRUE, sort = FALSE)
dmo <- fill_pqtl(dmo)
dmo[, in_qc_passed_matrix := as.integer(SAMPLE_ID %in% passed)]
# (M9) The recovered provider field has impossible values: 13 negative and 16 over
# a week. Flag them rather than silently shipping them as usable.
hh <- suppressWarnings(as.numeric(dmo$HOURS_FROM_COLLECTION_TO_FREEZING))
dmo[, hours_to_freezing_implausible := as.integer(!is.na(hh) & (hh < 0 | hh > 168))]
say("   flagged implausible hours-to-freezing: %d (negative %d, >168h %d)",
    sum(dmo$hours_to_freezing_implausible), sum(hh < 0, na.rm = TRUE), sum(hh > 168, na.rm = TRUE))
say("   %d -> %d columns", ncol(dm), ncol(dmo))
stopifnot(nrow(dmo) == 6191,
          sum(dmo$in_qc_passed_matrix) == 5928,
          sum(dmo$in_qc_passed_matrix == 0) == 263)
say("   in_qc_passed_matrix: %d passed / %d removed OK",
    sum(dmo$in_qc_passed_matrix), sum(dmo$in_qc_passed_matrix == 0))

# --------------------------------------------------------------------------
# 6b. The sample metadata schema file still carried the registry-derived
#     Disease_Group (populated for only 1,305 of 6,191) and no task force at
#     all, so an analyst reading it would take registry eligibility for the
#     disease group. Add the task force and rename the registry pair.
# --------------------------------------------------------------------------
say("== extending the sample metadata schema ==")
sch_name <- sprintf("FG3_batch03_sample_metadata_schema_%s.tsv", SUF)
sch <- fread(Rin(sch_name), colClasses = list(character = c("SAMPLE_ID", "FINNGENID")))
stopifnot(nrow(sch) == 6191)
say("   before: %d columns ; Disease_Group populated %d of %d",
    ncol(sch), sum(nz(sch$Disease_Group)), nrow(sch))
if ("Disease_Group" %in% names(sch))        setnames(sch, "Disease_Group", "Disease_Group_registry_derived")
if ("Disease_Group_source" %in% names(sch)) setnames(sch, "Disease_Group_source", "Disease_Group_registry_source")
sch[, Sample_set        := ss$Sample_set[match(SAMPLE_ID, ss$SAMPLE_ID)]]
sch[, Sample_set_source := ss$Sample_set_source[match(SAMPLE_ID, ss$SAMPLE_ID)]]
sch[, DISEASE_GROUP     := Sample_set]
sch <- merge(sch, onehot, by = "SAMPLE_ID", all.x = TRUE, sort = FALSE)
say("   after:  %d columns ; Sample_set populated %d of %d",
    ncol(sch), sum(nz(sch$Sample_set)), nrow(sch))

# --------------------------------------------------------------------------
# 7. Write. Every artefact keeps its filename; .tsv and .parquet stay in step.
# --------------------------------------------------------------------------
# (C2) The TSV convention is an empty field for missing, matching Batch 02. But the
# parquet twin must carry a TRUE NULL, not an empty string, or is.na() silently
# fails on it. Values read in from the pristine TSVs arrive as "", so normalise
# every character column before writing the parquet.
to_null <- function(dt) {
  out <- copy(dt)
  for (cn in names(out)) if (is.character(out[[cn]])) {
    v <- out[[cn]]; v[!is.na(v) & trimws(v) == ""] <- NA_character_
    out[[cn]] <- v
  }
  out
}
wr <- function(dt, stem) {
  f_tsv <- R(sprintf("%s_%s.tsv", stem, SUF)); f_pq <- R(sprintf("%s_%s.parquet", stem, SUF))
  fwrite(dt, f_tsv, sep = "\t", na = "", quote = FALSE)
  write_parquet(as.data.frame(to_null(dt)), f_pq)
  say("   wrote %s  (%d x %d)", basename(f_tsv), nrow(dt), ncol(dt))
  say("   wrote %s", basename(f_pq))
}
say("== writing ==")
wr(ol,  "comprehensive_outliers_list")
wr(qa,  "qc_annotated_metadata_all_6191_samples")
wr(dmo, "FG3_batch03_delivery_metadata")
f_sch <- R(sch_name)
fwrite(sch, f_sch, sep = "\t", na = "", quote = FALSE)
say("   wrote %s  (%d x %d)", basename(f_sch), nrow(sch), ncol(sch))

# --------------------------------------------------------------------------
# 8. Post-write verification
# --------------------------------------------------------------------------
say("== verification ==")
for (stem in c("comprehensive_outliers_list", "qc_annotated_metadata_all_6191_samples",
               "FG3_batch03_delivery_metadata")) {
  a <- fread(R(sprintf("%s_%s.tsv", stem, SUF)), colClasses = "character")
  b <- as.data.table(read_parquet(R(sprintf("%s_%s.parquet", stem, SUF))))
  say("   %-42s tsv %d x %d | parquet %d x %d | cols match %s",
      stem, nrow(a), ncol(a), nrow(b), ncol(b), identical(names(a), names(b)))
}
b2g <- fread(file.path(B02, "02.QCed_Batch02_Release_Dec2025",
                       "qc_annotated_metadata_all_2527_samples_fg3_batch_02.tsv"), nrows = 0)
missing_b02 <- setdiff(names(b2g), names(fread(R(sprintf("qc_annotated_metadata_all_6191_samples_%s.tsv", SUF)), nrows = 0)))
say("   Batch 02 qc_annotated columns still absent from Batch 03: %s",
    if (length(missing_b02)) paste(missing_b02, collapse = ", ") else "NONE")
# no empty-string-as-missing may survive in any parquet
for (stem in c("comprehensive_outliers_list", "qc_annotated_metadata_all_6191_samples",
               "FG3_batch03_delivery_metadata")) {
  b <- as.data.table(read_parquet(R(sprintf("%s_%s.parquet", stem, SUF))))
  ch <- names(b)[vapply(b, is.character, logical(1))]
  bad <- sum(vapply(ch, function(cn) sum(!is.na(b[[cn]]) & trimws(b[[cn]]) == ""), integer(1)))
  say("   %-42s parquet empty-string cells (want 0): %d", stem, bad)
}
# every row of every artefact must carry an identity
for (f in c(sprintf("comprehensive_outliers_list_%s.tsv", SUF),
            sprintf("qc_annotated_metadata_all_6191_samples_%s.tsv", SUF),
            sprintf("FG3_batch03_delivery_metadata_%s.tsv", SUF), sch_name)) {
  x <- fread(R(f), colClasses = "character")
  say("   %-52s blank FINNGENID: %d", f, sum(!nz(x$FINNGENID)))
}
pqn <- fread(R(sprintf("FG3_batch03_delivery_metadata_%s.tsv", SUF)), select = "QC_pqtl_mean_abs_z")
say("   QC_pqtl_mean_abs_z populated: %d of 6191", sum(nz(pqn[[1]])))
# no two columns may differ only by case, in any artefact
for (f in c(sprintf("comprehensive_outliers_list_%s.tsv", SUF),
            sprintf("qc_annotated_metadata_all_6191_samples_%s.tsv", SUF),
            sprintf("FG3_batch03_delivery_metadata_%s.tsv", SUF), sch_name)) {
  n <- names(fread(R(f), nrows = 0))
  dup <- names(which(table(tolower(n)) > 1))
  say("   %-52s case-only column collisions: %s", f,
      if (length(dup)) paste(dup, collapse = ", ") else "none")
}
say("done")
