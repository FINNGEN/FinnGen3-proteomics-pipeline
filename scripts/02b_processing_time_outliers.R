#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# 02b — processing-time outlier detection, re-run on the recovered dates
#       (closes the gap left by D-37)
#
# The processing-time sub-filter of step 02 never executed. The metadata the run
# was given had APPROX_TIMESTAMP_PROCESSING for 149 of 6,191 samples and
# APPROX_TIMESTAMP_FREEZING for none, so detect_processing_outliers() returned at
# its "No valid processing times found" guard. The published zero was therefore
# not a result.
#
# The provider selection record, recovered after the run (D-32), does carry
# freezing dates and an explicit HOURS_FROM_COLLECTION_TO_FREEZING. This step
# re-runs the filter against them, reproducing step 02's algorithm exactly
# (02_technical_outliers.R lines 196-233):
#
#   processing_hours = difftime(freezing, collection, units = "hours")
#   keep only  !is.na & > 0 & < 48          <- an a-priori validity window
#   center     = median(processing_hours)
#   mad        = mad(processing_hours, constant = 1.4826)
#   threshold  = 5 * mad
#   outlier    if |hours - center| > threshold
#
# Two things are checked that the original could not check, because both need the
# recovered field:
#   1. the provider's own HOURS_FROM_COLLECTION_TO_FREEZING is cross-validated
#      against hours recomputed from the freezing and collection timestamps. If
#      those disagree, one of them is wrong and neither should drive a QC verdict.
#   2. the effect of the hard (0, 48) window is reported explicitly, because it
#      discards rather than flags: a sample frozen 60 hours after collection is
#      dropped from the filter's view instead of being called an outlier, which
#      is the opposite of what an analyst would expect from a QC step.
#
# Read-only. Writes new artefacts; changes no existing QC verdict.
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({ library(data.table) })

OUT <- "/mnt/longGWAS_disk_100GB/long_gwas/Github_clones/fg3_olink_pipeline/output_batch03_phase1_qc"
REL <- "/mnt/longGWAS_disk_100GB/long_gwas/11.fg3_Proteomics/08.genewiz_batch3/01.QCed_Batch03_Release_Aug2026"
O <- function(...) file.path(OUT, ...)
say <- function(...) cat(sprintf(...), "\n", sep = "")
nz  <- function(x) !is.na(x) & trimws(as.character(x)) != ""

d <- fread(file.path(REL, "FG3_batch03_delivery_metadata_fg3_batch_03.tsv"),
           colClasses = list(character = c("SAMPLE_ID", "FINNGENID")))
say("== inputs ==")
say("   delivered samples ......................... %d", nrow(d))

# ---- what the run actually had, for the record -----------------------------
runmeta <- fread(O("qc/batch_03/00b_metadata_fg3_batch_03.tsv"),
                 colClasses = list(character = "SAMPLE_ID"))
say("   AS RUN: APPROX_TIMESTAMP_PROCESSING ....... %d of %d (%.1f%%)",
    sum(nz(runmeta$APPROX_TIMESTAMP_PROCESSING)), nrow(runmeta),
    100 * mean(nz(runmeta$APPROX_TIMESTAMP_PROCESSING)))
say("   AS RUN: APPROX_TIMESTAMP_FREEZING ......... %d of %d (%.1f%%)  <- why it never ran",
    sum(nz(runmeta$APPROX_TIMESTAMP_FREEZING)), nrow(runmeta),
    100 * mean(nz(runmeta$APPROX_TIMESTAMP_FREEZING)))

# ---- the recovered fields --------------------------------------------------
ts <- function(dt, tm) {
  dt <- trimws(as.character(dt)); tm <- trimws(as.character(tm))
  out <- ifelse(nz(dt), ifelse(nz(tm), paste(dt, tm), paste(dt, "00:00:00")), NA_character_)
  suppressWarnings(as.POSIXct(out, tz = "UTC"))
}
d[, coll_t := ts(collection_date, TIME_COLLECTION)]
d[, freez_t := ts(APPROX_DATE_FREEZING, TIME_FREEZING)]
d[, hours_recomputed := as.numeric(difftime(freez_t, coll_t, units = "hours"))]
d[, hours_provider := suppressWarnings(as.numeric(HOURS_FROM_COLLECTION_TO_FREEZING))]

say("\n== recovered coverage ==")
for (cn in c("collection_date", "APPROX_DATE_FREEZING", "TIME_FREEZING",
             "HOURS_FROM_COLLECTION_TO_FREEZING", "APPROX_DATE_PROCESSING")) {
  say("   %-34s %4d of %d (%5.1f%%)", cn, sum(nz(d[[cn]])), nrow(d), 100 * mean(nz(d[[cn]])))
}
say("   hours recomputable from timestamps  %4d of %d (%5.1f%%)",
    sum(!is.na(d$hours_recomputed)), nrow(d), 100 * mean(!is.na(d$hours_recomputed)))

# ---- (1) cross-validate the provider's hours against the timestamps --------
say("\n== cross-validating the provider's HOURS field against recomputed hours ==")
both <- d[!is.na(hours_provider) & !is.na(hours_recomputed)]
say("   both available ............................ %d", nrow(both))
if (nrow(both)) {
  dd <- both$hours_provider - both$hours_recomputed
  say("   agree to within 1 h ...................... %d (%.1f%%)", sum(abs(dd) < 1), 100 * mean(abs(dd) < 1))
  say("   agree to within 24 h ..................... %d (%.1f%%)", sum(abs(dd) < 24), 100 * mean(abs(dd) < 24))
  say("   median |difference| ...................... %.3f h ; max %.1f h",
      median(abs(dd)), max(abs(dd)))
  say("   Pearson r ................................ %.6f", cor(both$hours_provider, both$hours_recomputed))
}

# ---- implausible values ----------------------------------------------------
h <- d$hours_provider
say("\n== plausibility of the provider hours field (n = %d) ==", sum(!is.na(h)))
say("   negative (freezing BEFORE collection) .... %d  (min %.2f h)", sum(h < 0, na.rm = TRUE),
    suppressWarnings(min(h, na.rm = TRUE)))
say("   zero ..................................... %d", sum(h == 0, na.rm = TRUE))
say("   in the filter's (0, 48) window ........... %d", sum(h > 0 & h < 48, na.rm = TRUE))
say("   48-168 h ................................. %d", sum(h >= 48 & h <= 168, na.rm = TRUE))
say("   over 168 h (> 1 week) .................... %d  (max %.1f h)", sum(h > 168, na.rm = TRUE),
    suppressWarnings(max(h, na.rm = TRUE)))

# ---- (2) run the filter, exactly as step 02 defines it ---------------------
run_filter <- function(hours, ids, label, lo = 0, hi = 48) {
  m <- data.table(SAMPLE_ID = ids, processing_hours = hours)
  m <- m[!is.na(processing_hours) & processing_hours > lo & processing_hours < hi]
  say("\n== %s ==", label)
  if (!nrow(m)) { say("   no valid processing times -> filter inoperable"); return(NULL) }
  ctr <- median(m$processing_hours, na.rm = TRUE)
  md  <- mad(m$processing_hours, constant = 1.4826, na.rm = TRUE)
  thr <- 5 * md
  m[, processing_outlier := abs(processing_hours - ctr) > thr]
  say("   samples with valid processing time ....... %d (%.1f%% of 6,191)", nrow(m), 100 * nrow(m) / 6191)
  say("   center (median) .......................... %.2f h", ctr)
  say("   MAD ...................................... %.2f h", md)
  say("   threshold ................................ +/- %.2f h  (5 x MAD)", thr)
  say("   window ................................... (%g, %g) h", lo, hi)
  say("   OUTLIERS ................................. %d  (%.3f%% of those evaluated)",
      sum(m$processing_outlier), 100 * mean(m$processing_outlier))
  if (sum(m$processing_outlier))
    say("   flagged range ............................ %.2f - %.2f h",
        min(m[processing_outlier == TRUE, processing_hours]), max(m[processing_outlier == TRUE, processing_hours]))
  m[]
}

res <- run_filter(d$hours_provider, d$SAMPLE_ID, "FILTER AS SPECIFIED — provider hours, (0,48) window")

# sensitivity: the window is a scientific choice, so show what it hides
say("\n== sensitivity to the (0, 48) validity window ==")
say("   The window DISCARDS rather than flags. A sample frozen long after collection")
say("   leaves the filter's view entirely, so widening the window can only add flags.")
for (hi in c(48, 72, 168, Inf)) {
  m <- data.table(h = d$hours_provider)[!is.na(h) & h > 0 & h < hi]
  if (!nrow(m)) next
  ctr <- median(m$h); md <- mad(m$h, constant = 1.4826); thr <- 5 * md
  say("   window (0,%-4s): evaluated %4d | center %6.2f | MAD %5.2f | thr %6.2f | outliers %3d",
      ifelse(is.infinite(hi), "Inf", hi), nrow(m), ctr, md, thr, sum(abs(m$h - ctr) > thr))
}

# ---- who would be flagged, and are they in the release? --------------------
if (!is.null(res) && sum(res$processing_outlier)) {
  passed <- rownames(readRDS(file.path(REL, "npx_matrix_5928_qc_passed_fg3_batch_03.rds")))
  fl <- res[processing_outlier == TRUE]
  fl[, Sample_set := d$Sample_set_inferred[match(SAMPLE_ID, d$SAMPLE_ID)]]
  fl[, in_release := SAMPLE_ID %in% passed]
  fl[, already_flagged := d$QC_flag[match(SAMPLE_ID, d$SAMPLE_ID)]]
  say("\n== the processing-time outliers ==")
  print(as.data.frame(fl[order(-processing_hours)][, .(SAMPLE_ID, processing_hours = round(processing_hours, 2),
        Sample_set, in_release, already_flagged)]), row.names = FALSE)
  say("   of %d flagged: %d already excluded by another method, %d STILL IN THE RELEASE",
      nrow(fl), sum(fl$already_flagged == 1), sum(fl$in_release))
  fwrite(fl, O("outliers/02b_processing_time_outliers.tsv"), sep = "\t")
}
if (!is.null(res)) fwrite(res, O("outliers/02b_processing_time_stats.tsv"), sep = "\t")
saveRDS(res, O("outliers/02b_processing_time.rds"))
say("\n== wrote 02b_processing_time_stats.tsv / 02b_processing_time_outliers.tsv ==")
