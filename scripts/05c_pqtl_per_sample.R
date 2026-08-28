#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# 05c — pQTL provenance check, PER SAMPLE  (closes defect D-28)
#
# Step 05b keys the provenance check on FINNGENID. At scripts/05b_pqtl_outliers.R:371
# it overwrites the NPX matrix rownames with FINNGENID:
#
#     rownames(npx_matrix) <- dt_mapping$FINNGENID[idx[valid_idx]]
#
# For a repeat-sampled individual that makes the rownames non-unique, and the
# later `npx_matrix[common_samples, prot]` therefore returns only the FIRST
# matching row per individual. The consequence: of the 6,137 samples entering the
# step, only 4,611 (one aliquot per person) were ever scored, and 1,526 repeat
# measurements — in the one FG3 batch defined by serial sampling — were silently
# skipped. The step's own log line "Samples with pQTL stats: 6137" reports the
# merge total, not the number scored, which is why it went unnoticed.
#
# Genotype is a property of the individual, but NPX is a property of the tube, so
# the statistic is well defined per tube: use the individual's genotype with the
# tube's own protein values. This step does exactly that, reusing 05b's cached
# genotype extraction so nothing is re-derived.
#
# The arithmetic is byte-for-byte the 05b algorithm (05b lines 90-152, 1299-1420):
#   NPX    -> inverse rank normal per protein, over ALL 6,137 rows
#             (05b applies IRN at line 399, i.e. AFTER the rowname overwrite, so the
#             transform sees every tube even though only 4,611 are later scored;
#             ranks average ties, quantiles = (rank-0.5)/n, clamped to 1e-6, then qnorm)
#   Geno   = 2 - dosage                      (PLINK --export A flip)
#   stats  = per-genotype Mean and SD, groups with SD > 0 only
#   Z      = (Protein - Mean[Geno]) / SD[Geno]
#   per unit: MeanAbsZ, MedianAbsZ, MaxAbsZ, MedianAbsResidual, N_Prots
#   cutoff = mean(MeanAbsZ) + k_sd * sd(MeanAbsZ), rounded to 1 dp   (k_sd = 4)
#
# STAGE 1 validates the implementation by reproducing 05b's own MeanAbsZ on the
# tubes 05b actually read, then STAGE 2 extends to all 6,137. Without that
# ordering, a per-sample number that merely looks plausible could not be told
# apart from a reimplementation bug.
#
# Identifying the tubes 05b actually read takes care. Because `m[chr, ]` on
# duplicated rownames returns the FIRST match, the tube scored for each individual
# is that individual's first occurrence in the PCA-cleaned matrix row order --
# NOT the SampleID recorded in 05b's own output. Those disagree for 14 of the 4,611
# rows, and in all 14 the recorded SampleID is a tube PCA had already REMOVED, so
# it cannot be the tube whose NPX was read. 05b's SampleID column was evidently
# back-filled from FINNGENID through the mapping table rather than carried through
# the computation. Recorded as D-48: joining 05b_05b_pqtl_stats.tsv on SampleID
# mis-attributes 14 rows. Validation therefore keys on FINNGENID, which is sound.
#
# Read-only. Writes new artefacts; changes no existing QC verdict.
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
})

OUT <- "/mnt/longGWAS_disk_100GB/long_gwas/Github_clones/fg3_olink_pipeline/output_batch03_phase1_qc"
REL <- "/mnt/longGWAS_disk_100GB/long_gwas/11.fg3_Proteomics/08.genewiz_batch3/01.QCed_Batch03_Release_Aug2026"
K_SD <- 4                      # as used by the run: "MeanAbsZ: 1.6002 (mean + 4*SD)"
O <- function(...) file.path(OUT, ...)

say <- function(...) cat(sprintf(...), "\n", sep = "")
nz  <- function(x) !is.na(x) & trimws(as.character(x)) != ""

# --------------------------------------------------------------------------
# 1. Inputs — all reused from the 05b run
# --------------------------------------------------------------------------
say("== inputs ==")
npx <- readRDS(O("outliers/01_npx_matrix_pca_cleaned.rds"))
say("   NPX (PCA-cleaned, SampleID-indexed) : %d samples x %d proteins", nrow(npx), ncol(npx))
stopifnot(!anyDuplicated(rownames(npx)))

# 05b's inverse rank normalisation, per protein, verbatim (05b lines 105-138).
# Column-wise and therefore independent per protein, so restricting to the proteins
# actually used is exact. Crucially the transform runs over ALL rows of the
# PCA-cleaned matrix, matching 05b's call site at line 399.
irn_col <- function(x) {
  valid_idx <- which(!is.na(x)); n_valid <- length(valid_idx)
  if (n_valid == 0L) return(rep(NA_real_, length(x)))
  ranks <- rank(x[valid_idx], ties.method = "average", na.last = "keep")
  q <- (ranks - 0.5) / n_valid
  eps <- 1e-6
  q <- pmax(eps, pmin(1 - eps, q))
  out <- rep(NA_real_, length(x)); out[valid_idx] <- qnorm(q); out
}

map <- as.data.table(readRDS(file.path(REL, "00_sample_mapping_fg3_batch_03.rds")))
map <- map[nz(FINNGENID), .(SampleID = as.character(SampleID), FINNGENID = as.character(FINNGENID))]
say("   sample -> individual map            : %d tubes, %d individuals", nrow(map), uniqueN(map$FINNGENID))

raw_path <- Sys.glob(O("temp_work/pqtl_cache/batch_03/*/finngen_R14_exported.raw"))
stopifnot(length(raw_path) == 1)
geno <- fread(raw_path[1], colClasses = list(character = c("FID", "IID")))
say("   cached genotypes                    : %d individuals x %d variants", nrow(geno), ncol(geno) - 6L)

vars <- fread(O("pqtl/05b_05b_zscore_variants.tsv"))
say("   variant/protein pairs               : %d", nrow(vars))

prots <- intersect(unique(vars$trait), colnames(npx))
npx <- npx[, prots, drop = FALSE]
rn_keep <- rownames(npx)
npx <- apply(npx, 2, irn_col)
# apply() drops dimnames; 05b restores them explicitly at its lines 142-144
rownames(npx) <- rn_keep
colnames(npx) <- prots
stopifnot(!anyDuplicated(rownames(npx)), nrow(npx) == length(rn_keep))
say("   IRN applied over all %d rows x %d proteins: mean %.4f, sd %.4f (05b: 0 / 0.9999)",
    nrow(npx), ncol(npx), mean(npx, na.rm = TRUE), sd(as.vector(npx), na.rm = TRUE))

orig <- fread(O("outliers/05b_05b_pqtl_stats.tsv"),
              colClasses = list(character = c("SampleID", "FINNGENID")))
say("   05b stats (per individual)          : %d rows", nrow(orig))

geno_cols <- setdiff(names(geno), c("FID", "IID", "PAT", "MAT", "SEX", "PHENOTYPE"))
# .raw column names are "<rsid>_<countedAllele>"; strip the trailing allele to match
col_rsid <- sub("_[ACGT]+$", "", geno_cols)
rsid2col <- setNames(geno_cols, col_rsid)

# --------------------------------------------------------------------------
# 2. The 05b Z computation, parameterised by the unit of analysis
# --------------------------------------------------------------------------
# unit_ids : the tube ids to score. Genotype-group Mean/SD are computed over
#            exactly these rows, which is what 05b does over its own sample set.
compute_z <- function(unit_ids, label) {
  say("== computing Z over %s (%d units) ==", label, length(unit_ids))
  fg <- map$FINNGENID[match(unit_ids, map$SampleID)]
  keep <- !is.na(fg) & fg %in% geno$IID & unit_ids %in% rownames(npx)
  unit_ids <- unit_ids[keep]; fg <- fg[keep]
  gi <- match(fg, geno$IID)
  say("   scorable: %d tubes over %d individuals", length(unit_ids), uniqueN(fg))

  res <- vector("list", nrow(vars)); failed <- 0L
  for (i in seq_len(nrow(vars))) {
    prot <- vars$trait[i]; rs <- vars$rsid[i]
    gcol <- rsid2col[[rs]]
    if (is.null(gcol) || is.na(gcol) || !prot %in% colnames(npx)) { failed <- failed + 1L; next }
    df <- data.table(SampleID = unit_ids,
                     Geno     = 2 - as.numeric(geno[[gcol]][gi]),   # PLINK flip, as 05b
                     Protein  = as.numeric(npx[unit_ids, prot]))
    df <- df[!is.na(Protein) & !is.na(Geno)]
    if (!nrow(df)) { failed <- failed + 1L; next }
    st <- df[, .(Mean = mean(Protein), SD = sd(Protein)), by = Geno][!is.na(SD) & SD > 0]
    if (!nrow(st)) { failed <- failed + 1L; next }
    df <- merge(df, st, by = "Geno", all.x = FALSE)
    if (!nrow(df)) { failed <- failed + 1L; next }
    df[, Z := (Protein - Mean) / SD]

    b <- vars$beta[i]                                              # beta-weighted residual, as 05b
    if (!is.na(b) && is.finite(b) && nrow(df) > 3) {
      pm <- mean(df$Protein, na.rm = TRUE); mg <- mean(df$Geno, na.rm = TRUE)
      df[, Residual := Protein - (pm + b * (Geno - mg))]
      rsd <- sd(df$Residual, na.rm = TRUE)
      df[, StdResidual := if (!is.na(rsd) && rsd > 0) Residual / rsd else NA_real_]
    } else df[, StdResidual := NA_real_]
    res[[i]] <- df[, .(SampleID, rsid = rs, Z, StdResidual)]
  }
  res <- res[!vapply(res, is.null, logical(1))]
  say("   variants scored: %d ; skipped: %d", length(res), failed)
  dz <- rbindlist(res)
  s <- dz[, .(MeanAbsZ = mean(abs(Z), na.rm = TRUE),
              MedianAbsZ = median(abs(Z), na.rm = TRUE),
              MaxAbsZ = max(abs(Z), na.rm = TRUE),
              MedianAbsResidual = median(abs(StdResidual), na.rm = TRUE),
              N_Prots = .N,
              N_ValidResiduals = sum(!is.na(StdResidual))), by = SampleID]
  s[, FINNGENID := map$FINNGENID[match(SampleID, map$SampleID)]]
  s[]
}

# --------------------------------------------------------------------------
# 3. STAGE 1 — reproduce 05b on its own 4,611 tubes
# --------------------------------------------------------------------------
# the tube 05b actually read per individual = first occurrence in row order
rn0 <- rownames(npx)
fg0 <- map$FINNGENID[match(rn0, map$SampleID)]
first_tube <- data.table(SampleID = rn0, FINNGENID = fg0)[!is.na(FINNGENID)][, .SD[1], by = FINNGENID]
say("== identifying the tubes 05b actually read ==")
say("   first-occurrence tubes ... %d", nrow(first_tube))
say("   05b SampleID agrees ...... %d of %d  (%d name a PCA-REMOVED tube: D-48)",
    sum(orig$SampleID %in% first_tube$SampleID), nrow(orig),
    sum(!orig$SampleID %in% rn0))

v1 <- compute_z(first_tube$SampleID, "the tubes 05b actually read (first occurrence per individual)")
cmp <- merge(orig[, .(FINNGENID, MeanAbsZ_05b = MeanAbsZ, N_Prots_05b = N_Prots)],
             v1[, .(FINNGENID, MeanAbsZ_new = MeanAbsZ, N_Prots_new = N_Prots)], by = "FINNGENID")
say("== STAGE 1 validation (keyed on FINNGENID) ==")
say("   matched rows ............. %d of %d", nrow(cmp), nrow(orig))
say("   max |MeanAbsZ difference| . %.3e", max(abs(cmp$MeanAbsZ_new - cmp$MeanAbsZ_05b)))
say("   Pearson r ................ %.10f", cor(cmp$MeanAbsZ_new, cmp$MeanAbsZ_05b))
say("   N_Prots identical ........ %s", identical(cmp$N_Prots_05b, cmp$N_Prots_new))
cut1 <- round(mean(v1$MeanAbsZ) + K_SD * sd(v1$MeanAbsZ), 1)
say("   reproduced cutoff ........ %.4f -> %.1f   (05b: 1.6002 -> 1.6)",
    mean(v1$MeanAbsZ) + K_SD * sd(v1$MeanAbsZ), cut1)
say("   reproduced outliers ...... %d   (05b: 11)", sum(v1$MeanAbsZ > cut1))
VALID <- max(abs(cmp$MeanAbsZ_new - cmp$MeanAbsZ_05b)) < 1e-8 && nrow(cmp) == nrow(orig)
say("   VALIDATION: %s", if (VALID) "PASS" else "FAIL — stage 2 not trustworthy")
if (!VALID) say("   (stage 2 still computed, but treat its numbers as unverified)")

# --------------------------------------------------------------------------
# 4. STAGE 2 — every tube
# --------------------------------------------------------------------------
v2 <- compute_z(rownames(npx), "every tube in the PCA-cleaned matrix")
pop_mean <- mean(v2$MeanAbsZ, na.rm = TRUE); pop_sd <- sd(v2$MeanAbsZ, na.rm = TRUE)
cut2 <- round(pop_mean + K_SD * pop_sd, 1)
v2[, Outlier_MeanAbsZ := MeanAbsZ > cut2]
say("== STAGE 2 — per-sample provenance check ==")
say("   tubes scored .............. %d  (05b scored %d; +%d previously unchecked)",
    nrow(v2), nrow(orig), nrow(v2) - nrow(orig))
say("   individuals covered ....... %d", uniqueN(v2$FINNGENID))
say("   population MeanAbsZ ....... mean %.6f  sd %.6f", pop_mean, pop_sd)
say("   cutoff .................... %.4f -> %.1f", pop_mean + K_SD * pop_sd, cut2)
say("   outliers .................. %d  (%.3f%%)", sum(v2$Outlier_MeanAbsZ), 100 * mean(v2$Outlier_MeanAbsZ))

newly <- v2[Outlier_MeanAbsZ & !SampleID %in% orig$SampleID]
say("   of which NEWLY checked tubes: %d", nrow(newly))
say("   of which 05b already flagged: %d of its 11",
    v2[Outlier_MeanAbsZ & SampleID %in% orig$SampleID[orig$Outlier_MeanAbsZ], .N])

# --------------------------------------------------------------------------
# 5. Repeat-sampled individuals: within-person discordance
# --------------------------------------------------------------------------
v2[, n_tubes := .N, by = FINNGENID]
rep <- v2[n_tubes > 1]
say("== within-individual agreement across aliquots ==")
say("   individuals with >1 scored tube: %d covering %d tubes",
    uniqueN(rep$FINNGENID), nrow(rep))
disc <- rep[, .(n = .N, flagged = sum(Outlier_MeanAbsZ),
                spread = max(MeanAbsZ) - min(MeanAbsZ)), by = FINNGENID]
say("   individuals where aliquots DISAGREE on the flag: %d", disc[flagged > 0 & flagged < n, .N])
say("   median within-individual MeanAbsZ spread: %.4f (95th pct %.4f, max %.4f)",
    median(disc$spread), quantile(disc$spread, .95), max(disc$spread))

fwrite(v2, O("outliers/05c_pqtl_stats_per_sample.tsv"), sep = "\t")
fwrite(v2[Outlier_MeanAbsZ == TRUE][order(-MeanAbsZ)],
       O("outliers/05c_pqtl_outliers_per_sample.tsv"), sep = "\t")
fwrite(disc[order(-spread)], O("outliers/05c_pqtl_within_individual_spread.tsv"), sep = "\t")
saveRDS(list(stage1 = v1, stage2 = v2, cutoff = cut2, k_sd = K_SD, validated = VALID),
        O("outliers/05c_pqtl_per_sample.rds"))
say("== wrote 05c_pqtl_stats_per_sample.tsv, 05c_pqtl_outliers_per_sample.tsv, 05c_pqtl_within_individual_spread.tsv ==")
