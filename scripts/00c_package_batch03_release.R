#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------------------------------
# Step 00c: Package the Batch 03 Phase 1 QC release
#
# Takes the Step 05d outputs (which the pipeline writes under the single-batch layout, without the
# batch subdirectory or the _fg3_batch_03 suffix - defect D-14) and assembles the release package under
# the protocol section 8 naming convention, matching the Batch 02 December-2025 release.
#
# Emits three NPX matrices, each in rds / parquet / tsv:
#   npx_matrix_all_<N>_samples_fg3_batch_03            all delivered biological samples
#   npx_matrix_<M>_qc_passed_fg3_batch_03              after the union of all seven QC methods
#   npx_matrix_<K>_qc_passed_thl_stratum_fg3_batch_03  QC-passed restricted to the THL-metadata-complete
#                                                      stratum (the Batch 03 analogue of Batch 02's
#                                                      cohort-filtered 2,230 matrix; Batch 03 has no
#                                                      F64 or chromosomal-abnormality samples)
# plus above-LOD assay subsets, the annotated metadata, the comprehensive outlier list, the assay
# inventory, the per-assay LOD table, the tube-level sample mapping, the bridge candidate table and a
# counts report for the release note.
# ---------------------------------------------------------------------------------------------------

suppressPackageStartupMessages({ library(arrow); library(data.table); library(optparse) })

opt <- parse_args(OptionParser(option_list = list(
  make_option("--run_dir", type = "character", help = "pipeline output base_dir"),
  make_option("--staging", type = "character", help = "step 00a/00b output dir"),
  make_option("--outdir", type = "character", help = "release directory"),
  make_option("--batch", type = "character", default = "fg3_batch_03")
)))
say <- function(...) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), paste0(...)))
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)
B <- opt$batch

pick <- function(...) { for (p in c(...)) if (file.exists(p)) return(p); NA_character_ }
qc <- function(f) file.path(opt$run_dir, "qc", f)
ou <- function(f) file.path(opt$run_dir, "outliers", f)

write3 <- function(m, stem) {
  saveRDS(m, file.path(opt$outdir, paste0(stem, ".rds")))
  dt <- as.data.table(m, keep.rownames = "SampleID")
  write_parquet(dt, file.path(opt$outdir, paste0(stem, ".parquet")))
  fwrite(dt, file.path(opt$outdir, paste0(stem, ".tsv")), sep = "\t")
  say("  wrote ", stem, " (", nrow(m), " x ", ncol(m), ") in rds/parquet/tsv")
}

# ---- inputs ---------------------------------------------------------------------------------------
p_all   <- file.path(opt$staging, paste0("00a_npx_matrix_raw_", B, ".rds"))
p_meta  <- pick(qc("05d_qc_annotated_metadata.rds"), ou("05d_qc_annotated_metadata.rds"))
p_out   <- pick(qc("05d_comprehensive_outliers_list.tsv"), ou("05d_comprehensive_outliers_list.tsv"))
p_pass  <- pick(qc("05d_npx_matrix_all_qc_passed.rds"), ou("05d_npx_matrix_all_qc_passed.rds"))
p_map   <- qc("00_sample_mapping.rds")
p_b0meta<- file.path(opt$staging, paste0("00b_metadata_", B, ".tsv"))
p_assay <- file.path(opt$staging, paste0("00a_assay_inventory_", B, ".tsv"))
p_lod   <- file.path(opt$staging, paste0("00a_lod_from_neg_controls_", B, ".tsv"))
for (f in c(p_all, p_b0meta, p_assay, p_lod)) if (!file.exists(f)) stop("missing required input: ", f)
say("qc-annotated metadata : ", p_meta)
say("comprehensive outliers: ", p_out)
say("qc-passed matrix      : ", p_pass)

npx_all <- readRDS(p_all)
b0      <- fread(p_b0meta, colClasses = list(character = c("SAMPLE_ID", "FINNGENID")))
assays  <- fread(p_assay)
lod     <- fread(p_lod)

# ---- QC-passed sample set -------------------------------------------------------------------------
if (!is.na(p_pass)) {
  npx_pass <- readRDS(p_pass)
  pass_ids <- rownames(npx_pass)
} else if (!is.na(p_out)) {
  ol <- fread(p_out, colClasses = list(character = "SampleID"))
  pass_ids <- setdiff(rownames(npx_all), unique(ol$SampleID))
  npx_pass <- npx_all[pass_ids, , drop = FALSE]
} else stop("neither a QC-passed matrix nor a comprehensive outlier list was found")
say("QC-passed samples: ", length(pass_ids), " of ", nrow(npx_all))

# ---- assay sets -----------------------------------------------------------------------------------
excluded <- assays[vendor_excluded == TRUE]$assay_colname
keep_assay <- setdiff(colnames(npx_all), excluded)
say("assays: ", ncol(npx_all), " delivered biological, ", length(excluded),
    " vendor-excluded purged -> ", length(keep_assay))

npx_all_k  <- npx_all[, keep_assay, drop = FALSE]
npx_pass_k <- npx_pass[, intersect(keep_assay, colnames(npx_pass)), drop = FALSE]

# ---- three released matrices ----------------------------------------------------------------------
say("writing released matrices")
write3(npx_all_k,  sprintf("npx_matrix_all_%d_samples_%s", nrow(npx_all_k), B))
write3(npx_pass_k, sprintf("npx_matrix_%d_qc_passed_%s", nrow(npx_pass_k), B))

thl_ids <- b0[metadata_completeness == "thl_complete"]$SAMPLE_ID
thl_pass <- intersect(rownames(npx_pass_k), thl_ids)
write3(npx_pass_k[thl_pass, , drop = FALSE],
       sprintf("npx_matrix_%d_qc_passed_thl_stratum_%s", length(thl_pass), B))

# ---- above-LOD subsets ----------------------------------------------------------------------------
say("computing above-LOD assay subsets")
lodv <- lod[!is.na(LOD) & is_control_probe == FALSE & vendor_excluded == FALSE,
            .(assay_colname, LOD)]
common <- intersect(lodv$assay_colname, colnames(npx_pass_k))
lodv <- lodv[assay_colname %in% common]
sub <- npx_pass_k[, lodv$assay_colname, drop = FALSE]
below <- vapply(seq_along(lodv$assay_colname), function(j)
  mean(sub[, j] < lodv$LOD[j], na.rm = TRUE), numeric(1))
lodv[, frac_below_lod := below]
fwrite(lodv[order(frac_below_lod)], file.path(opt$outdir, sprintf("assay_lod_summary_%s.tsv", B)), sep = "\t")
for (thr in c(0.90, 0.50)) {
  keep <- lodv[frac_below_lod <= thr]$assay_colname
  say("  <= ", thr * 100, "% below LOD: ", length(keep), " assays")
  write3(npx_pass_k[, keep, drop = FALSE],
         sprintf("npx_matrix_above_lod_%d_samples_%d_proteins_%s", nrow(npx_pass_k), length(keep), B))
  writeLines(keep, file.path(opt$outdir, sprintf("assay_list_above_lod_%d_%s.txt", length(keep), B)))
}

# ---- annotated metadata, outlier list, ancillary tables -------------------------------------------
if (!is.na(p_meta)) {
  am <- as.data.table(readRDS(p_meta))
  idc <- intersect(c("SampleID", "SAMPLE_ID"), names(am))[1]
  am <- merge(am, b0, by.x = idc, by.y = "SAMPLE_ID", all.x = TRUE, suffixes = c("", "_b0"))
  fwrite(am, file.path(opt$outdir, sprintf("qc_annotated_metadata_all_%d_samples_%s.tsv", nrow(am), B)), sep = "\t")
  write_parquet(am, file.path(opt$outdir, sprintf("qc_annotated_metadata_all_%d_samples_%s.parquet", nrow(am), B)))
  say("annotated metadata: ", nrow(am), " x ", ncol(am))
} else {
  fwrite(b0, file.path(opt$outdir, sprintf("qc_annotated_metadata_all_%d_samples_%s.tsv", nrow(b0), B)), sep = "\t")
  say("WARNING: 05d annotated metadata absent; shipped the step-00b metadata instead")
}
if (!is.na(p_out)) {
  ol <- fread(p_out)
  fwrite(ol, file.path(opt$outdir, sprintf("comprehensive_outliers_list_%s.tsv", B)), sep = "\t")
  write_parquet(ol, file.path(opt$outdir, sprintf("comprehensive_outliers_list_%s.parquet", B)))
  say("outlier list: ", nrow(ol), " rows")
}
file.copy(p_b0meta, file.path(opt$outdir, sprintf("FG3_batch03_sample_metadata_schema_%s.tsv", B)), overwrite = TRUE)
file.copy(p_assay,  file.path(opt$outdir, sprintf("assay_qc_%s.tsv", B)), overwrite = TRUE)
file.copy(p_lod,    file.path(opt$outdir, sprintf("lod_from_negative_controls_%s.tsv", B)), overwrite = TRUE)
if (file.exists(p_map)) file.copy(p_map, file.path(opt$outdir, sprintf("00_sample_mapping_%s.rds", B)), overwrite = TRUE)
writeLines(keep_assay, file.path(opt$outdir, sprintf("protein_list_%d_%s.txt", length(keep_assay), B)))

# ---- counts report --------------------------------------------------------------------------------
flags <- grep("^QC_", names(if (!is.na(p_meta)) am else b0), value = TRUE)
flags <- flags[!grepl("_(pc|prob|rate|npx|z|zscore|prots|sex|residual)", flags)]
rep <- c("Batch 03 Phase 1 release packaging report",
         paste0("generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
         "",
         paste0("delivered biological samples (all matrix) : ", nrow(npx_all_k)),
         paste0("QC-passed samples                         : ", nrow(npx_pass_k)),
         paste0("THL-stratum QC-passed samples             : ", length(thl_pass)),
         paste0("assays released                           : ", length(keep_assay)),
         paste0("unique individuals (all)                  : ", uniqueN(b0$FINNGENID)),
         paste0("individuals with >1 sample                : ",
                uniqueN(b0[n_samples_this_individual > 1]$FINNGENID)),
         paste0("designated bridge rows                    : ", b0[is_designated_bridge == TRUE, .N]))
if (!is.na(p_meta)) for (f in flags)
  rep <- c(rep, sprintf("  %-40s : %s", f, sum(am[[f]] %in% c(TRUE, "TRUE", 1), na.rm = TRUE)))
writeLines(rep, file.path(opt$outdir, sprintf("packaging_report_%s.txt", B)))
cat(paste(rep, collapse = "\n"), "\n")
say("release packaged in ", opt$outdir)
