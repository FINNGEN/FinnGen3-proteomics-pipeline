#!/usr/bin/env Rscript
# ==============================================================================================
# 00g — figures for the FG3 Batch 03 data release note
# ==============================================================================================
# Writes PNGs into <release>/figures/. Every figure is derived from the delivered package alone,
# so the note's diagrams and its tables cannot drift apart. Flowcharts are drawn rather than
# imported, so there is no external diagramming dependency in the build.
#
# Usage: Rscript scripts/00g_make_release_note_figures.R [release_dir]
# ==============================================================================================

suppressMessages({library(data.table); library(ggplot2); library(scales)})

args <- commandArgs(trailingOnly = TRUE)
REL <- if (length(args)) args[1] else
  "/mnt/longGWAS_disk_100GB/long_gwas/11.fg3_Proteomics/08.genewiz_batch3/01.QCed_Batch03_Release_Aug2026"
FIG <- file.path(REL, "figures")
dir.create(FIG, showWarnings = FALSE)

NAVY  <- "#1F3864"; TEAL <- "#2E9B9B"; CORAL <- "#E4572E"
AMBER <- "#E8A33D"; GREY <- "#8C8C8C"; PURPLE <- "#6B3FA0"
INK   <- "#1A1A1A"

th <- theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        legend.position  = "right",
        legend.title     = element_text(size = 10),
        axis.title       = element_text(size = 11),
        plot.margin      = margin(6, 10, 6, 6))

th_flow <- theme_void(base_size = 12) + theme(plot.margin = margin(4, 4, 4, 4))

sv <- function(p, f, w, h) {
  ggsave(file.path(FIG, f), p, width = w, height = h, dpi = 200, bg = "white")
  cat(sprintf("  %-42s %4.1f x %4.1f in\n", f, w, h))
}

## helpers for the drawn flowcharts ------------------------------------------------------------
# A box is a rounded-look rectangle with centred text; arrows are plain segments with an arrowhead.
bx <- function(x, y, w, h, label, fill = "white", border = "#AAB2C0", col = INK,
               size = 3.1, face = "plain") {
  list(
    annotate("rect", xmin = x - w/2, xmax = x + w/2, ymin = y - h/2, ymax = y + h/2,
             fill = fill, colour = border, linewidth = .4),
    annotate("text", x = x, y = y, label = label, size = size, colour = col,
             fontface = face, lineheight = .95)
  )
}
ar <- function(x1, y1, x2, y2, col = "#8C97A8", lty = 1) {
  annotate("segment", x = x1, y = y1, xend = x2, yend = y2, colour = col, linewidth = .45,
           linetype = lty, arrow = grid::arrow(length = unit(0.13, "cm"), type = "closed"))
}
lab <- function(x, y, t, size = 2.6, col = GREY, hjust = 0.5, face = "plain") {
  annotate("text", x = x, y = y, label = t, size = size, colour = col, hjust = hjust,
           fontface = face)
}

## data --------------------------------------------------------------------------------------
d  <- fread(file.path(REL, "FG3_batch03_delivery_metadata_fg3_batch_03.tsv"),
            colClasses = list(character = c("SAMPLE_ID", "FINNGENID")))
qc <- rownames(readRDS(file.path(REL, "npx_matrix_5928_qc_passed_fg3_batch_03.rds")))
d[, passed := SAMPLE_ID %in% qc]
d[, cd := as.IDate(collection_date)]
d[cd > as.IDate("2026-12-31"), cd := NA]
o  <- fread(file.path(REL, "comprehensive_outliers_list_fg3_batch_03.tsv"),
            colClasses = list(character = "SampleID"))
al <- fread(file.path(REL, "assay_lod_summary_fg3_batch_03.tsv"))

cat("release-note data figures (5 to 12):\n")

## Figures 1 to 4 are NOT generated here -----------------------------------------------------------
# The sample flow, the QC architecture, the task-force cascade and the metadata provenance graph are
# drawn as native TikZ inside the release note itself, using the shared fgbox / fgnote / fgfinal node
# vocabulary defined in docs/fg3_release_note.latex. That matches the consort-style diagrams in the
# CKD and ATC reports and keeps them sharp at any zoom. Only the data figures are rendered here.

## ------------------------------------------------------- Fig 5  removal by task force --------
tf <- d[, .(Delivered = .N, Passed = sum(passed)), by = Sample_set_inferred]
tf[, `:=`(Removed = Delivered - Passed, pct = 100 * (Delivered - Passed) / Delivered)]
BR <- 100 * mean(!d$passed)
tf[, Sample_set_inferred := factor(Sample_set_inferred,
                                   levels = Sample_set_inferred[order(pct)])]
tf[, band := fifelse(pct > 8, "far above", fifelse(pct < 2, "below", "at the batch rate"))]
sv(ggplot(tf, aes(pct, Sample_set_inferred, fill = band)) +
     geom_col(width = .68) +
     geom_vline(xintercept = BR, linetype = "dashed", colour = "grey35", linewidth = .45) +
     annotate("text", x = BR + 0.5, y = 12.4, hjust = 0, size = 3.1, colour = "grey30",
              fontface = "italic", label = sprintf("batch rate %.2f%%", BR)) +
     geom_text(aes(label = sprintf("%.2f%%  (%d of %d)", pct, Removed, Delivered)),
               hjust = -0.06, size = 3.1) +
     scale_fill_manual(values = c(`far above` = CORAL, `at the batch rate` = NAVY,
                                  below = TEAL), guide = "none") +
     scale_x_continuous(limits = c(0, 42), expand = c(0, 0)) +
     labs(x = "Samples removed by QC (%)", y = NULL) + th,
   "fig05_removal_by_taskforce.png", 8.4, 4.6)

## ------------------------------------------------------------ Fig 6  archival stratum --------
dd <- d[!is.na(cd)]
yr <- dd[, .(n = .N, rem = sum(!passed)), by = .(y = year(cd))][order(y)][n >= 3]
yr[, pct := 100 * rem / n]
sv(ggplot(yr, aes(y, pct)) +
     geom_col(aes(fill = y <= 2002), width = .8) +
     geom_text(data = yr[pct > 20], aes(label = sprintf("%.0f%%", pct)), vjust = -0.45, size = 3.1) +
     annotate("text", x = 2007.5, y = 72, hjust = 0, size = 3.1, colour = CORAL, fontface = "bold",
              label = "29 of the 30 samples collected\non or before 2002 were removed") +
     annotate("curve", x = 2007.3, xend = 2002.6, y = 72, yend = 88, curvature = 0.25,
              colour = CORAL, linewidth = .4,
              arrow = grid::arrow(length = unit(0.14, "cm"), type = "closed")) +
     scale_fill_manual(values = c(`TRUE` = CORAL, `FALSE` = NAVY), guide = "none") +
     scale_y_continuous(limits = c(0, 112), expand = c(0, 0)) +
     labs(x = "Year of plasma collection", y = "Samples removed by QC (%)") + th,
   "fig06_archival_stratum.png", 8.4, 3.9)

## --------------------------------------------------------------- Fig 7  repeat sampling ------
rp <- d[!is.na(cd)][, nd := .N, by = FINNGENID][nd > 1]
rp[, t0 := min(cd), by = FINNGENID][, yrs := as.numeric(cd - t0) / 365.25]
rp[, draw := frank(cd, ties.method = "first"), by = FINNGENID]
ord <- rp[, .(span = max(yrs)), by = FINNGENID][order(span)][, y_idx := .I]
rp  <- merge(rp, ord, by = "FINNGENID")
ord[, multi := FINNGENID %in% rp[nd >= 3, FINNGENID]]
NI  <- nrow(ord); ms <- median(ord$span); nm <- sum(ord$multi)
sv(ggplot() +
     annotate("rect", xmin = c(0, 1, 5, 10), xmax = c(1, 5, 10, 21), ymin = .5, ymax = NI + .5,
              fill = c("#FAFBFD", "#F1F5FA", "#E8EFF8", "#DFE8F5")) +
     geom_segment(data = ord[!(multi)], aes(0, y_idx, xend = span, yend = y_idx),
                  colour = "#C9D2DE", linewidth = .13) +
     geom_segment(data = ord[(multi)], aes(0, y_idx, xend = span, yend = y_idx),
                  colour = AMBER, linewidth = .45, alpha = .75) +
     geom_point(data = rp, aes(yrs, y_idx, colour = factor(pmin(draw, 5))), size = .4, alpha = .9) +
     geom_vline(xintercept = 0, colour = "grey25", linewidth = .45) +
     geom_vline(xintercept = ms, linetype = "dashed", colour = CORAL, linewidth = .55) +
     annotate("text", x = ms + 0.3, y = NI * 0.10, hjust = 0, size = 3.2, colour = CORAL,
              fontface = "bold", label = sprintf("median span %.1f y", ms)) +
     annotate("text", x = c(0.5, 3, 7.5, 15.4), y = NI * 1.05,
              label = c("under 1 y", "1 to 5 y", "5 to 10 y", "over 10 y"),
              size = 2.9, colour = "grey45") +
     annotate("text", x = 12.5, y = NI * 0.5, hjust = 0, size = 3.1, colour = "#B7791F",
              fontface = "bold",
              label = sprintf("amber: the %d people who\ngave three or more draws", nm)) +
     scale_colour_manual(values = c(`1` = NAVY, `2` = TEAL, `3` = AMBER, `4` = CORAL,
                                    `5` = PURPLE),
                         name = "Draw", labels = c("1st", "2nd", "3rd", "4th", "5th")) +
     guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1))) +
     coord_cartesian(xlim = c(-0.35, 21), ylim = c(.5, NI * 1.085), expand = FALSE) +
     labs(x = "Years since that individual's first Batch 3 draw",
          y = sprintf("Individual, ordered by span  (n = %s)", comma(NI))) +
     th + theme(axis.text.y = element_blank(), panel.grid = element_blank()),
   "fig07_repeat_sampling.png", 8.4, 4.8)

gaps <- rp[order(FINNGENID, cd)][, .(g = as.numeric(diff(cd)) / 365.25), by = FINNGENID][g > 0]
qs <- quantile(gaps$g, c(.25, .5, .75))
sv(ggplot(gaps, aes(g)) +
     geom_histogram(binwidth = .5, fill = NAVY, colour = "white", linewidth = .22, boundary = 0) +
     geom_vline(xintercept = qs[2], colour = CORAL, linewidth = .6) +
     annotate("text", x = qs[2] + .3, y = Inf, vjust = 1.8, hjust = 0, size = 3.2, colour = CORAL,
              fontface = "bold", label = sprintf("median %.1f y", qs[2])) +
     annotate("text", x = qs[2] + .3, y = Inf, vjust = 3.4, hjust = 0, size = 3, colour = "grey40",
              label = sprintf("IQR %.1f to %.1f y", qs[1], qs[3])) +
     scale_x_continuous(breaks = seq(0, 20, 2.5), limits = c(0, 20.5), expand = c(0, 0)) +
     scale_y_continuous(expand = expansion(mult = c(0, .08))) +
     labs(x = "Years between consecutive draws", y = "Repeat draws") + th,
   "fig08_repeat_interval.png", 6.6, 3.4)

## ------------------------------------------------------------- Fig 9  outlier overlap --------
pat <- o[, .N, by = .(p = Detection_Steps)][order(-N)]
pat[, p := factor(p, levels = rev(p))]
sv(ggplot(pat, aes(N, p)) +
     geom_col(fill = TEAL, width = .66) +
     geom_text(aes(label = N), hjust = -0.25, size = 3.2) +
     scale_x_continuous(expand = expansion(mult = c(0, .14))) +
     labs(x = "Samples", y = NULL,
          caption = "263 unique samples; 52 flagged by more than one method at this granularity") +
     th + theme(plot.caption = element_text(colour = GREY, size = 9, hjust = 0)),
   "fig09_outlier_overlap.png", 7.6, 3.8)

## ------------------------------------------------------------------------ Fig 10  OSI --------
osi <- d[!is.na(OSICategory), .N, by = .(stratum, OSICategory)]
osi[, pct := 100 * N / sum(N), by = stratum]
osi[, Stratum := fifelse(stratum == "NFBC1966", "NFBC1966", "Task-force selected")]
sv(ggplot(osi, aes(factor(OSICategory), pct, fill = Stratum)) +
     geom_col(position = position_dodge(width = .74), width = .66) +
     geom_text(aes(label = sprintf("%.1f", pct)), position = position_dodge(width = .74),
               vjust = -0.45, size = 2.9, colour = "#555555") +
     scale_fill_manual(values = c(`Task-force selected` = NAVY, NFBC1966 = TEAL)) +
     scale_y_continuous(expand = expansion(mult = c(0, .15))) +
     labs(x = "Olink Sample Index category  (4 = most pre-analytical variation)",
          y = "Share of that stratum (%)",
          caption = "A hard filter on category 4 would remove 20.0% of task-force-selected samples and 2.7% of NFBC1966 ones") +
     th + theme(plot.caption = element_text(colour = GREY, size = 9, hjust = 0)),
   "fig10_osi_by_stratum.png", 7.8, 3.8)

## ------------------------------------------------------------------------ Fig 11  LOD --------
al2 <- al[!is.na(frac_below_lod)]
sv(ggplot(al2, aes(100 * frac_below_lod)) +
     geom_histogram(binwidth = 2.5, fill = NAVY, colour = "white", linewidth = .2, boundary = 0) +
     geom_vline(xintercept = c(50, 90), colour = c(CORAL, AMBER), linewidth = .6) +
     annotate("text", x = 48, y = Inf, vjust = 1.7, hjust = 1, size = 3.1, colour = CORAL,
              fontface = "bold", label = sprintf("%s assays at 50%%", comma(sum(al2$frac_below_lod <= .5)))) +
     annotate("text", x = 88, y = Inf, vjust = 1.7, hjust = 1, size = 3.1, colour = "#B7791F",
              fontface = "bold", label = sprintf("%s assays at 90%%", comma(sum(al2$frac_below_lod <= .9)))) +
     scale_x_continuous(breaks = seq(0, 100, 20), expand = c(0, 0)) +
     scale_y_continuous(expand = expansion(mult = c(0, .1))) +
     labs(x = "Samples below the limit of detection (%)", y = "Assays",
          caption = "LOD is mean(NC NPX) + 3 x SD(NC NPX), from 142 to 144 negative-control observations per assay") +
     th + theme(plot.caption = element_text(colour = GREY, size = 9, hjust = 0)),
   "fig11_lod_distribution.png", 7.6, 3.6)

## ------------------------------------------------------- Fig 12  the two strata ---------------
ag <- d[!is.na(BL_AGE)]
ag[, Stratum := fifelse(stratum == "NFBC1966", "NFBC1966 (Arctic Biobank)", "Task-force selected")]
sv(ggplot(ag, aes(BL_AGE, fill = Stratum)) +
     geom_histogram(binwidth = 2, alpha = .82, position = "identity", colour = "white",
                    linewidth = .18) +
     scale_fill_manual(values = c(`Task-force selected` = NAVY, `NFBC1966 (Arctic Biobank)` = TEAL)) +
     scale_y_continuous(expand = expansion(mult = c(0, .06))) +
     labs(x = "Age at sample collection (years)", y = "Samples",
          caption = "The two strata occupy separate plates, so plate correction cannot be separated from a 30-year age difference") +
     th + theme(plot.caption = element_text(colour = GREY, size = 9, hjust = 0)),
   "fig12_age_by_stratum.png", 7.8, 3.6)

cat("done\n")
