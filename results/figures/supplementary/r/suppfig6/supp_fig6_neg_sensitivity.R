#!/usr/bin/env Rscript
# Supplementary Fig. 6  -  Negative-mode sensitivity decomposition vs positive mode.
#   a  Five-kingdom composition (Archaea excluded): expected vs POS (adopted) vs
#      NEG (ceramide-anchored, phyla with >=5 reference samples). Bray-Curtis
#      distance NEG vs POS = 0.10.
#   b  Threshold robustness: Bray-Curtis distance (vs POS and vs expected) across
#      the minimum-atlas-samples-per-phylum cutoff; n>=5 is the adopted threshold.
#   c  Per-sample stability of the NEG composition (mean, SD, CV% per kingdom);
#      Animalia CV ~89%.
#   d  MS2 spectral-confirmation rate (cosine>=0.7, >=4 matched peaks) per kingdom
#      (~36% overall, near-uniform across kingdoms).
# Rebuilt on the locked release ncbi-phylum-2026-08-04-v1 (strict 16 analysis
# phyla) with the negative-mode SIMPER fingerprint atlas recomputed on those
# units. Stage 1 reproduced the published panels exactly (max |delta| = 0),
# so the differences here are attributable to the taxonomy correction.
# Producers: paper2_repro/scripts/suppfig6_rebuild_neg_simper.py,
#            paper2_repro/scripts/suppfig6_stage2_strict16.py
# Nature Communications: 183 mm wide, all text 5-7 pt, Helvetica.

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(scales)
})

.args <- commandArgs(trailingOnly = FALSE)
.f <- sub("--file=", "", .args[grep("--file=", .args)])
SCRIPT_DIR <- if (length(.f) > 0) dirname(normalizePath(.f)) else getwd()
source(file.path(SCRIPT_DIR, "soilmass_style.R"))

# -- Data ---------------------------------------------------------
data_dir <- file.path(SCRIPT_DIR, "data", "supp_neg_sensitivity")
out_dir  <- file.path(SCRIPT_DIR, "out")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Five-kingdom subset (Archaea excluded), order matches build_figure.py K5.
K5 <- c("Bacteria", "Fungi", "Plantae", "Animalia", "Protozoa")
K5_COL <- KINGDOM_COLOURS[K5]

EXP_C <- "#BFBFBF"                 # expected (renorm.) grey, matches matplotlib EXP_C
POS_C <- MODE_COLOURS["Positive"]  # POS adopted
NEG_C <- MODE_COLOURS["Negative"]  # NEG (Cer, n>=5)
ANNOT_GREY <- "#666666"
SPINE_GREY <- "#777777"

# ================================================================
#  Panel a: 5-kingdom composition (expected / POS / NEG)
# ================================================================
comp <- read.csv(file.path(data_dir, "pos_neg_5kingdom_comparison.csv"),
                 row.names = 1, check.names = FALSE)
bc_pos <- comp["BC_vs_POS", "NEG_noarch_n>=5"]   # 0.213

comp_long <- data.frame(
  kingdom = factor(rep(K5, 3), levels = K5),
  series  = factor(rep(c("Expected (renorm.)", "POS adopted", "NEG (Cer, n>=5)"),
                       each = length(K5)),
                   levels = c("Expected (renorm.)", "POS adopted", "NEG (Cer, n>=5)")),
  value   = c(comp[K5, "expected_5k"],
              comp[K5, "POS_adopted_5k"],
              comp[K5, "NEG_noarch_n>=5"])
)

SERIES_COL <- c("Expected (renorm.)" = EXP_C,
                "POS adopted"        = unname(POS_C),
                "NEG (Cer, n>=5)"    = unname(NEG_C))

pa <- ggplot(comp_long, aes(kingdom, value, fill = series)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.74) +
  scale_fill_manual(values = SERIES_COL, name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
  annotate("label", x = Inf, y = Inf,
           label = sprintf("NEG vs POS\nBray-Curtis = %.3f", bc_pos),
           hjust = 1.02, vjust = 1.05, size = pt2mm(5.5), colour = "black",
           fill = "white", label.size = 0.25, lineheight = 0.95) +
  annotate("text", x = Inf, y = Inf,
           label = "Archaea excluded\n(no validated NEG reference)",
           hjust = 1.02, vjust = 2.7, size = pt2mm(5), colour = ANNOT_GREY,
           fontface = "italic", lineheight = 0.95) +
  labs(x = NULL, y = "Composition (%)") +
  theme_nature() +
  theme(
    axis.text.x = element_text(angle = 20, hjust = 1),
    legend.position        = "inside",
    legend.position.inside = c(0.985, 0.62),
    legend.justification   = c(1, 1),
    legend.text     = element_text(size = 5),
    legend.key.size = unit(2.6, "mm"),
    legend.background = element_blank()
  )

# ================================================================
#  Panel b: threshold robustness (BC vs POS and vs expected)
# ================================================================
thr <- read.csv(file.path(data_dir, "robustness_C_threshold.csv"), check.names = FALSE)
thr_long <- thr |>
  select(min_n, BC_vs_POS, BC_vs_exp) |>
  pivot_longer(c(BC_vs_POS, BC_vs_exp), names_to = "metric", values_to = "bc") |>
  mutate(metric = factor(ifelse(metric == "BC_vs_POS", "BC vs POS", "BC vs expected"),
                         levels = c("BC vs POS", "BC vs expected")))

THR_COL <- c("BC vs POS" = unname(NEG_C), "BC vs expected" = SPINE_GREY)
THR_SHP <- c("BC vs POS" = 16,            "BC vs expected" = 15)
THR_LTY <- c("BC vs POS" = "solid",       "BC vs expected" = "dashed")

pb <- ggplot(thr_long, aes(min_n, bc, colour = metric, shape = metric, linetype = metric)) +
  geom_vline(xintercept = 5, colour = "black", linetype = "dotted", linewidth = 0.3) +
  geom_line(linewidth = 0.45) +
  geom_point(size = 1.4) +
  annotate("text", x = 5.12, y = max(thr_long$bc) * 0.985,
           label = "n>=5 (adopted)", hjust = 0, vjust = 1, size = pt2mm(5),
           colour = "black") +
  scale_colour_manual(values = THR_COL, name = NULL) +
  scale_shape_manual(values = THR_SHP, name = NULL) +
  scale_linetype_manual(values = THR_LTY, name = NULL) +
  scale_x_continuous(breaks = sort(unique(thr_long$min_n))) +
  labs(x = "min atlas samples per phylum", y = "Bray-Curtis distance") +
  theme_nature() +
  theme(
    legend.position        = "inside",
    legend.position.inside = c(0.985, 0.42),
    legend.justification   = c(1, 1),
    legend.text     = element_text(size = 5),
    legend.key.size = unit(3.2, "mm"),
    legend.background = element_blank()
  )

# ================================================================
#  Panel c: per-sample stability (mean +/- SD, CV% labels)
# ================================================================
ps <- read.csv(file.path(data_dir, "robustness_B_persample.csv"),
               row.names = 1, check.names = FALSE)
ps <- ps[K5, ]
psd <- data.frame(
  kingdom = factor(K5, levels = K5),
  mean    = ps[K5, "mean_%"],
  sd      = ps[K5, "SD"],
  cv      = ps[K5, "CV_%"]
)
psd$cv_lab <- sprintf("CV %.0f%%", psd$cv)

pc <- ggplot(psd, aes(kingdom, mean, fill = kingdom)) +
  geom_col(width = 0.74, colour = "white", linewidth = 0.2) +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.28, linewidth = 0.3) +
  geom_text(aes(y = mean + sd + 1.4, label = cv_lab), size = pt2mm(5), colour = "black") +
  scale_fill_manual(values = K5_COL, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(x = NULL, y = "NEG composition (%)") +
  theme_nature() +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))

# ================================================================
#  Panel d: MS2 spectral-confirmation rate per kingdom
# ================================================================
# Aggregate the per-feature confirmation table exactly as build_figure.py does:
#   by  = groupby(kingdom)['confirmed'].mean() * 100
#   cnt = groupby(kingdom)['confirmed'].count()
ms2 <- read.csv(file.path(data_dir, "neg_ms2_confirmation_per_feature.csv"),
                check.names = FALSE)
ms2$confirmed <- as.logical(ms2$confirmed)
agg <- ms2 |>
  filter(kingdom %in% K5) |>
  group_by(kingdom) |>
  summarise(pct = mean(confirmed) * 100, n = dplyr::n(), .groups = "drop")
agg <- agg[match(K5, agg$kingdom), ]
agg$kingdom <- factor(agg$kingdom, levels = K5)
agg$lab <- sprintf("%.0f%%\n(n=%d)", agg$pct, agg$n)
overall_pct <- mean(ms2$confirmed[ms2$kingdom %in% K5]) * 100

pd <- ggplot(agg, aes(kingdom, pct, fill = kingdom)) +
  geom_col(width = 0.74, colour = "white", linewidth = 0.2) +
  geom_text(aes(y = pct + 5.5, label = lab), size = pt2mm(5), colour = "black",
            lineheight = 0.9) +
  scale_fill_manual(values = K5_COL, guide = "none") +
  scale_y_continuous(limits = c(0, 100), expand = expansion(mult = c(0, 0))) +
  labs(x = NULL, y = "% features MS2-confirmed in atlas") +
  theme_nature() +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))

# ================================================================
#  Compose 2 x 2
# ================================================================
fig <- (pa | pb) / (pc | pd) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = NC_TAG_PT, face = "bold"))

# -- Save (pdf + png + svg) ---------------------------------------
save_all <- function(plot, stem, w_mm, h_mm) {
  w <- w_mm * MM_TO_IN; h <- h_mm * MM_TO_IN
  ggsave(paste0(stem, ".pdf"), plot, width = w, height = h, device = cairo_pdf, bg = "white")
  ggsave(paste0(stem, ".png"), plot, width = w, height = h, dpi = 600, bg = "white",
         device = ragg::agg_png)
  if (requireNamespace("svglite", quietly = TRUE)) {
    ggsave(paste0(stem, ".svg"), plot, width = w, height = h, device = svglite::svglite, bg = "white")
  } else {
    ggsave(paste0(stem, ".svg"), plot, width = w, height = h, device = grDevices::svg, bg = "white")
  }
}

save_all(fig, file.path(out_dir, "Supplementary_Fig6_neg_sensitivity"), NC_DOUBLE, 120)
cat(sprintf("Supplementary Fig 6 -> %s  (MS2 overall = %.1f%%)\n", out_dir, overall_pct))
