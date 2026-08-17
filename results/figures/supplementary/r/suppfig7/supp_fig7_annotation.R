#!/usr/bin/env Rscript
# Supplementary Fig. 7  -  Annotation pipeline contributions (a) + evidence tiers (b).
#   a  Per-source annotation contributions, stacked by confidence tier
#      (Gold/Silver/Bronze), with positive-mode quality-downgrade steps as
#      leftward bars. Positive and negative mode shown side by side.
#   b  Evidence-tier table: phyla grouped A-D by sampling depth in both modes.
# Port of figures_v5 panels/{annotation_waterfall,evidence_tiers}.py + compose/ed_annotation.py.
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
data_dir <- file.path(SCRIPT_DIR, "data", "supp_annotation")
out_dir  <- file.path(SCRIPT_DIR, "out")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

ANNOT_COLORS <- c(Gold = "#C7A535", Silver = "#8C8C8C", Bronze = "#B87333")
DOWNGRADE    <- "#C0504D"
COLOR_TEXT   <- "#222222"
COLOR_LIGHT  <- "#666666"
COLOR_SPINE  <- "#BBBBBB"

# phylum -> kingdom (ported from figures_v5 soilmass_style.py; shared R style file
# carries kingdoms only)
# Locked release ncbi-phylum-2026-08-04-v1: the 19 collection phyla. Retired
# labels (Euryarchaeota, Crenarchaeota, Amoebozoa, Bryophyta, Marchantiophyta,
# Trachaeophyta, Magnoliophyta, Charophyta) and the descriptive-only label
# Bicosoecida are no longer analysis-scheme units.
PHYLUM_TO_KINGDOM <- c(
  Bacillota = "Bacteria", Actinomycetota = "Bacteria", Pseudomonadota = "Bacteria", Cyanobacteriota = "Bacteria",
  Methanobacteriota = "Archaea", Thermoproteota = "Archaea",
  Ascomycota = "Fungi", Basidiomycota = "Fungi", Mucoromycota = "Fungi", Mortierellomycota = "Fungi",
  Streptophyta = "Plantae", Chlorophyta = "Plantae",
  Arthropoda = "Animalia", Nematoda = "Animalia", Mollusca = "Animalia",
  Discosea = "Protozoa", Evosea = "Protozoa", Heterolobosea = "Protozoa", Cercozoa = "Protozoa")

# Panel b rebuilt on the locked release ncbi-phylum-2026-08-04-v1:
# 19 collection phyla, corrected assignments, taxonomic breadth as
# distinct genera (the published species column is not reproducible).
# Producer: paper2_repro/scripts/suppfig7_evidence_tiers_release.py
# Panel a is taxonomy-independent and unchanged from the submission.
# ================================================================
#  Panel a: annotation waterfall (positive / negative mode)
# ================================================================
wf <- read.csv(file.path(data_dir, "waterfall_steps.csv"), check.names = FALSE)
wf$label <- gsub("[\r\n]+", " ", wf$label)   # frozen CSV wraps labels with newlines

FILL_LEVELS <- c("Gold", "Silver", "Bronze", "Quality downgrade")

build_waterfall <- function(mode, show_legend, x_breaks) {
  sub <- wf[wf$mode == mode, ]
  n   <- nrow(sub)
  sub$y <- (n - 1):0                          # first row at top (matches arange[::-1])
  maxamt  <- max(sub$amount)
  pos_end <- max(sub$amount[sub$kind == "positive"])
  neg_max <- if (any(sub$kind == "negative")) max(sub$amount[sub$kind == "negative"]) else 0
  off     <- maxamt * 0.02
  x_left  <- if (neg_max > 0) -(neg_max + maxamt * 0.14) else -(maxamt * 0.012)
  x_right <- pos_end + maxamt * 0.30

  # stacked positive tiers -> explicit left/right so Gold sits at the axis
  pos_long <- sub[sub$kind == "positive", c("y", "gold", "silver", "bronze")] |>
    pivot_longer(c(gold, silver, bronze), names_to = "tier", values_to = "val") |>
    mutate(tier = factor(tier, levels = c("gold", "silver", "bronze"))) |>
    arrange(y, tier) |>
    group_by(y) |>
    mutate(xmax = cumsum(val), xmin = xmax - val) |>
    ungroup() |>
    filter(val > 0) |>
    mutate(fill = factor(c(gold = "Gold", silver = "Silver", bronze = "Bronze")[as.character(tier)],
                         levels = FILL_LEVELS))

  neg_rows <- sub[sub$kind == "negative", ]
  pos_rows <- sub[sub$kind == "positive", ]

  p <- ggplot() +
    geom_rect(data = pos_long,
              aes(xmin = xmin, xmax = xmax, ymin = y - 0.31, ymax = y + 0.31, fill = fill),
              colour = "white", linewidth = 0.15)

  if (nrow(neg_rows) > 0) {
    neg_rows$fill <- factor("Quality downgrade", levels = FILL_LEVELS)
    p <- p +
      geom_rect(data = neg_rows,
                aes(xmin = -amount, xmax = 0, ymin = y - 0.31, ymax = y + 0.31, fill = fill),
                colour = "white", linewidth = 0.15) +
      geom_text(data = neg_rows,
                aes(x = -amount - off, y = y, label = paste0("−", comma(amount))),
                hjust = 1, size = pt2mm(5), colour = DOWNGRADE)
  }

  p <- p +
    geom_text(data = pos_rows,
              aes(x = amount + off, y = y, label = paste0("+", comma(amount))),
              hjust = 0, size = pt2mm(5), colour = COLOR_TEXT) +
    geom_vline(xintercept = 0, colour = COLOR_SPINE, linewidth = 0.4) +
    scale_fill_manual(values = c(ANNOT_COLORS, "Quality downgrade" = DOWNGRADE),
                      breaks = FILL_LEVELS, name = NULL, drop = FALSE) +
    scale_x_continuous(breaks = x_breaks) +
    scale_y_continuous(breaks = sub$y, labels = sub$label,
                       expand = expansion(add = 0.7)) +
    coord_cartesian(xlim = c(x_left, x_right), clip = "off") +
    labs(x = "Annotations contributed", y = NULL,
         title = if (mode == "POS") "Positive mode" else "Negative mode") +
    theme_nature() +
    theme(
      axis.text.y  = element_text(size = 5.5, colour = COLOR_TEXT),
      axis.ticks.y = element_blank(),
      plot.title   = element_text(size = 7, face = "bold", hjust = 0)
    )

  if (show_legend) {
    p <- p + theme(legend.position = "inside",
                   legend.position.inside = c(0.99, 0.03),
                   legend.justification = c(1, 0),
                   legend.text = element_text(size = 5),
                   legend.key.size = unit(2.6, "mm"),
                   legend.background = element_blank())
  } else {
    p <- p + theme(legend.position = "none")
  }
  p
}

pa_pos <- build_waterfall("POS", show_legend = TRUE,  x_breaks = c(0, 1000, 2000, 3000))
pa_neg <- build_waterfall("NEG", show_legend = FALSE, x_breaks = c(0, 500, 1000, 1500))

# ================================================================
#  Panel b: evidence-tier table
# ================================================================
ev <- read.csv(file.path(data_dir, "evidence_tiers.csv"), check.names = FALSE)
ev$tier_rank <- match(ev$tier, c("A", "B", "C", "D"))
ev <- ev[order(ev$tier_rank, -ev$pos_samples), ]
ev$kingdom <- PHYLUM_TO_KINGDOM[ev$phylum]
ev$kcol    <- KINGDOM_COLOURS[ev$kingdom]
n  <- nrow(ev)
ev$i <- 0:(n - 1)
ev$yr <- n - ev$i                              # data-row y: top row = n, bottom = 1

TIER_COLORS <- c(A = "#1B7837", B = "#A6D96A", C = "#FDAE61", D = "#D73027")

# column geometry (fractions of width), mirrors evidence_tiers.py
xw <- c(0.26, 0.12, 0.12, 0.12, 0.12, 0.12, 0.14)
xl <- c(0, head(cumsum(xw), -1))
xc <- xl + xw / 2

y_div   <- n + 0.5
y_chdr  <- n + 1.25
y_ghdr  <- n + 2.15
y_foot  <- 0.05

# alternating row shading (even 0-based rows)
shade <- ev[ev$i %% 2 == 0, ]
# value columns long form
val_long <- ev |>
  transmute(yr,
            c2 = pos_samples, c3 = pos_genera, c4 = pos_batches,
            c5 = neg_samples, c6 = neg_batches) |>
  pivot_longer(c(c2, c3, c4, c5, c6), names_to = "col", values_to = "v") |>
  mutate(xc = xc[as.integer(sub("c", "", col))])

ev$name_disp <- ifelse(!is.na(ev$note) & ev$note == "*", paste0(ev$phylum, " *"), ev$phylum)

foot <- paste0("Tier A robust (≥5 samples, ≥3 genera, ≥2 batches) · B moderate · ",
               "C preliminary (n=2) · D anecdotal (n=1).   * Arthropoda: single-batch exception.")

panel_b <- ggplot() +
  geom_rect(data = shade,
            aes(xmin = 0, xmax = 1, ymin = yr - 0.5, ymax = yr + 0.5),
            fill = "#F7F7F7") +
  # group headers
  annotate("text", x = 0.44, y = y_ghdr, label = "Positive mode",
           fontface = "bold", size = pt2mm(6), colour = COLOR_TEXT) +
  annotate("text", x = 0.74, y = y_ghdr, label = "Negative mode",
           fontface = "bold", size = pt2mm(6), colour = COLOR_TEXT) +
  # column headers
  annotate("text", x = xc, y = y_chdr,
           label = c("Phylum", "Samples", "Genera", "Batches", "Samples", "Batches", "Tier"),
           fontface = "bold", size = pt2mm(5.5), colour = COLOR_LIGHT) +
  geom_segment(aes(x = 0, xend = 1, y = y_div, yend = y_div),
               colour = COLOR_SPINE, linewidth = 0.4) +
  # phylum names (kingdom-coloured)
  geom_text(data = ev, aes(x = xl[1] + 0.01, y = yr, label = name_disp),
            hjust = 0, fontface = "bold", size = pt2mm(5.5), colour = ev$kcol) +
  # values
  geom_text(data = val_long, aes(x = xc, y = yr, label = v),
            size = pt2mm(5.5), colour = COLOR_TEXT) +
  # tier badges
  geom_rect(data = ev,
            aes(xmin = xl[7] + 0.012, xmax = xl[7] + xw[7] - 0.012,
                ymin = yr - 0.34, ymax = yr + 0.34),
            fill = TIER_COLORS[ev$tier]) +
  geom_text(data = ev, aes(x = xc[7], y = yr, label = tier),
            fontface = "bold", size = pt2mm(5.5), colour = "white") +
  # footnote
  annotate("text", x = 0, y = y_foot, label = foot, hjust = 0, vjust = 1,
           size = pt2mm(5), colour = COLOR_LIGHT) +
  coord_cartesian(xlim = c(-0.005, 1.005), ylim = c(-0.7, y_ghdr + 0.6), clip = "off") +
  theme_void() +
  theme(plot.margin = margin(2, 2, 2, 2, "mm"))

# ================================================================
#  Compose
# ================================================================
# Tag the waterfall pair (a) on its left sub-panel and the table (b) directly,
# rather than wrap_elements + auto-tags, so the waterfall fills its row instead
# of being vertically centred with a gap above the table.
pa_pos  <- pa_pos  + labs(tag = "a")
panel_b <- panel_b + labs(tag = "b")

ed1 <- (pa_pos | pa_neg) / panel_b +
  plot_layout(heights = c(0.52, 1.0)) &
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

save_all(ed1, file.path(out_dir, "Supplementary_Fig7_annotation"), NC_DOUBLE, 165)
cat("Supplementary Fig 7 ->", out_dir, "\n")
