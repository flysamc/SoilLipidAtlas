#!/usr/bin/env Rscript
# Fig 2 -- USER-DIRECTED LAYOUT VARIANT (2026-08-11), strict release
# ncbi-phylum-2026-08-04-v1.
#
# This departs from the submitted fig2_atlas.R composition on explicit author
# instruction: panel a spans the full figure height on the left with a broken
# count axis, panel b sits top-right, panel c is compact beneath it.
# Everything else is held to the submitted figure: soilmass_style.R theme_nature,
# Wong kingdom palette, annotation-tier colours, heatmap ramp, bar colours,
# 5-7 pt text, 183 mm double-column width.
#
# Axis break: the lower segment covers 0-1,150, which contains all 15 non-
# Streptophyta phyla complete (max Bacillota 1,021).  The upper segment covers
# 4,800-5,300 at THE SAME units-per-mm, so the two segments share one linear
# scale and only a range is omitted -- no compression.  Streptophyta (5,168)
# is the sole bar crossing the gap; its Gold|Bronze boundary (827) is visible
# in the lower segment and its Bronze|Unidentified boundary (3,642) falls
# inside the omitted range.  Stated in the caption.
#
#   Panel a: Biomarker counts per phylum (POS), annotation release 2026-08-06
#   Panel b: fastMASST kingdom x sample-type heatmap -- REIMPLEMENTATION
#   Panel c: Soil detection by selection method, 149 Pan-ReDU soil datasets
# Producers: paper2_repro/scripts/build_fig2bc_strict.py (b, c),
#            paper2_repro/annotation_summaries.py (a),
#            paper2_repro/scripts/build_fig2_strict16_render.py (this folder).

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(forcats)
  library(scales)
  library(ggtext)
})

.args <- commandArgs(trailingOnly = FALSE)
.f <- sub("--file=", "", .args[grep("--file=", .args)])
SCRIPT_DIR <- if (length(.f) > 0) dirname(normalizePath(.f)) else getwd()
source(file.path(SCRIPT_DIR, "soilmass_style.R"))

data_dir <- file.path(SCRIPT_DIR, "data")
out_dir  <- file.path(SCRIPT_DIR, "out")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

annot_colors <- c(
  Gold         = "#C7A535",
  Silver       = "#8C8C8C",
  Bronze       = "#B87333",
  Unidentified = "#ECECEC"
)
annot_order <- c("Gold", "Silver", "Bronze", "Unidentified")

# Axis-break geometry.  BREAK_LO must exceed every non-Streptophyta total and
# BREAK_HI must sit below the Streptophyta total; asserted after load.
# LO_PAD / HI_PAD are label room, not data range: each segment's physical width
# is set to its FULL span below, so both segments keep the same units per mm.
# Panels clip normally (clip = "on"); the crossing bar is meant to run off the
# edge of each segment, which is what makes the break read as a break.
BREAK_LO <- 1150
BREAK_HI <- 4800
XMAX     <- 5300
LO_MAX   <- 1400   # 1,150 of data + label room
HI_MAX   <- 5400   # 5,168 of data + label room

# -- Data ---------------------------------------------------------
tiers <- read.csv(file.path(data_dir, "tier_counts.csv")) |>
  filter(mode == "POS") |>
  mutate(
    tier    = factor(tier, levels = rev(annot_order)),
    kingdom = factor(kingdom, levels = KINGDOM_ORDER)
  )

phylum_totals <- tiers |>
  group_by(phylum, kingdom) |>
  summarise(total = sum(n_features), .groups = "drop")

phylum_order <- phylum_totals |>
  arrange(kingdom, desc(total)) |>
  pull(phylum)
tiers$phylum <- factor(tiers$phylum, levels = rev(phylum_order))
phylum_totals$phylum <- factor(phylum_totals$phylum, levels = rev(phylum_order))

# Guard the break: exactly one phylum may cross it, and nothing may exceed XMAX.
.above <- phylum_totals |> filter(total > BREAK_LO)
if (nrow(.above) != 1L) {
  stop(sprintf("axis break invalid: %d phyla exceed %d (%s); adjust BREAK_LO",
               nrow(.above), BREAK_LO, paste(.above$phylum, collapse = ", ")))
}
if (any(phylum_totals$total > XMAX) || .above$total < BREAK_HI) {
  stop("axis break invalid: upper segment does not contain the crossing bar")
}

stype <- read.csv(file.path(data_dir, "kingdom_sampletype_summary.csv"),
                  check.names = FALSE)

heatmap_types  <- c("Animal/Clinical", "Plant", "Bacterial", "Fungal", "Algae")
heatmap_labels <- c("Animal/\nClinical", "Plant", "Bacterial", "Fungal", "Algae")

hm_data <- stype |>
  filter(kingdom %in% KINGDOM_ORDER) |>
  select(kingdom, starts_with("pct_")) |>
  pivot_longer(-kingdom, names_to = "sample_type", values_to = "pct") |>
  mutate(
    sample_type = gsub("^pct_", "", sample_type),
    kingdom     = factor(kingdom, levels = rev(KINGDOM_ORDER))
  ) |>
  filter(sample_type %in% heatmap_types) |>
  mutate(sample_type = factor(sample_type, levels = heatmap_types))

soil_comp <- read.csv(file.path(data_dir, "shared_vs_exclusive_soil.csv"))
soil_comp$category <- factor(soil_comp$category,
                             levels = c("Indicator Value (IndVal)", "Composite scoring"))

# -- Panel a: split across the axis break -------------------------
phy_kingdom <- setNames(
  as.character(phylum_totals$kingdom[match(rev(phylum_order), phylum_totals$phylum)]),
  rev(phylum_order)
)
phy_html <- sapply(rev(phylum_order), function(p) {
  col <- KINGDOM_COLOURS[phy_kingdom[p]]
  paste0("<span style='color:", col, "'>", p, "</span>")
})

pa_base <- ggplot(tiers, aes(x = n_features, y = phylum, fill = tier)) +
  geom_bar(stat = "identity", width = 0.75, colour = "white", linewidth = 0.15) +
  scale_fill_manual(
    values = annot_colors[rev(annot_order)],
    breaks = annot_order,
    name   = "Annotation tier"
  ) +
  theme_nature() +
  theme(panel.grid.major.x = element_line(colour = "grey90", linewidth = 0.15))

# break marks: two short diagonals straddling the x-axis line at the cut,
# drawn just inside the panel edge so normal clipping keeps them intact
break_mark <- function(at, span) {
  list(
    annotate("segment", x = at - span, xend = at + span * 0.4,
             y = 0.42, yend = 0.86, linewidth = 0.3, colour = "black"),
    annotate("segment", x = at - span * 0.4, xend = at + span,
             y = 0.42, yend = 0.86, linewidth = 0.3, colour = "black")
  )
}

pa_left <- pa_base +
  geom_text(
    data = phylum_totals |> filter(total <= BREAK_LO),
    aes(x = total, y = phylum, label = format(total, big.mark = ",")),
    inherit.aes = FALSE, hjust = -0.12, size = pt2mm(5), colour = "grey45"
  ) +
  scale_x_continuous(breaks = seq(0, 1000, 250), expand = expansion(mult = c(0, 0))) +
  scale_y_discrete(labels = phy_html) +
  coord_cartesian(xlim = c(0, LO_MAX)) +
  break_mark(LO_MAX - 28, 22) +
  labs(x = "Number of phylum-enriched biomarkers", y = NULL, tag = "a") +
  theme(
    axis.text.y = element_markdown(size = 5),
    legend.position = "inside",
    legend.position.inside = c(0.99, 0.015),
    legend.justification = c(1, 0),
    legend.background = element_rect(fill = alpha("white", 0.9), colour = NA),
    legend.title = element_text(size = 5, face = "bold"),
    legend.text  = element_text(size = 5),
    legend.key.size = unit(2.5, "mm"),
    plot.margin = margin(2, 0, 2, 2, "mm")
  )

pa_right <- pa_base +
  geom_text(
    data = phylum_totals |> filter(total > BREAK_LO),
    aes(x = total, y = phylum, label = format(total, big.mark = ",")),
    inherit.aes = FALSE, hjust = -0.18, size = pt2mm(5),
    fontface = "bold", colour = "grey45"
  ) +
  scale_x_continuous(breaks = c(5000), expand = expansion(mult = c(0, 0))) +
  coord_cartesian(xlim = c(BREAK_HI, HI_MAX)) +
  break_mark(BREAK_HI + 28, 22) +
  labs(x = NULL, y = NULL) +
  theme(
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y  = element_blank(),
    legend.position = "none",
    plot.margin = margin(2, 3, 2, 1.5, "mm")
  )

# One linear scale across both segments: physical width tracks the full span
# each segment displays, so a millimetre means the same count on both sides.
panel_a <- pa_left + pa_right +
  plot_layout(widths = c(LO_MAX, HI_MAX - BREAK_HI))

# -- Panel b: MASST sample-type heatmap ---------------------------
kingdom_labels_df <- data.frame(
  kingdom = factor(KINGDOM_ORDER, levels = rev(KINGDOM_ORDER)),
  sample_type = factor(heatmap_types[1], levels = heatmap_types),
  label = KINGDOM_ORDER
)

pb <- ggplot(hm_data, aes(x = sample_type, y = kingdom, fill = pct)) +
  geom_tile(colour = "white", linewidth = 0.5) +
  geom_text(aes(label = paste0(round(pct), "%")),
            size = pt2mm(5.5), fontface = "bold",
            colour = ifelse(hm_data$pct > 25, "white", "#333333")) +
  geom_text(data = kingdom_labels_df,
            aes(y = kingdom, label = label), inherit.aes = FALSE,
            x = 0.35, size = pt2mm(5.5), fontface = "bold", hjust = 1,
            colour = KINGDOM_COLOURS[kingdom_labels_df$label]) +
  scale_fill_gradientn(
    colours = c("#FFFDE7", "#FFF9C4", "#FFE082", "#FFB300", "#E65100", "#B71C1C"),
    limits  = c(0, 60),
    name    = "Detection\nrate (%)"
  ) +
  scale_x_discrete(labels = heatmap_labels) +
  labs(x = "Public sample type (fastMASST)", y = NULL, tag = "b") +
  coord_cartesian(clip = "off") +
  theme_nature() +
  theme(
    axis.text.y  = element_blank(),
    axis.text.x  = element_text(size = 5, lineheight = 0.85),
    axis.ticks.y = element_blank(),
    axis.line    = element_blank(),
    panel.grid   = element_blank(),
    legend.position = "none",
    plot.margin  = margin(2, 2, 4, 20, "mm")
  )

# -- Panel c: soil detection by selection method ------------------
bar_colors <- c("Indicator Value (IndVal)" = "#2C5F8A",
                "Composite scoring"        = "#B0B0B0")

pc <- ggplot(soil_comp, aes(x = pct_soil, y = category, fill = category)) +
  geom_col(width = 0.5, show.legend = FALSE) +
  geom_text(
    aes(label = paste0(pct_soil, "%")),
    hjust = -0.12, size = pt2mm(6), fontface = "bold", colour = "#333333"
  ) +
  geom_text(
    aes(x = 0.6, label = paste0("n = ", format(n_total, big.mark = ","))),
    hjust = 0, size = pt2mm(5), colour = "white", fontface = "italic"
  ) +
  scale_fill_manual(values = bar_colors) +
  scale_x_continuous(
    expand = expansion(mult = c(0, 0.20)),
    breaks = seq(0, 30, 10),
    labels = function(x) paste0(x, "%")
  ) +
  scale_y_discrete(labels = function(x) sub(" \\(IndVal\\)", "\n(IndVal)", x)) +
  labs(x = "Biomarkers detected in soil datasets (Pan-ReDU)", y = NULL, tag = "c") +
  theme_nature() +
  theme(
    axis.text.y = element_text(size = 6, face = "bold", colour = "#333333"),
    panel.grid.major.x = element_line(colour = "grey90", linewidth = 0.15),
    plot.margin = margin(2, 2, 2, 2, "mm")
  )

# -- Compose: a full height left, b over c on the right -----------
right_col <- pb / pc + plot_layout(heights = c(2.5, 1))

fig2 <- panel_a | right_col
fig2 <- fig2 + plot_layout(widths = c(1.25, 1)) &
  theme(plot.tag = element_text(size = NC_TAG_PT, face = "bold"))

# -- Save ---------------------------------------------------------
w <- NC_DOUBLE * MM_TO_IN
h <- 160 * MM_TO_IN
ggsave(file.path(out_dir, "fig2_atlas.pdf"), fig2,
       width = w, height = h, device = cairo_pdf, bg = "white")
ggsave(file.path(out_dir, "fig2_atlas.png"), fig2,
       width = w, height = h, dpi = 600, bg = "white",
       device = ragg::agg_png)
cat("Fig 2 (wide-a layout) ->", out_dir, "\n")
