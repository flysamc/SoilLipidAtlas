#!/usr/bin/env Rscript
# Fig 5 — ClimGrass: two-estimator community composition (a) +
# phylum treatment responses (b).
#
# Panel a: pooled fc-weighted estimate (bars, 95% reference-sample
# bootstrap CI), marker-panel estimator (diamonds), literature biomass
# ranges (grey segments). Panel b: CLR fingerprint-set effect map.
# Display labels follow the locked ecological_group policy
# (Viridiplantae / Protists); CSV keys remain Plantae / Protozoa.
#
# This folder re-renders with only R: plotted CSVs live in r/data/.

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr)
  library(patchwork); library(scales)
})

.args <- commandArgs(trailingOnly = FALSE)
.f <- sub("--file=", "", .args[grep("--file=", .args)])
SCRIPT_DIR <- if (length(.f) > 0) dirname(normalizePath(.f)) else getwd()
source(file.path(SCRIPT_DIR, "soilmass_style.R"))
FONT <- "Arial"
data_dir <- file.path(SCRIPT_DIR, "data")

boot_df   <- read.csv(file.path(data_dir, "composition_fcweighted_kingdom_ci.csv"))
marker_df <- read.csv(file.path(data_dir, "kingdom_ci_marker_panel.csv"))
eff_df    <- read.csv(file.path(data_dir, "phylum_effects.csv"))

# Locked display labels (ecological_group policy).
DISPLAY <- c(Bacteria = "Bacteria", Fungi = "Fungi",
             Plantae = "Viridiplantae", Animalia = "Animalia",
             Protozoa = "Protists", Archaea = "Archaea")

# Panel-a bar order (top to bottom).
BAR_ORDER <- c("Bacteria", "Fungi", "Plantae", "Animalia", "Protozoa", "Archaea")

# Literature biomass ranges (%), Table S6; aligned to KINGDOM_ORDER
# (Bacteria, Archaea, Fungi, Plantae, Animalia, Protozoa).
EXPECTED <- data.frame(
  kingdom = KINGDOM_ORDER,
  lo = c(35, 2, 15, 15, 1, 0),
  hi = c(50, 8, 30, 30, 5, 1),
  stringsAsFactors = FALSE
)

# Strict-16 scheme (ncbi-phylum-2026-08-04-v1)
PHYLUM_KINGDOM <- c(
  Pseudomonadota    = "Bacteria",  Bacillota      = "Bacteria",
  Actinomycetota    = "Bacteria",
  Methanobacteriota = "Archaea",   Thermoproteota = "Archaea",
  Basidiomycota     = "Fungi",     Ascomycota     = "Fungi",
  Mucoromycota      = "Fungi",
  Streptophyta      = "Plantae",   Chlorophyta    = "Plantae",
  Arthropoda        = "Animalia",  Mollusca       = "Animalia",
  Nematoda          = "Animalia",
  Discosea          = "Protozoa",  Evosea         = "Protozoa",
  Heterolobosea     = "Protozoa"
)

stopifnot(all(BAR_ORDER %in% boot_df$kingdom),
          all(BAR_ORDER %in% marker_df$kingdom),
          all(eff_df$phylum %in% names(PHYLUM_KINGDOM)))

# ==================== Panel a ====================
pa_df <- boot_df |>
  mutate(mean_pct = mean * 100, lo_pct = ci_lo * 100, hi_pct = ci_hi * 100,
         kingdom = factor(kingdom, levels = rev(BAR_ORDER)))

pa_marker <- marker_df |>
  mutate(marker_pct = mean * 100,
         kingdom = factor(kingdom, levels = rev(BAR_ORDER)))

pa_expected <- EXPECTED |>
  mutate(kingdom = factor(kingdom, levels = rev(BAR_ORDER)))

# Compact visual key in the lower-right (Archaea = 1, Protists = 2).
# Numeric y is safe here because the discrete scale maps those levels to 1–2.
key_x0 <- 36.5
pa_key_bar <- data.frame(xmin = key_x0, xmax = key_x0 + 3.6,
                         ymin = 2.28, ymax = 2.58)
pa_key_dia <- data.frame(x = key_x0 + 1.8, y = 1.72)
pa_key_seg <- data.frame(x = key_x0, xend = key_x0 + 3.6, y = 1.18)

pa <- ggplot(pa_df, aes(y = kingdom, x = mean_pct)) +
  geom_col(aes(fill = as.character(kingdom)), width = 0.62,
           colour = "black", linewidth = 0.25, show.legend = FALSE) +
  geom_errorbarh(aes(xmin = lo_pct, xmax = hi_pct),
                 height = 0.18, linewidth = 0.45, colour = "black") +
  geom_segment(data = pa_expected,
               aes(y = as.numeric(kingdom) - 0.46,
                   yend = as.numeric(kingdom) - 0.46,
                   x = lo, xend = hi),
               inherit.aes = FALSE, colour = "grey35", linewidth = 1.8,
               lineend = "butt") +
  geom_point(data = pa_marker, aes(y = kingdom, x = marker_pct),
             shape = 23, size = 2.6, fill = "white", colour = "black",
             stroke = 0.5) +
  geom_rect(data = pa_key_bar,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            inherit.aes = FALSE, fill = KINGDOM_COLOURS[["Bacteria"]],
            colour = "black", linewidth = 0.25) +
  annotate("text", x = key_x0 + 4.3, y = 2.43, label = "fc-weighted",
           hjust = 0, vjust = 0.5, size = pt2mm(6), colour = "grey20",
           family = FONT) +
  geom_point(data = pa_key_dia, aes(x = x, y = y), inherit.aes = FALSE,
             shape = 23, size = 2.4, fill = "white", colour = "black",
             stroke = 0.5) +
  annotate("text", x = key_x0 + 4.3, y = 1.72, label = "marker panel",
           hjust = 0, vjust = 0.5, size = pt2mm(6), colour = "grey20",
           family = FONT) +
  geom_segment(data = pa_key_seg,
               aes(x = x, xend = xend, y = y, yend = y),
               inherit.aes = FALSE, colour = "grey35", linewidth = 1.8,
               lineend = "butt") +
  annotate("text", x = key_x0 + 4.3, y = 1.18, label = "literature range",
           hjust = 0, vjust = 0.5, size = pt2mm(6), colour = "grey20",
           family = FONT) +
  scale_fill_manual(values = KINGDOM_COLOURS) +
  scale_y_discrete(labels = DISPLAY[rev(BAR_ORDER)]) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.02)),
                     limits = c(0, 52)) +
  labs(x = "Mean fraction of matched soil lipid signal (%)", y = NULL) +
  theme_nature() +
  theme(
    text = element_text(family = FONT),
    axis.text.y = element_text(size = 7, family = FONT),
    panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.2)
  )

# ==================== Panel b ====================
eff <- eff_df |>
  mutate(kingdom = factor(PHYLUM_KINGDOM[phylum], levels = KINGDOM_ORDER),
         min_p = pmin(drought_p, climate_p))

xa <- max(abs(eff$drought_log2FC), 0.7) * 1.18
ya <- max(abs(eff$climate_log2FC), 0.4) * 1.28

# Call out the drought signal and the one nominal decline. No omnibus
# q < 0.05 exists; the unused FDR stroke encoding is not drawn.
label_specs <- tribble(
  ~phylum,          ~dx,   ~dy,   ~hjust, ~label_text,
  "Pseudomonadota", 0.08, -0.16,  0,      "Pseudomonadota\np = 0.005, q = 0.08",
  "Evosea",         0.07,  0.12,  0,      "Evosea\n(nominal)"
)

label_data <- label_specs |>
  left_join(eff, by = "phylum") |>
  mutate(x_end = drought_log2FC + dx,
         y_end = climate_log2FC + dy)

present_k <- KINGDOM_ORDER[KINGDOM_ORDER %in% levels(droplevels(eff$kingdom))]

stroke_key <- data.frame(
  x = -xa * 0.90,
  y = c(ya * 0.84, ya * 0.62),
  stroke = c(0.9, 0.4),
  lab = c("p < 0.05", "n.s.")
)

pb <- ggplot(eff, aes(x = drought_log2FC, y = climate_log2FC)) +
  geom_hline(yintercept = 0, colour = "#CCCCCC", linewidth = 0.4) +
  geom_vline(xintercept = 0, colour = "#CCCCCC", linewidth = 0.4) +
  geom_point(aes(fill = kingdom, size = mean_fraction),
             shape = 21, colour = "black", stroke = 0.4, alpha = 0.95) +
  geom_point(data = filter(eff, min_p < 0.05),
             aes(fill = kingdom, size = mean_fraction),
             shape = 21, colour = "black", stroke = 0.9, alpha = 0.95,
             show.legend = FALSE) +
  geom_segment(data = label_data,
               aes(x = drought_log2FC, y = climate_log2FC,
                   xend = x_end, yend = y_end),
               linewidth = 0.3, colour = "#666666", inherit.aes = FALSE) +
  geom_text(data = label_data,
            aes(x = x_end, y = y_end, label = label_text, hjust = hjust),
            vjust = 0.5, size = pt2mm(6), colour = "grey20",
            lineheight = 0.95, family = FONT, inherit.aes = FALSE) +
  geom_point(data = stroke_key[1, ], aes(x = x, y = y),
             inherit.aes = FALSE, shape = 21, size = 2.4,
             fill = "white", colour = "black", stroke = 0.9,
             show.legend = FALSE) +
  geom_point(data = stroke_key[2, ], aes(x = x, y = y),
             inherit.aes = FALSE, shape = 21, size = 2.4,
             fill = "white", colour = "black", stroke = 0.4,
             show.legend = FALSE) +
  annotate("text", x = -xa * 0.84, y = stroke_key$y[1],
           label = stroke_key$lab[1], hjust = 0, vjust = 0.5,
           size = pt2mm(6), colour = "grey20", family = FONT) +
  annotate("text", x = -xa * 0.84, y = stroke_key$y[2],
           label = stroke_key$lab[2], hjust = 0, vjust = 0.5,
           size = pt2mm(6), colour = "grey20", family = FONT) +
  annotate("text", x = xa * 0.96, y = -ya * 0.88,
           label = "drought-enriched", fontface = "italic",
           size = pt2mm(6), hjust = 1, vjust = 0.5, colour = "grey35",
           family = FONT) +
  annotate("text", x = -xa * 0.96, y = -ya * 0.88,
           label = "drought-depleted", fontface = "italic",
           size = pt2mm(6), hjust = 0, vjust = 0.5, colour = "grey35",
           family = FONT) +
  scale_fill_manual(values = KINGDOM_COLOURS, name = "Organism group",
                    limits = present_k, labels = DISPLAY[present_k]) +
  scale_size_area(max_size = 7, name = "Mean abundance",
                  breaks = c(0.01, 0.05, 0.12),
                  labels = c("1%", "5%", "12%")) +
  coord_cartesian(xlim = c(-xa, xa), ylim = c(-ya, ya)) +
  labs(
    x = expression("Drought effect ("*log[2]*" drought / no-drought)"),
    y = expression("Warming effect ("*log[2]*" future / ambient)")
  ) +
  theme_nature() +
  theme(
    text = element_text(family = FONT),
    legend.position = "right",
    legend.title = element_text(size = 6, face = "bold", family = FONT),
    legend.text = element_text(size = 6, family = FONT),
    legend.key.size = unit(3, "mm"),
    legend.spacing.y = unit(1, "mm"),
    legend.box.spacing = unit(2, "mm")
  ) +
  guides(
    fill = guide_legend(order = 1, override.aes = list(size = 3, stroke = 0.4)),
    size = guide_legend(order = 2, override.aes = list(fill = "#CCCCCC", stroke = 0.4))
  )

fig <- pa / pb +
  plot_layout(heights = c(0.62, 1.0)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = NC_TAG_PT, face = "bold", family = FONT))

outdir <- file.path(SCRIPT_DIR, "out")
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
w <- NC_DOUBLE * MM_TO_IN
h <- 150 * MM_TO_IN
ggsave(file.path(outdir, "Fig5_final.pdf"), fig,
       width = w, height = h, device = cairo_pdf, bg = "white")
ggsave(file.path(outdir, "Fig5_final.png"), fig,
       width = w, height = h, dpi = 600, bg = "white", device = ragg::agg_png)
cat("Fig5 final (R) ->", outdir, "\n")
