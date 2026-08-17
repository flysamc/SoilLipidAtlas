#!/usr/bin/env Rscript
# Fig 5 (final) -- ClimGrass: two-estimator community composition (a) +
# fingerprint-set treatment responses (b).
#
# R port of the locked final design (previously rendered by
# render_figure5_final.py): panel a shows the pooled fc-weighted estimate
# (bars, 95% reference-sample bootstrap CI) with the marker-panel estimator
# as diamonds and literature biomass ranges as grey segments. Panel b is the
# CLR fingerprint-set effect map (unchanged from fig6_climgrass_v2.R).
# Display labels follow the locked ecological_group policy
# (Viridiplantae / Protists); data keys remain Plantae / Protozoa.

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr)
  library(patchwork); library(scales)
})

.args <- commandArgs(trailingOnly = FALSE)
.f <- sub("--file=", "", .args[grep("--file=", .args)])
SCRIPT_DIR <- if (length(.f) > 0) dirname(normalizePath(.f)) else getwd()
source(file.path(SCRIPT_DIR, "soilmass_style.R"))
FONT <- "Arial"

WS <- dirname(SCRIPT_DIR)
boot_df   <- read.csv(file.path(WS, "composition_fcweighted_kingdom_ci.csv"))
marker_df <- read.csv(file.path(WS, "capped_marker_panel_test",
                                "kingdom_ci_marker_panel.csv"))
eff_df    <- read.csv(file.path(SCRIPT_DIR, "data", "phylum_effects.csv"))

# Locked display labels (ecological_group policy); dagger marks the
# ArchLips-only archaeal estimate.
DISPLAY <- c(Bacteria = "Bacteria", Fungi = "Fungi",
             Plantae = "Viridiplantae", Animalia = "Animalia",
             Protozoa = "Protists", Archaea = "Archaea†")

# Panel-a bar order (top to bottom), as in the locked final design.
BAR_ORDER <- c("Bacteria", "Fungi", "Plantae", "Animalia", "Protozoa", "Archaea")

# Literature biomass expectation ranges (%), identical to fig6_climgrass_v2.R.
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

# ==================== Panel a ====================
pa_df <- boot_df |>
  mutate(mean_pct = mean * 100, lo_pct = ci_lo * 100, hi_pct = ci_hi * 100,
         kingdom = factor(kingdom, levels = rev(BAR_ORDER)))

pa_marker <- marker_df |>
  mutate(marker_pct = mean * 100,
         kingdom = factor(kingdom, levels = rev(BAR_ORDER)))

pa_expected <- EXPECTED |>
  mutate(kingdom = factor(kingdom, levels = rev(BAR_ORDER)))

note <- paste0(
  "bars: lipid-derived community composition\n",
  "(fc-weighted, corrected; 95% CI, reference-sample\n",
  "bootstrap); diamonds: provenance of matched lipid\n",
  "signal (specific-marker panel; incl. plant necromass);\n",
  "grey segments: literature biomass ranges.\n",
  "† Archaea from ArchLips-validated ether lipids only\n",
  "(14 diagnostic markers); scale uncertain, no\n",
  "ether-lipid RIE standards")

pa <- ggplot(pa_df, aes(y = kingdom, x = mean_pct)) +
  geom_col(aes(fill = as.character(kingdom)), width = 0.62,
           colour = "black", linewidth = 0.25, show.legend = FALSE) +
  geom_errorbarh(aes(xmin = lo_pct, xmax = hi_pct),
                 height = 0.18, linewidth = 0.45, colour = "black") +
  geom_segment(data = pa_expected,
               aes(y = as.numeric(kingdom) - 0.46,
                   yend = as.numeric(kingdom) - 0.46,
                   x = lo, xend = hi),
               inherit.aes = FALSE, colour = "grey55", linewidth = 1.6,
               lineend = "butt") +
  geom_point(data = pa_marker, aes(y = kingdom, x = marker_pct),
             shape = 23, size = 2.6, fill = "white", colour = "black",
             stroke = 0.5) +
  scale_fill_manual(values = KINGDOM_COLOURS) +
  scale_y_discrete(labels = DISPLAY[rev(BAR_ORDER)]) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.04)),
                     limits = c(0, 55)) +
  labs(x = "Mean fraction of matched soil lipid signal (%)", y = NULL,
       title = "Organism-group composition (community estimate)") +
  annotate("text", x = 54.5, y = 1.05, label = note, hjust = 1, vjust = 0,
           size = pt2mm(5), colour = "grey40", family = FONT,
           lineheight = 1.05) +
  theme_nature() +
  theme(
    text = element_text(family = FONT),
    axis.text.y = element_text(size = 7, family = FONT),
    plot.title = element_text(size = 7, face = "bold", family = FONT,
                              hjust = 0),
    panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.2)
  )

# ==================== Panel b (unchanged from fig6_climgrass_v2.R,
# display labels conformed to the locked policy) ====================
eff <- eff_df |>
  mutate(kingdom = factor(PHYLUM_KINGDOM[phylum], levels = KINGDOM_ORDER),
         min_p = pmin(drought_p, climate_p),
         min_q = pmin(drought_q_bh, climate_q_bh))

xa <- max(abs(eff$drought_log2FC), 0.7) * 1.18
ya <- max(abs(eff$climate_log2FC), 0.4) * 1.28

label_specs <- tribble(
  ~phylum,            ~dx,    ~dy,    ~hjust, ~show_stats,
  "Pseudomonadota",   0.02,  -0.13,   0,      TRUE,
  "Evosea",           0.045,  0.09,   0,      TRUE,
  "Actinomycetota",  -0.045,  0.09,   1,      TRUE,
  "Streptophyta",     0.045,  0.06,   0,      FALSE,
  "Basidiomycota",    0.02,   0.14,   0,      FALSE
)

label_data <- label_specs |>
  left_join(eff, by = "phylum") |>
  mutate(x_end = drought_log2FC + dx,
         y_end = climate_log2FC + dy,
         label_text = ifelse(show_stats,
           sprintf("%s\np=%.3f; q=%.2f", phylum, drought_p, drought_q_bh),
           phylum))

lab_left   <- filter(label_data, hjust == 0)
lab_right  <- filter(label_data, hjust == 1)
present_k  <- KINGDOM_ORDER[KINGDOM_ORDER %in% levels(droplevels(eff$kingdom))]

pb <- ggplot(eff, aes(x = drought_log2FC, y = climate_log2FC)) +
  geom_hline(yintercept = 0, colour = "#CCCCCC", linewidth = 0.4) +
  geom_vline(xintercept = 0, colour = "#CCCCCC", linewidth = 0.4) +
  geom_point(aes(fill = kingdom, size = mean_fraction),
             shape = 21, colour = "black", stroke = 0.4, alpha = 0.95) +
  geom_point(data = filter(eff, min_p < 0.05),
             aes(fill = kingdom, size = mean_fraction),
             shape = 21, colour = "black", stroke = 0.9, alpha = 0.95,
             show.legend = FALSE) +
  geom_point(data = filter(eff, min_q < 0.05),
             aes(size = mean_fraction),
             shape = 21, fill = NA, colour = "black", stroke = 1.6,
             show.legend = FALSE) +
  geom_segment(data = label_data,
               aes(x = drought_log2FC, y = climate_log2FC,
                   xend = x_end, yend = y_end),
               linewidth = 0.3, colour = "#888888", inherit.aes = FALSE) +
  geom_text(data = lab_left,
            aes(x = x_end, y = y_end, label = label_text),
            hjust = 0, vjust = 0.5, size = pt2mm(5), colour = "grey30",
            lineheight = 0.9, family = FONT, inherit.aes = FALSE) +
  geom_text(data = lab_right,
            aes(x = x_end, y = y_end, label = label_text),
            hjust = 1, vjust = 0.5, size = pt2mm(5), colour = "grey30",
            lineheight = 0.9, family = FONT, inherit.aes = FALSE) +
  annotate("text", x = xa * 0.97, y = -ya * 0.97,
           label = "drought-enriched", fontface = "italic",
           size = pt2mm(5), hjust = 1, vjust = 0, colour = "#AAAAAA",
           family = FONT) +
  annotate("text", x = -xa * 0.97, y = -ya * 0.97,
           label = "drought-depleted", fontface = "italic",
           size = pt2mm(5), hjust = 0, vjust = 0, colour = "#AAAAAA",
           family = FONT) +
  scale_fill_manual(values = KINGDOM_COLOURS, name = "Organism group",
                    limits = present_k, labels = DISPLAY[present_k]) +
  scale_size_area(max_size = 7, name = "Mean abundance",
                  breaks = c(0.01, 0.05, 0.12),
                  labels = c("1%", "5%", "12%")) +
  coord_cartesian(xlim = c(-xa, xa), ylim = c(-ya, ya)) +
  labs(
    title = "Fingerprint-set treatment responses",
    x = expression("Drought effect  ("*log[2]*" drought / no-drought, fingerprint sets)"),
    y = expression("Warming effect  ("*log[2]*" future / ambient)"),
    caption = paste0(
      "Edge: thick = q < 0.05 (FDR), medium = p < 0.05, thin = n.s.  |  ",
      "CLR fingerprint-set statistic, exact stratified permutation test ",
      "(400 relabelings), n = 12 (6 vs 6).\n",
      "Pseudomonadota decline replicates the independent ClimGrass qSIP ",
      "prediction (pre-specified one-sided test, q = 0.005).\n",
      "† Archaea from ArchLips-validated ether lipids (14 diagnostic ",
      "markers); scale uncertain without ether-lipid RIE standards.")
  ) +
  theme_nature() +
  theme(
    text = element_text(family = FONT),
    plot.title = element_text(size = 7, face = "bold", family = FONT,
                              hjust = 0),
    legend.position = "right",
    legend.title = element_text(size = 6, face = "bold", family = FONT),
    legend.text = element_text(size = 5, family = FONT),
    legend.key.size = unit(3, "mm"),
    legend.spacing.y = unit(1, "mm"),
    legend.box.spacing = unit(2, "mm"),
    plot.caption = element_text(size = 5, colour = "#555555", hjust = 0,
                                family = FONT)
  ) +
  guides(
    fill = guide_legend(order = 1, override.aes = list(size = 3, stroke = 0.4)),
    size = guide_legend(order = 2, override.aes = list(fill = "#CCCCCC", stroke = 0.4))
  )

fig <- pa / pb +
  plot_layout(heights = c(0.72, 1.0)) +
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
