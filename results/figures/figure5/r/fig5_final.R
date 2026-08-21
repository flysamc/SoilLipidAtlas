#!/usr/bin/env Rscript
# Fig 5 — ClimGrass: biomass-converted community composition (a) + phylum responses (b)
# Public polished version: Panel a bars are biomass estimates (lipid signal / lipid content, Table S6)
# No internal datafile table in legend; methods pointer only.

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
eff_df    <- read.csv(file.path(data_dir, "phylum_effects.csv"))

DISPLAY <- c(Bacteria = "Bacteria", Fungi = "Fungi",
             Plantae = "Viridiplantae", Animalia = "Animalia",
             Protozoa = "Protists", Archaea = "Archaea")
BAR_ORDER <- c("Bacteria", "Fungi", "Plantae", "Animalia", "Protozoa", "Archaea")

# Lipid content as conversion factor (% dry mass, Table S6 v2)
LIPID_CONTENT <- c(Bacteria=7, Fungi=10, Plantae=2.5, Animalia=18, Protozoa=10, Archaea=6)

EXPECTED <- data.frame(
  kingdom = KINGDOM_ORDER,
  lo = c(35, 2, 15, 15, 1, 0),
  hi = c(50, 8, 30, 30, 5, 1),
  stringsAsFactors = FALSE
)

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
          all(eff_df$phylum %in% names(PHYLUM_KINGDOM)))

# Biomass conversion: biomass = (lipid_signal / lipid_content) / sum, renormalised
# Use single sum_mean denominator for mean and CI to preserve ordering (approx. from CI bounds)
bio <- boot_df %>%
  mutate(content = LIPID_CONTENT[kingdom],
         mean_raw = mean / content,
         lo_raw = ci_lo / content,
         hi_raw = ci_hi / content)
sum_mean <- sum(bio$mean_raw)
pa_df <- bio %>%
  mutate(mean_pct = mean_raw / sum_mean * 100,
         lo_pct = lo_raw / sum_mean * 100,
         hi_pct = hi_raw / sum_mean * 100,
         kingdom = factor(kingdom, levels = rev(BAR_ORDER)))

pa_expected <- EXPECTED %>%
  mutate(kingdom = factor(kingdom, levels = rev(BAR_ORDER)))

# No marker diamonds in polished biomass version — keep clean
key_x0 <- 36.5
pa_key_bar <- data.frame(xmin = key_x0, xmax = key_x0 + 3.6,
                         ymin = 1.72, ymax = 2.02)
pa_key_seg <- data.frame(x = key_x0, xend = key_x0 + 3.6, y = 1.18)

pa <- ggplot(pa_df, aes(y = kingdom, x = mean_pct)) +
  geom_col(aes(fill = as.character(kingdom)), width = 0.62,
           colour = "black", linewidth = 0.25, show.legend = FALSE) +
  geom_errorbar(aes(xmin = lo_pct, xmax = hi_pct),
                width = 0.18, linewidth = 0.45, colour = "black", orientation = "y") +
  geom_segment(data = pa_expected,
               aes(y = as.numeric(kingdom) - 0.46,
                   yend = as.numeric(kingdom) - 0.46,
                   x = lo, xend = hi),
               inherit.aes = FALSE, colour = "grey35", linewidth = 1.8,
               lineend = "butt") +
  geom_rect(data = pa_key_bar,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            inherit.aes = FALSE, fill = KINGDOM_COLOURS[["Bacteria"]],
            colour = "black", linewidth = 0.25) +
  annotate("text", x = key_x0 + 4.3, y = 1.87, label = "biomass estimate",
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
  labs(x = "Estimated biomass (%)", y = NULL) +
  theme_nature() +
  theme(
    text = element_text(family = FONT),
    axis.text.y = element_text(size = 7, family = FONT),
    panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.2)
  )

# Panel b unchanged
eff <- eff_df %>%
  mutate(kingdom = factor(PHYLUM_KINGDOM[phylum], levels = KINGDOM_ORDER),
         min_p = pmin(drought_p, climate_p))
xa <- max(abs(eff$drought_log2FC), 0.7) * 1.18
ya <- max(abs(eff$climate_log2FC), 0.4) * 1.28
label_specs <- tribble(
  ~phylum,          ~dx,   ~dy,   ~hjust, ~label_text,
  "Pseudomonadota", 0.08, -0.16,  0,      "Pseudomonadota\np = 0.005, FDR q = 0.08",
  "Evosea",        -0.09,  0.13,  1,      "Evosea\nnominal p = 0.050",
  "Actinomycetota", 0.08,  0.06,  0,      "Actinomycetota",
  "Bacillota",      0.08, -0.08,  0,      "Bacillota",
  "Ascomycota",    -0.12,  0.12,  1,      "Ascomycota",
  "Basidiomycota",  0.09,  0.09,  0,      "Basidiomycota",
  "Mucoromycota",   0.07,  0.07,  0,      "Mucoromycota"
)
label_data <- label_specs %>%
  left_join(eff, by = "phylum") %>%
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
cat("Fig5 biomass (public) ->", outdir, "\n")
# print biomass values for verification
print(pa_df[, c("kingdom","mean_pct","lo_pct","hi_pct")])
