#!/usr/bin/env Rscript
# Fig 5 (revision proposal) -- ClimGrass: kingdom composition (a) + phylum effect map (b)
# Adapted from the manuscript's fig6_climgrass.R: identical style/layout;
# changes are (1) strict-16 unit scheme in the phylum lookup, (2) effect
# statistics from the CLR fingerprint-set permutation test, (3) a below-
# expected marker (new data has one), (4) updated curated labels + caption.

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr)
  library(patchwork); library(scales); library(ggtext)
})

.args <- commandArgs(trailingOnly = FALSE)
.f <- sub("--file=", "", .args[grep("--file=", .args)])
SCRIPT_DIR <- if (length(.f) > 0) dirname(normalizePath(.f)) else getwd()
source(file.path(SCRIPT_DIR, "soilmass_style.R"))
FONT <- "Arial"

comp_df <- read.csv(file.path(SCRIPT_DIR, "data", "kingdom_composition.csv"))
eff_df  <- read.csv(file.path(SCRIPT_DIR, "data", "phylum_effects.csv"))

TREAT_ORDER  <- c("Ambient_Control", "Ambient_Drought",
                  "Future_Control",  "Future_Drought")
TREAT_COLORS <- c(Ambient_Control = "#8EC5E5", Ambient_Drought = "#0072B2",
                  Future_Control  = "#F0B67F", Future_Drought  = "#D55E00")
TREAT_LABELS <- c(Ambient_Control = "Ambient, control",
                  Ambient_Drought = "Ambient, drought",
                  Future_Control  = "Future, control",
                  Future_Drought  = "Future, drought")

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
comp_long <- comp_df |>
  pivot_longer(cols = all_of(KINGDOM_ORDER),
               names_to = "kingdom", values_to = "frac") |>
  mutate(pct = frac * 100,
         kingdom = factor(kingdom, levels = KINGDOM_ORDER),
         treatment = factor(treatment, levels = TREAT_ORDER))

summ <- comp_long |>
  group_by(kingdom, treatment) |>
  summarise(mean_pct = mean(pct), sd_pct = sd(pct), .groups = "drop")

gmean <- comp_long |>
  group_by(kingdom) |>
  summarise(grand_mean = mean(pct), .groups = "drop")

marker_df <- summ |>
  mutate(top = mean_pct + sd_pct) |>
  group_by(kingdom) |>
  summarise(y_top = max(top, na.rm = TRUE), .groups = "drop") |>
  left_join(gmean, by = "kingdom") |>
  left_join(EXPECTED, by = "kingdom") |>
  mutate(within = grand_mean >= lo & grand_mean <= hi,
         above  = grand_mean > hi,
         below  = grand_mean < lo,
         marker_y = pmax(y_top, hi) + 2.0,
         kingdom = factor(kingdom, levels = KINGDOM_ORDER))

expect_rects <- EXPECTED |>
  mutate(x = match(kingdom, KINGDOM_ORDER),
         xmin = x - 0.46, xmax = x + 0.46)

ymax_a <- max(c(summ$mean_pct + summ$sd_pct, EXPECTED$hi), na.rm = TRUE) * 1.18

x_html <- sapply(KINGDOM_ORDER, function(k) {
  lab <- if (k == "Archaea") "Archaea&#8224;" else k
  paste0("<span style='color:", KINGDOM_COLOURS[k],
         ";font-weight:bold'>", lab, "</span>")
})

key_y1 <- 44; key_y2 <- 40; key_y3 <- 36; key_y4 <- 32
key_x  <- 4.75

pa <- ggplot(summ, aes(x = kingdom, y = mean_pct, fill = treatment)) +
  geom_rect(data = expect_rects, inherit.aes = FALSE,
            aes(xmin = xmin, xmax = xmax, ymin = lo, ymax = hi),
            fill = "#DDDDDD", alpha = 0.55) +
  geom_col(position = position_dodge(width = 0.82), width = 0.19,
           colour = "black", linewidth = 0.2) +
  geom_errorbar(aes(ymin = mean_pct - sd_pct, ymax = mean_pct + sd_pct),
                position = position_dodge(width = 0.82), width = 0.12,
                linewidth = 0.4, colour = "black") +
  geom_point(data = filter(marker_df, within),
             aes(x = kingdom, y = marker_y), inherit.aes = FALSE,
             shape = 16, colour = "#3F8E3F", size = 2.2) +
  geom_point(data = filter(marker_df, above),
             aes(x = kingdom, y = marker_y), inherit.aes = FALSE,
             shape = 17, colour = "#C46210", size = 2.8) +
  geom_point(data = filter(marker_df, below),
             aes(x = kingdom, y = marker_y), inherit.aes = FALSE,
             shape = 25, colour = "#666666", fill = "#666666", size = 2.4) +
  scale_fill_manual(values = TREAT_COLORS, labels = TREAT_LABELS,
                    name = "Treatment  (climate x water)") +
  scale_x_discrete(labels = x_html) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.02)),
                     limits = c(0, ymax_a)) +
  labs(y = "Soil lipid signal attributed (%)", x = NULL) +
  annotate("rect", xmin = key_x - 0.18, xmax = key_x + 0.18,
           ymin = key_y1 - 1.6, ymax = key_y1 + 1.6,
           fill = "#DDDDDD", colour = NA) +
  annotate("text", x = key_x + 0.30, y = key_y1,
           label = "expected range", size = pt2mm(5.5), family = FONT,
           hjust = 0, vjust = 0.5, colour = "grey40") +
  annotate("point", x = key_x - 0.10, y = key_y2,
           shape = 16, colour = "#3F8E3F", size = 2) +
  annotate("text", x = key_x + 0.08, y = key_y2,
           label = "within", size = pt2mm(5), family = FONT,
           hjust = 0, vjust = 0.5, colour = "grey40") +
  annotate("point", x = key_x - 0.10, y = key_y3,
           shape = 17, colour = "#C46210", size = 2.4) +
  annotate("text", x = key_x + 0.08, y = key_y3,
           label = "above", size = pt2mm(5), family = FONT,
           hjust = 0, vjust = 0.5, colour = "grey40") +
  annotate("point", x = key_x - 0.10, y = key_y4,
           shape = 25, colour = "#666666", fill = "#666666", size = 2.2) +
  annotate("text", x = key_x + 0.08, y = key_y4,
           label = "below", size = pt2mm(5), family = FONT,
           hjust = 0, vjust = 0.5, colour = "grey40") +
  theme_nature() +
  theme(
    text = element_text(family = FONT),
    axis.text.x = element_markdown(size = 6, family = FONT),
    legend.position = "inside",
    legend.position.inside = c(0.99, 0.99),
    legend.justification = c(1, 1),
    legend.background = element_rect(fill = alpha("white", 0.92), colour = NA),
    legend.title = element_text(size = 5.5, face = "bold", family = FONT),
    legend.text = element_text(size = 5, family = FONT),
    legend.key.size = unit(3, "mm"),
    legend.key.height = unit(3, "mm"),
    legend.spacing.y = unit(0.5, "mm"),
    panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.2)
  ) +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE))

# ==================== Panel b ====================
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
lab_centre <- filter(label_data, hjust == 0.5)
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
  geom_text(data = lab_centre,
            aes(x = x_end, y = y_end, label = label_text),
            hjust = 0.5, vjust = 0.5, size = pt2mm(5), colour = "grey30",
            lineheight = 0.9, family = FONT, inherit.aes = FALSE) +
  annotate("text", x = xa * 0.97, y = -ya * 0.97,
           label = "drought-enriched", fontface = "italic",
           size = pt2mm(5), hjust = 1, vjust = 0, colour = "#AAAAAA", family = FONT) +
  annotate("text", x = -xa * 0.97, y = -ya * 0.97,
           label = "drought-depleted", fontface = "italic",
           size = pt2mm(5), hjust = 0, vjust = 0, colour = "#AAAAAA", family = FONT) +
  scale_fill_manual(values = KINGDOM_COLOURS, name = "Organism group",
                    limits = present_k) +
  scale_size_area(max_size = 7, name = "Mean abundance",
                  breaks = c(0.01, 0.05, 0.12),
                  labels = c("1%", "5%", "12%")) +
  coord_cartesian(xlim = c(-xa, xa), ylim = c(-ya, ya)) +
  labs(
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
  plot_layout(heights = c(0.8, 1.0)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = NC_TAG_PT, face = "bold", family = FONT))

outdir <- file.path(SCRIPT_DIR, "out")
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
w <- NC_DOUBLE * MM_TO_IN
h <- 145 * MM_TO_IN
ggsave(file.path(outdir, "Fig5_climgrass_v2.pdf"), fig,
       width = w, height = h, device = cairo_pdf, bg = "white")
ggsave(file.path(outdir, "Fig5_climgrass_v2.png"), fig,
       width = w, height = h, dpi = 600, bg = "white", device = ragg::agg_png)
cat("Fig5 v2 ->", outdir, "\n")
