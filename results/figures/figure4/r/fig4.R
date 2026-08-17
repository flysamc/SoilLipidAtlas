#!/usr/bin/env Rscript
# Figure 4 – Distributed fingerprints: SIMPER lollipop (a) + lipid-class heatmap (b)
# + subsampling curve (c). Nature Communications: 183 mm wide, 5-7 pt, Helvetica.
#
# Organism-group display names are defined in this script (locked ecological_group
# vocabulary: Viridiplantae, Protists). CSV kingdom columns are not used.

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(scales)
  library(ggtext)
  library(jsonlite)
})

.args <- commandArgs(trailingOnly = FALSE)
.f <- sub("--file=", "", .args[grep("--file=", .args)])
SCRIPT_DIR <- if (length(.f) > 0) dirname(normalizePath(.f)) else getwd()
source(file.path(SCRIPT_DIR, "soilmass_style.R"))

# ── Data paths ────────────────────────────────────────────
data_dir <- file.path(SCRIPT_DIR, "data")

simper_df   <- read.csv(file.path(data_dir, "simper_curves.csv"))
class_df    <- read.csv(file.path(data_dir, "class_composition.csv"))
class_order <- fromJSON(file.path(data_dir, "class_order.json"))
sub_df      <- read.csv(file.path(data_dir, "subsampling_curve.csv"))

# Locked organism-group labels (ecological_group). Edit here, not in the CSVs.
PHYLUM_KINGDOM <- c(
  Actinomycetota    = "Bacteria",
  Bacillota         = "Bacteria",
  Pseudomonadota    = "Bacteria",
  Methanobacteriota = "Archaea",
  Thermoproteota    = "Archaea",
  Ascomycota        = "Fungi",
  Basidiomycota     = "Fungi",
  Mucoromycota      = "Fungi",
  Chlorophyta       = "Viridiplantae",
  Streptophyta      = "Viridiplantae",
  Arthropoda        = "Animalia",
  Mollusca          = "Animalia",
  Nematoda          = "Animalia",
  Discosea          = "Protists",
  Evosea            = "Protists",
  Heterolobosea     = "Protists"
)
stopifnot(identical(unname(sort(unique(PHYLUM_KINGDOM))),
                    unname(sort(KINGDOM_ORDER))))
pk_map <- data.frame(
  phylum  = names(PHYLUM_KINGDOM),
  kingdom = unname(PHYLUM_KINGDOM),
  stringsAsFactors = FALSE
)
pk <- PHYLUM_KINGDOM
class_df$kingdom <- unname(PHYLUM_KINGDOM[as.character(class_df$phylum)])

# ══════════════════════════════════════════════════════════
# Panel a: SIMPER lollipop  –  features to reach 50% separation
# ══════════════════════════════════════════════════════════

fifty_pct <- simper_df |>
  group_by(phylum) |>
  filter(cumulative_pct >= 50) |>
  slice_min(simper_rank, n = 1) |>
  ungroup() |>
  mutate(
    kingdom = pk[phylum],
    kingdom = factor(kingdom, levels = KINGDOM_ORDER)
  ) |>
  arrange(simper_rank) |>
  mutate(phylum = factor(phylum, levels = phylum))

phy_html <- sapply(levels(fifty_pct$phylum), function(p) {
  col <- KINGDOM_COLOURS[pk[p]]
  paste0("<b><span style='color:", col, "'>", p, "</span></b>")
})

pa <- ggplot(fifty_pct, aes(x = simper_rank, y = phylum, colour = kingdom)) +
  geom_segment(aes(x = 0, xend = simper_rank, y = phylum, yend = phylum),
               linewidth = 0.3, colour = "grey78") +
  geom_point(size = 1.6) +
  geom_text(aes(label = simper_rank), hjust = -0.35, size = pt2mm(5),
            show.legend = FALSE) +
  scale_colour_manual(
    values = KINGDOM_COLOURS, name = NULL,
    limits = KINGDOM_ORDER,
    labels = function(x) paste0("<b>", x, "</b>"),
    guide = guide_legend(nrow = 2, override.aes = list(size = 1.8))
  ) +
  scale_x_continuous(
    limits = c(0, max(fifty_pct$simper_rank) * 1.12),
    expand = expansion(mult = c(0, 0))
  ) +
  scale_y_discrete(labels = phy_html) +
  labs(x = "SIMPER features to 50% separation", y = NULL) +
  theme_nature() +
  theme(
    legend.position  = "bottom",
    legend.margin     = margin(1, 0, 0, 0, "mm"),
    legend.box.margin = margin(0, 0, 0, 0, "mm"),
    legend.key.size   = unit(2.5, "mm"),
    legend.spacing.x  = unit(1, "mm"),
    legend.spacing.y  = unit(0.5, "mm"),
    legend.text       = element_markdown(size = 6),
    axis.text.y       = element_markdown(size = 5.5),
    panel.grid.major.x = element_line(linewidth = 0.15, colour = "grey90")
  )

# ══════════════════════════════════════════════════════════
# Panel b: Lipid-class heatmap  (POS top, NEG bottom)
# ══════════════════════════════════════════════════════════

phylum_order_heatmap <- pk_map |>
  mutate(kingdom = factor(kingdom, levels = KINGDOM_ORDER)) |>
  arrange(kingdom, phylum) |>
  pull(phylum)

neg_only_phyla <- setdiff(unique(class_df$phylum), phylum_order_heatmap)
phylum_order_heatmap <- c(phylum_order_heatmap, neg_only_phyla)

kingdom_cuts <- pk_map |>
  mutate(kingdom = factor(kingdom, levels = KINGDOM_ORDER)) |>
  arrange(kingdom, phylum) |>
  count(kingdom, name = "n") |>
  mutate(cut = cumsum(n) + 0.5) |>
  filter(cut < length(phylum_order_heatmap) + 0.4) |>
  pull(cut)

html_phylum_labs <- sapply(phylum_order_heatmap, function(p) {
  k <- if (p %in% names(pk)) pk[p] else class_df$kingdom[class_df$phylum == p][1]
  col <- if (!is.na(k) && k %in% names(KINGDOM_COLOURS)) KINGDOM_COLOURS[k] else "grey40"
  paste0("<b><span style='color:", col, "'>", p, "</span></b>")
})

heatmap_data <- class_df |>
  filter(phylum %in% phylum_order_heatmap) |>
  mutate(
    phylum = factor(phylum, levels = phylum_order_heatmap),
    class  = factor(class, levels = rev(class_order)),
    mode   = factor(mode, levels = c("POS", "NEG"),
                    labels = c("Positive mode", "Negative mode"))
  )

pb <- ggplot(heatmap_data, aes(x = phylum, y = class, fill = fraction)) +
  geom_tile(width = 0.92, height = 0.92, colour = NA) +
  geom_vline(xintercept = kingdom_cuts, colour = "black", linewidth = 0.55) +
  facet_wrap(~ mode, ncol = 1, strip.position = "top") +
  scale_fill_gradient(
    low = "#f7fbff", high = "#08306b",
    limits = c(0, 0.55), oob = squish,
    name = "Fraction of\nbiomarkers",
    labels = percent_format(accuracy = 1)
  ) +
  scale_x_discrete(labels = html_phylum_labs, expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(x = NULL, y = NULL) +
  theme_nature() +
  theme(
    axis.text.x       = element_markdown(angle = 45, hjust = 1, size = 5),
    axis.text.y       = element_text(size = 5),
    strip.text        = element_text(size = 6, face = "bold", hjust = 0),
    legend.key.height = unit(3.5, "mm"),
    legend.key.width  = unit(2, "mm"),
    legend.title      = element_text(size = 5, lineheight = 1.1),
    legend.text       = element_text(size = 5),
    legend.position   = "right",
    panel.spacing     = unit(2.5, "mm"),
    panel.background  = element_rect(fill = "black", colour = "black", linewidth = 0.4)
  )

# ══════════════════════════════════════════════════════════
# Panel c: Subsampling curve  –  r vs % features retained
# ══════════════════════════════════════════════════════════

sub_plot <- sub_df |>
  mutate(
    pct = fraction * 100,
    mode = factor(
      ifelse(mode == "POS", "Positive", "Negative"),
      levels = c("Positive", "Negative")
    )
  )

pos5 <- sub_plot |> filter(pct == 5, mode == "Positive") |> slice(1)
neg5 <- sub_plot |> filter(pct == 5, mode == "Negative") |> slice(1)

pc <- ggplot(sub_plot, aes(x = pct, y = mean_r, colour = mode, fill = mode)) +
  geom_ribbon(aes(ymin = p5_r, ymax = p95_r), alpha = 0.15, colour = NA) +
  geom_line(linewidth = 0.5) +
  geom_point(size = 1.2, shape = 21, stroke = 0.3) +
  geom_hline(yintercept = 0.95, linetype = "dashed", linewidth = 0.25,
             colour = "grey50") +
  annotate("text", x = 97, y = 0.956, label = "r = 0.95",
           hjust = 1, size = pt2mm(5), colour = "grey50") +
  annotate("text", x = 22, y = 0.82,
           label = sprintf("5%% (%s)\nr = %.3f",
                           format(pos5$n_features, big.mark = ","), pos5$mean_r),
           size = pt2mm(5), hjust = 0, colour = MODE_COLOURS["Positive"],
           lineheight = 1.0) +
  annotate("segment", x = 5, xend = 20,
           y = pos5$mean_r, yend = 0.835, linewidth = 0.25,
           colour = MODE_COLOURS["Positive"]) +
  annotate("text", x = 22, y = 0.72,
           label = sprintf("5%% (%s)\nr = %.3f",
                           format(neg5$n_features, big.mark = ","), neg5$mean_r),
           size = pt2mm(5), hjust = 0, colour = MODE_COLOURS["Negative"],
           lineheight = 1.0) +
  annotate("segment", x = 5, xend = 20,
           y = neg5$mean_r, yend = 0.735, linewidth = 0.25,
           colour = MODE_COLOURS["Negative"]) +
  scale_colour_manual(values = MODE_COLOURS, name = NULL) +
  scale_fill_manual(values = MODE_COLOURS, name = NULL) +
  scale_x_continuous(
    breaks = c(0, 5, 10, 20, 30, 50, 70, 90, 100),
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  scale_y_continuous(
    limits = c(0.63, 1.005),
    breaks = seq(0.65, 1.0, 0.05),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(x = "Features retained (% of quality-filtered)",
       y = "Distance matrix correlation (Pearson r)") +
  theme_nature() +
  theme(
    legend.position        = "inside",
    legend.position.inside = c(0.88, 0.22),
    legend.justification   = c(1, 0),
    legend.background      = element_rect(fill = alpha("white", 0.85), colour = NA),
    legend.key.size        = unit(2.5, "mm"),
    legend.spacing.y       = unit(0.3, "mm"),
    panel.grid.major       = element_line(linewidth = 0.1, colour = "grey92"),
    axis.title.y           = element_text(margin = margin(r = 1.5, unit = "mm")),
    plot.margin            = margin(2, 2, 2, 1, "mm")
  ) +
  guides(colour = guide_legend(override.aes = list(size = 1.5)),
         fill   = "none")

# ══════════════════════════════════════════════════════════
# Compose: a (left, full height) | b (right top) / c (right bottom)
# ══════════════════════════════════════════════════════════

right_col <- pb / pc + plot_layout(heights = c(1.15, 0.85))

fig4 <- (pa | right_col) +
  plot_layout(widths = c(0.85, 1.15)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = NC_TAG_PT, face = "bold"))

# ── Save ─────────────────────────────────────────────────
outdir <- file.path(SCRIPT_DIR, "out")
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

w <- NC_DOUBLE * MM_TO_IN
h <- 160 * MM_TO_IN

ggsave(file.path(outdir, "Figure_4.pdf"), fig4,
       width = w, height = h, device = cairo_pdf, bg = "white")
ggsave(file.path(outdir, "Figure_4.png"), fig4,
       width = w, height = h, dpi = 600, bg = "white",
       device = ragg::agg_png)
ggsave(file.path(outdir, "Figure_4b.pdf"), pb,
       width = 120 * MM_TO_IN, height = 90 * MM_TO_IN,
       device = cairo_pdf, bg = "white")
ggsave(file.path(outdir, "Figure_4b.png"), pb,
       width = 120 * MM_TO_IN, height = 90 * MM_TO_IN, dpi = 600,
       bg = "white", device = ragg::agg_png)
cat("Figure 4 →", outdir, "\n")
