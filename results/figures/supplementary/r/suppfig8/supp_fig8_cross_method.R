#!/usr/bin/env Rscript
# Supplementary Fig. 8  -  Cross-method fingerprint validation (4 panels).
#   a  Dendrogram-reconstruction Mantel r vs full 44k-feature dendrogram for
#      SIMPER, SCBD, CAP, L1 across per-phylum top-K, with a random-K null band.
#   b  Leave-one-out classification accuracy by method and K, full-substrate baseline.
#   c  Fraction of 769 ClimGrass verified soil features captured by each method's top-K union.
#   d  All-four-method consensus features per phylum at K=500 (bars coloured by kingdom;
#      light = at least 3 methods, dark = all 4).
# Port of method_validation_analysis20/06_figure/fig_cross_method_validation.py.
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
data_dir <- file.path(SCRIPT_DIR, "data", "supp_cross_method")
out_dir  <- file.path(SCRIPT_DIR, "out")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# -- Method ordering, display names, colours, markers (from the .py) --
METHOD_ORDER   <- c("simper", "scbd", "cap", "stability_l1")
METHOD_DISPLAY <- c(simper = "SIMPER", scbd = "SCBD", cap = "CAP",
                    stability_l1 = "L1 stability")
# Distinct colours, deliberately NOT overlapping the kingdom palette
METHOD_COLORS  <- c(simper = "#1F4E79", scbd = "#7E57C2",
                    cap = "#D55E00", stability_l1 = "#117733")
# pch (solid, coloured by line colour): o=16 (circle), s=15 (square),
# ^=17 (triangle), D=18 (diamond) -- matches the .py o/s/^/D markers
METHOD_SHAPES  <- c(simper = 16, scbd = 15, cap = 17, stability_l1 = 18)

K_BREAKS <- c(100, 250, 500, 1000, 2500)
GREY_DASH <- "#666666"

disp_levels <- unname(METHOD_DISPLAY[METHOD_ORDER])
col_disp    <- setNames(unname(METHOD_COLORS[METHOD_ORDER]), disp_levels)
shp_disp    <- setNames(unname(METHOD_SHAPES[METHOD_ORDER]), disp_levels)

# Shared method line+point layers (colour + shape only, so the two merge into a
# single legend key per method).
method_line_layers <- function(df, yvar) {
  df <- df[df$method %in% METHOD_ORDER, ]
  df$mlabel <- factor(METHOD_DISPLAY[df$method], levels = disp_levels)
  df <- df[order(df$mlabel, df$K), ]
  list(
    geom_line(data = df, aes(K, .data[[yvar]], colour = mlabel),
              linewidth = 0.5),
    geom_point(data = df, aes(K, .data[[yvar]], colour = mlabel, shape = mlabel),
               size = 1.5)
  )
}

# Method colour + shape scales (identical name/breaks -> one merged legend).
method_scales <- function() {
  list(
    scale_colour_manual(values = col_disp, breaks = disp_levels, name = NULL),
    scale_shape_manual(values = shp_disp, breaks = disp_levels, name = NULL),
    scale_x_log10(breaks = K_BREAKS, labels = as.character(K_BREAKS))
  )
}

inside_legend <- function(x, y, just) {
  theme(
    legend.position        = "inside",
    legend.position.inside = c(x, y),
    legend.justification   = just,
    legend.text            = element_text(size = 5),
    legend.key.size        = unit(2.6, "mm"),
    legend.spacing.y       = unit(0.3, "mm"),
    legend.background      = element_blank(),
    legend.key             = element_blank()
  )
}

# ================================================================
#  Panel a: Mantel r vs K, four methods + random-K null band
# ================================================================
mantel <- read.csv(file.path(data_dir, "panel_a_mantel_curves.csv"))
null_a <- read.csv(file.path(data_dir, "panel_a_random_null.csv"))

pa <- ggplot() +
  # null band (q05-q95) -> own one-level fill legend (methods no longer use fill)
  geom_ribbon(data = null_a,
              aes(x = K_per_phylum, ymin = null_q05, ymax = null_q95,
                  fill = "Random K-features (q05-q95)"),
              alpha = 0.18) +
  # null mean -> own one-level linetype legend
  geom_line(data = null_a,
            aes(x = K_per_phylum, y = null_mean,
                linetype = "Random K-features (mean)"),
            colour = GREY_DASH, linewidth = 0.4) +
  geom_hline(yintercept = 1.0, colour = "grey60", linewidth = 0.25,
             linetype = "dotted") +
  method_line_layers(mantel, "mantel_r_vs_full") +
  method_scales() +
  scale_fill_manual(values = c("Random K-features (q05-q95)" = "grey60"),
                    name = NULL) +
  scale_linetype_manual(values = c("Random K-features (mean)" = "dashed"),
                        name = NULL) +
  coord_cartesian(ylim = c(0.6, 1.01)) +
  labs(x = "Top-K features per phylum",
       y = "Mantel r vs full 44k-feature dendrogram",
       title = "Dendrogram reconstruction") +
  theme_nature() +
  theme(plot.title = element_text(size = 7, face = "bold", hjust = 0)) +
  inside_legend(0.99, 0.03, c(1, 0)) +
  guides(
    colour   = guide_legend(order = 1),
    shape    = guide_legend(order = 1),
    fill     = guide_legend(order = 2, override.aes = list(alpha = 0.18)),
    linetype = guide_legend(order = 3,
                            override.aes = list(colour = GREY_DASH))
  )

# ================================================================
#  Panel b: LOO classification accuracy vs K + full-substrate baseline
# ================================================================
loo <- read.csv(file.path(data_dir, "panel_b_loo_accuracy.csv"))
baseline <- loo$full_substrate_baseline[1]
loo$acc_pct <- loo$accuracy * 100

baseline_lab <- sprintf("Full 44k-substrate baseline (%.1f%%)", baseline * 100)

pb <- ggplot() +
  geom_hline(aes(yintercept = baseline * 100, linetype = baseline_lab),
             colour = "grey45", linewidth = 0.4) +
  method_line_layers(loo, "acc_pct") +
  method_scales() +
  scale_linetype_manual(values = setNames("dashed", baseline_lab), name = NULL) +
  coord_cartesian(ylim = c(28, 85)) +
  labs(x = "Top-K features per phylum",
       y = "Leave-one-out classification accuracy (%)",
       title = "Sample classification") +
  theme_nature() +
  theme(plot.title = element_text(size = 7, face = "bold", hjust = 0)) +
  inside_legend(0.99, 0.03, c(1, 0)) +
  guides(
    colour   = guide_legend(order = 1),
    shape    = guide_legend(order = 1),
    linetype = guide_legend(order = 2, override.aes = list(colour = "grey45"))
  )

# ================================================================
#  Panel c: Fraction of 769 ClimGrass features captured in top-K
# ================================================================
cg <- read.csv(file.path(data_dir, "panel_c_climgrass_overlap.csv"))
n_cg <- cg$n_climgrass_features[1]
cg$frac_pct <- cg$frac_of_climgrass_in_topK * 100

pc <- ggplot() +
  method_line_layers(cg, "frac_pct") +
  method_scales() +
  coord_cartesian(ylim = c(0, 100)) +
  labs(x = "Top-K features per phylum",
       y = sprintf("Fraction of ClimGrass-matched features\nin top-K (n=%d, %%)", n_cg),
       title = "Soil overlap (ClimGrass)") +
  theme_nature() +
  theme(plot.title = element_text(size = 7, face = "bold", hjust = 0)) +
  inside_legend(0.02, 0.97, c(0, 1)) +
  guides(colour = guide_legend(order = 1), shape = guide_legend(order = 1))

# ================================================================
#  Panel d: All-four-method consensus features per phylum (K=500)
# ================================================================
cons <- read.csv(file.path(data_dir, "panel_d_consensus_per_phylum.csv"))
# Sort all4 asc then geq3 asc (matches .py) -> largest all4 at TOP of horizontal bars
cons <- cons[order(cons$all4, cons$geq3), ]
cons$phylum <- factor(cons$phylum, levels = cons$phylum)
cons$kcol   <- KINGDOM_COLOURS[cons$kingdom]

# x-limit: max(geq3)+20, min 50 (matches .py)
x_max <- max(max(cons$geq3) + 20, 50)

# kingdom legend order = first appearance from bottom to top (dict.fromkeys on the
# already-sorted df), mirroring the .py "used" list
used_kingdoms <- unique(cons$kingdom)

pd <- ggplot(cons, aes(y = phylum)) +
  # background bar: >=3 methods (lighter)
  geom_col(aes(x = geq3, fill = kingdom), width = 0.78, alpha = 0.30) +
  # foreground bar: all 4 methods (full opacity)
  geom_col(aes(x = all4, fill = kingdom), width = 0.78, alpha = 1.0) +
  # count labels at the right edge of each all4 bar (only where all4 > 0)
  geom_text(data = cons[cons$all4 > 0, ],
            aes(x = all4 + 2, y = phylum, label = all4, colour = kingdom),
            hjust = 0, fontface = "bold", size = pt2mm(6)) +
  scale_fill_manual(values = KINGDOM_COLOURS, breaks = used_kingdoms,
                    name = "Organism group") +
  scale_colour_manual(values = KINGDOM_COLOURS, guide = "none") +
  scale_x_continuous(limits = c(0, x_max), expand = expansion(mult = c(0, 0.02))) +
  labs(x = "Features in top-500 of >=3 / all 4 methods", y = NULL,
       title = "Consensus features per phylum") +
  theme_nature() +
  theme(
    axis.text.y = element_text(size = 5.5),
    plot.title  = element_text(size = 7, face = "bold", hjust = 0)
  ) +
  inside_legend(0.99, 0.03, c(1, 0)) +
  theme(legend.title = element_text(size = 6, face = "bold")) +
  guides(fill = guide_legend(override.aes = list(alpha = 1)))

# ================================================================
#  Compose 2x2 with patchwork
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

save_all(fig, file.path(out_dir, "Supplementary_Fig8_cross_method"), NC_DOUBLE, 162)
cat("Supplementary Fig 8 ->", out_dir, "\n")
