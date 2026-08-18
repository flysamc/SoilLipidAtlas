#!/usr/bin/env Rscript
# Supplementary Fig. 2  -  Clade-conserved vs clade-exclusive features.
# Log-scale lollipop per clade: shared (filled) and exclusive (open) markers, so
# the orders-of-magnitude gap (large shared pool, negligible exclusive) is the
# visible message. Positive and negative mode shown as two panels.
# Port of figures_v5 panels/clade_features.py + compose/ed_clade.py.
# Nature Communications: 183 mm wide, all text 5-7 pt, Helvetica.

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  library(ggtext)
  library(scales)
})

.args <- commandArgs(trailingOnly = FALSE)
.f <- sub("--file=", "", .args[grep("--file=", .args)])
SCRIPT_DIR <- if (length(.f) > 0) dirname(normalizePath(.f)) else getwd()
source(file.path(SCRIPT_DIR, "soilmass_style.R"))

# -- Data ---------------------------------------------------------
data_dir <- file.path(SCRIPT_DIR, "data", "supp_clade")
out_dir  <- file.path(SCRIPT_DIR, "out")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Display names for clade keys that are opaque in the frozen data.
# 'Euryarchaeota_sg' = {Euryarchaeota + Methanobacteriota} = Euryarchaeota sensu lato.
CLADE_LABEL <- c("Euryarchaeota_sg" = "Euryarchaeota s.l.",
                 "All Plantae"      = "All Viridiplantae")  # locked display policy

clade <- read.csv(file.path(data_dir, "clade_counts.csv"), check.names = FALSE)
clade <- clade[order(clade$pos_shared), ]          # ascending: smallest at bottom
clade$y <- seq_len(nrow(clade)) - 1                 # y = 0 (bottom) .. n-1 (top)
clade$display <- ifelse(clade$clade %in% names(CLADE_LABEL),
                        CLADE_LABEL[clade$clade], clade$clade)
clade$col <- KINGDOM_COLOURS[clade$kingdom]

LEG_LEVELS  <- c("shared across clade", "clade-exclusive")
LINE_GREY   <- "#E6E6E6"   # COLOR_GRID
LIGHT_GREY  <- "#666666"   # COLOR_TEXT_LIGHT

# -- One lollipop panel -------------------------------------------
clade_panel <- function(shared_col, excl_col, title, show_y, show_legend) {
  pd <- clade
  pd$sh  <- pmax(pd[[shared_col]], 0.6)             # floor zeros onto the log axis
  pd$ex  <- pmax(pd[[excl_col]],   0.6)
  pd$shv <- pd[[shared_col]]
  pd$exv <- pd[[excl_col]]
  xmax   <- max(pd[[shared_col]]) * 6

  html_labs <- sprintf("<span style='color:%s'>%s</span>", pd$col, pd$display)

  p <- ggplot(pd) +
    geom_segment(aes(x = ex, xend = sh, y = y, yend = y),
                 colour = LINE_GREY, linewidth = 0.6) +
    geom_point(aes(x = sh, y = y, shape = "shared across clade"),
               fill = pd$col, colour = "white", size = 2.6, stroke = 0.4) +
    geom_point(aes(x = ex, y = y, shape = "clade-exclusive"),
               fill = "white", colour = pd$col, size = 2.3, stroke = 0.8) +
    geom_text(aes(x = sh * 1.25, y = y, label = comma(shv)),
              hjust = 0, colour = pd$col, size = pt2mm(5)) +
    geom_text(aes(x = ex * 0.8, y = y, label = as.character(exv)),
              hjust = 1, colour = LIGHT_GREY, size = pt2mm(5)) +
    scale_shape_manual(values = c("shared across clade" = 21, "clade-exclusive" = 21),
                       breaks = LEG_LEVELS, name = NULL) +
    scale_x_log10(breaks = 10^(0:5),
                  labels = scales::trans_format("log10", scales::math_format(10^.x))) +
    scale_y_continuous(breaks = pd$y, labels = html_labs,
                       expand = expansion(add = 0.6)) +
    coord_cartesian(xlim = c(0.4, xmax), clip = "off") +
    labs(x = "Features (log scale)", y = NULL, title = title) +
    theme_nature() +
    theme(
      axis.ticks.y = element_blank(),
      plot.title   = element_text(size = 7, face = "bold", hjust = 0),
      panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.15)
    ) +
    guides(shape = guide_legend(override.aes = list(
      fill = c("#777777", "white"), colour = c("white", "#777777"),
      stroke = c(0.4, 0.8), size = c(2.6, 2.3))))

  if (show_y) {
    p <- p + theme(axis.text.y = element_markdown(size = 6))
  } else {
    p <- p + theme(axis.text.y = element_blank())
  }
  if (show_legend) {
    p <- p + theme(legend.position = "inside",
                   legend.position.inside = c(0.99, 0.02),
                   legend.justification = c(1, 0),
                   legend.text = element_text(size = 5.5),
                   legend.key.size = unit(3, "mm"),
                   legend.background = element_blank())
  } else {
    p <- p + theme(legend.position = "none")
  }
  p
}

p_pos <- clade_panel("pos_shared", "pos_exclusive", "Positive mode",
                     show_y = TRUE,  show_legend = FALSE)
p_neg <- clade_panel("neg_shared", "neg_exclusive", "Negative mode",
                     show_y = FALSE, show_legend = TRUE)

# -- Pairwise phylum-level sharing heatmaps (panels c, d) ---------
pair <- read.csv(file.path(data_dir, "pairwise_sharing.csv"), check.names = FALSE)

# Locked display groups (Viridiplantae/Protists) -> legacy palette keys.
GROUP_PALETTE_KEY <- c(Viridiplantae = "Plantae", Protists = "Protozoa")

heat_panel <- function(mode_key, title, show_legend) {
  pd <- pair[pair$mode == mode_key, ]
  ord    <- unique(pd[order(pd$rank_1), "phylum_1"])
  ord    <- unique(c(ord, pd[order(pd$rank_2), "phylum_2"]))
  grp_of <- setNames(c(pd$group_1, pd$group_2), c(pd$phylum_1, pd$phylum_2))
  # lower triangle: x = earlier-ranked phylum (no last), y = later-ranked (no first)
  x_ord <- ord[-length(ord)]
  y_ord <- ord[-1]
  pd$xf <- factor(pd$phylum_1, levels = x_ord)
  pd$yf <- factor(pd$phylum_2, levels = rev(y_ord))
  pd$pct <- pd$frac_of_smaller * 100
  grp_key <- function(g) ifelse(g %in% names(GROUP_PALETTE_KEY), GROUP_PALETTE_KEY[g], g)
  lab_html <- function(phy) sprintf("<span style='color:%s'>%s</span>",
                                    KINGDOM_COLOURS[grp_key(grp_of[phy])], phy)
  x_labs <- lab_html(x_ord)
  y_labs <- lab_html(rev(y_ord))

  ggplot(pd, aes(x = xf, y = yf, fill = pct)) +
    geom_tile(colour = "white", linewidth = 0.25) +
    geom_text(aes(label = sprintf("%.0f", pct)),
              size = pt2mm(4.4), colour = "black") +
    scale_fill_gradient(low = "#FDF6EC", high = "#D9820F",
                        limits = c(0, 100), name = "Shared (%)") +
    scale_x_discrete(labels = x_labs, drop = FALSE) +
    scale_y_discrete(labels = y_labs, drop = FALSE) +
    coord_fixed(clip = "off") +
    labs(x = NULL, y = NULL, title = title) +
    theme_nature() +
    theme(
      axis.text.x  = element_markdown(size = 5.4, angle = 45, hjust = 1),
      axis.text.y  = element_markdown(size = 5.4),
      axis.ticks   = element_blank(),
      panel.grid   = element_blank(),
      plot.title   = element_text(size = 7, face = "bold", hjust = 0),
      legend.position = if (show_legend) "right" else "none",
      legend.title = element_text(size = 5.5),
      legend.text  = element_text(size = 5),
      legend.key.height = unit(6, "mm"),
      legend.key.width  = unit(2.5, "mm")
    )
}

p_heat_pos <- heat_panel("POS", "Positive mode", show_legend = FALSE)
p_heat_neg <- heat_panel("NEG", "Negative mode", show_legend = TRUE)

# -- Compose ------------------------------------------------------
ed3 <- ((p_pos | p_neg) + plot_layout(widths = c(1, 1))) /
  ((p_heat_pos | p_heat_neg) + plot_layout(widths = c(1, 1.12))) +
  plot_layout(heights = c(78, 100)) +
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

save_all(ed3, file.path(out_dir, "Supplementary_Fig2_clade"), NC_DOUBLE, 185)
cat("Supplementary Fig 2 ->", out_dir, "\n")
