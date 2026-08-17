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
CLADE_LABEL <- c("Euryarchaeota_sg" = "Euryarchaeota s.l.")

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

# -- Compose ------------------------------------------------------
ed3 <- (p_pos | p_neg) +
  plot_layout(widths = c(1, 1)) +
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

save_all(ed3, file.path(out_dir, "Supplementary_Fig2_clade"), NC_DOUBLE, 78)
cat("Supplementary Fig 2 ->", out_dir, "\n")
