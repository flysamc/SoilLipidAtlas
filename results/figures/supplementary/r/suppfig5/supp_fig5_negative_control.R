#!/usr/bin/env Rscript
# Supplementary Fig. 5  -  Negative control: pure isolates fed through the decomposition.
#   a  Kingdom confusion matrix, raw (uncorrected) method. Archaea acts as a
#      false-positive sink (Bacteria-to-Archaea 33%, Plantae-to-Archaea 21%).
#   b  Corrected pipeline (IS normalisation + 0.20 RIE floor + ArchLips archaeal
#      masking). The Archaea sink is removed (non-archaeal-to-Archaea at most 2%);
#      dominant kingdom correct for 80% of isolates (Bacteria 92, Plantae 94, Fungi 89).
#   c  Per-kingdom self-recovery, raw versus corrected (diagonal of a vs b).
# Rebuilt on the Figure 5 v2 (ArchLips-extended) decomposition under the
# locked strict release ncbi-phylum-2026-08-04-v1: 16 analysis phyla,
# 736-feature substrate, 164 isolates, leave-one-out reference.
# Producer: paper2_repro/scripts/suppfig5_negative_control_strict16.py
# The submitted analysis-19/16_negative_control/build_figure.py is not in
# the repository; published panel values are NOT reproduced here.
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
data_dir <- file.path(SCRIPT_DIR, "data", "supp_negative_control")
out_dir  <- file.path(SCRIPT_DIR, "out")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Six kingdoms in the canonical style order (Bacteria, Archaea, Fungi,
# Plantae, Animalia, Protozoa). build_figure.py used a Blues-then-Archaea
# ordering for display; we use KINGDOM_ORDER as the R-port convention.
K6 <- KINGDOM_ORDER

read_conf <- function(fn) {
  m <- read.csv(file.path(data_dir, fn), check.names = FALSE)
  rownames(m) <- m$true_kingdom
  as.matrix(m[K6, K6])            # reorder rows + cols to KINGDOM_ORDER
}
raw <- read_conf("negative_control_kingdom_confusion.csv")
cor <- read_conf("negative_control_corrected_kingdom_confusion.csv")

# Heatmap colour ramp matching matplotlib "Blues" (vmin 0, vmax 100).
BLUES <- c("#F7FBFF", "#DEEBF7", "#C6DBEF", "#9ECAE1", "#6BAED6",
           "#4292C6", "#2171B5", "#08519C", "#08306B")
DIAG_BOX <- "#D55E00"            # Wong vermillion, matches build_figure.py ec
RAW_FILL <- "#BBBBBB"
COR_FILL <- "#1F4E79"
COLOR_LIGHT <- "#777777"

# ================================================================
#  Panels a / b: confusion-matrix tile heatmaps
# ================================================================
conf_long <- function(M) {
  df <- as.data.frame(as.table(M))
  names(df) <- c("true", "assigned", "pct")
  df$true     <- factor(df$true,     levels = rev(K6))   # top row = Bacteria
  df$assigned <- factor(df$assigned, levels = K6)
  df$is_diag  <- as.character(df$true) == as.character(df$assigned)
  df
}

make_heat <- function(M, title) {
  df <- conf_long(M)
  ggplot(df, aes(assigned, true)) +
    geom_tile(aes(fill = pct), colour = "white", linewidth = 0.3) +
    # orange box around the diagonal (correct-kingdom) cells
    geom_tile(data = df[df$is_diag, ], fill = NA, colour = DIAG_BOX, linewidth = 0.6) +
    geom_text(aes(label = sprintf("%.0f", pct),
                  colour = pct > 55, fontface = ifelse(is_diag, "bold", "plain")),
              size = pt2mm(6)) +
    scale_fill_gradientn(colours = BLUES, limits = c(0, 100),
                         name = "Mean %\nassigned",
                         breaks = c(0, 25, 50, 75, 100)) +
    scale_colour_manual(values = c(`FALSE` = "black", `TRUE` = "white"),
                        guide = "none") +
    scale_x_discrete(position = "bottom", expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    labs(x = "Assigned group", y = "True group (pure isolate)", title = title) +
    coord_equal(clip = "off") +
    theme_nature() +
    theme(
      axis.line        = element_blank(),
      axis.ticks       = element_blank(),
      axis.text.x      = element_text(size = 6, angle = 35, hjust = 1),
      axis.text.y      = element_text(size = 6),
      plot.title       = element_text(size = 7, face = "bold", hjust = 0),
      legend.key.width = unit(2.4, "mm"),
      legend.key.height = unit(3.2, "mm"),
      legend.title     = element_text(size = 5.5, face = "bold"),
      legend.text      = element_text(size = 5)
    )
}

pa <- make_heat(raw, "Uncorrected pipeline")
pb <- make_heat(cor, "Corrected v2 (IS + RIE floor +\nArchLips reference + rule C)")

# ================================================================
#  Panel c: per-kingdom self-recovery, raw vs corrected
#  (diagonal of each confusion matrix, exactly as build_figure.py)
# ================================================================
own <- data.frame(
  kingdom = factor(K6, levels = rev(K6)),     # Bacteria at top
  Raw       = diag(raw),
  Corrected = diag(cor)
) |>
  pivot_longer(c(Raw, Corrected), names_to = "method", values_to = "own_pct") |>
  mutate(method = factor(method, levels = c("Raw", "Corrected")))

CHANCE <- 100 / 6                              # 1/6 = 16.7%

pc <- ggplot(own, aes(own_pct, kingdom, fill = method)) +
  geom_col(position = position_dodge(width = 0.72), width = 0.66) +
  geom_vline(xintercept = CHANCE, linetype = "dotted", linewidth = 0.4,
             colour = COLOR_LIGHT) +
  annotate("text", x = CHANCE + 2, y = 0.62, label = "chance (1/6)",
           hjust = 0, size = pt2mm(5), colour = COLOR_LIGHT) +
  scale_fill_manual(values = c(Raw = RAW_FILL, Corrected = COR_FILL), name = NULL) +
  scale_x_continuous(limits = c(0, 100), breaks = c(0, 25, 50, 75, 100),
                     expand = expansion(mult = c(0, 0.02))) +
  labs(x = "% assigned to own group", y = NULL,
       title = "Self-recovery per group") +
  coord_cartesian(clip = "off") +
  theme_nature() +
  theme(
    axis.text.y      = element_text(size = 6),
    plot.title       = element_text(size = 7, face = "bold", hjust = 0),
    legend.position        = "inside",
    legend.position.inside = c(0.99, 0.02),
    legend.justification   = c(1, 0),
    legend.text      = element_text(size = 5.5),
    legend.key.size  = unit(2.6, "mm"),
    legend.background = element_blank()
  )

# ================================================================
#  Compose
# ================================================================
fig <- (pa | pb | pc) +
  plot_layout(widths = c(1, 1, 0.95)) +
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

save_all(fig, file.path(out_dir, "Supplementary_Fig5_negative_control"), NC_DOUBLE, 72)
cat("Supplementary Fig 5 ->", out_dir, "\n")
