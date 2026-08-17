#!/usr/bin/env Rscript
# Supplementary Fig. 1 - submitted 4-panel layout, strict-16 data, R render.
#   a  cross-batch PCoA (168 samples), organism-group colours, silhouette 0.090
#   b  same ordination coloured by acquisition batch, silhouette 0.122
#   c  within-batch PCoA, batch 02 only, silhouette 0.106
#   d  within-batch Mantel r per batch (positive mode), permutation significance
# Layout mirrors the figure embedded in SLA-Supplementary-nature-comms.docx;
# styling follows soilmass_style.R (Wong palette, theme_nature) per the
# coauthor-package convention. Panels a-c are the submitted ordinations with
# relabelled source data (visuals unchanged); panel d carries the strict-16
# Mantel results (ncbi-phylum-2026-08-04-v1).
# Nature Communications: 183 mm double column, text 5-7 pt.

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  library(scales)
})

.args <- commandArgs(trailingOnly = FALSE)
.f <- sub("--file=", "", .args[grep("--file=", .args)])
SCRIPT_DIR <- if (length(.f) > 0) dirname(normalizePath(.f)) else getwd()
source(file.path(SCRIPT_DIR, "soilmass_style.R"))

data_dir <- file.path(SCRIPT_DIR, "data")
out_dir  <- file.path(SCRIPT_DIR, "out")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

coords <- read.csv(file.path(data_dir, "pcoa_coords.csv"))
b02    <- read.csv(file.path(data_dir, "pcoa_batch_02.csv"))
mant   <- read.csv(file.path(data_dir, "panel_d_mantel.csv"))

# ---- shared PCoA panel builder (legend inside, silhouette bottom-right) ----
pcoa_panel <- function(df, colour_var, palette, sil_label, xlab, ylab,
                       legend_cols = 2) {
  df$grp <- df[[colour_var]]
  if (colour_var == "kingdom") {
    df$grp <- factor(df$grp, levels = KINGDOM_ORDER)
    counts <- table(df$grp)
    df$n <- as.integer(counts[as.character(df$grp)])
    df <- df[order(-df$n), ]                 # big groups behind
  }
  ggplot(df, aes(pcoa1, pcoa2, fill = grp)) +
    geom_point(shape = 21, colour = "white", size = 1.6, stroke = 0.25,
               alpha = 0.9) +
    scale_fill_manual(values = palette, name = NULL) +
    annotate("text", x = Inf, y = -Inf, label = sil_label,
             hjust = 1.05, vjust = -0.9, size = pt2mm(6), colour = "#666666") +
    labs(x = xlab, y = ylab) +
    coord_cartesian(clip = "off") +
    theme_nature() +
    theme(
      legend.position        = "inside",
      legend.position.inside = c(0.01, 0.99),
      legend.justification   = c(0, 1),
      legend.text            = element_text(size = 5.5),
      legend.key.size        = unit(2.4, "mm"),
      legend.background      = element_blank()
    ) +
    guides(fill = guide_legend(ncol = legend_cols,
                               override.aes = list(size = 1.8)))
}

# ---- a: cross-batch, organism groups --------------------------------------
p_a <- pcoa_panel(coords, "kingdom", KINGDOM_COLOURS,
                  "group silhouette = 0.090",
                  "PCoA 1 (12.3%)", "PCoA 2 (7.5%)")

# ---- b: cross-batch, batches ----------------------------------------------
coords$batch_short <- sub("^([0-9]+)_.*$", "batch \\1",
                          sub("^batch_", "", coords$batch))
batch_levels <- sort(unique(coords$batch_short))
batch_pal <- setNames(
  colorRampPalette(c("#2A2D5E", "#5C6068", "#9C9265", "#E3CD4B", "#FDE725"))(
    length(batch_levels)), batch_levels)
p_b <- pcoa_panel(transform(coords, batch = batch_short), "batch", batch_pal,
                  "batch silhouette = 0.122",
                  "PCoA 1 (12.3%)", "PCoA 2 (7.5%)", legend_cols = 1)

# ---- c: within batch 02, organism groups ----------------------------------
p_c <- pcoa_panel(b02, "kingdom", KINGDOM_COLOURS,
                  "group silhouette = 0.106",
                  "PCoA 1 (21.7%)", "PCoA 2 (12.7%)", legend_cols = 1)

# ---- d: within-batch Mantel dot plot (strict-16) --------------------------
mant$label <- factor(mant$label, levels = rev(mant$label))
mant$sig   <- mant$p_perm < 0.05
mant$col   <- ifelse(mant$sig, "#0072B2", "#B0B0B0")
mant$mark  <- ifelse(mant$sig, "*", "n.s.")
r_lo <- sprintf("%.2f", min(mant$r))
r_hi <- sprintf("%.2f", max(mant$r))

p_d <- ggplot(mant, aes(x = r, y = label)) +
  geom_point(size = 2.6, colour = mant$col) +
  geom_text(aes(label = mark), vjust = -1.1, size = pt2mm(6.5),
            colour = mant$col, fontface = "bold") +
  annotate("text", x = 1.04, y = 1.62, hjust = 1, vjust = 0, size = pt2mm(5.5),
           colour = "#999999", lineheight = 1.1,
           label = sprintf(
             "r = %s-%s (positive mode)\nsignificance limited by permutation\nresolution in small batches",
             r_lo, r_hi)) +
  scale_x_continuous(limits = c(0, 1.05), breaks = seq(0, 1, 0.2)) +
  labs(x = "Within-batch Mantel r (lipidome vs phylogeny)", y = NULL) +
  coord_cartesian(clip = "off") +
  theme_nature() +
  theme(axis.text.y = element_text(size = 6))

# ---- compose --------------------------------------------------------------
fig <- (p_a | p_b) / (p_c | p_d) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = NC_TAG_PT, face = "bold"))

save_all <- function(plot, stem, w_mm, h_mm) {
  w <- w_mm * MM_TO_IN; h <- h_mm * MM_TO_IN
  ggsave(paste0(stem, ".pdf"), plot, width = w, height = h, device = cairo_pdf, bg = "white")
  ggsave(paste0(stem, ".png"), plot, width = w, height = h, dpi = 600, bg = "white",
         device = ragg::agg_png)
  if (requireNamespace("svglite", quietly = TRUE)) {
    ggsave(paste0(stem, ".svg"), plot, width = w, height = h, device = svglite::svglite, bg = "white")
  }
}
save_all(fig, file.path(out_dir, "Supplementary_Fig1_submitted_layout_strict16"), NC_DOUBLE, 160)
cat("Supplementary Fig 1 (submitted layout, strict-16) ->", out_dir, "\n")
