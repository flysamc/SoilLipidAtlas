#!/usr/bin/env Rscript
# Supplementary Fig. 3  -  MS2LDA motif enrichment heatmap (top motifs x phyla).
# Port of figures_v5 panels/ms2lda_heatmap.py + compose/ed_ms2lda.py.
# Nature Communications: 183 mm wide, all text 5-7 pt, Helvetica.

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  library(ggtext)
  library(reshape2)
  library(scales)
})

.args <- commandArgs(trailingOnly = FALSE)
.f <- sub("--file=", "", .args[grep("--file=", .args)])
SCRIPT_DIR <- if (length(.f) > 0) dirname(normalizePath(.f)) else getwd()
source(file.path(SCRIPT_DIR, "soilmass_style.R"))

# -- Data ---------------------------------------------------------
data_dir <- file.path(SCRIPT_DIR, "data", "supp_ms2lda")
out_dir  <- file.path(SCRIPT_DIR, "out")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

TOPN <- 15

# Phylum order + kingdom mapping (ported from figures_v5 soilmass_style.py; the
# shared R style file carries kingdoms only, so the phylum-level keys live here).
PHYLUM_ORDER <- c(
  "Bacillota", "Actinomycetota", "Pseudomonadota", "Cyanobacteriota",
  "Euryarchaeota", "Thermoproteota", "Methanobacteriota", "Crenarchaeota",
  "Ascomycota", "Basidiomycota", "Mucoromycota",
  "Trachaeophyta", "Streptophyta", "Chlorophyta", "Bryophyta", "Marchantiophyta", "Magnoliophyta", "Charophyta",
  "Arthropoda", "Nematoda", "Mollusca",
  "Amoebozoa", "Discosea", "Evosea", "Bicosoecida", "Cercozoa", "Heterolobosea",
  "Virus")
PHYLUM_TO_KINGDOM <- c(
  Bacillota = "Bacteria", Actinomycetota = "Bacteria", Pseudomonadota = "Bacteria", Cyanobacteriota = "Bacteria",
  Euryarchaeota = "Archaea", Thermoproteota = "Archaea", Methanobacteriota = "Archaea", Crenarchaeota = "Archaea",
  Ascomycota = "Fungi", Basidiomycota = "Fungi", Mucoromycota = "Fungi",
  Trachaeophyta = "Plantae", Chlorophyta = "Plantae", Bryophyta = "Plantae",
  Marchantiophyta = "Plantae", Magnoliophyta = "Plantae", Charophyta = "Plantae",
  Arthropoda = "Animalia", Nematoda = "Animalia", Mollusca = "Animalia",
  Amoebozoa = "Protozoa", Bicosoecida = "Protozoa", Cercozoa = "Protozoa", Heterolobosea = "Protozoa",
  # units of the 2026-08 strict release (ncbi-phylum-2026-08-04-v1)
  Streptophyta = "Plantae", Discosea = "Protozoa", Evosea = "Protozoa",
  Virus = "Viruses")

# Wong palette + grey for Viruses (R style file omits Viruses)
KCOL <- c(KINGDOM_COLOURS, Viruses = "#999999")

# Diverging colour map (ported verbatim from soilmass_style.py DIVERGING_COLORS)
DIVERGING_COLORS <- c("#08306B", "#2171B5", "#6BAED6", "#C6DBEF", "#FFFFFF",
                      "#FDD49E", "#FC8D59", "#E34A33", "#B30000")

# -- Enrichment matrix --------------------------------------------
# log2 enrichment per motif x phylum; pivot (missing combos -> 0); restrict to
# phyla present in PHYLUM_ORDER order; keep the TOPN motifs with the highest
# single-phylum enrichment, ordered by that maximum (descending).
enrich <- function(fname) {
  d <- read.csv(file.path(data_dir, fname))
  d$log2e <- log2((d$observed_frac + 1e-3) / (d$expected_frac + 1e-3))
  mat <- acast(d, motif ~ phylum, value.var = "log2e", fill = 0, fun.aggregate = mean)
  phyla <- PHYLUM_ORDER[PHYLUM_ORDER %in% colnames(mat)]
  mat <- mat[, phyla, drop = FALSE]
  rowmax <- apply(mat, 1, max)
  top <- names(sort(rowmax, decreasing = TRUE))[seq_len(TOPN)]
  list(mat = mat[top, , drop = FALSE], phyla = phyla, top = top)
}

# -- One heatmap --------------------------------------------------
heatmap_panel <- function(fname, title) {
  e <- enrich(fname)
  long <- melt(e$mat, varnames = c("motif", "phylum"), value.name = "log2e")
  long$motif  <- factor(long$motif,  levels = rev(e$top))     # top motif at top
  long$phylum <- factor(long$phylum, levels = e$phyla)
  vmax <- max(abs(e$mat))

  kof <- PHYLUM_TO_KINGDOM[e$phyla]
  pcol <- ifelse(is.na(KCOL[kof]), "#333333", KCOL[kof])
  x_html <- setNames(sprintf("<span style='color:%s'>%s</span>", pcol, e$phyla), e$phyla)

  ggplot(long, aes(phylum, motif, fill = log2e)) +
    geom_tile(colour = NA) +
    scale_fill_gradientn(colours = DIVERGING_COLORS, limits = c(-vmax, vmax),
                         name = expression(log[2] ~ "enrichment")) +
    scale_x_discrete(labels = x_html, expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    labs(x = NULL, y = NULL, title = sprintf("%s  (top %d motifs)", title, TOPN)) +
    theme_nature() +
    theme(
      axis.text.x  = element_markdown(angle = 90, hjust = 1, vjust = 0.5, size = 5),
      axis.text.y  = element_text(size = 5),
      axis.ticks   = element_blank(),
      axis.line    = element_blank(),
      plot.title   = element_text(size = 7, face = "bold", hjust = 0),
      legend.title = element_text(size = 5.5, face = "bold"),
      legend.text  = element_text(size = 5),
      legend.key.width  = unit(2, "mm"),
      legend.key.height = unit(4, "mm")
    )
}

p_pos <- heatmap_panel("motif_phylum_enrichment_pos.csv", "Positive mode")
p_neg <- heatmap_panel("motif_phylum_enrichment_neg.csv", "Negative mode")

# -- Compose (separate colour scales, so guides are NOT collected) ------
ed4 <- (p_pos | p_neg) +
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

save_all(ed4, file.path(out_dir, "Supplementary_Fig3_ms2lda_motifs"), NC_DOUBLE, 80)
cat("Supplementary Fig 3 ->", out_dir, "\n")
