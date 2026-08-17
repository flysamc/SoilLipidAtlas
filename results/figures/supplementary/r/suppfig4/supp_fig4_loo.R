#!/usr/bin/env Rscript
# Supplementary Fig. 4  -  NNLS leave-one-out classification confusion matrix.
# Port of analysis-17/positive/figures/fig5a_nnls_confusion/fig5a_nnls_confusion.py.
# Single-panel phylum x phylum heatmap: rows = true phylum, columns = predicted
# phylum, cells shaded by row-normalised fraction and labelled with sample counts.
# Nature Communications: single column (89 mm), all text 5-7 pt, Helvetica.

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(ggtext)
  library(scales)
})

.args <- commandArgs(trailingOnly = FALSE)
.f <- sub("--file=", "", .args[grep("--file=", .args)])
SCRIPT_DIR <- if (length(.f) > 0) dirname(normalizePath(.f)) else getwd()
source(file.path(SCRIPT_DIR, "soilmass_style.R"))

# -- Data ---------------------------------------------------------
data_dir <- file.path(SCRIPT_DIR, "data", "supp_loo")
out_dir  <- file.path(SCRIPT_DIR, "out")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# phylum -> kingdom (soilmass_style.R has no such map; define the 17 present
# plus the rest of the atlas for completeness).
PHYLUM_TO_KINGDOM <- c(
  Bacillota = "Bacteria", Actinomycetota = "Bacteria", Pseudomonadota = "Bacteria", Cyanobacteriota = "Bacteria",
  Euryarchaeota = "Archaea", Thermoproteota = "Archaea", Methanobacteriota = "Archaea", Crenarchaeota = "Archaea",
  Ascomycota = "Fungi", Basidiomycota = "Fungi", Mucoromycota = "Fungi",
  Trachaeophyta = "Plantae", Chlorophyta = "Plantae", Bryophyta = "Plantae",
  Marchantiophyta = "Plantae", Magnoliophyta = "Plantae", Charophyta = "Plantae",
  Arthropoda = "Animalia", Nematoda = "Animalia", Mollusca = "Animalia",
  Amoebozoa = "Protozoa", Bicosoecida = "Protozoa", Cercozoa = "Protozoa", Heterolobosea = "Protozoa",
  # units created by the 2026-08 taxonomy correction (ncbi-phylum-2026-08-04-v1)
  Discosea = "Protozoa", Evosea = "Protozoa", Halobacteriota = "Archaea",
  Tracheophyta = "Plantae", Mortierellomycota = "Fungi", Streptophyta = "Plantae")

cm <- read.csv(file.path(data_dir, "loo_confusion_matrix.csv"),
               row.names = 1, check.names = FALSE)
phyla <- rownames(cm)
raw   <- as.matrix(cm)
storage.mode(raw) <- "double"

# Row-normalise for shading (each row sums to 1); guard zero rows.
row_sums <- rowSums(raw)
row_sums[row_sums == 0] <- 1
norm <- raw / row_sums

# Per-phylum accuracy (diagonal / row sum) and overall.
diag_v        <- diag(raw)
row_sums_real <- rowSums(raw)
accuracy      <- diag_v / row_sums_real * 100
total_correct <- sum(diag_v)
total_samples <- sum(raw)
overall_acc   <- total_correct / total_samples * 100
n             <- length(phyla)

# -- Order phyla by kingdom, then alphabetically (matches matplotlib sort_key) --
KINGDOM_IDX <- setNames(seq_along(KINGDOM_ORDER), KINGDOM_ORDER)
sort_key <- function(p) {
  k <- PHYLUM_TO_KINGDOM[[p]]
  ki <- if (!is.null(k) && k %in% names(KINGDOM_IDX)) KINGDOM_IDX[[k]] else 99
  sprintf("%02d_%s", ki, p)
}
ord <- order(vapply(phyla, sort_key, character(1)))
sorted_phyla <- phyla[ord]

# Kingdom-coloured axis labels via ggtext markdown.
phylum_kingdom <- PHYLUM_TO_KINGDOM[sorted_phyla]
phylum_colour  <- KINGDOM_COLOURS[phylum_kingdom]
html_lab <- setNames(
  sprintf("<span style='color:%s'>%s</span>", phylum_colour, sorted_phyla),
  sorted_phyla)

# -- Long-format tile data ---------------------------------------
# x = predicted (columns), y = true (rows). y reversed so the first true
# phylum sits at the TOP, matching imshow's default origin.
df <- expand.grid(true = sorted_phyla, pred = sorted_phyla,
                  KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
df$count <- mapply(function(tr, pr) raw[tr, pr], df$true, df$pred)
df$frac  <- mapply(function(tr, pr) norm[tr, pr], df$true, df$pred)
df$true  <- factor(df$true, levels = rev(sorted_phyla))   # top-to-bottom
df$pred  <- factor(df$pred, levels = sorted_phyla)

# Cell-label data: only non-zero counts; white text on dark cells (frac > 0.55),
# dark text otherwise; bold on the diagonal (true == pred).
lab <- df[df$count > 0, ]
lab$txt_col <- ifelse(lab$frac > 0.55, "white", "#222222")   # COLOR_TEXT
lab$face    <- ifelse(as.character(lab$true) == as.character(lab$pred),
                      "bold", "plain")

# -- Per-phylum accuracy annotation (right of the matrix) --------
acc_df <- data.frame(
  phylum = sorted_phyla,
  acc    = accuracy[sorted_phyla],
  stringsAsFactors = FALSE)
acc_df$true  <- factor(acc_df$phylum, levels = rev(sorted_phyla))
acc_df$label <- ifelse(acc_df$acc > 0, sprintf("%.0f%%", acc_df$acc), "0%")
acc_df$col   <- ifelse(acc_df$acc >= 75, "#2E7D32",
                ifelse(acc_df$acc >= 40, "#F57F17", "#C62828"))
acc_x <- n + 1.2   # column position just right of the last predicted column

# -- Kingdom group separator lines -------------------------------
# Boundaries between kingdom blocks in the sorted ordering.
king_seq  <- PHYLUM_TO_KINGDOM[sorted_phyla]
breaks_at <- which(king_seq[-1] != king_seq[-n])  # index i means break after row i
# In x (predicted): vertical line between column i and i+1 -> at i + 0.5.
# In y (true, reversed): horizontal line between true-rank i and i+1.
vline_x <- breaks_at + 0.5
hline_y <- (n - breaks_at) + 0.5   # reversed-axis position

# -- Build the heatmap -------------------------------------------
p <- ggplot(df, aes(x = pred, y = true)) +
  geom_tile(aes(fill = frac), colour = "grey90", linewidth = 0.15) +
  scale_fill_gradient(low = "#FFFFFF", high = "#08306B",   # matplotlib Blues
                      limits = c(0, 1), name = "Fraction of\ntrue-class samples",
                      breaks = c(0, 0.25, 0.5, 0.75, 1.0)) +
  # kingdom group separators
  geom_vline(xintercept = vline_x, colour = "#999999", linewidth = 0.3) +
  geom_hline(yintercept = hline_y, colour = "#999999", linewidth = 0.3) +
  # cell counts
  geom_text(data = lab, aes(label = count, fontface = face),
            colour = lab$txt_col, size = pt2mm(5)) +
  # per-phylum accuracy column on the right
  geom_text(data = acc_df, aes(x = acc_x, y = true, label = label),
            colour = acc_df$col, hjust = 0, size = pt2mm(5), fontface = "bold",
            inherit.aes = FALSE) +
  annotate("text", x = acc_x, y = n + 0.85, label = "Acc",
           hjust = 0, vjust = 0, size = pt2mm(5), fontface = "bold",
           colour = "#222222") +
  scale_x_discrete(labels = html_lab, expand = expansion(add = c(0.5, 0.5))) +
  scale_y_discrete(labels = html_lab, expand = expansion(add = c(0.5, 1.0))) +
  coord_fixed(clip = "off") +
  labs(x = "Predicted phylum", y = "True phylum",
       title = sprintf("Leave-one-out classification: %d/%d (%.1f%%) across %d phyla",
                       total_correct, total_samples, overall_acc, n)) +
  theme_nature() +
  theme(
    axis.text.x       = element_markdown(size = 5, angle = 45, hjust = 1, vjust = 1),
    axis.text.y       = element_markdown(size = 5),
    axis.ticks        = element_blank(),
    axis.line         = element_blank(),
    plot.title        = element_text(size = 6, face = "plain", colour = "#666666",
                                     hjust = 0, margin = margin(b = 3)),
    legend.position   = "right",
    legend.title      = element_text(size = 5, face = "plain"),
    legend.text       = element_text(size = 5),
    legend.key.width  = unit(2.5, "mm"),
    legend.key.height = unit(4, "mm"),
    legend.box.spacing = unit(6, "mm"),
    plot.margin       = margin(2, 4, 4, 2, "mm")
  )

# -- Kingdom legend (compact, bottom) -----------------------------
# Build as a manual fill scale on a dummy layer so the kingdom swatches
# appear alongside the continuous fraction legend.
present_kingdoms <- KINGDOM_ORDER[KINGDOM_ORDER %in% unique(king_seq)]
king_leg <- data.frame(
  pred = factor(sorted_phyla[1], levels = sorted_phyla),
  true = factor(sorted_phyla[1], levels = rev(sorted_phyla)),
  kingdom = factor(present_kingdoms, levels = present_kingdoms))

p <- p +
  geom_point(data = king_leg,
             aes(x = pred, y = true, colour = kingdom),
             alpha = 0, size = 0, inherit.aes = FALSE) +
  scale_colour_manual(values = KINGDOM_COLOURS, breaks = present_kingdoms,
                      name = "Organism group") +
  guides(colour = guide_legend(
    order = 2, override.aes = list(alpha = 1, size = 1.8, shape = 15)),
    fill = guide_colourbar(order = 1))

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

# 1.5-column width: a square 17x17 heatmap plus an accuracy column and two
# legends need more room than a single 89 mm column allows.
NC_1P5 <- 120
save_all(p, file.path(out_dir, "Supplementary_Fig4_nnls_loo"), NC_1P5, 110)
cat(sprintf("Supplementary Fig 4 -> %s  (overall %d/%d = %.1f%%, %d phyla)\n",
            out_dir, total_correct, total_samples, overall_acc, n))
