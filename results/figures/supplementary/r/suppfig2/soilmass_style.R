library(ggplot2)
library(dplyr)
library(patchwork)
library(scales)

# ── Nature Communications figure spec ──────────────────────
# Text: Helvetica/Arial, 5–7 pt everywhere (axis, legend, annotation).
# Panel tags: 8 pt bold lowercase.  Lines: 0.25–1 pt.
# Widths: single 89 mm, 1.5-col 120 mm, double 183 mm.  Max height 247 mm.
# Legends: outside data area, compact titles (≤4 words or omit).
# ───────────────────────────────────────────────────────────

NC_SINGLE  <- 89
NC_DOUBLE  <- 183
NC_MAX_H   <- 247
MM_TO_IN   <- 1 / 25.4

# ggplot2 geom_text size is in mm; theme element_text size is in pt.
# 1 pt = 0.3528 mm.  Convenience converter:
pt2mm <- function(pt) pt * 0.3528

# Nature figure text range (pt)
NC_MIN_PT <- 5
NC_MAX_PT <- 7
NC_TAG_PT <- 8

# ── Wong colour-blind-safe palette ─────────────────────────
KINGDOM_COLOURS <- c(
  Bacteria = "#0072B2", Archaea  = "#E69F00", Fungi    = "#009E73",
  Plantae  = "#56B4E9", Animalia = "#D55E00", Protozoa = "#CC79A7"
)
KINGDOM_ORDER <- c("Bacteria", "Archaea", "Fungi", "Plantae", "Animalia", "Protozoa")

BATCH_COLOURS <- c(
  "Batch 01" = "#1b9e77", "Batch 02" = "#7570b3", "Batch 03" = "#666666",
  "Batch 04" = "#e7298a", "Batch 05" = "#e6ab02", "Batch 06" = "#a6761d"
)

MODE_COLOURS <- c(Positive = "#0072B2", Negative = "#E69F00")

# ── Nature-compliant ggplot2 theme ─────────────────────────
theme_nature <- function(base_size = 7) {
  theme_classic(base_size = base_size, base_family = "Helvetica") %+replace%
    theme(
      line               = element_line(linewidth = 0.3, colour = "black"),
      axis.line          = element_line(linewidth = 0.3),
      axis.ticks         = element_line(linewidth = 0.3),
      axis.ticks.length  = unit(1.2, "pt"),
      axis.title         = element_text(size = 7, colour = "black"),
      axis.text          = element_text(size = 6, colour = "black"),
      legend.title       = element_text(size = 6, face = "bold"),
      legend.text        = element_text(size = 5),
      legend.key.size    = unit(3, "mm"),
      legend.key.height  = unit(3, "mm"),
      legend.key.width   = unit(3, "mm"),
      legend.background  = element_blank(),
      legend.box.background = element_blank(),
      legend.margin      = margin(0, 0, 0, 0),
      legend.spacing     = unit(1, "mm"),
      plot.title         = element_text(size = 7, face = "bold", hjust = 0),
      plot.margin        = margin(2, 2, 2, 2, "mm"),
      strip.background   = element_blank(),
      strip.text         = element_text(size = 6),
      panel.grid         = element_blank()
    )
}

# ── Helpers ────────────────────────────────────────────────
# Wrap a phylum name in an HTML colour span (for ggtext axis labels)
colour_label <- function(phylum, pk_lookup) {
  kingdom <- pk_lookup[phylum]
  col <- KINGDOM_COLOURS[kingdom]
  paste0("<span style='color:", col, "'>", phylum, "</span>")
}
