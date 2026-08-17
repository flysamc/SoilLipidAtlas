#!/usr/bin/env Rscript
# Strict Figure 3c review-only producer.
#
# Target under the coauthor-organized main-figure layout:
#   Figure 3c = pairwise POS/NEG lipidome-similarity heatmaps.
#
# This producer intentionally does not use the legacy 17-phylum heatmap inputs
# or the earlier Figure 3c forest-plot producer. It consumes the strict-v1
# Figure 3a Bray-Curtis matrices and the strict Figure 3a display order.

suppressPackageStartupMessages({
  library(ggplot2)
  library(reshape2)
  library(scales)
  library(jsonlite)
  library(digest)
})

full_command <- commandArgs(trailingOnly = FALSE)
file_argument <- grep("^--file=", full_command, value = TRUE)
if (length(file_argument) != 1) stop("Cannot resolve the Figure 3c R producer path")
script_path <- normalizePath(sub("^--file=", "", file_argument), mustWork = TRUE)

parse_args <- function(values) {
  parsed <- list()
  for (value in values) {
    if (!startsWith(value, "--") || !grepl("=", value, fixed = TRUE)) {
      stop("Arguments must use --name=value syntax: ", value)
    }
    pieces <- strsplit(sub("^--", "", value), "=", fixed = TRUE)[[1]]
    parsed[[pieces[[1]]]] <- paste(pieces[-1], collapse = "=")
  }
  parsed
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
required_args <- c("source-panel-dir", "ordered-view-dir", "substrate-dir", "output-dir", "style")
missing_args <- setdiff(required_args, names(args))
if (length(missing_args) > 0) stop("Missing arguments: ", paste(missing_args, collapse = ", "))

source_panel_dir <- normalizePath(args[["source-panel-dir"]], mustWork = TRUE)
ordered_view_dir <- normalizePath(args[["ordered-view-dir"]], mustWork = TRUE)
substrate_dir <- normalizePath(args[["substrate-dir"]], mustWork = TRUE)
style_path <- normalizePath(args[["style"]], mustWork = TRUE)
output_dir <- args[["output-dir"]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
output_dir <- normalizePath(output_dir, mustWork = TRUE)
started <- Sys.time()

# The repository style file is recorded below as a provenance input, but it
# imports patchwork, which is not present in the pinned local R library. Keep
# the same numeric tokens, Helvetica base family, and Wong palette locally so
# this standalone panel remains runnable and style-compatible.
NC_DOUBLE <- 183
MM_TO_IN <- 1 / 25.4
NC_TAG_PT <- 8
KINGDOM_COLOURS <- c(
  Bacteria = "#0072B2", Archaea = "#E69F00", Fungi = "#009E73",
  Plantae = "#56B4E9", Animalia = "#D55E00", Protozoa = "#CC79A7"
)
theme_nature <- function(base_size = 7) {
  theme_classic(base_size = base_size, base_family = "Helvetica") %+replace%
    theme(
      line = element_line(linewidth = 0.3, colour = "black"),
      axis.line = element_line(linewidth = 0.3),
      axis.ticks = element_line(linewidth = 0.3),
      axis.ticks.length = grid::unit(1.2, "pt"),
      axis.title = element_text(size = 7, colour = "black"),
      axis.text = element_text(size = 6, colour = "black"),
      legend.title = element_text(size = 6, face = "bold"),
      legend.text = element_text(size = 5),
      legend.key.size = grid::unit(3, "mm"),
      legend.key.height = grid::unit(3, "mm"),
      legend.key.width = grid::unit(3, "mm"),
      legend.background = element_blank(),
      legend.box.background = element_blank(),
      legend.margin = margin(0, 0, 0, 0),
      legend.spacing = grid::unit(1, "mm"),
      plot.title = element_text(size = 7, face = "bold", hjust = 0),
      plot.margin = margin(2, 2, 2, 2, "mm"),
      strip.background = element_blank(),
      strip.text = element_text(size = 6),
      panel.grid = element_blank()
    )
}

TAXONOMY_RELEASE <- "ncbi-phylum-2026-08-04-v1"
STRICT_PHYLUMS <- c(
  "Actinomycetota", "Bacillota", "Pseudomonadota",
  "Methanobacteriota", "Thermoproteota",
  "Discosea", "Evosea", "Heterolobosea",
  "Mucoromycota", "Ascomycota", "Basidiomycota",
  "Chlorophyta", "Streptophyta",
  "Arthropoda", "Mollusca", "Nematoda"
)
DISPLAY_GROUP_ORDER <- c(
  "Bacteria", "Archaea", "Protists", "Fungi", "Viridiplantae", "Animalia"
)
DISPLAY_COLOURS <- c(
  Bacteria = unname(KINGDOM_COLOURS[["Bacteria"]]),
  Archaea = unname(KINGDOM_COLOURS[["Archaea"]]),
  Protists = unname(KINGDOM_COLOURS[["Protozoa"]]),
  Fungi = unname(KINGDOM_COLOURS[["Fungi"]]),
  Viridiplantae = unname(KINGDOM_COLOURS[["Plantae"]]),
  Animalia = unname(KINGDOM_COLOURS[["Animalia"]])
)

sha256 <- function(path) digest(file = path, algo = "sha256", serialize = FALSE)

file_record <- function(path) {
  absolute <- normalizePath(path, mustWork = TRUE)
  list(
    path = absolute,
    bytes = unname(file.info(absolute)$size),
    sha256 = sha256(absolute)
  )
}

read_matrix <- function(path) {
  tab <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  if (ncol(tab) != nrow(tab) + 1L) stop("Unexpected matrix shape: ", path)
  row_names <- tab[[1]]
  mat <- as.matrix(tab[-1])
  storage.mode(mat) <- "double"
  rownames(mat) <- row_names
  colnames(mat) <- colnames(tab)[-1]
  mat
}

validate_bray <- function(mat, mode) {
  if (!identical(dim(mat), c(16L, 16L))) stop(mode, ": expected a 16 x 16 matrix")
  if (!setequal(rownames(mat), STRICT_PHYLUMS) || !setequal(colnames(mat), STRICT_PHYLUMS)) {
    stop(mode, ": matrix phylum names differ from strict v1")
  }
  mat <- mat[STRICT_PHYLUMS, STRICT_PHYLUMS, drop = FALSE]
  if (max(abs(mat - t(mat))) > 1e-14) stop(mode, ": Bray-Curtis matrix is asymmetric")
  if (max(abs(diag(mat))) > 1e-14) stop(mode, ": Bray-Curtis diagonal is nonzero")
  if (min(mat) < -1e-14 || max(mat) > 1 + 1e-14) stop(mode, ": distance outside [0, 1]")
  mat
}

manifest_path <- file.path(source_panel_dir, "stage_manifest.json")
source_summary_path <- file.path(source_panel_dir, "figure3a_summary.json")
group_map_path <- file.path(source_panel_dir, "phylum_display_groups_and_counts.csv")
pos_bray_path <- file.path(source_panel_dir, "pos_braycurtis.csv")
neg_bray_path <- file.path(source_panel_dir, "neg_braycurtis.csv")
order_path <- file.path(ordered_view_dir, "biological_display_order.csv")
substrate_summary_path <- file.path(substrate_dir, "substrate_summary.json")
substrate_manifest_path <- file.path(substrate_dir, "stage_manifest.json")

manifest <- fromJSON(manifest_path, simplifyDataFrame = TRUE)
source_summary <- fromJSON(source_summary_path, simplifyVector = TRUE)
substrate_summary <- fromJSON(substrate_summary_path, simplifyVector = TRUE)
substrate_manifest <- fromJSON(substrate_manifest_path, simplifyDataFrame = TRUE)

if (!identical(manifest$taxonomy_release, TAXONOMY_RELEASE)) stop("Source panel A release mismatch")
if (!identical(source_summary$taxonomy_release, TAXONOMY_RELEASE)) stop("Source panel A summary release mismatch")
if (!identical(substrate_summary$taxonomy_release, TAXONOMY_RELEASE)) stop("Substrate release mismatch")
if (!identical(substrate_summary$status, "pass")) stop("Strict substrate is not passing")
if (!identical(substrate_manifest$status, "pass")) stop("Strict substrate manifest is not passing")
if (!identical(as.integer(substrate_summary$positive$n_phyla), 16L) ||
    !identical(as.integer(substrate_summary$negative$n_phyla), 16L)) {
  stop("Strict substrate does not have 16 phyla in both modes")
}
if (!identical(as.integer(substrate_summary$positive$n_quality_features), 45525L) ||
    !identical(as.integer(substrate_summary$negative$n_quality_features), 14896L)) {
  stop("Strict substrate quality-feature counts changed")
}
if (!identical(as.integer(substrate_summary$positive$n_mapped_samples), 164L) ||
    !identical(as.integer(substrate_summary$negative$n_mapped_samples), 192L)) {
  stop("Strict substrate sample denominators changed")
}

order_table <- read.csv(order_path, check.names = FALSE, stringsAsFactors = FALSE)
if (!identical(order_table$phylum, STRICT_PHYLUMS)) {
  stop("Figure 3a ordered-view order is not the declared strict display order")
}
if (!identical(as.integer(order_table$display_rank), seq_along(STRICT_PHYLUMS))) {
  stop("Figure 3a ordered-view ranks are not consecutive")
}
if (any(!order_table$display_group %in% DISPLAY_GROUP_ORDER)) stop("Unknown display group")
if (any(table(order_table$display_group) < 1L)) stop("Empty display group")

group_map <- read.csv(group_map_path, check.names = FALSE, stringsAsFactors = FALSE)
if (!setequal(group_map$phylum, STRICT_PHYLUMS)) stop("Source group map is not strict v1")
group_lookup <- setNames(group_map$ecological_group, group_map$phylum)
if (any(unname(group_lookup[order_table$phylum]) != order_table$display_group)) {
  stop("Ordered-view display groups differ from source Figure 3a groups")
}

pos_bray <- validate_bray(read_matrix(pos_bray_path), "POS")
neg_bray <- validate_bray(read_matrix(neg_bray_path), "NEG")
display_order <- order_table$phylum

make_heatmap_data <- function(bray, mode_label) {
  similarity <- 1 - bray[display_order, display_order, drop = FALSE]
  diag(similarity) <- NA_real_
  long <- reshape2::melt(similarity, varnames = c("row_phylum", "column_phylum"), value.name = "similarity")
  long <- long[!is.na(long$similarity), , drop = FALSE]
  long$mode <- mode_label
  long$x <- match(as.character(long$column_phylum), display_order)
  long$y <- length(display_order) - match(as.character(long$row_phylum), display_order) + 1
  long$row_group <- unname(group_lookup[as.character(long$row_phylum)])
  long$column_group <- unname(group_lookup[as.character(long$column_phylum)])
  long
}

build_panel <- function(bray, mode_label, show_legend) {
  panel_data <- make_heatmap_data(bray, mode_label)
  x_labels <- data.frame(
    x = seq_along(display_order),
    y = -0.55,
    phylum = display_order,
    colour = unname(DISPLAY_COLOURS[group_lookup[display_order]])
  )
  y_labels <- data.frame(
    x = 0.45,
    y = length(display_order) - seq_along(display_order) + 1,
    phylum = display_order,
    colour = unname(DISPLAY_COLOURS[group_lookup[display_order]])
  )

  ggplot(panel_data, aes(x = x, y = y, fill = similarity)) +
    geom_tile(colour = "white", linewidth = 0.12) +
    geom_text(
      data = x_labels,
      aes(x = x, y = y, label = phylum, colour = colour),
      angle = 45, hjust = 1, vjust = 1, size = 1.65, inherit.aes = FALSE
    ) +
    geom_text(
      data = y_labels,
      aes(x = x, y = y, label = phylum, colour = colour),
      hjust = 1, size = 1.65, inherit.aes = FALSE
    ) +
    scale_fill_gradientn(
      colours = c("#fff5f0", "#fb6a4a", "#99000d"),
      limits = c(0, 0.55),
      oob = scales::squish,
      name = "Lipidome similarity\n(1 - Bray-Curtis)",
      na.value = "grey92"
    ) +
    scale_colour_identity() +
    scale_x_continuous(
      breaks = seq_along(display_order),
      labels = NULL,
      limits = c(0.4, length(display_order) + 0.5),
      expand = c(0, 0)
    ) +
    scale_y_continuous(
      breaks = seq_along(display_order),
      labels = NULL,
      limits = c(-1.25, length(display_order) + 0.5),
      expand = c(0, 0)
    ) +
    coord_fixed(clip = "off") +
    labs(title = mode_label, x = NULL, y = NULL) +
    theme_nature() +
    theme(
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.line = element_blank(),
      legend.position = if (show_legend) "right" else "none",
      legend.key.height = grid::unit(4, "mm"),
      legend.key.width = grid::unit(2, "mm"),
      legend.title = element_text(size = 5),
      legend.text = element_text(size = 5),
      plot.title = element_text(size = 6, hjust = 0.5),
      plot.margin = margin(3, 3, 12, 14, "mm")
    )
}

p_pos <- build_panel(pos_bray, "Positive mode", show_legend = FALSE)
p_neg <- build_panel(neg_bray, "Negative mode", show_legend = FALSE)
legend_grob <- local({
  g <- ggplot2::ggplotGrob(build_panel(neg_bray, "Negative mode", show_legend = TRUE))
  idx <- which(grepl("^guide-box", g$layout$name))
  if (!length(idx)) stop("Could not extract the Figure 3c colour legend")
  g$grobs[[idx[[1]]]]
})

draw_figure <- function(path, kind) {
  width_in <- NC_DOUBLE * MM_TO_IN
  height_in <- 95 * MM_TO_IN
  if (kind == "png") {
    grDevices::png(path, width = width_in, height = height_in, units = "in", res = 600, bg = "white")
  } else if (kind == "pdf") {
    grDevices::cairo_pdf(path, width = width_in, height = height_in, bg = "white")
  } else if (kind == "svg") {
    grDevices::svg(path, width = width_in, height = height_in, bg = "white")
  } else {
    stop("Unknown output device: ", kind)
  }
  grid::grid.newpage()
  # Equal heatmap viewports; legend sits in its own column so coord_fixed
  # does not shrink the negative-mode tiles to make room for the colour bar.
  grid::pushViewport(grid::viewport(
    layout = grid::grid.layout(1, 3, widths = grid::unit(c(1, 1, 0.20), "null"))
  ))
  print(p_pos, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
  print(p_neg, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 2))
  grid::pushViewport(grid::viewport(layout.pos.row = 1, layout.pos.col = 3))
  grid::grid.draw(legend_grob)
  grid::popViewport()
  grid::popViewport()
  grid::grid.text("c", x = grid::unit(2, "mm"), y = grid::unit(height_in * 25.4 - 2, "mm"),
                 just = c("left", "top"), gp = grid::gpar(fontsize = NC_TAG_PT, fontface = "bold"))
  grDevices::dev.off()
}

figure_base <- file.path(output_dir, "Figure_3c_strict_heatmap_review_only")
figure_paths <- paste0(figure_base, c(".png", ".pdf", ".svg"))
draw_figure(figure_paths[[1]], "png")
draw_figure(figure_paths[[2]], "pdf")
draw_figure(figure_paths[[3]], "svg")

pair_data <- function(bray, mode) {
  similarity <- 1 - bray[display_order, display_order, drop = FALSE]
  diag(similarity) <- NA_real_
  long <- reshape2::melt(similarity, varnames = c("row_phylum", "column_phylum"), value.name = "similarity")
  long <- long[!is.na(long$similarity), , drop = FALSE]
  long$mode <- mode
  long$row_group <- unname(group_lookup[as.character(long$row_phylum)])
  long$column_group <- unname(group_lookup[as.character(long$column_phylum)])
  long$pair_class <- ifelse(long$row_group == long$column_group, "Within group", "Cross group")
  long$clipped_at_0_55 <- long$similarity > 0.55
  long[, c("mode", "row_phylum", "column_phylum", "row_group", "column_group", "pair_class", "similarity", "clipped_at_0_55")]
}

pair_table <- rbind(pair_data(pos_bray, "POS"), pair_data(neg_bray, "NEG"))
pair_path <- file.path(output_dir, "strict_figure3c_heatmap_pair_data_v1_review_only.csv")
write.csv(pair_table, pair_path, row.names = FALSE, quote = FALSE)

matrix_summary <- function(bray, mode) {
  sim <- 1 - bray
  diag(sim) <- NA_real_
  within <- pair_table$mode == mode & pair_table$pair_class == "Within group"
  cross <- pair_table$mode == mode & pair_table$pair_class == "Cross group"
  list(
    mode = mode,
    n_unique_phyla = length(display_order),
    n_off_diagonal_cells = sum(!is.na(sim)),
    n_unique_pairs = sum(!is.na(sim)) / 2,
    min_similarity = min(sim, na.rm = TRUE),
    max_similarity = max(sim, na.rm = TRUE),
    mean_within_group_similarity = mean(pair_table$similarity[within]),
    mean_cross_group_similarity = mean(pair_table$similarity[cross]),
    n_values_clipped_above_0_55 = sum(pair_table$mode == mode & pair_table$clipped_at_0_55)
  )
}

authority_files <- c(
  file.path("submission_source", "SLA-Main-article-nature-comms.docx"),
  file.path("submission_source", "SLA-Supplementary-nature-comms.docx"),
  file.path("submission_source", "SLA-figures-and-legends-nature-comms.docx"),
  file.path("COAUTHOR_PACKAGE_2026-08-03", "01_reviewer_comments", "reviewer_comments_tracker.md"),
  file.path("paper2_repro", "portable.py")
)
authority_files <- authority_files[file.exists(authority_files)]

summary <- list(
  schema_version = 1,
  stage_id = "figure3c_strict_heatmap_review_only",
  status = "review_only_complete_pending_submitted_S1_reconciliation",
  taxonomy_release = TAXONOMY_RELEASE,
  scope = "Strict Figure 3c POS/NEG pairwise lipidome-similarity heatmaps; disconnected from manuscript, tables and response.",
  figure_identity = "Coauthor-organized Main Figure 3c; not internal v4.9 Figure 4c and not the earlier Figure 3c forest plot.",
  method = list(
    similarity = "1 - Bray-Curtis",
    matrix_source = "strict Figure 3a phylum-centroid Bray-Curtis matrices",
    phylum_order_source = "strict Figure 3a ordered-view biological_display_order.csv",
    diagonal = "masked",
    shared_colour_scale = c(0, 0.55),
    out_of_range = "squished to the upper colour limit",
    markers = "none in heatmap; cross/within markers belong to Figure 3b only"
  ),
  population = list(
    strict_phyla = display_order,
    POS = list(features = 45525L, samples = 164L),
    NEG = list(features = 14896L, samples = 192L)
  ),
  display_groups = order_table,
  mode_summary = list(POS = matrix_summary(pos_bray, "POS"), NEG = matrix_summary(neg_bray, "NEG")),
  gates = list(
    taxonomy_release_locked = TRUE,
    strict_16_phyla = TRUE,
    strict_substrate_pass = TRUE,
    current_coauthor_figure_mapping_recorded = TRUE,
    legacy_17_phylum_inputs_used = FALSE,
    v4_9_figure_organization_used = FALSE,
    ssu_tree_used_for_heatmap = FALSE,
    figure3_connected = FALSE,
    submitted_S1_cell_reconciliation = FALSE
  ),
  inputs = list(
    source_panel_manifest = file_record(manifest_path),
    source_panel_summary = file_record(source_summary_path),
    source_group_map = file_record(group_map_path),
    source_POS_braycurtis = file_record(pos_bray_path),
    source_NEG_braycurtis = file_record(neg_bray_path),
    ordered_display_order = file_record(order_path),
    substrate_summary = file_record(substrate_summary_path),
    substrate_manifest = file_record(substrate_manifest_path),
    style = file_record(style_path),
    producer = file_record(script_path),
    authority = lapply(authority_files, file_record)
  ),
  outputs = lapply(c(figure_paths, pair_path), file_record),
  started_at_utc = format(started, tz = "UTC", usetz = TRUE),
  completed_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE)
)

summary_path <- file.path(output_dir, "figure3c_strict_heatmap_review_only_summary.json")
write_json(summary, summary_path, pretty = TRUE, auto_unbox = TRUE)

report_path <- file.path(output_dir, "FIGURE3C_STRICT_HEATMAP_REVIEW_ONLY_REPORT.md")
report_lines <- c(
  "# Strict Figure 3c heatmap - review only",
  "",
  "- Status: `review_only_complete_pending_submitted_S1_reconciliation`",
  paste0("- Taxonomy release: `", TAXONOMY_RELEASE, "`"),
  "- Figure identity: coauthor-organized Main Figure 3c.",
  "- Definition: POS and NEG pairwise lipidome similarity matrices, calculated as `1 - Bray-Curtis`.",
  paste0("- Strict phyla: ", length(display_order), "; POS features/samples: 45,525/164; NEG features/samples: 14,896/192."),
  "- Order: strict Figure 3a ordered-view biological display order.",
  "- Diagonal: masked. Shared colour scale: 0 to 0.55; values above 0.55 are clipped.",
  "- Cross-group triangles and within-group circles are not used in this heatmap; those markers belong to Figure 3b.",
  "- Legacy 17-phylum inputs, internal v4.9 figure organization, and SSU tree distances were not used.",
  "- Figure 3, manuscript, legends, tables, response tracker, and submitted source files remain disconnected.",
  "- Direct cell-level reconciliation to the authoritative submitted Supplementary Table S1 workbook remains pending.",
  "",
  "## Review outputs",
  "",
  paste0("- `", basename(figure_paths[[1]]), "`"),
  paste0("- `", basename(pair_path), "`"),
  paste0("- `", basename(summary_path), "`"),
  ""
)
writeLines(report_lines, report_path, useBytes = TRUE)

output_paths <- c(figure_paths, pair_path, summary_path, report_path)
output_records <- lapply(output_paths, file_record)
stage_manifest <- list(
  schema_version = 1,
  stage_id = "figure3c_strict_heatmap_review_only",
  status = "review_only_complete_pending_submitted_S1_reconciliation",
  taxonomy_release = TAXONOMY_RELEASE,
  generated_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  outputs = output_records,
  gates = summary$gates,
  downstream_not_run = c(
    "manuscript/legend/table/response update",
    "submitted-S1 cell-level reconciliation",
    "submission replacement"
  )
)
write_json(stage_manifest, file.path(output_dir, "stage_manifest.json"), pretty = TRUE, auto_unbox = TRUE)

cat("Strict Figure 3c heatmap review-only outputs written to:", output_dir, "\n")
