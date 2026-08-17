#!/usr/bin/env Rscript
# Assemble the complete coauthor-organized Main Figure 3 from the strict-v1
# R-rendered panel candidates: 3a dendrograms, 3b Mantel scatter, and 3c
# POS/NEG similarity heatmaps. Panel PNGs are reused only as rendered outputs;
# no legacy numerical input is read here.

suppressPackageStartupMessages({
  library(magick)
  library(jsonlite)
  library(digest)
})

full_command <- commandArgs(trailingOnly = FALSE)
file_argument <- grep("^--file=", full_command, value = TRUE)
if (length(file_argument) != 1) stop("Cannot resolve the Figure 3 composite R producer path")
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
required_args <- c("panel-a", "panel-b", "panel-c", "output-dir")
missing_args <- setdiff(required_args, names(args))
if (length(missing_args) > 0) stop("Missing arguments: ", paste(missing_args, collapse = ", "))

panel_a <- normalizePath(args[["panel-a"]], mustWork = TRUE)
panel_b <- normalizePath(args[["panel-b"]], mustWork = TRUE)
panel_c <- normalizePath(args[["panel-c"]], mustWork = TRUE)
output_dir <- args[["output-dir"]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
output_dir <- normalizePath(output_dir, mustWork = TRUE)
started <- Sys.time()

TAXONOMY_RELEASE <- "ncbi-phylum-2026-08-04-v1"
TREE_FREEZE_ID <- "figure3-ssu-curated-freeze-2026-08-04-v3"
PRIMARY_SSU_METRIC <- "anchor_set_mean__inventory_weighted"
EXPECTED_PANEL_B <- "Figure_3b_ssu_evolutionary_distance_review_only.png"
WIDTH_MM <- 183
HEIGHT_MM <- 195
MM_TO_IN <- 1 / 25.4

if (!identical(basename(panel_b), EXPECTED_PANEL_B)) {
  stop(
    "Figure 3b must be the finalized SSU evolutionary-distance panel: ",
    EXPECTED_PANEL_B
  )
}
ssu_summary_path <- file.path(
  dirname(panel_b), "figure3_ssu_distance_v3_review_only_summary.json"
)
if (!file.exists(ssu_summary_path)) stop("Missing SSU Figure 3b provenance summary")
ssu_summary <- fromJSON(ssu_summary_path, simplifyVector = TRUE)
if (!identical(ssu_summary$taxonomy_release, TAXONOMY_RELEASE) ||
    !identical(ssu_summary$tree_freeze_id, TREE_FREEZE_ID) ||
    !identical(ssu_summary$primary_metric, PRIMARY_SSU_METRIC)) {
  stop("Figure 3b SSU release, tree freeze, or primary metric does not match")
}

sha256 <- function(path) digest(file = path, algo = "sha256", serialize = FALSE)
file_record <- function(path) {
  absolute <- normalizePath(path, mustWork = TRUE)
  list(path = absolute, bytes = unname(file.info(absolute)$size), sha256 = sha256(absolute))
}

# The standalone strict 3a render includes a legend/caption below the two
# dendrograms. The manuscript composite puts the dendrograms in the upper-left
# slot, so trim only that trailing packaging region while preserving all branch
# and tip labels. The strict 3b and 3c renders are already composite-shaped.
img_a <- magick::image_read(panel_a) |>
  magick::image_crop(geometry = "2102x3100+0+0")
img_b <- magick::image_read(panel_b)
img_c <- magick::image_read(panel_c)

grob_a <- grid::rasterGrob(as.raster(img_a), interpolate = TRUE)
grob_b <- grid::rasterGrob(as.raster(img_b), interpolate = TRUE)
grob_c <- grid::rasterGrob(as.raster(img_c), interpolate = TRUE)

figure_base <- file.path(output_dir, "Figure_3_strict_review_only")
figure_paths <- paste0(figure_base, c(".png", ".pdf", ".svg"))

# Draw the three already-formatted R panel candidates on one manuscript-sized
# canvas. A raster viewport is used only for composition; the panel formatting
# itself remains produced by the existing strict R panel code.
draw_grob <- function(grob, vp) {
  grid::pushViewport(vp)
  grid::grid.draw(grob)
  grid::popViewport()
}

draw_one_page_clean <- function(path, kind) {
  width_in <- WIDTH_MM * MM_TO_IN
  height_in <- HEIGHT_MM * MM_TO_IN
  if (kind == "png") {
    grDevices::png(path, width = width_in, height = height_in, units = "in", res = 600, bg = "white")
  } else if (kind == "pdf") {
    grDevices::cairo_pdf(path, width = width_in, height = height_in, bg = "white")
  } else {
    grDevices::svg(path, width = width_in, height = height_in, bg = "white")
  }
  grid::grid.newpage()
  top_height <- 103 / HEIGHT_MM
  bottom_height <- 92 / HEIGHT_MM
  draw_grob(grob_a, grid::viewport(
    x = grid::unit(0, "npc"), y = grid::unit(bottom_height, "npc"),
    width = grid::unit(89 / WIDTH_MM, "npc"), height = grid::unit(top_height, "npc"),
    just = c("left", "bottom")
  ))
  draw_grob(grob_b, grid::viewport(
    x = grid::unit(89 / WIDTH_MM, "npc"), y = grid::unit(bottom_height, "npc"),
    width = grid::unit(94 / WIDTH_MM, "npc"), height = grid::unit(top_height, "npc"),
    just = c("left", "bottom")
  ))
  draw_grob(grob_c, grid::viewport(
    x = grid::unit(0, "npc"), y = grid::unit(0, "npc"),
    width = grid::unit(1, "npc"), height = grid::unit(bottom_height, "npc"),
    just = c("left", "bottom")
  ))
  grDevices::dev.off()
}

draw_one_page_clean(figure_paths[[1]], "png")
draw_one_page_clean(figure_paths[[2]], "pdf")
draw_one_page_clean(figure_paths[[3]], "svg")

authority_files <- c(
  file.path("submission_source", "SLA-Main-article-nature-comms.docx"),
  file.path("submission_source", "SLA-Supplementary-nature-comms.docx"),
  file.path("submission_source", "SLA-figures-and-legends-nature-comms.docx"),
  file.path("paper2_repro", "portable.py")
)
authority_files <- authority_files[file.exists(authority_files)]

summary <- list(
  schema_version = 1,
  stage_id = "figure3_strict_composite_review_only",
  status = "review_only_complete_pending_submitted_S1_reconciliation",
  taxonomy_release = TAXONOMY_RELEASE,
  tree_freeze_id = TREE_FREEZE_ID,
  panel_b_primary_metric = PRIMARY_SSU_METRIC,
  figure_identity = "Coauthor-organized Main Figure 3: 3a dendrogram, 3b Mantel scatter, 3c POS/NEG similarity heatmaps.",
  layout = list(width_mm = WIDTH_MM, height_mm = HEIGHT_MM, top_row_mm = 103, bottom_row_mm = 92),
  inputs = list(
    panel_a = file_record(panel_a),
    panel_b = file_record(panel_b),
    panel_b_ssu_summary = file_record(ssu_summary_path),
    panel_c = file_record(panel_c),
    producer = file_record(script_path),
    authority = lapply(authority_files, file_record)
  ),
  gates = list(
    taxonomy_release_locked = TRUE,
    strict_v1_panels_only = TRUE,
    coauthor_figure_mapping_used = TRUE,
    panel_b_ssu_v3_used = TRUE,
    ncbi_hierarchy_panel_used = FALSE,
    v4_9_figure_organization_used = FALSE,
    legacy_17_phylum_panel_used = FALSE,
    figure3_connected = FALSE,
    submitted_S1_cell_reconciliation = FALSE
  ),
  outputs = lapply(figure_paths, file_record),
  started_at_utc = format(started, tz = "UTC", usetz = TRUE),
  completed_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE)
)

summary_path <- file.path(output_dir, "figure3_strict_composite_review_only_summary.json")
write_json(summary, summary_path, pretty = TRUE, auto_unbox = TRUE)

report_path <- file.path(output_dir, "FIGURE3_STRICT_COMPOSITE_REVIEW_ONLY_REPORT.md")
writeLines(c(
  "# Strict Main Figure 3 composite - review only",
  "",
  "- 3a: strict POS/NEG UPGMA dendrograms.",
  "- 3b: finalized v3 SSU evolutionary-distance Mantel scatter; triangles are cross-group and circles are within-group.",
  "- 3c: strict POS/NEG pairwise similarity heatmaps.",
  "- Layout: 183 x 195 mm, matching the manuscript double-column composite.",
  "- No manuscript, supplementary, response, table, or submission file was changed.",
  "- Direct submitted-S1 cell-level reconciliation remains pending.",
  ""
), report_path, useBytes = TRUE)

output_paths <- c(figure_paths, summary_path, report_path)
stage_manifest <- list(
  schema_version = 1,
  stage_id = "figure3_strict_composite_review_only",
  status = "review_only_complete_pending_submitted_S1_reconciliation",
  taxonomy_release = TAXONOMY_RELEASE,
  tree_freeze_id = TREE_FREEZE_ID,
  inputs = list(
    panel_b = file_record(panel_b),
    panel_b_ssu_summary = file_record(ssu_summary_path)
  ),
  generated_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  outputs = lapply(output_paths, file_record),
  gates = summary$gates,
  downstream_not_run = c("manuscript/legend/table/response update", "submission replacement")
)
write_json(stage_manifest, file.path(output_dir, "stage_manifest.json"), pretty = TRUE, auto_unbox = TRUE)
cat("Strict Figure 3 composite written to:", output_dir, "\n")
