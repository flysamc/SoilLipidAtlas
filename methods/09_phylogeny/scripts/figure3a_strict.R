#!/usr/bin/env Rscript
# Strict Figure 3a: POS and NEG UPGMA dendrograms from the locked 16-phylum release.

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  library(scales)
  library(ggdendro)
  library(ggtext)
  library(ragg)
  library(jsonlite)
  library(digest)
})

full_command <- commandArgs(trailingOnly = FALSE)
file_argument <- grep("^--file=", full_command, value = TRUE)
if (length(file_argument) != 1) stop("Cannot resolve the Figure 3a R producer path")
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
required_args <- c("substrate-dir", "taxonomy-dir", "output-dir", "style")
missing_args <- setdiff(required_args, names(args))
if (length(missing_args) > 0) stop("Missing arguments: ", paste(missing_args, collapse = ", "))

substrate_dir <- normalizePath(args[["substrate-dir"]], mustWork = TRUE)
taxonomy_dir <- normalizePath(args[["taxonomy-dir"]], mustWork = TRUE)
style_path <- normalizePath(args[["style"]], mustWork = TRUE)
output_dir <- args[["output-dir"]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
output_dir <- normalizePath(output_dir, mustWork = TRUE)
started <- Sys.time()

source(style_path)

DISPLAY_COLOURS <- c(
  Bacteria = unname(KINGDOM_COLOURS[["Bacteria"]]),
  Archaea = unname(KINGDOM_COLOURS[["Archaea"]]),
  Fungi = unname(KINGDOM_COLOURS[["Fungi"]]),
  Viridiplantae = unname(KINGDOM_COLOURS[["Plantae"]]),
  Animalia = unname(KINGDOM_COLOURS[["Animalia"]]),
  Protists = unname(KINGDOM_COLOURS[["Protozoa"]])
)
DISPLAY_ORDER <- names(DISPLAY_COLOURS)
MIXED_COLOUR <- "grey55"

sha256 <- function(path) digest(file = path, algo = "sha256", serialize = FALSE)

file_record <- function(path) {
  absolute <- normalizePath(path, mustWork = TRUE)
  list(
    path = absolute,
    bytes = unname(file.info(absolute)$size),
    sha256 = sha256(absolute)
  )
}

taxonomy_path <- file.path(taxonomy_dir, "taxonomy_summary.json")
taxonomy <- fromJSON(taxonomy_path, simplifyVector = TRUE)
if (!identical(taxonomy$assertions$status, "pass")) stop("Taxonomy assertions did not pass")
release_id <- taxonomy$taxonomy_release
phyla <- as.character(taxonomy$analysis_phyla)
if (length(phyla) != 16 || anyDuplicated(phyla)) stop("Figure 3a requires exactly 16 unique phyla")

policy_path <- file.path(dirname(style_path), "..", "..", "..", "paper2_repro", "config", "taxonomy_policy.json")
policy_path <- normalizePath(policy_path, mustWork = TRUE)
policy <- fromJSON(policy_path, simplifyVector = TRUE)
if (!identical(policy$release_id, release_id) || !identical(policy$status, "locked")) {
  stop("Taxonomy policy is not the locked Figure 3 release")
}

substrate_manifest_path <- file.path(substrate_dir, "stage_manifest.json")
substrate_manifest <- fromJSON(substrate_manifest_path, simplifyDataFrame = TRUE)
if (!identical(substrate_manifest$status, "pass")) stop("Figure 3 substrate is not passing")
if (!identical(substrate_manifest$taxonomy_release, release_id)) stop("Substrate release mismatch")
for (index in seq_len(nrow(substrate_manifest$outputs))) {
  record <- substrate_manifest$outputs[index, ]
  observed <- file_record(record$path)
  if (observed$bytes != record$bytes || !identical(observed$sha256, record$sha256)) {
    stop("Substrate checksum mismatch: ", record$path)
  }
}

substrate_summary_path <- file.path(substrate_dir, "substrate_summary.json")
substrate_summary <- fromJSON(substrate_summary_path, simplifyVector = TRUE)
if (!identical(as.character(substrate_summary$analysis_phyla), phyla)) {
  stop("Substrate and taxonomy phylum order differ")
}

labels_path <- file.path(substrate_dir, "sample_labels.csv")
labels <- read.csv(labels_path, check.names = FALSE, stringsAsFactors = FALSE)
if (anyDuplicated(labels[c("mode", "sample_id")])) stop("Duplicate mode/sample label rows")
if (!setequal(unique(labels$phylum), phyla)) stop("Sample-label phyla differ from strict analysis set")

group_rows <- unique(labels[c("phylum", "ecological_group")])
if (anyDuplicated(group_rows$phylum)) stop("A phylum maps to multiple display groups")
pk <- setNames(group_rows$ecological_group, group_rows$phylum)
if (!setequal(names(pk), phyla)) stop("Display-group map does not cover the strict phyla")
if (length(setdiff(unique(pk), DISPLAY_ORDER)) > 0) stop("Unknown display group in sample labels")

bray_curtis <- function(profiles) {
  n <- nrow(profiles)
  result <- matrix(0, nrow = n, ncol = n, dimnames = list(rownames(profiles), rownames(profiles)))
  for (i in seq_len(n - 1)) {
    for (j in seq.int(i + 1, n)) {
      denominator <- sum(abs(profiles[i, ] + profiles[j, ]))
      if (!is.finite(denominator) || denominator <= 0) stop("Invalid Bray-Curtis denominator")
      value <- sum(abs(profiles[i, ] - profiles[j, ])) / denominator
      result[i, j] <- value
      result[j, i] <- value
    }
  }
  result
}

freeze_mode <- function(mode, filename) {
  source_path <- file.path(substrate_dir, filename)
  message("[", mode, "] reading ", source_path)
  table <- read.csv(source_path, check.names = FALSE, stringsAsFactors = FALSE)
  feature_id <- table[[1]]
  if (anyDuplicated(feature_id)) stop(mode, ": duplicate feature IDs")
  mode_labels <- labels[labels$mode == mode, , drop = FALSE]
  missing_samples <- setdiff(mode_labels$sample_id, names(table))
  if (length(missing_samples) > 0) stop(mode, ": substrate samples missing from table")

  sample_counts <- setNames(integer(length(phyla)), phyla)
  profiles <- matrix(NA_real_, nrow = length(phyla), ncol = nrow(table), dimnames = list(phyla, feature_id))
  for (phylum in phyla) {
    sample_ids <- mode_labels$sample_id[mode_labels$phylum == phylum]
    if (length(sample_ids) < 2) stop(mode, ": ", phylum, " has fewer than two samples")
    sample_counts[[phylum]] <- length(sample_ids)
    profiles[phylum, ] <- rowMeans(as.matrix(table[sample_ids]), na.rm = FALSE)
  }
  if (any(!is.finite(profiles))) stop(mode, ": non-finite centroid intensity")

  bc <- bray_curtis(profiles)
  if (!isTRUE(all.equal(bc, t(bc), tolerance = 1e-14))) stop(mode, ": asymmetric Bray-Curtis matrix")
  if (any(diag(bc) != 0)) stop(mode, ": non-zero Bray-Curtis diagonal")
  hc <- hclust(as.dist(bc), method = "average")
  coph <- as.matrix(cophenetic(hc))
  coph_r <- unname(cor(as.dist(bc), as.dist(coph), method = "pearson"))
  if (!is.finite(coph_r)) stop(mode, ": invalid cophenetic correlation")

  centroid_table <- data.frame(feature_id = feature_id, t(profiles), check.names = FALSE)
  centroid_path <- file.path(output_dir, paste0(tolower(mode), "_phylum_centroids.csv"))
  bc_path <- file.path(output_dir, paste0(tolower(mode), "_braycurtis.csv"))
  coph_path <- file.path(output_dir, paste0(tolower(mode), "_cophenetic.csv"))
  linkage_path <- file.path(output_dir, paste0(tolower(mode), "_upgma_linkage.csv"))
  write.csv(centroid_table, centroid_path, row.names = FALSE, quote = FALSE)
  write.csv(bc, bc_path, row.names = TRUE, quote = FALSE)
  write.csv(coph, coph_path, row.names = TRUE, quote = FALSE)

  resolve_child <- function(value) {
    if (value < 0) phyla[-value] else paste0("cluster_", length(phyla) + value)
  }
  linkage <- data.frame(
    cluster_id = paste0("cluster_", length(phyla) + seq_len(nrow(hc$merge))),
    left_child = vapply(hc$merge[, 1], resolve_child, character(1)),
    right_child = vapply(hc$merge[, 2], resolve_child, character(1)),
    bray_curtis_height = hc$height,
    n_phyla = NA_integer_,
    stringsAsFactors = FALSE
  )
  get_leaves <- function(node) {
    if (node < 0) return(-node)
    c(get_leaves(hc$merge[node, 1]), get_leaves(hc$merge[node, 2]))
  }
  linkage$n_phyla <- vapply(seq_len(nrow(hc$merge)), function(i) length(get_leaves(i)), integer(1))
  write.csv(linkage, linkage_path, row.names = FALSE, quote = FALSE)

  list(
    mode = mode,
    source_path = source_path,
    n_features = nrow(table),
    n_samples = sum(sample_counts),
    sample_counts = sample_counts,
    profiles = profiles,
    bc = bc,
    hc = hc,
    cophenetic = coph,
    cophenetic_r = coph_r,
    leaf_order = hc$labels[hc$order],
    outputs = c(centroid_path, bc_path, coph_path, linkage_path)
  )
}

pos <- freeze_mode("POS", "pos.csv")
neg <- freeze_mode("NEG", "neg.csv")

group_map <- data.frame(
  phylum = phyla,
  ecological_group = unname(pk[phyla]),
  POS_samples = unname(pos$sample_counts[phyla]),
  NEG_samples = unname(neg$sample_counts[phyla]),
  stringsAsFactors = FALSE
)
group_map_path <- file.path(output_dir, "phylum_display_groups_and_counts.csv")
write.csv(group_map, group_map_path, row.names = FALSE, quote = FALSE)

build_dendrogram <- function(result, title_text, tag = NULL, show_legend = FALSE) {
  hc <- result$hc
  dd <- dendro_data(hc, type = "rectangle")
  seg <- segment(dd)
  lab <- label(dd)
  lab$ecological_group <- unname(pk[as.character(lab$label)])
  lab$col <- unname(DISPLAY_COLOURS[lab$ecological_group])

  leaf_group <- pk[hc$labels]
  get_leaves <- function(node) {
    if (node < 0) return(-node)
    c(get_leaves(hc$merge[node, 1]), get_leaves(hc$merge[node, 2]))
  }
  merge_col <- vapply(seq_len(nrow(hc$merge)), function(i) {
    groups <- unique(leaf_group[get_leaves(i)])
    if (length(groups) == 1) unname(DISPLAY_COLOURS[[groups]]) else MIXED_COLOUR
  }, character(1))

  seg$col <- MIXED_COLOUR
  for (i in seq_len(nrow(seg))) {
    at_y <- which(abs(hc$height - seg$y[i]) < 1e-10)
    if (length(at_y) == 1) {
      seg$col[i] <- merge_col[at_y]
      next
    }
    at_yend <- which(abs(hc$height - seg$yend[i]) < 1e-10)
    if (length(at_yend) == 1) {
      seg$col[i] <- merge_col[at_yend]
      next
    }
    if (seg$yend[i] == 0) {
      leaf <- which(abs(lab$x - seg$xend[i]) < 0.5)
      if (length(leaf) == 1) seg$col[i] <- lab$col[leaf]
    }
  }

  plot <- ggplot() +
    geom_segment(
      data = seg,
      aes(x = x, y = y, xend = xend, yend = yend, colour = col),
      linewidth = 0.4
    ) +
    geom_text(
      data = lab,
      aes(x = x, y = -0.018, label = label, colour = col),
      size = pt2mm(5), angle = 90, hjust = 1,
      show.legend = FALSE
    ) +
    scale_colour_identity(
      name = NULL,
      breaks = unname(DISPLAY_COLOURS),
      labels = names(DISPLAY_COLOURS),
      guide = if (show_legend) "legend" else "none"
    ) +
    scale_y_continuous(
      name = "Bray-Curtis distance",
      breaks = seq(0, 1, 0.25),
      expand = expansion(mult = c(0, 0.03))
    ) +
    scale_x_continuous(expand = expansion(add = 0.6)) +
    coord_cartesian(ylim = c(-0.20, max(hc$height) * 1.05), clip = "off") +
    labs(
      title = sprintf("%s  (cophenetic r = %.3f)", title_text, result$cophenetic_r),
      tag = tag
    ) +
    theme_nature() +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.x = element_blank(),
      axis.line.x = element_blank(),
      plot.title = element_text(size = 6),
      plot.tag = element_text(size = NC_TAG_PT, face = "bold"),
      plot.tag.position = c(0.005, 1.0),
      plot.margin = margin(2, 2, 11, 2, "mm")
    ) +
    guides(colour = guide_legend(
      nrow = 2,
      byrow = TRUE,
      override.aes = list(linewidth = 1.2, size = 3, angle = 0)
    ))
  plot
}

pa_pos <- build_dendrogram(pos, "Positive mode", "a", show_legend = TRUE)
pa_neg <- build_dendrogram(neg, "Negative mode", show_legend = FALSE)
caption <- sprintf(
  paste0(
    "Strict NCBI phyla (n = 16).\n",
    "POS: %s features, %d samples; NEG: %s features, %d samples.\n",
    "UPGMA average linkage on phylum-centroid Bray-Curtis distances."
  ),
  format(pos$n_features, big.mark = ",", scientific = FALSE), pos$n_samples,
  format(neg$n_features, big.mark = ",", scientific = FALSE), neg$n_samples
)
figure <- (pa_pos / pa_neg / guide_area()) +
  plot_layout(guides = "collect", heights = c(1, 1, 0.12)) +
  plot_annotation(
    caption = caption,
    theme = theme(plot.caption = element_text(size = 5, hjust = 0.5, margin = margin(t = 2, unit = "mm")))
  ) &
  theme(legend.position = "bottom", legend.box.margin = margin(0, 0, 0, 0))

figure_base <- file.path(output_dir, "Figure_3a_strict_16phyla")
figure_paths <- paste0(figure_base, c(".png", ".pdf", ".svg"))
ggsave(figure_paths[[1]], figure, width = NC_SINGLE * MM_TO_IN, height = 150 * MM_TO_IN,
       dpi = 600, bg = "white", device = ragg::agg_png)
ggsave(figure_paths[[2]], figure, width = NC_SINGLE * MM_TO_IN, height = 150 * MM_TO_IN,
       bg = "white", device = grDevices::cairo_pdf)
ggsave(figure_paths[[3]], figure, width = NC_SINGLE * MM_TO_IN, height = 150 * MM_TO_IN,
       bg = "white", device = grDevices::svg)

software <- list(
  R = R.version.string,
  ggplot2 = as.character(packageVersion("ggplot2")),
  dplyr = as.character(packageVersion("dplyr")),
  patchwork = as.character(packageVersion("patchwork")),
  ggdendro = as.character(packageVersion("ggdendro")),
  ragg = as.character(packageVersion("ragg")),
  jsonlite = as.character(packageVersion("jsonlite")),
  digest = as.character(packageVersion("digest"))
)

summary <- list(
  taxonomy_release = release_id,
  stage_id = "figure3a",
  scope = "Figure 3a only; no taxonomic-tree distance, Mantel, partial Mantel, panel b or panel c",
  status = "candidate_visual_review",
  chart_contract = list(
    analytical_question = "How do strict-release phylum-centroid lipidomes cluster independently in POS and NEG?",
    family = "hierarchy",
    variant = "two stacked rectangular UPGMA dendrograms",
    renderer = "R hclust plus ggdendro/ggplot2/patchwork",
    palette_policy = "six display strata with named coloured tips; mixed branches neutral",
    non_colour_identification = "every terminal tip is named",
    surface = "89 mm standalone Figure 3a candidate"
  ),
  n_phyla = length(phyla),
  n_pairs_per_mode = choose(length(phyla), 2),
  analysis_phyla = phyla,
  display_groups_only = as.list(pk[phyla]),
  POS = list(
    n_features = pos$n_features,
    n_samples = pos$n_samples,
    sample_counts = as.list(pos$sample_counts),
    cophenetic_r = pos$cophenetic_r,
    leaf_order = pos$leaf_order
  ),
  NEG = list(
    n_features = neg$n_features,
    n_samples = neg$n_samples,
    sample_counts = as.list(neg$sample_counts),
    cophenetic_r = neg$cophenetic_r,
    leaf_order = neg$leaf_order
  ),
  font = list(
    family = "Helvetica as requested by soilmass_style.R",
    status = "candidate; final PDF font embedding remains a release gate"
  ),
  software = software,
  elapsed_seconds = unname(as.numeric(difftime(Sys.time(), started, units = "secs")))
)
summary_path <- file.path(output_dir, "figure3a_summary.json")
write_json(summary, summary_path, pretty = TRUE, auto_unbox = TRUE, digits = 16)

input_paths <- list(
  taxonomy_summary = taxonomy_path,
  taxonomy_policy = policy_path,
  substrate_manifest = substrate_manifest_path,
  substrate_summary = substrate_summary_path,
  sample_labels = labels_path,
  pos_substrate = file.path(substrate_dir, "pos.csv"),
  neg_substrate = file.path(substrate_dir, "neg.csv"),
  renderer = script_path,
  style = style_path
)
output_paths <- c(pos$outputs, neg$outputs, group_map_path, figure_paths, summary_path)
manifest <- list(
  schema_version = 1,
  stage_id = "figure3a",
  taxonomy_release = release_id,
  status = "candidate_visual_review",
  generated_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  inputs = lapply(input_paths, file_record),
  outputs = lapply(output_paths, file_record),
  substrate_stage_sha256 = sha256(substrate_manifest_path),
  substrate_stage_status = substrate_manifest$status,
  downstream_not_run = c(
    "taxonomic-tree distance",
    "Mantel and partial Mantel tests",
    "Figure 3b and Figure 3c",
    "manuscript, legend, supplementary table and response updates"
  )
)
manifest_path <- file.path(output_dir, "stage_manifest.json")
write_json(manifest, manifest_path, pretty = TRUE, auto_unbox = TRUE, digits = 16)

cat(toJSON(list(
  status = summary$status,
  n_phyla = summary$n_phyla,
  n_pairs_per_mode = summary$n_pairs_per_mode,
  POS = summary$POS,
  NEG = summary$NEG,
  figure = figure_paths[[1]]
), pretty = TRUE, auto_unbox = TRUE, digits = 16), "\n")
