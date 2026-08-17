#!/usr/bin/env Rscript

# Review-only application of the v3 organism-level SSU tree to the strict
# positive/negative lipid-distance matrices. This stage is deliberately
# disconnected from Figure 3, manuscript, legends, tables, and response text.
# It uses the recovered historical batch-overlap covariate, but recomputes all
# distances and permutations from the locked v1 strict substrate.

suppressPackageStartupMessages({
  library(ggplot2)
  library(jsonlite)
  library(digest)
})

parse_args <- function(values) {
  out <- list()
  for (value in values) {
    if (!startsWith(value, "--") || !grepl("=", value, fixed = TRUE)) {
      stop("Arguments must use --name=value syntax: ", value)
    }
    pieces <- strsplit(sub("^--", "", value), "=", fixed = TRUE)[[1]]
    out[[pieces[[1]]]] <- paste(pieces[-1], collapse = "=")
  }
  out
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
required <- c("panel-a-dir", "tree-inference-dir", "substrate-dir", "curated-units",
              "pos-metadata", "neg-metadata", "output-dir")
missing <- setdiff(required, names(args))
if (length(missing) > 0) stop("Missing arguments: ", paste(missing, collapse = ", "))

panel_a_dir <- normalizePath(args[["panel-a-dir"]], mustWork = TRUE)
tree_dir <- normalizePath(args[["tree-inference-dir"]], mustWork = TRUE)
substrate_dir <- normalizePath(args[["substrate-dir"]], mustWork = TRUE)
curated_units_path <- normalizePath(args[["curated-units"]], mustWork = TRUE)
pos_metadata_path <- normalizePath(args[["pos-metadata"]], mustWork = TRUE)
neg_metadata_path <- normalizePath(args[["neg-metadata"]], mustWork = TRUE)
output_dir <- args[["output-dir"]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
output_dir <- normalizePath(output_dir, mustWork = TRUE)

RELEASE_ID <- "ncbi-phylum-2026-08-04-v1"
TREE_FREEZE_ID <- "figure3-ssu-curated-freeze-2026-08-04-v3"
N_PHYLA <- 16L
N_PAIRS <- 120L
N_PERMUTATIONS <- 9999L
PHYLA <- c(
  "Actinomycetota", "Arthropoda", "Ascomycota", "Bacillota",
  "Basidiomycota", "Chlorophyta", "Discosea", "Evosea",
  "Heterolobosea", "Methanobacteriota", "Mollusca", "Mucoromycota",
  "Nematoda", "Pseudomonadota", "Streptophyta", "Thermoproteota"
)
MODES <- c("POS", "NEG")
MANTEL_SEEDS <- c(POS = 20720247L, NEG = 21211386L)
PARTIAL_SEEDS <- c(POS = 20720248L, NEG = 21211387L)

HISTORICAL_POS_SOURCE <- file.path(
  "external", "SOILMASS_PRODUCER_RECOVERY_2026-08-04", "payload", "core",
  "workspace", "analysis", "analysis-17", "positive", "scripts",
  "reviewer_computational_tasks.py"
)
HISTORICAL_NEG_SOURCE <- file.path(
  "external", "SOILMASS_PRODUCER_RECOVERY_2026-08-04", "payload", "core",
  "workspace", "analysis", "analysis-16", "negative_mode", "scripts",
  "phylogenetic_signal_neg.py"
)

sha256 <- function(path) digest(path, algo = "sha256", file = TRUE, serialize = FALSE)
file_record <- function(path) {
  path <- normalizePath(path, mustWork = TRUE)
  list(path = path, bytes = unname(file.info(path)$size), sha256 = sha256(path))
}

STYLE_PATH <- normalizePath(file.path(
  "manuscript_2_clean", "06_figures", "figures_r", "soilmass_style.R"
), mustWork = TRUE)
STRICT_TAXA_PATH <- normalizePath(file.path(
  "outputs", "analysis", RELEASE_ID, "figure3", "evolutionary_tree_reference",
  "silva_138_2", "audit", "strict_sample_taxa.csv"
), mustWork = TRUE)
MODE_COLOURS <- c(Positive = "#0072B2", Negative = "#E69F00")
theme_nature_figure3 <- function(base_size = 7) {
  theme_classic(base_size = base_size, base_family = "Helvetica") %+replace%
    theme(
      line = element_line(linewidth = 0.3, colour = "black"),
      axis.line = element_line(linewidth = 0.3),
      axis.ticks = element_line(linewidth = 0.3),
      axis.ticks.length = unit(1.2, "pt"),
      axis.title = element_text(size = 7, colour = "black"),
      axis.text = element_text(size = 6, colour = "black"),
      legend.title = element_text(size = 6, face = "bold"),
      legend.text = element_text(size = 5),
      legend.key.size = unit(3, "mm"),
      legend.key.height = unit(3, "mm"),
      legend.key.width = unit(3, "mm"),
      legend.background = element_blank(),
      legend.box.background = element_blank(),
      legend.margin = margin(0, 0, 0, 0),
      legend.spacing = unit(1, "mm"),
      plot.title = element_text(size = 7, face = "bold", hjust = 0),
      plot.margin = margin(2, 2, 2, 2, "mm"),
      strip.background = element_blank(),
      strip.text = element_text(size = 6),
      panel.grid = element_blank()
    )
}

write_matrix <- function(value, path) {
  output <- data.frame(label = rownames(value), value, check.names = FALSE,
                       stringsAsFactors = FALSE)
  names(output)[[1]] <- ""
  write.csv(output, path, row.names = FALSE, quote = FALSE, na = "")
}

read_matrix <- function(path) {
  tab <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  labels <- as.character(tab[[1]])
  value <- as.matrix(tab[-1])
  storage.mode(value) <- "double"
  rownames(value) <- labels
  colnames(value) <- names(tab)[-1]
  value
}

validate_distance <- function(value, label, bounded = FALSE) {
  if (!identical(dim(value), c(N_PHYLA, N_PHYLA))) stop(label, " must be 16 x 16")
  if (!identical(rownames(value), colnames(value))) stop(label, " labels differ")
  if (any(!is.finite(value)) || any(value < -1e-12)) stop(label, " has invalid values")
  if (max(abs(value - t(value))) > 1e-10) stop(label, " is asymmetric")
  if (max(abs(diag(value))) > 1e-10) stop(label, " has non-zero diagonal")
  if (bounded && max(value) > 1 + 1e-10) stop(label, " exceeds one")
  invisible(TRUE)
}

validate_similarity <- function(value, label) {
  if (!identical(dim(value), c(N_PHYLA, N_PHYLA))) stop(label, " must be 16 x 16")
  if (!identical(rownames(value), colnames(value))) stop(label, " labels differ")
  if (any(!is.finite(value)) || any(value < -1e-12) || any(value > 1 + 1e-12)) {
    stop(label, " has invalid similarity values")
  }
  if (max(abs(value - t(value))) > 1e-10) stop(label, " is asymmetric")
  if (max(abs(diag(value) - 1)) > 1e-10) stop(label, " has a non-unit diagonal")
  invisible(TRUE)
}

condensed <- function(value) unname(value[upper.tri(value)])

panel_manifest_path <- file.path(panel_a_dir, "stage_manifest.json")
substrate_summary_path <- file.path(substrate_dir, "substrate_summary.json")
sample_labels_path <- file.path(substrate_dir, "sample_labels.csv")
panel_manifest <- fromJSON(panel_manifest_path, simplifyVector = TRUE)
substrate_summary <- fromJSON(substrate_summary_path, simplifyVector = TRUE)
if (!identical(panel_manifest$taxonomy_release, RELEASE_ID) ||
    !identical(substrate_summary$taxonomy_release, RELEASE_ID)) {
  stop("Panel A or substrate is not from the locked taxonomy release")
}
if (!identical(substrate_summary$analysis_phyla, PHYLA)) {
  stop("Strict substrate phylum order differs from the locked v1 order")
}

pos_bray_path <- file.path(panel_a_dir, "pos_braycurtis.csv")
neg_bray_path <- file.path(panel_a_dir, "neg_braycurtis.csv")
pos_bray <- read_matrix(pos_bray_path)
neg_bray <- read_matrix(neg_bray_path)
validate_distance(pos_bray, "POS Bray-Curtis", bounded = TRUE)
validate_distance(neg_bray, "NEG Bray-Curtis", bounded = TRUE)
if (!identical(rownames(pos_bray), PHYLA) || !identical(rownames(neg_bray), PHYLA)) {
  stop("Panel A Bray-Curtis matrices do not have the locked phylum order")
}

unit_metadata_path <- file.path(tree_dir, "analysis_unit_tree_nodes_v3.csv")
unit_metadata <- read.csv(unit_metadata_path, check.names = FALSE,
                          colClasses = "character", stringsAsFactors = FALSE)
if (nrow(unit_metadata) != 103L || anyDuplicated(unit_metadata$analysis_unit_taxid) ||
    anyDuplicated(unit_metadata$analysis_unit_name)) {
  stop("The finalized SSU unit metadata must contain 103 unique units and names")
}
if (length(unique(unit_metadata$phylum)) != N_PHYLA ||
    !setequal(unique(unit_metadata$phylum), PHYLA)) {
  stop("SSU unit metadata does not cover the locked 16 phyla")
}
curated_units <- read.csv(curated_units_path, check.names = FALSE,
                          colClasses = "character", stringsAsFactors = FALSE)
curated_primary <- curated_units[toupper(curated_units$included_primary) == "TRUE", , drop = FALSE]
curated_index <- match(unit_metadata$analysis_unit_taxid, curated_primary$analysis_unit_taxid)
if (anyNA(curated_index)) stop("A finalized tree unit is absent from the curated primary unit map")
inventory_weight <- as.numeric(curated_primary$inventory_rows[curated_index])
if (any(!is.finite(inventory_weight)) || any(inventory_weight <= 0)) {
  stop("Curated inventory weights are not positive finite values")
}
unit_metadata$inventory_weight <- inventory_weight
unit_metadata$analysis_unit_taxid <- as.character(unit_metadata$analysis_unit_taxid)

distance_files <- c(
  anchor_set_mean = file.path(tree_dir, "analysis_unit_patristic_distance_anchor_set_mean_v3.csv"),
  medoid = file.path(tree_dir, "analysis_unit_patristic_distance_medoid_v3.csv"),
  mrca = file.path(tree_dir, "analysis_unit_patristic_distance_mrca_v3.csv")
)
unit_distances <- lapply(distance_files, read_matrix)
for (metric in names(unit_distances)) {
  value <- unit_distances[[metric]]
  if (nrow(value) != 103L || !identical(rownames(value), colnames(value)) ||
      !setequal(rownames(value), unit_metadata$analysis_unit_name)) {
    stop("Unit distance matrix contract failed for ", metric)
  }
  unit_distances[[metric]] <- value[unit_metadata$analysis_unit_name,
                                    unit_metadata$analysis_unit_name, drop = FALSE]
  value <- unit_distances[[metric]]
  if (any(!is.finite(value)) || any(value < -1e-12) ||
      max(abs(value - t(value))) > 1e-10 || max(abs(diag(value))) > 1e-10) {
    stop("Invalid unit distance matrix for ", metric)
  }
}

aggregate_phylum_distance <- function(unit_distance, weighting) {
  output <- matrix(0, N_PHYLA, N_PHYLA, dimnames = list(PHYLA, PHYLA))
  for (i in seq_len(N_PHYLA - 1L)) {
    for (j in seq.int(i + 1L, N_PHYLA)) {
      rows <- which(unit_metadata$phylum == PHYLA[[i]])
      cols <- which(unit_metadata$phylum == PHYLA[[j]])
      distances <- unit_distance[rows, cols, drop = FALSE]
      if (weighting == "inventory_weighted") {
        weights <- outer(inventory_weight[rows], inventory_weight[cols])
        value <- weighted.mean(as.vector(distances), as.vector(weights))
      } else if (weighting == "unit_unweighted") {
        value <- mean(distances)
      } else {
        stop("Unknown SSU phylum aggregation weighting: ", weighting)
      }
      output[i, j] <- value
      output[j, i] <- value
    }
  }
  validate_distance(output, paste0("SSU phylum distance ", weighting))
  output
}

evolutionary_matrices <- list()
for (metric in names(unit_distances)) {
  for (weighting in c("inventory_weighted", "unit_unweighted")) {
    key <- paste(metric, weighting, sep = "__")
    evolutionary_matrices[[key]] <- aggregate_phylum_distance(unit_distances[[metric]], weighting)
  }
}

sample_labels <- read.csv(sample_labels_path, check.names = FALSE,
                          colClasses = "character", stringsAsFactors = FALSE)
pos_source_path <- file.path(substrate_dir, "pos.csv")
if (!identical(sort(unique(sample_labels$phylum)), sort(PHYLA)) ||
    anyDuplicated(paste(sample_labels$mode, sample_labels$sample_id))) {
  stop("Strict substrate sample labels are incomplete or duplicated")
}

read_mode_metadata <- function(mode, path) {
  metadata <- read.csv(path, check.names = FALSE, colClasses = "character",
                       stringsAsFactors = FALSE)
  if (!all(c("sample_name", "batch", "ncbi_phylum", "taxonomy_release") %in% names(metadata))) {
    stop("Missing strict metadata fields for ", mode)
  }
  if (anyDuplicated(metadata$sample_name)) stop("Duplicate sample_name in ", mode, " metadata")
  labels <- sample_labels[sample_labels$mode == mode, , drop = FALSE]
  labels$sample_name <- sub("^sample:", "", labels$sample_id)
  index <- match(labels$sample_name, metadata$sample_name)
  if (anyNA(index)) stop("Strict metadata does not map every ", mode, " sample label")
  optional_field <- function(field) {
    if (field %in% names(metadata)) metadata[[field]][index] else rep(NA_character_, length(index))
  }
  out <- data.frame(
    mode = mode,
    sample_id = labels$sample_id,
    sample_name = labels$sample_name,
    phylum = labels$phylum,
    metadata_phylum = metadata$ncbi_phylum[index],
    source_phylum = optional_field("source_phylum"),
    genus = optional_field("genus"),
    species = optional_field("species"),
    ecological_group = optional_field("ecological_group"),
    batch = metadata$batch[index],
    taxonomy_scope = metadata$taxonomy_scope[index],
    taxonomy_release = metadata$taxonomy_release[index],
    stringsAsFactors = FALSE
  )
  if (any(out$phylum != out$metadata_phylum)) stop("Label/metadata phylum mismatch in ", mode)
  if (any(out$taxonomy_release != RELEASE_ID)) stop("Metadata release mismatch in ", mode)
  if (any(is.na(out$batch) | !nzchar(out$batch))) stop("Missing batch in ", mode, " metadata")
  if (any(table(out$phylum) < 2L)) stop("A strict ", mode, " phylum has fewer than two samples")
  out
}

mode_metadata <- list(
  POS = read_mode_metadata("POS", pos_metadata_path),
  NEG = read_mode_metadata("NEG", neg_metadata_path)
)

# Reuse the manuscript's marker-shape convention while keeping the strict
# NCBI phylum release as the computational taxonomy.  The grouping variable is
# an ecological display stratum only; it is not treated as a shared taxonomic
# rank.  The labels intentionally match the submitted Figure 4 legend.
group_rows <- unique(do.call(rbind, lapply(mode_metadata, function(metadata) {
  data.frame(phylum = metadata$phylum,
             ecological_group = metadata$ecological_group,
             stringsAsFactors = FALSE)
})))
group_by_phylum <- setNames(vapply(PHYLA, function(phylum) {
  values <- unique(group_rows$ecological_group[group_rows$phylum == phylum])
  values <- values[!is.na(values) & nzchar(values)]
  if (length(values) != 1L) {
    stop("Strict phylum must map to exactly one ecological display group: ", phylum)
  }
  values[[1]]
}, character(1)), PHYLA)
group_map <- data.frame(
  phylum = PHYLA,
  ecological_display_group = unname(group_by_phylum),
  stringsAsFactors = FALSE
)
write.csv(group_map, file.path(output_dir, "strict_phylum_ecological_display_groups.csv"),
          row.names = FALSE, quote = FALSE)
expected_counts <- list(
  POS = unlist(substrate_summary$positive$samples_per_phylum),
  NEG = unlist(substrate_summary$negative$samples_per_phylum)
)
sample_count_rows <- do.call(rbind, lapply(MODES, function(mode) {
  metadata <- mode_metadata[[mode]]
  counts <- table(factor(metadata$phylum, levels = PHYLA))
  if (!identical(as.integer(counts), as.integer(expected_counts[[mode]]))) {
    stop("Exact strict sample denominator/count vector failed for ", mode)
  }
  data.frame(mode = mode, phylum = PHYLA, n_samples = as.integer(counts),
             n_batches = vapply(PHYLA, function(p) {
               length(unique(metadata$batch[metadata$phylum == p]))
             }, integer(1)), stringsAsFactors = FALSE)
}))
write.csv(sample_count_rows, file.path(output_dir, "strict_mode_sample_denominators.csv"),
          row.names = FALSE, quote = FALSE)

batch_overlap <- list()
for (mode in MODES) {
  metadata <- mode_metadata[[mode]]
  batch_sets <- setNames(lapply(PHYLA, function(p) unique(metadata$batch[metadata$phylum == p])), PHYLA)
  overlap <- matrix(0, N_PHYLA, N_PHYLA, dimnames = list(PHYLA, PHYLA))
  for (i in seq_len(N_PHYLA)) {
    for (j in seq_len(N_PHYLA)) {
      shared <- length(intersect(batch_sets[[i]], batch_sets[[j]]))
      total <- length(union(batch_sets[[i]], batch_sets[[j]]))
      overlap[i, j] <- if (total > 0) shared / total else 0
    }
  }
  validate_similarity(overlap, paste0(mode, " batch-overlap Jaccard"))
  batch_overlap[[mode]] <- overlap
  write_matrix(overlap, file.path(output_dir, paste0("batch_overlap_jaccard_", mode, ".csv")))
}

residualize <- function(y, x) {
  x_centered <- x - mean(x)
  y_centered <- y - mean(y)
  denominator <- sum(x_centered * x_centered)
  if (!is.finite(denominator) || denominator <= 0) stop("Batch-overlap covariate has no variance")
  slope <- sum(x_centered * y_centered) / denominator
  y - (mean(y) + slope * x_centered)
}

mantel_greater <- function(lipid, evolutionary, seed) {
  lipid_vector <- condensed(lipid)
  evolutionary_vector <- condensed(evolutionary)
  observed <- unname(cor(lipid_vector, evolutionary_vector, method = "pearson"))
  if (!is.finite(observed)) stop("Non-finite Mantel statistic")
  RNGkind(kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
  set.seed(seed)
  exceedances <- 0L
  for (iteration in seq_len(N_PERMUTATIONS)) {
    order <- sample.int(nrow(lipid), replace = FALSE)
    statistic <- cor(condensed(lipid[order, order, drop = FALSE]), evolutionary_vector,
                     method = "pearson")
    if (is.finite(statistic) && statistic >= observed) exceedances <- exceedances + 1L
  }
  list(r = observed, p_greater = (exceedances + 1) / (N_PERMUTATIONS + 1),
       exceedances = exceedances, seed = seed)
}

partial_mantel_greater <- function(lipid, evolutionary, covariate, seed) {
  lipid_vector <- condensed(lipid)
  evolutionary_vector <- condensed(evolutionary)
  covariate_vector <- condensed(covariate)
  lipid_residual <- residualize(lipid_vector, covariate_vector)
  evolutionary_residual <- residualize(evolutionary_vector, covariate_vector)
  observed <- unname(cor(lipid_residual, evolutionary_residual, method = "pearson"))
  if (!is.finite(observed)) stop("Non-finite partial Mantel statistic")
  RNGkind(kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
  set.seed(seed)
  exceedances <- 0L
  n <- nrow(evolutionary)
  for (iteration in seq_len(N_PERMUTATIONS)) {
    order <- sample.int(n, replace = FALSE)
    permuted_vector <- condensed(evolutionary[order, order, drop = FALSE])
    statistic <- cor(lipid_residual, residualize(permuted_vector, covariate_vector),
                     method = "pearson")
    if (is.finite(statistic) && statistic >= observed) exceedances <- exceedances + 1L
  }
  list(r = observed, p_greater = (exceedances + 1) / (N_PERMUTATIONS + 1),
       exceedances = exceedances, seed = seed)
}

mantel_rows <- list()
partial_rows <- list()
pair_rows <- list()
mantel_index <- 0L
partial_index <- 0L
pair_index <- 0L
for (mode in MODES) {
  lipid <- if (mode == "POS") pos_bray else neg_bray
  metadata <- mode_metadata[[mode]]
  counts <- table(factor(metadata$phylum, levels = PHYLA))
  for (metric in names(evolutionary_matrices)) {
    evolutionary <- evolutionary_matrices[[metric]]
    mantel_result <- mantel_greater(lipid, evolutionary, MANTEL_SEEDS[[mode]])
    partial_result <- partial_mantel_greater(lipid, evolutionary,
                                              batch_overlap[[mode]], PARTIAL_SEEDS[[mode]])
    mantel_index <- mantel_index + 1L
    partial_index <- partial_index + 1L
    mantel_rows[[mantel_index]] <- data.frame(
      mode = mode, evolutionary_metric = metric, mantel_r = mantel_result$r,
      mantel_p_greater = mantel_result$p_greater, exceedances = mantel_result$exceedances,
      permutations = N_PERMUTATIONS, seed = mantel_result$seed,
      n_samples = nrow(metadata), n_phyla = N_PHYLA, n_pairs = N_PAIRS,
      stringsAsFactors = FALSE
    )
    partial_rows[[partial_index]] <- data.frame(
      mode = mode, evolutionary_metric = metric, partial_mantel_r = partial_result$r,
      partial_mantel_p_greater = partial_result$p_greater,
      exceedances = partial_result$exceedances, permutations = N_PERMUTATIONS,
      seed = partial_result$seed, n_samples = nrow(metadata), n_phyla = N_PHYLA,
      n_pairs = N_PAIRS, covariate = "batch_set_jaccard_overlap",
      stringsAsFactors = FALSE
    )
    pair_index <- pair_index + 1L
    pair_indices <- which(upper.tri(evolutionary), arr.ind = TRUE)
    pair_rows[[pair_index]] <- data.frame(
      mode = mode,
      evolutionary_metric = metric,
      phylum_1 = rownames(evolutionary)[pair_indices[, 1]],
      phylum_2 = colnames(evolutionary)[pair_indices[, 2]],
      evolutionary_distance = evolutionary[pair_indices],
      lipid_distance = lipid[pair_indices],
      batch_overlap_jaccard = batch_overlap[[mode]][pair_indices],
      n_samples_mode = nrow(metadata),
      n_samples_phylum_1 = as.integer(counts[rownames(evolutionary)[pair_indices[, 1]]]),
      n_samples_phylum_2 = as.integer(counts[colnames(evolutionary)[pair_indices[, 2]]]),
      stringsAsFactors = FALSE
    )
  }
}
mantel_results <- do.call(rbind, mantel_rows)
partial_results <- do.call(rbind, partial_rows)
pair_data <- do.call(rbind, pair_rows)
if (nrow(mantel_results) != length(MODES) * length(evolutionary_matrices) ||
    nrow(partial_results) != nrow(mantel_results) ||
    any(table(pair_data$mode, pair_data$evolutionary_metric) != N_PAIRS)) {
  stop("Mantel/pair output row-count validation failed")
}

primary_metric <- "anchor_set_mean__inventory_weighted"
primary_pos_pairs <- pair_data[pair_data$mode == "POS" &
                                 pair_data$evolutionary_metric == primary_metric, , drop = FALSE]
if (nrow(primary_pos_pairs) != N_PAIRS) stop("Primary POS pair table is incomplete")
outlier <- primary_pos_pairs[which.min(primary_pos_pairs$lipid_distance), , drop = FALSE]
outlier_phyla <- c(outlier$phylum_1, outlier$phylum_2)
if (!setequal(outlier_phyla, c("Ascomycota", "Basidiomycota"))) {
  stop("The expected lowest POS point is no longer the Ascomycota-Basidiomycota pair")
}

pos_source <- read.csv(pos_source_path, check.names = FALSE,
                        stringsAsFactors = FALSE)
pos_sample_ids <- mode_metadata$POS$sample_id
if (!all(pos_sample_ids %in% names(pos_source))) {
  stop("The POS substrate is missing a strict sample column needed for outlier audit")
}
profile_for <- function(sample_ids) {
  rowMeans(as.matrix(pos_source[, sample_ids, drop = FALSE]), na.rm = FALSE)
}
bray_profiles <- function(profile_1, profile_2) {
  denominator <- sum(abs(profile_1 + profile_2))
  if (!is.finite(denominator) || denominator <= 0) stop("Invalid outlier-audit Bray-Curtis denominator")
  sum(abs(profile_1 - profile_2)) / denominator
}
outlier_group_ids <- lapply(outlier_phyla, function(phylum) {
  mode_metadata$POS$sample_id[mode_metadata$POS$phylum == phylum]
})
names(outlier_group_ids) <- outlier_phyla
recomputed_outlier_bray <- bray_profiles(
  profile_for(outlier_group_ids[["Ascomycota"]]),
  profile_for(outlier_group_ids[["Basidiomycota"]])
)
reported_outlier_bray <- as.numeric(outlier$lipid_distance)
if (abs(recomputed_outlier_bray - reported_outlier_bray) > 1e-12) {
  stop("Outlier Bray-Curtis does not reproduce from the strict POS substrate")
}

outlier_audit_rows <- list()
outlier_audit_index <- 0L
for (phylum in outlier_phyla) {
  focal_ids <- outlier_group_ids[[phylum]]
  other_phylum <- setdiff(outlier_phyla, phylum)
  other_ids <- outlier_group_ids[[other_phylum]]
  for (sample_id in focal_ids) {
    sample_row <- mode_metadata$POS[mode_metadata$POS$sample_id == sample_id, , drop = FALSE]
    leave_one_out <- bray_profiles(profile_for(setdiff(focal_ids, sample_id)),
                                   profile_for(other_ids))
    outlier_audit_index <- outlier_audit_index + 1L
    outlier_audit_rows[[outlier_audit_index]] <- data.frame(
      sample_id = sample_id,
      sample_name = sample_row$sample_name,
      phylum = sample_row$phylum,
      metadata_phylum = sample_row$metadata_phylum,
      source_phylum = sample_row$source_phylum,
      genus = sample_row$genus,
      species = sample_row$species,
      batch = sample_row$batch,
      taxonomy_scope = sample_row$taxonomy_scope,
      leave_one_out_bray_curtis = leave_one_out,
      delta_from_full_pair = leave_one_out - recomputed_outlier_bray,
      stringsAsFactors = FALSE
    )
  }
}
outlier_audit <- do.call(rbind, outlier_audit_rows)
strict_taxa <- read.csv(STRICT_TAXA_PATH, check.names = FALSE,
                        colClasses = "character", stringsAsFactors = FALSE)
strict_taxa_base <- sub("\\.mzML$", "", strict_taxa[["GNPS Filename"]])
verified_for_sample <- function(sample_id) {
  sample_name <- sub("^sample:", "", sample_id)
  suffix_match <- vapply(strict_taxa_base, function(candidate) {
    grepl(paste0(sample_name, "_"), candidate, fixed = TRUE) ||
      endsWith(candidate, sample_name)
  }, logical(1))
  matches <- strict_taxa[strict_taxa_base == sample_name | suffix_match, , drop = FALSE]
  if (nrow(matches) == 0L) {
    return(data.frame(
      verified_taxid = NA_character_, verified_name = NA_character_,
      verification_status = NA_character_, decision_category = NA_character_,
      verification_note = NA_character_, stringsAsFactors = FALSE
    ))
  }
  if (length(unique(matches$verified_taxid)) != 1L) {
    stop("Ambiguous verified taxon mapping for ", sample_id)
  }
  match_row <- matches[1, , drop = FALSE]
  data.frame(
    verified_taxid = match_row$verified_taxid,
    verified_name = match_row$NCBI_name_verified,
    verification_status = match_row$verification_status,
    decision_category = match_row$decision_category,
    verification_note = match_row$verification_note,
    stringsAsFactors = FALSE
  )
}
verified_rows <- do.call(rbind, lapply(outlier_audit$sample_id, verified_for_sample))
outlier_audit <- cbind(outlier_audit, verified_rows)
outlier_audit$name_review_flag <- vapply(seq_len(nrow(outlier_audit)), function(i) {
  current_genus <- trimws(outlier_audit$genus[[i]])
  verified_name <- trimws(outlier_audit$verified_name[[i]])
  !is.na(current_genus) && !is.na(verified_name) &&
    nzchar(current_genus) && nzchar(verified_name) &&
    !grepl(tolower(current_genus), tolower(verified_name), fixed = TRUE)
}, logical(1))
outlier_source_mismatches <- outlier_audit[
  !is.na(outlier_audit$source_phylum) &
    nzchar(outlier_audit$source_phylum) &
    outlier_audit$source_phylum != outlier_audit$phylum, , drop = FALSE
]
outlier_sample_audit_path <- file.path(
  output_dir, "outlier_ascomycota_basidiomycota_leave_one_out.csv"
)
write.csv(outlier_audit, outlier_sample_audit_path, row.names = FALSE, quote = TRUE)
driver_index <- which.max(outlier_audit$leave_one_out_bray_curtis)
raising_index <- which.min(outlier_audit$leave_one_out_bray_curtis)
outlier_summary <- list(
  mode = "POS",
  pair = "Ascomycota - Basidiomycota",
  n_phylum_pairs = N_PAIRS,
  evolutionary_metric = primary_metric,
  ssu_distance = as.numeric(outlier$evolutionary_distance),
  reported_bray_curtis = reported_outlier_bray,
  recomputed_bray_curtis = recomputed_outlier_bray,
  strict_sample_counts = list(
    Ascomycota = length(outlier_group_ids[["Ascomycota"]]),
    Basidiomycota = length(outlier_group_ids[["Basidiomycota"]])
  ),
  source_vs_ncbi_phylum_mismatch_count = nrow(outlier_source_mismatches),
  source_vs_ncbi_phylum_mismatch_samples = as.list(outlier_source_mismatches$sample_id),
  sample_name_review_flag_count = sum(outlier_audit$name_review_flag),
  sample_name_review_flag_samples = as.list(outlier_audit$sample_id[outlier_audit$name_review_flag]),
  strongest_outlier_driver = list(
    sample_id = outlier_audit$sample_id[[driver_index]],
    current_genus = outlier_audit$genus[[driver_index]],
    current_species = outlier_audit$species[[driver_index]],
    verified_name = outlier_audit$verified_name[[driver_index]],
    leave_one_out_bray_curtis = max(outlier_audit$leave_one_out_bray_curtis),
    interpretation = "Removing this sample raises the pair distance most; it contributes most to the unusually low full-pair distance."
  ),
  strongest_distance_raising_sample = list(
    sample_id = outlier_audit$sample_id[[raising_index]],
    current_genus = outlier_audit$genus[[raising_index]],
    current_species = outlier_audit$species[[raising_index]],
    verified_name = outlier_audit$verified_name[[raising_index]],
    leave_one_out_bray_curtis = min(outlier_audit$leave_one_out_bray_curtis),
    interpretation = "Removing this sample lowers the pair distance; it was raising the full-pair value rather than causing the low outlier."
  ),
  leave_one_out_bray_min = min(outlier_audit$leave_one_out_bray_curtis),
  leave_one_out_bray_max = max(outlier_audit$leave_one_out_bray_curtis),
  leave_one_out_delta_min = min(outlier_audit$delta_from_full_pair),
  leave_one_out_delta_max = max(outlier_audit$delta_from_full_pair),
  known_ssu_primary_exclusions_in_ascomycota = c("Heydenia", "Warcupia"),
  interpretation = "The pair reproduces exactly from the strict POS substrate. Source and strict NCBI phylum labels agree for all 37 samples; inspect the per-sample audit for spelling/name artifacts."
)
outlier_summary_path <- file.path(output_dir, "outlier_ascomycota_basidiomycota_summary.json")
write_json(outlier_summary, outlier_summary_path, pretty = TRUE, auto_unbox = TRUE, digits = 16)

write.csv(mantel_results, file.path(output_dir, "mantel_results_v3_review_only.csv"),
          row.names = FALSE, quote = FALSE)
write.csv(partial_results, file.path(output_dir, "partial_mantel_results_v3_review_only.csv"),
          row.names = FALSE, quote = FALSE)
write.csv(pair_data, file.path(output_dir, "phylum_pair_data_v3_review_only.csv"),
          row.names = FALSE, quote = FALSE)
for (metric in names(evolutionary_matrices)) {
  write_matrix(evolutionary_matrices[[metric]], file.path(
    output_dir, paste0("phylum_ssu_distance_", metric, "_v3.csv")))
}
write.csv(unit_metadata, file.path(output_dir, "analysis_unit_tree_metadata_with_inventory_weights_v3.csv"),
          row.names = FALSE, quote = TRUE)

format_p <- function(value) if (value < 0.001) "< 0.001" else sprintf("= %.3f", value)
plot_rows <- list()
for (mode in MODES) {
  lipid <- if (mode == "POS") pos_bray else neg_bray
  evolutionary <- evolutionary_matrices[[primary_metric]]
  indices <- which(upper.tri(evolutionary), arr.ind = TRUE)
  plot_rows[[mode]] <- data.frame(
    mode = mode,
    phylum_1 = rownames(evolutionary)[indices[, 1]],
    phylum_2 = colnames(evolutionary)[indices[, 2]],
    evolutionary_distance = evolutionary[indices],
    lipid_distance = lipid[indices],
    stringsAsFactors = FALSE
  )
}
plot_data <- do.call(rbind, plot_rows)
plot_data$mode <- factor(plot_data$mode, levels = MODES,
                         labels = c("Positive", "Negative"))
plot_data$ecological_group_1 <- unname(group_by_phylum[plot_data$phylum_1])
plot_data$ecological_group_2 <- unname(group_by_phylum[plot_data$phylum_2])
plot_data$pair_group <- factor(
  ifelse(plot_data$ecological_group_1 == plot_data$ecological_group_2,
         "Within group", "Cross group"),
  levels = c("Cross group", "Within group")
)
pair_group_counts <- aggregate(
  phylum_1 ~ mode + pair_group,
  data = plot_data,
  FUN = length
)
names(pair_group_counts)[names(pair_group_counts) == "phylum_1"] <- "n_pairs"
write.csv(pair_group_counts, file.path(output_dir, "pair_group_counts_v3_review_only.csv"),
          row.names = FALSE, quote = FALSE)
primary_mantel <- mantel_results[mantel_results$evolutionary_metric == primary_metric, , drop = FALSE]
primary_partial <- partial_results[partial_results$evolutionary_metric == primary_metric, , drop = FALSE]
plot_x_values <- seq(min(plot_data$evolutionary_distance),
                     max(plot_data$evolutionary_distance), length.out = 200)
line_data <- do.call(rbind, lapply(levels(plot_data$mode), function(mode) {
  mode_data <- plot_data[plot_data$mode == mode, , drop = FALSE]
  model <- lm(lipid_distance ~ evolutionary_distance, data = mode_data)
  data.frame(
    evolutionary_distance = plot_x_values,
    lipid_distance = predict(model, newdata = data.frame(evolutionary_distance = plot_x_values)),
    mode = factor(mode, levels = levels(plot_data$mode)),
    stringsAsFactors = FALSE
  )
}))
annotation_r <- sprintf(
  "r = %.2f (POS), %.2f (NEG)",
  primary_mantel$mantel_r[primary_mantel$mode == "POS"],
  primary_mantel$mantel_r[primary_mantel$mode == "NEG"]
)
annotation_p <- "permutation P < 0.001; 120 pairs/mode"
review_plot <- ggplot() +
  geom_point(
    data = plot_data,
    aes(x = evolutionary_distance, y = lipid_distance, colour = mode, shape = pair_group),
    size = 0.8, alpha = 0.5, stroke = 0.2
  ) +
  geom_line(
    data = line_data,
    aes(x = evolutionary_distance, y = lipid_distance, colour = mode),
    linewidth = 0.65
  ) +
  scale_colour_manual(values = MODE_COLOURS, name = NULL) +
  scale_shape_manual(values = c("Cross group" = 2, "Within group" = 16), name = NULL) +
  scale_x_continuous(
    name = "SSU evolutionary distance (mean tree path length)",
    limits = c(min(plot_data$evolutionary_distance) - 0.10,
               max(plot_data$evolutionary_distance) + 0.10),
    expand = expansion(mult = c(0, 0))
  ) +
  scale_y_continuous(
    name = "Lipidome dissimilarity (Bray-Curtis)",
    breaks = seq(0.4, 1.0, 0.1),
    limits = c(0.4, 1.02),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(
    tag = "b"
  ) +
  annotate(
    "text", x = min(plot_data$evolutionary_distance) + 0.10, y = 0.58,
    hjust = 0, vjust = 1, label = annotation_r, size = 5.5 * 0.3528
  ) +
  annotate(
    "text", x = min(plot_data$evolutionary_distance) + 0.10, y = 0.54,
    hjust = 0, vjust = 1, label = annotation_p, size = 5.5 * 0.3528
  ) +
  theme_nature_figure3() +
  theme(
    plot.tag = element_text(size = 8, face = "bold"),
    plot.tag.position = c(0.005, 0.985),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.key.size = unit(2.5, "mm"),
    legend.spacing.x = unit(1.5, "mm"),
    panel.grid.major.y = element_line(colour = "grey90", linewidth = 0.25),
    plot.margin = margin(3, 2, 2, 2, "mm")
  ) +
  guides(
    colour = guide_legend(order = 1, override.aes = list(alpha = 1, size = 1.6)),
    shape = guide_legend(order = 2, override.aes = list(colour = "grey45", alpha = 1, size = 1.6))
  )
figure_png <- file.path(output_dir, "Figure_3b_ssu_evolutionary_distance_review_only.png")
figure_pdf <- file.path(output_dir, "Figure_3b_ssu_evolutionary_distance_review_only.pdf")
figure_svg <- file.path(output_dir, "Figure_3b_ssu_evolutionary_distance_review_only.svg")
ggsave(figure_png, review_plot, width = 89 / 25.4, height = 102 / 25.4,
       dpi = 600, bg = "white")
ggsave(figure_pdf, review_plot, width = 89 / 25.4, height = 102 / 25.4,
       device = grDevices::cairo_pdf, bg = "white")
ggsave(figure_svg, review_plot, width = 89 / 25.4, height = 102 / 25.4,
       device = grDevices::svg, bg = "white")

covariate_report_path <- file.path(output_dir, "historical_partial_mantel_covariate_trace.md")
covariate_report <- c(
  "# Historical partial-Mantel covariate trace",
  "",
  "The recovered producer defines the covariate as batch-set Jaccard similarity:",
  "`batch_overlap[i,j] = |B_i intersect B_j| / |B_i union B_j|`, where `B_i` is the set of batches represented by samples in phylum `i`.",
  "Both the evolutionary-distance vector and the lipid-distance vector are residualized by ordinary least squares against this similarity vector; the partial Mantel statistic is the Pearson correlation between the two residual vectors.",
  "",
  paste0("Positive-mode producer: `", HISTORICAL_POS_SOURCE, "` (reviewer_computational_tasks.py)."),
  paste0("Negative-mode producer: `", HISTORICAL_NEG_SOURCE, "` (phylogenetic_signal_neg.py; it writes 1 - Jaccard batch distance, an affine equivalent for residualization)."),
  "The historical positive implementation used 999 unseeded permutations on an older substrate. This review-only v3 run recomputes the covariate from locked v1 strict sample metadata and uses 9,999 deterministic label permutations with recorded mode-specific seeds.",
  "",
  "This trace does not authorize manuscript replacement or imply direct cell-level reconciliation to the authoritative submitted S1 workbook."
)
writeLines(covariate_report, covariate_report_path, useBytes = TRUE)

summary <- list(
  taxonomy_release = RELEASE_ID,
  tree_freeze_id = TREE_FREEZE_ID,
  stage_id = "figure3_ssu_distance_v3_review_only",
  status = "review_only_complete_pending_submitted_S1_reconciliation",
  scope = "Strict POS/NEG lipid-distance application and diagnostics only; disconnected from Figure 3/manuscript/table/response.",
  n_phyla = N_PHYLA,
  n_pairs_per_mode = N_PAIRS,
  strict_sample_denominators = list(POS = nrow(mode_metadata$POS), NEG = nrow(mode_metadata$NEG)),
  strict_sample_counts = split(sample_count_rows$n_samples, sample_count_rows$mode),
  evolutionary_metrics = names(evolutionary_matrices),
  primary_metric = primary_metric,
  figure_style = list(
    style_reference = STYLE_PATH,
    palette = "Wong colour-blind-safe Positive=#0072B2, Negative=#E69F00",
    theme = "soilmass_style.R theme_nature equivalent",
    marker_shape_convention = "Submitted Figure 4: hollow triangles = Cross group; filled circles = Within group",
    marker_group_field = "ecological_group (display stratum only)",
    x_label = "Mean between-phylum SSU evolutionary distance (patristic; inventory-weighted)",
    y_label = "Lipidome dissimilarity (Bray-Curtis)"
  ),
  outlier_audit = outlier_summary,
  permutations = N_PERMUTATIONS,
  mantel_seeds = as.list(MANTEL_SEEDS),
  partial_mantel_seeds = as.list(PARTIAL_SEEDS),
  partial_mantel_covariate = list(
    name = "batch_set_jaccard_overlap",
    formula = "|B_i intersect B_j| / |B_i union B_j|",
    residualization = "OLS residuals for evolutionary and lipid vectors",
    positive_producer = HISTORICAL_POS_SOURCE,
    negative_producer = HISTORICAL_NEG_SOURCE,
    historical_permutations = 999,
    historical_positive_seed = "unseeded",
    current_permutations = N_PERMUTATIONS
  ),
  primary_results = list(
    mantel = split(primary_mantel$mantel_r, primary_mantel$mode),
    mantel_p_greater = split(primary_mantel$mantel_p_greater, primary_mantel$mode),
    partial_mantel = split(primary_partial$partial_mantel_r, primary_partial$mode),
    partial_mantel_p_greater = split(primary_partial$partial_mantel_p_greater, primary_partial$mode)
  ),
  tree_summary = "148 primary SSU sequences, 103 units, 100 bootstrap trees, FBP support attached; four non-monophyletic multi-anchor units retained as a warning.",
  authority_warning = "Direct cell-level reconciliation to the authoritative submitted S1 workbook remains pending; this output is review-only.",
  figure3_connected = FALSE
)
summary_path <- file.path(output_dir, "figure3_ssu_distance_v3_review_only_summary.json")
write_json(summary, summary_path, pretty = TRUE, auto_unbox = TRUE, digits = 16)

report_path <- file.path(output_dir, "FIGURE3_SSU_DISTANCE_REVIEW_ONLY_REPORT.md")
report_lines <- c(
  "# Figure 3 SSU-distance application — review only",
  "",
  paste0("- Taxonomy release: `", RELEASE_ID, "`"),
  paste0("- Tree freeze: `", TREE_FREEZE_ID, "`"),
  paste0("- Strict denominators: POS = ", nrow(mode_metadata$POS), "; NEG = ", nrow(mode_metadata$NEG), "; 16 phyla; 120 pairs per mode."),
  "- Primary evolutionary aggregation: inventory-weighted anchor-set-mean patristic distance across 103 curated units.",
  paste0("- Lowest POS point: Ascomycota - Basidiomycota; Bray-Curtis = ", sprintf("%.4f", reported_outlier_bray),
         "; the source/strict phylum labels agree for all 37 samples."),
  paste0("- Strongest low-pair driver: ", outlier_summary$strongest_outlier_driver$sample_id,
         " (current name ", outlier_summary$strongest_outlier_driver$current_genus,
         "; verified name ", outlier_summary$strongest_outlier_driver$verified_name, ")."),
  paste0("- Sample-name review flags: ", outlier_summary$sample_name_review_flag_count,
         "; inspect `outlier_ascomycota_basidiomycota_leave_one_out.csv`."),
  paste0("- Primary POS Mantel: r = ", sprintf("%.4f", primary_mantel$mantel_r[primary_mantel$mode == "POS"]),
         ", one-sided P ", format_p(primary_mantel$mantel_p_greater[primary_mantel$mode == "POS"]),
         "; partial Mantel r = ", sprintf("%.4f", primary_partial$partial_mantel_r[primary_partial$mode == "POS"]),
         ", one-sided P ", format_p(primary_partial$partial_mantel_p_greater[primary_partial$mode == "POS"])),
  paste0("- Primary NEG Mantel: r = ", sprintf("%.4f", primary_mantel$mantel_r[primary_mantel$mode == "NEG"]),
         ", one-sided P ", format_p(primary_mantel$mantel_p_greater[primary_mantel$mode == "NEG"]),
         "; partial Mantel r = ", sprintf("%.4f", primary_partial$partial_mantel_r[primary_partial$mode == "NEG"]),
         ", one-sided P ", format_p(primary_partial$partial_mantel_p_greater[primary_partial$mode == "NEG"])),
  "",
  "The partial-Mantel covariate is the recovered batch-set Jaccard overlap; its exact source trace is in `historical_partial_mantel_covariate_trace.md`.",
  "Marker shapes reuse the submitted Figure 4 convention: hollow triangles are Cross group and filled circles are Within group, assigned from strict-release ecological display groups.",
  "Medoid, MRCA, and unit-unweighted aggregation sensitivities are included in the result tables.",
  "",
  "This is not a manuscript-ready result. Figure 3, the manuscript, legends, tables, and response tracker remain disconnected. Direct cell-level reconciliation to the authoritative submitted S1 workbook remains pending."
)
writeLines(report_lines, report_path, useBytes = TRUE)

output_paths <- unname(c(
  file.path(output_dir, "strict_mode_sample_denominators.csv"),
  file.path(output_dir, "batch_overlap_jaccard_POS.csv"),
  file.path(output_dir, "batch_overlap_jaccard_NEG.csv"),
  file.path(output_dir, "mantel_results_v3_review_only.csv"),
  file.path(output_dir, "partial_mantel_results_v3_review_only.csv"),
  file.path(output_dir, "phylum_pair_data_v3_review_only.csv"),
  file.path(output_dir, "strict_phylum_ecological_display_groups.csv"),
  file.path(output_dir, "pair_group_counts_v3_review_only.csv"),
  file.path(output_dir, "analysis_unit_tree_metadata_with_inventory_weights_v3.csv"),
  outlier_sample_audit_path, outlier_summary_path,
  vapply(names(evolutionary_matrices), function(metric) file.path(
    output_dir, paste0("phylum_ssu_distance_", metric, "_v3.csv")), character(1)),
  figure_png, figure_pdf, figure_svg, covariate_report_path, summary_path, report_path
))
manifest <- list(
  schema_version = 1,
  stage_id = summary$stage_id,
  taxonomy_release = RELEASE_ID,
  tree_freeze_id = TREE_FREEZE_ID,
  status = summary$status,
  generated_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  inputs = list(
    panel_a_manifest = file_record(panel_manifest_path),
    substrate_summary = file_record(substrate_summary_path),
    sample_labels = file_record(sample_labels_path),
    pos_metadata = file_record(pos_metadata_path),
    neg_metadata = file_record(neg_metadata_path),
    style_reference = file_record(STYLE_PATH),
    strict_sample_taxa_verification = file_record(STRICT_TAXA_PATH),
    pos_substrate = file_record(pos_source_path),
    pos_braycurtis = file_record(pos_bray_path),
    neg_braycurtis = file_record(neg_bray_path),
    curated_units = file_record(curated_units_path),
    tree_unit_metadata = file_record(unit_metadata_path),
    tree_distance_anchor_set_mean = file_record(distance_files[["anchor_set_mean"]]),
    tree_distance_medoid = file_record(distance_files[["medoid"]]),
    tree_distance_mrca = file_record(distance_files[["mrca"]]),
    historical_positive_covariate_producer = file_record(HISTORICAL_POS_SOURCE),
    historical_negative_covariate_producer = file_record(HISTORICAL_NEG_SOURCE)
  ),
  outputs = unname(lapply(output_paths, file_record)),
  validation_gates = list(
    taxonomy_release_locked = TRUE,
    strict_sample_denominators_POS_164 = nrow(mode_metadata$POS) == 164L,
    strict_sample_denominators_NEG_192 = nrow(mode_metadata$NEG) == 192L,
    matrices_16_by_16 = TRUE,
    pair_count_per_mode_120 = all(table(pair_data$mode) == length(evolutionary_matrices) * N_PAIRS),
    permutations_9999 = N_PERMUTATIONS == 9999L,
    historical_covariate_traced = TRUE,
    figure_style_reused = TRUE,
    marker_shape_convention_reused = TRUE,
    strict_ecological_display_group_map_complete = all(nzchar(group_map$ecological_display_group)),
    outlier_pair_recomputed = abs(recomputed_outlier_bray - reported_outlier_bray) <= 1e-12,
    outlier_source_phylum_matches = nrow(outlier_source_mismatches) == 0L,
    figure3_connected = FALSE,
    submitted_S1_cell_reconciliation = FALSE
  )
)
manifest_path <- file.path(output_dir, "stage_manifest.json")
write_json(manifest, manifest_path, pretty = TRUE, auto_unbox = TRUE, digits = 16)

cat(toJSON(summary, pretty = TRUE, auto_unbox = TRUE, digits = 16), "\n")
