#!/usr/bin/env Rscript

# Finalize the v3 organism-level SSU tree for review-only Figure 3 work.
# This stage combines the four deterministic bootstrap chunks, attaches FBP
# support to the checkpointed ML tree, and derives unit-level patristic
# sensitivities. It does not run lipid distances, Mantel tests, or rendering.

suppressPackageStartupMessages({
  library(phangorn)
  library(ape)
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
required <- c("alignment", "tip-map", "fit", "fit-summary", "bootstrap-dir", "output-dir")
missing <- setdiff(required, names(args))
if (length(missing) > 0) stop("Missing arguments: ", paste(missing, collapse = ", "))

alignment_path <- normalizePath(args[["alignment"]], mustWork = TRUE)
tip_map_path <- normalizePath(args[["tip-map"]], mustWork = TRUE)
fit_path <- normalizePath(args[["fit"]], mustWork = TRUE)
fit_summary_path <- normalizePath(args[["fit-summary"]], mustWork = TRUE)
bootstrap_dir <- normalizePath(args[["bootstrap-dir"]], mustWork = TRUE)
output_dir <- args[["output-dir"]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
output_dir <- normalizePath(output_dir, mustWork = TRUE)

RELEASE_ID <- "ncbi-phylum-2026-08-04-v1"
FREEZE_ID <- "figure3-ssu-curated-freeze-2026-08-04-v3"
N_TIPS <- 148L
N_BOOTSTRAP <- 100L

sha256 <- function(path) digest(path, algo = "sha256", file = TRUE, serialize = FALSE)
file_record <- function(path) {
  path <- normalizePath(path, mustWork = TRUE)
  list(path = path, bytes = unname(file.info(path)$size), sha256 = sha256(path))
}

write_matrix <- function(value, path) {
  out <- data.frame(label = rownames(value), value, check.names = FALSE,
                    stringsAsFactors = FALSE)
  names(out)[[1]] <- ""
  write.csv(out, path, row.names = FALSE, quote = FALSE, na = "")
}

validate_matrix <- function(value, name, expected_n) {
  if (!identical(dim(value), c(expected_n, expected_n))) stop(name, " dimensions invalid")
  if (is.null(rownames(value)) || !identical(rownames(value), colnames(value))) {
    stop(name, " row/column names differ")
  }
  if (any(!is.finite(value)) || any(value < -1e-12)) stop(name, " has invalid values")
  if (max(abs(value - t(value))) > 1e-10) stop(name, " is asymmetric")
  if (max(abs(diag(value))) > 1e-10) stop(name, " has non-zero diagonal")
  invisible(TRUE)
}

fit_summary <- fromJSON(fit_summary_path, simplifyVector = TRUE)
if (!identical(fit_summary$taxonomy_release, RELEASE_ID) ||
    !identical(fit_summary$freeze_id, FREEZE_ID) ||
    !identical(fit_summary$status, "checkpointed_gtr_gamma_nni_fit_complete")) {
  stop("The ML fit summary is not the locked v3 checkpoint")
}

tip_map <- read.csv(tip_map_path, check.names = FALSE, colClasses = "character",
                    stringsAsFactors = FALSE)
if (nrow(tip_map) != N_TIPS || anyDuplicated(tip_map$tip_id)) {
  stop("The v3 tip map must contain 148 unique primary tips")
}
if (!all(toupper(tip_map$included_primary) == "TRUE")) {
  stop("The final tip map contains a non-primary row")
}
if (length(unique(tip_map$taxonomy_release)) != 1 ||
    !identical(unique(tip_map$taxonomy_release), RELEASE_ID) ||
    length(unique(tip_map$freeze_id)) != 1 ||
    !identical(unique(tip_map$freeze_id), FREEZE_ID)) {
  stop("Tip-map release/freeze mismatch")
}

alignment <- read.dna(alignment_path, format = "fasta", as.matrix = TRUE)
if (nrow(alignment) != N_TIPS || ncol(alignment) != 1593L ||
    !identical(rownames(alignment), tip_map$tip_id)) {
  stop("Masked alignment does not match the v3 148-tip, 1,593-column contract")
}

fit <- readRDS(fit_path)
if (!inherits(fit, "pml")) stop("ML checkpoint is not a phangorn pml object")
ml_tree <- reorder.phylo(fit$tree, "cladewise")
if (Ntip(ml_tree) != N_TIPS || !setequal(ml_tree$tip.label, tip_map$tip_id) ||
    is.null(ml_tree$edge.length) ||
    any(!is.finite(ml_tree$edge.length)) || any(ml_tree$edge.length < 0)) {
  stop("ML tree validation failed")
}

chunk_ids <- 1:4
chunk_rds <- file.path(bootstrap_dir, sprintf("bootstrap_chunk_%02d.rds", chunk_ids))
chunk_json <- file.path(bootstrap_dir, sprintf("bootstrap_chunk_%02d.json", chunk_ids))
if (!all(file.exists(chunk_rds)) || !all(file.exists(chunk_json))) {
  stop("All four bootstrap chunks and their JSON summaries are required")
}

chunk_summaries <- lapply(chunk_json, function(path) fromJSON(path, simplifyVector = TRUE))
if (!all(vapply(seq_along(chunk_summaries), function(i) {
  s <- chunk_summaries[[i]]
  identical(as.integer(s$chunk), chunk_ids[[i]]) &&
    identical(as.integer(s$replicates), 25L) &&
    identical(s$status, "bootstrap_chunk_complete")
}, logical(1)))) stop("Bootstrap chunk summaries do not describe four 25-replicate chunks")

chunk_trees <- lapply(chunk_rds, readRDS)
if (!all(vapply(chunk_trees, length, integer(1)) == 25L)) {
  stop("A bootstrap chunk does not contain exactly 25 trees")
}
bootstrap_trees <- do.call(c, chunk_trees)
if (length(bootstrap_trees) != N_BOOTSTRAP) stop("Expected exactly 100 bootstrap trees")
# `do.call(c, list(multiPhylo, ...))` drops the multiPhylo class in some
# phangorn/ape combinations. Restore it before calling phangorn:::support;
# support treats a plain list as a single malformed tree collection.
class(bootstrap_trees) <- "multiPhylo"

bootstrap_trees <- lapply(bootstrap_trees, function(tree) {
  if (!inherits(tree, "phylo")) stop("Bootstrap object is not a phylo tree")
  if (is.rooted(tree)) tree <- unroot(tree)
  if (Ntip(tree) != N_TIPS || !setequal(tree$tip.label, ml_tree$tip.label)) {
    stop("Bootstrap tree tip set differs from the ML tree")
  }
  if (is.null(tree$edge.length) || any(!is.finite(tree$edge.length)) ||
      any(tree$edge.length < 0)) stop("Bootstrap tree has invalid branch lengths")
  tree
})
class(bootstrap_trees) <- "multiPhylo"

fbp_fraction <- phangorn:::support(ml_tree, bootstrap_trees, method = "FBP", scale = TRUE)
fbp_percent <- phangorn:::support(ml_tree, bootstrap_trees, method = "FBP", scale = FALSE)
if (length(fbp_fraction) != ml_tree$Nnode || any(!is.finite(fbp_fraction)) ||
    any(fbp_fraction < 0 | fbp_fraction > 1) ||
    any(!is.finite(fbp_percent)) || any(fbp_percent < 0 | fbp_percent > 100)) {
  stop("FBP support validation failed")
}

tree_with_fbp <- ml_tree
tree_with_fbp$node.label <- sprintf("%.2f", fbp_percent)
tree_with_fbp$comment <- paste0("FBP support percentage from ", N_BOOTSTRAP,
                                " deterministic bootstrap trees")

tip_distance <- cophenetic.phylo(ml_tree)
validate_matrix(tip_distance, "tip patristic distance", N_TIPS)
tip_index <- setNames(seq_along(ml_tree$tip.label), ml_tree$tip.label)

primary_units <- unique(tip_map$analysis_unit_taxid)
if (length(primary_units) != 103L) stop("Expected 103 primary analysis units")
unit_rows <- lapply(primary_units, function(unit_id) {
  rows <- tip_map[tip_map$analysis_unit_taxid == unit_id, , drop = FALSE]
  fields <- c("analysis_unit_name", "phylum", "representation", "proxy_flag")
  if (any(vapply(fields, function(field) length(unique(rows[[field]])) != 1, logical(1)))) {
    stop("Mixed unit metadata for analysis unit ", unit_id)
  }
  list(id = unit_id, rows = rows, tips = rows$tip_id)
})
unit_names <- vapply(unit_rows, function(x) unique(x$rows$analysis_unit_name), character(1))
if (anyDuplicated(unit_names)) stop("Analysis unit names are not unique")
unit_phyla <- vapply(unit_rows, function(x) unique(x$rows$phylum), character(1))
unit_representation <- vapply(unit_rows, function(x) unique(x$rows$representation), character(1))
unit_proxy <- vapply(unit_rows, function(x) unique(x$rows$proxy_flag), character(1))
unit_tip_sets <- lapply(unit_rows, function(x) x$tips)

within_mean <- vapply(unit_tip_sets, function(tips) {
  if (length(tips) < 2) return(0)
  mean(tip_distance[tips, tips][upper.tri(tip_distance[tips, tips])])
}, numeric(1))

medoid_tip <- vapply(unit_tip_sets, function(tips) {
  if (length(tips) == 1) return(tips[[1]])
  candidates <- tips[order(vapply(tips, function(tip) mean(tip_distance[tip, tips]), numeric(1)), tips)]
  candidates[[1]]
}, character(1))

node_for_unit <- vapply(seq_along(unit_tip_sets), function(i) {
  tip_ids <- unname(tip_index[unit_tip_sets[[i]]])
  if (length(tip_ids) == 1) return(as.integer(tip_ids[[1]]))
  node <- getMRCA(ml_tree, tip_ids)
  if (is.null(node) || is.na(node)) stop("Could not resolve MRCA for ", unit_names[[i]])
  as.integer(node)
}, integer(1))

descendant_count <- vapply(seq_along(node_for_unit), function(i) {
  node <- node_for_unit[[i]]
  if (node <= Ntip(ml_tree)) return(1L)
  length(phangorn::Descendants(ml_tree, node, type = "tips")[[1]])
}, integer(1))
anchors_monophyletic <- vapply(seq_along(node_for_unit), function(i) {
  if (length(unit_tip_sets[[i]]) == 1) return(NA)
  node <- node_for_unit[[i]]
  descendants <- phangorn::Descendants(ml_tree, node, type = "tips")[[1]]
  setequal(descendants, unname(tip_index[unit_tip_sets[[i]]]))
}, logical(1))

node_distance <- ape::dist.nodes(ml_tree)
if (any(!is.finite(node_distance))) stop("Non-finite node distances")
mrca_distance <- node_distance[as.character(node_for_unit), as.character(node_for_unit), drop = FALSE]
rownames(mrca_distance) <- unit_names
colnames(mrca_distance) <- unit_names

anchor_mean_distance <- matrix(0, nrow = length(primary_units), ncol = length(primary_units),
                               dimnames = list(unit_names, unit_names))
medoid_distance <- anchor_mean_distance
for (i in seq_along(unit_tip_sets)) {
  for (j in seq_along(unit_tip_sets)) {
    if (i == j) next
    anchor_mean_distance[i, j] <- mean(tip_distance[unit_tip_sets[[i]], unit_tip_sets[[j]], drop = FALSE])
    medoid_distance[i, j] <- tip_distance[medoid_tip[[i]], medoid_tip[[j]]]
  }
}
validate_matrix(anchor_mean_distance, "anchor-set-mean unit distance", length(primary_units))
validate_matrix(medoid_distance, "medoid unit distance", length(primary_units))
validate_matrix(mrca_distance, "MRCA unit distance", length(primary_units))

unit_metadata <- data.frame(
  unit_index = seq_along(primary_units),
  analysis_unit_taxid = primary_units,
  analysis_unit_name = unit_names,
  phylum = unit_phyla,
  representation = unit_representation,
  proxy = unit_proxy,
  n_anchor_sequences = lengths(unit_tip_sets),
  selected_tip_ids = vapply(unit_tip_sets, paste, collapse = " | ", FUN.VALUE = character(1)),
  medoid_tip_id = medoid_tip,
  within_anchor_mean_patristic = within_mean,
  mrca_node_id = node_for_unit,
  mrca_node_type = ifelse(node_for_unit <= Ntip(ml_tree), "tip", "internal_mrca"),
  mrca_descendant_tip_count = descendant_count,
  anchors_monophyletic = ifelse(is.na(anchors_monophyletic), "not_applicable",
                                as.character(anchors_monophyletic)),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

pair_rows <- list()
pair_index <- 1L
for (i in seq_len(length(primary_units) - 1L)) {
  for (j in seq.int(i + 1L, length(primary_units))) {
    pair_rows[[pair_index]] <- data.frame(
      unit_1 = primary_units[[i]],
      unit_2 = primary_units[[j]],
      name_1 = unit_names[[i]],
      name_2 = unit_names[[j]],
      phylum_1 = unit_phyla[[i]],
      phylum_2 = unit_phyla[[j]],
      anchor_set_mean = anchor_mean_distance[i, j],
      medoid = medoid_distance[i, j],
      mrca = mrca_distance[i, j],
      stringsAsFactors = FALSE
    )
    pair_index <- pair_index + 1L
  }
}
pair_sensitivities <- do.call(rbind, pair_rows)
if (nrow(pair_sensitivities) != choose(length(primary_units), 2)) {
  stop("Unit pair sensitivity table has the wrong number of rows")
}

fbp_table <- data.frame(
  node_id = Ntip(ml_tree) + seq_len(ml_tree$Nnode),
  fbp_support_fraction = unname(fbp_fraction),
  fbp_support_percent = unname(fbp_percent),
  stringsAsFactors = FALSE
)
fbp_table$descendant_tip_count <- vapply(fbp_table$node_id, function(node) {
  length(phangorn::Descendants(ml_tree, node, type = "tips")[[1]])
}, integer(1))
fbp_table$descendant_tip_hash <- vapply(fbp_table$node_id, function(node) {
  paste(sort(ml_tree$tip.label[phangorn::Descendants(ml_tree, node, type = "tips")[[1]]]),
        collapse = "|")
}, character(1))

combined_bootstrap_path <- file.path(output_dir, "bootstrap_trees_100_v3.rds")
tree_path <- file.path(output_dir, "organism_ssu_ml_gtr_gamma_nni_v3_fbp_supported.nwk")
tree_rds_path <- file.path(output_dir, "organism_ssu_ml_gtr_gamma_nni_v3_fbp_supported.rds")
fbp_path <- file.path(output_dir, "organism_ssu_ml_gtr_gamma_nni_v3_fbp_support.csv")
unit_meta_path <- file.path(output_dir, "analysis_unit_tree_nodes_v3.csv")
anchor_path <- file.path(output_dir, "analysis_unit_patristic_distance_anchor_set_mean_v3.csv")
medoid_path <- file.path(output_dir, "analysis_unit_patristic_distance_medoid_v3.csv")
mrca_path <- file.path(output_dir, "analysis_unit_patristic_distance_mrca_v3.csv")
pair_path <- file.path(output_dir, "analysis_unit_patristic_pair_sensitivities_v3.csv")
within_path <- file.path(output_dir, "analysis_unit_anchor_within_dispersion_v3.csv")
chunk_inventory_path <- file.path(output_dir, "bootstrap_chunk_inventory_v3.csv")

saveRDS(bootstrap_trees, combined_bootstrap_path)
write.tree(tree_with_fbp, tree_path)
saveRDS(tree_with_fbp, tree_rds_path)
write.csv(fbp_table, fbp_path, row.names = FALSE, quote = TRUE)
write.csv(unit_metadata, unit_meta_path, row.names = FALSE, quote = TRUE)
write_matrix(anchor_mean_distance, anchor_path)
write_matrix(medoid_distance, medoid_path)
write_matrix(mrca_distance, mrca_path)
write.csv(pair_sensitivities, pair_path, row.names = FALSE, quote = TRUE)
write.csv(data.frame(
  analysis_unit_taxid = primary_units,
  analysis_unit_name = unit_names,
  phylum = unit_phyla,
  n_anchor_sequences = lengths(unit_tip_sets),
  within_anchor_mean_patristic = within_mean,
  stringsAsFactors = FALSE
), within_path, row.names = FALSE, quote = TRUE)
chunk_inventory <- do.call(rbind, lapply(seq_along(chunk_summaries), function(i) {
  data.frame(
    chunk = as.integer(chunk_summaries[[i]]$chunk),
    replicates = as.integer(chunk_summaries[[i]]$replicates),
    seed = as.integer(chunk_summaries[[i]]$seed),
    rds_path = normalizePath(chunk_rds[[i]], mustWork = TRUE),
    rds_bytes = unname(file.info(chunk_rds[[i]])$size),
    rds_sha256 = sha256(chunk_rds[[i]]),
    json_path = normalizePath(chunk_json[[i]], mustWork = TRUE),
    json_bytes = unname(file.info(chunk_json[[i]])$size),
    json_sha256 = sha256(chunk_json[[i]]),
    stringsAsFactors = FALSE
  )
}))
write.csv(chunk_inventory, chunk_inventory_path, row.names = FALSE, quote = TRUE)

summary <- list(
  taxonomy_release = RELEASE_ID,
  freeze_id = FREEZE_ID,
  stage_id = "figure3_ssu_tree_finalize_v3",
  status = "complete_candidate_pending_submitted_S1_reconciliation",
  scope = paste(
    "Tree inference finalization only; no lipid matrix, Mantel/partial-Mantel,",
    "Figure 3, manuscript, legend, table, or response output was connected."
  ),
  fit = list(
    model = "GTR+Gamma (4 categories)",
    topology_search = "NNI",
    log_likelihood = unname(fit$logLik),
    sequences = Ntip(ml_tree),
    alignment_columns = ncol(alignment),
    branch_length_min = min(ml_tree$edge.length),
    branch_length_median = unname(median(ml_tree$edge.length)),
    branch_length_max = max(ml_tree$edge.length)
  ),
  bootstrap = list(
    chunks = 4,
    replicates_per_chunk = 25,
    total_replicates = length(bootstrap_trees),
    seeds = vapply(chunk_summaries, function(x) as.integer(x$seed), integer(1)),
    support_method = "FBP",
    support_scale = "fraction and percent; Newick labels are percent"
  ),
  unit_mapping = list(
    primary_analysis_units = length(primary_units),
    primary_sequences = nrow(tip_map),
    phyla = length(unique(unit_phyla)),
    proxy_units = sum(toupper(unit_proxy) == "TRUE"),
    multi_anchor_units = sum(lengths(unit_tip_sets) > 1),
    nonmonophyletic_multi_anchor_units = sum(anchors_monophyletic == FALSE, na.rm = TRUE),
    explicit_exclusions = c("931642 Heydenia", "352928 Warcupia")
  ),
  distances = list(
    primary = "mean pairwise patristic distance across all anchor-sequence pairs between units",
    sensitivity_1 = "patristic distance between within-unit medoid anchor tips",
    sensitivity_2 = "patristic distance between within-unit MRCA nodes",
    unit_pair_count = nrow(pair_sensitivities)
  ),
  authority_warning = paste(
    "Identity source is the verified annotated CSV; direct cell-level reconciliation",
    "to the authoritative submitted S1 workbook remains pending."
  ),
  software = list(
    R = R.version.string,
    ape = as.character(packageVersion("ape")),
    phangorn = as.character(packageVersion("phangorn")),
    jsonlite = as.character(packageVersion("jsonlite")),
    digest = as.character(packageVersion("digest"))
  )
)
summary_path <- file.path(output_dir, "figure3_ssu_tree_finalize_v3_summary.json")
write_json(summary, summary_path, pretty = TRUE, auto_unbox = TRUE, digits = 16)

input_paths <- c(alignment_path, tip_map_path, fit_path, fit_summary_path,
                 chunk_rds, chunk_json)
output_paths <- c(combined_bootstrap_path, tree_path, tree_rds_path, fbp_path,
                  unit_meta_path, anchor_path, medoid_path, mrca_path,
                  pair_path, within_path, chunk_inventory_path, summary_path)
manifest <- list(
  schema_version = 1,
  stage_id = summary$stage_id,
  taxonomy_release = RELEASE_ID,
  freeze_id = FREEZE_ID,
  status = summary$status,
  generated_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  inputs = lapply(input_paths, file_record),
  outputs = lapply(output_paths, file_record),
  validation_gates = list(
    alignment_148x1593 = nrow(alignment) == 148 && ncol(alignment) == 1593,
    ml_tips_148 = Ntip(ml_tree) == 148,
    bootstrap_trees_100 = length(bootstrap_trees) == 100,
    fbp_support_146 = length(fbp_fraction) == ml_tree$Nnode,
    unit_count_103 = length(primary_units) == 103,
    anchor_distance_symmetric = max(abs(anchor_mean_distance - t(anchor_mean_distance))) <= 1e-10,
    medoid_distance_symmetric = max(abs(medoid_distance - t(medoid_distance))) <= 1e-10,
    mrca_distance_symmetric = max(abs(mrca_distance - t(mrca_distance))) <= 1e-10,
    submitted_s1_reconciliation = FALSE,
    figure3_connected = FALSE
  )
)
manifest_path <- file.path(output_dir, "stage_manifest.json")
write_json(manifest, manifest_path, pretty = TRUE, auto_unbox = TRUE, digits = 16)

cat(toJSON(list(
  status = summary$status,
  sequences = Ntip(ml_tree),
  bootstrap_replicates = length(bootstrap_trees),
  fbp_nodes = length(fbp_fraction),
  analysis_units = length(primary_units),
  nonmonophyletic_multi_anchor_units = summary$unit_mapping$nonmonophyletic_multi_anchor_units,
  output_dir = output_dir
), pretty = TRUE, auto_unbox = TRUE, digits = 16), "\n")
