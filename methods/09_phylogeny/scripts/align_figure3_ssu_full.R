#!/usr/bin/env Rscript
# Structural SSU alignment and deterministic informative-column mask for the
# curated Figure 3 evolutionary-tree candidate.

suppressPackageStartupMessages({
  library(Biostrings)
  library(DECIPHER)
  library(jsonlite)
  library(digest)
})

parse_args <- function(values) {
  out <- list()
  for (value in values) {
    pieces <- strsplit(sub("^--", "", value), "=", fixed = TRUE)[[1]]
    out[[pieces[[1]]]] <- paste(pieces[-1], collapse = "=")
  }
  out
}
args <- parse_args(commandArgs(trailingOnly = TRUE))
required <- c("sequences", "tip-map", "output-dir")
missing <- setdiff(required, names(args))
if (length(missing)) stop("Missing arguments: ", paste(missing, collapse = ", "))

sequence_path <- normalizePath(args$sequences, mustWork = TRUE)
tip_map_path <- normalizePath(args$`tip-map`, mustWork = TRUE)
output_dir <- args$`output-dir`
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
output_dir <- normalizePath(output_dir, mustWork = TRUE)

RELEASE_ID <- "ncbi-phylum-2026-08-04-v1"
FREEZE_ID <- "figure3-ssu-curated-freeze-2026-08-04-v3"
set.seed(20260804)

sha256 <- function(path) digest(path, algo = "sha256", file = TRUE, serialize = FALSE)
record <- function(path) {
  path <- normalizePath(path, mustWork = TRUE)
  list(path = path, bytes = unname(file.info(path)$size), sha256 = sha256(path))
}

sequences <- readDNAStringSet(sequence_path, use.names = TRUE)
tip_map <- read.csv(tip_map_path, check.names = FALSE, colClasses = "character")
if (length(sequences) != nrow(tip_map) || !identical(names(sequences), tip_map$tip_id)) {
  stop("Curated sequence/tip-map mismatch")
}
if (length(unique(tip_map$analysis_unit_taxid)) != 103 ||
    length(unique(tip_map$phylum)) != 16) stop("Unexpected curated coverage")

message("Running DECIPHER structural SSU alignment...")
aligned <- AlignSeqs(sequences, iterations = 2, refinements = 1,
                     useStructures = TRUE, processors = 4, verbose = TRUE)
alignment_matrix <- toupper(as.matrix(aligned))
base_matrix <- alignment_matrix %in% c("A", "C", "G", "T")
dim(base_matrix) <- dim(alignment_matrix)
occupancy <- colMeans(base_matrix)
variable <- apply(alignment_matrix, 2, function(x) {
  length(unique(x[x %in% c("A", "C", "G", "T")])) >= 2
})
retained <- occupancy >= 0.50 & variable
if (sum(retained) < 500) stop("Fewer than 500 informative alignment columns")
masked_matrix <- alignment_matrix[, retained, drop = FALSE]
masked <- DNAStringSet(apply(masked_matrix, 1, paste0, collapse = ""))
names(masked) <- rownames(masked_matrix)

alignment_path <- file.path(output_dir, "ssu_structural_alignment_v3.fasta")
masked_path <- file.path(output_dir, "ssu_structural_alignment_masked_v3.fasta")
mask_path <- file.path(output_dir, "ssu_alignment_mask_v3.csv")
summary_path <- file.path(output_dir, "ssu_alignment_v3_summary.json")
manifest_path <- file.path(output_dir, "stage_manifest.json")
writeXStringSet(aligned, alignment_path, width = 80)
writeXStringSet(masked, masked_path, width = 80)
write.csv(data.frame(column = seq_along(occupancy), occupancy = occupancy,
                     variable = variable, retained = retained),
          mask_path, row.names = FALSE, quote = FALSE)

summary <- list(
  taxonomy_release = RELEASE_ID, freeze_id = FREEZE_ID,
  stage_id = "figure3_ssu_structural_alignment_v3",
  status = "structural_alignment_complete_candidate",
  sequences = length(sequences), analysis_units = length(unique(tip_map$analysis_unit_taxid)),
  phyla = length(unique(tip_map$phylum)), structural_alignment = TRUE,
  iterations = 2, refinements = 1, unmasked_columns = ncol(alignment_matrix),
  retained_variable_columns = sum(retained), occupancy_threshold = 0.50,
  retained_fraction = sum(retained) / length(retained)
)
write_json(summary, summary_path, pretty = TRUE, auto_unbox = TRUE, digits = 16)
outputs <- c(alignment_path, masked_path, mask_path, summary_path)
manifest <- list(
  schema_version = 1, stage_id = summary$stage_id, taxonomy_release = RELEASE_ID,
  freeze_id = FREEZE_ID, status = summary$status,
  generated_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE), random_seed = 20260804,
  software = list(R = R.version.string, Biostrings = as.character(packageVersion("Biostrings")),
                  DECIPHER = as.character(packageVersion("DECIPHER"))),
  inputs = list(sequences = record(sequence_path), tip_map = record(tip_map_path)),
  outputs = lapply(outputs, record),
  validation_gates = list(sequence_names_match_tip_map = identical(names(masked), tip_map$tip_id),
                          all_103_primary_units_present = length(unique(tip_map$analysis_unit_taxid)) == 103,
                          all_16_phyla_present = length(unique(tip_map$phylum)) == 16,
                          structural_alignment_used = TRUE,
                          informative_columns_at_least_500 = sum(retained) >= 500)
)
write_json(manifest, manifest_path, pretty = TRUE, auto_unbox = TRUE, digits = 16)
cat(toJSON(summary, pretty = TRUE, auto_unbox = TRUE, digits = 16), "\n")
