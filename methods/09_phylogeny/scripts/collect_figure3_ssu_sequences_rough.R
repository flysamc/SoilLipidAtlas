#!/usr/bin/env Rscript
# Collect the sequences in the v2 Figure 3 SSU candidate freeze for a rough
# end-to-end tree smoke test. SILVA records are already clipped to SSU rRNA.
# NCBI fallbacks are fetched as whole public accession records for this smoke
# test; annotated-feature extraction is a later publication-grade gate.

suppressPackageStartupMessages({
  library(Biostrings)
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
freeze_path <- normalizePath(args[["freeze"]], mustWork = TRUE)
silva_fasta_path <- normalizePath(args[["silva-fasta"]], mustWork = TRUE)
output_dir <- args[["output-dir"]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
output_dir <- normalizePath(output_dir, mustWork = TRUE)
cache_dir <- file.path(output_dir, "ncbi_efetch_cache")
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

sha256 <- function(path) digest(path, algo = "sha256", file = TRUE, serialize = FALSE)
file_record <- function(path) {
  path <- normalizePath(path, mustWork = TRUE)
  list(path = path, bytes = unname(file.info(path)$size), sha256 = sha256(path))
}

freeze <- read.csv(freeze_path, check.names = FALSE, colClasses = "character",
                   stringsAsFactors = FALSE)
if (nrow(freeze) != 145 || anyDuplicated(freeze$accession_version)) {
  stop("Expected 145 unique sequence records in v2 candidate freeze")
}
if (length(unique(freeze$analysis_unit_taxid)) != 105) stop("Freeze must cover 105 analysis units")
freeze$tip_id <- sprintf("T%03d_%s_A%s", seq_len(nrow(freeze)),
                         freeze$analysis_unit_taxid, freeze$anchor_index)

silva_rows <- which(freeze$source_database == "SILVA")
ncbi_rows <- which(freeze$source_database == "NCBI Nucleotide")
if (length(silva_rows) != 115 || length(ncbi_rows) != 30) {
  stop("Unexpected SILVA/NCBI split in v2 freeze")
}

message("Reading SILVA 138.2 NR99 truncated SSU FASTA...")
silva_all <- readRNAStringSet(silva_fasta_path, format = "fasta", use.names = TRUE)
silva_ids <- sub("[[:space:]].*$", "", names(silva_all))
silva_match <- match(freeze$silva_tree_tip[silva_rows], silva_ids)
if (any(is.na(silva_match))) {
  stop("Frozen SILVA tips missing from FASTA: ",
       paste(freeze$silva_tree_tip[silva_rows][is.na(silva_match)], collapse = ", "))
}
silva_selected_rna <- silva_all[silva_match]
silva_selected <- DNAStringSet(chartr("U", "T", as.character(silva_selected_rna)))
names(silva_selected) <- freeze$tip_id[silva_rows]
rm(silva_all, silva_selected_rna)
invisible(gc())

message("Fetching NCBI fallback accessions in two cached batches...")
ncbi_accessions <- freeze$accession_version[ncbi_rows]
batch_index <- split(seq_along(ncbi_accessions), ceiling(seq_along(ncbi_accessions) / 15))
batch_files <- character(length(batch_index))
for (index in seq_along(batch_index)) {
  accessions <- ncbi_accessions[batch_index[[index]]]
  batch_file <- file.path(cache_dir, sprintf("ncbi_ssu_batch_%02d.fasta", index))
  batch_files[[index]] <- batch_file
  if (!file.exists(batch_file) || file.info(batch_file)$size == 0) {
    query <- paste(accessions, collapse = ",")
    url_string <- paste0(
      "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?",
      "db=nuccore&rettype=fasta&retmode=text&tool=SoilMassFigure3&id=",
      URLencode(query, reserved = TRUE)
    )
    last_error <- ""
    complete <- FALSE
    for (attempt in seq_len(4)) {
      complete <- tryCatch({
        download.file(url_string, batch_file, mode = "wb", method = "libcurl", quiet = TRUE)
        file.exists(batch_file) && file.info(batch_file)$size > 0
      }, error = function(error) {
        last_error <<- conditionMessage(error)
        FALSE
      })
      if (complete) break
      Sys.sleep(2^attempt)
    }
    if (!complete) stop("NCBI EFetch failed for batch ", index, ": ", last_error)
    Sys.sleep(0.4)
  }
}

ncbi_sets <- lapply(batch_files, readDNAStringSet, format = "fasta", use.names = TRUE)
ncbi_all <- do.call(c, ncbi_sets)
ncbi_ids <- sub("[[:space:]].*$", "", names(ncbi_all))
ncbi_ids <- sub("^.*\\|", "", ncbi_ids)
ncbi_match <- match(ncbi_accessions, ncbi_ids)
if (any(is.na(ncbi_match))) {
  # NCBI headers can omit the requested version while retaining the accession.
  ncbi_match <- match(sub("\\.[0-9]+$", "", ncbi_accessions),
                      sub("\\.[0-9]+$", "", ncbi_ids))
}
if (any(is.na(ncbi_match))) {
  stop("Frozen NCBI accessions missing from cached FASTA: ",
       paste(ncbi_accessions[is.na(ncbi_match)], collapse = ", "))
}
ncbi_selected <- ncbi_all[ncbi_match]
names(ncbi_selected) <- freeze$tip_id[ncbi_rows]

all_sequences <- DNAStringSet(rep("N", nrow(freeze)))
all_sequences[silva_rows] <- silva_selected
all_sequences[ncbi_rows] <- ncbi_selected
names(all_sequences) <- freeze$tip_id

sequence_lengths <- width(all_sequences)
alphabet <- alphabetFrequency(all_sequences, baseOnly = FALSE)
ambiguous_count <- rowSums(alphabet[, !colnames(alphabet) %in% c("A", "C", "G", "T") , drop = FALSE])
tip_map <- freeze
tip_map$rough_sequence_length <- as.character(sequence_lengths)
tip_map$rough_ambiguous_bases <- as.character(ambiguous_count)
tip_map$rough_ambiguous_fraction <- sprintf("%.8f", ambiguous_count / sequence_lengths)
tip_map$rough_sequence_scope <- ifelse(
  tip_map$source_database == "SILVA",
  "SILVA clipped SSU rRNA",
  "NCBI whole accession record; annotated SSU extraction pending"
)

fasta_path <- file.path(output_dir, "figure3_ssu_sequences_rough.fasta")
tip_map_path <- file.path(output_dir, "sequence_tip_map.csv")
summary_path <- file.path(output_dir, "sequence_collection_summary.json")
manifest_path <- file.path(output_dir, "stage_manifest.json")
writeXStringSet(all_sequences, fasta_path, format = "fasta", width = 80)
write.csv(tip_map, tip_map_path, row.names = FALSE, quote = TRUE, na = "")

summary <- list(
  taxonomy_release = unique(freeze$taxonomy_release),
  freeze_id = unique(freeze$freeze_id),
  stage_id = "figure3_ssu_sequence_collection_rough",
  status = "rough_sequence_collection_complete",
  sequences = length(all_sequences),
  analysis_units = length(unique(freeze$analysis_unit_taxid)),
  source_counts = as.list(table(freeze$source_database)),
  length = list(min = min(sequence_lengths), median = unname(median(sequence_lengths)),
                max = max(sequence_lengths)),
  sequences_below_900_nt = sum(sequence_lengths < 900),
  sequences_above_3500_nt = sum(sequence_lengths > 3500),
  important_limitation = paste(
    "NCBI fallback records are whole public accessions in this rough pass;",
    "publication-grade analysis must extract and validate only annotated SSU features"
  )
)
write_json(summary, summary_path, pretty = TRUE, auto_unbox = TRUE, digits = 16)

outputs <- c(fasta_path, tip_map_path, summary_path, batch_files)
manifest <- list(
  schema_version = 1,
  stage_id = summary$stage_id,
  taxonomy_release = summary$taxonomy_release,
  freeze_id = summary$freeze_id,
  status = summary$status,
  generated_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  inputs = list(freeze = file_record(freeze_path), silva_fasta = file_record(silva_fasta_path)),
  outputs = lapply(outputs, file_record),
  validation_gates = list(
    expected_sequences_145 = length(all_sequences) == 145,
    expected_analysis_units_105 = length(unique(freeze$analysis_unit_taxid)) == 105,
    every_frozen_accession_matched = TRUE,
    ncbi_annotated_ssu_features_extracted = FALSE
  )
)
write_json(manifest, manifest_path, pretty = TRUE, auto_unbox = TRUE, digits = 16)

cat(toJSON(summary, pretty = TRUE, auto_unbox = TRUE, digits = 16), "\n")
