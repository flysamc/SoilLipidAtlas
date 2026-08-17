#!/usr/bin/env Rscript
# Curate the publication-candidate Figure 3 SSU sequence set from the preserved
# v2 candidate freeze. The primary set accepts only SSU sequences >=900 nt,
# replaces unusable exact records with declared same-genus SILVA proxies when
# possible, and excludes units with no defensible SSU route. Exact partials are
# retained only in a separate sensitivity set.

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
required <- c("freeze", "tip-map", "full-fasta", "alternative-audit",
              "alternative-fasta", "readiness", "silva-taxmap", "silva-fasta",
              "output-dir")
missing <- setdiff(required, names(args))
if (length(missing)) stop("Missing arguments: ", paste(missing, collapse = ", "))

paths <- lapply(args[required[required != "output-dir"]], normalizePath, mustWork = TRUE)
output_dir <- args[["output-dir"]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
output_dir <- normalizePath(output_dir, mustWork = TRUE)

RELEASE_ID <- "ncbi-phylum-2026-08-04-v1"
FREEZE_ID <- "figure3-ssu-curated-freeze-2026-08-04-v3"
MIN_PRIMARY_LENGTH <- 900L

sha256 <- function(path) digest(path, algo = "sha256", file = TRUE, serialize = FALSE)
file_record <- function(path) {
  path <- normalizePath(path, mustWork = TRUE)
  list(path = path, bytes = unname(file.info(path)$size), sha256 = sha256(path))
}
terminal_name <- function(path) {
  values <- strsplit(path, ";", fixed = TRUE)[[1]]
  values <- values[nzchar(values)]
  if (length(values)) tail(values, 1) else ""
}
valid_species <- function(value, genus) {
  grepl(paste0("^", genus, " [[:alpha:]][[:alnum:]_-]+$"), value) &&
    !grepl("(?i) (sp|cf|aff)\\.?$|uncultured|unidentified", value, perl = TRUE)
}

freeze <- read.csv(paths$freeze, check.names = FALSE, colClasses = "character")
tip_map <- read.csv(paths$`tip-map`, check.names = FALSE, colClasses = "character")
readiness <- read.csv(paths$readiness, check.names = FALSE, colClasses = "character")
alternatives <- read.csv(paths$`alternative-audit`, check.names = FALSE,
                         colClasses = "character")
full_sequences <- readDNAStringSet(paths$`full-fasta`)
alternative_sequences <- readDNAStringSet(paths$`alternative-fasta`)
if (nrow(freeze) != 145 || length(unique(freeze$analysis_unit_taxid)) != 105) {
  stop("Unexpected v2 freeze dimensions")
}
if (!identical(names(full_sequences), tip_map$tip_id)) stop("Full FASTA/tip map mismatch")

# Units with corrected explicit annotated SSU features.
replacement_units <- c("2482753", "4043")
replacement_rows <- alternatives[
  alternatives$analysis_unit_taxid %in% replacement_units &
    alternatives$recommended_replacement == "TRUE", , drop = FALSE]
if (nrow(replacement_rows) != 2) stop("Expected two approved exact replacements")

# Exact records below the primary length gate are represented by same-genus
# SILVA proxies. The requested anchor count is capped by available named taxa.
proxy_plan <- data.frame(
  analysis_unit_taxid = c("1331584", "186352", "44115", "67695", "50225", "47427"),
  proxy_genus = c("Auricularia", "Cephalotrichum", "Peziza", "Amanita", "Taraxacum", "Armillaria"),
  requested_anchors = c(3L, 2L, 3L, 3L, 1L, 1L),
  stringsAsFactors = FALSE
)
excluded_units <- c("931642", "352928") # Heydenia and Warcupia: no valid SSU/proxy route.
remove_units <- c(replacement_units, proxy_plan$analysis_unit_taxid, excluded_units)

kept <- freeze[!freeze$analysis_unit_taxid %in% remove_units, , drop = FALSE]
kept_old_tip <- tip_map$tip_id[match(paste(kept$analysis_unit_taxid, kept$accession_version),
                                     paste(tip_map$analysis_unit_taxid, tip_map$accession_version))]
if (anyNA(kept_old_tip)) stop("Could not map retained v2 sequences")
kept_sequences <- full_sequences[match(kept_old_tip, names(full_sequences))]

# Replace legacy MRCA labels: primary downstream distances aggregate all anchor
# pair distances, while MRCA/medoid variants are sensitivity analyses.
kept$representation[kept$representation == "genus_mrca"] <- "genus_anchor_set_mean"
kept$representation[kept$representation == "approved_genus_proxy_mrca"] <-
  "approved_genus_proxy_anchor_set_mean"

new_rows <- list()
new_sequences <- list()
new_index <- 0L
for (i in seq_len(nrow(replacement_rows))) {
  alt <- replacement_rows[i, , drop = FALSE]
  old <- freeze[freeze$analysis_unit_taxid == alt$analysis_unit_taxid, , drop = FALSE][1, , drop = FALSE]
  old$freeze_id <- FREEZE_ID
  old$representation <- "exact_taxon_tip"
  old$sequence_role <- "analysis_unit_tip"
  old$anchor_index <- "1"
  old$source_database <- "NCBI Nucleotide"
  old$source_release <- "GenBank record cached 2026-08-04"
  old$accession_version <- alt$accession_version
  old$source_organism <- alt$analysis_unit_name
  old$source_record_length <- alt$extracted_length
  old$silva_start <- old$silva_stop <- old$silva_tree_tip <- ""
  old$candidate_pool_size <- ""
  old$selection_score <- alt$improvement_nt
  old$selection_rule <- "Corrected parser: explicit annotated 18S rRNA feature; joined ranges concatenated; primary length >=900 nt"
  old$proxy_flag <- "FALSE"
  old$proxy_reason <- ""
  old$sequence_fetch_method <- "Annotated GenBank SSU feature extracted from cached record"
  old$approval_status <- "curated_primary_candidate_pending_submitted_S1_reconciliation"
  seq_index <- match(alt$accession_version, names(alternative_sequences))
  if (is.na(seq_index)) stop("Alternative sequence missing: ", alt$accession_version)
  new_index <- new_index + 1L
  new_rows[[new_index]] <- old
  new_sequences[[new_index]] <- alternative_sequences[seq_index]
}

taxmap <- read.delim(paths$`silva-taxmap`, header = TRUE, check.names = FALSE,
                     colClasses = "character", quote = "", comment.char = "")
names(taxmap) <- c("accession", "start", "stop", "ncbi_path", "submitted_name")
taxmap$length <- abs(as.numeric(taxmap$stop) - as.numeric(taxmap$start)) + 1
taxmap$terminal_name <- vapply(taxmap$ncbi_path, terminal_name, character(1))
taxmap$tree_tip <- paste(taxmap$accession, taxmap$start, taxmap$stop, sep = ".")

proxy_selected <- list()
for (i in seq_len(nrow(proxy_plan))) {
  plan <- proxy_plan[i, , drop = FALSE]
  marker <- paste0(";", plan$proxy_genus, ";")
  pool <- taxmap[grepl(marker, taxmap$ncbi_path, fixed = TRUE) & taxmap$length >= MIN_PRIMARY_LENGTH, , drop = FALSE]
  pool <- pool[vapply(pool$terminal_name, valid_species, logical(1), genus = plan$proxy_genus), , drop = FALSE]
  if (!nrow(pool)) stop("No SILVA proxy pool for ", plan$proxy_genus)
  pool$score <- 1000 - abs(pool$length - 1800)
  pool <- pool[order(-pool$score, pool$terminal_name, pool$tree_tip), , drop = FALSE]
  if (plan$proxy_genus != "Cephalotrichum") {
    pool <- pool[!duplicated(pool$terminal_name), , drop = FALSE]
  }
  selected <- head(pool, plan$requested_anchors)
  if (!nrow(selected)) stop("No selected proxy for ", plan$proxy_genus)
  proxy_selected[[plan$analysis_unit_taxid]] <- selected
}

message("Reading SILVA FASTA for curated proxy sequences...")
silva_all <- readRNAStringSet(paths$`silva-fasta`, use.names = TRUE)
silva_ids <- sub("[[:space:]].*$", "", names(silva_all))
for (taxid in names(proxy_selected)) {
  selected <- proxy_selected[[taxid]]
  unit <- readiness[readiness$verified_taxid == taxid, , drop = FALSE]
  if (nrow(unit) != 1) stop("Readiness unit missing: ", taxid)
  for (j in seq_len(nrow(selected))) {
    src <- selected[j, , drop = FALSE]
    fasta_index <- match(src$tree_tip, silva_ids)
    if (is.na(fasta_index)) stop("SILVA proxy missing from FASTA: ", src$tree_tip)
    template <- freeze[freeze$analysis_unit_taxid == taxid, , drop = FALSE][1, , drop = FALSE]
    template$freeze_id <- FREEZE_ID
    template$representation <- if (nrow(selected) > 1) "genus_proxy_anchor_set_mean" else "single_genus_proxy_tip"
    template$sequence_role <- "proxy_genus_anchor"
    template$anchor_index <- as.character(j)
    template$source_database <- "SILVA"
    template$source_release <- "138.2"
    template$accession_version <- src$accession
    template$source_record_taxid <- ""
    template$source_organism <- src$terminal_name
    template$source_record_length <- as.character(src$length)
    template$silva_start <- src$start
    template$silva_stop <- src$stop
    template$silva_tree_tip <- src$tree_tip
    template$candidate_pool_size <- as.character(nrow(selected))
    template$selection_score <- as.character(src$score)
    template$selection_rule <- "Same-genus named SILVA proxy; SSU >=900 nt; closest to 1800 nt; distinct species where available"
    template$proxy_flag <- "TRUE"
    template$proxy_reason <- paste0("Exact SSU failed the >=900 nt primary gate; explicit ",
                                    proxy_plan$proxy_genus[proxy_plan$analysis_unit_taxid == taxid],
                                    " genus proxy; exact partial retained only for sensitivity")
    template$sequence_fetch_method <- "SILVA 138.2 NR99 FASTA truncated SSU rRNA"
    template$approval_status <- "curated_primary_candidate_pending_submitted_S1_reconciliation"
    new_index <- new_index + 1L
    new_rows[[new_index]] <- template
    new_sequences[[new_index]] <- DNAStringSet(chartr("U", "T", as.character(silva_all[fasta_index])))
  }
}
rm(silva_all)
invisible(gc())

added <- do.call(rbind, new_rows)
curated <- rbind(kept, added)
curated$freeze_id <- FREEZE_ID
curated$approval_status <- "curated_primary_candidate_pending_submitted_S1_reconciliation"
order_index <- order(curated$phylum, curated$analysis_unit_name,
                     as.numeric(curated$anchor_index), curated$accession_version)
curated <- curated[order_index, , drop = FALSE]

# Reassemble sequences in the same row order using an accession-based map.
all_added_sequences <- do.call(c, new_sequences)
names(all_added_sequences) <- added$accession_version
names(kept_sequences) <- kept$accession_version
sequence_pool <- c(kept_sequences, all_added_sequences)
primary_sequences <- sequence_pool[match(curated$accession_version, names(sequence_pool))]
if (anyNA(match(curated$accession_version, names(sequence_pool)))) stop("Curated sequence map incomplete")
curated$tip_id <- sprintf("V3T%03d_%s_A%s", seq_len(nrow(curated)),
                          curated$analysis_unit_taxid, curated$anchor_index)
names(primary_sequences) <- curated$tip_id
curated$sequence_length <- as.character(width(primary_sequences))
curated$primary_length_gate <- as.character(width(primary_sequences) >= MIN_PRIMARY_LENGTH)
if (any(width(primary_sequences) < MIN_PRIMARY_LENGTH)) stop("Primary sequence below 900 nt")
if (anyDuplicated(curated$tip_id) || anyDuplicated(curated$accession_version)) stop("Duplicate curated tips/accessions")

# Preserve the rejected exact/partial sequences as a separate diagnostic set.
sensitivity_old <- tip_map[tip_map$analysis_unit_taxid %in% c(proxy_plan$analysis_unit_taxid, excluded_units), , drop = FALSE]
sensitivity_sequences <- full_sequences[match(sensitivity_old$tip_id, names(full_sequences))]
names(sensitivity_sequences) <- paste0("SENS_", sensitivity_old$tip_id)
sensitivity_old$sensitivity_tip_id <- names(sensitivity_sequences)
sensitivity_old$sensitivity_reason <- ifelse(
  sensitivity_old$analysis_unit_taxid %in% excluded_units,
  "No primary SSU/proxy route; retained only to quantify exclusion sensitivity",
  "Exact SSU below 900 nt or invalid broad record; genus proxy used in primary set"
)

unit_rows <- lapply(split(seq_len(nrow(curated)), curated$analysis_unit_taxid), function(idx) {
  x <- curated[idx, , drop = FALSE]
  data.frame(analysis_unit_taxid = x$analysis_unit_taxid[1],
             analysis_unit_name = x$analysis_unit_name[1], phylum = x$phylum[1],
             representation = x$representation[1], anchors = nrow(x),
             proxy = any(x$proxy_flag == "TRUE"), included_primary = TRUE,
             inventory_rows = x$analysis_unit_inventory_rows[1], stringsAsFactors = FALSE)
})
unit_table <- do.call(rbind, unit_rows)
excluded_table <- readiness[readiness$verified_taxid %in% excluded_units,
                            c("verified_taxid", "verified_name", "phylum", "n_inventory_rows"), drop = FALSE]
names(excluded_table) <- c("analysis_unit_taxid", "analysis_unit_name", "phylum", "inventory_rows")
excluded_table$representation <- "excluded_no_defensible_ssu_route"
excluded_table$anchors <- 0L
excluded_table$proxy <- FALSE
excluded_table$included_primary <- FALSE
unit_table <- rbind(unit_table, excluded_table[, names(unit_table)])
unit_table <- unit_table[order(unit_table$phylum, unit_table$analysis_unit_name), ]

freeze_path <- file.path(output_dir, "ssu_curated_freeze_v3.csv")
fasta_path <- file.path(output_dir, "ssu_curated_primary_v3.fasta")
unit_path <- file.path(output_dir, "ssu_curated_unit_representations_v3.csv")
sensitivity_fasta_path <- file.path(output_dir, "ssu_rejected_exact_partials_sensitivity_v3.fasta")
sensitivity_map_path <- file.path(output_dir, "ssu_rejected_exact_partials_sensitivity_v3.csv")
summary_path <- file.path(output_dir, "ssu_curated_freeze_v3_summary.json")
report_path <- file.path(output_dir, "SSU_CURATED_FREEZE_V3_REPORT.md")
manifest_path <- file.path(output_dir, "stage_manifest.json")
write.csv(curated, freeze_path, row.names = FALSE, quote = TRUE, na = "")
writeXStringSet(primary_sequences, fasta_path, width = 80)
write.csv(unit_table, unit_path, row.names = FALSE, quote = TRUE, na = "")
writeXStringSet(sensitivity_sequences, sensitivity_fasta_path, width = 80)
write.csv(sensitivity_old, sensitivity_map_path, row.names = FALSE, quote = TRUE, na = "")

summary <- list(
  taxonomy_release = RELEASE_ID, freeze_id = FREEZE_ID,
  status = "curated_primary_candidate_frozen_pending_submitted_S1_reconciliation",
  primary_sequences = length(primary_sequences),
  primary_analysis_units = length(unique(curated$analysis_unit_taxid)),
  excluded_analysis_units = nrow(excluded_table),
  excluded_units = as.list(setNames(excluded_table$analysis_unit_name, excluded_table$analysis_unit_taxid)),
  phyla = length(unique(curated$phylum)), minimum_primary_length = min(width(primary_sequences)),
  exact_replacements = as.list(setNames(replacement_rows$accession_version,
                                         replacement_rows$analysis_unit_name)),
  proxy_units = sum(vapply(split(curated$proxy_flag, curated$analysis_unit_taxid),
                           function(x) any(x == "TRUE"), logical(1))),
  primary_distance_rule = "Mean pairwise patristic distance across all anchor sequences per analysis-unit pair",
  sensitivity_rules = c("MRCA node", "medoid/single anchor", "rejected exact partials where alignable"),
  authority_warning = paste("Identity source is the verified annotated CSV; direct cell-level reconciliation",
                            "to the authoritative submitted S1 workbook is still required before release.")
)
write_json(summary, summary_path, pretty = TRUE, auto_unbox = TRUE, digits = 16)
report <- c(
  "# Figure 3 SSU curated freeze v3", "",
  paste0("- Taxonomy release: `", RELEASE_ID, "`"),
  paste0("- Freeze: `", FREEZE_ID, "`"),
  paste0("- Primary: ", summary$primary_sequences, " sequences / ", summary$primary_analysis_units,
         " analysis units / ", summary$phyla, " phyla"),
  paste0("- Primary length gate: >=", MIN_PRIMARY_LENGTH, " nt; observed minimum ",
         summary$minimum_primary_length, " nt"),
  paste0("- Excluded: ", paste(excluded_table$analysis_unit_name, collapse = ", ")),
  "- Multi-anchor primary representation: mean pairwise patristic distance, not MRCA.",
  "- Exact partials and invalid broad records are retained only in a separate sensitivity artifact.",
  "", "## Authority warning", "", summary$authority_warning
)
writeLines(report, report_path, useBytes = TRUE)

outputs <- c(freeze_path, fasta_path, unit_path, sensitivity_fasta_path,
             sensitivity_map_path, summary_path, report_path)
manifest <- list(
  schema_version = 1, stage_id = "figure3_ssu_curated_freeze_v3",
  taxonomy_release = RELEASE_ID, freeze_id = FREEZE_ID,
  status = summary$status, generated_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  inputs = lapply(paths, file_record), outputs = lapply(outputs, file_record),
  validation_gates = list(
    primary_sequences_all_at_least_900_nt = all(width(primary_sequences) >= MIN_PRIMARY_LENGTH),
    expected_two_exact_replacements = nrow(replacement_rows) == 2,
    expected_six_proxy_units = length(unique(proxy_plan$analysis_unit_taxid)) == 6,
    excluded_units_explicit = setequal(excluded_table$analysis_unit_name, c("Heydenia", "Warcupia")),
    all_16_phyla_retained = length(unique(curated$phylum)) == 16,
    accession_and_tip_ids_unique = !anyDuplicated(curated$accession_version) && !anyDuplicated(curated$tip_id),
    submitted_S1_cell_level_reconciliation_complete = FALSE
  )
)
write_json(manifest, manifest_path, pretty = TRUE, auto_unbox = TRUE, digits = 16)
cat(toJSON(summary, pretty = TRUE, auto_unbox = TRUE, digits = 16), "\n")
