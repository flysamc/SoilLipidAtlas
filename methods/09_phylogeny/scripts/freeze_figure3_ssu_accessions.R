#!/usr/bin/env Rscript
# Create a deterministic, versioned candidate SSU accession freeze for the
# organism-level evolutionary tree used downstream by Figure 3. SILVA 138.2
# supplies clipped SSU sequences whenever the verified taxon has coverage;
# NCBI Nucleotide candidates are used only for exact taxa absent from SILVA.
# Genus-labelled analysis units and the two user-approved genus proxies are
# represented by explicit multi-species anchor sets whose MRCA will represent
# that analysis unit after tree inference.

suppressPackageStartupMessages({
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
required_args <- c("readiness", "ncbi-ranks", "silva-coverage", "silva-taxmap",
                   "ncbi-candidates", "output-dir")
missing_args <- setdiff(required_args, names(args))
if (length(missing_args) > 0) stop("Missing arguments: ", paste(missing_args, collapse = ", "))

readiness_path <- normalizePath(args[["readiness"]], mustWork = TRUE)
ncbi_ranks_path <- normalizePath(args[["ncbi-ranks"]], mustWork = TRUE)
coverage_path <- normalizePath(args[["silva-coverage"]], mustWork = TRUE)
taxmap_path <- normalizePath(args[["silva-taxmap"]], mustWork = TRUE)
ncbi_candidates_path <- normalizePath(args[["ncbi-candidates"]], mustWork = TRUE)
output_dir <- args[["output-dir"]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
output_dir <- normalizePath(output_dir, mustWork = TRUE)

RELEASE_ID <- "ncbi-phylum-2026-08-04-v1"
FREEZE_ID <- "figure3-ssu-accession-candidate-freeze-2026-08-04-v2"
SILVA_RELEASE <- "138.2"

sha256 <- function(path) digest(path, algo = "sha256", file = TRUE, serialize = FALSE)
file_record <- function(path) {
  path <- normalizePath(path, mustWork = TRUE)
  list(path = path, bytes = unname(file.info(path)$size), sha256 = sha256(path))
}
normalise_taxid <- function(value) sub("\\.0$", "", trimws(as.character(value)))
terminal_name <- function(path) {
  pieces <- strsplit(path, ";", fixed = TRUE)[[1]]
  pieces <- pieces[nzchar(pieces)]
  if (length(pieces) == 0) "" else tail(pieces, 1)
}
infer_domain <- function(path) {
  if (grepl(";Bacteria <prokaryotes>;", path, fixed = TRUE)) "Bacteria" else
    if (grepl(";Archaea;", path, fixed = TRUE)) "Archaea" else "Eukaryota"
}
desired_length <- function(domain) ifelse(domain %in% c("Bacteria", "Archaea"), 1500, 1800)
valid_species_name <- function(name, genus) {
  escaped_genus <- gsub("([][{}()+*^$|\\?.])", "\\\\\\1", genus)
  ok <- grepl(paste0("^", escaped_genus, "[[:space:]]+[[:alpha:]][[:alnum:]_-]+"), name,
              ignore.case = TRUE, perl = TRUE)
  excluded <- grepl("(?i)[[:space:]](sp|cf|aff)\\.?([[:space:]]|$)|uncultured|unidentified|bacterium|archaeon|symbiont",
                    name, perl = TRUE)
  ok && !excluded
}

readiness <- read.csv(readiness_path, check.names = FALSE, colClasses = "character",
                      stringsAsFactors = FALSE)
ranks <- read.csv(ncbi_ranks_path, check.names = FALSE, colClasses = "character",
                  stringsAsFactors = FALSE)
coverage <- read.csv(coverage_path, check.names = FALSE, colClasses = "character",
                     stringsAsFactors = FALSE)
ncbi <- read.csv(ncbi_candidates_path, check.names = FALSE, colClasses = "character",
                 stringsAsFactors = FALSE)
ncbi$source_species <- ""
if (nrow(readiness) != 105 || anyDuplicated(readiness$verified_taxid)) {
  stop("Readiness input must contain 105 unique verified TaxIDs")
}
if (nrow(ranks) != 105 || anyDuplicated(ranks$verified_taxid)) {
  stop("NCBI rank input must contain 105 unique verified TaxIDs")
}
if (nrow(coverage) != 105 || anyDuplicated(coverage$verified_taxid)) {
  stop("SILVA coverage input must contain 105 unique verified TaxIDs")
}
if (length(unique(ncbi$verified_taxid)) != 103) {
  stop("NCBI candidate input must cover 103 verified TaxIDs")
}

taxmap <- read.delim(taxmap_path, header = TRUE, check.names = FALSE,
                     colClasses = "character", quote = "", comment.char = "",
                     stringsAsFactors = FALSE)
names(taxmap) <- c("primary_accession", "start", "stop", "ncbi_path", "submitted_name")
if (nrow(taxmap) != 510495) stop("Unexpected SILVA taxmap row count")
taxmap$start_num <- as.numeric(taxmap$start)
taxmap$stop_num <- as.numeric(taxmap$stop)
taxmap$region_length <- abs(taxmap$stop_num - taxmap$start_num) + 1
taxmap$tree_tip <- paste(taxmap$primary_accession, taxmap$start, taxmap$stop, sep = ".")
taxmap$terminal_name <- vapply(taxmap$ncbi_path, terminal_name, character(1))
taxmap$domain <- vapply(taxmap$ncbi_path, infer_domain, character(1))

score_silva <- function(rows) {
  target <- desired_length(rows$domain)
  in_acceptable_band <- ifelse(rows$domain %in% c("Bacteria", "Archaea"),
                               rows$region_length >= 1100 & rows$region_length <= 2200,
                               rows$region_length >= 1200 & rows$region_length <= 3500)
  in_acceptable_band * 1000 - abs(rows$region_length - target) +
    ifelse(grepl("^[A-Z]{1,4}[0-9]+$", rows$primary_accession), 1, 0)
}

silva_selection_row <- function(unit, source_row, representation, role, anchor_index,
                                pool_size, selection_rule, proxy_flag, proxy_reason) {
  data.frame(
    freeze_id = FREEZE_ID,
    taxonomy_release = RELEASE_ID,
    analysis_unit_taxid = unit$verified_taxid[[1]],
    analysis_unit_name = unit$verified_name[[1]],
    phylum = unit$phylum[[1]],
    analysis_unit_inventory_rows = unit$n_inventory_rows[[1]],
    representation = representation,
    sequence_role = role,
    anchor_index = as.character(anchor_index),
    source_database = "SILVA",
    source_release = SILVA_RELEASE,
    accession_version = source_row$primary_accession[[1]],
    source_record_taxid = "",
    source_organism = source_row$terminal_name[[1]],
    source_record_length = as.character(source_row$region_length[[1]]),
    silva_start = source_row$start[[1]],
    silva_stop = source_row$stop[[1]],
    silva_tree_tip = source_row$tree_tip[[1]],
    candidate_pool_size = as.character(pool_size),
    selection_score = as.character(source_row$selection_score[[1]]),
    selection_rule = selection_rule,
    proxy_flag = as.character(proxy_flag),
    proxy_reason = proxy_reason,
    sequence_fetch_method = "SILVA 138.2 NR99 FASTA already truncated to SSU rRNA",
    approval_status = "user_approved_candidate_freeze_2026-08-04",
    stringsAsFactors = FALSE
  )
}

ncbi_score <- function(rows) {
  length_num <- suppressWarnings(as.numeric(rows$sequence_length))
  target <- desired_length(rows$domain)
  title <- rows$title
  is_refseq <- grepl("^NR_", rows$accession_version)
  marker <- grepl("(?i)16S ribosomal|18S ribosomal|small subunit ribosomal", title, perl = TRUE)
  organelle <- grepl("(?i)mitochond|chloroplast|plastid", title, perl = TRUE)
  other_region <- grepl("(?i)internal transcribed spacer|large subunit|25S|26S|28S", title, perl = TRUE)
  complete <- grepl("(?i)complete sequence|complete gene", title, perl = TRUE)
  partial <- grepl("(?i)partial", title, perl = TRUE)
  in_band <- ifelse(rows$domain %in% c("Bacteria", "Archaea"),
                    length_num >= 1100 & length_num <= 2200,
                    length_num >= 1200 & length_num <= 3500)
  exact_record_taxid <- normalise_taxid(rows$record_taxid) == normalise_taxid(rows$verified_taxid)
  # Whole records containing ITS/LSU are never aligned as-is; the annotated SSU
  # feature is extracted later. Therefore record length is informative only for
  # standalone SSU records. A multi-region record with no "partial" qualifier is
  # preferred over a short partial amplicon.
  closeness <- ifelse(other_region, 0, pmax(-500, 200 - abs(length_num - target) / 5))
  in_band_score <- ifelse(other_region, 0, in_band * 300)
  complete_multiregion <- other_region & complete & !partial
  is_refseq * 2000 + marker * 500 + (!other_region) * 250 + complete * 120 -
    partial * 25 + in_band_score + complete_multiregion * 700 +
    exact_record_taxid * 25 + closeness -
    organelle * 10000 - as.numeric(rows$search_rank) / 1000
}

ncbi_selection_row <- function(unit, source_row, pool_size,
                               representation = "exact_taxon_tip",
                               role = "analysis_unit_tip", anchor_index = 1L,
                               selection_rule = NULL) {
  if (is.null(selection_rule)) {
    selection_rule <- paste(
      "Highest deterministic exact-taxon NCBI candidate score: RefSeq NR preferred;",
      "SSU marker required; clean standalone SSU and near-full length preferred;",
      "organelle records rejected; annotated SSU feature will be extracted"
    )
  }
  data.frame(
    freeze_id = FREEZE_ID,
    taxonomy_release = RELEASE_ID,
    analysis_unit_taxid = unit$verified_taxid[[1]],
    analysis_unit_name = unit$verified_name[[1]],
    phylum = unit$phylum[[1]],
    analysis_unit_inventory_rows = unit$n_inventory_rows[[1]],
    representation = representation,
    sequence_role = role,
    anchor_index = as.character(anchor_index),
    source_database = "NCBI Nucleotide",
    source_release = "live accession metadata retrieved 2026-08-04",
    accession_version = source_row$accession_version[[1]],
    source_record_taxid = source_row$record_taxid[[1]],
    source_organism = source_row$title[[1]],
    source_record_length = source_row$sequence_length[[1]],
    silva_start = "",
    silva_stop = "",
    silva_tree_tip = "",
    candidate_pool_size = as.character(pool_size),
    selection_score = as.character(source_row$selection_score[[1]]),
    selection_rule = selection_rule,
    proxy_flag = "FALSE",
    proxy_reason = "",
    sequence_fetch_method = "NCBI GenBank record; extract annotated 16S/18S rRNA feature only",
    approval_status = "user_approved_candidate_freeze_2026-08-04",
    stringsAsFactors = FALSE
  )
}

extract_title_species <- function(title, genus) {
  pattern <- paste0("^", genus, "[[:space:]]+([[:alpha:]][[:alnum:]_-]+)")
  hit <- regexec(pattern, title, ignore.case = TRUE, perl = TRUE)
  pieces <- regmatches(title, hit)[[1]]
  if (length(pieces) < 2) return(NA_character_)
  epithet <- pieces[[2]]
  if (tolower(gsub("\\.$", "", epithet)) %in% c("sp", "cf", "aff")) return(NA_character_)
  paste(genus, epithet)
}

selected <- list()
unit_map <- list()
ncbi_rejected <- list()
selection_index <- 0L
mapping_index <- 0L
rejection_index <- 0L

for (index in seq_len(nrow(readiness))) {
  unit <- readiness[index, , drop = FALSE]
  unit_rank <- ranks$ncbi_rank[match(unit$verified_taxid[[1]], ranks$verified_taxid)]
  if (length(unit_rank) != 1 || is.na(unit_rank)) stop("NCBI rank missing for analysis unit")
  coverage_row <- coverage[coverage$verified_taxid == unit$verified_taxid[[1]], , drop = FALSE]
  if (nrow(coverage_row) != 1) stop("Coverage row missing for ", unit$verified_taxid[[1]])
  status <- coverage_row$coverage_status[[1]]
  unit_selected <- NULL

  if (status %in% c("species_level_supported", "below_species_taxid_supported",
                    "verified_taxid_clade_supported", "species_name_path_supported_provisional")) {
    pool <- taxmap[startsWith(taxmap$ncbi_path, coverage_row$silva_path[[1]]), , drop = FALSE]
    pool$selection_score <- score_silva(pool)
    pool <- pool[pool$region_length >= 900, , drop = FALSE]
    if (nrow(pool) == 0) stop("No usable SILVA exact-taxon sequence for ", unit$verified_name[[1]])
    pool <- pool[order(-pool$selection_score, pool$tree_tip), , drop = FALSE]
    unit_selected <- silva_selection_row(
      unit, pool[1, , drop = FALSE], "exact_taxon_tip", "analysis_unit_tip", 1,
      nrow(pool),
      "Exact verified taxon SILVA clade; select near-full SSU sequence closest to domain target length",
      FALSE, ""
    )
  } else if (status %in% c("genus_level_only", "genus_name_path_only_provisional")) {
    genus <- unit$verified_name[[1]]
    pool <- taxmap[startsWith(taxmap$ncbi_path, coverage_row$silva_path[[1]]), , drop = FALSE]
    pool <- pool[vapply(pool$terminal_name, valid_species_name, logical(1), genus = genus), , drop = FALSE]
    pool$selection_score <- score_silva(pool)
    pool <- pool[pool$region_length >= 900, , drop = FALSE]
    if (nrow(pool) == 0) stop("No named descendant-species anchors for genus ", genus)
    best_by_species <- do.call(rbind, lapply(split(seq_len(nrow(pool)), pool$terminal_name), function(rows) {
      candidate <- pool[rows, , drop = FALSE]
      candidate[order(-candidate$selection_score, candidate$tree_tip), , drop = FALSE][1, , drop = FALSE]
    }))
    best_by_species <- best_by_species[order(-best_by_species$selection_score,
                                             best_by_species$terminal_name,
                                             best_by_species$tree_tip), , drop = FALSE]
    anchor_count <- min(3L, nrow(best_by_species))
    if (anchor_count < 2L) stop("Genus needs at least two named species anchors: ", genus)
    unit_selected <- do.call(rbind, lapply(seq_len(anchor_count), function(anchor_index) {
      silva_selection_row(
        unit, best_by_species[anchor_index, , drop = FALSE], "genus_mrca", "genus_anchor",
        anchor_index, nrow(pool),
        "Distinct named descendant species; best near-full SILVA SSU per species; up to three anchors",
        FALSE, ""
      )
    }))
  } else if (unit$sequence_source_tier[[1]] == "SILVA genus proxy required") {
    proxy_genus <- if (unit$verified_taxid[[1]] == "912681") "Micrasterias" else
      if (unit$verified_taxid[[1]] == "128215") "Sphagnum" else NA_character_
    if (is.na(proxy_genus)) stop("Unexpected genus-proxy taxon")
    marker <- paste0(";", proxy_genus, ";")
    pool <- taxmap[grepl(marker, taxmap$ncbi_path, fixed = TRUE), , drop = FALSE]
    pool <- pool[vapply(pool$terminal_name, valid_species_name, logical(1), genus = proxy_genus), , drop = FALSE]
    pool$selection_score <- score_silva(pool)
    pool <- pool[pool$region_length >= 900, , drop = FALSE]
    best_by_species <- do.call(rbind, lapply(split(seq_len(nrow(pool)), pool$terminal_name), function(rows) {
      candidate <- pool[rows, , drop = FALSE]
      candidate[order(-candidate$selection_score, candidate$tree_tip), , drop = FALSE][1, , drop = FALSE]
    }))
    best_by_species <- best_by_species[order(-best_by_species$selection_score,
                                             best_by_species$terminal_name,
                                             best_by_species$tree_tip), , drop = FALSE]
    if (nrow(best_by_species) < 3) stop("Proxy genus has fewer than three named species anchors")
    unit_selected <- do.call(rbind, lapply(seq_len(3), function(anchor_index) {
      silva_selection_row(
        unit, best_by_species[anchor_index, , drop = FALSE], "approved_genus_proxy_mrca",
        "proxy_genus_anchor", anchor_index, nrow(pool),
        "User-approved genus proxy; three distinct named SILVA descendant species selected deterministically",
        TRUE,
        paste0("No exact SSU route for verified species; represent by MRCA of ", proxy_genus,
               " anchors and require sensitivity analysis")
      )
    }))
  } else if (unit_rank == "genus") {
    genus <- unit$verified_name[[1]]
    pool <- ncbi[ncbi$verified_taxid == unit$verified_taxid[[1]], , drop = FALSE]
    if (nrow(pool) == 0) stop("No NCBI genus candidates for ", genus)
    pool$selection_score <- ncbi_score(pool)
    pool <- pool[!grepl("(?i)mitochond|chloroplast|plastid", pool$title, perl = TRUE), , drop = FALSE]
    pool$source_species <- vapply(pool$title, extract_title_species, character(1), genus = genus)
    named_pool <- pool[!is.na(pool$source_species), , drop = FALSE]
    if (nrow(named_pool) > 0) {
      best_by_species <- do.call(rbind, lapply(split(seq_len(nrow(named_pool)), named_pool$source_species),
                                              function(rows) {
        candidate <- named_pool[rows, , drop = FALSE]
        candidate[order(-candidate$selection_score, candidate$accession_version), , drop = FALSE][1, , drop = FALSE]
      }))
      best_by_species <- best_by_species[order(-best_by_species$selection_score,
                                               best_by_species$source_species,
                                               best_by_species$accession_version), , drop = FALSE]
    } else {
      best_by_species <- named_pool
    }
    if (nrow(best_by_species) >= 2) {
      anchor_count <- min(3L, nrow(best_by_species))
      unit_selected <- do.call(rbind, lapply(seq_len(anchor_count), function(anchor_index) {
        ncbi_selection_row(
          unit, best_by_species[anchor_index, , drop = FALSE], nrow(pool),
          representation = "genus_mrca", role = "genus_anchor",
          anchor_index = anchor_index,
          selection_rule = paste(
            "Verified NCBI genus absent from SILVA snapshot; select up to three distinct",
            "named descendant-species SSU records by deterministic quality score; use their MRCA"
          )
        )
      }))
    } else {
      pool <- pool[order(-pool$selection_score, pool$accession_version), , drop = FALSE]
      unit_selected <- ncbi_selection_row(
        unit, pool[1, , drop = FALSE], nrow(pool),
        representation = "genus_single_exemplar", role = "genus_exemplar",
        anchor_index = 1,
        selection_rule = paste(
          "Verified NCBI genus absent from SILVA snapshot and fewer than two named",
          "descendant species are available; retain one explicit SSU exemplar as a",
          "provisional limitation rather than calling it a genus MRCA"
        )
      )
    }
    selected_accessions <- unit_selected$accession_version
    rejected <- pool[!pool$accession_version %in% selected_accessions, , drop = FALSE]
    if (nrow(rejected) > 0) {
      rejection_index <- rejection_index + 1L
      rejected$selected_accession <- paste(selected_accessions, collapse = " | ")
      rejected$rejection_reason <- "Not selected for the explicit genus anchor/exemplar representation"
      ncbi_rejected[[rejection_index]] <- rejected
    }
  } else {
    pool <- ncbi[ncbi$verified_taxid == unit$verified_taxid[[1]], , drop = FALSE]
    if (nrow(pool) == 0) stop("No NCBI candidates for unresolved SILVA taxon ", unit$verified_name[[1]])
    pool$selection_score <- ncbi_score(pool)
    pool <- pool[!grepl("(?i)mitochond|chloroplast|plastid", pool$title, perl = TRUE), , drop = FALSE]
    pool <- pool[order(-pool$selection_score, pool$accession_version), , drop = FALSE]
    unit_selected <- ncbi_selection_row(unit, pool[1, , drop = FALSE], nrow(pool))
    if (nrow(pool) > 1) {
      rejection_index <- rejection_index + 1L
      rejected <- pool[-1, , drop = FALSE]
      rejected$selected_accession <- pool$accession_version[[1]]
      rejected$rejection_reason <- "Lower deterministic score than frozen candidate"
      ncbi_rejected[[rejection_index]] <- rejected
    }
  }

  selection_index <- selection_index + 1L
  selected[[selection_index]] <- unit_selected
  mapping_index <- mapping_index + 1L
  unit_map[[mapping_index]] <- data.frame(
    freeze_id = FREEZE_ID,
    taxonomy_release = RELEASE_ID,
    analysis_unit_taxid = unit$verified_taxid[[1]],
    analysis_unit_name = unit$verified_name[[1]],
    phylum = unit$phylum[[1]],
    representation = unique(unit_selected$representation),
    selected_sequences = nrow(unit_selected),
    accessions = paste(unit_selected$accession_version, collapse = " | "),
    source_databases = paste(unique(unit_selected$source_database), collapse = " | "),
    proxy_flag = any(unit_selected$proxy_flag == "TRUE"),
    application_node = if (unique(unit_selected$representation) %in% c(
      "exact_taxon_tip", "genus_single_exemplar"
    )) {
      "selected sequence tip"
    } else {
      "MRCA of selected anchor tips"
    },
    stringsAsFactors = FALSE
  )
}

selected <- do.call(rbind, selected)
unit_map <- do.call(rbind, unit_map)
if (nrow(unit_map) != 105 || anyDuplicated(unit_map$analysis_unit_taxid)) {
  stop("Freeze unit mapping does not contain 105 unique analysis units")
}
if (length(unique(selected$analysis_unit_taxid)) != 105 || anyDuplicated(selected$accession_version)) {
  stop("Freeze has missing analysis units or duplicated accessions")
}
if (sum(unit_map$proxy_flag) != 2) stop("Freeze must contain exactly two proxy analysis units")
if (any(as.numeric(selected$source_record_length) < 500)) stop("A frozen sequence record is unexpectedly short")

freeze_path <- file.path(output_dir, "ssu_accession_candidate_freeze.csv")
unit_map_path <- file.path(output_dir, "analysis_unit_node_mapping.csv")
rejected_path <- file.path(output_dir, "ncbi_rejected_candidates.csv")
summary_path <- file.path(output_dir, "accession_freeze_summary.json")
report_path <- file.path(output_dir, "ACCESSION_FREEZE_REPORT.md")
manifest_path <- file.path(output_dir, "stage_manifest.json")

write.csv(selected[order(selected$phylum, selected$analysis_unit_name,
                         as.numeric(selected$anchor_index)), ], freeze_path,
          row.names = FALSE, quote = TRUE, na = "")
write.csv(unit_map[order(unit_map$phylum, unit_map$analysis_unit_name), ], unit_map_path,
          row.names = FALSE, quote = TRUE, na = "")
if (length(ncbi_rejected) > 0) {
  rejected <- do.call(rbind, ncbi_rejected)
} else {
  rejected <- ncbi[0, , drop = FALSE]
}
write.csv(rejected, rejected_path, row.names = FALSE, quote = TRUE, na = "")

representation_counts <- as.list(table(unit_map$representation))
source_counts <- as.list(table(selected$source_database))
summary <- list(
  taxonomy_release = RELEASE_ID,
  freeze_id = FREEZE_ID,
  stage_id = "figure3_ssu_accession_candidate_freeze",
  status = "candidate_freeze_complete_pending_sequence_and_submitted_s1_validation",
  analysis_units = nrow(unit_map),
  frozen_sequences = nrow(selected),
  analysis_phyla = length(unique(unit_map$phylum)),
  representation_counts = representation_counts,
  sequence_source_counts = source_counts,
  proxy_units = unit_map$analysis_unit_name[unit_map$proxy_flag],
  frozen_submission_warning = paste(
    "Identity input is S1_verified_annotated.csv; direct cell-level reconciliation to the",
    "authoritative submitted SLA_Supplementary_Tables-nature-comms.xlsx remains pending"
  ),
  downstream_not_run = c("sequence extraction", "alignment", "tree inference",
                         "patristic distances", "Figure 3b", "manuscript changes")
)
write_json(summary, summary_path, pretty = TRUE, auto_unbox = TRUE, digits = 16)

report <- c(
  "# Figure 3 SSU accession candidate freeze",
  "",
  paste0("- Taxonomy release: `", RELEASE_ID, "`"),
  paste0("- Freeze ID: `", FREEZE_ID, "`"),
  paste0("- Analysis units: ", nrow(unit_map)),
  paste0("- Selected sequence records: ", nrow(selected)),
  paste0("- Analysis phyla: ", length(unique(unit_map$phylum))),
  paste0("- SILVA records: ", sum(selected$source_database == "SILVA")),
  paste0("- NCBI records: ", sum(selected$source_database == "NCBI Nucleotide")),
  paste0("- Exact taxon tips: ", sum(unit_map$representation == "exact_taxon_tip")),
  paste0("- Genus MRCA units: ", sum(unit_map$representation == "genus_mrca")),
  paste0("- Genus single-exemplar units: ", sum(unit_map$representation == "genus_single_exemplar")),
  paste0("- Approved proxy MRCA units: ", sum(unit_map$representation == "approved_genus_proxy_mrca")),
  "",
  "## Rules",
  "",
  "- NCBI taxonomy remains authoritative for analysis-unit identity and phylum.",
  "- SILVA 138.2 NR99 supplies clipped SSU sequences where the verified taxon is represented.",
  "- NCBI Nucleotide is used only for exact taxa absent from SILVA; an annotated 16S/18S feature must be extracted.",
  "- A genus unit is mapped to the MRCA of two or three explicit descendant-species anchors.",
  "- A genus with fewer than two named descendant-species SSU records remains an explicit provisional single-exemplar unit.",
  "- The two approved species-to-genus proxies remain explicit and require a downstream sensitivity analysis.",
  "",
  "## Authority warning",
  "",
  paste0("This is a candidate freeze, not a release-final freeze. ", summary$frozen_submission_warning, "."),
  "No submitted document, Figure 3 panel, legend, manuscript claim, table, or response was changed."
)
writeLines(report, report_path, useBytes = TRUE)

outputs <- c(freeze_path, unit_map_path, rejected_path, summary_path, report_path)
manifest <- list(
  schema_version = 1,
  stage_id = summary$stage_id,
  taxonomy_release = RELEASE_ID,
  freeze_id = FREEZE_ID,
  status = summary$status,
  generated_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  random_seed = NULL,
  deterministic_selection = TRUE,
  user_approval_context = "User approved the freeze and application tree on 2026-08-04",
  inputs = list(
    readiness = file_record(readiness_path),
    ncbi_ranks = file_record(ncbi_ranks_path),
    silva_coverage = file_record(coverage_path),
    silva_taxmap = file_record(taxmap_path),
    ncbi_candidates = file_record(ncbi_candidates_path)
  ),
  outputs = lapply(outputs, file_record),
  validation_gates = list(
    analysis_units_105 = nrow(unit_map) == 105,
    analysis_phyla_16 = length(unique(unit_map$phylum)) == 16,
    no_duplicate_accessions = !anyDuplicated(selected$accession_version),
    exactly_two_explicit_proxy_units = sum(unit_map$proxy_flag) == 2,
    submitted_s1_reconciliation = FALSE,
    sequences_extracted = FALSE,
    tree_inferred = FALSE
  )
)
write_json(manifest, manifest_path, pretty = TRUE, auto_unbox = TRUE, digits = 16,
           null = "null")

cat(toJSON(summary, pretty = TRUE, auto_unbox = TRUE, digits = 16), "\n")
