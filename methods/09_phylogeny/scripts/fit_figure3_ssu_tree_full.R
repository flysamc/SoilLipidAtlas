#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(Biostrings); library(phangorn); library(ape); library(jsonlite); library(digest)})
parse_args <- function(x) { out <- list(); for (v in x) { p <- strsplit(sub("^--", "", v), "=", fixed=TRUE)[[1]]; out[[p[1]]] <- paste(p[-1], collapse="=") }; out }
a <- parse_args(commandArgs(trailingOnly=TRUE))
alignment_path <- normalizePath(a$alignment, mustWork=TRUE); tip_map_path <- normalizePath(a$`tip-map`, mustWork=TRUE)
dir.create(a$`output-dir`, recursive=TRUE, showWarnings=FALSE); out <- normalizePath(a$`output-dir`, mustWork=TRUE)
seed <- 20260804L; set.seed(seed)
sha <- function(p) digest(p, algo="sha256", file=TRUE, serialize=FALSE)
rec <- function(p) { p <- normalizePath(p, mustWork=TRUE); list(path=p, bytes=unname(file.info(p)$size), sha256=sha(p)) }
aln <- readDNAStringSet(alignment_path); map <- read.csv(tip_map_path, check.names=FALSE, colClasses="character")
if (!identical(names(aln), map$tip_id) || length(aln) != 148) stop("Alignment/map mismatch")
dat <- phyDat(as.matrix(aln), type="DNA"); d <- dist.ml(dat, model="F81"); start <- NJ(d)
message("Optimizing checkpointed GTR+Gamma NNI tree...")
fit0 <- pml(start, dat, k=4)
fit <- optim.pml(fit0, model="GTR", optNni=TRUE, optBf=TRUE, optQ=TRUE,
                 optGamma=TRUE, optEdge=TRUE, rearrangement="NNI",
                 control=pml.control(epsilon=1e-7, maxit=20, trace=1))
tree <- reorder.phylo(fit$tree, "cladewise")
if (!setequal(tree$tip.label, map$tip_id) || any(!is.finite(tree$edge.length)) || any(tree$edge.length < 0)) stop("Tree validation failed")
fit_path <- file.path(out,"organism_ssu_ml_fit_nni_v3.rds"); tree_path <- file.path(out,"organism_ssu_ml_gtr_gamma_nni_v3_unrooted.nwk")
start_path <- file.path(out,"organism_ssu_start_nj_v3.nwk"); summary_path <- file.path(out,"ml_fit_nni_v3_summary.json")
saveRDS(fit, fit_path); write.tree(tree, tree_path); write.tree(start, start_path)
summary <- list(taxonomy_release="ncbi-phylum-2026-08-04-v1", freeze_id="figure3-ssu-curated-freeze-2026-08-04-v3",
 status="checkpointed_gtr_gamma_nni_fit_complete", model="GTR+Gamma (4 categories)", topology_search="NNI",
 log_likelihood=unname(fit$logLik), sequences=Ntip(tree), sites=attr(dat,"nr"), random_seed=seed,
 stochastic_search_note="Attempted stochastic ratchet was operationally unbounded in phangorn and was not promoted; NNI is the declared practical full-pass search.")
write_json(summary, summary_path, pretty=TRUE, auto_unbox=TRUE, digits=16)
outputs <- c(fit_path,tree_path,start_path,summary_path)
manifest <- list(schema_version=1, stage_id="figure3_ssu_ml_fit_nni_v3", taxonomy_release=summary$taxonomy_release,
 freeze_id=summary$freeze_id,status=summary$status,generated_at_utc=format(Sys.time(),tz="UTC",usetz=TRUE),random_seed=seed,
 inputs=list(alignment=rec(alignment_path),tip_map=rec(tip_map_path)),outputs=lapply(outputs,rec),
 validation_gates=list(tips_148=Ntip(tree)==148,finite_nonnegative_branches=all(is.finite(tree$edge.length)&tree$edge.length>=0),fit_checkpoint_written=file.exists(fit_path)))
write_json(manifest,file.path(out,"stage_manifest.json"),pretty=TRUE,auto_unbox=TRUE,digits=16); cat(toJSON(summary,pretty=TRUE,auto_unbox=TRUE,digits=16),"\n")
