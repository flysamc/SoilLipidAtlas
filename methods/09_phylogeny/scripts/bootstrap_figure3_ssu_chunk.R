#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(phangorn); library(ape); library(jsonlite); library(digest)})
parse_args <- function(x) { out <- list(); for (v in x) { p <- strsplit(sub("^--", "", v), "=", fixed=TRUE)[[1]]; out[[p[1]]] <- paste(p[-1], collapse="=") }; out }
a <- parse_args(commandArgs(trailingOnly=TRUE)); fit_path <- normalizePath(a$fit,mustWork=TRUE)
chunk <- as.integer(a$chunk); reps <- as.integer(a$reps); seed <- 20260804L + chunk * 1000L
dir.create(a$`output-dir`,recursive=TRUE,showWarnings=FALSE); out <- normalizePath(a$`output-dir`,mustWork=TRUE)
fit <- readRDS(fit_path); set.seed(seed); message("Bootstrap chunk ",chunk,": ",reps," replicates")
trees <- bootstrap.pml(fit,bs=reps,optNni=TRUE,multicore=FALSE,control=pml.control(trace=0))
if (length(trees)!=reps) stop("Bootstrap count mismatch")
path <- file.path(out,sprintf("bootstrap_chunk_%02d.rds",chunk)); saveRDS(trees,path)
sha <- function(p) digest(p,algo="sha256",file=TRUE,serialize=FALSE)
summary <- list(stage_id="figure3_ssu_bootstrap_chunk",chunk=chunk,replicates=reps,seed=seed,status="bootstrap_chunk_complete",path=normalizePath(path),sha256=sha(path))
write_json(summary,file.path(out,sprintf("bootstrap_chunk_%02d.json",chunk)),pretty=TRUE,auto_unbox=TRUE,digits=16); cat(toJSON(summary,pretty=TRUE,auto_unbox=TRUE,digits=16),"\n")
