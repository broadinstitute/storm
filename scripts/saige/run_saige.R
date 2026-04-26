#!/usr/bin/env Rscript
# Storm: driver for SAIGE Step 1 (null model) and Step 2 (single-variant test) on Storm DS VCFs.
#
# Requires:
#   - R package SAIGE installed, OR extdata on disk and SAIGE_EXTDATA set to that folder
#   - tabix index (.csi) next to each Step 2 VCF unless you set vcf_index_* in the config
#
# Pixi (SAIGE docs): run from the SAIGE checkout and set STORM_SAIGE_RSCRIPT to the absolute
# path of that environment's Rscript so nested Step 1/2 calls use the same interpreter, e.g.:
#   export STORM_SAIGE_RSCRIPT="$(cd \"$SAIGE_SRC\" && pixi run bash -lc 'command -v Rscript')"
#
# Usage:
#   Rscript scripts/saige/run_saige.R --config scripts/saige/example_chr22.config
#   Rscript scripts/saige/run_saige.R --config my.config --repo-root /path/to/storm
#   Rscript scripts/saige/run_saige.R --help
#
# You must supply Step 1 relatedness inputs: either plink_prefix (bed/bim/fam prefix)
# or sparse_grm_file + sparse_grm_sample_id_file (see SAIGE Step 1 docs).

`%||%` <- function(a, b) if (!is.null(a) && nzchar(as.character(a))) a else b

this_script_path <- function() {
  ca <- commandArgs(trailingOnly = FALSE)
  f <- ca[grep("^--file=", ca)]
  if (length(f) == 0L) return(NULL)
  normalizePath(sub("^--file=", "", f[1L]), mustWork = TRUE)
}

read_kv_config <- function(path) {
  lines <- readLines(path, warn = FALSE)
  lines <- trimws(sub("#.*$", "", lines))
  lines <- lines[nzchar(lines) & !startsWith(lines, "#")]
  n <- length(lines)
  keys <- character(n)
  vals <- character(n)
  for (i in seq_len(n)) {
    line <- lines[i]
    eq <- regexpr("=", line, fixed = TRUE)[1L]
    if (eq < 1L) stop("Bad config line (expected key=value): ", line)
    keys[i] <- trimws(substr(line, 1L, eq - 1L))
    vals[i] <- trimws(substr(line, eq + 1L, nchar(line)))
  }
  stats::setNames(as.list(vals), keys)
}

as_logical <- function(x, default = FALSE) {
  if (is.null(x) || !nzchar(x)) return(default)
  toupper(x) %in% c("1", "TRUE", "T", "YES", "Y")
}

as_int <- function(x, default) {
  if (is.null(x) || !nzchar(x)) return(default)
  as.integer(x)
}

resolve_path <- function(repo_root, p) {
  if (is.null(p)) return("")
  p <- trimws(as.character(p))
  if (!nzchar(p)) return("")
  if (startsWith(p, "/")) return(normalizePath(p, mustWork = FALSE))
  normalizePath(file.path(repo_root, p), mustWork = FALSE)
}

find_saige_extdata_script <- function(filename) {
  ext <- Sys.getenv("SAIGE_EXTDATA", "")
  if (nzchar(ext)) {
    alt <- file.path(ext, filename)
    if (file.exists(alt)) return(normalizePath(alt, mustWork = TRUE))
  }
  if (requireNamespace("SAIGE", quietly = TRUE)) {
    bundled <- system.file("extdata", filename, package = "SAIGE", mustWork = FALSE)
    if (nzchar(bundled) && file.exists(bundled)) return(bundled)
  }
  stop(
    "Cannot find ", filename, ". Set SAIGE_EXTDATA to the SAIGE extdata/ folder, or install the SAIGE R package.\n",
    "  Example: export SAIGE_EXTDATA=/path/to/SAIGE/extdata"
  )
}

run_rscript <- function(script_path, args) {
  rscript <- Sys.getenv("STORM_SAIGE_RSCRIPT", "Rscript")
  message("----\n", rscript, " ", script_path, " ...\n----")
  code <- system2(rscript, c(script_path, args), stdout = "", stderr = "")
  if (code != 0L) stop(rscript, " exited with code ", code)
  invisible(TRUE)
}

usage <- function() {
  cat("Usage:\n  Rscript run_saige.R --config FILE [--repo-root DIR]\n\n")
  cat("Optional: STORM_SAIGE_RSCRIPT=/abs/path/to/Rscript (Pixi SAIGE env) for nested Step 1/2.\n")
  cat("See scripts/saige/example_chr22.config for keys.\n")
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0L || any(args %in% c("-h", "--help"))) {
  usage()
  quit(status = 0L)
}

cfg_path <- NULL
repo_root_cli <- NULL
i <- 1L
while (i <= length(args)) {
  a <- args[[i]]
  if (a == "--config") {
    cfg_path <- args[[i + 1L]]
    i <- i + 2L
  } else if (startsWith(a, "--config=")) {
    cfg_path <- sub("^--config=", "", a)
    i <- i + 1L
  } else if (a == "--repo-root") {
    repo_root_cli <- args[[i + 1L]]
    i <- i + 2L
  } else if (startsWith(a, "--repo-root=")) {
    repo_root_cli <- sub("^--repo-root=", "", a)
    i <- i + 1L
  } else {
    stop("Unknown argument: ", a)
  }
}
if (is.null(cfg_path)) stop("Missing --config FILE")

cfg <- read_kv_config(cfg_path)

script_file <- this_script_path()
repo_root <- cfg$repo_root %||% repo_root_cli
if (is.null(repo_root) || !nzchar(repo_root)) {
  if (!is.null(script_file)) {
    repo_root <- normalizePath(file.path(dirname(script_file), "..", ".."), mustWork = TRUE)
  } else {
    repo_root <- normalizePath(getwd(), mustWork = TRUE)
  }
} else {
  repo_root <- normalizePath(repo_root, mustWork = TRUE)
}

pheno_file <- resolve_path(repo_root, cfg$pheno_file %||% "")
pheno_col <- cfg$pheno_col %||% "pheno_binary"
trait_type <- cfg$trait_type %||% "binary"
sample_id_col <- cfg$sample_id_col %||% "IID"
covar_col_list <- cfg$covar_col_list %||% ""
q_covar_col_list <- cfg$q_covar_col_list %||% ""
chrom <- cfg$chrom %||% "chr22"
n_threads <- as_int(cfg$n_threads, 4L)

output_dir <- resolve_path(repo_root, cfg$output_dir %||% "scratch/saige_runs")
run_prefix <- cfg$run_prefix %||% "chr22_smoke"
step1_prefix <- cfg$step1_output_prefix
if (is.null(step1_prefix) || !nzchar(step1_prefix)) {
  step1_prefix <- file.path(output_dir, run_prefix)
} else {
  step1_prefix <- resolve_path(repo_root, step1_prefix)
}
dir.create(dirname(step1_prefix), recursive = TRUE, showWarnings = FALSE)

plink_prefix <- resolve_path(repo_root, cfg$plink_prefix %||% "")
sparse_grm_file <- resolve_path(repo_root, cfg$sparse_grm_file %||% "")
sparse_grm_sample_id_file <- resolve_path(repo_root, cfg$sparse_grm_sample_id_file %||% "")

use_sparse_grm <- as_logical(cfg$use_sparse_grm, FALSE)
skip_variance_ratio <- as_logical(cfg$skip_variance_ratio_estimation, FALSE)
inv_normalize <- as_logical(cfg$inv_normalize, FALSE)
skip_step1 <- as_logical(cfg$skip_step1, FALSE)
skip_step2 <- as_logical(cfg$skip_step2, FALSE)

vcf_standard <- resolve_path(repo_root, cfg$vcf_standard_sv %||% "")
vcf_tr <- resolve_path(repo_root, cfg$vcf_tr_quantitative %||% "")
idx_std <- cfg$vcf_index_standard_sv
idx_tr <- cfg$vcf_index_tr_quantitative
if (!is.null(idx_std) && nzchar(idx_std)) idx_std <- resolve_path(repo_root, idx_std) else idx_std <- paste0(vcf_standard, ".csi")
if (!is.null(idx_tr) && nzchar(idx_tr)) idx_tr <- resolve_path(repo_root, idx_tr) else idx_tr <- paste0(vcf_tr, ".csi")

min_maf <- cfg$min_maf %||% "0"
min_mac <- cfg$min_mac %||% "20"
loco <- as_logical(cfg$loco, TRUE)
is_fast_test <- as_logical(cfg$is_fast_test, FALSE)

sparse_grm_step2 <- resolve_path(repo_root, cfg$sparse_grm_file_step2 %||% (cfg$sparse_grm_file %||% ""))
sparse_grm_sample_step2 <- resolve_path(
  repo_root,
  cfg$sparse_grm_sample_id_file_step2 %||% (cfg$sparse_grm_sample_id_file %||% "")
)

# --- validate inputs ---
stopifnot(file.exists(pheno_file))
if (!skip_step2) {
  stopifnot(file.exists(vcf_standard), file.exists(vcf_tr))
  stopifnot(file.exists(idx_std), file.exists(idx_tr))
}

has_plink <- nzchar(plink_prefix)
has_sparse <- nzchar(sparse_grm_file) && nzchar(sparse_grm_sample_id_file)

if (!skip_step1) {
  if (!has_plink && !has_sparse) {
    stop("Step 1 needs relatedness: set plink_prefix=... and/or sparse_grm_file= + sparse_grm_sample_id_file= in the config (see SAIGE Step 1 docs).")
  }
  if (has_plink) {
    for (suf in c(".bed", ".bim", ".fam")) {
      if (!file.exists(paste0(plink_prefix, suf))) {
        stop("Missing PLINK file: ", plink_prefix, suf)
      }
    }
  }
  if (has_sparse) {
    stopifnot(file.exists(sparse_grm_file), file.exists(sparse_grm_sample_id_file))
  }
  if (!skip_variance_ratio && !has_plink) {
    message("Note: skip_variance_ratio_estimation=FALSE usually requires plink_prefix for variance ratio markers (SAIGE docs).")
  }
}

step1_script <- find_saige_extdata_script("step1_fitNULLGLMM.R")
step2_script <- find_saige_extdata_script("step2_SPAtests.R")

model_rda <- paste0(step1_prefix, ".rda")
variance_ratio <- paste0(step1_prefix, ".varianceRatio.txt")

if (!skip_step1) {
  s1 <- c(
    paste0("--phenoFile=", pheno_file),
    paste0("--phenoCol=", pheno_col),
    paste0("--traitType=", trait_type),
    paste0("--sampleIDColinphenoFile=", sample_id_col),
    paste0("--outputPrefix=", step1_prefix),
    paste0("--nThreads=", n_threads),
    "--IsOverwriteVarianceRatioFile=TRUE"
  )
  if (nzchar(covar_col_list)) s1 <- c(s1, paste0("--covarColList=", covar_col_list))
  if (nzchar(q_covar_col_list)) s1 <- c(s1, paste0("--qCovarColList=", q_covar_col_list))
  if (inv_normalize) s1 <- c(s1, "--invNormalize=TRUE")
  if (has_plink) s1 <- c(s1, paste0("--plinkFile=", plink_prefix))
  if (has_sparse) {
    s1 <- c(
      s1,
      "--useSparseGRMtoFitNULL=TRUE",
      paste0("--sparseGRMFile=", sparse_grm_file),
      paste0("--sparseGRMSampleIDFile=", sparse_grm_sample_id_file)
    )
  }
  if (skip_variance_ratio) {
    s1 <- c(s1, "--skipVarianceRatioEstimation=TRUE")
  } else {
    s1 <- c(s1, "--skipVarianceRatioEstimation=FALSE")
  }
  run_rscript(step1_script, s1)
} else {
  message("Skipping Step 1 (--skip_step1=TRUE); expecting existing model: ", model_rda)
  stopifnot(file.exists(model_rda))
}

if (!file.exists(variance_ratio)) {
  message("Warning: variance ratio file not found: ", variance_ratio, " (Step 2 may still run if your Step 1 mode does not estimate it).")
}

if (!skip_step2) {
  step2_runs <- list(
    list(name = "standard_sv", vcf = vcf_standard, index = idx_std),
    list(name = "tr_quantitative", vcf = vcf_tr, index = idx_tr)
  )
  for (run in step2_runs) {
    out_txt <- paste0(step1_prefix, ".step2.", run$name, ".txt")
    s2 <- list(
      paste0("--vcfFile=", run$vcf),
      paste0("--vcfFileIndex=", run$index),
      "--vcfField=DS",
      paste0("--SAIGEOutputFile=", out_txt),
      paste0("--chrom=", chrom),
      paste0("--minMAF=", min_maf),
      paste0("--minMAC=", min_mac),
      paste0("--GMMATmodelFile=", model_rda),
      if (file.exists(variance_ratio)) paste0("--varianceRatioFile=", variance_ratio) else NULL,
      if (loco) "--LOCO=TRUE" else "--LOCO=FALSE",
      if (is_fast_test) "--is_fastTest=TRUE" else NULL
    )
    s2 <- unlist(Filter(Negate(is.null), s2))
    if (has_sparse && nzchar(sparse_grm_step2) && nzchar(sparse_grm_sample_step2) &&
        file.exists(sparse_grm_step2) && file.exists(sparse_grm_sample_step2)) {
      s2 <- c(
        s2,
        paste0("--sparseGRMFile=", sparse_grm_step2),
        paste0("--sparseGRMSampleIDFile=", sparse_grm_sample_step2)
      )
    }
    run_rscript(step2_script, s2)
  }
} else {
  message("Skipping Step 2 (--skip_step2=TRUE).")
}

message("Done.\n  Step 1 prefix: ", step1_prefix, "\n  Model: ", model_rda)
