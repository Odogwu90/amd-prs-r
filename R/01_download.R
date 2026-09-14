#!/usr/bin/env Rscript
# =============================================================================
# 01_download.R -- Download the public input data and verify checksums
# -----------------------------------------------------------------------------
# Purpose : Fetch every raw input for the AMD PRS project from its official
#           source, verify file integrity, and record what was downloaded
#           (version, date, checksum) in docs/data_sources.md.
#
# Run     : Rscript R/01_download.R        (from the repository root)
#
# Inputs  : Internet access only. No local inputs.
#
# Outputs : data/raw/1000G/
#             ALL.chr1.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz (+ .tbi)
#             ALL.chr10.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz (+ .tbi)
#             integrated_call_samples_v3.20130502.ALL.panel
#           data/raw/iamdgc/
#             26691988-GCST003219-EFO_0001365-build37.f.tsv.gz
#           data/raw/checksums.tsv        one row per file: url, date, md5, status
#           docs/data_sources.md          the block between the
#                                         <!-- download-log:start/end --> markers
#                                         is regenerated; hand-written text is kept
#
# Notes   : * Idempotent. A file that already exists with a matching md5 is
#             skipped, so an interrupted run can simply be re-run.
#           * Downloads go to <file>.part first and are renamed only after the
#             checksum passes, so a half-written file never looks complete.
#           * Total download is about 2.1 GB (the chr1 VCF alone is 1.1 GB).
#             Expect 10-60 min depending on connection.
# =============================================================================

suppressPackageStartupMessages({
  library(here)
  library(glue)
  library(dplyr)
  library(readr)
  library(tibble)
})

# download.file() gives up after 60 s by default; that is far too short for a
# 1 GB VCF, so allow up to 6 hours per file.
options(timeout = 6 * 60 * 60)

# ---- 1. Manifest of files to fetch -----------------------------------------
#
# Where the expected md5 values come from:
#   * 1000 Genomes: copied verbatim from the project's own file index
#     https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/current.tree
#     (index dated 2024-08-16; consulted 2026-09-14). This is the official
#     checksum source: the release directory itself ships no md5 files.
#   * GWAS Catalog GCST003219: the Catalog publishes an md5 only for the GRCh38
#     harmonised file (in <file>-meta.yaml). The GRCh37 formatted file used here
#     has no published checksum, so md5_expected is NA. On the first download the
#     script records the md5 it computes in data/raw/checksums.tsv and verifies
#     against that value on every later run. That detects later corruption but,
#     unlike an official checksum, cannot detect a bad *first* download.

kg_base <- "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502"
gc_base <- paste0("https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/",
                  "GCST003001-GCST004000/GCST003219")

manifest <- tribble(
  ~dataset, ~dest_dir, ~file,                                                                          ~md5_expected,
  "1000G",  "1000G",   "ALL.chr1.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz",      "76f1d3fef27c6c3f451cdfc515250a0e",
  "1000G",  "1000G",   "ALL.chr1.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz.tbi",  "8c40851a4e9d3c6c9a09461961c6e4c4",
  "1000G",  "1000G",   "ALL.chr10.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz",     "45814b06e651f8fb409364aa6d65257e",
  "1000G",  "1000G",   "ALL.chr10.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz.tbi", "266a7200332b971c741b1520e03bda80",
  "1000G",  "1000G",   "integrated_call_samples_v3.20130502.ALL.panel",                                   "7ee5675553088230530a7fe88c22f201",
  "IAMDGC", "iamdgc",  "harmonised/26691988-GCST003219-EFO_0001365-build37.f.tsv.gz",                    NA_character_
) |>
  mutate(
    url  = if_else(dataset == "1000G",
                   as.character(glue("{kg_base}/{file}")),
                   as.character(glue("{gc_base}/{file}"))),
    dest = here("data", "raw", dest_dir, basename(file))
  )

checksum_log <- here("data", "raw", "checksums.tsv")
data_sources <- here("docs", "data_sources.md")

# ---- 2. Helpers -------------------------------------------------------------

md5_of <- function(path) unname(tools::md5sum(path))

# md5 recorded on a previous run (used for files with no official checksum)
previous_md5 <- function(file_name) {
  if (!file.exists(checksum_log)) return(NA_character_)
  prev <- read_tsv(checksum_log, show_col_types = FALSE,
                   col_types = cols(.default = "c"))
  hit  <- prev |> filter(file == file_name, !is.na(md5))
  if (nrow(hit) == 0) NA_character_ else hit$md5[[nrow(hit)]]
}

row_result <- function(dataset, file, url, dest, md5, md5_source, status) {
  tibble(dataset = dataset, file = basename(file), url = url,
         date = Sys.Date(),
         size_bytes = if (file.exists(dest)) file.size(dest) else NA_real_,
         md5 = md5, md5_source = md5_source, status = status)
}

fetch_one <- function(dataset, dest_dir, file, md5_expected, url, dest) {
  dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)
  official   <- !is.na(md5_expected)
  reference  <- if (official) md5_expected else previous_md5(basename(file))
  ref_source <- if (official) "official" else "recorded on first download"

  # Already downloaded and intact? Then do nothing.
  if (file.exists(dest)) {
    md5 <- md5_of(dest)
    if (!is.na(reference) && identical(md5, reference)) {
      message(glue("[skip]  {basename(dest)}  md5 OK ({ref_source})"))
      return(row_result(dataset, file, url, dest, md5, ref_source, "present_verified"))
    }
    if (is.na(reference)) {
      message(glue("[keep]  {basename(dest)}  exists; no reference md5 -> recording it"))
      return(row_result(dataset, file, url, dest, md5, "recorded (existing file)",
                        "present_unverified"))
    }
    message(glue("[redo]  {basename(dest)}  md5 mismatch -> re-downloading"))
    unlink(dest)
  }

  # Download to a .part file, verify, then rename.
  part <- paste0(dest, ".part")
  message(glue("[get]   {url}"))
  t0 <- Sys.time()
  ok <- tryCatch(
    utils::download.file(url, part, mode = "wb", method = "libcurl", quiet = FALSE) == 0,
    error = function(e) { message("        ", conditionMessage(e)); FALSE }
  )
  if (!ok) {
    unlink(part)
    return(row_result(dataset, file, url, dest, NA_character_, ref_source,
                      "download_failed"))
  }
  md5 <- md5_of(part)
  if (!is.na(reference) && !identical(md5, reference)) {
    unlink(part)
    stop(glue("Checksum FAILED for {basename(dest)}\n",
              "  expected {reference}\n  got      {md5}"))
  }
  file.rename(part, dest)
  mins <- round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1)
  message(glue("[ok]    {basename(dest)}  {round(file.size(dest) / 1e6)} MB ",
               "in {mins} min; md5 {md5} ({ref_source})"))
  row_result(dataset, file, url, dest, md5, ref_source, "downloaded_verified")
}

# ---- 3. Download everything ------------------------------------------------

results <- purrr::pmap(manifest, fetch_one) |> bind_rows()

write_tsv(results, checksum_log)
message(glue("\nWrote {checksum_log}"))
print(select(results, dataset, file, status, md5), n = Inf, width = Inf)

# ---- 4. Update the download log in docs/data_sources.md --------------------
#
# Only the block between the two marker comments is regenerated, so the
# hand-written description of each source (licence, citation, build) is kept.

fmt_mb  <- function(x) ifelse(is.na(x), "", format(round(x / 1e6, 1), nsmall = 1))
fmt_md5 <- function(x) ifelse(is.na(x), "", paste0("`", x, "`"))

log_block <- c(
  "<!-- download-log:start -->",
  glue("_Generated by `R/01_download.R` on {format(Sys.time(), '%Y-%m-%d %H:%M')}; ",
       "R {getRversion()}. Do not edit this block by hand._"),
  "",
  "| Dataset | File | Size (MB) | md5 | md5 source | Status | Date |",
  "|---|---|---:|---|---|---|---|",
  glue_data(results,
            "| {dataset} | `{file}` | {fmt_mb(size_bytes)} | {fmt_md5(md5)} ",
            "| {md5_source} | {status} | {date} |"),
  "<!-- download-log:end -->"
)

if (file.exists(data_sources)) {
  lines <- read_lines(data_sources)
  s <- which(lines == "<!-- download-log:start -->")
  e <- which(lines == "<!-- download-log:end -->")
  if (length(s) == 1 && length(e) == 1 && e > s) {
    # seq_len() rather than seq(e + 1, n): when the end marker is the last
    # line, seq(n + 1, n) would count *backwards* and duplicate lines.
    after <- lines[e + seq_len(length(lines) - e)]
    lines <- c(lines[seq_len(s - 1)], log_block, after)
  } else {
    lines <- c(lines, "", "## Download log", "", log_block)
  }
} else {
  lines <- c("# Data sources", "", "## Download log", "", log_block)
}
write_lines(lines, data_sources)
message(glue("Updated {data_sources}"))

if (any(results$status == "download_failed")) {
  stop("One or more downloads failed; see messages above and re-run.")
}
