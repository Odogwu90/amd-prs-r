#!/usr/bin/env Rscript
# =============================================================================
# 02_qc.R -- Convert VCFs to PLINK 2 format and apply variant / sample QC
# -----------------------------------------------------------------------------
# Purpose : Turn the raw 1000 Genomes chr1 + chr10 VCFs into one QC-passed
#           PLINK fileset, applying decisions D-08 (variant QC) and D-09
#           (sample QC) from docs/decisions.md one filter at a time so the
#           number of variants and samples after every step is logged.
#
# Run     : Rscript R/02_qc.R                 (from the repository root)
#           Requires bin/plink2.exe (Windows) or bin/plink2; see bin/README.md.
#
# Inputs  : data/raw/1000G/ALL.chr1.*.vcf.gz, ALL.chr10.*.vcf.gz
#           data/raw/1000G/integrated_call_samples_v3.20130502.ALL.panel
#
# Outputs : data/processed/1kg_chr1_10_qc.{pgen,pvar,psam}  QC-passed genotypes
#           data/processed/1kg_chr1_10_qc.{bed,bim,fam}     same, PLINK 1 format
#                                                           (bigsnpr reads .bed)
#           data/processed/qc_log.csv                       counts after each step
#           data/processed/samples_qc.tsv                   kept samples + pop labels
#           data/processed/samples_removed_relatedness.tsv  samples dropped by KING
#           data/processed/logs/02_qc/*.log                 PLINK 2 logs per step
#           docs/decisions.md                               block between the
#                                                           <!-- qc-summary --> markers
#
# Filters (D-08, D-09), applied in this order:
#   per chromosome
#     s0  import VCF as-is (all records)               -> baseline count
#     s1  biallelic A/C/G/T SNPs only                  --snps-only just-acgt --max-alleles 2
#     s2  unique chr:pos:ref:alt IDs, drop duplicates  --set-all-var-ids @:#:$r:$a --rm-dup exclude-all
#     s3  minor allele frequency >= 0.01 (pooled)      --maf 0.01
#     s4  per-variant missingness < 0.02               --geno 0.02
#   combined
#     s5  concatenate chr1 + chr10                     --pmerge-list
#     s6  remove 2nd-degree or closer relatives        --king-cutoff 0.0884
#   No Hardy-Weinberg filter (mixed-ancestry pooled sample). No sex check.
#   Strand-ambiguous SNPs are handled at the summary-statistics matching step.
# =============================================================================

suppressPackageStartupMessages({
  library(here)
  library(glue)
  library(dplyr)
  library(readr)
  library(tibble)
  library(data.table)
})

# ---- 0. Settings ------------------------------------------------------------

MAF_MIN     <- 0.01     # D-08
GENO_MAX    <- 0.02     # D-08
KING_CUTOFF <- 0.0884   # D-09: KING-robust kinship threshold for 2nd degree
N_THREADS   <- max(1L, parallel::detectCores() - 1L)
# VCF import runs single-threaded: PLINK 2 (alpha 6.36 and 7.5, Windows
# builds of 10 Sep 2026) segfaults on the multi-threaded import of the
# 1000 Genomes chr1 VCF on this machine; the single-threaded import succeeds.
IMPORT_THREADS <- 1L
# Workspace cap in MiB. PLINK's default reserves half of physical RAM, which
# is more than is actually free on an 8 GB laptop with a browser open.
PLINK_MEMORY_MB <- 1500L

plink2 <- here("bin", if (.Platform$OS.type == "windows") "plink2.exe" else "plink2")
if (!file.exists(plink2)) stop(glue("PLINK 2 not found at {plink2}. See bin/README.md."))

raw_dir <- here("data", "raw", "1000G")
out_dir <- here("data", "processed")
log_dir <- file.path(out_dir, "logs", "02_qc")
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

vcf <- c(
  "1"  = file.path(raw_dir, "ALL.chr1.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"),
  "10" = file.path(raw_dir, "ALL.chr10.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz")
)
panel_file <- file.path(raw_dir, "integrated_call_samples_v3.20130502.ALL.panel")
missing_in <- c(vcf, panel_file)[!file.exists(c(vcf, panel_file))]
if (length(missing_in) > 0) {
  stop("Missing input(s); run `Rscript R/01_download.R` first:\n  ",
       paste(missing_in, collapse = "\n  "))
}

final_prefix <- file.path(out_dir, "1kg_chr1_10_qc")
tmp <- function(...) file.path(out_dir, paste0("tmp_", paste0(..., collapse = "")))

# ---- 1. Helpers -------------------------------------------------------------

# Run PLINK 2 with the given arguments; stop on failure; keep a copy of the log.
run_plink <- function(args, out, threads = N_THREADS) {
  message(glue("[plink2] {paste(args, collapse = ' ')} --out {basename(out)}"))
  res <- system2(plink2, c(args, "--threads", threads, "--memory", PLINK_MEMORY_MB,
                           "--silent", "--out", out),
                 stdout = TRUE, stderr = TRUE)
  status <- attr(res, "status")
  if (!is.null(status) && status != 0) {
    cat(res, sep = "\n")
    stop(glue("PLINK 2 failed at {basename(out)}; see {out}.log"))
  }
  file.copy(paste0(out, ".log"), file.path(log_dir, paste0(basename(out), ".log")),
            overwrite = TRUE)
  invisible(out)
}

# Count data rows in .pvar / .psam (both have '#'-prefixed header lines).
n_variants <- function(prefix) nrow(fread(paste0(prefix, ".pvar"), skip = "#CHROM", select = 1L))
n_samples  <- function(prefix) nrow(fread(paste0(prefix, ".psam"), skip = "#", select = 1L))

# Accumulate one row per QC step.
qc_log <- tibble()
log_step <- function(step, chromosome, filter, plink_options, prefix, prev_prefix = NULL) {
  nv <- n_variants(prefix); ns <- n_samples(prefix)
  pv <- if (is.null(prev_prefix)) NA_integer_ else n_variants(prev_prefix)
  ps <- if (is.null(prev_prefix)) NA_integer_ else n_samples(prev_prefix)
  row <- tibble(step, chromosome, filter, plink_options,
                n_variants = nv, n_samples = ns,
                variants_removed = pv - nv, samples_removed = ps - ns)
  qc_log <<- bind_rows(qc_log, row)
  message(glue("        -> {format(nv, big.mark = ',')} variants, {ns} samples",
               "{if (!is.na(pv)) glue(' (removed {format(pv - nv, big.mark = \",\")} variants, {ps - ns} samples)') else ''}"))
  invisible(row)
}

# ---- 2. Per-chromosome variant QC -----------------------------------------

for (chr in names(vcf)) {
  message(glue("\n=== chromosome {chr} ==="))
  p0 <- tmp("chr", chr, "_s0_import")
  run_plink(c("--vcf", vcf[[chr]], "--make-pgen"), p0, threads = IMPORT_THREADS)
  log_step("s0", chr, "import VCF (all records)", "--vcf --make-pgen", p0)

  p1 <- tmp("chr", chr, "_s1_biallelic_snp")
  run_plink(c("--pfile", p0, "--snps-only", "just-acgt", "--max-alleles", "2", "--make-pgen"), p1)
  log_step("s1", chr, "biallelic A/C/G/T SNPs only", "--snps-only just-acgt --max-alleles 2", p1, p0)

  p2 <- tmp("chr", chr, "_s2_unique_ids")
  run_plink(c("--pfile", p1, "--set-all-var-ids", "@:#:$r:$a", "--new-id-max-allele-len", "10", "truncate",
              "--rm-dup", "exclude-all", "--make-pgen"), p2)
  log_step("s2", chr, "unique chr:pos:ref:alt IDs; exact duplicates removed",
           "--set-all-var-ids @:#:$r:$a --rm-dup exclude-all", p2, p1)

  p3 <- tmp("chr", chr, "_s3_maf")
  run_plink(c("--pfile", p2, "--maf", MAF_MIN, "--make-pgen"), p3)
  log_step("s3", chr, glue("MAF >= {MAF_MIN} (pooled sample)"), glue("--maf {MAF_MIN}"), p3, p2)

  p4 <- tmp("chr", chr, "_s4_geno")
  run_plink(c("--pfile", p3, "--geno", GENO_MAX, "--make-pgen"), p4)
  log_step("s4", chr, glue("variant missingness < {GENO_MAX}"), glue("--geno {GENO_MAX}"), p4, p3)
}

# ---- 3. Concatenate chromosomes, then sample QC -----------------------------

message("\n=== chr1 + chr10 combined ===")
merge_list <- file.path(out_dir, "tmp_pmerge_list.txt")
writeLines(c(tmp("chr1_s4_geno"), tmp("chr10_s4_geno")), merge_list)
p5 <- tmp("all_s5_merged")
run_plink(c("--pmerge-list", merge_list, "--make-pgen"), p5)
log_step("s5", "1+10", "concatenate chr1 and chr10", "--pmerge-list", p5)

# KING-robust kinship; PLINK greedily drops one sample from each pair above
# the cutoff. The .king.cutoff.out.id file lists the removed samples.
run_plink(c("--pfile", p5, "--king-cutoff", KING_CUTOFF, "--make-pgen"), final_prefix)
log_step("s6", "1+10", glue("remove relatives, KING kinship > {KING_CUTOFF} (2nd degree)"),
         glue("--king-cutoff {KING_CUTOFF}"), final_prefix, p5)

# PLINK 1 binary copy for bigsnpr (snp_readBed).
run_plink(c("--pfile", final_prefix, "--make-bed"), final_prefix)

# ---- 4. Sample tables ------------------------------------------------------

panel <- fread(panel_file, select = 1:4, col.names = c("sample", "pop", "super_pop", "gender"))
stopifnot(nrow(panel) == 2504)

kept <- fread(paste0(final_prefix, ".psam"), skip = "#", select = 1L, col.names = "sample")
samples_qc <- kept |> left_join(panel, by = "sample")
if (anyNA(samples_qc$super_pop)) stop("Some kept samples are missing from the panel file.")
write_tsv(samples_qc, file.path(out_dir, "samples_qc.tsv"))

removed_file <- paste0(final_prefix, ".king.cutoff.out.id")
removed <- if (file.exists(removed_file) && length(readLines(removed_file)) > 1) {
  fread(removed_file, skip = "#", select = 1L, col.names = "sample") |> left_join(panel, by = "sample")
} else {
  tibble(sample = character(), pop = character(), super_pop = character(), gender = character())
}
write_tsv(removed, file.path(out_dir, "samples_removed_relatedness.tsv"))

# ---- 5. QC log -------------------------------------------------------------

write_csv(qc_log, file.path(out_dir, "qc_log.csv"))
message("\nQC log:")
print(qc_log |> select(-plink_options), n = Inf, width = Inf)

# ---- 6. Write the outcome block into docs/decisions.md ---------------------

fmt <- function(x) ifelse(is.na(x), "", format(x, big.mark = ","))
per_pop <- samples_qc |> count(super_pop, name = "n_kept") |>
  left_join(removed |> count(super_pop, name = "n_removed"), by = "super_pop") |>
  mutate(n_removed = coalesce(n_removed, 0L))

block <- c(
  "<!-- qc-summary:start -->",
  glue("_Generated by `R/02_qc.R` on {format(Sys.time(), '%Y-%m-%d %H:%M')}; ",
       "PLINK {trimws(sub('^PLINK ', '', system2(plink2, '--version', stdout = TRUE)[1]))}. ",
       "Do not edit this block by hand._"),
  "",
  "| Step | Chr | Filter | Variants | Samples | Variants removed | Samples removed |",
  "|---|---|---|---:|---:|---:|---:|",
  glue_data(qc_log, "| {step} | {chromosome} | {filter} | {fmt(n_variants)} | {n_samples} ",
            "| {fmt(variants_removed)} | {fmt(samples_removed)} |"),
  "",
  glue("Final fileset: `data/processed/1kg_chr1_10_qc` with ",
       "{fmt(n_variants(final_prefix))} variants and {n_samples(final_prefix)} samples."),
  "",
  "Samples kept / removed by super-population:",
  "",
  "| Super-population | Kept | Removed (relatedness) |",
  "|---|---:|---:|",
  glue_data(per_pop, "| {super_pop} | {n_kept} | {n_removed} |"),
  if (nrow(removed) > 0) c("", "Removed samples: ",
                           paste0("`", removed$sample, "` (", removed$pop, ")", collapse = ", ")) else NULL,
  "<!-- qc-summary:end -->"
)

decisions <- here("docs", "decisions.md")
lines <- read_lines(decisions)
s <- which(lines == "<!-- qc-summary:start -->"); e <- which(lines == "<!-- qc-summary:end -->")
if (length(s) == 1 && length(e) == 1 && e > s) {
  # seq_len() rather than seq(e + 1, n): when the end marker is the last line,
  # seq(n + 1, n) would count backwards and duplicate lines.
  after <- lines[e + seq_len(length(lines) - e)]
  lines <- c(lines[seq_len(s - 1)], block, after)
} else {
  lines <- c(lines, "", "### QC outcome (written by `R/02_qc.R`)", "", block)
}
write_lines(lines, decisions)
message(glue("Updated {decisions}"))

# ---- 7. Clean up intermediates (logs are kept in logs/02_qc) ---------------

unlink(list.files(out_dir, pattern = "^tmp_", full.names = TRUE))
message("Done. Intermediate tmp_* files removed.")
