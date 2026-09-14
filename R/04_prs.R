#!/usr/bin/env Rscript
# =============================================================================
# 04_prs.R -- Build the AMD polygenic risk scores (D-11)
# -----------------------------------------------------------------------------
# Purpose : Harmonise the IAMDGC 2016 summary statistics with the QC-passed
#           1000 Genomes genotypes and compute two scores per sample:
#             (a) "published" : the Fritsche 2016 independent variants on
#                 chr1/chr10 (GWAS Catalog GCST003219 curated list), weighted
#                 by the reported log odds ratio;
#             (b) "C+T"       : clumping + thresholding on all matched
#                 variants at p < 5e-8, 1e-5, 1e-3, 0.05.
#           Both are standardised (mean 0, SD 1) in the full 2,484-sample panel.
#
#           Effect sizes: the public IAMDGC release gives only a p-value and a
#           direction sign per variant (no beta/OR). For (b) each weight is an
#           approximate log-odds recovered from the signed z-score:
#             z    = sign * qnorm(p / 2, lower = FALSE)      (computed on log p,
#                                                             p goes below 1e-700)
#             beta = z * 2 / sqrt(n_eff * 2 f (1 - f)),  n_eff = 4 Nca Nco / (Nca + Nco)
#           with f = frequency of the effect allele in the 1000 Genomes EUR
#           samples (the GWAS was European). The same formula is checked
#           against the reported odds ratios of the published variants.
#
# Run     : Rscript R/04_prs.R                 (from the repository root)
#
# Inputs  : data/raw/iamdgc/Fritsche-26691988.txt.gz          author original
#           data/raw/iamdgc/gcst003219_associations.json      Catalog curated ORs
#           data/processed/1kg_chr1_10_qc.rds (+ .bk)         from 03_pca.R
#           data/processed/samples_qc.tsv
#
# Outputs : data/processed/sumstats_matched.tsv.gz   chr1/10 variants matched to genotypes
#           data/processed/prs_weights.tsv           variant, effect allele, weight, score(s)
#           data/processed/prs_scores.tsv            sample, pop, super_pop, raw + standardised scores
#           data/processed/prs_summary.tsv           variants per score, counting at each step
#           data/processed/prs_published_check.tsv   reported vs approximated log-OR
#           docs/decisions.md                        block between <!-- prs-summary --> markers
# =============================================================================

suppressPackageStartupMessages({
  library(here)
  library(glue)
  library(tidyverse)
  library(data.table)
  library(bigsnpr)
  library(bigstatsr)
  library(jsonlite)
})

# ---- 0. Settings (D-11, D-16) ---------------------------------------------

P_THRESHOLDS <- c(5e-8, 1e-5, 1e-3, 0.05)
CLUMP_R2     <- 0.1        # C+T clumping r2 (standard PRSice-style default)
CLUMP_KB     <- 250L       # C+T clumping window
EUR_MAF_MIN  <- 0.01       # variants rarer than this in EUR carry no usable GWAS signal
NCORES       <- max(1L, parallel::detectCores() - 1L)

proc_dir <- here("data", "processed")
ss_file  <- here("data", "raw", "iamdgc", "Fritsche-26691988.txt.gz")
api_file <- here("data", "raw", "iamdgc", "gcst003219_associations.json")
for (f in c(ss_file, api_file, file.path(proc_dir, "1kg_chr1_10_qc.rds"))) {
  if (!file.exists(f)) stop(glue("Missing {f}. Run 01_download.R / 03_pca.R first."))
}

# ---- 1. Genotypes ----------------------------------------------------------

obj <- snp_attach(file.path(proc_dir, "1kg_chr1_10_qc.rds"))
G   <- obj$genotypes
samples <- read_tsv(file.path(proc_dir, "samples_qc.tsv"), show_col_types = FALSE)
stopifnot(identical(obj$fam$sample.ID, samples$sample))
ind_eur <- which(samples$super_pop == "EUR")
n <- nrow(G)

# bigsnpr convention (as in its tutorials): map allele1 is the allele counted
# in G, so it is the effect allele `a1` for snp_match.
map <- obj$map |>
  transmute(chr = as.integer(chromosome), pos = physical.pos,
            a1 = allele1, a0 = allele2, id = marker.ID)

# ---- 2. Summary statistics -------------------------------------------------

message("Reading summary statistics...")
ss <- fread(ss_file, showProgress = FALSE,
            select = c("Marker", "Chrom", "Pos", "Allele1", "Allele2",
                       "Ncases", "Ncontrols", "GC.Pvalue", "Overall"),
            colClasses = list(character = c("GC.Pvalue", "Chrom")))
n_ss_all <- nrow(ss)
ss <- ss[Chrom %in% c("1", "10")]
n_ss_chr <- nrow(ss)

# log(p) from the text representation so that p = 6.5e-735 does not underflow.
logp_from_string <- function(p) {
  m <- regmatches(p, regexec("^([0-9.]+)[eE]([-+]?[0-9]+)$", p))
  vapply(seq_along(p), function(i) {
    if (length(m[[i]]) == 3) log(as.numeric(m[[i]][2])) + as.numeric(m[[i]][3]) * log(10)
    else log(as.numeric(p[i]))
  }, numeric(1))
}
ss[, logp := logp_from_string(GC.Pvalue)]
ss[, sign := fcase(Overall == "+", 1, Overall == "-", -1, default = NA_real_)]
n_ss_nodir <- sum(is.na(ss$sign) | is.na(ss$logp))
ss <- ss[!is.na(sign) & !is.na(logp)]
ss[, z := sign * qnorm(logp - log(2), lower.tail = FALSE, log.p = TRUE)]
ss[, lp10 := -logp / log(10)]                          # -log10(p)
n_eff <- with(ss[1], 4 * Ncases * Ncontrols / (Ncases + Ncontrols))

sumstats <- ss |>
  transmute(rsid = Marker, chr = as.integer(Chrom), pos = Pos,
            a1 = Allele1, a0 = Allele2, beta = z, lp10, sign) |>
  as_tibble()

# ---- 3. Match to genotypes (alleles, strand, ambiguous SNPs) ---------------

match_log <- character()
matched <- withCallingHandlers(
  snp_match(sumstats, map, strand_flip = TRUE, join_by_pos = TRUE, remove_dups = TRUE),
  message = function(m) { match_log <<- c(match_log, conditionMessage(m)); invokeRestart("muffleMessage") }
)
message(paste(trimws(match_log), collapse = "\n"))
# snp_match returns beta (here z) re-oriented to the genotype's a1 when the
# alleles were reversed; sign is kept in step with it.
matched <- matched |> mutate(sign = sign(beta)) |> as_tibble()
n_matched <- nrow(matched)

# ---- 4. EUR allele frequency and approximate log-odds ----------------------

cs  <- big_colstats(G, ind.row = ind_eur, ind.col = matched$`_NUM_ID_`, ncores = NCORES)
matched$f_eur <- cs$sum / (2 * length(ind_eur))
matched <- matched |>
  mutate(maf_eur = pmin(f_eur, 1 - f_eur),
         beta_hat = beta * 2 / sqrt(n_eff * 2 * f_eur * (1 - f_eur)),
         ambiguous = FALSE)
n_eur_rare <- sum(matched$maf_eur < EUR_MAF_MIN)
ct_set <- matched |> filter(maf_eur >= EUR_MAF_MIN)

write_tsv(ct_set |> select(rsid, chr, pos, a1, a0, z = beta, lp10, f_eur, beta_hat, col = `_NUM_ID_`),
          file.path(proc_dir, "sumstats_matched.tsv.gz"))

# ---- 5. Score (b): clumping + thresholding ---------------------------------

message(glue("Clumping {format(nrow(ct_set), big.mark = ',')} variants on EUR (r2 < {CLUMP_R2}, {CLUMP_KB} kb)..."))
t0 <- Sys.time()
# snp_clumping wants one statistic per genotype column; columns not in the
# C+T set get 0 and are excluded anyway.
S_full <- rep(0, ncol(G))
S_full[ct_set$`_NUM_ID_`] <- ct_set$lp10
keep_idx <- snp_clumping(G, infos.chr = map$chr, infos.pos = map$pos,
                         ind.row = ind_eur, S = S_full,
                         thr.r2 = CLUMP_R2, size = CLUMP_KB,
                         exclude = setdiff(cols_along(G), ct_set$`_NUM_ID_`),
                         ncores = NCORES)
# snp_clumping returns column indices of G; S was supplied in ct_set order, so
# translate back to ct_set rows.
keep_rows <- match(keep_idx, ct_set$`_NUM_ID_`)
stopifnot(!anyNA(keep_rows))
clumped <- ct_set[keep_rows, ]
message(glue("  kept {format(nrow(clumped), big.mark = ',')} index variants in ",
             "{round(as.numeric(difftime(Sys.time(), t0, units = 'mins')), 1)} min"))

thr_lp <- -log10(P_THRESHOLDS)
ct_scores <- snp_PRS(G, betas.keep = clumped$beta_hat, ind.keep = clumped$`_NUM_ID_`,
                     lpS.keep = clumped$lp10, thr.list = thr_lp)
colnames(ct_scores) <- paste0("prs_ct_", format(P_THRESHOLDS, scientific = TRUE, trim = TRUE))
n_per_thr <- vapply(thr_lp, function(t) sum(clumped$lp10 > t), integer(1))

# ---- 6. Score (a): published independent variants ---------------------------

assoc <- fromJSON(api_file)$`_embedded`$associations
pub <- tibble(
  rsid = vapply(assoc$snp_allele, function(x) x$rs_id[[1]], character(1)),
  or_catalog = as.numeric(assoc$or_per_copy_num),
  p_catalog = assoc$p_value,
  location_grch38 = vapply(assoc$locations, `[[`, character(1), 1),
  genes = vapply(assoc$mapped_genes, paste, character(1), collapse = ",")
) |>
  mutate(chr = as.integer(str_extract(location_grch38, "^[^:]+"))) |>
  filter(chr %in% c(1L, 10L))
n_pub_chr <- nrow(pub)

# Reported OR (Catalog stores it as >= 1 for the risk allele, allele unnamed);
# the original file's direction says whether Allele1 is that risk allele.
pub <- pub |>
  left_join(sumstats |> select(rsid, chr, pos, a1, a0, sign, lp10, beta_z = beta), by = c("rsid", "chr")) |>
  mutate(in_sumstats = !is.na(pos),
         log_or_reported_a1 = sign * abs(log(or_catalog)))

pub_ss <- pub |> filter(in_sumstats) |>
  transmute(rsid, chr, pos, a1, a0, beta = log_or_reported_a1)
pub_matched <- withCallingHandlers(
  snp_match(pub_ss, map, strand_flip = TRUE, join_by_pos = TRUE, match.min.prop = 0),
  message = function(m) invokeRestart("muffleMessage")
) |> as_tibble()
pub_cs <- big_colstats(G, ind.row = ind_eur, ind.col = pub_matched$`_NUM_ID_`)
pub_matched$f_eur <- pub_cs$sum / (2 * length(ind_eur))

pub_check <- pub |>
  left_join(pub_matched |> select(rsid, f_eur, log_or_reported_oriented = beta, col = `_NUM_ID_`), by = "rsid") |>
  mutate(in_genotypes = !is.na(col),
         log_or_reported = abs(log(or_catalog)),
         log_or_approx = abs(beta_z) * 2 / sqrt(n_eff * 2 * f_eur * (1 - f_eur)),
         ratio_approx_over_reported = log_or_approx / log_or_reported) |>
  select(rsid, genes, chr, pos, a1, a0, or_catalog, log_or_reported, direction_a1 = sign,
         p_catalog, lp10, f_eur, log_or_approx, ratio_approx_over_reported, in_sumstats, in_genotypes)
write_tsv(pub_check, file.path(proc_dir, "prs_published_check.tsv"))

pub_score <- big_prodVec(G, pub_matched$beta, ind.col = pub_matched$`_NUM_ID_`)

# ---- 7. Assemble, standardise, save ----------------------------------------

zscore <- function(x) (x - mean(x)) / sd(x)
scores <- samples |>
  mutate(prs_pub_raw = pub_score) |>
  bind_cols(as_tibble(ct_scores) |> rename_with(~ paste0(.x, "_raw"))) |>
  mutate(across(ends_with("_raw"), zscore, .names = "{str_remove(.col, '_raw')}"))
write_tsv(scores, file.path(proc_dir, "prs_scores.tsv"))

weights <- bind_rows(
  pub_matched |> transmute(score = "published", rsid, chr, pos, a1, a0, weight = beta,
                           lp10 = NA_real_, f_eur),
  clumped |> transmute(score = "ct", rsid, chr, pos, a1, a0, weight = beta_hat, lp10, f_eur,
                       min_threshold = P_THRESHOLDS[findInterval(-lp10, -thr_lp) + 1L] |>
                         (\(x) ifelse(lp10 > min(thr_lp), x, NA_real_))())
)
write_tsv(weights, file.path(proc_dir, "prs_weights.tsv"))

summary_tbl <- tribble(
  ~step, ~n,
  "summary-statistics rows, genome-wide",                           n_ss_all,
  "rows on chr1 + chr10",                                            n_ss_chr,
  "removed: missing direction or p-value",                           n_ss_nodir,
  "matched to QC'd genotypes (position + alleles, strand-flipped)",  n_matched,
  glue("removed: EUR MAF < {EUR_MAF_MIN}"),                          n_eur_rare,
  "C+T input variants",                                              nrow(ct_set),
  glue("C+T index variants after clumping (r2 < {CLUMP_R2}, {CLUMP_KB} kb, EUR LD)"), nrow(clumped),
  glue("C+T variants at p < {format(P_THRESHOLDS[1])}"),             n_per_thr[1],
  glue("C+T variants at p < {format(P_THRESHOLDS[2])}"),             n_per_thr[2],
  glue("C+T variants at p < {format(P_THRESHOLDS[3])}"),             n_per_thr[3],
  glue("C+T variants at p < {format(P_THRESHOLDS[4])}"),             n_per_thr[4],
  "published variants (Catalog GCST003219) on chr1 + chr10",         n_pub_chr,
  "published variants present in genotypes (MAF >= 0.01 pooled)",    nrow(pub_matched)
) |> mutate(step = as.character(step))
write_tsv(summary_tbl, file.path(proc_dir, "prs_summary.tsv"))
print(summary_tbl, n = Inf)

# ---- 8. Summary block in docs/decisions.md ---------------------------------

fmt  <- function(x) format(x, big.mark = ",")
fmt2 <- function(x) ifelse(is.na(x), "", formatC(x, format = "f", digits = 2))
fmt3 <- function(x) ifelse(is.na(x), "", formatC(x, format = "f", digits = 3))
block <- c(
  "<!-- prs-summary:start -->",
  glue("_Generated by `R/04_prs.R` on {format(Sys.time(), '%Y-%m-%d %H:%M')}; bigsnpr {packageVersion('bigsnpr')}; ",
       "n_eff = {round(n_eff)}. Do not edit this block by hand._"),
  "",
  "| Step | n |",
  "|---|---:|",
  glue_data(summary_tbl, "| {step} | {fmt(n)} |"),
  "",
  "snp_match report:",
  "",
  paste0("    ", trimws(match_log)),
  "",
  "Published variants on chr1/chr10 (GWAS Catalog curated list), reported vs approximated log-OR:",
  "",
  "| rsID | Genes | Chr:pos (GRCh37) | Allele1/Allele2 | Dir(A1) | OR (Catalog) | log-OR reported | f(A1) EUR | log-OR approx | approx/reported | in genotypes |",
  "|---|---|---|---|---:|---:|---:|---:|---:|---:|---|",
  glue_data(pub_check, "| {rsid} | {genes} | {chr}:{ifelse(is.na(pos), '', pos)} | {ifelse(is.na(a1), '', paste0(a1, '/', a0))} | ",
            "{ifelse(is.na(direction_a1), '', ifelse(direction_a1 > 0, '+', '-'))} | {fmt2(or_catalog)} | {fmt3(log_or_reported)} | ",
            "{fmt3(f_eur)} | {fmt3(log_or_approx)} | {fmt2(ratio_approx_over_reported)} | {ifelse(in_genotypes, 'yes', 'no')} |"),
  "<!-- prs-summary:end -->"
)

decisions <- here("docs", "decisions.md")
lines <- read_lines(decisions)
s <- which(lines == "<!-- prs-summary:start -->"); e <- which(lines == "<!-- prs-summary:end -->")
if (length(s) == 1 && length(e) == 1 && e > s) {
  after <- lines[e + seq_len(length(lines) - e)]
  lines <- c(lines[seq_len(s - 1)], block, after)
} else {
  lines <- c(lines, "", "### PRS outcome (written by `R/04_prs.R`)", "", block)
}
write_lines(lines, decisions)
message(glue("Updated {decisions}"))
