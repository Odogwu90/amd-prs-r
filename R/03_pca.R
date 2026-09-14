#!/usr/bin/env Rscript
# =============================================================================
# 03_pca.R -- Principal component analysis for genetic ancestry (D-10)
# -----------------------------------------------------------------------------
# Purpose : Compute 10 PCs on the QC-passed genotypes with bigsnpr::snp_autoSVD,
#           which (i) LD-prunes by clumping (r2 < 0.2 within 500 kb),
#           (ii) computes a randomised partial SVD, and (iii) iteratively
#           removes regions whose loadings are outliers (long-range LD) until
#           none remain. Attach 1000 Genomes labels and draw the report figures.
#
# Run     : Rscript R/03_pca.R                 (from the repository root)
#           Slow steps: reading the .bed into a file-backed matrix (~4 GB on
#           disk, memory-mapped) and clumping 1.55 M variants. Expect 15-40 min
#           on a laptop; uses all cores but one.
#
# Inputs  : data/processed/1kg_chr1_10_qc.{bed,bim,fam}
#           data/processed/samples_qc.tsv                (sample, pop, super_pop)
#
# Outputs : data/processed/1kg_chr1_10_qc.{bk,rds}        bigsnpr backing files
#           data/processed/pca_svd.rds                    full big_SVD object
#           data/processed/pca_scores.tsv                 sample, pop, super_pop, PC1..PC10
#           data/processed/pca_variance.tsv               per PC: singular value, variance share
#           data/processed/pca_variants_used.tsv          variants kept after pruning/outlier removal
#           data/processed/pca_lrld_regions.tsv           long-range LD regions removed
#           report/figures/03_pca_pc1-4.png               PC1 v PC2, PC3 v PC4 by super-population
#           report/figures/03_pca_facets.png              same, one panel per super-population
#           report/figures/03_pca_scree.png               variance explained per PC
#           docs/decisions.md                             block between <!-- pca-summary --> markers
# =============================================================================

suppressPackageStartupMessages({
  library(here)
  library(glue)
  library(tidyverse)
  library(bigsnpr)
  library(bigstatsr)
  library(patchwork)
})

# ---- 0. Settings (D-10) ----------------------------------------------------

K_PCS   <- 10L        # PCs to keep
THR_R2  <- 0.2        # clumping r2 threshold (bigsnpr default)
CLUMP_KB <- 500L      # window = 100 / THR_R2 kb (bigsnpr default expression)
NCORES  <- max(1L, parallel::detectCores() - 1L)
SEED    <- 20260914L  # randomised SVD; fixed for reproducibility

proc_dir <- here("data", "processed")
fig_dir  <- here("report", "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

bed_prefix <- file.path(proc_dir, "1kg_chr1_10_qc")
if (!file.exists(paste0(bed_prefix, ".bed"))) stop("Run R/02_qc.R first.")

super_levels <- c("AFR", "AMR", "EAS", "EUR", "SAS")   # fixed order everywhere

# Categorical colours: the first five slots of the project's reference
# palette, assigned in fixed order to the fixed super-population order.
# Five colours on one scatter exceed the palette's validated three, so shape
# is a secondary encoding and centroid labels give a third; see D-15.
pal_super <- c(AFR = "#2a78d6", AMR = "#eb6834", EAS = "#1baf7a",
               EUR = "#eda100", SAS = "#e87ba4")
shp_super <- c(AFR = 16, AMR = 17, EAS = 15, EUR = 18, SAS = 8)

# ---- 1. Read genotypes into a file-backed matrix ---------------------------
# snp_readBed writes <prefix>.bk (1 byte per genotype, memory-mapped) and
# <prefix>.rds. It is done once; later runs attach the existing files.

rds <- paste0(bed_prefix, ".rds")
if (!file.exists(rds)) {
  message("Reading .bed into a file-backed matrix (one-off, a few minutes)...")
  rds <- snp_readBed(paste0(bed_prefix, ".bed"))
}
obj <- snp_attach(rds)
G   <- obj$genotypes
CHR <- obj$map$chromosome
POS <- obj$map$physical.pos
message(glue("Genotype matrix: {nrow(G)} samples x {format(ncol(G), big.mark = ',')} variants"))

samples <- read_tsv(file.path(proc_dir, "samples_qc.tsv"), show_col_types = FALSE)
stopifnot(identical(obj$fam$sample.ID, samples$sample))   # same order as the .fam

# ---- 2. snp_autoSVD ---------------------------------------------------------
# Steps inside: clumping on MAF-ordered variants -> partial SVD (k = 10) ->
# rolling-window test on loadings -> drop outlier regions -> repeat.

set.seed(SEED)
t0 <- Sys.time()
svd <- snp_autoSVD(G, infos.chr = CHR, infos.pos = POS, k = K_PCS,
                   thr.r2 = THR_R2, size = CLUMP_KB, ncores = NCORES, verbose = TRUE)
elapsed_min <- round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1)
saveRDS(svd, file.path(proc_dir, "pca_svd.rds"))

kept  <- attr(svd, "subset")                # indices of variants used in the final SVD
lrldr <- attr(svd, "lrldr")                 # long-range LD regions removed (may be 0 rows)
message(glue("snp_autoSVD done in {elapsed_min} min: {format(length(kept), big.mark = ',')} variants used, ",
             "{nrow(lrldr)} long-range LD region(s) removed"))

# ---- 3. Variance explained --------------------------------------------------
# svd$d are singular values of the scaled, centred matrix restricted to `kept`.
# Sum of all squared singular values = sum over used columns of (n - 1) * var,
# so each PC's share of the total variance among the variants used is
# d^2 / that sum. big_colstats gives per-column variance of the raw 0/1/2
# genotypes; snp_scaleBinom divides each column by sqrt(2 p (1 - p)).

n <- nrow(G)
cs <- big_colstats(G, ind.col = kept, ncores = NCORES)
p  <- cs$sum / (2 * n)
scaled_var <- cs$var / (2 * p * (1 - p))
total_ss   <- (n - 1) * sum(scaled_var)

variance <- tibble(
  PC = seq_len(K_PCS),
  singular_value = svd$d,
  var_explained_pct = 100 * svd$d^2 / total_ss,
  share_of_top10_pct = 100 * svd$d^2 / sum(svd$d^2)
)
write_tsv(variance, file.path(proc_dir, "pca_variance.tsv"))

# ---- 4. Scores and variant tables ------------------------------------------

scores <- predict(svd) |>
  as_tibble(.name_repair = ~ paste0("PC", seq_len(K_PCS))) |>
  bind_cols(samples, y = _) |>
  mutate(super_pop = factor(super_pop, levels = super_levels))
write_tsv(scores, file.path(proc_dir, "pca_scores.tsv"))

obj$map |>
  slice(kept) |>
  transmute(chromosome, marker.ID, physical.pos, allele1, allele2) |>
  write_tsv(file.path(proc_dir, "pca_variants_used.tsv"))

lrldr_tbl <- if (nrow(lrldr) > 0) as_tibble(lrldr) else tibble(Chr = integer(), Start = integer(), Stop = integer(), Iter = integer())
write_tsv(lrldr_tbl, file.path(proc_dir, "pca_lrld_regions.tsv"))

# ---- 5. Figures -------------------------------------------------------------

centroids <- scores |>
  group_by(super_pop) |>
  summarise(across(PC1:PC4, median), .groups = "drop")

# k may arrive as "1" (from a column name) or 1L; index the table by integer.
axis_lab <- function(k) {
  k <- as.integer(k)
  glue("PC{k} ({round(variance$var_explained_pct[k], 2)} %)")
}

pc_plot <- function(x, y) {
  ggplot(scores, aes({{ x }}, {{ y }}, colour = super_pop, shape = super_pop)) +
    geom_point(size = 1.4, alpha = 0.6, stroke = 0.4) +
    geom_label(data = centroids, aes(label = super_pop), colour = "#0b0b0b",
               fill = "white", alpha = 0.85, size = 3, linewidth = 0, show.legend = FALSE) +
    scale_colour_manual(values = pal_super, name = "Super-population") +
    scale_shape_manual(values = shp_super, name = "Super-population") +
    labs(x = axis_lab(rlang::as_name(enquo(x)) |> str_remove("PC")),
         y = axis_lab(rlang::as_name(enquo(y)) |> str_remove("PC"))) +
    theme_minimal(base_size = 11) +
    theme(panel.grid.minor = element_blank(), legend.position = "bottom")
}

fig_pcs <- (pc_plot(PC1, PC2) | pc_plot(PC3, PC4)) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")
fig_pcs <- fig_pcs + plot_annotation(
  title = "1000 Genomes Phase 3, chr1 + chr10: principal components of genotype",
  subtitle = glue("{n} samples, {format(length(kept), big.mark = ',')} LD-pruned variants; ",
                  "axes show % of genotype variance explained"))
ggsave(file.path(fig_dir, "03_pca_pc1-4.png"), fig_pcs, width = 10, height = 5.5, dpi = 300, bg = "white")

bg <- scores |> select(-super_pop)
fig_facets <- ggplot(scores, aes(PC1, PC2)) +
  geom_point(data = bg, colour = "#c3c2b7", size = 0.8, alpha = 0.5) +
  geom_point(aes(colour = super_pop, shape = super_pop), size = 1.2, alpha = 0.8) +
  facet_wrap(~ super_pop, nrow = 1) +
  scale_colour_manual(values = pal_super, guide = "none") +
  scale_shape_manual(values = shp_super, guide = "none") +
  labs(x = axis_lab(1), y = axis_lab(2),
       title = "Each super-population against the full panel (grey)") +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank())
ggsave(file.path(fig_dir, "03_pca_facets.png"), fig_facets, width = 12, height = 3.2, dpi = 300, bg = "white")

fig_scree <- ggplot(variance, aes(factor(PC), var_explained_pct)) +
  geom_col(fill = "#2a78d6", width = 0.7) +
  geom_text(aes(label = sprintf("%.2f", var_explained_pct)), vjust = -0.4, size = 3, colour = "#52514e") +
  labs(x = "Principal component", y = "Variance explained (%)",
       title = "Scree plot: share of genotype variance per PC") +
  theme_minimal(base_size = 11) +
  theme(panel.grid.major.x = element_blank(), panel.grid.minor = element_blank())
ggsave(file.path(fig_dir, "03_pca_scree.png"), fig_scree, width = 6, height = 4, dpi = 300, bg = "white")

message("Figures written to report/figures/")

# ---- 6. Summary block in docs/decisions.md ---------------------------------

pop_spread <- scores |>
  group_by(super_pop) |>
  summarise(n = n(), PC1_median = median(PC1), PC2_median = median(PC2),
            PC1_sd = sd(PC1), PC2_sd = sd(PC2), .groups = "drop")

fmt2 <- function(x) formatC(x, format = "f", digits = 2)
block <- c(
  "<!-- pca-summary:start -->",
  glue("_Generated by `R/03_pca.R` on {format(Sys.time(), '%Y-%m-%d %H:%M')}; ",
       "bigsnpr {packageVersion('bigsnpr')}; snp_autoSVD ran {elapsed_min} min on {NCORES} cores. ",
       "Do not edit this block by hand._"),
  "",
  glue("* Input: {n} samples x {format(ncol(G), big.mark = ',')} QC-passed variants."),
  glue("* Clumping r2 < {THR_R2} within {CLUMP_KB} kb, then outlier-loading removal: ",
       "**{format(length(kept), big.mark = ',')} variants** used in the final SVD."),
  glue("* Long-range LD regions removed by the outlier test: {nrow(lrldr)}",
       if (nrow(lrldr) > 0) glue(" (", paste(glue_data(lrldr_tbl, "chr{Chr}:{Start}-{Stop}"), collapse = "; "), ")") else ""),
  "",
  "| PC | Singular value | Variance explained (%) | Share of top 10 (%) |",
  "|---:|---:|---:|---:|",
  glue_data(variance, "| {PC} | {fmt2(singular_value)} | {fmt2(var_explained_pct)} | {fmt2(share_of_top10_pct)} |"),
  "",
  "| Super-population | n | PC1 median | PC1 SD | PC2 median | PC2 SD |",
  "|---|---:|---:|---:|---:|---:|",
  glue_data(pop_spread, "| {super_pop} | {n} | {fmt2(PC1_median)} | {fmt2(PC1_sd)} | {fmt2(PC2_median)} | {fmt2(PC2_sd)} |"),
  "<!-- pca-summary:end -->"
)

decisions <- here("docs", "decisions.md")
lines <- read_lines(decisions)
s <- which(lines == "<!-- pca-summary:start -->"); e <- which(lines == "<!-- pca-summary:end -->")
if (length(s) == 1 && length(e) == 1 && e > s) {
  after <- lines[e + seq_len(length(lines) - e)]
  lines <- c(lines[seq_len(s - 1)], block, after)
} else {
  lines <- c(lines, "", "### PCA outcome (written by `R/03_pca.R`)", "", block)
}
write_lines(lines, decisions)
message(glue("Updated {decisions}"))
