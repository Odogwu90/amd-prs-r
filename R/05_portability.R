#!/usr/bin/env Rscript
# =============================================================================
# 05_portability.R -- Ancestry portability of the AMD PRS (primary aim, D-07)
# -----------------------------------------------------------------------------
# Purpose : Describe how each score built in 04_prs.R behaves across the five
#           1000 Genomes super-populations (AFR, AMR, EAS, EUR, SAS):
#             1. distribution per group: mean, SD, quantiles, and the mean
#                shift in EUR standard-deviation units (EUR is the discovery
#                ancestry, so it is the reference);
#             2. why the means differ: mean(raw score) = sum_j w_j * 2 f_gj,
#                so between-group mean differences are effect-allele
#                frequency differences; also the expected variance
#                sum_j w_j^2 * 2 f (1 - f) (Hardy-Weinberg, linkage
#                equilibrium) as a variance-ratio prediction;
#             3. variant coverage per group: proportion of score variants
#                polymorphic / with MAF >= 0.01, and agreement of effect-allele
#                frequencies with EUR;
#             4. relation to the PCs: R2 of each score on PC1-PC10 and group
#                means after PC adjustment (per D-10 covariates).
#           NOTE (D-07): 1000 Genomes has no AMD phenotype. Nothing here is
#           predictive accuracy; it is distribution shift and coverage.
#
# Run     : Rscript R/05_portability.R        (from the repository root)
#
# Inputs  : data/processed/prs_scores.tsv, prs_weights.tsv, prs_summary.tsv
#           data/processed/pca_scores.tsv
#           data/processed/1kg_chr1_10_qc.rds (+ .bk)   for per-group allele frequencies
#
# Outputs : data/processed/portability_summary.tsv      per score x super-population
#           data/processed/portability_by_pop.tsv       per score x 26 populations
#           data/processed/portability_coverage.tsv     variant coverage per score x group
#           data/processed/portability_pc.tsv           PC-adjustment results
#           data/processed/portability_contributions.tsv per-variant contribution to each
#                                                       group's mean shift (exact decomposition)
#           data/processed/portability_variant_freqs.tsv effect-allele frequency per group
#           report/figures/05_score_distributions.png
#           report/figures/05_mean_shift.png
#           report/figures/05_allele_freq_vs_eur.png
#           report/figures/05_score_vs_pc1.png
#           docs/decisions.md                           block between <!-- portability-summary --> markers
# =============================================================================

suppressPackageStartupMessages({
  library(here)
  library(glue)
  library(tidyverse)
  library(bigsnpr)
  library(bigstatsr)
  library(patchwork)
})

proc_dir <- here("data", "processed")
fig_dir  <- here("report", "figures")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

super_levels <- c("AFR", "AMR", "EAS", "EUR", "SAS")
pal_super <- c(AFR = "#2a78d6", AMR = "#eb6834", EAS = "#1baf7a", EUR = "#eda100", SAS = "#e87ba4")
shp_super <- c(AFR = 16, AMR = 17, EAS = 15, EUR = 18, SAS = 8)

# ---- 1. Inputs ---------------------------------------------------------------

scores  <- read_tsv(file.path(proc_dir, "prs_scores.tsv"), show_col_types = FALSE) |>
  mutate(super_pop = factor(super_pop, levels = super_levels))
weights <- read_tsv(file.path(proc_dir, "prs_weights.tsv"), show_col_types = FALSE)
summ    <- read_tsv(file.path(proc_dir, "prs_summary.tsv"), show_col_types = FALSE)
pcs     <- read_tsv(file.path(proc_dir, "pca_scores.tsv"), show_col_types = FALSE)
stopifnot(identical(scores$sample, pcs$sample))

obj <- snp_attach(file.path(proc_dir, "1kg_chr1_10_qc.rds"))
G   <- obj$genotypes
stopifnot(identical(obj$fam$sample.ID, scores$sample))
map <- obj$map |> transmute(chr = as.integer(chromosome), pos = physical.pos, a1 = allele1, a0 = allele2,
                            col = row_number())

# Score definitions: column name, label, variant set in `weights`.
thr <- c(5e-8, 1e-5, 1e-3, 0.05)
n_thr <- summ$n[str_detect(summ$step, "^C\\+T variants at p <")]
score_def <- tibble(
  score  = c("prs_pub", paste0("prs_ct_", format(thr, scientific = TRUE, trim = TRUE))),
  label  = c(glue("Published ({sum(weights$score == 'published')} variants)"),
             glue("C+T p<{format(thr, scientific = TRUE, trim = TRUE)} ({n_thr})")),
  set    = c("published", paste0("ct_", format(thr, scientific = TRUE, trim = TRUE))),
  p_thr  = c(NA, thr)
) |> mutate(label = factor(label, levels = label))

# Variant sets with genotype column indices.
variant_sets <- bind_rows(
  weights |> filter(score == "published") |> mutate(set = "published"),
  map_dfr(seq_along(thr), \(i) weights |> filter(score == "ct", lp10 > -log10(thr[i])) |>
            mutate(set = paste0("ct_", format(thr[i], scientific = TRUE, trim = TRUE))))
) |>
  inner_join(map, by = c("chr", "pos", "a1", "a0"))
stopifnot(all(table(variant_sets$set)[score_def$set[-1]] == n_thr))

# ---- 2. Distribution per group ---------------------------------------------

long <- scores |>
  select(sample, pop, super_pop, all_of(score_def$score)) |>
  pivot_longer(all_of(score_def$score), names_to = "score", values_to = "value") |>
  left_join(score_def |> select(score, label), by = "score")

eur_ref <- long |> filter(super_pop == "EUR") |>
  group_by(score) |> summarise(eur_mean = mean(value), eur_sd = sd(value), .groups = "drop")

summary_super <- long |>
  group_by(score, label, super_pop) |>
  summarise(n = n(), mean = mean(value), sd = sd(value), median = median(value),
            q10 = quantile(value, 0.10), q90 = quantile(value, 0.90), .groups = "drop") |>
  left_join(eur_ref, by = "score") |>
  mutate(shift_eur_sd = (mean - eur_mean) / eur_sd,          # mean shift in EUR SD units
         sd_ratio_vs_eur = sd / eur_sd,
         ci_lo = mean - 1.96 * sd / sqrt(n), ci_hi = mean + 1.96 * sd / sqrt(n),
         # share of this group above the EUR 90th percentile
         prop_above_eur_q90 = NA_real_)
eur_q90 <- long |> filter(super_pop == "EUR") |> group_by(score) |>
  summarise(q90_eur = quantile(value, 0.9), .groups = "drop")
summary_super <- summary_super |>
  left_join(long |> left_join(eur_q90, by = "score") |> group_by(score, super_pop) |>
              summarise(p90 = mean(value > q90_eur), .groups = "drop"), by = c("score", "super_pop")) |>
  mutate(prop_above_eur_q90 = p90) |> select(-p90)
write_tsv(summary_super, file.path(proc_dir, "portability_summary.tsv"))

summary_pop <- long |>
  group_by(score, label, super_pop, pop) |>
  summarise(n = n(), mean = mean(value), sd = sd(value), .groups = "drop") |>
  left_join(eur_ref, by = "score") |>
  mutate(shift_eur_sd = (mean - eur_mean) / eur_sd)
write_tsv(summary_pop, file.path(proc_dir, "portability_by_pop.tsv"))

# ---- 3. Allele frequencies, coverage, and the mean decomposition -----------

grp_rows <- split(seq_len(nrow(scores)), scores$super_pop)
freq_by_group <- function(cols) {
  map_dfr(names(grp_rows), \(g) {
    cs <- big_colstats(G, ind.row = grp_rows[[g]], ind.col = cols)
    tibble(col = cols, super_pop = g, f = cs$sum / (2 * length(grp_rows[[g]])))
  })
}
all_cols <- sort(unique(variant_sets$col))
freqs <- freq_by_group(all_cols) |> mutate(super_pop = factor(super_pop, levels = super_levels))

variant_freqs <- variant_sets |>
  select(set, rsid, chr, pos, a1, a0, weight, lp10, col) |>
  inner_join(freqs, by = "col", relationship = "many-to-many")
write_tsv(variant_freqs, file.path(proc_dir, "portability_variant_freqs.tsv"))

# Per set x group: coverage and the HWE/LE expectations.
coverage <- variant_freqs |>
  group_by(set, super_pop) |>
  summarise(n_variants = n(),
            prop_polymorphic = mean(f > 0 & f < 1),
            prop_maf_ge_0.01 = mean(pmin(f, 1 - f) >= 0.01),
            expected_mean_raw = sum(2 * f * weight),
            expected_var_raw  = sum(weight^2 * 2 * f * (1 - f)),
            .groups = "drop")
eur_f <- variant_freqs |> filter(super_pop == "EUR") |> select(set, col, f_eur = f)
freq_agreement <- variant_freqs |>
  inner_join(eur_f, by = c("set", "col")) |>
  group_by(set, super_pop) |>
  summarise(mean_abs_freq_diff_vs_eur = mean(abs(f - f_eur)),
            cor_freq_with_eur = if (n() > 2) cor(f, f_eur) else NA_real_,
            # weighted: |w| * |f - f_eur|, share of total |w|
            weighted_freq_diff = sum(abs(weight) * abs(f - f_eur)) / sum(abs(weight)),
            .groups = "drop")
coverage <- coverage |>
  left_join(freq_agreement, by = c("set", "super_pop")) |>
  left_join(score_def |> select(set, score, label), by = "set") |>
  group_by(set) |>
  # Under linkage equilibrium the score variance would be sum w^2 2f(1-f); the
  # ratio to EUR is a prediction to compare with the observed SD ratio. A gap
  # means LD between score variants (e.g. the CFH region) matters.
  mutate(expected_sd_ratio_le = sqrt(expected_var_raw / expected_var_raw[super_pop == "EUR"])) |>
  ungroup() |>
  select(score, label, set, super_pop, everything())
write_tsv(coverage, file.path(proc_dir, "portability_coverage.tsv"))

# Exact decomposition of the mean shift. mean(raw score in g) = sum_j 2 f_gj w_j
# holds without any assumption, so the shift of group g from EUR, in EUR SD
# units of the raw score, is the sum over variants of
#   c_gj = 2 w_j (f_gj - f_EUR,j) / sd_EUR(raw).
eur_sd_raw <- scores |> filter(super_pop == "EUR") |>
  summarise(across(all_of(paste0(score_def$score, "_raw")), sd)) |>
  pivot_longer(everything(), names_to = "score", values_to = "sd_eur_raw") |>
  mutate(score = str_remove(score, "_raw$"))
contributions <- variant_freqs |>
  inner_join(eur_f, by = c("set", "col")) |>
  left_join(score_def |> select(set, score, label), by = "set") |>
  left_join(eur_sd_raw, by = "score") |>
  mutate(contribution_eur_sd = 2 * weight * (f - f_eur) / sd_eur_raw) |>
  filter(super_pop != "EUR") |>
  select(score, label, set, super_pop, rsid, chr, pos, a1, weight, f, f_eur, contribution_eur_sd) |>
  arrange(label, super_pop, desc(abs(contribution_eur_sd)))
write_tsv(contributions, file.path(proc_dir, "portability_contributions.tsv"))
# Sanity check: the contributions must sum to the observed shift.
chk <- contributions |> group_by(score, super_pop) |> summarise(sum_c = sum(contribution_eur_sd), .groups = "drop") |>
  left_join(summary_super |> select(score, super_pop, shift_eur_sd), by = c("score", "super_pop"))
stopifnot(max(abs(chk$sum_c - chk$shift_eur_sd)) < 1e-6)

# ---- 4. PCs: how much of each score is ancestry ----------------------------

pc_cols <- paste0("PC", 1:10)
pc_res <- map_dfr(score_def$score, \(sc) {
  d <- bind_cols(scores |> select(super_pop, value = all_of(sc)), pcs |> select(all_of(pc_cols)))
  fit <- lm(reformulate(pc_cols, "value"), data = d)
  d$resid <- resid(fit)
  bind_rows(
    tibble(score = sc, stat = "r2_pc1_10", super_pop = NA, value = summary(fit)$r.squared),
    tibble(score = sc, stat = "cor_pc1", super_pop = NA, value = cor(d$value, d$PC1)),
    tibble(score = sc, stat = "cor_pc2", super_pop = NA, value = cor(d$value, d$PC2)),
    d |> group_by(super_pop) |> summarise(value = mean(resid), .groups = "drop") |>
      mutate(score = sc, stat = "mean_resid_after_pc", super_pop = as.character(super_pop)),
    d |> group_by(super_pop) |> summarise(value = sd(resid), .groups = "drop") |>
      mutate(score = sc, stat = "sd_resid_after_pc", super_pop = as.character(super_pop))
  )
}) |> left_join(score_def |> select(score, label), by = "score")
write_tsv(pc_res, file.path(proc_dir, "portability_pc.tsv"))

# ---- 5. Figures --------------------------------------------------------------

theme_set(theme_minimal(base_size = 11) + theme(panel.grid.minor = element_blank()))

# F1 distributions
fig_dist <- ggplot(long, aes(super_pop, value, fill = super_pop)) +
  geom_violin(colour = NA, alpha = 0.55, width = 0.9) +
  geom_boxplot(width = 0.18, outlier.size = 0.4, outlier.alpha = 0.4, fill = "white", colour = "#52514e") +
  geom_hline(yintercept = 0, colour = "#c3c2b7") +
  facet_wrap(~ label, nrow = 1) +
  scale_fill_manual(values = pal_super, guide = "none") +
  labs(x = NULL, y = "Standardised score (full-panel mean 0, SD 1)",
       title = "AMD PRS distributions by 1000 Genomes super-population",
       subtitle = "chr1 + chr10 scores from European-ancestry IAMDGC 2016 statistics; no phenotype is involved") +
  theme(strip.text = element_text(face = "bold"))
ggsave(file.path(fig_dir, "05_score_distributions.png"), fig_dist, width = 13, height = 4.6, dpi = 300, bg = "white")

# F2 mean shift in EUR SD units, with 95% CI, across scores
shift_df <- summary_super |>
  mutate(ci_lo_sd = (ci_lo - eur_mean) / eur_sd, ci_hi_sd = (ci_hi - eur_mean) / eur_sd)
fig_shift <- ggplot(shift_df, aes(label, shift_eur_sd, colour = super_pop, shape = super_pop, group = super_pop)) +
  geom_hline(yintercept = 0, colour = "#c3c2b7") +
  geom_line(linewidth = 0.6, alpha = 0.7) +
  geom_errorbar(aes(ymin = ci_lo_sd, ymax = ci_hi_sd), width = 0.15, linewidth = 0.5) +
  geom_point(size = 2.6) +
  geom_text(data = shift_df |> filter(label == levels(label)[nlevels(label)]),
            aes(label = super_pop), hjust = -0.35, size = 3.2, show.legend = FALSE) +
  scale_colour_manual(values = pal_super, name = "Super-population") +
  scale_shape_manual(values = shp_super, name = "Super-population") +
  scale_x_discrete(expand = expansion(add = c(0.4, 1.0))) +
  labs(x = NULL, y = "Mean shift from EUR (EUR SD units), 95% CI",
       title = "Mean score shift relative to the European reference",
       subtitle = "Zero = same mean as EUR. Shift reflects allele-frequency differences, not risk.") +
  theme(legend.position = "bottom")
ggsave(file.path(fig_dir, "05_mean_shift.png"), fig_shift, width = 9, height = 5.2, dpi = 300, bg = "white")

# F3 effect-allele frequency vs EUR, C+T 5e-8 set, published variants marked
set_5e8 <- score_def$set[2]
f3 <- variant_freqs |> filter(set == set_5e8, super_pop != "EUR") |>
  inner_join(eur_f |> filter(set == set_5e8) |> select(col, f_eur), by = "col") |>
  mutate(published = rsid %in% weights$rsid[weights$score == "published"],
         sign_lab = if_else(weight > 0, "risk-increasing allele", "risk-decreasing allele"))
fig_freq <- ggplot(f3, aes(f_eur, f)) +
  geom_abline(slope = 1, intercept = 0, colour = "#c3c2b7") +
  geom_point(aes(size = abs(weight), colour = super_pop, shape = super_pop), alpha = 0.6) +
  geom_point(data = f3 |> filter(published), shape = 5, size = 3.4, colour = "#0b0b0b", stroke = 0.7) +
  ggrepel::geom_text_repel(data = f3 |> filter(published), aes(label = rsid), size = 2.6,
                           colour = "#0b0b0b", max.overlaps = 20, seed = 1) +
  facet_wrap(~ super_pop, nrow = 1) +
  scale_colour_manual(values = pal_super, guide = "none") +
  scale_shape_manual(values = shp_super, guide = "none") +
  scale_size_continuous(range = c(0.6, 4), name = "|log-odds weight|") +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(x = "Effect-allele frequency in EUR", y = "Effect-allele frequency in group",
       title = glue("Effect-allele frequencies of the {n_thr[1]} C+T (p<5e-8) variants, each group vs EUR"),
       subtitle = "Diamonds: the five published variants. Points off the diagonal drive the mean shifts.") +
  theme(legend.position = "bottom")
ggsave(file.path(fig_dir, "05_allele_freq_vs_eur.png"), fig_freq, width = 12, height = 4.4, dpi = 300, bg = "white")

# F4 score vs PC1 (C+T 5e-8 and published)
f4 <- scores |> select(super_pop, prs_pub, `prs_ct_5e-08`) |> bind_cols(pcs |> select(PC1)) |>
  pivot_longer(c(prs_pub, `prs_ct_5e-08`), names_to = "score", values_to = "value") |>
  left_join(score_def |> select(score, label), by = "score")
fig_pc <- ggplot(f4, aes(PC1, value, colour = super_pop, shape = super_pop)) +
  geom_point(size = 1.1, alpha = 0.55) +
  geom_smooth(aes(PC1, value), inherit.aes = FALSE, method = "lm", formula = y ~ x,
              colour = "#0b0b0b", linewidth = 0.6, se = FALSE) +
  facet_wrap(~ label) +
  scale_colour_manual(values = pal_super, name = "Super-population") +
  scale_shape_manual(values = shp_super, name = "Super-population") +
  labs(x = "PC1 (AFR vs non-AFR axis)", y = "Standardised score",
       title = "Score against the first genotype principal component") +
  theme(legend.position = "bottom")
ggsave(file.path(fig_dir, "05_score_vs_pc1.png"), fig_pc, width = 10, height = 4.8, dpi = 300, bg = "white")
message("Figures written to report/figures/")

# ---- 6. Summary block in docs/decisions.md ---------------------------------

fmt2 <- function(x) ifelse(is.na(x), "", formatC(x, format = "f", digits = 2))
fmt3 <- function(x) ifelse(is.na(x), "", formatC(x, format = "f", digits = 3))
r2_tbl <- pc_res |> filter(stat == "r2_pc1_10") |> select(label, r2 = value)
block <- c(
  "<!-- portability-summary:start -->",
  glue("_Generated by `R/05_portability.R` on {format(Sys.time(), '%Y-%m-%d %H:%M')}. Do not edit this block by hand._"),
  "",
  "Mean shift from EUR in EUR-SD units (95% CI of the mean), SD ratio vs EUR, and share of each group above the EUR 90th percentile:",
  "",
  "| Score | Group | n | Mean | SD | Shift (EUR SD) | 95% CI | SD ratio | > EUR P90 |",
  "|---|---|---:|---:|---:|---:|---|---:|---:|",
  glue_data(shift_df |> arrange(label, super_pop),
            "| {label} | {super_pop} | {n} | {fmt2(mean)} | {fmt2(sd)} | {fmt2(shift_eur_sd)} | ",
            "{fmt2(ci_lo_sd)} to {fmt2(ci_hi_sd)} | {fmt2(sd_ratio_vs_eur)} | {fmt2(prop_above_eur_q90)} |"),
  "",
  "Variant coverage per group. 'SD ratio if LE' is the SD ratio the allele frequencies would give under linkage equilibrium; compare with the observed SD ratio above:",
  "",
  "| Score | Group | Variants | Polymorphic | MAF >= 0.01 | Mean |f - f_EUR| | Weighted |f - f_EUR| | r(f, f_EUR) | SD ratio if LE |",
  "|---|---|---:|---:|---:|---:|---:|---:|---:|",
  glue_data(coverage |> arrange(label, super_pop),
            "| {label} | {super_pop} | {n_variants} | {fmt2(prop_polymorphic)} | {fmt2(prop_maf_ge_0.01)} | ",
            "{fmt3(mean_abs_freq_diff_vs_eur)} | {fmt3(weighted_freq_diff)} | {fmt2(cor_freq_with_eur)} | ",
            "{fmt2(expected_sd_ratio_le)} |"),
  "",
  "Largest single-variant contributions to each group's mean shift (exact decomposition: shift = sum over variants of 2 w (f_group - f_EUR) / SD_EUR), C+T p<5e-8 and published scores:",
  "",
  "| Score | Group | rsID | Chr:pos | A1 | Weight | f group | f EUR | Contribution (EUR SD) |",
  "|---|---|---|---|---|---:|---:|---:|---:|",
  glue_data(contributions |> filter(set %in% c("published", score_def$set[2])) |>
              group_by(score, super_pop) |> slice_head(n = 3) |> ungroup() |> arrange(label, super_pop),
            "| {label} | {super_pop} | {rsid} | {chr}:{pos} | {a1} | {fmt2(weight)} | {fmt2(f)} | {fmt2(f_eur)} | {fmt2(contribution_eur_sd)} |"),
  "",
  "Share of score variance explained by PC1-PC10 (full panel), and group means of the PC-adjusted residual:",
  "",
  "| Score | R2 (PC1-10) | AFR | AMR | EAS | EUR | SAS |",
  "|---|---:|---:|---:|---:|---:|---:|",
  glue_data(pc_res |> filter(stat == "mean_resid_after_pc") |>
              select(label, super_pop, value) |> pivot_wider(names_from = super_pop, values_from = value) |>
              left_join(r2_tbl, by = "label"),
            "| {label} | {fmt2(r2)} | {fmt2(AFR)} | {fmt2(AMR)} | {fmt2(EAS)} | {fmt2(EUR)} | {fmt2(SAS)} |"),
  "<!-- portability-summary:end -->"
)

decisions <- here("docs", "decisions.md")
lines <- read_lines(decisions)
s <- which(lines == "<!-- portability-summary:start -->"); e <- which(lines == "<!-- portability-summary:end -->")
if (length(s) == 1 && length(e) == 1 && e > s) {
  after <- lines[e + seq_len(length(lines) - e)]
  lines <- c(lines[seq_len(s - 1)], block, after)
} else {
  lines <- c(lines, "", "### Portability outcome (written by `R/05_portability.R`)", "", block)
}
write_lines(lines, decisions)
message(glue("Updated {decisions}"))
print(shift_df |> select(label, super_pop, n, mean, sd, shift_eur_sd, sd_ratio_vs_eur), n = Inf)
