#!/usr/bin/env Rscript
# =============================================================================
# 06_gxe_simulation.R -- SIMULATED gene-environment interaction illustration
# -----------------------------------------------------------------------------
# Purpose : Teaching illustration only (D-12). 1000 Genomes has no phenotype
#           and no exposure data, so BOTH the outcome and the exposure below
#           are simulated. The only real quantity is the PRS from 04_prs.R.
#           The script shows how a PRS x exposure interaction analysis is set
#           up and reported: simulate smoking (binary, 20 % prevalence),
#           simulate AMD status from a logistic model with a PRS x smoking
#           term, then fit models with and without the interaction and with
#           ancestry PCs, and compare.
#
#           Everything this script writes is prefixed SIMULATED_ and the
#           report section is titled
#           "Illustrative simulation - no real phenotype data".
#
# Run     : Rscript R/06_gxe_simulation.R     (from the repository root)
#
# Inputs  : data/processed/prs_scores.tsv      (real PRS; column prs_ct_5e-08)
#           data/processed/pca_scores.tsv      (real PCs, used as covariates)
#
# Outputs : data/processed/SIMULATED_gxe_data.tsv     per-sample simulated exposure/outcome
#           data/processed/SIMULATED_gxe_models.tsv   tidy coefficients (OR, 95% CI) per model
#           data/processed/SIMULATED_gxe_tests.tsv    likelihood-ratio tests between models
#           report/figures/06_SIMULATED_gxe.png       predicted risk vs PRS by smoking status
#           docs/decisions.md                         block between <!-- gxe-summary --> markers
# =============================================================================

suppressPackageStartupMessages({
  library(here)
  library(glue)
  library(tidyverse)
  library(broom)
})

# ---- 0. Simulation design (D-12; effect sizes are illustrative, D-18) ------

SEED        <- 20260914L
PRS_COL     <- "prs_ct_5e-08"   # standardised C+T score at p < 5e-8 (real)
P_SMOKING   <- 0.20             # exposure prevalence
OR_PRS      <- 1.50             # per SD of PRS, in non-smokers
OR_SMOKING  <- 2.00             # smoking vs not, at PRS = 0
OR_INTERACT <- 1.30             # multiplicative interaction: PRS effect per SD is
                                # OR_PRS * OR_INTERACT in smokers
PREVALENCE  <- 0.15             # target overall outcome prevalence (sets the intercept)

proc_dir <- here("data", "processed")
fig_dir  <- here("report", "figures")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

scores <- read_tsv(file.path(proc_dir, "prs_scores.tsv"), show_col_types = FALSE)
pcs    <- read_tsv(file.path(proc_dir, "pca_scores.tsv"), show_col_types = FALSE)
stopifnot(identical(scores$sample, pcs$sample))

# ---- 1. Simulate exposure and outcome --------------------------------------

set.seed(SEED)
sim <- scores |>
  select(sample, pop, super_pop, prs = all_of(PRS_COL)) |>
  bind_cols(pcs |> select(PC1:PC10)) |>
  mutate(smoking_SIM = rbinom(n(), 1, P_SMOKING))

# Linear predictor without intercept; then solve for the intercept that gives
# the target prevalence (mean of the logistic over the sample).
lp0 <- with(sim, log(OR_PRS) * prs + log(OR_SMOKING) * smoking_SIM +
                 log(OR_INTERACT) * prs * smoking_SIM)
intercept <- uniroot(\(a) mean(plogis(a + lp0)) - PREVALENCE, c(-10, 10))$root
sim <- sim |>
  mutate(p_true_SIM = plogis(intercept + lp0),
         amd_SIM    = rbinom(n(), 1, p_true_SIM))
write_tsv(sim, file.path(proc_dir, "SIMULATED_gxe_data.tsv"))
message(glue("Simulated {sum(sim$amd_SIM)} cases among {nrow(sim)} samples ",
             "(prevalence {round(mean(sim$amd_SIM), 3)}); smokers {sum(sim$smoking_SIM)}"))

# ---- 2. Models -------------------------------------------------------------

m_main <- glm(amd_SIM ~ prs + smoking_SIM, family = binomial, data = sim)
m_int  <- glm(amd_SIM ~ prs * smoking_SIM, family = binomial, data = sim)
m_pcs  <- glm(amd_SIM ~ prs * smoking_SIM + PC1 + PC2 + PC3 + PC4 + PC5 +
                PC6 + PC7 + PC8 + PC9 + PC10, family = binomial, data = sim)
# Stratified: PRS effect within smokers and within non-smokers
m_s0 <- glm(amd_SIM ~ prs, family = binomial, data = sim |> filter(smoking_SIM == 0))
m_s1 <- glm(amd_SIM ~ prs, family = binomial, data = sim |> filter(smoking_SIM == 1))

tidy_or <- function(m, label) {
  tidy(m, exponentiate = TRUE, conf.int = TRUE) |>
    mutate(model = label, .before = 1)
}
models <- bind_rows(
  tidy_or(m_main, "1. PRS + smoking"),
  tidy_or(m_int,  "2. PRS x smoking"),
  tidy_or(m_pcs,  "3. PRS x smoking + PC1-10"),
  tidy_or(m_s0,   "4a. Non-smokers only"),
  tidy_or(m_s1,   "4b. Smokers only")
) |>
  mutate(truth = case_when(
    term == "prs" & model %in% c("1. PRS + smoking") ~ NA_real_,   # marginal, no single truth
    term == "prs"                                    ~ OR_PRS,
    term == "smoking_SIM"                            ~ OR_SMOKING,
    term == "prs:smoking_SIM"                        ~ OR_INTERACT,
    term == "prs" & model == "4b. Smokers only"      ~ OR_PRS * OR_INTERACT,
    TRUE ~ NA_real_)) |>
  mutate(truth = if_else(model == "4b. Smokers only" & term == "prs", OR_PRS * OR_INTERACT, truth))
write_tsv(models, file.path(proc_dir, "SIMULATED_gxe_models.tsv"))

lrt <- bind_rows(
  tidy(anova(m_main, m_int, test = "LRT")) |> mutate(comparison = "interaction vs main effects", .before = 1),
  tidy(anova(m_int, m_pcs, test = "LRT"))  |> mutate(comparison = "+ PC1-10 vs interaction", .before = 1)
)
write_tsv(lrt, file.path(proc_dir, "SIMULATED_gxe_tests.tsv"))

# ---- 3. Figure --------------------------------------------------------------

grid <- expand_grid(prs = seq(-3, 3, by = 0.05), smoking_SIM = c(0, 1))
grid <- grid |>
  mutate(fit = predict(m_int, newdata = grid, type = "link", se.fit = TRUE) |> (\(p) p$fit)(),
         se  = predict(m_int, newdata = grid, type = "link", se.fit = TRUE) |> (\(p) p$se.fit)(),
         p_hat = plogis(fit), lo = plogis(fit - 1.96 * se), hi = plogis(fit + 1.96 * se),
         p_true = plogis(intercept + log(OR_PRS) * prs + log(OR_SMOKING) * smoking_SIM +
                           log(OR_INTERACT) * prs * smoking_SIM),
         smoking = factor(if_else(smoking_SIM == 1, "Smoker (simulated)", "Non-smoker (simulated)"),
                          levels = c("Non-smoker (simulated)", "Smoker (simulated)")))
pal_smoke <- c("Non-smoker (simulated)" = "#2a78d6", "Smoker (simulated)" = "#eb6834")

fig <- ggplot(grid, aes(prs, p_hat, colour = smoking, fill = smoking)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15, colour = NA) +
  geom_line(aes(y = p_true), linetype = "dashed", linewidth = 0.6) +
  geom_line(linewidth = 1) +
  geom_rug(data = sim |> mutate(smoking = factor(if_else(smoking_SIM == 1, "Smoker (simulated)", "Non-smoker (simulated)"),
                                                 levels = levels(grid$smoking))),
           aes(prs, y = NULL), sides = "b", alpha = 0.15, length = unit(0.02, "npc")) +
  scale_colour_manual(values = pal_smoke, name = NULL) +
  scale_fill_manual(values = pal_smoke, name = NULL) +
  scale_y_continuous(labels = scales::percent) +
  labs(x = "AMD PRS (C+T p<5e-8, standardised; real)",
       y = "Simulated probability of AMD",
       title = "ILLUSTRATIVE SIMULATION - no real phenotype data",
       subtitle = glue("Outcome simulated with OR {OR_PRS} per SD (PRS), {OR_SMOKING} (smoking), ",
                       "{OR_INTERACT} (interaction)
Solid = fitted model, dashed = generating model, band = 95% CI")) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom",
        plot.title = element_text(face = "bold", colour = "#b3261e"))
ggsave(file.path(fig_dir, "06_SIMULATED_gxe.png"), fig, width = 8.5, height = 5.2, dpi = 300, bg = "white")
message("Figure written.")

# ---- 4. Summary block in docs/decisions.md ---------------------------------

fmt2 <- function(x) ifelse(is.na(x), "", formatC(x, format = "f", digits = 2))
fmtp <- function(x) ifelse(is.na(x), "", format(signif(x, 2)))
main_terms <- models |> filter(term != "(Intercept)", !str_detect(term, "^PC"))
block <- c(
  "<!-- gxe-summary:start -->",
  glue("_Generated by `R/06_gxe_simulation.R` on {format(Sys.time(), '%Y-%m-%d %H:%M')}; seed {SEED}. ",
       "**Everything in this block is simulated except the PRS.** Do not edit by hand._"),
  "",
  glue("* Design: PRS = `{PRS_COL}` (real); smoking ~ Bernoulli({P_SMOKING}); ",
       "logit P(AMD) = a + log({OR_PRS}) PRS + log({OR_SMOKING}) smoking + log({OR_INTERACT}) PRS x smoking; ",
       "a = {round(intercept, 3)} chosen for prevalence {PREVALENCE}."),
  glue("* Realised: {sum(sim$amd_SIM)} simulated cases / {nrow(sim)} ({round(100 * mean(sim$amd_SIM), 1)} %); ",
       "{sum(sim$smoking_SIM)} simulated smokers ({round(100 * mean(sim$smoking_SIM), 1)} %)."),
  "",
  "| Model | Term | OR | 95% CI | p | Generating value |",
  "|---|---|---:|---|---:|---:|",
  glue_data(main_terms, "| {model} | {term} | {fmt2(estimate)} | {fmt2(conf.low)} to {fmt2(conf.high)} | {fmtp(p.value)} | {fmt2(truth)} |"),
  "",
  "| Likelihood-ratio test | df | Deviance change | p |",
  "|---|---:|---:|---:|",
  glue_data(lrt |> filter(!is.na(p.value)), "| {comparison} | {df} | {fmt2(deviance)} | {fmtp(p.value)} |"),
  "<!-- gxe-summary:end -->"
)

decisions <- here("docs", "decisions.md")
lines <- read_lines(decisions)
s <- which(lines == "<!-- gxe-summary:start -->"); e <- which(lines == "<!-- gxe-summary:end -->")
if (length(s) == 1 && length(e) == 1 && e > s) {
  after <- lines[e + seq_len(length(lines) - e)]
  lines <- c(lines[seq_len(s - 1)], block, after)
} else {
  lines <- c(lines, "", "### SIMULATED gene-environment outcome (written by `R/06_gxe_simulation.R`)", "", block)
}
write_lines(lines, decisions)
message(glue("Updated {decisions}"))
print(main_terms |> select(model, term, estimate, conf.low, conf.high, p.value, truth), n = Inf)
