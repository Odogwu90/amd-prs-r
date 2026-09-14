#!/usr/bin/env Rscript
# =============================================================================
# 05_portability.R -- Ancestry portability of the AMD PRS (primary aim)
# -----------------------------------------------------------------------------
# Purpose : Compare the PRS distribution across the five 1000 Genomes
#           super-populations (AFR, AMR, EAS, EUR, SAS): mean and SD of the
#           score per group, standardisation against the EUR reference,
#           overlap with PCs, and the proportion of score variants that are
#           polymorphic / carry the same effect allele frequency in each group.
#           NOTE: 1000 Genomes has no AMD phenotype, so "portability" here
#           means score *distribution shift* and *variant coverage*, not
#           predictive accuracy. The report must state this clearly.
#
# Run     : Rscript R/05_portability.R        (from the repository root)
#
# Inputs  : data/processed/prs_scores.tsv
#           data/processed/prs_weights.tsv
#           data/processed/pca_scores.tsv
#           data/raw/1000G/integrated_call_samples_v3.20130502.ALL.panel
#
# Outputs : data/processed/portability_summary.tsv   per super-population stats
#           data/processed/portability_by_pop.tsv    per 26-population stats
#           figures are drawn by report/report.qmd from these files
#
# Status  : NOT YET IMPLEMENTED.
# =============================================================================

suppressPackageStartupMessages({
  library(here)
  library(glue)
  library(tidyverse)
  library(patchwork)
})

message("05_portability.R: not yet implemented (skeleton only).")
