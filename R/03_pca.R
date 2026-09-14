#!/usr/bin/env Rscript
# =============================================================================
# 03_pca.R -- Principal component analysis for genetic ancestry
# -----------------------------------------------------------------------------
# Purpose : Compute PCs on LD-pruned, QC-passed genotypes (bigsnpr::snp_autoSVD
#           on the .bed fileset), attach the 1000 Genomes population labels,
#           and produce PC plots coloured by super-population (AFR, AMR, EAS,
#           EUR, SAS). The PCs are used later to (a) confirm the labelled
#           ancestry groups separate as expected and (b) adjust the PRS
#           distributions in 05_portability.R.
#
# Run     : Rscript R/03_pca.R                (from the repository root)
#
# Inputs  : data/processed/1kg_chr1_10_qc.{bed,bim,fam}
#           data/raw/1000G/integrated_call_samples_v3.20130502.ALL.panel
#
# Outputs : data/processed/pca_scores.tsv        sample x PC1..PCk + pop labels
#           data/processed/pca_loadings.rds
#           data/processed/pca_variance.tsv       variance explained per PC
#           report figures are drawn from these files by report/report.qmd
#
# Status  : NOT YET IMPLEMENTED. Number of PCs to keep is an analyst decision;
#           see docs/decisions.md.
# =============================================================================

suppressPackageStartupMessages({
  library(here)
  library(glue)
  library(tidyverse)
  library(bigsnpr)
  library(bigstatsr)
})

message("03_pca.R: not yet implemented (skeleton only).")
