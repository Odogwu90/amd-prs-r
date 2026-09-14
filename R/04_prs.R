#!/usr/bin/env Rscript
# =============================================================================
# 04_prs.R -- Build the AMD polygenic risk score from IAMDGC summary statistics
# -----------------------------------------------------------------------------
# Purpose : Harmonise the IAMDGC 2016 (GCST003219, GRCh37) summary statistics
#           with the QC-passed 1000 Genomes genotypes (allele matching, strand
#           flips, removal of ambiguous A/T, C/G SNPs), select and weight
#           variants according to the method agreed in docs/decisions.md
#           (e.g. clumping + thresholding via bigsnpr::snp_clumping, or the
#           published 52-variant IAMDGC score), and compute a per-sample PRS.
#
# Run     : Rscript R/04_prs.R                (from the repository root)
#
# Inputs  : data/raw/iamdgc/26691988-GCST003219-EFO_0001365-build37.f.tsv.gz
#           data/processed/1kg_chr1_10_qc.{bed,bim,fam}
#
# Outputs : data/processed/sumstats_matched.tsv    variants after harmonisation
#           data/processed/prs_weights.tsv         variant, effect allele, weight
#           data/processed/prs_scores.tsv          sample id, PRS (raw)
#
# Status  : NOT YET IMPLEMENTED. Scoring method and p-value threshold are
#           analyst decisions; see docs/decisions.md.
# =============================================================================

suppressPackageStartupMessages({
  library(here)
  library(glue)
  library(tidyverse)
  library(data.table)
  library(bigsnpr)
  library(bigstatsr)
})

message("04_prs.R: not yet implemented (skeleton only).")
