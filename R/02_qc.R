#!/usr/bin/env Rscript
# =============================================================================
# 02_qc.R -- Convert VCFs to PLINK format and apply variant / sample QC
# -----------------------------------------------------------------------------
# Purpose : Turn the raw 1000 Genomes chr1 + chr10 VCFs into a single PLINK 2
#           fileset restricted to biallelic autosomal SNPs, then apply the QC
#           filters agreed in docs/decisions.md (MAF, genotype missingness,
#           HWE; sample missingness). Also produce a QC summary for the report.
#
# Run     : Rscript R/02_qc.R                 (from the repository root)
#           Requires bin/plink2.exe (Windows) or bin/plink2 (macOS/Linux);
#           see bin/README.md.
#
# Inputs  : data/raw/1000G/ALL.chr{1,10}.*.vcf.gz
#           data/raw/1000G/integrated_call_samples_v3.20130502.ALL.panel
#
# Outputs : data/processed/1kg_chr1_10_qc.{pgen,pvar,psam}   QC-passed genotypes
#           data/processed/1kg_chr1_10_qc.bed/.bim/.fam      same, PLINK 1 format
#                                                            (needed by bigsnpr)
#           data/processed/qc_summary.tsv                    counts before/after
#
# Status  : NOT YET IMPLEMENTED. QC thresholds are a scientific decision to be
#           made by the analyst; see docs/decisions.md.
# =============================================================================

suppressPackageStartupMessages({
  library(here)
  library(glue)
  library(tidyverse)
  library(data.table)
})

plink2 <- here("bin", if (.Platform$OS.type == "windows") "plink2.exe" else "plink2")
if (!file.exists(plink2)) {
  stop(glue("PLINK 2 binary not found at {plink2}. See bin/README.md."))
}

message("02_qc.R: not yet implemented (skeleton only).")
