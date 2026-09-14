# Analysis decisions log

One entry per decision, newest at the bottom. Each entry says what was decided,
why, and who decided it. Scientific choices (QC thresholds, number of PCs,
scoring method) are made by the analyst, not by tooling; entries marked
**PENDING** are waiting for that choice.

Format: `YYYY-MM-DD | D-nn | status | decision` followed by rationale.

---

### 2026-09-14 | D-01 | made | Project scope: chr1 + chr10 only

1000 Genomes Phase 3 chr1 and chr10 are used instead of the whole genome.
Rationale: keeps raw data near 2 GB and QC/PCA feasible on a laptop; the two
chromosomes carry the two largest AMD loci (*CFH* on 1q31, *ARMS2/HTRA1* on
10q26). Consequence: the PRS is a *partial* score and cannot be compared
numerically with published genome-wide AMD scores. The report must say so.

### 2026-09-14 | D-02 | made | Genome build: GRCh37 throughout

1000 Genomes Phase 3 is GRCh37, so the GRCh37 "formatted" IAMDGC file from the
GWAS Catalog is used rather than the GRCh38 "harmonised" file. Avoids lift-over
of either dataset. Trade-off: the GRCh38 file is the only one with a published
md5; the GRCh37 file's md5 is recorded locally on first download.

### 2026-09-14 | D-03 | made | Summary statistics source: GWAS Catalog, not amdgenetics.org

GCST003219 on the GWAS Catalog FTP is an official, stable, form-free copy with
documented terms (EMBL-EBI Terms of Use). The consortium portal could not be
verified over HTTPS at setup time.

### 2026-09-14 | D-04 | made | Tooling: PLINK 2 for VCF import and QC, bigsnpr for PCA and scoring

PLINK 2 reads `.vcf.gz` directly and is fast for filtering. bigsnpr/bigstatsr
work on PLINK 1 `.bed` files and provide `snp_autoSVD` (PCA with LD pruning
and long-range-LD region removal) and clumping tools in R, keeping the
statistics in R as the project convention requires. Binaries live in `bin/`
(gitignored) and are not installed system-wide.

### 2026-09-14 | D-05 | made | Package management: renv with an explicit manifest

`renv.lock` is snapshotted from the `Imports:` field of `DESCRIPTION`
(`renv::settings$snapshot.type("explicit")`) rather than from scanning
`library()` calls. Reason: the scripts are still stubs, so an implicit scan
would produce an incomplete lockfile. Adding a package = add it to
`DESCRIPTION`, `renv::install()`, `renv::snapshot()`.

### 2026-09-14 | D-06 | made | Quarto project rooted at the repository, execute-dir = project

`_quarto.yml` lives at the repo root with `execute-dir: project` so that code
chunks run with the same working directory as `Rscript R/*.R`, which means
`.Rprofile` (renv + here) is honoured during renders.

### 2026-09-14 | D-07 | made | Portability is assessed as distribution shift, not accuracy

1000 Genomes carries no AMD phenotype. "Portability" in this project therefore
means: how the PRS mean/variance, effect-allele frequencies and variant
coverage differ across super-populations, not how well the score predicts
disease outside Europeans. This limitation is stated in the report.

### 2026-09-14 | D-08 | final (analyst) | Variant QC

Applied in `02_qc.R`, in this order, on the pooled 2,504-sample panel:

1. Biallelic SNPs only (`--snps-only just-acgt --max-alleles 2`): indels,
   structural variants and multi-allelic sites are removed.
2. Minor allele frequency >= 0.01 (`--maf 0.01`), computed in the pooled sample.
3. Per-variant missingness < 0.02 (`--geno 0.02`).
4. Strand-ambiguous A/T and C/G SNPs are **not** removed here; they are removed
   at the summary-statistics matching step in `04_prs.R`, where strand
   orientation actually matters.
5. **No Hardy-Weinberg filter** on the pooled sample: the panel mixes five
   super-populations, so departures from HWE are expected from population
   structure alone (Wahlund effect) and would remove informative variants. If
   an HWE filter is ever applied, it is applied within super-population only.

### 2026-09-14 | D-09 | final (analyst) | Sample QC

Use all 2,504 samples of the Phase 3 integrated panel as the starting point.
Remove relatedness with PLINK 2 `--king-cutoff 0.0884` (KING-robust kinship
above the second-degree threshold; one sample of each flagged pair is dropped
by PLINK's greedy algorithm). No sex check: chrX is not in scope, so a sex
check is not possible and the 1000 Genomes panel sex is taken as given. No
per-sample missingness filter is needed because Phase 3 genotypes are
fully called (see the QC outcome below).

### 2026-09-14 | D-10 | final (analyst) | Principal components

Compute PCs with `bigsnpr::snp_autoSVD` on LD-pruned variants (the function
removes long-range LD regions automatically and iterates until no outlier
loadings remain). Keep 10 PCs. Plot PC1-PC4 coloured by super-population.

### 2026-09-14 | D-11 | final (analyst) | PRS construction: two scores side by side

(a) **Published-variant score**: the genome-wide-significant variants reported
by Fritsche et al. 2016 that lie on chr1 and chr10, weighted by the reported
log-odds. (b) **Clumping + thresholding**: `bigsnpr::snp_clumping` followed
by scoring at p < 5e-8, 1e-5, 1e-3 and 0.05. Both scores are standardised
(mean 0, SD 1) within the full 2,504-sample panel, so per-population means are
read as shifts relative to the pooled sample.

### 2026-09-14 | D-12 | final (analyst) | Gene-environment simulation design

Fixed seed. Simulated AMD liability = PRS + smoking + PRS x smoking, where
smoking is binary with 20 % prevalence. The report section is titled
"Illustrative simulation - no real phenotype data" and every output file is
prefixed `SIMULATED_`.

### 2026-09-14 | D-13 | made (assistant, technical; review welcome) | Variant identifiers and QC mechanics

* Variant IDs are rewritten as `chr:pos:ref:alt` (`--set-all-var-ids`) so every
  variant has a unique, self-describing ID; the handful of exact duplicates in
  the 1000 Genomes VCF are then removed with `--rm-dup exclude-all`. rsIDs are
  recovered from the summary-statistics file at the matching step.
* Filters are applied one at a time to separate PLINK 2 output files so that
  the count after each step can be logged (`data/processed/qc_log.csv`).
  Intermediate files are deleted at the end; PLINK logs are kept in
  `data/processed/logs/02_qc/`.
* KING kinship is estimated on chr1 + chr10 only (about 14 % of the autosomal
  genome). This is enough to detect second-degree relatives among 2,504
  samples but the estimates are noisier than a genome-wide KING run.
* **PLINK 2 alpha 6.36 instead of alpha 7.5.** The alpha 7.5 build (10 Sep
  2026) segfaulted about 131k variants into the unfiltered chr1 VCF import
  (reproducible with 1 or 3 threads; the file's md5 is verified). It succeeds
  only if biallelic SNPs are filtered at import, which would hide the baseline
  count. Alpha 6.36 (same date, maintenance line) imports the full file
  cleanly and supports every flag used here, so `bin/plink2.exe` is alpha 6.36.

### QC outcome (written by `R/02_qc.R`)

<!-- qc-summary:start -->
_Not yet run. Execute `Rscript R/02_qc.R` to populate this block._
<!-- qc-summary:end -->
