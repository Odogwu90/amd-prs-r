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

### PENDING | D-08 | analyst | Variant QC thresholds

To decide in `02_qc.R`: minimum MAF, maximum per-variant missingness, HWE
p-value threshold (and whether HWE is tested within populations or overall),
biallelic-SNP restriction, and exclusion of A/T, C/G ambiguous SNPs.

### PENDING | D-09 | analyst | Sample QC

To decide: per-sample missingness threshold; whether to keep all 2,504 or drop
related/duplicate samples (Phase 3 is nominally unrelated); sex checks are not
possible on chr1/chr10 alone.

### PENDING | D-10 | analyst | Number of PCs to retain

To decide after inspecting the scree plot from `03_pca.R`; also whether PCs are
computed on all samples jointly or within super-population.

### PENDING | D-11 | analyst | PRS construction method

Options: (a) the published IAMDGC 52-variant score restricted to chr1/chr10;
(b) clumping + p-value thresholding at one or more thresholds; (c) LDpred2-auto
via bigsnpr. Also: which LD reference (EUR subset of 1000 Genomes is the
natural choice since the GWAS is European).

### PENDING | D-12 | analyst | Simulation design for the G x E illustration

To decide: exposure prevalence, main-effect and interaction effect sizes,
sample size, and seed. All outputs prefixed `SIMULATED_`.
