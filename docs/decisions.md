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
  Even alpha 6.36 segfaults when the VCF import runs with 3 threads, so the
  import step uses `--threads 1` (about 7 min for chr1, 4 min for chr10); all
  PLINK calls are capped at `--memory 1500` because this 8 GB laptop had under
  1 GB free with a browser open and PLINK's default reserves half of RAM.

### 2026-09-14 | D-14 | made (assistant summary; analyst to review) | QC outcome: what was removed and why

Counts are in the auto-generated table below and in `data/processed/qc_log.csv`.

* **Import (s0).** 6,468,094 chr1 and 3,992,219 chr10 records, 2,504 samples,
  exactly as in the VCF headers; nothing is dropped at this stage.
* **Biallelic A/C/G/T SNPs (s1): 440,376 removed (4.2 %).** These are indels,
  multi-allelic sites and structural variants (`<CN0>`, `<INS:ME>` and similar
  symbolic alleles). Summary statistics and PRS tools assume one REF and one
  ALT allele, and indel/SV alleles cannot be matched reliably to a GWAS file.
* **Duplicate IDs (s2): 0 removed.** After rewriting IDs to `chr:pos:ref:alt`,
  no two records were identical on these chromosomes in the v5b call set. The
  step is kept because it is cheap and guarantees unique IDs downstream.
* **MAF >= 0.01 (s3): 8,467,121 removed (84.5 % of the SNPs).** Most variants
  found by sequencing 2,504 people are singletons or very rare; with 5,008
  chromosomes the threshold corresponds to roughly 50 minor-allele copies.
  Rare variants carry almost no weight in a PRS and are the least portable
  across ancestries, so they are excluded. Caveat: the filter is applied to
  the pooled sample, so a variant common in one small population but rare
  overall can be lost; this is worth remembering when interpreting portability.
* **Missingness < 0.02 (s4): 0 removed.** The Phase 3 integrated call set is
  fully called (every sample has a genotype at every site), so the filter is a
  no-op here. It stays in the pipeline so the same script is correct for
  genotype-array data, where it would matter.
* **Concatenation (s5): 1,552,816 variants.** chr1 943,790 + chr10 609,026.
* **Relatedness (s6): 20 samples removed (2,504 -> 2,484).** KING-robust
  kinship > 0.0884 flags second-degree or closer pairs; PLINK drops one member
  of each pair. Removed: 10 AFR (8 of them ASW), 4 EAS, 6 SAS, none in EUR or
  AMR. The 1000 Genomes "unrelated" set is known to contain cryptic
  relatives, and ASW in particular has documented ones, so this is plausible.
  Because kinship was estimated on two chromosomes only, pairs near the
  threshold could be misclassified; the removed IDs are listed in
  `data/processed/samples_removed_relatedness.tsv`. A future check is to
  compare them with the pedigree file
  `integrated_call_samples_v3.20200731.ALL.ped` on the 1000 Genomes FTP.

### QC outcome (written by `R/02_qc.R`)

<!-- qc-summary:start -->
_Generated by `R/02_qc.R` on 2026-09-14 19:58; PLINK v2.0.0-a.6.36 64-bit (10 Sep 2026). Do not edit this block by hand._

| Step | Chr | Filter | Variants | Samples | Variants removed | Samples removed |
|---|---|---|---:|---:|---:|---:|
| s0 | 1 | import VCF (all records) | 6,468,094 | 2504 |  |  |
| s1 | 1 | biallelic A/C/G/T SNPs only | 6,196,151 | 2504 |   271,943 |  0 |
| s2 | 1 | unique chr:pos:ref:alt IDs; exact duplicates removed | 6,196,151 | 2504 |         0 |  0 |
| s3 | 1 | MAF >= 0.01 (pooled sample) |   943,790 | 2504 | 5,252,361 |  0 |
| s4 | 1 | variant missingness < 0.02 |   943,790 | 2504 |         0 |  0 |
| s0 | 10 | import VCF (all records) | 3,992,219 | 2504 |  |  |
| s1 | 10 | biallelic A/C/G/T SNPs only | 3,823,786 | 2504 |   168,433 |  0 |
| s2 | 10 | unique chr:pos:ref:alt IDs; exact duplicates removed | 3,823,786 | 2504 |         0 |  0 |
| s3 | 10 | MAF >= 0.01 (pooled sample) |   609,026 | 2504 | 3,214,760 |  0 |
| s4 | 10 | variant missingness < 0.02 |   609,026 | 2504 |         0 |  0 |
| s5 | 1+10 | concatenate chr1 and chr10 | 1,552,816 | 2504 |  |  |
| s6 | 1+10 | remove relatives, KING kinship > 0.0884 (2nd degree) | 1,552,816 | 2484 |         0 | 20 |

Final fileset: `data/processed/1kg_chr1_10_qc` with 1,552,816 variants and 2484 samples.

Samples kept / removed by super-population:

| Super-population | Kept | Removed (relatedness) |
|---|---:|---:|
| AFR | 651 | 10 |
| AMR | 347 | 0 |
| EAS | 500 | 4 |
| EUR | 503 | 0 |
| SAS | 483 | 6 |

Removed samples: 
`HG00542` (CHS), `HG00595` (CHS), `HG00881` (CDX), `HG02179` (CDX), `HG02479` (ACB), `HG02682` (PJL), `HG03352` (ESN), `HG03754` (STU), `HG03873` (ITU), `HG03899` (STU), `NA19334` (LWK), `NA19913` (ASW), `NA19920` (ASW), `NA20274` (ASW), `NA20318` (ASW), `NA20321` (ASW), `NA20355` (ASW), `NA20362` (ASW), `NA20900` (GIH), `NA21135` (GIH)
<!-- qc-summary:end -->

### 2026-09-14 | D-15 | made (assistant summary; analyst to review) | PCA outcome and mechanics

* **bigsnpr's own rare-variant screen.** `snp_autoSVD` discards variants with
  minor allele count < 10 or MAF < 0.02 before clumping (its defaults), on top
  of the MAF >= 0.01 filter from D-08. 342,421 of 1,552,816 variants were set
  aside for the PCA only; they remain in the QC'd fileset for scoring.
* **Pruning.** Clumping at r2 < 0.2 within 500 kb kept 92,225 variants. The
  outlier-loading test then removed 117 variants in two regions:
  chr1:121.5-142.5 Mb and chr10:42.4-42.5 Mb. Both are pericentromeric and
  both appear on the standard long-range LD exclusion list (Price et al.
  2008), so the automatic detection behaved as expected. Final SVD on
  92,108 variants; converged at iteration 2. Seed 20260914.
* **What the PCs show.** PC1 (7.6 % of variance) separates AFR from everyone
  else; PC2 (2.7 %) separates EAS from EUR, with SAS between them and AMR
  spread along the EUR-AFR and EUR-EAS clines, as expected for admixed
  populations. PC3-PC4 (0.9 %, 0.7 %) mainly resolve AMR internal structure
  and SAS. PCs 5-10 each explain < 0.2 % and carry little population signal;
  they are kept as agreed in D-10 for use as covariates in 05.
* **Variance explained** is computed against the total variance of the
  92,108 scaled variants actually used in the SVD, not the full genome, so
  the percentages are comparable within this project only.
* **Figure encoding.** Five super-populations on one scatter exceed the three
  colours the project palette validates for all-pairs comparison, so the
  figures add point shape and centroid labels as secondary encodings and a
  faceted version with one population per panel. The palette validator
  (Node.js) is not installed on this machine, so the five hues were taken in
  the palette's documented order rather than re-validated locally.

### PCA outcome (written by `R/03_pca.R`)

<!-- pca-summary:start -->
_Generated by `R/03_pca.R` on 2026-09-14 20:38; bigsnpr 1.12.21; snp_autoSVD ran 11.4 min on 3 cores. Do not edit this block by hand._

* Input: 2484 samples x 1,552,816 QC-passed variants.
* Clumping r2 < 0.2 within 500 kb, then outlier-loading removal: **92,108 variants** used in the final SVD.
* Long-range LD regions removed by the outlier test: 2 (chr1:121467098-142536573; chr10:42378137-42527857)

| PC | Singular value | Variance explained (%) | Share of top 10 (%) |
|---:|---:|---:|---:|
| 1 | 4276.62 | 7.62 | 59.77 |
| 2 | 2553.33 | 2.72 | 21.31 |
| 3 | 1450.63 | 0.88 | 6.88 |
| 4 | 1271.82 | 0.67 | 5.29 |
| 5 | 682.01 | 0.19 | 1.52 |
| 6 | 640.80 | 0.17 | 1.34 |
| 7 | 587.22 | 0.14 | 1.13 |
| 8 | 579.96 | 0.14 | 1.10 |
| 9 | 543.32 | 0.12 | 0.96 |
| 10 | 463.83 | 0.09 | 0.70 |

| Super-population | n | PC1 median | PC1 SD | PC2 median | PC2 SD |
|---|---:|---:|---:|---:|---:|
| AFR | 651 | 149.03 | 19.80 | -7.07 | 6.21 |
| AMR | 347 | -44.33 | 22.33 | 27.89 | 27.16 |
| EAS | 500 | -61.94 | 1.20 | -87.22 | 2.80 |
| EUR | 503 | -50.63 | 3.17 | 66.51 | 4.42 |
| SAS | 483 | -46.94 | 1.46 | 15.31 | 8.11 |
<!-- pca-summary:end -->

### 2026-09-14 | D-16 | made (assistant, technical; analyst to review) | How the PRS weights were obtained, and what 04_prs.R found

**Why weights had to be derived.** The public IAMDGC 2016 release (author file
`Fritsche-26691988.txt.gz`) contains, per variant, only a genomic-control
corrected p-value and a direction sign for Allele1; no beta, odds ratio,
standard error or allele frequency. The GWAS Catalog's re-formatted copy
additionally turned every "-" sign into a missing value, so the author
original is used. Reported odds ratios exist only for the paper's 52
independent variants, which the Catalog curates (fetched through its API into
`gcst003219_associations.json`; it does not name the risk allele, so the
direction sign from the author file decides which allele the OR belongs to).

**Score (a), published variants.** Weight for Allele1 = sign x |log OR|.
Ten of the 52 variants lie on chr1/chr10 (eight independent CFH-region
signals, ARMS2/HTRA1, ARHGAP21). Five survive QC; the other five are rare
CFH-region variants (MAF < 1 %, including the R1210C variant rs121913059 with
OR 20) and were removed by the D-08 MAF filter. The score is therefore built
from 5 variants: rs10922109, rs570618, rs61818925 (CFH region), rs3750846
(ARMS2/HTRA1), rs12357257 (ARHGAP21).

**Score (b), clumping + thresholding.** Each weight is an approximate log-odds
recovered from the signed z-score, z = sign x qnorm(p/2, upper tail), computed
on log(p) because p goes below 1e-700, using the standard logistic-regression
relation beta ~ z x 2 / sqrt(n_eff x 2 f (1 - f)), n_eff = 4 Nca Nco / (Nca + Nco)
= 33,892, with f the Allele1 frequency in the 503 EUR samples (the discovery
ancestry). Validation against the five reported ORs: ratio approximated /
reported log-OR between 0.83 and 1.09 (table below), so the weights are on the
right scale and no rescaling is applied. Variants with EUR MAF < 0.01 (175,727)
are excluded because the European GWAS carries no usable signal for them.
Clumping uses the EUR samples as LD reference, r2 < 0.1 within 250 kb, ordered
by -log10 p (PRSice-style defaults). Thresholds p < 5e-8, 1e-5, 1e-3, 0.05
give 102, 139, 293 and 3,883 variants.

**Matching.** 1,528,142 chr1/chr10 rows; 165,323 strand-ambiguous A/T and C/G
SNPs removed (D-08 item 4); 949,179 matched by position and alleles with no
strand flips or allele reversals needed (both sources are GRCh37 forward
strand); the rest are indels or variants outside the QC'd set.

**Caveats for interpretation.** (i) The CFH/CFHR region has long-range LD, so
r2 < 0.1 within 250 kb still leaves many "independent" index variants there;
most of the 102 variants at p < 5e-8 sit in that region and partly tag the same
haplotypes. This is a known weakness of C+T at complex loci and one reason the
published-variant score is reported alongside it. (ii) Genomic-control
corrected p-values make z, and hence the derived weights, slightly
conservative. (iii) Both scores are standardised in the full panel, so a
super-population mean is a shift relative to the pooled sample.

### PRS outcome (written by `R/04_prs.R`)

<!-- prs-summary:start -->
_Generated by `R/04_prs.R` on 2026-09-14 21:30; bigsnpr 1.12.21; n_eff = 33892. Do not edit this block by hand._

| Step | n |
|---|---:|
| summary-statistics rows, genome-wide | 12,023,830 |
| rows on chr1 + chr10 |  1,528,142 |
| removed: missing direction or p-value |          0 |
| matched to QC'd genotypes (position + alleles, strand-flipped) |    949,179 |
| removed: EUR MAF < 0.01 |    175,727 |
| C+T input variants |    773,452 |
| C+T index variants after clumping (r2 < 0.1, 250 kb, EUR LD) |     31,151 |
| C+T variants at p < 5e-08 |        102 |
| C+T variants at p < 1e-05 |        139 |
| C+T variants at p < 0.001 |        293 |
| C+T variants at p < 0.05 |      3,883 |
| published variants (Catalog GCST003219) on chr1 + chr10 |         10 |
| published variants present in genotypes (MAF >= 0.01 pooled) |          5 |

snp_match report:

    1,528,142 variants to be matched.
    165,323 ambiguous SNPs have been removed.
    949,179 variants have been matched; 0 were flipped and 0 were reversed.

Published variants on chr1/chr10 (GWAS Catalog curated list), reported vs approximated log-OR:

| rsID | Genes | Chr:pos (GRCh37) | Allele1/Allele2 | Dir(A1) | OR (Catalog) | log-OR reported | f(A1) EUR | log-OR approx | approx/reported | in genotypes |
|---|---|---|---|---:|---:|---:|---:|---:|---:|---|
| rs3750846 | ARMS2,HTRA1-AS1 | 10:124215565 | C/T | + | 2.81 | 1.033 | 0.195 | 1.126 | 1.09 | yes |
| rs12357257 | ARHGAP21 | 10:24999593 | A/G | + | 1.11 | 0.104 | 0.251 | 0.097 | 0.93 | yes |
| rs191281603 | CFHR5 | 1: |  |  | 2.44 | 0.892 |  |  |  | no |
| rs35292876 | CFH | 1:196706642 | T/C | + | 2.42 | 0.884 |  |  |  | no |
| rs61818925 | CFHR1,CFHR4 | 1:196815450 | G/T | + | 1.67 | 0.511 | 0.579 | 0.426 | 0.83 | yes |
| rs187328863 | KCNT2 | 1:196380158 | T/C | + | 2.27 | 0.820 |  |  |  | no |
| rs148553336 | CFH,KCNT2 | 1:196613173 | C/T | - | 3.45 | 1.238 |  |  |  | no |
| rs121913059 | CFH | 1:196716375 | T/C | + | 20.28 | 3.010 |  |  |  | no |
| rs570618 | CFH | 1:196657064 | G/T | - | 2.38 | 0.867 | 0.640 | 0.833 | 0.96 | yes |
| rs10922109 | CFH | 1:196704632 | A/C | - | 2.63 | 0.968 | 0.425 | 0.827 | 0.85 | yes |
<!-- prs-summary:end -->
