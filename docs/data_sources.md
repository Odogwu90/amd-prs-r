# Data sources

All inputs are public. Nothing under `data/` is committed to git; run
`Rscript R/01_download.R` from the repository root to fetch everything and
regenerate the download log at the bottom of this file.

## 1. 1000 Genomes Project, Phase 3 (genotypes and sample panel)

| Item | Value |
|---|---|
| Release | Phase 3, release date 2013-05-02 (directory `release/20130502`) |
| Call set | `phase3_shapeit2_mvncall_integrated_v5b` (v5b re-release; file timestamps on the FTP are 2021-03-16) |
| Genome build | **GRCh37 / hg19** |
| Samples | 2,504 unrelated individuals, 26 populations, 5 super-populations (AFR, AMR, EAS, EUR, SAS) |
| Files used | `ALL.chr1.*.genotypes.vcf.gz` (1.17 GB), `ALL.chr10.*.genotypes.vcf.gz` (0.74 GB), their `.tbi` indexes, and `integrated_call_samples_v3.20130502.ALL.panel` (sample, pop, super_pop, gender) |
| Official location | <https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/> (EBI mirror of the official release) |
| Checksums | md5 values taken from the project file index `https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/current.tree` (index dated 2024-08-16, consulted 2026-09-14). The release directory itself has no md5 files. |
| Documentation | `README_phase3_callset_20150220`, `README_known_issues_20200731`, `README_vcf_info_annotation.20141104` in the same directory |
| Access terms | Open access, no registration. Fort Lauderdale / Toronto principles apply (data producers may publish first on global analyses; not relevant to a learning project). |
| Citation | The 1000 Genomes Project Consortium. A global reference for human genetic variation. *Nature* 526, 68-74 (2015). doi:10.1038/nature15393 |

Why only chr1 and chr10: keeps the download near 2 GB and the QC/PCA steps fast
on a laptop, while still giving a realistic mix of common variation. chr1 also
carries *CFH*, the largest-effect AMD locus, and chr10 carries *ARMS2/HTRA1*,
the second largest, so a score built from these two chromosomes captures a
meaningful share of the known AMD signal. This choice is logged in
`docs/decisions.md`.

## 2. IAMDGC 2016 GWAS summary statistics (advanced AMD)

| Item | Value |
|---|---|
| Study | Fritsche LG *et al.* A large genome-wide association study of age-related macular degeneration highlights contributions of rare and common variants. *Nat Genet* 48, 134-143 (2016). PMID 26691988. doi:10.1038/ng.3448 |
| Consortium | International AMD Genomics Consortium (IAMDGC) |
| Discovery sample | 16,144 advanced AMD cases, 17,832 controls, European ancestry; 12,023,830 variants, 1000 Genomes Phase 1 imputed |
| GWAS Catalog accession | **GCST003219** (<https://www.ebi.ac.uk/gwas/studies/GCST003219>) |
| Official download | GWAS Catalog summary-statistics FTP: <https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST003001-GCST004000/GCST003219/> |
| Files available there | `Fritsche-26691988.txt.gz` (author's original, 147 MB, GRCh37); `harmonised/26691988-GCST003219-EFO_0001365-build37.f.tsv.gz` (Catalog-formatted, standard column names, values unchanged, 157 MB, **GRCh37**); `harmonised/26691988-GCST003219-EFO_0001365.h.tsv.gz` (harmonised and lifted to **GRCh38**, 258 MB, md5 `13f682ec8a40c03d976abe217a1ad232`) |
| File used | `harmonised/...-build37.f.tsv.gz`, because the 1000 Genomes Phase 3 genotypes are GRCh37 and the formatted file has consistent column names without lift-over. No official md5 is published for this file; `01_download.R` records the md5 it computes on first download. |
| Access terms | No registration or form. Downloads are governed by the EMBL-EBI Terms of Use (<https://www.ebi.ac.uk/about/terms-of-use/>): data are freely available, attribution is expected (cite the study and the GWAS Catalog), and the user is responsible for respecting any third-party rights. The GWAS Catalog API lists `terms_of_license` as that page (checked 2026-09-14). |
| Consortium site | The original IAMDGC portal (<http://amdgenetics.org/>) also distributed these statistics. Its certificate could not be verified on 2026-09-14, so its current terms were not checked; the GWAS Catalog copy is used instead. |
| Citation for the Catalog | Cerezo M *et al.* The NHGRI-EBI GWAS Catalog: standards for reusability, sustainability and diversity. *Nucleic Acids Res* 53, D998-D1005 (2025). |

Column layout of the formatted file (per the Catalog `readme.txt`; confirm on
first inspection in `04_prs.R`): `variant_id`, `p_value`, `chromosome`,
`base_pair_location`, `odds_ratio` or `beta`, `standard_error`,
`effect_allele`, `other_allele`, plus any author-specific columns.

## 3. Software binaries (not data, but downloaded)

| Tool | Version installed | Source |
|---|---|---|
| PLINK 1.9 | v1.9.0-rc2 64-bit (13 Sep 2026) | <https://www.cog-genomics.org/plink/1.9/> -> `plink_win64_20260913.zip` |
| PLINK 2 | v2.0.0-a.7.5 64-bit (10 Sep 2026) | <https://www.cog-genomics.org/plink/2.0/> -> `alpha7/plink2_win64_20260910.zip` |

See `bin/README.md` for how to re-download them on a fresh clone.

## Download log

<!-- download-log:start -->
_Generated by `R/01_download.R` on 2026-09-14 12:18; R 4.6.1. Do not edit this block by hand._

| Dataset | File | Size (MB) | md5 | md5 source | Status | Date |
|---|---|---:|---|---|---|---|
| 1000G | `ALL.chr1.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz` | 1165.0 | `76f1d3fef27c6c3f451cdfc515250a0e` | official | present_verified | 2026-09-14 |
| 1000G | `ALL.chr1.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz.tbi` |    0.2 | `8c40851a4e9d3c6c9a09461961c6e4c4` | official | present_verified | 2026-09-14 |
| 1000G | `ALL.chr10.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz` |  741.7 | `45814b06e651f8fb409364aa6d65257e` | official | present_verified | 2026-09-14 |
| 1000G | `ALL.chr10.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz.tbi` |    0.1 | `266a7200332b971c741b1520e03bda80` | official | present_verified | 2026-09-14 |
| 1000G | `integrated_call_samples_v3.20130502.ALL.panel` |    0.1 | `7ee5675553088230530a7fe88c22f201` | official | present_verified | 2026-09-14 |
| IAMDGC | `26691988-GCST003219-EFO_0001365-build37.f.tsv.gz` |  156.6 | `fd6d8db958f96fb63770824f22920ea6` | recorded on first download | present_verified | 2026-09-14 |
<!-- download-log:end -->
