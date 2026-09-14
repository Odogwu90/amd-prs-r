# bin/

Local PLINK binaries. This directory is gitignored (except this file) because
the binaries are ~37 MB each and platform-specific.

Download the current stable builds from the official site and place the
executables here:

| Tool | Page | File expected here |
|---|---|---|
| PLINK 1.9 | <https://www.cog-genomics.org/plink/1.9/> | `plink.exe` (Windows) or `plink` |
| PLINK 2 (**alpha 6** line) | <https://www.cog-genomics.org/plink/2.0/> | `plink2.exe` (Windows) or `plink2` |

Versions used during development (2026-09-14): PLINK v1.9.0-rc2 (13 Sep 2026)
and PLINK v2.0.0-a.6.36 (10 Sep 2026), both 64-bit Windows builds. Use the
alpha 6 column of the PLINK 2 download table: the alpha 7.5 build crashed
(segmentation fault) while importing the unfiltered 1000 Genomes chr1 VCF,
whereas alpha 6.36 imports it cleanly (decision D-13). On macOS/Linux,
`chmod +x bin/plink bin/plink2` after downloading.

Check they run:

```sh
bin/plink --version
bin/plink2 --version
```
