# .Rprofile -- sourced automatically whenever R starts in this directory
# (interactive R, Rscript from the repo root, and Quarto renders because
# _quarto.yml sets execute-dir: project).

# 1. Activate renv so this project uses its own package library
#    (renv/library, gitignored) pinned by renv.lock.
source("renv/activate.R")

# 2. Anchor here::here() to the repository root. The empty `.here` file marks
#    the root, so paths like here("data", "raw") resolve correctly no matter
#    which sub-directory R was launched from. Guarded so the profile still
#    loads before `here` has been installed.
if (requireNamespace("here", quietly = TRUE)) {
  invisible(here::i_am(".here"))
}
