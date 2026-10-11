# Paths for every R step (mirrors config.py). Source from a script run with Rscript:
#   source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))), "config.R"))
.this <- normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1]))
CODE <- dirname(.this)
ROOT <- dirname(CODE)
REF  <- file.path(ROOT, "data", "original", "reference")
WORK <- file.path(ROOT, "data", "processed", "intermediate")   # contains data/, figures/, logs/, palettes.json
PARAMS <- file.path(CODE, "params.json")
# CellChat 2.2.0.9001 is installed in a separate R library on the analysis machine; override with ACC_CCC_RLIB
CCC_RLIB <- Sys.getenv("ACC_CCC_RLIB", "/Users/joezhao/Software/miniforge3/envs/ccc-rlib")
setwd(WORK)
