"""Paths and global settings for every Python step. All paths are relative to the project root
(the folder that contains data/, code/, figures/, tables/, doc/), so the project can be moved as a whole."""
import os
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CODE = ROOT / "code"
ORIG = ROOT / "data" / "original"                # raw inputs
FILTERED_MAT = ORIG / "filtered_mat"              # Cell Ranger filtered_feature_bc_matrix per sample (step 0)
REF = ORIG / "reference"                          # gene list, GWAS gene sets, SynGO release
PROC = ROOT / "data" / "processed"
H5AD = PROC / "h5ad"                              # final analysis objects (ACC_all_cells.h5ad, ACC_neurons.h5ad)
WORK = PROC / "intermediate"                      # intermediate objects and result tables
PIPE = WORK / "pipeline"                          # steps 1-4: split / clustering / annotation objects
DATA = WORK / "data"                              # steps 4-16: result tables, pseudobulk, CellChat, NMF ...
FIGS = WORK / "figures"                           # figure files as written by the analysis scripts
QC = WORK / "qc"                                  # split-threshold plots
FIG_OUT, TAB_OUT, DOC_OUT = ROOT / "figures", ROOT / "tables", ROOT / "doc"
# Step 0 output = input of step 1. Default: inside data/original; override with ACC_PREFILTER.
PREFILTER = Path(os.environ.get("ACC_PREFILTER", ORIG / "all_cell_adata_combined_pre_filter.h5ad"))
CODING_GENES = REF / "coding_genes_gencode_vM36.txt"
for p in (H5AD, WORK, PIPE, DATA, FIGS, QC, FIG_OUT, TAB_OUT, DOC_OUT, WORK / "logs"):
    p.mkdir(parents=True, exist_ok=True)

SAMPLES = {  # raw SampleID (Cell Ranger folder) -> (pipeline name, display name, group)
    "24090566_AC_Control_5": ("AC_Control_3", "CON_1", "CON"), "24090566_AC_Control_2": ("AC_Control_1", "CON_2", "CON"),
    "24090566_AC_Control_3": ("AC_Control_2", "CON_3", "CON"), "24090566_AC_ES_948": ("AC_ES_3", "CES_1", "CES"),
    "24090566_AC_ES_938": ("AC_ES_1", "CES_2", "CES"), "24090566_AC_ES_940": ("AC_ES_2", "CES_3", "CES"),
    "24090566_AC_PS_966": ("AC_PS_1", "CPS_1", "CPS"), "24090566_AC_PS_968": ("AC_PS_2", "CPS_2", "CPS"),
    "24090566_AC_PS_973": ("AC_PS_3", "CPS_3", "CPS")}
SEED = 0
