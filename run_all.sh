#!/usr/bin/env bash
# Run every step in order from the project's code/ folder:  bash run_all.sh
# Step 0 (matrices -> pre-filter object) needs data/original/filtered_mat/* and is run separately:
#   $PY 0_build_prefilter_from_matrices.py
# Set ACC_PREFILTER if the pre-filter object is not at data/original/all_cell_adata_combined_pre_filter.h5ad.
set -euo pipefail
cd "$(dirname "$0")"
export PYTHONPATH="$PWD${PYTHONPATH:+:$PYTHONPATH}"   # helper modules: config, common, pipeline_lib, manifest
export PYTHONDONTWRITEBYTECODE=1                       # keep code/ free of __pycache__
PY=${ACC_PY:-/Users/joezhao/Software/miniforge3/envs/sc-pipeline-python/bin/python}        # scanpy 1.11.5, decoupler 2.2.0
RE=${ACC_R_ENRICH:-/Users/joezhao/.claude-science/conda/envs/brain-enrich-r/bin/Rscript}    # R 4.5.3, clusterProfiler, limma, progeny
RC=${ACC_R_CCC:-/Users/joezhao/.claude-science/conda/envs/ccc/bin/Rscript}                  # R 4.5.3 + ccc-rlib CellChat 2.2.0.9001
PD=${ACC_PY_DOC:-/Users/joezhao/.claude-science/conda/envs/python/bin/python}               # any Python with pandas, anndata, openpyxl, python-docx
mkdir -p logs
run() { local log="logs/$1.log"; shift; echo "[$(date +%T)] $*"; "$@" > "$log" 2>&1; }
# ---- nuclei: QC, split, clustering, annotation, analysis objects (steps 1-4) ----
run 1_qc_split               $PY 1_qc_split.py
run 2_cluster_annotate_major $PY 2_cluster_annotate_major.py
run 3_neuron_subtypes        $PY 3_neuron_subtypes.py
run common                   $PY common.py                      # writes palettes.json for the R figures
run 4_build_objects          $PY 4_build_objects.py
# ---- analyses (steps 5-11) ----
run 5_de_pseudobulk_export   $PY 5_de_pseudobulk_export.py
run 6_stats                  $RE 6_stats.R                      # PROGENy network, limma-voom, GO, GSEA
run 7_cellchat_export        $PY 7_cellchat_export.py
run 8_cellchat               $RC 8_cellchat.R                   # coarse and fine groupings
run 9_progeny                $PY 9_progeny.py
OMP_NUM_THREADS=2 run 10_nmf_programs $PY 10_nmf_programs.py
run 11_nmf_go                $RE 11_nmf_go.R
# ---- figures (steps 12-16) ----
run 12_figure1               $PY 12_figure1.py
run 13_figure1EG             $RE 13_figure1EG.R
run 14_figure2               $PY 14_figure2.py
run 15_figure3               $PY 15_figure3.py
run 16_nmf_scores_figureS2   $PY 16_nmf_scores_figureS2.py
# ---- final layout and documents (steps 17-18) ----
run 17_assemble_outputs      $PY 17_assemble_outputs.py
run 18_build_docs            $PD 18_build_docs.py
echo "[$(date +%T)] done"
