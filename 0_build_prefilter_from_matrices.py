#!/usr/bin/env python3
"""Step 0 - Cell Ranger filtered matrices -> one pre-filter AnnData (QC metrics + Scrublet), as
process_sample() in script/reproducilibility.py. Needs data/original/filtered_mat/<SampleID>/ with
barcodes/features/matrix (see data/original/filtered_matrix_paths.csv). Not run on the analysis Mac:
the matrices are on the cluster; steps 1+ start from the resulting object (config.PREFILTER).
Output: PREFILTER (all_cell_adata_combined_pre_filter.h5ad)."""
import anndata as ad, scanpy as sc
from config import *

RANDOM_STATE = 0
objs = []
for sid, (pipe_name, disp, grp) in SAMPLES.items():
    m = FILTERED_MAT / sid
    m = m / "filtered_feature_bc_matrix" if (m / "filtered_feature_bc_matrix").exists() else m
    a = sc.read_10x_mtx(m, var_names="gene_symbols", cache=False)
    a.var_names_make_unique()
    a.obs["SampleName"] = sid
    genes = a.var_names.str.lower()
    a.var["mt"] = genes.str.startswith("mt-"); a.var["ribo"] = genes.str.startswith(("rps", "rpl"))
    sc.pp.calculate_qc_metrics(a, qc_vars=["mt", "ribo"], percent_top=None, log1p=False, inplace=True)
    sc.external.pp.scrublet(a, random_state=RANDOM_STATE)
    objs.append(a); print(sid, a.n_obs, int(a.obs.predicted_doublet.sum()), "doublets", flush=True)
pre = ad.concat(objs, join="outer", label="SampleID", keys=list(SAMPLES), index_unique="-")
pre.obs["SampleName"] = pre.obs["SampleID"].map({k: v[0] for k, v in SAMPLES.items()})
pre.obs["BrainRegion"] = pre.obs["SampleID"].str.split("_").str[1]
pre.obs["Experiment"] = pre.obs["SampleID"].str.split("_").str[2]
pre.write_h5ad(PREFILTER)
print(pre)
