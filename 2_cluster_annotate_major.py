#!/usr/bin/env python3
"""Step 2 - cluster the low/high-feature objects, assign major cell types, lineage check, ExNeuron pass 1."""


# ====================================================================================================
# 2a clustering  (Step 3 - cluster the low- and high-feature objects separately (protein-coding genes, GENCODE vM36:)
# ====================================================================================================
import time
import scanpy as sc
from pipeline_lib import *

t0 = time.time(); log = lambda m: print(f"[{time.time() - t0:7.1f}s] {m}", flush=True)
coding = set(CODING_GENES.read_text().split())
for name in ["high_feature", "low_feature"]:
    a = sc.read_h5ad(OUT / f"{name}_cell_adata_combined.h5ad")
    a = a[a.obs["BrainRegion"] == "AC", a.var_names.isin(coding)].copy()
    a = cluster(a, 1.0)
    a.write_h5ad(OUT / f"adata_{name}_cluster_unannotated.h5ad")
    marker_table(a, "leiden", OUT / f"markers_{name}.csv")
    log(f"{name}: {a.obs.leiden.nunique()} clusters")


# ====================================================================================================
# 2b annotation  (Step 4 - major cell-type annotation of the low/high objects (as in reproducilibility.py: one)
# ====================================================================================================
import time
import numpy as np, pandas as pd, scanpy as sc, anndata as ad
from pipeline_lib import *

t0 = time.time(); log = lambda m: print(f"[{time.time() - t0:7.1f}s] {m}", flush=True)

# cluster -> label maps written after reviewing markers_{high,low}_feature.csv (same style and label
# vocabulary as reproducilibility.py: X_ExNeuron = glia/neuron mixed clusters, excluded downstream)
MAPS = {
    "high_feature": {"0": "ExNeuron", "1": "ExNeuron", "2": "ExNeuron", "3": "ExNeuron", "4": "ExNeuron",
                     "5": "ExNeuron", "6": "ExNeuron", "7": "ExNeuron", "8": "Astro_ExNeuron", "9": "Oligo_ExNeuron",
                     "10": "ExNeuron", "11": "InNeuron", "12": "ExNeuron", "13": "InNeuron", "14": "OPC_ExNeuron",
                     "15": "ExNeuron", "16": "InNeuron", "17": "Micro_ExNeuron", "18": "InNeuron", "19": "InNeuron",
                     "20": "OPC_ExNeuron", "21": "Endo", "22": "ExNeuron", "23": "OPC", "24": "Mural_ExNeuron",
                     "25": "InNeuron"},
    "low_feature": {"0": "Micro", "1": "Astro", "2": "Oligo", "3": "Oligo", "4": "ExNeuron", "5": "Astro",
                    "6": "OPC", "7": "Endo", "8": "Mural", "9": "InNeuron", "10": "InNeuron", "11": "ExNeuron",
                    "12": "Oligo", "13": "Mural", "14": "Oligo_Astro", "15": "ExNeuron"},
}
removed = []
annotated = {}
for name in ["high_feature", "low_feature"]:
    a = sc.read_h5ad(OUT / f"adata_{name}_cluster_unannotated.h5ad")
    a.obs["major_celltype"] = a.obs["leiden"].astype(str).map(MAPS[name])
    assert a.obs.major_celltype.notna().all(), "unmapped cluster"
    ls = lineage_scores(a)
    lab = a.obs.major_celltype.values
    own = np.array([ls.iloc[i][t] if t in GLIA else np.nan for i, t in enumerate(lab)])
    flag = np.isin(lab, GLIA) & (own == 0) & (ls["Neuron"].values >= 0.5)
    rem = a.obs.loc[flag, ["SampleID", "leiden", "major_celltype", "n_genes_by_counts"]].assign(object=name)
    removed.append(rem)
    log(f"{name}: lineage check removes {flag.sum()} nuclei: " +
        str(rem.groupby(["major_celltype"]).size().to_dict()))
    a = a[~flag].copy()
    a.obs["major_celltype"] = a.obs["major_celltype"].astype("category")
    a.write_h5ad(OUT / f"adata_{name}_cluster.h5ad")
    annotated[name] = a
rem = pd.concat(removed)
rem.to_csv(OUT / "lineage_check_removed.csv")
pd.crosstab([rem.object, rem.major_celltype], rem.SampleID.astype(str).str[9:]).to_csv(OUT / "lineage_check_removed_summary.csv")

# ExNeuron pass 1 (reproducilibility.py: concat raw of low+high, ExNeuron only, res 0.5, 20 PCs)
comb = ad.concat([annotated["low_feature"].raw.to_adata(), annotated["high_feature"].raw.to_adata()])
comb.write_h5ad(OUT / "_combined_raw_tmp.h5ad")
kw = dict(regress=("total_counts", "total_counts_mt"), normalise=False)
ex = cluster(comb[comb.obs["major_celltype"] == "ExNeuron"].copy(), 0.5, n_pcs=20, **kw)
ex.write_h5ad(OUT / "adata_ExNeuron_cluster.h5ad")
marker_table(ex, "leiden", OUT / "markers_ExNeuron_pass1.csv")
log(f"ExNeuron pass 1: {ex.n_obs} nuclei, {ex.obs.leiden.nunique()} clusters")
