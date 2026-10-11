#!/usr/bin/env python3
"""Step 4 - build ACC_all_cells.h5ad / ACC_neurons.h5ad and the per-type, per-animal annotation check."""


# ====================================================================================================
# 4a objects  (Step 7 - build the two analysis objects every later step reads.)
# ====================================================================================================
import random, time, sys
from pathlib import Path
import numpy as np, pandas as pd, scanpy as sc, anndata as ad
sys.path.insert(0, str(Path(__file__).parent))
from common import DATA, EX_RENAME, IN_RENAME, SAMPLE_RENAME, EXP_RENAME
from config import *  # project paths

N_HVG, N_PCS = 2000, 30
t0 = time.time()
log = lambda m: print(f"[{time.time()-t0:7.1f}s] {m}", flush=True)

adata_l = sc.read_h5ad(PIPE / "adata_low_feature_cluster.h5ad")
adata_h = sc.read_h5ad(PIPE / "adata_high_feature_cluster.h5ad")
comb = ad.concat([adata_l.raw.to_adata(), adata_h.raw.to_adata()])
del adata_l, adata_h
ex = sc.read_h5ad(PIPE / "adata_ExNeuron_cluster2.h5ad")
inn = sc.read_h5ad(PIPE / "adata_InNeuron_cluster.h5ad")
ex_sub = ex.obs["minor_celltype"].astype(str).replace(EX_RENAME)
in_sub = inn.obs["minor_celltype"].astype(str).replace(IN_RENAME)
ex_umap = pd.DataFrame(ex.obsm["X_umap"], index=ex.obs_names)
ex, inn = ex.raw.to_adata(), inn.raw.to_adata()
glia = comb[comb.obs["major_celltype"].isin(["Micro", "Astro", "Oligo", "OPC", "Endo", "Mural"])]
adata = ad.concat([glia, ex, inn])
adata.obs_names_make_unique()
del comb, glia
log(f"combined {adata.shape}")

# ---- metadata -----------------------------------------------------------
adata.obs["Experiment"] = adata.obs["Experiment"].astype(str).map(EXP_RENAME)
adata.obs["Experiment2"] = adata.obs["Experiment"].replace({"CES": "SDS", "CPS": "SDS"})
adata.obs["sample"] = adata.obs["SampleID"].astype(str).map(SAMPLE_RENAME)
assert adata.obs["sample"].notna().all()
sub = pd.concat([ex_sub, in_sub])
adata.obs["subtype"] = sub.reindex(adata.obs_names).fillna(adata.obs["major_celltype"].astype(str)).values
adata.obs["major_celltype"] = adata.obs["major_celltype"].astype(str)
for c in ["Experiment", "Experiment2", "sample", "subtype", "major_celltype"]:
    adata.obs[c] = adata.obs[c].astype("category")

# ---- raw counts from the post-filter object (same barcodes) ---------------
cnt = sc.read_h5ad(PIPE / "all_cell_adata_combined_post_filter.h5ad")
missing = adata.obs_names.difference(cnt.obs_names)
log(f"counts object {cnt.shape}; barcodes missing from counts: {len(missing)}")
assert len(missing) == 0, missing[:5]
cnt = cnt[adata.obs_names, adata.var_names]
adata.layers["counts"] = cnt.X.copy()
del cnt

# ---- all-cell embedding (figure1_major_celltype_umap.py recipe) -----------
def embed(a, n_pcs=N_PCS):
    b = a.copy()
    sc.pp.highly_variable_genes(b, n_top_genes=N_HVG, batch_key="SampleID")
    b = b[:, b.var.highly_variable].copy()
    sc.pp.regress_out(b, ["total_counts", "total_counts_mt"])
    sc.pp.scale(b, max_value=10)
    sc.tl.pca(b, n_comps=50, random_state=7)
    random.seed(7)
    sc.external.pp.harmony_integrate(b, key=["SampleID"], adjusted_basis="X_pca", max_iter_harmony=50)
    sc.pp.neighbors(b, n_neighbors=20, n_pcs=n_pcs, method="umap")
    sc.tl.umap(b, random_state=0)
    return b.obsm["X_umap"], b.obsm["X_pca"][:, :n_pcs]

adata.obsm["X_umap"], adata.obsm["X_pca_harmony"] = embed(adata)
log("all-cell UMAP done")
adata.write_h5ad(H5AD / "ACC_all_cells.h5ad", compression="gzip")
md = adata.obs[["sample", "SampleID", "Experiment", "Experiment2", "major_celltype", "subtype",
                "n_genes_by_counts", "total_counts", "pct_counts_mt"]].copy()
md[["UMAP1", "UMAP2"]] = adata.obsm["X_umap"]
md.to_csv(DATA / "ACC_cell_metadata.csv")

# ---- neuron-only embedding ------------------------------------------------
neu = adata[adata.obs["major_celltype"].isin(["ExNeuron", "InNeuron"])].copy()
neu.obsm["X_umap"], _ = embed(neu)
neu.obsm["X_umap_ExNeuron_rerun"] = ex_umap.reindex(neu.obs_names).values
neu.write_h5ad(H5AD / "ACC_neurons.h5ad", compression="gzip")
log(f"neurons {neu.shape} done")
print(pd.crosstab(adata.obs["subtype"], adata.obs["Experiment"]))


# ====================================================================================================
# 4b annotation check  (Step 8 - annotation sanity check after the CPS_3 fix.)
# ====================================================================================================
import sys
from pathlib import Path
import numpy as np, pandas as pd, anndata as ad
from sklearn.neighbors import NearestNeighbors
sys.path.insert(0, str(Path(__file__).parent))
from common import DATA, ROOT, SAMPLE_ORDER
from config import *  # project paths

MARK = {"ExNeuron": ["Slc17a7", "Slc17a6", "Snap25"], "InNeuron": ["Gad1", "Gad2"], "Astro": ["Slc1a2", "Slc1a3", "Gpc5"],
        "Oligo": ["Mbp", "St18", "Prr5l"], "OPC": ["Vcan", "Pdgfra"], "Micro": ["Hexb", "P2ry12", "Cx3cr1", "Csf1r", "Tmem119", "C1qa", "Inpp5d"],
        "Endo": ["Flt1", "Cldn5", "Mecom"], "Mural": ["Rgs5", "Pdgfrb", "Colec12"]}
a = ad.read_h5ad(H5AD / "ACC_all_cells.h5ad", backed="r")
obs = a.obs[["sample", "major_celltype"]].astype(str).copy()
xy = np.asarray(a.obsm["X_umap"]); lab = obs.major_celltype.values
_, ind = NearestNeighbors(n_neighbors=16).fit(xy).kneighbors(xy)
obs["misplaced"] = (lab[ind[:, 1:]] == lab[:, None]).mean(1) < 0.5
obs["no_lineage_marker"] = False
for t, g in MARK.items():
    m = np.where(lab == t)[0]
    g = [x for x in g if x in a.var_names]
    X = a[m, g].to_memory().X
    obs.iloc[m, obs.columns.get_loc("no_lineage_marker")] = np.asarray((X > 0).sum(1)).ravel() == 0
r = obs.groupby(["major_celltype", "sample"]).agg(n=("misplaced", "size"), misplaced=("misplaced", "sum"),
                                                  no_lineage_marker=("no_lineage_marker", "sum")).reset_index()
r.to_csv(WORK / "logs" / "misplaced_check.csv", index=False)
print(r.pivot(index="major_celltype", columns="sample", values="misplaced").reindex(columns=SAMPLE_ORDER).to_string())
print(r[r.major_celltype == "Micro"].set_index("sample").reindex(SAMPLE_ORDER)[["n", "misplaced", "no_lineage_marker"]].to_string())
