#!/usr/bin/env python3
"""Step 3 - ExNeuron pass 2 and InNeuron subclustering, naming tables, subtype names."""


# ====================================================================================================
# 3a subclustering  (Step 5 - ExNeuron pass 2 and InNeuron subclustering, with naming tables for review.)
# ====================================================================================================
import time
import numpy as np, pandas as pd, scanpy as sc
from pipeline_lib import *

t0 = time.time(); log = lambda m: print(f"[{time.time() - t0:7.1f}s] {m}", flush=True)
PASS1_REMOVE = ["11"]
kw = dict(regress=("total_counts", "total_counts_mt"), normalise=False)
EX_NAME_GENES = ["Foxp2", "Cux2", "Col19a1", "Rorb", "Cdh9", "Trpc6", "Ptprm", "Tshz2", "Nell1", "Ntng1", "Fth1",
                 "Kcnip1", "Htr2c", "Arhgap6", "Slc17a7"]
IN_NAME_GENES = ["Pvalb", "Sst", "Vip", "Lamp5", "Ndnf", "Cnr1", "Gad1", "Gad2"]


def naming_table(a, genes, path):
    r = a.raw
    g = [x for x in genes if x in r.var_names]
    X = r[:, g].X; X = X.toarray() if hasattr(X, "toarray") else np.asarray(X)
    t = pd.DataFrame(X, columns=g).groupby(a.obs["leiden"].astype(str).values).mean().round(2)
    sc.tl.rank_genes_groups(a, groupby="leiden", method="wilcoxon", use_raw=True)
    top = sc.get.rank_genes_groups_df(a, group=None)
    t["top_markers"] = top.groupby("group").head(8).groupby("group")["names"].agg(", ".join)
    t["n"] = a.obs["leiden"].value_counts()
    t = t.loc[sorted(t.index, key=int)]
    t.to_csv(path)
    return a, top


ex = sc.read_h5ad(OUT / "adata_ExNeuron_cluster.h5ad")
ex2 = cluster(ex[~ex.obs["leiden"].astype(str).isin(PASS1_REMOVE)].raw.to_adata().copy(), 0.5, n_pcs=20, **kw)
ex2, top = naming_table(ex2, EX_NAME_GENES, OUT / "naming_ExNeuron_pass2.csv")
top.to_csv(OUT / "ExNeuron_wilcoxon_markers.csv")
ex2.write_h5ad(OUT / "adata_ExNeuron_cluster2_unnamed.h5ad")
log(f"ExNeuron pass 2: {ex2.n_obs} nuclei, {ex2.obs.leiden.nunique()} clusters")

comb = sc.read_h5ad(OUT / "_combined_raw_tmp.h5ad")
inn = cluster(comb[comb.obs["major_celltype"] == "InNeuron"].copy(), 0.3, n_pcs=10, **kw)
inn, top = naming_table(inn, IN_NAME_GENES, OUT / "naming_InNeuron.csv")
top.to_csv(OUT / "InNeuron_wilcoxon_markers.csv")
inn.write_h5ad(OUT / "adata_InNeuron_cluster_unnamed.h5ad")
log(f"InNeuron: {inn.n_obs} nuclei, {inn.obs.leiden.nunique()} clusters")


# ====================================================================================================
# 3b subtype names  (Step 6 - name ExNeuron / InNeuron subclusters after reviewing naming_ExNeuron_pass2.csv /)
# ====================================================================================================
import numpy as np, pandas as pd, scanpy as sc
from pipeline_lib import *

ex = sc.read_h5ad(OUT / "adata_ExNeuron_cluster2_unnamed.h5ad")
ex.obs["leiden_sub"] = ex.obs["leiden"].astype(str).astype("category")
sc.tl.leiden(ex, resolution=0.5, restrict_to=("leiden_sub", ["0"]), key_added="leiden_sub", n_iterations=-1)
g = ex.raw.to_adata()[:, ["Trpc6", "Col19a1"]].X.toarray()
chk = pd.DataFrame(g, columns=["Trpc6", "Col19a1"]).groupby(ex.obs["leiden_sub"].astype(str).values).mean()
chk.round(2).to_csv(OUT / "naming_ExNeuron_cluster0_subclusters.csv")
EX_MAP = {"1": "Foxp2", "2": "Cux2+Rorb", "3": "Cdh9", "4": "Ptprm+Nell1", "5": "Ptprm+Tshz2", "6": "Ntng1",
          "7": "Cux2+Fth1", "8": "Kcnip1", "9": "Cux2+Htr2c", "10": "Arhgap6"}
for k in chk.index[chk.index.str.startswith("0,")]:
    EX_MAP[k] = "Cux2+Trpc6" if chk.loc[k, "Trpc6"] > chk.loc[k, "Col19a1"] else "Cux2+Col19a1"
print("cluster 0 subclusters:", {k: EX_MAP[k] for k in chk.index if k.startswith("0,")})
ex.obs["minor_celltype"] = ex.obs["leiden_sub"].astype(str).map(EX_MAP)
assert ex.obs.minor_celltype.notna().all(), ex.obs.leiden_sub[ex.obs.minor_celltype.isna()].unique()
ex.obs["minor_celltype"] = ex.obs["minor_celltype"].astype("category")
ex.write_h5ad(OUT / "adata_ExNeuron_cluster2.h5ad")

inn = sc.read_h5ad(OUT / "adata_InNeuron_cluster_unnamed.h5ad")
IN_MAP = {"0": "Pvalb", "1": "Vip", "2": "Sst", "3": "Sst", "4": "Lamp5", "5": "Pvalb", "6": "Cnr1"}
inn = inn[inn.obs["leiden"].astype(str) != "7"].copy()
inn.obs["minor_celltype"] = inn.obs["leiden"].astype(str).map(IN_MAP)
assert inn.obs.minor_celltype.notna().all()
inn.obs["minor_celltype"] = inn.obs["minor_celltype"].astype("category")
inn.write_h5ad(OUT / "adata_InNeuron_cluster.h5ad")

for name, a in [("ExNeuron", ex), ("InNeuron", inn)]:
    print(name, a.obs.minor_celltype.value_counts().to_dict())
for f in ["_combined_raw_tmp.h5ad", "adata_ExNeuron_cluster2_unnamed.h5ad", "adata_InNeuron_cluster_unnamed.h5ad"]:
    (OUT / f).unlink(missing_ok=True)
