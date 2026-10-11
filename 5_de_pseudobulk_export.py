#!/usr/bin/env python3
"""Step 5 - cell-level differential expression and pseudobulk export."""


# ====================================================================================================
# 5a cell-level DE  (Step 9 - cell-level differential expression (method of the original notebooks).)
# ====================================================================================================
import sys, json, warnings
from pathlib import Path
import numpy as np, pandas as pd, scanpy as sc
sys.path.insert(0, str(Path(__file__).parent))
from common import DATA
from config import *  # project paths
warnings.simplefilter("ignore")

adata = sc.read_h5ad(H5AD / "ACC_all_cells.h5ad")
del adata.layers["counts"]
CONTRASTS = [("Experiment", "CES", "CON"), ("Experiment", "CPS", "CON"),
             ("Experiment", "CES", "CPS"), ("Experiment2", "SDS", "CON")]


def run(group_key, groups):
    out = []
    for col, cond, ref in CONTRASTS:
        for g in groups:
            a = adata[(adata.obs[group_key] == g) & adata.obs[col].isin([cond, ref])].copy()
            n1, n2 = (a.obs[col] == cond).sum(), (a.obs[col] == ref).sum()
            if min(n1, n2) < 3:
                continue
            a.obs[col] = a.obs[col].astype(str).astype("category")
            sc.tl.rank_genes_groups(a, groupby=col, groups=[cond], reference=ref,
                                    method="t-test", use_raw=False, pts=True)
            d = sc.get.rank_genes_groups_df(a, group=cond)
            d = d.rename(columns={"names": "gene", "logfoldchanges": "log2FC"})
            d = d.drop(columns=[c for c in d.columns if c.startswith("pct_nz")])
            m1 = (a.obs[col] == cond).values
            det = lambda X: pd.Series(np.asarray((X > 0).mean(0)).ravel(), index=a.var_names)
            d["pct_cond"] = det(a.X[m1]).reindex(d["gene"]).values
            d["pct_ref"] = det(a.X[~m1]).reindex(d["gene"]).values
            d.insert(0, "cluster", g); d.insert(0, "reference", ref); d.insert(0, "condition", cond)
            d["n_cells_cond"], d["n_cells_ref"] = n1, n2
            out.append(d)
            print(group_key, cond, ref, g, n1, n2, flush=True)
    return pd.concat(out, ignore_index=True)


major = run("major_celltype", list(adata.obs["major_celltype"].cat.categories))
major.to_csv(DATA / "DE_celllevel_major.csv.gz", index=False)
neu_types = [t for t in adata.obs["subtype"].cat.categories
             if t not in ("Astro", "Oligo", "OPC", "Micro", "Endo", "Mural")]
sub = run("subtype", neu_types)
sub.to_csv(DATA / "DE_celllevel_subtype.csv.gz", index=False)
print("done", major.shape, sub.shape)


# ====================================================================================================
# 5b pseudobulk export  (Step 10 - sum UMI counts per biological sample x cell group (pseudobulk) for limma-voom,)
# ====================================================================================================
import sys
from pathlib import Path
import numpy as np, pandas as pd, scanpy as sc, scipy.sparse as sp
sys.path.insert(0, str(Path(__file__).parent))
from common import DATA
from config import *  # project paths

adata = sc.read_h5ad(H5AD / "ACC_all_cells.h5ad")
for key, name in [("major_celltype", "major"), ("subtype", "subtype")]:
    cols, sums, means, meta = [], [], [], []
    for (g, s), idx in adata.obs.groupby([key, "sample"], observed=True).indices.items():
        cid = f"{g}|||{s}"
        cols.append(cid)
        sums.append(np.asarray(adata.layers["counts"][idx].sum(0)).ravel())
        means.append(np.asarray(adata.X[idx].mean(0)).ravel())
        meta.append(dict(id=cid, group=g, sample=s, condition=s.split("_")[0], n_cells=len(idx)))
    pd.DataFrame(np.vstack(sums).T.astype(np.int64), index=adata.var_names, columns=cols) \
        .to_csv(DATA / f"pseudobulk_counts_{name}.csv.gz")
    pd.DataFrame(np.vstack(means), index=cols, columns=adata.var_names) \
        .to_csv(DATA / f"sample_mean_logexpr_{name}.csv.gz")
    pd.DataFrame(meta).to_csv(DATA / f"pseudobulk_meta_{name}.csv", index=False)
    print(name, len(cols))
