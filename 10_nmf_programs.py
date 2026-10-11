#!/usr/bin/env python3
"""Step 10 - ExNeuron NMF: per-animal factors, robust factors and meta-programs."""


# ====================================================================================================
# 10a per-animal NMF  (NMF gene programs in ExNeuron (step 1/3): per-animal NMF.)
# ====================================================================================================
import sys, time
from pathlib import Path
import numpy as np, pandas as pd, scanpy as sc
from joblib import Parallel, delayed
from sklearn.decomposition import NMF
sys.path.insert(0, str(Path(__file__).parent))
from common import DATA, SEED
from config import *  # project paths

KS = range(4, 10); NTOP = 50; N_HVG = 2000
OUT = DATA / "nmf"; OUT.mkdir(exist_ok=True)
t0 = time.time(); log = lambda m: print(f"[{time.time() - t0:6.0f}s] {m}", flush=True)

ad = sc.read_h5ad(H5AD / "ACC_neurons.h5ad")
ex = ad[ad.obs.major_celltype == "ExNeuron"].copy(); del ad
ex.obs["sample"] = ex.obs["sample"].astype(str)
sc.pp.highly_variable_genes(ex, n_top_genes=N_HVG, batch_key="sample", flavor="seurat")
hvg = ex.var_names[ex.var.highly_variable].tolist()
(OUT / "nmf_hvg.txt").write_text("\n".join(hvg))
log(f"ExNeuron {ex.n_obs} nuclei, {len(hvg)} HVGs")


def prep(X):
    # Gavish et al. 2023: centre each gene within the animal and set negative values to zero, so that
    # factors capture relative (above-average) expression instead of ubiquitous highly expressed genes.
    X = X.toarray() if hasattr(X, "toarray") else np.asarray(X)
    return np.maximum(X - X.mean(0), 0).astype(np.float32)


def fit(sample, X, cells, k):
    nmf = NMF(n_components=k, init="nndsvda", beta_loss="frobenius", solver="cd", max_iter=1000, tol=1e-4, random_state=SEED)
    W = nmf.fit_transform(X); H = nmf.components_
    ids = [f"{sample}_k{k}_f{i + 1}" for i in range(k)]
    usage = pd.DataFrame(W / np.maximum(W.sum(1, keepdims=True), 1e-12), index=cells, columns=ids)
    Hn = H / np.maximum(H.sum(1, keepdims=True), 1e-12)          # loadings normalised per factor
    return pd.DataFrame(Hn.T, index=hvg, columns=ids), usage, nmf.reconstruction_err_


mats = {}
for s_ in sorted(ex.obs["sample"].unique()):
    m = (ex.obs["sample"] == s_).to_numpy()
    mats[s_] = (prep(ex[m, hvg].X), ex.obs_names[m])
jobs = [(s, k) for s in sorted(ex.obs["sample"].unique()) for k in KS]
res = Parallel(n_jobs=8, backend="threading")(delayed(fit)(s, mats[s][0], mats[s][1], k) for s, k in jobs)   # process pools are blocked in this sandbox
log(f"{len(res)} NMF fits done")
H = pd.concat([r[0] for r in res], axis=1)
H.to_csv(OUT / "nmf_factors_genes.csv.gz")
top = pd.DataFrame({c: H[c].sort_values(ascending=False).index[:NTOP].tolist() for c in H.columns})
top.to_csv(OUT / "nmf_factor_top50.csv", index=False)
U = pd.concat([r[1] for r in res], axis=1)
U.to_parquet(OUT / "nmf_cell_usage.parquet")
pd.DataFrame([dict(sample=s, k=k, reconstruction_err=r[2]) for (s, k), r in zip(jobs, res)]).to_csv(OUT / "nmf_fit_errors.csv", index=False)
log(f"{H.shape[1]} factors written")


# ====================================================================================================
# 10b meta-programs  (NMF gene programs in ExNeuron (step 2/3): robust factors -> clustering -> meta-programs (MPs).)
# ====================================================================================================
import sys
from collections import Counter
from pathlib import Path
import numpy as np, pandas as pd, anndata as ad
from scipy.cluster.hierarchy import linkage, fcluster, leaves_list
from scipy.spatial.distance import squareform
sys.path.insert(0, str(Path(__file__).parent))
from common import DATA
from config import *  # project paths

D = DATA / "nmf"; NTOP, NMP, CUT = 50, 30, 0.75
top = pd.read_csv(D / "nmf_factor_top50.csv"); H = pd.read_csv(D / "nmf_factors_genes.csv.gz", index_col=0)
f = np.array(top.columns); S = [set(top[c]) for c in f]; n = len(f)
samp = np.array([c.split("_k")[0] for c in f]); kk = np.array([int(c.split("_k")[1].split("_")[0]) for c in f])
O = np.array([[len(S[i] & S[j]) for j in range(n)] for i in range(n)])
pd.DataFrame(O, index=f, columns=f).to_csv(D / "factor_overlap.csv")

within = np.array([np.any(O[i, (samp == samp[i]) & (kk != kk[i])] >= 35) for i in range(n)])
across = np.array([np.any(O[i, samp != samp[i]] >= 10) for i in range(n)])
keep = within & across
score = np.array([O[i, (samp != samp[i]) & keep].sum() for i in range(n)])
sel = []
for i in np.argsort(-score):
    if keep[i] and not any(samp[j] == samp[i] and O[i, j] >= 10 for j in sel):
        sel.append(i)
sel = np.array(sorted(sel))
Os = O[np.ix_(sel, sel)]; Dm = 1 - Os / NTOP; np.fill_diagonal(Dm, 0)
Z = linkage(squareform(Dm, checks=False), "average"); cl = fcluster(Z, CUT, "distance")

# MPs ordered by number of animals, then size; MP ids 1..N
info = pd.DataFrame(dict(factor=f[sel], sample=samp[sel], k=kk[sel], cluster=cl))
stats = info.groupby("cluster").agg(n_factors=("factor", "size"), n_animals=("sample", "nunique"))
good = stats[stats.n_animals >= 3].sort_values(["n_animals", "n_factors"], ascending=False).index
mpid = {c: f"MP{i + 1}" for i, c in enumerate(good)}
info["MP"] = info.cluster.map(mpid).fillna("unassigned")
info["leaf_order"] = np.argsort(leaves_list(Z))
info.to_csv(D / "robust_factors.csv", index=False)

rows, genes = [], {}
for c in good:
    mem = info.loc[info.cluster == c, "factor"].tolist()
    r = pd.DataFrame({"n_factors_top50": pd.Series(Counter(g for m in mem for g in top[m]))})
    r["mean_loading"] = H.loc[r.index, mem].mean(1)
    r = r.sort_values(["n_factors_top50", "mean_loading"], ascending=False).head(NMP)
    genes[mpid[c]] = r.index.tolist()
    rows += [dict(MP=mpid[c], rank=i + 1, gene=g, n_factors_top50=int(v.n_factors_top50), n_factors=len(mem),
                  mean_loading=v.mean_loading) for i, (g, v) in enumerate(r.iterrows())]
pd.DataFrame(genes).to_csv(D / "metaprogram_genes.csv", index=False)
pd.DataFrame(rows).to_csv(D / "metaprogram_gene_table.csv", index=False)

# GO universes: the HVGs NMF could select (main) and genes detected in >= 1% of ExNeuron (sensitivity)
(D / "go_universe_hvg.txt").write_text((D / "nmf_hvg.txt").read_text())
a = ad.read_h5ad(H5AD / "ACC_neurons.h5ad", backed="r")
m = (a.obs.major_celltype == "ExNeuron").to_numpy()
X = a[m].to_memory().X
det = np.asarray((X > 0).mean(0)).ravel()
(D / "go_universe_expressed.txt").write_text("\n".join(a.var_names[det >= 0.01]))
print(f"factors {n}; recurrent within {within.sum()}, across {across.sum()}; selected {len(sel)}; MPs {len(good)}")
print(stats.loc[good].rename(index=mpid).to_string())
for k, v in genes.items():
    print(k, " ".join(v[:15]))
