"""Shared settings and functions of pipeline steps 1-6 (QC split, clustering, annotation).
Parameters follow script/reproducilibility.py; the nucleus-level lineage check and the KDE valley rule are additions."""
import random
from pathlib import Path
import numpy as np, pandas as pd, scanpy as sc

from config import *
OUT = PIPE                     # all pipeline objects and tables of steps 1-6
MIN_GENES, MAX_GENES, MIN_COUNTS, MAX_COUNTS, MAX_MT_PCT = 200, 4000, 200, 15000, 1.0
N_FEATURE = 2000
SAMPLE_RENAME = {"24090566_AC_Control_5": "AC_Control_3", "24090566_AC_ES_948": "AC_ES_3",
                 "24090566_AC_Control_2": "AC_Control_1", "24090566_AC_ES_938": "AC_ES_1",
                 "24090566_AC_PS_966": "AC_PS_1", "24090566_AC_ES_940": "AC_ES_2",
                 "24090566_AC_Control_3": "AC_Control_2", "24090566_AC_PS_973": "AC_PS_3",
                 "24090566_AC_PS_968": "AC_PS_2"}
# marker panel of reproducilibility.py (cell_markers) ...
CELL_MARKERS = {"Neuron": ["Syn1"], "ExNeuron": ["Slc17a7", "Slc17a6", "Slc17a8"], "InNeuron": ["Gad1", "Gad2"],
                "Oligo": ["St18", "Prr5l"], "Astro": ["Slc1a2", "Gpc5", "Slc1a3"],
                "Epen": ["Cfap299", "Dnah12", "Rgs22", "Wdr49"], "Mural": ["Slc6a20a", "Colec12"],
                "Micro": ["Inpp5d", "Tgfbr1", "Hexb"], "OPC": ["Vcan"], "Endo": ["Slco1a4", "Mecom", "Flt1"]}
# ... plus neuron_journal_gene_list of reproducilibility.py
JOURNAL = {"Neuron": "Snap25", "Astro": "Gja1", "Oligo": "Mog", "Micro": "Selplg", "OPC": "Pdgfra",
           "Epen": "Fam216b", "Endo": "Flt1", "Mural": "Rgs5"}
# NEW: lineage panel used only for the nucleus-level purity check (pipeline markers + journal markers)
LINEAGE = {k: sorted(set(v + ([JOURNAL[k]] if k in JOURNAL else []))) for k, v in CELL_MARKERS.items()}
LINEAGE["Neuron"] = ["Syn1", "Snap25", "Rbfox3", "Syt1"]
GLIA = ["Astro", "Oligo", "OPC", "Micro", "Endo", "Mural", "Epen"]


def cluster(a, res, n_pcs=None, regress=("n_genes_by_counts", "total_counts", "total_counts_mt"), normalise=True):
    """Clustering recipe of reproducilibility.py (normalise/log1p, HVG 2000 by SampleID, regress,
    scale 10, PCA 50 rs 7, Harmony on SampleID, 20 neighbours, UMAP, Leiden)."""
    if normalise:
        sc.pp.normalize_total(a, target_sum=1e4); sc.pp.log1p(a)
    sc.pp.highly_variable_genes(a, n_top_genes=N_FEATURE, batch_key="SampleID")
    a.raw = a
    a = a[:, a.var.highly_variable].copy()
    sc.pp.regress_out(a, list(regress)); sc.pp.scale(a, max_value=10)
    sc.tl.pca(a, n_comps=50, random_state=7)
    random.seed(7)
    sc.external.pp.harmony_integrate(a, key=["SampleID"], adjusted_basis="X_pca", max_iter_harmony=50)
    if n_pcs is None:
        sc.pp.neighbors(a, n_neighbors=20, method="umap")
    else:
        sc.pp.neighbors(a, n_neighbors=20, n_pcs=n_pcs, method="umap")
    sc.tl.umap(a, method="umap")
    sc.tl.leiden(a, resolution=res, n_iterations=-1)
    sc.tl.dendrogram(a, groupby="leiden")
    return a


def lineage_scores(a):
    """Per nucleus: fraction of each lineage's markers detected (from .raw, log-normalised)."""
    r = a.raw
    out = {}
    for k, genes in LINEAGE.items():
        g = [x for x in genes if x in r.var_names]
        X = r[:, g].X
        out[k] = np.asarray((X > 0).sum(1)).ravel() / max(len(g), 1)
    return pd.DataFrame(out, index=a.obs_names)


def marker_table(a, key, path):
    """Dot-plot numbers for the pipeline marker panel per cluster (fraction expressing, mean),
    plus cluster size, dominant-animal share and median genes. Written to CSV for review."""
    r = a.raw
    genes = [g for v in CELL_MARKERS.values() for g in v if g in r.var_names] + \
            [g for g in JOURNAL.values() if g in r.var_names]
    genes = list(dict.fromkeys(genes))
    X = r[:, genes].X
    X = X.toarray() if hasattr(X, "toarray") else np.asarray(X)
    lab = a.obs[key].astype(str).values
    rows = []
    ls = lineage_scores(a)
    for c in pd.unique(lab):
        m = lab == c
        row = {"cluster": c, "n": int(m.sum()),
               "top_sample": a.obs.loc[m, "SampleID"].astype(str).value_counts(normalize=True).index[0][9:],
               "top_sample_frac": round(a.obs.loc[m, "SampleID"].astype(str).value_counts(normalize=True).iloc[0], 2),
               "median_genes": int(np.median(a.obs.loc[m, "n_genes_by_counts"]))}
        for k in LINEAGE:
            row[f"lin_{k}"] = round(float(ls.loc[m, k].mean()), 2)
        for j, g in enumerate(genes):
            row[f"frac_{g}"] = round(float((X[m, j] > 0).mean()), 2)
        rows.append(row)
    t = pd.DataFrame(rows).set_index("cluster")
    t = t.loc[sorted(t.index, key=lambda s: [int(p) if p.isdigit() else p for p in s.replace(",", " ").split()])]
    t.to_csv(path)
    return t


# ---- KDE split threshold -------------------------------------------------------------
from scipy.signal import argrelextrema
from scipy.stats import gaussian_kde


def kde(values, grid_size=1000):          # Gaussian KDE of genes per nucleus (as find_kde_valley in reproducilibility.py)
    values = np.asarray(values, dtype=float); values = values[np.isfinite(values)]
    k = gaussian_kde(values); x = np.linspace(values.min(), values.max(), grid_size); return x, k(x)


def valley_fixed(x, d, min_density_frac=0.05):
    mx, mn = argrelextrema(d, np.greater)[0], argrelextrema(d, np.less)[0]
    mx = mx[d[mx] >= min_density_frac * d.max()]          # ignore tail wiggles
    glia = mx.min(); right = mx[mx > glia][np.argmax(d[mx[mx > glia]])]
    v = mn[(mn > glia) & (mn < right) & (x[mn] < 2 * x[glia])]
    return (float(x[v[np.argmin(d[v])]]) if len(v) else np.nan), float(x[glia]), float(x[right])
