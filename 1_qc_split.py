#!/usr/bin/env python3
"""Step 1 - nucleus QC, per-animal low/high-feature split, and the split plots (Supplementary Fig. S1)."""


# ====================================================================================================
# 1a QC and split  (Step 1 - nucleus QC and low/high-feature split.)
# ====================================================================================================
import time
import numpy as np, pandas as pd, scanpy as sc, anndata as ad
from pipeline_lib import *

t0 = time.time(); log = lambda m: print(f"[{time.time() - t0:7.1f}s] {m}", flush=True)
pre = sc.read_h5ad(PREFILTER)
pre.obs["SampleID"] = pre.obs["SampleID"].astype(str)
o = pre.obs
kd = []
for s_ in SAMPLES:
    x, d = kde(o.loc[o.SampleID == s_, "n_genes_by_counts"])
    vf, gp, rp = valley_fixed(x, d)
    kd.append(dict(SampleName=s_, valley=vf, glial_peak=gp, neuronal_peak=rp))
kd = pd.DataFrame(kd).set_index("SampleName")
used = kd.valley.fillna(kd.valley.median())
kd["threshold_used"] = used
kd["threshold_rule"] = np.where(kd.valley.notna(), "KDE valley", "median of other animals")
log("\n" + kd.round(1).to_string())
singlet = ~o["predicted_doublet"].fillna(False).astype(bool)
common = ((o.total_counts >= MIN_COUNTS) & (o.total_counts <= MAX_COUNTS) & (o.n_genes_by_counts >= MIN_GENES)
          & (o.n_genes_by_counts <= MAX_GENES) & (o.pct_counts_mt <= MAX_MT_PCT) & singlet)
v = o.SampleID.map(used)
high_m, low_m = common & (o.n_genes_by_counts > v), common & (o.n_genes_by_counts <= v)
keep_genes = ~pre.var_names.str.lower().str.startswith("mt-")
rec = [dict(SampleName=s, glial_peak=kd.loc[s, "glial_peak"], neuronal_peak=kd.loc[s, "neuronal_peak"], threshold_used=used[s],
            threshold_rule=kd.loc[s, "threshold_rule"], n_before=int((o.SampleID == s).sum()), n_high=int((high_m & (o.SampleID == s)).sum()),
            n_low=int((low_m & (o.SampleID == s)).sum()), n_after=int((common & (o.SampleID == s)).sum()),
            n_predicted_doublets=int((~singlet & (o.SampleID == s)).sum())) for s in SAMPLES]
pd.DataFrame(rec).to_csv(OUT / "sample_qc_thresholds.csv", index=False)


def finish(a):
    a = a[:, keep_genes].copy()
    a.obs["BrainRegion"] = a.obs["SampleID"].str.split("_").str[1]
    a.obs["Experiment"] = a.obs["SampleID"].str.split("_").str[2]
    a.obs["SampleName"] = a.obs["SampleID"].map(SAMPLE_RENAME)
    for c in ["SampleID", "SampleName", "BrainRegion", "Experiment"]:
        a.obs[c] = a.obs[c].astype("category")
    return a


adata_high, adata_low = finish(pre[high_m.values]), finish(pre[low_m.values])
del pre
adata_high.write_h5ad(OUT / "high_feature_cell_adata_combined.h5ad")
adata_low.write_h5ad(OUT / "low_feature_cell_adata_combined.h5ad")
ad.concat([adata_high, adata_low], join="outer", merge="same").write_h5ad(OUT / "all_cell_adata_combined_post_filter.h5ad")
log(f"high {adata_high.shape}, low {adata_low.shape}")


# ====================================================================================================
# 1b split plots  (Step 2 - per-animal KDE plots of genes per nucleus with the split threshold of step 1 (Supplementary Fig. S2).)
# ====================================================================================================
import numpy as np, pandas as pd, anndata as ad
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy.signal import argrelextrema
from pipeline_lib import *

plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Arial", "DejaVu Sans"], "font.size": 7,
                     "pdf.fonttype": 42, "axes.spines.top": False, "axes.spines.right": False, "axes.linewidth": 0.6})
DISPLAY = {k: v[1] for k, v in SAMPLES.items()}
thr = pd.read_csv(OUT / "sample_qc_thresholds.csv").set_index("SampleName")
obs = ad.read_h5ad(PREFILTER, backed="r").obs[["SampleID", "n_genes_by_counts"]]
obs["SampleID"] = obs["SampleID"].astype(str)


def panel(ax, s, small=False):
    g = obs.loc[obs.SampleID == s, "n_genes_by_counts"].to_numpy(float)
    x, d = kde(g); v_new, gp, rp = valley_fixed(x, d)
    used = thr.loc[s, "threshold_used"]
    assert np.isfinite(used) and (np.isnan(v_new) or np.isclose(v_new, used)), s
    peaks = argrelextrema(d, np.greater)[0]; peaks = peaks[d[peaks] >= 0.05 * d.max()]
    ax.fill_between(x, d, color="0.88", lw=0); ax.plot(x, d, color="0.2", lw=0.8)
    ax.plot(x[peaks], d[peaks], "v", color="0.35", ms=3 if small else 4)
    ax.axvline(used, color="#C0392B", lw=1.1, label=f"threshold {used:,.0f}" + ("" if np.isfinite(v_new) else " (median)"))
    lo = (g <= used).mean()
    ax.text(0.98, 0.97, f"low {lo:.0%} | high {1 - lo:.0%}", transform=ax.transAxes, ha="right", va="top", fontsize=5.5 if small else 6.5, color="#C0392B")
    ax.set_xlim(0, np.percentile(g, 99.5)); ax.set_ylim(0, d.max() * 1.15); ax.set_yticks([])
    ax.set_title(f"{DISPLAY[s]} ({s.replace('24090566_AC_', '')})", fontsize=7, color="0.15", loc="left")
    ax.legend(loc="upper right", bbox_to_anchor=(1.0, 0.9), fontsize=5.5 if small else 6.5, handlelength=1.4, frameon=False)
    return dict(SampleName=s, display=DISPLAY[s], glial_peak=gp, neuronal_peak=rp, threshold_used=used,
                threshold_rule=thr.loc[s, "threshold_rule"], peaks=";".join(f"{p:.0f}" for p in x[peaks]), frac_low=lo)


rows = []
for s in DISPLAY:
    fig, ax = plt.subplots(figsize=(4.2, 2.8))
    rows.append(panel(ax, s)); ax.set_xlabel("Genes per nucleus (n_genes_by_counts)"); ax.set_ylabel("Density")
    fig.tight_layout(); fig.savefig(QC / f"{s}_kde_split.pdf"); plt.close(fig)
fig, axs = plt.subplots(3, 3, figsize=(7.2, 5.6))
for ax, s in zip(axs.ravel(), DISPLAY):
    panel(ax, s, small=True)
for ax in axs[-1]:
    ax.set_xlabel("Genes per nucleus")
for ax in axs[:, 0]:
    ax.set_ylabel("Density")
fig.suptitle("Low/high-feature split per animal: KDE valley threshold (red)", fontsize=7.5)
fig.tight_layout(rect=(0, 0, 1, 0.96))
fig.savefig(QC / "kde_split_all_animals.pdf"); fig.savefig(QC / "kde_split_all_animals.png", dpi=300); plt.close(fig)
pd.DataFrame(rows).to_csv(QC / "kde_split_summary.csv", index=False)
print(pd.DataFrame(rows).round(1).to_string(index=False))
