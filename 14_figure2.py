#!/usr/bin/env python3
"""Step 14 - Figure 2: neuron subtypes (UMAP, composition), DEGs per subtype, CES/CPS DEG overlap, Cdh9+ DEG and GO overlap."""
import sys, json
from pathlib import Path
import numpy as np, pandas as pd, anndata as ad
import matplotlib.pyplot as plt
from scipy import stats
sys.path.insert(0, str(Path(__file__).parent))
from common import *  # noqa
from config import *  # project paths

P = json.load(open(CODE / "params.json"))
LFC, PADJ = P["deg_log2fc"], P["deg_padj"]
neu = ad.read_h5ad(H5AD / "ACC_neurons.h5ad", backed="r")
obs = neu.obs[["sample", "Experiment", "major_celltype", "subtype"]].copy()
obs["subtype"] = obs["subtype"].astype(str)
xy = np.asarray(neu.obsm["X_umap"])
types = [t for t in EX_TYPES + IN_TYPES if (obs.subtype == t).any()]
n_by = obs.subtype.value_counts()

# ---- 2A neuron subtype UMAP --------------------------------------------------
fig, ax = plt.subplots(figsize=(70 * MM, 62 * MM))
umap_scatter(ax, xy, obs.subtype.map(SUBTYPE), s=0.3)
label_centroids(ax, xy, obs.subtype, types, fontsize=5)
umap_axes_arrows(ax)
ax.set_title(f"Neurons, {len(obs):,} nuclei", fontsize=7)
ex = [t for t in types if t in EX_TYPES]; inn = [t for t in types if t in IN_TYPES]
h = lambda ts: [plt.Line2D([], [], marker="o", ls="", ms=4, mfc=SUBTYPE[t], mec="none", label=f"{t} ({n_by[t]:,})") for t in ts]
l1 = ax.legend(handles=h(ex), title="ExNeuron", loc="upper left", bbox_to_anchor=(1.0, 1.02), handletextpad=0.3, labelspacing=0.3)
ax.add_artist(l1)
ax.legend(handles=h(inn), title="InNeuron", loc="upper left", bbox_to_anchor=(1.62, 1.02), handletextpad=0.3, labelspacing=0.3)
src = obs.copy(); src[["UMAP1", "UMAP2"]] = xy
save(fig, "Fig2A_umap_neuron_subtypes", src.reset_index(names="cell"))

# ---- 2B neuron subtype composition (within neurons) --------------------------
ct = pd.crosstab(obs["sample"], obs.subtype).reindex(index=SAMPLE_ORDER, columns=types)
grp = pd.Series([s.split("_")[0] for s in ct.index], index=ct.index)
fig, a1 = plt.subplots(figsize=(88 * MM, 55 * MM))
frac = ct.div(ct.sum(1), axis=0); bottom = np.zeros(len(frac))
for t in types:
    a1.bar(range(len(frac)), frac[t], bottom=bottom, color=SUBTYPE[t], width=0.8, lw=0)
    bottom += frac[t].values
a1.set_xticks(range(len(frac)), frac.index, rotation=45, ha="right")
for tl in a1.get_xticklabels():
    tl.set_color(GROUP[tl.get_text().split("_")[0]])
a1.set_ylim(0, 1); a1.set_ylabel("Fraction of neurons per animal")
a1.yaxis.set_major_formatter(plt.FuncFormatter(lambda v, _: f"{v:.0%}"))
# replicate-level test (not drawn; reported in table2 subtype_test)
rows = []
tot = ct.sum(1); prop = (ct + 0.5).div(tot + 0.5 * ct.shape[1], axis=0); lg = np.log(prop / (1 - prop))
for cond in ["CES", "CPS"]:
    r = logit_prop_test(ct, grp, cond)
    for i, t in enumerate(types):
        a, b = lg.loc[grp == "CON", t], lg.loc[grp == cond, t]
        va, vb = a.var(ddof=1) / 3, b.var(ddof=1) / 3; se = np.sqrt(va + vb)
        dfw = (va + vb) ** 2 / (va ** 2 / 2 + vb ** 2 / 2)
        r.loc[i, "ci_low"] = r.loc[i, "logit_diff"] - stats.t.ppf(0.975, dfw) * se
        r.loc[i, "ci_high"] = r.loc[i, "logit_diff"] + stats.t.ppf(0.975, dfw) * se
    rows.append(r)
forest = pd.concat(rows, ignore_index=True)
a1.legend(handles=[plt.Rectangle((0, 0), 1, 1, fc=SUBTYPE[t]) for t in types[::-1]], labels=list(types)[::-1],   # same order as the stack
          loc="upper left", bbox_to_anchor=(1.02, 1.0), ncol=1, handlelength=0.8, labelspacing=0.25, fontsize=5.5, frameon=False)
save(fig, "Fig2B_composition_neuron_subtypes", {"per_sample_counts": ct.reset_index(), "logit_test": forest})

# ---- 2C mirrored DEG counts per subtype --------------------------------------
de = pd.read_csv(DATA / "DE_celllevel_subtype.csv.gz")
de = de[de.cluster.isin(types)]
sig = de[(de.pvals_adj < PADJ) & (de.log2FC.abs() > LFC) & (de.reference == "CON")]
cnt = sig.assign(dir=np.where(sig.log2FC > 0, "up", "down")).groupby(["condition", "cluster", "dir"]).size() \
    .unstack(["condition", "dir"]).reindex(types).fillna(0).astype(int)
for c in [("CES", "up"), ("CES", "down"), ("CPS", "up"), ("CPS", "down")]:
    if c not in cnt.columns:
        cnt[c] = 0
order = cnt.sum(1).sort_values(ascending=False).index.tolist()
cnt = cnt.loc[order]
ncell = obs.groupby(["subtype", "Experiment"], observed=True).size().unstack().reindex(order)
fig = plt.figure(figsize=(120 * MM, 118 * MM))
gs = fig.add_gridspec(5, 1, height_ratios=[0.42, 0.16, 1, 0.075, 1], hspace=0)
an = fig.add_subplot(gs[0]); at = fig.add_subplot(gs[2], sharex=an); ad = fig.add_subplot(gs[3], sharex=an); ab = fig.add_subplot(gs[4], sharex=an)
x = np.arange(len(order)); w = 0.3
# top annotation: nuclei per subtype and group (log10 axis)
for off, g in [(-0.24, "CON"), (0.0, "CES"), (0.24, "CPS")]:
    an.scatter(x + off, ncell[g].values, s=9, color=GROUP[g], lw=0, zorder=3, label=g)
an.set_yscale("log"); an.set_ylabel("Nuclei"); an.grid(axis="y", color="0.92", lw=0.5, zorder=0)
an.set_ylim(ncell.values.min() / 1.8, ncell.values.max() * 1.8)
an.tick_params(axis="x", bottom=False, labelbottom=False); an.spines["bottom"].set_visible(False)
an.legend(loc="upper right", ncol=3, handletextpad=0.1, columnspacing=0.8, borderaxespad=0.1, fontsize=5.5)
an.set_title(f"DEGs per neuron subtype (adj. P < {PADJ}, |log2FC| > {LFC})", fontsize=7)
# mirrored DEG counts: boxed panels, CES on top, CPS below (negative axis), subtype dots on the shared zero line
kw = dict(alpha=0.75, lw=0, zorder=2)
at.bar(x - w / 2 - 0.02, cnt[("CES", "up")], w, color=UP, label="Up", **kw); at.bar(x + w / 2 + 0.02, cnt[("CES", "down")], w, color=DOWN, label="Down", **kw)
ab.bar(x - w / 2 - 0.02, -cnt[("CPS", "up")], w, color=UP, **kw); ab.bar(x + w / 2 + 0.02, -cnt[("CPS", "down")], w, color=DOWN, **kw)
for a_ in (at, ab):
    for sp in a_.spines.values():
        sp.set_visible(True); sp.set_linewidth(1.0); sp.set_color("0.15")
    a_.set_ylabel("DEGs")
at.set_ylim(0, cnt[[("CES", "up"), ("CES", "down")]].values.max() * 1.08)
ab.set_ylim(-cnt[[("CPS", "up"), ("CPS", "down")]].values.max() * 1.08, 0)
at.set_title("CES vs CON", fontsize=7, pad=3)
at.tick_params(axis="x", bottom=False, labelbottom=False)
at.legend(loc="upper right", frameon=False, handlelength=0.6, handletextpad=0.4, fontsize=6)
ad.scatter(x, np.zeros(len(x)), s=26, c=[SUBTYPE[t] for t in order], edgecolors="white", linewidths=0.4, zorder=3, clip_on=False)
ad.set_ylim(-1, 1); ad.axis("off")
ab.set_xticks(x, order, rotation=60, ha="right"); ab.tick_params(axis="x", length=0)
ab.set_xlabel("CPS vs CON", fontsize=7, labelpad=4)
ab.yaxis.set_major_formatter(plt.FuncFormatter(lambda v, _: f"{abs(v):g}"))   # bars point down, counts shown positive
ab.set_xlim(-0.6, len(order) - 0.4)
fig.align_ylabels([an, at, ab])
src = cnt.copy(); src.columns = [f"{a}_{b}" for a, b in src.columns]; src = src.join(ncell.add_prefix("n_cells_"))
save(fig, "Fig2C_DEG_counts_neuron_subtypes", src.reset_index(names="subtype"))

# ---- 2D Jaccard similarity of CES and CPS DEG sets ---------------------------
jac = []
for t in types:
    A = set(sig[(sig.cluster == t) & (sig.condition == "CES")].gene); B = set(sig[(sig.cluster == t) & (sig.condition == "CPS")].gene)
    jac.append(dict(subtype=t, n_CES=len(A), n_CPS=len(B), n_shared=len(A & B),
                    jaccard=len(A & B) / len(A | B) if A | B else np.nan,
                    lineage="ExNeuron" if t in EX_TYPES else "InNeuron"))
jac = pd.DataFrame(jac).sort_values(["lineage", "jaccard"], ascending=[True, False])
fig, ax = plt.subplots(figsize=(95 * MM, 55 * MM))
ax.bar(range(len(jac)), jac.jaccard, color=[SUBTYPE[t] for t in jac.subtype], width=0.7)
for i, r in enumerate(jac.itertuples()):
    ax.text(i, r.jaccard, f"{r.n_shared}", ha="center", va="bottom", fontsize=4.5, color="0.3")
ax.set_xticks(range(len(jac)), jac.subtype, rotation=60, ha="right")
ax.set_ylabel("Jaccard index\n(CES vs CON DEGs, CPS vs CON DEGs)")
ne = (jac.lineage == "ExNeuron").sum()
ax.axvline(ne - 0.5, color="0.6", lw=0.6, ls=":")
ax.text(ne / 2 - 0.5, ax.get_ylim()[1] * 0.97, "ExNeuron", ha="center", fontsize=6)
ax.text(ne + (len(jac) - ne) / 2 - 0.5, ax.get_ylim()[1] * 0.97, "InNeuron", ha="center", fontsize=6)
ax.set_title("CES/CPS DEG overlap per subtype (bar labels = shared DEGs)", fontsize=6.5)
save(fig, "Fig2D_DEG_jaccard_CES_CPS", jac)

# ---- 2E Cdh9+ DEG Venn -------------------------------------------------------
c = sig[sig.cluster == "Cdh9+"]
A = set(c[c.condition == "CES"].gene); B = set(c[c.condition == "CPS"].gene)
fig, ax = plt.subplots(figsize=(60 * MM, 45 * MM))
venn2(ax, A, B, ["CES vs CON", "CPS vs CON"], [GROUP["CES"], GROUP["CPS"]], title="Cdh9+ DEGs")
memb = pd.DataFrame({"gene": sorted(A | B)}); memb["CES"] = memb.gene.isin(A); memb["CPS"] = memb.gene.isin(B)
save(fig, "Fig2E_Cdh9_DEG_venn", memb)

# ---- 2F Cdh9+ GO overlap (up, down, all DEGs) -----------------------------------------------------
go = pd.read_csv(DATA / "GO_per_contrast_ExNeuron_Cdh9.csv")
go_overlap_panel(go, "Cdh9+", "Fig2F_Cdh9_GO_overlap", dirs=("up", "down", "all"))   # bottom row: up + down combined

print("Fig2 done")

# ---- Supplementary Fig. S3: marker matrixplot of neuron subtypes -------------
# Markers = the genes in the subtype names (in subtype order) plus lineage markers; matrix plot of mean
# log-normalised expression per subtype, min-max scaled per gene (standard_scale="var").
def _name_genes(ts):
    out = []
    for t in ts:
        for g in t.replace("+", " ").split():
            if g not in out:
                out.append(g)
    return out
mk = {"Lineage": ["Slc17a7", "Gad1", "Gad2"], "ExNeuron subtypes": _name_genes(ex), "InNeuron subtypes": _name_genes(inn)}
mk = {k: [g for g in v if g in neu.var_names] for k, v in mk.items()}
sub = neu[:, [g for v in mk.values() for g in v]].to_memory()
sub.obs["subtype"] = pd.Categorical(sub.obs["subtype"].astype(str), categories=types)
genes_all = [g for v in mk.values() for g in v]
X = sub.X; mean = pd.DataFrame(np.vstack([np.asarray(X[(sub.obs.subtype == t).values].mean(0)).ravel() for t in types]), index=types, columns=genes_all)
scaled = (mean - mean.min()) / (mean.max() - mean.min())
frac = pd.DataFrame(np.vstack([np.asarray((X[(sub.obs.subtype == t).values] > 0).mean(0)).ravel() for t in types]), index=types, columns=genes_all)
grp = pd.Series({g: k for k, v in mk.items() for g in v})
src = pd.concat({"mean_logexpr": mean, "scaled": scaled, "frac_expressing": frac}, names=["value", "subtype"]).reset_index()
from matplotlib.colors import ListedColormap
# matrix plot (same values as sc.pl.matrixplot(standard_scale="var")) drawn as blocks so that ExNeuron and
# InNeuron rows, and the lineage / ExNeuron / InNeuron marker columns, are separated by gaps
row_groups = [("ExNeuron", ex), ("InNeuron", inn)]
col_groups = list(mk.items())
fig = plt.figure(figsize=(150 * MM, 70 * MM))
gs = fig.add_gridspec(2, len(col_groups) + 1, height_ratios=[len(r) for _, r in row_groups],
                      width_ratios=[len(c) for _, c in col_groups] + [1.2], hspace=0.06, wspace=0.04,
                      left=0.13, right=0.97, top=0.88, bottom=0.17)
for i, (rname, rows_) in enumerate(row_groups):
    for j, (cname, cols_) in enumerate(col_groups):
        ax_ = fig.add_subplot(gs[i, j])
        cross = (rname == "ExNeuron" and cname == "InNeuron subtypes") or (rname == "InNeuron" and cname == "ExNeuron subtypes")
        if cross:   # other-lineage subtype markers: not informative for this lineage, shown in grey
            ax_.pcolormesh(np.zeros((len(rows_), len(cols_))), cmap=ListedColormap(["#D9D9D9"]), edgecolors="white", linewidth=0.4)
        else:
            im = ax_.pcolormesh(scaled.loc[rows_, cols_].values, cmap=CMAP_SEQ, vmin=0, vmax=1, edgecolors="white", linewidth=0.4)
        ax_.set_xlim(0, len(cols_)); ax_.set_ylim(len(rows_), 0)
        ax_.set_xticks(np.arange(len(cols_)) + 0.5); ax_.set_yticks(np.arange(len(rows_)) + 0.5)
        ax_.tick_params(length=0, pad=2)
        for sp in ax_.spines.values():
            sp.set_visible(False)
        if j == 0:
            ax_.set_yticklabels(rows_, fontsize=6)
            ax_.text(-0.42 * max(len(c) for _, c in col_groups[:1]) - 3.2, len(rows_) / 2, rname, rotation=90, va="center", ha="center", fontsize=6.5)
        else:
            ax_.set_yticklabels([])
        if i == len(row_groups) - 1:
            ax_.set_xticklabels(cols_, rotation=90, fontsize=6, fontstyle="italic")
        else:
            ax_.set_xticklabels([])
            ax_.set_title(cname, fontsize=6.5, pad=3)
cax = fig.add_subplot(gs[:, -1]); cax.axis("off")
cb = fig.colorbar(im, ax=cax, fraction=0.9, aspect=12, location="right")
cb.set_label("Mean expression\n(scaled per gene)", fontsize=6); cb.ax.tick_params(labelsize=5.5)
save(fig, "FigS3_neuron_subtype_markers", src)
print("FigS3 markers:", {k: len(v) for k, v in mk.items()})
