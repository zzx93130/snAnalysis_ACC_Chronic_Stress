#!/usr/bin/env python3
"""Step 12 - Figure 1A-D and 1F (Python): atlas UMAP, marker heatmap, UMAP by group, DEG UpSet, ExNeuron GO.
Figure 1E and 1G are drawn in R by 13_figure1EG.R."""
import sys, json
from pathlib import Path
import numpy as np, pandas as pd, scanpy as sc
import matplotlib.pyplot as plt
from scipy import stats
sys.path.insert(0, str(Path(__file__).parent))
from common import *  # noqa
from config import *  # project paths

P = json.load(open(CODE / "params.json"))
LFC, PADJ = P["deg_log2fc"], P["deg_padj"]
md = pd.read_csv(DATA / "ACC_cell_metadata.csv", index_col=0)
xy = md[["UMAP1", "UMAP2"]].values
present = [k for k in MAJOR_ORDER if (md.major_celltype == k).any()]

# ---- 1A all-cell UMAP --------------------------------------------------------
fig, ax = plt.subplots(figsize=(70 * MM, 60 * MM))
umap_scatter(ax, xy, md.major_celltype.map(MAJOR), s=0.25)
label_centroids(ax, xy, md.major_celltype, present)
umap_axes_arrows(ax)
ax.set_title(f"ACC, {len(md):,} nuclei, 9 animals", fontsize=7)
counts = md.major_celltype.value_counts().reindex(present)
h = [plt.Line2D([], [], marker="o", ls="", ms=4, mfc=MAJOR[k], mec="none", label=f"{k} ({counts[k]:,})") for k in present]
ax.legend(handles=h, loc="center left", bbox_to_anchor=(1.0, 0.5), handletextpad=0.3, labelspacing=0.4)
save(fig, "Fig1A_umap_major_celltypes",
     md[["sample", "Experiment", "major_celltype", "UMAP1", "UMAP2"]].reset_index(names="cell"))

# ---- 1B marker heatmap -------------------------------------------------------
MARKERS = dict(ExNeuron=["Slc17a7", "Slc17a6"], InNeuron=["Gad1", "Gad2"], Astro=["Slc1a2", "Slc1a3"],
               Oligo=["Mbp", "St18"], OPC=["Vcan"], Micro=["Hexb", "Tgfbr1", "Inpp5d"],
               Endo=["Cldn5", "Pecam1", "Kdr"], Mural=["Rgs5", "Pdgfrb", "Kcnj8"])
ad = sc.read_h5ad(H5AD / "ACC_all_cells.h5ad", backed="r")
genes = [g for k in present for g in MARKERS[k] if g in ad.var_names]
X = ad[:, genes].to_memory().X
X = X.toarray() if hasattr(X, "toarray") else X
lab = ad.obs["major_celltype"].astype(str).values
mean = pd.DataFrame([X[lab == k].mean(0) for k in present], index=present, columns=genes)
frac = pd.DataFrame([(X[lab == k] > 0).mean(0) for k in present], index=present, columns=genes)
scaled = (mean - mean.min()) / (mean.max() - mean.min())
fig, ax = plt.subplots(figsize=(4.2 * len(genes) * MM + 25 * MM, 4.2 * len(present) * MM + 18 * MM))
im = ax.imshow(scaled.values, cmap="viridis", vmin=0, vmax=1, aspect="auto", interpolation="nearest")
ax.set_xticks(range(len(genes)), genes, rotation=90, style="italic")
ax.set_yticks(range(len(present)), present)
ax.tick_params(length=0)
for sp in ax.spines.values():
    sp.set_visible(False)
ax.set_xticks(np.arange(len(genes) + 1) - .5, minor=True); ax.set_yticks(np.arange(len(present) + 1) - .5, minor=True)
ax.grid(which="minor", color="white", lw=0.6); ax.tick_params(which="minor", length=0)
st = 0
for k in present:
    n = len([g for g in MARKERS[k] if g in genes])
    ax.plot([st - .4, st + n - .6], [-0.75, -0.75], color=MAJOR[k], lw=3, solid_capstyle="butt", clip_on=False)
    st += n
cb = fig.colorbar(im, ax=ax, fraction=0.04, pad=0.03, ticks=[0, .5, 1])
cb.set_label("Mean expression\n(min-max scaled per gene)", fontsize=6); cb.outline.set_visible(False)
src = mean.stack().rename("mean_logexpr").to_frame().join(frac.stack().rename("frac_expressing")) \
    .join(scaled.stack().rename("scaled")).reset_index(names=["celltype", "gene"])
save(fig, "Fig1B_marker_heatmap", src)

# ---- 1C UMAP coloured by group (scanpy sc.pl.umap, random draw order) ----------------------------------
def fig1c():
    md = pd.read_csv(DATA / "ACC_cell_metadata.csv", index_col=0)
    md = md.iloc[np.random.default_rng(SEED).permutation(len(md))]          # random draw order
    a = sc.AnnData(obs=md[["sample", "Experiment", "Experiment2"]].copy())
    a.obsm["X_umap"] = md[["UMAP1", "UMAP2"]].to_numpy()
    a.obs["Experiment"] = pd.Categorical(a.obs.Experiment, categories=GROUP_ORDER)
    a.obs["Experiment2"] = pd.Categorical(a.obs.Experiment2, categories=["CON", "SDS"])
    n1 = a.obs.Experiment.value_counts(); n2 = a.obs.Experiment2.value_counts()
    a.obs["group_lab"] = a.obs.Experiment.map({g: f"{g} (n = {n1[g]:,})" for g in GROUP_ORDER}).astype("category")
    a.obs["group_lab"] = a.obs.group_lab.cat.reorder_categories([f"{g} (n = {n1[g]:,})" for g in GROUP_ORDER])
    a.obs["sds_lab"] = a.obs.Experiment2.map({g: f"{g} (n = {n2[g]:,})" for g in ["CON", "SDS"]}).astype("category")
    a.obs["sds_lab"] = a.obs.sds_lab.cat.reorder_categories([f"{g} (n = {n2[g]:,})" for g in ["CON", "SDS"]])

    fig, axs = plt.subplots(1, 2, figsize=(150 * MM, 60 * MM), gridspec_kw=dict(wspace=0.55))
    pt = 120000 / len(a)                                       # scanpy default point size
    for ax, key, pal, title in [(axs[0], "group_lab", [GROUP[g] for g in GROUP_ORDER], "CON / CES / CPS"),
                                (axs[1], "sds_lab", [GROUP["CON"], GROUP["SDS"]], "CON / SDS (CES + CPS)")]:
        sc.pl.umap(a, color=key, palette=pal, ax=ax, show=False, frameon=False, size=pt * 0.3,
                   title=title, legend_loc="right margin", legend_fontsize=6)
        for c in ax.collections:
            c.set_rasterized(True); c.set_alpha(0.55)
        ax.title.set_fontsize(7)
        lg = ax.get_legend()
        for h in lg.legend_handles:
            h.set_alpha(1); h.set_sizes([12])
    umap_axes_arrows(axs[0])
    src = md[["sample", "Experiment", "Experiment2", "UMAP1", "UMAP2"]].sort_index().reset_index(names="cell")
    save(fig, "Fig1C_umap_combined", src)
    print(n1.to_dict(), n2.to_dict())


fig1c()

# ---- 1D UpSet of DEGs per major type: CES and CPS as paired bars over one dot matrix --------------------
def fig1d():
    TOP = 12
    de = pd.read_csv(DATA / "DE_celllevel_major.csv.gz")
    sig = de[(de.pvals_adj < PADJ) & (de.log2FC.abs() > LFC) & (de.reference == "CON")]
    order = [k for k in MAJOR_ORDER if k in set(de.cluster)]
    conds = ["CES", "CPS"]

    combo, sizes = {}, {}
    for c in conds:
        s = sig[sig.condition == c]
        sets = {k: set(s[s.cluster == k].gene) for k in order}
        sizes[c] = [len(sets[k]) for k in order]
        allg = sorted(set().union(*sets.values()))
        memb = pd.DataFrame({k: [g in sets[k] for g in allg] for k in order}, index=allg)
        combo[c] = memb.groupby(order).size().rename(c)
    tab = pd.concat(combo.values(), axis=1).fillna(0).astype(int)
    tab = tab[np.asarray(tab.index.to_frame().any(axis=1))]
    keep = set(tab.sort_values("CES", ascending=False).head(TOP).index) | set(tab.sort_values("CPS", ascending=False).head(TOP).index)
    tab = tab.loc[[i for i in tab.index if i in keep]]
    tab = tab.loc[tab.sum(1).sort_values(ascending=False).index]

    nb = len(tab); x = np.arange(nb); w = 0.38
    fig = plt.figure(figsize=(max(120, 40 + nb * 6.2) * MM, 64 * MM))
    gs = fig.add_gridspec(2, 3, height_ratios=[1.35, 1], width_ratios=[1.0, 0.62, nb * 0.3], hspace=0.05, wspace=0.02)
    axb = fig.add_subplot(gs[0, 2]); axm = fig.add_subplot(gs[1, 2], sharex=axb)
    axs = fig.add_subplot(gs[1, 0], sharey=axm); axn = fig.add_subplot(gs[1, 1], sharey=axm); axn.axis("off")
    for off, c in [(-w / 2, "CES"), (w / 2, "CPS")]:
        axb.bar(x + off, tab[c], width=w, color=GROUP[c], label=f"{c} vs CON")
        for xi, n in zip(x + off, tab[c]):
            axb.text(xi, n, str(n), ha="center", va="bottom", fontsize=4.2, rotation=90 if n >= 1000 else 0)
    axb.set_ylabel("Intersection size (DEGs)"); axb.tick_params(axis="x", bottom=False, labelbottom=False)
    axb.spines["bottom"].set_visible(False)
    axb.legend(loc="upper right", handlelength=1, borderaxespad=0)
    axb.set_title(f"DEGs vs CON per major cell type (adj. P < {PADJ}, |log2FC| > {LFC})", fontsize=7, loc="left")
    ys = np.arange(len(order))[::-1]
    for i in range(len(order)):
        if i % 2 == 0:
            axm.axhspan(ys[i] - 0.5, ys[i] + 0.5, color="0.96", lw=0, zorder=0)
    for xi, key in zip(x, tab.index):
        on = [ys[i] for i, k in enumerate(order) if key[i]]
        axm.scatter([xi] * len(order), ys, s=12, c="0.86", lw=0, zorder=2)
        axm.scatter([xi] * len(on), on, s=12, c="0.2", lw=0, zorder=3)
        if len(on) > 1:
            axm.plot([xi, xi], [min(on), max(on)], c="0.2", lw=0.9, zorder=3)
    axm.set_yticks(ys, order); axm.tick_params(axis="both", left=False, labelleft=False, bottom=False, labelbottom=False)
    for sp in axm.spines.values():
        sp.set_visible(False)
    axm.set_ylim(-0.6, len(order) - 0.4); axm.set_xlim(-0.6, nb - 0.4)
    h = 0.4
    for off, c in [(h / 2, "CES"), (-h / 2, "CPS")]:
        axs.barh(ys + off, sizes[c], height=h, color=GROUP[c])
        for y, s in zip(ys + off, sizes[c]):
            axs.text(s, y, f"{s:,} ", ha="right", va="center", fontsize=4.2, color=GROUP[c])
    axs.set_xlim(max(max(v) for v in sizes.values()) * 1.35, 0); axs.tick_params(axis="y", left=False, labelleft=False); axs.set_xlabel("DEGs per cell type")
    axs.spines["left"].set_visible(False); axs.spines["right"].set_visible(False)
    for y, k in zip(ys, order):
        axn.add_patch(plt.Rectangle((0.02, y - 0.22), 0.07, 0.44, color=MAJOR[k], transform=axn.transData, clip_on=False))
        axn.text(0.14, y, k, ha="left", va="center", fontsize=6)
    axn.set_xlim(0, 1)
    src = tab.reset_index()
    save(fig, "Fig1D_upset_DEGs_merged_CES_CPS", src)
    print(tab.to_string())


fig1d()
de = pd.read_csv(DATA / "DE_celllevel_major.csv.gz")
sig = de[(de.pvals_adj < PADJ) & (de.log2FC.abs() > LFC) & (de.reference == "CON")]
cnt = sig.groupby(["condition", "cluster"]).size().unstack(0).reindex(present).fillna(0).astype(int)
ncell = md.groupby(["major_celltype", "Experiment"]).size().unstack().reindex(present)
cnt.join(ncell.add_prefix("n_cells_")).to_csv(FIGS / "Fig1D_DEG_counts_per_celltype.csv")


# ---- 1F ExNeuron GO overlap (CES vs CPS) + shared terms ----------------------
go = pd.read_csv(DATA / "GO_per_contrast_ExNeuron_Cdh9.csv")
go_overlap_panel(go, "ExNeuron", "Fig1F_ExNeuron_GO_overlap", dirs=("up", "down", "all"))   # bottom row: up + down combined
print("Fig1 done")
