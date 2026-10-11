#!/usr/bin/env python3
"""Step 15 - Figure 3: shared CES/CPS response of Cdh9+ ExNeurons (GO of shared DEGs, log2FC concordance,
CellChat, PROGENy, gene violins)."""
import sys, json, re
from pathlib import Path
import numpy as np, pandas as pd, anndata as ad
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from scipy import stats
sys.path.insert(0, str(Path(__file__).parent))
from common import *  # noqa
from config import *  # project paths

P = json.load(open(CODE / "params.json"))
LFC, PADJ = P["deg_log2fc"], P["deg_padj"]

# ---- 3A GO of shared concordant Cdh9+ DEGs: representative terms (step 6), fold enrichment, gene matrix ----
D = DATA / "fig3a_followup"; N_DOWN = 10
reps = pd.read_csv(D / "fig3a_representative_terms.csv")
ora = pd.read_csv(DATA / "GO_Cdh9_shared_DEGs.csv")
gw = pd.read_csv(D / "gsea_pseudobulk_Cdh9_CES_CPS_wide.csv")
conc = set(gw.loc[gw.concordant == True, "Description"])
sh = pd.read_csv(DATA / "Cdh9_shared_DEGs.csv").set_index("gene")
sh["mean_lfc"] = sh[["CES", "CPS"]].mean(1)
nu = sh["class"].value_counts()


def select(st):
    r = reps[reps.set == st].sort_values("p.adjust")
    if st == "Common down":
        r = r.head(N_DOWN)
    r = r.copy()
    r["gsea_concordant"] = [any(t in conc for t in [d] + (str(rt).split("; ") if isinstance(rt, str) else []))
                            for d, rt in zip(r.Description, r.related_terms)]
    r["label"] = [d if len(d) <= 38 else d[:36] + "..." for d in r.Description]
    r["label"] = [f"{l}  (+{n})" if n > 0 else l for l, n in zip(r.label, r.n_related)]
    return r.reset_index(drop=True)


def gene_cols(r, nmax):
    genes = pd.Series([g for s in r.geneID for g in s.split("/")]).value_counts()
    keep = genes[genes >= 2].index.tolist()
    if len(keep) < 6:
        keep = genes.index[:min(nmax, len(genes))].tolist()
    keep = keep[:nmax]
    return sorted(keep, key=lambda g: (-genes[g], g))


sel = {st: select(st) for st in ["Common up", "Common down"]}
nrow = max(len(v) for v in sel.values())
xmax = 5 * np.ceil(max(v.FoldEnrichment.max() for v in sel.values()) / 5)
pmax = max((-np.log10(v["p.adjust"])).max() for v in sel.values())
lmax = np.abs(sh.loc[sh["class"] != "Opposite direction", "mean_lfc"]).quantile(0.98)
genes = {st: gene_cols(sel[st], 16 if st == "Common up" else 22) for st in sel}

ROW = {"Common up": 5.0, "Common down": 3.6}          # row height (mm): the down part is drawn more compactly
TOPM, BOTM = 14, 17                                   # mm above (titles) and below (gene names, legend) the panels
size = lambda c: 4 + 2.6 * np.asarray(c)            # dot area per gene count
area = max(ROW[st] * (len(sel[st]) + 0.4) for st in sel)
H = TOPM + area + BOTM
fig = plt.figure(figsize=(183 * MM, H * MM))
top = 1 - TOPM / H
outer = fig.add_gridspec(1, 2, wspace=0.85, left=0.2, right=0.93, top=top, bottom=top - area / H)
src = []
for k, (st, col) in enumerate([("Common up", UP), ("Common down", DOWN)]):
    r = sel[st]; g = genes[st]
    inner = outer[k].subgridspec(1, 2, width_ratios=[1, 0.12 * len(g) + 0.1], wspace=0.1)
    ax = fig.add_subplot(inner[0]); axm = fig.add_subplot(inner[1], sharey=ax)
    y = np.arange(len(r))[::-1]
    sca = ax.scatter(r.FoldEnrichment, y, s=size(r.Count), c=-np.log10(r["p.adjust"]), cmap=CMAP_SEQ, vmin=0, vmax=pmax, lw=0, zorder=3)
    ax.set_yticks(y, r.label, fontsize=5.3); ax.set_ylim(-0.7, len(r) - 0.3)
    ax.set_xlim(0, xmax); ax.set_xlabel("Fold enrichment", fontsize=6); ax.tick_params(axis="x", labelsize=5.5)
    ax.grid(axis="x", color="0.92", lw=0.5, zorder=0)
    ax.set_title(f"{st} in CES and CPS (n = {nu.get(st, 0)})", fontsize=7, color=col, weight="bold", loc="right")
    M = np.full((len(r), len(g)), np.nan)
    for i, gs in enumerate(r.geneID):
        s_ = set(gs.split("/"))
        for j, gene in enumerate(g):
            if gene in s_:
                M[i, j] = sh.loc[gene, "mean_lfc"]
    im = axm.imshow(M[::-1], cmap="RdBu_r", norm=TwoSlopeNorm(0, -lmax, lmax), aspect="auto",
                    extent=(-0.5, len(g) - 0.5, -0.5, len(r) - 0.5), interpolation="nearest")
    axm.set_xticks(range(len(g)), g, rotation=90, fontsize=4.8, style="italic")
    axm.tick_params(axis="y", left=False, labelleft=False); axm.tick_params(axis="x", length=0)
    axm.set_xticks(np.arange(len(g) + 1) - 0.5, minor=True); axm.set_yticks(np.arange(len(r) + 1) - 0.5, minor=True)
    axm.grid(which="minor", color="0.88", lw=0.3); axm.tick_params(which="minor", length=0)
    for sp in axm.spines.values():
        sp.set_visible(False)
    hk = ROW[st] * (len(r) + 0.4) / H                    # top-aligned panel of its own height
    for a in (ax, axm):
        b = a.get_position(); a.set_position([b.x0, top - hk, b.width, hk])
    src.append(r.assign(genes_shown=";".join(g)))
cb_h = min(0.3, area / H * 0.42)
cax = fig.add_axes([0.955, top - cb_h, 0.012, cb_h]); cb = fig.colorbar(sca, cax=cax); cb.set_label("-log10 adj. P", fontsize=5.5); cb.ax.tick_params(labelsize=5); cb.outline.set_visible(False)
cax2 = fig.add_axes([0.955, top - 2.25 * cb_h, 0.012, cb_h]); cb2 = fig.colorbar(im, cax=cax2); cb2.set_label("mean log2FC\n(CES, CPS)", fontsize=5.5); cb2.ax.tick_params(labelsize=5); cb2.outline.set_visible(False)
h = [plt.scatter([], [], s=size(c), c="0.5") for c in [5, 15, 30]]
fig.legend(h, ["5 genes", "15 genes", "30 genes"], loc="upper left", bbox_to_anchor=(0.2, top - area / H - 11 / H), ncol=3, fontsize=5.5, frameon=False)
fig.suptitle("GO enrichment of Cdh9+ DEGs shared by CES and CPS (redundant terms merged, +n = merged terms)", fontsize=7.5, y=1 - 4 / H)
save(fig, "Fig3A_Cdh9_shared_DEG_GO", pd.concat(src))


# ---- 3B log2FC concordance scatter ------------------------------------------
w = pd.read_csv(DATA / "Cdh9_shared_DEGs.csv")
col = w["class"].map({"Common up": UP, "Common down": DOWN, "Opposite direction": NS})
fig, ax = plt.subplots(figsize=(75 * MM, 68 * MM))
ax.scatter(w.CES, w.CPS, c=col, s=7, lw=0, alpha=0.85)
lim = np.abs(w[["CES", "CPS"]].values).max() * 1.08
ax.plot([-lim, lim], [-lim, lim], ls="--", color="0.5", lw=0.6)
ax.axhline(0, color="0.8", lw=0.5); ax.axvline(0, color="0.8", lw=0.5)
w["dist"] = np.hypot(w.CES, w.CPS)
lab = pd.concat([w[w["class"] == c].nlargest(5, "dist") for c in ["Common up", "Common down"]] +
                [w[w["class"] == "Opposite direction"].nlargest(3, "dist")])
from adjustText import adjust_text
texts = [ax.text(r.CES, r.CPS, r.gene, fontsize=5, style="italic") for r in lab.itertuples()]
adjust_text(texts, ax=ax, arrowprops=dict(arrowstyle="-", lw=0.4, color="0.4"))
r, p = stats.pearsonr(w.CES, w.CPS)
ax.text(0.03, 0.97, f"Pearson r = {r:.2f}\nn = {len(w)} shared DEGs", transform=ax.transAxes, va="top", fontsize=6)
ax.set_xlabel("CES vs CON log2FC"); ax.set_ylabel("CPS vs CON log2FC")
cnts = w["class"].value_counts()
ax.legend(handles=[plt.Line2D([], [], marker="o", ls="", ms=4, mfc=c, mec="none", label=f"{k} ({cnts.get(k, 0)})")
                   for k, c in [("Common up", UP), ("Common down", DOWN), ("Opposite direction", NS)]],
          loc="lower right", fontsize=5.5)
ax.set_title("Cdh9+: DEGs shared by CES and CPS", fontsize=7)
save(fig, "Fig3B_Cdh9_shared_DEG_log2FC_scatter", w.drop(columns="dist"))

# ---- 3C CellChat differential number / strength: only edges to/from Cdh9+ ----
# Partners with few interactions with Cdh9+ are de-emphasised: node area and opacity follow the absolute
# number (or strength) of Cdh9+ <-> partner interactions, mean of CON and the stress group (hub_network(total=)).
# Two versions: fine groups (all subtypes) and merged groups (other ExNeuron, InNeuron, glia; CellChat run on
# the coarse grouping, data/cellchat_coarse).
def cdh9_hub_figure(cc, nodes_all, cols, name, figsize):
    M = {(c, m): pd.read_csv(cc / f"{c}_net_{m}.csv", index_col=0) for c in GROUP_ORDER for m in ["count", "weight"]}
    nodes = [t for t in nodes_all if t in M[("CON", "count")].index]
    fig, axs = plt.subplots(2, 2, figsize=figsize)
    srcs = []
    for j, (m, lab) in enumerate([("count", "number of interactions"), ("weight", "interaction strength")]):
        diff = {c: M[(c, m)] - M[("CON", m)] for c in ["CES", "CPS"]}
        vmax = max(max(np.abs(d.loc["Cdh9+"]).max(), np.abs(d["Cdh9+"]).max()) for d in diff.values())   # same scale CES/CPS
        for i, cond in enumerate(["CES", "CPS"]):
            tot = {k: np.mean([M[(c, m)].loc["Cdh9+", k] + M[(c, m)].loc[k, "Cdh9+"] for c in ["CON", cond]]) for k in nodes if k != "Cdh9+"}
            e = hub_network(axs[i, j], diff[cond], cols, "Cdh9+", [k for k in nodes if k != "Cdh9+"], vmax=vmax, total=tot,
                            title=f"{cond} vs CON: differential {lab}\n(Cdh9+ outgoing and incoming)")
            srcs.append(e.assign(comparison=f"{cond}-CON", measure=m,
                                 partner_total=[tot.get(r if s == "Cdh9+" else s, np.nan) for s, r in zip(e.source, e.target)]))
    fig.legend(handles=[plt.Line2D([], [], color=UP, lw=2, label="Increased vs CON"),
                        plt.Line2D([], [], color=DOWN, lw=2, label="Decreased vs CON"),
                        plt.Line2D([], [], ls="", marker="o", ms=6, mfc="0.5", mec="none", label="node size/opacity = interactions with Cdh9+")],
               loc="lower center", ncol=3)
    save(fig, name, pd.concat(srcs))
    return M


M = cdh9_hub_figure(DATA / "cellchat_fine", EX_TYPES + IN_TYPES + ["Astro", "Micro", "Oligo", "OPC", "Endo", "Mural"],
                    {**SUBTYPE, **MAJOR}, "Fig3C_cellchat_Cdh9_diff_count_weight", (165 * MM, 150 * MM))
cdh9_hub_figure(DATA / "cellchat_coarse", ["Cdh9+", "OtherExNeuron", "InNeuron", "Astro", "Micro", "Oligo", "OPC", "Endo", "Mural"],
                {**MAJOR, "Cdh9+": SUBTYPE["Cdh9+"], "OtherExNeuron": MAJOR["ExNeuron"]},
                "Fig3C_merged_cellchat_Cdh9_diff_count_weight", (150 * MM, 140 * MM))
tot = pd.DataFrame({c: {m: M[(c, m)].values.sum() for m in ["count", "weight"]} for c in GROUP_ORDER}).T
tot.to_csv(FIGS / "Fig3C_cellchat_total_count_weight.csv")

# ---- 3D/E PROGENy heatmaps ---------------------------------------------------
act = pd.read_csv(DATA / "progeny_activity_long.csv")
dif = pd.read_csv(DATA / "progeny_differential.csv")
samples = [s for g in ["CON", "CES", "CPS"] for s in SAMPLE_ORDER if s.startswith(g)]


def progeny_panel(feats, name, title):
    mat = act.pivot_table(index="feature", columns="sample", values="activity").reindex(index=feats, columns=samples)
    z = mat.sub(mat.mean(1), axis=0).div(mat.std(1, ddof=1), axis=0)
    dd = dif.set_index("feature").reindex(feats)
    h = max(3.2 * len(feats) * MM + 25 * MM, 45 * MM)
    fig = plt.figure(figsize=(130 * MM, h))
    gs = fig.add_gridspec(2, 6, height_ratios=[0.05, 1], width_ratios=[9, 2.4, 0.5, 0.3, 1.6, 0.3], wspace=0.12, hspace=0.03)
    axg = fig.add_subplot(gs[0, 0]); ax = fig.add_subplot(gs[1, 0]); axd = fig.add_subplot(gs[1, 1])
    cax1 = fig.add_subplot(gs[1, 3]); cax2 = fig.add_subplot(gs[1, 5])
    axg.imshow([[list(plt.matplotlib.colors.to_rgb(GROUP[s.split("_")[0]])) for s in samples]], aspect="auto")
    axg.set_xticks([]); axg.set_yticks([]); [sp.set_visible(False) for sp in axg.spines.values()]
    for gname, xc in [("CON", 1), ("CES", 4), ("CPS", 7)]:
        axg.text(xc, -0.9, gname, ha="center", va="bottom", fontsize=6, color=GROUP[gname])
    im = ax.imshow(z.values, cmap=CMAP_DIV, norm=TwoSlopeNorm(0, -2, 2), aspect="auto", interpolation="nearest")
    ax.set_yticks(range(len(feats)), [("↑ " if dd.loc[f, "direction"] == "up_in_both" else "↓ ") + f for f in feats], fontsize=5.5)
    ax.set_xticks(range(len(samples)), samples, rotation=90, fontsize=5)
    for x in [2.5, 5.5]:
        ax.axvline(x, color="k", lw=0.6)
    nup = (dd.direction == "up_in_both").sum()
    if 0 < nup < len(feats):
        ax.axhline(nup - 0.5, color="k", lw=0.8); axd.axhline(nup - 0.5, color="k", lw=0.8)
    cb = fig.colorbar(im, cax=cax1); cb.set_label("Row z-score (heatmap)", fontsize=5.5); cb.outline.set_visible(False)
    cax1.yaxis.set_label_position("left")
    B = dd[["beta_CES", "beta_CPS"]].values; m = np.nanmax(np.abs(B))
    im2 = axd.imshow(B, cmap=CMAP_DIV, norm=TwoSlopeNorm(0, -m, m), aspect="auto")
    for i in range(len(feats)):
        for j, c in enumerate(["CES", "CPS"]):
            axd.text(j, i, f"{B[i, j]:.2f}{sig_symbol(dd.iloc[i]['fdr_' + c])}", ha="center", va="center", fontsize=4.8)
    axd.set_xticks([0, 1], ["CES\n−CON", "CPS\n−CON"], fontsize=5.5); axd.xaxis.tick_top()
    axd.set_yticks([])
    cb2 = fig.colorbar(im2, cax=cax2); cb2.set_label("Mean activity difference (right panel)", fontsize=5.5); cb2.outline.set_visible(False)
    fig.suptitle(title, fontsize=7, y=1.0 + 0.25 / (h / MM / 10))
    ax.set_xlabel("Welch t-test per feature (n = 3 vs 3 animals), BH-FDR: · <0.1, * <0.05, ** <0.01, *** <0.001", fontsize=5)
    src = z.add_prefix("z_").join(mat.add_prefix("activity_")).join(dd[["beta_CES", "beta_CPS", "pvalue_CES", "pvalue_CPS",
                                                                       "fdr_CES", "fdr_CPS", "direction"]])
    save(fig, name, src.reset_index(), index=False)


TOPN = 20
conc = dif[dif.direction != "discordant"]
feats = list(conc[conc.direction == "up_in_both"].head(TOPN).feature) + list(conc[conc.direction == "down_in_both"].head(TOPN).feature)
progeny_panel(feats, "Fig3D_progeny_common_all_celltypes",
              f"PROGENy: top {TOPN} common increases and decreases (CES and CPS vs CON), all cell types")
c9 = conc[conc.group == "Cdh9+"]
feats9 = list(c9[c9.direction == "up_in_both"].feature) + list(c9[c9.direction == "down_in_both"].feature)
progeny_panel(feats9, "Fig3E_progeny_Cdh9", "Cdh9+: PROGENy pathways changed in the same direction in CES and CPS")

# ---- 3F violins of selected genes in Cdh9+ -----------------------------------
GENES = ["Camk1d", "Grik3", "Iqsec1", "Shank2"]
a = ad.read_h5ad(H5AD / "ACC_all_cells.h5ad", backed="r")
m = (a.obs["subtype"] == "Cdh9+").values
sub = a[m, GENES].to_memory()
X = sub.X.toarray() if hasattr(sub.X, "toarray") else np.asarray(sub.X)
ob = sub.obs[["sample", "Experiment"]].astype(str)
de = pd.read_csv(DATA / "DE_celllevel_subtype.csv.gz"); de = de[(de.cluster == "Cdh9+") & (de.reference == "CON") & de.gene.isin(GENES)]
stars = lambda p: "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else "ns"
fig, axs = plt.subplots(1, 4, figsize=(150 * MM, 42 * MM), gridspec_kw=dict(wspace=0.45))
rng = np.random.default_rng(SEED); rows = []
for ax, (k, g) in zip(axs, enumerate(GENES)):
    for i, cond in enumerate(GROUP_ORDER):
        v = X[(ob.Experiment == cond).values, k]
        parts = ax.violinplot(v, positions=[i], widths=0.8, showextrema=False)
        for b in parts["bodies"]:
            b.set_facecolor(GROUP[cond]); b.set_edgecolor("0.3"); b.set_linewidth(0.4); b.set_alpha(0.75)
        sm = pd.Series(X[:, k]).groupby(ob["sample"].values).mean()
        ss = [s for s in SAMPLE_ORDER if s.startswith(cond)]
        ax.scatter(i + rng.uniform(-0.12, 0.12, len(ss)), sm.reindex(ss), s=10, c="white", edgecolors="k", lw=0.6, zorder=3)
        rows += [dict(gene=g, sample=s, condition=cond, mean_logexpr=sm[s], n_cells=int((ob["sample"] == s).sum())) for s in ss]
    top = X[:, k].max()
    for j, cond in enumerate(["CES", "CPS"]):                # brackets CON-CES and CON-CPS
        r1 = de[(de.gene == g) & (de.condition == cond)]
        p = r1.pvals_adj.iloc[0] if len(r1) else np.nan
        yb = top * (1.06 + 0.13 * j)
        ax.plot([0, 0, j + 1, j + 1], [yb - top * 0.02, yb, yb, yb - top * 0.02], color="0.2", lw=0.6)
        ax.text((j + 1) / 2, yb, stars(p) if np.isfinite(p) else "n.t.", ha="center", va="bottom", fontsize=6 if np.isfinite(p) and p < 0.05 else 5)
        rows.append(dict(gene=g, sample=f"{cond} vs CON", condition=cond, log2FC=r1.log2FC.iloc[0] if len(r1) else np.nan, adj_P=p, label=stars(p)))
    ax.set_ylim(0, top * 1.32)
    ax.set_xticks(range(3), GROUP_ORDER); ax.set_title(g, style="italic", fontsize=7)
axs[0].set_ylabel("Log-normalised expression")
fig.text(0.5, 1.03, "Cdh9+ ExNeurons; violins = nuclei, white dots = animal means. Cell-level Welch t-test vs CON, BH: * < 0.05, ** < 0.01, *** < 0.001",
         ha="center", fontsize=5.5)
save(fig, "Fig3F_Cdh9_gene_violins", pd.DataFrame(rows))
print("Fig3 done")
