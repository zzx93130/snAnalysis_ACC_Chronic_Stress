"""Step 16 - ExNeuron NMF meta-programs: annotation, scoring, group tests and Supplementary Fig. S2.

Annotation: each MP is named from the GO terms enriched in its top-30 genes (expressed-gene universe,
08c_nmf_go.R) and, in brackets, the ExNeuron subtype with the highest MP score. The names are curated
below (ANNOT) from data/nmf/metaprogram_GO_expressed_universe.csv; the terms they rest on are printed
in the annotation table. Terms that also pass with the stricter HVG universe are flagged.
Scores: scanpy score_genes with the 30 MP genes on log-normalised data, per nucleus.
Group test: per animal x subtype mean score, Welch t-test CES or CPS vs CON (n = 3 vs 3), BH within contrast.
Output: figures/FigS2_ExNeuron_NMF_metaprograms.{pdf,png} + source data; data/nmf/metaprogram_annotation.csv,
data/nmf/metaprogram_scores_per_animal.csv, data/nmf/metaprogram_group_test.csv.
"""
import sys
from pathlib import Path
import numpy as np, pandas as pd, scanpy as sc
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from scipy import stats
from statsmodels.stats.multitest import multipletests
sys.path.insert(0, str(Path(__file__).parent))
from common import *  # noqa
from config import *  # project paths

D = DATA / "nmf"
ANNOT = {  # MP: (functional name from GO of top-30 genes, terms used)
    "MP1": ("Glycosylation & K+ channels", "heparan sulfate sulfotransferase; GalNAc-transferase; K+ transmembrane transporter"),
    "MP2": ("Synapse / cell-junction assembly", "cell junction assembly; postsynaptic density membrane; synapse assembly"),
    "MP3": ("Neuronal differentiation TFs", "CNS neuron differentiation (Foxp2, Nfia, Nfib, Sox5)"),
    "MP4": ("Ribosome & oxidative phosphorylation", "cytosolic ribosome; structural constituent of ribosome; respiratory chain complex IV"),
    "MP5": ("Postsynaptic density & cadherin adhesion", "postsynaptic density membrane; homophilic cell adhesion; catenin complex"),
    "MP6": ("Ligand-gated cation channels & axon guidance", "intracellular ligand-gated cation channel; cell adhesion molecule binding; axon guidance"),
    "MP7": ("Ionotropic glutamate receptor complex", "ionotropic glutamate receptor complex; neurotransmitter receptor complex"),
    "MP8": ("Cadherin-mediated adhesion", "catenin complex; cadherin binding; cell-cell adhesion mediated by cadherin"),
}
mpg = pd.read_csv(D / "metaprogram_genes.csv"); MPS = list(mpg.columns)
go = pd.read_csv(D / "metaprogram_GO_expressed_universe.csv")
goh = pd.read_csv(D / "metaprogram_GO_hvg_universe.csv")
info = pd.read_csv(D / "robust_factors.csv")

ad = sc.read_h5ad(H5AD / "ACC_neurons.h5ad"); ex = ad[ad.obs.major_celltype == "ExNeuron"].copy(); del ad
for m in MPS:
    sc.tl.score_genes(ex, mpg[m].tolist(), score_name=m, random_state=SEED)
o = ex.obs; o["subtype"] = o["subtype"].astype(str); o["sample"] = o["sample"].astype(str)
sub_order = [s for s in SUBTYPE if s in set(o.subtype)]
sub_mean = o.groupby("subtype", observed=True)[MPS].mean().loc[sub_order]
z = (sub_mean - sub_mean.mean()) / sub_mean.std()
top_sub = z.idxmax()

ann = []
for m in MPS:
    g = go[go.MP == m].sort_values("p.adjust")
    ann.append(dict(MP=m, name=f"{ANNOT[m][0]} ({top_sub[m]})", function=ANNOT[m][0], top_subtype=top_sub[m],
                    n_animals=info.loc[info.MP == m, "sample"].nunique(), n_factors=int((info.MP == m).sum()),
                    go_terms_used=ANNOT[m][1], top_go_term=g.Description.iloc[0] if len(g) else "",
                    top_go_padj=g["p.adjust"].iloc[0] if len(g) else np.nan, n_go_terms_expressed_universe=len(g),
                    n_go_terms_hvg_universe=int((goh.MP == m).sum()), top30_genes=" ".join(mpg[m])))
ann = pd.DataFrame(ann); ann.to_csv(D / "metaprogram_annotation.csv", index=False)

ps = o.groupby(["subtype", "sample", "Experiment"], observed=True)[MPS].mean().reset_index()
ps.to_csv(D / "metaprogram_scores_per_animal.csv", index=False)
rows = []
for st in sub_order + ["All ExNeuron"]:
    g = ps[ps.subtype == st] if st != "All ExNeuron" else o.groupby(["sample", "Experiment"], observed=True)[MPS].mean().reset_index()
    for c in ["CES", "CPS"]:
        for m in MPS:
            x, y = g.loc[g.Experiment == c, m], g.loc[g.Experiment == "CON", m]
            if len(x) >= 2 and len(y) >= 2:
                t, p = stats.ttest_ind(x, y, equal_var=False)
                rows.append(dict(subtype=st, contrast=f"{c} vs CON", MP=m, diff=x.mean() - y.mean(), t=t, p=p, n=len(x) + len(y)))
gt = pd.DataFrame(rows)
gt["fdr"] = gt.groupby("contrast")["p"].transform(lambda p: multipletests(p, method="fdr_bh")[1])
gt.to_csv(D / "metaprogram_group_test.csv", index=False)

# ------------------------------------------------------------------ figure (Supplementary Fig. S2)
names = dict(zip(ann.MP, ann.name)); mpcol = dict(zip(MPS, plt.cm.tab10(np.arange(len(MPS)))))
fig = plt.figure(figsize=(183 * MM, 190 * MM))
# A: colour key (left) | factor heatmap | lines to MP function and genes (right)
gA = fig.add_gridspec(1, 3, top=0.9674, bottom=0.7332, left=0.075, right=0.985, width_ratios=[0.035, 1, 1.75], wspace=0.08)
gC = fig.add_gridspec(1, 1, top=0.6674, bottom=0.5363, left=0.255, right=0.9)                # B (+ term labels below)
gB = fig.add_gridspec(1, 2, top=0.2579, bottom=0.1079, left=0.075, right=0.87, width_ratios=[1, 1.9], wspace=0.3)   # C, D

axA = fig.add_subplot(gA[0, 1])
ov = pd.read_csv(D / "factor_overlap.csv", index_col=0)
fi = info.sort_values("leaf_order"); M = ov.loc[fi.factor, fi.factor].to_numpy()
im = axA.imshow(M, cmap="Reds", vmin=0, vmax=30, interpolation="nearest", aspect="auto")
axA.set_xticks([]); axA.set_yticks([])
n = len(fi)
for i, (mp, s_) in enumerate(zip(fi.MP, fi["sample"])):
    axA.add_patch(plt.Rectangle((-3.4, i - 0.5), 1.5, 1, color=mpcol.get(mp, "0.85"), clip_on=False))
    axA.add_patch(plt.Rectangle((-1.7, i - 0.5), 1.1, 1, color=GROUP[s_.split("_")[0]], clip_on=False))
axA.text(-2.65, n + 0.2, "MP", fontsize=4.4, ha="center", va="top", rotation=90); axA.text(-1.15, n + 0.2, "group", fontsize=4.4, ha="center", va="top", rotation=90)
axA.set_title(f"A  {n} robust NMF factors -> {len(MPS)} meta-programs: function (GO of top 30 genes) and top genes", fontsize=6.5, loc="left", x=-0.1)
cax = fig.add_subplot(gA[0, 0]); cb = fig.colorbar(im, cax=cax); cb.outline.set_visible(False)
cax.yaxis.set_ticks_position("left"); cax.yaxis.set_label_position("left"); cb.set_label("Shared top-50 genes", fontsize=5.5); cb.ax.tick_params(labelsize=5)
pos = cax.get_position(); cax.set_position([pos.x0 - 0.035, pos.y0 + pos.height * 0.25, pos.width * 0.6, pos.height * 0.5])
axT = fig.add_subplot(gA[0, 2]); axT.set_xlim(0, 1); axT.set_ylim(n - 0.5, -0.5); axT.axis("off")
ty = np.linspace(-0.5 + n / len(MPS) / 2, n - 0.5 - n / len(MPS) / 2, len(MPS))
for k, r in ann.reset_index(drop=True).iterrows():
    idx = np.where(fi.MP.to_numpy() == r.MP)[0]
    yb, yt = idx.mean(), ty[k]
    axT.plot([0.0, 0.025, 0.065, 0.085], [yb, yb, yt, yt], color=mpcol[r.MP], lw=0.8, solid_capstyle="round", clip_on=False)
    axT.plot([0.0, 0.0], [idx.min() - 0.45, idx.max() + 0.45], color=mpcol[r.MP], lw=1.6, clip_on=False)
    axT.text(0.095, yt - 0.25, r.MP, fontsize=6, weight="bold", color=mpcol[r.MP], va="bottom")
    axT.text(0.165, yt - 0.25, f"{r.function}  [{r.top_subtype}; {r.n_animals} animals]", fontsize=5.6, va="bottom")
    axT.text(0.165, yt + 0.25, ", ".join(r.top30_genes.split()[:10]), fontsize=4.8, style="italic", color="0.3", va="top")

# B GO dot plot: top 3 terms per MP (compact)
axC = fig.add_subplot(gC[0, 0])
sel = go.sort_values("p.adjust").groupby("MP").head(3)
terms = []
for m in MPS:
    terms += [t for t in sel[sel.MP == m].sort_values("p.adjust").Description if t not in terms]
xs = {m: i for i, m in enumerate(MPS)}; ys = {t: j for j, t in enumerate(terms)}
gsub = go[go.Description.isin(terms)]
size = lambda c: 4 + 6 * np.asarray(c)
sca = axC.scatter([ys[t] for t in gsub.Description], [xs[m] for m in gsub.MP], s=size(gsub.Count),
                  c=-np.log10(gsub["p.adjust"]), cmap="viridis", vmin=1.3, vmax=6, lw=0)
hv = set(zip(goh.MP, goh.Description))
for t, m, c_ in zip(gsub.Description, gsub.MP, gsub.Count):
    if (m, t) in hv:
        axC.scatter(ys[t], xs[m], s=size(c_), facecolor="none", edgecolor="k", lw=0.5)
axC.set_yticks(range(len(MPS)), [f"{m}  {ANNOT[m][0]}" for m in MPS], fontsize=5.2)
axC.set_xticks(range(len(terms)), terms, rotation=50, ha="right", rotation_mode="anchor", fontsize=4.8)   # full term names
axC.invert_yaxis(); axC.set_xlim(-0.6, len(terms) - 0.4); axC.set_ylim(len(MPS) - 0.5, -0.5)
for lab, m in zip(axC.get_yticklabels(), MPS):
    lab.set_color(mpcol[m])
cb = fig.colorbar(sca, ax=axC, fraction=0.015, pad=0.01); cb.set_label("-log10 adj. P", fontsize=5.2); cb.outline.set_visible(False); cb.ax.tick_params(labelsize=4.8)
for c in [3, 6, 9]:
    axC.scatter([], [], s=size(c), c="0.5", label=f"{c} genes")
axC.legend(loc="upper left", bbox_to_anchor=(1.07, 1.0), fontsize=4.8, handletextpad=0.2, labelspacing=0.6)
axC.set_title("B  GO enrichment of MP top-30 genes, top 3 terms per MP (expressed-gene universe; ring = also FDR < 0.05 with HVG universe)", fontsize=6.5, loc="left", x=-0.28)

# C MP score per subtype (z-score across subtypes)
axD = fig.add_subplot(gB[0, 0])
im = axD.imshow(z.T.to_numpy(), cmap="RdBu_r", norm=TwoSlopeNorm(0, -2.5, 2.5), aspect="auto", interpolation="nearest")
axD.set_xticks(range(len(sub_order)), sub_order, rotation=60, ha="right", fontsize=5.0)
axD.set_yticks(range(len(MPS)), MPS, fontsize=5.2)
for lab, m in zip(axD.get_yticklabels(), MPS):
    lab.set_color(mpcol[m])
for lab in axD.get_xticklabels():
    if lab.get_text() == "Cdh9+":
        lab.set_weight("bold")
axD.tick_params(length=0)
cb = fig.colorbar(im, ax=axD, fraction=0.05, pad=0.03); cb.set_label("Mean MP score\n(z across subtypes)", fontsize=5.2); cb.outline.set_visible(False); cb.ax.tick_params(labelsize=4.8)
axD.set_title("C  MP activity per ExNeuron subtype", fontsize=6.5, loc="left")

# D CES / CPS vs CON per subtype (wider)
axE = fig.add_subplot(gB[0, 1])
cols = [(st, c) for st in sub_order + ["All ExNeuron"] for c in ["CES vs CON", "CPS vs CON"]]
lim = np.nanpercentile(np.abs(gt["diff"]), 98)
psz = lambda p: 3 + 16 * np.clip(-np.log10(p), 0, 3)
for xi, (st, c) in enumerate(cols):
    g = gt[(gt.subtype == st) & (gt.contrast == c)].set_index("MP").reindex(MPS)
    axE.scatter([xi] * len(MPS), range(len(MPS)), s=psz(g.p), c=g["diff"], cmap="RdBu_r", norm=TwoSlopeNorm(0, -lim, lim), lw=0)
    for yi, (p, q) in enumerate(zip(g.p, g.fdr)):
        if q < 0.1:
            axE.scatter(xi, yi, s=psz(p), facecolor="none", edgecolor="k", lw=0.6)
axE.set_xticks([i * 2 + 0.5 for i in range(len(sub_order) + 1)], sub_order + ["All ExNeuron"], rotation=60, ha="right", fontsize=5.0)
for lab in axE.get_xticklabels():
    if lab.get_text() == "Cdh9+":
        lab.set_weight("bold")
for i in range(len(sub_order) + 1):
    axE.axvspan(i * 2 - 0.5, i * 2 + 1.5, color="0.95" if i % 2 else "white", lw=0, zorder=0)
axE.set_yticks(range(len(MPS)), MPS, fontsize=5.2); axE.set_xlim(-0.6, len(cols) - 0.4); axE.set_ylim(len(MPS) - 0.5, -0.5)
for lab, m in zip(axE.get_yticklabels(), MPS):
    lab.set_color(mpcol[m])
axE.tick_params(length=0)
for sp in axE.spines.values():
    sp.set_visible(False)
sm = plt.cm.ScalarMappable(cmap="RdBu_r", norm=TwoSlopeNorm(0, -lim, lim))
cb = fig.colorbar(sm, ax=axE, fraction=0.03, pad=0.02); cb.set_label("Score difference vs CON\n(animal means)", fontsize=5.2); cb.outline.set_visible(False); cb.ax.tick_params(labelsize=4.8)
for pv in [0.05, 0.01, 0.001]:
    axE.scatter([], [], s=psz(pv), c="0.5", label=f"P = {pv:g}")
axE.legend(loc="upper left", bbox_to_anchor=(1.2, 1.0), fontsize=4.8, labelspacing=0.7, handletextpad=0.2)
axE.set_title("D  Stress vs CON per subtype (left CES, right CPS; ring = FDR < 0.1)", fontsize=6.5, loc="left")

save(fig, "FigS2_ExNeuron_NMF_metaprograms",
     {"annotation": ann, "GO_top3": sel, "subtype_score_z": z.reset_index(names="subtype"), "group_test": gt})
print(ann[["MP", "name", "n_animals", "top_go_term", "top_go_padj", "n_go_terms_expressed_universe", "n_go_terms_hvg_universe"]].to_string(index=False))
print(gt[gt.fdr < 0.1].sort_values("p").to_string(index=False))
print(gt[gt.subtype == "Cdh9+"].pivot(index="MP", columns="contrast", values=["diff", "p"]).round(4).to_string())
