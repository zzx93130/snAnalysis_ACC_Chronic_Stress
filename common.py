"""Shared style, palettes and export helpers for the ACC CES/CPS paper figures.

Every figure script imports this module so fonts, colours and file naming stay
identical across figures. Palettes are also written to WORK/palettes.json (read by the R scripts).
"""
from pathlib import Path
import json
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

from config import *   # ROOT, DATA (result tables), FIGS (figure files), H5AD, PIPE, WORK ...

SEED = 0

mpl.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans"],
    "font.size": 7, "axes.titlesize": 8, "axes.labelsize": 7,
    "xtick.labelsize": 6, "ytick.labelsize": 6, "legend.fontsize": 6,
    "legend.title_fontsize": 7, "legend.frameon": False,
    "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
    "axes.spines.top": False, "axes.spines.right": False, "axes.linewidth": 0.6,
    "xtick.major.width": 0.6, "ytick.major.width": 0.6,
    "xtick.major.size": 2.5, "ytick.major.size": 2.5,
    "savefig.dpi": 300, "savefig.bbox": "tight", "savefig.pad_inches": 0.02,
    "figure.dpi": 150,
})
MM = 1 / 25.4  # inches per mm

# ---- palettes (one per categorical variable, reused everywhere) ----------
GROUP = {"CON": "#8C8C8C", "CES": "#1B9E77", "CPS": "#D95F02", "SDS": "#7570B3"}
GROUP_ORDER = ["CON", "CES", "CPS"]
# major cell-type palette supplied by the user (celltype_colors); Endo is not in it -> muted teal added
MAJOR = {"ExNeuron": "#A3C6DB", "InNeuron": "#BF7A7A", "Astro": "#3E74AB",
         "Oligo": "#A3B746", "OPC": "#C6B745", "Micro": "#D2AF90",
         "Endo": "#5E8C7A", "Mural": "#C3BFBF"}
MAJOR_ORDER = list(MAJOR)
# rerun neuron subtype labels (reproducilibility.py); '+' appended to single-gene names
EX_RENAME = {"Foxp2": "Foxp2+", "Cdh9": "Cdh9+", "Ntng1": "Ntng1+", "Kcnip1": "Kcnip1+",
             "Arhgap6": "Arhgap6+"}
IN_RENAME = {"Vip": "Vip+", "Pvalb": "Pvalb+", "Sst": "Sst+", "Lamp5": "Lamp5+", "Cnr1": "Cnr1+"}
SUBTYPE = {
    # excitatory
    "Cdh9+": "#6A3D9A", "Foxp2+": "#7FB2D9", "Cux2+Col19a1": "#2F5F9E", "Cux2+Rorb": "#8FA9C4",
    "Cux2+Trpc6": "#1F78B4", "Cux2+Fth1": "#B7B7DD", "Cux2+Htr2c": "#9E7BB5",
    "Ptprm+Tshz2": "#9DD3E8", "Ptprm+Nell1": "#4FA3A0", "Ntng1+": "#2E7D4F",
    "Kcnip1+": "#86A83C", "Arhgap6+": "#4B4B7F",
    # inhibitory
    "Pvalb+": "#D9534F", "Sst+": "#F28E2B", "Vip+": "#E6C34A", "Lamp5+": "#C77CB0", "Cnr1+": "#9C6644",
}
EX_TYPES = [k for k in SUBTYPE if k not in ("Pvalb+", "Sst+", "Vip+", "Lamp5+", "Cnr1+")]
IN_TYPES = ["Pvalb+", "Sst+", "Vip+", "Lamp5+", "Cnr1+"]
UP, DOWN, NS = "#C0392B", "#2C6FB7", "#BDBDBD"
CMAP_DIV = "RdBu_r"      # centred quantities only
CMAP_SEQ = "viridis"     # non-negative quantities

SAMPLE_RENAME = {  # SampleID (raw folder) -> display name, from figure1_major_celltype_umap.py
    "24090566_AC_Control_5": "CON_1", "24090566_AC_Control_2": "CON_2", "24090566_AC_Control_3": "CON_3",
    "24090566_AC_ES_948": "CES_1", "24090566_AC_ES_938": "CES_2", "24090566_AC_ES_940": "CES_3",
    "24090566_AC_PS_966": "CPS_1", "24090566_AC_PS_968": "CPS_2", "24090566_AC_PS_973": "CPS_3"}
SAMPLE_ORDER = ["CON_1", "CON_2", "CON_3", "CES_1", "CES_2", "CES_3", "CPS_1", "CPS_2", "CPS_3"]
EXP_RENAME = {"Control": "CON", "ES": "CES", "PS": "CPS", "CON": "CON", "CES": "CES", "CPS": "CPS"}


def save(fig, name, source=None, index=False):
    """Write FIGS/name.pdf + .png (+ name_source_data.csv)."""
    fig.savefig(FIGS / f"{name}.pdf")
    fig.savefig(FIGS / f"{name}.png", dpi=300)
    if source is not None:
        if isinstance(source, dict):
            for k, v in source.items():
                v.to_csv(FIGS / f"{name}_source_data_{k}.csv", index=index)
        else:
            source.to_csv(FIGS / f"{name}_source_data.csv", index=index)
    plt.close(fig)


def write_palettes():
    json.dump({"group": GROUP, "major_celltype": MAJOR, "neuron_subtype": SUBTYPE,
               "direction": {"up": UP, "down": DOWN, "ns": NS}},
              open(WORK / "palettes.json", "w"), indent=2)


if __name__ == "__main__":
    write_palettes()
    print("palettes written", WORK / "palettes.json")


# ====================================================================================================
# Plotting helpers (UMAP, Venn, UpSet, circle network, statistics)
# ====================================================================================================
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyArrowPatch
from matplotlib.lines import Line2D
from scipy import stats
from statsmodels.stats.multitest import multipletests


def umap_scatter(ax, xy, colors, s=0.6, alpha=0.9, order_seed=SEED):
    order = np.random.default_rng(order_seed).permutation(len(xy))
    ax.scatter(xy[order, 0], xy[order, 1], c=np.asarray(colors)[order], s=s, lw=0, alpha=alpha, rasterized=True)
    ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_visible(False)


def umap_axes_arrows(ax, size=0.14):
    for end, lab, pos, rot in [((size, 0.02), "UMAP1", (size / 2 + 0.02, -0.03), 0),
                               ((0.02, size), "UMAP2", (-0.03, size / 2 + 0.02), 90)]:
        ax.annotate("", xy=end, xytext=(0.02, 0.02), xycoords="axes fraction",
                    arrowprops=dict(arrowstyle="-|>", lw=0.6, color="0.2", mutation_scale=5))
        ax.text(*pos, lab, transform=ax.transAxes, rotation=rot, ha="center", va="center", fontsize=5.5)


def label_centroids(ax, xy, labels, names, fontsize=6):
    import matplotlib.patheffects as pe
    for n in names:
        p = xy[np.asarray(labels) == n]
        if len(p) == 0:
            continue
        c = p[np.argmin(((p - np.median(p, 0)) ** 2).sum(1))]
        ax.text(*c, n, fontsize=fontsize, ha="center", va="center", weight="bold",
                path_effects=[pe.withStroke(linewidth=1.8, foreground="white")])


def legend_dots(ax_or_fig, names, colors, title=None, ncol=1, **kw):
    h = [Line2D([], [], marker="o", ls="", ms=4, mfc=colors[n], mec="none", label=n) for n in names]
    return ax_or_fig.legend(handles=h, title=title, ncol=ncol, handletextpad=0.3, labelspacing=0.35,
                            columnspacing=0.8, **kw)


def venn2(ax, a, b, labels, colors, title=None):
    """Area-free two-set Venn (equal circles), counts and percentages printed."""
    a, b = set(a), set(b)
    n_a, n_b, n_ab = len(a - b), len(b - a), len(a & b)
    tot = max(n_a + n_b + n_ab, 1)
    for x, c in [(-0.5, colors[0]), (0.5, colors[1])]:
        ax.add_patch(Circle((x, 0), 1.0, fc=c, ec="0.2", lw=0.8, alpha=0.45))
    for x, n in [(-0.95, n_a), (0, n_ab), (0.95, n_b)]:
        ax.text(x, 0, f"{n}\n({100 * n / tot:.0f}%)", ha="center", va="center", fontsize=7)
    ax.text(-0.85, 1.12, labels[0], ha="center", va="bottom", fontsize=7, color="0.1")
    ax.text(0.85, 1.12, labels[1], ha="center", va="bottom", fontsize=7, color="0.1")
    ax.set_xlim(-1.65, 1.65); ax.set_ylim(-1.2, 1.45); ax.set_aspect("equal"); ax.axis("off")
    if title:
        ax.set_title(title, fontsize=8)
    return n_a, n_ab, n_b


def upset(fig, sets, order, color="#3B6FB6", min_size=1, max_bars=12, title=None, set_colors=None):
    """Minimal UpSet: exclusive intersections, sorted by size. Returns the intersection table."""
    allg = sorted(set().union(*sets.values()))
    memb = pd.DataFrame({k: [g in sets[k] for g in allg] for k in order}, index=allg)
    combos = memb.groupby(order).size().reset_index(name="n")
    combos = combos[combos[order].any(axis=1) & (combos.n >= min_size)].sort_values("n", ascending=False).head(max_bars)
    gs = fig.add_gridspec(2, 3, height_ratios=[1.3, 1], width_ratios=[0.55, 0.32, max(len(combos), 1) * 0.28],
                          hspace=0.05, wspace=0.02)
    axb = fig.add_subplot(gs[0, 2]); axm = fig.add_subplot(gs[1, 2], sharex=axb); axs = fig.add_subplot(gs[1, 0], sharey=axm)
    axn = fig.add_subplot(gs[1, 1], sharey=axm); axn.axis("off")
    x = np.arange(len(combos))
    axb.bar(x, combos.n, color=color, width=0.65)
    for xi, n in zip(x, combos.n):
        axb.text(xi, n, str(n), ha="center", va="bottom", fontsize=5)
    axb.set_ylabel("Intersection size"); axb.tick_params(axis="x", bottom=False, labelbottom=False)
    axb.spines["bottom"].set_visible(False)
    ys = np.arange(len(order))[::-1]
    for xi, (_, r) in zip(x, combos.iterrows()):
        on = [ys[i] for i, k in enumerate(order) if r[k]]
        axm.scatter([xi] * len(order), ys, s=14, c="0.85", lw=0)
        axm.scatter([xi] * len(on), on, s=14, c=color, lw=0)
        if len(on) > 1:
            axm.plot([xi, xi], [min(on), max(on)], c=color, lw=1)
    axm.set_yticks(ys, order); axm.tick_params(axis="y", left=False, labelleft=False)
    axm.tick_params(axis="x", bottom=False, labelbottom=False)
    for sp in axm.spines.values():
        sp.set_visible(False)
    sizes = [len(sets[k]) for k in order]
    axs.barh(ys, sizes, color=[set_colors[k] for k in order] if set_colors else "0.5", height=0.6)
    axs.invert_xaxis(); axs.tick_params(axis="y", left=False, labelleft=False); axs.set_xlabel("Set size")
    for y, k in zip(ys, order):
        axn.text(0.04, y, k, ha="left", va="center", fontsize=6)
    axs.spines["left"].set_visible(False); axs.spines["right"].set_visible(False)
    for y, s in zip(ys, sizes):
        axs.text(s, y, f"{s} ", ha="right", va="center", fontsize=5)
    axm.set_ylim(-0.6, len(order) - 0.4)
    if title:
        axb.set_title(title, fontsize=8)
    return combos


def circle_network(ax, mat, colors, order, title=None, max_lw=4, node_size=None, top=None,
                   pos_color="#C0392B", neg_color="#2C6FB7", self_loops=True):
    """CellChat-style circle plot of a (signed) sender x receiver matrix.
    Edge width proportional to |value|; red = increase, blue = decrease (diff matrices).
    `top` keeps edges whose |value| is in the top fraction (as netVisual_diffInteraction(top=))."""
    n = len(order)
    ang = np.pi / 2 - 2 * np.pi * np.arange(n) / n
    P = {k: (np.cos(a), np.sin(a)) for k, a in zip(order, ang)}
    m = mat.reindex(index=order, columns=order).fillna(0)
    vals = m.values.copy()
    if top is not None and np.abs(vals).max() > 0:
        thr = np.quantile(np.abs(vals[vals != 0]), 1 - top) if (vals != 0).sum() else 0
        vals[np.abs(vals) < thr] = 0
    vmax = np.abs(vals).max() if np.abs(vals).max() > 0 else 1
    kept = []
    for i, s in enumerate(order):
        for j, r in enumerate(order):
            v = vals[i, j]
            if v == 0:
                continue
            kept.append((s, r, v))
            col = pos_color if v > 0 else neg_color
            lw = 0.3 + max_lw * abs(v) / vmax
            if s == r:
                if not self_loops:
                    continue
                x, y = P[s]
                ax.add_patch(Circle((x * 1.13, y * 1.13), 0.09, fc="none", ec=col, lw=lw, alpha=0.8))
                continue
            (x1, y1), (x2, y2) = P[s], P[r]
            ax.add_patch(FancyArrowPatch((x1 * 0.93, y1 * 0.93), (x2 * 0.93, y2 * 0.93),
                                         connectionstyle="arc3,rad=0.18", arrowstyle="-|>",
                                         mutation_scale=4 + lw * 1.5, lw=lw, color=col, alpha=0.7,
                                         shrinkA=4, shrinkB=6))
    for k in order:
        x, y = P[k]
        ax.add_patch(Circle((x, y), 0.075 if node_size is None else node_size[k], fc=colors.get(k, "0.6"),
                            ec="white", lw=0.5, zorder=5))
        ax.text(x * 1.14 + (0.05 if x > 0.2 else -0.05 if x < -0.2 else 0), y * 1.16 + (0.04 if abs(x) <= 0.2 else 0), k,
                ha="left" if x > 0.2 else "right" if x < -0.2 else "center", va="center", fontsize=5.5)
    ax.set_xlim(-1.7, 1.7); ax.set_ylim(-1.45, 1.45); ax.set_aspect("equal"); ax.axis("off")
    if title:
        ax.set_title(title, fontsize=7)
    return pd.DataFrame(kept, columns=["source", "target", "value"])


def logit_prop_test(counts, groups, cond, ref="CON"):
    """Replicate-aware composition test (propeller-style): per-sample proportions, logit
    transform with +0.5 pseudocount, Welch t-test cond vs ref per cell type, BH across types."""
    tot = counts.sum(1)
    prop = (counts + 0.5).div(tot + 0.5 * counts.shape[1], axis=0)
    lg = np.log(prop / (1 - prop))
    a, b = lg[groups == ref], lg[groups == cond]
    rows = []
    for ct in counts.columns:
        t = stats.ttest_ind(b[ct], a[ct], equal_var=False)
        pa, pb = (counts[ct] / tot)[groups == ref], (counts[ct] / tot)[groups == cond]
        rows.append(dict(celltype=ct, condition=cond, reference=ref, mean_prop_ref=pa.mean(),
                         mean_prop_cond=pb.mean(), prop_ratio=pb.mean() / pa.mean() if pa.mean() > 0 else np.nan,
                         logit_diff=b[ct].mean() - a[ct].mean(), t=t.statistic, pvalue=t.pvalue))
    r = pd.DataFrame(rows)
    r["fdr"] = multipletests(r.pvalue.fillna(1), method="fdr_bh")[1]
    return r


def sig_symbol(p):
    return "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else "·" if p < 0.1 else ""


MIN_ROWS = 8   # a GO bar panel is laid out for at least this many rows


def go_overlap_panel(go, cluster, name, top_n=12, dirs=("up", "down")):
    """One row per DEG direction: Venn of CES/CPS GO terms + bars of shared terms.
    dirs may include "all" = GO of up- and down-regulated DEGs combined (one gene list)."""
    fig, axs = plt.subplots(len(dirs), 2, figsize=(150 * MM, 60 * MM * len(dirs)), gridspec_kw=dict(width_ratios=[1, 2.4]))
    srcs = []
    dlab = {"up": "up-regulated DEGs", "down": "down-regulated DEGs", "all": "all DEGs (up + down)"}
    for r, dirn in enumerate(dirs):
        g = go[(go.cluster == cluster) & (go.direction == dirn) & (go["p.adjust"] < 0.05)]
        g = g.assign(key=g.ID + " | " + g.Description)
        A = set(g[g.condition == "CES"].key); B = set(g[g.condition == "CPS"].key)
        venn2(axs[r, 0], A, B, ["CES vs CON", "CPS vs CON"], [GROUP["CES"], GROUP["CPS"]],
              title=f"{cluster}: GO terms, {dlab[dirn]}")
        shared = A & B
        ax = axs[r, 1]
        if not shared:
            ax.text(0.5, 0.5, "No shared significant GO terms", ha="center", va="center", transform=ax.transAxes)
            ax.axis("off"); continue
        s = g[g.key.isin(shared)]
        best = s.groupby("key")["p.adjust"].min().sort_values().head(top_n)
        s = s[s.key.isin(best.index)].copy()
        s["score"] = -np.log10(s["p.adjust"]) * np.where(s.condition == "CES", -1, 1)
        order = list(best.index[::-1])
        for cond in ["CES", "CPS"]:
            t = s[s.condition == cond].set_index("key").reindex(order)
            ax.barh(range(len(order)), t.score, color=GROUP[cond], height=0.7, label=f"{cond} vs CON")
        lab = [k.split(" | ", 1)[1] for k in order]
        lab = [l if len(l) < 48 else l[:46] + "…" for l in lab]
        ax.set_yticks(range(len(order)), lab, fontsize=5.5)
        pad = max(0.0, (MIN_ROWS - len(order)) / 2)          # few terms: keep bars thin, centred in the row
        ax.set_ylim(-0.6 - pad, len(order) - 0.4 + pad)
        m = np.abs(s.score).max() * 1.1
        ax.set_xlim(-m, m); ax.axvline(0, color="0.2", lw=0.6)
        ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda v, _: f"{abs(v):g}"))
        ax.set_xlabel("-log10(adj. P)"); ax.set_title(f"Shared terms ({len(shared)}), top {len(order)} by adj. P", fontsize=7)
        ax.legend(loc="lower right", fontsize=5)
        srcs.append(s.assign(direction=dirn))
    fig.tight_layout()
    save(fig, name, pd.concat(srcs) if srcs else pd.DataFrame())


def hub_network(ax, mat, colors, hub, partners, title=None, max_lw=3.5, vmax=None,
                pos_color="#C0392B", neg_color="#2C6FB7", total=None, alpha_floor=0.25):
    """Hub-centred CellChat diff plot: `hub` in the centre, partners on a circle.
    Only hub->partner, partner->hub and the hub self-loop are drawn. Edge width ~ |value|;
    red = increase, blue = decrease. Outgoing and incoming edges of the same pair bend to
    opposite sides (arc3 with the same rad on reversed direction), so both stay visible.
    `total` (partner -> absolute interactions with the hub) de-emphasises sparse partners: node area
    scales with total, and node, label and edge opacity = clip(total / median(total), alpha_floor, 1)."""
    n = len(partners)
    ang = np.pi / 2 - 2 * np.pi * (np.arange(n) + 0.5) / n
    P = {k: (np.cos(a), np.sin(a)) for k, a in zip(partners, ang)}
    m = mat.reindex(index=[hub] + partners, columns=[hub] + partners).fillna(0)
    vals = [(hub, p, m.loc[hub, p]) for p in partners] + [(p, hub, m.loc[p, hub]) for p in partners] + [(hub, hub, m.loc[hub, hub])]
    vmax = vmax or max(abs(v) for *_, v in vals) or 1
    if total is not None:
        t = pd.Series(total, dtype=float).reindex(partners).fillna(0)
        alpha = (t / t.median()).clip(alpha_floor, 1).to_dict()
        rad = (0.03 + 0.055 * np.sqrt(t / t.max())).to_dict()
    else:
        alpha, rad = {k: 1.0 for k in partners}, {k: 0.065 for k in partners}
    kept = []
    for s, r, v in vals:
        if v == 0:
            continue
        kept.append((s, r, v))
        col = pos_color if v > 0 else neg_color
        lw = 0.3 + max_lw * abs(v) / vmax
        if s == r:
            ax.add_patch(Circle((0, 0.2), 0.11, fc="none", ec=col, lw=lw, alpha=0.8, zorder=3))
            continue
        a, b = ((0, 0), P[r]) if s == hub else (P[s], (0, 0))
        ax.add_patch(FancyArrowPatch(a, b, connectionstyle="arc3,rad=0.22", arrowstyle="-|>",
                                     mutation_scale=4 + lw * 1.5, lw=lw, color=col, alpha=0.75 * alpha[r if s == hub else s],
                                     shrinkA=9 if s == hub else 5, shrinkB=5 if s == hub else 11, zorder=2))
    ax.add_patch(Circle((0, 0), 0.13, fc=colors.get(hub, "0.4"), ec="white", lw=0.8, zorder=6))
    ax.text(0, -0.21, hub, ha="center", va="top", fontsize=6.5, weight="bold", zorder=7,
            path_effects=[__import__("matplotlib.patheffects", fromlist=["withStroke"]).withStroke(linewidth=2, foreground="white")])
    for k in partners:
        x, y = P[k]
        ax.add_patch(Circle((x, y), rad[k], fc=colors.get(k, "0.6"), ec="white", lw=0.5, zorder=5, alpha=alpha[k]))
        ax.text(x * 1.12 + (0.03 if x > 0.05 else -0.03 if x < -0.05 else 0), y * 1.13, k,
                ha="left" if x > 0.05 else "right" if x < -0.05 else "center", va="center", fontsize=5.5,
                color="0.15" if alpha[k] >= 0.6 else "0.6")
    ax.set_xlim(-1.75, 1.75); ax.set_ylim(-1.35, 1.35); ax.set_aspect("equal"); ax.axis("off")
    if title:
        ax.set_title(title, fontsize=7)
    return pd.DataFrame(kept, columns=["source", "target", "value"])
