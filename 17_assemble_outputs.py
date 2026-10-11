#!/usr/bin/env python3
"""Step 17 - final layout: figures/ (figureN panels + source data), tables/ (numbered Excel tables + full DE
results), data/processed/ (cell metadata, CellChat objects, md5 checksums). Specification: manifest.py."""
import hashlib, shutil
import numpy as np, pandas as pd, anndata as ad
from config import *
from manifest import FIGURES, EXTRA_SOURCE, TABLES

LFC, PADJ = 0.25, 0.05
DIR = {"FIGS": FIGS, "QC": QC}
def md5(p, chunk=1 << 24):
    h = hashlib.md5()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(chunk), b""):
            h.update(b)
    return h.hexdigest()

# ---- figures --------------------------------------------------------------------------------------
fig_rows = []
for panel, name, d, stem, scripts, inputs, legend in FIGURES:
    src = DIR[d]
    for ext in (".pdf", ".png"):
        shutil.copy2(src / f"{stem}{ext}", FIG_OUT / f"{name}{ext}")
    srcs = [(p, p.stem.split("_source_data", 1)[1]) for p in sorted(src.glob(f"{stem}_source_data*.csv"))]
    for e, suf in EXTRA_SOURCE.get(panel, []):
        srcs.append((QC / e[3:] if e.startswith("QC:") else FIGS / e, f"_{suf}" if suf else ""))
    for s, suffix in srcs:
        shutil.copy2(s, FIG_OUT / f"{name}_source_data{suffix}.csv")
    fig_rows.append(dict(panel=panel, name=name, n_source=len(srcs)))
print(pd.DataFrame(fig_rows).to_string(index=False))

# ---- cell metadata (all QC-passing nuclei with final status) --------------------------------------
R = lambda f: ad.read_h5ad(PIPE / f, backed="r").obs
post = R("all_cell_adata_combined_post_filter.h5ad")
hi_ids = set(ad.read_h5ad(PIPE / "high_feature_cell_adata_combined.h5ad", backed="r").obs_names)
hi, lo = R("adata_high_feature_cluster.h5ad"), R("adata_low_feature_cluster.h5ad")
lin = pd.read_csv(PIPE / "lineage_check_removed.csv", index_col=0)
ex1, ex2, inn = R("adata_ExNeuron_cluster.h5ad"), R("adata_ExNeuron_cluster2.h5ad"), R("adata_InNeuron_cluster.h5ad")
allc = ad.read_h5ad(H5AD / "ACC_all_cells.h5ad", backed="r"); neu = ad.read_h5ad(H5AD / "ACC_neurons.h5ad", backed="r")
c = pd.DataFrame(index=post.index)
c["barcode"] = [i.split("-")[0] for i in c.index]
c["raw_SampleID"] = post.SampleID.astype(str)
c["sample_display"] = c.raw_SampleID.map({k: v[1] for k, v in SAMPLES.items()})
c["sample_pipeline"] = c.raw_SampleID.map({k: v[0] for k, v in SAMPLES.items()})
c["group"] = c.raw_SampleID.map({k: v[2] for k, v in SAMPLES.items()}); c["pooled_group"] = np.where(c.group == "CON", "CON", "SDS")
for k in ["n_genes_by_counts", "total_counts", "pct_counts_mt", "pct_counts_ribo", "doublet_score"]:
    c[k] = post[k].values
c["split_object"] = np.where(c.index.isin(hi_ids), "high_feature", "low_feature")
sp = pd.concat([hi[["leiden", "major_celltype"]], lo[["leiden", "major_celltype"]]])
c["leiden_split_object"] = sp.leiden.astype(str).reindex(c.index)
c["major_celltype_split_object"] = sp.major_celltype.astype(str).reindex(c.index)
c.loc[lin.index, "leiden_split_object"] = lin.leiden.astype(str); c.loc[lin.index, "major_celltype_split_object"] = lin.major_celltype
c["leiden_ExNeuron_pass1"] = ex1.leiden.astype(str).reindex(c.index)
c["leiden_ExNeuron_pass2"] = ex2.leiden.astype(str).reindex(c.index)
c["leiden_ExNeuron_sub"] = ex2.leiden_sub.astype(str).reindex(c.index)
c["leiden_InNeuron"] = inn.leiden.astype(str).reindex(c.index)
c["major_celltype"] = allc.obs.major_celltype.astype(str).reindex(c.index)
c["subtype"] = allc.obs.subtype.astype(str).reindex(c.index)
st = pd.Series("final", index=c.index)
st[c.index.isin(lin.index)] = "excluded: lineage check"
mixed = c.major_celltype_split_object.str.contains("_", na=False) & ~c.index.isin(lin.index)
st[mixed] = "excluded: mixed cluster (" + c.loc[mixed, "major_celltype_split_object"] + ")"
st[(c.leiden_ExNeuron_pass1 == "11") & c.major_celltype.isna()] = "excluded: ExNeuron pass-1 cluster 11 (low quality)"
in_all = set(hi.index[hi.major_celltype == "InNeuron"]) | set(lo.index[lo.major_celltype == "InNeuron"])
st[c.index.isin(in_all) & c.major_celltype.isna()] = "excluded: InNeuron cluster 7 (mixed Ex/In)"
c["status"] = st
assert (c.status == "final").sum() == allc.n_obs and c.loc[c.status == "final", "major_celltype"].notna().all(), c.status.value_counts()
c = c.join(pd.DataFrame(allc.obsm["X_umap"], index=allc.obs_names, columns=["UMAP1_allcells", "UMAP2_allcells"]))
c = c.join(pd.DataFrame(neu.obsm["X_umap"], index=neu.obs_names, columns=["UMAP1_neurons", "UMAP2_neurons"]))
c.index.name = "cell_id"
c.to_csv(PROC / "cell_metadata.csv.gz")
print(c.status.value_counts().to_string())

# ---- CellChat objects -------------------------------------------------------------------------------
for cfg in ("fine", "coarse"):
    for f in sorted((DATA / f"cellchat_{cfg}").glob("*.rds")):
        shutil.copy2(f, PROC / "cellchat" / f"cellchat_{cfg}_{f.name}")

# ---- tables -------------------------------------------------------------------------------------------
thr = pd.read_csv(PIPE / "sample_qc_thresholds.csv").set_index("SampleName")
pre = ad.read_h5ad(PREFILTER, backed="r").obs[["SampleID", "n_genes_by_counts", "total_counts", "pct_counts_mt"]].copy()
pre["SampleID"] = pre.SampleID.astype(str)
med = lambda o: o.groupby("SampleID")[["total_counts", "n_genes_by_counts", "pct_counts_mt"]].median()
qm = med(pre).add_prefix("median_").add_suffix("_before_qc").join(med(c.rename(columns={"raw_SampleID": "SampleID"})).add_prefix("median_").add_suffix("_after_qc"))
disp = {k: v[1] for k, v in SAMPLES.items()}
order = [v[1] for v in SAMPLES.values()]
fin = c[c.status == "final"]
X = {}
s_ = pd.DataFrame({"sample_display": [v[1] for v in SAMPLES.values()], "raw_SampleID": list(SAMPLES), "sample_pipeline": [v[0] for v in SAMPLES.values()],
                   "group": [v[2] for v in SAMPLES.values()]}).set_index("raw_SampleID")
s_ = s_.join(thr[["n_before", "n_predicted_doublets", "n_after", "glial_peak", "neuronal_peak", "threshold_used", "threshold_rule", "n_low", "n_high"]])
s_["scrublet_doublet_pct"] = (100 * s_.n_predicted_doublets / s_.n_before).round(2)
s_["nuclei_final"] = fin.groupby("raw_SampleID").size()
s_ = s_.rename(columns={"n_before": "nuclei_in_matrix", "n_after": "nuclei_after_qc", "threshold_used": "split_threshold", "threshold_rule": "split_rule",
                        "n_low": "n_low_feature", "n_high": "n_high_feature"})
X["samples"] = s_.reset_index()
qm.insert(0, "sample_display", qm.index.map(disp)); X["qc_medians"] = qm.reset_index()
ex = c.groupby(["sample_display", "status"]).size().unstack(fill_value=0).reindex(order); ex.loc["all animals"] = ex.sum(); X["exclusions"] = ex.reset_index()
g = pre.groupby("SampleID")
cap = pd.DataFrame({"nuclei_pre_QC": g.size(), "median_genes": g.n_genes_by_counts.median(),
                    "n_gt_4000_genes": g.n_genes_by_counts.apply(lambda x: int((x > 4000).sum())), "n_gt_15000_UMIs": g.total_counts.apply(lambda x: int((x > 15000).sum()))})
cap["frac_gt_4000_genes"] = (cap.n_gt_4000_genes / cap.nuclei_pre_QC).round(3); cap.insert(0, "sample_display", cap.index.map(disp)); X["gene_cap"] = cap.reset_index()
mx = []
for t in ["Astro", "Oligo", "OPC", "Micro", "Mural"]:
    for s in order:
        n_mixed = int(((c.status == f"excluded: mixed cluster ({t}_ExNeuron)") & (c.sample_display == s)).sum())
        kept = int(((fin.major_celltype == t) & (fin.sample_display == s)).sum())
        mx.append(dict(celltype=t, sample=s, in_mixed_X_cluster=n_mixed, kept=kept, frac_in_mixed=n_mixed / max(n_mixed + kept, 1)))
X["mixed"] = pd.DataFrame(mx)
X["major_counts"] = fin.groupby(["sample_display", "major_celltype"]).size().unstack(fill_value=0).reindex(order).reset_index()
X["subtype_counts"] = fin[fin.major_celltype.isin(["ExNeuron", "InNeuron"])].groupby(["sample_display", "subtype"]).size().unstack(fill_value=0).reindex(order).reset_index()
summ, sig = [], {}
for lvl in ("major", "subtype"):
    cl = pd.read_csv(DATA / f"DE_celllevel_{lvl}.csv.gz"); pb = pd.read_csv(DATA / f"DE_pseudobulk_{lvl}.csv.gz")
    cl["sig"] = (cl.pvals_adj < PADJ) & (cl.log2FC.abs() > LFC)
    s2 = cl.groupby(["cluster", "condition", "reference"]).agg(n_nuclei_cond=("n_cells_cond", "first"), n_nuclei_ref=("n_cells_ref", "first"),
            genes_tested=("gene", "size"), DEG_up=("sig", lambda x: int((x & (cl.loc[x.index, "log2FC"] > 0)).sum())),
            DEG_down=("sig", lambda x: int((x & (cl.loc[x.index, "log2FC"] < 0)).sum()))).reset_index()
    p2 = pb.assign(s=pb["adj.P.Val"] < PADJ).groupby(["cluster", "condition"]).agg(pb_genes_tested=("gene", "size"), pb_genes_adjP_lt_0_05=("s", "sum"),
            pb_n_animals=("n_samples", "first")).reset_index()
    s2 = s2.merge(p2, on=["cluster", "condition"], how="left"); s2.insert(0, "level", lvl); summ.append(s2)
    sig[lvl] = cl[cl.sig & (cl.reference == "CON")].drop(columns="sig").merge(
        pb[["cluster", "condition", "gene", "logFC", "P.Value", "adj.P.Val"]].rename(columns={"logFC": "pb_logFC", "P.Value": "pb_P", "adj.P.Val": "pb_adjP"}),
        on=["cluster", "condition", "gene"], how="left")
    shutil.copy2(DATA / f"DE_celllevel_{lvl}.csv.gz", TAB_OUT / f"table3_full_celllevel_{lvl}.csv.gz")
    shutil.copy2(DATA / f"DE_pseudobulk_{lvl}.csv.gz", TAB_OUT / f"table3_full_pseudobulk_{lvl}.csv.gz")
X["deg_summary"] = pd.concat(summ); X["sig_major"] = sig["major"]; X["sig_subtype"] = sig["subtype"]
for cfg in ("fine", "coarse"):
    dd = DATA / f"cellchat_{cfg}"; nets = []
    for cond in ("CON", "CES", "CPS"):
        for kind in ("count", "weight"):
            m = pd.read_csv(dd / f"{cond}_net_{kind}.csv", index_col=0).stack().rename("value").reset_index()
            m.columns = ["source", "target", "value"]; m.insert(0, "measure", kind); m.insert(0, "condition", cond); nets.append(m)
    X[f"cc_net_{cfg}"] = pd.concat(nets)
    if cfg == "fine":
        X["cc_lr_fine"] = pd.concat([pd.read_csv(dd / f"{cond}_communications.csv").assign(condition=cond) for cond in ("CON", "CES", "CPS")]).query("source == 'Cdh9+' or target == 'Cdh9+'")
X["cc_balance"] = pd.concat([pd.read_csv(DATA / f"cellchat_{cfg}" / "balancing.csv").assign(grouping=cfg) for cfg in ("fine", "coarse")])
ls_ = pd.read_csv(PIPE / "lineage_check_removed.csv", index_col=0); ls_["animal"] = ls_.SampleID.astype(str).map(disp)
X["lineage"] = ls_.groupby(["object", "major_celltype", "animal"]).size().unstack(fill_value=0).reset_index()
PFX = {"D": DATA, "P": PIPE, "F": FIGS, "Q": QC, "W": WORK}
tab_rows = []
for tid, (title, unit, scripts, sheets) in TABLES.items():
    about = pd.DataFrame([dict(sheet=s, contents=n, source=src.split(":", 1)[1] if not src.startswith("X:") else "built in 17_assemble_outputs.py") for s, src, n in sheets])
    with pd.ExcelWriter(TAB_OUT / f"{tid}.xlsx", engine="openpyxl") as w:
        pd.DataFrame({"field": ["table", "title", "statistical unit / model", "scripts", "DEG threshold"],
                      "value": [tid, title, unit, scripts, f"adj. P < {PADJ}, |log2FC| > {LFC}"]}).to_excel(w, sheet_name="about", index=False)
        about.to_excel(w, sheet_name="about", index=False, startrow=7)
        for sheet, src, note in sheets:
            k, v = src.split(":", 1)
            df = X[v] if k == "X" else pd.read_csv(PFX[k] / v)
            assert len(df) < 1_048_000, (tid, sheet, len(df))
            df.to_excel(w, sheet_name=sheet[:31], index=False)
            tab_rows.append(dict(table=tid, sheet=sheet, rows=len(df)))
print(pd.DataFrame(tab_rows).groupby("table").rows.agg(["count", "sum"]).to_string())

# ---- checksums of processed data ------------------------------------------------------------------------
ck = [dict(file=str(p.relative_to(ROOT)), bytes=p.stat().st_size, md5=md5(p)) for p in sorted(PROC.rglob("*"))
      if p.is_file() and "intermediate" not in p.parts]
pd.DataFrame(ck).to_csv(PROC / "checksums_md5.csv", index=False)
print(pd.DataFrame(ck).to_string(index=False))
