#!/usr/bin/env python3
"""Step 18 - Word documents in doc/: 1_methods, 2_code_figure_table_map, 3_data_and_objects, 4_qc_and_annotation,
5_figure_legends. Every number in the text is read from the outputs of steps 1-17. Needs python-docx."""
import re, ast
import numpy as np, pandas as pd, anndata as ad
from config import *
from manifest import FIGURES, EXTRA_SOURCE, TABLES

# ---- Word helpers (python-docx) ----------------------------------------------------------------------
from pathlib import Path
import hashlib, datetime
import pandas as pd
from docx import Document
from docx.shared import Pt, Cm
from docx.enum.section import WD_ORIENT

TODAY = datetime.date.today().isoformat()
NULL = "null"


def new_doc(title, subtitle=None, landscape=False):
    d = Document()
    st = d.styles["Normal"]; st.font.name = "Arial"; st.font.size = Pt(9.5)
    sec = d.sections[0]
    if landscape:
        sec.orientation = WD_ORIENT.LANDSCAPE; sec.page_width, sec.page_height = sec.page_height, sec.page_width
    for side in ("left_margin", "right_margin", "top_margin", "bottom_margin"):
        setattr(sec, side, Cm(1.8))
    d.add_heading(title, 0)
    p = d.add_paragraph(f"ACC snRNA-seq, chronic stress (CES / CPS vs CON). Document version {TODAY}.")
    if subtitle:
        d.add_paragraph(subtitle)
    return d


def h(d, text, level=1):
    d.add_heading(text, level)


def para(d, text, bold_prefix=None, italic=False):
    p = d.add_paragraph()
    if bold_prefix:
        r = p.add_run(bold_prefix); r.bold = True
    r = p.add_run(text); r.italic = italic
    return p


def bullets(d, items):
    for it in items:
        if isinstance(it, tuple):
            p = d.add_paragraph(style="List Bullet"); r = p.add_run(it[0]); r.bold = True; p.add_run(it[1])
        else:
            d.add_paragraph(str(it), style="List Bullet")


def table(d, df, font=7.5, max_rows=None, col_widths=None):
    df = df.copy()
    if max_rows is not None and len(df) > max_rows:
        df = df.head(max_rows)
    t = d.add_table(rows=1, cols=len(df.columns)); t.style = "Light Grid Accent 1"
    for j, c in enumerate(df.columns):
        cell = t.rows[0].cells[j]; cell.text = str(c)
        for r in cell.paragraphs[0].runs:
            r.bold = True; r.font.size = Pt(font)
    for _, row in df.iterrows():
        cells = t.add_row().cells
        for j, v in enumerate(row):
            cells[j].text = NULL if (v is None or (isinstance(v, float) and pd.isna(v))) else str(v)
            for r in cells[j].paragraphs[0].runs:
                r.font.size = Pt(font)
    if col_widths:
        for j, w in enumerate(col_widths):
            for row in t.rows:
                row.cells[j].width = Cm(w)
    d.add_paragraph()
    return t


def save(d, path):
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True); d.save(path); return path


def md5(path, chunk=1 << 24):
    m = hashlib.md5()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(chunk), b""):
            m.update(b)
    return m.hexdigest()


def human(n):
    for u in ["B", "KB", "MB", "GB"]:
        if n < 1024:
            return f"{n:.1f} {u}"
        n /= 1024
    return f"{n:.1f} TB"


LFC, PADJ = 0.25, 0.05
fmt = lambda x: f"{x:,.0f}"
# ---------------------------------------------------------------- values used in the text
cm = pd.read_csv(PROC / "cell_metadata.csv.gz", index_col=0, low_memory=False)
fin = cm[cm.status == "final"]
thr = pd.read_csv(PIPE / "sample_qc_thresholds.csv").set_index("SampleName")
thr.insert(0, "animal", thr.index.map({k: v[1] for k, v in SAMPLES.items()}))
cap = pd.read_excel(TAB_OUT / "table1_samples_QC.xlsx", sheet_name="gene_cap_loss").set_index("sample_display")
mixed = pd.read_excel(TAB_OUT / "table1_samples_QC.xlsx", sheet_name="glia_in_mixed_clusters")
excl = cm.status.value_counts()
sh = pd.read_csv(DATA / "Cdh9_shared_DEGs.csv")
jac = pd.read_csv(FIGS / "Fig2D_DEG_jaccard_CES_CPS_source_data.csv").set_index("subtype")
gw = pd.read_csv(DATA / "fig3a_followup" / "gsea_pseudobulk_Cdh9_CES_CPS_wide.csv")
pe = pd.read_csv(FIGS / "Fig3E_progeny_Cdh9_source_data.csv")
rf = pd.read_csv(DATA / "nmf" / "robust_factors.csv"); mpg = pd.read_csv(DATA / "nmf" / "metaprogram_genes.csv")
r_shared = np.corrcoef(sh.CES, sh.CPS)[0, 1]
V = dict(n_qc=len(cm), n_final=fmt(len(fin)), n_neurons=fmt(fin.major_celltype.isin(["ExNeuron", "InNeuron"]).sum()),
         n_ex=fmt((fin.major_celltype == "ExNeuron").sum()), n_in=fmt((fin.major_celltype == "InNeuron").sum()),
         cdh9_jaccard=f"{jac.loc['Cdh9+', 'jaccard']:.3f}", cdh9_shared=int(jac.loc["Cdh9+", "n_shared"]),
         cdh9_ces_only=int(jac.loc["Cdh9+", "n_CES"] - jac.loc["Cdh9+", "n_shared"]), cdh9_cps_only=int(jac.loc["Cdh9+", "n_CPS"] - jac.loc["Cdh9+", "n_shared"]),
         n_up=int((sh["class"] == "Common up").sum()), n_down=int((sh["class"] == "Common down").sum()), n_opp=int(sh["class"].str.startswith("Opposite").sum()),
         r_shared=f"{r_shared:.2f}", n_conc=int(sh["class"].isin(["Common up", "Common down"]).sum()),
         cps3_thr=fmt(thr.loc["24090566_AC_PS_973", "threshold_used"]), n_robust=len(rf), n_mp=mpg.shape[1],
         gsea_conc=int(gw.concordant.fillna(False).astype(bool).sum()))
up_p = [f.split(" | ")[1] for f in pe.loc[pe.direction == "up_in_both", "feature"]]
dn_p = [f.split(" | ")[1] for f in pe.loc[pe.direction == "down_in_both", "feature"]]
n_mix = int(cm.status.str.startswith("excluded: mixed").sum()); n_lin = int(excl.get("excluded: lineage check", 0))
n_p11 = int(excl.get("excluded: ExNeuron pass-1 cluster 11 (low quality)", 0)); n_in7 = int(excl.get("excluded: InNeuron cluster 7 (mixed Ex/In)", 0))
valid = thr[thr.threshold_rule == "KDE valley"].threshold_used
scr = thr.set_index("animal").n_predicted_doublets
low_scr = scr[scr < 50]; hi_scr = scr[scr >= 50]
ver = pd.read_csv(CODE / "environments" / "key_package_versions.csv")
VV = {(r.environment, r.package): r.version for r in ver.itertuples()}
py = lambda p: VV[("sc-pipeline-python", p)]; rr = lambda p: VV[("brain-enrich-r", p)]
print({k: V[k] for k in list(V)[:12]})

# ================================================================ 1_methods.docx
P = [
 ("QC", "genes per nucleus", "200-4,000", "pipeline_lib.py"), ("QC", "UMIs per nucleus", "200-15,000", "pipeline_lib.py"),
 ("QC", "mitochondrial UMI fraction", "<= 1%", "pipeline_lib.py"),
 ("QC", "doublets", "scanpy.external.pp.scrublet, default parameters, random_state 0, per animal; predicted doublets removed", "0_build_prefilter_from_matrices.py"),
 ("QC", "genes", "mt- genes removed; protein-coding genes (GENCODE vM36, 19,074) used for clustering and analysis", "1_qc_split.py; 2_cluster_annotate_major.py"),
 ("Split", "low/high-feature threshold", "Gaussian KDE (1,000-point grid) of genes per nucleus; deepest minimum between the lowest (glial) peak and the tallest peak above it, below 2x the glial peak; peaks < 5% of max density ignored; median of the other animals if none", "1_qc_split.py"),
 ("Clustering", "normalisation", "normalize_total 1e4, log1p", "pipeline_lib.cluster"),
 ("Clustering", "HVG", "2,000, batch_key = animal", "pipeline_lib.cluster"),
 ("Clustering", "regression", "n_genes_by_counts, total_counts, total_counts_mt (major types); total_counts, total_counts_mt (neurons)", "pipeline_lib.cluster"),
 ("Clustering", "scaling / PCA", "scale max 10; PCA 50 components, random_state 7", "pipeline_lib.cluster"),
 ("Clustering", "integration", "Harmony (harmonypy), key = animal, max_iter 50", "pipeline_lib.cluster"),
 ("Clustering", "graph / UMAP", "20 neighbours (50 PCs major types; 20 PCs ExNeuron; 10 PCs InNeuron); UMAP", "2; 3"),
 ("Clustering", "Leiden resolution", "1.0 major types; 0.5 ExNeuron pass 1 and 2; 0.5 restricted re-clustering of ExNeuron cluster 0; 0.3 InNeuron", "2; 3"),
 ("Annotation", "lineage check", "glial/vascular-labelled nucleus removed if 0 of its lineage markers detected and >= 2 of 4 pan-neuronal genes (Syn1, Snap25, Rbfox3, Syt1)", "2_cluster_annotate_major.py"),
 ("Embedding", "final UMAP", "HVG 2,000 by animal, regress total_counts + total_counts_mt, scale 10, PCA 50 (seed 7), Harmony, 20 neighbours on 30 PCs", "4_build_objects.py"),
 ("DE", "cell-level test", "scanpy rank_genes_groups method='t-test' (Welch), log-normalised expression, BH", "5_de_pseudobulk_export.py"),
 ("DE", "DEG threshold", f"BH adj. P < {PADJ} and |log2FC| > {LFC} (all analyses)", "params.json"),
 ("DE", "contrasts", "CES vs CON, CPS vs CON, CES vs CPS, SDS (CES + CPS) vs CON; per major type and per subtype", "5_de_pseudobulk_export.py"),
 ("Pseudobulk", "limma-voom", "summed UMIs per animal x group; >= 10 nuclei per animal and >= 2 animals per arm; filterByExpr, TMM, voom, ~0 + condition", "5; 6"),
 ("GO ORA", "enrichGO", "org.Mm.eg.db SYMBOL, ont ALL, BH, P < 0.05, q < 0.2, gene sets 10-500, simplify 0.7 (Wang)", "6_stats.R"),
 ("GO ORA", "universe", "genes detected in >= 1% of nuclei of the group in either arm (main); all tested genes (sensitivity)", "6_stats.R"),
 ("GO ORA", "gene exclusion", "^(mt-|Rpl|Rps|Hsp) removed from per-contrast DEG lists (not from the Cdh9+ shared sets)", "6_stats.R"),
 ("GSEA", "gseGO", "Cdh9+ pseudobulk limma t ranks, ont ALL, gene sets 15-500, BH, eps 0, fgsea multilevel, seed TRUE; concordant = BH < 0.05 in both contrasts, same NES sign", "6_stats.R"),
 ("Composition", "test", "neuron subtypes: per-animal proportions, +0.5 pseudocount, logit; Welch t-test vs CON (n = 3 vs 3), BH across subtypes", "common.logit_prop_test; 14_figure2.py"),
 ("Disease overlap", "test", "7 GWAS traits (human -> mouse, homologene); one-sided hypergeometric, BH; background = genes detected in >= 1% of the type; neuron-enriched gene sensitivity", "13_figure1EG.R"),
 ("CellChat", "settings", "CellChatDB.mouse, triMean, population.size FALSE, nboot 100, seed 20260724, min.cells 30, thresholds 0.05; L-R set = union of over-expressed interactions per animal", "8_cellchat.R"),
 ("CellChat", "balancing", "groups with >= 30 nuclei in every animal; downsampled to min(smallest animal, 500) nuclei per animal (seed 0)", "7_cellchat_export.py"),
 ("PROGENy", "settings", "mouse model, top 500 genes; decoupler ULM (tmin 5) on per-animal mean log-expression of groups with >= 20 nuclei; Welch t-test vs CON; BH within contrast", "6; 9"),
 ("NMF", "factorisation", "ExNeuron per animal; 2,000 HVGs (seurat flavour, batch = animal); log-expression centred per gene within animal, negatives set to 0; sklearn NMF k = 4-9, nndsvda, Frobenius, cd, max_iter 1000, tol 1e-4, random_state 0", "10_nmf_programs.py"),
 ("NMF", "meta-programs", "top-50 genes; within-animal recurrence >= 35/50, cross-animal >= 10/50, redundancy < 10/50; average linkage on 1 - overlap/50 cut 0.75; >= 3 animals; top 30 genes per MP", "10_nmf_programs.py"),
 ("NMF", "scoring / test", "scanpy score_genes; per animal x subtype mean; Welch t-test vs CON; BH within contrast", "16_nmf_scores_figureS2.py"),
 ("Global", "random seeds", "0 (analysis), 7 (PCA/Harmony), 20260724 (CellChat), 1 (step 6, GSEA)", "config.py; scripts")]
par = pd.DataFrame(P, columns=["step", "parameter", "value", "script"])

d = new_doc("Methods (draft)", "Written from the code in code/. Fields marked null must be completed by the lab.")
h(d, "Animals, stress paradigms and tissue (to be completed)")
para(d, "Mice: strain null, sex null, age null, supplier null, housing null. Chronic stress models: CES (raw label ES) protocol null; CPS (raw label PS) protocol null; controls (CON) null. "
        "Ethics approval: committee null, protocol number null. ACC dissection: coordinates null. n = 3 animals per group (9 libraries), one library per animal.")
h(d, "Single-nucleus RNA sequencing (to be completed)")
para(d, "Nuclei isolation null; 10x Genomics chemistry null; target nuclei per library null; sequencer null; read configuration null. Reads were processed with Cell Ranger (version null; reference null), "
        "and filtered feature-barcode matrices were imported with scanpy read_10x_mtx (gene symbols). Protein-coding genes were defined by GENCODE vM36 (gene_type protein_coding).")
h(d, "Quality control and low/high-feature split")
para(d, f"Per animal, QC metrics were computed with scanpy {py('scanpy')} (mitochondrial genes = mt- prefix). Nuclei with 200-4,000 detected genes, 200-15,000 UMIs and <= 1% mitochondrial UMIs were kept, "
        f"and Scrublet doublets (scanpy.external.pp.scrublet, default parameters, random_state 0) were removed, leaving {fmt(V['n_qc'])} nuclei (per animal {fmt(thr.n_after.min())}-{fmt(thr.n_after.max())}). "
        "Because neuronal nuclei carry about twice as many genes as glial nuclei, each animal was split into a low-feature and a high-feature object at the valley of its genes-per-nucleus density (Gaussian KDE): "
        "the deepest minimum between the lowest-gene (glial) peak and the tallest peak above it, required to lie below twice the glial peak, ignoring peaks below 5% of the maximum density. "
        f"Thresholds were {fmt(valid.min())}-{fmt(valid.max())} genes; CPS_3 had no valid valley (neuronal peak at {fmt(thr.loc['24090566_AC_PS_973', 'neuronal_peak'])} genes), so the median of the other animals ({V['cps3_thr']} genes) was used (Supplementary Fig. S1; table1_samples_QC).")
para(d, f"Scrublet called almost no doublets in {len(low_scr)} animals ({', '.join(f'{a}: {n}' for a, n in low_scr.items())}, against {fmt(hi_scr.min())}-{fmt(hi_scr.max())} in the others). "
        "Doublets in these animals were removed only by the cluster- and nucleus-level steps below; this limits the composition analysis.", italic=True)
h(d, "Clustering and annotation")
para(d, f"Low- and high-feature objects were clustered separately: normalisation to 10,000 UMIs and log1p, 2,000 highly variable genes (batch = animal), regression of genes, UMIs and mitochondrial UMIs per nucleus, scaling (max 10), PCA (50 components, seed 7), "
        f"Harmony integration by animal (harmonypy {py('harmonypy')}), a 20-nearest-neighbour graph, UMAP and Leiden clustering (resolution 1.0; leidenalg {py('leidenalg')}). "
        "Each cluster received one major label from canonical markers (ExNeuron Slc17a7/Slc17a6/Slc17a8; InNeuron Gad1/Gad2; Astro Slc1a2/Gpc5/Slc1a3/Gja1; Oligo St18/Prr5l/Mog; OPC Vcan/Pdgfra; Micro Inpp5d/Tgfbr1/Hexb/Selplg; Endo Slco1a4/Mecom/Flt1; Mural Slc6a20a/Colec12/Rgs5). "
        f"Clusters co-expressing glial and neuronal markers ({fmt(n_mix)} nuclei) were excluded, and a nucleus-level lineage check removed glial/vascular-labelled nuclei expressing none of their lineage markers but >= 2 of Syn1, Snap25, Rbfox3 and Syt1 ({fmt(n_lin)} nuclei).")
para(d, f"Excitatory neurons were re-clustered (regression of UMIs and mitochondrial UMIs, PCA, Harmony, 20 neighbours on 20 PCs, Leiden 0.5). One pass-1 cluster of low-quality nuclei ({fmt(n_p11)} nuclei; median 864 genes, lowest pan-neuronal score) was removed and the rest re-clustered with the same settings. "
        "Pass-2 cluster 0, which co-expressed Cux2, Col19a1 and Trpc6, was re-clustered (Leiden 0.5, restricted) and split into Cux2+Col19a1 and Cux2+Trpc6 by mean Trpc6 versus Col19a1. "
        f"Inhibitory neurons were clustered with 10 PCs and resolution 0.3; one cluster with low Gad1/Gad2 and excitatory markers ({fmt(n_in7)} nuclei) was removed. "
        "Subtypes were named from marker expression and Wilcoxon top markers: 12 ExNeuron (Foxp2+, Cux2+Rorb, Cux2+Col19a1, Cdh9+, Cux2+Trpc6, Ptprm+Nell1, Ptprm+Tshz2, Ntng1+, Cux2+Fth1, Kcnip1+, Cux2+Htr2c, Arhgap6+) and 5 InNeuron subtypes (Pvalb+, Sst+, Vip+, Lamp5+, Cnr1+). "
        f"The final dataset has {V['n_final']} nuclei ({V['n_neurons']} neurons). A common UMAP of all final nuclei was computed for figures (2,000 HVGs by animal, regression of UMIs and mitochondrial UMIs, PCA 50, Harmony, 20 neighbours on 30 PCs).")
h(d, "Differential expression")
para(d, "For each major type and subtype, CES vs CON, CPS vs CON, CES vs CPS and SDS (CES + CPS) vs CON were compared at the nucleus level with scanpy rank_genes_groups (method = 't-test', Welch; log-normalised expression; Benjamini-Hochberg). "
        f"DEGs were defined as BH-adjusted P < {PADJ} and |log2FC| > {LFC} throughout. Because nuclei from the same animal are not independent, key results were checked with a replicate-level pseudobulk model: UMIs summed per animal and cell group (>= 10 nuclei per animal), "
        f"edgeR {rr('edgeR')} filterByExpr and TMM, limma {rr('limma')} voom, design ~0 + condition, contrasts CES-CON and CPS-CON (n = 3 vs 3 animals). "
        f"The CES/CPS overlap per subtype was summarised with the Jaccard index of DEG sets; Cdh9+ DEGs shared by both models with the same sign define the common-up (n = {V['n_up']}) and common-down (n = {V['n_down']}) sets.")
h(d, "Functional enrichment")
para(d, f"GO over-representation used clusterProfiler {rr('clusterProfiler')} enrichGO with org.Mm.eg.db {rr('org.Mm.eg.db')} (all ontologies, BH, P < 0.05, q < 0.2, gene-set size 10-500) followed by simplify (cutoff 0.7, Wang similarity). "
        "The universe was the set of genes detected in >= 1% of nuclei of the group in either arm; an all-genes universe was run as a sensitivity analysis. Mitochondrial, ribosomal-protein and Hsp genes were removed from per-contrast DEG lists, but not from the Cdh9+ shared sets. "
        f"Threshold-free GSEA (gseGO, fgsea {rr('fgsea')} multilevel, set size 15-500) was run on Cdh9+ pseudobulk limma t statistics; terms with BH < 0.05 in both contrasts and the same NES sign were called concordant. "
        "For Fig. 3A, significant terms of each shared set were clustered by gene overlap (Jaccard distance, average linkage, cut 0.7) and the term with the lowest adjusted P represents each cluster.")
h(d, "Cell-type composition")
para(d, "Major cell-type proportions are shown pooled per group (Fig. 1E, descriptive). Neuron subtype proportions were computed per animal (+0.5 pseudocount), logit-transformed and compared with CON by Welch t-test (n = 3 vs 3 animals; BH across subtypes), with 95% confidence intervals of the logit difference (table2; per-animal composition in Fig. 2B).")
h(d, "Overlap with psychiatric GWAS gene sets")
para(d, f"Gene sets for ADHD, bipolar disorder, depression, PTSD, anxiety, insomnia and neuroticism (data/original/reference/GWAS_Neuron_Disease_Gene_Dataset.csv; source null) were mapped from human to mouse with homologene {rr('homologene')}. "
        "SDS-vs-CON DEGs of each major type were tested with a one-sided hypergeometric test against genes detected in >= 1% of that type (BH). For glial types, neuron-enriched genes (CON ExNeuron mean log-expression >= 0.5 and >= 2x the glial mean) were flagged and the test repeated without them. "
        f"Chord diagrams: circlize {rr('circlize')}.")
h(d, "Cell-cell communication")
para(d, f"CellChat {VV[('ccc', 'CellChat')]} (R {VV[('ccc', 'R')]}) was run per condition on balanced nuclei: groups with >= 30 nuclei in every animal, downsampled to min(smallest animal count, 500) nuclei per animal. "
        "The ligand-receptor set was the union of over-expressed interactions found in each animal (CellChatDB.mouse). Communication probabilities used triMean, population.size = FALSE, 100 permutations, filterCommunication(min.cells = 30) and thresholds of 0.05. "
        "Two groupings were used: fine (neuron subtypes plus glial/vascular types; Fig. 3C) and merged (Cdh9+, other ExNeurons, InNeurons plus glial/vascular types; Fig. 3C merged). Differential interaction number and strength (stress minus CON) are shown for edges to and from Cdh9+.")
h(d, "Pathway activity (PROGENy)")
para(d, f"Mean log-expression per animal and cell group (>= 20 nuclei) was scored with the PROGENy mouse model (top 500 genes; progeny {rr('progeny')}) using the univariate linear model of decoupler {py('decoupler')} (tmin 5), and compared with CON by Welch t-test (n = 3 vs 3; BH within contrast).")
h(d, "NMF gene programs")
para(d, f"ExNeuron nuclei were factorised per animal (scikit-learn {py('scikit-learn')} NMF, k = 4-9, nndsvda, Frobenius loss, seed 0) on 2,000 HVGs, using log-expression centred per gene within each animal with negative values set to zero (Gavish et al., Nature 618, 598-606, 2023; doi:10.1038/s41586-023-06130-4). "
        f"Robust factors (>= 35/50 top genes shared with another k in the same animal, >= 10/50 with another animal, non-redundant within animal; {V['n_robust']} factors) were clustered (average linkage on 1 - overlap/50, cut 0.75; >= 3 animals) into {V['n_mp']} meta-programs defined by their 30 most recurrent genes and annotated from GO enrichment. "
        "MP scores were averaged per animal and subtype and compared with CON by Welch t-test (BH within contrast).")
h(d, "Statistics and software")
para(d, "Tests are two-sided unless stated otherwise; multiple testing is controlled with Benjamini-Hochberg. n = 3 animals per group. Cell-level statistics treat nuclei as observations and are reported with animal-level (pseudobulk) checks. "
        f"Software: Python {py('python')}, scanpy {py('scanpy')}, anndata {py('anndata')}, numpy {py('numpy')}, scipy {py('scipy')}, pandas {py('pandas')}, scikit-learn {py('scikit-learn')}, statsmodels {py('statsmodels')}, decoupler {py('decoupler')}, harmonypy {py('harmonypy')}, leidenalg {py('leidenalg')}, matplotlib {py('matplotlib')}; "
        f"R {rr('R')}, clusterProfiler {rr('clusterProfiler')}, org.Mm.eg.db {rr('org.Mm.eg.db')}, GO.db {rr('GO.db')}, limma {rr('limma')}, edgeR {rr('edgeR')}, fgsea {rr('fgsea')}, progeny {rr('progeny')}, homologene {rr('homologene')}, circlize {rr('circlize')}, ggalluvial {rr('ggalluvial')}, ggplot2 {rr('ggplot2')}; CellChat {VV[('ccc', 'CellChat')]}. Full lists: code/environments/.")
h(d, "Parameters"); table(d, par, font=6.5)
h(d, "Per-animal QC and split thresholds")
table(d, thr.reset_index()[["animal", "SampleName", "n_before", "n_predicted_doublets", "n_after", "glial_peak", "neuronal_peak", "threshold_used", "threshold_rule", "n_low", "n_high"]].round(0), font=6.5)
save(d, DOC_OUT / "1_methods.docx")

# ================================================================ 2_code_figure_table_map.docx
OUTS = {"0": "data/original/all_cell_adata_combined_pre_filter.h5ad", "1": "pipeline/sample_qc_thresholds.csv, {high,low}_feature objects; qc/kde_split_* (Fig. S1)",
        "2": "pipeline/adata_*_cluster.h5ad, markers_*.csv, lineage_check_removed.csv, ExNeuron pass 1", "3": "pipeline/adata_ExNeuron_cluster2.h5ad, adata_InNeuron_cluster.h5ad, naming_*.csv",
        "4": "h5ad/ACC_all_cells.h5ad, ACC_neurons.h5ad, data/ACC_cell_metadata.csv, logs/misplaced_check.csv", "5": "data/DE_celllevel_*.csv.gz, pseudobulk_*.csv.gz, sample_mean_logexpr_*.csv.gz",
        "6": "data/progeny_mouse_top500.csv.gz, DE_pseudobulk_*.csv.gz, GO_*.csv, Cdh9_shared_DEGs.csv, fig3a_followup/ (GSEA, representative terms)",
        "7": "data/cellchat_{fine,coarse}/ inputs", "8": "data/cellchat_*/*.rds, net_*.csv, communications.csv", "9": "data/progeny_*.csv",
        "10": "data/nmf/ factors and meta-programs", "11": "data/nmf/metaprogram_GO_*.csv", "12": "Fig. 1A-D, 1F", "13": "Fig. 1E, 1G", "14": "Fig. 2A-F",
        "15": "Fig. 3A-F", "16": "data/nmf/ scores and tests; Fig. S2", "17": "figures/, tables/, data/processed/cell_metadata.csv.gz, cellchat/, checksums_md5.csv", "18": "doc/*.docx"}
def purpose(p):
    t = p.read_text(errors="ignore")
    m = re.search(r'"""(.*?)(\n\n|""")', t, re.S)
    if p.suffix == ".py" and m:
        return " ".join(m.group(1).split())[:230]
    lines = [l.lstrip("# ").strip() for l in t.splitlines()[:6] if l.startswith("#") and not l.startswith("#!")]
    return " ".join(lines)[:230]
srows = []
for p in sorted(CODE.glob("*"), key=lambda p: (int(re.match(r"(\d+)", p.name).group(1)) if re.match(r"\d+", p.name) else 99, p.name)):
    if p.suffix not in (".py", ".R", ".sh") or p.name.startswith("_"):
        continue
    step = re.match(r"(\d+[ab]?)_", p.name)
    srows.append(dict(step=step.group(1) if step else "helper", script=p.name, language={"py": "Python", "R": "R", "sh": "bash"}[p.suffix[1:]],
                      purpose=purpose(p), main_outputs=OUTS.get(step.group(1), "") if step else ""))
scr = pd.DataFrame(srows)
frows = []
for panel, name, dd, stem, scripts, inputs, legend in FIGURES:
    srcs = sorted(FIG_OUT.glob(f"{name}_source_data*.csv"))
    frows.append(dict(panel=panel, figure_file=f"figures/{name}.pdf|png", source_data="; ".join(f"figures/{s.name}" for s in srcs), scripts=scripts, inputs=inputs))
fmap = pd.DataFrame(frows)
trows = [dict(table=f"tables/{tid}.xlsx", title=t[0], sheets=", ".join(s[0] for s in t[3]), scripts=t[2], statistical_unit=t[1]) for tid, t in TABLES.items()]
trows += [dict(table=f"tables/{p.name}", title="full DE results (all genes)", sheets="-", scripts="5_de_pseudobulk_export.py; 6_stats.R", statistical_unit="nucleus / animal") for p in sorted(TAB_OUT.glob("table3_full_*"))]
tmap = pd.DataFrame(trows)
ck = pd.read_csv(PROC / "checksums_md5.csv")
d = new_doc("Code, figure and table map", "Which script produces which figure, table and object; how to run the code.", landscape=True)
h(d, "1. Project layout")
bullets(d, ["data/original/ - filtered_mat/ (one placeholder folder per sample; copy the Cell Ranger matrices here), filtered_matrix_paths.csv (cluster locations), reference/ (protein-coding gene list, GWAS gene sets, SynGO 1.3).",
            "data/processed/ - h5ad/ (final objects), cellchat/ (CellChat objects), cell_metadata.csv.gz, checksums_md5.csv; intermediate/ = working files of steps 1-16, regenerated by run_all.sh (pipeline objects, result tables, raw figure files, logs).",
            "code/ - numbered steps 0-18, helpers (config.py, config.R, common.py, pipeline_lib.py, manifest.py), params.json, run_all.sh, environments/, logs/.",
            "figures/ - figure1A ... figure3F, figureS1, figureS2 (pdf, png, source-data csv). tables/ - table1 ... table10 (xlsx) and full DE results (csv.gz). doc/ - these Word documents."])
h(d, "2. How to run")
bullets(d, ["Step 0 (cluster or after copying the matrices into data/original/filtered_mat/): python 0_build_prefilter_from_matrices.py -> data/original/all_cell_adata_combined_pre_filter.h5ad.",
            "Steps 1-18: bash code/run_all.sh (set ACC_PREFILTER if the pre-filter object is elsewhere; ACC_PY, ACC_R_ENRICH, ACC_R_CCC, ACC_CCC_RLIB select the interpreters). Run time on the analysis Mac (16 cores, 48 GB): about 45 min.",
            "All paths are relative to the project folder (config.py / config.R), so the folder can be moved as a whole.",
            "Environments: sc-pipeline-python (Python steps), brain-enrich-r (R steps except CellChat), ccc + ccc-rlib (CellChat). Package lists: code/environments/. Step 18 needs python-docx."])
h(d, "3. Scripts"); table(d, scr, font=6.5)
h(d, "4. Figures"); table(d, fmap, font=6)
h(d, "5. Tables"); table(d, tmap, font=6.5)
h(d, "6. Processed objects (md5)"); table(d, ck, font=6.5)
vf = WORK / "rerun_verification.csv"
if vf.exists():
    h(d, "7. Rerun verification")
    para(d, "The consolidated pipeline was rerun and its results compared with those of the previous version of this folder (same inputs and software).")
    table(d, pd.read_csv(vf), font=6.5)
save(d, DOC_OUT / "2_code_figure_table_map.docx")

# ================================================================ 3_data_and_objects.docx
fm = pd.read_csv(ORIG / "filtered_matrix_paths.csv")
def obj_desc(p):
    a = ad.read_h5ad(p, backed="r")
    lay = ", ".join(a.layers.keys()) or "none"
    return dict(file=str(p.relative_to(ROOT)), nuclei=a.n_obs, genes=a.n_vars, X="log-normalised (ln(1 + counts per 10,000))", layers=lay,
                obsm=", ".join(a.obsm.keys()), obs_columns=", ".join(a.obs.columns))
objs = pd.DataFrame([obj_desc(H5AD / f) for f in ("ACC_all_cells.h5ad", "ACC_neurons.h5ad")])
CM_DICT = {"cell_id": "barcode-raw_SampleID", "barcode": "10x cell barcode", "raw_SampleID": "Cell Ranger folder name", "sample_display": "animal name used in figures (CON_1 ... CPS_3)",
           "sample_pipeline": "animal name inside the original pipeline", "group": "CON / CES / CPS", "pooled_group": "CON / SDS (CES + CPS)", "n_genes_by_counts": "detected genes",
           "total_counts": "UMIs", "pct_counts_mt": "mitochondrial UMI %", "pct_counts_ribo": "ribosomal-protein UMI %", "doublet_score": "Scrublet score",
           "split_object": "low_feature / high_feature (step 1)", "leiden_split_object": "Leiden cluster in its split object (step 3)", "major_celltype_split_object": "cluster label incl. mixed X_ExNeuron (step 4)",
           "leiden_ExNeuron_pass1": "ExNeuron pass-1 cluster", "leiden_ExNeuron_pass2": "ExNeuron pass-2 cluster", "leiden_ExNeuron_sub": "pass-2 cluster after splitting cluster 0",
           "leiden_InNeuron": "InNeuron cluster", "major_celltype": "final major type (final nuclei only)", "subtype": "final subtype / type", "status": "final or exclusion reason",
           "UMAP1_allcells": "all-cell UMAP", "UMAP2_allcells": "all-cell UMAP", "UMAP1_neurons": "neuron UMAP", "UMAP2_neurons": "neuron UMAP"}
nulls = pd.DataFrame([("data/original/filtered_mat", "Cell Ranger matrices (copy from cluster paths)"), ("FASTQ", "paths and md5 per library"), ("sequencing", "platform, read layout, run/batch"),
                      ("10x chemistry", "kit version"), ("Cell Ranger", "version and reference genome"), ("animals", "strain, sex, age, supplier, housing"), ("stress protocols", "CES and CPS procedures, duration"),
                      ("tissue", "collection date, dissection coordinates, nuclei isolation"), ("ethics", "committee and protocol number"), ("GWAS gene sets", "origin of GWAS_Neuron_Disease_Gene_Dataset.csv"),
                      ("deposition", "GEO accession, reviewer token, repository DOI, code URL")], columns=["field", "missing information (null)"])
dep = pd.DataFrame([dict(item="raw FASTQ + filtered matrices (9 libraries)", repository="GEO/SRA", accession=NULL, reviewer_token=NULL),
                    dict(item="data/processed (h5ad, cell metadata, CellChat)", repository="GEO supplementary or Zenodo", accession=NULL, reviewer_token=NULL),
                    dict(item="figures/ source data and tables/", repository="journal supplementary + Zenodo", accession=NULL, reviewer_token=NULL),
                    dict(item="code/", repository="GitHub + Zenodo release", accession=NULL, reviewer_token=NULL)])
d = new_doc("Data and objects", "Raw data, processed objects, metadata dictionary and data availability.", landscape=True)
h(d, "1. Raw data (data/original)")
para(d, "Nine single-nucleus RNA-seq libraries from mouse anterior cingulate cortex, one per animal: 3 CON, 3 CES (raw label ES), 3 CPS (raw label PS); SDS = CES + CPS. "
        "The Cell Ranger filtered matrices are on the compute cluster; data/original/filtered_mat/ holds one placeholder folder per sample to copy them into.")
table(d, fm.drop(columns=["local_folder"]), font=6.5)
para(d, "Reference files (data/original/reference): coding_genes_gencode_vM36.txt (19,074 protein-coding genes used in all objects), GWAS_Neuron_Disease_Gene_Dataset.csv (7 traits), syngo/ (SynGO 1.3 release).")
h(d, "2. Processed objects (data/processed)")
table(d, objs, font=6.5)
bullets(d, ["X = log-normalised expression; layers['counts'] = raw UMI counts (integers); obsm['X_umap'] = UMAP used in the figures; obsm['X_pca_harmony'] = Harmony-corrected PCs.",
            "cellchat/cellchat_{fine,coarse}_{CON,CES,CPS}_CellChat.rds: CellChat objects per condition; *_merged_CON_{CES,CPS}.rds: merged objects for comparisons (readRDS in R; net$count / net$weight).",
            "cell_metadata.csv.gz: all QC-passing nuclei with their final status (dictionary below); intermediate/: working files of the pipeline, regenerated by run_all.sh."])
h(d, "3. Cell metadata dictionary"); table(d, pd.DataFrame(list(CM_DICT.items()), columns=["column", "meaning"]), font=7)
h(d, "4. Samples"); table(d, pd.read_excel(TAB_OUT / "table1_samples_QC.xlsx", sheet_name="samples"), font=6)
h(d, "5. Fields still null"); table(d, nulls, font=7)
h(d, "6. Data and code availability (draft)")
para(d, f"The snRNA-seq data (raw FASTQ and Cell Ranger filtered matrices for nine libraries) have been deposited in GEO under accession {NULL} (reviewer token {NULL}). "
        f"Processed AnnData objects ({V['n_final']} nuclei; log-normalised expression with raw counts in layers['counts']), cell metadata and CellChat objects are available at {NULL} (DOI {NULL}). Source data are provided with this paper. "
        f"Code is available at {NULL} and archived at Zenodo (DOI {NULL}).")
table(d, dep, font=7)
save(d, DOC_OUT / "3_data_and_objects.docx")

# ================================================================ 4_qc_and_annotation.docx
exw = pd.read_excel(TAB_OUT / "table1_samples_QC.xlsx", sheet_name="exclusions")
mxp = mixed.pivot(index="celltype", columns="sample", values="frac_in_mixed").reindex(columns=[v[1] for v in SAMPLES.values()]).round(2)
cps = mixed[mixed["sample"].str.startswith("CPS")].groupby("celltype").frac_in_mixed.agg(["min", "max"])
d = new_doc("QC, exclusions and annotation", "What was removed, why, and the evidence behind the cell-type labels.", landscape=True)
h(d, "1. Nucleus flow")
para(d, f"QC-passing nuclei ({fmt(V['n_qc'])}) -> low/high-feature split -> clustering -> exclusions below -> {V['n_final']} final nuclei.")
table(d, exw, font=6.5)
h(d, "2. Exclusion criteria")
table(d, pd.DataFrame([
  dict(step="Mixed clusters (X_ExNeuron, Oligo_Astro)", nuclei=n_mix, criterion="Leiden cluster co-expresses glial and neuronal markers (table10 markers_* sheets)", rationale="neuron-glia doublets / ambient RNA"),
  dict(step="Lineage check", nuclei=n_lin, criterion="glial/vascular label, 0 own-lineage markers, >= 2 of 4 pan-neuronal genes", rationale="neuronal nuclei in glial clusters (table10 lineage_check)"),
  dict(step="ExNeuron pass-1 cluster 11", nuclei=n_p11, criterion="median 864 genes (others 1,869-2,681), lowest pan-neuronal score 0.62", rationale="low-quality neuronal nuclei"),
  dict(step="InNeuron cluster 7", nuclei=n_in7, criterion="Gad1/Gad2 0.39/0.38 (others 0.6-2.0), no subtype marker, excitatory top markers", rationale="mixed Ex/In cluster")]), font=7)
h(d, "3. Low/high-feature split (Fig. S1)")
para(d, f"Each animal was split at the valley of its genes-per-nucleus density (Gaussian KDE): the deepest minimum between the glial peak and the tallest peak above it, below twice the glial peak; peaks below 5% of the maximum density are ignored. "
        f"Thresholds were {fmt(valid.min())}-{fmt(valid.max())} genes in {len(valid)} animals; CPS_3 has no valid valley (glial peak {fmt(thr.loc['24090566_AC_PS_973', 'glial_peak'])}, neuronal peak {fmt(thr.loc['24090566_AC_PS_973', 'neuronal_peak'])} genes) and uses the median of the others ({V['cps3_thr']}). Per animal: table1 split_kde.")
h(d, "4. Doublet detection and composition caveats")
para(d, f"Scrublet's automatic threshold called {', '.join(f'{a}: {n}' for a, n in low_scr.items())} doublets, against {fmt(hi_scr.min())}-{fmt(hi_scr.max())} in the other animals, so doublets in these animals were removed only by the mixed-cluster and lineage steps. "
        f"The fraction of glia-labelled nuclei falling in mixed clusters is highest in CPS (Astro {cps.loc['Astro', 'min']:.2f}-{cps.loc['Astro', 'max']:.2f}, Oligo {cps.loc['Oligo', 'min']:.2f}-{cps.loc['Oligo', 'max']:.2f}); "
        "removing mixed clusters therefore removes glia unevenly, and glial proportions, particularly in CPS, are under-estimated.")
table(d, mxp.reset_index(), font=7)
para(d, f"Upper QC cap: CPS_3 was sequenced more deeply (median {fmt(cap.loc['CPS_3', 'median_genes'])} genes per nucleus before QC). The 4,000-gene cap removed {cap.loc['CPS_3', 'frac_gt_4000_genes']:.1%} of its nuclei "
        f"({fmt(cap.loc['CPS_3', 'n_gt_4000_genes'])} of {fmt(cap.loc['CPS_3', 'nuclei_pre_QC'])}) against {cap.drop('CPS_3').frac_gt_4000_genes.min():.1%}-{cap.drop('CPS_3').frac_gt_4000_genes.max():.1%} in the other animals; the cap was applied identically to all animals.")
table(d, cap.reset_index(), font=7)
h(d, "5. Annotation evidence (table10_annotation)")
bullets(d, ["Major types: one label per Leiden cluster from the marker panel; sheets markers_low_feature / markers_high_feature (fraction expressing each marker, lineage scores, cluster size, dominant-animal share, median genes).",
            "Neuron subtypes: naming_ExNeuron / naming_InNeuron (mean expression of naming genes + top 8 Wilcoxon markers); naming_ExNeuron_c0 (Trpc6 vs Col19a1 split of cluster 0).",
            "Post-hoc check: annotation_check (nuclei whose 15 UMAP neighbours are mostly another type, and nuclei with no lineage marker, per type and animal)."])
save(d, DOC_OUT / "4_qc_and_annotation.docx")

# ================================================================ 5_figure_legends.docx
V2 = dict(V, **{})
d = new_doc("Figure legends", "Main and supplementary figures in figures/. Values are read from the results.")
for panel, name, dd, stem, scripts, inputs, legend in FIGURES:
    leg = legend.format(**V)
    if panel == "3E":
        leg += f" Up in both: {', '.join(up_p)}; down in both: {', '.join(dn_p)}."
    para(d, f" {leg}", bold_prefix=f"Fig. {panel.replace('_merged', ' (merged groups)')} ({name}).")
para(d, "DEG threshold everywhere: BH-adjusted P < 0.05 and |log2FC| > 0.25 (cell-level Welch t-test); replicate-level checks use animals as the unit (pseudobulk limma-voom; logit-proportion tests; per-animal scores).", italic=True)
save(d, DOC_OUT / "5_figure_legends.docx")
print("docs written:", sorted(p.name for p in DOC_OUT.glob("*.docx")))
