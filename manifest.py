"""What goes into figures/ and tables/: one entry per figure panel and per table, with the scripts and inputs
behind it. Used by 17_assemble_outputs.py (copy/rename/build) and 18_build_docs.py (map + legends)."""

# panel, output name, source dir ("FIGS" or "QC"), source stem, scripts, inputs, legend
FIGURES = [
 ("1A", "figure1A_umap_major_celltypes", "FIGS", "Fig1A_umap_major_celltypes", "12_figure1.py; 4_build_objects.py", "ACC_all_cells.h5ad",
  "UMAP of {n_final} ACC nuclei from 9 animals (3 CON, 3 CES, 3 CPS), coloured by major cell type; numbers = nuclei per type. Embedding: 2,000 HVGs (batch = animal), regression of total and mitochondrial counts, scaling, PCA (50 components), Harmony (batch = animal), 20 neighbours on 30 PCs, UMAP."),
 ("1B", "figure1B_marker_heatmap", "FIGS", "Fig1B_marker_heatmap", "12_figure1.py", "ACC_all_cells.h5ad",
  "Mean log-normalised expression of canonical markers per major type, min-max scaled per gene (viridis). Bars above columns indicate the cell type each marker defines."),
 ("1C", "figure1C_umap_by_group", "FIGS", "Fig1C_umap_combined", "12_figure1.py", "ACC_all_cells.h5ad",
  "UMAP coloured by group (left: CON, CES, CPS; right: CON vs SDS = CES + CPS pooled), nuclei drawn in random order; nuclei per group in the legend."),
 ("1D", "figure1D_upset_DEGs", "FIGS", "Fig1D_upset_DEGs_merged_CES_CPS", "12_figure1.py; 5_de_pseudobulk_export.py", "table3_DEGs",
  "UpSet of DEGs per major cell type (cell-level Welch t-test, BH adj. P < 0.05, |log2FC| > 0.25). Exclusive cell-type combinations are computed separately for CES vs CON (green) and CPS vs CON (orange) and shown as paired bars over one dot matrix; left bars = DEGs per cell type."),
 ("1E", "figure1E_composition_alluvial", "FIGS", "Fig1E_composition_alluvial", "13_figure1EG.R", "cell_metadata.csv.gz",
  "Major cell-type composition of nuclei pooled per group (CON | SDS and CON | CES | CPS); ribbons connect the same cell type, numbers = nuclei per group. Pooled fractions are descriptive. Glial proportions are under-estimated in CPS (see doc/4_qc_and_annotation)."),
 ("1F", "figure1F_ExNeuron_GO", "FIGS", "Fig1F_ExNeuron_GO_overlap", "12_figure1.py; 6_stats.R", "table5_GO_enrichment",
  "GO enrichment (clusterProfiler enrichGO, all ontologies, BH, simplify 0.7; universe = genes detected in >= 1% of ExNeuron) of ExNeuron DEGs: up-regulated, down-regulated and all DEGs. Venn = significant terms in CES and CPS; bars = -log10 adj. P of shared terms (CES left, CPS right)."),
 ("1G", "figure1G_disease_gene_overlap", "FIGS", "Fig1G_disease_chord", "13_figure1EG.R", "table6_disease_overlap; reference/GWAS_Neuron_Disease_Gene_Dataset.csv",
  "Overlap of SDS vs CON DEGs (adj. P < 0.05, |log2FC| > 0.25) of each major type with psychiatric GWAS gene sets (7 traits; human genes mapped to mouse with homologene). Numbers after cell types = DEGs; ribbon width = overlapping genes. ExNeuron links are highlighted (dark outline); observed and expected overlaps per trait are in table6. Other ribbons: dark = enriched (one-sided hypergeometric, expressed-gene background, BH FDR < 0.05) and still enriched without neuron-enriched genes; light = enriched only through neuron-enriched genes; grey = not enriched."),
 ("2A", "figure2A_umap_neuron_subtypes", "FIGS", "Fig2A_umap_neuron_subtypes", "14_figure2.py; 4_build_objects.py", "ACC_neurons.h5ad",
  "UMAP of {n_neurons} neurons coloured by subtype: 12 ExNeuron and 5 InNeuron subtypes (cluster-level annotation from marker expression)."),
 ("2B", "figure2B_composition_neuron_subtypes", "FIGS", "Fig2B_composition_neuron_subtypes", "14_figure2.py", "table2_cell_composition",
  "Neuron subtype composition per animal (fraction of all neurons in each animal; x-axis labels coloured by group). The replicate-level test (logit proportion, Welch t-test vs CON, n = 3 vs 3, BH) is in table2."),
 ("2C", "figure2C_DEG_counts_neuron_subtypes", "FIGS", "Fig2C_DEG_counts_neuron_subtypes", "14_figure2.py; 5_de_pseudobulk_export.py", "table3_DEGs",
  "Top: nuclei per subtype in CON, CES and CPS (log10 scale). Below: number of up- and down-regulated DEGs per neuron subtype for CES vs CON (upper panel) and CPS vs CON (lower panel, bars plotted downwards, axis shows counts); dots on the zero line = subtype colours as in Fig. 2A."),
 ("2D", "figure2D_DEG_jaccard_CES_CPS", "FIGS", "Fig2D_DEG_jaccard_CES_CPS", "14_figure2.py", "table3_DEGs",
  "Jaccard index of CES and CPS DEG sets per subtype (bar labels = shared DEGs). Cdh9+ ExNeurons have the highest overlap ({cdh9_jaccard}; {cdh9_shared} shared DEGs)."),
 ("2E", "figure2E_Cdh9_DEG_venn", "FIGS", "Fig2E_Cdh9_DEG_venn", "14_figure2.py", "table3_DEGs",
  "Cdh9+ DEGs in CES vs CON and CPS vs CON ({cdh9_ces_only} CES-only, {cdh9_shared} shared, {cdh9_cps_only} CPS-only)."),
 ("2F", "figure2F_Cdh9_GO_overlap", "FIGS", "Fig2F_Cdh9_GO_overlap", "14_figure2.py; 6_stats.R", "table5_GO_enrichment",
  "GO enrichment of Cdh9+ DEGs (as Fig. 1F): up-regulated, down-regulated and all DEGs (up + down). Venn = significant terms in CES and CPS; bars = -log10 adj. P of shared terms."),
 ("3A", "figure3A_Cdh9_shared_DEG_GO", "FIGS", "Fig3A_Cdh9_shared_DEG_GO", "6_stats.R; 15_figure3.py", "table4_Cdh9_shared_DEGs",
  "GO enrichment of Cdh9+ DEGs shared by CES and CPS with the same sign (left: common up, n = {n_up}; right: common down, n = {n_down}). Significant terms are clustered by gene overlap (Jaccard, average linkage, cut 0.7); one representative per cluster (+n = merged terms). x = fold enrichment, size = genes, colour = -log10 adj. P. Matrix: genes behind the terms, colour = mean log2FC (CES, CPS)."),
 ("3B", "figure3B_Cdh9_shared_DEG_log2FC", "FIGS", "Fig3B_Cdh9_shared_DEG_log2FC_scatter", "15_figure3.py", "table4_Cdh9_shared_DEGs",
  "log2FC of the {cdh9_shared} shared Cdh9+ DEGs in CES vs CON and CPS vs CON ({n_up} common up, {n_down} common down, {n_opp} opposite; Pearson r = {r_shared})."),
 ("3C", "figure3C_cellchat_Cdh9", "FIGS", "Fig3C_cellchat_Cdh9_diff_count_weight", "15_figure3.py; 7_cellchat_export.py; 8_cellchat.R", "cellchat/cellchat_fine_*.rds; table7_CellChat",
  "CellChat v2 differential number (left) and strength (right) of interactions, CES vs CON (top) and CPS vs CON (bottom), edges to and from Cdh9+ (fine groups). Edge width = |difference| (same scale for CES and CPS); red = increase, blue = decrease. Node area and opacity = interactions of the partner with Cdh9+ (mean of CON and stress group)."),
 ("3C_merged", "figure3C_merged_cellchat_Cdh9_merged_groups", "FIGS", "Fig3C_merged_cellchat_Cdh9_diff_count_weight", "15_figure3.py; 7_cellchat_export.py; 8_cellchat.R", "cellchat/cellchat_coarse_*.rds; table7_CellChat",
  "As Fig. 3C, with CellChat run on merged groups: Cdh9+, other ExNeuron subtypes (OtherExNeuron), all InNeuron and the glial/vascular types."),
 ("3D", "figure3D_progeny_all_celltypes", "FIGS", "Fig3D_progeny_common_all_celltypes", "15_figure3.py; 9_progeny.py", "table8_PROGENy",
  "PROGENy pathway activities (mouse model, top 500 genes; univariate linear model on per-animal mean expression of groups with >= 20 nuclei): top 20 features changed in the same direction in CES and CPS. Heatmap = row z-score per animal; right = mean activity difference vs CON; Welch t-test, BH."),
 ("3E", "figure3E_progeny_Cdh9", "FIGS", "Fig3E_progeny_Cdh9", "15_figure3.py; 9_progeny.py", "table8_PROGENy",
  "PROGENy pathways changed in the same direction in CES and CPS in Cdh9+ ExNeurons."),
 ("3F", "figure3F_Cdh9_gene_violins", "FIGS", "Fig3F_Cdh9_gene_violins", "15_figure3.py", "ACC_neurons.h5ad; table3_DEGs",
  "Expression of Camk1d, Grik3, Iqsec1 and Shank2 in Cdh9+ nuclei per group (violins = nuclei, white dots = animal means). Brackets: CES and CPS vs CON, cell-level Welch t-test, BH adj. P: * < 0.05, ** < 0.01, *** < 0.001, ns = not significant."),
 ("S1", "figureS1_split_threshold_per_animal", "QC", "kde_split_all_animals", "1_qc_split.py", "pre-filter object; table1_samples_QC",
  "Low/high-feature split per animal: density of genes per nucleus (Gaussian KDE) with the detected peaks and the split threshold (red; deepest minimum between the glial peak and the tallest peak above it, below 2x the glial peak). CPS_3 has no valid valley and uses the median of the other animals ({cps3_thr} genes). Text = fraction of nuclei below and above the threshold."),
 ("S2", "figureS2_ExNeuron_NMF_metaprograms", "FIGS", "FigS2_ExNeuron_NMF_metaprograms", "10_nmf_programs.py; 11_nmf_go.R; 16_nmf_scores_figureS2.py", "ACC_neurons.h5ad; table9_NMF",
  "NMF gene programs in ExNeurons. (A) {n_robust} robust factors (per-animal NMF, k = 4-9, Gavish et al. filters; heatmap = shared top-50 genes) clustered into {n_mp} meta-programs (MPs); lines connect each MP block to its function (GO of the top 30 genes; [subtype with the highest score; animals contributing]) and top 10 genes. (B) GO enrichment of MP genes (ring = also significant with the HVG universe). (C) MP activity per subtype (z-score). (D) CES/CPS vs CON per subtype, animal means, Welch t-test (ring = FDR < 0.1)."), ("S3", "figureS3_neuron_subtype_markers", "FIGS", "FigS3_neuron_subtype_markers", "14_figure2.py", "ACC_neurons.h5ad",
  "Marker expression of neuron subtypes (matrix plot). Rows = 12 ExNeuron (top) and 5 InNeuron (bottom) subtypes; columns = lineage markers (Slc17a7, Gad1, Gad2) and the genes that name each subtype. Colour = mean log-normalised expression per subtype, min-max scaled per gene across all neurons; grey = ExNeuron subtype markers in InNeuron subtypes and InNeuron subtype markers in ExNeuron subtypes (not shown; values in the source data)."),
]
# extra source-data files (beyond <stem>_source_data*.csv) per panel
EXTRA_SOURCE = {"1D": [("Fig1D_DEG_counts_per_celltype.csv", "DEG_counts")], "3C": [("Fig3C_cellchat_total_count_weight.csv", "total_count_weight")],
                "S1": [("QC:kde_split_summary.csv", "")]}

# table id -> (title, statistical unit / model, scripts, sheets [(sheet, source, note)])
# source prefixes: D: = DATA, P: = PIPE, F: = FIGS, Q: = QC, W: = WORK, X: = built in 17_assemble_outputs.py
TABLES = {
 "table1_samples_QC": ("Samples, QC and split thresholds", "animal", "1_qc_split.py; 17_assemble_outputs.py", [
    ("samples", "X:samples", "one row per animal: IDs, group, nuclei before/after QC, Scrublet doublets, split threshold, final nuclei"),
    ("qc_medians", "X:qc_medians", "median UMIs, genes and mitochondrial % per nucleus before and after QC"),
    ("exclusions", "X:exclusions", "nuclei per animal and status (final or exclusion reason)"),
    ("gene_cap_loss", "X:gene_cap", "nuclei above the 4,000-gene and 15,000-UMI caps per animal"),
    ("glia_in_mixed_clusters", "X:mixed", "glia-labelled nuclei in mixed (X_ExNeuron) clusters vs kept, per animal"),
    ("split_kde", "Q:kde_split_summary.csv", "KDE peaks, threshold and rule per animal")]),
 "table2_cell_composition": ("Cell-type composition", "animal; neuron subtypes: logit(proportion + 0.5) Welch t-test (n = 3 vs 3), BH across subtypes", "12_figure1.py; 13_figure1EG.R; 14_figure2.py", [
    ("major_counts_per_animal", "X:major_counts", "final nuclei per animal and major type"),
    ("subtype_counts_per_animal", "X:subtype_counts", "final neurons per animal and subtype"),
    ("subtype_test", "F:Fig2B_composition_neuron_subtypes_source_data_logit_test.csv", "Fig. 2B"),
    ("pooled_fractions", "F:Fig1E_composition_alluvial_source_data.csv", "Fig. 1E (descriptive)")]),
 "table3_DEGs": ("Differential expression", "nucleus: scanpy t-test (Welch), BH per group x contrast; animal: limma-voom pseudobulk (n = 3 vs 3), BH", "5_de_pseudobulk_export.py; 6_stats.R", [
    ("DEG_summary", "X:deg_summary", "DEGs (adj. P < 0.05, |log2FC| > 0.25) per group and contrast, nuclei per arm, pseudobulk genes adj. P < 0.05"),
    ("DEGs_major", "X:sig_major", "significant cell-level DEGs, major types, CES/CPS/SDS vs CON, with pseudobulk logFC/P of the same gene"),
    ("DEGs_subtype", "X:sig_subtype", "as above, neuron subtypes and glial types")]),
 "table4_Cdh9_shared_DEGs": ("Cdh9+ DEGs shared by CES and CPS and their GO terms", "gene; GO ORA: hypergeometric, BH; GSEA: Cdh9+ pseudobulk t ranks, BH", "6_stats.R; 15_figure3.py", [
    ("shared_DEGs", "D:Cdh9_shared_DEGs.csv", "class = Common up / Common down / Opposite direction"),
    ("GO_shared", "D:GO_Cdh9_shared_DEGs.csv", "Fig. 3A ORA, expressed-gene universe"),
    ("GO_shared_allgene_bg", "D:GO_Cdh9_shared_DEGs_allgene_background.csv", "sensitivity: all tested genes as universe"),
    ("Fig3A_terms", "F:Fig3A_Cdh9_shared_DEG_GO_source_data.csv", "Fig. 3A representative terms"),
    ("GSEA_CES_CPS", "D:fig3a_followup/gsea_pseudobulk_Cdh9_CES_CPS_wide.csv", "gseGO on Cdh9+ pseudobulk t; concordant = BH < 0.05 in both, same sign"),
    ("ORA_terms_in_GSEA", "D:fig3a_followup/fig3a_ora_terms_in_gsea.csv", "Fig. 3A ORA terms with their GSEA NES")]),
 "table5_GO_enrichment": ("GO enrichment per contrast", "gene; hypergeometric (enrichGO), BH, simplify 0.7", "6_stats.R; 12_figure1.py; 14_figure2.py", [
    ("ExNeuron_Cdh9_per_contrast", "D:GO_per_contrast_ExNeuron_Cdh9.csv", "up / down / all DEGs, CES and CPS vs CON, expressed-gene universe"),
    ("allgene_background", "D:GO_per_contrast_ExNeuron_Cdh9_allgene_background.csv", "sensitivity: all tested genes as universe"),
    ("Fig1F", "F:Fig1F_ExNeuron_GO_overlap_source_data.csv", "Fig. 1F"),
    ("Fig2F", "F:Fig2F_Cdh9_GO_overlap_source_data.csv", "Fig. 2F")]),
 "table6_disease_overlap": ("Overlap of SDS vs CON DEGs with psychiatric GWAS gene sets", "gene; one-sided hypergeometric, BH", "13_figure1EG.R", [
    ("SDS_vs_CON", "F:Fig1G_disease_chord_source_data.csv", "Fig. 1G")]),
 "table7_CellChat": ("Cell-cell communication (CellChat v2)", "balanced nuclei per animal; permutation test (100), P < 0.05", "7_cellchat_export.py; 8_cellchat.R; 15_figure3.py", [
    ("Fig3C_Cdh9_edges", "F:Fig3C_cellchat_Cdh9_diff_count_weight_source_data.csv", "Fig. 3C"),
    ("Fig3C_merged_edges", "F:Fig3C_merged_cellchat_Cdh9_diff_count_weight_source_data.csv", "Fig. 3C merged groups"),
    ("fine_net", "X:cc_net_fine", "interaction count and weight per source x target, CON/CES/CPS (fine groups)"),
    ("coarse_net", "X:cc_net_coarse", "as above, merged groups"),
    ("fine_Cdh9_LR", "X:cc_lr_fine", "significant L-R interactions with Cdh9+ as source or target, per condition"),
    ("cells_used", "X:cc_balance", "nuclei per group and animal used")]),
 "table8_PROGENy": ("Pathway activity (PROGENy)", "animal (n = 3 vs 3); Welch t-test, BH within contrast", "6_stats.R; 9_progeny.py; 15_figure3.py", [
    ("differential", "D:progeny_differential.csv", "CES/CPS vs CON per group x pathway"),
    ("activity_per_animal", "D:progeny_activity_long.csv", "ULM activity per animal, group, pathway"),
    ("Fig3D", "F:Fig3D_progeny_common_all_celltypes_source_data.csv", "Fig. 3D"),
    ("Fig3E", "F:Fig3E_progeny_Cdh9_source_data.csv", "Fig. 3E")]),
 "table9_NMF": ("ExNeuron NMF meta-programs", "factors per animal; MP scores per animal x subtype, Welch t-test, BH within contrast", "10_nmf_programs.py; 11_nmf_go.R; 16_nmf_scores_figureS2.py", [
    ("MP_genes", "D:nmf/metaprogram_genes.csv", "top 30 genes per MP"),
    ("MP_annotation", "D:nmf/metaprogram_annotation.csv", "functional names"),
    ("MP_GO_expressed", "D:nmf/metaprogram_GO_expressed_universe.csv", "GO of MP genes, expressed-gene universe"),
    ("MP_GO_hvg", "D:nmf/metaprogram_GO_hvg_universe.csv", "GO of MP genes, HVG universe"),
    ("MP_group_test", "D:nmf/metaprogram_group_test.csv", "CES/CPS vs CON per subtype"),
    ("MP_scores_animal", "D:nmf/metaprogram_scores_per_animal.csv", "mean score per animal x subtype"),
    ("robust_factors", "D:nmf/robust_factors.csv", "factors kept and their MP")]),
 "table10_annotation": ("Clustering and annotation evidence", "cluster", "2_cluster_annotate_major.py; 3_neuron_subtypes.py; 4_build_objects.py", [
    ("markers_low_feature", "P:markers_low_feature.csv", "per Leiden cluster: size, dominant animal, median genes, lineage scores, marker detection"),
    ("markers_high_feature", "P:markers_high_feature.csv", "as above, high-feature object"),
    ("markers_ExNeuron_pass1", "P:markers_ExNeuron_pass1.csv", "ExNeuron pass-1 clusters"),
    ("naming_ExNeuron", "P:naming_ExNeuron_pass2.csv", "name-gene expression + top Wilcoxon markers"),
    ("naming_ExNeuron_c0", "P:naming_ExNeuron_cluster0_subclusters.csv", "Trpc6 vs Col19a1 in cluster-0 subclusters"),
    ("naming_InNeuron", "P:naming_InNeuron.csv", "name-gene expression + top Wilcoxon markers"),
    ("lineage_check", "X:lineage", "nuclei removed by the lineage check, per animal and type"),
    ("annotation_check", "W:logs/misplaced_check.csv", "UMAP-neighbour and lineage-marker check per type and animal")]),
}
