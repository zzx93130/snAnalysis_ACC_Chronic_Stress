# Step 6 - R statistics: PROGENy network export, pseudobulk limma-voom, GO over-representation (expressed-gene
# and all-gene universes), GSEA of Cdh9+ pseudobulk t statistics and the representative Fig. 3A GO terms.
source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1])), "config.R"))


# ====================================================================================================
# 6a PROGENy mouse model export
# ====================================================================================================
# Export the PROGENy mouse model (top 500 genes per pathway by p-value) in decoupler long format.
suppressMessages(library(progeny))
m <- progeny::getModel(organism = "Mouse", top = 500)
long <- do.call(rbind, lapply(colnames(m), function(p) {
  w <- setNames(m[, p], rownames(m)); w <- w[w != 0]
  data.frame(source = p, target = names(w), weight = unname(w))
}))
write.csv(long, gzfile("data/progeny_mouse_top500.csv.gz"), row.names = FALSE)
cat("progeny", as.character(packageVersion("progeny")), nrow(long), length(unique(long$source)), "\n")


# ====================================================================================================
# 6b pseudobulk limma-voom
# ====================================================================================================
# Step 12 - replicate-aware check: pseudobulk limma-voom (biological sample = unit, n = 3 per group).
# For each cell group with >= 10 cells in every sample of the two arms, sum counts per sample,
# filterByExpr, TMM, voom, lmFit(~0 + condition), contrasts CES-CON and CPS-CON.
# Output: data/DE_pseudobulk_{major,subtype}.csv.gz
suppressMessages({library(limma); library(edgeR); library(data.table)})
MIN_CELLS <- 10
for (lvl in c("major", "subtype")) {
  cnt <- as.matrix(data.frame(fread(sprintf("data/pseudobulk_counts_%s.csv.gz", lvl)), row.names = 1, check.names = FALSE))
  meta <- fread(sprintf("data/pseudobulk_meta_%s.csv", lvl))
  out <- list()
  for (g in unique(meta$group)) {
    m <- meta[group == g & n_cells >= MIN_CELLS]
    if (any(table(factor(m$condition, levels = c("CON", "CES", "CPS"))) < 2)) next
    y <- DGEList(cnt[, m$id, drop = FALSE])
    cond <- factor(m$condition, levels = c("CON", "CES", "CPS"))
    design <- model.matrix(~0 + cond); colnames(design) <- levels(cond)
    keep <- filterByExpr(y, design); y <- calcNormFactors(y[keep, , keep.lib.sizes = FALSE])
    v <- voom(y, design)
    fit <- eBayes(contrasts.fit(lmFit(v, design), makeContrasts(CES - CON, CPS - CON, levels = design)))
    for (k in 1:2) {
      tt <- topTable(fit, coef = k, number = Inf, sort.by = "none")
      out[[paste(g, k)]] <- data.frame(cluster = g, condition = c("CES", "CPS")[k], reference = "CON",
        gene = rownames(tt), logFC = tt$logFC, AveExpr = tt$AveExpr, t = tt$t, P.Value = tt$P.Value,
        adj.P.Val = tt$adj.P.Val, n_samples = nrow(m))
    }
  }
  res <- rbindlist(out)
  fwrite(res, sprintf("data/DE_pseudobulk_%s.csv.gz", lvl))
  cat(lvl, nrow(res), length(unique(res$cluster)), "groups\n")
}


# ====================================================================================================
# 6c GO over-representation (run for both gene universes)
# ====================================================================================================
# Step 13 - GO over-representation (clusterProfiler::enrichGO, org.Mm.eg.db, ont = "ALL").
# Parameters follow script/function_enrichment20260325.R (BH, simplify cutoff 0.7 Wang) and
# script/Cdh9_neuron_related_analysis.R (shared concordant Cdh9+ DEGs).
#  A. Per-contrast GO for ExNeuron (major) and Cdh9+ (subtype): up, down, and all DEGs (up + down combined);
#     DEG lists exclude ^(mt-|Rpl|Rps|Hsp) as in function_enrichment20260325.R.
#  B. GO of Cdh9+ DEGs shared by CES and CPS with concordant direction (common up / common down);
#     no gene exclusion, as in Cdh9_neuron_related_analysis.R.
# Usage: Rscript 13_go_enrichment.R [expressed|all]
# Background = genes detected in >= 1% of cells of that group in either arm (both contrasts for B).
suppressMessages({library(clusterProfiler); library(org.Mm.eg.db); library(data.table); library(jsonlite)})
P <- fromJSON(PARAMS)
LFC <- P$deg_log2fc; PADJ <- P$deg_padj; MINPCT <- P$background_min_pct
excl <- "^(mt-|Rpl|Rps|Hsp)"

run_go <- function(genes, universe) {
  genes <- intersect(unique(genes), universe)
  if (length(genes) < 5) return(NULL)
  ego <- tryCatch(enrichGO(genes, OrgDb = org.Mm.eg.db, keyType = "SYMBOL", ont = "ALL",
                           universe = universe, pAdjustMethod = "BH", pvalueCutoff = 0.05,
                           qvalueCutoff = 0.2, minGSSize = 10, maxGSSize = 500, readable = FALSE),
                  error = function(e) NULL)
  if (is.null(ego) || nrow(as.data.frame(ego)) == 0) return(NULL)
  ego <- tryCatch(simplify(ego, cutoff = 0.7, by = "p.adjust", select_fun = min, measure = "Wang"),
                  error = function(e) ego)
  d <- as.data.frame(ego); d$n_input <- length(genes); d$n_universe <- length(universe); d
}
for (BG_MODE in c("expressed", "all")) {
# expressed (main): genes detected in >= 1% of cells; all: every gene tested (original scripts' choice)
bg_of <- function(d) if (BG_MODE == "all") unique(d$gene) else unique(d[pct_cond >= MINPCT | pct_ref >= MINPCT, gene])
SUF <- if (BG_MODE == "all") "_allgene_background" else ""

major <- fread("data/DE_celllevel_major.csv.gz"); sub <- fread("data/DE_celllevel_subtype.csv.gz")
res <- list()
for (spec in list(list(d = major, cl = "ExNeuron"), list(d = sub, cl = "Cdh9+"))) {
  for (cond in c("CES", "CPS")) {
    d <- spec$d[cluster == spec$cl & condition == cond & reference == "CON"]
    uni <- bg_of(d)
    for (dir in c("up", "down", "all")) {
      s <- d[pvals_adj < PADJ & (if (dir == "up") log2FC > LFC else if (dir == "down") log2FC < -LFC else abs(log2FC) > LFC)]
      s <- s[!grepl(excl, gene, ignore.case = TRUE)]
      g <- run_go(s$gene, uni)
      if (!is.null(g)) res[[paste(spec$cl, cond, dir)]] <- data.table(cluster = spec$cl, condition = cond,
                                                                    reference = "CON", direction = dir, g)
      cat(spec$cl, cond, dir, nrow(s), "DEGs ->", ifelse(is.null(g), 0, nrow(g)), "terms\n")
    }
  }
}
fwrite(rbindlist(res, fill = TRUE), paste0("data/GO_per_contrast_ExNeuron_Cdh9", SUF, ".csv"))

# ---- B. shared concordant Cdh9+ DEGs ---------------------------------------
d <- sub[cluster == "Cdh9+" & reference == "CON" & condition %in% c("CES", "CPS")]
uni <- Reduce(intersect, lapply(split(d, d$condition), bg_of))
sig <- d[pvals_adj < PADJ & abs(log2FC) > LFC]
w <- dcast(sig, gene ~ condition, value.var = "log2FC")
w <- w[!is.na(CES) & !is.na(CPS)]
w[, class := fifelse(CES > 0 & CPS > 0, "Common up", fifelse(CES < 0 & CPS < 0, "Common down", "Opposite direction"))]
if (SUF == "") fwrite(w, "data/Cdh9_shared_DEGs.csv")
resB <- list()
for (cl in c("Common up", "Common down")) {
  g <- run_go(w[class == cl, gene], uni)
  if (!is.null(g)) resB[[cl]] <- data.table(set = cl, g)
  cat(cl, sum(w$class == cl), "genes ->", ifelse(is.null(g), 0, nrow(g)), "terms\n")
}
fwrite(rbindlist(resB, fill = TRUE), paste0("data/GO_Cdh9_shared_DEGs", SUF, ".csv"))

}


# ====================================================================================================
# 6d GSEA (Cdh9+ pseudobulk t) and representative Fig. 3A terms
# ====================================================================================================
dir.create("data/fig3a_followup", showWarnings = FALSE)
suppressMessages({library(clusterProfiler); library(org.Mm.eg.db); library(data.table); library(limma); library(edgeR)
  library(homologene); library(readxl); library(jsonlite)})
OUT <- "data/fig3a_followup"; P <- fromJSON(PARAMS); LFC <- P$deg_log2fc; PADJ <- P$deg_padj; MINPCT <- P$background_min_pct
set.seed(1)
sh <- fread("data/Cdh9_shared_DEGs.csv")

# ---- A. GSEA on pseudobulk t ------------------------------------------------
pb <- fread("data/DE_pseudobulk_subtype.csv.gz")[cluster == "Cdh9+"]
gsea_file <- file.path(OUT, "gsea_pseudobulk_Cdh9_all.csv")
gs <- rbindlist(lapply(c("CES", "CPS"), function(cond) {
  d <- pb[condition == cond]; r <- sort(setNames(d$t, d$gene), decreasing = TRUE)
  g <- gseGO(r, OrgDb = org.Mm.eg.db, keyType = "SYMBOL", ont = "ALL", minGSSize = 15, maxGSSize = 500,
             pvalueCutoff = 1, pAdjustMethod = "BH", eps = 0, seed = TRUE, verbose = FALSE)
  data.table(condition = cond, as.data.frame(g))
}), fill = TRUE)
fwrite(gs, gsea_file)
w <- dcast(gs[, .(condition, ID, Description, ONTOLOGY, NES, p.adjust, setSize)], ID + Description + ONTOLOGY + setSize ~ condition,
           value.var = c("NES", "p.adjust"))
w[, concordant := p.adjust_CES < 0.05 & p.adjust_CPS < 0.05 & sign(NES_CES) == sign(NES_CPS)]
fwrite(w, file.path(OUT, "gsea_pseudobulk_Cdh9_CES_CPS_wide.csv"))
ora <- fread("data/GO_Cdh9_shared_DEGs.csv")[p.adjust < 0.05]
ora_g <- merge(ora[, .(set, ID, Description, ONTOLOGY, ora_padj = p.adjust)], w, by = c("ID", "Description", "ONTOLOGY"), all.x = TRUE)
fwrite(ora_g, file.path(OUT, "fig3a_ora_terms_in_gsea.csv"))
cat("GSEA terms:", nrow(w), "concordant:", sum(w$concordant, na.rm = TRUE),
    "(up", sum(w$concordant & w$NES_CES > 0, na.rm = TRUE), "/ down", sum(w$concordant & w$NES_CES < 0, na.rm = TRUE), ")\n")
print(ora_g[, .(n = .N, in_gsea = sum(!is.na(NES_CES)), concordant = sum(concordant, na.rm = TRUE)), by = set])

# representative terms: significant ORA terms clustered by gene overlap (Jaccard of term genes, average
# linkage, cut at distance 0.7, i.e. terms sharing on average >= 30% of their genes are merged); representative = lowest adj. P; n_related = other terms in the cluster.
reps <- rbindlist(lapply(c("Common up", "Common down"), function(cl) {
  o <- ora[set == cl][order(p.adjust)]; gl <- strsplit(o$geneID, "/")
  J <- outer(seq_along(gl), seq_along(gl), Vectorize(function(i, j) length(intersect(gl[[i]], gl[[j]])) / length(union(gl[[i]], gl[[j]]))))
  k <- if (nrow(o) > 1) cutree(hclust(as.dist(1 - J), "average"), h = 0.7) else 1
  o[, cluster := k]; o[, rep := !duplicated(cluster)]
  o[, `:=`(n_related = .N - 1L, related_terms = paste(Description[-1], collapse = "; ")), by = cluster]
  o[rep == TRUE]
}))
fwrite(reps, file.path(OUT, "fig3a_representative_terms.csv"))
cat("representative terms:", paste(reps[, .N, by = set][, paste(set, N)], collapse = ", "), "\n")
writeLines(capture.output(sessionInfo()), "logs/6_stats_sessionInfo.txt")
