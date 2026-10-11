# NMF gene programs in ExNeuron (step 3a): GO enrichment of the top-30 genes of every meta-program.
# clusterProfiler enrichGO, ont = ALL, BH, pvalueCutoff 0.05, qvalueCutoff 0.2, minGSSize 10, simplify 0.7.
# Universe (main): the 2,000 HVGs the NMF could select from (an MP gene can only come from them);
# sensitivity: genes detected in >= 1% of ExNeuron nuclei (the universe used elsewhere in this figure set).
suppressMessages({library(clusterProfiler); library(org.Mm.eg.db); library(data.table)})
source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1])), "config.R")); setwd(file.path(WORK, "data", "nmf"))
mp <- fread("metaprogram_genes.csv")
run <- function(univ_file, tag) {
  univ <- readLines(univ_file)
  res <- rbindlist(lapply(names(mp), function(m) {
    e <- enrichGO(mp[[m]], OrgDb = org.Mm.eg.db, keyType = "SYMBOL", ont = "ALL", universe = univ,
                  pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2, minGSSize = 10, maxGSSize = 500)
    if (is.null(e) || nrow(as.data.frame(e)) == 0) return(NULL)
    e <- tryCatch(simplify(e, cutoff = 0.7, by = "p.adjust", select_fun = min), error = function(err) e)
    d <- as.data.table(as.data.frame(e)); d[, MP := m]; d
  }), fill = TRUE)
  # also keep the best nominal terms for MPs without significant terms (annotation aid, flagged)
  nom <- rbindlist(lapply(names(mp), function(m) {
    e <- enrichGO(mp[[m]], OrgDb = org.Mm.eg.db, keyType = "SYMBOL", ont = "ALL", universe = univ,
                  pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1, minGSSize = 10, maxGSSize = 500)
    d <- as.data.table(as.data.frame(e))[order(pvalue)][1:5]; d[, MP := m]; d }), fill = TRUE)
  fwrite(res, sprintf("metaprogram_GO_%s.csv", tag)); fwrite(nom, sprintf("metaprogram_GO_top5_nominal_%s.csv", tag))
  cat("==", tag, "\n"); print(res[, .N, by = MP])
  print(res[order(MP, p.adjust)][, head(.SD, 4), by = MP][, .(MP, ONTOLOGY, Description, GeneRatio, p.adjust = signif(p.adjust, 2))])
}
run("go_universe_hvg.txt", "hvg_universe")
run("go_universe_expressed.txt", "expressed_universe")
