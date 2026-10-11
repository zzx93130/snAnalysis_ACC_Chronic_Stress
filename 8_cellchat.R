# Step 8 - CellChat v2 per condition on balanced cells (parameters of script/cellchat_split3.R):
# CellChatDB.mouse, triMean, population.size = FALSE, nboot = 100, filterCommunication(min.cells = 30).
# The ligand-receptor set is the union of over-expressed interactions identified in each
# biological sample (as in cellchat_split3.R section 7); networks are then inferred on cells
# pooled per condition. Runs the coarse and fine groupings.
# Run with the 'ccc' R 4.5 environment and .libPaths() including envs/ccc-rlib (CellChat 2.2.0.9001).

source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1])), "config.R"))
.libPaths(c(CCC_RLIB, .libPaths()))
suppressMessages({library(CellChat); library(Matrix); library(future)})
for (cfg in c("coarse", "fine")) {
dir <- file.path("data", paste0("cellchat_", cfg))
seed_use <- 20260724; set.seed(seed_use)
options(future.globals.maxSize = 8 * 1024^3)
options(parallelly.maxWorkers.localhost = c(16, 16)); plan("multicore", workers = 6)

X <- readMM(file.path(dir, "matrix.mtx")); X <- as(X, "CsparseMatrix")
rownames(X) <- readLines(file.path(dir, "genes.txt"))
meta <- read.csv(file.path(dir, "meta.csv"), row.names = 1); colnames(X) <- rownames(meta)
levels_use <- sort(unique(meta$cell_group))
meta$cell_group <- factor(meta$cell_group, levels = levels_use)
DB <- CellChatDB.mouse

prep <- function(cells) {
  m <- meta[cells, , drop = FALSE]; m$cell_group <- droplevels(m$cell_group)
  cc <- createCellChat(object = X[, cells], meta = m, group.by = "cell_group")
  cc@DB <- DB; cc <- subsetData(cc)
  cc <- identifyOverExpressedGenes(cc); identifyOverExpressedInteractions(cc)
}
lr <- unique(unlist(lapply(split(rownames(meta), meta$sample), function(cells) rownames(prep(cells)@LR$LRsig))))
lr <- intersect(lr, rownames(DB$interaction)); common_lr <- DB$interaction[lr, , drop = FALSE]
write.csv(data.frame(interaction_name = lr), file.path(dir, "shared_LR_pairs.csv"), row.names = FALSE)
message("shared LR pairs: ", length(lr))

objs <- list()
for (cond in c("CON", "CES", "CPS")) {
  cells <- rownames(meta)[meta$condition == cond]
  m <- meta[cells, , drop = FALSE]
  cc <- createCellChat(object = X[, cells], meta = m, group.by = "cell_group")
  cc@DB <- DB; cc <- subsetData(cc)
  cc <- identifyOverExpressedGenes(cc); cc <- identifyOverExpressedInteractions(cc)
  cc@LR$LRsig <- common_lr
  cc <- computeCommunProb(cc, type = "triMean", LR.use = common_lr, population.size = FALSE,
                          nboot = 100, seed.use = seed_use)
  cc <- filterCommunication(cc, min.cells = 30)
  cc <- computeCommunProbPathway(cc, thresh = 0.05)
  cc <- aggregateNet(cc, thresh = 0.05)
  saveRDS(cc, file.path(dir, paste0(cond, "_CellChat.rds")))
  write.csv(cc@net$count, file.path(dir, paste0(cond, "_net_count.csv")))
  write.csv(cc@net$weight, file.path(dir, paste0(cond, "_net_weight.csv")))
  write.csv(subsetCommunication(cc, thresh = 0.05), file.path(dir, paste0(cond, "_communications.csv")), row.names = FALSE)
  objs[[cond]] <- cc; message(cond, " done")
}
for (cmp in list(c("CON", "CES"), c("CON", "CPS"))) {
  mg <- mergeCellChat(objs[cmp], add.names = cmp, cell.prefix = TRUE)
  saveRDS(mg, file.path(dir, sprintf("merged_%s_%s.rds", cmp[1], cmp[2])))
}
writeLines(capture.output(sessionInfo()), file.path("logs", paste0("8_cellchat_", cfg, "_sessionInfo.txt")))
}
