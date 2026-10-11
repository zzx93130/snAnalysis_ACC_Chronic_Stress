# Step 13 - Figure 1E (composition alluvial) and Figure 1G (disease gene-set chord), R.
source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1])), "config.R"))


# ====================================================================================================
# 13a Figure 1E
# ====================================================================================================
# Figure 1E - major cell-type composition as a 5-bar alluvial plot (style of the user's plot_5bar_flow):
#   CON | SDS   (gap)   CON | CES | CPS
# Bars = nuclei pooled per group (fraction of all nuclei in that group); ribbons connect the same cell type.
# Colours = project palette (palettes.json, same as Fig1A). Pooled fractions are descriptive.
suppressMessages({library(dplyr); library(ggplot2); library(ggalluvial); library(scales); library(jsonlite); library(readr)})
pal <- unlist(fromJSON("palettes.json")$major_celltype)
md <- read_csv("data/ACC_cell_metadata.csv", show_col_types = FALSE)
ct_levels <- names(pal)[names(pal) %in% md$major_celltype]
md$major_celltype <- factor(md$major_celltype, levels = ct_levels)

axis_levels <- c("CON_1", "SDS", "gap", "CON_2", "CES", "CPS")
frac <- function(group, levels, axes) md %>% filter(.data[[group]] %in% levels) %>%
  count(Group = .data[[group]], Celltype = major_celltype, name = "Freq", .drop = FALSE) %>%
  group_by(Group) %>% mutate(Frac = Freq / sum(Freq), n_group = sum(Freq)) %>% ungroup() %>%
  mutate(Axis = factor(unname(setNames(axes, levels)[Group]), levels = axis_levels))
d1 <- frac("Experiment2", c("CON", "SDS"), c("CON_1", "SDS"))
d2 <- frac("Experiment", c("CON", "CES", "CPS"), c("CON_2", "CES", "CPS"))
dat <- bind_rows(mutate(d1, comparison = "CON vs SDS"), mutate(d2, comparison = "CON vs CES vs CPS"))
write_csv(dat %>% select(comparison, Group, Axis, Celltype, Freq, n_group, Frac), "figures/Fig1E_composition_alluvial_source_data.csv")

nlab <- dat %>% distinct(Axis, n_group) %>% mutate(lab = comma(n_group))
p <- ggplot() +
  geom_alluvium(data = d1, aes(x = Axis, y = Frac, alluvium = Celltype, fill = Celltype), alpha = 0.55, width = 0.45, color = NA) +
  geom_alluvium(data = d2, aes(x = Axis, y = Frac, alluvium = Celltype, fill = Celltype), alpha = 0.55, width = 0.45, color = NA) +
  geom_stratum(data = dat, aes(x = Axis, y = Frac, stratum = Celltype, fill = Celltype), width = 0.45, color = "white", linewidth = 0.3) +
  geom_text(data = nlab, aes(x = Axis, y = 1.03, label = lab), size = 1.9, vjust = 0, color = "grey30") +
  scale_x_discrete(limits = axis_levels, drop = FALSE,
                   labels = c(CON_1 = "CON", SDS = "SDS", gap = "", CON_2 = "CON", CES = "CES", CPS = "CPS")) +
  scale_y_continuous(labels = percent_format(), breaks = seq(0, 1, 0.25), expand = c(0, 0)) +
  scale_fill_manual(values = pal[ct_levels], breaks = ct_levels) +
  labs(x = NULL, y = "Cell proportion", fill = NULL) +
  annotate("segment", x = 0.7, xend = 2.3, y = -0.105, yend = -0.105, linewidth = 0.3) +
  annotate("segment", x = 3.7, xend = 6.3, y = -0.105, yend = -0.105, linewidth = 0.3) +
  annotate("text", x = 1.5, y = -0.16, label = "CON vs SDS", size = 2.3) +
  annotate("text", x = 5, y = -0.16, label = "CON vs CES vs CPS", size = 2.3) +
  coord_cartesian(ylim = c(0, 1), clip = "off") +
  theme_classic(base_size = 7, base_family = "Arial") +
  theme(axis.title.y = element_text(size = 7), axis.text.x = element_text(size = 7, color = "black"),
        axis.text.y = element_text(size = 6, color = "black"), axis.ticks.x = element_blank(),
        axis.line = element_line(linewidth = 0.3), axis.ticks.y = element_line(linewidth = 0.3),
        legend.position = "right", legend.text = element_text(size = 6), legend.key.height = unit(3.2, "mm"),
        legend.key.width = unit(3.2, "mm"), plot.margin = margin(12, 4, 26, 4))
ggsave("figures/Fig1E_composition_alluvial.pdf", p, width = 100, height = 68, units = "mm", device = cairo_pdf)
ggsave("figures/Fig1E_composition_alluvial.png", p, width = 100, height = 68, units = "mm", dpi = 300)
print(as.data.frame(dat %>% select(Axis, Celltype, Frac) %>% tidyr::pivot_wider(names_from = Axis, values_from = Frac) %>% mutate(across(-Celltype, ~ round(100 * .x, 1)))))


# ====================================================================================================
# 13b Figure 1G
# ====================================================================================================
# Figure 1G (revised) - overlap of stress DEGs per major cell type with psychiatric-disease gene sets
# (logic of script/GWAS_Psychdisease_chord.R).
# Changes vs v1 (archived in figures/archive_Fig1G_v1_padj0.1_lfc0.5/):
#   * DEG threshold = the paper-wide one (adj. P < 0.05, |log2FC| > 0.25) instead of adj. P < 0.1,
#     |log2FC| > 0.5. Neuronal disease genes change by ~0.25-0.45 log2FC and were cut by the 0.5 filter.
#   * Ribbon shading: dark = enriched (hypergeometric, BH-FDR < 0.05) and still enriched after removing
#     neuron-enriched genes; light = enriched only with neuron-enriched genes; grey = not enriched.
#   * Each glial overlap gene is flagged when it is neuron-enriched (CON mean log-expression in ExNeuron
#     >= 2x that in the glial type and >= 0.5): such up-regulated genes in glia point to ambient RNA or
#     neuron-glia doublets. A sensitivity test without them is in the source data.
# Gene sets: Brain/data/GWAS_Neuron_Disease_Gene_Dataset.csv, 7 traits, human -> mouse with homologene.
# Background: genes detected in >= 1% of cells of that cell type (either group). One-sided hypergeometric, BH.
suppressMessages({library(data.table); library(homologene); library(circlize); library(jsonlite)})
src <- file.path(REF, "GWAS_Neuron_Disease_Gene_Dataset.csv")
if (!file.exists(src)) { message("Disease gene table not found: ", src); quit(status = 2) }
PADJ <- 0.05; LFC <- 0.25; FDR_SIG <- 0.05
pal <- fromJSON("palettes.json")
traits <- c("ADHD", "Bipolar.disorder", "Depression", "PTSD", "Anxiety", "Insomnia", "Neuroticism")
gw <- read.csv(src, check.names = TRUE)[, traits]
sets <- lapply(gw, function(g) { g <- unique(trimws(g)); g <- g[g != "" & !is.na(g)]
  unique(homologene(g, inTax = 9606, outTax = 10090)[["10090"]]) })
names(sets) <- gsub(".", " ", names(sets), fixed = TRUE)
cts <- c("ExNeuron", "InNeuron", "Astro", "Micro", "Oligo")
glia <- c("Astro", "Micro", "Oligo")

# neuron-enriched genes per glial type (CON animals only)
me <- fread("data/sample_mean_logexpr_major.csv.gz"); setnames(me, 1, "id")
me[, c("ct", "sample") := tstrsplit(id, "|||", fixed = TRUE)]
mcon <- me[grepl("^CON", sample)]
cmean <- function(x) colMeans(as.matrix(mcon[mcon$ct == x, !c("id", "ct", "sample")]))
ex_mean <- cmean("ExNeuron")
neuro <- lapply(setNames(glia, glia), function(g) { gm <- cmean(g); names(ex_mean)[ex_mean >= 0.5 & ex_mean >= 2 * gm] })

hyper <- function(deg, s, bg) { s <- intersect(s, bg); k <- length(intersect(deg, s))
  list(k = k, n_set = length(s), exp = length(s) * length(deg) / length(bg),
       p = phyper(k - 1, length(s), length(bg) - length(s), length(deg), lower.tail = FALSE)) }

overlap_table <- function(cond) {
  de <- fread("data/DE_celllevel_major.csv.gz")[condition == cond & reference == "CON"]
  rows <- list()
  for (ct in cts) {
    d <- de[cluster == ct]; bg <- unique(d[pct_cond >= 0.01 | pct_ref >= 0.01, gene])
    sig <- d[pvals_adj < PADJ & abs(log2FC) > LFC]; deg <- unique(sig$gene); up <- sig[log2FC > 0, gene]
    nn <- if (ct %in% glia) neuro[[ct]] else character(0)
    for (tr in names(sets)) {
      h <- hyper(deg, sets[[tr]], bg); g <- sort(intersect(deg, sets[[tr]]))
      h2 <- hyper(setdiff(deg, nn), setdiff(sets[[tr]], nn), setdiff(bg, nn))
      rows[[paste(ct, tr)]] <- data.table(contrast = paste(cond, "vs CON"), celltype = ct, trait = tr,
        overlap = h$k, expected = round(h$exp, 2), n_deg = length(deg), n_set_in_bg = h$n_set, n_bg = length(bg), p = h$p,
        n_up = sum(g %in% up), n_neuron_enriched = sum(g %in% nn),
        genes = paste(g, collapse = ";"), neuron_enriched_genes = paste(intersect(g, nn), collapse = ";"),
        overlap_without_neuron_enriched = h2$k, p_without_neuron_enriched = h2$p)
    }
  }
  tab <- rbindlist(rows)
  tab[, fdr := p.adjust(p, "BH")][, fdr_without_neuron_enriched := p.adjust(p_without_neuron_enriched, "BH")]
  setcolorder(tab, c(setdiff(names(tab), c("genes", "neuron_enriched_genes")), "genes", "neuron_enriched_genes"))
  tab
}

ct_cols <- unlist(pal$major_celltype)[cts]
tr_cols <- setNames(colorRampPalette(c("#D9D9D9", "#4D4D4D"))(length(sets)), names(sets))   # grey traits so cell-type colours stand out
draw_chord <- function(tab, title, cex = 0.55, highlight = NULL) {
  links <- tab[overlap > 0, .(from = celltype, to = trait, value = overlap,
                               cls = ifelse(fdr >= FDR_SIG, 0L, ifelse(fdr_without_neuron_enriched < FDR_SIG, 2L, 1L)))]
  lc <- ifelse(links$cls == 2L, adjustcolor(ct_cols[links$from], alpha.f = 0.85),
        ifelse(links$cls == 1L, adjustcolor(ct_cols[links$from], alpha.f = 0.28), adjustcolor("grey80", alpha.f = 0.5)))
  hl <- links$from %in% highlight
  if (any(hl)) {                                # highlighted cell type: full colour + dark outline, others muted
    lc[!hl] <- adjustcolor(lc[!hl], alpha.f = 0.55)
    lc[hl] <- adjustcolor(ct_cols[links$from[hl]], alpha.f = 0.9)
  }
  lb <- ifelse(hl, adjustcolor("#08306B", 0.9), NA)
  ord <- order(hl, links$cls)                   # highlighted, then enriched ribbons drawn last (on top)
  circos.clear()
  circos.par(start.degree = 90, canvas.xlim = c(-1.42, 1.42), canvas.ylim = c(-1.42, 1.42),
             gap.after = c(rep(2, length(cts) - 1), 8, rep(1.5, length(sets) - 1), 8), points.overflow.warning = FALSE)
  chordDiagram(as.data.frame(links[ord, .(from, to, value)]), col = lc[ord], order = c(cts, names(sets)),
               link.border = lb[ord], link.lwd = 0.6,
               grid.col = c(ct_cols, tr_cols), annotationTrack = "grid", preAllocateTracks = list(track.height = 0.08),
               link.sort = TRUE, link.decreasing = FALSE)
  ndeg <- tab[, .(n = n_deg[1]), by = celltype]
  circos.trackPlotRegion(track.index = 1, bg.border = NA, panel.fun = function(x, y) {
    s <- get.cell.meta.data("sector.index")
    lab <- if (s %in% cts) sprintf("%s (%d)", s, ndeg[celltype == s, n]) else s
    circos.text(mean(get.cell.meta.data("xlim")), get.cell.meta.data("ylim")[1] + 0.2, lab,
                facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), cex = cex) })
  title(title, cex.main = 0.7, font.main = 1, line = 0)
}
# ---------- main Fig1G: SDS vs CON ----------
sds <- overlap_table("SDS")
fwrite(sds, "figures/Fig1G_disease_chord_source_data.csv")
mainfig <- function() {
  par(mar = c(0.2, 0.2, 1.4, 0.2))
  draw_chord(sds, sprintf("SDS vs CON DEGs (adj. P < %.2g, |log2FC| > %.2g)\nvs psychiatric GWAS genes", PADJ, LFC), highlight = "ExNeuron")
}
cairo_pdf("figures/Fig1G_disease_chord.pdf", 4.0, 4.1, family = "Arial"); mainfig(); invisible(dev.off())
png("figures/Fig1G_disease_chord.png", 4.0, 4.1, units = "in", res = 300, family = "Arial"); mainfig(); invisible(dev.off())

cat("Main (SDS) FDR<0.05:\n"); print(sds[fdr < FDR_SIG, .(celltype, trait, overlap, expected, n_up, n_neuron_enriched, fdr = signif(fdr, 2), fdr_wo = signif(fdr_without_neuron_enriched, 2))])
cat("ExNeuron/InNeuron SDS:\n"); print(sds[celltype %in% c("ExNeuron", "InNeuron"), .(celltype, trait, overlap, expected, fdr = signif(fdr, 2))])
