###############################################################################
# Mouse GIST gamma-delta T-cell analysis: Imatinib vs Vehicle
###############################################################################

library(DESeq2)
library(dplyr)
library(ggplot2)
library(matrixStats)
library(fgsea)
library(org.Mm.eg.db)
library(msigdbr)
library(pheatmap)
library(EnhancedVolcano)
library(KEGGREST)
library(BiocParallel)
library(gridExtra)

setwd("C:/GIST")
dir.create("Figures", showWarnings = FALSE)
dir.create("GDT_Revision/Supplemental_Tables", recursive = TRUE, showWarnings = FALSE)

###############################################################################
# Load and prepare count data
###############################################################################

countData <- read.csv(
  "110723_DV ydT raw_counts_reduced.csv",
  header = TRUE,
  check.names = FALSE
)

summed_data <- aggregate(. ~ refGene, data = countData, FUN = sum)
rownames(summed_data) <- summed_data$refGene
summed_data$refGene <- NULL

count_matrix <- as.matrix(summed_data)
storage.mode(count_matrix) <- "integer"

write.table(
  count_matrix,
  file = "RNA_Mixture.txt",
  sep = "\t",
  quote = FALSE,
  col.names = NA
)

metaData <- read.delim(
  "GDT_Mouse_Metadata_Reduced.txt",
  header = TRUE,
  stringsAsFactors = FALSE
)

rownames(metaData) <- metaData$SampleID
metaData$Condition <- factor(metaData$Condition, levels = c("Vehicle", "Imatinib"))

stopifnot(all(colnames(count_matrix) %in% rownames(metaData)))
metaData <- metaData[colnames(count_matrix), , drop = FALSE]
stopifnot(identical(colnames(count_matrix), rownames(metaData)))

###############################################################################
# Differential expression with DESeq2
###############################################################################

dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = metaData,
  design = ~ Condition
)

dds$Condition <- relevel(dds$Condition, ref = "Vehicle")
dds <- DESeq(dds)

rld <- rlog(dds, blind = FALSE)

###############################################################################
# Plot selected genes
###############################################################################

genes_to_plot <- c(
  "Trgv1", "Trgv2", "Trgv3", "Trgv4", "Trgv5", "Trgv6",
  "Trgv7", "Trdc", "Il17a", "Ifng", "Nkg7", "Cd3e"
)

genes_to_plot <- intersect(genes_to_plot, rownames(dds))

gene_plots <- lapply(genes_to_plot, function(g) {
  gene_data <- plotCounts(
    dds,
    gene = g,
    intgroup = "Condition",
    returnData = TRUE
  )
  
  ggplot(gene_data, aes(x = Condition, y = count, color = Condition)) +
    geom_point(size = 3) +
    stat_summary(
      fun = mean,
      geom = "crossbar",
      aes(ymin = after_stat(y), ymax = after_stat(y)),
      width = 0.5
    ) +
    scale_y_log10() +
    theme_bw() +
    labs(title = g, x = NULL, y = "Normalized counts") +
    guides(color = "none") +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
})

names(gene_plots) <- genes_to_plot

if (length(gene_plots) > 0) {
  gene_plot_grid <- marrangeGrob(grobs = gene_plots, ncol = 6, nrow = 2)
  ggsave(
    "Figures/Selected_Gene_Expression.pdf",
    gene_plot_grid,
    width = 14,
    height = 7
  )
}

###############################################################################
# Heatmap of selected genes
###############################################################################

sampleInfo <- as.data.frame(colData(dds)[, "Condition", drop = FALSE])
rownames(sampleInfo) <- rownames(colData(dds))

selected_genes <- c(
  "Trdc", "Trgv1", "Trgv2", "Trgv4", "Trgv5", "Trgv6", "Trgv7",
  "Trdv4", "Cxcr6", "Ccr2", "Il17a", "Cd44", "Ccr6", "Il23r",
  "Rorc", "Nr1d1", "Maf", "Tnfrsf25", "Ccl5", "Klrd1", "Klrc2",
  "Xaf1", "Runx3", "Gzmb", "Pdcd1", "Havcr2", "Lag3", "Tnfrsf9",
  "Pgk1", "Sdc2", "Igsf9"
)

selected_genes <- intersect(selected_genes, rownames(assay(rld)))

if (length(selected_genes) > 0) {
  pheatmap(
    assay(rld)[selected_genes, , drop = FALSE],
    scale = "row",
    show_rownames = TRUE,
    show_colnames = FALSE,
    annotation_col = sampleInfo,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    color = colorRampPalette(c("navy", "white", "red"))(50),
    filename = "Figures/Genes_Heatmap.pdf",
    width = 4.5,
    height = 6
  )
}

###############################################################################
# PCA and normalized counts
###############################################################################

pca_plot <- plotPCA(
  rld,
  intgroup = "Condition",
  ntop = nrow(rld)
) + theme_bw()

ggsave(
  "Figures/PCA.pdf",
  pca_plot,
  width = 6,
  height = 5
)

normalized_counts <- counts(dds, normalized = TRUE)

write.table(
  normalized_counts,
  file = "Normalized_Counts_Reduced.txt",
  sep = "\t",
  quote = FALSE,
  col.names = NA
)

###############################################################################
# Differential expression results
###############################################################################

resultsNames(dds)

res_unshrunk <- results(
  dds,
  name = "Condition_Imatinib_vs_Vehicle"
)

res_shrunk <- lfcShrink(
  dds,
  coef = "Condition_Imatinib_vs_Vehicle",
  type = "apeglm"
)

write.table(
  as.data.frame(res_shrunk),
  file = "GDT_Revision/Supplemental_Tables/Condition_Imatinib_vs_Vehicle_DESeq2_Second_Sort.txt",
  sep = "\t",
  quote = FALSE,
  col.names = NA
)

###############################################################################
# Prepare gene rankings for GSEA
###############################################################################

res_rank <- as.data.frame(res_unshrunk) %>%
  tibble::rownames_to_column("gene") %>%
  filter(!is.na(stat)) %>%
  group_by(gene) %>%
  summarize(stat = mean(stat), .groups = "drop")

ranks <- res_rank$stat
names(ranks) <- res_rank$gene
ranks <- sort(ranks, decreasing = TRUE)

###############################################################################
# Prepare Hallmark, Reactome, and KEGG pathway collections
###############################################################################

hallmark_msig <- msigdbr(
  species = "Mus musculus",
  category = "H"
)

reactome_msig <- msigdbr(
  species = "Mus musculus",
  category = "C2",
  subcategory = "CP:REACTOME"
)

pathways_hallmark <- split(hallmark_msig$gene_symbol, hallmark_msig$gs_name)
pathways_reactome <- split(reactome_msig$gene_symbol, reactome_msig$gs_name)

kegg_raw <- keggLink("pathway", "mmu")
pathway_ids <- gsub("\\D", "", kegg_raw)
entrez_ids <- sub("mmu:", "", names(kegg_raw))
kegg_to_entrez <- split(entrez_ids, pathway_ids)

kegg_titles <- keggList("pathway", "mmu")
title_ids <- gsub("\\D", "", names(kegg_titles))
clean_titles <- sub(" - Mus musculus \\(mouse\\)", "", kegg_titles)
names(clean_titles) <- title_ids

pathways_kegg <- lapply(kegg_to_entrez, function(ids) {
  symbols <- mapIds(
    org.Mm.eg.db,
    keys = ids,
    column = "SYMBOL",
    keytype = "ENTREZID",
    multiVals = "first"
  )
  unname(na.omit(symbols))
})

matched_titles <- clean_titles[names(pathways_kegg)]
formatted_names <- gsub("[^A-Za-z0-9]+", "_", toupper(matched_titles))
final_names <- ifelse(
  is.na(formatted_names),
  paste0("PATHWAY_", names(pathways_kegg)),
  formatted_names
)

names(pathways_kegg) <- paste0("KEGG_", final_names)

combined_pathways <- c(
  pathways_hallmark,
  pathways_reactome,
  pathways_kegg
)

###############################################################################
# GSEA
###############################################################################

fgsea_res <- fgsea(
  pathways = combined_pathways,
  stats = ranks,
  minSize = 15,
  maxSize = 500,
  BPPARAM = SerialParam()
) %>%
  as.data.frame() %>%
  arrange(padj)

fgsea_export <- fgsea_res %>%
  mutate(leadingEdge = vapply(leadingEdge, paste, collapse = ";", character(1)))

write.table(
  fgsea_export,
  file = "GDT_Revision/Supplemental_Tables/Mouse_GSEA_Imatinib_vs_Control_UpdatedPathways.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

###############################################################################
# Plot top significantly enriched pathways
###############################################################################

sig_pathways <- fgsea_res %>%
  filter(!is.na(padj), padj < 0.05)

top_up <- sig_pathways %>%
  filter(NES > 0) %>%
  arrange(desc(NES)) %>%
  head(50)

top_down <- sig_pathways %>%
  filter(NES < 0) %>%
  arrange(NES) %>%
  head(50)

plot_df <- bind_rows(top_up, top_down) %>%
  distinct(pathway, .keep_all = TRUE) %>%
  arrange(NES) %>%
  mutate(pathway = factor(pathway, levels = unique(pathway)))

if (nrow(plot_df) > 0) {
  p_gsea_top <- ggplot(plot_df, aes(x = reorder(pathway, NES), y = NES, fill = -log10(padj))) +
    geom_bar(stat = "identity", position = "identity", width = 0.7) +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 10) +
    coord_flip() +
    labs(x = NULL, y = "Normalized Enrichment Score", fill = "-Log10(FDR)") +
    theme_bw() +
    theme(
      axis.text.x = element_text(size = 8, color = "black"),
      axis.text.y = element_text(size = 12, color = "black"),
      axis.line = element_blank(),
      plot.margin = margin(5, 5, 5, 5, unit = "pt")
    )
  
  ggsave(
    "Figures/GSEA_Top_Pathways.pdf",
    p_gsea_top,
    width = 11,
    height = 8
  )
}

###############################################################################
# Plot selected pathways
###############################################################################

selected_pathways <- c(
  "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
  "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
  "KEGG_PYRUVATE_METABOLISM_MUS_MUSCULUS_HOUSE_MOUSE_",
  "HALLMARK_FATTY_ACID_METABOLISM",
  "KEGG_GLUTATHIONE_METABOLISM_MUS_MUSCULUS_HOUSE_MOUSE_",
  "HALLMARK_XENOBIOTIC_METABOLISM",
  "REACTOME_INTERLEUKIN_12_SIGNALING",
  "HALLMARK_GLYCOLYSIS",
  "KEGG_DRUG_METABOLISM_OTHER_ENZYMES_MUS_MUSCULUS_HOUSE_MOUSE_",
  "HALLMARK_MTORC1_SIGNALING",
  "HALLMARK_HYPOXIA",
  "KEGG_TH17_CELL_DIFFERENTIATION_MUS_MUSCULUS_HOUSE_MOUSE_",
  "HALLMARK_INFLAMMATORY_RESPONSE",
  "KEGG_NATURAL_KILLER_CELL_MEDIATED_CYTOTOXICITY_MUS_MUSCULUS_HOUSE_MOUSE_",
  "KEGG_FC_GAMMA_R_MEDIATED_PHAGOSOME_FORMATION_MUS_MUSCULUS_HOUSE_MOUSE_",
  "KEGG_ALLOGRAFT_REJECTION_MUS_MUSCULUS_HOUSE_MOUSE_",
  "REACTOME_CDC42_GTPASE_CYCLE",
  "KEGG_ANTIGEN_PROCESSING_AND_PRESENTATION_MUS_MUSCULUS_HOUSE_MOUSE_",
  "REACTOME_RAC1_GTPASE_CYCLE",
  "KEGG_GRAFT_VERSUS_HOST_DISEASE_MUS_MUSCULUS_HOUSE_MOUSE_",
  "KEGG_JAK_STAT_SIGNALING_PATHWAY_MUS_MUSCULUS_HOUSE_MOUSE_",
  "HALLMARK_MITOTIC_SPINDLE",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "KEGG_CYTOSOLIC_DNA_SENSING_PATHWAY_MUS_MUSCULUS_HOUSE_MOUSE_",
  "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "HALLMARK_INTERFERON_ALPHA_RESPONSE"
)

selected_pathways_clean <- sub(
  "_MUS_MUSCULUS_HOUSE_MOUSE_$",
  "",
  selected_pathways
)

sig_pathways_selected <- sig_pathways %>%
  filter(pathway %in% selected_pathways_clean) %>%
  mutate(pathway = factor(pathway, levels = pathway[order(NES)]))

if (nrow(sig_pathways_selected) > 0) {
  p_gsea_selected <- ggplot(
    sig_pathways_selected,
    aes(x = pathway, y = NES, fill = -log10(padj))
  ) +
    geom_bar(stat = "identity", position = "identity", width = 0.7) +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 12) +
    coord_flip() +
    labs(x = NULL, y = "Normalized Enrichment Score", fill = "-Log10(FDR)") +
    theme_bw() +
    theme(
      axis.text.x = element_text(size = 8, color = "black"),
      axis.text.y = element_text(size = 12, color = "black"),
      axis.line = element_blank(),
      plot.margin = margin(5, 5, 5, 5, unit = "pt")
    )
  
  ggsave(
    "Figures/GSEA_Barplot.pdf",
    p_gsea_selected,
    width = 11,
    height = 3.5
  )
}

###############################################################################
# Volcano plot
###############################################################################

res1 <- as.data.frame(res_unshrunk) %>%
  tibble::rownames_to_column("genename") %>%
  filter(!is.na(padj), !is.na(log2FoldChange))

genes_to_label <- c(
  "Thbs2", "A2m", "Adgrd1", "Hmcn2", "Spock1", "Prkg1", "Kcnk3",
  "Igfbp5", "S100a8", "Gja1", "Cgnl1", "Fn1", "Ror2", "Cxcr2",
  "Enpep", "S100a9", "Atp1a2", "Tcf21", "Islr", "Apol7c", "Ccdc80",
  "Trpm2", "Gm20559", "Il1a", "Ccl12", "Ccl7", "Ifi213", "Igsf9",
  "Oas2", "Lilrb4a", "Gm43197", "Mx1", "Epha2"
)

p_volcano <- EnhancedVolcano(
  res1,
  pCutoff = 0.05,
  FCcutoff = 1,
  lab = res1$genename,
  selectLab = intersect(genes_to_label, res1$genename),
  x = "log2FoldChange",
  y = "padj",
  raster = TRUE,
  drawConnectors = TRUE,
  xlim = c(-8, 8),
  ylim = c(0, 40),
  labSize = 4,
  max.overlaps = 20,
  title = "Imatinib versus vehicle"
)

ggsave(
  "Figures/Volcano_Plot.pdf",
  p_volcano,
  width = 7,
  height = 7
)

###############################################################################
# Export significantly downregulated genes
###############################################################################

filtered_rows <- res1 %>%
  filter(log2FoldChange < -1, padj < 0.05)

write.table(
  filtered_rows,
  file = "GDT_Revision/Supplemental_Tables/Downregulated_Genes_Imatinib_vs_Vehicle.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

nrow(filtered_rows)
