library(Seurat)
library(ggplot2)
library(matrixStats)
library(DESeq2)
library(tidyverse)
library(S4Vectors)
library(fgsea)
library(org.Hs.eg.db)
library(msigdbr)
library(data.table)
library(pheatmap)
library(stringr)
library(Seurat)
library(gridExtra)
library(EnhancedVolcano)
library(ggpubr)
library(scCustomize)

setwd("C:/GIST")
gdt = readRDS("Murine gdt.rds")
Idents(gdt) = "seurat_clusters"
markers = FindAllMarkers(gdt)
write.table(markers, file = "Mouse_scRNA_Markers.tsv", quote = F, sep = '\t')

gdt$celltype = Idents(gdt)
table(gdt$celltype)
p1 <- DimPlot(gdt, group.by = "celltype", raster = F, pt.size = 0.15, split.by = "treatment") + coord_fixed()
ggsave(p1, file = "Figures/Mouse_GDT_UMAP.pdf", width = 3, height = 4.5)

DefaultAssay(gdt) = "RNA"
p2 <- FeaturePlot(gdt, features = "Tcrg-V6", raster = F, pt.size = 0.15) + coord_fixed()
p2
ggsave(p2, file = "Figures6/Mouse_GDT_UMAP_FeaturePlot.pdf", width = 3, height = 4.5)

p2 <- FeaturePlot(gdt, features = "Tcrg-V7", raster = F, pt.size = 0.15) + coord_fixed()
p2
ggsave(p2, file = "Figures6/Mouse_GDT_UMAP_FeaturePlot2.pdf", width = 3, height = 4.5)

p2 <- FeaturePlot(gdt, features = "Trgv1", raster = F, pt.size = 0.15) + coord_fixed()
p2
ggsave(p2, file = "Figures6/Mouse_GDT_UMAP_FeaturePlot2.pdf", width = 3, height = 4.5)


Stacked_VlnPlot(gdt, features = c("Il17a", "Ifng", "Klrd1", "Gzmb"), group.by = "seurat_clusters", pt.size = 0.5)

levels(Idents(gdt))
levels(as.numeric(Idents(gdt)))

counts = as.matrix(GetAssayData(gdt, assay = "RNA", slot = "counts"))
length(colnames(counts))
length(gdt$seurat_clusters)

counts_cs = rbind(Idents(gdt), counts)
write.table(counts_cs, sep = '\t', quote = F, "Murine_CIBERSORTx_Ref2.txt")

##############
countData <- read.csv("110723_DV ydT raw_counts_reduced.csv", sep=",", header = TRUE)
countData = as.data.frame(countData)
summed_data <- aggregate(. ~ refGene, data = countData, FUN = sum)

rownames(summed_data) = summed_data$refGene
summed_data$refGene = NULL
write.table(summed_data, sep = '\t', quote = F, "RNA_Mixture.txt")
metaData <- read.csv("GDT_Mouse_Metadata_Reduced.txt", sep="\t", header = TRUE)
metaData

dds <- DESeqDataSetFromMatrix(countData=summed_data, 
                              colData=metaData, 
                              design = ~ Condition, tidy = FALSE)
dds
dds$Condition <- relevel(dds$Condition, "Vehicle")
dds <- DESeq(dds)
resultsNames(dds)
rld <- rlog(dds, blind=FALSE)
head(rld)

gene = c("Trgv1", "Trgv2", "Trgv3", "Trgv4", "Trgv5", "Trgv6", "Trgv7", "Trdc", 
         "Il17a", "Ifng", "Nkg7", "Cd3e")
length(gene)
dds$Condition <- relevel(dds$Condition, "Vehicle")
gene_plots <- list()
for(g in gene){
  gene_data <- plotCounts(dds, gene = g, intgroup = "Condition", returnData = TRUE)
  gene_plots[[g]] <- ggplot(gene_data, aes(factor(Condition), count)) +
    geom_point(aes(color = Condition), size = 3) +  
    stat_summary(fun.y="mean",geom="crossbar", 
                 mapping=aes(ymin=..y.., ymax=..y..), width=0.5) +
    scale_y_log10() + theme_bw() + ggtitle(g) + RotatedAxis() + xlab("") + 
    ylab("Normalized Counts") + guides(fill = FALSE, colour = F, linetype = FALSE, shape = FALSE) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold")) 
}
grid.arrange(grobs = gene_plots, ncol = 6, nrow = 2)

dds$Condition <- relevel(dds$Condition, "Imatinib")

sampleInfo <- as.data.frame(colData(dds)[,c("Condition")])
rownames(sampleInfo) <- colData(dds)[,c("SampleID")]
colnames(sampleInfo) <- c("Condition")
sampleInfo[] <- lapply(sampleInfo, function(x) str_replace_all(x, "-", "."))
rownames(sampleInfo) <- gsub("-", ".", rownames(sampleInfo))

topVarGenes <- head(order(rowVars(assay(rld)), decreasing = TRUE), 50)

selected_genes = c("Trdc", "Trgv1", "Trgv2", "Trgv4", "Trgv5", "Trgv6", 
                   "Trgv7", "Trdv4", "Cxcr6", "Ccr2", "Il17a", "Cd44", 
                   "Ccr6", "Il23r", "Rorc", "Nr1d1", "Maf", 
                   "Tnfrsf25", "Ccl5", "Klrd1", "Klrc2", "Xaf1", "Runx3", "Gzmb", 
                   "Pdcd1", "Havcr2", "Lag3", "Tnfrsf9", "Pgk1", "Sdc2", "Igsf9")
indices_selected_genes <- match(selected_genes, rownames(assay(rld)))

p2 = pheatmap(assay(rld)[indices_selected_genes, ], scale = 'row', show_rownames = T, 
         annotation_col = sampleInfo, show_colnames = F, cluster_rows = F,
         cluster_cols = F, 
         color=colorRampPalette(c("navy", "white", "red"))(50),)
ggsave(p2, file = "Figures/Genes_Heatmap.pdf", width = 4.5, height = 6)


#PCA with top n variable genes 
DESeq2::plotPCA(rld, intgroup = "Condition", ntop = 300000) + theme_bw()
normalized_counts <- counts(dds, normalized=TRUE)
write.table(normalized_counts, file="Normalized_Counts_Reduced.txt", sep="\t", quote=F, col.names=NA)


#Gene Set Enrichment Analysis 
#Line below is key code to generate results, ensure tidy = TRUE
resultsNames(dds)
res <- results(dds, name="Condition_Imatinib_vs_Vehicle")
res <- lfcShrink(dds, res = res, coef = "Condition_Imatinib_vs_Vehicle", type = "apeglm")
write.table(res, file = "Condition_Imatinib_vs_Vehicle_DESeq2_Second_Sort.txt", sep = '\t', quote = F)

res <- results(dds, name="Condition_Imatinib_vs_Vehicle", tidy = TRUE)
res2 <- res %>% 
  dplyr::select(row, stat) %>% 
  na.omit() %>%
  distinct() %>%
  group_by(row) %>%
  summarize(stat=mean(stat))
res2
ranks <- deframe(res2)

pathways.kegg <- msigdbr(species = "mouse", category = "C2", subcategory = "KEGG")
pathways.reactome <- msigdbr(species = "mouse", category = "C2", subcategory = "REACTOME")
pathways.hallmark <- msigdbr(species = "mouse", category = "H")

pathway_list_kegg = split(x = pathways.kegg$gene_symbol, f = pathways.kegg$gs_name)
pathway_list_hallmark = split(x = pathways.hallmark$gene_symbol, f = pathways.hallmark$gs_name)
pathway_list <- c(pathway_list_kegg, pathway_list_hallmark)

GSEA_result <-  fgsea(pathways = pathway_list, stats = ranks)

fgseaResTidy <- GSEA_result %>%
  as_tibble() %>%
  arrange(desc(NES))

head(GSEA_result[order(pval), ])

fgseaResTidy %>% 
  dplyr::select(-ES, -leadingEdge) %>% 
  arrange(padj) %>% 
  DT::datatable()
write.table(as.matrix(GSEA_result), quote = F, sep = '\t', file = "GSEA_Kegg_Hallmark.tsv")

sig_pathways <- GSEA_result[GSEA_result$padj<0.05, ]
sig_pathways <- sig_pathways[order(sig_pathways$NES), ]
sig_pathways = na.omit(sig_pathways)

top_50_pathways <- head(sig_pathways, 20)
bottom_50_pathways <- tail(sig_pathways, 20)
sig_pathways_filter <- rbind(top_50_pathways, bottom_50_pathways)
sig_pathways_filter <- sig_pathways_filter[!duplicated(sig_pathways_filter[["pathway"]]), ]
sig_pathways_filter$pathway <- factor(
  sig_pathways_filter$pathway,
  levels = sig_pathways_filter$pathway[order(sig_pathways_filter$NES)]
)

selected_pathways = c("HALLMARK_OXIDATIVE_PHOSPHORYLATION", "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
                      "HALLMARK_FATTY_ACID_METABOLISM", "KEGG_PYRUVATE_METABOLISM", 
                      "KEGG_GLUTATHIONE_METABOLISM", "KEGG_GLYCOLYSIS_GLUCONEOGENESIS", 
                      "HALLMARK_MYC_TARGETS_V1", "KEGG_DRUG_METABOLISM_CYTOCHROME_P450", 
                      "KEGG_BUTANOATE_METABOLISM", "KEGG_ARGININE_AND_PROLINE_METABOLISM", 
                      "HALLMARK_INFLAMMATORY_RESPONSE", "KEGG_NATURAL_KILLER_CELL_MEDIATED_CYTOTOXICITY", 
                      "KEGG_ALLOGRAFT_REJECTION", "KEGG_ANTIGEN_PROCESSING_AND_PRESENTATION", 
                      "KEGG_GRAFT_VERSUS_HOST_DISEASE", "HALLMARK_TNFA_SIGNALING_VIA_NFKB", 
                      "KEGG_CYTOSOLIC_DNA_SENSING_PATHWAY", "KEGG_JAK_STAT_SIGNALING_PATHWAY", 
                      "KEGG_NOD_LIKE_RECEPTORd_SIGNALING_PATHWAY", "HALLMARK_INTERFERON_GAMMA_RESPONSE",
                      "HALLMARK_INTERFERON_ALPHA_RESPONSE")
sig_pathways_selected = sig_pathways[sig_pathways$pathway %in% selected_pathways, ]
sig_pathways_selected$pathway <- factor(
  sig_pathways_selected$pathway,
  levels = sig_pathways_selected$pathway[order(sig_pathways_selected$NES)]
)

p3 <- ggplot(sig_pathways_selected, aes(x = pathway, y = NES, fill = -log10(padj))) +
  geom_bar(stat = "identity", position = "identity") +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint =12) +
  labs(x = "", y = "Normalized Enrichment Score", fill = "-log10(FDR)") +  # Update axis and legend labels
  coord_flip()+ ggtitle("") + 
  theme(
    axis.text.x = element_text(angle = 0, vjust = 1, hjust = 1, family = "sans", size = 8, color = "black"),
    axis.text.y = element_text(family = "sans", size = 12, color = "black"),
    axis.line = element_blank(), 
    plot.margin = margin(l = 5, r = 5, b = 5, t = 5, unit = "pt")
  ) + theme_bw()
p3
ggsave(p3, file = "Figures/GSEA_Barplot.pdf", width = 7, height = 4.5)

#Leading edge plots 
p3a <- plotEnrichment(pathway_list[["HALLMARK_INTERFERON_GAMMA_RESPONSE"]],
               ranks) + labs(title="Interferon Gamma Response") + theme_bw()
p3b <- plotEnrichment(pathway_list[["HALLMARK_TNFA_SIGNALING_VIA_NFKB"]],
                      ranks) + labs(title="TNA Signaling via NFkB") + theme_bw()
p3c <- plotEnrichment(pathway_list[["KEGG_ANTIGEN_PROCESSING_AND_PRESENTATION"]],
                      ranks) + labs(title="Antigen Presentation/Processing") + theme_bw()
p3d <- plotEnrichment(pathway_list[["KEGG_NATURAL_KILLER_CELL_MEDIATED_CYTOTOXICITY"]],
                      ranks) + labs(title="NK Cell Mediated Cytotoxicity") + theme_bw()
p3e <- plotEnrichment(pathway_list[["HALLMARK_OXIDATIVE_PHOSPHORYLATION"]],
                      ranks) + labs(title="Oxidative Phosphorylation") + theme_bw()
p3f <- plotEnrichment(pathway_list[["HALLMARK_FATTY_ACID_METABOLISM"]],
                      ranks) + labs(title="Fatty Acid Metabolism") + theme_bw()
p3g <- plotEnrichment(pathway_list[["KEGG_GLYCOLYSIS_GLUCONEOGENESIS"]],
                      ranks) + labs(title="Glycolysis/Gluconeogenesis") + theme_bw()
p6 <- grid_plot <- grid.arrange(p3a, p3b, p3c, p3d, p3e, p3f, p3g, nrow = 2)
p6
ggsave(p6, file = "Figures/GSEA_Individual_Plots.pdf", width = 11.33, height = 4)

#Volcano 
res1 = as.data.frame(na.omit(res))
res1$genename = rownames(res1) 
p4 <- EnhancedVolcano(res1, pCutoff = 0.05, FCcutoff = 1, 
                lab = res1$genename,
                selectLab = c("Thbs2", "A2m", "Adgrd1", "Hmcn2", "Spock1", 
                              "Prkg1", "Kcnk3", "Igfbp5", "S100a8", 
                              "Gja1", "Cgnl1", "Fn1", "Ror2", 
                              "Cxcr2", "Enpep", "S100a9", "Atp1a2", "Tcf21", 
                              "Islr", "Apol7c", "Ccdc80", "Trpm2", "Gm20559", 
                              "Il1a", "Ccl12", "Ccl7", "Ifi213", "Igsf9", 
                              "Oas2", "Lilrb4a", "Gm43197", "Mx1", "Epha2"),
                x = 'log2FoldChange',
                y = 'padj', raster = T, 
                drawConnectors = TRUE,
                xlim = c(-8,8),
                ylim = c(0, 40),
                labSize = 4.0, max.overlaps = 20,
                title = "Imatinib versus vehicle")
ggsave(p4, file = "Figures/Volcano_Plot.pdf", width = 7, height = 7)
filtered_rows <- subset(res1, log2FoldChange  < -1 & padj < 0.05)
nrow(filtered_rows)

#Deconvolution results 
results = read.csv("CIBERSORTx_SMode_Reduced.csv", header = TRUE)
rownames(results) <- results$Mixture
results$Mixture <- NULL
metaData
results$phenotype <- metaData$Condition 
results$phenotype <- factor(results$phenotype)
results <- subset(results, !grepl("TISSUE", phenotype))
dim(results)

p0 <- ggplot(results, aes(x = phenotype, y = X1, fill = phenotype)) +
  geom_boxplot() + geom_jitter() + scale_fill_brewer(palette = "Set3") +
  labs(
    x = "Phenotype",
    y = "Proportion of Cluster 1"
  ) + 
  theme_bw() + RotatedAxis() + 
  theme(
    legend.position = "none"
  ) +  stat_compare_means(method = "t.test" 
  )
p1 <- ggplot(results, aes(x = phenotype, y = X2, fill = phenotype)) +
  geom_boxplot() + geom_jitter() +scale_fill_brewer(palette = "Set3") +
  labs(
    x = "Phenotype",
    y = "Proportion of Cluster 2"
  ) +
  theme_bw() + RotatedAxis() + 
  theme(
    legend.position = "none"
  )+  stat_compare_means(method = "t.test" 
  )
p2 <- ggplot(results, aes(x = phenotype, y = X3, fill = phenotype)) +
  geom_boxplot() + geom_jitter() + scale_fill_brewer(palette = "Set3") +
  labs(
    x = "Phenotype",
    y = "Proportion of Cluster 3"
  ) +
  theme_bw() + RotatedAxis() + 
  theme(
    legend.position = "none"
  )+  stat_compare_means(method = "t.test" 
  )
p3 <- ggplot(results, aes(x = phenotype, y = X4, fill = phenotype)) +
  geom_boxplot() + geom_jitter() + scale_fill_brewer(palette = "Set3") +
  labs(
    x = "Phenotype",
    y = "Proportion of Cluster 4"
  ) +
  theme_bw() + RotatedAxis() + 
  theme(
    legend.position = "none"
  )+  stat_compare_means(method = "t.test",  
  )
p5 <- grid_plot <- grid.arrange(p0, p1, p2, p3, ncol = 2)
p5
ggsave(p5, file = "Figures/CIBERSORT_Boxplots.pdf", width = 4, height = 6)

