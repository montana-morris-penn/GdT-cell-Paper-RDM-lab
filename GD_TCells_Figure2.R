library(Seurat)
library(survival)
library(dplyr)
library(limma)
library(ggpubr)
library(SeuratDisk)
library(dittoSeq)
library(S4Vectors)
library(fgsea)
library(org.Hs.eg.db)
library(msigdbr)
library(gridExtra)
library(ggrepel)
library(tidyverse)
library(cowplot)
library(RColorBrewer)
library(clusterProfiler)
library(CytoTRACE)
library(ggpubr)
library(AUCell)
library(pheatmap)

setwd("/mnt/isilon/tan_lab/sussmanj/GIST_Multiome/GD_TCells")

gdt_seurat <- readRDS("/mnt/isilon/tan_lab/sussmanj/GIST_Multiome/GD_TCells/GDT_Seurat_Annotated_05262024.RDS")

cti <- GetAssayData(gdt_seurat, slot = 'data', assay = 'RNA')
cti <- as.data.frame(cti) # cytotrace error does not mention but need to convert matrix to DF to run
cytotrace <- CytoTRACE(cti, ncores = 8) 
gdt_seurat <- AddMetaData(gdt_seurat, cytotrace$CytoTRACE, col.name = "CytoTRACE_score")
#saveRDS(gdt_seurat, "/mnt/isilon/tan_lab/sussmanj/GIST_Multiome/GD_TCells/GDT_Seurat_Annotated_CytoTRACE_05262024.RDS")

gdt_seurat = readRDS("/mnt/isilon/tan_lab/sussmanj/GIST_Multiome/GD_TCells/GDT_Seurat_Annotated_CytoTRACE_05262024.RDS")
DimPlot(gdt_seurat, raster = T, raster.dpi = c(3000,3000), pt.size = 4.5, group.by = "samples", cols = colors4) + coord_fixed()
ggsave(p3, filename = "Figures/UMAP_Sample.pdf", width = 7, height = 7)

p1 = FeaturePlot(gdt_seurat, reduction = "umap.rpca", features = "CytoTRACE_score", 
                 raster = T, raster.dpi = c(3000,3000), pt.size = 4.5,)+ coord_fixed()+
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
ggsave(p1, filename = "Figures/UMAP_CytoTRACE.pdf", width = 7, height = 7)

VlnPlot(gdt_seurat, features = "CytoTRACE_score", group.by = "cell_type", raster = F, pt.size = 0) +
  geom_boxplot(width=0.1, fill="white", outlier.shape = NA)

VlnPlot(gdt_seurat, features = "CytoTRACE_score", group.by = "broad_class", split.by = "condition", raster = F, pt.size = 0) +
  geom_boxplot(width=0.1, fill="white", outlier.shape = NA)


gdt_seurat$broad_class = factor(gdt_seurat$broad_class, levels = c("VD1+", "VD2+", "VD3+", "Mixed/Unknown"))
p2 = VlnPlot(gdt_seurat, features = "CytoTRACE_score", group.by = "broad_class", split.by = "condition", raster = F, pt.size = 0) 
p2 
ggsave(p2, filename = "Figures/Violin_CytoTRACE.pdf", width = 4, height = 3)

#SCENIC figures 
gdt_rna <- readRDS("GDT_Seurat_Annotated_CytoTRACE_SCENIC_05292024.RDS")

Idents(gdt_rna) = "cell_type"
combined_averages <- as.data.frame(AverageExpression(gdt_rna, return.seurat = FALSE, slot = "data", assays = "SCENIC"))
regulons <- c("POLR3G(+)","EOMES(+)","ASCL2(+)","CEBPA(+)","TBX21(+)","ZBTB7A(+)",
              "RUNX2(+)","TCF7(+)","RARG(+)","NF1(+)",
              "NFAT5(+)","PRDM1(+)", "IRF6(+)","IRF8(+)", "STAT1(+)","KLF2(+)", "ZNF700(+)","PPARA(+)","TCF7L1(+)","GLIS1(+)",
              "E2F7(+)","PATZ1(+)","TFAP4(+)",
              "RORC(+)","CEBPD(+)","BCL3(+)","TCF3(+)","RORA(+)","HOXA9(+)",
              "NR3C2(+)","GRHL1(+)","FOSL2(+)","TCF7L2(+)","NR1H2(+)","ZNF768(+)")
regulons = unique(regulons)

p1 = pheatmap(combined_averages[regulons, ], scale = "row", cluster_rows = F, cluster_col = F, angle_col = 45, 
              color = colorRampPalette(c("blue", "white", "red"))(100)) 
ggsave(p1, file = "Figures/SCENIC_Heamtap.pdf", width = 5, height = 8.25)

#Volcano plot 
##Differential Gene analysis
DefaultAssay(gdt_rna) = "SCENIC"
Idents(gdt_rna) = "condition"
detfs = FindAllMarkers(gdt_rna)
write.table(detfs, "SCENIC_Regulons_Condition_Differential.txt", sep = '\t', quote = F)

Idents(gdt_rna) = "broad_class"
detfs = FindAllMarkers(gdt_rna)
write.table(detfs, "SCENIC_Regulons_Broad_Class_Differential.txt", sep = '\t', quote = F)

groupby = "condition"
compartment = "All"
gdt_subset = gdt_rna
comparisons = c("T", "P")

groupby = "broad_class"
compartment = "T"
gdt_subset = subset(gdt_rna, subset=condition==compartment)
comparisons = c("VD1+", "VD2+")

group1 = comparisons[1]
group2= comparisons[2]

print(compartment)
print(paste0(group1, "_vs_", group2))
de.genes = FindMarkers(gdt_subset, 
                       verbose = T,
                       assay = "SCENIC", 
                       logfc.threshold = 0, 
                       min.pct = 0.001, 
                       ident.1 = group1, 
                       ident.2 = group2, 
                       group.by = groupby,
                       max.cells.per.ident = 1500000)
de.genes$group = "de.genes"
de.genes$gene = rownames(de.genes)
de.genes$enriched = ifelse(de.genes$avg_log2FC < 0, yes = group2, no = group1)

tfs <- c(
  "CEBPD(+)", "CEBPB(+)", "BCL3(+)", "RORC(+)", "RUNX3((+))(+)", "RORA(+)", "TCF7(+)", 
  "TCF3(+)", "ASCL2(+)", "GATA6(+)", "ZBTB7A(+)", "KLF2(+)", "MYC(+)",
  "EOMES(+)", "TBX21(+)", "POLR3R(+)", "RUNX3(+)", "TCF7L1(+)", "TCF4(+)",
  "HIVEP1(+)", "BATF3(+)", "FOXO3(+)", "MAFB(+)", "NFIC(+)", 
  "E2F7(+)", "NFYB(+)", "IRF6(+)", "NR3C2(+)", "NFAT5(+)",
  "ZNF148(+)", "NF1(+)", "HMBOX1(+)", "GRHL1(+)", "GLIS1(+)", "ZFP64(+)"
)
label_points <- de.genes %>% filter(gene %in% unique(tfs))  # Replace "gene1", "gene2", "gene3" with your gene names

colors3 = readRDS("Figures/Condition_Colors.rds")

gene.stringent = ggplot(de.genes %>% filter(abs(avg_log2FC) > 0.01, p_val_adj < 1e-4), 
                        aes(x = avg_log2FC, y = -log(p_val_adj+1e-300, base = 10),
                            color = enriched,
                            label = gene)) + 
  geom_point(data = de.genes, color = "gray")+
  geom_point()+
  #geom_text_repel(data = label_points, aes(label = gene), size = 3, max.overlaps = 30) + 
  geom_text_repel(size = 3, max.overlaps = 20) + 
  theme_bw() + 
  labs(x = "Log2FC", 
       y = "-log(pv_adjust)") + #xlim(c(-3, 3))+ ylim(c(2, 300)) +
  scale_color_manual(values =  c("#BD20BF","#9DA52F"))+
  #scale_color_manual(values =  c("#F6222E","#66B0FF"))+
  theme(legend.position = "none")
gene.stringent
ggsave(gene.stringent, filename = paste0("Figures/",group1,"_vs_",group2,"_",compartment,"_SCENIC_Volcano_Labels.pdf"), width = 6, height = 6)

p1 <- FeaturePlot(gdt_rna, features = c("CEBPD(+)"), raster = TRUE, cols = c('grey', 'red'),
                  raster.dpi = c(800,800), pt.size = 2.5, max.cutoff = 'q99', min.cutoff = 'q1'
)  + coord_fixed() + ggtitle("CEBPD Regulon") +
  theme(axis.title.y=element_blank(), axis.text.y=element_blank(),axis.ticks.y=element_blank()) + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank()) 
p2 <- FeaturePlot(gdt_rna, features = c("PRDM1(+)"), raster = TRUE, cols = c('grey', 'red'),
                  raster.dpi = c(800,800), pt.size = 2.5, max.cutoff = 'q99', min.cutoff = 'q1'
) + coord_fixed()+ ggtitle("PRDM1 Regulon") +
  theme(axis.title.y=element_blank(), axis.text.y=element_blank(),axis.ticks.y=element_blank()) + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank()) 
p3 <- FeaturePlot(gdt_rna, features = c("RORC(+)"), raster = TRUE, cols = c('grey', 'red'),
                  raster.dpi = c(800,800), pt.size = 2.5, max.cutoff = 'q99', min.cutoff = 'q1'
)  + coord_fixed()+ggtitle("RORC Regulon") +
  theme(axis.title.y=element_blank(), axis.text.y=element_blank(),axis.ticks.y=element_blank()) + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank()) 
p4 <- FeaturePlot(gdt_rna, features = c("EOMES(+)"), raster = TRUE, cols = c('grey', 'red'),
                  raster.dpi = c(800,800), pt.size = 2.5, max.cutoff = 'q99', min.cutoff = 'q1'
)  + coord_fixed()+ggtitle("EOMES Regulon") +
  theme(axis.title.y=element_blank(), axis.text.y=element_blank(),axis.ticks.y=element_blank()) + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank()) 
p5 <- FeaturePlot(gdt_rna, features = c("TBX21(+)"), raster = TRUE, cols = c('grey', 'red'),
                  raster.dpi = c(800,800), pt.size = 2.5, max.cutoff = 'q99', min.cutoff = 'q1'
)  + coord_fixed()+ggtitle("TBX21 Regulon") +
  theme(axis.title.y=element_blank(), axis.text.y=element_blank(),axis.ticks.y=element_blank()) + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank()) 
p6 <- FeaturePlot(gdt_rna, features = c("NFAT5(+)"), raster = TRUE, cols = c('grey', 'red'),
                  raster.dpi = c(800,800), pt.size = 2.5, max.cutoff = 'q99', min.cutoff = 'q1'
)  + coord_fixed()+ggtitle("NFAT5 Regulon") +
  theme(axis.title.y=element_blank(), axis.text.y=element_blank(),axis.ticks.y=element_blank()) + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank()) 

p_combined <- plot_grid(p1, p2, p3, p4, p5, p6, nrow = 2, ncol=3)
ggsave(p_combined, filename = "Figures/FeaturePlot_Grid_SCENIC.pdf", device = "pdf", width = 10, height = 7)

#Violin plots of regulons 
gdt_seurat = gdt_rna
gdt_seurat$broad_class = factor(gdt_seurat$broad_class, levels = c("VD1+", "VD2+", "VD3+", "Mixed/Unknown"))
p1 = VlnPlot(gdt_seurat, features = "RORA(+)", group.by = "broad_class", split.by = "condition", 
             raster = F, pt.size = 0, cols = colors3) + stat_compare_means(method = "wilcox.test")
p2 = VlnPlot(gdt_seurat, features = "POLR3G(+)", group.by = "broad_class", split.by = "condition", 
             raster = F, pt.size = 0, cols = colors3) + stat_compare_means(method = "wilcox.test")
p3 = VlnPlot(gdt_seurat, features = "RUNX2(+)", group.by = "broad_class", split.by = "condition", 
             raster = F, pt.size = 0, cols = colors3) + stat_compare_means(method = "wilcox.test")
p4 = VlnPlot(gdt_seurat, features = "TBX21(+)", group.by = "broad_class", split.by = "condition", 
             raster = F, pt.size = 0, cols = colors3) + stat_compare_means(method = "wilcox.test")
p5 = VlnPlot(gdt_seurat, features = "RUNX3(+)", group.by = "broad_class", split.by = "condition", 
             raster = F, pt.size = 0, cols = colors3) + stat_compare_means(method = "wilcox.test")
p_combined <- plot_grid(p1, p2, p3, p4, p5, ncol=1)
ggsave(p_combined, filename = "Figures/ViolinPlots_Regulons_SCENIC_2.pdf", width = 7, height = 18)

#Module scores 
Memory <- c("IL7R", "CCR7", "ZNF683", "TCF7", "LEF1", "CXCR3")
Exhaustion <- c("PDCD1", "HAVCR2", "TIGIT", "LAG3", "CD244", "CD160")
Proliferation <- c("MKI67", "STMN1", "TYMS", "TOP2A", "PCNA", "CD74")
Activation <- c("IL2RB", "IL2RA", "HLA-DRA", "HLA-DRB1", "CD28", "CD27", "CD226", "TNFSF9")
Cytotoxicity <- c("GZMH", "NKG7", "KLRD1", "GZMK", "GNLY", "PRF1")

gene_modules = list(Memory, Exhaustion, Proliferation, Activation, Cytotoxicity)
names(gene_modules) = c("Memory", "Exhaustion", "Proliferation", "Activation", "Cytotoxicity")

DefaultAssay(gdt_seurat) = "RNA"
gdt_seurat = AddModuleScore(gdt_seurat, features = gene_modules, name = "ActivityModule")
p1 = VlnPlot(gdt_seurat, features = "ActivityModule1", group.by = "broad_class", split.by = "condition", 
             raster = F, pt.size = 0, cols = colors3) + stat_compare_means(method = "wilcox.test") + ggtitle("Memory Score")
p2 = VlnPlot(gdt_seurat, features = "ActivityModule2", group.by = "broad_class", split.by = "condition", 
             raster = F, pt.size = 0, cols = colors3) + stat_compare_means(method = "wilcox.test")+ ggtitle("Exhaustion Score")
p3 = VlnPlot(gdt_seurat, features = "ActivityModule3", group.by = "broad_class", split.by = "condition", 
             raster = F, pt.size = 0, cols = colors3) + stat_compare_means(method = "wilcox.test")+ ggtitle("Proliferation Score")
p4 = VlnPlot(gdt_seurat, features = "ActivityModule4", group.by = "broad_class", split.by = "condition", 
             raster = F, pt.size = 0, cols = colors3) + stat_compare_means(method = "wilcox.test")+ ggtitle("Activation Score")
p5 = VlnPlot(gdt_seurat, features = "ActivityModule5", group.by = "broad_class", split.by = "condition", 
             raster = F, pt.size = 0, cols = colors3) + stat_compare_means(method = "wilcox.test")+ ggtitle("Cytotoxicity Score")
p_combined <- plot_grid(p1, p2, p3, p4, p5, ncol=1)
ggsave(p_combined, filename = "Figures/ViolinPlots_GeneModules.pdf", width = 7, height = 18)

wilcox.test(subset(gdt_subset, subset=broad_class=="VD1+")$ActivityModule1, subset(gdt_subset, subset=broad_class=="VD2+")$ActivityModule1)
wilcox.test(subset(gdt_subset, subset=broad_class=="VD1+")$ActivityModule2, subset(gdt_subset, subset=broad_class=="VD2+")$ActivityModule2)
wilcox.test(subset(gdt_subset, subset=broad_class=="VD1+")$ActivityModule3, subset(gdt_subset, subset=broad_class=="VD2+")$ActivityModule3)
wilcox.test(subset(gdt_subset, subset=broad_class=="VD1+")$ActivityModule4, subset(gdt_subset, subset=broad_class=="VD2+")$ActivityModule4)
wilcox.test(subset(gdt_subset, subset=broad_class=="VD1+")$ActivityModule5, subset(gdt_subset, subset=broad_class=="VD2+")$ActivityModule5)

v#Try AUCell
expr.matrix = GetAssayData(gdt_seurat, assay = "RNA", layer = "counts")
cells_rankings <- AUCell_buildRankings(expr.matrix, nCores = 12)
#saveRDS(cells_rankings, file="Cell_Rankings_AUCell.rds")
cells_AUC <- AUCell_calcAUC(gene_modules, cells_rankings, 
                                     aucMaxRank=nrow(cells_rankings)*0.05, verbose = T)
gdt_seurat[["AUCell"]] <- CreateAssayObject(cells_AUC@assays@data$AUC)
DefaultAssay(gdt_seurat) = "AUCell"
rownames(gdt_seurat)
p1 = VlnPlot(gdt_seurat, features = "Memory", group.by = "broad_class", split.by = "condition", 
             raster = F, pt.size = 0, cols = colors3) + stat_compare_means(method = "wilcox.test") + ggtitle("Memory Score")
p2 = VlnPlot(gdt_seurat, features = "Exhaustion", group.by = "broad_class", split.by = "condition", 
             raster = F, pt.size = 0, cols = colors3) + stat_compare_means(method = "wilcox.test")+ ggtitle("Exhaustion Score")
p3 = VlnPlot(gdt_seurat, features = "Proliferation", group.by = "broad_class", split.by = "condition", 
             raster = F, pt.size = 0, cols = colors3) + stat_compare_means(method = "wilcox.test")+ ggtitle("Proliferation Score")
p4 = VlnPlot(gdt_seurat, features = "Activation", group.by = "broad_class", split.by = "condition", 
             raster = F, pt.size = 0, cols = colors3) + stat_compare_means(method = "wilcox.test")+ ggtitle("Activation Score")
p5 = VlnPlot(gdt_seurat, features = "Cytotoxicity", group.by = "broad_class", split.by = "condition", 
             raster = F, pt.size = 0, cols = colors3) + stat_compare_means(method = "wilcox.test")+ ggtitle("Cytotoxicity Score")
p_combined <- plot_grid(p1, p2, p3, p4, p5, ncol=1)
ggsave(p_combined, filename = "Figures/ViolinPlots_AUCell_Sets.pdf", width = 7, height = 18)


#GO enrichment of regulon targets
DefaultAssay(gdt_seurat) = "SCENIC"
regulon_targets = readRDS("pySCENIC_Output/Regulon_Target_Genes.rds")
symbols = regulon_targets$`POLR3G(+)`
entrezids = mapIds(org.Hs.eg.db, key = symbols, keytype = "SYMBOL", column = "ENTREZID")
entrezids = unique(entrezids)
ego <- enrichGO(gene = entrezids,
                OrgDb = org.Hs.eg.db,
                ont = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.05,
                readable = TRUE)
head(ego)
result = as.data.frame(ego)
head(result)

go_bp_list <- result[order(result$qvalue), ]
go_bp_list <- go_bp_list[go_bp_list$qvalue < 0.05, ]
go_bp_list$logqvalue <- -log10(go_bp_list$qvalue)
go_bp_list$Description <- factor(
  go_bp_list$Description, 
  levels = go_bp_list$Description[order(go_bp_list$logqvalue)]
)
#write.table(go_bp_list, file = "POLR3G_Target_Pathways.txt", sep = '\t', quote = F)

eomes = c('regulation of T cell activation', 'positive regulation of cell activation', 
          'leukocyte cell-cell adhesion', 'natural killer cell mediated cytotoxicity', 'alpha-beta T cell activation')
go_bp_list <- go_bp_list[go_bp_list$Description %in% eomes, ]
cebpd = c('regulation of leukocyte differentiation', 'positive regulation of cell activation', 'leukocyte proliferation', 
          'positive regulation of lymphocyte activation', 'leukocyte cell-cell adhesion')
go_bp_list <- go_bp_list[go_bp_list$Description %in% cebpd, ]
runx2 = c('lymphocyte differentiation', 'T cell differentiation', 'regulation of T cell activation', 
          'leukocyte cell-cell adhesion', 'positive regulation of lymphocyte proliferation', 
          'positive regulation of cytokine production')
go_bp_list <- go_bp_list[go_bp_list$Description %in% runx2, ]
polr3g = c('natural killer cell chemotaxis', 'positive regulation of lymphocyte chemotaxis', 
          'positive regulation of lymphocyte migration', 'chemokine-mediated signaling pathway', 
          'cellular response to type II interferon')
go_bp_list <- go_bp_list[go_bp_list$Description %in% polr3g, ]


go_bp_list$Description <- factor(
  go_bp_list$Description, 
  levels = go_bp_list$Description[order(go_bp_list$logqvalue)]
)
p1 <- ggplot(go_bp_list, aes(x = Description, y = logqvalue)) +
  geom_bar(stat = "identity", fill = "orange") +
  labs(x = "", y = "-log10(FDR)") +
  theme_bw() + coord_flip()
p1
ggsave("Figures/POLR3G_Selected_Targets.pdf", plot = p1, device = "pdf", height = 1.5, width = 4)

#Trying diffusion components 
library(destiny)
library(SingleCellExperiment)

dm_data <- Embeddings(gdt_seurat, reduction = 'integrated.rpca')
all.dm <- DiffusionMap(dm_data, verbose = TRUE)
saveRDS(all.dm, "All_Diffusion.RDS")
gdt_seurat[["DM"]] <- CreateDimReducObject(embeddings = all.dm@eigenvectors, key = "DM_")
DimPlot(gdt_seurat, group.by = "broad_class", reduction = "DM") + coord_fixed()
FeaturePlot(gdt_seurat, features = "CytoTRACE_score", reduction = "DM") + coord_fixed() +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

#Correlations with cytotrace
DefaultAssay(gdt_seurat) = "SCENIC"
Idents(gdt_seurat) = "broad_class"
scenic_matrix = as.data.frame(t(GetAssayData(gdt_seurat, assay = "SCENIC", slot = "data")))

tf_names <- c()
correlation_r <- c()
p_values <- c()
for(tf in rownames(gdt_seurat)){
  correlation <- cor.test(scenic_matrix[[tf]], gdt_seurat$CytoTRACE_score)
  tf_names <- c(tf_names, tf)
  correlation_r <- c(correlation_r, correlation$estimate)
  p_values <- c(p_values, correlation$p.value)
}
correlation_results <- data.frame(
  TF = tf_names,
  Correlation_r = correlation_r,
  P_value = p_values
)

correlation_results <- correlation_results %>%
  arrange(Correlation_r)
top_bottom_5 <- correlation_results %>%
  filter(row_number() %in% c(1:5, (n()-4):n()))
waterfall_plot <- ggplot(correlation_results, aes(x = reorder(TF, -Correlation_r), y = Correlation_r)) +
  geom_point(size = 3, color = "steelblue") +
  geom_text_repel(data = top_bottom_5, aes(label = TF), size = 3.5, max.overlaps = 20) +
  theme_bw() +
  labs(title = "Correlations with CytoTRACE", x = "Transcription Factors", y = "Correlation Coefficient (r)") +
  theme(axis.text.x = element_blank(),  
        axis.ticks.x = element_blank()) 
print(waterfall_plot)



write.table()
