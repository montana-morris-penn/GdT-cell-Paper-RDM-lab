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

setwd("/mnt/isilon/tan_lab/sussmanj/GIST_Multiome/GD_TCells")
gdt_seurat <- readRDS("/mnt/isilon/tan_lab/sussmanj/GIST_Multiome/GD_TCells/GDT_Seurat_Annotated_CytoTRACE_05262024.RDS")
gdt_seurat <- readRDS("/mnt/isilon/tan_lab/sussmanj/GIST_Multiome/GD_TCells/merged_gdt_res1.rds")

#UMAPs
colors = DiscretePalette(14, palette = "stepped", shuffle = FALSE)
names = levels(Idents(gdt_seurat))
names(colors) = names

p1 = DimPlot(gdt_seurat, raster = T, raster.dpi = c(3000,3000), pt.size = 4.5, cols = colors, label = T, label.box = T, repel = T) + coord_fixed()
ggsave(p1, filename = "Figures/UMAP_Clusters.pdf", width = 7, height = 7)

p2 = DimPlot(gdt_seurat, raster = T, raster.dpi = c(3000,3000), pt.size = 4.5, group.by = "condition", cols = colors3) + coord_fixed()
ggsave(p2, filename = "Figures/UMAP_Condition.pdf", width = 7, height = 7)

p3 = DimPlot(gdt_seurat, raster = T, raster.dpi = c(3000,3000), pt.size = 4.5, group.by = "samples", cols = colors4) + coord_fixed()
ggsave(p3, filename = "Figures/UMAP_Sample.pdf", width = 7, height = 7)

gdt_seurat$cell_type = Idents(gdt_seurat)
table(gdt_seurat$cell_type)
Idents(gdt_seurat) = gdt_seurat$cell_type
gdt_seurat = RenameIdents(gdt_seurat, "Vd2+ 1" = "VD2+", 
                          "Vd2+ 2" = "VD2+", 
                          "Vd2+ 3" = "VD3+", 
                          "Vd2+ 4" = "VD2+",
                          "Vd2+ 5" = "VD2+",
                          "Vd1+ 1" = "VD1+",
                          "Vd1+ 2" = "VD1+",
                          "Vd1+ 3" = "VD1+",
                          "Vd1+ 4" = "VD1+",
                          "Vd1+ 5" = "VD1+",
                          "Vd1+ 6" = "VD1+",
                          "Vd1+ 7" = "VD1+", 
                          "mixed gdt" = "Mixed/Unknown", 
                          "Unknown" = "Mixed/Unknown")
gdt_seurat$broad_class = Idents(gdt_seurat)

colors2 = DiscretePalette(4, palette = "polychrome", shuffle = TRUE)
names(colors2) = c("VD2+", "VD1+", "VD3+", "Mixed/Unknown")
saveRDS(colors2, "Figures/Broad_Class_Colors.rds")
p4 = DimPlot(gdt_seurat, raster = T, raster.dpi = c(3000,3000), pt.size = 4.5, group.by = "broad_class", cols = colors2, reduction = "umap.rpca") + coord_fixed()
ggsave(p4, filename = "Figures/UMAP_VD1_VD2.pdf", width = 7, height = 7)

source("/mnt/isilon/tan_lab/sussmanj/Temp/GBM_Immune/Figure_Functions.R")
gdt_seurat$condition = factor(gdt_seurat$condition)
gdt_seurat$broad_class = factor(gdt_seurat$broad_class, levels = c("VD1+", "VD2+", "VD3+", "Mixed/Unknown"))
colors2 = colors2[levels(gdt_seurat$broad_class)]
barplots(gdt_seurat@meta.data, x="condition", y="broad_class", width=8, height=10,
         colors = colors2, filename = "Figures/Barplot_Condition_Broad_Class.pdf", xlabs = levels(gdt_seurat$condition))

gdt_seurat$cell_type = factor(gdt_seurat$cell_type, levels = c("Vd1+ 1", "Vd1+ 2", "Vd1+ 3", "Vd1+ 4", "Vd1+ 5", "Vd1+ 6", "Vd1+ 7", 
                              "Vd2+ 1", "Vd2+ 2", "Vd2+ 4", "Vd2+ 5", "Vd2+ 3", "mixed gdt", "Unknown"))
DimPlot(gdt_seurat, cols = colors, group.by = "cell_type", 
        label = T, label.box = T, repel = T) + coord_fixed()

gdt_seurat$condition = factor(gdt_seurat$condition)
colors3 = DiscretePalette(3, palette = "polychrome", shuffle = TRUE)
names(colors3) = c("P", "T")
colors3 = colors3[1:2]
#saveRDS(colors3, "Figures/Condition_Colors.rds")
colors3 = readRDS("Figures/Condition_Colors.rds")
barplots(gdt_seurat@meta.data, x="cell_type", y="condition", width=14, height=8,
         colors = colors3, filename = "Figures/Barplot_Condition_CellType.pdf", xlabs = levels(gdt_seurat$cell_type))

gdt_seurat$samples = factor(gdt_seurat$samples)
colors4 = DiscretePalette(6, palette = "polychrome", shuffle = TRUE)
names(colors4) = levels(gdt_seurat$samples)
#saveRDS(colors4, "Figures/Sample_Colors.rds")
colors4 = readRDS("Figures/Sample_Colors.rds")
barplots(gdt_seurat@meta.data, x="cell_type", y="samples", width=14, height=8,
         colors = colors4, filename = "Figures/Barplot_Sample_CellType.pdf", xlabs = levels(gdt_seurat$cell_type))

#Featureplots 
p1 <- FeaturePlot(gdt_seurat, features = c("TRDV1"), raster = TRUE, 
                   raster.dpi = c(800,800), pt.size = 2.5, max.cutoff = 'q99'
                  )  + coord_fixed() +
  theme(axis.title.y=element_blank(), axis.text.y=element_blank(),axis.ticks.y=element_blank()) + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank()) 
p2 <- FeaturePlot(gdt_seurat, features = c("TRDV2"), raster = TRUE, 
                   raster.dpi = c(800,800), pt.size = 2.5, max.cutoff = 'q99'
                  ) + coord_fixed()+
  theme(axis.title.y=element_blank(), axis.text.y=element_blank(),axis.ticks.y=element_blank()) + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank()) 
p3 <- FeaturePlot(gdt_seurat, features = c("TRGV9"), raster = TRUE, 
                   raster.dpi = c(800,800), pt.size = 2.5, max.cutoff = 'q99'
                  )  + coord_fixed()+
  theme(axis.title.y=element_blank(), axis.text.y=element_blank(),axis.ticks.y=element_blank()) + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank()) 
p4 <- FeaturePlot(gdt_seurat, features = c("TRDV3"), raster = TRUE, 
                   raster.dpi = c(800,800), pt.size = 2.5, max.cutoff = 'q99'
                  )  + coord_fixed()+
  theme(axis.title.y=element_blank(), axis.text.y=element_blank(),axis.ticks.y=element_blank()) + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank()) 
p5 <- FeaturePlot(gdt_seurat, features = c("TRGV8"), raster = TRUE, 
                   raster.dpi = c(800,800), pt.size = 2.5, max.cutoff = 'q99'
                  )  + coord_fixed()+
  theme(axis.title.y=element_blank(), axis.text.y=element_blank(),axis.ticks.y=element_blank()) + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank()) 
p6 <- FeaturePlot(gdt_seurat, features = c("TRDC"), raster = TRUE, 
                   raster.dpi = c(800,800), pt.size = 2.5, max.cutoff = 'q99'
                  )  + coord_fixed()+
  theme(axis.title.y=element_blank(), axis.text.y=element_blank(),axis.ticks.y=element_blank()) + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank()) 
p7 <- FeaturePlot(gdt_seurat, features = c("KLRC1"), raster = TRUE, 
                   raster.dpi = c(800,800), pt.size = 2.5, max.cutoff = 'q99'
                  )  + coord_fixed()+
  theme(axis.title.y=element_blank(), axis.text.y=element_blank(),axis.ticks.y=element_blank()) + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank()) 
p8 <- FeaturePlot(gdt_seurat, features = c("KLRC2"), raster = TRUE, 
                   raster.dpi = c(800,800), pt.size = 2.5, max.cutoff = 'q99'
                  )  + coord_fixed()+
  theme(axis.title.y=element_blank(), axis.text.y=element_blank(),axis.ticks.y=element_blank()) + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank()) 
p9 <- FeaturePlot(gdt_seurat, features = c("KLRB1"), raster = TRUE, 
                   raster.dpi = c(800,800), pt.size = 2.5, max.cutoff = 'q99'
                  )  + coord_fixed()+
  theme(axis.title.y=element_blank(), axis.text.y=element_blank(),axis.ticks.y=element_blank()) + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank()) 
p10 <- FeaturePlot(gdt_seurat, features = c("GZMK"), raster = TRUE, 
                    raster.dpi = c(800,800), pt.size = 2.5, max.cutoff = 'q99'
                   )  + coord_fixed()+
  theme(axis.title.y=element_blank(), axis.text.y=element_blank(),axis.ticks.y=element_blank()) + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank()) 
p11 <- FeaturePlot(gdt_seurat, features = c("GZMH"), raster = TRUE, 
                    raster.dpi = c(800,800), pt.size = 2.5, max.cutoff = 'q99'
                   )  + coord_fixed()+
  theme(axis.title.y=element_blank(), axis.text.y=element_blank(),axis.ticks.y=element_blank()) + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank()) 
p12 <- FeaturePlot(gdt_seurat, features = c("CD8A"), raster = TRUE, 
                    raster.dpi = c(800,800), pt.size = 2.5, max.cutoff = 'q99'
                   ) + coord_fixed()+
  theme(axis.title.y=element_blank(), axis.text.y=element_blank(),axis.ticks.y=element_blank()) + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank()) 
p13 <- FeaturePlot(gdt_seurat, features = c("TBX21"), raster = TRUE, 
                    raster.dpi = c(800,800), pt.size = 2.5, max.cutoff = 'q99'
                   ) + coord_fixed()+
  theme(axis.title.y=element_blank(), axis.text.y=element_blank(),axis.ticks.y=element_blank()) + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank()) 
p14 <- FeaturePlot(gdt_seurat, features = c("RORA"), raster = TRUE, 
                    raster.dpi = c(800,800), pt.size = 2.5, max.cutoff = 'q99'
                   )  + coord_fixed()+
  theme(axis.title.y=element_blank(), axis.text.y=element_blank(),axis.ticks.y=element_blank()) + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank()) 
p15 <- FeaturePlot(gdt_seurat, features = c("CEBPD"), raster = TRUE, 
                    raster.dpi = c(800,800), pt.size = 2.5, max.cutoff = 'q99'
                   )  + coord_fixed()+
  theme(axis.title.y=element_blank(), axis.text.y=element_blank(),axis.ticks.y=element_blank()) + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank()) 
p16 <- FeaturePlot(gdt_seurat, features = c("ZEB2"), raster = TRUE, 
                    raster.dpi = c(800,800), pt.size = 2.5, max.cutoff = 'q99'
                   )  + coord_fixed()+
  theme(axis.title.y=element_blank(), axis.text.y=element_blank(),axis.ticks.y=element_blank()) + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank()) 
p_combined <- plot_grid(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, 
                        p11, p12, p13, p14, p15, p16, nrow = 4, ncol=4)
ggsave(p_combined, filename = "Figures/FeaturePlot_Grid_2.pdf", device = "pdf", width = 16, height = 16)


##Differential Gene analysis
DefaultAssay(gdt_seurat) = "RNA"
gdt_seurat = JoinLayers(gdt_seurat)
gdt_seurat = NormalizeData(gdt_seurat)
#saveRDS(gdt_seurat, "GDT_Seurat_Annotated_05262024.RDS")

groupby = "broad_class"
compartment = "T"
gdt_subset = subset(gdt_seurat, subset=condition==compartment)
comparisons = c("VD1+", "VD2+")

groupby = "condition"
compartment = "VD1+"
gdt_subset = subset(gdt_seurat, subset=broad_class==compartment)
comparisons = c("P", "T")

groupby = "condition"
compartment = "All"
gdt_subset = gdt_seurat
comparisons = c("T", "P")

group1 = comparisons[1]
group2= comparisons[2]

print(compartment)
print(paste0(group1, "_vs_", group2))
de.genes = FindMarkers(gdt_subset, 
                       verbose = T,
                       assay = "RNA", 
                       logfc.threshold = 0, 
                       min.pct = 0.001, 
                       ident.1 = group1, 
                       ident.2 = group2, 
                       group.by = groupby,
                       max.cells.per.ident = 1500)
de.genes$group = "de.genes"
de.genes$gene = rownames(de.genes)
de.genes$enriched = ifelse(de.genes$avg_log2FC < 0, yes = group2, no = group1)

gene.stringent = ggplot(de.genes %>% filter(abs(avg_log2FC) > 0.5, p_val_adj < 1e-10), 
                                 aes(x = avg_log2FC, y = -log(p_val_adj+1e-300, base = 10),
                                     color = enriched,
                                     label = gene)) + 
  geom_point(data = de.genes, color = "gray")+
  geom_point()+
  geom_text_repel(size = 3, max.overlaps = 50) + 
  theme_bw() + 
  labs(x = "Log2FC", 
       y = "-log(pv_adjust)") + #xlim(c(-3, 3))+ ylim(c(2, 300)) +
  scale_color_manual(values =  c("#BD20BF","#9DA52F"))+
  theme(legend.position = "none")
gene.stringent
ggsave(gene.stringent, filename = paste0("Figures/",group1,"_vs_",group2,"_",compartment,"_Volcano.pdf"), width = 6, height = 6)

#GSEA
res <- de.genes %>% dplyr::select(gene, avg_log2FC) 

# de.genes$rank = sign(de.genes$avg_log2FC)*-log10(de.genes$p_val_adj)
# res <- de.genes %>% dplyr::select(gene, rank)

ranks <- deframe(res)

pathways.kegg <- msigdbr(species = "human", category = "C2", subcategory = "KEGG")
pathways.hallmark <- msigdbr(species = "human", category = "H")

pathway_list_kegg = split(x = pathways.kegg$gene_symbol, f = pathways.kegg$gs_name)
pathway_list_hallmark = split(x = pathways.hallmark$gene_symbol, f = pathways.hallmark$gs_name)
pathway_list <- c(pathway_list_kegg, pathway_list_hallmark)

GSEA_result <-  fgsea(pathways = pathway_list, stats = ranks)
write.table(as.data.frame(GSEA_result[,-"leadingEdge"]), file = paste0("Figures/",group1,"_vs_",group2,"_",compartment,"_GSEA_Barplot.txt"), sep = '\t', quote = F)

sig_pathways <- GSEA_result[GSEA_result$padj<0.10, ]
sig_pathways <- sig_pathways[order(sig_pathways$NES), ]
sig_pathways$pathway <- factor(
  sig_pathways$pathway, 
  levels = sig_pathways$pathway[order(sig_pathways$NES)]
)

selected_pathways = c("KEGG_CYTOKINE_CYTOKINE_RECEPTOR_INTERACTION", 
                      "HALLMARK_INFLAMMATORY_RESPONSE",
                      "HALLMARK_IL6_JAK_STAT3_SIGNALING",
                      "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", 
                      "KEGG_COMPLEMENT_AND_COAGULATION_CASCADES", 
                      "HALLMARK_IL2_STAT5_SIGNALING", 
                      "KEGG_CHEMOKINE_SIGNALING_PATHWAY",
                      "HALLMARK_TNFA_SIGNALING_VIA_NFKB", 
                      "HALLMARK_ADIPOGENESIS", 
                      "HALLMARK_INTERFERON_GAMMA_RESPONSE")
sig_pathways_filter = sig_pathways[sig_pathways$pathway %in% selected_pathways, ]

top_pathways <- head(sig_pathways[sig_pathways$NES<0, ], 20)
bottom_pathways <- tail(sig_pathways[sig_pathways$NES>0, ], 20)
sig_pathways_filter <- rbind(top_pathways, bottom_pathways)

sig_pathways_filter$pathway <- factor(
  sig_pathways_filter$pathway,
  levels = sig_pathways_filter$pathway[order(sig_pathways_filter$NES)]
)
p1 <- ggplot(sig_pathways_filter, aes(x = pathway, y = NES, fill = -log10(padj))) +
  geom_bar(stat = "identity", position = "identity") +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 4.5) +
  labs(x = "", y = "Normalized Enrichment Score", fill = "FDR") +  # Update axis and legend labels
  theme_bw() + coord_flip()+
  theme(
    axis.text.x = element_text(vjust = 0, hjust = 0.5, family = "sans", size = 8, color = "black"),
    axis.text.y = element_text(family = "sans", size = 8, color = "black"),
    axis.line = element_blank(),
    plot.margin = margin(l = 5, r = 5, b = 5, t = 5, unit = "pt"), 
    panel.grid = element_blank())
p1
ggsave(paste0("Figures/",group1,"_vs_",group2,"_",compartment,"_GSEA_Barplot.pdf"), 
       plot = p1, device = 'pdf', height = 2, width = 6)

#GO enrichment 
symbols = de.genes[de.genes$enriched == "T", ]$gene
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
dotplot(ego, showCategory=20) + ggtitle("Tumor Enrichment")

#Dotplot of genes
genes <- c("TRDV1", "TRDV3", "TRGV8", "TRGV3", "TRGC2", "GZMH", "NKG7", "KLRD1", "KLRK1", "GNLY",
           "KLRF1", "FCGR3A", "FCRL6", "NCR1", "CD244", "CD8A", "KLRC2", "CX3CR1", "HLA-DRB1", "IFNG",
           "TIGIT", "TBX21", "EOMES", "TYROBP", "ZEB2", "ZNF683", "IL2RB", "IKZF2", "AOAH", "PRKCH",
           "PIK3AP1", "BCL2", "RUNX1", "TRDV2", "TRGV9", "CEBPD", "KLRB1", "GZMK", "KLRC1", "IL7R",
           "IL23R", "CD44", "CD69", "ZBTB16", "RORA", "CAMK4", "BTN3A1", "ANXA1", "IL17RA", "IFNGR1",
           "CXCR6", "AREG", "THEMIS", "TCF7", "LEF1", "LTB", "SYNE2", "XCL1", "SELL", "CCL4", "CCR7",
           "KIR2DL1", "MAL", "CDC20", "AURKB", "MKI67", "CDK1")

features = unique(genes)
#gdt_seurat$cell_type = Idents(gdt_seurat)
gdt_seurat$cell_type = factor(gdt_seurat$cell_type, levels = rev(c("Vd1+ 1", "Vd1+ 2", "Vd1+ 3", "Vd1+ 4", "Vd1+ 5", "Vd1+ 6", "Vd1+ 7", 
                                                               "Vd2+ 1", "Vd2+ 2", "Vd2+ 4", "Vd2+ 5", "Vd2+ 3", "mixed gdt", "Unknown")))
DefaultAssay(gdt_seurat) = "RNA"
gdt_seurat = JoinLayers(gdt_seurat)
gdt_seurat = NormalizeData(gdt_seurat)

pdot = DotPlot(gdt_seurat, features = features, group.by = "cell_type", 
        assay = "RNA") + scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) + 
        RotatedAxis()
pdot
ggsave(pdot, filename = "Figures/Dotplot.pdf", width = 18, height = 3.5)



