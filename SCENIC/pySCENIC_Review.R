library(Seurat)
library(tidyverse)
library(data.table)
library(SCopeLoomR)
library(SeuratData)
library(SeuratDisk)
library(tidyverse)
library(ggridges)
library(heatmaply)
library(AUCell)
library(ggpubr)
library(gridExtra)
library(pheatmap)
library(Matrix)
library(SCENIC)

setwd("/mnt/isilon/tan_lab/sussmanj/GIST_Multiome/GD_TCells")

out_dir <- "pySCENIC_Output"
out_files <- list.files(out_dir, full.names = TRUE)
seurat_path <- "/mnt/isilon/tan_lab/sussmanj/GIST_Multiome/GD_TCells/GDT_All_Processed_Revised13_GDT_QC_Integrated2.RDS"

analysis <- "GDT"
path2loom <- out_files[grep("pySCENIC_AUCell_Revision.loom", out_files)]
loom <- open_loom(path2loom, mode = "r")
regulonsMat <- get_regulons(loom, column.attr.name = "Regulons")
dim(regulonsMat)
write.csv(regulonsMat,file = paste0(out_dir, "/", analysis, "_pySCENIC_regulons_incidence_matrix_revision.csv"))
regulons <- SCENIC::regulonsToGeneLists(regulonsMat) 
head(regulons)
saveRDS(regulons, file = "pySCENIC_Output/Regulon_Target_Genes_Revision.rds")
regulonsAUC <- get_regulons_AUC(loom, column.attr.name = 'RegulonsAUC')
saveRDS(regulonsAUC, file = "regulonsAUC_revision.rds")
close_loom(loom)

### Add desired data to original Seurat object (RNA only)
gdt_seurat <- readRDS(seurat_path)
AUCmat <- AUCell::getAUC(regulonsAUC)
AUCmat[1:10, 1:10]
rownames(AUCmat)
gdt_seurat[['SCENIC']] <- CreateAssayObject(data = AUCmat)
DefaultAssay(gdt_seurat) <- "SCENIC"
gdt_seurat <- ScaleData(gdt_seurat, assay = 'SCENIC', features = rownames(AUCmat))
saveRDS(gdt_seurat, file = paste0("/mnt/isilon/tan_lab/sussmanj/GIST_Multiome/GD_TCells/GDT_All_Processed_Revised13_GDT_QC_Integrated2_SCENIC.RDS"))

###
#START FROM HERE#
gdt_rna <- readRDS("/mnt/isilon/tan_lab/sussmanj/GIST_Multiome/GD_TCells/GDT_All_Processed_Revised13_GDT_QC_Integrated2_SCENIC.RDS")

DefaultAssay(gdt_rna) = "SCENIC"
FeaturePlot(gdt_rna, features = c("TCF7(+)"), raster = FALSE, 
            min.cutoff = "q1", max.cutoff = "q99", reduction = "umap.rpca") + coord_fixed()
VlnPlot(gdt_rna, features = c("TCF7(+)"), group.by = 'cell_type', split.by = "condition", pt.size = 0)

Idents(gdt_rna) <- "cell_type"
deg <- FindAllMarkers(gdt_rna, logfc.threshold = 0.01, only.pos = TRUE)
write.table(deg, file = "SCENIC_Clusters_Differential_TF_Revision.tsv", sep = '\t', quote = F)

Idents(gdt_rna) <- "Source"
deg2 <- FindAllMarkers(gdt_rna, logfc.threshold = 0.01, only.pos = TRUE)
write.table(deg2, file = "SCENIC_Source_Differential_TF_Revision.tsv", sep = '\t', quote = F)

gdt_subset = subset(gdt_rna, subset = Source == "Tumor")
Idents(gdt_subset) = "broad_class"
detfs = FindAllMarkers(gdt_subset, logfc.threshold = 0.01, only.pos = TRUE)
write.table(detfs, "SCENIC_Regulons_Tumor_Broad_Class_Differential.txt", sep = '\t', quote = F)

gdt_subset = subset(gdt_rna, subset = Source == "PBMC")
Idents(gdt_subset) = "broad_class"
detfs = FindAllMarkers(gdt_subset, logfc.threshold = 0.01, only.pos = TRUE)
write.table(detfs, "SCENIC_Regulons_PBMC_Broad_Class_Differential.txt", sep = '\t', quote = F)

