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
seurat_path <- "GDT_Seurat_Annotated_CytoTRACE_05262024.RDS"

analysis <- "GDT"
path2loom <- out_files[grep("pySCENIC_AUCell.loom", out_files)]
loom <- open_loom(path2loom, mode = "r")
regulonsMat <- get_regulons(loom, column.attr.name = "Regulons")
dim(regulonsMat)
write.csv(regulonsMat,file = paste0(out_dir, "/", analysis, "_pySCENIC_regulons_incidence_matrix.csv"))
regulons <- SCENIC::regulonsToGeneLists(regulonsMat) 
head(regulons)
regulonsAUC <- get_regulons_AUC(loom, column.attr.name = 'RegulonsAUC')
saveRDS(regulonsAUC, file = "regulonsAUC.rds")
close_loom(loom)

### Add desired data to original Seurat object (RNA only)
gdt_rna <- readRDS(seurat_path)
AUCmat <- AUCell::getAUC(regulonsAUC)
AUCmat[1:10, 1:10]
rownames(AUCmat)
gdt_rna[['SCENIC']] <- CreateAssayObject(data = AUCmat)
DefaultAssay(gdt_rna) <- "SCENIC"
gdt_rna <- ScaleData(gdt_rna, assay = 'SCENIC', features = rownames(AUCmat))
saveRDS(gdt_rna, file = paste0("GDT_Seurat_Annotated_CytoTRACE_SCENIC_05292024.RDS"))

###
#START FROM HERE#
gdt_rna <- readRDS("GDT_Seurat_Annotated_CytoTRACE_SCENIC_05292024.RDS")

DefaultAssay(gdt_rna) = "SCENIC"
FeaturePlot(gdt_rna, features = c("TCF7(+)"), raster = FALSE, 
            min.cutoff = "q1", max.cutoff = "q99", reduction = "umap.rpca") + coord_fixed()
VlnPlot(gdt_rna, features = c("TCF7(+)"), group.by = 'cell_type', split.by = "condition", pt.size = 0)

Idents(gdt_rna) <- "cell_type"
deg <- FindAllMarkers(gdt_rna, logfc.threshold = 0.01, only.pos = TRUE)
write.table(deg, file = "SCENIC_Clusters_Differential_TF.tsv", sep = '\t', quote = F)

Idents(gdt_rna) <- "condition"
deg <- FindAllMarkers(gdt_rna, logfc.threshold = 0.01, only.pos = TRUE)
write.table(deg, file = "SCENIC_Condition_Differential_TF.tsv", sep = '\t', quote = F)

