library(Seurat)
library(survival)
library(dplyr)
library(limma)
library(edgeR)
library(ggpubr)
library(SeuratDisk)
library(dittoSeq)
library(DoubletFinder)

setwd("/mnt/isilon/tan_lab/sussmanj/GIST_Multiome/GD_TCells")
gdt_seurat <- readRDS("/mnt/isilon/tan_lab/sussmanj/GIST_Multiome/GD_TCells/Gamma_Delta_TCells.rds")

DimPlot(gdt_seurat) + coord_fixed()

#Doublet finder 
FindDoublets <- function(seurat.rna, PCs = 1:50, exp_rate = 0.075, sct = FALSE){
  
  ## pK identification
  sweep.res.list <- paramSweep(seurat.rna, PCs = PCs, sct = sct)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  
  ## Homotypic Doublet proportion Estimate
  annotations <- seurat.rna@meta.data$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
  nExp_poi <- round(exp_rate * length(seurat.rna$seurat_clusters))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  
  ## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
  seurat.rna <- doubletFinder(seurat.rna, PCs = PCs, pN = 0.25,
                              pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, 
                              sct = sct)
  
  seurat.rna <- doubletFinder(seurat.rna, PCs = PCs, pN = 0.25, 
                              pK = 0.09, nExp = nExp_poi.adj,
                              reuse.pANN = paste0("pANN_0.25_0.09_", nExp_poi), 
                              sct = sct)
  doublet_var = paste0('DF.classifications_0.25_0.09_', nExp_poi.adj)
  seurat.rna[['Doublet_Singlet']] = seurat.rna[[doublet_var]]
  
  mnames = names(seurat.rna@meta.data)
  seurat.rna@meta.data[, grep(mnames, pattern = '0.25_0.09')] <- NULL
  return(seurat.rna)
}

#Reprocess
DefaultAssay(gdt_seurat) = "RNA"
gdt_seurat[["integrated"]] <- NULL
gdt_seurat[["prediction.score.celltype.l1"]] <- NULL
gdt_seurat[["prediction.score.celltype.l2"]] <- NULL
gdt_seurat[["predicted_ADT"]] <- NULL
gdt_seurat[["SCT"]] <- NULL

ob.list <- SplitObject(gdt_seurat, split.by = "samples")
for (i in 1:length(ob.list)) {
  DefaultAssay(ob.list[[i]]) <- "RNA"
  ob.list[[i]][["percent.mt"]] <- PercentageFeatureSet(ob.list[[i]], pattern = "^MT-")
  ob.list[[i]] <- subset(ob.list[[i]], subset = nFeature_RNA > 100 & nFeature_RNA < 10000 & percent.mt < 10)
  ob.list[[i]] <- SCTransform(ob.list[[i]])
  ob.list[[i]] <- RunPCA(ob.list[[i]], features = VariableFeatures(object = ob.list[[i]]))
  ob.list[[i]] <- FindNeighbors(ob.list[[i]], dims = 1:30)
  ob.list[[i]] <- FindClusters(ob.list[[i]], resolution = 1)
  ob.list[[i]] <- FindDoublets(ob.list[[i]], PCs = 1:30, sct = TRUE, exp_rate = (length(colnames(ob.list[[i]]))/125000))
}
saveRDS(ob.list, file = "Object_List_Temp.RDS")

merged_gdt = merge(ob.list[[1]], y = ob.list[-1])
saveRDS(merged_gdt, file = "GDT_All_Processed.RDS")

#Process
merged_gdt <- readRDS("GDT_All_Processed.RDS")
merged_gdt <- subset(merged_gdt, nFeature_RNA > 300 & nCount_RNA > 1000 & nFeature_RNA < 10000 & percent.mt < 10)
merged_gdt <- subset(merged_gdt, subset = Doublet_Singlet == "Singlet")

DefaultAssay(merged_gdt) = "RNA"
merged_gdt[["SCT"]] = NULL
counts = GetAssayData(merged_gdt, assay = "RNA", slot = "counts")
merged_gdt[["RNA"]] = CreateAssay5Object(counts = counts)
merged_gdt

merged_gdt <- JoinLayers(merged_gdt)
merged_gdt <- SCTransform(merged_gdt, verbose = TRUE)
merged_gdt <- RunPCA(merged_gdt)
merged_gdt <- RunUMAP(merged_gdt, dims = 1:30, reduction.name = "umap.unintegrated")
DimPlot(merged_gdt, group.by = "orig.ident") + coord_fixed()
DimPlot(merged_gdt, group.by = "samples", label = TRUE) + coord_fixed()
FeaturePlot(merged_gdt, features = "percent.mt") + coord_fixed()

merged_gdt[["RNA"]] <- split(merged_gdt[["RNA"]], f = merged_gdt$samples)
merged_gdt <- SCTransform(merged_gdt, vars.to.regress = c("percent.mt"))
merged_gdt <- RunPCA(merged_gdt)
merged_gdt <- IntegrateLayers(object = merged_gdt, 
                                   method = RPCAIntegration, new.reduction = "integrated.rpca",
                                   normalization.method = "SCT")
merged_gdt <- RunUMAP(merged_gdt, dims = 1:30, reduction = "integrated.rpca", 
                           reduction.name = "umap.rpca")
merged_gdt <- CellCycleScoring(merged_gdt, s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes, set.ident = FALSE)
DimPlot(merged_gdt, group.by = "orig.ident", reduction = "umap.rpca") + coord_fixed()
DimPlot(merged_gdt, group.by = "samples", reduction = "umap.rpca") + coord_fixed()

DefaultAssay(merged_gdt) = "RNA"
merged_gdt = NormalizeData(merged_gdt)
FeaturePlot(merged_gdt, features = c("TRDV1", "TRDV2"), reduction = "umap.rpca") + coord_fixed()
VlnPlot(merged_gdt, features = c("TRDV1", "TRDV2", "CD3E", "TRGV9"), group.by = "SCT_snn_res.0.3", pt.size = 0, ncol = 2) 
Idents(merged_gdt) = "SCT_snn_res.0.3"
markers = FindAllMarkers(merged_gdt)
saveRDS(markers, "Markers_All_0.3.RDS")

#saveRDS(merged_gdt, "GDT_All_Integrated.RDS")

#Go back and check the cell cycle phase 
FeaturePlot(merged_gdt, features = "percent.mt", reduction = "umap.rpca") + coord_fixed()
DimPlot(merged_gdt, group.by = "Phase", reduction = "umap.rpca") + coord_fixed()

DefaultAssay(merged_gdt) <- "SCT"
merged_gdt <- FindNeighbors(merged_gdt, dims = 1:30, reduction = "integrated.rpca")
merged_gdt <- FindClusters(merged_gdt, resolution = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8))
DimPlot(merged_gdt, group.by = "SCT_snn_res.0.3", reduction = "umap.rpca") + coord_fixed()
DimPlot(merged_gdt, group.by = "Doublet_Singlet", reduction = "umap.rpca") + coord_fixed()
DimPlot(merged_gdt, group.by = "Phase", reduction = "umap.rpca") + coord_fixed()


#Project to PBMCs
#reference <- LoadH5Seurat("pbmc_multimodal.h5seurat")
#saveRDS(reference, file = "PBMC_Azimuth_Reference.rds")

reference = readRDS("PBMC_Azimuth_Reference.rds")
DimPlot(object = reference, reduction = "wnn.umap", group.by = "celltype.l2",
        label = TRUE, label.size = 3, repel = TRUE, raster = F, pt.size = 0.15, label.box = T) + NoLegend()

DefaultAssay(merged_gdt) = "RNA"
anchors <- FindTransferAnchors(
  reference = reference,
  query = merged_gdt,
  normalization.method = "SCT",
  reference.reduction = "spca",
  dims = 1:50
)

merged_gdt <- MapQuery(
  anchorset = anchors,
  query = merged_gdt,
  reference = reference,
  refdata = list(
    celltype.l1 = "celltype.l1",
    celltype.l2 = "celltype.l2",
    predicted_ADT = "ADT"
  ),
  reference.reduction = "spca", 
  reduction.model = "wnn.umap"
)
#saveRDS(merged_gdt, "GDT_All_Integrated_Projected.RDS")


DimPlot(merged_gdt, group.by = "predicted.celltype.l2", reduction = "umap.rpca", 
        label = T, label.box = T, repel = T, pt.size = 0.5) + coord_fixed()
table(merged_gdt$predicted.celltype.l2)

FeaturePlot(merged_gdt, features = "predicted.celltype.l2.score", reduction = "umap") + coord_fixed()
VlnPlot(merged_gdt, features = "predicted.celltype.l2.score", group.by = "predicted.celltype.l2", pt.size = 0)

FeaturePlot(merged_gdt, features = c("CD3D", "CD3E", "TRDC"), reduction = "umap.rpca")
dittoBarPlot(merged_gdt, var = "predicted.celltype.l2", group.by = "cell_type")



