library(Seurat)
library(tidyverse)
library(magrittr)
library(SingleCellExperiment)
library(SeuratObject)
library(SeuratDisk)
library(SCopeLoomR)


seurat.rna <- readRDS('/mnt/isilon/tan_lab/sussmanj/GIST_Multiome/GD_TCells/merged_gdt_res1.rds')
DefaultAssay(seurat.rna) <- "RNA"
seurat.rna = JoinLayers(seurat.rna)
seurat.rna = NormalizeData(seurat.rna)


DimPlot(seurat.rna, group.by = "condition") + coord_fixed()
FeaturePlot(seurat.rna, features = c("TRDV1", "TRDV2"), reduction = "umap")

dgem <- GetAssayData(seurat.rna, assay = "RNA", layer = "counts")
dim(dgem)
cellnames = colnames(dgem)

DimPlot(seurat.rna, cells = cellnames) + coord_fixed()
DimPlot(seurat.rna, cells = cellnames, reduction = "umap.rpca") + coord_fixed()

#############################
# Compute logical vector for every gene reflecting whether
# it has more than zero count per cell & and is detected in at 1 per thousand of cells in the global dataset
c <- floor(dim(dgem)[2]*.001); c  
nonzero <- dgem > 0L
keep_genes <- rowSums(as.matrix(nonzero)) >= c
filtered_dgem <- dgem[keep_genes, ]

# Further exclude mitochondrial genes
idx <- grep("^MT-", rownames(filtered_dgem)) #Change for mouse 
rownames(filtered_dgem)[idx]
filtered_dgem <- filtered_dgem[-idx, ]
dim(filtered_dgem) 

####Using filtered genes from above
file.name <- "Combined_pySCENIC_RNA_Filtered.loom"
project.title <- "GDT"
build_loom(
  file.name = file.name,
  dgem = filtered_dgem,
  title = project.title,
  genome = "human",
)

loom <- open_loom(file.path = file.name, mode = "r+")
finalize(loom = loom)
