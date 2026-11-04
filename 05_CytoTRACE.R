###############################################################################
# CytoTRACE analysis of paired PBMC and Tumor samples
###############################################################################

library(Seurat)
library(CytoTRACE)
library(dplyr)
library(rstatix)

setwd("/mnt/isilon/tan_lab/sussmanj/GIST_Multiome/GD_TCells")

dir.create("Figures_Revision", showWarnings = FALSE, recursive = TRUE)


###############################################################################
# 1. Load Seurat object
###############################################################################

gdt_seurat <- readRDS(
  "GDT_All_Processed_Revised13_GDT_QC_Integrated2_SCENIC.RDS"
)

DefaultAssay(gdt_seurat) <- "RNA"
gdt_seurat <- JoinLayers(gdt_seurat)


###############################################################################
# 2. Calculate CytoTRACE separately for each patient
###############################################################################

calculate_cytotrace <- function(object, samples) {
  
  subset_obj <- subset(object, subset = orig.ident %in% samples)
  
  expression_matrix <- GetAssayData(
    subset_obj,
    assay = "RNA",
    layer = "data"
  )
  
  cytotrace_results <- CytoTRACE(
    as.data.frame(expression_matrix),
    ncores = 1
  )
  
  scores <- cytotrace_results$CytoTRACE[colnames(subset_obj)]
  
  subset_obj$CytoTRACE_score <- scores
  
  return(subset_obj)
}

HG610subset <- calculate_cytotrace(
  gdt_seurat,
  c("HG610_P", "HG610_T")
)

HG618subset <- calculate_cytotrace(
  gdt_seurat,
  c("HG618_P", "HG618_T")
)


###############################################################################
# 3. Add CytoTRACE scores to original Seurat object
###############################################################################

cytotrace_scores <- c(
  HG610subset$CytoTRACE_score,
  HG618subset$CytoTRACE_score
)

names(cytotrace_scores) <- c(
  colnames(HG610subset),
  colnames(HG618subset)
)

gdt_seurat <- AddMetaData(
  gdt_seurat,
  metadata = cytotrace_scores,
  col.name = "CytoTRACE_score"
)


###############################################################################
# 4. Violin plots: CytoTRACE by broad class and Source
###############################################################################

broad_class_order <- c("VD1+", "VD2+", "VD3+", "Mixed")

HG618subset$broad_class <- factor(
  HG618subset$broad_class,
  levels = broad_class_order
)

p_HG618 <- VlnPlot(
  HG618subset,
  features = "CytoTRACE_score",
  group.by = "broad_class",
  split.by = "Source",
  raster = FALSE,
  pt.size = 0
)

p_HG618

ggsave(
  p_HG618,
  filename = "Figures_Revision/Violin_CytoTRACE_618.pdf",
  width = 4,
  height = 3
)


HG610subset$broad_class <- factor(
  HG610subset$broad_class,
  levels = broad_class_order
)

p_HG610 <- VlnPlot(
  HG610subset,
  features = "CytoTRACE_score",
  group.by = "broad_class",
  split.by = "Source",
  raster = FALSE,
  pt.size = 0
)

p_HG610

ggsave(
  p_HG610,
  filename = "Figures_Revision/Violin_CytoTRACE_610.pdf",
  width = 4,
  height = 3
)


###############################################################################
# 5. Statistical comparisons: PBMC vs Tumor within each broad class
###############################################################################

cytotrace_stats <- function(object) {
  object@meta.data %>%
    group_by(broad_class) %>%
    wilcox_test(CytoTRACE_score ~ Source) %>%
    adjust_pvalue(method = "BH") %>%
    add_significance()
}

stats_HG618 <- cytotrace_stats(HG618subset)
stats_HG610 <- cytotrace_stats(HG610subset)

stats_HG618
stats_HG610


###############################################################################
# 6. Save annotated Seurat object
###############################################################################

saveRDS(
  gdt_seurat,
  "GDT_All_Processed_Revised13_GDT_QC_Integrated2_SCENIC_CytoTRACE.RDS"
)

