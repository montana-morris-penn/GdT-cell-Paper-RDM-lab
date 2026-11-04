###############################################################################
# GIST γδ T-cell scRNA-seq Analysis
# Single-cell RNA sequencing processing, integration, and TCR lineage analysis
###############################################################################

#------------------------------------------------------------------------------
# 1. Load libraries
#------------------------------------------------------------------------------

library(Seurat)
library(dplyr)
library(ggplot2)
library(scales)
library(RColorBrewer)
library(cowplot)
library(scCustomize)


#------------------------------------------------------------------------------
# 2. Set directories
#------------------------------------------------------------------------------

project_dir <- "/mnt/isilon/tan_lab/sussmanj/GIST_Multiome/GD_TCells"

setwd(project_dir)

base_new <- file.path(project_dir, "scRNA_Data/New_Data_Cellranger")
base_orig <- file.path(project_dir, "scRNA_Data/Original_Data_Cellranger")
trust4_dir <- file.path(project_dir, "scRNA_Data/TRUST4_Output")

figure_dir <- file.path(project_dir, "Figures_Revision")

dir.create(figure_dir, showWarnings = FALSE, recursive = TRUE)


#------------------------------------------------------------------------------
# 3. Sample metadata
#------------------------------------------------------------------------------

meta_df <- data.frame(
  Sample = c(
    "HG612_T", "HG719_T", "HG722_T", "HG674_T",
    "HG679_T", "HG683_T", "HG703_T",
    "HG618_T", "HG618_P", "HG589_P", "HG610_P",
    "HG642_P", "HG610_T"
  ),
  
  Folder = c(
    "RonDeMatteo-20260722-GEX-HG612-GDT",
    "RonDeMatteo-20260722-GEX-HG719-GDT",
    "RonDeMatteo-20260722-GEX-HG722-GDT",
    "RonDeMatteo-20260804-GEX-HG674-GDT",
    "RonDeMatteo-20260804-GEX-HG679-GDT",
    "RonDeMatteo-20260804-GEX-HG683-GDT",
    "RonDeMatteo-20260804-GEX-HG703-GDT",
    "RonDematteo-11062023-HG618-T-GEX",
    "RonDematteo-11062023-HG618-P-GEX",
    "RonDematteo-11062023-HG589-P-GEX",
    "RonDematteo-11032023-HG610-P-GEX",
    "RonDematteo-10112023-HG642-P-GEX",
    "RonDematteo-09252023-HG610-T-YDT-GEX"
  ),
  
  BaseDir = c(
    rep(base_new, 7),
    rep(base_orig, 6)
  ),
  
  Assay = c(
    rep("V4", 7),
    rep("V3", 6)
  ),
  
  Source = c(
    rep("Tumor", 8),
    rep("PBMC", 4),
    "Tumor"
  ),
  
  Status = c(
    "Resistant", "Untreated", "Untreated",
    "Resistant", "Untreated", "Untreated",
    "Resistant", "Untreated", "Untreated",
    "Resistant", "Untreated", "Untreated",
    "Untreated"
  ),
  
  Location = c(
    "SmallBowel", "Stomach", "Stomach",
    "SmallBowel", "Stomach", "Stomach",
    "SmallBowel", "Stomach", "Stomach",
    "Peritoneum", "Stomach", "Stomach",
    "Stomach"
  ),
  
  stringsAsFactors = FALSE
)


#------------------------------------------------------------------------------
# 4. Load individual samples
#------------------------------------------------------------------------------

ob.list <- list()

for (i in seq_len(nrow(meta_df))) {
  
  sample_name <- meta_df$Sample[i]
  
  message("Loading sample: ", sample_name)
  
  counts_path <- file.path(
    meta_df$BaseDir[i],
    meta_df$Folder[i],
    "outs/filtered_feature_bc_matrix"
  )
  
  counts <- Read10X(counts_path)
  
  obj <- CreateSeuratObject(
    counts = counts,
    project = sample_name,
    assay = "RNA"
  )
  
  # Add sample metadata
  obj$Sample <- sample_name
  obj$Assay <- meta_df$Assay[i]
  obj$Source <- meta_df$Source[i]
  obj$Status <- meta_df$Status[i]
  obj$Location <- meta_df$Location[i]
  
  # Calculate mitochondrial gene expression
  obj$percent.mt <- PercentageFeatureSet(
    obj,
    pattern = "^MT-"
  )
  
  ob.list[[sample_name]] <- obj
}

#------------------------------------------------------------------------------
# 5. Initial sample-level quality control and preprocessing
#------------------------------------------------------------------------------

for (i in seq_along(ob.list)) {
  
  DefaultAssay(ob.list[[i]]) <- "RNA"
  
  # Recalculate mitochondrial percentage
  ob.list[[i]][["percent.mt"]] <- PercentageFeatureSet(
    ob.list[[i]],
    pattern = "^MT-"
  )
  
  # Filter low-quality cells
  ob.list[[i]] <- subset(
    ob.list[[i]],
    subset =
      nFeature_RNA > 300 &
      nCount_RNA > 1000 &
      nFeature_RNA < 10000 &
      percent.mt < 10
  )
  
  # Normalize and perform initial dimensionality reduction
  ob.list[[i]] <- NormalizeData(ob.list[[i]])
  ob.list[[i]] <- FindVariableFeatures(ob.list[[i]])
  ob.list[[i]] <- ScaleData(ob.list[[i]])
  ob.list[[i]] <- RunPCA(
    ob.list[[i]],
    features = VariableFeatures(ob.list[[i]])
  )
}


#------------------------------------------------------------------------------
# 6. Merge samples and perform initial analysis
#------------------------------------------------------------------------------

merged_gdt <- merge(
  ob.list[[1]],
  y = ob.list[-1]
)

merged_gdt <- JoinLayers(merged_gdt)

merged_gdt <- NormalizeData(merged_gdt)
merged_gdt <- FindVariableFeatures(merged_gdt)
merged_gdt <- ScaleData(merged_gdt)
merged_gdt <- RunPCA(merged_gdt)

# Unintegrated UMAP
merged_gdt <- RunUMAP(
  merged_gdt,
  dims = 1:30,
  reduction.name = "umap.unintegrated"
)

# Initial QC and batch assessment
DimPlot(merged_gdt, group.by = "Assay", raster = FALSE) +
  coord_fixed()

FeaturePlot(merged_gdt, features = "percent.mt", raster = FALSE) +
  coord_fixed()


#------------------------------------------------------------------------------
# 7. Cell cycle scoring and variable feature selection
#------------------------------------------------------------------------------

merged_gdt <- CellCycleScoring(
  merged_gdt,
  s.features = cc.genes$s.genes,
  g2m.features = cc.genes$g2m.genes,
  set.ident = FALSE
)

DimPlot(merged_gdt, group.by = "Phase", raster = FALSE) +
  coord_fixed()

VlnPlot(
  merged_gdt,
  features = c("percent.mt", "nCount_RNA", "nFeature_RNA"),
  group.by = "Sample",
  pt.size = 0
)

VlnPlot(
  merged_gdt,
  features = c("percent.mt", "nCount_RNA", "nFeature_RNA"),
  group.by = "Assay",
  pt.size = 0
)


# Split RNA assay by sample for integration
merged_gdt[["RNA"]] <- split(
  merged_gdt[["RNA"]],
  f = merged_gdt$Sample
)

merged_gdt <- NormalizeData(merged_gdt)

merged_gdt <- FindVariableFeatures(
  merged_gdt,
  selection.method = "vst",
  nfeatures = 2000
)

# Remove cell-cycle and proliferation-associated genes from HVGs
hvgs <- VariableFeatures(merged_gdt)

cell_cycle_genes <- unique(c(
  cc.genes$s.genes,
  cc.genes$g2m.genes
))

extra_proliferation_genes <- c(
  "MKI67", "TOP2A", "CENPF", "ASPM",
  "BUB1", "NUSAP1", "CDK1", "STMN1"
)

hvgs_clean <- setdiff(
  hvgs,
  c(cell_cycle_genes, extra_proliferation_genes)
)

VariableFeatures(merged_gdt) <- hvgs_clean

# Scale data while regressing mitochondrial gene expression
merged_gdt <- ScaleData(
  merged_gdt,
  vars.to.regress = "percent.mt"
)

merged_gdt <- RunPCA(merged_gdt)


#------------------------------------------------------------------------------
# 8. RPCA integration
#------------------------------------------------------------------------------

merged_gdt <- IntegrateLayers(
  object = merged_gdt,
  method = RPCAIntegration,
  new.reduction = "integrated.rpca",
  normalization.method = "LogNormalize"
)

merged_gdt <- RunUMAP(
  merged_gdt,
  dims = 1:30,
  reduction = "integrated.rpca",
  reduction.name = "umap.rpca"
)

# Assess integration
DimPlot(merged_gdt, group.by = "Sample",
        reduction = "umap.rpca", raster = FALSE) +
  coord_fixed()

DimPlot(merged_gdt, group.by = "Assay",
        reduction = "umap.rpca", raster = FALSE) +
  coord_fixed()

DimPlot(merged_gdt, group.by = "Phase",
        reduction = "umap.rpca", raster = FALSE) +
  coord_fixed()

DimPlot(merged_gdt, group.by = "Sample",
        reduction = "umap.unintegrated", raster = FALSE) +
  coord_fixed()


#------------------------------------------------------------------------------
# 9. Initial clustering and contaminant identification
#------------------------------------------------------------------------------

merged_gdt <- FindNeighbors(
  merged_gdt,
  dims = 1:30,
  reduction = "integrated.rpca"
)

merged_gdt <- FindClusters(
  merged_gdt,
  resolution = c(0.1, 0.3, 0.4, 0.5)
)

DimPlot(
  merged_gdt,
  group.by = "RNA_snn_res.0.5",
  reduction = "umap.rpca",
  label = TRUE,
  raster = FALSE
) + coord_fixed()


# Evaluate canonical contaminant markers
FeaturePlot(
  merged_gdt,
  features = c("KIT", "PDGFRA"),
  reduction = "umap.rpca",
  raster = FALSE
) & coord_fixed()

FeaturePlot(
  merged_gdt,
  features = c("CD79A", "CD19"),
  reduction = "umap.rpca",
  raster = FALSE
) & coord_fixed()

FeaturePlot(
  merged_gdt,
  features = c("PECAM1", "VWF"),
  reduction = "umap.rpca",
  raster = FALSE
) & coord_fixed()

FeaturePlot(
  merged_gdt,
  features = c("CD14", "CD68"),
  reduction = "umap.rpca",
  raster = FALSE
) & coord_fixed()

FeaturePlot(
  merged_gdt,
  features = c("CD4", "CD8A"),
  reduction = "umap.rpca",
  raster = FALSE
) & coord_fixed()


# Calculate lineage module scores
gd_genes <- list(c("TRDC", "TRDV1", "TRDV2", "TRGV9"))
ab_genes <- list(c("CD4", "CD8A", "CD8B", "TRBC1", "TRBC2"))
nk_genes <- list(c("NCAM1", "NKG7", "KLRD1", "GNLY", "FCGR3A"))

merged_gdt <- AddModuleScore(
  merged_gdt,
  features = gd_genes,
  name = "GD_T_Score"
)

merged_gdt <- AddModuleScore(
  merged_gdt,
  features = ab_genes,
  name = "AB_T_Score"
)

merged_gdt <- AddModuleScore(
  merged_gdt,
  features = nk_genes,
  name = "NK_Score"
)

FeaturePlot(
  merged_gdt,
  features = c("GD_T_Score1", "AB_T_Score1", "NK_Score1"),
  reduction = "umap.rpca",
  ncol = 3,
  raster = FALSE,
  min.cutoff = "q1"
) & coord_fixed()

saveRDS(
  merged_gdt,
  file = "GDT_All_Processed_Revised13_Integrated_NoCC.RDS"
)

###############################################################################
# 10. TRUST4 TCR Reconstruction
###############################################################################

# TRUST4 was performed on the original V3 samples (samples 8–13).

samples_5prime <- c(
  "RonDematteo-11062023-HG618-T-GEX",
  "RonDematteo-11062023-HG618-P-GEX",
  "RonDematteo-11062023-HG589-P-GEX",
  "RonDematteo-11032023-HG610-P-GEX",
  "RonDematteo-10112023-HG642-P-GEX",
  "RonDematteo-09252023-HG610-T-YDT-GEX"
)

tcr_list <- list()

for (i in seq_along(samples_5prime)) {
  
  sample_name <- samples_5prime[i]
  
  # Corresponding Seurat sample index
  seurat_index <- i + 7
  
  report_file <- file.path(
    trust4_dir,
    paste0(sample_name, "_trust4_barcode_report.tsv")
  )
  
  if (!file.exists(report_file)) {
    warning("TRUST4 report not found: ", report_file)
    next
  }
  
  df <- read.delim(
    report_file,
    header = TRUE,
    sep = "\t",
    comment.char = "",
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  
  # Standardize barcodes to match merged Seurat object
  df$`#barcode` <- ifelse(
    grepl("-1$", df$`#barcode`),
    paste0(df$`#barcode`, "_", seurat_index),
    paste0(df$`#barcode`, "-1_", seurat_index)
  )
  
  tcr_list[[sample_name]] <- df
}

tcr_df <- bind_rows(tcr_list)


#------------------------------------------------------------------------------
# 11. Classify TCR lineages
#------------------------------------------------------------------------------

tcr_classified <- tcr_df %>%
  group_by(`#barcode`) %>%
  summarise(
    has_ab = any(
      grepl("TRA|TRB", chain1) |
        grepl("TRA|TRB", chain2),
      na.rm = TRUE
    ),
    
    has_gd = any(
      grepl("TRG|TRD", chain1) |
        grepl("TRG|TRD", chain2),
      na.rm = TRUE
    ),
    
    .groups = "drop"
  ) %>%
  mutate(
    TCR_Lineage = case_when(
      has_ab & !has_gd ~ "ab_T_cell",
      has_gd & !has_ab ~ "gd_T_cell",
      has_ab & has_gd ~ "Dual_TCR / Ambiguous",
      TRUE ~ "Unassigned"
    )
  )


#------------------------------------------------------------------------------
# 12. Map TRUST4 classifications onto Seurat metadata
#------------------------------------------------------------------------------

merged_gdt$TRUST4_TCR_Type <- "No TCR Data (New Assay / Unassigned)"

matched_barcodes <- intersect(
  colnames(merged_gdt),
  tcr_classified$`#barcode`
)

lineage_vec <- setNames(
  tcr_classified$TCR_Lineage,
  tcr_classified$`#barcode`
)

merged_gdt$TRUST4_TCR_Type[matched_barcodes] <-
  lineage_vec[matched_barcodes]


# Restrict visualization to V3 samples with TRUST4 data
v3_only <- subset(
  merged_gdt,
  subset = Assay == "V3"
)


#------------------------------------------------------------------------------
# 13. TRUST4 lineage visualization
#------------------------------------------------------------------------------

tcr_colors <- c(
  "ab_T_cell" = "firebrick2",
  "gd_T_cell" = "forestgreen",
  "Dual_TCR / Ambiguous" = "purple",
  "No TCR Data (New Assay / Unassigned)" = "grey85"
)

p1 <- DimPlot(
  v3_only,
  group.by = "TRUST4_TCR_Type",
  cols = tcr_colors,
  raster = TRUE,
  raster.dpi = c(3000, 3000),
  pt.size = 6.5,
  order = FALSE,
  reduction = "umap.rpca"
) +
  ggtitle("TRUST4 Lineage Overlay (Samples 8–13)") +
  coord_fixed()

p1

ggsave(
  filename = file.path(figure_dir, "TRUST4_Overlay_UMAP.pdf"),
  plot = p1,
  width = 15,
  height = 10
)


# Absolute cell counts per cluster
p_counts <- ggplot(
  v3_only@meta.data,
  aes(
    x = factor(RNA_snn_res.0.4),
    fill = TRUST4_TCR_Type
  )
) +
  geom_bar(
    position = "stack",
    color = "black",
    linewidth = 0.2
  ) +
  scale_fill_manual(values = tcr_colors) +
  labs(
    x = "Seurat Cluster",
    y = "Cell Count",
    fill = "TCR Lineage",
    title = "TCR Lineage Cell Counts per Cluster"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(
      size = 14,
      face = "bold",
      hjust = 0.5
    ),
    axis.text = element_text(size = 10, color = "black"),
    axis.title = element_text(size = 12, face = "bold"),
    legend.title = element_text(face = "bold")
  )

p_counts


# Relative lineage composition per cluster
p_proportions <- ggplot(
  v3_only@meta.data,
  aes(
    x = factor(RNA_snn_res.0.4),
    fill = TRUST4_TCR_Type
  )
) +
  geom_bar(
    position = "fill",
    color = "black",
    linewidth = 0.2
  ) +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(values = tcr_colors) +
  labs(
    x = "Seurat Cluster",
    y = "Proportion of Cells",
    fill = "TCR Lineage",
    title = "TCR Lineage Composition per Cluster"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(
      size = 14,
      face = "bold",
      hjust = 0.5
    ),
    axis.text = element_text(size = 10, color = "black"),
    axis.title = element_text(size = 12, face = "bold"),
    legend.title = element_text(face = "bold")
  )

p_proportions

###############################################################################
# 14. Remove non-T-cell populations
###############################################################################

# Preserve the full integrated object
merged_gdt_full <- merged_gdt

# Remove clusters previously identified as non-T-cell contaminants
merged_gdt <- subset(
  merged_gdt_full,
  subset = RNA_snn_res.0.5 %in% c(14, 16, 15, 13, 8),
  invert = TRUE
)

DimPlot(
  merged_gdt,
  group.by = "RNA_snn_res.0.5",
  reduction = "umap.rpca",
  label = TRUE,
  raster = FALSE
) + coord_fixed()


#------------------------------------------------------------------------------
# 15. Reprocess and reintegrate the T-cell-enriched population
#------------------------------------------------------------------------------

merged_gdt <- JoinLayers(merged_gdt)

merged_gdt <- CellCycleScoring(
  merged_gdt,
  s.features = cc.genes$s.genes,
  g2m.features = cc.genes$g2m.genes,
  set.ident = FALSE
)

merged_gdt[["RNA"]] <- split(
  merged_gdt[["RNA"]],
  f = merged_gdt$Sample
)

merged_gdt <- NormalizeData(merged_gdt)
merged_gdt <- FindVariableFeatures(merged_gdt)

# Exclude cell-cycle-associated genes
hvgs <- VariableFeatures(merged_gdt)

cell_cycle_genes <- unique(c(
  cc.genes$s.genes,
  cc.genes$g2m.genes
))

extra_proliferation_genes <- c(
  "MKI67", "TOP2A", "CENPF", "ASPM",
  "BUB1", "NUSAP1", "CDK1", "STMN1"
)

hvgs_clean <- setdiff(
  hvgs,
  c(cell_cycle_genes, extra_proliferation_genes)
)

VariableFeatures(merged_gdt) <- hvgs_clean

merged_gdt <- ScaleData(
  merged_gdt,
  vars.to.regress = "percent.mt"
)

merged_gdt <- RunPCA(merged_gdt)

# RPCA integration
merged_gdt <- IntegrateLayers(
  object = merged_gdt,
  method = RPCAIntegration,
  new.reduction = "integrated.rpca",
  normalization.method = "LogNormalize"
)

# Dimensionality reduction and clustering
merged_gdt <- RunUMAP(
  merged_gdt,
  dims = 1:30,
  reduction = "integrated.rpca",
  reduction.name = "umap.rpca"
)

merged_gdt <- FindNeighbors(
  merged_gdt,
  dims = 1:30,
  reduction = "integrated.rpca"
)

merged_gdt <- FindClusters(
  merged_gdt,
  resolution = c(0.1, 0.2, 0.3, 0.4, 0.5)
)


# Assess clustering and sample composition
DimPlot(merged_gdt, group.by = "Phase",
        reduction = "umap.rpca", raster = FALSE) +
  coord_fixed()

DimPlot(merged_gdt, group.by = "Assay",
        reduction = "umap.rpca", raster = FALSE) +
  coord_fixed()

DimPlot(merged_gdt, group.by = "Source",
        reduction = "umap.rpca", raster = FALSE) +
  coord_fixed()

DimPlot(
  merged_gdt,
  group.by = "RNA_snn_res.0.5",
  reduction = "umap.rpca",
  label = TRUE,
  raster = FALSE
) + coord_fixed()


###############################################################################
# 16. Identify and remove putative alpha-beta T-cell contaminants
###############################################################################

all_genes <- rownames(merged_gdt)

# Define lineage-associated marker sets
gd_genes <- list(intersect(
  c("TRDC", "TRDV1", "TRDV2", "TRGV9", "TRGC1", "TRGC2"),
  all_genes
))

ab_genes <- list(intersect(
  c("TRAC", "TRBC1", "TRBC2", "CD4", "CD8A", "CD8B"),
  all_genes
))

nk_genes <- list(intersect(
  c("NCAM1", "FCGR3A", "KLRD1", "NKG7", "GNLY"),
  all_genes
))

myeloid_genes <- list(intersect(
  c("CD14", "LYZ", "CD68", "AIF1", "CSF1R"),
  all_genes
))

b_genes <- list(intersect(
  c("MS4A1", "CD79A", "CD79B", "CD19"),
  all_genes
))

# Marker panel for visualization
special_genes <- c(
  "CD3E",
  gd_genes[[1]],
  ab_genes[[1]],
  nk_genes[[1]],
  myeloid_genes[[1]],
  b_genes[[1]]
)


# Calculate module scores
merged_gdt <- AddModuleScore(
  merged_gdt,
  features = gd_genes,
  name = "Score_gdT"
)

merged_gdt <- AddModuleScore(
  merged_gdt,
  features = ab_genes,
  name = "Score_abT"
)

merged_gdt <- AddModuleScore(
  merged_gdt,
  features = nk_genes,
  name = "Score_NK"
)

merged_gdt <- AddModuleScore(
  merged_gdt,
  features = myeloid_genes,
  name = "Score_Myeloid"
)

merged_gdt <- AddModuleScore(
  merged_gdt,
  features = b_genes,
  name = "Score_Bcell"
)


# Rename module score metadata columns
merged_gdt$Score_gdT <- merged_gdt$Score_gdT1
merged_gdt$Score_abT <- merged_gdt$Score_abT1
merged_gdt$Score_NK <- merged_gdt$Score_NK1
merged_gdt$Score_Myeloid <- merged_gdt$Score_Myeloid1
merged_gdt$Score_Bcell <- merged_gdt$Score_Bcell1

# Calculate relative alpha-beta versus gamma-delta T-cell score
merged_gdt$abT_vs_gdT_Index <-
  merged_gdt$Score_abT - merged_gdt$Score_gdT


# Define putative alpha-beta T-cell contaminants
ab_cutoff <- quantile(merged_gdt$Score_abT, 0.90)
gd_cutoff <- quantile(merged_gdt$Score_gdT, 0.25)

merged_gdt$contaminant_flag <- ifelse(
  merged_gdt$Score_abT > ab_cutoff &
    merged_gdt$Score_gdT < gd_cutoff,
  "Putative_abT_Contaminant",
  "Pass"
)


#------------------------------------------------------------------------------
# 17. Evaluate lineage module scores
#------------------------------------------------------------------------------

FeaturePlot(
  merged_gdt,
  features = c(
    "Score_gdT",
    "Score_abT",
    "Score_NK",
    "abT_vs_gdT_Index"
  ),
  ncol = 2,
  reduction = "umap.rpca",
  raster = FALSE,
  cols = c("grey85", "red")
) & coord_fixed()


VlnPlot(
  merged_gdt,
  features = c(
    "Score_gdT",
    "Score_abT",
    "Score_NK",
    "Score_Myeloid",
    "Score_Bcell",
    "abT_vs_gdT_Index"
  ),
  group.by = "RNA_snn_res.0.4",
  pt.size = 0,
  ncol = 3
)


VlnPlot(
  merged_gdt,
  features = special_genes,
  group.by = "RNA_snn_res.0.4",
  pt.size = 0,
  ncol = 3
)


DotPlot(
  merged_gdt,
  features = special_genes,
  group.by = "RNA_snn_res.0.4"
) +
  RotatedAxis()


DimPlot(
  merged_gdt,
  group.by = "contaminant_flag",
  reduction = "umap.rpca",
  raster = FALSE,
  order = TRUE,
  cols = c(
    "Pass" = "grey80",
    "Putative_abT_Contaminant" = "firebrick"
  )
) & coord_fixed()


# Summarize putative contaminants by cluster
contam_summary <- merged_gdt@meta.data %>%
  group_by(RNA_snn_res.0.4, contaminant_flag) %>%
  summarise(
    cell_count = n(),
    .groups = "drop"
  ) %>%
  group_by(RNA_snn_res.0.4) %>%
  mutate(
    percent_of_cluster = cell_count / sum(cell_count) * 100
  )

print(
  contam_summary %>%
    filter(contaminant_flag == "Putative_abT_Contaminant")
)

###############################################################################
# 18. Remove putative contaminating clusters
###############################################################################

# Preserve the T-cell-enriched object before final filtering
merged_gdt_Tall <- merged_gdt

# Exclude clusters 4, 10, and 12 at resolution 0.4
merged_gdt <- subset(
  merged_gdt_Tall,
  subset = RNA_snn_res.0.4 %in% c(4, 10, 12),
  invert = TRUE
)


###############################################################################
# 19. Reintegrate the γδ T-cell-enriched population
###############################################################################

merged_gdt <- JoinLayers(merged_gdt)

merged_gdt <- CellCycleScoring(
  merged_gdt,
  s.features = cc.genes$s.genes,
  g2m.features = cc.genes$g2m.genes,
  set.ident = FALSE
)

merged_gdt[["RNA"]] <- split(
  merged_gdt[["RNA"]],
  f = merged_gdt$Sample
)

merged_gdt <- NormalizeData(merged_gdt)
merged_gdt <- FindVariableFeatures(merged_gdt)

merged_gdt <- ScaleData(
  merged_gdt,
  vars.to.regress = "percent.mt"
)

merged_gdt <- RunPCA(merged_gdt)

merged_gdt <- IntegrateLayers(
  object = merged_gdt,
  method = RPCAIntegration,
  new.reduction = "integrated.rpca",
  normalization.method = "LogNormalize"
)

merged_gdt <- RunUMAP(
  merged_gdt,
  dims = 1:30,
  reduction = "integrated.rpca",
  reduction.name = "umap.rpca"
)

merged_gdt <- FindNeighbors(
  merged_gdt,
  dims = 1:30,
  reduction = "integrated.rpca"
)

merged_gdt <- FindClusters(
  merged_gdt,
  resolution = c(0.1, 0.2, 0.3, 0.4, 0.5)
)


# Assess final clustering
DimPlot(merged_gdt, group.by = "Phase",
        reduction = "umap.rpca", raster = FALSE) +
  coord_fixed()

DimPlot(merged_gdt, group.by = "Assay",
        reduction = "umap.rpca", raster = FALSE) +
  coord_fixed()

DimPlot(merged_gdt, group.by = "Source",
        reduction = "umap.rpca", raster = FALSE) +
  coord_fixed()

DimPlot(
  merged_gdt,
  group.by = "RNA_snn_res.0.5",
  reduction = "umap.rpca",
  label = TRUE,
  raster = FALSE
) + coord_fixed()


###############################################################################
# 20. Cluster marker identification
###############################################################################

Idents(merged_gdt) <- "RNA_snn_res.0.5"

markers2 <- FindAllMarkers(merged_gdt)

markers2$pct.diff <- markers2$pct.1 - markers2$pct.2

write.table(
  markers2,
  file = "GDT_14Cluster_Markers.tsv",
  sep = "\t",
  quote = FALSE
)


###############################################################################
# 21. Cluster annotation
###############################################################################

annotation_genes <- c(
  # Gamma-delta T-cell markers
  "TRDC", "TRGC1", "TRGC2",
  "TRDV1", "TRDV2", "TRDV3",
  "TRGV9", "TRGV2", "TRGV4",
  
  # Alpha-beta T-cell markers
  "TRAC", "TRBC1", "TRBC2",
  "CD3D", "CD3E", "CD3G",
  "CD4", "CD8A", "CD8B",
  
  # NK-associated markers
  "GNLY", "NKG7", "KLRD1",
  "GZMB", "GZMK", "FCGR3A", "NCAM1",
  
  # Naive/memory-associated markers
  "CCR7", "IL7R"
)


VlnPlot(
  merged_gdt,
  features = annotation_genes,
  group.by = "RNA_snn_res.0.4",
  pt.size = 0,
  ncol = 3
)


DotPlot(
  merged_gdt,
  features = annotation_genes,
  group.by = "RNA_snn_res.0.4"
) +
  RotatedAxis()


FeaturePlot(
  merged_gdt,
  features = c("CD3E", "TRDC"),
  reduction = "umap.rpca",
  raster = FALSE
) & coord_fixed()


FeaturePlot(
  merged_gdt,
  features = c("nCount_RNA", "nFeature_RNA"),
  reduction = "umap.rpca",
  raster = FALSE
) & coord_fixed()


FeaturePlot(
  merged_gdt,
  features = "percent.mt",
  reduction = "umap.rpca",
  raster = FALSE
) & coord_fixed()


DimPlot(
  merged_gdt,
  group.by = "Source",
  reduction = "umap.rpca",
  label = TRUE,
  raster = FALSE
) + coord_fixed()


###############################################################################
# 22. Additional quality control
###############################################################################

merged_gdt[["percent.ribo"]] <- PercentageFeatureSet(
  merged_gdt,
  pattern = "^RP[SL]"
)

VlnPlot(
  merged_gdt,
  features = c(
    "percent.mt",
    "percent.ribo",
    "nCount_RNA",
    "nFeature_RNA"
  ),
  group.by = "RNA_snn_res.0.5",
  pt.size = 0,
  ncol = 2
)


###############################################################################
# 23. Filter low-quality cells and reintegrate
###############################################################################

# Retain cells with adequate feature complexity and ribosomal content
merged_gdt <- subset(
  merged_gdt,
  subset = nFeature_RNA > 500 & percent.ribo > 5
)


DimPlot(
  merged_gdt,
  group.by = "RNA_snn_res.0.4",
  reduction = "umap.rpca",
  label = TRUE,
  raster = FALSE
) + coord_fixed()


VlnPlot(
  merged_gdt,
  features = annotation_genes,
  group.by = "RNA_snn_res.0.4",
  pt.size = 0,
  ncol = 3
)


# Reprocess and reintegrate following QC
merged_gdt <- JoinLayers(merged_gdt)

merged_gdt <- CellCycleScoring(
  merged_gdt,
  s.features = cc.genes$s.genes,
  g2m.features = cc.genes$g2m.genes,
  set.ident = FALSE
)

merged_gdt[["RNA"]] <- split(
  merged_gdt[["RNA"]],
  f = merged_gdt$Sample
)

merged_gdt <- NormalizeData(merged_gdt)
merged_gdt <- FindVariableFeatures(merged_gdt)

merged_gdt <- ScaleData(
  merged_gdt,
  vars.to.regress = "percent.mt"
)

merged_gdt <- RunPCA(merged_gdt)

merged_gdt <- IntegrateLayers(
  object = merged_gdt,
  method = RPCAIntegration,
  new.reduction = "integrated.rpca",
  normalization.method = "LogNormalize"
)

merged_gdt <- RunUMAP(
  merged_gdt,
  dims = 1:30,
  reduction = "integrated.rpca",
  reduction.name = "umap.rpca"
)

merged_gdt <- FindNeighbors(
  merged_gdt,
  dims = 1:30,
  reduction = "integrated.rpca"
)

merged_gdt <- FindClusters(
  merged_gdt,
  resolution = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7)
)


# Assess final QC and clustering
DimPlot(merged_gdt, group.by = "Phase",
        reduction = "umap.rpca", raster = FALSE) +
  coord_fixed()

DimPlot(merged_gdt, group.by = "Assay",
        reduction = "umap.rpca", raster = FALSE) +
  coord_fixed()

DimPlot(merged_gdt, group.by = "Source",
        reduction = "umap.rpca", raster = FALSE) +
  coord_fixed()

DimPlot(
  merged_gdt,
  group.by = "RNA_snn_res.0.6",
  reduction = "umap.rpca",
  label = TRUE,
  label.box = TRUE,
  raster = FALSE
) + coord_fixed()


# Cluster sizes
table(merged_gdt$RNA_snn_res.0.6)


# Assess TRDV1 and TRDV2 expression
FeatureScatter(
  merged_gdt,
  feature1 = "TRDV1",
  feature2 = "TRDV2"
)

FeaturePlot(
  merged_gdt,
  features = c("TRDV1", "TRDV2"),
  blend = TRUE,
  raster = FALSE,
  reduction = "umap.rpca",
  ncol = 2,
  order = TRUE,
  blend.threshold = 0.5,
  cols = c("grey", "red", "green")
) + coord_fixed()


# Explore higher clustering resolutions
merged_gdt <- FindClusters(
  merged_gdt,
  resolution = seq(0.8, 1.2, by = 0.1)
)

merged_gdt <- FindClusters(
  merged_gdt,
  resolution = seq(1.3, 1.7, by = 0.1)
)

merged_gdt <- JoinLayers(merged_gdt)


# Final higher-resolution clustering visualization
DimPlot(
  merged_gdt,
  group.by = "RNA_snn_res.1.6",
  reduction = "umap.rpca",
  raster = FALSE,
  label = TRUE,
  label.box = TRUE
) + coord_fixed()


VlnPlot(
  merged_gdt,
  features = annotation_genes,
  group.by = "RNA_snn_res.1.6",
  pt.size = 0,
  ncol = 3
)


DotPlot(
  merged_gdt,
  features = annotation_genes,
  group.by = "RNA_snn_res.0.7"
) +
  scale_colour_gradientn(
    colours = rev(RColorBrewer::brewer.pal(11, "RdBu"))
  ) +
  RotatedAxis()

###############################################################################
# 24. Figure generation
###############################################################################

#------------------------------------------------------------------------------
# 24.1. Initial dataset: contaminant marker visualization
#------------------------------------------------------------------------------

seurat1 <- readRDS(
  "GDT_All_Processed_Revised13_Integrated_NoCC.RDS"
)

# Helper function for standardized FeaturePlots
make_feature_plot <- function(
    object,
    feature,
    reduction = "umap.rpca",
    pt.size = 2.5,
    raster.dpi = c(800, 800),
    max.cutoff = "q99") {
  
  FeaturePlot(
    object,
    features = feature,
    reduction = reduction,
    raster = TRUE,
    raster.dpi = raster.dpi,
    pt.size = pt.size,
    max.cutoff = max.cutoff
  ) +
    coord_fixed() +
    theme(
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank()
    )
}


# Generate individual marker plots
contaminant_markers <- c(
  "KIT", "PDGFRA",
  "CD79A", "CD19",
  "CD14", "CD68",
  "PECAM1", "VWF"
)

contaminant_plots <- lapply(
  contaminant_markers,
  function(gene) make_feature_plot(seurat1, gene)
)

# Combine into a single figure
p_combined <- cowplot::plot_grid(
  plotlist = contaminant_plots,
  nrow = 2,
  ncol = 4
)

ggsave(
  filename = file.path(
    figure_dir,
    "FeaturePlot_Grid_Contaminant.pdf"
  ),
  plot = p_combined,
  device = "pdf",
  width = 16,
  height = 8
)


#------------------------------------------------------------------------------
# 24.2. Initial dataset: UMAP
#------------------------------------------------------------------------------

p1 <- DimPlot(
  seurat1,
  reduction = "umap.rpca",
  group.by = "RNA_snn_res.0.5",
  raster = TRUE,
  raster.dpi = c(3000, 3000),
  pt.size = 4.5,
  label = TRUE
) +
  coord_fixed()

ggsave(
  filename = file.path(figure_dir, "UMAP_All_1.pdf"),
  plot = p1,
  width = 10,
  height = 10
)


#------------------------------------------------------------------------------
# 24.3. T-cell-enriched dataset: UMAP and marker visualization
#------------------------------------------------------------------------------

seurat2 <- readRDS(
  "GDT_All_Processed_Revised13_TOnly_NoCC_Integrated.RDS"
)

p1 <- DimPlot(
  seurat2,
  reduction = "umap.rpca",
  group.by = "RNA_snn_res.0.4",
  raster = TRUE,
  raster.dpi = c(3000, 3000),
  pt.size = 4.5,
  label = TRUE
) +
  coord_fixed()

ggsave(
  filename = file.path(figure_dir, "UMAP_All_2.pdf"),
  plot = p1,
  width = 10,
  height = 10
)


# Order clusters for DotPlot visualization
seurat2$RNA_snn_res.0.4 <- factor(
  seurat2$RNA_snn_res.0.4,
  levels = rev(0:12)
)

DefaultAssay(seurat2) <- "RNA"

pdot <- DotPlot(
  seurat2,
  features = special_genes,
  group.by = "RNA_snn_res.0.4",
  assay = "RNA"
) +
  scale_colour_gradientn(
    colours = rev(RColorBrewer::brewer.pal(11, "RdBu"))
  ) +
  RotatedAxis()

ggsave(
  filename = file.path(
    figure_dir,
    "Dotplot_Remove_Contaminant_T.pdf"
  ),
  plot = pdot,
  width = 7.75,
  height = 4
)


#------------------------------------------------------------------------------
# 24.4. γδ T-cell-enriched dataset: UMAP
#------------------------------------------------------------------------------

seurat3 <- readRDS(
  "GDT_All_Processed_Revised13_GDT_Integrated.RDS"
)

p1 <- DimPlot(
  seurat3,
  reduction = "umap.rpca",
  group.by = "RNA_snn_res.0.4",
  raster = TRUE,
  raster.dpi = c(3000, 3000),
  pt.size = 4.5,
  label = TRUE
) +
  coord_fixed()

ggsave(
  filename = file.path(figure_dir, "UMAP_All_3.pdf"),
  plot = p1,
  width = 10,
  height = 10
)


#------------------------------------------------------------------------------
# 24.5. γδ T-cell dataset: quality control violin plots
#------------------------------------------------------------------------------

seurat3[["percent.ribo"]] <- PercentageFeatureSet(
  seurat3,
  pattern = "^RP[SL]"
)

n_clusters <- length(levels(Idents(seurat3)))

default_seurat_cols <- scales::hue_pal()(n_clusters)

p1 <- scCustomize::Stacked_VlnPlot(
  seurat_object = seurat3,
  group.by = "RNA_snn_res.0.4",
  colors_use = default_seurat_cols,
  features = c(
    "CD3E",
    "TRDC",
    "percent.mt",
    "nCount_RNA",
    "percent.ribo"
  )
)

ggsave(
  filename = file.path(
    figure_dir,
    "ViolinQuality_First.pdf"
  ),
  plot = p1,
  width = 6,
  height = 4.5
)


#------------------------------------------------------------------------------
# 24.6. Sample-level quality control
#------------------------------------------------------------------------------

# NOTE: gdt_seurat must be loaded or defined before running this section.
# Sample_Colors.rds contains the manually defined sample color palette.

colors4 <- readRDS("Figures/Sample_Colors.rds")

p1 <- scCustomize::Stacked_VlnPlot(
  seurat_object = gdt_seurat,
  group.by = "orig.ident",
  colors_use = colors4,
  features = c(
    "percent.mt",
    "nCount_RNA",
    "nFeature_RNA",
    "percent.ribo"
  ),
  x_lab_rotate = 45
) +
  RotatedAxis()

ggsave(
  filename = file.path(
    figure_dir,
    "ViolinQuality_Second.pdf"
  ),
  plot = p1,
  width = 7,
  height = 4.5
)

#------------------------------------------------------------------------------
# 24.7. TRDV expression summary and cluster composition
#------------------------------------------------------------------------------
# 11. TRDV expression summary and cluster composition
###############################################################################
# TRDV1/2/3 expression by cluster
trdv_colors <- c(
  "TRDV1_only" = "#377EB8",
  "TRDV2_only" = "#E41A1C",
  "TRDV3_only" = "#984EA3",
  "Multiple" = "#4DAF4A",
  "None" = "#D3D3D3"
)

plot_data <- FetchData(
  merged_gdt,
  vars = c("TRDV1", "TRDV2", "TRDV3", "RNA_snn_res.0.8")
) %>%
  mutate(
    n_expressed = (TRDV1 > 0) + (TRDV2 > 0) + (TRDV3 > 0),
    trdv_status = case_when(
      n_expressed == 0 ~ "None",
      n_expressed > 1 ~ "Multiple",
      TRDV1 > 0 ~ "TRDV1_only",
      TRDV2 > 0 ~ "TRDV2_only",
      TRUE ~ "TRDV3_only"
    ),
    trdv_status = factor(trdv_status, levels = names(trdv_colors))
  ) %>%
  count(RNA_snn_res.0.8, trdv_status) %>%
  group_by(RNA_snn_res.0.8) %>%
  mutate(Percentage = 100 * n / sum(n))

p <- ggplot(
  plot_data,
  aes(x = RNA_snn_res.0.8, y = Percentage, fill = trdv_status)
) +
  geom_col(width = 0.7, color = "black", linewidth = 0.2) +
  scale_fill_manual(values = trdv_colors, name = "TRDV Status") +
  scale_y_continuous(
    labels = function(x) paste0(x, "%"),
    expand = c(0, 0),
    limits = c(0, 100.5)
  ) +
  labs(
    title = "TRDV Chain Expression by Cluster",
    x = "Cell Type",
    y = "Percentage of Cells"
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

print(p)

ggsave(
  "Figures_Revision/TRDV12_barplot_clusters_original.pdf",
  p,
  width = 5,
  height = 5
)
