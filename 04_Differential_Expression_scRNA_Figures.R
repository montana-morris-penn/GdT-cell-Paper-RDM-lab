###############################################################################
# GIST gamma-delta T-cell analysis
###############################################################################

library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggrepel)
library(cowplot)
library(RColorBrewer)
library(scales)
library(speckle)
library(writexl)

setwd("/mnt/isilon/tan_lab/sussmanj/GIST_Multiome/GD_T")

gdt_seurat <- readRDS("GDT_All_Processed_Revised13_GDT_QC_Integrated2.RDS")

dir.create("Figures_Revision", showWarnings = FALSE, recursive = TRUE)
dir.create("Tables", showWarnings = FALSE, recursive = TRUE)


###############################################################################
# 1. Colors
###############################################################################

# Random colors, with a fixed seed to keep colors consistent across runs.
set.seed(123)

random_palette <- function(n) {
  sample(grDevices::hcl.colors(max(n, 30), palette = "Dynamic"), n)
}


###############################################################################
# 2. Cell-type and broad-class annotations
###############################################################################

cell_type_map <- c(
  "2" = "VD1+ CTL response", "3" = "VD1+ CTL response",
  "4" = "VD1+ Naive",
  "10" = "VD1+ Inflammatory", "11" = "VD1+ Inflammatory",
  "5" = "VD1+ MHC-II+", "18" = "VD1+ Interferon+",
  "0" = "VD2+ innate-like", "1" = "VD2+ innate-like",
  "9" = "VD2+ Type 3 Immunity", "12" = "VD2+ Type 3 Immunity",
  "6" = "VD2+ CTL response", "15" = "VD2+ MHC-II positive",
  "17" = "VD2+ Interferon+", "7" = "VD2+ Inflammatory",
  "13" = "VD3+", "8" = "Mixed",
  "14" = "Proliferating", "16" = "Proliferating"
)

broad_class_map <- c(
  "2" = "VD1+", "3" = "VD1+", "4" = "VD1+",
  "10" = "VD1+", "11" = "VD1+", "5" = "VD1+", "18" = "VD1+",
  "0" = "VD2+", "1" = "VD2+", "9" = "VD2+", "12" = "VD2+",
  "6" = "VD2+", "15" = "VD2+", "17" = "VD2+", "7" = "VD2+",
  "13" = "VD3+", "8" = "Mixed", "14" = "Mixed", "16" = "Mixed"
)

cluster_ids <- as.character(gdt_seurat$RNA_snn_res.0.8)

gdt_seurat$cell_type <- unname(cell_type_map[cluster_ids])
gdt_seurat$broad_class <- unname(broad_class_map[cluster_ids])

if (anyNA(gdt_seurat$cell_type) || anyNA(gdt_seurat$broad_class)) {
  stop("Some clusters are missing annotations.")
}

cell_type_order <- c(
  "VD1+ CTL response", "VD1+ Naive", "VD1+ Inflammatory",
  "VD1+ MHC-II+", "VD1+ Interferon+", "VD2+ innate-like",
  "VD2+ Type 3 Immunity", "VD2+ CTL response", "VD2+ Inflammatory",
  "VD2+ Interferon+", "VD2+ MHC-II positive", "VD3+", "Mixed",
  "Proliferating"
)

broad_class_order <- c("VD1+", "VD2+", "VD3+", "Mixed")
source_order <- c("PBMC", "Tumor")

gdt_seurat$cell_type <- factor(gdt_seurat$cell_type, levels = cell_type_order)
gdt_seurat$broad_class <- factor(gdt_seurat$broad_class, levels = broad_class_order)

cell_type_cols <- setNames(random_palette(length(cell_type_order)), cell_type_order)
broad_class_cols <- setNames(random_palette(length(broad_class_order)), broad_class_order)
source_cols <- setNames(random_palette(length(source_order)), source_order)

sample_order <- sort(unique(as.character(gdt_seurat$orig.ident)))
sample_cols <- setNames(random_palette(length(sample_order)), sample_order)


###############################################################################
# 3. UMAP plots
###############################################################################

# Generate UMAPs using the same settings and save each plot.
make_umap <- function(object, group, colors = NULL, label = FALSE, filename,
                      width = 8, height = 7) {
  p <- DimPlot(
    object, reduction = "umap.rpca", group.by = group, cols = colors,
    raster = TRUE, raster.dpi = c(3000, 3000), pt.size = 4.5,
    label = label, repel = TRUE
  ) + coord_fixed() + theme_classic()
  
  ggsave(filename, p, width = width, height = height)
  return(p)
}

p_umap_clusters <- make_umap(
  gdt_seurat, "RNA_snn_res.0.8",
  filename = "Figures_Revision/UMAP_Clusters.pdf"
)

p_umap_celltype <- make_umap(
  gdt_seurat, "cell_type", cell_type_cols, TRUE,
  "Figures_Revision/UMAP_Cell_Types.pdf", 10, 8
)

p_umap_source <- make_umap(
  gdt_seurat, "Source", source_cols,
  filename = "Figures_Revision/UMAP_Source.pdf"
)

p_umap_samples <- make_umap(
  gdt_seurat, "orig.ident", sample_cols,
  filename = "Figures_Revision/UMAP_Samples.pdf", width = 11, height = 8
)

p_umap_broad <- make_umap(
  gdt_seurat, "broad_class", broad_class_cols, TRUE,
  "Figures_Revision/UMAP_Broad_Class.pdf"
)

if ("Assay" %in% colnames(gdt_seurat@meta.data)) {
  p_umap_assay <- make_umap(
    gdt_seurat, "Assay",
    filename = "Figures_Revision/UMAP_Assay.pdf"
  )
}


###############################################################################
# 4. Composition barplots
###############################################################################

barplots <- function(data, x, y, colors, as_percent = TRUE,
                     filename = NULL, width = 8, height = 6) {
  data[[y]] <- factor(data[[y]], levels = names(colors))
  data[[x]] <- factor(data[[x]], levels = unique(data[[x]]))
  
  p <- ggplot(data, aes(x = .data[[x]], y = n, fill = .data[[y]])) +
    geom_bar(
      stat = "identity",
      position = if (as_percent) "fill" else "stack",
      color = "black", linewidth = 0.2
    ) +
    scale_fill_manual(values = colors, drop = FALSE) +
    theme_classic() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.title = element_blank()
    ) +
    labs(x = NULL, y = if (as_percent) "Proportion" else "Cell count")
  
  if (as_percent) p <- p + scale_y_continuous(labels = percent_format())
  if (!is.null(filename)) ggsave(filename, p, width = width, height = height)
  p
}

meta <- gdt_seurat@meta.data %>% filter(Source %in% source_order)

broad_counts <- meta %>% count(Source, broad_class, name = "n")
celltype_counts <- meta %>% count(Source, cell_type, name = "n")
sample_celltype_counts <- meta %>% count(orig.ident, cell_type, name = "n")

p_broad_composition <- barplots(
  broad_counts, "Source", "broad_class", broad_class_cols,
  filename = "Figures_Revision/Broad_Class_Composition.pdf", width = 6
)

p_celltype_composition <- barplots(
  celltype_counts, "Source", "cell_type", cell_type_cols,
  filename = "Figures_Revision/Cell_Type_Composition.pdf"
)

p_sample_celltype_composition <- barplots(
  sample_celltype_counts, "orig.ident", "cell_type", cell_type_cols,
  filename = "Figures_Revision/Cell_Type_Composition_By_Sample.pdf",
  width = 12, height = 7
)


###############################################################################
# 5. Feature plots
###############################################################################

DefaultAssay(gdt_seurat) <- "RNA"

feature_genes <- c(
  "TRDV1", "TRDV2", "TRGV9", "TRDV3", "TRGV8", "TRDC", "KLRC1", "KLRC2",
  "KLRB1", "GZMK", "GZMH", "CD8A", "TBX21", "RORA", "CEBPD", "ZEB2"
)

feature_genes <- intersect(feature_genes, rownames(gdt_seurat[["RNA"]]))

feature_plots <- lapply(feature_genes, function(gene) {
  FeaturePlot(
    gdt_seurat, features = gene, reduction = "umap.rpca",
    raster = TRUE, raster.dpi = c(800, 800),
    pt.size = 2.5, max.cutoff = "q99"
  ) +
    ggtitle(gene) + coord_fixed() +
    theme(
      axis.title = element_blank(), axis.text = element_blank(),
      axis.ticks = element_blank(), axis.line = element_blank()
    )
})

p_feature_grid <- plot_grid(plotlist = feature_plots, ncol = 4)

ggsave(
  "Figures_Revision/FeaturePlots.pdf",
  p_feature_grid, width = 16, height = 16
)



###############################################################################
# 6. DotPlot of cell-type marker genes
###############################################################################

dotplot_genes <- c(
  "TRDV1", "TRDV3", "TRGV8", "TRGV3", "TRGC2", "GZMH", "NKG7", "PRF1",
  "GNLY", "KLRD1", "TBX21", "EOMES", "TYROBP", "LAMP1", "KLRF1", "CX3CR1",
  "FCGR3A", "FCRL6", "NCR1", "CD244", "CD8A", "PDCD1", "KLRC2", "TCF7",
  "LEF1", "MAL", "ZNF683", "KLRK1", "IL2RB", "IKZF2", "AOAH", "PRKCH",
  "THEMIS", "STIM1", "CARD11", "BCL2", "TRDV2", "TRGV9", "CEBPD", "KLRB1",
  "GZMK", "IL7R", "KLRC1", "IL12RB2", "IL23R", "DPP4", "CCR6", "RORC",
  "ZBTB16", "RORA", "CD44", "NFAT5", "NFKB1", "ANXA1", "IFNG", "CD69",
  "ZEB2", "TIGIT", "CCL4", "CCR7", "HLA-DRB1", "HLA-DRA", "HLA-DQA2",
  "KIR2DL1", "CDC20", "AURKB", "MKI67", "CDK1"
)

dotplot_genes <- intersect(dotplot_genes, rownames(gdt_seurat[["RNA"]]))

p_dotplot <- DotPlot(
  gdt_seurat, features = dotplot_genes,
  group.by = "cell_type", assay = "RNA"
) +
  scale_color_gradientn(colors = rev(brewer.pal(11, "RdBu"))) +
  RotatedAxis() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_blank())

ggsave(
  "Figures_Revision/DotPlot_Cell_Type_Markers.pdf",
  p_dotplot, width = 18, height = 4
)


###############################################################################
# 7. PBMC vs Tumor cell-type frequencies
###############################################################################

# Calculate sample-level counts and percentages.
# Each sample is treated as the biological replicate.

frequency_data <- meta %>%
  mutate(
    Source = factor(Source, levels = source_order),
    cell_type = factor(cell_type, levels = cell_type_order)
  ) %>%
  count(orig.ident, Source, cell_type, name = "n") %>%
  complete(
    nesting(orig.ident, Source), cell_type,
    fill = list(n = 0)
  ) %>%
  group_by(orig.ident, Source) %>%
  mutate(total_cells = sum(n), percent = 100 * n / total_cells) %>%
  ungroup()


# Descriptive plot only; statistical testing is performed using propeller.

p_frequency <- ggplot(
  frequency_data,
  aes(x = cell_type, y = percent, fill = Source)
) +
  geom_boxplot(
    position = position_dodge(width = 0.8),
    outlier.shape = NA, width = 0.65
  ) +
  geom_point(
    aes(color = Source),
    position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8),
    size = 2, alpha = 0.8
  ) +
  scale_fill_manual(values = source_cols) +
  scale_color_manual(values = source_cols) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.title = element_blank()
  ) +
  labs(x = NULL, y = "Cells per sample (%)")

ggsave(
  "Figures_Revision/Frequency_Clusters_PBMC_Tumor.pdf",
  p_frequency, width = 12, height = 7
)


###############################################################################
# 8. propeller statistical testing
###############################################################################

# Sample-level compositional testing for PBMC vs Tumor.
# No additional Wilcoxon tests are performed.

run_propeller <- function(metadata, cluster_col,
                          sample_col = "orig.ident",
                          group_col = "Source") {
  metadata <- metadata %>%
    filter(
      .data[[group_col]] %in% source_order,
      !is.na(.data[[cluster_col]]),
      !is.na(.data[[sample_col]])
    )
  
  propeller_results <- speckle::propeller(
    clusters = metadata[[cluster_col]],
    sample = metadata[[sample_col]],
    group = factor(metadata[[group_col]], levels = source_order),
    transform = "asin"
  )
  
  result_table <- as.data.frame(propeller_results$tt)
  result_table$cluster <- rownames(result_table)
  result_table <- result_table %>% relocate(cluster)
  
  list(results = result_table, full_output = propeller_results)
}


# Cell-type and broad-class comparisons.

propeller_celltype <- run_propeller(meta, "cell_type")
propeller_broad <- run_propeller(meta, "broad_class")

write.table(
  propeller_celltype$results,
  "Tables/Propeller_Cell_Type_PBMC_vs_Tumor.tsv",
  sep = "\t", quote = FALSE, row.names = FALSE
)

write.table(
  propeller_broad$results,
  "Tables/Propeller_Broad_Class_PBMC_vs_Tumor.tsv",
  sep = "\t", quote = FALSE, row.names = FALSE
)


###############################################################################
# 9. Frequency source data export
###############################################################################

source_data_frequency <- frequency_data %>%
  select(cell_type, orig.ident, Source, n, percent) %>%
  pivot_wider(
    names_from = Source,
    values_from = c(n, percent),
    names_glue = "{.value}_{Source}",
    values_fill = 0
  ) %>%
  rename(
    Cell_Type = cell_type,
    Sample = orig.ident,
    Cell_Count_PBMC = n_PBMC,
    Percentage_PBMC = percent_PBMC,
    Cell_Count_Tumor = n_Tumor,
    Percentage_Tumor = percent_Tumor
  )

write_xlsx(
  source_data_frequency,
  "Figures_Revision/SourceData_Figure1c_and_1e.xlsx"
)


###############################################################################
# 10. Differential expression markers by cell type
###############################################################################

DefaultAssay(gdt_seurat) <- "RNA"
Idents(gdt_seurat) <- "cell_type"

markers2 <- FindAllMarkers(
  gdt_seurat, assay = "RNA", only.pos = TRUE
)

markers2$pct.diff <- markers2$pct.1 - markers2$pct.2

write.table(
  markers2,
  "Tables/Supplemental_Table_DEGs_scRNA.tsv",
  sep = "\t", quote = FALSE, row.names = FALSE
)


###############################################################################
# 11. Save annotated Seurat object
###############################################################################

saveRDS(
  gdt_seurat,
  "GDT_All_Processed_Revised13_GDT_QC_Integrated2_Annotated.RDS"
)

