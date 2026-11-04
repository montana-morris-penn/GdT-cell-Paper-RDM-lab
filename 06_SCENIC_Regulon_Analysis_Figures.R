###############################################################################
# SCENIC regulon analysis and figures
###############################################################################

library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(patchwork)
library(ggpubr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(writexl)

setwd("/mnt/isilon/tan_lab/sussmanj/GIST_Multiome/GD_TCells")

dir.create("Figures_Revision", showWarnings = FALSE, recursive = TRUE)

gdt_seurat <- readRDS(
  "GDT_All_Processed_Revised13_GDT_QC_Integrated2_SCENIC.RDS"
)

DefaultAssay(gdt_seurat) <- "SCENIC"


###############################################################################
# 1. Average regulon activity heatmap
###############################################################################

combined_averages <- AverageExpression(
  gdt_seurat,
  assays = "SCENIC",
  return.seurat = FALSE,
  layer = "data"
)$SCENIC

regulons <- c(
  "POLR3G(+)", "EOMES(+)", "ASCL2(+)", "TBX21(+)", "RUNX3(+)",
  "PRDM1(+)", "KLF2(+)", "ZBTB7A(+)", "PATZ1(+)", "LEF1(+)",
  "SOX4(+)", "TCF4(+)", "NFAT5(+)", "NR3C2(+)", "STAT1(+)",
  "ATF2(+)", "SMAD4(+)", "IKZF1(+)", "FOXO1(+)", "ZBTB46(+)",
  "TBL1XR1(+)", "ZNF493(+)", "MAF(+)", "STAT3(+)", "BATF(+)",
  "PPARA(+)", "CEBPD(+)", "BHLHE40(+)", "NR1D2(+)", "TFAP4(+)",
  "RORC(+)", "NR1D1(+)", "RUNX2(+)", "PPARG(+)", "SMAD3(+)",
  "NFKB1(+)", "FOSL2(+)", "IRF4(+)", "ZNF384(+)", "E2F2(+)",
  "E2F7(+)"
)

regulons <- intersect(unique(regulons), rownames(combined_averages))

pheatmap(
  combined_averages[regulons, ],
  scale = "row",
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  angle_col = 45,
  color = colorRampPalette(c("blue", "white", "red"))(100),
  filename = "Figures_Revision/SCENIC_Heatmap_Cluster.pdf",
  width = 5,
  height = 8.25
)


###############################################################################
# 2. Differential regulon activity
###############################################################################

DefaultAssay(gdt_seurat) <- "SCENIC"

# Tumor versus PBMC
Idents(gdt_seurat) <- "Source"

de_source <- FindMarkers(
  gdt_seurat,
  ident.1 = "Tumor",
  ident.2 = "PBMC",
  logfc.threshold = 0,
  min.pct = 0.001,
  max.cells.per.ident = 1000
)

de_source$gene <- rownames(de_source)

write.table(
  de_source,
  "SCENIC_Regulons_Source_Differential.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# Broad-class comparisons within PBMC
pbmc_subset <- subset(gdt_seurat, subset = Source == "PBMC")
Idents(pbmc_subset) <- "broad_class"

de_pbmc <- FindMarkers(
  pbmc_subset,
  ident.1 = "VD1+",
  ident.2 = "VD2+",
  logfc.threshold = 0,
  min.pct = 0.001,
  max.cells.per.ident = 1000
)

de_pbmc$gene <- rownames(de_pbmc)

write.table(
  de_pbmc,
  "SCENIC_Regulons_PBMC_VD1_vs_VD2_Differential.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


###############################################################################
# 3. Volcano plot: PBMC VD1+ versus VD2+
###############################################################################

red_regulons <- c(
  "TFEB(+)", "PATZ1(+)", "ZNF384(+)", "RFX1(+)", "HOXA9(+)",
  "SIRT6(+)", "ASCL2(+)", "POLR3A(+)", "BCL3(+)", "CEBPB(+)",
  "CEBPD(+)", "E2F4(+)", "KLF2(+)", "NR1D1(+)", "PRDM1(+)",
  "RORC(+)", "RUNX3(+)", "TBX21(+)", "TCF7(+)", "ZBTB7A(+)",
  "TFAP4(+)", "RUNX2(+)"
)

blue_regulons <- c(
  "NR3C2(+)", "RARA(+)", "GADD45A(+)", "ZBTB46(+)", "ATF2(+)",
  "STAT3(+)", "PPARA(+)", "NFAT5(+)", "SMAD3(+)", "FOSL2(+)",
  "FOXO1(+)", "NFKB1(+)", "NR3C1(+)", "REL(+)", "STAT1(+)",
  "TBL1XR1(+)", "ZEB1(+)", "ZNF493(+)", "NFIA(+)", "JDP2(+)"
)

tfs <- unique(c(red_regulons, blue_regulons))

de_pbmc <- de_pbmc %>%
  mutate(
    enriched = ifelse(avg_log2FC < 0, "VD2+", "VD1+")
  )

label_points <- de_pbmc %>%
  filter(gene %in% tfs)

gene_volcano <- ggplot(
  de_pbmc,
  aes(
    x = avg_log2FC,
    y = -log10(p_val_adj + 1e-300),
    color = enriched
  )
) +
  geom_point(color = "grey", alpha = 0.7) +
  geom_point(
    data = de_pbmc %>% filter(abs(avg_log2FC) > 0.01, p_val_adj < 1e-4)
  ) +
  geom_text_repel(
    data = label_points,
    aes(label = gene),
    size = 3,
    max.overlaps = 30
  ) +
  scale_color_manual(values = c("VD1+" = "#D85FF7", "VD2+" = "#AAF400")) +
  theme_bw() +
  labs(
    x = "Log2 fold change",
    y = "-log10 adjusted P value"
  ) +
  theme(legend.position = "none")

ggsave(
  "Figures_Revision/PBMC_VD1_vs_VD2_SCENIC_Volcano_Labels.pdf",
  gene_volcano,
  width = 6,
  height = 6
)


###############################################################################
# 4. FeaturePlots of selected regulons
###############################################################################

feature_regulons <- c(
  "CEBPD(+)", "PRDM1(+)", "RORC(+)",
  "EOMES(+)", "TBX21(+)", "NFAT5(+)"
)

feature_plots <- lapply(feature_regulons, function(regulon) {
  FeaturePlot(
    gdt_seurat,
    features = regulon,
    reduction = "umap.rpca",
    raster = TRUE,
    raster.dpi = c(800, 800),
    pt.size = 2.5,
    cols = c("grey", "red"),
    max.cutoff = "q99",
    min.cutoff = "q1"
  ) +
    coord_fixed() +
    ggtitle(gsub("\\(\\+\\)", "", regulon)) +
    theme(
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank()
    )
})

feature_grid <- wrap_plots(feature_plots, ncol = 3, nrow = 2)

ggsave(
  "Figures_Revision/FeaturePlot_Grid_SCENIC.pdf",
  feature_grid,
  width = 10,
  height = 7
)


###############################################################################
# 5. Sample-level analysis of eight selected regulons
###############################################################################

tfs <- c("RORC", "KLF10", "CEBPD", "PPARG", "PBX3", "TBX21", "MNT", "PPARA")

all_regs <- rownames(gdt_seurat[["SCENIC"]])

regulons <- setNames(
  vapply(tfs, function(tf) {
    hits <- grep(
      paste0("^", tf, "(\\(\\+\\))?$"),
      all_regs,
      value = TRUE
    )
    if (length(hits) > 0) hits[1] else NA_character_
  }, character(1)),
  tfs
)

if (anyNA(regulons)) {
  warning(
    "Missing regulons: ",
    paste(names(regulons)[is.na(regulons)], collapse = ", ")
  )
}

regulons <- regulons[!is.na(regulons)]

df_summary <- FetchData(
  gdt_seurat,
  vars = c("Sample", "broad_class", "Source", unname(regulons))
) %>%
  filter(
    broad_class %in% c("VD1+", "VD2+"),
    Source %in% c("PBMC", "Tumor")
  ) %>%
  group_by(Sample, broad_class, Source) %>%
  summarise(
    across(all_of(unname(regulons)), ~ median(.x, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  mutate(
    Group = factor(
      paste(broad_class, Source, sep = "_"),
      levels = c(
        "VD1+_PBMC", "VD2+_PBMC",
        "VD1+_Tumor", "VD2+_Tumor"
      )
    ),
    pair_id = paste(Sample, Source, sep = "_")
  )


###############################################################################
# 6. Statistical tests and plotting function
###############################################################################

# Paired tests: match VD1+ and VD2+ by Sample within each tissue
paired_p <- function(data, tissue) {
  wide <- data %>%
    filter(Source == tissue) %>%
    select(Sample, broad_class, AUC) %>%
    pivot_wider(names_from = broad_class, values_from = AUC) %>%
    drop_na(`VD1+`, `VD2+`)
  
  if (nrow(wide) < 2) return(NA_real_)
  
  t.test(wide$`VD1+`, wide$`VD2+`, paired = TRUE)$p.value
}


# Unpaired Welch tests: PBMC versus Tumor within each broad class
welch_p <- function(data, cell_class) {
  subset_data <- data %>% filter(broad_class == cell_class)
  
  if (length(unique(subset_data$Source)) < 2) return(NA_real_)
  
  t.test(AUC ~ Source, data = subset_data, var.equal = FALSE)$p.value
}


plot_regulon <- function(tf_name, reg_col) {
  
  df_tf <- df_summary %>%
    select(Sample, Source, broad_class, Group, pair_id, AUC = all_of(reg_col))
  
  raw_p <- c(
    paired_p(df_tf, "PBMC"),
    paired_p(df_tf, "Tumor"),
    welch_p(df_tf, "VD1+"),
    welch_p(df_tf, "VD2+")
  )
  
  adj_p <- p.adjust(raw_p, method = "BH")
  
  y_max <- max(df_tf$AUC, na.rm = TRUE)
  y_min <- min(df_tf$AUC, na.rm = TRUE)
  y_range <- ifelse(y_max == y_min, 1, y_max - y_min)
  
  stat_df <- data.frame(
    group1 = c(
      "VD1+_PBMC", "VD1+_Tumor",
      "VD1+_PBMC", "VD2+_PBMC"
    ),
    group2 = c(
      "VD2+_PBMC", "VD2+_Tumor",
      "VD1+_Tumor", "VD2+_Tumor"
    ),
    p.adj.format = ifelse(
      is.na(adj_p),
      "NA",
      paste0("q = ", signif(adj_p, 3))
    ),
    y.position = y_max + y_range * c(0.08, 0.08, 0.25, 0.42)
  )
  
  ggplot(df_tf, aes(x = Group, y = AUC, fill = broad_class)) +
    geom_line(aes(group = pair_id), color = "grey60", alpha = 0.6) +
    geom_boxplot(outlier.shape = NA, width = 0.5, alpha = 0.7) +
    geom_jitter(width = 0, size = 1.5) +
    scale_fill_manual(values = c("VD1+" = "#BD20BF", "VD2+" = "#9DA52F")) +
    scale_x_discrete(
      labels = c(
        "VD1+\nPBMC", "VD2+\nPBMC",
        "VD1+\nTumor", "VD2+\nTumor"
      )
    ) +
    stat_pvalue_manual(
      stat_df,
      label = "p.adj.format",
      size = 2.8,
      tip.length = 0.01
    ) +
    coord_cartesian(
      ylim = c(y_min - y_range * 0.05, y_max + y_range * 0.58)
    ) +
    theme_bw() +
    labs(
      title = tf_name,
      x = NULL,
      y = "Sample Median AUC"
    ) +
    theme(
      legend.position = "none",
      plot.title = element_text(face = "bold", hjust = 0.5, size = 12),
      axis.text = element_text(color = "black", size = 8),
      panel.grid = element_blank()
    )
}


###############################################################################
# 7. Generate eight-regulon figure
###############################################################################

plot_list <- lapply(names(regulons), function(tf) {
  plot_regulon(tf, regulons[[tf]])
})

names(plot_list) <- names(regulons)

final_grid <- wrap_plots(plot_list, ncol = 4)

ggsave(
  "Figures_Revision/SCENIC_8Regulons_Grid_FDR_Rearranged.pdf",
  final_grid,
  width = 9,
  height = 8
)


###############################################################################
# 8. Export source data
###############################################################################

source_data_scenic <- df_summary %>%
  select(Sample, Source, broad_class, all_of(unname(regulons))) %>%
  rename_with(
    ~ names(regulons)[match(.x, unname(regulons))],
    .cols = all_of(unname(regulons))
  ) %>%
  pivot_longer(
    cols = all_of(names(regulons)),
    names_to = "Regulon",
    values_to = "Sample_Median_AUC"
  ) %>%
  pivot_wider(
    names_from = c(Source, broad_class),
    values_from = Sample_Median_AUC,
    names_sep = "_"
  ) %>%
  arrange(Regulon, Sample) %>%
  select(
    Regulon,
    Sample,
    PBMC_VD1 = "PBMC_VD1+",
    PBMC_VD2 = "PBMC_VD2+",
    Tumor_VD1 = "Tumor_VD1+",
    Tumor_VD2 = "Tumor_VD2+"
  )

write_xlsx(
  list("Source Data" = source_data_scenic),
  "Figures_Revision/SourceData_Figure_SCENIC_Regulons.xlsx"
)


###############################################################################
# 9. GO enrichment example: TBX21 regulon targets
###############################################################################

regulon_targets <- readRDS(
  "pySCENIC_Output/Regulon_Target_Genes_Revision.rds"
)

write_xlsx(
  as.data.frame(
    lapply(regulon_targets, function(x) {
      length(x) <- max(lengths(regulon_targets))
      x
    }),
    check.names = FALSE
  ),
  "Figures_Revision/Regulon_targets_full.xlsx"
)

# Example: TBX21
symbols <- regulon_targets[["TBX21(+)"]]

entrezids <- mapIds(
  org.Hs.eg.db,
  keys = symbols,
  keytype = "SYMBOL",
  column = "ENTREZID",
  multiVals = "first"
)

entrezids <- unique(na.omit(unname(entrezids)))

ego <- enrichGO(
  gene = entrezids,
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.25,
  readable = TRUE
)

go_bp_list <- as.data.frame(ego) %>%
  filter(qvalue < 0.25) %>%
  mutate(logFDR = -log10(p.adjust))

tbx21_terms <- c(
  "immune response-activating cell surface receptor signaling pathway",
  "phagocytosis",
  "T cell receptor signaling pathway",
  "natural killer cell mediated immunity",
  "regulation of T cell activation"
)

go_bp_list <- go_bp_list %>%
  filter(Description %in% tbx21_terms) %>%
  arrange(logFDR) %>%
  mutate(Description = factor(Description, levels = Description))

p_go <- ggplot(go_bp_list, aes(x = Description, y = logFDR)) +
  geom_col(fill = "#E6E6E6") +
  coord_flip() +
  theme_bw() +
  labs(x = NULL, y = "-log10(FDR)")

ggsave(
  "Figures_Revision/TBX21_Selected_Targets_Revision.pdf",
  p_go,
  width = 7,
  height = 1.5
)
