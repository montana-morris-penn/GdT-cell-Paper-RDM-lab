###############################################################################
# GIST pseudobulk differential expression, PCA, and GSEA
###############################################################################

library(Seurat)
library(DESeq2)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(EnhancedVolcano)
library(msigdbr)
library(fgsea)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(KEGGREST)
library(BiocParallel)

setwd("/mnt/isilon/tan_lab/sussmanj/GIST_Multiome/GD_TCells")
dir.create("Figures_Revision", showWarnings = FALSE)
dir.create("Tables", showWarnings = FALSE)

###############################################################################
# 1. Load Seurat object and prepare metadata
###############################################################################

gdt <- readRDS("GDT_All_Processed_Revised13_GDT_QC_Integrated2.RDS")

if (!"Treatment" %in% colnames(gdt@meta.data)) gdt$Treatment <- gdt$Status

gdt$Treatment <- factor(gdt$Treatment, levels = c("Untreated", "Resistant"))
gdt$Source <- factor(gdt$Source, levels = c("PBMC", "Tumor"))
DefaultAssay(gdt) <- "RNA"

###############################################################################
# 2. Define pseudobulk aggregation and DESeq2 functions
###############################################################################

make_pseudobulk <- function(object, group_vars, metadata_vars) {
  counts <- AggregateExpression(
    object, group.by = group_vars, assays = "RNA",
    slot = "counts", return.seurat = FALSE
  )$RNA
  
  metadata <- object@meta.data %>%
    select(all_of(metadata_vars)) %>%
    distinct()
  
  colnames(counts) <- gsub("-", "_", colnames(counts))
  
  if (length(group_vars) == 1) {
    metadata$Sample <- gsub("-", "_", metadata$Sample)
    rownames(metadata) <- metadata$Sample
  } else {
    metadata$sample_cluster <- gsub("-", "_", metadata$sample_cluster)
    rownames(metadata) <- metadata$sample_cluster
  }
  
  metadata <- metadata[colnames(counts), , drop = FALSE]
  stopifnot(identical(rownames(metadata), colnames(counts)))
  
  list(counts = counts, metadata = metadata)
}

run_deseq <- function(pb, design) {
  dds <- DESeqDataSetFromMatrix(
    countData = round(pb$counts), colData = pb$metadata, design = design
  )
  DESeq(dds)
}

###############################################################################
# 3. PBMC versus Tumor pseudobulk analysis
#    Design: ~ Assay + Source + Treatment
###############################################################################

pb_all <- make_pseudobulk(
  gdt, group_vars = "Sample",
  metadata_vars = c("Sample", "Source", "Assay", "Treatment")
)

print(with(pb_all$metadata, table(Source, Treatment)))
print(table(pb_all$metadata$Assay))
stopifnot(!anyNA(pb_all$metadata[, c("Source", "Assay", "Treatment")]))

dds_source <- run_deseq(pb_all, design = ~ Assay + Source + Treatment)

# Tumor versus PBMC, adjusted for Assay and Treatment
res_source <- results(dds_source, contrast = c("Source", "Tumor", "PBMC")) %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Gene") %>%
  filter(!is.na(padj)) %>%
  arrange(padj)

write.table(
  res_source, "Tables/Pseudobulk_Tumor_vs_PBMC.tsv",
  sep = "\t", quote = FALSE, row.names = FALSE
)

###############################################################################
# 4. Tumor-only Resistant versus Untreated pseudobulk analysis
#    Design: ~ Assay + Treatment
###############################################################################

tumor <- subset(gdt, subset = Source == "Tumor")

pb_tumor <- make_pseudobulk(
  tumor, group_vars = "Sample",
  metadata_vars = c("Sample", "Assay", "Treatment")
)

print(table(pb_tumor$metadata$Treatment, pb_tumor$metadata$Assay))
stopifnot(!anyNA(pb_tumor$metadata[, c("Assay", "Treatment")]))

dds_tumor <- run_deseq(pb_tumor, design = ~ Assay + Treatment)

# Resistant versus Untreated, adjusted for Assay
res_tumor <- results(dds_tumor, contrast = c("Treatment", "Resistant", "Untreated")) %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Gene") %>%
  filter(!is.na(padj)) %>%
  arrange(padj)

write.table(
  res_tumor, "Tables/Pseudobulk_Tumor_Resistant_vs_Untreated.tsv",
  sep = "\t", quote = FALSE, row.names = FALSE
)

###############################################################################
# 5. PCA plotting function
###############################################################################

plot_pca <- function(dds, color_by, filename, title) {
  vst_data <- vst(dds, blind = TRUE)
  expression <- assay(vst_data)
  
  # Select the 3,000 most variable genes
  gene_variance <- apply(expression, 1, var)
  selected <- head(order(gene_variance, decreasing = TRUE), min(3000, length(gene_variance)))
  pca <- prcomp(t(expression[selected, ]))
  variance <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)
  
  pca_df <- data.frame(
    pca$x[, 1:3], Sample = rownames(pca$x),
    colData(vst_data)[, color_by, drop = FALSE]
  )
  
  plots <- lapply(list(c("PC1", "PC2"), c("PC2", "PC3")), function(axes) {
    ggplot(pca_df, aes(
      x = .data[[axes[1]]], y = .data[[axes[2]]],
      color = .data[[color_by]], label = Sample
    )) +
      geom_point(size = 3) +
      geom_text_repel(max.overlaps = Inf) +
      coord_fixed() +
      theme_bw() +
      labs(
        title = title,
        x = paste0(axes[1], ": ", variance[as.integer(sub("PC", "", axes[1]))], "% variance"),
        y = paste0(axes[2], ": ", variance[as.integer(sub("PC", "", axes[2]))], "% variance"),
        color = color_by
      )
  })
  
  p <- patchwork::wrap_plots(plots, ncol = 2)
  ggsave(filename, p, width = 10, height = 5)
  return(p)
}

###############################################################################
# 6. PCA: primary pseudobulk comparisons
###############################################################################

pca_source <- plot_pca(
  dds_source, color_by = "Source",
  filename = "Figures_Revision/PCA_Tumor_vs_PBMC.pdf",
  title = "Tumor versus PBMC"
)

pca_tumor <- plot_pca(
  dds_tumor, color_by = "Treatment",
  filename = "Figures_Revision/PCA_Tumor_Resistant_vs_Untreated.pdf",
  title = "Tumor: Resistant versus Untreated"
)

###############################################################################
# 7. PCA only: PBMC VD1+ versus VD2+
#    No differential expression or GSEA for this comparison
###############################################################################

pbmc_vd <- subset(gdt, subset = Source == "PBMC" & broad_class %in% c("VD1+", "VD2+"))

pb_vd <- make_pseudobulk(
  pbmc_vd, group_vars = c("Sample", "broad_class"),
  metadata_vars = c("Sample", "Source", "Assay", "Treatment", "broad_class")
)

dds_vd <- DESeqDataSetFromMatrix(
  countData = round(pb_vd$counts), colData = pb_vd$metadata, design = ~ 1
)

vst_vd <- vst(dds_vd, blind = TRUE)
expr_vd <- assay(vst_vd)
selected_vd <- head(order(apply(expr_vd, 1, var), decreasing = TRUE), 3000)

pca_vd <- prcomp(t(expr_vd[selected_vd, ]))

pca_vd_df <- data.frame(
  pca_vd$x[, 1:3], Sample = rownames(pca_vd$x),
  broad_class = pb_vd$metadata$broad_class
)

pca_vd_var <- round(100 * pca_vd$sdev^2 / sum(pca_vd$sdev^2), 1)

p_vd <- ggplot(pca_vd_df, aes(PC1, PC2, color = broad_class, label = Sample)) +
  geom_point(size = 3) +
  geom_text_repel(max.overlaps = Inf) +
  coord_fixed() +
  theme_bw() +
  labs(
    x = paste0("PC1: ", pca_vd_var[1], "% variance"),
    y = paste0("PC2: ", pca_vd_var[2], "% variance"),
    color = "Broad class"
  )

ggsave("Figures_Revision/PCA_PBMC_VD1_vs_VD2.pdf", p_vd, width = 6, height = 5)

###############################################################################
# 8. Volcano plots for the two primary comparisons
###############################################################################

plot_volcano <- function(results_df, title, filename, fc_cutoff = 0.5, p_cutoff = 0.05) {
  genes_to_label <- results_df %>%
    filter(padj < p_cutoff, abs(log2FoldChange) > fc_cutoff) %>%
    arrange(padj) %>%
    slice_head(n = 30) %>%
    pull(Gene)
  
  n_up <- sum(results_df$padj < p_cutoff & results_df$log2FoldChange > fc_cutoff, na.rm = TRUE)
  n_down <- sum(results_df$padj < p_cutoff & results_df$log2FoldChange < -fc_cutoff, na.rm = TRUE)
  
  p <- EnhancedVolcano(
    results_df, lab = results_df$Gene, selectLab = genes_to_label,
    x = "log2FoldChange", y = "padj", pCutoff = p_cutoff, FCcutoff = fc_cutoff,
    colAlpha = 0.85, pointSize = 0.75, labSize = 3.5, max.overlaps = 30,
    drawConnectors = TRUE, arrowheads = FALSE, colConnectors = "black",
    hline = NULL, vline = NULL, gridlines.major = FALSE, gridlines.minor = FALSE,
    raster = TRUE, title = title,
    subtitle = paste0("Upregulated: ", n_up, " | Downregulated: ", n_down)
  )
  
  ggsave(filename, p, width = 9.2, height = 9, dpi = 600)
  return(p)
}

volcano_source <- plot_volcano(
  res_source, title = "Tumor versus PBMC",
  filename = "Figures_Revision/Volcano_Tumor_vs_PBMC.pdf"
)

volcano_tumor <- plot_volcano(
  res_tumor, title = "Tumor: Resistant versus Untreated",
  filename = "Figures_Revision/Volcano_Tumor_Resistant_vs_Untreated.pdf"
)

###############################################################################
# 9. GSEA pathway preparation
###############################################################################

hallmark <- msigdbr(species = "Homo sapiens", category = "H")
reactome <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME")

pathways_hallmark <- split(hallmark$gene_symbol, hallmark$gs_name)
pathways_reactome <- split(reactome$gene_symbol, reactome$gs_name)

# KEGG pathways
kegg_raw <- keggLink("pathway", "hsa")
pathway_ids <- gsub("\\D", "", unname(kegg_raw))
entrez_ids <- sub("hsa:", "", names(kegg_raw))
kegg_to_entrez <- split(entrez_ids, pathway_ids)

pathways_kegg <- lapply(kegg_to_entrez, function(ids) {
  symbols <- mapIds(org.Hs.eg.db, keys = ids, column = "SYMBOL",
                    keytype = "ENTREZID", multiVals = "first")
  unique(na.omit(unname(symbols)))
})

kegg_titles <- keggList("pathway", "hsa")
title_ids <- gsub("\\D", "", names(kegg_titles))
clean_titles <- sub(" - Homo sapiens \\(human\\)", "", kegg_titles)
names(clean_titles) <- title_ids

matched_titles <- clean_titles[names(pathways_kegg)]
formatted_names <- gsub("[^A-Za-z0-9]+", "_", toupper(matched_titles))

names(pathways_kegg) <- paste0(
  "KEGG_",
  ifelse(is.na(formatted_names), paste0("PATHWAY_", names(pathways_kegg)), formatted_names)
)

combined_pathways <- c(pathways_hallmark, pathways_reactome, pathways_kegg)

###############################################################################
# 10. GSEA function and customized pathway plots
###############################################################################

run_gsea <- function(results_df, comparison, pathways, midpoint, width, height) {
  gene_list <- results_df$stat
  names(gene_list) <- results_df$Gene
  gene_list <- sort(gene_list[!is.na(gene_list) & !duplicated(names(gene_list))], decreasing = TRUE)
  
  fgsea_res <- fgsea(
    pathways = combined_pathways, stats = gene_list,
    minSize = 15, maxSize = 500, BPPARAM = SerialParam()
  ) %>%
    as.data.frame() %>%
    arrange(padj)
  
  write.table(
    fgsea_res, paste0("Tables/GSEA_", comparison, ".tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  
  plot_df <- fgsea_res %>%
    filter(pathway %in% pathways) %>%
    distinct(pathway, .keep_all = TRUE) %>%
    arrange(NES) %>%
    mutate(pathway = factor(pathway, levels = pathway),
           neg_log10_padj = -log10(padj))
  
  p <- ggplot(plot_df, aes(x = reorder(pathway, NES), y = NES, fill = -log10(padj))) +
    geom_bar(stat = "identity", position = "identity", width = 0.7) +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = midpoint) +
    coord_flip() +
    labs(x = "", y = "Normalized Enrichment Score", fill = "-log10(FDR)", title = "") +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 0, vjust = 1, hjust = 1, family = "sans", size = 8, color = "black"),
      axis.text.y = element_text(family = "sans", size = 12, color = "black"),
      axis.line = element_blank(),
      plot.margin = margin(l = 5, r = 5, b = 5, t = 5, unit = "pt")
    )
  
  print(p)
  ggsave(paste0("Figures_Revision/", comparison, "_GSEA_Barplot.pdf"),
         plot = p, device = "pdf", height = height, width = width)
  
  return(list(results = fgsea_res, plot = p))
}

###############################################################################
# 11. GSEA pathway selections
###############################################################################

pathways_source <- c(
  "HALLMARK_E2F_TARGETS",
  "KEGG_CYTOKINE_CYTOKINE_RECEPTOR_INTERACTION",
  "REACTOME_CHEMOKINE_RECEPTORS_BIND_CHEMOKINES",
  "HALLMARK_IL6_JAK_STAT3_SIGNALING",
  "HALLMARK_INFLAMMATORY_RESPONSE",
  "HALLMARK_IL2_STAT5_SIGNALING",
  "REACTOME_CDC42_GTPASE_CYCLE",
  "REACTOME_RAC1_GTPASE_CYCLE",
  "HALLMARK_ALLOGRAFT_REJECTION",
  "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "KEGG_IL_17_SIGNALING_PATHWAY",
  "REACTOME_FCGR_ACTIVATION",
  "HALLMARK_GLYCOLYSIS"
)

pathways_tumor <- c(
  "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "HALLMARK_MTORC1_SIGNALING",
  "KEGG_FATTY_ACID_METABOLISM",
  "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
  "REACTOME_HSF1_ACTIVATION",
  "KEGG_AUTOPHAGY_OTHER",
  "REACTOME_CELLULAR_RESPONSE_TO_HYPOXIA",
  "HALLMARK_FATTY_ACID_METABOLISM",
  "KEGG_INTEGRATED_STRESS_RESPONSE_ISR_SIGNALING_PATHWAY",
  "HALLMARK_UNFOLDED_PROTEIN_RESPONSE",
  "HALLMARK_GLYCOLYSIS",
  "KEGG_JAK_STAT_SIGNALING_PATHWAY",
  "REACTOME_EPH_EPHRIN_SIGNALING",
  "HALLMARK_IL2_STAT5_SIGNALING",
  "REACTOME_FCGAMMA_RECEPTOR_FCGR_DEPENDENT_PHAGOCYTOSIS",
  "KEGG_TNF_SIGNALING_PATHWAY",
  "REACTOME_EPHB_MEDIATED_FORWARD_SIGNALING",
  "KEGG_IL_17_SIGNALING_PATHWAY",
  "REACTOME_CDC42_GTPASE_CYCLE",
  "REACTOME_RAC1_GTPASE_CYCLE",
  "KEGG_CYTOKINE_CYTOKINE_RECEPTOR_INTERACTION",
  "KEGG_NF_KAPPA_B_SIGNALING_PATHWAY"
)

###############################################################################
# 12. GSEA: primary comparisons
###############################################################################

gsea_source <- run_gsea(
  res_source, comparison = "Tumor_vs_PBMC",
  pathways = pathways_source, midpoint = 10, width = 9, height = 5
)

gsea_tumor <- run_gsea(
  res_tumor, comparison = "Tumor_Resistant_vs_Untreated",
  pathways = pathways_tumor, midpoint = 5, width = 13, height = 5.5
)