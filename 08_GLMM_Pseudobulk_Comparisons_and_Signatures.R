###############################################################################
# GIST gamma-delta T-cell pseudobulk differential expression and signature scores
###############################################################################

library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(EnhancedVolcano)
library(edgeR)
library(limma)
library(variancePartition)
library(DESeq2)
library(msigdbr)
library(fgsea)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(KEGGREST)
library(BiocParallel)
library(ggpubr)
library(patchwork)
library(writexl)

setwd("/mnt/isilon/tan_lab/sussmanj/GIST_Multiome/GD_TCells")
dir.create("Figures_Revision", showWarnings = FALSE)
dir.create("Tables", showWarnings = FALSE)

###############################################################################
# Load data and prepare pseudobulk counts and metadata
###############################################################################

gdt <- readRDS("GDT_All_Processed_Revised13_GDT_QC_Integrated2.RDS")

gdt@meta.data <- gdt@meta.data %>%
  left_join(gdt_metadata, by = "Sample")

pb_counts <- AggregateExpression(
  gdt,
  assays = "RNA",
  group.by = c("Sample", "broad_class"),
  slot = "counts",
  return.seurat = FALSE
)$RNA

sample_cluster <- gdt@meta.data %>%
  distinct(Sample, Source, Assay, Status, Location, broad_class) %>%
  mutate(sample_cluster = paste0(Sample, "_", broad_class)) %>%
  mutate(sample_cluster = gsub("-", "_", sample_cluster))

colnames(pb_counts) <- gsub("-", "_", colnames(pb_counts))
rownames(sample_cluster) <- sample_cluster$sample_cluster
pb_metadata <- sample_cluster[colnames(pb_counts), , drop = FALSE]

stopifnot(identical(rownames(pb_metadata), colnames(pb_counts)))

pb_metadata$Source <- factor(pb_metadata$Source)
pb_metadata$Status <- factor(pb_metadata$Status)
pb_metadata$broad_class <- factor(pb_metadata$broad_class, levels = c("VD2+", "VD1+"))

###############################################################################
# Prepare pathway collections
###############################################################################

pathways_hallmark <- msigdbr(species = "Homo sapiens", category = "H") %>%
  split(x = .$gene_symbol, f = .$gs_name)

pathways_reactome <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME") %>%
  split(x = .$gene_symbol, f = .$gs_name)

kegg_genes <- keggLink("pathway", "hsa")
kegg_names <- keggList("pathway", "hsa")

kegg_df <- data.frame(
  pathway = sub("path:", "", unname(kegg_genes)),
  entrez = sub("hsa:", "", names(kegg_genes)),
  stringsAsFactors = FALSE
)

kegg_df$symbol <- mapIds(
  org.Hs.eg.db,
  keys = kegg_df$entrez,
  column = "SYMBOL",
  keytype = "ENTREZID",
  multiVals = "first"
)

kegg_df <- kegg_df %>%
  filter(!is.na(symbol)) %>%
  mutate(pathway_name = unname(kegg_names[pathway]))

pathways_kegg <- split(kegg_df$symbol, kegg_df$pathway_name)
names(pathways_kegg) <- paste0(
  "KEGG_",
  toupper(gsub("[^A-Za-z0-9]+", "_", names(pathways_kegg)))
)

combined_pathways <- c(pathways_hallmark, pathways_reactome, pathways_kegg)

###############################################################################
# Part 1: Differential expression using limma-voom and dream
###############################################################################

run_dream_de <- function(counts, metadata, tissue, formula, output_prefix) {
  keep_samples <- metadata$Source == tissue
  counts_sub <- counts[, keep_samples, drop = FALSE]
  metadata_sub <- droplevels(metadata[keep_samples, , drop = FALSE])
  
  metadata_sub$broad_class <- relevel(metadata_sub$broad_class, ref = "VD2+")
  
  dge <- DGEList(counts = counts_sub)
  keep_genes <- rowSums(cpm(dge) > 1) >= 3
  dge <- dge[keep_genes, , keep.lib.sizes = FALSE]
  dge <- calcNormFactors(dge)
  
  vobj <- voomWithDreamWeights(dge, formula, metadata_sub)
  fit <- dream(vobj, formula, metadata_sub)
  fit <- eBayes(fit)
  
  results <- topTable(
    fit,
    coef = "broad_classVD1+",
    number = Inf,
    sort.by = "P"
  ) %>%
    tibble::rownames_to_column("gene") %>%
    arrange(P.Value)
  
  write.table(
    results,
    file = paste0("Tables/", output_prefix, "_DE.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  
  results
}

tumor_de <- run_dream_de(
  counts = pb_counts,
  metadata = pb_metadata,
  tissue = "Tumor",
  formula = ~ broad_class + Assay + Status + (1 | Sample),
  output_prefix = "Tumor_VD1_vs_VD2"
)

pbmc_de <- run_dream_de(
  counts = pb_counts,
  metadata = pb_metadata,
  tissue = "PBMC",
  formula = ~ broad_class + Status + (1 | Sample),
  output_prefix = "PBMC_VD1_vs_VD2"
)

###############################################################################
# Volcano plots
###############################################################################

make_volcano <- function(de_results, title, filename, xlim, ylim) {
  de_results <- de_results %>%
    mutate(
      significance = case_when(
        adj.P.Val < 0.05 & logFC > 0.5 ~ "VD1+ higher",
        adj.P.Val < 0.05 & logFC < -0.5 ~ "VD2+ higher",
        TRUE ~ "Not significant"
      )
    )
  
  p <- ggplot(de_results, aes(x = logFC, y = -log10(adj.P.Val), color = significance)) +
    geom_point(size = 1.5, alpha = 0.8) +
    scale_color_manual(values = c(
      "VD2+ higher" = "#AAF400",
      "VD1+ higher" = "#D85FF7",
      "Not significant" = "grey75"
    )) +
    geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "grey50") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50") +
    coord_cartesian(xlim = xlim, ylim = ylim) +
    labs(title = title, x = "Log2 fold change (VD1+ vs VD2+)", y = "-Log10 adjusted P value", color = NULL) +
    theme_bw() +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10, color = "black"),
      legend.position = "top"
    )
  
  ggsave(filename, p, width = 6, height = 5)
  p
}

volcano_tumor <- make_volcano(
  tumor_de,
  "Tumor: VD1+ vs VD2+",
  "Figures_Revision/Volcano_Tumor_VD1_vs_VD2.pdf",
  xlim = c(-6, 4.5),
  ylim = c(0, 4)
)

volcano_pbmc <- make_volcano(
  pbmc_de,
  "PBMC: VD1+ vs VD2+",
  "Figures_Revision/Volcano_PBMC_VD1_vs_VD2.pdf",
  xlim = c(-5, 5),
  ylim = c(0, 3.5)
)

###############################################################################
# GSEA
###############################################################################

tumor_gsea_pathways <- c(
  "REACTOME_EFFECTS_OF_PIP2_HYDROLYSIS",
  "REACTOME_IMMUNOREGULATORY_INTERACTIONS_BETWEEN_A_LYMPHOID_AND_A_NON_LYMPHOID_CELL",
  "REACTOME_FCGR_ACTIVATION",
  "KEGG_ANTIGEN_PROCESSING_AND_PRESENTATION",
  "KEGG_NATURAL_KILLER_CELL_MEDIATED_CYTOTOXICITY",
  "HALLMARK_GLYCOLYSIS",
  "KEGG_IL_17_SIGNALING_PATHWAY",
  "KEGG_INTEGRATED_STRESS_RESPONSE_ISR_SIGNALING_PATHWAY",
  "HALLMARK_IL2_STAT5_SIGNALING",
  "HALLMARK_UNFOLDED_PROTEIN_RESPONSE",
  "HALLMARK_APOPTOSIS",
  "REACTOME_CO_INHIBITION_BY_PD_1",
  "HALLMARK_HYPOXIA",
  "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "KEGG_NECROPTOSIS",
  "REACTOME_INTERLEUKIN_12_SIGNALING",
  "REACTOME_INTERLEUKIN_1_SIGNALING"
)

pbmc_gsea_pathways <- c(
  "KEGG_FC_GAMMA_R_MEDIATED_PHAGOSOME_FORMATION",
  "KEGG_NATURAL_KILLER_CELL_MEDIATED_CYTOTOXICITY",
  "KEGG_CHEMOKINE_SIGNALING_PATHWAY",
  "KEGG_OXIDATIVE_PHOSPHORYLATION",
  "REACTOME_MITOCHONDRIAL_TRANSLATION_ELONGATION",
  "KEGG_RIBOSOME",
  "KEGG_LEUKOCYTE_TRANSENDOTHELIAL_MIGRATION"
)

run_gsea <- function(de_results, selected_pathways, output_prefix, plot_width = 13, plot_height = 3.5) {
  ranks <- de_results$t
  names(ranks) <- de_results$gene
  ranks <- sort(ranks[!is.na(ranks)], decreasing = TRUE)
  
  pathways <- combined_pathways[names(combined_pathways) %in% selected_pathways]
  
  gsea <- fgsea(
    pathways = pathways,
    stats = ranks,
    minSize = 15,
    maxSize = 500,
    BPPARAM = SerialParam()
  ) %>%
    arrange(padj)
  
  gsea_export <- gsea %>%
    mutate(leadingEdge = vapply(leadingEdge, paste, collapse = ";", character(1)))
  
  write.table(
    gsea_export,
    file = paste0("Tables/", output_prefix, "_GSEA.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  
  plot_df <- gsea %>%
    filter(!is.na(padj), padj < 0.25) %>%
    mutate(pathway = factor(pathway, levels = pathway[order(NES)]))
  
  p <- ggplot(plot_df, aes(x = reorder(pathway, NES), y = NES, fill = -log10(padj))) +
    geom_bar(stat = "identity", position = "identity", width = 0.7) +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 10) +
    coord_flip() +
    labs(x = NULL, y = "Normalized enrichment score", fill = "-Log10 FDR") +
    theme_bw() +
    theme(
      axis.text.x = element_text(size = 8, color = "black"),
      axis.text.y = element_text(size = 12, color = "black"),
      axis.line = element_blank(),
      plot.margin = margin(5, 5, 5, 5)
    )
  
  ggsave(
    paste0("Figures_Revision/", output_prefix, "_GSEA.pdf"),
    p,
    width = plot_width,
    height = plot_height
  )
  
  gsea
}

tumor_gsea <- run_gsea(tumor_de, tumor_gsea_pathways, "Tumor_VD1_vs_VD2")
pbmc_gsea <- run_gsea(pbmc_de, pbmc_gsea_pathways, "PBMC_VD1_vs_VD2")

###############################################################################
# Part 2: Pseudobulk signature scores
###############################################################################

signatureScore_zscore_robust <- function(matrix, pfeatures, nfeatures = NULL) {
  pfeatures <- intersect(pfeatures, rownames(matrix))
  
  rzscore <- function(x) {
    sd.robust <- sd(x) / 2 + mad(x) / 2
    if (is.na(sd.robust) || sd.robust <= 0) return(x - median(x))
    (x - median(x)) / sd.robust
  }
  
  pscore <- matrix[pfeatures, , drop = FALSE]
  pscore <- t(apply(pscore, 1, rzscore))
  pscore <- Matrix::colSums(pscore)
  
  if (!is.null(nfeatures) && length(nfeatures) > 0) {
    nfeatures <- intersect(nfeatures, rownames(matrix))
    nscore <- matrix[nfeatures, , drop = FALSE]
    nscore <- t(apply(nscore, 1, rzscore))
    nscore <- Matrix::colSums(nscore)
    return((pscore - nscore) / length(c(pfeatures, nfeatures)))
  }
  
  pscore / length(pfeatures)
}

dds <- DESeqDataSetFromMatrix(
  countData = round(pb_counts),
  colData = pb_metadata,
  design = ~ 1
)

vsd <- vst(dds, blind = TRUE)
vst_mat <- assay(vsd)

signature_scores <- as.data.frame(vapply(
  combined_pathways,
  function(genes) signatureScore_zscore_robust(vst_mat, genes),
  FUN.VALUE = numeric(ncol(vst_mat))
))

signature_meta <- pb_metadata %>%
  mutate(sample_cluster = rownames(pb_metadata)) %>%
  bind_cols(signature_scores)

###############################################################################
# Signature comparison helper functions
###############################################################################

test_signature <- function(df, group1, group2, paired = FALSE) {
  if (paired) {
    paired_df <- df %>%
      filter(Group %in% c(group1, group2)) %>%
      select(Sample, Group, Score) %>%
      pivot_wider(names_from = Group, values_from = Score) %>%
      drop_na(all_of(c(group1, group2)))
    
    if (nrow(paired_df) < 2) return(NA_real_)
    
    return(t.test(paired_df[[group1]], paired_df[[group2]], paired = TRUE)$p.value)
  }
  
  x <- df %>% filter(Group == group1) %>% pull(Score)
  y <- df %>% filter(Group == group2) %>% pull(Score)
  
  if (length(x) < 2 || length(y) < 2) return(NA_real_)
  t.test(x, y)$p.value
}

make_signature_plot <- function(data, pathways, comparison) {
  plot_list <- list()
  
  if (comparison == "PBMC_vs_Tumor") {
    data <- data %>%
      filter(Source %in% c("PBMC", "Tumor")) %>%
      mutate(
        Group = paste0(broad_class, "_", Source),
        Group = factor(Group, levels = c("VD1+_PBMC", "VD2+_PBMC", "VD1+_Tumor", "VD2+_Tumor"))
      )
    
    comparisons <- list(
      c("VD1+_PBMC", "VD2+_PBMC"),
      c("VD1+_Tumor", "VD2+_Tumor"),
      c("VD1+_PBMC", "VD1+_Tumor"),
      c("VD2+_PBMC", "VD2+_Tumor")
    )
    
    paired_flags <- c(TRUE, TRUE, FALSE, FALSE)
    fill_values <- c("VD1+" = "#BD20BF", "VD2+" = "#9DA52F")
    data$Class <- factor(as.character(data$broad_class), levels = c("VD1+", "VD2+"))
    data$Source <- factor(data$Source, levels = c("PBMC", "Tumor"))
    x_mapping <- "Group"
  } else {
    data <- data %>%
      filter(Source == "Tumor", Status %in% c("Untreated", "Resistant")) %>%
      mutate(
        Group = paste0(broad_class, "_", Status),
        Group = factor(Group, levels = c("VD1+_Untreated", "VD1+_Resistant", "VD2+_Untreated", "VD2+_Resistant"))
      )
    
    comparisons <- list(
      c("VD1+_Untreated", "VD1+_Resistant"),
      c("VD2+_Untreated", "VD2+_Resistant"),
      c("VD1+_Untreated", "VD2+_Untreated"),
      c("VD1+_Resistant", "VD2+_Resistant")
    )
    
    paired_flags <- c(FALSE, FALSE, TRUE, TRUE)
    fill_values <- c("VD1+" = "#BD20BF", "VD2+" = "#9DA52F")
    data$Class <- factor(as.character(data$broad_class), levels = c("VD1+", "VD2+"))
    data$Status <- factor(data$Status, levels = c("Untreated", "Resistant"))
    x_mapping <- "Group"
  }
  
  for (pathway in pathways) {
    plot_df <- data %>%
      select(Sample, Source, Status, broad_class, Class, Group, all_of(pathway)) %>%
      rename(Score = all_of(pathway)) %>%
      filter(!is.na(Score))
    
    pvals <- lapply(seq_along(comparisons), function(i) {
      p <- test_signature(plot_df, comparisons[[i]][1], comparisons[[i]][2], paired_flags[i])
      data.frame(
        group1 = comparisons[[i]][1],
        group2 = comparisons[[i]][2],
        p = p
      )
    }) %>%
      bind_rows() %>%
      mutate(padj = p.adjust(p, method = "BH"))
    
    y_max <- max(plot_df$Score, na.rm = TRUE)
    y_min <- min(plot_df$Score, na.rm = TRUE)
    y_range <- max(y_max - y_min, 0.5)
    
    pvals <- pvals %>%
      mutate(
        y.position = y_max + y_range * (0.12 + 0.14 * (row_number() - 1)),
        label = ifelse(is.na(padj), "NA", paste0("FDR = ", signif(padj, 2)))
      )
    
    p <- ggplot(plot_df, aes(x = Group, y = Score, fill = Class)) +
      geom_boxplot(outlier.shape = NA, alpha = 0.75, width = 0.65) +
      geom_jitter(aes(color = Class), width = 0.12, size = 1.8, alpha = 0.8) +
      scale_fill_manual(values = fill_values) +
      scale_color_manual(values = fill_values) +
      stat_pvalue_manual(pvals, label = "label", y.position = "y.position", tip.length = 0.01, size = 3) +
      labs(title = pathway, x = NULL, y = "Robust z-score signature", fill = NULL, color = NULL) +
      theme_bw() +
      theme(
        plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
        axis.text.x = element_text(size = 8, angle = 45, hjust = 1, color = "black"),
        axis.text.y = element_text(size = 8, color = "black"),
        axis.title.y = element_text(size = 9),
        legend.position = "none"
      )
    
    if (comparison == "PBMC_vs_Tumor") {
      p <- p + geom_line(
        data = plot_df,
        aes(x = Group, y = Score, group = interaction(Sample, Source)),
        inherit.aes = FALSE,
        color = "grey75",
        alpha = 0.6
      )
    } else {
      p <- p + geom_line(
        data = plot_df,
        aes(x = Group, y = Score, group = interaction(Sample, Status)),
        inherit.aes = FALSE,
        color = "grey75",
        alpha = 0.6
      )
    }
    
    plot_list[[pathway]] <- p
  }
  
  wrap_plots(plot_list, ncol = 5)
}

###############################################################################
# PBMC vs Tumor signature comparisons
###############################################################################

pbmc_tumor_pathways <- c(
  "KEGG_CYTOKINE_CYTOKINE_RECEPTOR_INTERACTION",
  "HALLMARK_INFLAMMATORY_RESPONSE",
  "REACTOME_RAC1_GTPASE_CYCLE",
  "KEGG_IL_17_SIGNALING_PATHWAY",
  "REACTOME_EFFECTS_OF_PIP2_HYDROLYSIS",
  "REACTOME_IMMUNOREGULATORY_INTERACTIONS_BETWEEN_A_LYMPHOID_AND_A_NON_LYMPHOID_CELL",
  "KEGG_NATURAL_KILLER_CELL_MEDIATED_CYTOTOXICITY",
  "HALLMARK_APOPTOSIS",
  "HALLMARK_HYPOXIA",
  "HALLMARK_OXIDATIVE_PHOSPHORYLATION"
)

pbmc_tumor_plot <- make_signature_plot(
  signature_meta,
  pbmc_tumor_pathways,
  "PBMC_vs_Tumor"
)

ggsave(
  "Figures_Revision/Signature_Zscore_10Pathways_PBMCTumor_FDR.pdf",
  pbmc_tumor_plot,
  width = 14,
  height = 7.5
)

pbmc_tumor_source <- signature_meta %>%
  filter(Source %in% c("PBMC", "Tumor")) %>%
  select(Sample, Source, Status, broad_class, all_of(pbmc_tumor_pathways))

write.table(
  pbmc_tumor_source,
  "Tables/SourceData_Signatures_PBMC_vs_Tumor.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

###############################################################################
# Tumor Resistant vs Untreated signature comparisons
###############################################################################

treatment_pathways <- c(
  "KEGG_FATTY_ACID_METABOLISM",
  "KEGG_INTEGRATED_STRESS_RESPONSE_ISR_SIGNALING_PATHWAY",
  "HALLMARK_GLYCOLYSIS",
  "HALLMARK_HYPOXIA",
  "REACTOME_EPH_EPHRIN_SIGNALING",
  "KEGG_TNF_SIGNALING_PATHWAY",
  "KEGG_IL_17_SIGNALING_PATHWAY",
  "REACTOME_RAC1_GTPASE_CYCLE",
  "KEGG_CYTOKINE_CYTOKINE_RECEPTOR_INTERACTION"
)

treatment_plot <- make_signature_plot(
  signature_meta,
  treatment_pathways,
  "Tumor_Resistant_vs_Untreated"
)

ggsave(
  "Figures_Revision/Signature_Zscore_9Pathways_Status_FDR_5x2.pdf",
  treatment_plot,
  width = 14,
  height = 7.5
)

treatment_source <- signature_meta %>%
  filter(Source == "Tumor", Status %in% c("Untreated", "Resistant")) %>%
  select(Sample, Source, Status, broad_class, all_of(treatment_pathways))

write.table(
  treatment_source,
  "Tables/SourceData_Signatures_Tumor_Resistant_vs_Untreated.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

###############################################################################
# Export signature source data to Excel
###############################################################################

write_xlsx(
  list(
    PBMC_vs_Tumor = pbmc_tumor_source,
    Tumor_Resistant_vs_Untreated = treatment_source
  ),
  "Figures_Revision/SourceData_Figure_Signatures.xlsx"
)