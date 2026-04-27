#Set WD
setwd("~/Library/CloudStorage/Box-Box/TAM_sc/")
options(future.globals.maxSize = 12 * 1024^3)
library(ggpubr)
library(patchwork)

#Load Libraries ----
library(BiocManager)
library(tidyverse) 
library(Seurat)
library(SeuratWrappers)
library(SeuratData)
library(sctransform)
library(harmony)
library(ggpubr)
library(ggstatsplot) 
library(patchwork)
library(cowplot) #organize plots
library(DT) #for tables
library(RColorBrewer)
library(dittoSeq)
#library(CytoTRACE)
#library(monocle3)
library(HDF5Array)
library(Matrix)
library(Matrix.utils)
library(DESeq2)
library(pheatmap)
library(EnhancedVolcano)
library(msigdbr)
library(fgsea)
library(viridis)
library(writexl)
library(SingleCellExperiment)
library(ashr)
library(gprofiler2)
library(dplyr)
library(tidyr)
library(SCENIC)

#load environment/objects
setwd("/mnt/isilon/tan_lab/sussmanj/GIST_Multiome/GD_TCells")
options(future.globals.maxSize= 9891289600)

se_HG_RNA.Harmony <- readRDS("/mnt/isilon/tan_lab/sussmanj/GIST_Multiome/Reference_Mapping/se_HG_RNA.Harmony_annot.rds")
se_HG_RNA.Harmony$cell_type = Idents(se_HG_RNA.Harmony)

se_SC <- subset(se_HG_RNA.Harmony, subset = assay == "scRNAseq")
target_idents <- c("HG590_T", "HG607_T", "HG610_T", "HG589_T", "HG605_T")
selected_subset <- subset(se_SC, subset = orig.ident %in% target_idents)

DefaultAssay(selected_subset) = "RNA"
selected_subset[["SCT"]] = NULL
selected_subset[["RNA"]] <- split(selected_subset[["RNA"]], f = selected_subset$orig.ident)
selected_subset <- SCTransform(selected_subset)
selected_subset <- RunPCA(selected_subset)

selected_subset <- IntegrateLayers(object = selected_subset, method = RPCAIntegration, new.reduction = "integrated.rpca", 
                                   normalization.method = "SCT", verbose = T)
selected_subset <- RunUMAP(selected_subset, dims = 1:30, reduction = "integrated.rpca", reduction.name = "rpca.umap")
FeaturePlot(selected_subset, features = c("CD3E"), reduction = "rpca.umap", order = T) + coord_fixed()
selected_subset <- FindNeighbors(selected_subset, dims = 1:30, reduction = "integrated.rpca")
selected_subset <- FindClusters(selected_subset, resolution = c(0.2, 0.3, 0.5))
DimPlot(selected_subset, group.by = "SCT_snn_res.0.2", reduction = "rpca.umap", label = T) + coord_fixed()
DimPlot(selected_subset, group.by = "orig.ident", reduction = "rpca.umap", label = T) + coord_fixed()


t_cell_subset = subset(selected_subset, subset = SCT_snn_res.0.2 %in% c(1,5,3,7,8))
DefaultAssay(t_cell_subset) = "RNA"
t_cell_subset = JoinLayers(t_cell_subset)
t_cell_subset[["SCT"]] = NULL
t_cell_subset[["RNA"]] <- split(t_cell_subset[["RNA"]], f = t_cell_subset$orig.ident)
t_cell_subset <- SCTransform(t_cell_subset)
t_cell_subset <- RunPCA(t_cell_subset)

t_cell_subset <- IntegrateLayers(object = t_cell_subset, method = RPCAIntegration, new.reduction = "integrated.rpca", 
                                 normalization.method = "SCT", verbose = T)
t_cell_subset <- RunUMAP(t_cell_subset, dims = 1:30, reduction = "integrated.rpca", reduction.name = "rpca.umap")
FeaturePlot(t_cell_subset, features = c("TRDC", "TRGC1", "TRGC2", "TRDV1", "TRDV2", "TRGV9"), 
            reduction = "rpca.umap", order = T, ncol = 3) & coord_fixed()
tmarkers = c("CD4", "IL7R", "FOXP3", "CD8A", "CD8B", "GZMB", "PRF1", "NKG7", "CCR7")
FeaturePlot(t_cell_subset, features = tmarkers, 
            reduction = "rpca.umap", order = T, ncol = 3) & coord_fixed()
t_cell_subset <- FindNeighbors(t_cell_subset, dims = 1:30, reduction = "integrated.rpca")
DefaultAssay(t_cell_subset) = "SCT"
t_cell_subset <- FindClusters(t_cell_subset, resolution = c(0.2, 0.3, 0.5))
t_cell_subset <- FindClusters(t_cell_subset, resolution = c(0.25))

DimPlot(t_cell_subset, group.by = "SCT_snn_res.0.2", reduction = "rpca.umap", label = T) + coord_fixed()
DimPlot(t_cell_subset, group.by = "orig.ident", reduction = "rpca.umap", label = T) + coord_fixed()


reference<-readRDS("GDT_Seurat_Annotated_CytoTRACE_SCENIC_05292024.RDS")
DefaultAssay(reference) = "RNA"

DimPlot(t_cell_subset, group.by = "SCT_snn_res.0.2", reduction = "rpca.umap", label = T) + coord_fixed()
DimPlot(t_cell_subset, group.by = "orig.ident", reduction = "rpca.umap", label = T) + coord_fixed()
DimPlot(t_cell_subset, group.by = "condition", reduction = "rpca.umap", label = T) + coord_fixed()

#Gamma delta score 
t_cell_subset = JoinLayers(t_cell_subset)
DefaultAssay(t_cell_subset) = "RNA"
gd_t_genes <- list(c("TRDC", "TRGC1", "TRGC2", "TRDV1", "TRDV2", "TRGV9"))
t_cell_subset <- AddModuleScore(object = t_cell_subset,features = gd_t_genes, name = "gdT_Score")
p1 <- FeaturePlot(t_cell_subset, features = "gdT_Score1", max.cutoff = 'q99',reduction = "rpca.umap", order = T) +
  ggtitle("Gamma Delta T Cell Module Score") + coord_fixed()
p2 <- VlnPlot(t_cell_subset, features = "gdT_Score1",  
              group.by = "SCT_snn_res.0.2", pt.size = 0) +
  geom_boxplot(width = 0.2, fill = "white", outlier.shape = NA) + 
  theme(legend.position = "none") +
  ggtitle("Gamma Delta Score by Cluster")
p1 | p2

#subset gdt cells----
library(ggrepel)
DefaultAssay(t_cell_subset) = "RNA"
gdsubset <- subset(t_cell_subset, subset = SCT_snn_res.0.2 %in% c(1, 3))
DimPlot(gdsubset, reduction = "rpca.umap")
gdsubset = NormalizeData(gdsubset)


markers = FindMarkers(gdsubset, ident.1 = "resistant", ident.2 = "untreated", 
                      group.by = "condition", max.cells.per.ident = 1000)
markers$gene <- rownames(markers)

#write_xlsx(markers, "untx_v_res_gdts_DEG.xlsx")

#volcano plot----
library(ggplot2)
library(dplyr)
library(ggrepel)

plot_volcano_pub <- function(
    markers,
    gene_col = "gene",
    fc_col = "avg_log2FC",
    p_col = "p_val_adj",
    
    fc_cutoff = 0.5,
    p_cutoff = 0.05,
    
    select_genes_1 = NULL,
    select_genes_2 = NULL,
    label_color_1 = "black",
    label_color_2 = "black",
    
    auto_label_n = 15,
    label_only_sig = TRUE,
    label_alpha_sig = 1,
    label_alpha_nonsig = 0.4,
    
    color_up = "#D55E00",
    color_down = "#0072B2",
    color_ns = "gray80",
    
    point_size = 2,
    point_alpha_ns = 0.5,
    point_alpha_sig = 0.9,
    
    label_size = 4,
    label_fontface = "italic",
    max_overlaps = Inf,
    box_padding = 0.5,
    point_padding = 0.35,
    force = 2,
    force_pull = 0.5,
    max_time = 5,
    max_iter = 20000,
    min_segment_length = 0,
    
    show_connectors = TRUE,
    connector_color = "black",
    connector_size = 0.4,
    connector_alpha = 0.8,
    connector_linetype = "solid",
    
    label_box = TRUE,
    label_box_fill = "white",
    label_box_color = "black",
    label_box_label_size = 0.2,
    label_r = unit(0.12, "lines"),
    label_padding = unit(0.15, "lines"),
    
    show_cutoff_lines = TRUE,
    cutoff_line_color = "black",
    cutoff_line_type = "dashed",
    cutoff_line_size = 0.6,
    
    title = NULL,
    subtitle = NULL,
    xlab = "Log2FC",
    ylab = "-log10(FDR)",
    
    title_size = 16,
    subtitle_size = 13,
    axis_title_size = 14,
    axis_text_size = 12,
    axis_line_size = 0.8,
    border_size = 0.8,
    
    x_breaks = NULL,
    x_limits = NULL,
    y_limits = NULL,
    
    show_legend = FALSE,
    seed = 123
) {
  df <- markers %>%
    as.data.frame() %>%
    mutate(
      .gene = .data[[gene_col]],
      .fc   = .data[[fc_col]],
      .p    = .data[[p_col]]
    ) %>%
    filter(!is.na(.gene), !is.na(.fc), !is.na(.p)) %>%
    mutate(
      .neglog10p = -log10(pmax(.p, 1e-300)),
      sig_dir = case_when(
        .p < p_cutoff & .fc >=  fc_cutoff ~ "Up",
        .p < p_cutoff & .fc <= -fc_cutoff ~ "Down",
        TRUE ~ "NS"
      )
    )
  
  df$sig_dir <- factor(df$sig_dir, levels = c("Up", "Down", "NS"))
  
  p <- ggplot(df, aes(x = .fc, y = .neglog10p)) +
    geom_point(
      data = df %>% filter(sig_dir == "NS"),
      aes(color = sig_dir),
      size = point_size,
      alpha = point_alpha_ns
    ) +
    geom_point(
      data = df %>% filter(sig_dir != "NS"),
      aes(color = sig_dir),
      size = point_size,
      alpha = point_alpha_sig
    ) +
    scale_color_manual(
      values = c(Up = color_up, Down = color_down, NS = color_ns),
      drop = FALSE
    ) +
    labs(
      title = title,
      subtitle = subtitle,
      x = xlab,
      y = ylab,
      color = NULL
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(size = title_size, face = "bold"),
      plot.subtitle = element_text(size = subtitle_size),
      axis.title = element_text(size = axis_title_size, face = "bold"),
      axis.text = element_text(size = axis_text_size, color = "black"),
      panel.border = element_rect(color = "black", fill = NA, linewidth = border_size),
      axis.line = element_line(color = "black", linewidth = axis_line_size),
      panel.grid = element_blank(),
      legend.position = if (show_legend) "right" else "none",
      plot.margin = margin(10, 30, 10, 10)
    )
  
  if (show_cutoff_lines) {
    p <- p +
      geom_vline(
        xintercept = c(-fc_cutoff, fc_cutoff),
        color = cutoff_line_color,
        linetype = cutoff_line_type,
        linewidth = cutoff_line_size
      ) +
      geom_hline(
        yintercept = -log10(p_cutoff),
        color = cutoff_line_color,
        linetype = cutoff_line_type,
        linewidth = cutoff_line_size
      )
  }
  
  if (!is.null(x_breaks)) {
    p <- p + scale_x_continuous(
      breaks = x_breaks,
      limits = x_limits,
      expand = expansion(mult = c(0.02, 0.08))
    )
  }
  
  if (!is.null(x_limits) || !is.null(y_limits)) {
    p <- p + coord_cartesian(xlim = x_limits, ylim = y_limits, clip = "off")
  }
  
  make_label_df <- function(gene_set, lab_color) {
    if (is.null(gene_set) || length(gene_set) == 0) return(NULL)
    
    out <- df %>% filter(.gene %in% gene_set)
    
    if (label_only_sig) {
      out <- out %>% filter(sig_dir != "NS")
    }
    
    if (nrow(out) == 0) return(NULL)
    
    out %>%
      mutate(
        label_color = lab_color,
        label_alpha = ifelse(sig_dir == "NS", label_alpha_nonsig, label_alpha_sig)
      )
  }
  
  label_df <- bind_rows(
    make_label_df(select_genes_1, label_color_1),
    make_label_df(select_genes_2, label_color_2)
  )
  
  if (nrow(label_df) == 0) {
    label_df <- df %>%
      filter(sig_dir != "NS") %>%
      mutate(rank_score = abs(.fc) * .neglog10p) %>%
      arrange(desc(rank_score)) %>%
      slice_head(n = auto_label_n) %>%
      mutate(
        label_color = "black",
        label_alpha = label_alpha_sig
      )
  } else {
    label_df <- label_df %>%
      distinct(.gene, .keep_all = TRUE)
  }
  
  if (nrow(label_df) > 0) {
    if (label_box) {
      p <- p +
        geom_label_repel(
          data = label_df,
          aes(label = .gene, alpha = label_alpha),
          color = label_df$label_color,
          fill = label_box_fill,
          label.colour = label_box_color,
          label.size = label_box_label_size,
          label.r = label_r,
          label.padding = label_padding,
          box.padding = box_padding,
          point.padding = point_padding,
          size = label_size,
          fontface = label_fontface,
          max.overlaps = max_overlaps,
          force = force,
          force_pull = force_pull,
          max.time = max_time,
          max.iter = max_iter,
          min.segment.length = min_segment_length,
          segment.color = if (show_connectors) connector_color else NA,
          segment.size = connector_size,
          segment.alpha = connector_alpha,
          segment.linetype = connector_linetype,
          seed = seed,
          show.legend = FALSE
        )
    } else {
      p <- p +
        geom_text_repel(
          data = label_df,
          aes(label = .gene, alpha = label_alpha),
          color = label_df$label_color,
          box.padding = box_padding,
          point.padding = point_padding,
          size = label_size,
          fontface = label_fontface,
          max.overlaps = max_overlaps,
          force = force,
          force_pull = force_pull,
          max.time = max_time,
          max.iter = max_iter,
          min.segment.length = min_segment_length,
          segment.color = if (show_connectors) connector_color else NA,
          segment.size = connector_size,
          segment.alpha = connector_alpha,
          segment.linetype = connector_linetype,
          seed = seed,
          show.legend = FALSE
        )
    }
  }
  
  p + scale_alpha_identity()
}


markers$gene <- rownames(markers)
#VD1
genes_set1 <- c("TRGC2","TRDV1","FCGR3A","CX3CR1","GNLY","PRF1","GZMH","ZEB2","IKZF2","TBX21","KIR2DL1","KLRD1","TRDC","XCL2","FGFBP2","CD8A","FCRL6", "KLRC2", "ITGAD", "XCL1", "NCR1", "KLRF1")
#VD2
genes_set2 <- c("DPP4","CXCR6","ANXA3","CD69","TCF4","GZMK","ANXA1","IL7R","LTB","THEMIS","CAMK4","DUSP16","TNFAIP3", "CD28", "TRGV9", "RORA", "RUNX2", "ANK3","IL17RA", "CD44")

combo_set<-c("DPP4","CXCR6","CD44","TRGC2","TRDV1","FCGR3A","CX3CR1","GNLY","PRF1","GZMH","ZEB2","IKZF2","TBX21","KIR2DL1","KLRD1","TRDC","XCL2","FGFBP2","CD8A","FCRL6", "KLRC2", "ITGAD", "XCL1", "NCR1", "KLRF1", "CD69","TCF4","GZMK","ANXA1","IL7R","LTB","THEMIS","CAMK4","DUSP16","TNFAIP3", "CD28", "TRGV9", "RORA", "RUNX2", "ANK3","IL17RA")
p <- plot_volcano_pub(
  markers = markers,
  select_genes_1 = genes_set1,
  select_genes_2 = genes_set2,
  label_color_1 = "mediumpurple4",
  label_color_2 = "chartreuse4",
  fc_cutoff = 0.5,
  p_cutoff = 0.01,
  color_up = "dodgerblue3",
  color_down = "tomato2",
  color_ns = "gray85",
  show_cutoff_lines = FALSE,
  label_box = TRUE,
  label_box_fill = "white",
  label_box_color = "black",
  label_box_label_size = 0.2,
  show_connectors = TRUE,
  connector_color = "black",
  connector_size = 0.9,
  label_only_sig = TRUE,
  label_alpha_sig = 1,
  label_alpha_nonsig = 0.7,
  label_size = 6,
  box_padding = 0.6,
  point_padding = 0.45,
  force = 3,
  force_pull = 0.3,
  max_overlaps = Inf,
  max_time = 10,
  max_iter = 50000,
  title = NULL,
  subtitle = NULL,
  x_breaks = seq(-10, 10, by = 5),
  x_limits = c(-10, 10),
  y_limits = c(0,200),
  label_fontface = "plain",
  axis_title_size = 22,
  axis_text_size = 18
)

print(p)

#gprofiler----
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(gprofiler2)

filter_genes_for_ora <- function(
    markers,
    direction = c("up", "down"),
    gene_col = "gene",
    padj_col = "p_val_adj",
    logfc_col = "avg_log2FC",
    pct1_col = "pct.1",
    pct2_col = "pct.2",
    padj_max = 0.05,
    logfc_min = 0.25,
    min_pct_high = 0.10,
    only_positive_fc_threshold = TRUE,
    drop_na_genes = TRUE,
    return = c("table", "genes")
) {
  direction <- match.arg(direction)
  return <- match.arg(return)
  
  df <- markers %>%
    as.data.frame()
  
  needed_cols <- c(gene_col, padj_col, logfc_col, pct1_col, pct2_col)
  missing_cols <- setdiff(needed_cols, colnames(df))
  if (length(missing_cols) > 0) {
    stop("Missing required columns in markers: ", paste(missing_cols, collapse = ", "))
  }
  
  df <- df %>%
    mutate(
      .gene  = .data[[gene_col]],
      .padj  = as.numeric(.data[[padj_col]]),
      .logfc = as.numeric(.data[[logfc_col]]),
      .pct1  = as.numeric(.data[[pct1_col]]),
      .pct2  = as.numeric(.data[[pct2_col]])
    )
  
  if (drop_na_genes) {
    df <- df %>% filter(!is.na(.gene), .gene != "")
  }
  
  df <- df %>% filter(!is.na(.padj), !is.na(.logfc), !is.na(.pct1), !is.na(.pct2))
  
  if (direction == "up") {
    df_filt <- df %>%
      filter(
        .padj <= padj_max,
        .logfc >= logfc_min,
        .pct1 >= min_pct_high
      ) %>%
      arrange(.padj, desc(.logfc), desc(.pct1))
  } else {
    df_filt <- df %>%
      filter(
        .padj <= padj_max,
        .logfc <= -logfc_min,
        .pct2 >= min_pct_high
      ) %>%
      arrange(.padj, .logfc, desc(.pct2))
  }
  
  df_filt <- df_filt %>%
    distinct(.gene, .keep_all = TRUE)
  
  if (return == "genes") {
    return(df_filt$.gene)
  } else {
    df_filt %>%
      transmute(
        gene = .gene,
        p_val_adj = .padj,
        avg_log2FC = .logfc,
        pct.1 = .pct1,
        pct.2 = .pct2,
        high_group_pct = ifelse(direction == "up", .pct1, .pct2),
        direction = direction
      ) %>%
      return()
  }
}


run_gost_ora <- function(
    genes,
    organism = "hsapiens",
    sources = c("GO:MF", "GO:CC", "GO:BP", "KEGG", "REAC"),
    user_threshold = 0.05,
    correction_method = "g_SCS",
    domain_scope = "annotated",
    ordered_query = FALSE,
    multi_query = FALSE,
    significant = TRUE,
    exclude_iea = FALSE,
    measure_underrepresentation = FALSE,
    evcodes = FALSE,
    numeric_ns = NULL,
    custom_bg = NULL
) {
  genes <- unique(genes)
  genes <- genes[!is.na(genes) & genes != ""]
  
  if (length(genes) == 0) {
    stop("No genes supplied to run_gost_ora().")
  }
  
  gost(
    query = genes,
    organism = organism,
    ordered_query = ordered_query,
    multi_query = multi_query,
    significant = significant,
    exclude_iea = exclude_iea,
    measure_underrepresentation = measure_underrepresentation,
    evcodes = evcodes,
    user_threshold = user_threshold,
    correction_method = correction_method,
    domain_scope = domain_scope,
    custom_bg = custom_bg,
    numeric_ns = numeric_ns,
    sources = sources
  )
}


tidy_gost_results <- function(gost_res) {
  if (is.null(gost_res) || is.null(gost_res$result) || nrow(gost_res$result) == 0) {
    return(tibble())
  }
  
  gost_res$result %>%
    as_tibble() %>%
    mutate(
      pathway_number = row_number(),
      FDR = p_value,
      neg_log10_FDR = -log10(FDR)
    ) %>%
    select(
      pathway_number,
      source,
      term_id,
      term_name,
      FDR,
      neg_log10_FDR,
      intersection_size,
      term_size,
      query_size,
      effective_domain_size,
      precision,
      recall,
      everything()
    )
}


print_gost_summary <- function(
    gost_tbl,
    n = 30,
    source_order = c("GO:BP", "GO:CC", "GO:MF", "KEGG", "REAC", "WP")
) {
  if (nrow(gost_tbl) == 0) {
    message("No enriched terms found.")
    return(invisible(gost_tbl))
  }
  
  source_order <- unique(as.character(source_order))
  source_present <- unique(as.character(gost_tbl$source))
  source_levels <- unique(c(source_order, setdiff(source_present, source_order)))
  
  out <- gost_tbl %>%
    mutate(source = factor(as.character(source), levels = source_levels)) %>%
    arrange(source, FDR, desc(intersection_size)) %>%
    select(
      pathway_number,
      source,
      term_name,
      term_id,
      intersection_size,
      term_size,
      FDR,
      neg_log10_FDR
    ) %>%
    slice_head(n = n)
  
  print(out, n = n)
  invisible(out)
}


plot_gost_terms <- function(
    gost_tbl,
    pathways,
    title = "ORA enrichment",
    subtitle = NULL,
    xlab = expression(-log[10](FDR)),
    base_size = 13,
    sort_desc = TRUE,
    
    bar_color = "grey85",
    bar_outline = "grey60",
    bar_linewidth = 0.3,
    bar_height = 0.82,
    
    text_color = "black",
    text_size = 5,
    text_x = 0.06,
    text_fontface = "plain",
    
    title_size = 18,
    title_face = "bold",
    
    axis_title_size = 14,
    axis_text_size = 12,
    
    axis_line_x_color = "black",
    axis_line_y_color = "black",
    axis_line_x_width = 0.8,
    axis_line_y_width = 0.8,
    
    expand_mult = c(0, 0.05),
    x_limits = NULL,
    x_breaks = waiver(),
    x_decimals = NULL
) {
  if (nrow(gost_tbl) == 0) {
    stop("gost_tbl is empty.")
  }
  
  if (is.numeric(pathways)) {
    plot_df <- gost_tbl %>%
      dplyr::filter(pathway_number %in% pathways)
  } else {
    plot_df <- gost_tbl %>%
      dplyr::filter(term_name %in% pathways | term_id %in% pathways)
  }
  
  if (nrow(plot_df) == 0) {
    stop("None of the selected pathways were found.")
  }
  
  source_abbrev <- c(
    "GO:BP" = "GO",
    "GO:CC" = "GO",
    "GO:MF" = "GO",
    "KEGG"  = "K",
    "REAC"  = "R",
    "WP"    = "WP",
    "HPA"   = "H"
  )
  
  plot_df <- plot_df %>%
    dplyr::mutate(
      source_short = dplyr::recode(source, !!!source_abbrev, .default = as.character(source)),
      term_label = paste0(term_name, " (", source_short, ")")
    )
  
  if (sort_desc) {
    plot_df <- plot_df %>%
      dplyr::arrange(dplyr::desc(neg_log10_FDR)) %>%
      dplyr::mutate(term_label = factor(term_label, levels = rev(term_label)))
  } else {
    plot_df <- plot_df %>%
      dplyr::arrange(neg_log10_FDR) %>%
      dplyr::mutate(term_label = factor(term_label, levels = rev(term_label)))
  }
  
  if (is.null(x_decimals)) {
    x_labels_fun <- waiver()
  } else {
    x_labels_fun <- \(x) formatC(x, format = "f", digits = x_decimals)
  }
  
  p <- ggplot(plot_df, aes(x = neg_log10_FDR, y = term_label)) +
    geom_col(
      width = bar_height,
      fill = bar_color,
      color = bar_outline,
      linewidth = bar_linewidth
    ) +
    geom_text(
      aes(x = text_x, label = term_label),
      hjust = 0,
      vjust = 0.5,
      color = text_color,
      size = text_size,
      fontface = text_fontface,
      lineheight = 0.95
    ) +
    labs(
      title = title,
      subtitle = subtitle,
      x = xlab,
      y = NULL
    ) +
    coord_cartesian(clip = "off", xlim = x_limits) +
    theme_classic(base_size = base_size) +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.line.x = element_line(color = axis_line_x_color, linewidth = axis_line_x_width),
      axis.line.y = element_line(color = axis_line_y_color, linewidth = axis_line_y_width),
      
      axis.title.x = element_text(size = axis_title_size, color = "black"),
      axis.text.x = element_text(size = axis_text_size, color = "black"),
      
      plot.title = element_text(size = title_size, face = title_face, hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5),
      
      plot.margin = margin(10, 20, 10, 20)
    )
  
  if (is.null(x_limits)) {
    p <- p + scale_x_continuous(
      breaks = x_breaks,
      labels = x_labels_fun,
      expand = expansion(mult = expand_mult)
    )
  } else {
    p <- p + scale_x_continuous(
      breaks = x_breaks,
      labels = x_labels_fun,
      expand = expansion(mult = c(0, 0))
    )
  }
  
  p
}

ora_up_tbl <- filter_genes_for_ora(
  markers = markers,
  direction = "up",
  gene_col = "gene",
  padj_col = "p_val_adj",
  logfc_col = "avg_log2FC",
  pct1_col = "pct.1",
  pct2_col = "pct.2",
  padj_max = 0.01,
  logfc_min = 1,
  min_pct_high = 0.20,
  return = "table"
)

ora_up_genes <- ora_up_tbl$gene
length(ora_up_genes)
head(ora_up_tbl)

ora_down_tbl <- filter_genes_for_ora(
  markers = markers,
  direction = "down",
  gene_col = "gene",
  padj_col = "p_val_adj",
  logfc_col = "avg_log2FC",
  pct1_col = "pct.1",
  pct2_col = "pct.2",
  padj_max = 0.01,
  logfc_min = 1,
  min_pct_high = 0.20,
  return = "table"
)

ora_down_genes <- ora_down_tbl$gene
length(ora_down_genes)
head(ora_down_tbl)

gost_up <- run_gost_ora(
  genes = ora_up_genes,
  organism = "hsapiens",
  sources = c("GO:MF", "GO:CC", "GO:BP", "KEGG", "REAC", "WP"),
  user_threshold = 0.05,
  correction_method = "g_SCS",
  domain_scope = "annotated",
  ordered_query = FALSE,
  multi_query = FALSE,
  significant = TRUE,
  exclude_iea = FALSE,
  measure_underrepresentation = FALSE,
  evcodes = FALSE,
  numeric_ns = NULL
)

gost_down <- run_gost_ora(
  genes = ora_down_genes,
  organism = "hsapiens",
  sources = c("GO:MF", "GO:CC", "GO:BP", "KEGG", "REAC", "WP"),
  user_threshold = 0.05,
  correction_method = "g_SCS",
  domain_scope = "annotated",
  ordered_query = FALSE,
  multi_query = FALSE,
  significant = TRUE,
  exclude_iea = FALSE,
  measure_underrepresentation = FALSE,
  evcodes = FALSE,
  numeric_ns = NULL
)

gost_up_tbl <- tidy_gost_results(gost_up)
gost_down_tbl <- tidy_gost_results(gost_down)

print_gost_summary(gost_up_tbl, n = 40)
print_gost_summary(gost_down_tbl, n = 40)

# writexl::write_xlsx(
#   list(
#     ora_up_genes = ora_up_tbl,
#     ora_down_genes = ora_down_tbl,
#     gost_up = gost_up_tbl,
#     gost_down = gost_down_tbl
#   ),
#   path = "gost_ora_all_gdts_cond.xlsx"
# )

plot_gost_terms(
  gost_tbl = gost_up_tbl,
  pathways = c(4, 186, 189, 194, 196, 202, 209),
  title = "Resistant",
  bar_color = "steelblue1",
  bar_outline = "dodgerblue4",
  text_size = 7,
  text_x = 0.3,
  axis_title_size = 22,
  axis_text_size = 18,
  title_size = 28,
  bar_linewidth = 1,
  x_limits = c(0,10)
)

plot_gost_terms(
  gost_tbl = gost_down_tbl,
  pathways = c(4, 7, 12, 19, 61),
  title = "Untreated",
  bar_color = "tomato1",
  bar_outline = "tomato4",
  text_size = 7,
  text_x = 0.3,
  axis_title_size = 22,
  axis_text_size = 18,
  title_size = 28,
  bar_linewidth = 1,
  x_decimals = 1
  
)



#GSEA----
msigdbr_species()
print(msigdbr_collections(), n =26)

kegg.medicus <- msigdbr(
  species    = "Homo sapiens",
  category   = "C2",
  subcategory= "CP:KEGG_MEDICUS"
)

kegg.legacy <- msigdbr(
  species    = "Homo sapiens",
  category   = "C2",
  subcategory= "CP:KEGG_LEGACY"
)
reactome<-msigdbr(
  species    = "Homo sapiens",
  category   = "C2",
  subcategory= "CP:REACTOME"
)


pathway_list_react = split(x = reactome$gene_symbol, f = reactome$gs_name)
pathways.hallmark <- msigdbr(species = "Homo sapiens", category = "H")
#pathway_list_immune <- split(x = pathways.immune$gene_symbol, f = pathways.immune$gs_name)
pathway_list_hallmark = split(x = pathways.hallmark$gene_symbol, f = pathways.hallmark$gs_name)
pathway_list_kegg_med = split(x = kegg.medicus$gene_symbol, f = kegg.medicus$gs_name)
pathway_list_kegg_leg = split(x = kegg.legacy$gene_symbol, f = kegg.legacy$gs_name)
pathway_list <- c(pathway_list_kegg_leg, pathway_list_hallmark)


library(dplyr)
library(tibble)

markers2 = FindMarkers(gdsubset, ident.1 = "resistant", ident.2 = "untreated", 
                       group.by = "condition", min.pct = 0.001, logfc.threshold = 0)
markers2$gene = rownames(markers2)
res <- markers2 %>% dplyr::select(gene, avg_log2FC) 
ranks <- deframe(res)

GSEA_result <-  fgsea(pathways = pathway_list, stats = ranks)



sig_pathways <- GSEA_result[GSEA_result$padj<0.1, ]
sig_pathways <- sig_pathways[order(sig_pathways$NES), ]
sig_pathways$pathway <- factor(
  sig_pathways$pathway, 
  levels = sig_pathways$pathway[order(sig_pathways$NES)]
)
top_pathways <- head(sig_pathways[sig_pathways$NES<0, ], 50)
bottom_pathways <- tail(sig_pathways[sig_pathways$NES>0, ], 50)
sig_pathways_filter <- rbind(top_pathways, bottom_pathways)
sig_pathways_filter$pathway <- factor(
  sig_pathways_filter$pathway,
  levels = sig_pathways_filter$pathway[order(sig_pathways_filter$NES)]
)
p1 <- ggplot(sig_pathways_filter, aes(x = pathway, y = NES, fill = -log10(padj))) +
  geom_bar(stat = "identity", position = "identity") +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 3) +
  labs(x = "", y = "Normalized Enrichment Score", fill = "FDR") +  # Update axis and legend labels
  theme_bw() + coord_flip()+
  theme(
    axis.text.x = element_text(vjust = 0, hjust = 0.5, family = "sans", size = 8, color = "black"),
    axis.text.y = element_text(family = "sans", size = 8, color = "black"),
    axis.line = element_blank(),
    plot.margin = margin(l = 5, r = 5, b = 5, t = 5, unit = "pt"), 
    panel.grid = element_blank())
p1

#module scoring_dotplot----


Vd1_dot<-DotPlot(gdsubset, features = genes_set1, cluster.idents = T, 
        group.by = "condition", scale = F) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + coord_flip() +
#  ggtitle("Vd1_signature_genes") +
  theme(plot.title = element_text(hjust = 0.5))

Vd2_dot<-DotPlot(gdsubset, features = genes_set2, cluster.idents = T, 
                 group.by = "condition", scale = F) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + coord_flip() +
#  ggtitle("Vd2_signature_genes") +
  theme(plot.title = element_text(hjust = 0.5))

combo_dot<-DotPlot(gdsubset, features = combo_set, cluster.idents = T, 
                 group.by = "condition", scale = F) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + coord_flip() +
  #  ggtitle("combo_signature_genes") +
  theme(plot.title = element_text(hjust = 0.5))
combo_dot


library(Seurat)
library(ggplot2)
library(dplyr)
library(ggtext)
library(grid)
library(cowplot)

dotplot_combo_seurat_style <- function(
    object,
    genes_set1,
    genes_set2,
    group.by = "condition",
    group.order = NULL,
    group.labels = NULL,
    cluster.idents = TRUE,
    scale_expr = FALSE,
    flip = TRUE,
    
    set1_name = "Vd1 signature",
    set2_name = "Vd2 signature",
    set1_color = "mediumpurple4",
    set2_color = "chartreuse4",
    
    expr_colors = NULL,
    expr_limits = NULL,
    
    x_title = NULL,
    y_title = NULL,
    
    gene_fontface = "plain",
    gene_text_size = 12,
    gene_text_color_default = "black",
    
    condition_text_size = 12,
    condition_text_color = "black",
    condition_angle = 45,
    
    axis_title_size = 14,
    axis_title_face = "plain",
    axis_title_color = "black",
    
    axis_line_color = "black",
    axis_line_width = 0.8,
    
    tick_color = "black",
    tick_width = 0.8,
    tick_length = 0.2,
    
    panel_background = "white",
    plot_background = "white",
    panel_border = FALSE,
    panel_border_color = "black",
    panel_border_width = 0.8,
    
    legend_title_size = 12,
    legend_text_size = 11,
    
    plot_title = NULL,
    plot_title_size = 14,
    plot_title_face = "bold",
    
    add_gene_set_legend = TRUE,
    gene_legend_title = "Gene set",
    
    expr_legend_title = "Average Expression",
    pct_legend_title = "Percent Expressed",
    
    legend_rel_widths = c(1, 0.38),
    legend_stack_heights = c(1, 0.42),
    
    manual_gene_legend_text_size = 11,
    manual_gene_legend_title_size = 12,
    manual_gene_legend_point_size = 5
) {
  combo_set <- unique(c(genes_set1, genes_set2))
  
  gene_df <- tibble(
    feature = combo_set,
    gene_set = case_when(
      feature %in% genes_set1 ~ set1_name,
      feature %in% genes_set2 ~ set2_name,
      TRUE ~ "Other"
    )
  ) %>%
    mutate(
      gene_label = case_when(
        gene_set == set1_name & gene_fontface == "italic" ~
          paste0("<span style='color:", set1_color, ";'><i>", feature, "</i></span>"),
        gene_set == set2_name & gene_fontface == "italic" ~
          paste0("<span style='color:", set2_color, ";'><i>", feature, "</i></span>"),
        gene_set == set1_name & gene_fontface == "plain" ~
          paste0("<span style='color:", set1_color, ";'>", feature, "</span>"),
        gene_set == set2_name & gene_fontface == "plain" ~
          paste0("<span style='color:", set2_color, ";'>", feature, "</span>"),
        TRUE ~ paste0("<span style='color:", gene_text_color_default, ";'>", feature, "</span>")
      )
    )
  
  label_map <- gene_df$gene_label
  names(label_map) <- gene_df$feature
  
  if (!is.null(group.order)) {
    object[[group.by]] <- factor(
      object[[group.by]][, 1],
      levels = group.order
    )
  }
  
  if (is.null(group.labels)) {
    group_label_map <- NULL
  } else {
    if (is.null(names(group.labels))) {
      if (!is.null(group.order) && length(group.labels) == length(group.order)) {
        group_label_map <- group.labels
        names(group_label_map) <- group.order
      } else {
        stop("group.labels must be a named vector, or match the length of group.order.")
      }
    } else {
      group_label_map <- group.labels
    }
  }
  
  p <- DotPlot(
    object = object,
    features = combo_set,
    cluster.idents = cluster.idents,
    group.by = group.by,
    scale = scale_expr
  )
  
  p <- p + scale_x_discrete(labels = label_map)
  
  if (!is.null(group_label_map)) {
    p <- p + scale_y_discrete(labels = group_label_map)
  }
  
  if (!is.null(expr_colors)) {
    p <- p + scale_color_gradientn(
      colours = expr_colors,
      limits = expr_limits,
      oob = scales::squish,
      name = expr_legend_title
    )
  }
  
  p <- p + guides(
    size = guide_legend(title = pct_legend_title),
    color = guide_colorbar(title = expr_legend_title)
  )
  
  if (is.null(x_title)) x_title <- if (flip) group.by else "Features"
  if (is.null(y_title)) y_title <- if (flip) "Features" else group.by
  
  internal_x_title <- if (flip) y_title else x_title
  internal_y_title <- if (flip) x_title else y_title
  
  p <- p +
    labs(
      title = plot_title,
      x = internal_x_title,
      y = internal_y_title
    ) +
    theme(
      legend.position = "right",
      legend.title = element_text(size = legend_title_size, face = "bold"),
      legend.text = element_text(size = legend_text_size),
      
      axis.title.x = element_text(
        size = axis_title_size,
        face = axis_title_face,
        color = axis_title_color
      ),
      axis.title.y = element_text(
        size = axis_title_size,
        face = axis_title_face,
        color = axis_title_color
      ),
      
      axis.line = element_line(
        color = axis_line_color,
        linewidth = axis_line_width
      ),
      axis.ticks = element_line(
        color = tick_color,
        linewidth = tick_width
      ),
      axis.ticks.length = unit(tick_length, "cm"),
      
      panel.background = element_rect(fill = panel_background, color = NA),
      plot.background = element_rect(fill = plot_background, color = NA),
      
      panel.border = if (panel_border) {
        element_rect(fill = NA, color = panel_border_color, linewidth = panel_border_width)
      } else {
        element_blank()
      },
      
      plot.title = element_text(
        size = plot_title_size,
        face = plot_title_face,
        hjust = 0.5
      )
    )
  
  if (flip) {
    p <- p +
      coord_flip() +
      theme(
        axis.text.y = ggtext::element_markdown(size = gene_text_size),
        axis.text.x = element_text(
          size = condition_text_size,
          color = condition_text_color,
          angle = condition_angle,
          hjust = 0.5,
          vjust = 1
        )
      )
  } else {
    p <- p +
      theme(
        axis.text.x = ggtext::element_markdown(
          size = gene_text_size,
          angle = condition_angle,
          hjust = 1,
          vjust = 1
        ),
        axis.text.y = element_text(
          size = condition_text_size,
          color = condition_text_color
        )
      )
  }
  
  seurat_legend <- cowplot::get_legend(p)
  p_no_legend <- p + theme(legend.position = "none")
  
  if (!add_gene_set_legend) {
    return(
      cowplot::plot_grid(
        p_no_legend,
        cowplot::ggdraw(seurat_legend),
        nrow = 1,
        rel_widths = legend_rel_widths
      )
    )
  }
  
  seurat_legend_plot <- cowplot::ggdraw(seurat_legend) +
    theme(
      plot.margin = margin(t = 0, r = 0, b = -30, l = 0)
    )
  
  gene_legend_panel <- ggplot() +
    xlim(0, 1) + ylim(0.40, 1.00) +
    annotate(
      "text", x = 0.00, y = 0.97, label = gene_legend_title,
      hjust = 0.07,
      fontface = "bold",
      size = manual_gene_legend_title_size / 3
    ) +
    annotate(
      "point", x = 0.00, y = 0.74, shape = 15,
      size = manual_gene_legend_point_size, colour = set1_color
    ) +
    annotate(
      "text", x = 0.12, y = 0.74, label = set1_name,
      hjust = 0, vjust = 0.5, size = manual_gene_legend_text_size / 3
    ) +
    annotate(
      "point", x = 0.00, y = 0.56, shape = 15,
      size = manual_gene_legend_point_size, colour = set2_color
    ) +
    annotate(
      "text", x = 0.12, y = 0.56, label = set2_name,
      hjust = 0, vjust = 0.5, size = manual_gene_legend_text_size / 3
    ) +
    coord_cartesian(clip = "off") +
    theme_void() +
    theme(
      plot.margin = margin(t = -80, r = 0, b = 110, l = 0)
    )
  
  right_legends <- cowplot::plot_grid(
    seurat_legend_plot,
    gene_legend_panel,
    ncol = 1,
    rel_heights = legend_stack_heights,
    align = "v"
  )
  
  cowplot::plot_grid(
    p_no_legend,
    right_legends,
    nrow = 1,
    rel_widths = legend_rel_widths,
    align = "h"
  )
}

combo_dot <- dotplot_combo_seurat_style(
  object = gdsubset,
  genes_set1 = genes_set1,
  genes_set2 = genes_set2,
  group.by = "condition",
  group.order = c("untreated", "resistant"),
  group.labels = c(
    untreated = "UT",
    resistant = "Res"
  ),
  cluster.idents = TRUE,
  scale_expr = FALSE,
  flip = TRUE,
  
  set1_name = "Vd1",
  set2_name = "Vd2",
  set1_color = "mediumpurple4",
  set2_color = "chartreuse4",
  
  expr_colors = c("grey95", "blue", "darkblue"),
  expr_limits = NULL,
  
  x_title = "",
  y_title = "Genes",
  
  gene_text_size = 16,
  condition_text_size = 22,
  condition_angle = 0,
  
  axis_title_size = 24,
  axis_title_face = "bold",
  
  axis_line_color = "black",
  axis_line_width = 0.9,
  tick_color = "black",
  tick_width = 0.9,
  tick_length = 0.22,
  
  legend_title_size = 16,
  legend_text_size = 14,
  
  expr_legend_title = "Average Expression",
  pct_legend_title = "Percent Expressed",
  gene_legend_title = "Signature Genes",
  
  add_gene_set_legend = TRUE,
  legend_rel_widths = c(0.8, 0.4),
  legend_stack_heights = c(1, 0.2),
  manual_gene_legend_text_size = 16,
  manual_gene_legend_title_size = 18,
  manual_gene_legend_point_size = 5
)

combo_dot

#reference transfer, broad----

DimPlot(object = reference, reduction = "umap.rpca", group.by = "broad_class",
        label = TRUE, label.size = 3, repel = TRUE, raster = F, pt.size = 0.15, label.box = T) + NoLegend()

# Reference transfer from reference -> gdsubset ----
library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)

DefaultAssay(reference) <- "RNA"
DefaultAssay(gdsubset)  <- "RNA"

reference <- NormalizeData(reference, verbose = FALSE)
gdsubset  <- NormalizeData(gdsubset, verbose = FALSE)

shared_features <- intersect(rownames(reference), rownames(gdsubset))

reference <- FindVariableFeatures(reference, verbose = FALSE)
transfer_features <- intersect(VariableFeatures(reference), shared_features)

reference <- ScaleData(reference, features = transfer_features, verbose = FALSE)
reference <- RunPCA(reference, features = transfer_features, npcs = 50, verbose = FALSE)

anchors <- FindTransferAnchors(
  reference = reference,
  query = gdsubset,
  normalization.method = "LogNormalize",
  reference.assay = "RNA",
  query.assay = "RNA",
  dims = 1:50,
  features = transfer_features,
  reduction = "pcaproject"
)

pred_broad <- TransferData(
  anchorset = anchors,
  refdata = reference$broad_class,
  dims = 1:50
)

gdsubset <- AddMetaData(gdsubset, metadata = pred_broad)

gdsubset$broad_class_pred_rna  <- gdsubset$predicted.id
gdsubset$broad_class_score_rna <- gdsubset$prediction.score.max

summary(gdsubset$broad_class_score_rna)

DimPlot(
  gdsubset,
  reduction = "rpca.umap",
  group.by = "broad_class_pred_rna",
  label = TRUE
) + coord_fixed()


#prop bar plots----
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(tibble)

meta_gd <- gdsubset@meta.data %>%
  as.data.frame() %>%
  rownames_to_column("cell")

prop_broad_sample <- meta_gd %>%
  dplyr::count(orig.ident, condition, broad_class_pred_rna, name = "n_cells") %>%
  dplyr::group_by(orig.ident, condition) %>%
  dplyr::mutate(prop = n_cells / sum(n_cells)) %>%
  dplyr::ungroup()

prop_broad_condition <- meta_gd %>%
  dplyr::count(condition, broad_class_pred_rna, name = "n_cells") %>%
  dplyr::group_by(condition) %>%
  dplyr::mutate(prop = n_cells / sum(n_cells)) %>%
  dplyr::ungroup()

plot_broad_stacked <- function(
    data,
    plot_type = c("sample", "condition"),
    
    sample_col = "orig.ident",
    condition_col = "condition",
    class_col = "broad_class_pred_rna",
    prop_col = "prop",
    
    sample_order = NULL,
    condition_order = NULL,
    class_order = NULL,
    
    sample_labels = NULL,
    condition_labels = NULL,
    class_labels = NULL,
    
    fill_colors = NULL,
    drop_unused_levels = TRUE,
    
    bar_width = NULL,
    bar_color = "black",
    bar_linewidth = 0.2,
    
    facet_scales = "free_x",
    facet_space = "free_x",
    panel_spacing_x = 1,
    
    xlab = NULL,
    ylab = "Proportion",
    legend_title = "Broad class",
    plot_title = NULL,
    plot_subtitle = NULL,
    
    title_size = 16,
    title_face = "bold",
    subtitle_size = 13,
    
    axis_title_size = 14,
    axis_title_face = "plain",
    axis_title_color = "black",
    
    axis_text_size_x = 12,
    axis_text_size_y = 12,
    axis_text_color = "black",
    axis_text_angle_x = NULL,
    axis_text_hjust_x = NULL,
    axis_text_vjust_x = 1,
    
    y_labels = scales::percent_format(),
    y_limits = c(0, 1.05),
    y_breaks = seq(0, 1, 0.25),
    y_expand = c(0, 0),
    
    axis_line_color = "black",
    axis_line_linewidth = 0.8,
    axis_ticks_color = "black",
    axis_ticks_linewidth = 0.8,
    tick_length = 0.2,
    
    strip_text_size = 13,
    strip_text_face = "bold",
    strip_text_color = "black",
    strip_background_fill = NA,
    strip_background_color = NA,
    strip_background_linewidth = 0,
    strip_text_margin_b = 4,
    
    legend_position = "right",
    legend_title_size = 12,
    legend_title_face = "bold",
    legend_text_size = 11,
    
    show_legend = TRUE,
    
    show_panel_grid_major_y = FALSE,
    panel_grid_major_y_color = "grey85",
    panel_grid_major_y_linewidth = 0.3,
    
    panel_background_fill = "white",
    plot_background_fill = "white"
) {
  plot_type <- match.arg(plot_type)
  df <- as.data.frame(data)
  
  if (!sample_col %in% colnames(df) && plot_type == "sample") {
    stop(sample_col, " not found in data.")
  }
  if (!condition_col %in% colnames(df)) {
    stop(condition_col, " not found in data.")
  }
  if (!class_col %in% colnames(df)) {
    stop(class_col, " not found in data.")
  }
  if (!prop_col %in% colnames(df)) {
    stop(prop_col, " not found in data.")
  }
  
  if (!is.null(sample_order) && plot_type == "sample") {
    df[[sample_col]] <- factor(df[[sample_col]], levels = sample_order)
  }
  if (!is.null(condition_order)) {
    df[[condition_col]] <- factor(df[[condition_col]], levels = condition_order)
  }
  if (!is.null(class_order)) {
    df[[class_col]] <- factor(df[[class_col]], levels = class_order)
  }
  
  make_label_map <- function(labels, order_vec = NULL, what = "labels") {
    if (is.null(labels)) return(NULL)
    
    if (is.null(names(labels))) {
      if (!is.null(order_vec) && length(labels) == length(order_vec)) {
        out <- labels
        names(out) <- order_vec
        return(out)
      } else {
        stop(what, " must be a named vector, or have the same length as its corresponding order vector.")
      }
    }
    
    labels
  }
  
  sample_label_map <- make_label_map(sample_labels, sample_order, "sample_labels")
  condition_label_map <- make_label_map(condition_labels, condition_order, "condition_labels")
  class_label_map <- make_label_map(class_labels, class_order, "class_labels")
  
  if (plot_type == "sample") {
    if (is.null(bar_width)) bar_width <- 0.8
    if (is.null(xlab)) xlab <- "Sample"
    if (is.null(axis_text_angle_x)) axis_text_angle_x <- 45
    if (is.null(axis_text_hjust_x)) axis_text_hjust_x <- 1
    x_var <- sample_col
  } else {
    if (is.null(bar_width)) bar_width <- 0.7
    if (is.null(xlab)) xlab <- NULL
    if (is.null(axis_text_angle_x)) axis_text_angle_x <- 0
    if (is.null(axis_text_hjust_x)) axis_text_hjust_x <- 0.5
    x_var <- condition_col
  }
  
  p <- ggplot(
    df,
    aes_string(x = x_var, y = prop_col, fill = class_col)
  ) +
    geom_col(
      width = bar_width,
      color = bar_color,
      linewidth = bar_linewidth
    ) +
    scale_y_continuous(
      labels = y_labels,
      limits = y_limits,
      breaks = y_breaks,
      expand = y_expand
    ) +
    labs(
      title = plot_title,
      subtitle = plot_subtitle,
      x = xlab,
      y = ylab,
      fill = legend_title
    ) +
    theme_classic() +
    theme(
      plot.title = element_text(size = title_size, face = title_face),
      plot.subtitle = element_text(size = subtitle_size),
      
      axis.title.x = element_text(
        size = axis_title_size,
        face = axis_title_face,
        color = axis_title_color
      ),
      axis.title.y = element_text(
        size = axis_title_size,
        face = axis_title_face,
        color = axis_title_color
      ),
      
      axis.text.x = element_text(
        size = axis_text_size_x,
        angle = axis_text_angle_x,
        hjust = axis_text_hjust_x,
        vjust = axis_text_vjust_x,
        color = axis_text_color
      ),
      axis.text.y = element_text(
        size = axis_text_size_y,
        color = axis_text_color
      ),
      
      axis.line = element_line(
        color = axis_line_color,
        linewidth = axis_line_linewidth
      ),
      axis.ticks = element_line(
        color = axis_ticks_color,
        linewidth = axis_ticks_linewidth
      ),
      axis.ticks.length = unit(tick_length, "cm"),
      
      legend.position = if (show_legend) legend_position else "none",
      legend.title = element_text(size = legend_title_size, face = legend_title_face),
      legend.text = element_text(size = legend_text_size),
      
      panel.background = element_rect(fill = panel_background_fill, color = NA),
      plot.background = element_rect(fill = plot_background_fill, color = NA),
      
      panel.grid.major.y = if (show_panel_grid_major_y) {
        element_line(color = panel_grid_major_y_color, linewidth = panel_grid_major_y_linewidth)
      } else {
        element_blank()
      },
      panel.grid.minor = element_blank()
    )
  
  if (plot_type == "sample") {
    if (!is.null(sample_label_map)) {
      p <- p + scale_x_discrete(labels = sample_label_map)
    }
    
    p <- p +
      facet_grid(
        stats::as.formula(paste(". ~", condition_col)),
        scales = facet_scales,
        space = facet_space,
        labeller = if (!is.null(condition_label_map)) {
          labeller(.cols = as_labeller(condition_label_map))
        } else {
          "label_value"
        }
      ) +
      theme(
        strip.text = element_text(
          size = strip_text_size,
          face = strip_text_face,
          color = strip_text_color,
          margin = margin(b = strip_text_margin_b)
        ),
        strip.background = element_rect(
          fill = strip_background_fill,
          color = strip_background_color,
          linewidth = strip_background_linewidth
        ),
        panel.spacing.x = unit(panel_spacing_x, "lines")
      )
  } else {
    if (!is.null(condition_label_map)) {
      p <- p + scale_x_discrete(labels = condition_label_map)
    }
  }
  
  if (!is.null(fill_colors) && !is.null(class_label_map)) {
    p <- p + scale_fill_manual(
      values = fill_colors,
      labels = class_label_map,
      drop = drop_unused_levels
    )
  } else if (!is.null(fill_colors) && is.null(class_label_map)) {
    p <- p + scale_fill_manual(
      values = fill_colors,
      drop = drop_unused_levels
    )
  } else if (is.null(fill_colors) && !is.null(class_label_map)) {
    p <- p + scale_fill_discrete(
      labels = class_label_map,
      drop = drop_unused_levels
    )
  }
  
  p
}

sample_order <- c("HG589_T", "HG590_T", "HG605_T", "HG607_T", "HG610_T")
condition_order <- c("untreated", "resistant")
class_order <- c("VD1+", "VD2+", "VD3+", "Mixed/Unknown")

broad_cols <- c(
  "Mixed/Unknown" = "#556EA4",
  "VD1+" = "#976898",
  "VD2+" = "#ADC954",
  "VD3+" = "#A64A34"
)

p_sample <- plot_broad_stacked(
  data = prop_broad_sample,
  plot_type = "sample",
  
  sample_order = sample_order,
  condition_order = condition_order,
  class_order = class_order,
  
  condition_labels = c(
    untreated = "UT",
    resistant = "Res"
  ),
  
  fill_colors = broad_cols,
  
  bar_width = 0.8,
  bar_color = "black",
  bar_linewidth = 0.25,
  
  xlab = "Sample",
  ylab = "Cell type percentage",
  legend_title = "Cell type",
  plot_title = NULL,
  
  title_size = 16,
  axis_title_size = 18,
  axis_title_face = "bold",
  axis_text_size_x = 14,
  axis_text_size_y = 14,
  axis_text_angle_x = 45,
  axis_text_hjust_x = 1,
  
  strip_text_size = 16,
  strip_text_face = "bold",
  panel_spacing_x = 1.2,
  strip_background_fill = "white",
  strip_background_color = "black",
  strip_background_linewidth = 0.8,
  strip_text_margin_b = 2,
  
  legend_position = "right",
  legend_title_size = 16,
  legend_text_size = 14,
  
  y_breaks = seq(0, 1, 0.25),
  y_limits = c(0, 1.01),
  show_panel_grid_major_y = FALSE,
  show_legend = TRUE
)

p_sample

p_condition <- plot_broad_stacked(
  data = prop_broad_condition,
  plot_type = "condition",
  
  condition_order = condition_order,
  class_order = class_order,
  condition_labels = c(
    untreated = "UT",
    resistant = "Res"
  ),
  
  fill_colors = broad_cols,
  
  bar_width = 0.7,
  bar_color = "black",
  bar_linewidth = 0.25,
  
  xlab = "",
  ylab = "Cell type percentage",
  legend_title = "Cell type",
  plot_title = NULL,
  
  title_size = 16,
  axis_title_size = 22,
  axis_title_face = "bold",
  axis_text_size_x = 22,
  axis_text_size_y = 16,
  
  legend_position = "right",
  legend_title_size = 22,
  legend_text_size = 18,
  
  y_limits = c(0, 1),
  y_breaks = seq(0, 1, 0.25),
  show_panel_grid_major_y = FALSE,
  show_legend = TRUE,
  axis_text_angle_x = 0,
  axis_text_vjust_x = 1,
  axis_text_hjust_x = 0.5
)

p_condition
