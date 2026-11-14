# GSE8167 | gdT signature | RFS Kaplan-Meier and univariate Cox

library(GEOquery)
library(Biobase)
library(data.table)
library(readxl)
library(survival)
library(survminer)
library(ggplot2)

# Package  Version
# 1    GEOquery   2.76.0
# 2     Biobase   2.68.0
# 3       dplyr    1.2.1
# 4       tidyr    1.3.2
# 5     stringr    1.6.0
# 6  data.table 1.18.6.1
# 7      readxl    1.6.0
# 8    survival   3.8.12
# 9   survminer    0.5.2
# 10    ggplot2    4.0.3


# Expression matrix ----
eset <- getGEO("GSE8167", GSEMatrix = TRUE, getGPL = TRUE)[[1]]
expr <- exprs(eset)  # Already normalized and log2-transformed
features <- fData(eset)
features <- features[match(rownames(expr), rownames(features)), , drop = FALSE]
stopifnot("Gene Symbol" %in% colnames(features))

probe_map <- data.table(
  ProbeID = rownames(expr),
  Annotation = trimws(as.character(features[["Gene Symbol"]]))
)[!is.na(Annotation) & Annotation != ""]
probe_map[, Gene := strsplit(Annotation, "\\s*///\\s*")]
probe_map <- probe_map[, .(Gene = trimws(unlist(Gene))), by = ProbeID]
probe_map <- probe_map[!is.na(Gene) & Gene != ""]

probe_expr <- merge(
  probe_map,
  as.data.table(expr, keep.rownames = "ProbeID"),
  by = "ProbeID", sort = FALSE
)
expr_long <- melt(
  probe_expr,
  id.vars = c("ProbeID", "Gene"),
  measure.vars = colnames(expr),
  variable.name = "Sample", value.name = "Expression"
)
gene_expr <- expr_long[, .(Expression = median(Expression, na.rm = TRUE)),
                       by = .(Gene, Sample)]
gene_expr <- dcast(gene_expr, Gene ~ Sample, value.var = "Expression")
expr_mat <- as.matrix(gene_expr[, -1L, with = FALSE])
rownames(expr_mat) <- gene_expr$Gene
storage.mode(expr_mat) <- "numeric"

# gdT signature ----
sig_genes <- c("TRDV3", "TRGC2", "TRDC")
missing_genes <- setdiff(sig_genes, rownames(expr_mat))
if (length(missing_genes)) stop("Missing genes: ", paste(missing_genes, collapse = ", "))

sig_expr <- expr_mat[sig_genes, , drop = FALSE]
if (any(!is.finite(sig_expr))) stop("Non-finite signature expression values.")
sig_center <- apply(sig_expr, 1, median, na.rm = TRUE)
sig_scale <- (apply(sig_expr, 1, sd, na.rm = TRUE) +
                apply(sig_expr, 1, mad, na.rm = TRUE)) / 2
sig_scale[!is.finite(sig_scale) | sig_scale <= 0] <- 1
sig_z <- sweep(sweep(sig_expr, 1, sig_center, "-"), 1, sig_scale, "/")
gdt_score <- colMeans(sig_z)
names(gdt_score) <- colnames(expr_mat)

# RFS data ----
meta <- as.data.frame(read_excel("japan_meta_clean.xlsx"), check.names = FALSE) #metadata contained in source data file
stopifnot(all(c("sample", "RFS", "RFS_Months") %in% colnames(meta)))
meta$sample <- trimws(as.character(meta$sample))
if (anyNA(meta$sample) || any(meta$sample == "") || anyDuplicated(meta$sample)) {
  stop("Missing or duplicated sample IDs in metadata.")
}

rfs_event <- tolower(trimws(as.character(meta$RFS)))
meta$RFS_event <- NA_integer_
meta$RFS_event[rfs_event %in% c("1")] <- 1L
meta$RFS_event[rfs_event %in% c("0")] <- 0L
meta$RFS_Months <- suppressWarnings(as.numeric(meta$RFS_Months))

common <- intersect(names(gdt_score), meta$sample)
if (length(common) < 5) stop("Too few matched samples.")
survival_data <- meta[match(common, meta$sample), , drop = FALSE]
survival_data$gdT_score <- unname(gdt_score[common])
survival_data <- survival_data[
  is.finite(survival_data$gdT_score) &
    is.finite(survival_data$RFS_Months) & survival_data$RFS_Months >= 0 &
    survival_data$RFS_event %in% c(0L, 1L), , drop = FALSE
]
if (nrow(survival_data) < 5 || sum(survival_data$RFS_event) < 2) {
  stop("Insufficient RFS observations or events.")
}

# Kaplan-Meier ----
gdt_cutpoint <- median(survival_data$gdT_score)
survival_data$gdT_group <- ifelse(
  survival_data$gdT_score > gdt_cutpoint, "High",
  ifelse(survival_data$gdT_score < gdt_cutpoint, "Low", NA_character_)
)
survival_data$gdT_group <- factor(survival_data$gdT_group, levels = c("Low", "High"))
km_data <- droplevels(survival_data[!is.na(survival_data$gdT_group), , drop = FALSE])
if (nlevels(km_data$gdT_group) != 2) stop("KM requires two groups.")

km_fit <- survfit(Surv(RFS_Months, RFS_event) ~ gdT_group, data = km_data)
logrank <- survdiff(Surv(RFS_Months, RFS_event) ~ gdT_group, data = km_data)
logrank_p <- pchisq(logrank$chisq, df = 1, lower.tail = FALSE)

km_theme <- theme_bw(base_size = 20, base_family = "Helvetica") +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(linewidth = 1.8),
    panel.border = element_rect(fill = NA, linewidth = 1),
    axis.title = element_text(size = 22),
    axis.text = element_text(size = 22),
    plot.title = element_text(size = 20, face = "bold"),
    legend.title = element_blank(),
    legend.text = element_text(size = 20)
  )

km_plot <- ggsurvplot(
  km_fit, data = km_data,
  palette = c("#2C7BB6", "#D7191C"),
  legend.labs = c("Low", "High"), legend.title = "gdT score", legend = "right",
  conf.int = FALSE, risk.table = FALSE, pval = TRUE, pval.size = 6,
  censor = TRUE, censor.size = 4, censor.shape = 3,
  size = 1.5, break.time.by = 25,
  xlab = "Time (months)", ylab = "Survival Probability",
  title = "gdT score - RFS", ggtheme = km_theme
)

km_stats <- data.frame(
  n = nrow(km_data), Events = sum(km_data$RFS_event),
  n_Low = sum(km_data$gdT_group == "Low"),
  n_High = sum(km_data$gdT_group == "High"),
  Events_Low = sum(km_data$RFS_event == 1 & km_data$gdT_group == "Low"),
  Events_High = sum(km_data$RFS_event == 1 & km_data$gdT_group == "High"),
  Logrank_p = logrank_p
)

# Univariate Cox ----
cox_data <- survival_data
cox_data$gdT_per_SD <- as.numeric(scale(cox_data$gdT_score))
if (any(!is.finite(cox_data$gdT_per_SD))) stop("gdT score has zero variance.")
cox_fit <- coxph(Surv(RFS_Months, RFS_event) ~ gdT_per_SD,
                 data = cox_data, ties = "efron", x = TRUE)
cox_summary <- summary(cox_fit)
cox_results <- data.frame(
  n = nrow(cox_data), Events = sum(cox_data$RFS_event),
  HR_per_SD = unname(cox_summary$conf.int[1, "exp(coef)"]),
  Lower_95CI = unname(cox_summary$conf.int[1, "lower .95"]),
  Upper_95CI = unname(cox_summary$conf.int[1, "upper .95"]),
  Wald_p = unname(cox_summary$coefficients[1, "Pr(>|z|)"])
)

# Results ----
print(km_plot)
print(km_stats)
print(cox_results)
print(summary(cox_fit))
# View(survival_data)

