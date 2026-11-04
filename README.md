# Single-cell profiling reveals γδ T-cell heterogeneity associated with prognosis and imatinib resistance in human GIST

Analysis code accompanying the manuscript:

> Zeng S\*, Sussman JH\*, Morris M\*, Fu Y, Perez JE, Rong J, Tardy K, Kwak H, Hartlein T, Liu S, Do K, Rossi F, Zhang JQ, Xu X, DeMatteo RP. **Single-cell profiling reveals γδ T-cell heterogeneity associated with prognosis and imatinib resistance in human GIST.** *(Under revision, Nature Communications.)*
>
> \*Equal contribution

**Lab:** DeMatteo Lab (RDM lab), Department of Surgery, Perelman School of Medicine, University of Pennsylvania
**Corresponding author:** Ronald P. DeMatteo, ronald.dematteo@pennmedicine.upenn.edu

---

## Abstract

Gamma-delta (γδ) T cells are emerging effectors of anti-tumor immunity, but their role in gastrointestinal stromal tumor (GIST), the most common human sarcoma, remains poorly defined. We integrated single-cell RNA-sequencing (scRNA-seq), T cell receptor profiling, immunophenotyping, and functional assays to characterize circulating and intratumoral γδ T cells. We constructed a scRNA-seq atlas of γδ T cells from 11 patients and identified 14 transcriptomic clusters spanning heterogeneous differentiation and activation states. Intratumoral γδ T cells, and Vδ1+ cells, were enriched in immunomodulatory signatures and effector phenotypes. Vδ1+ cells predominated within tumors and their frequency within the γδ T cell compartment correlated with favorable prognostic features. Expanded γδ T cells from healthy donors and GIST tumors were cytotoxic *ex vivo* via the MICA-NKG2D axis. Compared with untreated tumors, imatinib-resistant tumors exhibited a shift toward Vδ2+ cells, altered metabolic signatures, and reduced effector-associated features. In a genetically engineered murine GIST model, PD-L1 blockade enhanced imatinib efficacy and restored tumor γδ T cell interferon and IL17 activity. Together, these findings characterize γδ T cell states in human GIST and identify associations with tumor prognosis and imatinib resistance, supporting further investigation of γδ T cells as potential immunotherapy targets.

**Keywords:** γδ T cells, Vδ1 T cells, single-cell transcriptomics, T cell receptor (TCR) repertoire, gastrointestinal stromal tumor, tumor immunity

---

## Repository overview

This repository contains the code used for the computational analyses in the paper: processing of the human γδ T cell scRNA-seq data, atlas construction and annotation, differential expression and pathway analyses, transcription factor regulon inference, TCR repertoire analysis, and analysis of the murine bulk RNA-seq data.

Scripts are numbered in the order they are intended to be run.

```
.
├── 01_Run_CellRanger_Example.sh                     # Example Cell Ranger run (FASTQ -> count matrices)
├── 02_Run_TRUST4_TCR.sh                             # TCR inference from 5' scRNA-seq with TRUST4
├── 03_Process_Integrate_GD_TCell_scRNA_Seurat.R     # QC, filtering, RPCA integration, clustering, contaminant removal
├── 04_Differential_Expression_scRNA_Figures.R       # Marker genes, cluster annotation, scRNA-seq figures
├── 05_CytoTRACE.R                                   # CytoTRACE differentiation-potential analysis
├── 06_SCENIC_Regulon_Analysis_Figures.R             # Downstream analysis and plots of pySCENIC regulon activity
├── 07_DESeq2_Pseudobulk_Comparisons.R               # Pseudobulk DESeq2 (+ fgsea) comparisons
├── 08_GLMM_Pseudobulk_Comparisons_and_Signatures.R  # limma-voom/dream mixed models (Vδ1 vs Vδ2) and gene signature scores
├── 09_Murine_GDT_RNAseq.R                           # Bulk RNA-seq of sorted murine γδ T cells (vehicle vs imatinib)
├── SCENIC/                                          # pySCENIC workflow files
└── TCR_Analysis/                                    # TCRδ repertoire (ImmunoSEQ) analysis
```

---

## Data availability

The raw and processed data are **not** stored in this repository. They are available from the following sources:

| Data | Location |
|---|---|
| Human γδ T cell scRNA-seq (raw reads, Cell Ranger outputs, and the Seurat object with metadata used in this work) | NCBI GEO: [GSE301351](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE301351) |
| Bulk RNA-seq of sorted murine γδ T cells (*Kit*^V558del/+^ mice, vehicle vs imatinib) | NCBI GEO: [GSE300283](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE300283) |
| TCRαδ sequencing of sorted human γδ T cells (ImmunoSEQ) | Zenodo: <https://zenodo.org/records/23366688> |
| Remaining data | Within the article and its supplementary materials, or from the corresponding author on reasonable request |

The survival analysis in the paper uses a previously published RNA microarray dataset of primary GIST (see the manuscript for the reference).

---

## Analysis workflow

1. **Preprocessing (`01`, `02`).** Reads from the 10x Genomics platform (NextGEM v3.1 5′ or v4 3′ assays) were processed with Cell Ranger v9.0.1 against the GRCh38-2024-A reference. For the 5′ samples, TCR identity was inferred from the scRNA-seq reads with TRUST4 v1.1.10 and used to help identify contaminating αβ T cells.
2. **Atlas construction (`03`).** Cells were filtered (nFeature_RNA, nCount_RNA, percent mitochondrial and ribosomal reads), normalized, and integrated across samples with reciprocal PCA (RPCA) in Seurat. Clusters of non-T cells (tumor, B, endothelial, myeloid), αβ T cells, and doublets were removed iteratively, and the data were re-integrated and re-clustered to produce the final atlas of 84,291 γδ T cells from 4 blood and 9 tumor samples.
3. **Annotation and differential expression (`04`).** Cluster markers (Wilcoxon rank-sum test) were used to annotate 14 clusters by predominant Vδ class (Vδ1, Vδ2, Vδ3, or mixed). Differential abundance between groups was tested with *propeller*.
4. **Differentiation state (`05`).** CytoTRACE was used to assess differentiation potential across clusters.
5. **Transcription factor regulons (`SCENIC/`, `06`).** Regulon inference was run with pySCENIC (GRNBoost2, cisTarget v9 motif collection, hg38 RefSeq r80 databases, AUCell scoring). The resulting regulon activities were analyzed and plotted in R.
6. **Pseudobulk comparisons (`07`, `08`).** Counts were aggregated per sample and analyzed with DESeq2 (PBMC vs tumor; untreated vs imatinib-resistant), followed by pre-ranked GSEA with fgsea using MSigDB Hallmark and Reactome collections and KEGG. Vδ1-predominant vs Vδ2-predominant clusters were compared with a mixed-model framework (limma-voom with dream, sample as a random effect). Gene signature scores were computed from robust z-scores.
7. **TCR repertoire (`TCR_Analysis/`).** TCRδ CDR3 repertoires from sorted PBMC and tumor γδ T cells were analyzed for clonal expansion, repertoire diversity, and V-gene usage, including repeated downsampling to equal depth.
8. **Murine RNA-seq (`09`).** Bulk RNA-seq from sorted tumor γδ T cells of *Kit*^V558del/+^ mice treated with vehicle or imatinib was analyzed with DESeq2 and fgsea (reads aligned to mm10 with STAR and quantified with featureCounts).

---

## Software and versions

Versions reported in the manuscript:

| Tool | Version |
|---|---|
| Cell Ranger | 9.0.1 (reference: GRCh38-2024-A) |
| TRUST4 | 1.1.10 |
| R | 4.5.1 |
| Seurat | 5.3.1 |
| DESeq2 | 1.50.2 |
| pySCENIC | see `SCENIC/` (cisTarget human v9 motif collection) |

Other packages used include `apeglm`, `fgsea`, `propeller`, `limma`/`variancePartition` (dream), `clusterProfiler`, `SCopeLoomR`, CytoTRACE, and `survival`. Please check the individual scripts for the full list of `library()` calls and any version-specific notes.

---

## Usage

```bash
# 1. Clone the repository
git clone https://github.com/montana-morris-penn/GdT-cell-Paper-RDM-lab.git
cd GdT-cell-Paper-RDM-lab

# 2. Download the data (see "Data availability" above)
#    - GEO GSE301351 (scRNA-seq) and GSE300283 (murine bulk RNA-seq)
#    - Zenodo record 23366688 (TCRαδ sequencing)

# 3. Run the scripts in numerical order, editing the input/output paths at the top of each script
```

To skip the preprocessing steps (`01`–`03`), you can start from the Seurat object with metadata deposited at GEO under GSE301351 and run the downstream scripts (`04` onward).

> **Note:** Several scripts contain paths to local or cluster directories. Update these to match your own environment before running.

---

## Citation

If you use this code or data, please cite:

> Zeng S, Sussman JH, Morris M, et al. Single-cell profiling reveals γδ T-cell heterogeneity associated with prognosis and imatinib resistance in human GIST. *Nature Communications* (in revision). \[DOI to be added upon publication\]

---

## Contact

For questions about the code, please open a GitHub issue or contact the authors. For correspondence about the paper, contact Ronald P. DeMatteo (ronald.dematteo@pennmedicine.upenn.edu), Department of Surgery, Perelman School of Medicine, University of Pennsylvania, 3400 Spruce Street, Philadelphia, PA 19104.

## License

\[Add a license here, e.g. MIT or GPL-3.0, and include a `LICENSE` file in the repository.\]
