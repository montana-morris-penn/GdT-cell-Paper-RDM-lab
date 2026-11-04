#!/bin/bash

CELLRANGER_BIN="/mnt/isilon/tan_lab/sussmanj/Single_Cell_Tools/cellranger-9.0.1/cellranger"
FASTQ_DIR="/mnt/isilon/tan_lab/sussmanj/GIST_Multiome/GD_TCells/scRNA_Data/Original_FASTQ"
REF_DIR="/mnt/isilon/tan_lab/sussmanj/Single_Cell_Tools/refdata-gex-GRCh38-2024-A"

SAMPLES=(
    "RonDematteo-11062023-HG618-T-GEX"
    "RonDematteo-11062023-HG618-P-GEX"
    "RonDematteo-11062023-HG589-P-GEX"
    "RonDematteo-11032023-HG610-P-GEX"
    "RonDematteo-10112023-HG642-P-GEX"
    "RonDematteo-09252023-HG610-T-YDT-GEX"
)

for SAMPLE in "${SAMPLES[@]}"; do
    sbatch << EOF
#!/bin/bash
#SBATCH --job-name=${SAMPLE}
#SBATCH -c 8
#SBATCH --mem=64G
#SBATCH -t 99:00:00
#SBATCH --output=${SAMPLE}_%j.out
#SBATCH --error=${SAMPLE}_%j.err

${CELLRANGER_BIN} count \
    --id="${SAMPLE}" \
    --sample="${SAMPLE}" \
    --fastqs="${FASTQ_DIR}" \
    --transcriptome="${REF_DIR}" \
    --create-bam=true \
    --localcores=8 \
    --localmem=60
EOF
done
