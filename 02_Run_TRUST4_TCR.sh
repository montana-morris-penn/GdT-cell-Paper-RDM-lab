#!/bin/bash

BASE_DIR="/mnt/isilon/tan_lab/sussmanj/GIST_Multiome/GD_TCells/scRNA_Data"
TRUST4_BIN="${BASE_DIR}/TRUST4/run-trust4"
TRUST4_FA="${BASE_DIR}/TRUST4/hg38_bcrtcr.fa"
TRUST4_REF="${BASE_DIR}/TRUST4/human_IMGT+C.fa"

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
#SBATCH --job-name=TRUST4_${SAMPLE}
#SBATCH -c 8
#SBATCH --mem=72G
#SBATCH -t 12:00:00
#SBATCH --output=${BASE_DIR}/TRUST4_${SAMPLE}_%j.out
#SBATCH --error=${BASE_DIR}/TRUST4_${SAMPLE}_%j.err

cd ${BASE_DIR}

${TRUST4_BIN} \
  -b ${BASE_DIR}/${SAMPLE}/outs/possorted_genome_bam.bam \
  --barcode CB \
  --UMI UB \
  -f ${TRUST4_FA} \
  --ref ${TRUST4_REF} \
  -t 8 \
  -o ${BASE_DIR}/${SAMPLE}_trust4

EOF
done
