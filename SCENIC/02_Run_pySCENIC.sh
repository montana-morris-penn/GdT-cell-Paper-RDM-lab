#!/bin/bash
#SBATCH -c 75                   # Request N cores [>12 cores recommended for pyscenic2 steps]
#SBATCH -t 10-0            
#SBATCH --mem 700G              # total amount of memory requested [>100G recommended]
#SBATCH --job-name pyscenic    # Job name
#SBATCH -o %j.out               # File to which standard out will be written
#SBATCH -e %j.err               # File to which standard err will be written

mkdir pySCENIC_Output

#Merged all TFs together for maximum TFs
singularity exec --bind /mnt/isilon/tan_lab/sussmanj/ \
    -H /mnt/isilon/tan_lab/sussmanj/GIST_Multiome/GD_TCells \
    --cleanenv /mnt/isilon/tan_lab/sussmanj/GIST_Multiome/GD_TCells/aertslab-pyscenic-0.12.1.sif \
arboreto_with_multiprocessing.py \
        /mnt/isilon/tan_lab/sussmanj/GIST_Multiome/GD_TCells/Combined_pySCENIC_RNA_Filtered_Revision.loom \
        cisTarget_databases/hs_hgnc_tfs.txt \
        --num_workers 75 \
        -o pySCENIC_Output/pySCENIC_GRN_adjacencies_revision.csv \
        --method grnboost2 \
        --seed 69420

date
echo "regulatory networks inferred"
sleep 5


# Run pySCENIC get regulatory network from command line interface
singularity exec --bind /mnt/isilon/tan_lab/sussmanj/ \
    -H /mnt/isilon/tan_lab/sussmanj/GIST_Multiome/GD_TCells \
    --cleanenv /mnt/isilon/tan_lab/sussmanj/GIST_Multiome/GD_TCells/aertslab-pyscenic-0.12.1.sif \
    pyscenic ctx pySCENIC_Output/pySCENIC_GRN_adjacencies_revision.csv \
    /mnt/isilon/tan_lab/sussmanj/GIST_Multiome/GD_TCells/cisTarget_databases/*.mc9nr.genes_vs_motifs.rankings.feather \
    --annotations_fname /mnt/isilon/tan_lab/sussmanj/GIST_Multiome/GD_TCells/cisTarget_databases/motifs-v9-nr.hgnc-m0.001-o0.0.tbl \
    --expression_mtx_fname /mnt/isilon/tan_lab/sussmanj/GIST_Multiome/GD_TCells/Combined_pySCENIC_RNA_Filtered_Revision.loom \
    --output pySCENIC_Output/pySCENIC_CTX_regulons_revision.csv \
    --mask_dropouts \
    --num_workers 60

date
echo "regulon modules defined"
sleep 5

singularity exec --bind /mnt/isilon/tan_lab/sussmanj/ \
    -H /mnt/isilon/tan_lab/sussmanj/GIST_Multiome/GD_TCells \
    --cleanenv /mnt/isilon/tan_lab/sussmanj/GIST_Multiome/GD_TCells/aertslab-pyscenic-0.12.1.sif \
    pyscenic aucell \
    /mnt/isilon/tan_lab/sussmanj/GIST_Multiome/GD_TCells/Combined_pySCENIC_RNA_Filtered_Revision.loom \
    pySCENIC_Output/pySCENIC_CTX_regulons_revision.csv \
    --output pySCENIC_Output/pySCENIC_AUCell_Revision.loom \
    --num_workers 60

echo "cell scoring completed"

