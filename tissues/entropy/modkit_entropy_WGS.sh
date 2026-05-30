#!/bin/bash

# samtools 1.21
# modkit v0.5.1-rc1
THREADS=24

REF_DIR="/scratch/eger/projects/MethSmoothEval/datasets/reference_genomes"
REF="GCA_000001405.15_GRCh38_no_alt_analysis_set.fa"

# Input directories and BAM files
MAINDIR="/scratch/eger/projects/MethSmoothEval/tissues"
OUTDIR1="${MAINDIR}/modkit/entropy"

SAMPLE="Bulk_FC_Control_02"
BAM="${MAINDIR}/data_transfer/rt38291/Bulk_FC_Control_02.phased.bam"

SAMPLE="hg002_blood"
BAM="/scratch/eger/projects/MethSmoothEval/tissues/clair3/hg002_blood/hg002_blood.phased.bam"

# modkit entropy --force --threads 24 
# --in-bam PPMI_3404/PPMI_3404_merged.sorted_minimap2_bq10_filtered.bam 
# -o ENTROPY/PPMI_3404/PPMI_3404.haplotagged.modkit_entropy.bedgraph 
# --cpg --ref /data/CARDPB/resources/hg38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa 
# --log-filepath ENTROPY/PPMI_3404/PPMI_3404.haplotagged.modkit_entropy.log

modkit entropy \
    --force \
    --cpg \
    --in-bam "${BAM}" \
    --ref "${REF_DIR}/${REF}" \
    --threads "${THREADS}" \
    --log-filepath "${OUTDIR1}/${SAMPLE}.phased.modkit_entropy.log" \
    -o "${OUTDIR1}/${SAMPLE}.phased.modkit_entropy.bedgraph"