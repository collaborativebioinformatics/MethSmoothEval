#!/bin/bash

## INSTALL
# conda activate methyl
# git clone https://github.com/nloyfer/wgbs_tools.git
# git clone https://github.com/nloyfer/UXM_deconv.git

# cd wgbs_tools
# python setup.py

## RUN
THREADS=30
ATLAS="Atlas.U25.l4.hg38"
# ATLAS="Atlas.U250.l4.hg38"
ATLAS_FILE="/home/eger/software/UXM_deconv/supplemental/${ATLAS}.full.tsv"
ATLAS_BED="/home/eger/software/UXM_deconv/supplemental/${ATLAS}.bed"

# create atlas bed
awk 'BEGIN {OFS="\t"} NR>1 {print $1, $2-1, $3, $6, ".", "."}' ${ATLAS_FILE} > ${ATLAS_BED} 

# Inputs
MAINDIR="/scratch/eger/projects/MethSmoothEval/tissues"
INDIR1="${MAINDIR}/clair3"
OUTDIR="/scratch/eger/projects/MethSmoothEval/tissues/UXM_deconv/${ATLAS}"
mkdir -p ${OUTDIR}
cd ${OUTDIR}

declare -a SAMPLES=("Bulk_FC_Control_02" "hg002_blood" "colo829bl")

# make pat files
for SAMPLE in "${SAMPLES[@]}"; do
    if [ $SAMPLE == "Bulk_FC_Control_02" ]; then
        BAM="${MAINDIR}/data_transfer/rt38291/Bulk_FC_Control_02.phased.bam"
    else
        BAM="${MAINDIR}/clair3/${SAMPLE}/${SAMPLE}.phased.bam"
    fi    
    # all reads
    wgbstools bam2pat -np -L \
        ${ATLAS_BED} \
        -@ ${THREADS} \
        --genome hg38 \
        --force \
        --no_beta \
        --out_dir ${OUTDIR} \
        ${BAM}

    # HP1 reads
    wgbstools bam2pat -np -L \
        ${ATLAS_BED} \
        -@ ${THREADS} \
        --genome hg38 \
        --force \
        --no_beta \
        --out_dir ${OUTDIR} \
        "${MAINDIR}/clair3/${SAMPLE}/${SAMPLE}.HP1.bam"

    # HP2 reads
    wgbstools bam2pat -np -L \
        ${ATLAS_BED} \
        -@ ${THREADS} \
        --genome hg38 \
        --force \
        --no_beta \
        --out_dir ${OUTDIR} \
        "${MAINDIR}/clair3/${SAMPLE}/${SAMPLE}.HP2.bam"
done

# deconvolve brain
uxm deconv \
    --include Neuron Oligodend \
    --atlas ${ATLAS_FILE} \
    --output Bulk_FC_Control_02.${ATLAS}.UMX_deconv.csv \
    --threads ${THREADS} \
    Bulk_FC_Control_02.*.pat.gz

uxm plot Bulk_FC_Control_02.${ATLAS}.UMX_deconv.csv \
    -o ~/scripts/MethSmoothEval/tissues/Bulk_FC_Control_02.${ATLAS}.UMX_deconv.pdf

# deconvolve blood
uxm deconv \
    --include Blood-B Blood-Granul Blood-Mono+Macro Blood-NK Blood-T \
    --atlas ${ATLAS_FILE} \
    --output Blood.${ATLAS}.UMX_deconv.csv \
    --threads ${THREADS} \
    *bl*.pat.gz

uxm plot Blood.${ATLAS}.UMX_deconv.csv \
    -o ~/scripts/MethSmoothEval/tissues/Blood.${ATLAS}.UMX_deconv.pdf
