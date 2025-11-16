#!/bin/bash

# singularity pull docker://hkubal/clair3:latest
# samtools 1.21

THREADS=30
PLATFORM="ont"
MODEL="r1041_e82_400bps_sup_v430"
REF="GCA_000001405.15_GRCh38_no_alt_analysis_set.fa"
REF_DIR="/scratch/eger/projects/MethSmoothEval/datasets/reference_genomes"
MODEL_DIR="/home/eger/references/clair3"

## hg002 blood
BAM1="hg002_blood.sorted.bam"
INDIR1="/scratch/eger/projects/MethSmoothEval/tissues/data_transfer/rt38291/hg002_blood_combined"
OUTDIR1="/scratch/eger/projects/MethSmoothEval/tissues/clair3/hg002_blood"
mkdir -p ${OUTDIR1}
cd ${OUTDIR1}

# variant calling & phasing
singularity run \
    --bind /usr/lib/locale/ \
    --bind /home/eger/software/:/home/eger/software/ \
    --bind /scratch/:/scratch/ \
    /home/eger/software/clair3_latest.sif \
    /opt/bin/run_clair3.sh \
    --bam_fn=${INDIR1}/${BAM1} \
    --ref_fn=${REF_DIR}/${REF} \
    --threads=${THREADS} \
    --platform=${PLATFORM} \
    --model_path=${MODEL_DIR}/${MODEL} \
    --output=${OUTDIR1} \
    --enable_phasing

# haplotagging
singularity run \
    --bind /usr/lib/locale/ \
    --bind /home/eger/software/:/home/eger/software/ \
    --bind /scratch/:/scratch/ \
    /home/eger/software/clair3_latest.sif \
    whatshap haplotag \
    -r ${REF_DIR}/${REF} \
    --ignore-read-groups \
    --output-threads ${THREADS} \
    -o ${OUTDIR1}/hg002_blood.phased.bam \
    ${OUTDIR1}/phased_merge_output.vcf.gz \
    ${INDIR1}/${BAM1}
# Total alignments processed:                  29488091
# Alignments that could be tagged:             13643164
# Alignments spanning multiple phase sets:            4
samtools index -@ ${THREADS} ${OUTDIR1}/hg002_blood.phased.bam

## colo829bl
BAM2="colo829bl.merged.bam"
INDIR2="/scratch/eger/projects/MethSmoothEval/tissues/data_transfer/colo829bl"
OUTDIR2="/scratch/eger/projects/MethSmoothEval/tissues/clair3/colo829bl"
mkdir -p ${OUTDIR2}
cd ${OUTDIR2}

# merging
samtools merge \
    -@ ${THREADS} \
    -f ${OUTDIR2}/${BAM2} \
    ${INDIR2}/PAU59807.d052sup4305mCG_5hmCGvHg38.bam \
    ${INDIR2}/PAU61427.d052sup4305mCG_5hmCGvHg38.bam
samtools index -@ ${THREADS} ${OUTDIR2}/${BAM2}

# variant calling & phasing
singularity run \
    --bind /usr/lib/locale/ \
    --bind /home/eger/software/:/home/eger/software/ \
    --bind /scratch/:/scratch/ \
    /home/eger/software/clair3_latest.sif \
    /opt/bin/run_clair3.sh \
    --bam_fn=${OUTDIR2}/${BAM2} \
    --ref_fn=${REF_DIR}/${REF} \
    --threads=${THREADS} \
    --platform=${PLATFORM} \
    --model_path=${MODEL_DIR}/${MODEL} \
    --output=${OUTDIR2} \
    --enable_phasing

# haplotagging
singularity run \
    --bind /usr/lib/locale/ \
    --bind /home/eger/software/:/home/eger/software/ \
    --bind /scratch/:/scratch/ \
    /home/eger/software/clair3_latest.sif \
    whatshap haplotag \
    -r ${REF_DIR}/${REF} \
    --ignore-read-groups \
    --output-threads ${THREADS} \
    -o ${OUTDIR2}/colo829bl.phased.bam \
    ${OUTDIR2}/phased_merge_output.vcf.gz \
    ${OUTDIR2}/${BAM2}
# Total alignments processed:                  21603193
# Alignments that could be tagged:             11969319
# Alignments spanning multiple phase sets:           11
samtools index -@ ${THREADS} ${OUTDIR2}/colo829bl.phased.bam
    
