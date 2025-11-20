#!/bin/bash
#SBATCH --cpus-per-task=40
#SBATCH --mem=100g
#SBATCH --mail-type=BEGIN,TIME_LIMIT_90,END
#SBATCH --time=24:00:00

module load samtools clair3 whatshap

# MODIFY
SAMPLE_ID="hg002_blood"
BAM="${SAMPLE_ID}.sorted.bam" # minimap2 aligned BAM
REF_DIR="/scratch/eger/projects/MethSmoothEval/datasets/reference_genomes"
MODEL_DIR="/home/eger/references/clair3"
INDIR="/scratch/eger/projects/MethSmoothEval/tissues/data_transfer/rt38291/hg002_blood_combined"
OUTDIR="/scratch/eger/projects/MethSmoothEval/tissues/clair3/${SAMPLE_ID}"

# presets
PLATFORM="ont"
REF="GCA_000001405.15_GRCh38_no_alt_analysis_set.fa"
# MODEL="r1041_e82_400bps_sup_v430" # 5kHz
MODEL="r1041_e82_400bps_sup_v410" # 4kHz
# already be available in # $CLAIR3_MODELS on biowulf?
# download: https://github.com/nanoporetech/rerio/blob/master/clair3_models/r1041_e82_400bps_sup_v410_model

mkdir -p ${OUTDIR}
cd ${OUTDIR}

# variant calling & phasing
clair3 \
    --bam_fn=${INDIR}/${BAM} \
    --ref_fn=${REF_DIR}/${REF} \
    --threads=${SLURM_CPUS_PER_TASK} \
    --platform=${PLATFORM} \
    --model_path=${MODEL_DIR}/${MODEL} \
    --output=${OUTDIR} \
    --enable_phasing

# haplotagging
whatshap haplotag \
    -r ${REF_DIR}/${REF} \
    --ignore-read-groups \
    --output-threads ${SLURM_CPUS_PER_TASK} \
    -o ${OUTDIR}/${SAMPLE_ID}.phased.bam \
    ${OUTDIR}/phased_merge_output.vcf.gz \
    ${INDIR}/${BAM}
samtools index -@ ${SLURM_CPUS_PER_TASK} ${OUTDIR}/${SAMPLE_ID}.phased.bam
