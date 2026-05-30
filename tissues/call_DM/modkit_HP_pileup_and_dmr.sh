#!/bin/bash

# samtools 1.21
# modkit v0.5.1-rc1
THREADS=30

REF_DIR="/scratch/eger/projects/MethSmoothEval/datasets/reference_genomes"
REF="GCA_000001405.15_GRCh38_no_alt_analysis_set.fa"

# Input directories and BAM files
MAINDIR="/scratch/eger/projects/MethSmoothEval/tissues"
declare -a SAMPLES=("Bulk_FC_Control_02" "hg002_blood" "colo829bl")

# Output directories
OUTDIR1="${MAINDIR}/modkit/pileup"
OUTDIR2="${MAINDIR}/modkit/dmr"
mkdir -p "${OUTDIR1}" "${OUTDIR2}"

# Loop over all samples
for SAMPLE in "${SAMPLES[@]}"; do
    if [ ${SAMPLE} == "Bulk_FC_Control_02" ]; then
        BAM="${MAINDIR}/data_transfer/rt38291/Bulk_FC_Control_02.phased.bam"
    else
        BAM="${MAINDIR}/clair3/${SAMPLE}/${SAMPLE}.phased.bam"
    fi    
    BASE="${OUTDIR1}/${SAMPLE}.cpg"
    OUTFILE="${BASE}.bed.gz"

    echo "======================================"
    echo "Processing ${SAMPLE} ..."
    echo "======================================"

    ### --- Standard CpG pileup ---
    if [[ -s "${OUTFILE}" ]]; then
        echo "Skipping ${SAMPLE}: ${OUTFILE} already exists."
    else
        echo "Running standard CpG pileup for ${SAMPLE} ..."
        modkit pileup \
            --cpg \
            --ignore h \
            --combine-strands \
            --ref "${REF_DIR}/${REF}" \
            --threads "${THREADS}" \
            --log-filepath "${BASE}.log" \
            "${BAM}" - | \
            bgzip -c > "${OUTFILE}"
        tabix -p bed "${OUTFILE}"
        echo "Finished CpG pileup for ${SAMPLE}."
    fi

    ### --- HP-partitioned pileup ---
    HP_PREFIX="${BASE}.hp"
    HP1="${HP_PREFIX}_1.bed.gz"
    HP2="${HP_PREFIX}_2.bed.gz"
    HPUNG="${HP_PREFIX}_ungrouped.bed.gz"

    if [[ -s "${HP1}" && -s "${HP2}" && -s "${HPUNG}" ]]; then
        echo "Skipping ${SAMPLE}: HP bed files already exist."
    else
        echo "Running HP-partitioned CpG pileup for ${SAMPLE} ..."
        modkit pileup --cpg \
            --partition-tag HP \
            --ignore h \
            --combine-strands \
            --ref "${REF_DIR}/${REF}" \
            --threads "${THREADS}" \
            --prefix "${HP_PREFIX}" \
            --log-filepath "${HP_PREFIX}.log" \
            "${BAM}" \
            "${OUTDIR1}"

        # Compress and index all HP bed files
        for HPBED in "${HP_PREFIX}"_1.bed "${HP_PREFIX}"_2.bed "${HP_PREFIX}"_ungrouped.bed; do
            if [[ -s "${HPBED}" ]]; then
                bgzip "${HPBED}"
                tabix -p bed "${HPBED}.gz"
            fi
        done
        echo "Finished HP-partitioned CpG pileup for ${SAMPLE}."
    fi

    ### --- DMR analysis between haplotypes ---
    DMR_OUT="${OUTDIR2}/${SAMPLE}.hp_dmr.bed"
    DMR_SEG="${OUTDIR2}/${SAMPLE}.hp_dmr.segments.txt"
    DMR_LOG="${OUTDIR2}/${SAMPLE}.hp_dmr.log"

    if [[ -s "${DMR_OUT}" && -s "${DMR_SEG}" ]]; then
        echo "Skipping ${SAMPLE}: DMR output already exists."
    else
        echo "Running DMR analysis for ${SAMPLE} ..."
        modkit dmr pair \
            -a "${HP1}" \
            -b "${HP2}" \
            -o "${DMR_OUT}" \
            --segment "${DMR_SEG}" \
            --min-valid-coverage 5 \
            --fine-grained \
            --ref "${REF_DIR}/${REF}" \
            --base C \
            --threads "${THREADS}" \
            --log-filepath "${DMR_LOG}"
        bgzip "${DMR_OUT}"
        echo "Finished DMR analysis for ${SAMPLE}."
    fi

    echo
done

### --- DMR analysis between samples ---
SAMP1="Bulk_FC_Control_02"
SAMP2="hg002_blood"
modkit dmr pair \
    -a ${OUTDIR1}/${SAMP1}.cpg.bed.gz \
    -b ${OUTDIR1}/${SAMP2}.cpg.bed.gz \
    -o ${OUTDIR2}/${SAMP1}_v_${SAMP2}.dmr.bed \
    --segment ${OUTDIR2}/${SAMP1}_v_${SAMP2}.dmr.segments.txt \
    --min-valid-coverage 5 \
    --fine-grained \
    --ref "${REF_DIR}/${REF}" \
    --base C \
    --threads "${THREADS}" \
    --log-filepath ${OUTDIR2}/${SAMP1}_v_${SAMP2}.dmr.log
bgzip ${OUTDIR2}/${SAMP1}_v_${SAMP2}.dmr.bed

SAMP1="Bulk_FC_Control_02"
SAMP2="colo829bl"
modkit dmr pair \
    -a ${OUTDIR1}/${SAMP1}.cpg.bed.gz \
    -b ${OUTDIR1}/${SAMP2}.cpg.bed.gz \
    -o ${OUTDIR2}/${SAMP1}_v_${SAMP2}.dmr.bed \
    --segment ${OUTDIR2}/${SAMP1}_v_${SAMP2}.dmr.segments.txt \
    --min-valid-coverage 5 \
    --fine-grained \
    --ref "${REF_DIR}/${REF}" \
    --base C \
    --threads "${THREADS}" \
    --log-filepath ${OUTDIR2}/${SAMP1}_v_${SAMP2}.dmr.log
bgzip ${OUTDIR2}/${SAMP1}_v_${SAMP2}.dmr.bed

SAMP1="hg002_blood"
SAMP2="colo829bl"
modkit dmr pair \
    -a ${OUTDIR1}/${SAMP1}.cpg.bed.gz \
    -b ${OUTDIR1}/${SAMP2}.cpg.bed.gz \
    -o ${OUTDIR2}/${SAMP1}_v_${SAMP2}.dmr.bed \
    --segment ${OUTDIR2}/${SAMP1}_v_${SAMP2}.dmr.segments.txt \
    --min-valid-coverage 5 \
    --fine-grained \
    --ref "${REF_DIR}/${REF}" \
    --base C \
    --threads "${THREADS}" \
    --log-filepath ${OUTDIR2}/${SAMP1}_v_${SAMP2}.dmr.log
bgzip ${OUTDIR2}/${SAMP1}_v_${SAMP2}.dmr.bed
