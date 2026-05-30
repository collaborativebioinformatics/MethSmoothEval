#!/bin/bash

# 5kHz samples from baylor

eval "$(conda shell.bash hook)"
conda activate aplanat
# samtools 1.21
# modkit v0.5.1-rc1

THREADS=60

SCRIPT_DIR="/home/eger/scripts/MethSmoothEval/tissues"
SOFT_DIR="/home/eger/software/CARDlongread-ONT-meth-benchmark"
REF_DIR="/scratch/eger/projects/MethSmoothEval/datasets/reference_genomes"

# Input directories and BAM files
MAINDIR="/scratch/eger/projects/MethSmoothEval/tissues"
OUTDIR1="${MAINDIR}/modkit/sample-probs"
OUTDIR2="${MAINDIR}/modkit/entropy"
OUTDIR3="${MAINDIR}/modkit/dmr"
OUTDIR4="${MAINDIR}/DM_analysis/DMRs_for_entropy"
OUTDIR5="${MAINDIR}/DSS"

REF="${REF_DIR}/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa"

# SAMPLE1="Bulk_FC_Control_02"
# SAMPLE2="hg002_blood"
SAMPLE1="hg002_blood"
SAMPLE2="Bulk_FC_Control_02"

SAMPLES=(${SAMPLE1} ${SAMPLE2})
# BAM1="${MAINDIR}/data_transfer/rt38291/Bulk_FC_Control_02.phased.bam"
# BAM2="${MAINDIR}/clair3/hg002_blood/hg002_blood.phased.bam"
BAM1="${MAINDIR}/clair3/hg002_blood/hg002_blood.phased.bam"
BAM2="${MAINDIR}/data_transfer/rt38291/Bulk_FC_Control_02.phased.bam"

# # compare modkit sample-probs methylation benchmarks across 
# # different samples by drawing methylation likelihood lineplots.
# mkdir -p ${OUTDIR1}/probs_comparison
# cd ${OUTDIR1}

# modkit sample-probs \
#     --force \
#     --hist \
#     -t ${THREADS} \
#     ${BAM1} \
#     -o ${OUTDIR1}/${SAMPLE1} \
#     --prefix ${SAMPLE1}

# modkit sample-probs \
#     --force \
#     --hist \
#     -t ${THREADS} \
#     ${BAM2} \
#     -o ${OUTDIR1}/${SAMPLE2} \
#     --prefix ${SAMPLE2}

# python ${SOFT_DIR}/modkit_sample_probs_comparison.py \
#     --input ${SAMPLE1}/${SAMPLE1}_probabilities.tsv ${SAMPLE2}/${SAMPLE2}_probabilities.tsv \
#     --names ${SAMPLE1} ${SAMPLE2}\
#     --output ${OUTDIR1}/probs_comparison/${SAMPLE1}_${SAMPLE2}


## get entropy for DMRs
mkdir -p ${OUTDIR2}/DMR_entropy
cd ${OUTDIR2}/DMR_entropy

# unsmoothed DMRs
modkit entropy \
    --force \
    --threads ${THREADS} \
    --in-bam ${BAM1} \
	-o ${SAMPLE1}/${SAMPLE1}_v_${SAMPLE2}.dss_unsmoothed_dmrs \
    --cpg \
	--regions ${OUTDIR4}/${SAMPLE1}_v_${SAMPLE2}.DMLtest_DMRs.bed \
	--ref ${REF} \
	--log-filepath ${SAMPLE1}/${SAMPLE1}_v_${SAMPLE2}.dss_unsmoothed_dmrs.log

modkit entropy \
    --force \
    --threads ${THREADS} \
    --in-bam ${BAM2} \
	-o ${SAMPLE2}/${SAMPLE1}_v_${SAMPLE2}.dss_unsmoothed_dmrs \
    --cpg \
	--regions ${OUTDIR4}/${SAMPLE1}_v_${SAMPLE2}.DMLtest_DMRs.bed \
	--ref ${REF} \
	--log-filepath ${SAMPLE2}/${SAMPLE1}_v_${SAMPLE2}.dss_unsmoothed_dmrs.log

# smoothed DMRs
modkit entropy \
    --force \
    --threads ${THREADS} \
    --in-bam ${BAM1} \
	-o ${SAMPLE1}/${SAMPLE1}_v_${SAMPLE2}.dss_smoothed_dmrs \
    --cpg \
	--regions ${OUTDIR4}/${SAMPLE1}_v_${SAMPLE2}.DMLtestwSmoothing_DMRs.bed \
	--ref ${REF} \
	--log-filepath ${SAMPLE1}/${SAMPLE1}_v_${SAMPLE2}.dss_smoothed_dmrs.log

modkit entropy \
    --force \
    --threads ${THREADS} \
    --in-bam ${BAM2} \
	-o ${SAMPLE2}/${SAMPLE1}_v_${SAMPLE2}.dss_smoothed_dmrs \
    --cpg \
	--regions ${OUTDIR4}/${SAMPLE1}_v_${SAMPLE2}.DMLtestwSmoothing_DMRs.bed \
	--ref ${REF} \
	--log-filepath ${SAMPLE2}/${SAMPLE1}_v_${SAMPLE2}.dss_smoothed_dmrs.log

# segments
modkit entropy \
    --force \
    --threads ${THREADS} \
    --in-bam ${BAM1} \
	-o ${SAMPLE1}/${SAMPLE1}_v_${SAMPLE2}.modkit_dmr_segments \
    --cpg \
	--regions ${OUTDIR3}/${SAMPLE1}_v_${SAMPLE2}.dmr.segments.txt \
	--ref ${REF} \
	--log-filepath ${SAMPLE1}/${SAMPLE1}_v_${SAMPLE2}.modkit_dmr_segments.log

modkit entropy \
    --force \
    --threads ${THREADS} \
    --in-bam ${BAM2} \
	-o ${SAMPLE2}/${SAMPLE1}_v_${SAMPLE2}.modkit_dmr_segments \
    --cpg \
	--regions ${OUTDIR3}/${SAMPLE1}_v_${SAMPLE2}.dmr.segments.txt \
	--ref ${REF} \
	--log-filepath ${SAMPLE2}/${SAMPLE1}_v_${SAMPLE2}.modkit_dmr_segments.log

# visualize pairwise comparisons between methylation entropies 
# of two samples for which differentially methylated regions (DMRs) 
# were identified with several different methods (modkit dmr pair 
# fine-grained, DSS/bsseq with smoothing, and DSS/bsseq without smoothing)
mkdir -p ${OUTDIR2}/DMR_entropy/pairwise_comparison

python ${SOFT_DIR}/CARDlongread_methylation_entropy_pairwise_comparison.py \
    --sample_name_1 ${SAMPLE1} \
    --sample_name_2 ${SAMPLE2} \
    --sample_1_bulk_entropy ${OUTDIR2}/${SAMPLE1}.phased.modkit_entropy.bedgraph \
    --sample_2_bulk_entropy ${OUTDIR2}/${SAMPLE2}.phased.modkit_entropy.bedgraph \
	--sample_1_modkit_dmr_entropy ${SAMPLE1}/${SAMPLE1}_v_${SAMPLE2}.modkit_dmr_segments/regions.bed \
	--sample_2_modkit_dmr_entropy ${SAMPLE2}/${SAMPLE1}_v_${SAMPLE2}.modkit_dmr_segments/regions.bed \
	--sample_1_dss_unsmoothed_dmr_entropy ${SAMPLE1}/${SAMPLE1}_v_${SAMPLE2}.dss_unsmoothed_dmrs/regions.bed \
	--sample_2_dss_unsmoothed_dmr_entropy ${SAMPLE2}/${SAMPLE1}_v_${SAMPLE2}.dss_unsmoothed_dmrs/regions.bed \
	--sample_1_dss_smoothed_dmr_entropy ${SAMPLE1}/${SAMPLE1}_v_${SAMPLE2}.dss_smoothed_dmrs/regions.bed \
	--sample_2_dss_smoothed_dmr_entropy ${SAMPLE2}/${SAMPLE1}_v_${SAMPLE2}.dss_smoothed_dmrs/regions.bed \
	--modkit_dmr_segments ${OUTDIR3}/${SAMPLE1}_v_${SAMPLE2}.dmr.segments.txt \
	--dss_unsmoothed_dmrs ${OUTDIR5}/${SAMPLE1}_v_${SAMPLE2}_DMLtest_DMRs.tsv \
	--dss_smoothed_dmrs ${OUTDIR5}/${SAMPLE1}_v_${SAMPLE2}_DMLtestwSmoothing_DMRs.tsv \
	--output_prefix pairwise_comparison/${SAMPLE1}_v_${SAMPLE2}_entropy_comparison \
	--plot_title "${SAMPLE1} v. ${SAMPLE2} entropy comparison"




