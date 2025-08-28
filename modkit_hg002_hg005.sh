#!/bin/bash

# modkit 0.5.0
# ONT human data

DIR="/scratch/eger/projects/MethSmoothEval"
REF="${DIR}/datasets/reference_genomes/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa"

# create pileups
modkit pileup \
    --cpg \
    --ignore h \
    --combine-strands \
    --ref ${REF} \
    --threads 30 \
    "${DIR}/datasets/human/hg002/HG002.ont.chr22.bam" - | \
bgzip -c > "${DIR}/analysis/Sarah_analysis/human_ONT/HG002.ont.chr22.cpg.bed.gz"
tabix -p bed "${DIR}/analysis/Sarah_analysis/human_ONT/HG002.ont.chr22.cpg.bed.gz"

modkit pileup \
    --cpg \
    --ignore h \
    --combine-strands \
    --ref ${REF} \
    --threads 30 \
    "${DIR}/datasets/human/hg005/HG005.ont.chr22.bam" - | \
bgzip -c > "${DIR}/analysis/Sarah_analysis/human_ONT/HG005.ont.chr22.cpg.bed.gz"
tabix -p bed "${DIR}/analysis/Sarah_analysis/human_ONT/HG005.ont.chr22.cpg.bed.gz"

# get dmr output
modkit dmr pair \
    -a "${DIR}/analysis/Sarah_analysis/human_ONT/HG002.ont.chr22.cpg.bed.gz" \
    -b "${DIR}/analysis/Sarah_analysis/human_ONT/HG005.ont.chr22.cpg.bed.gz" \
    -o "${DIR}/analysis/Sarah_analysis/human_ONT/single_base_dmr_HG002_and_HG005.bed" \
    --ref ${REF} \
    --base C \
    --threads 30 \
    --log-filepath "${DIR}/analysis/Sarah_analysis/human_ONT/single_base_dmr_HG002_and_HG005.log"
