#!/usr/bin/bash

set -e

mkdir -p outputs && cd outputs ;

find ../datasets/ -type f -name "*.bam.bai" | fgrep '.pacbio.' | while read bamindex1; do
  filename1=$(basename "${bamindex1}") ;
  echo "$filename1" ;
  aligned_bam_to_cpg_scores --bam ${bamindex1/.bai/} --output-prefix ${filename1/.bam.bai/} --threads 16
done

