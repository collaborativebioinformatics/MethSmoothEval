#!/usr/bin/bash

set -e

mkdir -p outputs && cd outputs ;

find ../datasets/ -type f -name "*.bam.bai" | fgrep '.ont.' | while read bamindex1; do
  filename1=$(basename "${bamindex1}") ;
  echo "$filename1" ;
  modkit pileup ${bamindex1/.bai/} ${filename1/.bam.bai/.pileup.bed} --log-filepath ${filename1/.bam.bai/.pileup.log}
done

