#!/usr/bin/bash

### Creating Bedgraph files for Metilene output

set -e

mkdir -p outputs/filt_bedgraph && cd outputs/filt_bedgraph ;

find /mnt/project/analysis/ -type f -name "*bed*" | fgrep -v -e '.tbi' -e 'single_base' | grep pacbio | while read infile1; do
  outfile1=$(basename "${infile1}" | sed -e 's/.bed.gz$/.filt.bedgraph/') ; # Output file name
  echo "Processing - ${outfile1}"
  if [[ ! -f "${outfile1}" ]]; then
    zcat ${infile1} | tail -n +9 | awk -F'\t' '$6>=5' | awk -F'\t' 'BEGIN{OFS="\t"} {print $1,$2,$3,$4}' >tmp.${outfile1} ; # Atleast 5X coverage and get coords
    sortBed -i tmp.${outfile1} > ${outfile1} ; # Sorting the bedfile 
    rm -f tmp.${outfile1} ;
  fi  
done  

