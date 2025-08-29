#!/usr/bin/bash

set -e

mkdir -p outputs/metilene/ && cd outputs/metilene/ ;
mkdir -p inputs/ results/ ;

cat lst_pairs_bedgraph.txt | while read line1; do 
  file1=$(echo $line1 | sed -e 's/ .*//'); 
  file2=$(echo $line1 | sed -e 's/.* //');
  name1=$(basename "$file1" | sed -e 's/.chr22.*//' -e 's/.pacbio.*//') ;
  name2=$(basename "$file2" | sed -e 's/.chr22.*//' -e 's/.pacbio.*//') ; 
  printf " PROCESSing PacBio ${name1}__${name2}\n\tSTART\n" ;
  metilene_input.pl -in1 $file1 -in2 $file2 -o inputs/${name1}__${name2}.pacbio.metilene.input ; # Creating Metiline input file format

  if [[ ! -f "results/${name1}__${name2}.pacbio.metilene.out.tsv" ]]; then # Avoid rerun if output exists
    metilene -t 16 -a g1 -b g2 inputs/${name1}__${name2}.pacbio.metilene.input >results/${name1}__${name2}.pacbio.metilene.out.tsv 2>results/${name1}__${name2}.pacbio.metilene.log ; # Running Metilene
  fi
done
  
