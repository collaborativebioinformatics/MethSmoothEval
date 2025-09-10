#!/usr/bin/bash

set -e

###

mkdir -p outputs/metilene/ && cd outputs/metilene/ ;
mkdir -p inputs/ results/ ;

#
cat lst_pairs_bedgraph.sr.filt.txt | grep -e 'smoothed' | while read line1; do 
  file1=$(echo $line1 | sed -e 's/ .*//'); 
  file2=$(echo $line1 | sed -e 's/.* //');
  name1=$(basename "$file1" | sed -e 's/.chr22.*//' ) ;
  name2=$(basename "$file2" | sed -e 's/.chr22.*//' ) ; 
  printf " PROCESSing smoothed SR ${name1}__${name2} : $file1 $file2\n\tSTART\n" ;
  metilene_input.pl -in1 $file1 -in2 $file2 -o inputs/${name1}__${name2}.smoothed.sr.metilene.input ; # Creating Metiline input file format

  if [[ ! -f "results/${name1}__${name2}.smoothed.sr.metilene.out.tsv" ]]; then # Avoid rerun if output exists
    metilene -t 16 -a g1 -b g2 inputs/${name1}__${name2}.smoothed.sr.metilene.input >results/${name1}__${name2}.smoothed.sr.metilene.out.tsv 2>results/${name1}__${name2}.smoothed.sr.metilene.log ; # Running Metilene
  fi
done


#
cat lst_pairs_bedgraph.sr.filt.txt | grep -v -e 'smoothed' | while read line1; do
  file1=$(echo $line1 | sed -e 's/ .*//'); 
  file2=$(echo $line1 | sed -e 's/.* //');
  name1=$(basename "$file1" | sed -e 's/.chr22.*//' ) ;
  name2=$(basename "$file2" | sed -e 's/.chr22.*//' ) ; 
  printf " PROCESSing raw SR ${name1}__${name2} : $file1 $file2\n\tSTART\n" ;
  metilene_input.pl -in1 $file1 -in2 $file2 -o inputs/${name1}__${name2}.raw.sr.metilene.input ; # Creating Metiline input file format

  if [[ ! -f "results/${name1}__${name2}.raw.sr.metilene.out.tsv" ]]; then # Avoid rerun if output exists
    metilene -t 16 -a g1 -b g2 inputs/${name1}__${name2}.raw.sr.metilene.input >results/${name1}__${name2}.raw.sr.metilene.out.tsv 2>results/${name1}__${name2}.raw.sr.metilene.log ; # Running Metilene
  fi
done
  
