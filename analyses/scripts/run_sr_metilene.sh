#!/usr/bin/bash

set -e

mkdir -p outputs/metilene/ && cd outputs/metilene/
mkdir -p inputs/ results/

cat lst_pairs_bedgraph.txt | while read line1; do 
    file1=$(echo $line1 | awk '{print $1}')
    file2=$(echo $line1 | awk '{print $2}')
    
    name1=$(basename "$file1" | sed -e 's/.chr22.*//' -e 's/.cpg.*//')
    name2=$(basename "$file2" | sed -e 's/.chr22.*//' -e 's/.cpg.*//')
    
    printf " PROCESSing ${name1}__${name2}\n\tSTART\n"
    
    # Create Metilene input file for SR
    metilene_input.pl -in1 "$file1" -in2 "$file2" -o inputs/${name1}__${name2}.metilene.input
    
    # Run Metilene if output doesn't exist
    if [[ ! -f "results/${name1}__${name2}.metilene.out.tsv" ]]; then
        metilene -t 16 -a g1 -b g2 inputs/${name1}__${name2}.metilene.input \
            > results/${name1}__${name2}.metilene.out.tsv \
            2> results/${name1}__${name2}.metilene.log
    fi
done
