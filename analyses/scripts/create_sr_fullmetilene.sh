#!/usr/bin/bash

set -e

mkdir -p outputs/filt_bedgraph && cd outputs/filt_bedgraph ;

# Filtering and creating sorted bedgraph
find /mnt/project/analysis/Senthil_analysis/ -type f -name "*bedGraph" | grep -v 'MinCov5' | while read infile1; do
  outfile1=$(basename "${infile1}" | sed -e 's/.bedGraph$/.filt.bedgraph/') ; # Output file name
  echo "Processing - ${outfile1}"
  if [[ ! -f "${outfile1}" ]]; then
    cat ${infile1} | awk -F'\t' '$5+$6 >= 5' | awk -F'\t' 'BEGIN{OFS="\t"} {print $1,$2,$3,$4}' >tmp.${outfile1} ; # Atleast 5X coverage and get coords
    sortBed -i tmp.${outfile1} > ${outfile1} ; # Sorting the bedfile
    rm -f tmp.${outfile1} ;
  fi  
done  


find /mnt/project/analysis/Senthil_analysis/ -type f -name "*bedGraph" | grep -e 'MinCov5' | while read infile1; do
  outfile1=$(basename "${infile1}" | sed -e 's/.bedGraph$/.filt.bedgraph/') ; # Output file name
  echo "Processing - ${outfile1}"
  if [[ ! -f "${outfile1}" ]]; then
    cat ${infile1} | awk -F'\t' '$6>=5' | awk -F'\t' 'BEGIN{OFS="\t"} {print $1,$2,$3,$4}' >tmp.${outfile1} ; # Atleast 5X coverage and get coords
    sortBed -i tmp.${outfile1} > ${outfile1} ; # Sorting the bedfile
    rm -f tmp.${outfile1} ;
  fi
done

cd ../../;

###

mkdir -p outputs/metilene && cd outputs/metilene ;

\ls -d /opt/notebooks/outputs/filt_bedgraph/* | grep -e '_TruSeq_' -e '_EMSeq_' -e 'TrueMethyl' >lst_files_bedgraph.sr.txt ;

# Read items into an array
mapfile -t items < lst_files_bedgraph.sr.txt ;

num_items=${#items[@]}

# Iterate through the array to create pairs
printf "" >lst_pairs_bedgraph.sr.txt ;
for ((i = 0; i < num_items; i++)); do
  for ((j = i + 1; j < num_items; j++)); do
    echo "${items[i]} ${items[j]}" >>lst_pairs_bedgraph.sr.txt ;
  done
done

cd ../../ ;


# Filtering pairs for proper comparison
grep -e 'HG005.*HG002' -e 'HG002.*HG005' outputs/metilene/lst_pairs_bedgraph.sr.txt | grep -e '_TruSeq_.*_TruSeq_' -e '_EMSeq_.*_EMSeq_' -e '_TrueMethylOX_.*_TrueMethylOX_' -e '_TrueMethylBS_.*_TrueMethylBS_' >outputs/metilene/lst_pairs_bedgraph.sr.filt.txt
### MANUAL FILTERING NEEDED 
printf "PLEASE FILTER MANUALLY for appropriate pairs : outputs/metilene/lst_pairs_bedgraph.sr.filt.txt\n\n"

