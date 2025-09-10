#!/usr/bin/bash

set -e

mkdir -p outputs/metilene && cd outputs/metilene ;

\ls -d /opt/notebooks/outputs/filt_bedgraph/* | grep 'ont' >lst_files_bedgraph.ont.txt ;

# Read items into an array
mapfile -t items < lst_files_bedgraph.ont.txt ;

num_items=${#items[@]}

# Iterate through the array to create pairs
printf "" >lst_pairs_bedgraph.ont.txt ;
for ((i = 0; i < num_items; i++)); do
  for ((j = i + 1; j < num_items; j++)); do
    echo "${items[i]} ${items[j]}" >>lst_pairs_bedgraph.ont.txt ;
  done
done
