#!/usr/bin/bash

set -e

mkdir -p output/metilene && cd output/metilene ;

\ls -d /opt/notebooks/outputs/filt_bedgraph/* >lst_files_bedgraph.txt ;

# Read items into an array
mapfile -t items < lst_files_bedgraph.txt ;

num_items=${#items[@]}

# Iterate through the array to create pairs
touch lst_pairs_bedgraph.txt ;
for ((i = 0; i < num_items; i++)); do
  for ((j = i + 1; j < num_items; j++)); do
    echo "${items[i]} ${items[j]}" >>lst_pairs_bedgraph.txt ;
  done
done
