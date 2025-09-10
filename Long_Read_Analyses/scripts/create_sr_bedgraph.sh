#!/usr/bin/env bash
set -e

mkdir -p outputs/filt_bedgraph
cd outputs/filt_bedgraph

# Find .bed or .bed.gz, exclude ONT/Pacbio/single_base/.tbi
find /mnt/project/analysis/ -type f \( -iname "*.bed" -o -iname "*.bed.gz" \) \
    | grep -vi -e '.tbi' -e 'single_base' -e 'ont' -e 'pacbio' \
    | while read infile; do

    outfile=$(basename "${infile}" | sed -e 's/.bed.gz$/.filt.bedgraph/' -e 's/.bed$/.filt.bedgraph/')

    echo "Processing - ${outfile}"

    if [[ ! -f "${outfile}" ]]; then
        # Use zcat for gzipped files, cat otherwise
        if file "$infile" | grep -q gzip; then
            reader="zcat $infile"
        else
            reader="cat $infile"
        fi

        $reader | awk 'BEGIN{OFS="\t"}
        /^#/ {next} 
        {
            # If column 4 is "m" or "u", take column 11 for methylation %
            if ($4=="m" || $4=="u") {
                print $1, $2, $3, $11
            } else {
                # Otherwise take column 4 (already numeric)
                print $1, $2, $3, $4
            }
        }' > tmp.${outfile}

        sortBed -i tmp.${outfile} > "${outfile}"
        rm -f tmp.${outfile}
    fi
done
