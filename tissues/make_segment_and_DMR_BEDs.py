#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# makes universe files
# gets beds for upset plots

import os 
import pandas as pd
import numpy as np
from pybedtools import BedTool

def build_regions(df: pd.DataFrame,
                  chrom_col: str = "chr",
                  pos_col: str = "pos",
                  max_gap_bp: int = 100,
                  min_cpgs: int = 3,
                  min_len_bp: int = 50) -> pd.DataFrame:
    """
    Create regions by grouping consecutive CpGs where distance <= max_gap_bp.
    BED coordinates: start is 0-based, end is 0-based exclusive.
    We'll treat each CpG at position `pos` (1-based) as [pos-1, pos).
    Region spans from first CpG start to last CpG end.
    """
    # Keep only needed cols, drop missing
    d = df[[chrom_col, pos_col]].dropna().copy()
    d[pos_col] = d[pos_col].astype(int)

    # Sort
    d.sort_values([chrom_col, pos_col], inplace=True)

    # Determine breaks between CpGs per chromosome
    # A new region starts when:
    #   - chromosome changes, OR
    #   - gap between current and previous CpG > max_gap_bp
    prev_pos = d.groupby(chrom_col)[pos_col].shift(1)
    gap = d[pos_col] - prev_pos
    new_region = (prev_pos.isna()) | (gap > max_gap_bp)

    # Assign region IDs per chromosome
    d["region_id"] = new_region.groupby(d[chrom_col]).cumsum()

    # Aggregate to regions
    agg = (d.groupby([chrom_col, "region_id"])
             .agg(first_pos=(pos_col, "min"),
                  last_pos=(pos_col, "max"),
                  n_cpgs=(pos_col, "size"))
             .reset_index(drop=False))

    # Convert to BED coords:
    # first CpG at pos => start = pos-1
    # last CpG at pos  => end   = pos (exclusive for that CpG)
    agg["start"] = agg["first_pos"] - 1
    agg["end"] = agg["last_pos"]  # exclusive
    agg["length_bp"] = agg["end"] - agg["start"]

    # Filters
    keep = (agg["n_cpgs"] >= min_cpgs) & (agg["length_bp"] >= min_len_bp)
    out = agg.loc[keep, [chrom_col, "start", "end", "n_cpgs", "length_bp"]].copy()
    out.rename(columns={chrom_col: "chrom"}, inplace=True)

    # Optional name column
    out.insert(3, "name", [f"bg_{i}" for i in range(len(out))])

    return out

# Convert overlap pairs into true clipped overlap intervals
def clipped_overlap(row):
    a_start, a_end = int(row[1]), int(row[2])
    b_start, b_end = int(row[4]), int(row[5])
    return [
        row[0],
        max(a_start, b_start),
        min(a_end, b_end),
    ]


if __name__ == '__main__':
    # reciprocal overlap threshold
    F = 0.2
    
    DIR1 = '/scratch/eger/projects/MethSmoothEval/tissues/modkit/dmr'
    DIR2 = '/scratch/eger/projects/MethSmoothEval/tissues/DSS'
    OUTDIR = '/scratch/eger/projects/MethSmoothEval/tissues/DM_analysis'

    sample_pairs = ['Bulk_FC_Control_02_v_hg002_blood', 
                   'PPMI_3404_v_HBCC_81951_FTX',
                   'PPMI_3404_v_HBCC_82044_FTX']

    seg_cols = [
        'chrom',
        'start',
        'end',
        'state_name',
        'score',
        'n_sites',
        'a_counts',
        'b_counts',
        'a_mod_percentages',
        'b_mod_percentages',
        'a_frac_modified',
        'b_frac_modified',
        'effect_size',
        'cohen_h', 
        'cohen_h_low',
        'cohen_h_high'
    ]

    # no effect size filtering
    for pair in sample_pairs:
        fname1 = os.path.join(DIR1, pair+'.dmr.segments.txt')
        fname2 = os.path.join(DIR2, pair+'_DMLtestwSmoothing_DMRs.tsv')
        fname3 = os.path.join(DIR2, pair+'_DMLtestwSmoothing.tsv.gz')

        ## create the universe
        out_bed1 = os.path.join(OUTDIR, pair+'.background.50bp_3CpG.bed')
        if not os.path.exists(out_bed1):
            # load the DMPs
            df3 = pd.read_csv(fname3, sep='\t')
            
            background = build_regions(df3)
            background[['chrom', 'start', 'end']].to_csv(out_bed1, sep='\t', index=False, header=False)
        
        # save the DMRs as beds
        out_bed2 = os.path.join(OUTDIR, pair+'.different_segments.50bp_3CpG.bed')
        out_bed3 = os.path.join(OUTDIR, pair+'.DMRs.50bp_3CpG.bed')

        if not os.path.exists(out_bed2):
            # load the segments
            df1 = pd.read_csv(fname1, sep='\t', header=None, names=seg_cols)
            df1['length'] = df1['end'] - df1['start']
            df1 = df1[(df1['n_sites'] >= 3) & (df1['length'] >= 50) & (df3['state_name'] == 'different')]

            df1[['chrom', 'start', 'end']].to_csv(out_bed2, sep='\t', index=False, header=False)
        
        if not os.path.exists(out_bed3):
            # load the DMRs
            df2 = pd.read_csv(fname2, sep='\t')
            df2['start'] = df2['start'].astype(int)
            df2['end'] = df2['end'].astype(int)
            
            df2[['chr', 'start', 'end']].to_csv(out_bed3, sep='\t', index=False, header=False)

        # convert to bedtools
        bt1 = BedTool(out_bed2)
        bt2 = BedTool(out_bed3)

        # Intersect bt1 and bt2 reciprocally, returning overlap *pairs*
        overlap = bt1.intersect(bt2, wa=True, wb=True, f=F, r=True)

        df = overlap.to_dataframe(disable_auto_names=True, header=None)
        overlap_df = df.apply(clipped_overlap, axis=1, result_type="expand")
        overlap_df.columns = ["chrom", "start", "end"]
        
        shared_regions = BedTool.from_dataframe(overlap_df).sort().merge()

        modkit_unique = bt1.subtract(shared_regions)
        dss_unique = bt2.subtract(shared_regions)

        print(pair)
        print('Total modkit:', len(bt1))
        print('Total DSS:', len(bt2))
        print('Total shared:', len(shared_regions))
        print('Unique modkit:', len(modkit_unique))
        print('Unique DSS:', len(dss_unique))
        print('\n')

        # save the beds
        out_bed4 = os.path.join(OUTDIR, pair+'.modkit.50bp_3CpG.unique.bed')
        out_bed5 = os.path.join(OUTDIR, pair+'.DSS.50bp_3CpG.unique.bed')
        out_bed6 = os.path.join(OUTDIR, pair+'.50bp_3CpG.shared.bed')

        modkit_unique.saveas(out_bed4)
        dss_unique.saveas(out_bed5)
        shared_regions.saveas(out_bed6)

        print(out_bed4, out_bed5, out_bed6)
        print('\n')

        """
        # make beds of regions unique to DSS and modkit (< 20% overlap)
        # -v      Only report those entries in A that have _no overlaps_ with B
        unique_bt1 = bt1.intersect(bt2, v=True, f=0.2)
        unique_bt2 = bt2.intersect(bt1, v=True, f=0.2)

        shared_bt1 = bt1.intersect(bt2, u=True, f=0.2)
        shared_bt2 = bt2.intersect(bt1, u=True, f=0.2)

        out_bed4 = os.path.join(OUTDIR, pair+'.different_segments.50bp_3CpG.unique.bed')
        out_bed5 = os.path.join(OUTDIR, pair+'.DMRs.50bp_3CpG.unique.bed')
        out_bed6 = os.path.join(OUTDIR, pair+'.different_segments.50bp_3CpG.shared.bed')
        out_bed7 = os.path.join(OUTDIR, pair+'.DMRs.50bp_3CpG.shared.bed')
        
        unique_bt1.saveas(out_bed4)
        unique_bt2.saveas(out_bed5)
        shared_bt1.saveas(out_bed6)
        shared_bt2.saveas(out_bed7)
        """

    """
    # with effect size filtering
    for pair in sample_pairs:
        fname1 = os.path.join(DIR1, pair+'.dmr.segments.txt')
        fname2 = os.path.join(DIR2, pair+'_DMLtestwSmoothing_DMRs.tsv')

        # load the segments
        df1 = pd.read_csv(fname1, sep='\t', header=None, names=seg_cols)
        df1['length'] = df1['end'] - df1['start']

        # load the DMRs
        df2 = pd.read_csv(fname2, sep='\t')

        # subset segments to match DMR settings
        df3 = df1[(df1['n_sites'] >= 3) & (df1['length'] >= 50)]
        
        # save the beds
        out_bed1 = os.path.join(OUTDIR, pair+'.all_segments.50bp_3CpG.bed')
        out_bed2 = os.path.join(OUTDIR, pair+'.different_segments.50bp_3CpG.diff_01.bed')
        out_bed3 = os.path.join(OUTDIR, pair+'.DMRs.50bp_3CpG.diff_01.bed')

        if not os.path.exists(out_bed1):
            df3[['chrom', 'start', 'end']].to_csv(out_bed1, sep='\t', index=False, header=False)

        if not os.path.exists(out_bed2):
            df3[df3['state_name'] == 'different'][
            ['chrom', 'start', 'end']].to_csv(out_bed2, sep='\t', index=False, header=False)
        
        if not os.path.exists(out_bed3):
            df2[df2['diff.Methy'].abs() > 0.1][
            ['chr', 'start', 'end']].to_csv(out_bed3, sep='\t', index=False, header=False)

        # convert to bedtools
        bt1 = BedTool(out_bed2)
        bt2 = BedTool(out_bed3)

        # make beds of regions unique to DSS and modkit (< 20% overlap)
        # -v      Only report those entries in A that have _no overlaps_ with B
        unique_bt1 = bt1.intersect(bt2, v=True, f=0.2)
        unique_bt2 = bt2.intersect(bt1, v=True, f=0.2)

        out_bed4 = os.path.join(OUTDIR, pair+'.different_segments.50bp_3CpG.diff_01.unique.bed')
        out_bed5 = os.path.join(OUTDIR, pair+'.DMRs.50bp_3CpG.diff_01.unique.bed')
        
        unique_bt1.saveas(out_bed4)
        unique_bt2.saveas(out_bed5)

    """
        

