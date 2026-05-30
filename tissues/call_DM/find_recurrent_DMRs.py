#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from upsetplot import UpSet
from pybedtools import BedTool


def _prepare_bed_df(df, value_col, pair_name, chrom_col='chrom'):
    """Return a BED-like dataframe with chr/start/end/value for one sample pair."""
    out = df[[chrom_col, 'start', 'end', value_col]].copy()
    out.columns = ['chr', 'start', 'end', pair_name]
    out['start'] = out['start'].astype(int)
    out['end'] = out['end'].astype(int)
    out = out.dropna(subset=['chr', 'start', 'end', pair_name])
    out = out[out['end'] > out['start']]
    return out.sort_values(['chr', 'start', 'end']).reset_index(drop=True)


def _bedtool_from_df(df):
    """Create a pybedtools BedTool from a dataframe with chr/start/end/value."""
    return BedTool.from_dataframe(df[['chr', 'start', 'end', df.columns[3]]])


def _pick_largest_per_locus(df):
    """
    Keep one entry per locus.

    Rows are grouped by exact output locus chr/start/end. If duplicate rows are
    produced by overlapping input intervals, keep the row whose contributing
    original regions span the largest total length.
    """
    if df.empty:
        return df

    df = df.copy()
    df['span_sum'] = (
        (df['a_end'] - df['a_start'])
        + (df['b_end'] - df['b_start'])
        + (df['c_end'] - df['c_start'])
    )
    df = (
        df.sort_values(['chr', 'start', 'end', 'span_sum'], ascending=[True, True, True, False])
          .drop_duplicates(['chr', 'start', 'end'], keep='first')
    )
    return df.drop(columns=['span_sum']).reset_index(drop=True)


def get_region_membership_df(df_a, df_b, df_c, pair_names):
    """
    Build an UpSet membership dataframe using original regions as the counting
    unit.

    The output has exactly one row per original input region across the three
    sample pairs, so its number of rows should equal:
        len(df_a) + len(df_b) + len(df_c)

    The boolean columns indicate whether each original region overlaps each of
    the three pair-specific region sets.
    """
    dfs = [
        df_a.reset_index(drop=True).copy(),
        df_b.reset_index(drop=True).copy(),
        df_c.reset_index(drop=True).copy(),
    ]
    bts = [_bedtool_from_df(df) for df in dfs]

    out_dfs = []

    for i, query_df in enumerate(dfs):
        query_bt = bts[i]

        out = query_df[['chr', 'start', 'end']].copy().reset_index(drop=True)
        out.insert(0, 'source_pair', pair_names[i])
        out.insert(1, 'region_id', [f'{pair_names[i]}__{k}' for k in range(len(out))])

        for j, target_bt in enumerate(bts):
            if i == j:
                out[pair_names[j]] = True
                continue

            # -c preserves one row per original query interval and appends the
            # number of overlaps with target_bt.
            counts = [int(feature.fields[-1]) for feature in query_bt.intersect(target_bt, c=True)]

            if len(counts) != len(out):
                raise RuntimeError(
                    f'pybedtools intersect returned {len(counts)} rows, but expected '
                    f'{len(out)} rows for {pair_names[i]} vs {pair_names[j]}'
                )

            out[pair_names[j]] = [count > 0 for count in counts]

        out_dfs.append(out)

    membership_df = pd.concat(out_dfs, ignore_index=True)

    expected_n = sum(len(df) for df in dfs)
    if len(membership_df) != expected_n:
        raise RuntimeError(
            f'Membership dataframe has {len(membership_df)} rows, but expected '
            f'{expected_n} rows from the three input dataframes.'
        )

    return membership_df


def _bedtool_from_indexed_df(df):
    """Create a BedTool with chr/start/end/value/query_idx fields."""
    return BedTool.from_dataframe(df[['chr', 'start', 'end', 'value', 'query_idx']])


def get_same_direction_region_membership_df(df_a, df_b, df_c, pair_names):
    """
    Build an UpSet membership dataframe using original regions as the counting
    unit, but only count overlaps when the effect-size direction matches.

    For each original DMR, membership in its source comparison is True when its
    own effect size has a non-zero sign. Membership in another comparison is
    True only when at least one overlapping DMR in that comparison has the same
    non-zero sign as the query DMR.
    """
    dfs = []
    for df in (df_a, df_b, df_c):
        tmp = df.reset_index(drop=True).copy()
        value_col = tmp.columns[3]
        tmp = tmp.rename(columns={value_col: 'value'})
        tmp['query_idx'] = np.arange(len(tmp), dtype=int)
        tmp['direction'] = np.sign(tmp['value'].astype(float))
        dfs.append(tmp)

    bts = [_bedtool_from_indexed_df(df) for df in dfs]
    out_dfs = []

    for i, query_df in enumerate(dfs):
        query_bt = bts[i]

        out = query_df[['chr', 'start', 'end']].copy().reset_index(drop=True)
        out.insert(0, 'source_pair', pair_names[i])
        out.insert(1, 'region_id', [f'{pair_names[i]}__{k}' for k in range(len(out))])

        query_directions = query_df['direction'].to_numpy()

        for j, target_bt in enumerate(bts):
            if i == j:
                out[pair_names[j]] = query_directions != 0
                continue

            has_same_direction_overlap = np.zeros(len(query_df), dtype=bool)

            # Retain both query and target records so we can compare effect-size
            # signs for each overlapping interval pair.
            for feature in query_bt.intersect(target_bt, wa=True, wb=True):
                fields = feature.fields
                query_value = float(fields[3])
                query_idx = int(fields[4])
                target_value = float(fields[8])

                query_sign = np.sign(query_value)
                target_sign = np.sign(target_value)
                if query_sign != 0 and query_sign == target_sign:
                    has_same_direction_overlap[query_idx] = True

            out[pair_names[j]] = has_same_direction_overlap

        out_dfs.append(out)

    membership_df = pd.concat(out_dfs, ignore_index=True)

    expected_n = sum(len(df) for df in dfs)
    if len(membership_df) != expected_n:
        raise RuntimeError(
            f'Same-direction membership dataframe has {len(membership_df)} rows, '
            f'but expected {expected_n} rows from the three input dataframes.'
        )

    return membership_df

def intersect_three_pairs(df_a, df_b, df_c, pair_names):
    """
    Intersect regions across three sample pairs using pybedtools.

    The output locus is the actual 3-way overlap:
        max(starts), min(ends)

    If multiple combinations produce the same output locus, retain the row whose
    three contributing input regions have the largest total length.

    Output columns:
        chr, start, end, <pair1 value>, <pair2 value>, <pair3 value>
    """
    a_bt = _bedtool_from_df(df_a)
    b_bt = _bedtool_from_df(df_b)
    c_bt = _bedtool_from_df(df_c)

    # First intersect pair 1 with pair 2, retaining both full original intervals.
    ab = a_bt.intersect(b_bt, wa=True, wb=True)

    # Then intersect that pairwise result with pair 3, retaining the pairwise row
    # and the full interval from pair 3.
    abc = ab.intersect(c_bt, wa=True, wb=True)

    rows = []
    for feature in abc:
        fields = feature.fields

        a_chr = fields[0]
        a_start = int(fields[1])
        a_end = int(fields[2])
        a_val = float(fields[3])

        b_chr = fields[4]
        b_start = int(fields[5])
        b_end = int(fields[6])
        b_val = float(fields[7])

        c_chr = fields[8]
        c_start = int(fields[9])
        c_end = int(fields[10])
        c_val = float(fields[11])

        start = max(a_start, b_start, c_start)
        end = min(a_end, b_end, c_end)

        if end <= start:
            continue

        rows.append({
            'chr': a_chr,
            'start': start,
            'end': end,
            pair_names[0]: a_val,
            pair_names[1]: b_val,
            pair_names[2]: c_val,
            'a_start': a_start,
            'a_end': a_end,
            'b_start': b_start,
            'b_end': b_end,
            'c_start': c_start,
            'c_end': c_end,
        })

    overlap_df = pd.DataFrame(rows)
    if overlap_df.empty:
        return pd.DataFrame(columns=['chr', 'start', 'end', *pair_names])

    overlap_df = _pick_largest_per_locus(overlap_df)
    overlap_df = overlap_df[['chr', 'start', 'end', *pair_names]]
    return overlap_df.sort_values(['chr', 'start', 'end']).reset_index(drop=True)


def filter_same_effect_direction(df, pair_names):
    """
    Keep only 3-way overlapping loci where the effect-size direction is the
    same across all three sample comparisons.

    A locus is retained when all three comparison values are positive or all
    three are negative. Loci with any zero effect size are excluded because
    their direction is undefined.
    """
    if df.empty:
        return df.copy()

    out = df.copy()
    signs = np.sign(out[pair_names].astype(float))
    same_direction = signs.nunique(axis=1).eq(1) & signs.iloc[:, 0].ne(0)
    return out.loc[same_direction].reset_index(drop=True)


if __name__ == '__main__':
    DIR1 = '/scratch/eger/projects/MethSmoothEval/tissues/modkit/dmr'
    DIR2 = '/scratch/eger/projects/MethSmoothEval/tissues/DSS'
    OUTDIR = '/scratch/eger/projects/MethSmoothEval/tissues/DM_analysis/DMRs_recurrent'
    os.makedirs(OUTDIR, exist_ok=True)

    sample_pairs = [
        'hg002_blood_v_Bulk_FC_Control_02',
        'PPMI_3404_v_HBCC_81951_FTX',
        'PPMI_3404_v_HBCC_82044_FTX'
    ]

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

    fname1 = os.path.join(DIR1, sample_pairs[0] + '.dmr.segments.txt')
    fname2 = os.path.join(DIR2, sample_pairs[0] + '_DMLtestwSmoothing_DMRs.tsv')
    fname3 = os.path.join(DIR2, sample_pairs[0] + '_DMLtest_DMRs.tsv')

    fname4 = os.path.join(DIR1, sample_pairs[1] + '.dmr.segments.txt')
    fname5 = os.path.join(DIR2, sample_pairs[1] + '_DMLtestwSmoothing_DMRs.tsv')
    fname6 = os.path.join(DIR2, sample_pairs[1] + '_DMLtest_DMRs.tsv')

    fname7 = os.path.join(DIR1, sample_pairs[2] + '.dmr.segments.txt')
    fname8 = os.path.join(DIR2, sample_pairs[2] + '_DMLtestwSmoothing_DMRs.tsv')
    fname9 = os.path.join(DIR2, sample_pairs[2] + '_DMLtest_DMRs.tsv')

    df1 = pd.read_csv(fname1, sep='\t', header=None, names=seg_cols)
    df1 = df1[df1['state_name'] == 'different']

    df2 = pd.read_csv(fname2, sep='\t')
    df2['start'] = df2['start'].astype(int)
    df2['end'] = df2['end'].astype(int)

    df3 = pd.read_csv(fname3, sep='\t')
    df3['start'] = df3['start'].astype(int)
    df3['end'] = df3['end'].astype(int)

    df4 = pd.read_csv(fname4, sep='\t', header=None, names=seg_cols)
    df4 = df4[df4['state_name'] == 'different']

    df5 = pd.read_csv(fname5, sep='\t')
    df5['start'] = df5['start'].astype(int)
    df5['end'] = df5['end'].astype(int)

    df6 = pd.read_csv(fname6, sep='\t')
    df6['start'] = df6['start'].astype(int)
    df6['end'] = df6['end'].astype(int)

    df7 = pd.read_csv(fname7, sep='\t', header=None, names=seg_cols)
    df7 = df7[df7['state_name'] == 'different']

    df8 = pd.read_csv(fname8, sep='\t')
    df8['start'] = df8['start'].astype(int)
    df8['end'] = df8['end'].astype(int)

    df9 = pd.read_csv(fname9, sep='\t')
    df9['start'] = df9['start'].astype(int)
    df9['end'] = df9['end'].astype(int)

    # Prepare BED-like dataframes for each region type.
    # Segments use effect_size; DSS DMRs use diff.Methy.
    seg_pair1 = _prepare_bed_df(df1, 'effect_size', sample_pairs[0], chrom_col='chrom')
    seg_pair2 = _prepare_bed_df(df4, 'effect_size', sample_pairs[1], chrom_col='chrom')
    seg_pair3 = _prepare_bed_df(df7, 'effect_size', sample_pairs[2], chrom_col='chrom')

    smooth_pair1 = _prepare_bed_df(df2, 'diff.Methy', sample_pairs[0], chrom_col='chr')
    smooth_pair2 = _prepare_bed_df(df5, 'diff.Methy', sample_pairs[1], chrom_col='chr')
    smooth_pair3 = _prepare_bed_df(df8, 'diff.Methy', sample_pairs[2], chrom_col='chr')

    unsmooth_pair1 = _prepare_bed_df(df3, 'diff.Methy', sample_pairs[0], chrom_col='chr')
    unsmooth_pair2 = _prepare_bed_df(df6, 'diff.Methy', sample_pairs[1], chrom_col='chr')
    unsmooth_pair3 = _prepare_bed_df(df9, 'diff.Methy', sample_pairs[2], chrom_col='chr')

    # New output dataframes: overlapping loci across the 3 sample pairs.
    segments_overlap_df = intersect_three_pairs(seg_pair1, seg_pair2, seg_pair3, sample_pairs)
    smooth_dmrs_overlap_df = intersect_three_pairs(smooth_pair1, smooth_pair2, smooth_pair3, sample_pairs)
    unsmooth_dmrs_overlap_df = intersect_three_pairs(unsmooth_pair1, unsmooth_pair2, unsmooth_pair3, sample_pairs)

    # Same-direction subsets of the 3-way overlaps.
    # These retain only loci where all three sample-pair effect sizes have the
    # same sign: all positive or all negative.
    segments_same_direction_df = filter_same_effect_direction(segments_overlap_df, sample_pairs)
    smooth_dmrs_same_direction_df = filter_same_effect_direction(smooth_dmrs_overlap_df, sample_pairs)
    unsmooth_dmrs_same_direction_df = filter_same_effect_direction(unsmooth_dmrs_overlap_df, sample_pairs)

    # Membership dataframes for UpSet plots.
    # These use original DMRs/segments as the counting unit.
    smooth_dmrs_membership_df = get_region_membership_df(smooth_pair1, smooth_pair2, smooth_pair3, sample_pairs)
    unsmooth_dmrs_membership_df = get_region_membership_df(unsmooth_pair1, unsmooth_pair2, unsmooth_pair3, sample_pairs)

    # Same-direction membership dataframes for UpSet plots.
    # These still use original DMRs as the counting unit, but an overlap only
    # counts if the overlapping DMR has the same non-zero effect-size sign as
    # the query DMR.
    smooth_dmrs_same_direction_membership_df = get_same_direction_region_membership_df(
        smooth_pair1, smooth_pair2, smooth_pair3, sample_pairs
    )
    unsmooth_dmrs_same_direction_membership_df = get_same_direction_region_membership_df(
        unsmooth_pair1, unsmooth_pair2, unsmooth_pair3, sample_pairs
    )

    # Optional file outputs.
    segments_overlap_df.to_csv(
        os.path.join(OUTDIR, 'segments_overlap_across_3_pairs.tsv'),
        sep='\t',
        index=False
    )
    smooth_dmrs_overlap_df.to_csv(
        os.path.join(OUTDIR, 'smooth_DMRs_overlap_across_3_pairs.tsv'),
        sep='\t',
        index=False
    )
    unsmooth_dmrs_overlap_df.to_csv(
        os.path.join(OUTDIR, 'unsmooth_DMRs_overlap_across_3_pairs.tsv'),
        sep='	',
        index=False
    )

    segments_same_direction_df.to_csv(
        os.path.join(OUTDIR, 'segments_overlap_across_3_pairs_same_direction.tsv'),
        sep='	',
        index=False
    )
    smooth_dmrs_same_direction_df.to_csv(
        os.path.join(OUTDIR, 'smooth_DMRs_overlap_across_3_pairs_same_direction.tsv'),
        sep='	',
        index=False
    )
    unsmooth_dmrs_same_direction_df.to_csv(
        os.path.join(OUTDIR, 'unsmooth_DMRs_overlap_across_3_pairs_same_direction.tsv'),
        sep='	',
        index=False
    )

    smooth_dmrs_membership_df.to_csv(
        os.path.join(OUTDIR, 'smooth_DMRs_upset_membership_loci.tsv'),
        sep='	',
        index=False
    )
    unsmooth_dmrs_membership_df.to_csv(
        os.path.join(OUTDIR, 'unsmooth_DMRs_upset_membership_loci.tsv'),
        sep='	',
        index=False
    )

    smooth_dmrs_same_direction_membership_df.to_csv(
        os.path.join(OUTDIR, 'smooth_DMRs_upset_membership_loci_same_direction.tsv'),
        sep='	',
        index=False
    )
    unsmooth_dmrs_same_direction_membership_df.to_csv(
        os.path.join(OUTDIR, 'unsmooth_DMRs_upset_membership_loci_same_direction.tsv'),
        sep='	',
        index=False
    )
    print('segments_overlap_df:', segments_overlap_df.shape)
    print('smooth_dmrs_overlap_df:', smooth_dmrs_overlap_df.shape)
    print('unsmooth_dmrs_overlap_df:', unsmooth_dmrs_overlap_df.shape)
    print('smooth_dmrs_membership_df:', smooth_dmrs_membership_df.shape)
    print('unsmooth_dmrs_membership_df:', unsmooth_dmrs_membership_df.shape)