#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# use output from ~/scripts/MethSmoothEval/tissues/call_DM/find_recurrent_DMRs.py

import os
import warnings
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from pybedtools import BedTool

warnings.filterwarnings('ignore', category=UserWarning)


def clean_bed_df(df):
    """Return a clean 3-column BED dataframe."""
    out = df.iloc[:, :3].copy()
    out.columns = ['chrom', 'start', 'end']
    out = out.dropna().copy()
    out['chrom'] = out['chrom'].astype(str)
    out['start'] = out['start'].astype(int)
    out['end'] = out['end'].astype(int)
    out = out[out['end'] > out['start']].drop_duplicates().reset_index(drop=True)
    return out


def df_to_bt(df):
    """Convert dataframe to sorted BedTool."""
    if df.empty:
        return BedTool('', from_string=True)
    return BedTool.from_dataframe(clean_bed_df(df)).sort()


def bt_to_df(bt):
    """Convert BedTool to clean dataframe."""
    rows = [x.fields[:3] for x in bt]
    if len(rows) == 0:
        return pd.DataFrame(columns=['chrom', 'start', 'end'])
    df = pd.DataFrame(rows, columns=['chrom', 'start', 'end'])
    return clean_bed_df(df)


def count_regions(bt):
    return len(bt_to_df(bt))


def intersect_keep_larger(bt_a, bt_b):
    """
    Intersect two BED sets and, for each overlapping pair, keep the larger
    of the two original intervals.
    """
    a_df = bt_to_df(bt_a)
    b_df = bt_to_df(bt_b)
    if a_df.empty or b_df.empty:
        return BedTool('', from_string=True)

    a_df = a_df.copy().reset_index(drop=True)
    b_df = b_df.copy().reset_index(drop=True)
    a_df['a_id'] = ['A{}'.format(i) for i in range(len(a_df))]
    b_df['b_id'] = ['B{}'.format(i) for i in range(len(b_df))]
    a_df['a_len'] = a_df['end'] - a_df['start']
    b_df['b_len'] = b_df['end'] - b_df['start']

    a_bt = BedTool.from_dataframe(a_df[['chrom', 'start', 'end', 'a_id', 'a_len']]).sort()
    b_bt = BedTool.from_dataframe(b_df[['chrom', 'start', 'end', 'b_id', 'b_len']]).sort()
    inter = a_bt.intersect(b_bt, wa=True, wb=True)

    kept = []
    seen = set()
    for rec in inter:
        f = rec.fields
        a_key = (f[0], int(f[1]), int(f[2]), f[3], int(f[4]))
        b_key = (f[5], int(f[6]), int(f[7]), f[8], int(f[9]))

        if a_key[4] >= b_key[4]:
            chosen = a_key[:3]
        else:
            chosen = b_key[:3]

        if chosen not in seen:
            seen.add(chosen)
            kept.append(chosen)

    if len(kept) == 0:
        return BedTool('', from_string=True)

    kept_df = pd.DataFrame(kept, columns=['chrom', 'start', 'end'])
    return df_to_bt(kept_df)


def recurrent_bed_from_pairs(bed_paths):
    """Recurrent DMRs = intervals with overlap support across all 3 pairwise DMR files."""
    current = BedTool(bed_paths[0]).sort()
    for p in bed_paths[1:]:
        current = intersect_keep_larger(current, BedTool(p).sort())
    return df_to_bt(bt_to_df(current))

def subset_overlap(bt_a, bt_b):
    """Regions from A that overlap B at least once."""
    a_df = bt_to_df(bt_a)
    b_df = bt_to_df(bt_b)
    if a_df.empty or b_df.empty:
        return BedTool('', from_string=True)
    return df_to_bt(bt_to_df(bt_a.intersect(bt_b, u=True)))


def subset_nonoverlap(bt_a, bt_b):
    """Regions from A that do not overlap B."""
    a_df = bt_to_df(bt_a)
    if a_df.empty:
        return BedTool('', from_string=True)
    b_df = bt_to_df(bt_b)
    if b_df.empty:
        return df_to_bt(a_df)

    # get the intersection and then select all not in it
    inter_bt = bt_a.intersect(bt_b, u=True)
        
    return df_to_bt(bt_to_df(bt_a.intersect(inter_bt, v=True)))


def percent_overlap_by_class(region_bt, ccre_df, class_name):
    """Percent of regions in region_bt overlapping >=1 cCRE of a given class."""
    region_df = bt_to_df(region_bt)
    total = len(region_df)
    if total == 0:
        return np.nan

    class_df = ccre_df[ccre_df['class'] == class_name][['chrom', 'start', 'end']].copy()
    class_bt = df_to_bt(class_df)
    n_overlap = count_regions(region_bt.intersect(class_bt, u=True))
    return 100.0 * n_overlap / total


def build_matrix(region_sets, ccre_df, classes_of_interest):
    matrix = pd.DataFrame(index=classes_of_interest, columns=list(region_sets.keys()), dtype=float)
    for col_name, bt in region_sets.items():
        for class_name in classes_of_interest:
            matrix.loc[class_name, col_name] = percent_overlap_by_class(bt, ccre_df, class_name)
    return matrix


def columnwise_minmax_norm(matrix):
    """Normalize each column independently to [0, 1] for coloring only."""
    normed = pd.DataFrame(index=matrix.index, columns=matrix.columns, dtype=float)
    for col in matrix.columns:
        vals = matrix[col].astype(float)
        valid = vals.dropna()
        if valid.empty:
            normed[col] = np.nan
            continue
        vmin = valid.min()
        vmax = valid.max()
        if np.isclose(vmin, vmax):
            normed[col] = 0.5
            continue
        norm = Normalize(vmin=vmin, vmax=vmax)
        normed[col] = vals.map(lambda x: norm(x) if pd.notna(x) else np.nan)
    return normed


def plot_matrix(matrix, out_png, title):
    display_vals = matrix.copy()
    color_vals = columnwise_minmax_norm(display_vals)

    nrows, ncols = display_vals.shape
    fig, ax = plt.subplots(figsize=(6, 4), dpi=200)

    im = ax.imshow(color_vals.values, cmap='RdBu_r', vmin=0, vmax=1, aspect='auto')

    ax.set_xticks(np.arange(ncols))
    ax.set_yticks(np.arange(nrows))
    ax.set_xticklabels(display_vals.columns, rotation=35, ha='right', fontsize=10)
    ax.set_yticklabels(display_vals.index, fontsize=11)
    ax.set_title(title, fontsize=14, pad=16)

    ax.set_xticks(np.arange(-0.5, ncols, 1), minor=True)
    ax.set_yticks(np.arange(-0.5, nrows, 1), minor=True)
    ax.grid(which='minor', color='black', linestyle='-', linewidth=1)
    ax.tick_params(which='minor', bottom=False, left=False)

    for i in range(nrows):
        for j in range(ncols):
            val = display_vals.iloc[i, j]
            txt = 'NA' if pd.isna(val) else '{:.1f}%'.format(val)

            cval = color_vals.iloc[i, j]
            # white text at extremes, black in middle
            if pd.isna(cval):
                text_color = 'black'
            elif cval <= 0.2 or cval >= 0.8:
                text_color = 'white'
            else:
                text_color = 'black'

            ax.text(j, i, txt, ha='center', va='center', color=text_color, fontsize=10, fontweight='bold')

    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label('% regions overlapping cCRE class\n (per column normalized)', fontsize=10)

    plt.tight_layout()
    plt.savefig(out_png, bbox_inches='tight')
    plt.close(fig)

def print_set_sizes(label, region_sets):
    print('{} region counts:'.format(label))
    for name, bt in region_sets.items():
        print('  {:<22s} {:,}'.format(name, count_regions(bt)))

if __name__ == '__main__':
    REF_DIR1 = '/scratch/eger/projects/MethSmoothEval/tissues/ref_data/ENCODE_cCREs'
    DIR1 = '/scratch/eger/projects/MethSmoothEval/tissues/DSS'
    DIR2 = '/scratch/eger/projects/MethSmoothEval/tissues/DM_analysis/DMRs_recurrent'
    DIR3 = '/scratch/eger/projects/MethSmoothEval/tissues/modkit/entropy'
    OUT_DIR = os.getcwd()

    ### recurrent DMRs ###
    unsmooth_file = os.path.join(DIR2, 'unsmooth_DMRs_overlap_across_3_pairs_same_direction.tsv')
    smooth_file = os.path.join(DIR2, 'smooth_DMRs_overlap_across_3_pairs_same_direction.tsv')

    # make bed files
    df_unsmooth = pd.read_csv(unsmooth_file, sep='\t')
    df_smooth = pd.read_csv(smooth_file, sep='\t')
    
    unsmooth_bed = os.path.join(DIR2, 'unsmooth_DMRs_overlap_across_3_pairs_same_direction.bed')
    smooth_bed = os.path.join(DIR2, 'smooth_DMRs_overlap_across_3_pairs_same_direction.bed')

    df_unsmooth.to_csv(unsmooth_bed, sep='\t', header=False, index=False)
    df_smooth.to_csv(smooth_bed, sep='\t', header=False, index=False)
    
    unsmooth_bt = BedTool(unsmooth_bed)  # 4123
    smooth_bt = BedTool(smooth_bed) # 64465

    DMR_region_sets = {
        'unsmoothed ∩ smoothed': subset_overlap(unsmooth_bt, smooth_bt),
        'smoothed ∩ unsmoothed': subset_overlap(smooth_bt, unsmooth_bt),
        'smoothed only': subset_nonoverlap(smooth_bt, unsmooth_bt),
        'unsmoothed only': subset_nonoverlap(unsmooth_bt, smooth_bt)}
    
    for name, bt in DMR_region_sets.items():
        print('  {:<22s} {:,}'.format(name, len(bt)))
      # unsmoothed ∩ smoothed  4,065
      # smoothed ∩ unsmoothed  3,507
      # smoothed only          60,958
      # unsmoothed only        58
    
    ### tissue-specific high entropy regions ###
    ent_tsv = os.path.join(DIR3, 'tissue_entropy_regions_25_75_DMR_overlaps.tsv')
    ent_df = pd.read_csv(ent_tsv, sep='\t')
    ent_df['class'] = ent_df['high_tissue']
    
    for name, count in dict(ent_df['class'].value_counts()).items():
        print('  {:<22s} {:,}'.format(name, count))
      # brain                  3,837
      # blood                  364

    
    ### Define sets of cCREs ###
    classes_of_interest = ['PLS', 'pELS', 'dELS', 'CTCF-bound']
    ccre_bed = os.path.join(REF_DIR1, 'GRCh38-cCREs.bed')
    ccre_df = pd.read_csv(ccre_bed, sep='\t', header=None, 
                          names=['chrom', 'start', 'end', 'id1', 'id2', 'class'])

    
    # keep only ccres in those classes
    # ordered by priority, make sure sets are unique
    ccre_dict = {}
    for i, ccre_class in enumerate(classes_of_interest):
        class_bed = os.path.join(REF_DIR1, 'GRCh38-cCREs.'+ccre_class+'.bed')
        df = pd.read_csv(class_bed, sep='\t', header=None, 
                          names=['chrom', 'start', 'end', 'id1', 'id2', 'class'])
        if i == 0:
            ccre_dict[ccre_class] = df['id1'].to_list()
        elif i == 1:
            ccre_dict[ccre_class] = list(set(df['id1']) - set(ccre_dict['PLS']))
        elif i == 2:
            ccre_dict[ccre_class] = list(set(df['id1']) - set(ccre_dict['PLS']) - set(ccre_dict['pELS']))
        else:
            ccre_dict[ccre_class] = list(set(df['id1']) - set(ccre_dict['PLS']) - 
                                         set(ccre_dict['pELS']) - set(ccre_dict['dELS']))
    all_ccres = []
    for ccre_class in ccre_dict.keys():
        all_ccres += ccre_dict[ccre_class]
        
    ccre_df = ccre_df[ccre_df['id1'].isin(all_ccres)].reset_index(drop=True)

    # change to CTCF-bound
    ccre_df.loc[ccre_df['class'].isin(['CA-CTCF', 'TF', 'CA-H3K4me3']), 'class'] = 'CTCF-bound'
    print(' '.join(classes_of_interest), 'cCRE regions: {:,}'.format(len(ccre_df)))
    # PLS pELS dELS CTCF-bound cCRE regions: 1,958,107

    for name, count in dict(ccre_df['class'].value_counts()).items():
        print('  {:<22s} {:,}'.format(name, count))

    # dELS                   1,469,205
    # pELS                   249,464
    # CTCF-bound             191,906
    # PLS                    47,532

    ### recurrent DMRs & cCREs overlap ###
    # create the matrix
    mat = build_matrix(DMR_region_sets, ccre_df, classes_of_interest)

    png = os.path.join(OUT_DIR, 'recurrent_DMR_cCRE_matrix.png')
    plot_matrix(mat, png, 'recurrent DMRs overlap with cCRE classes')

    print('Saved:')
    print('  {}'.format(png))

    print('Raw percent matrices:')
    print('{}'.format(mat.round(2)))


    ### recurrent DMRs & tissue-specific high entropy regions overlap ###
    # create the matrix
    mat = build_matrix(DMR_region_sets, ent_df, ['blood', 'brain'])

    png = os.path.join(OUT_DIR, 'recurrent_DMR_high_entropy_matrix.png')
    plot_matrix(mat, png, 'recurrent DMRs overlap with high entropy regions')

    print('Saved:')
    print('  {}'.format(png))

    print('Raw percent matrices:')
    print('{}'.format(mat.round(2)))
