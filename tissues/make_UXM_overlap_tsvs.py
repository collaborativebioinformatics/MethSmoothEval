#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# DMR files are outputs from `DSS_sample_callDMR_UXM_size.R`

# take all DMP and DMR files and extract only sites/regions overlapping with atlas regions
# covert all files into BEDs

# last two columns in all BEDs should be:
# (1) location of the overlapping atlas region 
# (2) True/False for whether the effect size direction matches the atlas region direction

import os 
import pandas as pd
import numpy as np
from typing import Optional, Mapping
from pybedtools import BedTool

def swap_modkit_dmr_pair_samples(df: pd.DataFrame, *, inplace: bool = False) -> pd.DataFrame:
    """
    Swap sample A/B fields in `modkit dmr pair` output.

    Expected columns (typical):
      chrom, start, end, name, score, strand,
      a_counts, a_total, b_counts, b_total,
      a_mod_percentages, b_mod_percentages,
      a_pct_modified, b_pct_modified,
      map_pvalue,
      effect_size, cohen_h, cohen_h_low, cohen_h_high

    Swapping A<->B implies:
      - a_* <-> b_* for counts/totals/percentages
      - effect_size flips sign (A - B becomes B - A)
      - cohen_h and its CI flip sign (Cohen's h is directional here)
      - map_pvalue unchanged (two-sided / symmetric test)

    Parameters
    ----------
    df : pd.DataFrame
        modkit dmr pair results.
    inplace : bool
        If True, mutate df in place and return it.

    Returns
    -------
    pd.DataFrame
    """
    out = df if inplace else df.copy(deep=True)

    def _swap_cols(c1: str, c2: str):
        if c1 in out.columns and c2 in out.columns:
            tmp = out[c1].copy()
            out[c1] = out[c2]
            out[c2] = tmp

    # Swap A/B columns
    _swap_cols("a_counts", "b_counts")
    _swap_cols("a_total", "b_total")
    _swap_cols("a_mod_percentages", "b_mod_percentages")
    _swap_cols("a_pct_modified", "b_pct_modified")

    # Flip signed effect columns
    if "effect_size" in out.columns:
        out["effect_size"] = -out["effect_size"]

    # Cohen's h (and CI) are directional and should flip sign.
    # Also swap low/high after sign flip so the interval remains ordered.
    if "cohen_h" in out.columns:
        out["cohen_h"] = -out["cohen_h"]

    have_ci = ("cohen_h_low" in out.columns) and ("cohen_h_high" in out.columns)
    if have_ci:
        low = out["cohen_h_low"].copy()
        high = out["cohen_h_high"].copy()
        out["cohen_h_low"] = -high
        out["cohen_h_high"] = -low

    # map_pvalue, score, strand, etc. unchanged
    return out

def swap_dss_DMLtest_samples(df: pd.DataFrame, *, inplace: bool = False) -> pd.DataFrame:
    """
    Swap group/condition 1 and 2 in a DSS::DMLtest() output dataframe.

    DMLtest() output columns typically include:
      chr, pos, mu1, mu2, diff, diff.se, stat, phi1, phi2, pval, fdr

    Where:
      - mu1, mu2 are mean methylation levels for group 1 and group 2
      - diff is (mu1 - mu2)
      - stat is the Wald test statistic, typically diff / diff.se (signed)
      - phi1, phi2 are dispersion parameters per group

    Swapping groups therefore:
      - mu1 <-> mu2
      - phi1 <-> phi2
      - diff -> -diff
      - stat -> -stat
      - diff.se stays the same (a standard error, non-directional)
      - pval and fdr stay the same (two-sided test)

    Parameters
    ----------
    df : pd.DataFrame
        DMLtest() results as a pandas dataframe.
    inplace : bool
        If True, mutate df in place and return it. Otherwise return a copy.

    Returns
    -------
    pd.DataFrame
    """
    out = df if inplace else df.copy(deep=True)

    def _swap_cols(c1: str, c2: str):
        if c1 in out.columns and c2 in out.columns:
            tmp = out[c1].copy()
            out[c1] = out[c2]
            out[c2] = tmp

    # Swap group-specific estimates
    _swap_cols("mu1", "mu2")
    _swap_cols("phi1", "phi2")

    # Flip signed contrasts/statistics
    if "diff" in out.columns:
        out["diff"] = -out["diff"]

    if "stat" in out.columns:
        out["stat"] = -out["stat"]

    # diff.se, pval, fdr unchanged (if present)
    return out

def swap_modkit_segments_samples(
    df: pd.DataFrame,
    *,
    inplace: bool = False,
    state_name_map: Optional[Mapping[str, str]] = None,
) -> pd.DataFrame:
    """
    Swap sample A and B columns in a modkit --segment segments dataframe.

    Expected columns (subset ok):
      - a_counts <-> b_counts
      - a_mod_percentages <-> b_mod_percentages
      - a_frac_modified <-> b_frac_modified
      - effect_size sign flips (assumed A-B or directionally tied to A vs B)
      - cohen_h sign flips; cohen_h_low/high flip + sign-adjust so low <= high

    Other columns (chrom/start/end/state_name/score/n_sites, etc.) are left unchanged
    unless state_name_map is provided.

    Parameters
    ----------
    df : pd.DataFrame
        Input segments dataframe.
    inplace : bool
        If True, mutate df in place and return it. Otherwise return a copy.
    state_name_map : dict-like, optional
        If your `state_name` encodes direction (e.g., "A_hyper" vs "B_hyper"),
        provide a mapping to swap labels. Unmapped values are left as-is.

    Returns
    -------
    pd.DataFrame
    """
    out = df if inplace else df.copy(deep=True)

    def _swap_cols(c1: str, c2: str):
        if c1 in out.columns and c2 in out.columns:
            tmp = out[c1].copy()
            out[c1] = out[c2]
            out[c2] = tmp

    # Swap A/B payload columns
    _swap_cols("a_counts", "b_counts")
    _swap_cols("a_mod_percentages", "b_mod_percentages")
    _swap_cols("a_frac_modified", "b_frac_modified")

    # Flip directional summary stats
    if "effect_size" in out.columns:
        out["effect_size"] = -out["effect_size"]

    if "cohen_h" in out.columns:
        out["cohen_h"] = -out["cohen_h"]

    if "cohen_h_low" in out.columns and "cohen_h_high" in out.columns:
        # When you negate an interval [low, high], it becomes [-high, -low]
        low = out["cohen_h_low"]
        high = out["cohen_h_high"]
        out["cohen_h_low"] = -high
        out["cohen_h_high"] = -low

        # Optional: enforce ordering if any rows ended up inverted due to weird inputs
        bad = out["cohen_h_low"] > out["cohen_h_high"]
        if bad.any():
            tmp = out.loc[bad, "cohen_h_low"].copy()
            out.loc[bad, "cohen_h_low"] = out.loc[bad, "cohen_h_high"]
            out.loc[bad, "cohen_h_high"] = tmp

    # Optionally swap directional state labels
    if state_name_map is not None and "state_name" in out.columns:
        out["state_name"] = out["state_name"].map(lambda x: state_name_map.get(x, x))

    return out

def swap_dss_callDMR_samples(df: pd.DataFrame, *, inplace: bool = False) -> pd.DataFrame:
    """
    Swap group/condition 1 and 2 in a DSS::callDMR() output dataframe.

    callDMR() typically returns columns:
      chr, start, end, length, nCG, meanMethy1, meanMethy2, diff.Methy, areaStat
    where diff.Methy is (group1 - group2) and areaStat is the (signed) sum of locus test stats.
    Swapping groups therefore:
      - meanMethy1 <-> meanMethy2
      - diff.Methy  -> -diff.Methy
      - areaStat    -> -areaStat

    Notes:
      - If a column like 'postprob.overThreshold' is present (can appear in some DSS outputs
        when delta > 0 at the DML level), we *cannot* transform it correctly without
        recomputing from the underlying per-locus results, so we leave it unchanged.

    Parameters
    ----------
    df : pd.DataFrame
        callDMR() output as a pandas dataframe.
    inplace : bool
        If True, mutate df in place and return it. Otherwise return a copy.

    Returns
    -------
    pd.DataFrame
    """
    out = df if inplace else df.copy(deep=True)

    def _swap_cols(c1: str, c2: str):
        if c1 in out.columns and c2 in out.columns:
            tmp = out[c1].copy()
            out[c1] = out[c2]
            out[c2] = tmp

    # Swap group means
    _swap_cols("meanMethy1", "meanMethy2")

    # Flip direction of group1-group2 quantities
    if "diff.Methy" in out.columns:
        out["diff.Methy"] = -out["diff.Methy"]

    if "areaStat" in out.columns:
        out["areaStat"] = -out["areaStat"]

    return out

def direction_matches_atlas(effect_size, best_group, best_dir, eps=0.0):
    """
    Returns True/False if the sign of (brain - blood) effect_size matches the atlas direction,
    accounting for whether best_group is a brain or blood cell type.

    best_dir:
      - "hypo"  => best_group < rest
      - "hyper" => best_group > rest

    eps: treat |effect_size| <= eps as ambiguous (returns np.nan)
    """
    if pd.isna(effect_size) or pd.isna(best_group) or pd.isna(best_dir):
        return np.nan

    if abs(effect_size) <= eps:
        return np.nan

    # +1 means we expect effect_size > 0, -1 means expect effect_size < 0
    if best_group in BRAIN_TYPES:
        expected = -1 if best_dir == "hypo" else +1
    elif best_group in BLOOD_TYPES:
        expected = +1 if best_dir == "hypo" else -1
    else:
        # unknown / other cell type label
        return np.nan

    observed = +1 if effect_size > 0 else -1
    return observed == expected

def subset_atlas_overlap(df):
    # adds 3 additional columns from atlas file
    bt = BedTool.from_dataframe(df[['chr', 'start', 'end']])
    
    bt_atlas = bt.intersect(atlas_bt, wa=True, wb=True)
    df_atlas = bt_atlas.to_dataframe(disable_auto_names=True, header=None)
    df_atlas = df_atlas[[0, 1, 2, 6, 7, 8]]
    df_atlas.columns = ['chr', 'start', 'end', 'best_group', 'best_dir', 'location']

    # merge all single site info into single dataframe
    df_atlas = df.merge(df_atlas, on=['chr', 'start', 'end'], how='inner')

    # determine if direction matches
    if "effect_size" in df_atlas.columns:
        df_atlas["dir_match_atlas"] = df_atlas.apply(
            lambda r: direction_matches_atlas(r["effect_size"], r["best_group"], r["best_dir"]),
            axis=1
        )
    else:
        df_atlas["dir_match_atlas"] = df_atlas.apply(
            lambda r: direction_matches_atlas(r["diff.Methy"], r["best_group"], r["best_dir"]),
            axis=1
        )
    return df_atlas

if __name__ == '__main__':

    dmr_cols = [
        'chr',
        'start',
        'end',
        'name',
        'score',
        'strand',
        'a_counts',
        'a_total',
        'b_counts',
        'b_total',
        'a_mod_percentages',
        'b_mod_percentages',
        'a_pct_modified',
        'b_pct_modified',
        'map_pvalue',
        'effect_size',
        'cohen_h', 
        'cohen_h_low',
        'cohen_h_high'
    ]
    
    seg_cols = [
        'chr',
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
    
    REF_DIR = '/home/eger/software/UXM_deconv/supplemental'
    DIR1 = '/scratch/eger/projects/MethSmoothEval/tissues/modkit/dmr'
    DIR2 = '/scratch/eger/projects/MethSmoothEval/tissues/DSS'
    DIR3 = '/scratch/eger/projects/MethSmoothEval/tissues/DM_analysis/UXM/DSS_DMRs'
    
    OUTDIR = '/scratch/eger/projects/MethSmoothEval/tissues/DM_analysis/UXM'

    # load atlas & convert to BED
    BRAIN_TYPES = {"Neuron", "Oligodend", "Neuron:Oligodend"}
    BLOOD_TYPES = {"Blood-B", "Blood-T", "Blood-NK", "Blood-Granul", "Blood-Mono+Macro"}

    # atlas_name = 'U25'
    # fname0 = os.path.join(REF_DIR, 'Atlas.'+atlas_name+'.Brain_Blood_lifted_to_hg38.tsv') # 199 regions
    atlas_name = 'U250'
    fname0 = os.path.join(REF_DIR, 'Atlas.'+atlas_name+'.Brain_Blood_lifted_to_hg38.tsv') # 199 regions
    
    df0 = pd.read_csv(fname0, sep='\t')
    atlas_bt = BedTool.from_dataframe(df0[['chr', 'start', 'end', 'best_group', 'best_dir', 'location']])

    # 3 comparisons
    sample_pairs = [
        "Bulk_FC_Control_02_v_hg002_blood",
        "HBCC_81951_FTX_v_PPMI_3404",
        "HBCC_82044_FTX_v_PPMI_3404"]

    ##############################################################################
    # DMPs and DMRs
    for pair in sample_pairs:        
        if pair == "Bulk_FC_Control_02_v_hg002_blood":
            # single-sites
            fname1 = os.path.join(DIR1, pair+'.dmr.bed.gz')
            fname2 = os.path.join(DIR2, pair+'_DMLtest.tsv.gz')
            fname3 = os.path.join(DIR2, pair+'_DMLtestwSmoothing.tsv.gz')
            
            # regions
            fname4 = os.path.join(DIR1, pair+'.dmr.segments.txt')
            fname5 = os.path.join(DIR2, pair+'_DMLtest_DMRs.tsv')
            fname6 = os.path.join(DIR2, pair+'_DMLtestwSmoothing_DMRs.tsv')

            # load the files
            df1 = pd.read_csv(fname1, sep='\t', header=None, names=dmr_cols)
            df2 = pd.read_csv(fname2, sep='\t')
            df3 = pd.read_csv(fname3, sep='\t')

            df4 = pd.read_csv(fname4, sep='\t', header=None, names=seg_cols)
            df5 = pd.read_csv(fname5, sep='\t')
            df6 = pd.read_csv(fname6, sep='\t')
        
        else:
            # get the original names
            samples = pair.split('_v_')
            orig_pair = samples[1]+'_v_'+samples[0]
            
            # single-sites
            fname1 = os.path.join(DIR1, orig_pair+'.dmr.bed.gz')
            fname2 = os.path.join(DIR2, orig_pair+'_DMLtest.tsv.gz')
            fname3 = os.path.join(DIR2, orig_pair+'_DMLtestwSmoothing.tsv.gz')
            
            # regions
            fname4 = os.path.join(DIR1, orig_pair+'.dmr.segments.txt')
            fname5 = os.path.join(DIR2, orig_pair+'_DMLtest_DMRs.tsv')
            fname6 = os.path.join(DIR2, orig_pair+'_DMLtestwSmoothing_DMRs.tsv')

            # load the files
            df1 = swap_modkit_dmr_pair_samples(pd.read_csv(fname1, sep='\t', header=None, names=dmr_cols))
            df2 = swap_dss_DMLtest_samples(pd.read_csv(fname2, sep='\t'))
            df3 = swap_dss_DMLtest_samples(pd.read_csv(fname3, sep='\t'))

            df4 = swap_modkit_segments_samples(pd.read_csv(fname4, sep='\t', header=None, names=seg_cols))
            df5 = swap_dss_callDMR_samples(pd.read_csv(fname5, sep='\t'))
            df6 = swap_dss_callDMR_samples(pd.read_csv(fname6, sep='\t'))

            """
            # swap the sample order for all files & save them
            # single-sites
            new_fname1 = os.path.join(DIR1, pair+'.dmr.bed.gz')
            new_fname2 = os.path.join(DIR2, pair+'_DMLtest.tsv.gz')
            new_fname3 = os.path.join(DIR2, pair+'_DMLtestwSmoothing.tsv.gz')
            
            # regions
            new_fname4 = os.path.join(DIR1, pair+'.dmr.segments.txt')
            new_fname5 = os.path.join(DIR2, pair+'_DMLtest_DMRs.tsv')
            new_fname6 = os.path.join(DIR2, pair+'_DMLtestwSmoothing_DMRs.tsv')

            df1.to_csv(new_fname1, sep='\t', header=False, index=False, compression='gzip')
            df2.to_csv(new_fname2, sep='\t', index=False, compression='gzip')
            df3.to_csv(new_fname3, sep='\t', index=False, compression='gzip')

            df4.to_csv(new_fname4, sep='\t', header=False, index=False)
            df5.to_csv(new_fname5, sep='\t', index=False)
            df6.to_csv(new_fname6, sep='\t', index=False)
            """
        
        ## single sites        
        df1_atlas = subset_atlas_overlap(df1)

        # change col names for DMLtest outputs for merging
        df2.columns = ['chr', 'end', 'mu1', 'mu2', 'diff', 'diff.se', 'DML_stat', 
                       'phi1', 'phi2', 'DML_pval', 'DML_fdr']
        df3.columns = ['chr', 'end', 'mu1', 'mu2', 'diff', 'diff.se', 'Smooth_DML_stat', 
                       'phi1', 'phi2', 'Smooth_DML_pval', 'Smooth_DML_fdr']
        df1_atlas = df1_atlas.merge(df2[['chr', 'end', 'DML_stat', 'DML_pval', 'DML_fdr']],
                                    on=['chr', 'end'], how='inner').merge(
                                    df3[['chr', 'end', 'Smooth_DML_stat', 'Smooth_DML_pval', 'Smooth_DML_fdr']],
                                    on=['chr', 'end'], how='inner')
        
        # save as .tsv
        out1 = os.path.join(OUTDIR, pair+'.all_CpGs_in_'+atlas_name+'.tsv')
        df1_atlas.to_csv(out1, sep='\t', index=False)

        ## Segments & DMRs
        # DMRs - subtract 1 from start to make the same as segments
        df5['start'] = df5['start'] - 1
        df6['start'] = df6['start'] - 1
        df5['end'] = df5['end'].astype(int)
        df6['end'] = df6['end'].astype(int)

        df4_atlas = subset_atlas_overlap(df4)
        df5_atlas = subset_atlas_overlap(df5)
        df6_atlas = subset_atlas_overlap(df6)
        
        out4 = os.path.join(OUTDIR, pair+'.segments_in_'+atlas_name+'.tsv')
        out5 = os.path.join(OUTDIR, pair+'.DMLtest_DMRs_in_'+atlas_name+'.tsv')
        out6 = os.path.join(OUTDIR, pair+'.DMLtestwSmoothing_DMRs_in_'+atlas_name+'.tsv')
        
        df4_atlas.to_csv(out4, sep='\t', index=False)
        df5_atlas.to_csv(out5, sep='\t', index=False)
        df6_atlas.to_csv(out6, sep='\t', index=False)

        print('\n'.join([out1, out4, out5, out6]))

    """
    ##############################################################################
    # DMRs only
    for pair in sample_pairs:        
        if pair == "Bulk_FC_Control_02_v_hg002_blood":
            # regions
            fname4 = os.path.join(DIR1, pair+'.dmr.segments.txt')
            fname5 = os.path.join(DIR3, pair+'_DMLtest_DMRs.tsv')
            fname6 = os.path.join(DIR3, pair+'_DMLtestwSmoothing_DMRs.tsv')

            df4 = pd.read_csv(fname4, sep='\t', header=None, names=seg_cols)
            df5 = pd.read_csv(fname5, sep='\t')
            df6 = pd.read_csv(fname6, sep='\t')
        
        else:
            # get the original names
            samples = pair.split('_v_')
            orig_pair = samples[1]+'_v_'+samples[0]
            
            # regions
            fname4 = os.path.join(DIR1, orig_pair+'.dmr.segments.txt')
            fname5 = os.path.join(DIR3, orig_pair+'_DMLtest_DMRs.tsv')
            fname6 = os.path.join(DIR3, orig_pair+'_DMLtestwSmoothing_DMRs.tsv')

            # load the files
            df4 = swap_modkit_segments_samples(pd.read_csv(fname4, sep='\t', header=None, names=seg_cols))
            df5 = swap_dss_callDMR_samples(pd.read_csv(fname5, sep='\t'))
            df6 = swap_dss_callDMR_samples(pd.read_csv(fname6, sep='\t'))

        ## Segments & DMRs
        # DMRs - subtract 1 from start to make the same as segments
        df5['start'] = df5['start'] - 1
        df6['start'] = df6['start'] - 1
        df5['end'] = df5['end'].astype(int)
        df6['end'] = df6['end'].astype(int)

        df4_atlas = subset_atlas_overlap(df4)
        df5_atlas = subset_atlas_overlap(df5)
        df6_atlas = subset_atlas_overlap(df6)
        
        out4 = os.path.join(DIR3, pair+'.segments_in_atlas.tsv')
        out5 = os.path.join(DIR3, pair+'.DMLtest_DMRs_in_atlas.tsv')
        out6 = os.path.join(DIR3, pair+'.DMLtestwSmoothing_DMRs_in_atlas.tsv')
        
        df4_atlas.to_csv(out4, sep='\t', index=False)
        df5_atlas.to_csv(out5, sep='\t', index=False)
        df6_atlas.to_csv(out6, sep='\t', index=False)
        """


