#!/usr/bin/env python3
"""
BigWig comparison script that handles sparse data properly
"""

import pyBigWig
import matplotlib.pyplot as plt
import numpy as np
import argparse

def compare_pileup_modes_fixed(model_bw, count_bw, chrom, start, end, output_prefix):
    """Compare model vs count modes with proper sparse data handling"""
    
    print(f"Comparing {chrom}:{start:,}-{end:,}")
    
    # Open BigWig files
    bw_model = pyBigWig.open(model_bw)
    bw_count = pyBigWig.open(count_bw)
    
    # Get intervals instead of values for sparse data
    print("Extracting intervals...")
    model_intervals = bw_model.intervals(chrom, start, end)
    count_intervals = bw_count.intervals(chrom, start, end)
    
    print(f"Model intervals: {len(model_intervals)}")
    print(f"Count intervals: {len(count_intervals)}")
    
    # Convert intervals to positions and values
    model_positions = []
    model_values = []
    for interval in model_intervals:
        pos, end_pos, value = interval
        model_positions.append(pos)
        model_values.append(value)
    
    count_positions = []
    count_values = []
    for interval in count_intervals:
        pos, end_pos, value = interval
        count_positions.append(pos)
        count_values.append(value)
    
    # Create the comparison plot
    plt.figure(figsize=(15, 12))
    
    # Top plot: Line comparison using scatter for sparse data
    plt.subplot(3, 1, 1)
    plt.scatter(model_positions, model_values, label='Model (Smoothed)', alpha=0.7, s=1, color='blue')
    plt.scatter(count_positions, count_values, label='Count (Unsmoothed)', alpha=0.5, s=1, color='orange')
    plt.ylabel('Methylation Frequency (%)')
    plt.title(f'Pileup Mode Comparison: {chrom}:{start:,}-{end:,}')
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.xlim(start, end)
    
    # Create matching datasets for difference calculation
    # Find common positions (approximately)
    model_dict = {pos: val for pos, val in zip(model_positions, model_values)}
    count_dict = {pos: val for pos, val in zip(count_positions, count_values)}
    
    # Find positions that exist in both datasets
    common_positions = set(model_positions) & set(count_positions)
    
    if common_positions:
        common_pos_list = sorted(list(common_positions))
        common_model_vals = [model_dict[pos] for pos in common_pos_list]
        common_count_vals = [count_dict[pos] for pos in common_pos_list]
        
        # Middle plot: Difference at common positions
        plt.subplot(3, 1, 2)
        differences = np.array(common_model_vals) - np.array(common_count_vals)
        plt.scatter(common_pos_list, differences, alpha=0.6, s=1, color='red')
        plt.axhline(y=0, color='black', linestyle='--', alpha=0.5)
        plt.ylabel('Difference\n(Model - Count)')
        plt.title(f'Smoothing Effect at {len(common_positions)} Common Positions')
        plt.grid(True, alpha=0.3)
        plt.xlim(start, end)
        
        # Add difference statistics
        plt.text(0.02, 0.98, 
                f'Mean diff: {np.mean(differences):.2f}\nStd diff: {np.std(differences):.2f}', 
                transform=plt.gca().transAxes, 
                verticalalignment='top',
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
        
        # Bottom plot: Correlation using common positions
        plt.subplot(3, 1, 3)
        plt.scatter(common_count_vals, common_model_vals, alpha=0.5, s=1)
        plt.xlabel('Count Mode (Unsmoothed) %')
        plt.ylabel('Model Mode (Smoothed) %')
        plt.title('Model vs Count Mode Correlation')
        
        # Add correlation coefficient
        corr = np.corrcoef(common_model_vals, common_count_vals)[0,1]
        plt.text(0.05, 0.95, f'r = {corr:.3f}', transform=plt.gca().transAxes, 
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
        
        # Add diagonal line
        min_val = min(min(common_model_vals), min(common_count_vals))
        max_val = max(max(common_model_vals), max(common_count_vals))
        plt.plot([min_val, max_val], [min_val, max_val], 'r--', alpha=0.7)
        
        print(f"Common positions: {len(common_positions)}")
        print(f"Correlation: {corr:.3f}")
        print(f"Mean difference: {np.mean(differences):.2f}")
        
    else:
        print("Warning: No exactly matching positions found between model and count data")
        
        # Still do correlation with all data
        plt.subplot(3, 1, 3)
        plt.scatter(count_values, model_values, alpha=0.5, s=1)
        plt.xlabel('Count Mode (Unsmoothed) %')
        plt.ylabel('Model Mode (Smoothed) %')
        plt.title('Model vs Count Mode Correlation (All Data)')
        
        # Correlation using all data
        corr = np.corrcoef(model_values, count_values)[0,1]
        plt.text(0.05, 0.95, f'r = {corr:.3f}', transform=plt.gca().transAxes, 
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    plt.tight_layout()
    plt.savefig(f'{output_prefix}_mode_comparison.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # Create a summary plot showing data density
    plt.figure(figsize=(15, 8))
    
    # Data density plot
    plt.subplot(2, 1, 1)
    # Create histogram of positions to show data density
    bin_size = (end - start) // 200  # 200 bins across region
    bins = range(start, end + bin_size, bin_size)
    
    model_hist, _ = np.histogram(model_positions, bins=bins)
    count_hist, _ = np.histogram(count_positions, bins=bins)
    
    bin_centers = [(bins[i] + bins[i+1]) / 2 for i in range(len(bins)-1)]
    
    plt.plot(bin_centers, model_hist, label='Model data density', alpha=0.8, color='blue')
    plt.plot(bin_centers, count_hist, label='Count data density', alpha=0.8, color='orange')
    plt.ylabel('Data Points per Bin')
    plt.title('Data Density Across Region')
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    # Value distribution comparison
    plt.subplot(2, 1, 2)
    plt.hist(model_values, bins=50, alpha=0.6, label='Model values', density=True, color='blue')
    plt.hist(count_values, bins=50, alpha=0.6, label='Count values', density=True, color='orange')
    plt.xlabel('Methylation Frequency (%)')
    plt.ylabel('Density')
    plt.title('Value Distribution Comparison')
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(f'{output_prefix}_data_summary.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # Close BigWig files
    bw_model.close()
    bw_count.close()
    
    print(f"Plots saved:")
    print(f"  - {output_prefix}_mode_comparison.png")
    print(f"  - {output_prefix}_data_summary.png")

def main():
    parser = argparse.ArgumentParser(description='pb-cpg-tools pileup mode comparison')
    parser.add_argument('--model', required=True, help='Model mode BigWig file')
    parser.add_argument('--count', required=True, help='Count mode BigWig file')
    parser.add_argument('--region', required=True, help='Region to analyze (chr:start-end)')
    parser.add_argument('--output', required=True, help='Output prefix')
    
    args = parser.parse_args()
    
    # Parse region
    chrom, coords = args.region.split(':')
    start, end = map(int, coords.split('-'))
    
    compare_pileup_modes_fixed(args.model, args.count, chrom, start, end, args.output)

if __name__ == '__main__':
    main()