#!/usr/bin/env python3
"""
CpG Density Analysis Tool
Analyzes CpG methylation density from bed files with configurable parameters
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import argparse
import gzip
import sys
import os

def load_bed_file(filepath, chromosome=None):
    """Load bed file, optionally filtering by chromosome"""
    print(f"Loading {filepath}...")
    
    # Handle compressed files
    if filepath.endswith('.gz'):
        df = pd.read_csv(filepath, sep='\t', header=None, compression='gzip')
    else:
        df = pd.read_csv(filepath, sep='\t', header=None)
    
    # Set standard column names
    df.columns = ['chr', 'start', 'end'] + [f'col{i}' for i in range(3, len(df.columns))]
    
    # Filter by chromosome if specified
    if chromosome:
        if not chromosome.startswith('chr'):
            chromosome = f'chr{chromosome}'
        df = df[df['chr'] == chromosome]
        print(f"Filtered to {chromosome}: {len(df):,} CpG sites")
    else:
        print(f"Total CpG sites: {len(df):,}")
    
    if len(df) == 0:
        print(f"Error: No data found for chromosome {chromosome}")
        sys.exit(1)
        
    return df

def calculate_basic_stats(df):
    """Calculate and display basic statistics"""
    print("\n=== Basic Statistics ===")
    
    if len(df['chr'].unique()) == 1:
        chrom = df['chr'].iloc[0]
        span_start = df['start'].min()
        span_end = df['end'].max()
        span_length = span_end - span_start
        
        print(f"Chromosome: {chrom}")
        print(f"CpG sites: {len(df):,}")
        print(f"Genomic span: {span_start:,} - {span_end:,}")
        print(f"Total length: {span_length/1e6:.2f} Mb")
        print(f"Overall density: {len(df)/(span_length/1000):.2f} CpGs per kb")
    else:
        cpg_per_chr = df.groupby('chr').size()
        print("CpGs per chromosome:")
        for chrom, count in cpg_per_chr.items():
            print(f"  {chrom}: {count:,}")

def plot_chromosome_distribution(df, output_dir):
    """Plot CpG distribution across chromosomes"""
    if len(df['chr'].unique()) <= 1:
        return
        
    cpg_per_chr = df.groupby('chr').size().sort_index()
    
    plt.figure(figsize=(12, 6))
    bars = plt.bar(range(len(cpg_per_chr)), cpg_per_chr.values)
    plt.title('CpG Sites per Chromosome', fontsize=14)
    plt.xlabel('Chromosome')
    plt.ylabel('Number of CpG Sites')
    plt.xticks(range(len(cpg_per_chr)), cpg_per_chr.index, rotation=45)
    
    # Add value labels on bars
    for bar, value in zip(bars, cpg_per_chr.values):
        plt.text(bar.get_x() + bar.get_width()/2, bar.get_height() + max(cpg_per_chr.values)*0.01,
                f'{value:,}', ha='center', va='bottom', fontsize=8)
    
    plt.tight_layout()
    output_path = os.path.join(output_dir, 'cpg_per_chromosome.png')
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"Saved: {output_path}")
    plt.show()

def calculate_density_windows(df, window_size, chromosome=None):
    """Calculate CpG density in genomic windows"""
    if len(df['chr'].unique()) > 1 and chromosome is None:
        print("Error: Multiple chromosomes found. Please specify --chromosome for window analysis.")
        return None
        
    chrom_data = df.copy()
    chrom = chrom_data['chr'].iloc[0]
    
    print(f"\n=== Window Analysis ({window_size/1000:.0f}kb windows) ===")
    
    min_pos = chrom_data['start'].min()
    max_pos = chrom_data['end'].max()
    
    print(f"Analyzing {chrom}: {min_pos:,} - {max_pos:,}")
    
    densities = []
    window_count = 0
    
    for window_start in range(int(min_pos), int(max_pos), window_size):
        window_end = window_start + window_size
        cpgs_in_window = chrom_data[
            (chrom_data['start'] >= window_start) & 
            (chrom_data['start'] < window_end)
        ]
        
        density = len(cpgs_in_window) / (window_size / 1000)  # CpGs per kb
        densities.append({
            'chr': chrom,
            'window_start': window_start,
            'window_end': window_end,
            'window_center': window_start + window_size // 2,
            'cpg_count': len(cpgs_in_window),
            'density_per_kb': density
        })
        window_count += 1
    
    density_df = pd.DataFrame(densities)
    
    print(f"Created {window_count:,} windows")
    print(f"Mean density: {density_df['density_per_kb'].mean():.2f} CpGs/kb")
    print(f"Std density: {density_df['density_per_kb'].std():.2f} CpGs/kb")
    print(f"Max density: {density_df['density_per_kb'].max():.2f} CpGs/kb")
    
    return density_df

def plot_density_profile(density_df, window_size, output_dir, chromosome=None):
    """Plot CpG density along chromosome"""
    if density_df is None or len(density_df) == 0:
        return
        
    chrom = density_df['chr'].iloc[0]
    
    plt.figure(figsize=(15, 8))
    
    # Main density plot
    plt.subplot(2, 1, 1)
    plt.plot(density_df['window_center']/1e6, density_df['density_per_kb'], 
             'b-', alpha=0.7, linewidth=1)
    plt.title(f'CpG Density along {chrom} ({window_size/1000:.0f}kb windows)', fontsize=14)
    plt.ylabel('CpG Density (per kb)')
    plt.grid(True, alpha=0.3)
    
    # Add mean line
    mean_density = density_df['density_per_kb'].mean()
    plt.axhline(y=mean_density, color='r', linestyle='--', alpha=0.7, 
                label=f'Mean: {mean_density:.2f} CpGs/kb')
    plt.legend()
    
    # Density histogram
    plt.subplot(2, 1, 2)
    plt.hist(density_df['density_per_kb'], bins=50, alpha=0.7, color='skyblue', edgecolor='black')
    plt.axvline(x=mean_density, color='r', linestyle='--', alpha=0.7)
    plt.xlabel('CpG Density (per kb)')
    plt.ylabel('Number of Windows')
    plt.title('Distribution of CpG Density Values')
    plt.grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    # Save plot
    chrom_clean = chrom.replace('chr', '')
    output_path = os.path.join(output_dir, f'{chrom_clean}_cpg_density_{window_size//1000}kb.png')
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"Saved: {output_path}")
    plt.show()

def save_density_data(density_df, output_dir, chromosome, window_size):
    """Save density data to CSV"""
    if density_df is None or len(density_df) == 0:
        return
        
    chrom_clean = chromosome.replace('chr', '') if chromosome else 'all'
    output_path = os.path.join(output_dir, f'{chrom_clean}_density_{window_size//1000}kb.csv')
    density_df.to_csv(output_path, index=False)
    print(f"Saved density data: {output_path}")

def main():
    parser = argparse.ArgumentParser(
        description='Analyze CpG density from methylation bed files',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic analysis of chr22
  python cpg_density_analyzer.py data.bed.gz --chromosome 22
  
  # Full genome with custom window size
  python cpg_density_analyzer.py data.bed.gz --window-size 50000
  
  # Quick analysis with plots
  python cpg_density_analyzer.py data.bed.gz --chromosome 22 --window-size 10000 --output results/
        """
    )
    
    parser.add_argument('input_file', 
                       help='Input bed file (can be gzipped)')
    
    parser.add_argument('--chromosome', '-c', 
                       help='Analyze specific chromosome (e.g., 22, chr22)')
    
    parser.add_argument('--window-size', '-w', type=int, default=100000,
                       help='Window size for density calculation in bp (default: 100000)')
    
    parser.add_argument('--output', '-o', default='.',
                       help='Output directory for plots and data (default: current directory)')
    
    parser.add_argument('--no-plots', action='store_true',
                       help='Skip generating plots')
    
    parser.add_argument('--save-data', action='store_true',
                       help='Save density data to CSV file')
    
    args = parser.parse_args()
    
    # Validate inputs
    if not os.path.exists(args.input_file):
        print(f"Error: Input file {args.input_file} not found")
        sys.exit(1)
    
    if not os.path.exists(args.output):
        os.makedirs(args.output)
        print(f"Created output directory: {args.output}")
    
    print(f"CpG Density Analyzer")
    print(f"Input file: {args.input_file}")
    print(f"Window size: {args.window_size:,} bp ({args.window_size/1000:.0f}kb)")
    print(f"Output directory: {args.output}")
    
    # Load data
    df = load_bed_file(args.input_file, args.chromosome)
    
    # Basic statistics
    calculate_basic_stats(df)
    
    # Chromosome distribution plot (if multiple chromosomes)
    if not args.no_plots and len(df['chr'].unique()) > 1:
        plot_chromosome_distribution(df, args.output)
    
    # Window-based density analysis
    if len(df['chr'].unique()) == 1 or args.chromosome:
        density_df = calculate_density_windows(df, args.window_size, args.chromosome)
        
        if not args.no_plots:
            plot_density_profile(density_df, args.window_size, args.output, args.chromosome)
        
        if args.save_data:
            save_density_data(density_df, args.output, args.chromosome, args.window_size)
    else:
        print("\nSkipping window analysis - multiple chromosomes detected.")
        print("Use --chromosome option to analyze a specific chromosome.")

if __name__ == '__main__':
    main()
