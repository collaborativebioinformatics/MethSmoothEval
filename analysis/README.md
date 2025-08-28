# CpG Density Analyzer

Analyzes CpG methylation density from BED files using sliding genomic windows.

## How it works

1. **Loads BED file** containing CpG site positions (chr, start, end)
2. **Divides chromosome** into fixed-size windows (default: 100kb)
3. **Counts CpG sites** in each window
4. **Calculates density** as CpGs per kilobase for each window
5. **Generates plots** showing density along the chromosome

## Usage

### Basic analysis
```bash
python cpg_density.py data.bed.gz --chromosome 22
```

### Custom parameters
```bash
# Different window size
python cpg_density.py data.bed.gz --chromosome 22 --window-size 50000

# Save data and plots
python cpg_density.py data.bed.gz --chromosome 22 --output results/ --save-data

# Analyze all chromosomes (no window analysis)
python cpg_density.py data.bed.gz
```

## Command line options

| Option | Description | Default |
|--------|-------------|---------|
| `input_file` | Input BED file (can be gzipped) | Required |
| `--chromosome`, `-c` | Analyze specific chromosome (e.g., 22, chr22) | All chromosomes |
| `--window-size`, `-w` | Window size in bp | 100000 |
| `--output`, `-o` | Output directory | Current directory |
| `--no-plots` | Skip generating plots | Generate plots |
| `--save-data` | Save density data to CSV | Don't save |

## Calculations

### Basic statistics
- **Total CpG sites**: Count of all CpG positions
- **Genomic span**: Distance from first to last CpG
- **Overall density**: `total_cpgs / (genomic_span / 1000)` CpGs per kb

### Window analysis
For each window:
- **Window density** = `cpgs_in_window / (window_size / 1000)`
- **Statistics**: Mean, standard deviation, and maximum density across all windows

### Example calculation
- Window: 100,000 bp (100 kb)
- CpG sites in window: 150
- **Density = 150 ÷ 100 = 1.5 CpGs per kb**

## Output files

### Plots (PNG format)
- `cpg_per_chromosome.png` - CpG distribution across chromosomes (if multiple)
- `{chr}_cpg_density_{window_size}kb.png` - Two-panel plot:
  - **Top**: Density profile along chromosome (with mean line)
  - **Bottom**: Histogram of density values

### Data files (if `--save-data` used)
- `{chr}_density_{window_size}kb.csv` - Raw density data

**CSV columns:**
```
chr,window_start,window_end,window_center,cpg_count,density_per_kb
chr22,16050000,16150000,16100000,45,0.45
chr22,16150000,16250000,16200000,120,1.20
```

## Console output example

```
Loading data.bed.gz...
Filtered to chr22: 45,678 CpG sites

=== Basic Statistics ===
Chromosome: chr22
CpG sites: 45,678
Genomic span: 16,050,000 - 51,244,566
Total length: 35.19 Mb
Overall density: 1.30 CpGs per kb

=== Window Analysis (100kb windows) ===
Analyzing chr22: 16,050,000 - 51,244,566
Created 352 windows
Mean density: 1.30 CpGs/kb
Std density: 1.85 CpGs/kb
Max density: 12.45 CpGs/kb
```

## Interpretation

- **High density regions** (>3 CpGs/kb): Potential CpG islands
- **Low density regions** (<1 CpG/kb): Typical genomic background
- **Peaks in density plot**: Regions of interest for further analysis

## Requirements

```python
pandas
matplotlib
numpy
argparse
gzip
```

## Input file format

Standard BED format (tab-separated):
```
chr1    1000    1001
chr1    1050    1051
chr1    1200    1201
...
```

Supports both compressed (.gz) and uncompressed files.