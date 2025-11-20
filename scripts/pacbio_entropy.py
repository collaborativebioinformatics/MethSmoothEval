#!/usr/bin/env python3
"""
PacBio Methylation Entropy Analysis Script
Adapts ONT modkit entropy methodology for PacBio BAM files using MM/ML tags.

This script calculates methylation entropy from PacBio SMRT sequencing data
by extracting base modification calls from MM/ML tags and applying Shannon 
entropy calculations following the modkit approach.

Usage:
    python pacbio_entropy.py --bam input.bam --ref reference.fasta --output entropy.bed [options]
"""

import argparse
import pysam
import numpy as np
from collections import defaultdict, Counter
import logging
from pathlib import Path
import re
from typing import Dict, List, Tuple, Optional, Set

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class MMMLParser:
    """
    Parser for SAM MM/ML base modification tags.
    
    Follows the SAM specification v1.7+ for base modification encoding.
    """
    
    @staticmethod
    def parse_mm_tag(mm_string: str, seq: str, is_reverse: bool) -> Dict[str, List[Tuple[int, str]]]:
        """
        Parse MM tag to extract modification positions and types.
        
        Args:
            mm_string: MM tag value (e.g., "C+m,5,12,0;")
            seq: Read sequence
            is_reverse: Whether read is reverse complemented (FLAG 0x10)
            
        Returns:
            Dictionary mapping base type to list of (position, mod_code) tuples
            Positions are 0-based in SEQ orientation
        """
        modifications = defaultdict(list)
        
        # Split by semicolon to get each modification type
        for mod_spec in mm_string.rstrip(';').split(';'):
            if not mod_spec:
                continue
                
            # Parse: base[+-]mod[.?]?,positions
            match = re.match(r'([ACGTUN])([-+])([a-z]+|\d+)([.?]?),(.+)', mod_spec)
            if not match:
                # Handle case with no positions: "C+m;"
                match = re.match(r'([ACGTUN])([-+])([a-z]+|\d+)([.?]?)', mod_spec)
                if match:
                    continue  # No modifications present
                logger.warning(f"Could not parse MM specification: {mod_spec}")
                continue
            
            base, strand, mod_codes, skip_flag, positions_str = match.groups()
            positions = [int(x) for x in positions_str.split(',')]
            
            # Handle multi-modification codes (e.g., "mh" means could be m or h)
            if len(mod_codes) > 1 and not mod_codes.isdigit():
                # Multi-modification: each position could be any of these mods
                mod_list = list(mod_codes)
            else:
                mod_list = [mod_codes]
            
            # Find all bases of the specified type in the sequence
            # MM tags always refer to the as-sequenced (forward) orientation
            base_positions = []
            for i, b in enumerate(seq):
                if base == 'N' or b.upper() == base:
                    base_positions.append(i)
            
            # Convert delta encoding to absolute positions
            # The deltas are relative to the PREVIOUS modified position, not absolute index
            current_base_idx = 0  # Index into base_positions array
            
            for skip_count in positions:
                # Skip this many bases of the specified type
                current_base_idx += skip_count
                
                if current_base_idx < len(base_positions):
                    seq_pos = base_positions[current_base_idx]
                    # Store modification info
                    # For multi-mods, we'll track that multiple types are possible
                    for mod_code in mod_list:
                        modifications[base].append((seq_pos, mod_code, strand))
                    # Move to next position for the next delta
                    current_base_idx += 1
                else:
                    # This can happen if MM tag is out of sync with sequence
                    # Log only at debug level to avoid spam
                    logger.debug(f"Position index {current_base_idx} exceeds {len(base_positions)} "
                               f"base positions for {base} in sequence of length {len(seq)}")
                    break
        
        return modifications
    
    @staticmethod
    def parse_ml_tag(ml_array: List[int], mm_mods: Dict[str, List[Tuple[int, str]]], 
                     multi_mod_format: bool = False) -> Dict[Tuple[int, str], float]:
        """
        Parse ML tag to extract modification probabilities.
        
        Args:
            ml_array: ML tag byte array (0-255 values)
            mm_mods: Parsed MM modifications
            multi_mod_format: Whether MM uses multi-modification format
            
        Returns:
            Dictionary mapping (position, mod_code) to probability (0.0-1.0)
        """
        probabilities = {}
        ml_idx = 0
        
        # Flatten modifications in the order they appear
        all_mods = []
        for base_type in sorted(mm_mods.keys()):
            for pos, mod_code, strand in mm_mods[base_type]:
                all_mods.append((pos, mod_code, strand))
        
        # Map probabilities to positions
        for pos, mod_code, strand in all_mods:
            if ml_idx < len(ml_array):
                # Convert 0-255 to 0.0-1.0 probability
                prob = (ml_array[ml_idx] + 0.5) / 256.0
                probabilities[(pos, mod_code)] = prob
                ml_idx += 1
            else:
                logger.warning(f"ML array shorter than expected modifications")
                break
        
        return probabilities


class PacBioEntropyCalculator:
    """
    Calculate methylation entropy from PacBio BAM files using MM/ML tags.
    
    Implements the entropy calculation methodology from ONT modkit
    adapted for PacBio's SAM-compliant modification tags.
    """
    
    def __init__(self, reference_fasta: str, motif: str = 'CG', 
                 min_coverage: int = 10, filter_threshold: float = 0.66,
                 window_size: int = 50, num_positions: int = 4,
                 max_filtered_frac: float = 0.5, combine_strands: bool = False,
                 mod_codes: Optional[List[str]] = None):
        """
        Initialize the entropy calculator.
        
        Args:
            reference_fasta: Path to reference genome FASTA file
            motif: DNA motif to analyze (default: CG for CpG sites)
            min_coverage: Minimum coverage threshold for analysis (default: 10)
                         Note: modkit uses 3, but 10-30 recommended for reliable estimates
            filter_threshold: Minimum probability threshold for confident calls (default: 0.66)
                            Calls above this are 'm', below (1-threshold) are 'u', between are '*'
            window_size: Maximum bp distance for num_positions motifs (default: 50)
                        This ensures analyzed sites are in the same regulatory context
            num_positions: Number of modification sites per entropy calculation (default: 4)
            max_filtered_frac: Maximum fraction of filtered positions allowed per read (default: 0.5)
                              Reads with more wildcards than this are discarded
            combine_strands: Combine modification counts from both strands (default: False)
                           When True with CG motif, behaves like modkit's --cpg flag
            mod_codes: Modification codes to analyze (default: ['m'] for 5mC)
        """
        self.reference = pysam.FastaFile(reference_fasta)
        self.motif = motif.upper()
        self.min_coverage = min_coverage
        self.filter_threshold = filter_threshold
        self.window_size = window_size
        self.num_positions = num_positions
        self.max_filtered_frac = max_filtered_frac
        self.combine_strands = combine_strands
        self.mod_codes = mod_codes or ['m']  # Default to 5mC
        self.parser = MMMLParser()
        
    def find_motif_positions(self, chromosome: str, start: int, end: int) -> List[int]:
        """
        Find all motif positions in the specified genomic region.
        
        When combine_strands is True and motif is palindromic (e.g., CG),
        returns positions from both strands mapped to the positive strand.
        
        Args:
            chromosome: Chromosome name
            start: Start position (0-based)
            end: End position (0-based)
            
        Returns:
            List of motif positions within the region
        """
        try:
            sequence = self.reference.fetch(chromosome, start, end).upper()
            positions = []
            
            # Find all occurrences of the motif on forward strand
            for i in range(len(sequence) - len(self.motif) + 1):
                if sequence[i:i+len(self.motif)] == self.motif:
                    positions.append(start + i)
            
            # If combining strands and motif is self-complementary (like CG)
            if self.combine_strands:
                # Get reverse complement of motif
                complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
                rev_motif = ''.join(complement.get(b, b) for b in reversed(self.motif))
                
                # For CG motif, reverse complement is also CG (palindrome)
                # For other motifs, find reverse complement occurrences
                if rev_motif != self.motif:
                    for i in range(len(sequence) - len(rev_motif) + 1):
                        if sequence[i:i+len(rev_motif)] == rev_motif:
                            # Map to positive strand position
                            positions.append(start + i)
            
            return sorted(set(positions))  # Remove duplicates and sort
            
        except Exception as e:
            logger.warning(f"Could not fetch sequence for {chromosome}:{start}-{end}: {e}")
            return []
    
    def extract_modification_data(self, bamfile: pysam.AlignmentFile, 
                                 chromosome: str, start: int, end: int,
                                 motif_positions: Set[int]) -> Dict[str, Dict[int, Tuple[str, float]]]:
        """
        Extract modification data from BAM reads using MM/ML tags.
        
        Args:
            bamfile: Opened BAM file handle
            chromosome: Chromosome name
            start: Start position
            end: End position
            motif_positions: Set of reference positions to analyze
            
        Returns:
            Dictionary mapping read IDs to position -> (call, probability) mappings
        """
        read_data = defaultdict(dict)
        
        for read in bamfile.fetch(chromosome, start, end):
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            
            # Check for MM/ML tags
            if not read.has_tag('MM'):
                continue
            
            read_id = read.query_name
            is_reverse = read.is_reverse
            
            try:
                # Validate MN tag if present (sequence length check)
                if read.has_tag('MN'):
                    expected_len = read.get_tag('MN')
                    if expected_len != len(read.query_sequence):
                        logger.debug(f"Read {read_id}: MN tag ({expected_len}) doesn't match "
                                   f"current SEQ length ({len(read.query_sequence)}). "
                                   f"MM/ML tags may be stale. Skipping read.")
                        continue
                
                # Parse MM tag
                mm_string = read.get_tag('MM')
                mm_mods = self.parser.parse_mm_tag(mm_string, read.query_sequence, is_reverse)
                
                # Parse ML tag if present
                probabilities = {}
                if read.has_tag('ML'):
                    ml_array = read.get_tag('ML')
                    probabilities = self.parser.parse_ml_tag(ml_array, mm_mods)
                
                # Get reference positions from read
                query_to_ref = {}
                ref_positions = read.get_reference_positions(full_length=True)
                for query_pos, ref_pos in enumerate(ref_positions):
                    if ref_pos is not None:
                        query_to_ref[query_pos] = ref_pos
                
                # Extract modifications at motif positions
                for base_type in mm_mods:
                    # Only process relevant base types (e.g., 'C' for CpG)
                    if self.motif[0] != base_type:
                        continue
                    
                    for query_pos, mod_code, strand in mm_mods[base_type]:
                        # Filter by modification code
                        if mod_code not in self.mod_codes:
                            continue
                        
                        # If combining strands, accept both + and - strand calls
                        # Otherwise, only accept calls from the appropriate strand
                        if not self.combine_strands and strand == '-':
                            continue
                        
                        # Convert query position to reference position
                        if query_pos in query_to_ref:
                            ref_pos = query_to_ref[query_pos]
                            
                            # Only include positions in our motif set
                            if ref_pos in motif_positions:
                                # Get probability
                                prob = probabilities.get((query_pos, mod_code), 0.5)
                                
                                # Convert to call: m (methylated), u (unmethylated), * (uncertain)
                                if prob >= self.filter_threshold:
                                    call = 'm'
                                elif prob <= (1.0 - self.filter_threshold):
                                    call = 'u'
                                else:
                                    call = '*'
                                
                                read_data[read_id][ref_pos] = (call, prob)
                
            except Exception as e:
                logger.debug(f"Could not process read {read_id}: {e}")
                continue
        
        return read_data
    
    def calculate_entropy_window(self, read_data: Dict[str, Dict[int, Tuple[str, float]]], 
                               motif_positions: List[int]) -> Optional[Tuple[float, int]]:
        """
        Calculate methylation entropy for a set of motif positions.
        
        Implements the Shannon entropy formula from modkit:
        ME = (-1/N) × Σ Pr(ni) × log₂(Pr(ni))
        
        Args:
            read_data: Modification data for reads in the region
            motif_positions: List of motif positions to analyze
            
        Returns:
            Tuple of (entropy value, coverage) or None if insufficient data
        """
        if len(motif_positions) < self.num_positions:
            return None
        
        # Use first num_positions motifs in the window
        positions = sorted(motif_positions)[:self.num_positions]
        
        # Calculate max_filtered threshold
        max_filtered_positions = int(self.max_filtered_frac * self.num_positions)
        
        # Collect methylation patterns across reads
        patterns = []
        
        for read_id, read_mods in read_data.items():
            pattern = []
            valid_positions = 0
            filtered_positions = 0
            
            for pos in positions:
                if pos in read_mods:
                    call, prob = read_mods[pos]
                    pattern.append(call)
                    if call == '*':
                        filtered_positions += 1
                    else:
                        valid_positions += 1
                else:
                    pattern.append('*')  # No coverage
                    filtered_positions += 1
            
            # Apply max_filtered_positions filter (like modkit)
            if filtered_positions > max_filtered_positions:
                continue  # Discard read with too many uncertain calls
            
            # Only include patterns with sufficient valid positions
            min_valid = max(1, int(len(positions) * 0.5))
            if valid_positions >= min_valid:
                patterns.append(''.join(pattern))
        
        if len(patterns) < self.min_coverage:
            return None
        
        # Calculate entropy using prefix trie approach
        entropy = self._calculate_shannon_entropy(patterns)
        coverage = len(patterns)
        
        return (entropy, coverage)
    
    def _calculate_shannon_entropy(self, patterns: List[str]) -> float:
        """
        Calculate Shannon entropy from methylation patterns.
        
        Implements the modkit entropy algorithm with wildcard handling.
        
        Args:
            patterns: List of methylation patterns (strings with 'm', 'u', '*')
            
        Returns:
            Shannon entropy value
        """
        # Expand patterns with wildcards using prefix trie logic
        expanded_patterns = []
        
        for pattern in patterns:
            if '*' in pattern:
                # For wildcards, contribute fractionally to possible patterns
                wildcard_count = pattern.count('*')
                
                # Generate all possible expansions (limit for computational efficiency)
                if wildcard_count <= 4:
                    expansions = self._expand_wildcards(pattern)
                    weight = 1.0 / len(expansions)
                    
                    # Add weighted contributions
                    for expansion in expansions:
                        # Add multiple times to maintain integer counts for frequency
                        expanded_patterns.extend([expansion] * max(1, int(weight * 100)))
                else:
                    # Too many wildcards - use most likely pattern
                    # Replace wildcards with unmethylated (conservative)
                    solid_pattern = pattern.replace('*', 'u')
                    expanded_patterns.append(solid_pattern)
            else:
                expanded_patterns.append(pattern)
        
        if not expanded_patterns:
            return 0.0
        
        # Count pattern frequencies
        pattern_counts = Counter(expanded_patterns)
        total_count = sum(pattern_counts.values())
        
        if total_count == 0:
            return 0.0
        
        # Calculate Shannon entropy: H = -Σ p(i) * log2(p(i))
        entropy = 0.0
        for count in pattern_counts.values():
            probability = count / total_count
            if probability > 0:
                entropy -= probability * np.log2(probability)
        
        return entropy
    
    def _expand_wildcards(self, pattern: str) -> List[str]:
        """
        Expand a pattern with wildcards (*) to all possible concrete patterns.
        
        Args:
            pattern: Pattern string with wildcards
            
        Returns:
            List of all possible concrete patterns
        """
        if '*' not in pattern:
            return [pattern]
        
        # Find first wildcard and recursively expand
        idx = pattern.index('*')
        prefix = pattern[:idx]
        suffix = pattern[idx+1:]
        
        results = []
        for replacement in ['m', 'u']:
            new_pattern = prefix + replacement + suffix
            results.extend(self._expand_wildcards(new_pattern))
        
        return results
    
    def process_region(self, bamfile: pysam.AlignmentFile, chromosome: str, 
                      start: int, end: int) -> List[Tuple[str, int, int, float, int]]:
        """
        Process a genomic region and calculate entropy values.
        
        Uses a sliding window approach to find sets of num_positions motifs
        within window_size base pairs of each other.
        
        Args:
            bamfile: Opened BAM file handle
            chromosome: Chromosome name
            start: Region start position
            end: Region end position
            
        Returns:
            List of tuples: (chromosome, start, end, entropy, coverage)
        """
        results = []
        
        # Find motif positions in the region
        motif_positions = self.find_motif_positions(chromosome, start, end)
        
        if len(motif_positions) < self.num_positions:
            return results
        
        # Convert to set for faster lookup
        motif_set = set(motif_positions)
        
        # Extract modification data for the region
        read_data = self.extract_modification_data(bamfile, chromosome, start, end, motif_set)
        
        # Calculate entropy for windows where num_positions motifs fit within window_size bp
        i = 0
        while i <= len(motif_positions) - self.num_positions:
            # Get potential window of num_positions consecutive motifs
            window_positions = motif_positions[i:i + self.num_positions]
            
            # Check if all positions fit within window_size
            window_span = window_positions[-1] - window_positions[0]
            
            if window_span <= self.window_size:
                # Valid window - calculate entropy
                window_start = window_positions[0]
                window_end = window_positions[-1] + len(self.motif)
                
                result = self.calculate_entropy_window(read_data, window_positions)
                
                if result is not None:
                    entropy, coverage = result
                    results.append((chromosome, window_start, window_end, entropy, coverage))
                
                # Move to next position (creates overlapping windows)
                i += 1
            else:
                # Window too large - skip first position and try next
                i += 1
        
        return results
    
    def run_analysis(self, bam_path: str, output_path: str, 
                    regions: Optional[List[Tuple[str, int, int]]] = None,
                    output_format: str = 'bed'):
        """
        Run complete entropy analysis on BAM file.
        
        Args:
            bam_path: Path to input BAM file
            output_path: Path for output file
            regions: Optional list of specific regions to analyze
            output_format: Output format ('bed' or 'bedgraph')
        """
        logger.info(f"Starting entropy analysis of {bam_path}")
        logger.info(f"Parameters: motif={self.motif}, min_coverage={self.min_coverage}, "
                   f"filter_threshold={self.filter_threshold}, num_positions={self.num_positions}")
        
        with pysam.AlignmentFile(bam_path, 'rb') as bamfile:
            all_results = []
            windows_processed = 0
            
            if regions:
                # Process specified regions
                for chrom, start, end in regions:
                    logger.info(f"Processing region {chrom}:{start}-{end}")
                    results = self.process_region(bamfile, chrom, start, end)
                    all_results.extend(results)
                    windows_processed += 1
            else:
                # Process entire genome in windows
                for ref_info in bamfile.header.as_dict()['SQ']:
                    chrom = ref_info['SN']
                    length = ref_info['LN']
                    
                    logger.info(f"Processing chromosome {chrom} (length: {length:,} bp)")
                    
                    chrom_windows = 0
                    for start in range(0, length, self.window_size):
                        end = min(start + self.window_size * 2, length)  # Overlap windows
                        results = self.process_region(bamfile, chrom, start, end)
                        all_results.extend(results)
                        chrom_windows += 1
                        windows_processed += 1
                        
                        if chrom_windows % 50 == 0:
                            logger.info(f"  Processed {chrom_windows} windows on {chrom}, "
                                      f"{len(all_results)} total entropy windows calculated")
                    
                    logger.info(f"Completed {chrom}: {chrom_windows} windows processed")
        
        # Write results
        if output_format == 'bedgraph':
            self._write_bedgraph_output(all_results, output_path)
        else:
            self._write_bed_output(all_results, output_path)
            
        logger.info(f"="*60)
        logger.info(f"Analysis complete!")
        logger.info(f"Entropy windows calculated: {len(all_results)}")
        if len(all_results) > 0:
            entropies = [e for _, _, _, e, _ in all_results]
            coverages = [c for _, _, _, _, c in all_results]
            logger.info(f"Entropy range: {min(entropies):.3f} - {max(entropies):.3f}")
            logger.info(f"Mean entropy: {np.mean(entropies):.3f} (median: {np.median(entropies):.3f})")
            logger.info(f"Coverage range: {min(coverages)} - {max(coverages)} reads")
            logger.info(f"Mean coverage: {np.mean(coverages):.1f} reads")
        else:
            logger.warning(f"No entropy windows calculated! Check:")
            logger.warning(f"  - Do you have {self.num_positions} {self.motif} motifs within {self.window_size} bp?")
            logger.warning(f"  - Is coverage >= {self.min_coverage}?")
            logger.warning(f"  - Are MM/ML tags present in the BAM?")
        logger.info(f"Results written to {output_path}")
        logger.info(f"="*60)
    
    def _write_bed_output(self, results: List[Tuple[str, int, int, float, int]], 
                         output_path: str):
        """
        Write results to BED format file.
        
        Args:
            results: List of entropy results
            output_path: Output file path
        """
        with open(output_path, 'w') as f:
            # Write header
            f.write(f"# PacBio Methylation Entropy Analysis\n")
            f.write(f"# Motif: {self.motif}, Modification codes: {','.join(self.mod_codes)}\n")
            f.write(f"# Parameters: num_positions={self.num_positions}, window_size={self.window_size}, "
                   f"min_coverage={self.min_coverage}, filter_threshold={self.filter_threshold:.2f}\n")
            f.write(f"# Max filtered fraction: {self.max_filtered_frac:.1%}, "
                   f"combine_strands={self.combine_strands}\n")
            f.write(f"# Columns: chromosome, start, end, entropy, coverage\n")
            
            # Write data
            for chrom, start, end, entropy, coverage in sorted(results):
                f.write(f"{chrom}\t{start}\t{end}\t{entropy:.6f}\t{coverage}\n")
    
    def _write_bedgraph_output(self, results: List[Tuple[str, int, int, float, int]], 
                              output_path: str):
        """
        Write results to bedGraph format file.
        
        Args:
            results: List of entropy results
            output_path: Output file path
        """
        with open(output_path, 'w') as f:
            # Write track header
            f.write(f"track type=bedGraph name=\"PacBio_Methylation_Entropy\" "
                   f"description=\"Methylation entropy from PacBio MM/ML tags\"\n")
            
            # Write data (bedGraph uses 4 columns: chr, start, end, value)
            for chrom, start, end, entropy, coverage in sorted(results):
                f.write(f"{chrom}\t{start}\t{end}\t{entropy:.6f}\n")


def main():
    """Main function to run the PacBio entropy analysis."""
    parser = argparse.ArgumentParser(
        description="Calculate methylation entropy from PacBio BAM files with MM/ML tags",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    # Basic 5mC entropy analysis on CpG sites
    python pacbio_entropy.py --bam sample.bam --ref genome.fa --output entropy.bed
    
    # Analyze 5hmC modifications
    python pacbio_entropy.py --bam sample.bam --ref genome.fa --output entropy.bed \\
        --mod-codes h
    
    # Multiple modification types
    python pacbio_entropy.py --bam sample.bam --ref genome.fa --output entropy.bed \\
        --mod-codes m h
    
    # Custom motif and parameters
    python pacbio_entropy.py --bam sample.bam --ref genome.fa --output entropy.bed \\
        --motif CCWGG --min-coverage 15 --window-size 5000
    
    # Analyze specific regions
    python pacbio_entropy.py --bam sample.bam --ref genome.fa --output entropy.bed \\
        --regions chr1:1000000-2000000 chr2:500000-1500000
    
    # Output as bedGraph for visualization
    python pacbio_entropy.py --bam sample.bam --ref genome.fa --output entropy.bedgraph \\
        --format bedgraph
        """
    )
    
    parser.add_argument('--bam', required=True, 
                       help='Input PacBio BAM file with MM/ML tags')
    parser.add_argument('--ref', required=True,
                       help='Reference genome FASTA file')
    parser.add_argument('--output', required=True,
                       help='Output file for entropy values')
    parser.add_argument('--motif', default='CG',
                       help='DNA motif to analyze (default: CG for CpG sites)')
    parser.add_argument('--mod-codes', nargs='+', default=['m'],
                       help='Modification codes to analyze (default: m for 5mC). '
                            'Options: m (5mC), h (5hmC), f (5fC), c (5caC), a (6mA), etc.')
    parser.add_argument('--min-coverage', type=int, default=5,
                       help='Minimum coverage threshold (default: 5). '
                            'Modkit uses 3, but 10-30 recommended for reliable estimates')
    parser.add_argument('--filter-threshold', type=float, default=0.66,
                       help='Minimum probability threshold for confident calls (default: 0.66). '
                            'Modkit uses ~0.1, but 0.66 recommended for PacBio kinetic data')
    parser.add_argument('--window-size', type=int, default=50,
                       help='Maximum bp distance for num_positions motifs (default: 50, matches modkit). '
                            'This ensures analyzed sites are in the same regulatory context')
    parser.add_argument('--num-positions', type=int, default=4,
                       help='Number of motif positions per entropy calculation (default: 4, matches modkit)')
    parser.add_argument('--max-filtered-frac', type=float, default=0.5,
                       help='Maximum fraction of uncertain calls allowed per read (default: 0.5, matches modkit). '
                            'Reads with more wildcards are discarded')
    parser.add_argument('--combine-strands', action='store_true',
                       help='Combine modification counts from both strands (like modkit --cpg). '
                            'Use with --motif CG for standard CpG analysis')
    parser.add_argument('--regions', nargs='*',
                       help='Specific regions to analyze (format: chr:start-end)')
    parser.add_argument('--format', choices=['bed', 'bedgraph'], default='bed',
                       help='Output format (default: bed)')
    parser.add_argument('--verbose', action='store_true',
                       help='Enable verbose logging')
    
    args = parser.parse_args()
    
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)
    
    # Parse regions if provided
    regions = None
    if args.regions:
        regions = []
        for region_str in args.regions:
            match = re.match(r'(.+):(\d+)-(\d+)', region_str)
            if match:
                chrom, start, end = match.groups()
                regions.append((chrom, int(start), int(end)))
            else:
                logger.error(f"Invalid region format: {region_str}")
                return 1
    
    # Initialize calculator
    calculator = PacBioEntropyCalculator(
        reference_fasta=args.ref,
        motif=args.motif,
        min_coverage=args.min_coverage,
        filter_threshold=args.filter_threshold,
        window_size=args.window_size,
        num_positions=args.num_positions,
        max_filtered_frac=args.max_filtered_frac,
        combine_strands=args.combine_strands,
        mod_codes=args.mod_codes
    )
    
    # Run analysis
    try:
        calculator.run_analysis(args.bam, args.output, regions, args.format)
        return 0
    except Exception as e:
        logger.error(f"Analysis failed: {e}", exc_info=args.verbose)
        return 1


if __name__ == '__main__':
    exit(main())