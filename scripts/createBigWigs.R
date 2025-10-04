


create_dmr_granges <- function(dmr_df, score_column, name_prefix = "") {
  if(nrow(dmr_df) == 0) {
    return(GRanges())
  }
  
  gr <- GRanges(
    seqnames = dmr_df$chr,
    ranges = IRanges(start = dmr_df$start, end = dmr_df$end),
    score = dmr_df[[score_column]],
    name = paste0(name_prefix, "_", 1:nrow(dmr_df)),
    length = dmr_df$end - dmr_df$start + 1,
    nCG = dmr_df$nCG
  )
  
  # Get chromosomes present in the data
  data_chrs <- unique(as.character(seqnames(gr)))
  
  # Get chromosomes that exist in reference genome
  ref_chrs <- seqnames(seqinfo_genome)
  available_chrs <- data_chrs[data_chrs %in% ref_chrs]
  
  # Start with available chromosomes from reference
  if(length(available_chrs) > 0) {
    available_seqinfo <- seqinfo_genome[available_chrs]
  } else {
    available_seqinfo <- Seqinfo()
  }
  
  # Handle chromosomes not in reference genome
  missing_chrs <- data_chrs[!data_chrs %in% ref_chrs]
  if(length(missing_chrs) > 0) {
    cat("Warning: DMR chromosomes not found in reference genome:", paste(missing_chrs, collapse = ", "), "\n")
    
    # Create seqinfo for missing chromosomes
    max_length <- max(seqlengths(seqinfo_genome), na.rm = TRUE)
    missing_seqinfo <- Seqinfo(seqnames = missing_chrs,
                              seqlengths = rep(max_length, length(missing_chrs)),
                              genome = genome_build)
    
    # Combine seqinfos
    if(length(available_seqinfo) > 0) {
      available_seqinfo <- merge(available_seqinfo, missing_seqinfo)
    } else {
      available_seqinfo <- missing_seqinfo
    }
  }
  
  # Assign seqinfo to GRanges
  seqinfo(gr) <- available_seqinfo
  
  return(gr)
}


create_granges <- function(df, score_column, name_prefix = "") {
  if(nrow(df) == 0) {
    return(GRanges())
  }
  
  gr <- GRanges(
    seqnames = df$chr,
    ranges = IRanges(start = df$pos, end = df$pos),
    score = df[[score_column]],
    name = paste0(name_prefix, "_", 1:nrow(df))
  )
  
  # Get chromosomes present in the data
  data_chrs <- unique(as.character(seqnames(gr)))
  
  # Get chromosomes that exist in reference genome
  ref_chrs <- seqnames(seqinfo_genome)
  available_chrs <- data_chrs[data_chrs %in% ref_chrs]
  
  # Start with available chromosomes from reference
  if(length(available_chrs) > 0) {
    available_seqinfo <- seqinfo_genome[available_chrs]
  } else {
    available_seqinfo <- Seqinfo()
  }
  
  # Handle chromosomes not in reference genome
  missing_chrs <- data_chrs[!data_chrs %in% ref_chrs]
  if(length(missing_chrs) > 0) {
    cat("Warning: Chromosomes not found in reference genome:", paste(missing_chrs, collapse = ", "), "\n")
    cat("These will be assigned maximum chromosome length.\n")
    
    # Create seqinfo for missing chromosomes
    max_length <- max(seqlengths(seqinfo_genome), na.rm = TRUE)
    missing_seqinfo <- Seqinfo(seqnames = missing_chrs,
                              seqlengths = rep(max_length, length(missing_chrs)),
                              genome = genome_build)
    
    # Combine seqinfos
    if(length(available_seqinfo) > 0) {
      available_seqinfo <- merge(available_seqinfo, missing_seqinfo)
    } else {
      available_seqinfo <- missing_seqinfo
    }
  }
  
  # Assign seqinfo to GRanges
  seqinfo(gr) <- available_seqinfo
  
  return(gr)
}



genome_build <- "hg38"  # Change this if using a different genome build
cat("Using genome build:", genome_build, "\n")

seqinfo_genome <- seqinfo(BSgenome.Hsapiens.UCSC.hg38)

create_bigwig <- function(methyl_results, 
                          output_dir = "./bigwig",
                          verbose = FALSE) {
  # Create BigWig files from methylation analysis results
  #
  # Args:
  #   methyl_results: The list object returned by run_methylation_analysis function
  #   output_dir: Directory to save BigWig files (default: "./bigwig")
  #   verbose: Whether to print progress messages (default: TRUE)
  
  # Extract objects from the methylation analysis results
  dmlTest_smoothed <- methyl_results$dmlTest_smoothed
  dmlTest_unsmoothed <- methyl_results$dmlTest_unsmoothed
  dmls_smoothed <- methyl_results$dmls_smoothed
  dmls_unsmoothed <- methyl_results$dmls_unsmoothed
  dmrs_smoothed <- methyl_results$dmrs_smoothed
  dmrs_unsmoothed <- methyl_results$dmrs_unsmoothed
  BSobj <- methyl_results$BSobj
  p_threshold <- methyl_results$p_threshold
  delta_threshold <- methyl_results$delta_threshold

    # Extract sample names from BSobj
  sample_names <- sampleNames(BSobj)
  sample1_name <- sample_names[1]
  sample2_name <- sample_names[2]

  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

# ===== 1. Methylation Difference Tracks =====
  if(verbose) cat("\n=== Creating Methylation Difference Tracks ===\n")
  
  # Smoothed methylation differences
  smoothed_bw_path <- NULL
  if(nrow(dmlTest_smoothed) > 0) {
    filename <- paste0(sample1_name, "_vs_", sample2_name, "_methylation_difference_smoothed.bw")
    smoothed_diff_gr <- create_granges(dmlTest_smoothed, "diff", "smoothed_diff")
    smoothed_bw_path <- export.bw(smoothed_diff_gr, filename)
    if(verbose) cat("Created:", filename, "\n")
  }
  
  # Unsmoothed methylation differences
  unsmoothed_bw_path <- NULL
  if(nrow(dmlTest_unsmoothed) > 0) {
    filename <- paste0(sample1_name, "_vs_", sample2_name, "_methylation_difference_unsmoothed.bw")
    unsmoothed_diff_gr <- create_granges(dmlTest_unsmoothed, "diff", "unsmoothed_diff")
    unsmoothed_bw_path <- export.bw(unsmoothed_diff_gr, filename)
    if(verbose) cat("Created:", filename, "\n")
  }
  

  # ===== 2. Statistical Significance Tracks (-log10 p-values) =====
  if(verbose) cat("\n=== Creating Statistical Significance Tracks ===\n")

  # Smoothed -log10 p-values
  smoothed_pval_bw_path <- NULL
  if(nrow(dmlTest_smoothed) > 0) {
    filename <- paste0(sample1_name, "_vs_", sample2_name, "_significance_smoothed.bw")
    
    # Handle zero and very small p-values
    smoothed_pval_clean <- dmlTest_smoothed
    smoothed_pval_clean$log10_pval <- ifelse(
      smoothed_pval_clean$pval > 0 & is.finite(smoothed_pval_clean$pval),
      -log10(smoothed_pval_clean$pval),
      -log10(1e-300)  # Cap very small p-values
    )
    
    smoothed_pval_gr <- create_granges(smoothed_pval_clean, "log10_pval", "smoothed_pval")
    smoothed_pval_bw_path <- export.bw(smoothed_pval_gr, filename)
    if(verbose) cat("Created:", filename, "\n")
  }

  # Unsmoothed -log10 p-values
  unsmoothed_pval_bw_path <- NULL
  if(nrow(dmlTest_unsmoothed) > 0) {
    filename <- paste0(sample1_name, "_vs_", sample2_name, "_significance_unsmoothed.bw")
    
    unsmoothed_pval_clean <- dmlTest_unsmoothed
    unsmoothed_pval_clean$log10_pval <- ifelse(
      unsmoothed_pval_clean$pval > 0 & is.finite(unsmoothed_pval_clean$pval),
      -log10(unsmoothed_pval_clean$pval),
      -log10(1e-300)
    )
    
    unsmoothed_pval_gr <- create_granges(unsmoothed_pval_clean, "log10_pval", "unsmoothed_pval")
    unsmoothed_pval_bw_path <- export.bw(unsmoothed_pval_gr, filename)
    if(verbose) cat("Created:", filename, "\n")
  }

  # ===== 3. Significant DML Tracks (Binary and Weighted) =====
  if(verbose) cat("\n=== Creating Significant DML Tracks ===\n")

  # Track paths
  smoothed_binary_bw_path <- NULL
  smoothed_weighted_bw_path <- NULL
  unsmoothed_binary_bw_path <- NULL
  unsmoothed_weighted_bw_path <- NULL

  # Smoothed significant DMLs - Binary track (1 for significant, 0 for not)
  if(nrow(dmls_smoothed) > 0) {
    # Binary track
    binary_filename <- paste0(sample1_name, "_vs_", sample2_name, "_significant_DMLs_smoothed_binary.bw")
    dmls_smoothed_binary <- dmls_smoothed
    dmls_smoothed_binary$binary_score <- 1
    
    dmls_smoothed_binary_gr <- create_granges(dmls_smoothed_binary, "binary_score", "smoothed_sig_binary")
    smoothed_binary_bw_path <- export.bw(dmls_smoothed_binary_gr, binary_filename)
    if(verbose) cat("Created:", binary_filename, "\n")
    
    # Weighted track
    weighted_filename <- paste0(sample1_name, "_vs_", sample2_name, "_significant_DMLs_smoothed_weighted.bw")
    dmls_smoothed_weighted_gr <- create_granges(dmls_smoothed, "diff", "smoothed_sig_weighted")
    smoothed_weighted_bw_path <- export.bw(dmls_smoothed_weighted_gr, weighted_filename)
    if(verbose) cat("Created:", weighted_filename, "\n")
  }

  # Unsmoothed significant DMLs
  if(nrow(dmls_unsmoothed) > 0) {
    # Binary track
    binary_filename <- paste0(sample1_name, "_vs_", sample2_name, "_significant_DMLs_unsmoothed_binary.bw")
    dmls_unsmoothed_binary <- dmls_unsmoothed
    dmls_unsmoothed_binary$binary_score <- 1
    
    dmls_unsmoothed_binary_gr <- create_granges(dmls_unsmoothed_binary, "binary_score", "unsmoothed_sig_binary")
    unsmoothed_binary_bw_path <- export.bw(dmls_unsmoothed_binary_gr, binary_filename)
    if(verbose) cat("Created:", binary_filename, "\n")
    
    # Weighted track
    weighted_filename <- paste0(sample1_name, "_vs_", sample2_name, "_significant_DMLs_unsmoothed_weighted.bw")
    dmls_unsmoothed_weighted_gr <- create_granges(dmls_unsmoothed, "diff", "unsmoothed_sig_weighted")
    unsmoothed_weighted_bw_path <- export.bw(dmls_unsmoothed_weighted_gr, weighted_filename)
    if(verbose) cat("Created:", weighted_filename, "\n")
  }

  # Track paths
  smoothed_dmr_diff_path <- NULL
  smoothed_dmr_sig_path <- NULL
  smoothed_dmr_cpg_path <- NULL
  unsmoothed_dmr_diff_path <- NULL
  unsmoothed_dmr_sig_path <- NULL
  unsmoothed_dmr_cpg_path <- NULL

  # Smoothed DMRs - Multiple scoring schemes
  if(nrow(dmrs_smoothed) > 0) {
    # DMR track weighted by methylation difference
    diff_filename <- paste0(sample1_name, "_vs_", sample2_name, "_DMRs_smoothed_methylation_difference.bw")
    dmrs_smoothed_diff_gr <- create_dmr_granges(dmrs_smoothed, "diff.Methy", "smoothed_dmr_diff")
    smoothed_dmr_diff_path <- export.bw(dmrs_smoothed_diff_gr, diff_filename)
    if(verbose) cat("Created:", diff_filename, "\n")
    
    # DMR track weighted by area statistic (DSS significance measure)
    if(verbose) cat("Processing smoothed DMR significance track using areaStat...\n")
    
    if("areaStat" %in% colnames(dmrs_smoothed)) {
      sig_filename <- paste0(sample1_name, "_vs_", sample2_name, "_DMRs_smoothed_significance.bw")
      dmrs_smoothed_sig_gr <- create_dmr_granges(dmrs_smoothed, "areaStat", "smoothed_dmr_sig")
      smoothed_dmr_sig_path <- export.bw(dmrs_smoothed_sig_gr, sig_filename)
      if(verbose) cat("Created:", sig_filename, " (using areaStat)\n")
    } else if(verbose) {
      cat("Warning: No 'areaStat' column found in smoothed DMR data. Available columns:\n")
      print(colnames(dmrs_smoothed))
    }
    
    # DMR track weighted by number of CpGs
    cpg_filename <- paste0(sample1_name, "_vs_", sample2_name, "_DMRs_smoothed_CpG_count.bw")
    dmrs_smoothed_cpg_gr <- create_dmr_granges(dmrs_smoothed, "nCG", "smoothed_dmr_cpg")
    smoothed_dmr_cpg_path <- export.bw(dmrs_smoothed_cpg_gr, cpg_filename)
    if(verbose) cat("Created:", cpg_filename, "\n")
  }

  # Unsmoothed DMRs
  if(nrow(dmrs_unsmoothed) > 0) {
    # DMR track weighted by methylation difference
    diff_filename <- paste0(sample1_name, "_vs_", sample2_name, "_DMRs_unsmoothed_methylation_difference.bw")
    dmrs_unsmoothed_diff_gr <- create_dmr_granges(dmrs_unsmoothed, "diff.Methy", "unsmoothed_dmr_diff")
    unsmoothed_dmr_diff_path <- export.bw(dmrs_unsmoothed_diff_gr, diff_filename)
    if(verbose) cat("Created:", diff_filename, "\n")
    
    # DMR track weighted by area statistic (DSS significance measure)
    if(verbose) cat("Processing unsmoothed DMR significance track using areaStat...\n")
    
    if("areaStat" %in% colnames(dmrs_unsmoothed)) {
      sig_filename <- paste0(sample1_name, "_vs_", sample2_name, "_DMRs_unsmoothed_significance.bw")
      dmrs_unsmoothed_sig_gr <- create_dmr_granges(dmrs_unsmoothed, "areaStat", "unsmoothed_dmr_sig")
      unsmoothed_dmr_sig_path <- export.bw(dmrs_unsmoothed_sig_gr, sig_filename)
      if(verbose) cat("Created:", sig_filename, " (using areaStat)\n")
    } else if(verbose) {
      cat("Warning: No 'areaStat' column found in unsmoothed DMR data. Available columns:\n")
      print(colnames(dmrs_unsmoothed))
    }
    
    # DMR track weighted by number of CpGs
    cpg_filename <- paste0(sample1_name, "_vs_", sample2_name, "_DMRs_unsmoothed_CpG_count.bw")
    dmrs_unsmoothed_cpg_gr <- create_dmr_granges(dmrs_unsmoothed, "nCG", "unsmoothed_dmr_cpg")
    unsmoothed_dmr_cpg_path <- export.bw(dmrs_unsmoothed_cpg_gr, cpg_filename)
    if(verbose) cat("Created:", cpg_filename, "\n")
  }



  # ===== 5. Comparison Tracks =====
  if(verbose) cat("\n=== Creating Comparison Tracks ===\n")

  # Track paths
  consensus_binary_path <- NULL
  consensus_weighted_path <- NULL
  smoothed_specific_path <- NULL
  unsmoothed_specific_path <- NULL

  # Method agreement track (where both methods find significant DMLs)
  if(nrow(dmls_smoothed) > 0 && nrow(dmls_unsmoothed) > 0) {
    
    # Find overlapping positions
    smoothed_positions <- paste(dmls_smoothed$chr, dmls_smoothed$pos, sep = "_")
    unsmoothed_positions <- paste(dmls_unsmoothed$chr, dmls_unsmoothed$pos, sep = "_")
    overlapping_positions <- intersect(smoothed_positions, unsmoothed_positions)
    
    if(length(overlapping_positions) > 0) {
      # Create consensus track
      consensus_filename <- paste0(sample1_name, "_vs_", sample2_name, "_DMLs_consensus_both_methods.bw")
      consensus_dmls <- dmls_smoothed[paste(dmls_smoothed$chr, dmls_smoothed$pos, sep = "_") %in% overlapping_positions, ]
      consensus_dmls$consensus_score <- 1
      
      consensus_gr <- create_granges(consensus_dmls, "consensus_score", "consensus")
      consensus_binary_path <- export.bw(consensus_gr, consensus_filename)
      if(verbose) cat("Created:", consensus_filename, "\n")
      
      # Weighted consensus by average methylation difference
      weighted_filename <- paste0(sample1_name, "_vs_", sample2_name, "_DMLs_consensus_weighted_average.bw")
      consensus_weighted <- consensus_dmls
      unsmoothed_subset <- dmls_unsmoothed[paste(dmls_unsmoothed$chr, dmls_unsmoothed$pos, sep = "_") %in% overlapping_positions, ]
      
      # Match order for averaging
      match_idx <- match(paste(consensus_weighted$chr, consensus_weighted$pos, sep = "_"),
                        paste(unsmoothed_subset$chr, unsmoothed_subset$pos, sep = "_"))
      consensus_weighted$avg_diff <- (consensus_weighted$diff + unsmoothed_subset$diff[match_idx]) / 2
      
      consensus_weighted_gr <- create_granges(consensus_weighted, "avg_diff", "consensus_weighted")
      consensus_weighted_path <- export.bw(consensus_weighted_gr, weighted_filename)
      if(verbose) cat("Created:", weighted_filename, "\n")
    }
    
    # Method-specific tracks (unique to each method)
    # Smoothed-specific DMLs
    smoothed_specific <- dmls_smoothed[!smoothed_positions %in% unsmoothed_positions, ]
    if(nrow(smoothed_specific) > 0) {
      smoothed_specific_filename <- paste0(sample1_name, "_vs_", sample2_name, "_DMLs_smoothed_specific.bw")
      smoothed_specific_gr <- create_granges(smoothed_specific, "diff", "smoothed_specific")
      smoothed_specific_path <- export.bw(smoothed_specific_gr, smoothed_specific_filename)
      if(verbose) cat("Created:", smoothed_specific_filename, "\n")
    }
    
    # Unsmoothed-specific DMLs
    unsmoothed_specific <- dmls_unsmoothed[!unsmoothed_positions %in% smoothed_positions, ]
    if(nrow(unsmoothed_specific) > 0) {
      unsmoothed_specific_filename <- paste0(sample1_name, "_vs_", sample2_name, "_DMLs_unsmoothed_specific.bw")
      unsmoothed_specific_gr <- create_granges(unsmoothed_specific, "diff", "unsmoothed_specific")
      unsmoothed_specific_path <- export.bw(unsmoothed_specific_gr, unsmoothed_specific_filename)
      if(verbose) cat("Created:", unsmoothed_specific_filename, "\n")
    }
  }

  # Final return statement with ALL tracks
  return(list(
    methylation_diff_smoothed = smoothed_bw_path,
    methylation_diff_unsmoothed = unsmoothed_bw_path,
    significance_smoothed = smoothed_pval_bw_path,
    significance_unsmoothed = unsmoothed_pval_bw_path,
    significant_dml_smoothed_binary = smoothed_binary_bw_path,
    significant_dml_smoothed_weighted = smoothed_weighted_bw_path,
    significant_dml_unsmoothed_binary = unsmoothed_binary_bw_path,
    significant_dml_unsmoothed_weighted = unsmoothed_weighted_bw_path,
    dmr_smoothed_diff = smoothed_dmr_diff_path,
    dmr_smoothed_significance = smoothed_dmr_sig_path,
    dmr_smoothed_cpg_count = smoothed_dmr_cpg_path,
    dmr_unsmoothed_diff = unsmoothed_dmr_diff_path,
    dmr_unsmoothed_significance = unsmoothed_dmr_sig_path,
    dmr_unsmoothed_cpg_count = unsmoothed_dmr_cpg_path,
    consensus_binary = consensus_binary_path,
    consensus_weighted = consensus_weighted_path,
    smoothed_specific = smoothed_specific_path,
    unsmoothed_specific = unsmoothed_specific_path
  ))


}
