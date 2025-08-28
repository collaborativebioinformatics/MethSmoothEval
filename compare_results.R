# compare_results.R
# Compare DMR results between smoothed and raw data

# Load both result files
dmr_raw <- read.table("dmr_results_raw.txt", header=TRUE, sep="\t")
dmr_smoothed <- read.table("dmr_results_smoothed.txt", header=TRUE, sep="\t")

dml_raw <- read.table("significant_dmls_raw.txt", header=TRUE, sep="\t")
dml_smoothed <- read.table("significant_dmls_smoothed.txt", header=TRUE, sep="\t")

cat("=============================================================================\n")
cat("COMPARISON: RAW vs SMOOTHED DATA ANALYSIS\n")
cat("=============================================================================\n\n")

# Compare DMR counts
cat("DMR COMPARISON:\n")
cat("---------------\n")
cat(paste("Raw data DMRs:", nrow(dmr_raw), "\n"))
cat(paste("Smoothed data DMRs:", nrow(dmr_smoothed), "\n"))
cat(paste("Difference:", nrow(dmr_smoothed) - nrow(dmr_raw), "\n\n"))

# Compare DML counts
cat("DML COMPARISON:\n")
cat("---------------\n")
cat(paste("Raw data significant DMLs:", nrow(dml_raw), "\n"))
cat(paste("Smoothed data significant DMLs:", nrow(dml_smoothed), "\n"))
cat(paste("Difference:", nrow(dml_smoothed) - nrow(dml_raw), "\n\n"))

if(nrow(dmr_raw) > 0 && nrow(dmr_smoothed) > 0) {
  
  # Compare DMR characteristics
  cat("DMR CHARACTERISTICS COMPARISON:\n")
  cat("-------------------------------\n")
  
  comparison_df <- data.frame(
    Metric = c("Mean Length (bp)", "Median Length (bp)", 
               "Mean CpGs per DMR", "Median CpGs per DMR",
               "Mean |Methylation Difference|", "Median |Methylation Difference|"),
    Raw = c(round(mean(dmr_raw$length), 1), round(median(dmr_raw$length), 1),
            round(mean(dmr_raw$nCG), 1), round(median(dmr_raw$nCG), 1),
            round(mean(abs(dmr_raw$meanMethy1 - dmr_raw$meanMethy2)), 3),
            round(median(abs(dmr_raw$meanMethy1 - dmr_raw$meanMethy2)), 3)),
    Smoothed = c(round(mean(dmr_smoothed$length), 1), round(median(dmr_smoothed$length), 1),
                 round(mean(dmr_smoothed$nCG), 1), round(median(dmr_smoothed$nCG), 1),
                 round(mean(abs(dmr_smoothed$meanMethy1 - dmr_smoothed$meanMethy2)), 3),
                 round(median(abs(dmr_smoothed$meanMethy1 - dmr_smoothed$meanMethy2)), 3))
  )
  
  print(comparison_df)
  cat("\n")
  
  # Create comparison plots
  cat("Generating comparison plots...\n")
  if(!dir.exists("comparison_plots")) dir.create("comparison_plots")
  
  # Side-by-side histograms for DMR length
  png("comparison_plots/dmr_length_comparison.png", width=1200, height=600)
  par(mfrow=c(1,2))
  hist(dmr_raw$length, breaks=30, main="DMR Length - Raw Data", 
       xlab="Length (bp)", col="lightblue", xlim=c(0, max(c(dmr_raw$length, dmr_smoothed$length))))
  hist(dmr_smoothed$length, breaks=30, main="DMR Length - Smoothed Data", 
       xlab="Length (bp)", col="lightgreen", xlim=c(0, max(c(dmr_raw$length, dmr_smoothed$length))))
  dev.off()
  
  # Side-by-side histograms for CpG count
  png("comparison_plots/dmr_cpg_comparison.png", width=1200, height=600)
  par(mfrow=c(1,2))
  hist(dmr_raw$nCG, breaks=20, main="CpGs per DMR - Raw Data", 
       xlab="Number of CpGs", col="lightblue", xlim=c(0, max(c(dmr_raw$nCG, dmr_smoothed$nCG))))
  hist(dmr_smoothed$nCG, breaks=20, main="CpGs per DMR - Smoothed Data", 
       xlab="Number of CpGs", col="lightgreen", xlim=c(0, max(c(dmr_raw$nCG, dmr_smoothed$nCG))))
  dev.off()
  
  # Methylation difference comparison
  png("comparison_plots/methylation_diff_comparison.png", width=1200, height=600)
  par(mfrow=c(1,2))
  meth_diff_raw <- abs(dmr_raw$meanMethy1 - dmr_raw$meanMethy2)
  meth_diff_smooth <- abs(dmr_smoothed$meanMethy1 - dmr_smoothed$meanMethy2)
  hist(meth_diff_raw, breaks=20, main="Methylation Difference - Raw Data", 
       xlab="|Methylation Difference|", col="lightblue", 
       xlim=c(0, max(c(meth_diff_raw, meth_diff_smooth))))
  hist(meth_diff_smooth, breaks=20, main="Methylation Difference - Smoothed Data", 
       xlab="|Methylation Difference|", col="lightgreen", 
       xlim=c(0, max(c(meth_diff_raw, meth_diff_smooth))))
  dev.off()
  
  cat("Comparison plots saved to comparison_plots/ directory\n\n")
}

# Check for overlapping DMRs (basic genomic overlap)
if(nrow(dmr_raw) > 0 && nrow(dmr_smoothed) > 0) {
  cat("GENOMIC OVERLAP ANALYSIS:\n")
  cat("-------------------------\n")
  
  # Simple overlap check (same chromosome and overlapping coordinates)
  overlaps <- 0
  for(i in 1:nrow(dmr_raw)) {
    raw_chr <- dmr_raw$chr[i]
    raw_start <- dmr_raw$start[i]
    raw_end <- dmr_raw$end[i]
    
    # Check if any smoothed DMR overlaps with this raw DMR
    overlap_found <- any(
      dmr_smoothed$chr == raw_chr &
      dmr_smoothed$start <= raw_end &
      dmr_smoothed$end >= raw_start
    )
    
    if(overlap_found) overlaps <- overlaps + 1
  }
  
  cat(paste("Raw DMRs with genomic overlap in smoothed results:", overlaps, "\n"))
  cat(paste("Percentage overlap:", round(100 * overlaps / nrow(dmr_raw), 1), "%\n\n"))
}

cat("=============================================================================\n")
cat("Comparison complete!\n")
cat("=============================================================================\n")