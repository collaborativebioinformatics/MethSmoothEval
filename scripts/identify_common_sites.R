#!/usr/bin/env Rscript

# identify_common_sites.R
# Identifies common CpG sites across multiple methylation technologies
# and generates coverage statistics

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
})

# ============================================================================
# Command-line argument parsing
# ============================================================================

option_list <- list(
  make_option(c("-p", "--pacbio"), type="character", default=NULL,
              help="PacBio methylation file (BED format: chr, start, end, cov, mod_count)",
              metavar="FILE"),
  
  make_option(c("-e", "--emseq"), type="character", default=NULL,
              help="EM-seq methylation file (Bismark cov format: chr, start, end, meth%, count_meth, count_unmeth)",
              metavar="FILE"),
  
  make_option(c("-t", "--trueseq"), type="character", default=NULL,
              help="TruSeq WGBS methylation file (Bismark cov format)",
              metavar="FILE"),
  
  make_option(c("-x", "--trueseqiox"), type="character", default=NULL,
              help="TruSeq IOX methylation file (Bismark cov format)",
              metavar="FILE"),
  
  make_option(c("-o", "--output"), type="character", default="common_sites.rds",
              help="Output file for common sites (RDS format) [default=%default]",
              metavar="FILE"),
  
  make_option(c("-s", "--stats"), type="character", default="coverage_stats.txt",
              help="Output file for coverage statistics [default=%default]",
              metavar="FILE"),
  
  make_option(c("-c", "--min_coverage"), type="integer", default=5,
              help="Minimum coverage threshold to include a site [default=%default]",
              metavar="INT"),
  
  make_option(c("-v", "--verbose"), action="store_true", default=FALSE,
              help="Print verbose output")
)

opt_parser <- OptionParser(option_list=option_list,
                          description="Identify common CpG sites across methylation technologies")
opt <- parse_args(opt_parser)

# Check that at least 2 input files are provided
input_files <- c(opt$pacbio, opt$emseq, opt$trueseq, opt$trueseqiox)
input_names <- c("PacBio", "EMseq", "TruSeq", "TruSeqIOX")
provided <- !sapply(input_files, is.null)

if(sum(provided) < 2) {
  stop("At least 2 input files must be provided for comparison")
}

if(opt$verbose) {
  cat("=== Identify Common CpG Sites ===\n")
  cat("Input files:\n")
  for(i in which(provided)) {
    cat(sprintf("  %s: %s\n", input_names[i], input_files[i]))
  }
  cat(sprintf("Minimum coverage: %d\n", opt$min_coverage))
  cat(sprintf("Output: %s\n", opt$output))
  cat(sprintf("Stats: %s\n", opt$stats))
  cat("\n")
}

# ============================================================================
# Data loading and processing functions
# ============================================================================

#' Read PacBio format (aligned_bam_to_cpg_scores output)
#' Format: chrom, begin, end, cov, mod_count
#' Coordinates: 0-based start, 1-based end
read_pacbio <- function(file, min_cov = 5) {
  if(opt$verbose) cat(sprintf("Reading PacBio file: %s\n", file))
  

  #chrom  begin   end     mod_score       type    cov     mod_count       unmod_count     avg_mod_score   avg_unmod_score
  #chr22   10514791        10514792        60.0    Total   5       3       2       0.864   0.007

  dt <- fread(file, col.names = c("chr", "start", "end", "mod_score", "type", "N", "X", "unmod_count", "avg_mod_score", "avg_unmod_score"))
  
  # Filter by coverage
  dt <- dt[N >= min_cov]
  
  # Convert to 1-based coordinates for consistency
  dt[, pos := start + 1]
  
  # Calculate beta value
  dt[, beta := X / N]
  
  # Return standardized format
  return(dt[, .(chr, pos, N, X, beta)])
}

#' Read Bismark coverage format (EM-seq, TruSeq, TruSeqIOX)
#' Format: chr, start, end, methylation%, count_methylated, count_unmethylated
#' Coordinates: 1-based start and end
read_bismark <- function(file, min_cov = 5) {
  if(opt$verbose) cat(sprintf("Reading Bismark format file: %s\n", file))
  
  dt <- fread(file, col.names = c("chr", "start", "end", "meth_pct", "X", "unmeth"))
  
  # Calculate total coverage
  dt[, N := X + unmeth]
  
  # Filter by coverage
  dt <- dt[N >= min_cov]
  
  # Use start position as the CpG position (1-based)
  dt[, pos := start]
  
  # Calculate beta value (methylation proportion)
  dt[, beta := X / N]
  
  # Return standardized format
  return(dt[, .(chr, pos, N, X, beta)])
}

#' Standardize chromosome names
standardize_chr <- function(dt) {
  # Remove "chr" prefix if present for consistency
  dt[, chr := gsub("^chr", "", chr)]
  
  # Filter to standard chromosomes (1-22, X, Y, MT/M)
  standard_chr <- c(as.character(1:22), "X", "Y", "M", "MT")
  dt <- dt[chr %in% standard_chr]
  
  # Standardize MT to M
  dt[chr == "MT", chr := "M"]
  
  return(dt)
}

# ============================================================================
# Load all provided datasets
# ============================================================================

datasets <- list()
tech_names <- character()

if(!is.null(opt$pacbio)) {
  datasets[["PacBio"]] <- standardize_chr(read_pacbio(opt$pacbio, opt$min_coverage))
  tech_names <- c(tech_names, "PacBio")
}

if(!is.null(opt$emseq)) {
  datasets[["EMseq"]] <- standardize_chr(read_bismark(opt$emseq, opt$min_coverage))
  tech_names <- c(tech_names, "EMseq")
}

if(!is.null(opt$trueseq)) {
  datasets[["TruSeq"]] <- standardize_chr(read_bismark(opt$trueseq, opt$min_coverage))
  tech_names <- c(tech_names, "TruSeq")
}

if(!is.null(opt$trueseqiox)) {
  datasets[["TruSeqIOX"]] <- standardize_chr(read_bismark(opt$trueseqiox, opt$min_coverage))
  tech_names <- c(tech_names, "TruSeqIOX")
}

if(opt$verbose) {
  cat("\nData loaded:\n")
  for(tech in names(datasets)) {
    cat(sprintf("  %s: %s CpG sites\n", tech, format(nrow(datasets[[tech]]), big.mark=",")))
  }
  cat("\n")
}

# ============================================================================
# Calculate coverage statistics for each technology
# ============================================================================

coverage_stats <- rbindlist(lapply(tech_names, function(tech) {
  dt <- datasets[[tech]]
  
  data.table(
    Technology = tech,
    Total_Sites = nrow(dt),
    Mean_Coverage = round(mean(dt$N), 2),
    Median_Coverage = median(dt$N),
    Q25_Coverage = quantile(dt$N, 0.25),
    Q75_Coverage = quantile(dt$N, 0.75),
    Min_Coverage = min(dt$N),
    Max_Coverage = max(dt$N),
    Mean_Beta = round(mean(dt$beta), 3),
    Median_Beta = round(median(dt$beta), 3),
    Sites_5x = sum(dt$N >= 5),
    Sites_10x = sum(dt$N >= 10),
    Sites_20x = sum(dt$N >= 20),
    Sites_50x = sum(dt$N >= 50)
  )
}))

# ============================================================================
# Identify common sites across all technologies
# ============================================================================

if(opt$verbose) cat("Identifying common CpG sites...\n")

# Set keys for efficient merging
for(tech in names(datasets)) {
  setkey(datasets[[tech]], chr, pos)
}

# Merge all datasets iteratively
common_sites <- datasets[[1]]
setnames(common_sites, 
         c("N", "X", "beta"), 
         paste0(c("N", "X", "beta"), "_", names(datasets)[1]))

for(i in 2:length(datasets)) {
  tech_name <- names(datasets)[i]
  tech_data <- copy(datasets[[i]])
  
  setnames(tech_data,
           c("N", "X", "beta"),
           paste0(c("N", "X", "beta"), "_", tech_name))
  
  common_sites <- merge(common_sites, tech_data, by = c("chr", "pos"))
  
  if(opt$verbose) {
    cat(sprintf("  After merging %s: %s common sites\n", 
                tech_name, format(nrow(common_sites), big.mark=",")))
  }
}

# Add common site count to coverage stats
coverage_stats[, Common_Sites := nrow(common_sites)]

# Calculate pairwise overlap statistics
if(length(datasets) > 1) {
  pairwise_overlaps <- rbindlist(lapply(1:(length(datasets)-1), function(i) {
    rbindlist(lapply((i+1):length(datasets), function(j) {
      tech1 <- names(datasets)[i]
      tech2 <- names(datasets)[j]
      
      dt1 <- datasets[[tech1]]
      dt2 <- datasets[[tech2]]
      
      merged <- merge(dt1, dt2, by = c("chr", "pos"))
      
      data.table(
        Tech1 = tech1,
        Tech2 = tech2,
        Sites_Tech1 = nrow(dt1),
        Sites_Tech2 = nrow(dt2),
        Overlap = nrow(merged),
        Jaccard = round(nrow(merged) / (nrow(dt1) + nrow(dt2) - nrow(merged)), 3),
        Overlap_Pct_Tech1 = round(100 * nrow(merged) / nrow(dt1), 2),
        Overlap_Pct_Tech2 = round(100 * nrow(merged) / nrow(dt2), 2)
      )
    }))
  }))
}

# ============================================================================
# Save outputs
# ============================================================================

if(opt$verbose) cat("\nSaving outputs...\n")

# Save common sites as RDS
output_data <- list(
  common_sites = common_sites,
  individual_datasets = datasets,
  tech_names = tech_names,
  min_coverage = opt$min_coverage
)

saveRDS(output_data, opt$output)
if(opt$verbose) cat(sprintf("  Common sites saved to: %s\n", opt$output))

# Save coverage statistics
fwrite(coverage_stats, opt$stats, sep = "\t")
if(opt$verbose) cat(sprintf("  Coverage stats saved to: %s\n", opt$stats))

# Save pairwise overlaps if applicable
if(exists("pairwise_overlaps")) {
  overlap_file <- sub("\\.(txt|tsv)$", "_pairwise.txt", opt$stats)
  fwrite(pairwise_overlaps, overlap_file, sep = "\t")
  if(opt$verbose) cat(sprintf("  Pairwise overlaps saved to: %s\n", overlap_file))
}

# ============================================================================
# Print summary
# ============================================================================

cat("\n=== Summary ===\n")
cat(sprintf("Common sites across all %d technologies: %s\n\n", 
            length(datasets), format(nrow(common_sites), big.mark=",")))

cat("Coverage Statistics:\n")
print(coverage_stats[, .(Technology, Total_Sites, Mean_Coverage, Median_Coverage, 
                        Mean_Beta, Common_Sites)])

if(exists("pairwise_overlaps")) {
  cat("\nPairwise Overlaps:\n")
  print(pairwise_overlaps[, .(Tech1, Tech2, Overlap, Jaccard, 
                              Overlap_Pct_Tech1, Overlap_Pct_Tech2)])
}

cat("\nDone!\n")