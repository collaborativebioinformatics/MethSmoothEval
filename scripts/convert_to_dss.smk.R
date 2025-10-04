#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(data.table)
})





# Parse command line arguments
option_list <- list(
  make_option("--input", type="character", help="Input masked CpG BED file"),
  make_option("--output", type="character", help="Output DSS bed file"),
  make_option("--sample", type="character", help="Sample name"),
  make_option(c("-c", "--min_coverage"), type="integer", default=5,
              help="Minimum coverage threshold [default %default]"),
  make_option(c("-v", "--verbose"), action="store_true", default=FALSE,
              help="Print verbose output")
)

opt <- parse_args(OptionParser(option_list=option_list))

# Validate arguments
if(is.null(opt$input) || is.null(opt$output) || is.null(opt$sample)) {
  stop("Required arguments missing. Usage: Rscript convert_to_dss.R --input INPUT --output OUTPUT --sample SAMPLE")
}


if (opt$verbose) {
  message("Converting ", toupper(opt$sample), " methylation data to DSS format...")
  message("Input file: ", opt$input)
  message("Output file: ", opt$output)
  message("Minimum coverage: ", opt$min_coverage)
  }



convert_pacbio_count <- function(file, outfile, min_coverage = 5, verbose = FALSE) {
  if (verbose) message("Reading PacBio file...")
  
  df <- fread(file, header=TRUE)
  setnames(df, gsub("^#", "", colnames(df))) # rename '#chrom' → 'chrom'
  
  required_cols <- c("chrom","begin","end","cov","mod_count")
  if (!all(required_cols %in% colnames(df))) {
    stop("PacBio file must contain columns: ", paste(required_cols, collapse=", "))
  }
  
  if (verbose) {
    message("Initial sites: ", nrow(df))
    message("Mean coverage: ", round(mean(df$cov), 2))
  }
  
  dss_df <- df %>%
    mutate(chr = chrom,
           pos = begin + 1,       # convert 0-based to 1-based
           N   = cov,
           X   = mod_count) %>%
    filter(N >= min_coverage) %>%  # Apply coverage filter
    select(chr, pos, N, X)
  
  if (verbose) {
    message("Sites after filtering (coverage >= ", min_coverage, "): ", nrow(dss_df))
    message("Mean methylation rate: ", round(mean(dss_df$X / dss_df$N), 3))
  }
  
  write.table(dss_df, file=outfile, sep="\t", quote=FALSE,
              row.names=FALSE, col.names=TRUE)
  
  return(nrow(dss_df))
}


converted_bed <- convert_pacbio_count(opt$input, opt$output, opt$min_coverage, opt$verbose)



message("Conversion complete: ", converted_bed)