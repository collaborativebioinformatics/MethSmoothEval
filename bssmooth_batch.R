#!/usr/bin/env Rscript

library(bsseq)
library(rtracklayer)
library(optparse)

# Parse command line arguments
option_list <- list(
  make_option(c("-f", "--file-list"), type="character", help="Text file with list of input files"),
  make_option(c("-o", "--output-dir"), type="character", default="smoothed_output", help="Output directory"),
  make_option(c("--ns"), type="integer", default=70, help="Number of sites for smoothing"),
  make_option(c("--h"), type="integer", default=1000, help="Bandwidth in bp"),
  make_option(c("--cores"), type="integer", default=4, help="Number of cores")
)

opt <- parse_args(OptionParser(option_list=option_list))

# Create output directory
dir.create(opt$output_dir, showWarnings = FALSE)

# Read file list
files <- readLines(opt$file_list)
files <- files[files != ""]  # Remove empty lines

cat("Processing", length(files), "files...\n")

# Process each file individually
for(i in 1:length(files)) {
  cat("Processing file", i, "of", length(files), ":", files[i], "\n")
  
  tryCatch({
    # Read pb-cpg-tools format (skip header lines starting with ##)
    if(grepl("\\.gz$", files[i])) {
      con <- gzfile(files[i])
      cat("  -> Reading gzipped file\n")
    } else {
      con <- file(files[i])
      cat("  -> Reading uncompressed file\n")
    }
    
    # Read all lines and filter out header
    all_lines <- readLines(con)
    close(con)
    
    # Remove comment lines starting with ##
    data_lines <- all_lines[!grepl("^##", all_lines)]
    
    # Parse the actual data
    data <- read.table(text = paste(data_lines, collapse = "\n"), 
                      header = TRUE, stringsAsFactors = FALSE)
    
    # Convert to BSseq format using pb-cpg-tools columns
    M <- matrix(data$est_mod_count, ncol=1)
    Cov <- matrix(data$cov, ncol=1)
    pos <- data$begin
    chr <- data$chrom
    
    # Create raw BSseq object
    sample_name <- tools::file_path_sans_ext(basename(files[i]))
    if(grepl("\\.gz$", files[i])) {
      sample_name <- tools::file_path_sans_ext(sample_name)
    }
    
    BSobj_raw <- BSseq(chr = chr, pos = pos, M = M, Cov = Cov,
                       sampleNames = sample_name)
    
    # Run BSsmooth
    cat("  -> Running BSsmooth...\n")
    BSobj_smooth <- BSmooth(BSobj_raw, ns=opt$ns, h=opt$h, verbose=FALSE)
    
    # Save RDS files
    base_name <- sample_name
    raw_output <- file.path(opt$output_dir, paste0(base_name, "_raw.rds"))
    smooth_output <- file.path(opt$output_dir, paste0(base_name, "_smoothed.rds"))
    
    saveRDS(BSobj_raw, raw_output)
    saveRDS(BSobj_smooth, smooth_output)
    
    # Export BigWig files
    cat("  -> Exporting BigWig files...\n")
    
    # Raw BigWig
    raw_meth_values <- getMeth(BSobj_raw, type = "raw")
    raw_gr_with_scores <- granges(BSobj_raw)
    raw_gr_with_scores$score <- raw_meth_values[,1]
    #raw_bw_output <- file.path(opt$output_dir, paste0(base_name, "_raw.bw"))
    #export.bw(raw_gr_with_scores, raw_bw_output)
    
    # Smoothed BigWig  
    smooth_meth_values <- getMeth(BSobj_smooth, type = "smooth")
    smooth_gr_with_scores <- granges(BSobj_smooth)
    smooth_gr_with_scores$score <- smooth_meth_values[,1]
    smooth_bw_output <- file.path(opt$output_dir, paste0(base_name, "_smoothed.bw"))
    export.bw(smooth_gr_with_scores, smooth_bw_output)
    
    cat("  -> Raw RDS:", raw_output, "\n")
    cat("  -> Smooth RDS:", smooth_output, "\n") 
    #cat("  -> Raw BigWig:", raw_bw_output, "\n")
    cat("  -> Smooth BigWig:", smooth_bw_output, "\n")
    
  }, error = function(e) {
    cat("  ERROR processing", files[i], ":", e$message, "\n")
  })
}

cat("Done! All files in:", opt$output_dir, "\n")
