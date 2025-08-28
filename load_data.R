

# load_data.R
# Data Loading Script for DMR Analysis

# Load required libraries
library(DSS)
library(BiocParallel)



source("sample_config.R")



load_bs_data <- function(sample_files, sample_names, coverage_threshold = 5) {
  
  cat("Loading bisulfite sequencing data files...\n")
  
  # Validate inputs
  if(length(sample_files) != length(sample_names)) {
    stop("Number of sample files must match number of sample names")
  }
  
  # Load all data files
  dat_list <- list()
  for(i in 1:length(sample_files)) {
    cat(paste("Loading", sample_files[i], "\n"))
    
    if(!file.exists(sample_files[i])) {
      stop(paste("File not found:", sample_files[i]))
    }
    
    dat_list[[i]] <- read.table(sample_files[i], header=TRUE, stringsAsFactors=FALSE)
    
    # Check data format
    required_cols <- c("chr", "pos", "N", "X")
    if(!all(required_cols %in% colnames(dat_list[[i]]))) {
      stop(paste("File", sample_files[i], "missing required columns:", 
                 paste(required_cols[!required_cols %in% colnames(dat_list[[i]])], collapse=", ")))
    }
    
    # Basic data validation
    cat(paste("  Rows:", nrow(dat_list[[i]]), "\n"))
    cat(paste("  Chromosomes:", length(unique(dat_list[[i]]$chr)), "\n"))
    cat(paste("  Mean coverage:", round(mean(dat_list[[i]]$N, na.rm=TRUE), 2), "\n"))
  }
  
  # Create BSseq object
  cat("Creating BSseq object...\n")
  BSobj <- makeBSseqData(dat_list, sampleNames = sample_names)
  
  # Filter out CpG sites with low coverage
  cat(paste("Filtering sites with coverage <", coverage_threshold, "\n"))
  initial_sites <- nrow(BSobj)
  BSobj <- BSobj[rowSums(getCoverage(BSobj, type="Cov")) >= coverage_threshold]
  final_sites <- nrow(BSobj)
  
  cat(paste("CpG sites before filtering:", initial_sites, "\n"))
  cat(paste("CpG sites after filtering:", final_sites, "\n"))
  cat(paste("Sites removed:", initial_sites - final_sites, "\n"))
  
  cat("BSseq object created successfully!\n")
  print(BSobj)
  
  return(BSobj)
}

BSobj_raw <- load_bs_data(sample_files, sample_names, coverage_threshold = 5)

cat("\nData loading complete! BSobj_raw is now available.\n")
cat("Use source('dmr_analysis.R') to run the analysis.\n")