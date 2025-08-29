################################################################################
## The bsseq User’s Guide
# https://www.bioconductor.org/packages/devel/bioc/vignettes/bsseq/inst/doc/bsseq.html

################################################################################
## conda activate bs
# conda install bioconda::bioconductor-bsseq
# conda install bioconda::bioconductor-dss
# conda install bioconda::bioconductor-annotationhub
# conda install bioconda::bioconductor-annotatr
# conda install bioconda::bioconductor-txdb.hsapiens.ucsc.hg38.knowngene
# conda install bioconda::bioconductor-org.hs.eg.db

################################################################################
# modified script from Senthil
# Batch BSmooth Analysis for ONT modkit pileups
# Filters: chromosome 22, coverage >= 5
# Output: BEDGraph files with smoothed methylation values

# Load required libraries
library(bsseq)
library(GenomicRanges)
library(data.table)

# 1. Define Paths and Parameters
# -----------------------------------------------------------
base_dir <- "/scratch/eger/projects/MethSmoothEval"
manifest_file <- "/scratch/eger/projects/MethSmoothEval/analysis/Sarah_analysis/human_ONT/file_list_modkit_pileups.txt"
output_dir <- paste0(base_dir, "/analysis/Sarah_analysis/human_ONT/BSmooth_results/")  # Output directory

# Create the output directory if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  cat("Created output directory:", output_dir, "\n")
}

# 2. Read the Manifest File
# -----------------------------------------------------------
cat("Reading manifest file:", manifest_file, "\n")
file_list <- readLines(manifest_file)
cat("Found", length(file_list), "samples to process:\n")
print(file_list)

# 3. Process Each Sample
# -----------------------------------------------------------
for (file_path in file_list) {
  
    # Extract sample name from file path for output naming
    file_name <- basename(file_path)
    parts <- strsplit(file_name, "\\.")[[1]]
    sample_name <- paste(parts[1:2], collapse = ".")
    output_bedgraph <- paste0(output_dir, sample_name, ".chr22.MinCov5_smoothed.bedGraph")
    
    cat("Output will be saved to:", output_bedgraph, "\n")
    
    # 4. Load and Filter Data for Chr22
    # -----------------------------------------------------------
    cat("Importing and filtering for chr22...\n")

    bs <- read.modkit(file_path, rmZeroCov = FALSE, strandCollapse = FALSE)
    
    # give the BSseq object clean, unique rownames if needed
    if (is.null(rownames(bs)) || anyNA(rownames(bs))) {
      rr <- rowRanges(bs)
      rn <- paste0(as.character(seqnames(rr)), ":", start(rr))
      rownames(bs) <- make.unique(rn)
    }
    bg_data <- bs
    cat("Number of CpG sites on chr22:", nrow(bg_data), "\n")
    
    # 5. Apply Coverage Filter
    # -----------------------------------------------------------
    cat("Applying coverage filter (>= 5 reads)...\n")
    
    min_cov <- 5
    
    # keep loci that meet thresholds
    cov <- getCoverage(bg_data)
    keepLoci <- which(cov[, 1] >= min_cov)
    bg_data_filtered <- bg_data[keepLoci, ]
    
    cat("Sites after coverage filter:", nrow(bg_data_filtered), "\n")
    cat("Percentage retained:", round(nrow(bg_data_filtered)/nrow(bg_data) * 100, 2), "%\n")
    
    if (nrow(bg_data_filtered) < 50) {
      cat("Warning: Too few CpG sites (< 50) on chr22 after filtering. Skipping", file_path, "\n")
      next
    }
    
    # 8. Sort and Smooth
    # -----------------------------------------------------------
    cat("Sorting and running BSmooth...\n")
    bs_obj_smoothed <- BSmooth(bg_data_filtered,
                               ns = 70,
                               h = 1000,
                               maxGap = 10^8,
                               verbose = TRUE)
    
    # 9. Extract and Save Results as BEDGraph
    # -----------------------------------------------------------
    cat("Extracting smoothed values for BEDGraph format...\n")
    smoothed_meth <- getMeth(bs_obj_smoothed, type = "smooth", what = "perBase")
    smoothed_gr <- granges(bs_obj_smoothed)
    
    # Create BEDGraph format data frame
    # BEDGraph format: chrom, start, end, value
    # Note: BEDGraph uses 0-based start coordinates, so we subtract 1 from start
    output_df <- data.frame(
      chrom = as.character(seqnames(smoothed_gr)),
      start = start(smoothed_gr) - 1,  # Convert to 0-based for BED format
      end = end(smoothed_gr),
      value = round(smoothed_meth[, 1] * 100, 4)  # Convert to percentage 
    )
    
    # 10. Write BEDGraph file with proper header
    # -----------------------------------------------------------
    cat("Writing BEDGraph file:", output_bedgraph, "\n")
    
    # Write track header line (optional but useful for genome browsers)
    track_line <- paste0('track type=bedGraph name="', sample_name, 
                         '" description="BSmooth smoothed methylation values for ', sample_name,
                         ' (chr22, cov>=5)" visibility=full color=0,100,200 altColor=200,100,0 priority=20')
    
    # Write to file
    writeLines(track_line, output_bedgraph)
    fwrite(output_df, file = output_bedgraph, sep = "\t", quote = FALSE, 
           col.names = FALSE, append = TRUE)
    
    cat("✓ Successfully processed", sample_name, "and saved as BEDGraph\n")
    
}

cat("\n" , rep("=", 70), "\n", sep = "")
cat("Batch processing completed!\n")
cat("BEDGraph results saved in:", output_dir, "\n")
cat(rep("=", 70), "\n", sep = "")
