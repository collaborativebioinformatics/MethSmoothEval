# Load required libraries
library(bsseq)
library(GenomicRanges)
library(data.table)

# 1. Define Paths and Parameters
# -----------------------------------------------------------
manifest_file <- "list_bedgraph"  # File containing list of BEDGraph paths
base_dir <- "/lustre07/scratch/senthil/Hackthon_2025/"  # Base directory for paths
output_dir <- paste0(base_dir, "analysis/BSmooth_results/")  # Output directory

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
  
  # Construct full path to input file
  full_path <- paste0(base_dir, file_path)
  cat("\n" , rep("=", 70), "\n", sep = "")
  cat("Processing sample:", file_path, "\n")
  cat(rep("=", 70), "\n\n", sep = "")
  
  # Extract sample name from file path for output naming
  file_name <- basename(file_path)
  sample_name <- gsub("\\.bedGraph$", "", file_name)
  output_bedgraph <- paste0(output_dir, sample_name, "_chr22_MinCov5_smoothed.bedGraph")
  
  cat("Output will be saved to:", output_bedgraph, "\n")
  
  # 4. Load and Filter Data for Chr22
  # -----------------------------------------------------------
  cat("Importing and filtering for chr22...\n")
  
  tryCatch({
    # Method 1: Use grep for efficient loading (recommended)
    bg_data <- fread(cmd = paste("grep -w '^chr22'", full_path), 
                     sep = "\t", 
                     header = FALSE,
                     col.names = c("chr", "start", "end", "percent_meth", "count_meth", "count_unmeth"))
    
    # Method 2: Alternative if grep fails - read all then subset
    if (nrow(bg_data) == 0) {
      cat("Grep method failed, trying alternative method...\n")
      bg_data <- fread(full_path, sep = "\t", header = FALSE)
      colnames(bg_data) <- c("chr", "start", "end", "percent_meth", "count_meth", "count_unmeth")
      bg_data <- bg_data[chr == "chr22", ]
    }
    
    cat("Number of CpG sites on chr22:", nrow(bg_data), "\n")
    
    if (nrow(bg_data) == 0) {
      cat("Warning: No data found for chr22 in", file_path, "Skipping.\n")
      next
    }
    
    # 5. Apply Coverage Filter
    # -----------------------------------------------------------
    cat("Applying coverage filter (>= 5 reads)...\n")
    bg_data$coverage <- bg_data$count_meth + bg_data$count_unmeth
    bg_data_filtered <- bg_data[coverage >= 5, ]
    
    cat("Sites after coverage filter:", nrow(bg_data_filtered), "\n")
    cat("Percentage retained:", round(nrow(bg_data_filtered)/nrow(bg_data) * 100, 2), "%\n")
    
    if (nrow(bg_data_filtered) < 50) {
      cat("Warning: Too few CpG sites (< 50) on chr22 after filtering. Skipping", file_path, "\n")
      next
    }
    
    # 6. Prepare Data for BSseq Object
    # -----------------------------------------------------------
    cat("Preparing data matrices...\n")
    gr <- GRanges(seqnames = bg_data_filtered$chr,
                  ranges = IRanges(start = bg_data_filtered$start + 1,
                                   end = bg_data_filtered$end),
                  strand = "*")
    
    M_matrix <- as.matrix(bg_data_filtered$count_meth)
    Cov_matrix <- as.matrix(bg_data_filtered$coverage)
    
    # 7. Create BSseq Object
    # -----------------------------------------------------------
    cat("Creating BSseq object...\n")
    bs_obj <- BSseq(chr = as.character(seqnames(gr)),
                    pos = start(gr),
                    M = M_matrix,
                    Cov = Cov_matrix,
                    sampleNames = sample_name)
    
    cat("Sample coverage summary:\n")
    print(summary(getCoverage(bs_obj)))
    
    # 8. Sort and Smooth
    # -----------------------------------------------------------
    cat("Sorting and running BSmooth...\n")
    bs_obj_sorted <- sort(bs_obj)
    
    bs_obj_smoothed <- BSmooth(bs_obj_sorted,
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
      value = round(smoothed_meth[, 1] * 100, 2)  # Convert to percentage 0-100
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
    
  }, error = function(e) {
    cat("❌ Error processing", file_path, ":\n")
    cat("Error message:", e$message, "\n")
    cat("Skipping to next sample.\n")
  })
}

cat("\n" , rep("=", 70), "\n", sep = "")
cat("Batch processing completed!\n")
cat("BEDGraph results saved in:", output_dir, "\n")
cat(rep("=", 70), "\n", sep = "")
