


run_methSmooth_analysis <- function(sample1_file, sample2_file, 
                                     sample1_name, sample2_name,
                                     p_threshold = 0.05, 
                                     delta_threshold = 0.1,
                                     verbose = TRUE) {
  # Function to perform both smoothed and unsmoothed methylation analysis
  # 
  # Args:
  #   sample1_file: Path to bed file for sample 1
  #   sample2_file: Path to bed file for sample 2
  #   sample1_name: Name for sample 1 in the analysis
  #   sample2_name: Name for sample 2 in the analysis
  #   p_threshold: P-value threshold for significance (default: 0.05)
  #   delta_threshold: Methylation difference threshold (default: 0.1 or 10%)
  #   verbose: Whether to print analysis progress and results (default: TRUE)
  #
  # Returns:
  #   A list containing all relevant objects:
  #     - BSobj: The BSseq object created from input files
  #     - dmlTest_smoothed: DMLtest results with smoothing
  #     - dmlTest_unsmoothed: DMLtest results without smoothing
  #     - dmls_smoothed: Significant DMLs from smoothed analysis
  #     - dmls_unsmoothed: Significant DMLs from unsmoothed analysis
  #     - overlap_results: Information about overlapping DMLs between methods
  
  # Read data files
  if (verbose) cat("Reading methylation data files...\n")
  sample1_data <- read.table(sample1_file, header = TRUE, stringsAsFactors = FALSE)
  sample2_data <- read.table(sample2_file, header = TRUE, stringsAsFactors = FALSE)
  
  # Create BSseq object
  BSobj <- makeBSseqData(
    list(sample1_data, sample2_data),
    sampleNames = c(sample1_name, sample2_name)
  )
  
  if (verbose) {
    cat("BSseq object created:\n")
    print(BSobj)
  }
  
  # Run DML test with smoothing
  if (verbose) cat("\nRunning DML test with smoothing...\n")
  dmlTest_smoothed <- DMLtest(BSobj, 
                     group1 = c(sample1_name), 
                     group2 = c(sample2_name),
                     smoothing = TRUE)
  
  if (verbose) cat("Smoothed DML test completed. Total tests performed:", nrow(dmlTest_smoothed), "\n")
  
  # Run DML test without smoothing
  if (verbose) cat("\nRunning DML test without smoothing...\n")
  dmlTest_unsmoothed <- DMLtest(BSobj, 
                       group1 = c(sample1_name), 
                       group2 = c(sample2_name),
                       smoothing = FALSE,
                       equal.disp = TRUE)
  
  if (verbose) {
    cat("Unsmoothed DML test completed. Total tests performed:", nrow(dmlTest_unsmoothed), "\n")
    cat("\nP-value summaries:\n")
    cat("Smoothed analysis:\n")
    print(summary(dmlTest_smoothed$pval))
    cat("\nUnsmoothed analysis:\n")
    print(summary(dmlTest_unsmoothed$pval))
  }
  
  # Call significant DMLs
  if (verbose) cat("\nCalling significant DMLs...\n")
  
  dmls_smoothed <- callDML(dmlTest_smoothed, 
                          p.threshold = p_threshold,
                          delta = delta_threshold)
  
  dmls_unsmoothed <- callDML(dmlTest_unsmoothed, 
                            p.threshold = p_threshold,
                            delta = delta_threshold)
  
  # Analyze results
  overlap_results$dmls <- list()
  
  if (verbose) {
    cat("\n--- DML Results Comparison ---\n")
    cat("Smoothed analysis - Significant DMLs:", nrow(dmls_smoothed), "\n")
    cat("Unsmoothed analysis - Significant DMLs:", nrow(dmls_unsmoothed), "\n")
    cat("Criteria: p-value <", p_threshold, "and |methylation difference| >", delta_threshold, "\n")
  }
  
  # Find overlapping DMLs
  if (nrow(dmls_smoothed) > 0 && nrow(dmls_unsmoothed) > 0) {
    smoothed_positions <- paste(dmls_smoothed$chr, dmls_smoothed$pos, sep = "_")
    unsmoothed_positions <- paste(dmls_unsmoothed$chr, dmls_unsmoothed$pos, sep = "_")
    overlapping_positions <- intersect(smoothed_positions, unsmoothed_positions)
    
    # Store overlap information
    overlap_results$dmls$positions <- overlapping_positions
    overlap_results$dmls$count <- length(overlapping_positions)
    overlap_results$dmls$pct_of_smoothed <- round(length(overlapping_positions) / nrow(dmls_smoothed) * 100, 1)
    overlap_results$dmls$pct_of_unsmoothed <- round(length(overlapping_positions) / nrow(dmls_unsmoothed) * 100, 1)
    
    if (verbose) {
      cat("Overlapping significant DMLs between methods:", overlap_results$dmls$count, "\n")
      cat("Overlap as % of smoothed DMLs:", overlap_results$dmls$pct_of_smoothed, "%\n")
      cat("Overlap as % of unsmoothed DMLs:", overlap_results$dmls$pct_of_unsmoothed, "%\n")
    }
  }
  
  # Print top results if available
  if (verbose) {
    if (nrow(dmls_smoothed) > 0) {
      cat("\nTop 10 Smoothed DMLs:\n")
      print(head(dmls_smoothed, 10))
    }
    
    if (nrow(dmls_unsmoothed) > 0) {
      cat("\nTop 10 Unsmoothed DMLs:\n")
      print(head(dmls_unsmoothed, 10))
    }
  }
  
  # Save top DMLs to CSV files
  if(nrow(dmls_smoothed) > 0) {
  smoothed_csv_file <- paste0(sample1_name, "_vs_", sample2_name, "_smoothed_DMLs.csv")
  write.csv(dmls_smoothed, file = smoothed_csv_file, row.names = FALSE)
  if(verbose) cat("\nSmoothed DMLs saved to:", smoothed_csv_file, "\n")
  }

  if(nrow(dmls_unsmoothed) > 0) {
  unsmoothed_csv_file <- paste0(sample1_name, "_vs_", sample2_name, "_unsmoothed_DMLs.csv")
  write.csv(dmls_unsmoothed, file = unsmoothed_csv_file, row.names = FALSE)
  if(verbose) cat("Unsmoothed DMLs saved to:", unsmoothed_csv_file, "\n")
  }


  # Call DMRs for both methods
  if (verbose) cat("\nCalling DMRs...\n")
  dmrs_smoothed <- callDMR(dmlTest_smoothed, 
                        p.threshold = p_threshold,     
                        delta = delta_threshold)
  dmrs_unsmoothed <- callDMR(dmlTest_unsmoothed, 
                          p.threshold = p_threshold,     
                          delta = delta_threshold)

	# Store DMR information in results
	overlap_results$dmrs <- list()

	if (verbose) {
	  if(nrow(dmrs_smoothed) > 0) {
	    cat("\nTop 10 Smoothed DMRs:\n")
	    print(head(dmrs_smoothed, 10))
	  }
	  if(nrow(dmrs_unsmoothed) > 0) {
	    cat("\nTop 10 Unsmoothed DMRs:\n")
	    print(head(dmrs_unsmoothed, 10))
	  }
	}

	# Save top DMRs to CSV files
	smoothed_dmr_csv_file <- paste0(sample1_name, "_vs_", sample2_name, "_smoothed_DMRs.csv")
	unsmoothed_dmr_csv_file <- paste0(sample1_name, "_vs_", sample2_name, "_unsmoothed_DMRs.csv")

	if(nrow(dmrs_smoothed) > 0) {
	  write.csv(dmrs_smoothed, file = smoothed_dmr_csv_file, row.names = FALSE)
	  if(verbose) cat("\nSmoothed DMRs saved to:", smoothed_dmr_csv_file, "\n")
	}

	if(nrow(dmrs_unsmoothed) > 0) {
	  write.csv(dmrs_unsmoothed, file = unsmoothed_dmr_csv_file, row.names = FALSE)
	  if(verbose) cat("Unsmoothed DMRs saved to:", unsmoothed_dmr_csv_file, "\n")
	}

	if(nrow(dmrs_smoothed) > 0 && nrow(dmrs_unsmoothed) > 0) {
	  if (verbose) {
	    cat("\n--- DMR Comparison ---\n")
	    cat("Method comparison:\n")
	    cat("Smoothed DMRs: ", nrow(dmrs_smoothed), "\n")
	    cat("Unsmoothed DMRs: ", nrow(dmrs_unsmoothed), "\n")
	  }
	  
	  # Check for overlapping DMRs (simplified overlap check)
	  smoothed_regions <- paste(dmrs_smoothed$chr, dmrs_smoothed$start, dmrs_smoothed$end, sep = "_")
	  unsmoothed_regions <- paste(dmrs_unsmoothed$chr, dmrs_unsmoothed$start, dmrs_unsmoothed$end, sep = "_")
	  overlapping_dmrs <- intersect(smoothed_regions, unsmoothed_regions)
	  
	  overlap_results$dmrs$count <- length(overlapping_dmrs)
	  overlap_results$dmrs$regions <- overlapping_dmrs
	  
	  if (verbose) {
	    cat("Identical DMRs between methods:", length(overlapping_dmrs), "\n")
	  }
	}



  # Return all objects and tables for downstream analysis
  return(list(
    BSobj = BSobj,
    dmlTest_smoothed = dmlTest_smoothed,
    dmlTest_unsmoothed = dmlTest_unsmoothed,
    dmls_smoothed = dmls_smoothed,
    dmls_unsmoothed = dmls_unsmoothed,
    dmrs_smoothed = dmrs_smoothed,
    dmrs_unsmoothed = dmrs_unsmoothed,
    overlap_results = overlap_results
  ))
}




create_plots <- function(methyl_results, 
                         output_dir = "./plots",
                         verbose = TRUE) {
  # Create plots from methylation analysis results
  #
  # Args:
  #   methyl_results: The list object returned by run_methylation_analysis function
  #   output_dir: Directory to save plot images (default: "./plots")
  #   verbose: Whether to print progress messages (default: TRUE)
  #
  # Returns:
  #   A list containing plot functions for reuse
  
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
  
  # Load necessary libraries
  if (!requireNamespace("scales", quietly = TRUE)) {
    install.packages("scales")
  }
  library(scales)  # For alpha function
  
  # List to store all plots
  plots <- list()
  
  if (verbose) cat("Creating combined axis visualizations...\n")

	  # ===== 1. DML P-value Distribution - Overlaid =====
	if (verbose) cat("Creating overlaid DML p-value distribution...\n")

	# Create file paths with sample names
	grid_plot_path <- file.path(output_dir, paste0(sample1_name, "_vs_", sample2_name, "_DML_pvalue_overlaid_grid.png"))
	wide_plot_path <- file.path(output_dir, paste0(sample1_name, "_vs_", sample2_name, "_DML_pvalue_overlaid_wide.png"))

	# Create histogram data for both analyses
	hist_smooth <- hist(dmlTest_smoothed$pval, breaks = 50, plot = FALSE)
	hist_unsmooth <- hist(dmlTest_unsmoothed$pval, breaks = 50, plot = FALSE)

	# Find common y-axis limits
	max_count <- max(c(hist_smooth$counts, hist_unsmooth$counts))

	# Function to create the plot (for reuse)
	create_pval_dist_plot <- function() {
	  # White background
	  par(bg = "white", mar = c(5, 5, 4, 2))
	  
	  # Plot smoothed first (background)
	  hist(dmlTest_smoothed$pval, breaks = 50,
	       main = paste0("DML P-value Distribution\n", 
	                    sample1_name, " vs ", sample2_name, 
	                    "\nSmoothed vs Unsmoothed Analysis"),
	       xlab = "P-value",
	       ylab = "Frequency",
	       col = alpha("blue", 0.5),
	       border = "blue",
	       ylim = c(0, max_count * 1.1),
	       cex.main = 1.3,
	       cex.lab = 1.2)
	  
	  # Overlay unsmoothed (foreground)
	  hist(dmlTest_unsmoothed$pval, breaks = 50,
	       col = alpha("red", 0.5),
	       border = "red",
	       add = TRUE)
	  
	  # Add significance threshold line
	  abline(v = p_threshold, col = "black", lty = 2, lwd = 2)
	  
	  # Add legend
	  legend("topright", 
	         legend = c("Smoothed", "Unsmoothed", paste("p =", p_threshold)),
	         col = c("blue", "red", "black"),
	         lty = c(1, 1, 2),
	         lwd = c(8, 8, 2),
	         cex = 1.1)
	}


	# ===== 2. DML Methylation Difference Distribution - Overlaid =====
	if (verbose) cat("Creating overlaid DML methylation difference distribution...\n")

	# Create file paths with sample names
	grid_plot_path_diff <- file.path(output_dir, paste0(sample1_name, "_vs_", sample2_name, "_DML_methylation_diff_overlaid_grid.png"))
	wide_plot_path_diff <- file.path(output_dir, paste0(sample1_name, "_vs_", sample2_name, "_DML_methylation_diff_overlaid_wide.png"))

	# Create histogram data for methylation differences
	hist_smooth_diff <- hist(dmlTest_smoothed$diff, breaks = 50, plot = FALSE)
	hist_unsmooth_diff <- hist(dmlTest_unsmoothed$diff, breaks = 50, plot = FALSE)

	# Find common y-axis limits
	max_count_diff <- max(c(hist_smooth_diff$counts, hist_unsmooth_diff$counts))

	# Function to create the plot (for reuse)
	create_meth_diff_plot <- function() {
	  # White background
	  par(bg = "white", mar = c(5, 5, 4, 2))
	  
	  # Plot smoothed first
	  hist(dmlTest_smoothed$diff, breaks = 50,
	       main = paste0("DML Methylation Difference Distribution\n", 
	                     sample1_name, " vs ", sample2_name, 
	                     "\nSmoothed vs Unsmoothed Analysis"),
	       xlab = paste0("Methylation Difference (", sample1_name, " - ", sample2_name, ")"),
	       ylab = "Frequency",
	       col = alpha("blue", 0.5),
	       border = "blue",
	       ylim = c(0, max_count_diff * 1.1),
	       cex.main = 1.3,
	       cex.lab = 1.2)
	  
	  # Overlay unsmoothed
	  hist(dmlTest_unsmoothed$diff, breaks = 50,
	       col = alpha("red", 0.5),
	       border = "red",
	       add = TRUE)
	  
	  # Add threshold lines
	  abline(v = 0, col = "black", lty = 1, lwd = 1)
	  abline(v = c(-delta_threshold, delta_threshold), col = "darkgreen", lty = 2, lwd = 2)
	  
	  # Add legend
	  legend("topright", 
	         legend = c("Smoothed", "Unsmoothed", "No difference", 
	                    paste("±", delta_threshold*100, "% threshold")),
	         col = c("blue", "red", "black", "darkgreen"),
	         lty = c(1, 1, 1, 2),
	         lwd = c(8, 8, 1, 2),
	         cex = 1.1)
	}

	# Create grid version (square-ish)
	png(grid_plot_path_diff, width = 800, height = 800, res = 120)
	create_meth_diff_plot()
	dev.off()

	# Create wide version (for presentations)
	png(wide_plot_path_diff, width = 1200, height = 700, res = 120)
	create_meth_diff_plot()
	dev.off()

	# Save the plot in the list
	plots$meth_diff <- create_meth_diff_plot

	if (verbose) cat("Methylation difference distribution plots saved to:\n", 
	                 grid_plot_path_diff, "\n", wide_plot_path_diff, "\n")



	# ===== 3. DML Volcano Plot - Overlaid =====
	if (verbose) cat("Creating overlaid volcano plot...\n")

	# Create file paths with sample names
	grid_plot_path_volcano <- file.path(output_dir, paste0(sample1_name, "_vs_", sample2_name, "_volcano_plot_overlaid_grid.png"))
	wide_plot_path_volcano <- file.path(output_dir, paste0(sample1_name, "_vs_", sample2_name, "_volcano_plot_overlaid_wide.png"))

	# Create function to generate volcano plot
	create_volcano_plot <- function() {
	  # White background
	  par(bg = "white", mar = c(5, 5, 4, 2))
	  
	  # Clean data for plotting (remove NA, infinite, and zero p-values)
	  smoothed_clean <- data.frame(
	    pval = dmlTest_smoothed$pval,
	    diff = dmlTest_smoothed$diff
	  )
	  smoothed_clean <- smoothed_clean[!is.na(smoothed_clean$pval) & 
	                                  smoothed_clean$pval > 0 & 
	                                  is.finite(smoothed_clean$pval) &
	                                  !is.na(smoothed_clean$diff) &
	                                  is.finite(smoothed_clean$diff), ]
	  
	  unsmoothed_clean <- data.frame(
	    pval = dmlTest_unsmoothed$pval,
	    diff = dmlTest_unsmoothed$diff
	  )
	  unsmoothed_clean <- unsmoothed_clean[!is.na(unsmoothed_clean$pval) & 
	                                      unsmoothed_clean$pval > 0 & 
	                                      is.finite(unsmoothed_clean$pval) &
	                                      !is.na(unsmoothed_clean$diff) &
	                                      is.finite(unsmoothed_clean$diff), ]
	  
	  # Calculate plot limits with cleaned data
	  if(nrow(smoothed_clean) > 0 && nrow(unsmoothed_clean) > 0) {
	    x_max <- max(c(-log10(smoothed_clean$pval), -log10(unsmoothed_clean$pval)), na.rm = TRUE)
	    y_max <- max(c(abs(smoothed_clean$diff), abs(unsmoothed_clean$diff)), na.rm = TRUE)
	    y_min <- min(c(smoothed_clean$diff, unsmoothed_clean$diff), na.rm = TRUE)
	    
	    # Handle edge cases
	    if(!is.finite(x_max) || x_max <= 0) x_max <- 10
	    if(!is.finite(y_max) || y_max <= 0) y_max <- 1
	    if(!is.finite(y_min)) y_min <- -1
	    
	    # Plot background points (all data)
	    plot(-log10(smoothed_clean$pval), smoothed_clean$diff,
	         xlim = c(0, x_max * 1.05),
	         ylim = c(y_min * 1.1, y_max * 1.1),
	         xlab = "-log10(P-value)",
	         ylab = paste0("Methylation Difference (", sample1_name, " - ", sample2_name, ")"),
	         main = paste0("Volcano Plot: DML Results\n", 
	                      sample1_name, " vs ", sample2_name, 
	                      "\nSmoothed vs Unsmoothed Analysis"),
	         pch = 16, cex = 0.3, col = alpha("blue", 0.3),
	         cex.main = 1.3,
	         cex.lab = 1.2)
	    
	    # Overlay unsmoothed background points
	    points(-log10(unsmoothed_clean$pval), unsmoothed_clean$diff,
	           pch = 16, cex = 0.3, col = alpha("red", 0.3))
	    
	    # Highlight significant points if they exist
	    if(nrow(dmls_smoothed) > 0) {
	      # Create position identifier for matching
	      smoothed_positions <- paste(dmlTest_smoothed$chr, dmlTest_smoothed$pos, sep = "_")
	      dmls_smoothed_positions <- paste(dmls_smoothed$chr, dmls_smoothed$pos, sep = "_")
	      
	      # Find indices of significant DMLs in the full dataset
	      sig_smooth_idx <- which(smoothed_positions %in% dmls_smoothed_positions &
	                             !is.na(dmlTest_smoothed$pval) & 
	                             dmlTest_smoothed$pval > 0 &
	                             is.finite(dmlTest_smoothed$pval) &
	                             !is.na(dmlTest_smoothed$diff) &
	                             is.finite(dmlTest_smoothed$diff))
	      
	      if(length(sig_smooth_idx) > 0) {
	        points(-log10(dmlTest_smoothed$pval[sig_smooth_idx]), 
	               dmlTest_smoothed$diff[sig_smooth_idx],
	               pch = 16, cex = 0.8, col = "blue")
	      }
	    }
	    
	    if(nrow(dmls_unsmoothed) > 0) {
	      # Create position identifier for matching
	      unsmoothed_positions <- paste(dmlTest_unsmoothed$chr, dmlTest_unsmoothed$pos, sep = "_")
	      dmls_unsmoothed_positions <- paste(dmls_unsmoothed$chr, dmls_unsmoothed$pos, sep = "_")
	      
	      # Find indices of significant DMLs in the full dataset
	      sig_unsmooth_idx <- which(unsmoothed_positions %in% dmls_unsmoothed_positions &
	                               !is.na(dmlTest_unsmoothed$pval) & 
	                               dmlTest_unsmoothed$pval > 0 &
	                               is.finite(dmlTest_unsmoothed$pval) &
	                               !is.na(dmlTest_unsmoothed$diff) &
	                               is.finite(dmlTest_unsmoothed$diff))
	      
	      if(length(sig_unsmooth_idx) > 0) {
	        points(-log10(dmlTest_unsmoothed$pval[sig_unsmooth_idx]), 
	               dmlTest_unsmoothed$diff[sig_unsmooth_idx],
	               pch = 17, cex = 0.8, col = "red")  # Different shape for distinction
	      }
	    }
	    
	    # Add threshold lines
	    if(is.finite(-log10(p_threshold))) {
	      abline(v = -log10(p_threshold), col = "darkgreen", lty = 2, lwd = 2)
	    }
	    abline(h = c(-delta_threshold, delta_threshold), col = "darkgreen", lty = 2, lwd = 2)
	    abline(h = 0, col = "black", lty = 1, lwd = 1)
	    
	    # Add legend
	    legend_items <- c("Smoothed (all)", "Unsmoothed (all)")
	    legend_colors <- c(alpha("blue", 0.3), alpha("red", 0.3))
	    legend_pch <- c(16, 16)
	    
	    if(nrow(dmls_smoothed) > 0) {
	      legend_items <- c(legend_items, "Smoothed (sig)")
	      legend_colors <- c(legend_colors, "blue")
	      legend_pch <- c(legend_pch, 16)
	    }
	    
	    if(nrow(dmls_unsmoothed) > 0) {
	      legend_items <- c(legend_items, "Unsmoothed (sig)")
	      legend_colors <- c(legend_colors, "red")
	      legend_pch <- c(legend_pch, 17)
	    }
	    
	    legend("topleft", 
	           legend = legend_items,
	           col = legend_colors,
	           pch = legend_pch,
	           cex = 1.0,
	           pt.cex = 1.2)
	    
	    # Add data summary text
	    text(x = x_max * 0.7, y = y_max * 0.9,
	         labels = paste("Smoothed points:", nrow(smoothed_clean),
	                       "\nUnsmoothed points:", nrow(unsmoothed_clean)),
	         cex = 1.0, adj = 0)
	  } else {
	    # Create an empty plot with a message if no valid data
	    plot(0, 0, type = "n", 
	         xlim = c(0, 1), ylim = c(0, 1),
	         xlab = "", ylab = "",
	         main = "No valid data for volcano plot")
	    text(0.5, 0.5, "No valid data available for volcano plot", cex = 1.5)
	    
	    if (verbose) cat("Warning: No valid data for volcano plot. Created empty plot.\n")
	  }
	}

	# Create grid version (square-ish)
	png(grid_plot_path_volcano, width = 800, height = 800, res = 120)
	create_volcano_plot()
	dev.off()

	# Create wide version (for presentations)
	png(wide_plot_path_volcano, width = 1200, height = 700, res = 120)
	create_volcano_plot()
	dev.off()

	# Save the plot in the list
	plots$volcano <- create_volcano_plot

	if (verbose) cat("Volcano plot saved to:\n", 
	                grid_plot_path_volcano, "\n", wide_plot_path_volcano, "\n")

	# ===== 6. DML Count per Chromosome - Bar Chart =====
	if (verbose) cat("Creating DML count per chromosome bar chart...\n")

	# Create file paths with sample names
	grid_plot_path_chr <- file.path(output_dir, paste0(sample1_name, "_vs_", sample2_name, "_DML_count_per_chromosome_grid.png"))
	wide_plot_path_chr <- file.path(output_dir, paste0(sample1_name, "_vs_", sample2_name, "_DML_count_per_chromosome_wide.png"))

	# Function to create the DML per chromosome plot
	create_chr_count_plot <- function() {
	  # White background
	  par(bg = "white", mar = c(5, 5, 4, 2))
	  
	  # Check if both DML results have data
	  if(nrow(dmls_smoothed) > 0 && nrow(dmls_unsmoothed) > 0) {
	    # Count DMLs per chromosome for both methods
	    smooth_chr_counts <- table(dmls_smoothed$chr)
	    unsmooth_chr_counts <- table(dmls_unsmoothed$chr)
	    
	    # Get all unique chromosomes
	    all_chrs <- sort(unique(c(names(smooth_chr_counts), names(unsmooth_chr_counts))))
	    
	    # Create data frame with counts for both methods
	    chr_counts <- data.frame(
	      chr = all_chrs,
	      smoothed = 0,
	      unsmoothed = 0
	    )
	    
	    # Fill in counts
	    for(i in 1:length(all_chrs)) {
	      chr <- all_chrs[i]
	      if(chr %in% names(smooth_chr_counts)) {
	        chr_counts$smoothed[i] <- smooth_chr_counts[chr]
	      }
	      if(chr %in% names(unsmooth_chr_counts)) {
	        chr_counts$unsmoothed[i] <- unsmooth_chr_counts[chr]
	      }
	    }
	    
	    # Order chromosomes naturally (if possible)
	    if(all(grepl("^chr[0-9]+$", chr_counts$chr)) || all(grepl("^[0-9]+$", chr_counts$chr))) {
	      # Extract numbers
	      chr_nums <- as.numeric(gsub("chr", "", chr_counts$chr))
	      chr_counts <- chr_counts[order(chr_nums), ]
	    }
	    
	    # Set up the barplot
	    max_count <- max(c(chr_counts$smoothed, chr_counts$unsmoothed))
	    
	    # Calculate bar positions
	    bar_count <- nrow(chr_counts)
	    bar_width <- 0.35
	    positions <- barplot(chr_counts$smoothed, plot = FALSE)
	    x_smooth <- positions - bar_width/2
	    x_unsmooth <- positions + bar_width/2
	    
	    # Create the plot
	    barplot(height = rep(0, bar_count), 
	            names.arg = chr_counts$chr, 
	            ylim = c(0, max_count * 1.1),
	            main = paste0("DML Count per Chromosome\n", 
	                         sample1_name, " vs ", sample2_name),
	            xlab = "Chromosome",
	            ylab = "DML Count",
	            cex.main = 1.3,
	            cex.lab = 1.2,
	            border = NA)
	    
	    # Add bars for each method
	    for(i in 1:bar_count) {
	      rect(xleft = x_smooth[i] - bar_width/2, 
	           ybottom = 0, 
	           xright = x_smooth[i] + bar_width/2, 
	           ytop = chr_counts$smoothed[i],
	           col = alpha("blue", 0.7), 
	           border = "blue")
	      
	      rect(xleft = x_unsmooth[i] - bar_width/2, 
	           ybottom = 0, 
	           xright = x_unsmooth[i] + bar_width/2, 
	           ytop = chr_counts$unsmoothed[i],
	           col = alpha("red", 0.7), 
	           border = "red")
	    }
	    
	    # Add grid lines
	    grid_lines <- axTicks(2)
	    for(y in grid_lines) {
	      abline(h = y, col = "lightgray", lty = 2)
	    }
	    
	    # Add legend
	    legend("topright", 
	           legend = c("Smoothed", "Unsmoothed"),
	           fill = c(alpha("blue", 0.7), alpha("red", 0.7)),
	           border = c("blue", "red"),
	           cex = 1.1)
	    
	    # Add total count
	    text(x = positions[length(positions)] * 0.9, 
	         y = max_count * 0.9,
	         labels = paste("Total DMLs:\nSmoothed:", sum(chr_counts$smoothed),
	                       "\nUnsmoothed:", sum(chr_counts$unsmoothed)),
	         adj = 1,
	         cex = 1.1)
	    
	  } else {
	    # Create an empty plot with a message if no DMLs
	    plot(0, 0, type = "n", 
	         xlim = c(0, 1), ylim = c(0, 1),
	         xlab = "", ylab = "",
	         main = "No DML data available")
	    text(0.5, 0.5, "No DML data available for chromosome distribution", cex = 1.5)
	    
	    if (verbose) cat("Warning: No DML data for chromosome distribution. Created empty plot.\n")
	  }
	}

	# Create grid version (square-ish)
	png(grid_plot_path_chr, width = 800, height = 800, res = 120)
	create_chr_count_plot()
	dev.off()

	# Create wide version (for presentations)
	png(wide_plot_path_chr, width = 1200, height = 700, res = 120)
	create_chr_count_plot()
	dev.off()

	# Save the plot in the list
	plots$chr_count <- create_chr_count_plot
	plot_paths$chr_count <- list(grid = grid_plot_path_chr, wide = wide_plot_path_chr)

	if (verbose) cat("DML per chromosome plots saved to:\n", 
	                grid_plot_path_chr, "\n", wide_plot_path_chr, "\n")
	# ===== 4. DMR Length Distribution - Overlaid =====
	if (verbose) cat("Creating overlaid DMR length distribution...\n")

	# Create file paths with sample names
	grid_plot_path_dmr_len <- file.path(output_dir, paste0(sample1_name, "_vs_", sample2_name, "_DMR_length_overlaid_grid.png"))
	wide_plot_path_dmr_len <- file.path(output_dir, paste0(sample1_name, "_vs_", sample2_name, "_DMR_length_overlaid_wide.png"))

	# Function to create the DMR length plot
	create_dmr_length_plot <- function() {
	  # White background
	  par(bg = "white", mar = c(5, 5, 4, 2))
	  
	  # Check if both DMR results have data
	  if(nrow(dmrs_smoothed) > 0 && nrow(dmrs_unsmoothed) > 0) {
	    # Find the range of lengths across both datasets
	    min_len <- min(c(dmrs_smoothed$length, dmrs_unsmoothed$length), na.rm = TRUE)
	    max_len <- max(c(dmrs_smoothed$length, dmrs_unsmoothed$length), na.rm = TRUE)
	    
	    # Create common breaks for both histograms to ensure equal-sized bins
	    common_breaks <- seq(from = min_len, to = max_len, length.out = 31)  # 30 bins
	    
	    # Create histogram data using the common breaks
	    hist_smooth_dmr <- hist(dmrs_smoothed$length, breaks = common_breaks, plot = FALSE)
	    hist_unsmooth_dmr <- hist(dmrs_unsmoothed$length, breaks = common_breaks, plot = FALSE)
	    
	    # Find common y-axis limits
	    max_count_dmr <- max(c(hist_smooth_dmr$counts, hist_unsmooth_dmr$counts))
	    
	    # Plot smoothed DMR lengths
	    hist(dmrs_smoothed$length, breaks = common_breaks,
	         main = paste0("DMR Length Distribution\n", 
	                      sample1_name, " vs ", sample2_name, 
	                      "\nSmoothed vs Unsmoothed Analysis"),
	         xlab = "DMR Length (bp)",
	         ylab = "Frequency",
	         col = alpha("blue", 0.5),
	         border = "blue",
	         ylim = c(0, max_count_dmr * 1.1),
	         cex.main = 1.3,
	         cex.lab = 1.2)
	    
	    # Overlay unsmoothed DMR lengths using the same breaks
	    hist(dmrs_unsmoothed$length, breaks = common_breaks,
	         col = alpha("red", 0.5),
	         border = "red",
	         add = TRUE)
	    
	    # Add statistics text
	    text(x = max_len * 0.7,
	         y = max_count_dmr * 0.9,
	         labels = paste("Smoothed DMRs:", nrow(dmrs_smoothed), 
	                       "\nMean length:", round(mean(dmrs_smoothed$length), 1), "bp",
	                       "\n\nUnsmoothed DMRs:", nrow(dmrs_unsmoothed),
	                       "\nMean length:", round(mean(dmrs_unsmoothed$length), 1), "bp"),
	         cex = 1.1,
	         adj = 0)
	    
	    # Add legend
	    legend("topright", 
	           legend = c("Smoothed DMRs", "Unsmoothed DMRs"),
	           col = c("blue", "red"),
	           lwd = 8,
	           cex = 1.1)
	  } else {
	    # Create an empty plot with a message if no DMRs
	    plot(0, 0, type = "n", 
	         xlim = c(0, 1), ylim = c(0, 1),
	         xlab = "", ylab = "",
	         main = "No DMR data available")
	    text(0.5, 0.5, "No DMR data available for length distribution", cex = 1.5)
	    
	    if (verbose) cat("Warning: No DMR data for length distribution. Created empty plot.\n")
	  }
	}

	# Create grid version (square-ish)
	png(grid_plot_path_dmr_len, width = 800, height = 800, res = 120)
	create_dmr_length_plot()
	dev.off()

	# Create wide version (for presentations)
	png(wide_plot_path_dmr_len, width = 1200, height = 700, res = 120)
	create_dmr_length_plot()
	dev.off()

	# Save the plot in the list
	plots$dmr_length <- create_dmr_length_plot

	if (verbose) cat("DMR length distribution plots saved to:\n", 
	                grid_plot_path_dmr_len, "\n", wide_plot_path_dmr_len, "\n")


	# ===== 5. DMR Methylation Difference Distribution - Overlaid =====
	if (verbose) cat("Creating overlaid DMR methylation difference distribution...\n")

	# Create file paths with sample names
	grid_plot_path_dmr_diff <- file.path(output_dir, paste0(sample1_name, "_vs_", sample2_name, "_DMR_methylation_diff_overlaid_grid.png"))
	wide_plot_path_dmr_diff <- file.path(output_dir, paste0(sample1_name, "_vs_", sample2_name, "_DMR_methylation_diff_overlaid_wide.png"))

	# Function to create the DMR methylation difference plot
	create_dmr_diff_plot <- function() {
	  # White background
	  par(bg = "white", mar = c(5, 5, 4, 2))
	  
	  # Check if both DMR results have data
	  if(nrow(dmrs_smoothed) > 0 && nrow(dmrs_unsmoothed) > 0) {
	    # Find the range of differences across both datasets
	    min_diff <- min(c(dmrs_smoothed$diff.Methy, dmrs_unsmoothed$diff.Methy), na.rm = TRUE)
	    max_diff <- max(c(dmrs_smoothed$diff.Methy, dmrs_unsmoothed$diff.Methy), na.rm = TRUE)
	    
	    # Create common breaks for both histograms to ensure equal-sized bins
	    common_breaks <- seq(from = min_diff, to = max_diff, length.out = 21)  # 20 bins
	    
	    # Create histogram data using the common breaks
	    hist_smooth_dmr_diff <- hist(dmrs_smoothed$diff.Methy, breaks = common_breaks, plot = FALSE)
	    hist_unsmooth_dmr_diff <- hist(dmrs_unsmoothed$diff.Methy, breaks = common_breaks, plot = FALSE)
	    
	    # Find common y-axis limits
	    max_count_dmr_diff <- max(c(hist_smooth_dmr_diff$counts, hist_unsmooth_dmr_diff$counts))
	    
	    # Plot smoothed DMR methylation differences
	    hist(dmrs_smoothed$diff.Methy, breaks = common_breaks,
	         main = paste0("DMR Methylation Difference Distribution\n", 
	                      sample1_name, " vs ", sample2_name, 
	                      "\nSmoothed vs Unsmoothed Analysis"),
	         xlab = paste0("DMR Methylation Difference (", sample1_name, " - ", sample2_name, ")"),
	         ylab = "Frequency",
	         col = alpha("blue", 0.5),
	         border = "blue",
	         ylim = c(0, max_count_dmr_diff * 1.1),
	         cex.main = 1.3,
	         cex.lab = 1.2)
	    
	    # Overlay unsmoothed DMR methylation differences using the same breaks
	    hist(dmrs_unsmoothed$diff.Methy, breaks = common_breaks,
	         col = alpha("red", 0.5),
	         border = "red",
	         add = TRUE)
	    
	    # Add reference line at zero
	    abline(v = 0, col = "black", lty = 1, lwd = 1)
	    
	    # Count hyper/hypo DMRs
	    hyper_smooth <- sum(dmrs_smoothed$diff.Methy > 0)
	    hypo_smooth <- sum(dmrs_smoothed$diff.Methy < 0)
	    hyper_unsmooth <- sum(dmrs_unsmoothed$diff.Methy > 0)
	    hypo_unsmooth <- sum(dmrs_unsmoothed$diff.Methy < 0)
	    
	    # Add statistics text
	    text(x = min_diff * 0.8,
	         y = max_count_dmr_diff * 0.8,
	         labels = paste("Smoothed:", hyper_smooth, "hyper,", hypo_smooth, "hypo",
	                       "\nUnsmoothed:", hyper_unsmooth, "hyper,", hypo_unsmooth, "hypo"),
	         cex = 1.1,
	         adj = 0)
	    
	    # Add legend
	    legend("topright", 
	           legend = c("Smoothed DMRs", "Unsmoothed DMRs", "No difference"),
	           col = c("blue", "red", "black"),
	           lty = c(1, 1, 1),
	           lwd = c(8, 8, 1),
	           cex = 1.1)
	  } else {
	    # Create an empty plot with a message if no DMRs
	    plot(0, 0, type = "n", 
	         xlim = c(0, 1), ylim = c(0, 1),
	         xlab = "", ylab = "",
	         main = "No DMR data available")
	    text(0.5, 0.5, "No DMR data available for methylation difference distribution", cex = 1.5)
	    
	    if (verbose) cat("Warning: No DMR data for methylation difference distribution. Created empty plot.\n")
	  }
	}

	# Create grid version (square-ish)
	png(grid_plot_path_dmr_diff, width = 800, height = 800, res = 120)
	create_dmr_diff_plot()
	dev.off()

	# Create wide version (for presentations)
	png(wide_plot_path_dmr_diff, width = 1200, height = 700, res = 120)
	create_dmr_diff_plot()
	dev.off()

	# Save the plot in the list
	plots$dmr_diff <- create_dmr_diff_plot

	if (verbose) cat("DMR methylation difference distribution plots saved to:\n", 
	                grid_plot_path_dmr_diff, "\n", wide_plot_path_dmr_diff, "\n")

	# ===== Create Combined Grid of All Plots =====
	if (verbose) cat("Creating combined grid of all plots...\n")

	# Load required libraries for grid layout
	required_packages <- c("grid", "gridExtra", "png")
	for (pkg in required_packages) {
	  if (!requireNamespace(pkg, quietly = TRUE)) {
	    install.packages(pkg)
	  }
	  library(pkg, character.only = TRUE)
	}

	# Create file path for the grid plot - landscape orientation
	grid_plot_path_combined <- file.path(output_dir, paste0(sample1_name, "_vs_", sample2_name, "_methylation_analysis_grid.png"))

	# List of plot file paths (grid versions)
	grid_plots <- list(
	  grid_plot_path,           # p-value dist
	  grid_plot_path_diff,      # methylation diff
	  grid_plot_path_volcano,   # volcano
 	  grid_plot_path_chr        # NEW: dml per chromosome
	  grid_plot_path_dmr_len,   # dmr length
	  grid_plot_path_dmr_diff,  # dmr methylation diff
	)

	# Read all plot PNGs into a list of raster grobs
	plot_grobs <- list()
	for (i in 1:length(grid_plots)) {
	  if (file.exists(grid_plots[[i]])) {
	    img <- png::readPNG(grid_plots[[i]])
	    plot_grobs[[i]] <- grid::rasterGrob(img)
	  }
	}

	# Create a 2x3 grid (landscape orientation)
	png(grid_plot_path_combined, width = 2400, height = 1600, res = 150, bg = "white")

	grid::grid.newpage()
	grid::pushViewport(grid::viewport(layout = grid::grid.layout(2, 3)))

	# First row
	for (i in 1:3) {
	  if (i <= length(plot_grobs)) {
	    grid::pushViewport(grid::viewport(layout.pos.row = 1, layout.pos.col = i))
	    grid::grid.draw(plot_grobs[[i]])
	    grid::popViewport()
	  }
	}

	# Second row
	for (i in 4:6) {
	  if (i <= length(plot_grobs)) {
	    grid::pushViewport(grid::viewport(layout.pos.row = 2, layout.pos.col = i-3))
	    grid::grid.draw(plot_grobs[[i]])
	    grid::popViewport()
	  }
	}

	dev.off()

	if (verbose) cat("Combined plot grid saved to:", grid_plot_path_combined, "\n")

	# Add grid plot path to the returned list
	plot_paths$grid_combined <- grid_plot_path_combined


	  # Return plot functions and file paths
  return(list(
    plots = plots,
    paths = plot_paths
  ))

}