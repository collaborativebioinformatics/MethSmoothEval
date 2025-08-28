# runner_updated.R
# Run the complete DMR analysis pipeline with smoothing comparison

cat("=============================================================================\n")
cat("DMR ANALYSIS PIPELINE - SMOOTHING COMPARISON\n")
cat("=============================================================================\n\n")

# Load data and configuration
cat("Loading data and configuration...\n")
source("load_data.R")
source("samples_config.R")

# =============================================================================
# RUN ANALYSIS WITHOUT SMOOTHING
# =============================================================================

cat("\n=== RUNNING ANALYSIS WITHOUT SMOOTHING ===\n")

# Temporarily modify the analysis script to disable smoothing
analysis_script <- readLines("dmr_analysis_fixed.R")
analysis_script <- gsub("USE_SMOOTHING <- TRUE", "USE_SMOOTHING <- FALSE", analysis_script)
writeLines(analysis_script, "temp_analysis_raw.R")

# Run the raw analysis
source("temp_analysis_raw.R")

# =============================================================================
# RUN ANALYSIS WITH SMOOTHING
# =============================================================================

cat("\n=== RUNNING ANALYSIS WITH SMOOTHING ===\n")

# Reset and run with smoothing enabled
analysis_script <- readLines("dmr_analysis_fixed.R")
analysis_script <- gsub("USE_SMOOTHING <- TRUE", "USE_SMOOTHING <- FALSE", analysis_script)
analysis_script <- gsub("USE_SMOOTHING <- FALSE", "USE_SMOOTHING <- TRUE", analysis_script)
writeLines(analysis_script, "temp_analysis_smooth.R")

# Run the smoothed analysis
source("temp_analysis_smooth.R")

# =============================================================================
# COMPARE RESULTS
# =============================================================================

cat("\n=== COMPARING RAW vs SMOOTHED RESULTS ===\n")
source("compare_longread_results.R")

# =============================================================================
# CLEANUP
# =============================================================================

# Remove temporary files
if(file.exists("temp_analysis_raw.R")) file.remove("temp_analysis_raw.R")
if(file.exists("temp_analysis_smooth.R")) file.remove("temp_analysis_smooth.R")

cat("\n=============================================================================\n")
cat("PIPELINE COMPLETE\n")
cat("=============================================================================\n")
cat("Check the following outputs:\n")
cat("- Individual comparison files (*_dmls_raw.txt, *_dmrs_smoothed.txt, etc.)\n")
cat("- comparison_summary.txt (overall summary table)\n")
cat("- comparison_plots/ directory (all visualization files)\n")
cat("- plots_raw/ and plots_smoothed/ directories (analysis-specific plots)\n")
cat("=============================================================================\n")