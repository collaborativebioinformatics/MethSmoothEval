

source("load_data.R")

# Test WITHOUT smoothing
USE_SMOOTHING <- FALSE
source("dmr_analysis.R")

# Test WITH smoothing  
USE_SMOOTHING <- TRUE
source("dmr_analysis.R")

# Compare both results
source("compare_results.R")