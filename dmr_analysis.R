# dmr_analysis.R  
# DMR Analysis - Focus on smoothing effects in long-read data

# Load additional libraries for BigWig export
library(rtracklayer)
library(GenomicRanges)

# =============================================================================
# ANALYSIS PARAMETERS
# =============================================================================

USE_SMOOTHING <- TRUE  # Toggle this for comparison runs
GENERATE_BIGWIGS <- FALSE  # Skip for focused analysis
P_THRESHOLD <- 0.001
MIN_LENGTH <- 50
MIN_CG <- 3
DIS_MERGE <- 100
N_CORES <- 2

cat("=============================================================================\n")
cat("SMOOTHING EFFECTS ANALYSIS\n")
cat("=============================================================================\n")
cat(paste("Smoothing:", ifelse(USE_SMOOTHING, "ENABLED", "DISABLED"), "\n"))
cat(paste("Focus: Long-read platform smoothing artifacts\n"))
cat("=============================================================================\n\n")

# =============================================================================
# CHECK DATA
# =============================================================================

if(!exists("BSobj_raw")) stop("BSobj_raw not found. Run source('load_data.R') first.")
if(!exists("sample_info")) stop("sample_info not found. Run source('samples_config.R') first.")

# Extract factors
tissue <- factor(sample_info$tissue)
individual <- factor(sample_info$individual)
platform <- factor(sample_info$platform)
tech_category <- factor(sample_info$tech_category)

cat("SAMPLE OVERVIEW:\n")
print(sample_info[,c("sample", "individual", "tissue", "platform", "tech_category")])
cat("\n")

# =============================================================================
# DATA PREPROCESSING
# =============================================================================

BSobj <- BSobj_raw

if(USE_SMOOTHING) {
  cat("Applying smoothing...\n")
  start_time <- Sys.time()
  BSobj <- BSmooth(BSobj, BPPARAM = MulticoreParam(workers = N_CORES))
  end_time <- Sys.time()
  cat(paste("Completed in", round(difftime(end_time, start_time, units="mins"), 2), "minutes\n\n"))
} else {
  cat("Using raw data\n\n")
}

# =============================================================================
# FOCUSED ANALYSIS: LONG-READ SAMPLES ONLY
# =============================================================================

# Identify long-read samples (ONT + PacBio)
longread_idx <- which(tech_category == "Long")
cat(paste("Long-read samples identified:", length(longread_idx), "\n"))
cat("Long-read samples:\n")
print(sample_info[longread_idx, c("sample", "platform", "individual", "tissue")])

if(length(longread_idx) < 3) {
  stop("Insufficient long-read samples for analysis")
}

# Subset to long-read samples
BSobj_lr <- BSobj[, longread_idx]
lr_sample_info <- sample_info[longread_idx, ]

# Design matrix for long-read samples only
lr_tissue <- factor(lr_sample_info$tissue)
lr_individual <- factor(lr_sample_info$individual)
lr_platform <- factor(lr_sample_info$platform)

design_lr <- model.matrix(~lr_tissue + lr_individual + lr_platform)

cat("\nLong-read design matrix:\n")
print(design_lr)
cat("Testing coefficients:\n")
print(colnames(design_lr))

# =============================================================================
# DMR ANALYSIS ON LONG-READ DATA
# =============================================================================

cat("\nRunning DML analysis on long-read samples...\n")
dmlTest_lr <- DMLfit.multiFactor(BSobj_lr, design = design_lr)

# Key contrasts for long-read platform effects
lr_contrasts <- list(
  ONT_vs_PacBio = "lr_platformPacBio",
  tissue_lung = "lr_tissuelung", 
  tissue_liver = "lr_tissueliver",
  tissue_colon = "lr_tissuecolon"
)

all_results <- list()

for(contrast_name in names(lr_contrasts)) {
  coef <- lr_contrasts[[contrast_name]]
  
  # Skip if coefficient doesn't exist
  if(!coef %in% colnames(design_lr)) {
    cat(paste("Skipping", contrast_name, "- coefficient not in design\n"))
    next
  }
  
  cat(paste("Testing:", contrast_name, "\n"))
  
  dmls <- DMLtest.multiFactor(dmlTest_lr, coef = coef)
  dmls.sig <- dmls[dmls$pvals < P_THRESHOLD, ]
  
  # Call DMRs
  if(nrow(dmls.sig) >= MIN_CG) {
    dmrs <- callDMR(dmlTest_lr, p.threshold=P_THRESHOLD, minlen=MIN_LENGTH, 
                    minCG=MIN_CG, dis.merge=DIS_MERGE, coef=coef)
  } else {
    dmrs <- data.frame()
  }
  
  all_results[[contrast_name]] <- list(
    dmls = dmls,
    dmls.sig = dmls.sig,
    dmrs = dmrs
  )
  
  # Save results
  smooth_suffix <- ifelse(USE_SMOOTHING, "_smoothed", "_raw")
  write.table(dmls.sig, paste0("longread_", contrast_name, "_dmls", smooth_suffix, ".txt"), 
              sep="\t", row.names=FALSE, quote=FALSE)
  if(nrow(dmrs) > 0) {
    write.table(dmrs, paste0("longread_", contrast_name, "_dmrs", smooth_suffix, ".txt"), 
                sep="\t", row.names=FALSE, quote=FALSE)
  }
  
  cat(paste("  DMLs:", nrow(dmls.sig), "| DMRs:", nrow(dmrs), "\n"))
}

# =============================================================================
# SAVE RESULTS
# =============================================================================

# Save processed BSobj and results for reuse
smooth_suffix <- ifelse(USE_SMOOTHING, "_smoothed", "_raw")
save(BSobj, BSobj_lr, dmlTest_lr, all_results, sample_info, 
     file = paste0("longread_analysis", smooth_suffix, ".RData"))
cat(paste("Saved analysis to longread_analysis", smooth_suffix, ".RData\n\n", sep=""))

# =============================================================================
# HG002 VALIDATION COMPARISON
# =============================================================================

cat("\n=== HG002 VALIDATION ===\n")
hg002_idx <- which(individual == "HG002")

if(length(hg002_idx) >= 3) {
  BSobj_hg002 <- BSobj[, hg002_idx]
  hg002_info <- sample_info[hg002_idx, ]
  
  design_hg002 <- model.matrix(~platform, data = hg002_info)
  dmlTest_hg002 <- DMLfit.multiFactor(BSobj_hg002, design = design_hg002)
  
  # Test long-read vs short-read in HG002
  hg002_contrasts <- c("ONT_vs_EMSeq", "PacBio_vs_EMSeq")
  
  for(hg_contrast in hg002_contrasts) {
    platform_coef <- switch(hg_contrast,
                           "ONT_vs_EMSeq" = "platformONT",
                           "PacBio_vs_EMSeq" = "platformPacBio")
    
    if(platform_coef %in% colnames(design_hg002)) {
      dmls_hg <- DMLtest.multiFactor(dmlTest_hg002, coef = platform_coef)
      dmls_hg.sig <- dmls_hg[dmls_hg$pvals < P_THRESHOLD, ]
      
      if(nrow(dmls_hg.sig) >= MIN_CG) {
        dmrs_hg <- callDMR(dmlTest_hg002, p.threshold=P_THRESHOLD, minlen=MIN_LENGTH,
                          minCG=MIN_CG, dis.merge=DIS_MERGE, coef=platform_coef)
      } else {
        dmrs_hg <- data.frame()
      }
      
      all_results[[hg_contrast]] <- list(dmls.sig = dmls_hg.sig, dmrs = dmrs_hg)
      
      smooth_suffix <- ifelse(USE_SMOOTHING, "_smoothed", "_raw")
      write.table(dmls_hg.sig, paste0("HG002_", hg_contrast, "_dmls", smooth_suffix, ".txt"), 
                  sep="\t", row.names=FALSE, quote=FALSE)
      if(nrow(dmrs_hg) > 0) {
        write.table(dmrs_hg, paste0("HG002_", hg_contrast, "_dmrs", smooth_suffix, ".txt"), 
                    sep="\t", row.names=FALSE, quote=FALSE)
      }
      
      cat(paste(hg_contrast, "- DMLs:", nrow(dmls_hg.sig), "| DMRs:", nrow(dmrs_hg), "\n"))
    }
  }
}

# =============================================================================
# SUMMARY
# =============================================================================

cat("\n=============================================================================\n")
cat("ANALYSIS SUMMARY\n") 
cat("=============================================================================\n")
cat(paste("Smoothing:", ifelse(USE_SMOOTHING, "Applied", "Not applied"), "\n"))
cat(paste("Long-read samples analyzed:", length(longread_idx), "\n"))

summary_df <- data.frame(
  Contrast = names(all_results),
  DMLs = sapply(all_results, function(x) nrow(x$dmls.sig)),
  DMRs = sapply(all_results, function(x) nrow(x$dmrs))
)
print(summary_df)

# Save summary
smooth_suffix <- ifelse(USE_SMOOTHING, "_smoothed", "_raw")
write.table(summary_df, paste0("longread_analysis_summary", smooth_suffix, ".txt"), 
            sep="\t", row.names=FALSE, quote=FALSE)

cat("\nKey outputs:\n")
cat("- longread_*_dmls_*.txt (DML results)\n") 
cat("- longread_*_dmrs_*.txt (DMR results)\n")
cat("- HG002_*_*.txt (validation results)\n")
cat(paste("- longread_analysis_summary", smooth_suffix, ".txt\n"))
cat("=============================================================================\n")