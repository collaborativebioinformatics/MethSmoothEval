# compare_longread_results.R
# Compare smoothing effects on long-read DMR detection

cat("=============================================================================\n")
cat("LONG-READ SMOOTHING EFFECTS COMPARISON\n")
cat("=============================================================================\n\n")

# Define analysis types
analysis_types <- c(
  "longread_ONT_vs_PacBio",
  "longread_tissue_lung", 
  "longread_tissue_liver",
  "longread_tissue_colon",
  "HG002_ONT_vs_EMSeq",
  "HG002_PacBio_vs_EMSeq"
)

# Initialize storage
comparison_results <- list()
summary_df <- data.frame(
  Analysis = character(),
  Raw_DMLs = numeric(),
  Smoothed_DMLs = numeric(),
  Raw_DMRs = numeric(),
  Smoothed_DMRs = numeric(),
  DML_Fold = numeric(),
  DMR_Fold = numeric(),
  stringsAsFactors = FALSE
)

# =============================================================================
# LOAD AND COMPARE RESULTS
# =============================================================================

for(analysis in analysis_types) {
  cat(paste("Processing", analysis, "..."))
  
  # File paths
  dml_raw_file <- paste0(analysis, "_dmls_raw.txt")
  dml_smooth_file <- paste0(analysis, "_dmls_smoothed.txt")
  dmr_raw_file <- paste0(analysis, "_dmrs_raw.txt")
  dmr_smooth_file <- paste0(analysis, "_dmrs_smoothed.txt")
  
  # Load files
  dml_raw <- if(file.exists(dml_raw_file)) read.table(dml_raw_file, header=TRUE, sep="\t") else data.frame()
  dml_smooth <- if(file.exists(dml_smooth_file)) read.table(dml_smooth_file, header=TRUE, sep="\t") else data.frame()
  dmr_raw <- if(file.exists(dmr_raw_file)) read.table(dmr_raw_file, header=TRUE, sep="\t") else data.frame()
  dmr_smooth <- if(file.exists(dmr_smooth_file)) read.table(dmr_smooth_file, header=TRUE, sep="\t") else data.frame()
  
  # Calculate fold changes
  dml_fold <- ifelse(nrow(dml_raw) > 0, nrow(dml_smooth) / nrow(dml_raw), 
                     ifelse(nrow(dml_smooth) > 0, Inf, 1))
  dmr_fold <- ifelse(nrow(dmr_raw) > 0, nrow(dmr_smooth) / nrow(dmr_raw), 
                     ifelse(nrow(dmr_smooth) > 0, Inf, 1))
  
  # Store results
  comparison_results[[analysis]] <- list(
    dml_raw = dml_raw, dml_smooth = dml_smooth,
    dmr_raw = dmr_raw, dmr_smooth = dmr_smooth
  )
  
  # Add to summary
  summary_df <- rbind(summary_df, data.frame(
    Analysis = analysis,
    Raw_DMLs = nrow(dml_raw),
    Smoothed_DMLs = nrow(dml_smooth),
    Raw_DMRs = nrow(dmr_raw),
    Smoothed_DMRs = nrow(dmr_smooth),
    DML_Fold = round(dml_fold, 2),
    DMR_Fold = round(dmr_fold, 2)
  ))
  
  cat(paste(" DMLs:", nrow(dml_raw), "→", nrow(dml_smooth), 
            "| DMRs:", nrow(dmr_raw), "→", nrow(dmr_smooth), "\n"))
}

# =============================================================================
# SMOOTHING EFFECTS SUMMARY
# =============================================================================

cat("\nSMOOTHING EFFECTS SUMMARY:\n")
print(summary_df)

# Platform-specific effects
longread_platform <- summary_df[summary_df$Analysis == "longread_ONT_vs_PacBio", ]
hg002_effects <- summary_df[grepl("HG002", summary_df$Analysis), ]

cat("\n=== KEY FINDINGS ===\n")
if(nrow(longread_platform) > 0) {
  cat("ONT vs PacBio smoothing effect:\n")
  cat(paste("- DML fold change:", longread_platform$DML_Fold, "\n"))
  cat(paste("- DMR fold change:", longread_platform$DMR_Fold, "\n"))
}

if(nrow(hg002_effects) > 0) {
  cat("HG002 long-read vs short-read:\n")
  for(i in 1:nrow(hg002_effects)) {
    row <- hg002_effects[i, ]
    cat(paste("- ", row$Analysis, ": DMR fold =", row$DMR_Fold, "\n"))
  }
}

# =============================================================================
# OVERLAP ANALYSIS
# =============================================================================

cat("\n=== DMR OVERLAP ANALYSIS ===\n")
for(analysis in names(comparison_results)) {
  result <- comparison_results[[analysis]]
  dmr_raw <- result$dmr_raw
  dmr_smooth <- result$dmr_smooth
  
  if(nrow(dmr_raw) > 0 && nrow(dmr_smooth) > 0) {
    # Simple overlap check
    overlaps <- 0
    for(i in 1:nrow(dmr_raw)) {
      overlap_found <- any(
        dmr_smooth$chr == dmr_raw$chr[i] &
        dmr_smooth$start <= dmr_raw$end[i] &
        dmr_smooth$end >= dmr_raw$start[i]
      )
      if(overlap_found) overlaps <- overlaps + 1
    }
    
    overlap_pct <- round(100 * overlaps / nrow(dmr_raw), 1)
    cat(paste(analysis, "- Overlap:", overlap_pct, "%\n"))
  }
}

# =============================================================================
# VISUALIZATION
# =============================================================================

if(!dir.exists("smoothing_plots")) dir.create("smoothing_plots")

# Main comparison plot
png("smoothing_plots/smoothing_effects_summary.png", width=1200, height=800)
par(mfrow=c(2,2))

# DML effects
barplot(summary_df$DML_Fold, names.arg=summary_df$Analysis,
        main="DML Smoothing Effects", col="lightblue", las=2, cex.names=0.7)
abline(h=1, lty=2, col="red")

# DMR effects
barplot(summary_df$DMR_Fold, names.arg=summary_df$Analysis,
        main="DMR Smoothing Effects", col="lightgreen", las=2, cex.names=0.7)
abline(h=1, lty=2, col="red")

# Raw vs smoothed counts
barplot(t(as.matrix(summary_df[,c("Raw_DMRs", "Smoothed_DMRs")])), 
        beside=TRUE, names.arg=summary_df$Analysis,
        main="DMR Counts: Raw vs Smoothed", col=c("orange", "purple"),
        las=2, legend.text=c("Raw", "Smoothed"), cex.names=0.7)

# Focus on platform comparison
platform_idx <- which(summary_df$Analysis == "longread_ONT_vs_PacBio")
if(length(platform_idx) > 0) {
  barplot(c(summary_df$DML_Fold[platform_idx], summary_df$DMR_Fold[platform_idx]),
          names.arg=c("DML Fold", "DMR Fold"),
          main="ONT vs PacBio Smoothing Effect", col="cyan")
  abline(h=1, lty=2, col="red")
}

dev.off()

# =============================================================================
# SAVE RESULTS
# =============================================================================

write.table(summary_df, "smoothing_effects_summary.txt", sep="\t", row.names=FALSE, quote=FALSE)

cat("\n=============================================================================\n")
cat("SMOOTHING COMPARISON COMPLETE\n")
cat("=============================================================================\n")
cat("Generated:\n")
cat("- smoothing_effects_summary.txt\n")
cat("- smoothing_plots/smoothing_effects_summary.png\n")
cat("=============================================================================\n")