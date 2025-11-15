#!/usr/bin/env Rscript

# site_level_agreement.R
# Calculates pairwise correlations, concordance metrics, and generates
# scatter plots comparing methylation levels across technologies

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(ggplot2)
  library(cowplot)
  library(corrplot)
  library(RColorBrewer)
})

# ============================================================================
# Command-line argument parsing
# ============================================================================

option_list <- list(
  make_option(c("-i", "--input"), type="character", default=NULL,
              help="Input RDS file from identify_common_sites.R (required)",
              metavar="FILE"),
  
  make_option(c("-o", "--output-dir"), type="character", default=".",
              help="Output directory for plots and results [default=%default]",
              metavar="DIR"),
  
  make_option(c("-p", "--prefix"), type="character", default="site_agreement",
              help="Prefix for output files [default=%default]",
              metavar="STRING"),
  
  make_option(c("--scatter-width"), type="numeric", default=12,
              help="Width of scatter plot PDF in inches [default=%default]",
              metavar="NUM"),
  
  make_option(c("--scatter-height"), type="numeric", default=10,
              help="Height of scatter plot PDF in inches [default=%default]",
              metavar="NUM"),
  
  make_option(c("--cor-width"), type="numeric", default=8,
              help="Width of correlation plot PDF in inches [default=%default]",
              metavar="NUM"),
  
  make_option(c("--cor-height"), type="numeric", default=8,
              help="Height of correlation plot PDF in inches [default=%default]",
              metavar="NUM"),
  
  make_option(c("--hex-bins"), type="integer", default=100,
              help="Number of bins for hexbin plots [default=%default]",
              metavar="INT"),
  
  make_option(c("--min-sites"), type="integer", default=1000,
              help="Minimum common sites required for analysis [default=%default]",
              metavar="INT"),
  
  make_option(c("-v", "--verbose"), action="store_true", default=FALSE,
              help="Print verbose output")
)

opt_parser <- OptionParser(option_list=option_list,
                          description="Site-level agreement analysis for methylation data")
opt <- parse_args(opt_parser)

# Check required arguments
if(is.null(opt$input)) {
  stop("Input RDS file (--input) is required")
}

if(!file.exists(opt$input)) {
  stop(sprintf("Input file does not exist: %s", opt$input))
}

# Create output directory if needed
if(!dir.exists(opt$`output-dir`)) {
  dir.create(opt$`output-dir`, recursive = TRUE)
  if(opt$verbose) cat(sprintf("Created output directory: %s\n", opt$`output-dir`))
}

if(opt$verbose) {
  cat("=== Site-Level Agreement Analysis ===\n")
  cat(sprintf("Input: %s\n", opt$input))
  cat(sprintf("Output directory: %s\n", opt$`output-dir`))
  cat(sprintf("Prefix: %s\n", opt$prefix))
  cat("\n")
}

# ============================================================================
# Load data
# ============================================================================

if(opt$verbose) cat("Loading common sites data...\n")

input_data <- readRDS(opt$input)
common_sites <- input_data$common_sites
tech_names <- input_data$tech_names

if(nrow(common_sites) < opt$`min-sites`) {
  stop(sprintf("Insufficient common sites (%d) for analysis. Minimum required: %d",
               nrow(common_sites), opt$`min-sites`))
}

if(opt$verbose) {
  cat(sprintf("  Technologies: %s\n", paste(tech_names, collapse=", ")))
  cat(sprintf("  Common sites: %s\n", format(nrow(common_sites), big.mark=",")))
  cat("\n")
}

# ============================================================================
# Extract beta values for each technology
# ============================================================================

# Build a data.table with beta values for each technology
beta_cols <- paste0("beta_", tech_names)
betas <- common_sites[, c("chr", "pos", ..beta_cols), with=FALSE]
setnames(betas, beta_cols, tech_names)

if(opt$verbose) {
  cat("Beta value summary:\n")
  for(tech in tech_names) {
    cat(sprintf("  %s: mean=%.3f, median=%.3f, range=[%.3f, %.3f]\n",
                tech, 
                mean(betas[[tech]], na.rm=TRUE),
                median(betas[[tech]], na.rm=TRUE),
                min(betas[[tech]], na.rm=TRUE),
                max(betas[[tech]], na.rm=TRUE)))
  }
  cat("\n")
}

# ============================================================================
# A. Calculate pairwise correlations and metrics
# ============================================================================

if(opt$verbose) cat("Calculating pairwise correlations...\n")

calc_correlations <- function(x, y, name_x, name_y) {
  # Remove any NA values
  valid <- complete.cases(x, y)
  x_valid <- x[valid]
  y_valid <- y[valid]
  
  if(length(x_valid) < 10) {
    warning(sprintf("Insufficient valid data points for %s vs %s", name_x, name_y))
    return(NULL)
  }
  
  data.table(
    Comparison = paste(name_x, "vs", name_y),
    N_sites = length(x_valid),
    Pearson = cor(x_valid, y_valid, method = "pearson"),
    Spearman = cor(x_valid, y_valid, method = "spearman"),
    MAD = mean(abs(x_valid - y_valid)),
    RMSD = sqrt(mean((x_valid - y_valid)^2)),
    Mean_Diff = mean(x_valid - y_valid),
    Median_Diff = median(x_valid - y_valid)
  )
}

# Calculate all pairwise comparisons
comparisons_list <- list()
idx <- 1

for(i in 1:(length(tech_names)-1)) {
  for(j in (i+1):length(tech_names)) {
    comp <- calc_correlations(
      betas[[tech_names[i]]], 
      betas[[tech_names[j]]], 
      tech_names[i], 
      tech_names[j]
    )
    if(!is.null(comp)) {
      comparisons_list[[idx]] <- comp
      idx <- idx + 1
    }
  }
}

comparisons <- rbindlist(comparisons_list)

# ============================================================================
# B. Concordance at methylation thresholds
# ============================================================================

if(opt$verbose) cat("Calculating concordance at methylation thresholds...\n")

calc_concordance <- function(x, y, name_x, name_y) {
  # Remove any NA values
  valid <- complete.cases(x, y)
  x_valid <- x[valid]
  y_valid <- y[valid]
  
  if(length(x_valid) < 10) {
    return(NULL)
  }
  
  # Low methylation: 0-20%
  low_x <- x_valid <= 0.2
  low_y <- y_valid <= 0.2
  low_union <- sum(low_x | low_y)
  low_concordance <- if(low_union > 0) sum(low_x & low_y) / low_union else NA
  
  # Intermediate: 20-80%
  int_x <- x_valid > 0.2 & x_valid < 0.8
  int_y <- y_valid > 0.2 & y_valid < 0.8
  int_union <- sum(int_x | int_y)
  int_concordance <- if(int_union > 0) sum(int_x & int_y) / int_union else NA
  
  # High: 80-100%
  high_x <- x_valid >= 0.8
  high_y <- y_valid >= 0.8
  high_union <- sum(high_x | high_y)
  high_concordance <- if(high_union > 0) sum(high_x & high_y) / high_union else NA
  
  data.table(
    Comparison = paste(name_x, "vs", name_y),
    Low_0_20_pct = low_concordance,
    Low_0_20_n = low_union,
    Intermediate_20_80_pct = int_concordance,
    Intermediate_20_80_n = int_union,
    High_80_100_pct = high_concordance,
    High_80_100_n = high_union
  )
}

concordance_list <- list()
idx <- 1

for(i in 1:(length(tech_names)-1)) {
  for(j in (i+1):length(tech_names)) {
    conc <- calc_concordance(
      betas[[tech_names[i]]], 
      betas[[tech_names[j]]], 
      tech_names[i], 
      tech_names[j]
    )
    if(!is.null(conc)) {
      concordance_list[[idx]] <- conc
      idx <- idx + 1
    }
  }
}

concordance <- rbindlist(concordance_list)

# ============================================================================
# C. Combine metrics and save
# ============================================================================

all_metrics <- merge(comparisons, concordance, by = "Comparison")

metrics_file <- file.path(opt$`output-dir`, paste0(opt$prefix, "_metrics.csv"))
fwrite(all_metrics, metrics_file)

if(opt$verbose) {
  cat("\nPairwise Correlation Metrics:\n")
  print(all_metrics[, .(Comparison, N_sites, Pearson, Spearman, MAD)])
  cat("\n")
}

# ============================================================================
# D. Create scatter plots
# ============================================================================

if(opt$verbose) cat("Generating scatter plots...\n")

create_scatter <- function(x, y, name_x, name_y, metrics_row) {
  # Remove NA values
  valid <- complete.cases(x, y)
  df <- data.frame(x = x[valid], y = y[valid])
  
  # Calculate correlation for subtitle
  pearson_r <- metrics_row$Pearson
  spearman_rho <- metrics_row$Spearman
  mad <- metrics_row$MAD
  
  p <- ggplot(df, aes(x = x, y = y)) +
    geom_hex(bins = opt$`hex-bins`) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", 
                color = "red", linewidth = 1) +
    scale_fill_viridis_c(trans = "log10", name = "Count",
                        labels = scales::comma) +
    labs(
      x = paste(name_x, "methylation (β)"),
      y = paste(name_y, "methylation (β)"),
      title = paste(name_x, "vs", name_y),
      subtitle = sprintf(
        "Pearson r = %.3f | Spearman ρ = %.3f | MAD = %.3f | n = %s", 
        pearson_r, spearman_rho, mad,
        format(nrow(df), big.mark=",")
      )
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 10),
      axis.text = element_text(size = 10),
      legend.position = "right"
    ) +
    coord_fixed(xlim = c(0, 1), ylim = c(0, 1))
  
  return(p)
}

# Generate all pairwise scatter plots
scatter_plots <- list()
idx <- 1

for(i in 1:(length(tech_names)-1)) {
  for(j in (i+1):length(tech_names)) {
    metrics_row <- all_metrics[Comparison == paste(tech_names[i], "vs", tech_names[j])]
    
    if(nrow(metrics_row) > 0) {
      p <- create_scatter(
        betas[[tech_names[i]]], 
        betas[[tech_names[j]]], 
        tech_names[i], 
        tech_names[j],
        metrics_row
      )
      scatter_plots[[idx]] <- p
      idx <- idx + 1
    }
  }
}

# Save scatter plots
scatter_file <- file.path(opt$`output-dir`, paste0(opt$prefix, "_scatter.pdf"))

pdf(scatter_file, width = opt$`scatter-width`, height = opt$`scatter-height`)

# Arrange plots in grid
n_plots <- length(scatter_plots)

if(opt$verbose) cat(sprintf("  Number of plots: %s\n", n_plots))


if(n_plots == 1) {
  print(scatter_plots[[1]])
} else if(n_plots <= 4) {
  print(plot_grid(plotlist = scatter_plots, ncol = 2))
} else if(n_plots <= 6) {
  print(plot_grid(plotlist = scatter_plots, ncol = 3))
} else {
  # Print in batches of 6
  for(i in seq(1, n_plots, by = 6)) {
    batch <- scatter_plots[i:min(i+5, n_plots)]
    print(plot_grid(plotlist = batch, ncol = 3))
  }
}

dev.off()

if(opt$verbose) cat(sprintf("  Scatter plots saved: %s\n", scatter_file))

# ============================================================================
# E. Correlation matrix plot
# ============================================================================

if(opt$verbose) cat("Generating correlation matrix...\n")

cor_matrix_pearson <- cor(betas[, ..tech_names], 
                          use = "complete.obs", 
                          method = "pearson")

cor_matrix_spearman <- cor(betas[, ..tech_names],
                           use = "complete.obs",
                           method = "spearman")

cor_file <- file.path(opt$`output-dir`, paste0(opt$prefix, "_correlations.pdf"))

pdf(cor_file, width = opt$`cor-width`, height = opt$`cor-height`)

# Pearson correlation matrix
corrplot(cor_matrix_pearson, 
         method = "color", 
         type = "upper",
         addCoef.col = "black", 
         number.cex = 1.2,
         tl.col = "black", 
         tl.srt = 45,
         tl.cex = 1.2,
         col = colorRampPalette(c("#67001F", "#B2182B", "#D6604D", 
                                  "#F4A582", "#FDDBC7", "#FFFFFF",
                                  "#D1E5F0", "#92C5DE", "#4393C3",
                                  "#2166AC", "#053061"))(200),
         title = "Pearson Correlation of Methylation β-values",
         mar = c(0,0,2,0),
         cl.cex = 1)

# Spearman correlation matrix
corrplot(cor_matrix_spearman,
         method = "color",
         type = "upper",
         addCoef.col = "black",
         number.cex = 1.2,
         tl.col = "black",
         tl.srt = 45,
         tl.cex = 1.2,
         col = colorRampPalette(c("#67001F", "#B2182B", "#D6604D",
                                  "#F4A582", "#FDDBC7", "#FFFFFF",
                                  "#D1E5F0", "#92C5DE", "#4393C3",
                                  "#2166AC", "#053061"))(200),
         title = "Spearman Correlation of Methylation β-values",
         mar = c(0,0,2,0),
         cl.cex = 1)

dev.off()

if(opt$verbose) cat(sprintf("  Correlation plots saved: %s\n", cor_file))

# ============================================================================
# F. Bland-Altman plots (Difference vs Average)
# ============================================================================

if(opt$verbose) cat("Generating Bland-Altman plots...\n")

create_bland_altman <- function(x, y, name_x, name_y, metrics_row) {
  valid <- complete.cases(x, y)
  x_valid <- x[valid]
  y_valid <- y[valid]
  
  avg <- (x_valid + y_valid) / 2
  diff <- x_valid - y_valid
  
  mean_diff <- mean(diff)
  sd_diff <- sd(diff)
  
  df <- data.frame(avg = avg, diff = diff)
  
  p <- ggplot(df, aes(x = avg, y = diff)) +
    geom_hex(bins = opt$`hex-bins`) +
    geom_hline(yintercept = mean_diff, color = "red", linewidth = 1) +
    geom_hline(yintercept = mean_diff + 1.96*sd_diff, 
               color = "red", linetype = "dashed", linewidth = 0.8) +
    geom_hline(yintercept = mean_diff - 1.96*sd_diff,
               color = "red", linetype = "dashed", linewidth = 0.8) +
    geom_hline(yintercept = 0, color = "black", linetype = "dotted") +
    scale_fill_viridis_c(trans = "log10", name = "Count") +
    labs(
      x = "Average methylation (β)",
      y = sprintf("Difference (%s - %s)", name_x, name_y),
      title = paste("Bland-Altman:", name_x, "vs", name_y),
      subtitle = sprintf(
        "Mean diff = %.4f | SD = %.4f | 95%% limits: [%.3f, %.3f]",
        mean_diff, sd_diff,
        mean_diff - 1.96*sd_diff,
        mean_diff + 1.96*sd_diff
      )
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 9)
    )
  
  return(p)
}

ba_plots <- list()
idx <- 1

for(i in 1:(length(tech_names)-1)) {
  for(j in (i+1):length(tech_names)) {
    metrics_row <- all_metrics[Comparison == paste(tech_names[i], "vs", tech_names[j])]
    
    if(nrow(metrics_row) > 0) {
      p <- create_bland_altman(
        betas[[tech_names[i]]],
        betas[[tech_names[j]]],
        tech_names[i],
        tech_names[j],
        metrics_row
      )
      ba_plots[[idx]] <- p
      idx <- idx + 1
    }
  }
}

ba_file <- file.path(opt$`output-dir`, paste0(opt$prefix, "_bland_altman.pdf"))

pdf(ba_file, width = opt$`scatter-width`, height = opt$`scatter-height`)

n_plots <- length(ba_plots)
if(n_plots == 1) {
  print(ba_plots[[1]])
} else if(n_plots <= 4) {
  print(plot_grid(plotlist = ba_plots, ncol = 2))
} else if(n_plots <= 6) {
  print(plot_grid(plotlist = ba_plots, ncol = 3))
} else {
  for(i in seq(1, n_plots, by = 6)) {
    batch <- ba_plots[i:min(i+5, n_plots)]
    print(plot_grid(plotlist = batch, ncol = 3))
  }
}

dev.off()

if(opt$verbose) cat(sprintf("  Bland-Altman plots saved: %s\n", ba_file))

# ============================================================================
# Summary output
# ============================================================================

cat("\n=== Site-Level Agreement Analysis Complete ===\n")
cat(sprintf("\nResults saved to: %s\n", opt$`output-dir`))
cat(sprintf("  Metrics: %s\n", basename(metrics_file)))
cat(sprintf("  Scatter plots: %s\n", basename(scatter_file)))
cat(sprintf("  Correlation matrix: %s\n", basename(cor_file)))
cat(sprintf("  Bland-Altman plots: %s\n", basename(ba_file)))

cat("\nKey Findings:\n")
cat(sprintf("  Technologies compared: %d\n", length(tech_names)))
cat(sprintf("  Common CpG sites: %s\n", format(nrow(common_sites), big.mark=",")))
cat("\nPairwise Correlations:\n")
print(all_metrics[, .(Comparison, Pearson, Spearman, MAD)])

cat("\nDone!\n")