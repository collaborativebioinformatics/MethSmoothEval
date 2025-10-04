#!/usr/bin/env Rscript


#install.packages("viridis")
#if (!require("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("BSgenome.Hsapiens.UCSC.hg38", force=TRUE)


suppressPackageStartupMessages({
    library(DSS)
    library(bsseq)
    library(dplyr)
    library(data.table)
    library(BiocParallel)
    library(ggplot2)
    library(gridExtra)
    library(scales)
    library(viridis)
    library(rtracklayer)
    library(GenomicRanges)
    library(IRanges)
    library(BSgenome.Hsapiens.UCSC.hg38)
    library(optparse)

})


# Parse command line arguments
option_list <- list(
  make_option("--input1", type="character", help="Input 1st masked CpG DSS BED file"),
  make_option("--input2", type="character", help="Input 2nd masked CpG DSS BED file"),
  make_option("--outputDir", type="character", help="output dir"),
  make_option("--sample1", type="character", help="Sample 1 name"),
  make_option("--sample2", type="character", help="Sample 2 name"),
  make_option("--util", type="character", help="analysis utilities"),
  make_option("--bw", type="character", help="bigWig utilities"),
  make_option("--p_threshold", type="numeric", default=0.05, help="p-value threshold [default %default]"),
  make_option("--delta_threshold", type="numeric", default=0.1, help="delta of methylation difference threshold [default %default]"),
  make_option(c("-v", "--verbose"), action="store_true", default=FALSE,
              help="Print verbose output")
)

opt <- parse_args(OptionParser(option_list=option_list))

# Validate arguments
if(is.null(opt$input1) || is.null(opt$input2) || is.null(opt$sample1) || is.null(opt$sample2) || is.null(opt$util) || is.null(opt$bw) || is.null(opt$outputDir)) {
  stop("Required arguments missing.")
}


source(opt$util)


methylOut <- run_methSmooth_analysis(opt$input1, opt$input2, opt$sample1, opt$sample2, opt$p_threshold, opt$delta_threshold, opt$verbose)


plotss <- create_plots(methylOut, opt$outputDir, opt$verbose)

source(opt$bw)

bigwigs <- create_bigwig(methylOut, opt$outputDir, opt$verbose)


message("Analysis done.")







