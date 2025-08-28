#!/usr/bin/env Rscript

# Convert PacBio or ONT methylation output to DSS-ready format
# Output: chr pos N X   (no header, tab-delimited)
# Compatible with DSS::makeBSseqData()

suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(data.table)
})

# -------- ARGUMENT PARSER --------
option_list <- list(
  make_option(c("-i", "--input"), type="character", help="Input file (PacBio or ONT, can be .gz)"),
  make_option(c("-o", "--output"), type="character", help="Output DSS-ready file"),
  make_option(c("-t", "--type"), type="character", help="Input type: 'pacbio' or 'ont'"),
  make_option(c("-p", "--prob_threshold"), type="double", default=0.5,
              help="ONT probability threshold for methylation call [default %default]")
)

opt <- parse_args(OptionParser(option_list=option_list))

if (is.null(opt$input) || is.null(opt$output) || is.null(opt$type)) {
  stop("Error: must provide --input, --output, and --type (pacbio or ont)")
}

# -------- PACBIO HANDLER --------
convert_pacbio <- function(file, outfile) {
  df <- fread(file, header=TRUE)             # keep the header line
  setnames(df, gsub("^#", "", colnames(df))) # rename '#chrom' → 'chrom'
  
  required_cols <- c("chrom","begin","end","cov","est_mod_count")
  if (!all(required_cols %in% colnames(df))) {
    stop("PacBio file must contain columns: ", paste(required_cols, collapse=", "))
  }
  
  dss_df <- df %>%
    mutate(chr = chrom,
           pos = begin + 1,       # convert 0-based to 1-based
           N   = cov,
           X   = est_mod_count) %>%
    select(chr, pos, N, X)
  
  write.table(dss_df, file=outfile, sep="\t", quote=FALSE,
              row.names=FALSE, col.names=FALSE)
}


# -------- ONT HANDLER --------
convert_ont <- function(file, outfile, prob_threshold=0.5) {
  # fread handles .gz automatically
  df <- fread(file, header=FALSE)#TRUE, comment.char="##")
  
  #required_cols <- c("chr","pos","probability")
  #if (!all(required_cols %in% colnames(df))) {
  #  stop("ONT file must contain columns: ", paste(required_cols, collapse=", "))
  #}
  

  colnames(df) <- c(
    "chrom", "start", "end", "mod_code", "score", "strand",
    "start2", "end2", "color",
    "Nvalid_cov", "fraction_modified", "Nmod", "Ncanonical",
    "Nother_mod", "Ndelete", "Nfail", "Ndiff", "Nnocall"
  )

  
  dss_df <- df %>%
    mutate(chr = chrom,
           pos = start + 1,        # convert 0-based to 1-based
           N   = Nvalid_cov,
           X   = Nmod) %>%
    select(chr, pos, N, X)
  
  write.table(dss_df, file=outfile, sep="\t", quote=FALSE,
              row.names=FALSE, col.names=FALSE)
}

# -------- DISPATCH --------
if (tolower(opt$type) == "pacbio") {
  convert_pacbio(opt$input, opt$output)
  message("Converted PacBio file to DSS-ready format: ", opt$output)
} else if (tolower(opt$type) == "ont") {
  convert_ont(opt$input, opt$output, opt$prob_threshold)
  message("Converted ONT file to DSS-ready format: ", opt$output)
} else {
  stop("Unknown type: must be 'pacbio' or 'ont'")
}

