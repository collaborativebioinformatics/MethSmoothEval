#!/usr/bin/env Rscript

# Convert PacBio, ONT, or Short Read methylation output to DSS-ready format
# Output: chr pos N X   (no header, tab-delimited)
# Compatible with DSS::makeBSseqData()

suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(data.table)
})

# -------- ARGUMENT PARSER --------
option_list <- list(
  make_option(c("-i", "--input"), type="character", help="Input file (PacBio, ONT, or Short Read, can be .gz)"),
  make_option(c("-o", "--output"), type="character", help="Output DSS-ready file"),
  make_option(c("-t", "--type"), type="character", help="Input type: 'pacbio', 'ont', or 'sr' (short read)"),
  make_option(c("-p", "--prob_threshold"), type="double", default=0.5,
              help="ONT probability threshold for methylation call [default %default]"),
  make_option(c("-c", "--min_coverage"), type="integer", default=5,
              help="Minimum coverage threshold [default %default]"),
  make_option(c("-v", "--verbose"), action="store_true", default=FALSE,
              help="Print verbose output")
)

opt <- parse_args(OptionParser(option_list=option_list))

if (is.null(opt$input) || is.null(opt$output) || is.null(opt$type)) {
  stop("Error: must provide --input, --output, and --type (pacbio, ont, or sr)")
}

if (opt$verbose) {
  message("Converting ", toupper(opt$type), " methylation data to DSS format...")
  message("Input file: ", opt$input)
  message("Output file: ", opt$output)
  if (tolower(opt$type) == "ont") {
    message("ONT probability threshold: ", opt$prob_threshold)
  }
  message("Minimum coverage: ", opt$min_coverage)
}

# -------- PACBIO HANDLER --------
convert_pacbio <- function(file, outfile, min_coverage = 5, verbose = FALSE) {
  if (verbose) message("Reading PacBio file...")
  
  df <- fread(file, header=TRUE)
  setnames(df, gsub("^#", "", colnames(df))) # rename '#chrom' → 'chrom'
  
  required_cols <- c("chrom","begin","end","cov","est_mod_count")
  if (!all(required_cols %in% colnames(df))) {
    stop("PacBio file must contain columns: ", paste(required_cols, collapse=", "))
  }
  
  if (verbose) {
    message("Initial sites: ", nrow(df))
    message("Mean coverage: ", round(mean(df$cov), 2))
  }
  
  dss_df <- df %>%
    mutate(chr = chrom,
           pos = begin + 1,       # convert 0-based to 1-based
           N   = cov,
           X   = est_mod_count) %>%
    filter(N >= min_coverage) %>%  # Apply coverage filter
    select(chr, pos, N, X)
  
  if (verbose) {
    message("Sites after filtering (coverage >= ", min_coverage, "): ", nrow(dss_df))
    message("Mean methylation rate: ", round(mean(dss_df$X / dss_df$N), 3))
  }
  
  write.table(dss_df, file=outfile, sep="\t", quote=FALSE,
              row.names=FALSE, col.names=FALSE)
  
  return(nrow(dss_df))
}

# -------- ONT HANDLER --------
convert_ont <- function(file, outfile, prob_threshold = 0.5, min_coverage = 3, verbose = FALSE) {
  if (verbose) message("Reading ONT file...")
  
  # fread handles .gz automatically
  df <- fread(file, header=FALSE)
  
  colnames(df) <- c(
    "chrom", "start", "end", "mod_code", "score", "strand",
    "start2", "end2", "color",
    "Nvalid_cov", "fraction_modified", "Nmod", "Ncanonical",
    "Nother_mod", "Ndelete", "Nfail", "Ndiff", "Nnocall"
  )
  
  if (verbose) {
    message("Initial sites: ", nrow(df))
    message("Mean coverage: ", round(mean(df$Nvalid_cov), 2))
  }
  
  dss_df <- df %>%
    mutate(chr = chrom,
           pos = start + 1,        # convert 0-based to 1-based
           N   = Nvalid_cov,
           X   = Nmod) %>%
    filter(N >= min_coverage) %>%  # Apply coverage filter
    select(chr, pos, N, X)
  
  if (verbose) {
    message("Sites after filtering (coverage >= ", min_coverage, "): ", nrow(dss_df))
    message("Mean methylation rate: ", round(mean(dss_df$X / dss_df$N), 3))
  }
  
  write.table(dss_df, file=outfile, sep="\t", quote=FALSE,
              row.names=FALSE, col.names=FALSE)
  
  return(nrow(dss_df))
}

# -------- SHORT READ HANDLER --------
convert_short_read <- function(file, outfile, min_coverage = 5, verbose = FALSE) {
  if (verbose) message("Reading Short Read bedgraph file...")
  
  # Try to read the file and detect format
  df <- tryCatch({
    fread(file, header = FALSE)
  }, error = function(e) {
    stop("Error reading file: ", e$message)
  })
  
  if (verbose) message("Read ", nrow(df), " regions")
  
  # Check format and assign column names
  if (ncol(df) == 6) {
    # Standard format: chr, start, end, methylation_pct, methylated_reads, unmethylated_reads
    setnames(df, c("chr", "start", "end", "methylation_pct", "methylated", "unmethylated"))
    
    if (verbose) {
      message("Detected 6-column format: chr, start, end, meth_pct, meth_reads, unmeth_reads")
      message("Mean coverage: ", round(mean(df$methylated + df$unmethylated), 2))
    }
    
    # Convert to DSS format
    dss_df <- df %>%
      mutate(chr = chr,
             pos = round((start + end) / 2),  # midpoint
             N = methylated + unmethylated,
             X = methylated) %>%
      filter(N >= min_coverage,
             X >= 0, X <= N,  # Valid methylation counts
             !is.na(N), !is.na(X)) %>%
      select(chr, pos, N, X)
      
  } else if (ncol(df) == 5) {
    # Alternative format: chr, start, end, methylated_reads, unmethylated_reads
    setnames(df, c("chr", "start", "end", "methylated", "unmethylated"))
    
    if (verbose) {
      message("Detected 5-column format: chr, start, end, meth_reads, unmeth_reads")
      message("Mean coverage: ", round(mean(df$methylated + df$unmethylated), 2))
    }
    
    dss_df <- df %>%
      mutate(chr = chr,
             pos = round((start + end) / 2),
             N = methylated + unmethylated,
             X = methylated) %>%
      filter(N >= min_coverage,
             X >= 0, X <= N,
             !is.na(N), !is.na(X)) %>%
      select(chr, pos, N, X)
      
  } else if (ncol(df) == 4) {
    # Minimal format: chr, start, end, methylation_pct
    setnames(df, c("chr", "start", "end", "methylation_pct"))
    
    if (verbose) {
      message("Detected 4-column format: chr, start, end, meth_pct")
      message("Warning: Estimating coverage from methylation percentage")
    }
    
    # Estimate coverage (not ideal but sometimes necessary)
    default_coverage <- 10
    dss_df <- df %>%
      mutate(chr = chr,
             pos = round((start + end) / 2),
             N = default_coverage,
             X = round(methylation_pct * default_coverage / 100)) %>%
      filter(N >= min_coverage) %>%
      select(chr, pos, N, X)
      
  } else {
    stop("Unexpected number of columns: ", ncol(df), ". Expected 4, 5, or 6 columns")
  }
  
  # Sort by genomic coordinates
  if (verbose) message("Sorting by genomic coordinates...")
  
  # Convert chromosome to factor for proper sorting
  chr_levels <- paste0("chr", c(1:22, "X", "Y", "M"))
  chr_levels <- c(chr_levels, setdiff(unique(dss_df$chr), chr_levels))
  
  dss_df$chr <- factor(dss_df$chr, levels = chr_levels)
  dss_df <- dss_df[order(dss_df$chr, dss_df$pos)]
  dss_df$chr <- as.character(dss_df$chr)
  
  # Remove duplicates (same position)
  n_before <- nrow(dss_df)
  dss_df <- dss_df[!duplicated(paste(dss_df$chr, dss_df$pos))]
  n_after <- nrow(dss_df)
  
  if (verbose && n_before != n_after) {
    message("Removed ", n_before - n_after, " duplicate positions")
  }
  
  if (verbose) {
    message("Final sites: ", nrow(dss_df))
    message("Chromosomes: ", length(unique(dss_df$chr)))
    message("Mean coverage: ", round(mean(dss_df$N), 2))
    message("Mean methylation rate: ", round(mean(dss_df$X / dss_df$N), 3))
  }
  
  write.table(dss_df, file=outfile, sep="\t", quote=FALSE,
              row.names=FALSE, col.names=FALSE)
  
  return(nrow(dss_df))
}

# -------- DISPATCH --------
sites_converted <- 0

if (tolower(opt$type) == "pacbio") {
  sites_converted <- convert_pacbio(opt$input, opt$output, opt$min_coverage, opt$verbose)
  message("Converted PacBio file to DSS-ready format: ", opt$output)
  
} else if (tolower(opt$type) == "ont") {
  sites_converted <- convert_ont(opt$input, opt$output, opt$prob_threshold, opt$min_coverage, opt$verbose)
  message("Converted ONT file to DSS-ready format: ", opt$output)
  
} else if (tolower(opt$type) == "sr") {
  sites_converted <- convert_short_read(opt$input, opt$output, opt$min_coverage, opt$verbose)
  message("Converted Short Read file to DSS-ready format: ", opt$output)
  
} else {
  stop("Unknown type: must be 'pacbio', 'ont', or 'sr'")
}

if (opt$verbose) {
  message("Conversion complete! ", sites_converted, " sites written to output file.")
}

# -------- R FUNCTION VERSIONS (for interactive use) --------
if (interactive()) {
  
  #' Convert PacBio methylation data to DSS format
  convert_pacbio_to_dss <- function(input_file, output_file, min_coverage = 5, verbose = TRUE) {
    convert_pacbio(input_file, output_file, min_coverage, verbose)
  }
  
  #' Convert ONT methylation data to DSS format  
  convert_ont_to_dss <- function(input_file, output_file, prob_threshold = 0.5, min_coverage = 5, verbose = TRUE) {
    convert_ont(input_file, output_file, prob_threshold, min_coverage, verbose)
  }
  
  #' Convert Short Read bedgraph to DSS format
  convert_sr_to_dss <- function(input_file, output_file, min_coverage = 5, verbose = TRUE) {
    convert_short_read(input_file, output_file, min_coverage, verbose)
  }
  
  #' Batch convert multiple files to DSS format
  batch_convert_to_dss <- function(input_files, input_types, output_dir = NULL, min_coverage = 5, verbose = TRUE) {
    
    if (length(input_files) != length(input_types)) {
      stop("input_files and input_types must have the same length")
    }
    
    if (!is.null(output_dir) && !dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }
    
    results <- data.frame(
      input_file = input_files,
      input_type = input_types,
      output_file = character(length(input_files)),
      sites_converted = integer(length(input_files)),
      status = character(length(input_files)),
      stringsAsFactors = FALSE
    )
    
    for (i in seq_along(input_files)) {
      input_file <- input_files[i]
      input_type <- tolower(input_types[i])
      
      # Generate output filename
      if (!is.null(output_dir)) {
        output_file <- file.path(output_dir, gsub("\\.(bed|bedgraph|txt|gz)$", ".dss.txt", basename(input_file)))
      } else {
        output_file <- gsub("\\.(bed|bedgraph|txt|gz)$", ".dss.txt", input_file)
      }
      
      results$output_file[i] <- output_file
      
      if (verbose) {
        message("Processing file ", i, " of ", length(input_files), ": ", basename(input_file))
      }
      
      # Convert based on type
      tryCatch({
        if (input_type == "pacbio") {
          sites <- convert_pacbio(input_file, output_file, min_coverage, verbose)
        } else if (input_type == "ont") {
          sites <- convert_ont(input_file, output_file, 0.5, min_coverage, verbose)
        } else if (input_type == "sr") {
          sites <- convert_short_read(input_file, output_file, min_coverage, verbose)
        } else {
          stop("Unknown input type: ", input_type)
        }
        
        results$sites_converted[i] <- sites
        results$status[i] <- "SUCCESS"
        
      }, error = function(e) {
        results$sites_converted[i] <- 0
        results$status[i] <- paste("ERROR:", e$message)
        if (verbose) message("ERROR: ", e$message)
      })
    }
    
    if (verbose) {
      successful <- sum(results$status == "SUCCESS")
      message("Batch conversion complete: ", successful, " of ", length(input_files), " files converted successfully")
      print(results)
    }
    
    return(results)
  }
  
  message("Interactive functions loaded:")
  message("  convert_pacbio_to_dss(input, output)")
  message("  convert_ont_to_dss(input, output)")
  message("  convert_sr_to_dss(input, output)")
  message("  batch_convert_to_dss(files, types, output_dir)")
}
