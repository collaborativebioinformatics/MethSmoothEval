# smooth_bsobj.R
# Apply smoothing to any BSobj input

smooth_bsobj <- function(bsobj, cores = 2, save_output = TRUE, output_name = NULL) {
  
  if(class(bsobj)[1] != "BSseq") {
    stop("Input must be a BSseq object")
  }
  
  cat(paste("Smoothing", ncol(bsobj), "samples with", cores, "cores...\n"))
  start_time <- Sys.time()
  
  # Apply smoothing
  bsobj_smooth <- BSmooth(bsobj, BPPARAM = MulticoreParam(workers = cores))
  
  end_time <- Sys.time()
  duration <- round(difftime(end_time, start_time, units="mins"), 2)
  cat(paste("Completed in", duration, "minutes\n"))
  
  # Save if requested
  if(save_output) {
    if(is.null(output_name)) {
      output_name <- paste0("BSobj_smoothed_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".RData")
    }
    save(bsobj_smooth, file = output_name)
    cat(paste("Saved to", output_name, "\n"))
  }
  
  return(bsobj_smooth)
}

# Usage examples:
# smoothed <- smooth_bsobj(BSobj_raw)
# smoothed <- smooth_bsobj(my_bsobj, cores = 4, save_output = FALSE)