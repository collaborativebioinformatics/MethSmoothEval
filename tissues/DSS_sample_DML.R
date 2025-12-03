################################################################################
## The bsseq User’s Guide
# https://www.bioconductor.org/packages/devel/bioc/vignettes/bsseq/inst/doc/bsseq.html

## The DSS User’s Guide
# https://bioconductor.org/packages/devel/bioc/vignettes/DSS/inst/doc/DSS.html

################################################################################
## conda activate bs
# conda install bioconda::bioconductor-bsseq
# conda install bioconda::bioconductor-dss
# conda install bioconda::bioconductor-annotationhub
# conda install bioconda::bioconductor-annotatr
# conda install bioconda::bioconductor-txdb.hsapiens.ucsc.hg38.knowngene
# conda install bioconda::bioconductor-org.hs.eg.db

# DSS_2.54.0                 
# bsseq_1.42.0  

################################################################################
library(bsseq)
library(DSS)
library(BiocParallel)
library(dplyr)
library(tidyr)
library(data.table) 

# functions
run_DSS_pair <- function(
    obj1,
    obj2,
    samplenames,
    min_cov = 5,
    OUTDIR = "./",
    cores = 4
) {
    # make object
    BSobj1 <- bsseq::combine(obj1, obj2)
    sampleNames(BSobj1) <- samplenames
    pData(BSobj1)$sample <- samplenames
    
    sample1 <- samplenames[1]
    sample2 <- samplenames[2]
    
    cov <- getCoverage(BSobj1)
    
    # --- Filter loci by coverage ---
    keepLoci <- which(
        cov[, sample1] >= min_cov &
        cov[, sample2] >= min_cov
    )
    message("Loci kept: ", length(keepLoci))
    
    BSobj2 <- BSobj1[keepLoci, ]

    # --- Run DSS with smoothing ---
    message("Running DMLtest with smoothing...")
    dml_smoothed <- DMLtest(
        BSobj2,
        group1 = c(sample1),
        group2 = c(sample2),
        smoothing = TRUE,
        ncores = cores
    )

    out1 <- file.path(
        OUTDIR,
        paste0(sample1, "_v_", sample2, "_DMLtestwSmoothing.tsv.gz")
    )
    data.table::fwrite(dml_smoothed, out1, sep = "\t",
                       quote = FALSE, na = "")
    message("Saved smoothed DML results to: ", out1)

    # --- Run DSS without smoothing ---
    message("Running DMLtest without smoothing...")
    dml_nosmooth <- DMLtest(
        BSobj2,
        group1 = c(sample1),
        group2 = c(sample2),
        equal.disp = TRUE,
        smoothing = FALSE,
        ncores = cores
    )

    out2 <- file.path(
        OUTDIR,
        paste0(sample1, "_v_", sample2, "_DMLtest.tsv.gz")
    )
    data.table::fwrite(dml_nosmooth, out2, sep = "\t",
                       quote = FALSE, na = "")
    message("Saved non-smoothed DML results to: ", out2)

    return(list(
        smoothed = dml_smoothed,
        nosmooth = dml_nosmooth,
        keepLoci = keepLoci
    ))
}

# settings
n_cores <- 40
mincov <- 5

# dirs
DIR <- "/scratch/eger/projects/MethSmoothEval/tissues/"
oDIR <- paste0(DIR, "DSS/")

################################################################################
## Create bsseq object

sample_names <- c("Bulk_FC_Control_02", "hg002_blood", "colo829bl")

# get bed file paths
bed1 <- paste0(DIR, "modkit/pileup/", sample_names[1], ".cpg.bed.gz") 
bed2 <- paste0(DIR, "modkit/pileup/", sample_names[2], ".cpg.bed.gz") 
bed3 <- paste0(DIR, "modkit/pileup/", sample_names[3], ".cpg.bed.gz") 
    
# load the modkit pileups
bs1 <- read.modkit(bed1, strandCollapse = FALSE)
bs2 <- read.modkit(bed2, strandCollapse = FALSE)
bs3 <- read.modkit(bed3, strandCollapse = FALSE)

################################################################################
## pairwise

res <- run_DSS_pair(bs1, bs2,
    sample_names[1:2],
    min_cov = mincov,
    OUTDIR = oDIR,
    cores = n_cores
)
# Loci kept: 28422002

res <- run_DSS_pair(bs1, bs3,
    c(sample_names[1], sample_names[3]),
    min_cov = mincov,
    OUTDIR = oDIR,
    cores = n_cores
)
# Loci kept: 28370678

res <- run_DSS_pair(bs2, bs3,
    sample_names[2:3],
    min_cov = mincov,
    OUTDIR = oDIR,
    cores = n_cores
)
# Loci kept: 28541221

################################################################################
