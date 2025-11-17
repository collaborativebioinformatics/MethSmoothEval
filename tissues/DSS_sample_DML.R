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
    BSobj,
    cov,
    sample1,
    sample2,
    min_cov = 5,
    OUTDIR = "./",
    cores = 4
) {
    # --- Filter loci by coverage ---
    keepLoci <- which(
        cov[, sample1] >= min_cov &
        cov[, sample2] >= min_cov
    )
    message("Loci kept: ", length(keepLoci))
    
    BSobj2 <- BSobj[keepLoci, ]

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
n_cores <- 30
mincov <- 5

# dirs
DIR <- "/scratch/eger/projects/MethSmoothEval/tissues/"
oDIR <- paste0(DIR, "DSS/")

################################################################################
## Create bsseq object

sampleList <- c("Bulk_FC_Control_02", "hg002_blood", "colo829bl")

# get bed file paths
bed1 <- paste0(DIR, "modkit/pileup/", sampleList[1], ".cpg.bed.gz") 
bed2 <- paste0(DIR, "modkit/pileup/", sampleList[2], ".cpg.bed.gz") 
bed3 <- paste0(DIR, "modkit/pileup/", sampleList[3], ".cpg.bed.gz") 
    
# load the modkit pileups
bs1 <- read.modkit(bed1, rmZeroCov = FALSE, strandCollapse = FALSE)
bs2 <- read.modkit(bed2, rmZeroCov = FALSE, strandCollapse = FALSE)
bs3 <- read.modkit(bed3, rmZeroCov = FALSE, strandCollapse = FALSE)
    
# combine into one object
BSobj1 <- bsseq::combine(bs1, bs2, bs3)
sampleNames(BSobj1) <- sampleList
pData(BSobj1)$sample <- sampleList
  # 29191643 methylation loci
  # 3 samples

# save the object
saveRDS(BSobj1, file = paste0(OUTDIR, "ALL.unphased.rds"))

################################################################################
## Load object
BSobj1 <- readRDS(paste0(OUTDIR, "ALL.unphased.rds"))
all_cov <- getCoverage(BSobj1)

################################################################################
## pairwise

res <- run_DSS_pair(
    BSobj = BSobj1,
    cov = all_cov,
    sample1 = "Bulk_FC_Control_02",
    sample2 = "hg002_blood",
    min_cov = mincov,
    OUTDIR = oDIR,
    cores = n_cores
)
# Loci kept: 28435808

res <- run_DSS_pair(
    BSobj = BSobj1,
    cov = all_cov,
    sample1 = "Bulk_FC_Control_02",
    sample2 = "colo829bl",
    min_cov = mincov,
    OUTDIR = oDIR,
    cores = n_cores
)
# Loci kept: 28370678

res <- run_DSS_pair(
    BSobj = BSobj1,
    cov = all_cov,
    sample1 = "hg002_blood",
    sample2 = "colo829bl",
    min_cov = mincov,
    OUTDIR = oDIR,
    cores = n_cores
)
# Loci kept: 28541221

################################################################################
## 1 Brain v. 2 Blood

# keep loci that meet cov threshold for all samples
keepLoci <- which(all_cov[, sampleNames(BSobj1)[1]] >= min_cov & 
                 all_cov[, sampleNames(BSobj1)[2]] >= min_cov &
                all_cov[, sampleNames(BSobj1)[3]] >= min_cov)   
length(keepLoci) 
BSobj2 <- BSobj1[keepLoci, ]
# 28286171

# DSS with smoothing
dmlTest.sm1 <- DMLtest(BSobj2, 
                      group1=c(sampleNames(BSobj1)[1]), 
                      group2=c(sampleNames(BSobj1)[2], sampleNames(BSobj1)[3]), 
                      smoothing=TRUE,
                      ncores=cores)

# compressed output filename
out1 <- paste0(OUTDIR, "BrainvBlood_DMLtestwSmoothing.tsv.gz")
fwrite(dmlTest.sm1, out1, sep = "\t", quote = FALSE, na = "")

# DSS without smoothing
dmlTest.sm2 <- DMLtest(BSobj2, 
                      group1=c(sampleNames(BSobj1)[1]), 
                      group2=c(sampleNames(BSobj1)[2], sampleNames(BSobj1)[3]), 
                      equal.disp=TRUE,
                      smoothing=FALSE,
                      ncores=cores)

# compressed output filename
out2 <- paste0(OUTDIR, "BrainvBlood_DMLtest.tsv.gz")
fwrite(dmlTest.sm2, out2, sep = "\t", quote = FALSE, na = "")
################################################################################