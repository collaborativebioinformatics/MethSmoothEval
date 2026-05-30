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

# settings
cores <- 30
min_cov <- 5

# dirs
DIR <- "/scratch/eger/projects/MethSmoothEval/tissues/"
OUTDIR <- paste0(DIR, "DSS/")

sampleList <- c("Bulk_FC_Control_02", "hg002_blood", "colo829bl")

for (sample_id in sampleList) {
    # get bed file paths
    bed1 <- paste0(DIR, "modkit/pileup/", sample_id, ".cpg.hp_1.bed.gz") 
    bed2 <- paste0(DIR, "modkit/pileup/", sample_id, ".cpg.hp_2.bed.gz") 
    
    # load the modkit pileups
    bs1 <- read.modkit(bed1, rmZeroCov = FALSE, strandCollapse = FALSE)
    bs2 <- read.modkit(bed2, rmZeroCov = FALSE, strandCollapse = FALSE)
    
    # combine into one object
    BSobj1 <- bsseq::combine(bs1, bs2)
    sampleNames(BSobj1) <- c(paste0(sample_id, ".hp_1"), paste0(sample_id, ".hp_2"))
    pData(BSobj1)$sample  <- c(paste0(sample_id, ".hp_1"), paste0(sample_id, ".hp_2"))

    # save the object
    saveRDS(BSobj1, file = paste0(OUTDIR, sample_id, ".hp.rds"))
    
    # keep loci that meet cov threshold for both HPs
    cov <- getCoverage(BSobj1)
    keepLoci <- which(cov[, sampleNames(BSobj1)[1]] >= min_cov & cov[, sampleNames(BSobj1)[2]] >= min_cov)
    length(keepLoci) 
    BSobj1 <- BSobj1[keepLoci, ]
    
    # DSS with smoothing
    dmlTest.sm1 <- DMLtest(BSobj1, 
                          group1=c(sampleNames(BSobj1)[1]), 
                          group2=c(sampleNames(BSobj1)[2]), 
                          smoothing=TRUE,
                          ncores=cores)

    # compressed output filename
    out1 <- paste0(OUTDIR, sample_id, ".hp_DMLtestwSmoothing.tsv.gz")
    fwrite(dmlTest.sm1, out1, sep = "\t", quote = FALSE, na = "")

    # DSS without smoothing
    dmlTest.sm2 <- DMLtest(BSobj1, 
                          group1=c(sampleNames(BSobj1)[1]), 
                          group2=c(sampleNames(BSobj1)[2]), 
                          equal.disp=TRUE,
                          smoothing=FALSE,
                          ncores=cores)
    
    # compressed output filename
    out2 <- paste0(OUTDIR, sample_id, ".hp_DMLtest.tsv.gz")
    fwrite(dmlTest.sm2, out2, sep = "\t", quote = FALSE, na = "")
}
