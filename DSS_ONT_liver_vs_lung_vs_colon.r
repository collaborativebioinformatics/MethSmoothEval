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

################################################################################
library(bsseq)
library(DSS)
library(BiocParallel)
library(dplyr)
library(tidyr)

DIR <- "/scratch/eger/projects/MethSmoothEval"

# can only have 'm' modifications
bed1 <- paste0(DIR, "/analysis/Jon_analysis/human/st001/st001.liver.chr22.cpg.bed.gz") 
bed2 <- paste0(DIR, "/analysis/Jon_analysis/human/st001/st001.lung.chr22.cpg.bed.gz") 
bed3 <- paste0(DIR, "/analysis/Jon_analysis/human/st002/st002.colon.chr22.cpg.bed.gz") 
bed4 <- paste0(DIR, "/analysis/Jon_analysis/human/st002/st002.lung.chr22.cpg.bed.gz") 

# load the modkit pileups
bs1 <- read.modkit(bed1, rmZeroCov = FALSE, strandCollapse = FALSE)
bs2 <- read.modkit(bed2, rmZeroCov = FALSE, strandCollapse = FALSE)
bs3 <- read.modkit(bed3, rmZeroCov = FALSE, strandCollapse = FALSE)
bs4 <- read.modkit(bed4, rmZeroCov = FALSE, strandCollapse = FALSE)

# combine samples
BSobj1 <- bsseq::combine(bs1, bs2)
sampleNames(BSobj1) <- c("st001.liver", "st001.lung")
pData(BSobj1)$sample  <- c("st001.liver", "st001.lung")

BSobj2 <- bsseq::combine(bs3, bs4)
sampleNames(BSobj2) <- c("st002.colon", "st002.lung")
pData(BSobj2)$sample <- c("st002.colon", "st002.lung")

################################################################################
## Filter low coverage regions
min_cov <- 5

# keep loci that meet thresholds
cov <- getCoverage(BSobj1)
keepLoci <- which(cov[, "st001.liver"] >= min_cov & cov[, "st001.lung"] >= min_cov)
length(keepLoci) 
BSobj1 <- BSobj1[keepLoci, ]

cov <- getCoverage(BSobj2)
keepLoci <- which(cov[, "st002.colon"] >= min_cov & cov[, "st002.lung"] >= min_cov)
length(keepLoci) 
BSobj2 <- BSobj2[keepLoci, ]

################################################################################
# DSS with smoothing
dmlTest.sm1 <- DMLtest(BSobj1, 
                      group1=c("st001.liver"), 
                      group2=c("st001.lung"), 
                      smoothing=TRUE)

out <- paste0(DIR, "/analysis/Jon_analysis/human/st001/DSS_DMLtestwSmoothing_st001_liver_vs_lung_mincov5.tsv")
write.table(dmlTest.sm1, out, quote=F, na="", row.names=F, sep="\t")

dmlTest.sm2 <- DMLtest(BSobj2, 
                      group1=c("st002.colon"), 
                      group2=c("st002.lung"), 
                      smoothing=TRUE)

out <- paste0(DIR, "/analysis/Jon_analysis/human/st002/DSS_DMLtestwSmoothing_st002_colon_vs_lung_mincov5.tsv")
write.table(dmlTest.sm2, out, quote=F, na="", row.names=F, sep="\t")
