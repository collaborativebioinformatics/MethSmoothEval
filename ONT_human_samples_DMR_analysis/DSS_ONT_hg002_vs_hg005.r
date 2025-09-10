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
library(AnnotationHub)
library(GenomicRanges)
library(annotatr)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

DIR <- "/scratch/eger/projects/MethSmoothEval"

# can only have 'm' modifications
bed1 <- paste0(DIR, "/analysis/Sarah_analysis/human_ONT/HG002.ont.chr22.cpg.mONLY.bed.gz") 
bed2 <- paste0(DIR, "/analysis/Sarah_analysis/human_ONT/HG005.ont.chr22.cpg.mONLY.bed.gz") 

# load the modkit pileups
bs1 <- read.modkit(bed1, rmZeroCov = FALSE, strandCollapse = FALSE)
bs2 <- read.modkit(bed2, rmZeroCov = FALSE, strandCollapse = FALSE)

# combine the two samples
BSobj <- bsseq::combine(bs1, bs2)
sampleNames(BSobj) <- c("HG002", "HG005")
pData(BSobj)$sample <- c("HG002", "HG005")

################################################################################
## Filter low coverage regions
min_cov <- 5

# keep loci that meet thresholds
cov <- getCoverage(BSobj)
keepLoci <- which(cov[, "HG002"] >= min_cov & cov[, "HG005"] >= min_cov)
length(keepLoci) 
BSobj <- BSobj[keepLoci, ]

################################################################################
# DSS with smoothing
dmlTest.sm <- DMLtest(BSobj, 
                      group1=c("HG002"), 
                      group2=c("HG005"), 
                      smoothing=TRUE)

out <- paste0(DIR, "/analysis/Sarah_analysis/human_ONT/DSS_DMLtestwSmoothing_HG002_vs_HG005_mincov5.tsv")
write.table(dmlTest.sm, out, quote=F, na="", row.names=F, sep="\t")
