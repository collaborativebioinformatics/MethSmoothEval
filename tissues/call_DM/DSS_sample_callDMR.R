################################################################################
library(bsseq)
library(DSS)
library(dplyr)
library(tidyr)
library(data.table) 

# dirs
DIR <- "/scratch/eger/projects/MethSmoothEval/tissues/"
oDIR <- paste0(DIR, "DSS/")

sample_names <- c("Bulk_FC_Control_02", "hg002_blood", "colo829bl")

################################################################################
## Call DMRs
# callDMR(DMLresult, delta=0, p.threshold=1e-5,
#         minlen=50, minCG=3, dis.merge=100, pct.sig=0.5)

# load the DML results
sample1 <- sample_names[1]
sample2 <- sample_names[2]

# load the smoothed results
fname1 <- file.path(oDIR,
        paste0(sample1, "_v_", sample2, "_DMLtestwSmoothing.tsv.gz"))
dml_smoothed <- data.table::fread(fname1)

# default settings
dmr1 <- callDMR(dml_smoothed) # 259161

# sort DMRs
chrom_order <- c(paste0("chr", 1:22), "chrX", "chrY")
dmr1_sorted <- dmr1 %>%
  mutate(chr_factor = factor(chr, levels = chrom_order)) %>%
  arrange(chr_factor, start) %>%
  select(-chr_factor)

out1 <- file.path(
    oDIR,
    paste0(sample1, "_v_", sample2, "_DMLtestwSmoothing_DMRs.tsv")
)
data.table::fwrite(dmr1_sorted, out1, sep = "\t",
                   quote = FALSE, na = "")

################################################################################
# load the DML results
DIR <- "/scratch/eger/projects/MethSmoothEval/tissues/data_transfer/CARD_dmr_outputs/PPMI/DSS/PPMI_3404"
oDIR <- "/scratch/eger/projects/MethSmoothEval/tissues/DSS/"

sample1 <- "PPMI_3404"
sample2 <- "HBCC_81951_FTX"

# load the smoothed results
fname1 <- file.path(DIR,
        paste0(sample1, "_v_", sample2, "_DMLtestwSmoothing.tsv.gz"))
dml_smoothed <- data.table::fread(fname1)

# default settings
dmr1 <- callDMR(dml_smoothed, delta=0, p.threshold=1e-5,
             minlen=50, minCG=3, dis.merge=100, pct.sig=0.5) # 136827

# sort DMRs
dmr1 <- dmr1[order(as.numeric(rownames(dmr1))),]

out1 <- file.path(
    oDIR,
    paste0(sample1, "_v_", sample2, "_DMLtestwSmoothing_DMRs.tsv")
)
data.table::fwrite(dmr1, out1, sep = "\t",
                   quote = FALSE, na = "")


sample1 <- "PPMI_3404"
sample2 <- "HBCC_82044_FTX"

# load the smoothed results
fname1 <- file.path(DIR,
        paste0(sample1, "_v_", sample2, "_DMLtestwSmoothing.tsv.gz"))
dml_smoothed <- data.table::fread(fname1)

# default settings
dmr1 <- callDMR(dml_smoothed) # 129703

# sort DMRs
dmr1 <- dmr1[order(as.numeric(rownames(dmr1))),]

out1 <- file.path(
    oDIR,
    paste0(sample1, "_v_", sample2, "_DMLtestwSmoothing_DMRs.tsv")
)
data.table::fwrite(dmr1, out1, sep = "\t",
                   quote = FALSE, na = "")

