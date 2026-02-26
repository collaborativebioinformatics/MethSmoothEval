################################################################################
## UXM regions have a range of 5-54 CpGs (see histogram in `UXM_markers.ipynb`).

# conda activate bs
################################################################################
library(bsseq)
library(DSS)
library(dplyr)
library(tidyr)
library(data.table) 

################################################################################
run_callDMR_atlas_settings <- function(
    sample1,
    sample2,
    input_dir,
    output_dir) {

    # No smoothing & smoothing file
    fname1 <- file.path(input_dir,
            paste0(sample1, "_v_", sample2, "_DMLtest.tsv.gz"))
    
    fname2 <- file.path(input_dir,
            paste0(sample1, "_v_", sample2, "_DMLtestwSmoothing.tsv.gz"))
    
    dml1 <- data.table::fread(fname1)
    dml2 <- data.table::fread(fname2)

    # settings determined by atlas region sizes
    dmr1 <- callDMR(dml1, 
                    delta=0, 
                    p.threshold=1e-5,
                    minlen=0, # default = 50
                    minCG=4, # default = 3
                    dis.merge=100, 
                    pct.sig=0.5) 
    
    dmr2 <- callDMR(dml2, 
                    delta=0, 
                    p.threshold=1e-5,
                    minlen=0, # default = 50
                    minCG=4, # default = 3
                    dis.merge=100, 
                    pct.sig=0.5)  
    
    print(paste0(sample1, "_v_", sample2))
    print(paste0("No smoothing: ", dim(dmr1)[1]))
    print(paste0("With smoothing: ", dim(dmr2)[1]))

    # sort DMRs
    dmr1 <- dmr1[order(as.numeric(rownames(dmr1))),]
    dmr2 <- dmr2[order(as.numeric(rownames(dmr2))),]

    # save files
    out1 <- file.path(output_dir,
        paste0(sample1, "_v_", sample2, "_DMLtest_DMRs.tsv"))
    out2 <- file.path(output_dir,
        paste0(sample1, "_v_", sample2, "_DMLtestwSmoothing_DMRs.tsv"))
    
    data.table::fwrite(dmr1, out1, sep = "\t", quote = FALSE, na = "")
    data.table::fwrite(dmr2, out2, sep = "\t", quote = FALSE, na = "")
    
}
       
# dirs
iDIR <- "/scratch/eger/projects/MethSmoothEval/tissues/DSS/"
oDIR <- "/scratch/eger/projects/MethSmoothEval/tissues/DM_analysis/UXM/DSS_DMRs/"
if (!dir.exists(oDIR)) dir.create(oDIR)

################################################################################
# Brain v. Blood comparisons

run_callDMR_atlas_settings("Bulk_FC_Control_02",
               "hg002_blood",
               iDIR,
               oDIR)
# [1] "Bulk_FC_Control_02_v_hg002_blood"
# [1] "No smoothing: 120434"
# [1] "With smoothing: 233178"

run_callDMR_atlas_settings("PPMI_3404",
               "HBCC_81951_FTX",
               iDIR,
               oDIR)
# [1] "PPMI_3404_v_HBCC_81951_FTX"
# [1] "No smoothing: 11141"
# [1] "With smoothing: 117956"

run_callDMR_atlas_settings("PPMI_3404",
               "HBCC_82044_FTX",
               iDIR,
               oDIR)
# [1] "PPMI_3404_v_HBCC_82044_FTX"
# [1] "No smoothing: 9761"
# [1] "With smoothing: 111657"