# expanded_sample_config.R
# Sample configuration for expanded dataset

# =============================================================================
# SAMPLE FILES AND NAMES
# =============================================================================

sample_files <- c(
  "HG002.ont.chr22.cpg.filter.dss.txt",
  "HG002.pacbio.chr22.cpg.filter.dss.txt", 
  "HG002_EMSeq_LAB01_REP01.chr22.dss.filter.txt",
  "HG002_TruSeq_LAB01_REP02.chr22.dss.filter.txt",
  "HG002_TrueMethylOX_LAB01_REP02.chr22.dss.filter.txt",
  "st001.liver.chr22.ont.dss.filter.txt",
  "st001.liver.chr22.pacbio.dss.filter.txt",
  "st001.lung.chr22.ont.dss.filter.txt",
  "st001.lung.chr22.pacbio.dss.filter.txt",
  "st002.colon.chr22.ont.dss.filter.txt",
  "st002.colon.chr22.pacbio.dss.filter.txt",
  "st002.lung.chr22.ont.dss.filter.txt",
  "st002.lung.chr22.pacbio.dss.filter.txt"
)

sample_names <- c(
  "HG002_ONT",
  "HG002_PacBio",
  "HG002_EMSeq_REP01", 
  "HG002_TruSeq_REP02",
  "HG002_TrueMethylOX_REP02",
  "st001_liver_ONT",
  "st001_liver_PacBio",
  "st001_lung_ONT",
  "st001_lung_PacBio", 
  "st002_colon_ONT",
  "st002_colon_PacBio",
  "st002_lung_ONT",
  "st002_lung_PacBio"
)

# =============================================================================
# EXPERIMENTAL FACTORS
# =============================================================================

individual <- factor(c(
  "HG002", "HG002", "HG002", "HG002", "HG002",  # HG002 samples
  "st001", "st001", "st001", "st001",           # st001 samples
  "st002", "st002", "st002", "st002"            # st002 samples
))

tissue <- factor(c(
  "blood", "blood", "blood", "blood", "blood",  # HG002 (blood/reference)
  "liver", "liver", "lung", "lung",             # st001 samples
  "colon", "colon", "lung", "lung"              # st002 samples
))

platform <- factor(c(
  "ONT", "PacBio", "EMSeq", "TruSeq", "TrueMethylOX",  # HG002 platforms
  "ONT", "PacBio", "ONT", "PacBio",                    # st001 platforms
  "ONT", "PacBio", "ONT", "PacBio"                     # st002 platforms
))

# Technology category (Long-read vs Short-read)
tech_category <- factor(c(
  "Long", "Long", "Short", "Short", "Short",    # HG002
  "Long", "Long", "Long", "Long",               # st001
  "Long", "Long", "Long", "Long"                # st002
))

# =============================================================================
# SAMPLE METADATA
# =============================================================================

sample_info <- data.frame(
  file = sample_files,
  sample = sample_names,
  individual = individual,
  tissue = tissue,
  platform = platform,
  tech_category = tech_category,
  stringsAsFactors = FALSE
)

cat("EXPANDED SAMPLE CONFIGURATION:\n")
cat("==============================\n")
print(sample_info)
cat("\n")

# Summary statistics
cat("DATASET SUMMARY:\n")
cat("================\n")
cat("Total samples:", nrow(sample_info), "\n")
cat("Individuals:", length(unique(individual)), "(", paste(unique(individual), collapse=", "), ")\n")
cat("Tissues:", length(unique(tissue)), "(", paste(unique(tissue), collapse=", "), ")\n")
cat("Platforms:", length(unique(platform)), "(", paste(unique(platform), collapse=", "), ")\n")
cat("Tech categories:", length(unique(tech_category)), "(", paste(unique(tech_category), collapse=", "), ")\n\n")

# Cross-tabulations
cat("PLATFORM x INDIVIDUAL:\n")
print(table(platform, individual))
cat("\nTISSUE x INDIVIDUAL:\n")
print(table(tissue, individual))
cat("\nTECH CATEGORY x TISSUE:\n")
print(table(tech_category, tissue))

cat("Configuration complete! Variables ready for model.matrix().\n")