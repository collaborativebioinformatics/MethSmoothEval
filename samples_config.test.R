# expanded_sample_config.R
# Sample configuration for expanded dataset

# =============================================================================
# SAMPLE FILES AND NAMES
# =============================================================================

sample_files <- c(
  "st001.liver.chr22.ont.dss.filter.txt",
  "st001.liver.chr22.pacbio.dss.filter.txt",
  "st001.lung.chr22.ont.dss.filter.txt",
  "st001.lung.chr22.pacbio.dss.filter.txt"
)

sample_names <- c(
  "st001_liver_ONT",
  "st001_liver_PacBio",
  "st001_lung_ONT",
  "st001_lung_PacBio" 
)

# =============================================================================
# EXPERIMENTAL FACTORS
# =============================================================================

individual <- factor(c(
  # HG002 samples
  "st001", "st001", "st001", "st001"          # st001 samples
          # st002 samples
))

tissue <- factor(c(
  # HG002 (blood/reference)
  "liver", "liver", "lung", "lung"           # st001 samples
         # st002 samples
))

platform <- factor(c(
 # HG002 platforms
  "ONT", "PacBio", "ONT", "PacBio"                  # st001 platforms
               # st002 platforms
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