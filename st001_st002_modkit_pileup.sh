# repeating pileup making on st001 and st002
# original script written by Sarah Eger
# create pileups

# st001 liver and then lung
/home/dnanexus/dist_modkit_v0.5.0_5120ef7/modkit pileup \
    --cpg \
    --ignore h \
    --combine-strands \
    --ref ~/ref/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa \
    --threads 30 \
    "/home/dnanexus/human/st001/st001.liver.chr22.bam" - | \
bgzip -c > "/home/dnanexus/human/st001/st001.liver.chr22.cpg.bed.gz"
tabix -p bed "/home/dnanexus/human/st001/st001.liver.chr22.cpg.bed.gz"

/home/dnanexus/dist_modkit_v0.5.0_5120ef7/modkit pileup \
    --cpg \
    --ignore h \
    --combine-strands \
    --ref ~/ref/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa \
    --threads 30 \
    "/home/dnanexus/human/st001/st001.lung.chr22.bam" - | \
bgzip -c > "/home/dnanexus/human/st001/st001.lung.chr22.cpg.bed.gz"
tabix -p bed "/home/dnanexus/human/st001/st001.lung.chr22.cpg.bed.gz"

# st002 colon and then lung
/home/dnanexus/dist_modkit_v0.5.0_5120ef7/modkit pileup \
    --cpg \
    --ignore h \
    --combine-strands \
    --ref ~/ref/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa \
    --threads 30 \
    "/home/dnanexus/human/st002/st002.colon.chr22.bam" - | \
bgzip -c > "/home/dnanexus/human/st002/st002.colon.chr22.cpg.bed.gz"
tabix -p bed "/home/dnanexus/human/st002/st002.colon.chr22.cpg.bed.gz"

/home/dnanexus/dist_modkit_v0.5.0_5120ef7/modkit pileup \
    --cpg \
    --ignore h \
    --combine-strands \
    --ref ~/ref/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa \
    --threads 30 \
    "/home/dnanexus/human/st002/st002.lung.chr22.bam" - | \
bgzip -c > "/home/dnanexus/human/st002/st002.lung.chr22.cpg.bed.gz"
tabix -p bed "/home/dnanexus/human/st002/st002.lung.chr22.cpg.bed.gz"

# get dmr output for st001 and st002
# st001 liver vs. lung
/home/dnanexus/dist_modkit_v0.5.0_5120ef7/modkit dmr pair \
    -a "/home/dnanexus/human/st001/st001.liver.chr22.cpg.bed.gz" \
    -b "/home/dnanexus/human/st001/st001.lung.chr22.cpg.bed.gz" \
    -o "/home/dnanexus/human/st001/single_base_dmr_st001_liver_vs_lung.bed" \
	--segment st001_liver_vs_lung_segments.txt \
	--fine-grained \
    --ref ~/ref/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa \
    --base C \
    --threads 30 \
    --log-filepath "/home/dnanexus/human/st001/single_base_dmr_st001_liver_vs_lung.log"
	
# st002 colon vs. lung
/home/dnanexus/dist_modkit_v0.5.0_5120ef7/modkit dmr pair \
    -a "/home/dnanexus/human/st002/st002.colon.chr22.cpg.bed.gz" \
    -b "/home/dnanexus/human/st002/st002.lung.chr22.cpg.bed.gz" \
    -o "/home/dnanexus/human/st002/single_base_dmr_st002_colon_vs_lung.bed" \
	--segment st002_colon_vs_lung_segments.txt \
	--fine-grained \
    --ref ~/ref/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa \
    --base C \
    --threads 30 \
    --log-filepath "/home/dnanexus/human/st002/single_base_dmr_st002_colon_vs_lung.log"