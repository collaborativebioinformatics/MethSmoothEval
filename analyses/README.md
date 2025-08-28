### Running the tools 
1) Activate respective environments
2) Check for directory structure - the dataset directory is one dir prior to the bash scripts (restricting to human data)
scripts/
├── run_modkit.sh
└── run_pbcpgtools.sh
datasets/
├── human
│   ├── hg002
│   │   ├── GSM5649436_TruSeq_HG002_LAB01_REP02.bedGraph.gz
│   │   ├── GSM5649520_EMSeq_HG002_LAB01_REP01.bedGraph.gz
│   │   ├── GSM5649543_TrueMethylOX_HG002_LAB01_REP02.bedGraph.gz
│   │   ├── HG002.ont.chr22.bam
│   │   ├── HG002.ont.chr22.bam.bai
│   │   ├── HG002.pacbio.chr22.bam
│   │   └── HG002.pacbio.chr22.bam.bai
│   └── hg005
│       ├── GSM5649427_TruSeq_HG005_LAB01_REP02.bedGraph.gz
│       ├── GSM5649556_TrueMethylBS_HG005_LAB01_REP02.bedGraph.gz
│       ├── GSM5649571_EMSeq_HG005_LAB02_REP01.bedGraph.gz
│       ├── HG005.ont.chr22.bam
│       └── HG005.ont.chr22.bam.bai
└── reference_genomes
    ├── GCA_000001405.15_GRCh38_no_alt_analysis_set.fa
    ├── GCA_000001405.15_GRCh38_no_alt_analysis_set.fa.fai
    ├── GCA_000001405.15_GRCh38_no_alt_analysis_set.mmi
3) Execute the bash scripts

Example :
conda env create -f envs/pbcpgtools.yml
cd scripts/ && bash scripts/run_pbcpgtools.sh

4) Workflow for PacBio data analysis
- Creating Methylation calls : Activate conda env pbcpgtools.yml and execute cpgtools script
  run_pbcpgtools.sh
- Running Metilene for DMR calls : Activate conda env metilene.yml and execute in this order 
  create_pacbio_bedgraph.sh -> create_pacbio_pairs.sh -> run_pacbio_metilene.sh
