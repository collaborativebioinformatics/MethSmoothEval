# Snakefile

import os
import json
from snakemake.utils import min_version

# Enforce minimum Snakemake version
min_version("7.0.0")

# Load configuration
configfile: config["config_file"] #"config.json"


# Container settings - applied to all rules by default
container: config["container"]



# Helper function to extract information from input JSON
def get_input_files():
    with open(config["input_json"]) as f:
        data = json.load(f)
    return data

# Get input data
input_data = get_input_files()

SAMPLE_NAMES = list(input_data["bam_files"].keys())



BED_FILE = input_data["bed_file"]
R_SCRIPT = input_data["r_script"]
DSS_CONVERT_SCRIPT = input_data.get("dss_convert_script", "scripts/convert_to_dss.R")

# Generate all possible pairs of BAM files for R analysis
SAMPLE_PAIRS = []
for i in range(len(SAMPLE_NAMES)):
    for j in range(i+1, len(SAMPLE_NAMES)):
        SAMPLE_PAIRS.append((SAMPLE_NAMES[i], SAMPLE_NAMES[j]))

#print(SAMPLE_NAMES)

#print(SAMPLE_PAIRS)


# Define output directory structure
RESULTS_DIR = config.get("results_dir", "results")
CpG_DIR = os.path.join(RESULTS_DIR, "cpg_scores")
MASKED_DIR = os.path.join(RESULTS_DIR, "masked_cpg")
DSS_DIR = os.path.join(RESULTS_DIR, "dss_format")
ANALYSIS_DIR = os.path.join(RESULTS_DIR, "r_analysis")

# Target rule
rule all:
    input:
        # CpG score files
        expand(os.path.join(CpG_DIR, "{sample}.cpg.bed.gz"), sample=SAMPLE_NAMES),
        # Masked CpG files
        expand(os.path.join(MASKED_DIR, "{sample}.masked.cpg.bed"), sample=SAMPLE_NAMES),
        #DSS format
        expand(os.path.join(DSS_DIR, "{sample}.dss.bed"), sample=SAMPLE_NAMES),
        # R analysis results for all pairs
        expand(os.path.join(ANALYSIS_DIR, "{sample1}_vs_{sample2}_pb_methylation_analysis_grid.png"), zip, 
               sample1=[pair[0] for pair in SAMPLE_PAIRS],
               sample2=[pair[1] for pair in SAMPLE_PAIRS])

# Rule to generate CpG scores from BAM files
rule generate_cpg_scores:
    input:
        bam = lambda wildcards: input_data["bam_files"][wildcards.sample]
    output:
        cpg = os.path.join(CpG_DIR, "{sample}.cpg.bed.gz")
    log:
        os.path.join("logs", "cpg", "{sample}.log")
    threads: config.get("threads", 4)
    resources:
        mem_mb = config.get("memory_mb", 8000)
    params:
        min_coverage = config.get("min_coverage", 5),
        output_prefix = "{sample}",
        modes = config.get("pileup", "count")
    shell:
        """
        mkdir -p {CpG_DIR}
        aligned_bam_to_cpg_scores \
            --bam {input.bam} \
            --output-prefix {params.output_prefix} \
            --pileup-mode {params.modes} \
            --threads {threads} \
            --min-coverage {params.min_coverage} \
            2> {log} && \
        mv {wildcards.sample}.combined.bed.gz {output.cpg}
        mv {wildcards.sample}.combined.bed.gz.tbi {output.cpg}.tbi
        mv {wildcards.sample}* {CpG_DIR}
        """

# Rule to mask CpG scores with BED file
rule repeat_mask_cpg_scores:
    input:
        cpg = os.path.join(CpG_DIR, "{sample}.cpg.bed.gz"),
        bed = BED_FILE
    output:
        masked = os.path.join(MASKED_DIR, "{sample}.masked.cpg.bed") ## naming
    log:
        os.path.join("logs", "masking", "{sample}.log")
    shell:
        """
        mkdir -p {MASKED_DIR}

        zcat {input.cpg} | grep "^#" > {output.masked}.header
        bedtools intersect -v \
            -a <( zcat {input.cpg}| grep -v "^#") \
            -b <( zcat {input.bed} |grep -v "^#"|cut -f 1-4) \
            > {output.masked}.content \
            2> {log}

        cat {output.masked}.header {output.masked}.content | grep -v "^##" > {output.masked}
        
        rm {output.masked}.header {output.masked}.content

        """

rule convert_to_dss:
    input:
        masked_cpg = os.path.join(MASKED_DIR, "{sample}.masked.cpg.bed"),
        script = DSS_CONVERT_SCRIPT
    output:
        dss = os.path.join(DSS_DIR, "{sample}.dss.bed")
    log:
        os.path.join("logs", "dss_conversion", "{sample}.log")
    resources:
        mem_mb = config.get("r_memory_mb", 8000)
    params:
        min_coverage = config.get("min_coverage", 5)        
    shell:
        """
        mkdir -p {DSS_DIR}
        Rscript {input.script} \
            --input {input.masked_cpg} \
            --output {output.dss} \
            --sample {wildcards.sample} \
            -c {params.min_coverage} -v \
            2> {log}
        """


# Rule to run R analysis on each pair of masked CpG files
rule analyze_cpg_pair:
    input:
        dss1 = os.path.join(DSS_DIR, "{sample1}.dss.bed"),
        dss2 = os.path.join(DSS_DIR, "{sample2}.dss.bed"),
        script = R_SCRIPT
    output:
        result = os.path.join(ANALYSIS_DIR, "{sample1}_vs_{sample2}_pb_methylation_analysis_grid.png")
    log:
        os.path.join("logs", "analysis", "{sample1}_{sample2}.log")
    resources:
        mem_mb = config.get("r_memory_mb", 8000),
        threads = workflow.cores #config.get("threads")
    params:
        p_threshold = config.get("p_threshold", 0.05),
        delta_threshold = config.get("delta_threshold", 0.1),
        util = input_data["analysis_util"],
        bw = input_data["bigWigs"]
    shell:
        """
        mkdir -p {ANALYSIS_DIR}
        Rscript {input.script} \
            --input1 {input.dss1} \
            --input2 {input.dss2} \
            --sample1 {wildcards.sample1} \
            --sample2 {wildcards.sample2} \
            --outputDir {ANALYSIS_DIR} \
            --p_threshold {params.p_threshold} \
            --delta_threshold {params.delta_threshold} \
            --util {params.util} --bw {params.bw} --verbose \
            > {log} 2>&1
        
        # Move all pb_* files from R script output location to analysis dir
        mv {wildcards.sample1}_vs_{wildcards.sample2}_* {ANALYSIS_DIR}/ 2>/dev/null || true
        
        # Ensure the required output file exists
        touch {output.result}
        """