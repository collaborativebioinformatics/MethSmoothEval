# MethSmoothEval (Group 6 - 2025 BCM HGSC SVcrying Hackathon)
## Group members
[Yilei Fu](https://github.com/Fu-Yilei),
[Robel Dagnew](https://github.com/redndgreen8), 
[Jon Moller](https://github.com/molleraj),
[Rohit Kolora](https://github.com/evolgen),
[Alex Leonard](https://github.com/ASLeonard),
[Halimat Atanda](https://github.com/chisomgold), 
[Arijita Sarkar](https://github.com/arijita88),
[Senthilkumar Kailasam](https://github.com/dksenthil),
[Zahra Seyfollahi](https://github.com/ZSeyfollahi),
Nilabja Bhattacharjee,
[Sarah Eger](https://github.com/saraheger),
[Michael Olufemi](https://github.com/michael-olufemi),
[Ben Braun](https://github.com/bennibraun)

## Introduction
While the gold standard for DNA methylation detection is bisulfite sequencing, the advent of long-read sequencing has made it possible to detect DNA methylation and additional modifications in native nucleic acid. The goal of this project is to evaluating effects of smoothing on long read methylation calling.

### Key Takeaway
**Use smoothing with caution**, especially where single CpG effect have biological relevance like tissue-specific or cell-type specific methylation analyses. Some specific examples of "lost" CpG signals are highlighted in the [slides](https://docs.google.com/presentation/d/1xntJzpwD1QStqhe31eTOESBfSB2EzlydXQ2MT2NmyDY/edit?usp=sharing).

Pipeline/flowchart of MethSmoothEval
![Flowchart](img/flowchart.svg)


## Method
This project investigates whether smoothing of DNA methylation signals—historically essential in short-read bisulfite sequencing (WGBS)—is still necessary in the era of high-accuracy long-read sequencing (ONT, PacBio).  

We aimed to benchmark **single-CpG vs smoothed-CpG analysis** to assess:  
- When smoothing improves detection of differentially methylated CpGs/regions (DMCpGs/DMRs).  
- When smoothing introduces artifacts or masks true single-CpG variation.  
- Whether per-read accuracy and sequence context in ONT impact CpG methylation detection.  

1. **Data Sources**  
   - Short-read: SEQC2 EpiQC data (WGBS / EM-Seq / oxidative bisulfite sequencing).  
   - Long-read: ONT datasets (HG002, HG005).  

2. **Benchmarking**  
   - Run smoothing-based tools: `BSmooth`, `DSS`.  
   - Run long-read tools: `Modkit`, `pb-CpG-tools`.  
   - Compare per-CpG and DMR calls across technologies.  

3. **Key Analyses**  
   - CpG density profiling across hg38 (focus on chr22).  
   - Evaluate whether smoothing disrupts single-CpG interpretation (e.g., IGV inspection of low vs high DNAm CpGs in proximity).  
   - Identify counter-examples where ONT single-CpGs are biologically informative but missed in short-read smoothing pipelines.  
   - Assess tissue-specific effects (blood vs brain vs tumor).  

4. **Validation**  
   - Cross-check ONT read-level calls with EM-Seq / TrueMethyl.  
   - Ask: Are ONT-only CpG calls false positives, or do they reflect true biological signal?  

## Intended output
- **Catalog of CpG regions on hg38**  
  - Where single-CpG analysis is reliable.  
  - Where smoothing remains necessary.  
- **Benchmarking pipeline**  
  - Reproducible and fast for new long-read chemistries.  
  - Extensible to tissue- and cell-type–specific studies.  
- **Case Studies**  
  - Examples of ONT CpGs not recovered by EM-Seq smoothing.  
  - Markers distinguishing HG002 vs HG005.  
- **Follow-up Questions**  
  - Sequence context dependencies for CpG calling accuracy.  
  - Reliability of methylation calling in repetitive / TR regions.  
  - Tissue-specific CpG resolution requirements.  

## Results
![robel-density](img/22_cpg_density_100kb.png)
![EM-bias](img/EM-Seq-bias.png)

