# Somatic_Variant_Calling_Mutect2
Somatic Variant Calling with Mutect2: A Beginner-Friendly Workflow Following GATK Best Practices

## Files and Folders
- **README.md:** Project documentation.
- **script.sh:** Shell script containing the commands for each step in the workflow.
- **datalinks.txt:** Contains URLs for datasets used in the project.
- **results/:** Contains results of variant calling workflow.

## Overview
This project performs somatic variant calling on **Illumina WGS data** using **Mutect2** from **GATK**. The analysis follows the best practices as outlined by GATK, with pre-processing, variant calling, and filtering steps. The datasets used in this project are publicly available from the Genome in a Bottle (GIAB) consortium and include samples from HG008. The goal is to identify somatic mutations in cancer samples through comprehensive variant analysis.

## Installation and Setup
### System requirements
- **UNIX OS**
- **Java 17.0.12**
### Tools required
- **GATK** (Genome Analysis Toolkit) 
- **BWA** (Burrows-Wheeler Aligner) 
- **SAMtools** 
- **FASTQC** 
- **MULTIQC** 

## Data
### Key Details:
-**Tumor Sample:** HG008-T
-	**Source:** Pancreatic Ductal Adenocarcinoma Cell Line (PDAC)
-	**Sequencing Platform:** Illumina NovaSeq 6000
-	**Read Length:** Paired-end 150bp

-**Matched Normal Sample:** HG008-N
-	**Source:** Duodenal Solid Tissue Normal (STN)
-	**Purpose:** Used as a control for somatic variant calling to distinguish somatic mutations from germline variants
-	**Sequencing Platform:** Illumina NovaSeq 6000
-	**Read Length:** Paired-end 150bp

## Workflow
The project follows a step-by-step process:  
### **1. Pre-processing**  
The raw sequencing reads are processed through several steps to ensure data quality and proper alignment before variant calling.  

- **Quality Control using FastQC**  
  The raw sequencing reads (*FASTQ files*) are assessed for quality metrics such as base quality scores, GC content, adapter contamination, and sequence duplication levels. This step helps to identify any issues in the data before proceeding.  

- **Alignment using BWA-MEM**  
  The high-quality reads are aligned to the reference genome (hg38) using the **BWA-MEM** algorithm. This step maps the sequencing reads to their corresponding genomic locations, producing a **SAM file** (Sequence Alignment/Map).  

- **Mark Duplicates and Sort using GATK**  
  Duplicate reads, which may have arisen from PCR amplification during library preparation, are identified and flagged using **GATK MarkDuplicates**. Sorting the reads ensure they are arranged in genomic order, improving processing efficiency in downstream analyses.  

- **Base Quality Recalibration (BQSR) using GATK**  
  Systematic errors in base quality scores introduced by sequencing instruments are corrected using **GATK BaseRecalibrator**. This step relies on known variation databases ( **dbSNP**) to refine quality scores, increasing variant calling accuracy.  

- **Collect Alignment & Insert Size Metrics using GATK**  
  Alignment statistics, such as read mapping percentages, coverage depth, and insert size distribution, are gathered using **GATK CollectAlignmentSummaryMetrics** and **CollectInsertSizeMetrics**. These metrics provide insights into sequencing performance and alignment quality.  


### **2. Variant Calling**

- **Call Variants using Gatk Mutect2**
  Somatic variants are identified using Mutect2, with the analysis limited to 10 chromosomes to reduce computational
  demands. The workflow remains the same for all chromosomes analysis, and users can modify the -L option accordingly. The
  output is an unfiltered VCF file containing raw variant calls. The --f1r2-tar-gz option is used
  during this step to collect strand orientation bias data, which is later required for filtering artifacts caused by
  sequencing errors.

- **Estimating Cross Sample Contamination**
  This step helps estimate the fraction of cross-sample contamination. The **GATK GetPileupSummaries** tool summarizes read support for known population variants. The **GATK
  CalculateContamination** tool estimates the contamination fraction. The output, HG008_pair_calculatecontamination.table, is used in **FilterMutectCalls**.

- **Read Orientation Artifacts Estimation**
  To filter out strand bias artifacts, a read orientation model is created using data from Mutect2. This helps identify sequencing errors caused by strand-specific bias. The
  model is generated from F1R2 counts collected during the variant calling step. It is later used to refine variant filtering and improve accuracy.

- **Filtering Variants**
  Raw variant calls are filtered using **FilterMutectCalls** to remove false positives based on multiple factors, including:
    - Sequencing artifacts (using the read orientation model). 
    - Low-quality variants and insufficient evidence supporting somatic mutations.
    - Cross-sample contamination (using estimated contamination levels).
  The final output is a high-confidence set of somatic variants ready for further annotation and interpretation.
   



