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
