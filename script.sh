#!/bin/bash

#-------------------------------------------------------------------
# Script to call somatic variants from tumor and normal sample pairs
#-------------------------------------------------------------------


# *****Download sample datasets*****

# Normal sample 
wget -P /home/aishamehmood/Project_Mutect2/reads https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data_somatic/HG008/Liss_lab/NYGC_Illumina-WGS_20231023/HG008-N-D_CGGACAAC-AATCCGGA_H3LLJDSXC_L001_001.R1.fastq.gz

wget -P /home/aishamehmood/Project_Mutect2/reads https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data_somatic/HG008/Liss_lab/NYGC_Illumina-WGS_20231023/HG008-N-D_CGGACAAC-AATCCGGA_H3LLJDSXC_L001_001.R2.fastq.gz

# Tumor sample
wget -P /home/aishamehmood/Project_Mutect2/reads https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data_somatic/HG008/Liss_lab/NYGC_Illumina-WGS_20231023/HG008-T_TTCCTGTT-AAGATACT_HJVY2DSX7_L001_001.R1.fastq.gz

wget -P /home/aishamehmood/Project_Mutect2/reads https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data_somatic/HG008/Liss_lab/NYGC_Illumina-WGS_20231023/HG008-T_TTCCTGTT-AAGATACT_HJVY2DSX7_L001_001.R2.fastq.gz

# Subset 1,000,000 reads from each FASTQ file (this step is done to reduce data size due to system limitations) 
cd ../reads/

seqtk sample -s100 HG008-N-D_CGGACAAC-AATCCGGA_H3LLJDSXC_L001_001.R1.fastq.gz 1000000 | gzip > HG008-N-D_CGGACAAC-AATCCGGA_subset_H3LLJDSXC_L001_001.R1.fastq.gz
seqtk sample -s100 HG008-N-D_CGGACAAC-AATCCGGA_H3LLJDSXC_L001_001.R2.fastq.gz 1000000 | gzip > HG008-N-D_CGGACAAC-AATCCGGA_subset_H3LLJDSXC_L001_001.R2.fastq.gz
seqtk sample -s100 HG008-T_TTCCTGTT-AAGATACT_HJVY2DSX7_L001_001.R1.fastq.gz 1000000 | gzip > HG008-T_TTCCTGTT-AAGATACT_subset_HJVY2DSX7_L001_001.R1.fastq.gz
seqtk sample -s100 HG008-T_TTCCTGTT-AAGATACT_HJVY2DSX7_L001_001.R2.fastq.gz 1000000 | gzip > HG008-T_TTCCTGTT-AAGATACT_subset_HJVY2DSX7_L001_001.R2.fastq.gz
cd

# *****Download additional files*****

# Reference files
wget -P /home/aishamehmood/Project_Mutect2/additional_files/hg38/ https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
gunzip /home/aishamehmood/Project_Mutect2/additional_files/hg38/hg38.fa.gz

# Index reference file (required by Mutect2)
samtools faidx /home/aishamehmood/Project_Mutect2/additional_files/hg38/hg38.fa

# Create a sequence dictionary for the reference file (required by Mutect2)
gatk CreateSequenceDictionary R=/home/aishamehmood/Project_Mutect2/additional_files/hg38/hg38.fa O=/home/aishamehmood/Project_Mutect2/additional_files/hg38/hg38.dict

# Download known sites files from GATK resource bundle for BQSR
wget -P /home/aishamehmood/Project_Mutect2/additional_files/hg38/ https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf

wget -P /home/aishamehmood/Project_Mutect2/additional_files/hg38/ https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.idx

# *****Set directories path*****

ref=/home/aishamehmood/Project_Mutect2/additional_files/hg38/hg38.fa
known_sites=/home/aishamehmood/Project_Mutect2/additional_files/hg38/Homo_sapiens_assembly38.dbsnp138.vcf
project_dir=/home/aishamehmood/Project_Mutect2
aligned_reads=$project_dir/aligned
reads=$project_dir/reads
results=$project_dir/results


# -----------------------------------------
# Step 1: Pre-processing
# -----------------------------------------

# 1.1. Quality control - Run FastQC
for file in ${reads}/*.gz; do            #Loop through all .gz files in the reads folder and run FastQC
    fastqc $file -o ${results}/
done
# No trimming required, quality looks okay.



