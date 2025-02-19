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

# 1.2. Map to reference using BWA-MEM
# Index reference using BWA
bwa index ${ref}

# BWA alignment
bwa mem -t 4 -R "@RG\tID:HG008-N-D\tPL:ILLUMINA\tSM:HG008-N-D" ${ref} ${reads}/HG008-N-D_CGGACAAC-AATCCGGA_subset_H3LLJDSXC_L001_001.R1.fastq.gz ${reads}/HG008-N-D_CGGACAAC-AATCCGGA_subset_H3LLJDSXC_L001_001.R2.fastq.gz > ${aligned_reads}/HG008-N-D.paired.sam

bwa mem -t 4 -R "@RG\tID:HG008-T\tPL:ILLUMINA\tSM:HG008-T" ${ref} ${reads}/HG008-T_TTCCTGTT-AAGATACT_subset_HJVY2DSX7_L001_001.R1.fastq.gz ${reads}/HG008-T_TTCCTGTT-AAGATACT_subset_HJVY2DSX7_L001_001.R2.fastq.gz > ${aligned_reads}/HG008-T.paired.sam

# 1.3. Mark duplicates and sort - GATK
/home/aishamehmood/apps/gatk-4.6.0.0/gatk MarkDuplicatesSpark -I ${aligned_reads}/HG008-T.paired.sam -O ${aligned_reads}/HG008-T_sorted_dedup_reads.bam
/home/aishamehmood/apps/gatk-4.6.0.0/gatk MarkDuplicatesSpark -I ${aligned_reads}/HG008-N-D.paired.sam -O ${aligned_reads}/HG008-N-D_sorted_dedup_reads.bam

# 1.4. Base quality recalibration
# Build the model and then apply the model to adjust base quality scores
/home/aishamehmood/apps/gatk-4.6.0.0/gatk BaseRecalibrator -I ${aligned_reads}/HG008-N-D_sorted_dedup_reads.bam -R ${ref} --known-sites ${known_sites} -O ${aligned_reads}/HG008-N-D_recal_data.table

/home/aishamehmood/apps/gatk-4.6.0.0/gatk BaseRecalibrator -I ${aligned_reads}/HG008-T_sorted_dedup_reads.bam -R ${ref} --known-sites ${known_sites} -O ${aligned_reads}/HG008-T_recal_data.table

/home/aishamehmood/apps/gatk-4.6.0.0/gatk ApplyBQSR -I ${aligned_reads}/HG008-N-D_sorted_dedup_reads.bam -R ${ref} --bqsr-recal-file ${aligned_reads}/HG008-N-D_recal_data.table -O ${aligned_reads}/HG008-N-D_sorted_dedup_bqsr_reads.bam

/home/aishamehmood/apps/gatk-4.6.0.0/gatk ApplyBQSR -I ${aligned_reads}/HG008-T_sorted_dedup_reads.bam -R ${ref} --bqsr-recal-file ${aligned_reads}/HG008-T_recal_data.table -O ${aligned_reads}/HG008-T_sorted_dedup_bqsr_reads.bam

# 1.5. Collect Alignment & Insert Size Metrics
/home/aishamehmood/apps/gatk-4.6.0.0/gatk CollectAlignmentSummaryMetrics R=${ref} I=${aligned_reads}/HG008-N-D_sorted_dedup_bqsr_reads.bam  O=${aligned_reads}/HG008-N-D_alignment_metrics.txt

/home/aishamehmood/apps/gatk-4.6.0.0/gatk CollectAlignmentSummaryMetrics R=${ref} I=${aligned_reads}/HG008-T_sorted_dedup_bqsr_reads.bam  O=${aligned_reads}/HG008-T_alignment_metrics.txt

/home/aishamehmood/apps/gatk-4.6.0.0/gatk CollectInsertSizeMetrics INPUT=${aligned_reads}/HG008-N-D_sorted_dedup_bqsr_reads.bam OUTPUT=${aligned_reads}/HG008-N-D_insert_size_metrics.txt HISTOGRAM_FILE=${aligned_reads}/HG008-N-D_insert_size_histogram.pdf

/home/aishamehmood/apps/gatk-4.6.0.0/gatk CollectInsertSizeMetrics INPUT=${aligned_reads}/HG008-T_sorted_dedup_bqsr_reads.bam OUTPUT=${aligned_reads}/HG008-T_insert_size_metrics.txt HISTOGRAM_FILE=${aligned_reads}/HG008-T_insert_size_histogram.pdf


# *****Download Mutect2 additional files*****

# gnomAD
wget -P /home/aishamehmood/Project_Mutect2/additional_files/mutect2_additional_files https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/somatic-hg38/af-only-gnomad.hg38.vcf.gz                                            
wget -P /home/aishamehmood/Project_Mutect2/additional_files/mutect2_additional_files https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/somatic-hg38/af-only-gnomad.hg38.vcf.gz.tbi                                              

# PoN
wget -P /home/aishamehmood/Project_Mutect2/additional_files/mutect2_additional_files https://storage.googleapis.com/gatk-best-practices/somatic-hg38/1000g_pon.hg38.vcf.gz                                                
wget -P /home/aishamehmood/Project_Mutect2/additional_files/mutect2_additional_files https://storage.googleapis.com/gatk-best-practices/somatic-hg38/1000g_pon.hg38.vcf.gz.tbi 


# Set directories path
gatk_path=/home/aishamehmood/apps/gatk-4.6.0.0/gatk
mutect2_supporting_files=/home/aishamehmood/Project_Mutect2/additional_files/mutect2_additional_files


# -----------------------------------------
# Step 2: Variant Calling with Mutect2
# -----------------------------------------

# 2.1. Call variants using Gatk Mutect2
# Note: Only a subset of chromosomes was used for variant calling due to memory constraints (8GB RAM)
# This is for practice/demo purposes
${gatk_path} Mutect2 -R ${ref} \
     -I ${aligned_reads}/HG008-T_sorted_dedup_bqsr_reads.bam \
     -I ${aligned_reads}/HG008-N-D_sorted_dedup_bqsr_reads.bam \
     -tumor HG008-T \
     -normal HG008-N-D \
     --germline-resource ${mutect2_supporting_files}/af-only-gnomad.hg38.vcf.gz \
     --panel-of-normals ${mutect2_supporting_files}/1000g_pon.hg38.vcf.gz \
     -L chr1 -L chr3 -L chr5 -L chr7 -L chr9 -L chr11 -L chr13 -L chr15 -L chr17 -L chr19 \
     -O ${results}/HG008_somatic_variants_mutect2_selected_chroms.vcf.gz \
     --f1r2-tar-gz ${results}/HG008_f1r2_selected_chroms.tar.gz

# The below code snippet handles full WGS data for adequate computational resources
#${gatk_path} Mutect2 \
#    -R ${ref} \
#    -I ${aligned_reads}/HG008-T_sorted_dedup_bqsr_reads.bam \
#    -I ${aligned_reads}/HG008-N-D_sorted_dedup_bqsr_reads.bam \
#    -tumor HG008-T \
#    -normal HG008-N-D \
#    --germline-resource ${mutect2_supporting_files}/af-only-gnomad.hg38.vcf.gz \
#    --panel-of-normals ${mutect2_supporting_files}/1000g_pon.hg38.vcf.gz \
#    -O ${results}/HG008_somatic_variants_mutect2.vcf.gz \
#    --f1r2-tar-gz ${results}/HG008_f1r2_selected_chroms.tar.gz


# 2.2. Estimate cross sample contamination                                                                                                                                                                                                              # Download intervals list                                                                                                   
wget -P /home/aishamehmood/Project_Mutect2/additional_files/mutect2_additional_files https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/exome_calling_regions.v1.1.interval_list

# GetPileupSummaries
# Tumor
${gatk_path} GetPileupSummaries \
    --java-options '-Xmx3G' --tmp-dir ${project_dir}/tmp/ \
    -I ${aligned_reads}/HG008-T_sorted_dedup_bqsr_reads.bam \
    -V ${mutect2_supporting_files}/af-only-gnomad.hg38.vcf.gz \
    -L ${mutect2_supporting_files}/exome_calling_regions.v1.1.interval_list \
    -O ${results}/HG008_T_getpileupsummaries.table

# Normal
${gatk_path} GetPileupSummaries \
    --java-options '-Xmx3G' --tmp-dir ${project_dir}/tmp/ \
    -I ${aligned_reads}/HG008-N-D_sorted_dedup_bqsr_reads.bam  \
    -V ${mutect2_supporting_files}/af-only-gnomad.hg38.vcf.gz \
    -L ${mutect2_supporting_files}/exome_calling_regions.v1.1.interval_list \
    -O ${results}/HG008-N-D_getpileupsummaries.table

# Calculate contamination
${gatk_path} CalculateContamination \
    -I ${results}/HG008_T_getpileupsummaries.table \
    -matched ${results}/HG008-N-D_getpileupsummaries.table \
    -O ${results}/HG008_pair_calculatecontamination.table

                                              