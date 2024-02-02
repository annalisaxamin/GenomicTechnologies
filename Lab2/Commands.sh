#!/bin/sh
# STEP 1: fastqc
fastqc -o output_fastqc/ -f fastq *.fastq.gz

# STEP 2: trimming
java -jar /usr/local/Trimmomatic-0.39/trimmomatic-0.39.jar PE RPE-ROI_S6_L001_R1_001.fastq.gz RPE-ROI_S6_L001_R2_001.fastq.gz RPE-ROI_paired_R1.fastq.gz RPE-ROI_unpaired_R1.fastq.gz RPE-ROI_paired_R2.fastq.gz RPE-ROI_unpaired_R2.fastq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:True LEADING:3 TRAILING:30 MINLEN:36

java -jar /usr/local/Trimmomatic-0.39/trimmomatic-0.39.jar PE RPE-MITO-NXT_S1_L001_R1_001.fastq.gz RPE-MITO-NXT_S1_L001_R2_001.fastq.gz RPE-MITO_paired_R1.fastq.gz RPE-MITO_unpaired_R1.fastq.gz RPE-MITO_paired_R2.fastq.gz RPE-MITO_unpaired_R2.fastq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:True LEADING:3 TRAILING:30 MINLEN:36

# STEP 3: fastqc after trimming
fastqc -o output_fastqc_after_trimm/ -f fastq RPE-MITO_paired_R1.fastq.gz

fastqc -o output_fastqc_after_trimm/ -f fastq RPE-MITO_paired_R2.fastq.gz

fastqc -o output_fastqc_after_trimm/ -f fastq RPE-ROI_paired_R1.fastq.gz

fastqc -o output_fastqc_after_trimm/ -f fastq RPE-ROI_paired_R2.fastq.gz
