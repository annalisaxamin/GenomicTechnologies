#!/bin/sh

# STEP 1: building the reference
bwa index -p mitoIndex GRCh38.primary_assembly.genome.fa.gz

# STEP 2: aligning the reads
bwa mem -t 2 -R "@RG\tID:RPE-MITO-NXT\tSM:RPE-MITO-NXT" mitoIndex RPE-MITO_paired_R1.fastq.gz RPE-MITO_paired_R2.fastq.gz > RPE-MITO-NXT.sam

samtools view -b RPE-MITO-NXT.sam > RPE-MITO-NXT.bam


bwa mem -t 2 -R "@RG\tID:RPE-ROI\tSM:RPE-ROI" mitoIndex RPE-ROI_paired_R1.fastq.gz RPE-ROI_paired_R2.fastq.gz > RPE-ROI.sam

samtools view -b RPE-ROI.sam > RPE-ROI.bam

# STEP 3: BAM indexing
samtools sort RPE-MITO-NXT.bam > RPE-MITO_sorted.bam
samtools index RPE-MITO_sorted.bam

samtools sort RPE-ROI.bam > RPE-ROI_sorted.bam
samtools index RPE-ROI_sorted.bam


# STEP 4: visualizing the alignment 
# use IGV