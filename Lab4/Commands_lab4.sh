#!/bin/sh

# STEP 1: Mutect 2
# source /usr/local/GenomicsTechnologies/setup_env.sh
gatk CreateSequenceDictionary -R GRCh38.primary_assembly.genome.fa

samtools faidx GRCh38.primary_assembly.genome.fa
gatk Mutect2 --reference GRCh38.primary_assembly.genome.fa --input RPE-MITO_sorted.bam --output RPE-MITO.unfiltered.vcf.gz
gatk Mutect2 --reference GRCh38.primary_assembly.genome.fa --input RPE-ROI_sorted.bam --output RPE-ROI.unfiltered.vcf.gz

# STEP 2: filter variants
gatk FilterMutectCalls --reference GRCh38.primary_assembly.genome.fa --variant RPE-ROI.unfiltered.vcf.gz --output RPE-ROI.filtered.vcf.gz

gatk FilterMutectCalls --reference GRCh38.primary_assembly.genome.fa --variant RPE-MITO.unfiltered.vcf.gz --output RPE-MITO.filtered.vcf.gz
