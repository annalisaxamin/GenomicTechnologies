#!/bin/sh

# STEP 1: filter VCF to keep only chrM
# source /usr/local/GenomicsTechnologies/setup_env.sh
bcftools view -r chrM RPE-ROI.filtered.vcf.gz > RPE-ROI.filtered.chrM.vcf.gz

bcftools view -r chrM RPE-MITO.filtered.vcf.gz > RPE-MITO.filtered.chrM.vcf.gz

# STEP 2: VEP
# use interface https://www.ensembl.org/Tools/VEP