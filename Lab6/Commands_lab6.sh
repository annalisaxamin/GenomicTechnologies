#!/bin/sh

# Correct ALL VCF files in folder
files=*.vcf;

for f in $files; do
    echo $f;
    cat $f | python3 correct_vcf.py | bgzip > ${f%.vcf}.corrected.vcf.gz;
    tabix -p vcf ${f%.vcf}.corrected.vcf.gz;
done