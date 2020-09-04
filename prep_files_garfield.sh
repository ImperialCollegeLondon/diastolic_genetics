#!/bin/bash

# for *_full.bgen.stats.gz
# grep chr and save file
# awk pos, p-value and sort

# test at this point and see if results are obtainable with current package annotation

chrs=($(zcat /gpfs01/bhcbio/projects/UK_Biobank/20190102_UK_Biobank_Imaging/Results/GWAS/GWAS_diastolic_BOLT/Results/bolt_LAV_repl.bgen.stats.gz | tail -n +2 | cut -f 2 | sort | uniq))
echo "${chrs[*]}"

for chr in "${chrs[@]}"; do
  echo $chr
  zcat /gpfs01/bhcbio/projects/UK_Biobank/20190102_UK_Biobank_Imaging/Results/GWAS/GWAS_diastolic_BOLT/Results/bolt_LAV_repl.bgen.stats.gz | awk -v chr="$chr" '{ if($2 == chr) { print $3 " " $16}}' > garfield/garfield-data/pval/lav/chr$chr
done

