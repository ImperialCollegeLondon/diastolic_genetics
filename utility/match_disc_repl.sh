#!/bin/bash

myarr=($(awk '{ print $1 }' data/res_bolt_radial_pdsr_full.txt))

# for i in "${myarr[@]}"
# do
#   :
#   echo $i
# done

zcat /gpfs01/bhcbio/projects/UK_Biobank/20190102_UK_Biobank_Imaging/Results/GWAS/GWAS_diastolic_BOLT/Results/bolt_radial_PDSR_repl.bgen.stats.gz | ( IFS='|'; grep -E "^(${myarr[*]})\>" ) > data/res_bolt_radial_pdsr_repl_match.txt
zcat /gpfs01/bhcbio/projects/UK_Biobank/20190102_UK_Biobank_Imaging/Results/GWAS/GWAS_diastolic_BOLT/Results/bolt_radial_PDSR_disc.bgen.stats.gz | ( IFS='|'; grep -E "^(${myarr[*]})\>" ) > data/res_bolt_radial_pdsr_disc_match.txt
