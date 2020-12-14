#!/bin/bash

myarr=($(awk '{ print $1 }' data/res_bolt_lav_disc.txt))

# for i in "${myarr[@]}"
# do
#   :
#   echo $i
# done

# zcat /gpfs01/bhcbio/projects/UK_Biobank/20190102_UK_Biobank_Imaging/Results/GWAS/GWAS_diastolic_BOLT/Results/bolt_radial_PDSR_repl.bgen.stats.gz | ( IFS='|'; grep -E "^(${myarr[*]})\>" ) > data/res_bolt_radial_pdsr_repl_disc_match.txt
cat /gpfs01/bhcbio/projects/UK_Biobank/20190102_UK_Biobank_Imaging/Results/GWAS/GWAS_diastolic_BOLT/Results/bolt_LAV_repl.bgen.stats | ( IFS='|'; grep -E "^(${myarr[*]})\>" ) > data/res_bolt_lav_repl_disc_match.txt
