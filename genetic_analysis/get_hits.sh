#!/bin/bash

# zcat /gpfs01/bhcbio/projects/UK_Biobank/20190102_UK_Biobank_Imaging/Results/GWAS/GWAS_diastolic_BOLT/Results/bolt_LAV_full.bgen.stats | awk '$16<=5e-8' - > /gpfs01/home/glkwj/projects/diastolic_genetics/genetic_analysis/data/res_bolt_lav_full.txt
# zcat /gpfs01/bhcbio/projects/UK_Biobank/20190102_UK_Biobank_Imaging/Results/GWAS/GWAS_diastolic_BOLT/Results/bolt_long_PDSR_full.bgen.stats | awk '$16<=5e-8' - > /gpfs01/home/glkwj/projects/diastolic_genetics/genetic_analysis/data/res_bolt_long_pdsr_full.txt
# zcat /gpfs01/bhcbio/projects/UK_Biobank/20190102_UK_Biobank_Imaging/Results/GWAS/GWAS_diastolic_BOLT/Results/bolt_radial_PDSR_full.bgen.stats | awk '$16<=5e-8' - > /gpfs01/home/glkwj/projects/diastolic_genetics/genetic_analysis/data/res_bolt_radial_pdsr_full.txt

zcat /gpfs01/bhcbio/projects/UK_Biobank/20190102_UK_Biobank_Imaging/Results/GWAS/GWAS_diastolic_BOLT/Results/bolt_LAV_disc.bgen.stats | awk '$16<=5e-8' - > /gpfs01/home/glkwj/projects/diastolic_genetics/genetic_analysis/data/res_bolt_lav_disc.txt
zcat /gpfs01/bhcbio/projects/UK_Biobank/20190102_UK_Biobank_Imaging/Results/GWAS/GWAS_diastolic_BOLT/Results/bolt_long_PDSR_disc.bgen.stats | awk '$16<=5e-8' - > /gpfs01/home/glkwj/projects/diastolic_genetics/genetic_analysis/data/res_bolt_long_pdsr_disc.txt
zcat /gpfs01/bhcbio/projects/UK_Biobank/20190102_UK_Biobank_Imaging/Results/GWAS/GWAS_diastolic_BOLT/Results/bolt_radial_PDSR_disc.bgen.stats | awk '$16<=5e-8' - > /gpfs01/home/glkwj/projects/diastolic_genetics/genetic_analysis/data/res_bolt_radial_pdsr_disc.txt