~/projects/ldsc/munge_sumstats.py \
--sumstats /gpfs01/bhcbio/projects/UK_Biobank/20190102_UK_Biobank_Imaging/Results/GWAS/GWAS_diastolic_BOLT/Results/bolt_LAV_full.bgen.stats \
--out /gpfs01/bhcbio/projects/UK_Biobank/20190102_UK_Biobank_Imaging/Results/GWAS/GWAS_diastolic_BOLT/Results/bolt_LAV_full_ldsc \
--p P_BOLT_LMM \
--N 36541 \
--a1 Allele0 --a2 Allele1

~/projects/ldsc/ldsc.py \
--h2 /gpfs01/bhcbio/projects/UK_Biobank/20190102_UK_Biobank_Imaging/Results/GWAS/GWAS_diastolic_BOLT/Results/bolt_LAV_full_ldsc.sumstats.gz \
--ref-ld-chr data/eur_w_ld_chr/ \
--w-ld-chr data/eur_w_ld_chr/ \
--out /gpfs01/bhcbio/projects/UK_Biobank/20190102_UK_Biobank_Imaging/Results/GWAS/GWAS_diastolic_BOLT/Results/bolt_LAV_full_ldsc_h2
