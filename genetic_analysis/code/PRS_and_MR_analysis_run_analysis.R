## Assume that the GWAS results from radial PDSR are stored in a data frame called radialTab which contains at least the
## following columns:
## SNP: RSIDs (character)
## CHR: number chromosome of SNP
## MAF: minor allele frequency
## ESTIMATE: point estimate of effect size
## SE: standard error
## P: P-value of BOLT-LMM

## In the same way: long PDSR is stored in a data frame called longTab and LAVmax (indexed to BSA) is stored
## in a data frame called lavTab.

## bgen_file_for_chrom: list with path to bgen files from UK Biobank (e.g., _004_ukb_imp_chr1_v3.bgen)
## sample_file: path to sample file (e.g., ukb40616_imp_chr17_v3_s487297.sample)
## tmpdir/plink.tmpdir: path to directories for storing temporary results

library(tidyverse)
library(ggplot2)
library(data.table)


###################### Step 1: perform clumping to identify independent signals #############################

#### radial PDSR

## obtain list of SNPs via clumping
radial <- gwas_clumping(score_input=radialTab, score_threshold=10^-6, clump_r2=0.1, tmpdir=tmpdir, 
                        plink.tmpdir=plink.tmpdir, snpset=NA,
                        minMaf=0.01, bgen_file_for_chrom=bgen_file_for_chrom, sample_file=sample_file)
## remove list structure and exclude empty entries (are generated if no SNP is selected for a specific CHR)
radial <- unlist(radial)
radialSNPs <- radial[which(radial!="")]

#### long PDSR

## obtain list of SNPs via clumping
long <- gwas_clumping(score_input=longTab, score_threshold=10^-6, clump_r2=0.1, tmpdir=tmpdir, 
                      plink.tmpdir=plink.tmpdir, snpset=NA,
                      minMaf=0.01, bgen_file_for_chrom=bgen_file_for_chrom, sample_file=sample_file)
## remove list structure and exclude empty entries (are generated if no SNP is selected for a specific CHR)
long <- unlist(long)
longSNPs <- long[which(long!="")]

#### LAVmax (indexed to BSA)

## obtain list of SNPs via clumping
lav <- gwas_clumping(score_input=lavTab, score_threshold=10^-6, clump_r2=0.1, tmpdir=tmpdir, 
                     plink.tmpdir=plink.tmpdir, snpset=NA,
                     minMaf=0.01, bgen_file_for_chrom=bgen_file_for_chrom, sample_file=sample_file)
## remove list structure and exclude empty entries (are generated if no SNP is selected for a specific CHR)
lav <- unlist(lav)
lavSNPs <- lav[which(lav!="")]


###################### Step 2: compute and evaluate genetic score #############################################

## Assume that:
## - geno is a data frame containing the genetic data (columns SID plus additional columns for all genotypes, column
##   names of geno_data match entries in radialSNPs). Data is included for all QCed European subjects.
## - pheno is a data frame containing the measured values for radial PDSR/longitudinal PDSR/LAVmax, a column with the SIDs plus
##   all covariates necessary for the model
## - covariates is a data frame containing age, sex (first column SID, afterwards sex and age)
## - imaging is a vector containing SIDs of all QCed European samples who were part of the imaging study (releases 1-3)
## - non-imaging is a vector containing SIDs of all QCed European samples who were NOT part of the imaging study (release (1-3))
## - pca is a data frame containing SID plus the first ten principal components
## - phenotypes_bin is a data frame containing all binary traits which we want to consider for the PheWAS
## - phenotypes_quant is a data frame containing all quantitative traits which we want to consider for the PheWAS
## - meta_quant/meta_bin: dataframe containing two columns, first column called 'endpoint' with names of phenotypes as used in the 
##  code (e.g., standing_height), second column is the display name (e.g., Standing height) and called 'description'

#### radial PDSR

## compute genetic score
radialScore <- computeGeneticScoreForSnps(radialSNPs, trait="radial_PDSR", geno_data = geno, 
                                                  pheno_data=pheno)

## split into imaging/non-imaging subset
radialScoreMRI <- radialScore$Score %>% dplyr::filter(SID %in% imaging) %>% dplyr::select(SID, Score)
radialScoreNonMRI <- radialScore$Score %>% dplyr::filter((SID %in% non_imaging)) %>% dplyr::select(SID, Score)

## calculate explained variability with score
explained_radial <- dplyr::inner_join(radialScoreMRI, pheno, by="SID")
summary(lm(radial_PDSR~Score, data=explained_radial))
# explained variability = R^2 of this model

## run PheWAS
QuantRadial <- getQuantPheWASResults(radialScoreNonMRI, covariates, pca, phenotypes=phenotypes_quant, 
                                 scale_phenotypes = TRUE)
BinRadial <- getPheWASResults(radialScoreNonMRI, covariates, pca, phenotypes=phenotypes_bin, cores = 1)

## combine results and apply multiplicity correction
Radial <- rbind(QuantRadial %>% dplyr::select(phenotype, estimate, conf.low, conf.high, p.value, statistic), 
              BinRadial %>% dplyr::select(phenotype, estimate, conf.low, conf.high, p.value, statistic)) 
RadialMT <- Radial %>% dplyr::mutate(p.adjusted=p.adjust(p.value, method="fdr")) %>% dplyr::filter(p.adjusted < 0.05)

#### long PDSR

## compute genetic score
longScore <- computeGeneticScoreForSnps(longSNPs, trait="long_PDSR", geno_data = geno, 
                                                  pheno_data=pheno)

## split into imaging/non-imaging subset
longScoreMRI <- longScore$Score %>% dplyr::filter(SID %in% imaging) %>% dplyr::select(SID, Score)
longScoreNonMRI <- longScore$Score %>% dplyr::filter((SID %in% non_imaging)) %>% dplyr::select(SID, Score)

## calculate explained variability with score
explained_long <- dplyr::inner_join(longScoreMRI, pheno, by="SID")
summary(lm(long_PDSR~Score, data=explained_long))
# explained variability = R^2 of this model

## run PheWAS
QuantLong <- getQuantPheWASResults(longScoreNonMRI, covariates, pca, phenotypes=phenotypes_quant, 
                                     scale_phenotypes = FALSE)
BinLong <- getPheWASResults(longScoreNonMRI, covariates, pca, phenotypes=phenotypes_bin, cores = 1)

## combine results and apply multiplicity correction
Long <- rbind(QuantLong %>% dplyr::select(phenotype, estimate, conf.low, conf.high, p.value, statistic), 
             BinLong %>% dplyr::select(phenotype, estimate, conf.low, conf.high, p.value, statistic)) 
LongMT <- Long %>% dplyr::mutate(p.adjusted=p.adjust(p.value, method="fdr")) %>% dplyr::filter(p.adjusted < 0.05)

#### LAVmax (indexed to BSA)

## compute genetic score
lavScore <- computeGeneticScoreForSnps(lavSNPs, trait="lav", geno_data = geno, 
                                                pheno_data=pheno)

## split into imaging/non-imaging subset
lavScoreMRI <- lavScore$Score %>% dplyr::filter(SID %in% imaging) %>% dplyr::select(SID, Score)
lavScoreNonMRI <- lavScore$Score %>% dplyr::filter((SID %in% non_imaging)) %>% dplyr::select(SID, Score)

## calculate explained variability with score
explained_lav <- dplyr::inner_join(lavScoreMRI, pheno, by="SID")
summary(lm(lav~Score, data=explained_lav))
# explained variability = R^2 of this model


## run PheWAS
QuantLav <- getQuantPheWASResults(lavScoreNonMRI, covariates, pca, phenotypes=phenotypes_quant, 
                                   scale_phenotypes = FALSE)
BinLav <- getPheWASResults(lavScoreNonMRI, covariates, pca, phenotypes=phenotypes_bin, cores = 1)

## combine results and apply multiplicity correction
Lav <- rbind(QuantLav %>% dplyr::select(phenotype, estimate, conf.low, conf.high, p.value, statistic), 
             BinLav %>% dplyr::select(phenotype, estimate, conf.low, conf.high, p.value, statistic)) 
LavMT <- Lav %>% dplyr::mutate(p.adjusted=p.adjust(p.value, method="fdr")) %>% dplyr::filter(p.adjusted < 0.05)

#### do cross-validation: which results stay if any of the SNPs is excluded from the score at P<0.05
## radial PDSR
resRadial <- list()
for (i in 1:length(RadialMT$phenotype)) {
  resRadial[[i]] <- matrix(data=NA, ncol=7, nrow=length(radialSNPs))
  for (j in 1:length(radialSNPs)) {
    resRadial[[i]][j,1:6] <- unlist(leave_one_out_cv(snps_include=radialSNPs[-j], trait=RadialMT$phenotype[i], imaging_trait="radial_PDSR", geno=geno,
                                                  covariates=covariates, pca=pca, phenotypes_quant=phenotypes_quant,
                                                  phenotypes_bin=phenotypes_bin, phenotypes_score=pheno, imaging=imaging))
    resRadial[[i]][j,7] <- radialSNPs[j]
  }
}
names(resRadial) <- RadialMT$phenotype

finalRadial <- c()
for (i in 1:length(resRadial)) {
  if (all(p.adjust(as.numeric(resRadial[[i]][,5], method="fdr"))  < 0.05)) {
    finalRadial <- c(finalRadial, resRadial[[i]][1,1])
  }
}

## long PDSR
resLong <- list()
for (i in 1:length(LongMT$phenotype)) {
  resLong[[i]] <- matrix(data=NA, ncol=7, nrow=length(longSNPs))
  for (j in 1:length(longSNPs)) {
    resLong[[i]][j,1:6] <- unlist(leave_one_out_cv(snps_include=longSNPs[-j], trait=LongMT$phenotype[i], imaging_trait="long_PDSR", geno=geno,
                                                     covariates=covariates, pca=pca, phenotypes_quant=phenotypes_quant,
                                                     phenotypes_bin=phenotypes_bin, phenotypes_score=pheno, imaging=imaging))
    resLong[[i]][j,7] <- longSNPs[j]
  }
}
names(resLong) <- LongMT$phenotype

finalLong <- c()
for (i in 1:length(resLong)) {
  if (all(p.adjust(as.numeric(resLong[[i]][,5], method="fdr"))  < 0.05)) {
    finalLong <- c(finalLong, resLong[[i]][1,1])
  }
}

## LAVmax(i)
resLav <- list()
for (i in 1:length(LavMT$phenotype)) {
  resLav[[i]] <- matrix(data=NA, ncol=7, nrow=length(lavSNPs))
  for (j in 1:length(lavSNPs)) {
    resLav[[i]][j,1:6] <- unlist(leave_one_out_cv(snps_include=lavSNPs[-j], trait=LavMT$phenotype[i], imaging_trait="lav", geno=geno,
                                                  covariates=covariates, pca=pca, phenotypes_quant=phenotypes_quant,
                                                  phenotypes_bin=phenotypes_bin, phenotypes_score=pheno, imaging=imaging))
    resLav[[i]][j,7] <- lavSNPs[j]
  }
}
names(resLav) <- LavMT$phenotype

finalLav <- c()
for (i in 1:length(resLav)) {
  if (all(p.adjust(as.numeric(resLav[[i]][,5], method="fdr"))  < 0.05)) {
    finalLav <- c(finalLav, resLav[[i]][1,1])
  }
}

## get final list in the same way for long PDSR (finalLav) and radial PDSR (finalRadial)
traitsSel <- union(union(finalLong, finalRadial), finalLav)
traitsSelQuant <- traitsSel[which(traitsSel %in% names(phenotypes_quant))]
traitsSelBin <- traitsSel[which(traitsSel %in% names(phenotypes_bin))]

#### display results
## quantitative hits
figure_data <- bind_rows(Lav %>% dplyr::filter(phenotype %in% traitsSelQuant),
                         Radial %>% dplyr::filter(phenotype %in% traitsSelQuant),
                         Long %>% dplyr::filter(phenotype %in% traitsSelQuant)) %>% 
  dplyr::mutate(prs_trait=rep(c("LAVmax_i", "PDSR_rr","PDSR_ll"), each=length(traitsSelQuant))) %>% dplyr::mutate(p.adjusted=p.adjust(p.value,
                                                                                                                  method="fdr"))
figure_data <- dplyr::left_join(figure_data %>% dplyr::rename(endpoint=phenotype), meta_quant, by="endpoint")

p1 <- plot_PRS(figure_data, title="Associations to quantitative traits", xlim=c(-5,5), xlab="Change in trait per 1 SD increase  in diastolic function trait", zero=0, vjust=1.5, ylim=3, topmargin=0.05)

## binary hits
figure_data_bin <- bind_rows(Lav %>% dplyr::filter(phenotype %in% traitsSelBin),
                             Radial  %>% dplyr::filter(phenotype %in% traitsSelBin),
                             Long  %>% dplyr::filter(phenotype %in% traitsSelBin)) %>% 
  dplyr::mutate(prs_trait=rep(c("LAVmax_i", "PDSR_rr","PDSR_ll"), each=length(traitsSelBin)))  %>% dplyr::mutate(p.adjusted=p.adjust(p.value,
                                                                                                                                     method="fdr"))
figure_data_bin <- dplyr::left_join(figure_data_bin %>% dplyr::rename(endpoint=phenotype), meta_bin, by="endpoint")

p2 <- plot_PRS(figure_data_bin, title="Associations to binary traits", xlim=c(-5,5), 
                  xlab="Log(Odds ratio) per 1 SD increase in diastolic function trait", zero=0, vjust=0.5, topmargin=0.1, ylim=2)

cowplot::plot_grid(NULL, p1, NULL, p2, ncol=1, labels = c("A", "", "", "B"), rel_heights = c(-0.1, 4.3,0,1), vjust = 1)


###################### Step 3: MR analysis ###################################################################

## Assumption: GWAS summary stats for diastolic and non-diastolic traits are stored within a data frame with 
## the following columns:
## * P: containing the P-values
## * SNP: containing the RSIDs
## * CHR: chromosome
## * MAF: minor allele frequency
## * ESTIMATE: point estimates
## * SE: standard error
##
## We assume that the alleles of the GWAS results for diastolic and non-diastolic traits are matched, i.e., the
## directionality of the alleles are aligned.

## example run
mr_analysis(diastolic_gwas=diastolic_gwas, nondiastolic_gwas=nondiastolic_gwas, name_diastolic="PDSR_ll", name_nondiastolic="diastolic_blood_pressure", 
            score_threshold=10^-6, clump_kb=1000, clump_r2=0.1, tmpdir=tmpdir, bgen_file_for_chrom=bgen_file_for_chrom, sample_file=sample_file, minMaf=0)

####################### Step 4: comparison to observational dataset #############################################

## We assume, for example, the following data structure for the data frame'obs_data_in':
## - Heart_Failure.Event: binary variable where 1 reflects a case and 0 reflects a control
## - radial_PDSR: value for PDSR_{rr}
## - Age: age at CMR visit
## - Sex coded as binary variables
## - Diabetes_mellitus: 1 if a subject had a diabetes diagnosis prior to the CMR visit, 1 otherwise
## - dbp: diastolic blood pressure measurement
## - bmi: measurement for bmi

## example for Heart failure vs. radial PDSR
confint.default(glm(Heart_Failure.Event~radial_PDSR+Age+Sex+Diabetes_mellitus+dbp+bmi, 
                              data=obs_data_in, family="binomial"), level = 0.95)[2,] 
