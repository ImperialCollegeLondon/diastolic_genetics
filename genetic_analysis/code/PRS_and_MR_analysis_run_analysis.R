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
##

## del start
radial <- read.table("/gpfs01/bhcbio/projects/UK_Biobank/20190102_UK_Biobank_Imaging/Results/GWAS/GWAS_diastolic_BOLT/Results/bolt_radial_PDSR_full.bgen.stats", header=TRUE)
radialTab <- radial %>% dplyr::rename(P=P_BOLT_LMM, ESTIMATE=BETA)
radialTab$MAF <- apply(cbind(radialTab$A1FREQ, 1-radialTab$A1FREQ), 1, min)
rm(radial)

long <- read.table("/gpfs01/bhcbio/projects/UK_Biobank/20190102_UK_Biobank_Imaging/Results/GWAS/GWAS_diastolic_BOLT/Results/bolt_long_PDSR_full.bgen.stats", header=TRUE)
longTab <- long %>% dplyr::rename(P=P_BOLT_LMM, ESTIMATE=BETA)
longTab$MAF <- apply(cbind(longTab$A1FREQ, 1-longTab$A1FREQ), 1, min)
rm(long)

lav <- read.table("/gpfs01/bhcbio/projects/UK_Biobank/20190102_UK_Biobank_Imaging/Results/GWAS/GWAS_diastolic_BOLT/Results/bolt_LAV_full.bgen.stats", header=TRUE)
lavTab <- lav %>% dplyr::rename(P=P_BOLT_LMM, ESTIMATE=BETA)
lavTab$MAF <- apply(cbind(lavTab$A1FREQ, 1-lavTab$A1FREQ), 1, min)
rm(lav)

bgen_file_for_chrom <- list()
for (i in 1:22) {
  bgen_file_for_chrom[[i]] = file.path(config::get("imputation_dir"), paste0("_004_ukb_imp_chr", 
                                                                             i, "_v3.bgen"))
}
sample_file <- config::get("sample_autosome_file")
tmpdir <- "/gpfs01/bhcbio/projects/UK_Biobank/20190102_UK_Biobank_Imaging/test"
plink.tmpdir <- "/gpfs01/bhcbio/projects/UK_Biobank/20190102_UK_Biobank_Imaging/test"
withdrawn_samples <- "/gpfs01/bhcbio/projects/UK_Biobank/20190102_UK_Biobank_Imaging/Data/Original/WithdrawnSamples/w40616_20200820.csv"
tmpdir="/gpfs01/bhcbio/projects/UK_Biobank/20190102_UK_Biobank_Imaging/test"
plink.tmpdir="/gpfs01/bhcbio/projects/UK_Biobank/20190102_UK_Biobank_Imaging/test"

# get MRI subjects
mri_subjects <- data.table::fread("/gpfs01/bhcbio/projects/UK_Biobank/20190102_UK_Biobank_Imaging/Data/Original/Phenotypes/ImperialPhenos/clinical_measures_39k_collated_11052020.csv") %>% dplyr::pull(eid_40616)
mri_subjects_2 <- data.table::fread("/gpfs01/bhcbio/projects/UK_Biobank/20190102_UK_Biobank_Imaging/Data/Original/Phenotypes/ImperialPhenos/clinical_measures_27k_collated_qced_04032020.csv")%>% dplyr::pull(eid_40616)
phenotypes <- readSexAndAge() 
imaging <- mri_subjects
non_imaging <- phenotypes$SID[which(!(phenotypes$SID %in% mri_subjects))] 
## del end


###################### Step 1: perform clumping to identify independent signals #############################

#### radial PDSR

## obtain list of SNPs via clumping
radial <- gwas_clumping(score_input=radialTab, score_threshold=10^-6, clump_r2=0.2, tmpdir=tmpdir, 
                        plink.tmpdir=plink.tmpdir, snpset=NA,
                        minMaf=0.0001, bgen_file_for_chrom=bgen_file_for_chrom, sample_file=sample_file)
## remove list structure and exclude empty entries (are generated if no SNP is selected for a specific CHR)
radial <- unlist(radial)
radialSNPs <- radial[which(radial!="")]

#### long PDSR

## obtain list of SNPs via clumping
long <- gwas_clumping(score_input=longTab, score_threshold=10^-6, clump_r2=0.2, tmpdir=tmpdir, 
                      plink.tmpdir=plink.tmpdir, snpset=NA,
                      minMaf=0.0001, bgen_file_for_chrom=bgen_file_for_chrom, sample_file=sample_file)
## remove list structure and exclude empty entries (are generated if no SNP is selected for a specific CHR)
long <- unlist(long)
longSNPs <- long[which(long!="")]

#### LAVmax (indexed to BSA)

## obtain list of SNPs via clumping
lav <- gwas_clumping(score_input=lavTab, score_threshold=10^-6, clump_r2=0.2, tmpdir=tmpdir, 
                     plink.tmpdir=plink.tmpdir, snpset=NA,
                     minMaf=0.0001, bgen_file_for_chrom=bgen_file_for_chrom, sample_file=sample_file)
## remove list structure and exclude empty entries (are generated if no SNP is selected for a specific CHR)
lav <- unlist(lav)
lavSNPs <- lav[which(lav!="")]


###################### Step 2: compute and evaluate genetic score #############################################

## Assume that:
## - geno is a data frame containing the genetic data (columns SID plus additional columns for all genotypes, column
##   names of geno_data match entries in radialSNPs). Data is included for all QCed European subjects.
## - pheno is a data frame containing the measured values for radial PDSR, a column with the SIDs plus
##   all covariates necessary for the model
## - imaging is a vector containing SIDs of all QCed European samples who were part of the imaging study (releases 1-3)
## - non-imaging is a vector containing SIDs of all QCed European samples who were NOT part of the imaging study (release (1-3))
## - pca is a data frame containing SID plus the first ten principal components
## - phenotypes_bin is a data frame containing all binary traits which we want to consider for the PheWAS
## - phenotypes_quant is a data frame containing all quantitative traits which we want to consider for the PheWAS

##del start
geno <- getTibbleFromSnpMatrix(getImputedV3GenotypesForVariants(c(radialSNPs, longSNPs, lavSNPs)))
pca <- readPCA()
age <- readSexAndAge(2)
quant <- readProcessedQuantTraits(2) %>% dplyr::rename(lav=LAV.max_log)
quant2 <- readProcessedQuantTraits(0) 
covariates <- age

tmp <- dplyr::left_join(pca, age, by="SID")
pheno <- dplyr::left_join(tmp, quant, by="SID")
phenotypes_quant <- quant2
phenotypes_bin <- readProcessedBinaryTraits()[,1:3]

names(geno)[-1] <- paste0("x_", names(geno)[-1])
names(geno)[-1] <- gsub(":", "_", names(geno)[-1])
radialSNPs <-  paste0("x_", radialSNPs)
radialSNPs <- gsub(":", "_", radialSNPs)
longSNPs <-  paste0("x_", longSNPs)
longSNPs <- gsub(":", "_", longSNPs)
lavSNPs <-  paste0("x_", lavSNPs)
lavSNPs <- gsub(":", "_", lavSNPs)
## end del

#### radial PDSR

## compute genetic score
radialScore <- computeGeneticScoreForSnps(radialSNPs, trait="radial_PDSR", geno_data = geno, 
                                                  pheno_data=pheno)

## split into imaging/non-imaging subset
radialScoreMRI <- radialScore$Score %>% dplyr::filter(SID %in% imaging) %>% dplyr::select(SID, Score)
radialScoreNonMRI <- radialScore$Score %>% dplyr::filter((SID %in% non_imaging)) %>% dplyr::select(SID, Score)

## calculate explained variability with score
explained_radial <- dplyr::inner_join(radialScoreMRI, phenotypes, by="SID")
summary(lm(radial_PDSR~Score, data=explained_radial))
# explained variability = R^2 of this model

## run PheWAS
QuantRadial <- getQuantPheWASResults(radialScoreNonMRI, covariates, pca, phenotypes=phenotypes_quant, 
                                 scale_phenotypes = TRUE)
BinRadial <- getPheWASResults(radialScoreNonMRI, covariates, pca, phenotypes=phenotypes_bin, cores = 4)

## combine results and apply multiplicity correction
Radial <- rbind(QuantRadial %>% dplyr::select(phenotype, estimate, conf.low, conf.high, p.value, statistic), 
              BinRadial %>% dplyr::select(phenotype, estimate, conf.low, conf.high, p.value, statistic)) 
RadialMT <- Radial %>% dplyr::mutate(p.adjusted=p.adjust(p.value, method="fdr")) %>% dplyr::filter(p.adjusted < 0.05)

#### long PDSR

## compute genetic score
longScore <- computeGeneticScoreForSnpsModified(longSNPs, trait="long_PDSR", geno_data = geno, 
                                                  pheno_data=phenotypes)

## split into imaging/non-imaging subset
longScoreMRI <- longScore$Score %>% dplyr::filter(SID %in% imaging) %>% dplyr::select(SID, Score)
longScoreNonMRI <- longScore$Score %>% dplyr::filter((SID %in% non_imaging)) %>% dplyr::select(SID, Score)

## calculate explained variability with score
explained_long <- dplyr::inner_join(longScoreMRI, phenotypes, by="SID")
summary(lm(long_PDSR~Score, data=explained_long))
# explained variability = R^2 of this model

## run PheWAS
QuantLong <- getQuantPheWASResults(longScoreNonMRI, covariates, pca, phenotypes=quant, 
                                     scale_phenotypes = FALSE)
BinLong <- getPheWASResults(longScoreNonMRI, covariates, pca, phenotypes=bin, cores = 4)

## combine results and apply multiplicity correction
Long <- rbind(QuantLong %>% dplyr::select(phenotype, estimate, conf.low, conf.high, p.value, statistic), 
             BinLong %>% dplyr::select(phenotype, estimate, conf.low, conf.high, p.value, statistic)) 
LongMT <- Long %>% dplyr::mutate(p.adjusted=p.adjust(p.value, method="fdr")) %>% dplyr::filter(p.adjusted < 0.05)

#### LAVmax (indexed to BSA)

## compute genetic score
lavScore <- computeGeneticScoreForSnpsModified(lavSNPs, trait="lav", geno_data = geno, 
                                                pheno_data=phenotypes)

## split into imaging/non-imaging subset
lavScoreMRI <- lavScore$Score %>% dplyr::filter(SID %in% imaging) %>% dplyr::select(SID, Score)
lavScoreNonMRI <- lavScore$Score %>% dplyr::filter((SID %in% non_imaging)) %>% dplyr::select(SID, Score)

## calculate explained variability with score
explained_lav <- dplyr::inner_join(lavScoreMRI, phenotypes, by="SID")
summary(lm(lav_PDSR~Score, data=explained_lav))
# explained variability = R^2 of this model


## run PheWAS
QuantLav <- getQuantPheWASResults(lavScoreNonMRI, covariates, pca, phenotypes=quant, 
                                   scale_phenotypes = FALSE)
BinLav <- getPheWASResults(lavScoreNonMRI, covariates, pca, phenotypes=bin, cores = 4)

## combine results and apply multiplicity correction
Lav <- rbind(QuantLav %>% dplyr::select(phenotype, estimate, conf.low, conf.high, p.value, statistic), 
             BinLav %>% dplyr::select(phenotype, estimate, conf.low, conf.high, p.value, statistic)) 
LavMT <- Lav %>% dplyr::mutate(p.adjusted=p.adjust(p.value, method="fdr")) %>% dplyr::filter(p.adjusted < 0.05)

#### Selection of results for displaying 

## which traits are associated to all PDSR
intersect(intersect(LAVMT$phenotype, LongMT$phenotype), RadialMT$phenotype)

## which traits are associated to at least one trait?
in_any <- unique(union(union(LavMT$phenotype, LongMT$phenotype), RadialMT$phenotype))
Lav_inany <- Lav %>% dplyr::filter(phenotype %in% in_any)
Radial_inany <- Radial %>% dplyr::filter(phenotype %in% in_any)
Long_inany <- Long %>% dplyr::filter(phenotype %in% in_any)

## which binary traits are associated to at least one trait
Lav_bin <- Lav %>% dplyr::filter((phenotype %in% in_any) & (phenotype %in% BinLav$phenotype))
Radial_bin <- Radial %>% dplyr::filter((phenotype %in% in_any) & (phenotype %in% BinRadial$phenotype))
Long_bin <- Long %>% dplyr::filter((phenotype %in% in_any) & (phenotype %in% BinLong$phenotype))
top_bin <- Long_bin$phenotypes

## which quantitative traits are associated to at least one trait and which ones are most significant
LAV_quant <- LAV %>% dplyr::filter((phenotype %in% in_any) & (phenotype %in% QuantLav$phenotype))
Radial_quant <- Radial %>% dplyr::filter((phenotype %in% in_any) & (phenotype %in% QuantRadial$phenotype))
Long_quant <- Long %>% dplyr::filter((phenotype %in% in_any) & (phenotype %in% QuantLong$phenotype))
# identify top 10 associations across quantitative traits to be used for display in main paper
pe <- apply(cbind(Lav_quant$p.value, Radial_quant$p.value, Long_quant$p.value), 1, min)
top_traits <- Long_quant$phenotype[order(pe)[1:10]]

###################### Step 3: MR analysis ###################################################################

## del start
diastolic_gwas <- radialTab
nondiastolic_gwas <- radialTab
## del stop

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
mr_analysis(diastolic_gwas=diastolic_gwas, nondiastolic_gwas=nondiastolic_gwas, name_diastolic="PDSR_ll", name_nondiastolic="diastolic_blood_pressure", score_threshold=10^-6,
clump_kb=1000, clump_r2=0.2, tmpdir=tmpdir, bgen_file_for_chrom=bgen_file_for_chrom, sample_file=sample_file, minMaf=0)
