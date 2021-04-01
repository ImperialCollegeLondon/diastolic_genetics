###################### Step 1: perform clumping to identify independent signals #############################

## Clumping is performed with PLINK

## However, clumping is performed with PLINK1.9 which can, on the other hand, not work with
## bgen files (UKBB data format for imputed data that we want to use).
## That is why we first transform the bgen files to bed files for the subset of all SNPs with P-values below the defined 
## threshold (filter_bgen_to_bed) with PLINK2.0 and afterwards perform clumping based on the bed files with PLINK1.9
## (plink_clumping). The function gwas_clumping is combining both steps and returning the SNP list.

## Function to filter bgen files to selected SNPs and export data as bed file
##
## Input:
## - rs_file: path to file were SNP RSIDs are stored (generated within gwas_clumping)
## - outdir: a directory to which the bed file is written to
## - bgen_file_for_chrom: the path to the input bgen-file (for one CHR as provided by UKBB)
## - sample_file: file path for sample file from UKBB (belonging to bgen-file)
## - withdrawn_samples: path to file with withdrawn samples SID (will be excluded prior to analysis)
##
## Output:
## - File path with bed file (string)

filter_bgen_to_bed <- function (rs_file, outdir, bgen_file_for_chrom, sample_file, withdrawn_samples) {
  # extract SNPs from SNP list (rs_file)
  filtered_bgen = tempfile(paste0("filtered_", stringr::str_replace(base::basename(rs_file), 
                                                                    "\\.", "_")), outdir, fileext = ".bgen")
  cmd = paste0("bgen/bgen-112018/bin/bgenix -g ", bgen_file_for_chrom, " -incl-rsids ", rs_file, 
               " > ", filtered_bgen)
  system(cmd)
  print(cmd)
  
  # export data as bed
  cmd = paste0("Plink_2.00_20190716/plink2 --out \"", plinkfile, "\" ")
  cmd = paste0(cmd, " --out \"", plinkfile, "\" ")
  cmd = paste0(cmd, "--bgen \"", filtered_bgen, "\" ")
  cmd = paste0(cmd, "--sample \"", sample_file, 
               "\" ")
  cmd = paste0(cmd, "--remove-fam \"", withdrawn_samples, 
               "\" ")
  cmd = paste0(cmd, "--real-ref-alleles ")
  cmd = paste0(cmd, "--rm-dup exclude-all ")
  cmd = paste0(cmd, "--make-bed ")
  system(cmd, ignore.stdout = TRUE)
  
  return(plinkfile)

  }

## Function to perform the actual clumping
##
## Input
## - plink_basename: path to bfiles
## - score_file: path to file where GWAS summary stats for selected SNPs are stored (see gwas_clumping)
## - outdir: path to a directory where results will be stored
## - clump_p1: value for p1 threshold (see PLINK documentation for details)
## - clump_kb: value for clumping window (see PLINK documentation for details)
## - clump_r2: value for R2 threshold (see PLINK documentation for details)
## - withdrawn_samples: path to file with withdrawn samples SID (will be excluded prior to analysis)

plink_clumping <- function (plink_basename, score_file, outdir, clump_p1, clump_kb, 
          clump_r2, withdrawn_samples)  {
  plinkfile = tempfile("plink", outdir)
  
  ## run clumping with PLINK
  cmd <- "plink-1.90/plink"
  cmd = paste0(cmd, " --out \"", plinkfile, "\" ")
  cmd = paste0(cmd, " --bfile \"", plink_basename, "\" ")
  cmd = paste0(cmd, "--remove-fam \"", withdrawn_samples, 
               "\" ")
  cmd = paste0(cmd, "--clump \"", score_file, "\" ")
  cmd = paste0(cmd, "--clump-snp-field ID ")
  cmd = paste0(cmd, "--clump-p1 ", clump_p1, " ")
  cmd = paste0(cmd, "--clump-p2 0.001 ")
  cmd = paste0(cmd, "--clump-kb ", clump_kb, " ")
  cmd = paste0(cmd, "--clump-r2 ", clump_r2, " ")
  system(cmd, ignore.stdout = TRUE)
  
  ## import clumping files and extract relevant information (SNP, P-value)
  ret = data.table::fread((list.files(path = outdir, pattern = paste0(basename(plinkfile), 
                                                                      ".clumped"), full.names = T))[1]) %>% 
                                                                      dplyr::select(SNP, P)
  ## delete unnecessary files
  purrr::map(list.files(path = dirname(plink_basename), pattern = paste0(basename(plink_basename), 
                                                                         ".*"), full.names = T), file.remove)
  return(ret)
}

## Function to select SNPs with clumping based on GWAS summary statistics
##
## Input:
## - score_input: data frame with GWAS result summary statistics 
## - pval: name of column with the P-values
## - rsid: name of column with RSIDs
## - chr: name of column with information on chromosome
## - score_threshold: cut-off for P-value - only SNPs with smaller P-values are considered for clumping
## - clump_kb: windows for clumping (see  PLINK documentation for details)
## - clump_r2: R2 threshold (see PLINK documention for details)
## - snpset: vector containing RSIDs that clumping should be restricted to
## - tmpdir: a directory to which the filtered temporrary GWAS results are written to
## - plink.tmpdir: a directory to which temporary PLINK files are written to
## - minMaf: minimal MAF so that we consider the SNP
## - bgen_file_for_chrom: the path to the input bgen-file (for one CHR as provided by UKBB)
## - sample_file: location of sample file corresponding to bgen genotypes
##
## Output: a list (one entry per chromsome) containing data frames with selected SNPs and P-values

gwas_clumping <- function(score_input, pval="P", rsid="SNP", chr="CHR", MAFfreq="MAF", score_threshold=5*10^-8, 
                          clump_kb=1000, clump_r2=0.1, 
                          snpset=NA, tmpdir, plink.tmpdir, minMaf=0, bgen_file_for_chrom, sample_file) {
  
  ## reformat GWAS to "standard" column names and export
  names(score_input)[which(names(score_input)==pval)] <- "P"
  names(score_input)[which(names(score_input)==rsid)] <- "ID"
  names(score_input)[which(names(score_input)==chr)] <- "Chr"
  names(score_input)[which(names(score_input)==MAFfreq)] <- "MAF"
  score_file <- tempfile(tmpdir = tmpdir)
  
  if (any(!is.na(snpset))) {
    score_input <- score_input %>% dplyr::filter(ID %in% snpset)
  }
  
  ## select snps below threshold and with reasonable MAF
  scores_selected <- score_input %>% dplyr::mutate(P = as.numeric(P)) %>% dplyr::filter(P < score_threshold, MAF > minMaf)
  
  ## write out data
  data.table::fwrite(scores_selected, file = score_file, sep = "\t")
  filtered_rs_file = tempfile(tmpdir = plink.tmpdir)
  scores_selected %>% dplyr::select(ID) %>% readr::write_delim(path = filtered_rs_file)
  
  ## count number of snps per CHR below threshold (only do following calculations if there is at least a snp on the CHR)
  freq.file <- scores_selected
  freqs <- numeric(22)
  for (i in 1:22) {
    freqs[i] <- sum(freq.file$Chr==i)
  }
  
  ## do clumping
  ret <- list()
  for (i in 1:22) {
    if (freqs[i]>0) {
      plink_file = filter_bgen_to_bed(filtered_rs_file, outdir=tmpdir, bgen_file_for_chrom[[i]], sample_file, withdrawn_samples)
      ret[[i]] = unlist(plink_clumping(plink_file, score_file=score_file, outdir=tmpdir, clump_p1=score_threshold, 
                                       clump_kb=clump_kb, clump_r2=clump_r2, withdrawn_samples)[,1])
      purrr::map(list.files(path = dirname(plink_file), pattern = paste0(basename(plink_file),   ".*"),
                            full.names = T),   file.remove)
    }
    if (freqs[i]==0) {
      ret[[i]] <- c("")
    }
    names(ret)[i] <- paste0("Chr", i)
  }
  file.remove(filtered_rs_file)
  return(ret)
}


###################### Step 2: compute and evaluate genetic score #############################################


## Function to compute a PRS based on a list of SNPs. Model is estimated using UK Biobank data using
## age (at imaging visit), sex, array, PC1-PC10 and the imaging center as covariates
## 
## Input:
## - snps: The vector of SNPs to be included in the score.
## - trait 	The trait for which the score is computed. Must be found in pheno_data.
## - pheno_data: A data frame with all phenotype data (SID column, plus trait column, also
##   including covariates like pca, sex, age, array)
## - geno_data: A data frame with the genotypes, one column per SNP and an SID column. Genotypes are coded 0/1/2.
##   Column names for genotypes need to match SNPs in snps input vector.
## - mean_impute_genotypes: Should missing genotype values be imputed (with the mean value for that genotype)?
## 
## Output: A list containing the geno_data with an additional Score column, the score weights, and p-values 

computeGeneticScoreForSnps <- function(snps, trait, pheno_data, geno_data, 
                                               mean_impute_genotypes = T){
  
  score_formula = formula(paste0(trait, " ~ ", paste(snps, 
                                                     collapse = " + "), " + Age_Imaging + Sex + Array + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + 
                                                     PC7 + PC8 + PC9 + PC10+Center1_Imaging+Center2_Imaging"))
  if (mean_impute_genotypes) {
    geno_data = tidyimpute::impute(geno_data, na.tools::na.mean)
  }
  fit =  lm(score_formula, data = pheno_data %>% dplyr::inner_join(geno_data))
            
  list(Score = geno_data %>% dplyr::mutate(Score = as.matrix(geno_data[,snps]) %*% fit$coefficients[snps]), 
       Weights = fit$coefficients)
}


## Run PheWAS for association of PRS with quantitative traits
## 
## Input:
## - score: A data frame containing two columns, the first called "SID" containing the subject IDs, and 
##  the second containing the (numeric) score (e.g. genotypes of a SNP), for which associations are calculated.
## - covariates: A data frame containing all covariates that are used in the linear regression (e.g. Sex, Age), and a 
##   column named SID specifying the sample ids. Default is to read all pre-processed phenotypes using readSexAndAge()
## - pca: A data frame containing the principal components of the genotype matrix (used as covariates in the linear regression, analogous to the covariates table) Default is to read all pre-processed phenotypes using readPca()
## - phenotypes: A data frame with all quantitative traits, for which associations should be computed.
## - scale_phenotypes: If TRUE, the phenotypes are scaled to mean 0 and standard deviation 1.
##
## Output: A data.frame with columns term, estimate, std.error, statistic, p.value, phenotype, conf.low, conf.high, 
##         n specifying the results of the linear regressions. term is the score name (i.e. the name of the second column of 
##         the score argument), phenotype is the name of the trait of the association. The other columns contain the results 
##         of the linear regression, i.e. the estimate for the score with p-value, confidence intervals etc. 
getQuantPheWASResults <- function (score, covariates, pca, phenotypes, 
          scale_phenotypes = FALSE) {
  if (scale_phenotypes) {
    phenotypes <- phenotypes %>% dplyr::mutate_at(dplyr::vars(-SID), ~as.numeric(scale(.)))
  }
  score_name = names(score)[2]
  phenocov <- dplyr::left_join(phenotypes, covariates, by = c("SID"))
  data = dplyr::left_join(score, phenocov, by = "SID") %>%  dplyr::left_join(pca, by = "SID")
  result <- data.frame()
  formula_covariates = "Age + Sex + Array + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10"
  for (var_name in names(phenotypes %>% dplyr::select(-SID))) {
    my_formula <- paste0("`", var_name, "` ~ ", score_name, 
                         " + ", formula_covariates)
    fit <- tryCatch(lm(as.formula(my_formula), data = data), 
                    error = function(e) NULL)
    if (length(fit) == 0) {
      warning("Error in model fitting, please check input.")
      next
    }
    tmp <- broom::tidy(fit) %>% dplyr::filter(term == names(score)[2])
    if (nrow(tmp) == 0) {
      tmp <- data.frame(term = NA, estimate = NA, std.error = NA, 
                        statistic = NA, p.value = NA, phenotype = var_name, 
                        conf.low = NA, conf.high = NA, n = NA)
    }
    else {
      tmp <- cbind((broom::tidy(fit) %>% dplyr::mutate(phenotype = var_name)), 
                   confint.default(fit) %>% as.data.frame() %>% 
                     dplyr::rename(conf.low = 1, conf.high = 2))
      tmp <- tmp %>% dplyr::filter(term == names(score)[2])
      counts <- data.frame(n = length(fit$residuals))
      tmp <- cbind(tmp, counts)
    }
    result <- rbind(result, tmp)
  }
  result
}

## Compute PheWAS for association of PRS with binary phenotypes using logistic regression
##
## Input:
## - score: A data frame containing two columns, the first called "SID" containing the subject IDs, and the second 
##   containing the (numeric) score (e.g. genotypes of a SNP), for which associations are calculated
## - covariates: A data frame containing all covariates that are used in the linear regression (e.g. Sex, Age). 
## - pca: A data frame containing the principal components of the genotype matrix (used as covariates in the  
##   regression, analogous to the covariates table) 
## - phenotypes: A data frame with all binary traits (coded as 0,1), for which associations should be computed. 
## - cores: Number of cores used for parallelization, default is 4.
##
## Output: A data.frame with columns term,estimate,std.error,statistic,p.value,phenotype,conf.low,conf.high,cases, 
##         controls specifying the results of the logistic regressions. term is the score name (i.e. the name of the second 
##         column of the score argument), phenotype is the name of the trait of the association. The other columns contain 
##         the results of the linear regression, i.e. the estimate for the score with p-value, confidence intervals etc. 
getPheWASResults <- function (score, covariates, pca, phenotypes, cores = 4) {
  score_name = names(score)[2]
  data = dplyr::inner_join(score, covariates, by = "SID") %>% 
    dplyr::inner_join(pca, by = "SID")

  data <- data %>% dplyr::inner_join(phenotypes, by = "SID")
  formula_covariates = "Age + Sex + Array + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10"
  inner_phewas <- function(var_name, data, score_name,
                           formula_covariates = "Age + Sex + Array + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10") {
    
    my_formula <- paste0("`", var_name, "` ~ ", score_name, " + ", formula_covariates)
    
    reg_data <- data %>% filter(!is.na(!!as.name(score_name))) %>% filter(!is.na(!!as.name(var_name)))
    
    safe_glm <- purrr::safely(glm)
    fit <- safe_glm(as.formula(my_formula), family='binomial', data = reg_data)
    counts <- data.frame(cases=sum(fit$result$y==1), controls=sum(fit$result$y==0))
    
    if (is.null(fit$result)) {
      return (data.frame(term=NA,estimate=NA,std.error=NA,statistic=NA,p.value=NA,phenotype=var_name,conf.low=NA,conf.high=NA, cases=NA, controls=NA))
    }
    tmp <- broom::tidy(fit$result) %>% filter(term==score_name)
    
    if (nrow(tmp)==0) {
      return (data.frame(term=NA,estimate=NA,std.error=NA,statistic=NA,p.value=NA,phenotype=var_name,conf.low=NA,conf.high=NA, cases=NA, controls=NA))
    } else {
      tmp <- cbind((broom::tidy(fit$result) %>% mutate(phenotype = var_name)), broom::confint_tidy(fit$result, func = stats::confint.default))
      tmp <- tmp %>% filter(term==score_name)
      tmp <- cbind(tmp, counts)
    }
    tmp
  }
  

  ret <- bind_rows(parallel::mclapply(as.list(names(phenotypes[2:ncol(phenotypes)])), 
                                        FUN = inner_phewas, data, score_name, formula_covariates, 
                                        mc.cores = cores))
  return(ret)
}

## Compute PheWAS for association of PRS with binary phenotypes using logistic regression
##
## Input:
## - snps_included: RSIDs (or similar identifier) - should correspond to genetic data stored in 'geno'
## - trait: name of trait of interest (association to PRS as called in 'phenotypes')
## - imaging_trait: name of imaging trait (trait that PRS will be constructed for, need to be named in 'phenotypes_score')
## - geno: tibble with genetic data (first column SID, afterwards SNP matrix)
## - restrict_bin: list of binary traits (trait needs to be within restrict_quant or restrict_bin)
## - restrict_quant: list of quantitative traits (trait needs to be within restrict_quant or restrict_bin)
## - covariates: data set which contains covariates (e.g., age, sex) for PRS construction (first column SID, afterwards phenotypes)
## - pca: dataset containing the PCA as used for the PheWAS (first column SID, afterwards phenotypes)
## - phenotypes_quant: dataset containing the phenotype (named trait) (first column SID, afterwards phenotype - if quantitative)
## - phenotypes_bin: dataset containing the phenotype (named trait) (first column SID, afterwards phenotype - if binary)
## - phenotypes_score: dataset containing the imaging trait data (first column SID, afterwards phenotype)
## - imaging: vector with subject IDs belonging to the imaging set
##
## Output: A vector with columns  phenotype, estimate, lower and upper CI, p-value and the statistics for the specific PRS and the trait
leave_one_out_cv <- function(snps_include, trait, imaging_trait, geno, covariates, pca, phenotypes_quant, phenotypes_bin, phenotypes_score) {
  score <- computeGeneticScoreForSnps(snps_include, trait=imaging_trait, geno_data = geno, pheno_data=phenotypes_score, imaging=imaging)
  score_nonmri <- score$Score %>% dplyr::filter(!(SID %in% imaging)) %>% dplyr::select(SID, Score)
  if (trait %in% names(phenotypes_bin)) {
    res <- getPheWASResults(score=score_nonmri, covariates=covariates, pca=pca, phenotypes=phenotypes_bin %>% dplyr::select(SID, trait))
  }
  if (trait %in% names(phenotypes_quant)) {
    res <- getQuantPheWASResults(score=score_nonmri, phenotypes = phenotypes_quant  %>% dplyr::select(SID, trait), 
                                 covariates=covariates, pca=pca, scale_phenotypes = TRUE)
  }
  return(res %>% dplyr::select(phenotype, estimate, conf.low, conf.high, p.value, statistic))
}

## Code to generate PRS forestplot as displayed in the paper
##
## Input: 
## - figure_data: a data frame with the following entries:
## * endpoint: names of phenotypes (only used for internal purposes, not displayed)
## * estimate: point estimate of associations phenotype vs. PRS 
## * conf.low: lower CI limit of associations phenotype vs. PRS 
## * conf.high: upper CI limit of associations phenotype vs. PRS 
## * prs_trait: name of diastolic function parameter:   "PDSR_ll", "PDSR_rr" or "LAVmax_i"
## * description: read name for trait (will be displayed)
## * p.adjusted: p-value- if >0.05: line will be displayed in grey
## - title: title to display in figure
## - xlim: limits on x-axis
## - xlab: label on x-axis
## - zero: position of vertical line
## - vjust, ylim, topmargin: features for fine-tuning display of results; please see in the code
plot_PRS <- function(figure_data, title, xlim, xlab, zero=1, vjust=0.5, ylim=0, topmargin=0.1) {
  quantiles <- figure_data %>% dplyr::arrange(endpoint)  %>% rowid_to_column() %>% 
    dplyr::mutate(pos = if_else((rowid %% 3) == 1, 1.5, 1)) %>% 
    mutate(ci = sprintf("%.2f [%.2f,%.2f]",estimate,conf.low, conf.high)) %>% 
    mutate(description_str = if_else(prs_trait!="PDSR_ll", "", as.character(description))) %>% 
    mutate(ff=1)
  quantiles <- quantiles %>% add_row(ci="Point estimate [95% CI]", pos=max(quantiles$pos+1), 
                                     description_str="", prs_trait="PRS", conf.low=-NA, ff=1) 
  quantiles$pos <- cumsum(quantiles$pos)
  quantiles$prs_trait[which(quantiles$prs_trait=="PDSR_ll")] <- "PDSR[ll]"
  quantiles$prs_trait[which(quantiles$prs_trait=="PDSR_rr")] <- "PDSR[rr]"
  quantiles$prs_trait[which(quantiles$prs_trait=="LAVmax_i")] <- "LAVmax[i]"
  
  # for lines in Plot
  seqline <- quantiles$pos[seq(4, length(quantiles$pos), by=3)]-1
  seqline <- seqline[-length(seqline)]
  colors <- c(rep(1:3, (length(unique(quantiles$description))-1)))
  for (i in 1:(nrow(quantiles)-1)) {
    if (quantiles$p.adjusted[i] >0.05) {
      colors[i]  <- "nonsig"
    }
  }
  man_col <- c("blueviolet", "deeppink4", "chartreuse4", "grey70")
  if (!(1 %in% colors)) {
    man_col <- c("deeppink4", "chartreuse4", "grey70")
  }
  if (!(2 %in% colors)) {
    man_col <- c("blueviolet",  "chartreuse4", "grey70")
  }
  if (!(3 %in% colors)) {
    man_col <- c("blueviolet", "deeppink4", "grey70")
  }
  
  
  lines <- ggplot(quantiles %>% dplyr::filter(!is.na(estimate)),
                  aes(x=estimate, y=pos - 0.5, color=colors, xmin=conf.low, xmax=conf.high)) +
    geom_pointrange(aes(col=colors))+
    geom_vline(aes(fill=colors),xintercept =zero, linetype=2) +
    geom_hline(aes(fill=colors),yintercept =seqline, linetype=3) +
    geom_errorbar(aes(xmin=conf.low, xmax=conf.high,col=colors),width=0.5,cex=1)+
    scale_y_continuous(limits=c(ylim-2,max(quantiles$pos)+1-ylim), expand=c(0.01,0.1))+
    coord_cartesian(xlim = c(min(quantiles$conf.low), max(quantiles$conf.high)))+    cowplot::theme_cowplot()  +
    theme(legend.position = "none",
          plot.margin = unit(c(10,10,10,0), "pt"), 
          axis.ticks.y = element_blank(),
          axis.line.y = element_blank(),
          axis.text.y = element_blank(),) +
    xlab(xlab) + ylab("") +
    scale_color_manual(values=man_col)
  
  table <- ggplot(quantiles  %>% dplyr::filter(!is.na(estimate)), aes(y=pos - 0.5, fontface=ff)) + 
    geom_text(hjust=0, aes(x=0, label=description_str)) +
    scale_fill_brewer(type = "seq", palette = "Greys") +   theme_void()+
    theme(
      panel.border = element_blank(),
      axis.ticks = element_blank(),
      axis.ticks.x = element_line(colour="white"),
      axis.line = element_blank(),
      axis.line.y = element_line(colour="white"),
      axis.text.y = element_blank(),
      axis.text.x = element_text(colour="white"),
      legend.position = "none",
      plot.margin = unit(c(10,0,10,0), "pt")) + 
    xlab(" ") + ylab("") +
    scale_y_continuous(limits=c(ylim-2,max(quantiles$pos)+1-ylim), expand=c(0.01,0.1))+coord_cartesian(xlim = c(0,4))
  
  # now add the title
  if (!is.null(title)) {
    plot_title <- cowplot::ggdraw() + 
      cowplot::draw_label(title,
                 fontface = 'bold',
                 size = 12,
                 x = 0,
                 hjust = 0, vjust=vjust
      ) +
      theme(
        # add margin on the left of the drawing canvas,
        # so title is aligned with left edge of first plot
        plot.margin = margin(0, 0, 0, 20)
      )
    
    plot_row <- cowplot::plot_grid(table, lines, ncol=2, rel_widths = c(0.25,1), label_size = 0, align = "h", axis = "b")
    cowplot::plot_grid(plot_title, plot_row, ncol = 1,
              # rel_heights values control vertical title margins
              rel_heights = c(topmargin, 1)
    )
  } else {
    cowplot::plot_grid(table, lines, ncol=2, rel_widths = c(1.5,1), label_size = 0, align = "h", axis = "b")
  }
}



###################### Step 3: MR analysis ###################################################################


## Calculate MR estimates for bi-directional MR analysis
##
## Input:
## - diastolic_gwas/nondiastolic_gwas: GWAS summary stats for the diastolic and nondiastolic trait with the following columns:
## * P: containing the P-values
## * SNP: containing the RSIDs
## * CHR: chromosome
## * MAF: minor allele frequency
## * ESTIMATE: point estimates
## * SE: standard error
## - name_diastolic/name_nondiastolic: name of diastolic/nondiastolic trait that will be used in result table
## - score_threshold: cut-off for P-value - only SNPs with smaller P-values are considered for clumping
## - clump_kb: windows for clumping (see  PLINK documentation for details)
## - clump_r2: R2 threshold (see PLINK documention for details)
## - tmpdir: a directory to which the filtered temporrary GWAS results are written to
## - bgen_file_for_chrom: the path to the input bgen-file (for one CHR as provided by UKBB)
## - sample_file: location of sample file corresponding to bgen genotypes
## - minMaf: minimal MAF so that we consider the SNP for clumping
##
## Result: a data frame with the results from both directionalities with the columns
##         Method (which MR-method), Estimate (MR estimate), Sd (standard deviation), P-value, Direction (diastolic parameter as cause or
##         consequence),  DirectionN (numeric indicator for direction), N_snps (Number of SNPs selected in clumping step),
##         ND (name non-diastolic trait), D (name diastolic trait)

mr_analysis <- function(diastolic_gwas, nondiastolic_gwas, name_diastolic, name_nondiastolic, score_threshold,
                        clump_kb, clump_r2, tmpdir, bgen_file_for_chrom, sample_file, minMaf=0) {
  
  ## find overlap between snpsets
  snpset <- intersect(diastolic_gwas$SNP, nondiastolic_gwas$SNP)
  
  ## select snps via clumping
  snps_diastolic <-  gwas_clumping(score_input=diastolic_gwas, score_threshold=score_threshold, clump_r2=clump_r2,clump_kb=clump_kb, tmpdir=tmpdir, 
                                   plink.tmpdir=tmpdir, snpset=snpset,
                                   minMaf=minMaf, bgen_file_for_chrom=bgen_file_for_chrom, sample_file=sample_file)
  snps_nondiastolic <-  gwas_clumping(score_input=nondiastolic_gwas, score_threshold=score_threshold, clump_r2=clump_r2,clump_kb=clump_kb, tmpdir=tmpdir, 
                                      plink.tmpdir=tmpdir, snpset=snpset,
                                      minMaf=minMaf, bgen_file_for_chrom=bgen_file_for_chrom, sample_file=sample_file)
  snps_nondiastolic <- unlist(snps_nondiastolic)
  if (any(snps_nondiastolic=="")) {
    snps_nondiastolic <- snps_nondiastolic[which(snps_nondiastolic!="")]
  }
  snps_diastolic <- unlist(snps_diastolic)
  if (any(snps_diastolic=="")) {
    snps_diastolic <- snps_diastolic[which(snps_diastolic!="")]
  }
  
  ## extract association results from GWAS
  res <- list()
  ## 1) diastolic => non diastolic
  mr_1a <- diastolic_gwas %>% dplyr::filter(SNP %in% snps_diastolic) %>% dplyr::select(SNP, ESTIMATE, SE)
  mr_1b <- nondiastolic_gwas %>% dplyr::filter(SNP %in% snps_diastolic) %>% dplyr::select(SNP, ESTIMATE, SE) 
  
  mr_1 <- dplyr::left_join(mr_1a %>% dplyr::rename(bx=ESTIMATE, bxse=SE), mr_1b  %>% dplyr::rename(by=ESTIMATE, byse=SE), by="SNP") 
  
  # calculate MR estimates
  input_1 <- MendelianRandomization::mr_input(bx=mr_1$bx, bxse=mr_1$bxse, 
                                              by=mr_1$by, byse=mr_1$byse, snps=mr_1$SNP)
  res[[1]] <- MendelianRandomization::mr_allmethods(input_1, method=c("main"))@Values %>% dplyr::filter(Method!="Simple median") %>%
    dplyr::rename(Sd=`Std Error`) %>% dplyr::select(Method, Estimate, Sd, `P-value`)
  
  presso1 <- data.frame(bx=mr_1$bx, bxse=mr_1$bxse, 
                        by=mr_1$by, byse=mr_1$byse, snps=mr_1$SNP)
  presso_tmp <-  try(MRPRESSO::mr_presso(BetaOutcome = "by", BetaExposure = "bx", SdOutcome = "byse", SdExposure = "bxse", 
                                         OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = presso1, NbDistribution = 1000,  
                                         SignifThreshold = 0.05))
  if (class(presso_tmp)=="try-error") {
    presso_tmp <- data.frame(Method="Outlier-corrected", Estimate=NA, Sd=NA, "P-value"=-1)
    names(presso_tmp)[4] <- "P-value"
  } else {
    presso_tmp <- presso_tmp$`Main MR results` %>% dplyr::rename(Method=`MR Analysis`, Estimate=`Causal Estimate`) %>% 
      dplyr::select(Method, Estimate, Sd, `P-value`) %>%  dplyr::filter(Method=="Outlier-corrected")
  }
  
  res[[1]] <- rbind(res[[1]], presso_tmp) %>%    dplyr::mutate(Direction = "Diastolic => Non-diastolic", 
                                                               DirectionN=1, N_snps=length(snps_diastolic)) 
  
  # 2) non diastolic =>  diastolic
  mr_2a <- diastolic_gwas %>% dplyr::filter(SNP %in% snps_nondiastolic) %>% dplyr::select(SNP, ESTIMATE, SE)
  mr_2b <- nondiastolic_gwas %>% dplyr::filter(SNP %in% snps_nondiastolic) %>% dplyr::select(SNP, ESTIMATE, SE) 
  
  mr_2 <- dplyr::left_join(mr_2b %>% dplyr::rename(bx=ESTIMATE, bxse=SE), mr_2a  %>% dplyr::rename(by=ESTIMATE, byse=SE), by="SNP") 
  
  # calculate MR estimates
  input_2 <- MendelianRandomization::mr_input(bx=mr_2$bx, bxse=mr_2$bxse, 
                                              by=mr_2$by, byse=mr_2$byse, snps=mr_2$SNP)
  res[[2]] <- MendelianRandomization::mr_allmethods(input_2, method=c("main"))@Values %>% dplyr::filter(Method!="Simple median") %>%
    dplyr::rename(Sd=`Std Error`) %>% dplyr::select(Method, Estimate, Sd, `P-value`)
  
  presso2 <- data.frame(bx=mr_2$bx, bxse=mr_2$bxse, 
                        by=mr_2$by, byse=mr_2$byse, snps=mr_2$SNP)
  presso_tmp <-  try(MRPRESSO::mr_presso(BetaOutcome = "by", BetaExposure = "bx", SdOutcome = "byse", SdExposure = "bxse", 
                                         OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = presso2, NbDistribution = 1000,  SignifThreshold = 0.05))
  if (class(presso_tmp)=="try-error") {
    presso_tmp <- data.frame(Method="Outlier-corrected", Estimate=NA, Sd=NA, "P-value"=-1)
    names(presso_tmp)[4] <- "P-value"
  } else {
    presso_tmp <- presso_tmp$`Main MR results` %>% dplyr::rename(Method=`MR Analysis`, Estimate=`Causal Estimate`) %>% 
      dplyr::select(Method, Estimate, Sd, `P-value`) %>%  dplyr::filter(Method=="Outlier-corrected")
  }
  res[[2]] <- rbind(res[[2]], presso_tmp) %>%    dplyr::mutate(Direction = "Non-diastolic => Diastolic", DirectionN=2, 
                                                               N_snps=length(snps_nondiastolic))
  
  ## combine both directions into one result data frame
  res <- dplyr::bind_rows(res[[1]], res[[2]]) %>% dplyr::mutate(ND=name_nondiastolic, D=name_diastolic)
  return(res)
}
