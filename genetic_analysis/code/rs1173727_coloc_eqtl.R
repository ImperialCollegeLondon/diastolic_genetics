# #### Script for NPR3 analyses ####
#==================================#

# generic libs
library(data.table)
library(dplyr)
library(tibble)
# for eqtls
require("org.Hs.eg.db")
require(ggplot2)
require(grid)
require(httr)
require(jsonlite)
require(xml2)

# load function set
source("./manuscript_code_for_git/coloc_eqtl_rs1173727_functions.R")

#### Coloc Analysis ####
gene_info <- data.table::data.table(external_name="NPR3", 
                                    start = 32689070,
                                    end = 32791724,
                                    chr = 5)
# fetch data to be analysed for colocalization
input_data <- coloc_get_my_data(gene_info=gene_info, flanking=200000, subset=1:5)

# coloc of all traits
coloc_out <- coloc_loop(input_data)



#### GWAS NPR3 SNP ####
snp <- "rs1173727"

#### examples for eqtl plots ####

# fetch eqtls
eqtls <- list()
for(i in 1:length(locus_genes)) {
  tf <- locus_genes[i]
  eqtls[[tf]] <- get_eqtl_ensembl(gene=tf)
}

# setup dataframe with gver 38 infor for npr3
gene_info_npr3 <- data.table::data.table(external_name="NPR3", start=32689070, end=32791724, seq_region_name=5)

flank <- 500000
for(i in 1:length(eqtls)) {
  eqt <- eqtls[[i]]
  mygene <- names(eqtls)[i]
  # adjust left and right borders as you wish
  flank_left <- flank_right <- flank
  plot_eqtl_wrap(mygene, tab=eqt, flank_left=flank_left, flank_right=flank_right, 
                 mark_region=NULL, #ifelse(mygene %in% c("NPR1", "NPR2"), NULL, gene_info_npr3),
                 mark_snps = locus_snps[1:3],
                 mark_snps_size = 3,
                 sig_snps_size = 2,
                 sig_snps_p_cut = 5)
}


# show NPR3 locus eQTL, marking LOC340113
tab <- eqtls[["NPR3"]]
locx_info <- data.table::data.table(start=32947549,
                                    end=32962573,
                                    external_name="LOC340113")
plot_eqtl_flanking_region(tab=tab, gene_info = gene_info_npr3, 
                          flank_left = 400000, flank_right=400000, mark_region = locx_info,
                          mark_snps = locus_snps,
                          mark_snps_size = 3,
                          sig_snps_size = 2,
                          sig_snps_p_cut = 5,
                          debug=F)







