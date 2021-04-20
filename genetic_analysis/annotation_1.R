
require(tidyverse)
require(epiChoose)
require(ggrepel)
require(VennDiagram)
require(devtools)
require(ghql)
require(graphql)
require(rtracklayer)
require(httr)
require(org.Hs.eg.db)
require(readxl)


# MANHATTAN PLOTS ANNOTATED -----------------------------------------------

# manhattan plots

# radial
# stats = summary stats as output by bolt-lmm

# apply maf > 0.005 and hwe filter
stats = stats %>% filter(A1FREQ <= 0.995)
stats = stats %>% filter(SNP %in% snps_in)
stats = stats %>% dplyr::select(SNP, CHR, BP, P_BOLT_LMM_INF)
stats = dplyr::rename(stats, p_value=P_BOLT_LMM_INF, chr=CHR, pos=BP)
stats$phenotype = NA
stats$gene = NA
stats$phenotype[match(c("rs2234962","rs2644262","rs11535974"), stats$SNP)] = ""
stats$gene[match(c("rs2234962","rs2644262","rs11535974"), stats$SNP)] = c("BAG3","FHOD3","AC023158.1")
png(filename="data/manhattan_disc_radial.png", width=1200, height=350)
manh_plot(stats, pcut_label=1e-5)
dev.off()

stats = read_tsv(paste0("/gpfs01/bhcbio/projects/UK_Biobank/20190102_UK_Biobank_Imaging/Results/GWAS/GWAS_diastolic_BOLT/Results/bolt_", "radial_PDSR", "_full.bgen.stats"))
# apply maf > 0.005 and hwe filter
stats = stats %>% filter(A1FREQ <= 0.995)
stats = stats %>% filter(SNP %in% snps_in)
stats = stats %>% dplyr::select(SNP, CHR, BP, P_BOLT_LMM_INF)
stats = dplyr::rename(stats, p_value=P_BOLT_LMM_INF, chr=CHR, pos=BP)
stats$phenotype = NA
stats$gene = NA
stats$phenotype[match(c("rs2234962","rs2644262","rs9388001","rs11535974","rs11170519"), stats$SNP)] = ""
stats$gene[match(c("rs2234962","rs2644262","rs9388001","rs11535974","rs11170519"), stats$SNP)] = c("BAG3","FHOD3","GJA1","AC023158.1","SP1")
png(filename="data/manhattan_full_radial.png", width=1200, height=350)
manh_plot(stats, pcut_label=1e-5)
dev.off()

stats = read_tsv(paste0("/gpfs01/bhcbio/projects/UK_Biobank/20190102_UK_Biobank_Imaging/Results/GWAS/GWAS_diastolic_BOLT/Results/bolt_", "radial_PDSR", "_full_apr_covar.bgen.stats"))
# apply maf > 0.005 and hwe filter
stats = stats %>% filter(A1FREQ <= 0.995)
stats = stats %>% filter(SNP %in% snps_in)
stats = stats %>% dplyr::select(SNP, CHR, BP, P_BOLT_LMM_INF)
stats = dplyr::rename(stats, p_value=P_BOLT_LMM_INF, chr=CHR, pos=BP)
stats$phenotype = NA
stats$gene = NA
stats$phenotype[match(c("rs2234962","rs369533272"), stats$SNP)] = "Radial PDSR"
stats$gene[match(c("rs2234962","rs369533272"), stats$SNP)] = c("BAG3","FHOD3")
png(filename="data/manhattan_full_radial_apr_covar.png", width=1200, height=350)
manh_plot(stats, pcut_label=1e-5)
dev.off()

# lav

stats = read_tsv(paste0("/gpfs01/bhcbio/projects/UK_Biobank/20190102_UK_Biobank_Imaging/Results/GWAS/GWAS_diastolic_BOLT/Results/bolt_", "LAV", "_disc.bgen.stats"))
# apply maf > 0.005 and hwe filter
stats = stats %>% filter(A1FREQ <= 0.995)
stats = stats %>% filter(SNP %in% snps_in)
stats = stats %>% dplyr::select(SNP, CHR, BP, P_BOLT_LMM_INF)
stats = dplyr::rename(stats, p_value=P_BOLT_LMM_INF, chr=CHR, pos=BP)
stats$phenotype = NA
stats$gene = NA
stats$phenotype[match(c("rs59985551"), stats$SNP)] = ""
stats$gene[match(c("rs59985551"), stats$SNP)] = c("EFEMP1")
png(filename="data/manhattan_disc_lav.png", width=1200, height=350)
manh_plot(stats, pcut_label=1e-5)
dev.off()

stats = read_tsv(paste0("/gpfs01/bhcbio/projects/UK_Biobank/20190102_UK_Biobank_Imaging/Results/GWAS/GWAS_diastolic_BOLT/Results/bolt_", "LAV", "_full.bgen.stats"))
# apply maf > 0.005 and hwe filter
stats = stats %>% filter(A1FREQ <= 0.995)
stats = stats %>% filter(SNP %in% snps_in)
stats = stats %>% dplyr::select(SNP, CHR, BP, P_BOLT_LMM_INF)
stats = dplyr::rename(stats, p_value=P_BOLT_LMM_INF, chr=CHR, pos=BP)
stats$phenotype = NA
stats$gene = NA
stats$phenotype[match(c("rs1173727","rs59985551"), stats$SNP)] = ""
stats$gene[match(c("rs1173727","rs59985551"), stats$SNP)] = c("NPR3","EFEMP1")
png(filename="data/manhattan_full_lav.png", width=1200, height=350)
manh_plot(stats, pcut_label=1e-5)
dev.off()

stats = read_tsv(paste0("/gpfs01/bhcbio/projects/UK_Biobank/20190102_UK_Biobank_Imaging/Results/GWAS/GWAS_diastolic_BOLT/Results/bolt_", "LAV", "_full_apr_covar.bgen.stats"))
stats = stats %>% dplyr::select(SNP, CHR, BP, P_BOLT_LMM_INF)
stats = dplyr::rename(stats, p_value=P_BOLT_LMM_INF, chr=CHR, pos=BP)
stats$phenotype = NA
stats$gene = NA
stats$phenotype[match(c("rs1173727","rs59985551","rs35489511"), stats$SNP)] = "LAV"
stats$gene[match(c("rs1173727","rs59985551","rs35489511"), stats$SNP)] = c("NPR3","EFEMP1","CDK6")
png(filename="data/manhattan_full_lav_apr_covar.png", width=1200, height=350)
manh_plot(stats, pcut_label=1e-5)
dev.off()

# long

stats = read_tsv(paste0("/gpfs01/bhcbio/projects/UK_Biobank/20190102_UK_Biobank_Imaging/Results/GWAS/GWAS_diastolic_BOLT/Results/bolt_", "long_PDSR", "_disc.bgen.stats"))
# apply maf > 0.005 and hwe filter
stats = stats %>% filter(A1FREQ <= 0.995)
stats = stats %>% filter(SNP %in% snps_in)
stats = stats %>% dplyr::select(SNP, CHR, BP, P_BOLT_LMM_INF)
stats = dplyr::rename(stats, p_value=P_BOLT_LMM_INF, chr=CHR, pos=BP)
stats$phenotype = NA
stats$gene = NA
stats$phenotype[match(c("rs11970286"), stats$SNP)] = ""
stats$gene[match(c("rs11970286"), stats$SNP)] = c("PLN")
png(filename="data/manhattan_disc_long.png", width=1200, height=350)
manh_plot(stats, pcut_label=1e-5)
dev.off()

stats = read_tsv(paste0("/gpfs01/bhcbio/projects/UK_Biobank/20190102_UK_Biobank_Imaging/Results/GWAS/GWAS_diastolic_BOLT/Results/bolt_", "long_PDSR", "_full.bgen.stats"))
# apply maf > 0.005 and hwe filter
stats = stats %>% filter(A1FREQ <= 0.995)
stats = stats %>% filter(SNP %in% snps_in)
stats = stats %>% dplyr::select(SNP, CHR, BP, P_BOLT_LMM_INF)
stats = dplyr::rename(stats, p_value=P_BOLT_LMM_INF, chr=CHR, pos=BP)
stats$phenotype = NA
stats$gene = NA
stats$phenotype[match(c("rs11970286","rs11535974","rs499715","rs10261575"), stats$SNP)] = ""
stats$gene[match(c("rs11970286","rs11535974","rs499715","rs10261575"), stats$SNP)] = c("PLN","AC023158.1","FHOD3","PHF14")
png(filename="data/manhattan_full_long.png", width=1200, height=350)
manh_plot(stats, pcut_label=1e-5)
dev.off()

stats = read_tsv(paste0("/gpfs01/bhcbio/projects/UK_Biobank/20190102_UK_Biobank_Imaging/Results/GWAS/GWAS_diastolic_BOLT/Results/bolt_", "long_PDSR", "_full_apr_covar.bgen.stats"))
stats = stats %>% dplyr::select(SNP, CHR, BP, P_BOLT_LMM_INF)
stats = dplyr::rename(stats, p_value=P_BOLT_LMM_INF, chr=CHR, pos=BP)
stats$phenotype = NA
stats$gene = NA
stats$phenotype[match(c("rs11535974"), stats$SNP)] = ""
stats$gene[match(c("rs11535974"), stats$SNP)] = c("AC023158.1")
png(filename="data/manhattan_full_long_apr_covar.png", width=1200, height=350)
manh_plot(stats, pcut_label=1e-3)
dev.off()


# HERITABILITY ------------------------------------------------------------

# heritability plot

h2 = data.frame(
  gwas = c("Radial PDSR","Longitdinal PDSR","LAV"),
  estimate = c(0.1416,0.1238,0.2138),
  se = c(0.0148,0.0148,0.0173)
)

ggplot(h2, aes(y=factor(gwas), x=estimate)) + geom_point(size=1) + geom_errorbarh(aes(xmax=estimate+se, xmin=estimate-se), size=0.5, height=0.1, color="gray50") + theme_thesis(15) + xlab("H2") + ylab("") + coord_flip() + scale_x_continuous(limits=c(0,0.5))


# LEAD SNP HEATMAP --------------------------------------------------------

# input = list of lead variants
# res = summary stats for the variants

to_plot = as.matrix(res[,c(3,6,9)])
to_plot = -to_plot
rownames(to_plot) = res$variant

# make_colors <- function(colors, cutoff_fraction, num_colors) {
#   stopifnot(length(colors) == 4)
#   ramp1 = colorRampPalette(colors[1:2])(num_colors * cutoff_fraction)
#   ramp2 = colorRampPalette(colors[3:4])(num_colors * (1-cutoff_fraction))
#   return(c(ramp1, ramp2))
# }

# cutoff_distance = 0  
# cols = make_colors(c("blue","white","white","red"), cutoff_distance/max(to_plot), 100)

p_annot = matrix(paste0("p = ", unlist(res[,c(2,5,8)]), ", beta = ", -unlist(res[,c(3,6,9)])), nrow=9, byrow=FALSE)
colnames(p_annot) = colnames(to_plot)

genes = c("BAG3","FHOD3","PLN","AC023158.1","EFEMP1","NPR3","GJA1","PHF14","SP1")
rownames(to_plot) = paste(rownames(to_plot), genes, sep=", ")
rownames(p_annot) = rownames(to_plot)

palette_length = 100
my_color = colorRampPalette(c("red","white","blue"))(palette_length)
# length(breaks) == length(palette_length) + 1
# use floor and ceiling to deal with even/odd length pallette lengths
my_breaks = c(seq(min(to_plot), 0, length.out=ceiling(palette_length/2)+1), 
              seq(max(to_plot)/palette_length, max(to_plot), length.out=floor(palette_length/2))
)

to_plot = to_plot[,c(1,3,2)]
p_annot = p_annot[,c(1,3,2)]

pheatmap(to_plot, cluster_rows=FALSE, cluster_cols=FALSE, display_numbers=p_annot, color=my_color, breaks=my_breaks, cex=1.2, fontsize_number=12, number_color="black")

