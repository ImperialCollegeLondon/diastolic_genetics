
require(tidyverse)
require(epiChoose)
require(ggrepel)
require(VennDiagram)
require(UKBRlib)
require(devtools)
require(ghql)
require(graphql)
require(rtracklayer)
require(httr)
require(org.Hs.eg.db)
require(readxl)
load_all("~/links/bullseye/")

Sys.setenv(R_CONFIG_ACTIVE="standard")


# SIGNIFICANT HITS --------------------------------------------------------

# /gpfs01/bhcbio/projects/UK_Biobank/20190102_UK_Biobank_Imaging/Results/GWAS/GWAS_diastolic_BOLT/Results
# sig hits
# non-infin
# zcat /gpfs01/bhcbio/projects/UK_Biobank/20190102_UK_Biobank_Imaging/Results/GWAS/GWAS_diastolic_BOLT/Results/bolt_LAV_full.bgen.stats | awk '$16<=5e-8' - > /gpfs01/home/glkwj/projects/diastolic_genetics/genetic_analysis/data/res_bolt_lav_full.txt
# zcat /gpfs01/bhcbio/projects/UK_Biobank/20190102_UK_Biobank_Imaging/Results/GWAS/GWAS_diastolic_BOLT/Results/bolt_long_PDSR_full.bgen.stats | awk '$16<=5e-8' - > /gpfs01/home/glkwj/projects/diastolic_genetics/genetic_analysis/data/res_bolt_long_pdsr_full.txt
# zcat /gpfs01/bhcbio/projects/UK_Biobank/20190102_UK_Biobank_Imaging/Results/GWAS/GWAS_diastolic_BOLT/Results/bolt_radial_PDSR_full.bgen.stats | awk '$16<=5e-8' - > /gpfs01/home/glkwj/projects/diastolic_genetics/genetic_analysis/data/res_bolt_radial_pdsr_full.txt

# zcat /gpfs01/bhcbio/projects/UK_Biobank/20190102_UK_Biobank_Imaging/Results/GWAS/GWAS_diastolic_BOLT/Results/bolt_LAV_disc.bgen.stats | awk '$16<=5e-8' - > /gpfs01/home/glkwj/projects/diastolic_genetics/genetic_analysis/data/res_bolt_lav_disc.txt
# zcat /gpfs01/bhcbio/projects/UK_Biobank/20190102_UK_Biobank_Imaging/Results/GWAS/GWAS_diastolic_BOLT/Results/bolt_long_PDSR_disc.bgen.stats | awk '$16<=5e-8' - > /gpfs01/home/glkwj/projects/diastolic_genetics/genetic_analysis/data/res_bolt_long_pdsr_disc.txt
# zcat /gpfs01/bhcbio/projects/UK_Biobank/20190102_UK_Biobank_Imaging/Results/GWAS/GWAS_diastolic_BOLT/Results/bolt_radial_PDSR_disc.bgen.stats | awk '$16<=5e-8' - > /gpfs01/home/glkwj/projects/diastolic_genetics/genetic_analysis/data/res_bolt_radial_pdsr_disc.txt

# infin
# awk '$14<=5e-8' bolt_lmm_res_bgen > output

header = c("SNP","CHR","BP","GENPOS","ALLELE1","ALLELE0","A1FREQ","INFO","CHISQ_LINREG","P_LINREG","BETA","SE","CHISQ_BOLT_LMM_INF","P_BOLT_LMM_INF","CHISQ_BOLT_LMM","P_BOLT_LMM")
measures = c("lav", "long", "radial")

input_dat = list(
  lav = read_tsv("data/res_bolt_lav_disc.txt", col_names=FALSE),
  lav_full = read_tsv("data/res_bolt_lav_full.txt", col_names=FALSE),
  lav_repl = read_tsv("data/res_bolt_lav_repl_disc_match.txt", col_names=FALSE),
  long = read_tsv("data/res_bolt_long_pdsr_disc.txt", col_names=FALSE),
  long_full = read_tsv("data/res_bolt_long_pdsr_full.txt", col_names=FALSE),
  long_repl = read_tsv("data/res_bolt_long_pdsr_repl_disc_match.txt", col_names=FALSE),
  radial = read_tsv("data/res_bolt_radial_pdsr_disc.txt", col_names=FALSE),
  radial_full = read_tsv("data/res_bolt_radial_pdsr_full.txt", col_names=FALSE),
  radial_repl = read_tsv("data/res_bolt_radial_pdsr_repl_disc_match.txt", col_names=FALSE)
)

for(i in 1:length(input_dat)) names(input_dat[[i]]) = header
lapply(input_dat, function(x) dim(x)[1])

lapply(input_dat[c(2,5,8)], function(x) x %>% group_by(CHR) %>% summarise(N=n(), lower=min(BP), upper=max(BP)))

x = bind_rows(input_dat[c(1,4,7)], .id="GWAS_NAME") # %>% View(title="disc") # disc
y = bind_rows(input_dat[c(3,6,9)], .id="GWAS_NAME") # %>% View(title="repl") # repl
z = bind_rows(input_dat[c(2,5,8)], .id="GWAS_NAME") # %>% View(title="full") # full


# there are 3 hits across the 3 gwas that come up in the disc but not the full
# check them here ...
# added to genetic_analysis/locuszoom/ - look like fps

hits_tc = c(
  "rs13184632",
  "1:178221156_CA_C",
  "rs149751856"
)
hits_tc_meta = getMetadataForImputedSnpsV2(by="snp", hits_tc)
hits_tc_meta$gwas = c("long_PDSR","LAV","radial_PDSR")
win = 1e5

for(i in 1:length(hits_tc_meta$RSID)) {
  
  stats = read_tsv(paste0("/gpfs01/bhcbio/projects/UK_Biobank/20190102_UK_Biobank_Imaging/Results/GWAS/GWAS_diastolic_BOLT/Results/bolt_", hits_tc_meta$gwas[i], "_disc.bgen.stats"))
  stats = stats %>% dplyr::select(SNP, CHR, BP, P_BOLT_LMM_INF, BETA, SE)
  stats$N = 2e4
  stats = dplyr::rename(stats, P=P_BOLT_LMM_INF, ESTIMATE=BETA)
  locus_zoom_region(gwas_res=stats, CHR=hits_tc_meta$Chromosome[i], start=hits_tc_meta$Pos[i]-win, stop=hits_tc_meta$Pos[i]+win, clean=TRUE)
  
}

# replication significant?

lead_disc_snps = bind_rows(input_dat[c(1,4,7)], .id="GWAS_NAME") %>% group_by(GWAS_NAME, X2) %>% dplyr::slice(which.min(X16)) %>% filter(X2!=1, X2!=3)
names(lead_disc_snps)[-1] = header

all_repl = bind_rows(
  list(
    radial = read_tsv("data/res_bolt_radial_pdsr_repl_match.txt", col_names=FALSE),
    long = read_tsv("data/res_bolt_long_pdsr_repl_match.txt", col_names=FALSE),
    lav = read_tsv("data/res_bolt_lav_repl_match.txt", col_names=FALSE)
  ),
  .id="GWAS_NAME"
)
names(all_repl)[-1] = header
all_repl %>% group_by(GWAS_NAME, CHR) %>% dplyr::slice(which.min(P_BOLT_LMM)) %>% filter(CHR!=1, CHR!=3) %>% summarise(p=P_BOLT_LMM<0.05/5, p_value=P_BOLT_LMM)

to_plot = lead_disc_snps %>% left_join(all_repl, by=c("GWAS_NAME","SNP","CHR","BP","ALLELE1","ALLELE0"))
to_plot %>% ggplot(aes(P_BOLT_LMM.x, P_BOLT_LMM.y)) + geom_point() + geom_vline(xintercept=5e-8) + geom_hline(yintercept=0.05/5)
to_plot %>% ggplot(aes(BETA.x, BETA.y)) + geom_point() + geom_hline(yintercept=0) + geom_vline(xintercept=0)

# plot the replication - these are the plots from the teamsite conversation

lead_disc_snps = bind_rows(data.frame(GWAS_NAME=NA), lead_disc_snps)
lead_disc_snps[1,c(1,3)] = lead_disc_snps[2,c(1,3)]
lead_disc_snps$SNP[1] = "rs1173727"
lead_disc_snps$BP[1] = 32830521
lead_disc_snps$gwas_root = c("LAV","long_PDSR","radial_PDSR","radial_PDSR","radial_PDSR")

gwas_tots = data.frame()

for(variant_ix in 1:length(lead_disc_snps$SNP)) {
  
  snp_pos = lead_disc_snps$BP[variant_ix]
  win = 2e5
  
  region_granges = GenomicRanges::GRanges(
    seqnames = lead_disc_snps$CHR[variant_ix], 
    ranges = IRanges::IRanges(start=lead_disc_snps$BP[variant_ix]-win, end=lead_disc_snps$BP[variant_ix]+win), 
    strand = "*")
  region_granges
  
  gwas_disc = read.table(pipe(paste0("awk 'NR==1 {print}; $2==", lead_disc_snps$CHR[variant_ix], " && $3>", lead_disc_snps$BP[variant_ix]-win, " && $3<", lead_disc_snps$BP[variant_ix]+win, " {print}' /gpfs01/bhcbio/projects/UK_Biobank/20190102_UK_Biobank_Imaging/Results/GWAS/GWAS_diastolic_BOLT/Results/bolt_", lead_disc_snps$gwas_root[variant_ix], "_disc.bgen.stats")), header=TRUE)
  gwas_repl = read.table(pipe(paste0("awk 'NR==1 {print}; $2==", lead_disc_snps$CHR[variant_ix], " && $3>", lead_disc_snps$BP[variant_ix]-win, " && $3<", lead_disc_snps$BP[variant_ix]+win, " {print}' /gpfs01/bhcbio/projects/UK_Biobank/20190102_UK_Biobank_Imaging/Results/GWAS/GWAS_diastolic_BOLT/Results/bolt_", lead_disc_snps$gwas_root[variant_ix], "_repl.bgen.stats")), header=TRUE)
  gwas_full = read.table(pipe(paste0("awk 'NR==1 {print}; $2==", lead_disc_snps$CHR[variant_ix], " && $3>", lead_disc_snps$BP[variant_ix]-win, " && $3<", lead_disc_snps$BP[variant_ix]+win, " {print}' /gpfs01/bhcbio/projects/UK_Biobank/20190102_UK_Biobank_Imaging/Results/GWAS/GWAS_diastolic_BOLT/Results/bolt_", lead_disc_snps$gwas_root[variant_ix], "_full.bgen.stats")), header=TRUE)
  
  gwas = rbind(
    cbind(gwas_disc, group="disc"),
    cbind(gwas_repl, group="repl"),
    cbind(gwas_full, group="full")
  )
  # gwas = makeGRangesFromDataFrame(gwas, keep.extra.columns=TRUE, start.field="BP", end.field="BP")
  # gwas = lift_over(gwas, dir="37_to_38")
  gwas_tots = rbind(gwas_tots, cbind(gwas, locus=paste(lead_disc_snps$CHR[variant_ix], lead_disc_snps$SNP[variant_ix], lead_disc_snps$gwas_root[variant_ix], sep="_")))
  
}

gwas_tots_label = gwas_tots %>% group_by(locus, group) %>% summarise(lead_snp=SNP[which.min(P_BOLT_LMM)], P_BOLT_LMM_MIN=P_BOLT_LMM[which.min(P_BOLT_LMM)], BP_MIN=BP[which.min(P_BOLT_LMM)])
gwas_tots_label$repl_sig = NA
gwas_tots_label$repl_sig[gwas_tots_label$group=="repl"] = gwas_tots_label$P_BOLT_LMM_MIN[gwas_tots_label$group=="repl"] < 0.05/5

gwas_tots %>% ggplot(aes(x=BP, y=-log(P_BOLT_LMM, base=10), color=group)) + geom_point(alpha=0.5, size=1) + theme_thesis(15) + ylab("-log10(P)") + xlab("") + geom_hline(yintercept = c(-log(5e-8, base=10), -log(0.05/5, base=10)), lty=2, color="red") + facet_wrap(~locus, scales="free") + geom_text_repel(data=gwas_tots_label, aes(x=BP_MIN, y=-log(P_BOLT_LMM_MIN, base=10), label=lead_snp), fontface="bold", size=3, force=0.5, box.padding=0.5)


# PULSE RATE CONDITIONAL ANALYSIS -----------------------------------------

# /gpfs01/bhcbio/projects/UK_Biobank/20190102_UK_Biobank_Imaging/Results/GWAS/GWAS_diastolic_BOLT/BOLT_*_full_apr_covar.sh

covar = read_delim("/gpfs01/bhcbio/projects/UK_Biobank/20190102_UK_Biobank_Imaging/Data/Derived/GWAS_covariates/UKB_covariates.plink", delim=" ")
quant_traits = readProcessedQuantTraits(instance=0:2, type="all")
quant_traits = quant_traits %>% dplyr::select(SID, pulse_rate_adj)
table(covar$FID %in% quant_traits$SID)
covar$pulse_rate_adj = quant_traits$pulse_rate_adj[match(covar$FID, quant_traits$SID)]
write_delim(covar, file="/gpfs01/bhcbio/projects/UK_Biobank/20190102_UK_Biobank_Imaging/Data/Derived/GWAS_covariates/UKB_covariates_with_apr.plink", delim=" ")

input_dat = list(
  lav_full = read_tsv("data/res_bolt_lav_full_apr_covar.txt", col_names=FALSE),
  long_full = read_tsv("data/res_bolt_long_pdsr_full_apr_covar.txt", col_names=FALSE),
  radial_full = read_tsv("data/res_bolt_radial_pdsr_full_apr_covar.txt", col_names=FALSE)
)
for(i in 1:length(input_dat)) names(input_dat[[i]]) = header
lapply(input_dat, function(x) dim(x)[1])

# lead_snps = lapply(input_dat, function(x) x %>% group_by(CHR) %>% dplyr::slice(which.min(P_BOLT_LMM)))
# lead_snps = bind_rows(lead_snps, .id="GWAS")

radial = read.table("/gpfs01/bhcbio/projects/UK_Biobank/20190102_UK_Biobank_Imaging/Results/GWAS/GWAS_diastolic_BOLT/Results/bolt_radial_PDSR_full_apr_covar.bgen.stats", header=TRUE)
radial_tab <- radial %>% dplyr::rename(P=P_BOLT_LMM)
pdf(file="manhattan_radial_apr_covar.pdf", width=15, height=4)
manhattan_plot(radial_tab, header="Radial Strain Rate (Conditioned on Pulse Rate)")
dev.off()
rm(radial); rm(radial_tab)

long = read.table("/gpfs01/bhcbio/projects/UK_Biobank/20190102_UK_Biobank_Imaging/Results/GWAS/GWAS_diastolic_BOLT/Results/bolt_long_PDSR_full_apr_covar.bgen.stats", header=TRUE)
long_tab <- long %>% dplyr::rename(P=P_BOLT_LMM)
pdf(file="manhattan_long_apr_covar.pdf", width=15, height=4)
manhattan_plot(long_tab, header="Longitudinal Strain Rate (Conditioned on Pulse Rate)")
dev.off()
rm(long); rm(long_tab)

lav = read.table("/gpfs01/bhcbio/projects/UK_Biobank/20190102_UK_Biobank_Imaging/Results/GWAS/GWAS_diastolic_BOLT/Results/bolt_LAV_full_apr_covar.bgen.stats", header=TRUE)
lav_tab <- lav %>% dplyr::rename(P=P_BOLT_LMM)
pdf(file="manhattan_lav_apr_covar.pdf", width=15, height=4)
manhattan_plot(lav_tab, header="LAV (Conditioned on Pulse Rate)")
dev.off()
rm(lav); rm(lav_tab)


# MANHATTAN PLOTS ANNOTATED -----------------------------------------------

# radial

stats = read_tsv(paste0("/gpfs01/bhcbio/projects/UK_Biobank/20190102_UK_Biobank_Imaging/Results/GWAS/GWAS_diastolic_BOLT/Results/bolt_", "radial_PDSR", "_disc.bgen.stats"))
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
stats$phenotype[match(c("rs11535974"), stats$SNP)] = "Longitudinal PDSR"
stats$gene[match(c("rs11535974"), stats$SNP)] = c("AC023158.1")
png(filename="data/manhattan_full_long_apr_covar.png", width=1200, height=350)
manh_plot(stats, pcut_label=1e-5)
dev.off()


# QQ PLOTS ----------------------------------------------------------------

radial = read.table("/gpfs01/bhcbio/projects/UK_Biobank/20190102_UK_Biobank_Imaging/Results/GWAS/GWAS_diastolic_BOLT/Results/bolt_radial_PDSR_full.bgen.stats", header=TRUE)
pdf(file="qq_radial.pdf", width=7, height=6)
qq(radial$P_BOLT_LMM)
dev.off()
rm(radial)

long = read.table("/gpfs01/bhcbio/projects/UK_Biobank/20190102_UK_Biobank_Imaging/Results/GWAS/GWAS_diastolic_BOLT/Results/bolt_long_PDSR_full.bgen.stats", header=TRUE)
pdf(file="qq_long.pdf", width=7, height=6)
qq(long$P_BOLT_LMM)
dev.off()
rm(long)

lav = read.table("/gpfs01/bhcbio/projects/UK_Biobank/20190102_UK_Biobank_Imaging/Results/GWAS/GWAS_diastolic_BOLT/Results/bolt_LAV_full.bgen.stats", header=TRUE)
pdf(file="qq_lav.pdf", width=7, height=6)
qq(lav$P_BOLT_LMM)
dev.off()
rm(lav)

# heritability

h2 = data.frame(
  gwas = c("Radial PDSR","Longitdinal PDSR","LAV"),
  estimate = c(0.1416,0.1238,0.2138),
  se = c(0.0148,0.0148,0.0173)
)

ggplot(h2, aes(y=factor(gwas), x=estimate)) + geom_point(size=1) + geom_errorbarh(aes(xmax=estimate+se, xmin=estimate-se), size=0.5, height=0.1, color="gray50") + theme_thesis(15) + xlab("H2") + ylab("") + coord_flip() + scale_x_continuous(limits=c(0,0.5))


# UK DIGITAL HEART REPLICATION --------------------------------------------

# get the 13 lead snps

lead_snps = lapply(input_dat[c(2,5,8)], function(x) x %>% group_by(CHR) %>% dplyr::slice(which.min(P_BOLT_LMM)))
lead_snps = bind_rows(lead_snps, .id="GWAS")

replace_snp = input_dat$long_full %>% filter(CHR=="18") %>% arrange(P_BOLT_LMM)
lead_snps[8,] = cbind("long_full", replace_snp[2,])
replace_snp = input_dat$radial_full %>% filter(CHR=="6") %>% arrange(P_BOLT_LMM)
lead_snps[10,] = cbind("radial_full", replace_snp[2,])

add_long = c("rs572786197") # 1 extra on chr 4
add_radial = c("rs9269887","rs1580396")

lead_snps = rbind(
  lead_snps,
  z %>% filter(GWAS_NAME=="long_full", SNP %in% add_long),
  z %>% filter(GWAS_NAME=="radial_full", SNP %in% add_radial)
)

lead_snps$GWAS[14:16] = c("long_full","radial_full","radial_full")
lead_snps = lead_snps %>% arrange(GWAS, CHR, BP)

lead_snps$root = c(rep("LAV",2), rep("long_PDSR",7), rep("radial_PDSR",7))
lead_snps_with_junk = lead_snps
lead_snps = lead_snps[c(13,16,5,14,1,2,12,6,15),] # paper order

# annotate for table 2
meta_data = getMetadataForUKBVariants(by="snpid", snpids=lead_snps$SNP)
meta_data_filt = dplyr::select(meta_data, -ukb_genosrc, -imputed_INFO) %>% group_by(across(-MAF_ukb)) %>% summarise(MAF_ukb=max(MAF_ukb))
lead_snps = left_join(lead_snps, meta_data_filt, by=c("SNP"="rsid"))
ft = lead_snps %>% dplyr::select(SNP, CHR, REF_ukb, ALT_ukb, MAF_ukb, BETA, SE, P_BOLT_LMM)
ft$SE = format(round(ft$SE,4), nsmall=4)
ft$BETA = -as.numeric(ft$BETA)
ft$BETA = format(round(ft$BETA,4), nsmall=4)
ft$MAF_ukb = format(round(ft$MAF_ukb,4), nsmall=4)
ft
write_csv(ft, "data/lead_snps_for_table.csv")

win = 1e4
stats_lav = read_tsv("/gpfs01/bhcbio/projects/UK_Biobank/20190102_UK_Biobank_Imaging/Results/GWAS/GWAS_diastolic_BOLT/Results/bolt_LAV_full.bgen.stats")
stats_long = read_tsv("/gpfs01/bhcbio/projects/UK_Biobank/20190102_UK_Biobank_Imaging/Results/GWAS/GWAS_diastolic_BOLT/Results/bolt_long_PDSR_full.bgen.stats")
stats_radial = read_tsv("/gpfs01/bhcbio/projects/UK_Biobank/20190102_UK_Biobank_Imaging/Results/GWAS/GWAS_diastolic_BOLT/Results/bolt_radial_PDSR_full.bgen.stats")

res_all = vector("list", length(lead_snps$GWAS))
names(res_all) = lead_snps$SNP

for(lead_ix in 10:length(lead_snps$GWAS)) {
  
  print(paste("Variant:", lead_ix))
  # lead_ix = 1
  
  # get the list of snps in ld > 0.8
  
  if(lead_snps$GWAS[lead_ix]=="lav_full") {
    background_snps = stats_lav %>% dplyr::filter(CHR==lead_snps$CHR[lead_ix], BP > lead_snps$BP[lead_ix]-win, BP < lead_snps$BP[lead_ix]+win) %>% dplyr::select(SNP)
  }
  if(lead_snps$GWAS[lead_ix]=="long_full") {
    background_snps = stats_long %>% dplyr::filter(CHR==lead_snps$CHR[lead_ix], BP > lead_snps$BP[lead_ix]-win, BP < lead_snps$BP[lead_ix]+win) %>% dplyr::select(SNP)
  }
  if(lead_snps$GWAS[lead_ix]=="radial_full") {
    background_snps = stats_radial %>% dplyr::filter(CHR==lead_snps$CHR[lead_ix], BP > lead_snps$BP[lead_ix]-win, BP < lead_snps$BP[lead_ix]+win) %>% dplyr::select(SNP)
  }
  
  res_all[[lead_ix]] = search_proxies(snp=lead_snps$SNP[lead_ix], background_snps=background_snps$SNP, ld_thres=0.8)
  print(res_all[[lead_ix]])
  
}

to_send = bind_rows(input_dat[c(2,5,8)], .id="GWAS_NAME")


# DISCOVERY VERSUS REPLICATION --------------------------------------------

# is the directionality repeated in the validation set?
# bash script for filtering *repl data by *disc hits
# match_disc_repl.sh

for(i in 1:length(measures)) {
  discovery = input_dat[[which(names(input_dat)==paste0(measures[i], "_full"))]]
  replication = input_dat[[which(names(input_dat)==paste0(measures[i], "_repl"))]]
  print(data.frame(Discovery=discovery$BETA, Replication=replication$BETA) %>% ggplot(aes(Discovery, Replication)) + geom_point() + theme_thesis(20) + ggtitle(measures[i]))
}


# DISC REPL CHECK ---------------------------------------------------------

roots = c("lav","long_pdsr","radial_pdsr")
options(stringsAsFactors=FALSE)

hits = list(
  lav = c("rs59985551","rs1173727","rs35489511"),
  long = c("rs2275950","rs11970286","rs10261575","rs11535974","rs499715"),
  radial = c("rs528236848","rs9388001","rs2234962","rs11170519","rs369533272")
)

hits_all = vector("list", length(hits))
names(hits_all) = roots

for(i in 1:length(hits_all)) {
  
  full = read_tsv(paste0("data/res_bolt_", roots[i], "_full.txt"), col_names=FALSE)
  disc = read_tsv(paste0("data/res_bolt_", roots[i], "_disc_match.txt"), col_names=FALSE)
  repl = read_tsv(paste0("data/res_bolt_", roots[i], "_repl_match.txt"), col_names=FALSE)
  
  hits_all[[i]] = data.frame(
    hits = hits[[i]],
    full_p = full$X16[match(hits[[i]], full$X1)],
    full_e = full$X11[match(hits[[i]], full$X1)],
    disc_p = disc$X16[match(hits[[i]], disc$X1)],
    disc_e = disc$X11[match(hits[[i]], disc$X1)],
    repl_p = repl$X16[match(hits[[i]], repl$X1)],
    repl_e = repl$X11[match(hits[[i]], repl$X1)]
  )
}

hits_all = bind_rows(hits_all, .id="id")
hits_all
hits_all %>% dplyr::select(c(1,2,3,5,7)) %>% gather("Group","P",3:5) %>% ggplot(aes(fill=Group, y=-log(P,base=10), x=hits)) + geom_bar(position="dodge", stat="identity") + theme_thesis(20) + facet_wrap(~id, scales="free") + geom_hline(yintercept=-log(5e-8, base=10), linetype="dashed", color="red", size=1) + ylab("P")

hits_all %>% dplyr::select(c(1,2,4,6,8)) %>% gather("Group","Beta",3:5) %>% ggplot(aes(fill=Group, y=Beta, x=hits)) + geom_bar(position="dodge", stat="identity") + theme_thesis(20) + facet_wrap(~id, scales="free")


# CONDITIONAL ANALYSIS ----------------------------------------------------

require(UKBRlib)
require(snpStats)
require(jsonlite)
require(httr)

Sys.setenv(R_CONFIG_ACTIVE="standard")

# need snp, a1, a2, freq, beta, se, p, n
# for full gwas

roots = c("radial_PDSR","long_PDSR","LAV")

for(i in 1:length(roots)) {
  
  stats = read_tsv(paste0("/gpfs01/bhcbio/projects/UK_Biobank/20190102_UK_Biobank_Imaging/Results/GWAS/GWAS_diastolic_BOLT/Results/bolt_", roots[i], "_full.bgen.stats"))
  stats = stats %>% dplyr::select(CHR, SNP, ALLELE0, ALLELE1, BETA, SE, P_BOLT_LMM) %>% dplyr::rename(A2=ALLELE1, A1=ALLELE0, b=BETA, se=SE, p=P_BOLT_LMM) %>% mutate(N=39559)
  
  # get snp stats by chromosome
  
  # chrs = unique(stats$CHR)
  # snp_chrs = vector("list", length(chrs))
  # names(snp_chrs) = chrs
  # for(j in 1:length(chrs)) {
  #   print(paste("Chr:", chrs[j]))
  #   snp_chrs[[j]] = getMetadataForImputedSnps(stats %>% filter(CHR==chrs[j]) %>% dplyr::select(SNP) %>% unlist() %>% as.character() %>% unique()) %>% dplyr::select(-SNP,-Info) %>% dplyr::rename(SNP=RSID, A1M=A1, A2M=A2)
  # }
  # snp_chrs_df = bind_rows(snp_chrs, .id="CHR")  
  
  # join freqs to df
  
  stats = stats %>% dplyr::left_join(snp_chrs_df %>% dplyr::select(SNP, MAF, MA), by="SNP") %>% dplyr::mutate(MAF=dplyr::if_else(A1==MA, MAF, 1-MAF)) 
  stats = stats %>% dplyr::select(SNP, A1, A2, MAF, b, se, p, N) %>% dplyr::rename(freq=MAF) %>% dplyr::filter(!is.na(freq))
  
  prior_snps = input_dat$lav_full %>% dplyr::select(SNP, CHR)
  
  dat_ca = conditional_analysis(ma_file=stats, prior_snps=prior_snps, dir="~/projects/diastolic_genetics/cojo", thresh=5e-8)
  
}


# LEAD SNP HEATMAP --------------------------------------------------------

lead_snps
roots = c("LAV","long_PDSR","radial_PDSR")
res = data.frame(variant=lead_snps$SNP)

for(i in 1:length(roots)) {
  
  stats = read_tsv(paste0("/gpfs01/bhcbio/projects/UK_Biobank/20190102_UK_Biobank_Imaging/Results/GWAS/GWAS_diastolic_BOLT/Results/bolt_", roots[i], "_full.bgen.stats"))
  dat = stats %>% filter(SNP %in% lead_snps$SNP)
  dat = dat[match(lead_snps$SNP, dat$SNP),]
  to_add = data.frame(dat$P_BOLT_LMM, dat$BETA, dat$SE)  
  names(to_add) = c(paste("p", roots[i], sep="_"), paste("beta", roots[i], sep="_"), paste("se", roots[i], sep="_"))
  res = cbind(res, to_add)
}

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

p_annot = matrix(paste0("p = ", unlist(res[,c(2,5,8)]), ", beta = ", -unlist(res[,c(3,6,9)])), nrow=8, byrow=FALSE)
colnames(p_annot) = colnames(to_plot)

genes = c("BAG3","FHOD3","PLN","AC023158.1","EFEMP1","NPR3","GJA1","PHF14")
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

png(filename="data/heatmap.png", height=785, width=1060)
pheatmap(to_plot, cluster_rows=FALSE, cluster_cols=FALSE, display_numbers=p_annot, color=my_color, breaks=my_breaks, cex=1.2, fontsize_number=12, number_color="black")
dev.off()


# LOCUS ZOOM --------------------------------------------------------------

source("genetic_analysis/code/locus_zoom_region_temp.R")
roots = c("LAV","long_PDSR","radial_PDSR")
win = 1e5

for(i in 1:length(roots)) {
  
  stats = read_tsv(paste0("/gpfs01/bhcbio/projects/UK_Biobank/20190102_UK_Biobank_Imaging/Results/GWAS/GWAS_diastolic_BOLT/Results/bolt_", roots[i], "_full.bgen.stats"))
  stats = stats %>% dplyr::select(SNP, CHR, BP, P_BOLT_LMM_INF, BETA, SE)
  stats$N = 2e4
  stats = dplyr::rename(stats, P=P_BOLT_LMM_INF, ESTIMATE=BETA)
  
  loci = input_dat[[which(names(input_dat)==paste0(measures[i], "_full"))]] %>% group_by(CHR) %>% summarise(N=n(), start=min(BP)-win, end=min(BP)+win)
  # cutoff = 5e7 # for > 1 peak on same chr
  # loci = rbind(
  #   input_dat[[which(names(input_dat)==paste0(measures[i], "_full"))]] %>% filter(CHR!=6) %>% group_by(CHR) %>% summarise(N=n(), start=min(BP)-win, end=max(BP)+win),
  #   input_dat[[which(names(input_dat)==paste0(measures[i], "_full"))]] %>% filter(CHR==6, BP<cutoff) %>% group_by(CHR) %>% summarise(N=n(), start=min(BP)-win, end=max(BP)+win),
  #   input_dat[[which(names(input_dat)==paste0(measures[i], "_full"))]] %>% filter(CHR==6, BP>cutoff) %>% group_by(CHR) %>% summarise(N=n(), start=min(BP)-win, end=max(BP)+win)
  # )
  
  ix = 4
  locus_zoom_region_temp(gwas_res=stats, CHR=lead_snps_with_junk$CHR[ix], start=lead_snps_with_junk$BP[ix]-win, stop=lead_snps_with_junk$BP[ix]+win)
  
}


# PASCAL ------------------------------------------------------------------

require(fgsea)
require(reshape2)
require(readxl)
require(epiChoose)

mapping = read_excel("~/projects/mapping/HumanGeneList_17Sep2018_workup_betterensembl_list.xlsx")
roots = c("lav","long_pdsr","radial_pdsr")

for(i in 1:length(roots)) {
  
  ps = read_tsv(paste0("~/projects/pascal_hpc/PASCAL/output/bolt_", roots[i], "_full.PathwaySet--msigBIOCARTA_KEGG_REACTOME--sum.txt"))
  gs = read_tsv(paste0("~/projects/pascal_hpc/PASCAL/output/bolt_", roots[i], "_full.sum.genescores.txt"))
  gs_fused = read_tsv(paste0("~/projects/pascal_hpc/PASCAL/output/bolt_", roots[i], "_full.sum.fusion.genescores.txt"))
  gs = arrange(gs, pvalue)
  gs
  
  # grab the pathway lists and merge all the data
  gene_sets = gmtPathways("~/projects/pascal_hpc/PASCAL/resources/genesets/msigdb/msigBIOCARTA_KEGG_REACTOME.gmt")
  gene_sets = melt(gene_sets) %>% tbl_df()
  names(gene_sets) = c("entrez","pathway")
  gene_sets$gene = mapping$Symbol[match(gene_sets$entrez, mapping$EntrezGeneID)]
  gene_sets_merged = left_join(gene_sets, gs, by=c("gene"="gene_symbol"))
  gene_sets_merged = left_join(gene_sets_merged, ps, by=c("pathway"="Name"))
  
  # the distribution of pathway p-values
  ps %>% ggplot(aes(x=chi2Pvalue)) + geom_histogram(binwidth=0.01) + xlab("Chi-squared P-value") + theme_thesis(20) + ylab("")
  
  # the top pathway hit
  # TODO is the p-value corrected?
  head(ps)
  gene_sets_merged %>% filter(pathway %in% ps$Name[1:4]) %>% group_by(gene) %>% summarise(N=n()) %>% arrange(desc(N))
  gene_sets_merged %>% filter(pathway==ps$Name[1]) %>% arrange(pvalue)
  
  # select the top/bottom 10 pathways
  my_ps = c(
    ps %>% head(10) %>% select(Name) %>% unlist() %>% as.character(),
    ps %>% tail(10) %>% select(Name) %>% unlist() %>% as.character()
  )
  
  # plot the gene scores for these pathways
  gene_sets_merged %>% filter(pathway %in% my_ps) %>% ggplot(aes(x=factor(pathway, levels=my_ps), y=pvalue)) + geom_boxplot() + theme_thesis(10) + xlab("") + ylab("")
  
  # pathways ranked by percentage of significant gene scores
  gene_sets_merged %>% group_by(pathway) %>% summarise(perc_sig_genes=sum(pvalue<0.05)/n(), pathway_n=n()) %>% arrange(desc(perc_sig_genes))
  
  # pathways ranked by total significant gene scores
  gene_sets_merged %>% group_by(pathway) %>% summarise(perc_sig_genes=sum(pvalue<0.05), pathway_n=n()) %>% arrange(desc(perc_sig_genes))
  
  my_pathway = ps$Name[1]
  
  # results for pathway
  gene_sets_merged %>% filter(pathway==my_pathway) %>% select(gene, pvalue)
  grep("HIST1H1", gs_fused$gene_symbol, value=TRUE)
  
  # what are the significant gene scores for this pathway?
  gene_sets_merged %>% filter(pathway==my_pathway, pvalue<0.05) %>% arrange(gene) %>% select(gene, pvalue)
  
  # for some of these pathways, what is the breakdown between sig and sub-sig genes?
  # could different gene scores be used?
  
}


# GARFIELD ----------------------------------------------------------------

my_measure = measures[1]

# example
setwd("garfield/")
garfield.run(
  paste0(my_measure, ".output"), 
  data.dir="garfield-data",
  trait=my_measure,
  run.option="prep",
  chrs=c(1:22),
  exclude=c(895,975,976,977,978,979,980)
)

garfield.run(
  paste0(my_measure, ".output"),
  data.dir="garfield-data",
  run.option="perm",
  nperm=100000,
  thresh=c(0.1,0.01,0.001,1e-04,1e-05,1e-06,1e-07,1e-08),
  pt_thresh = c(1e-05, 1e-06, 1e-07, 1e-08), maf.bins=5, tags.bins=5,
  tss.bins=5, prep.file=paste0(my_measure, ".output.prep"), optim_mode=TRUE, minit=100, thresh_perm=0.0001
)

garfield.plot(paste0(my_measure, ".output.perm"), num_perm=100000, output_prefix=paste0(my_measure, ".output"), plot_title=my_measure, filter=10, tr=-log10(0.05/498))


# HARDY WEINBERG CHECK ----------------------------------------------------

require(tidyverse)
require(UKBRlib)
Sys.setenv(R_CONFIG_ACTIVE="imaging")

# code from /gpfs01/bhcbio/projects/UK_Biobank/20190102_UK_Biobank_Imaging/Programs/OtherAnalysis/DiastolicFunction

check <- getTibbleFromSnpMatrix(getImputedV3GenotypesForVariants("rs11111", excluded=c()))

phenotypes <- readUKBiobankPheno(c(22006,22001,22027,22019,31,22005))

relateness <- data.table::fread(config::get("sample_relateness_file"))

# non-white british subset
getEuropean <- data.table::fread("/gpfs01/bhcbio/projects/UK_Biobank/20190102_UK_Biobank_Imaging/Data/Derived/European_samples_filtered_by_HapMapIII_CGRCh37_non_redacted_18545_40616.txt")$eid_40616

# genetic sex mismatch
sex_mismatch <- phenotypes$SID[which(phenotypes$x22001_0_0!=phenotypes$x31_0_0)]
length(sex_mismatch)

# heterogeneity outlier
het_outlier <- phenotypes$SID[which(phenotypes$x22027_0_0==1)]
length(het_outlier)

# sex_chrom_aneuploidy
sex_chrom <- phenotypes$SID[which(phenotypes$x22019_0_0==1)]
length(sex_chrom)

# relateness
related <- unique(relateness$ID1[which(relateness$Kinship>0.0884)])
length(related)

# total number of excluded subjects
nrow(check) - length(unique(c(sex_mismatch, het_outlier, sex_chrom, related)))

# number of europeans that passed qc
getEuropeanQC <- getEuropean[which(!(getEuropean) %in% c(sex_mismatch, het_outlier, sex_chrom, related))]
length(getEuropeanQC)

# compare with subjects included in MRI study
mri_subjects <- data.table::fread("/gpfs01/bhcbio/projects/UK_Biobank/20190102_UK_Biobank_Imaging/Data/Original/Phenotypes/ImperialPhenos/clinical_measures_39k_collated_11052020.csv") %>% dplyr::pull(eid_40616)
mri_subjects_2 <- data.table::fread("/gpfs01/bhcbio/projects/UK_Biobank/20190102_UK_Biobank_Imaging/Data/Original/Phenotypes/ImperialPhenos/clinical_measures_27k_collated_qced_04032020.csv")%>% dplyr::pull(eid_40616)

# exclude subjects with qc issues (in full release, but not in release 1-2, mail maria) plus subjects with genetic qc issues
# plus non-europeans
qc_issues_exclude <- phenotypes$SID[which((phenotypes$SID %in% mri_subjects_2) & !(phenotypes$SID %in% mri_subjects))]
mri_subjects_2 <- mri_subjects_2[which((mri_subjects_2 %in% getEuropeanQC))]
mri_subjects <- mri_subjects[which((mri_subjects %in% getEuropeanQC))]

# full imaging set
length(mri_subjects)

## first two releases imaging set
length(mri_subjects_2)

# last release imaging set
mri_subjects_lr <- mri_subjects[which(!(mri_subjects %in% mri_subjects_2))]
length(mri_subjects_lr)

# non imaging samples
length(which(!(phenotypes$SID %in% mri_subjects) & (phenotypes$SID %in% (unique(c(sex_mismatch, het_outlier, sex_chrom, related))))))

# non imaging samples wb
length(which(!(phenotypes$SID %in% mri_subjects) & (phenotypes$SID %in% getEuropeanQC)))

# define imaging and non-imaging subset
imaging <- mri_subjects
non_imaging <- phenotypes$SID[which(!(phenotypes$SID %in% mri_subjects) & (phenotypes$SID %in% getEuropeanQC))] 

# check the lead snps first
lead_snps_with_junk
hw_dat = data.frame()
for(i in 2:dim(lead_snps_with_junk)[1]) {
  my_chr = lead_snps_with_junk$CHR[i]
  my_snp = lead_snps_with_junk$SNP[i]
  dat = read.table(pipe(paste0("awk 'NR==1 {print}; $1==", my_chr, " && $2==\"", my_snp, "\" {print}' /gpfs01/bhcbio/projects/UK_Biobank/20170616_UK_Biobank_Genotyping/Data/Derived/Genotypes/freq_hardy/_001_ukb_cal_chr", my_chr, "_v2.hwe")), header=TRUE)
  if(dim(dat)[1]==0) next
  hw_dat = bind_rows(hw_dat, dat)
}

# this will only work for genotyped snps ...

geno_data = getImputedV3GenotypesForVariants(ids=lead_snps_with_junk$SNP)
hw_in = col.summary(geno_data)
geno_data_df = getTibbleFromSnpMatrix(geno_data)
table(geno_data_df$rs2234962)
apply(geno_data_df[,-1], 2, function(x) HWChisq(as.numeric(table(x))))

binary_traits = readProcessedBinaryTraits(input_subset="all", coreOnly=FALSE, renaming=TRUE)
binary_traits = binary_traits[,c(1,which(names(binary_traits) %in% ""))]

# binary_traits[,2:dim(binary_traits)[2]] = binary_traits[,2:dim(binary_traits)[2]] + 1
# apply(binary_traits[,-1], 2, table)

# for(i in c(6:8)) {
#   print(paste("Phenotype:", names(binary_traits)[i]))
#   binary_traits_edit = binary_traits[,c(1,i)]
#   print(table(binary_traits_edit[,2]))
#   runGWAS(data=binary_traits_edit, queue="medium_priority", explanation="MAP4K1 target assessment")
# }

# run plink again to get full hwe stats

template = readLines("/gpfs01/bhcbio/projects/UK_Biobank/20190102_UK_Biobank_Imaging/Results/GWAS/GWAS_radial_PDSR/cmds/gwas_chr4_hwe.sh")

for(i in c(2:3,5:22)) {
  
  print(paste("Writing File:", i))
  template_copy = template
  template_copy[2] = str_replace(template[2], "4", as.character(i))
  template_copy[11] = str_replace(template[11], "chr4", paste0("chr",as.character(i)))
  template_copy[24] = str_replace(template[24], "chr4", paste0("chr",as.character(i)))
  file_conn <- file(paste0("/gpfs01/bhcbio/projects/UK_Biobank/20190102_UK_Biobank_Imaging/Results/GWAS/GWAS_radial_PDSR/cmds/gwas_chr", i, "_hwe.sh"))
  writeLines(template_copy, file_conn)
  close(file_conn)
}

bash_script <- "/gpfs01/bhcbio/projects/UK_Biobank/20190102_UK_Biobank_Imaging/Results/GWAS/GWAS_radial_PDSR/cmds/run_hwe.sh"
cat("#!/bin/bash", file=bash_script, sep="\n")

for(i in c(2:3,5:22)) {
  cat(paste0("bash gwas_chr", i, "_hwe.sh"), file=bash_script, append=TRUE, sep="\n")
}

# check the lead snps first
lead_snps_with_junk
hw_dat_2 = data.frame()
for(i in 12:dim(lead_snps_with_junk)[1]) {
  print(i)
  my_chr = lead_snps_with_junk$CHR[i]
  my_snp = lead_snps_with_junk$SNP[i]
  dat = read.table(pipe(paste0("awk 'NR==1 {print}; $1==", my_chr, " && $2==\"", my_snp, "\" {print}' /gpfs01/bhcbio/projects/UK_Biobank/20190102_UK_Biobank_Imaging/Results/GWAS/GWAS_radial_PDSR/chunks/gwas_chr", my_chr, "_hwe.hardy")), header=FALSE)
  if(dim(dat)[1]==0) next
  dat$V3 = as.character(dat$V3)
  dat$V4 = as.character(dat$V4)
  hw_dat_2 = bind_rows(hw_dat_2, dat)
}

# get ld for the 4 snps

ix = c(3,4,7,10)
res_all = vector("list", length(ix))
names(res_all) = lead_snps_with_junk$SNP[ix]
mc = 1

for(i in ix) {
  print(mc)  
  res = getMetadataForUKBVariants(by="region", chr=lead_snps_with_junk$CHR[i], pos_start=lead_snps_with_junk$BP[i]-1000, pos_end=lead_snps_with_junk$BP[i]+1000)
  res_all[[mc]] = get_ld(query_snps=lead_snps_with_junk$SNP[i], target_snps=res$rsid, name_qtl="test")
  mc = mc+1
}

# read data
snps_in <- c()
for (i in 1:22) {
  print(i)
  dat <- data.table::fread(paste0("/gpfs01/bhcbio/projects/UK_Biobank/20190102_UK_Biobank_Imaging/Results/GWAS/GWAS_radial_PDSR/chunks/gwas_chr", i, "_hwe.hardy")) %>% dplyr::filter(P>0.05) %>% dplyr::pull(ID)
  snps_in <- c(snps_in, dat)
}


