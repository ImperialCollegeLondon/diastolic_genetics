
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

Sys.setenv(R_CONFIG_ACTIVE="imaging")


# UK DIGITAL HEART REPLICATION --------------------------------------------

# get the 13 lead snps

lead_snps = lapply(input_dat[c(2,5,8)], function(x) x %>% group_by(CHR) %>% dplyr::slice(which.min(P_BOLT_LMM)))
lead_snps = bind_rows(lead_snps, .id="GWAS")

replace_snp = input_dat$long_full %>% filter(CHR=="18") %>% arrange(P_BOLT_LMM)
lead_snps[8,] = cbind("lav_full", replace_snp[2,])

roots = c("LAV","long_PDSR","radial_PDSR")
win = 1e5
stats = read_tsv(paste0("/gpfs01/bhcbio/projects/UK_Biobank/20190102_UK_Biobank_Imaging/Results/GWAS/GWAS_diastolic_BOLT/Results/bolt_", roots[3], "_full.bgen.stats"))
lead_ix = 1

# get the list of snps in ld > 0.8
background_snps = stats %>% dplyr::filter(CHR==lead_snps$CHR[lead_ix], BP > lead_snps$BP[lead_ix]-win, BP < lead_snps$BP[lead_ix]+win) %>% dplyr::select(SNP)

search_proxies(lead_snps$SNP[1])


# SIGNIFICANT HITS --------------------------------------------------------

# /gpfs01/bhcbio/projects/UK_Biobank/20190102_UK_Biobank_Imaging/Results/GWAS/GWAS_diastolic_BOLT/Results
# sig hits
# non-infin
# awk '$16<=5e-8' bolt_lmm_res_bgen > output
# infin
# awk '$14<=5e-8' bolt_lmm_res_bgen > output

header = c("SNP","CHR","BP","GENPOS","ALLELE1","ALLELE0","A1FREQ","INFO","CHISQ_LINREG","P_LINREG","BETA","SE","CHISQ_BOLT_LMM_INF","P_BOLT_LMM_INF","CHISQ_BOLT_LMM","P_BOLT_LMM")
measures = c("lav", "long", "radial")

input_dat = list(
  lav = read_tsv("data/res_bolt_lav_disc.txt", col_names=FALSE),
  lav_full = read_tsv("data/res_bolt_lav_full.txt", col_names=FALSE),
  lav_repl_match = read_tsv("data/res_bolt_lav_repl_match.txt", col_names=FALSE),
  long = read_tsv("data/res_bolt_long_pdsr_disc.txt", col_names=FALSE),
  long_full = read_tsv("data/res_bolt_long_pdsr_full.txt", col_names=FALSE),
  long_repl_match = read_tsv("data/res_bolt_long_pdsr_repl_match.txt", col_names=FALSE),
  radial = read_tsv("data/res_bolt_radial_pdsr_disc.txt", col_names=FALSE),
  radial_full = read_tsv("data/res_bolt_radial_pdsr_full.txt", col_names=FALSE),
  radial_repl_match = read_tsv("data/res_bolt_radial_pdsr_repl_match.txt", col_names=FALSE)
)

for(i in 1:length(input_dat)) names(input_dat[[i]]) = header
lapply(input_dat, function(x) dim(x)[1])

lapply(input_dat[c(2,5,8)], function(x) x %>% group_by(CHR) %>% summarise(N=n(), lower=min(BP), upper=max(BP)))


# DISCOVERY VERSUS REPLICATION --------------------------------------------

# is the directionality repeated in the validation set?
# bash script for filtering *repl data by *disc hits
# match_disc_repl.sh

for(i in 1:length(measures)) {
  
  discovery = input_dat[[which(names(input_dat)==paste0(measures[i], "_full"))]]
  replication = input_dat[[which(names(input_dat)==paste0(measures[i], "_repl_match"))]]
  print(data.frame(Discovery=discovery$BETA, Replication=replication$BETA) %>% ggplot(aes(Discovery, Replication)) + geom_point() + theme_thesis(20) + ggtitle(measures[i]))
}


# CONDITIONAL ANALYSIS ----------------------------------------------------

require(UKBRlib)
require(snpStats)
require(jsonlite)
require(httr)

Sys.setenv(R_CONFIG_ACTIVE="imaging")

# need snp, a1, a2, freq, beta, se, p, n
# for all alleles in locus?

# pick a dataset

for(i in 1:length(measures)) {
  
  dat = input_dat[[which(names(input_dat)==paste0(measures[i], "_full"))]]
  
  locus_plots = dat %>% group_by(CHR) %>% do(plots = ggplot(data=., aes(x=BP, y=-log(P_BOLT_LMM, base=10))) + geom_point() + geom_text_repel(aes(label=SNP), fontface="bold", size=4, force=0.5) + ggtitle(.$CHR[1]))
  locus_plots$plots
  dat %>% group_by(CHR) %>% dplyr::slice(which.min(P_BOLT_LMM)) %>% View()
  
  dat_filt = dat %>% dplyr::select(SNP, ALLELE1, ALLELE0, BETA, SE, P_BOLT_LMM) %>% dplyr::rename(A1=ALLELE1, A2=ALLELE0, b=BETA, se=SE, p=P_BOLT_LMM) %>% mutate(N=2e4)
  head(dat_filt)
  
  # get snp stats
  dat_meta = getMetadataForImputedSnps(dat_filt$SNP) %>% dplyr::select(-SNP,-Info) %>% dplyr::rename(SNP=RSID, A1M=A1, A2M=A2)
  
  dat_comb = dplyr::left_join(dat_filt, dat_meta, by="SNP") %>% dplyr::mutate(MAF=dplyr::if_else(A1==MA, MAF, 1-MAF)) 
  dat_comb = dat_comb %>% dplyr::select(SNP, A1, A2, MAF, b, se, p, N) %>% dplyr::rename(freq=MAF) %>% dplyr::filter(!is.na(freq))
  
  source("R/conditional_analysis.R")
  
  chrs = unique(dat$CHR)
  for(j in 1:length(chrs)) {
    dat_ca = conditional_analysis(dat_comb, dir="/gpfs01/bhcbio/projects/UK_Biobank/20170616_UK_Biobank_Genotyping/Programs/TargetAnalyses/Tie2/test", thres=5e-8, chr=chrs[j])
  }
  
  # test_ma = read_tsv("data/madata_example")
  # dat_ca = conditional_analysis(test_ma, dir="/gpfs01/bhcbio/projects/UK_Biobank/20170616_UK_Biobank_Genotyping/Programs/TargetAnalyses/Tie2/test", thres=5e-8, chr=9)
  
}

lmm %>% group_by(X2) %>% dplyr::slice(order(X16)[1:2])
lmm %>% group_by(X2) %>% summarise(start=min(X3), end=max(X3))


# LOCUS ZOOM --------------------------------------------------------------

source("R/locus_zoom_region.R")
source("R/locus_zoom_region_temp.R")
roots = c("LAV","long_PDSR","radial_PDSR")
win = 1e5

for(i in 1:length(roots)) {
  
  stats = read_tsv(paste0("/gpfs01/bhcbio/projects/UK_Biobank/20190102_UK_Biobank_Imaging/Results/GWAS/GWAS_diastolic_BOLT/Results/bolt_", roots[i], "_full.bgen.stats"))
  stats = stats %>% select(SNP, CHR, BP, P_BOLT_LMM_INF, BETA, SE)
  stats$N = 2e4
  stats = dplyr::rename(stats, P=P_BOLT_LMM_INF, ESTIMATE=BETA)
  
  loci = input_dat[[which(names(input_dat)==paste0(measures[i], "_full"))]] %>% group_by(CHR) %>% summarise(N=n(), start=min(BP)-win, end=max(BP)+win)
  # cutoff = 5e7 # for > 1 peak on same chr
  # loci = rbind(
  #   input_dat[[which(names(input_dat)==paste0(measures[i], "_full"))]] %>% filter(CHR!=6) %>% group_by(CHR) %>% summarise(N=n(), start=min(BP)-win, end=max(BP)+win),
  #   input_dat[[which(names(input_dat)==paste0(measures[i], "_full"))]] %>% filter(CHR==6, BP<cutoff) %>% group_by(CHR) %>% summarise(N=n(), start=min(BP)-win, end=max(BP)+win),
  #   input_dat[[which(names(input_dat)==paste0(measures[i], "_full"))]] %>% filter(CHR==6, BP>cutoff) %>% group_by(CHR) %>% summarise(N=n(), start=min(BP)-win, end=max(BP)+win)
  # )
  
  locus_zoom_region(gwas_res=stats, CHR=loci$CHR, start=loci$start, stop=loci$end)
  
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


# DISC REPL CHECK ---------------------------------------------------------

roots = c("LAV","long_PDSR","radial_PDSR")
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
hits_all %>% select(c(1,2,3,5,7)) %>% gather("Group","P",3:5) %>% ggplot(aes(fill=Group, y=-log(P,base=10), x=hits)) + geom_bar(position="dodge", stat="identity") + theme_thesis(20) + facet_wrap(~id, scales="free") + geom_hline(yintercept=-log(5e-8, base=10), linetype="dashed", color="red", size=1) + ylab("P")

hits_all %>% select(c(1,2,4,6,8)) %>% gather("Group","Beta",3:5) %>% ggplot(aes(fill=Group, y=Beta, x=hits)) + geom_bar(position="dodge", stat="identity") + theme_thesis(20) + facet_wrap(~id, scales="free")


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


# V2G ---------------------------------------------------------------------

mapping = read_excel("data/HumanGeneList_17Sep2018_workup_betterensembl_list.xlsx")

hits = data.frame(
  variant = c("rs59985551","rs1173727","rs35489511","rs2275950","rs11970286","rs10261575","rs11535974","rs499715","rs528236848","rs9388001","rs2234962","rs11170519","rs369533272"),
  gwas = c(rep("lav",3), rep("long",5), rep("radial",5)),
  gwas_root = c(rep("LAV",3), rep("long_PDSR",5), rep("radial_PDSR",5))
)

hits_annot = getMetadataForGenotypedSnps(hits$variant)
hits_ot = bind_rows(lapply(hits$variant, get_coordinates_from_rsid))
hits = left_join(hits, hits_ot, by=c("variant"="rsId"))
hits$pos_38 = as.numeric(str_replace(hits$id, "^.*_([0-9]+)_.*$", "\\1"))
hits$chr = as.numeric(str_extract(hits$id, "^[0-9]+"))
hits$ref = str_replace(hits$id, "^.*([A-Z]+)_[A-Z]+$", "\\1")
hits$alt = str_replace(hits$id, "^.*[A-Z]+_([A-Z]+)$", "\\1")
hits$closest_gene = NA
hits_gr = makeGRangesFromDataFrame(hits, keep.extra.columns=FALSE, start.field="pos_38", end.field="pos_38")
hits_lo = lift_over(hits_gr, dir="38_to_37")
hits$pos_37 = hits_lo$start

# http://htmlpreview.github.io/?https://github.com/kauralasoo/eQTL-Catalogue-resources/blob/master/scripts/tabix_use_case.html
source("code/import_eqtl_catalog.R")
source("code/scan_tabix_df.R")

tabix_paths = read.delim("https://raw.githubusercontent.com/eQTL-Catalogue/eQTL-Catalogue-resources/master/tabix/tabix_ftp_paths.tsv", sep="\t", header=TRUE, stringsAsFactors=FALSE) %>% dplyr::as_tibble()
imported_tabix_paths = read.delim("https://raw.githubusercontent.com/eQTL-Catalogue/eQTL-Catalogue-resources/master/tabix/tabix_ftp_paths_imported.tsv", sep="\t", header=TRUE, stringsAsFactors=FALSE) %>% dplyr::as_tibble()

eqtlgen = read_tsv("2019-12-11-cis-eQTLsFDR0.05-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt")


# RADIAL - RS528236848 ----------------------------------------------------

variant_ix = 9

v2g = get_V2G_data(hits$id[variant_ix])
hits$closest_gene[variant_ix] = v2g$gene[which.min(unlist(lapply(v2g$distances, function(x) x$tissues))),]

# make a granges object for the variant

snp_pos = hits$pos_38[variant_ix]
win = 2e6
region_granges = GenomicRanges::GRanges(
  seqnames = hits$chr[variant_ix], 
  ranges = IRanges::IRanges(start=snp_pos-win, end=snp_pos+win), 
  strand = "*")
region_granges

# gwas dat

gwas = read.table(pipe(paste0("awk 'NR==1 {print}; $2==", hits$chr[variant_ix], " && $3>", hits$pos_37[variant_ix]-win, " && $3<", hits$pos_37[variant_ix]+win, " {print}' /gpfs01/bhcbio/projects/UK_Biobank/20190102_UK_Biobank_Imaging/Results/GWAS/GWAS_diastolic_BOLT/Results/bolt_", hits$gwas_root[variant_ix], "_full.bgen.stats")), header=TRUE)
gwas = makeGRangesFromDataFrame(gwas, keep.extra.columns=TRUE, start.field="BP", end.field="BP")
gwas = lift_over(gwas, dir="37_to_38")
write_tsv(gwas, paste0("data/", hits$variant[variant_ix], "/gwas.txt"))

# check intensity/cluster plots for genotyping

snp_mat = getImputedV3GenotypesForVariants(hits$variant[variant_ix])
attr(snp_mat, "metadata")
snpStats::col.summary(snp_mat)
geno_data = getTibbleFromSnpMatrix(snp_mat)

phewas = bind_rows(
  getPheWASResults(geno_data %>% dplyr::select(SID, !!as.name(hits$variant[variant_ix])), slurm=F, cores=8),
  getQuantPheWASResults(geno_data %>% dplyr::select(SID, !!as.name(hits$variant[variant_ix])))
) %>% tbl_df()

phewas %>% arrange(p.value) %>% filter(p.value < 0.05) %>% View()
pheWASForest(phewas %>% filter(p.value < 0.05) %>% dplyr::select(-n) %>% arrange(p.value))

traits = c("sbp_adj","pulse_rate","dbp_adj","cholesterol")
gwas_dat = vector("list", length(traits))
coloc_res_gwas = gwas_dat
names(gwas_dat) = traits
names(coloc_res_gwas) = names(gwas_dat)

for(i in 1:length(gwas_dat)) {
  
  gwas_comp = loadGWAS(trait=traits[i], type="linear") # not working
  gwas_comp = as.data.frame(gwas_comp) %>% dplyr::filter(CHR==hits$chr[variant_ix], (BP > hits$pos_37[variant_ix]-win & BP < hits$pos_37[variant_ix]+win)) %>% arrange(BP)
  gwas_comp = makeGRangesFromDataFrame(gwas_comp, keep.extra.columns=TRUE, start.field="BP", end.field="BP")
  gwas_comp = lift_over(gwas_comp, dir="37_to_38")
  
  common_variants = intersect(gwas$start, gwas_comp$start)
  coloc_input = data.frame(
    beta1 = as.numeric(gwas$BETA[match(common_variants, gwas$start)]),
    se1 = as.numeric(gwas$SE[match(common_variants, gwas$start)]),
    beta2 = as.numeric(gwas_comp$ESTIMATE[match(common_variants, gwas_comp$start)]),
    se2 = as.numeric(gwas_comp$SE[match(common_variants, gwas_comp$start)])
  )
  coloc_res_gwas[[i]] = coloc_wrapper(coloc_input)
  
  gwas_comp = gwas_comp %>% dplyr::select(P, start) %>% mutate(p_value = -log(P,10) / max(-log(P,10)))
  gwas_dat[[i]] = gwas_comp
}

coloc_plot_gwas = bind_rows(
  gwas %>% dplyr::select(P_BOLT_LMM, start) %>% mutate(p_value = -log(P_BOLT_LMM,10) / max(-log(P_BOLT_LMM,10)), Group="GWAS"),
  bind_rows(gwas_dat, .id="Group")
)
ggplot(coloc_plot_gwas, aes(x=start, y=p_value, color=Group)) + geom_point(alpha=0.5, size=1) + theme_thesis(15) + ylab("-log10(P)") + xlab("")


# RADIAL - RS9388001 ------------------------------------------------------

# get closest gene from open targets genetics
variant_ix = 10

v2g = get_V2G_data(hits$id[variant_ix])
hits$closest_gene[variant_ix] = v2g$gene[which.min(unlist(lapply(v2g$distances, function(x) x$tissues))),]

# make a granges object for the variant

snp_pos = hits$pos_38[variant_ix]
win = 2e6
region_granges = GenomicRanges::GRanges(
  seqnames = hits$chr[variant_ix], 
  ranges = IRanges::IRanges(start=snp_pos-win, end=snp_pos+win), 
  strand = "*")
region_granges

# "Alasoo_2018","BLUEPRINT","BrainSeq","CEDAR","Fairfax_2012","Fairfax_2014","GENCORD","GEUVADIS","HipSci","Kasela_2017","Lepik_2017","Naranbhai_2015","Nedelec_2016","Quach_2016","ROSMAP","Schmiedel_2018","Schwartzentruber_2018","TwinsUK","van_de_Bunt_2015"
# "macrophage","monocyte","neutrophil","CD4+,T,cell","DLPFC","CD8+,T,cell","transverse,colon","platelet","rectum","B,cell","ileum","LCL","fibroblast","T,cell","iPSC","blood","Tfh,cell","Th17,cell","Th1,cell","Th2,cell","Treg,naive","Treg,memory","CD16+,monocyte","NK,cell","sensory,neuron","adipose","skin","pancreatic,islet"

# pull in eqtl
to_pull = data.frame(
  ensembl_id = v2g$gene_id,
  study = c("eQTLGen", "eQTLGen", "GTEx_V8", NA),
  tissue = c("Blood", "Blood", "Nerve - Tibial", NA),
  source = c("eqtlgen", "eqtlgen", "api", "none"),
  gene = mapping$external_gene_name[match(v2g$gene_id, mapping$ensembl_gene_id)]

)

for(i in 1:dim(to_pull)[1]) {
  
  if(to_pull$source[i]=="api") {
    
    if(to_pull$study[i]=="GTEx_V8") {
      eqtl_df = dplyr::filter(imported_tabix_paths, study==to_pull$study[i], tissue_label==to_pull$tissue[i])  
    } else {
      eqtl_df = dplyr::filter(tabix_paths, study==to_pull$study[i], tissue_label==to_pull$tissue[i])  
    }
    column_names = colnames(readr::read_tsv(eqtl_df$ftp_path[1], n_max=1))
    summary_stats = import_eqtl_catalog(eqtl_df$ftp_path[1], region_granges, selected_gene_id=to_pull$ensembl_id[i], column_names)
    print(ggplot(summary_stats, aes(x=position, y=-log(pvalue,10))) + geom_point() + geom_vline(xintercept=hits$pos_38[variant_ix]) + theme_thesis(15) + ggtitle(paste(unlist(to_pull[i,]), collapse=", ")))
    write_tsv(summary_stats, paste0("data/", hits$variant[variant_ix] ,"/eqtl_", to_pull$ensembl_id[i], "_", to_pull$study[i], "_", to_pull$tissue[i], ".txt"))
    
  }
  
  if(to_pull$source[i]=="eqtlgen") {
    
    summary_stats = eqtlgen %>% filter(Gene==to_pull$ensembl_id[i])
    summary_stats = makeGRangesFromDataFrame(summary_stats, keep.extra.columns=TRUE, start.field="SNPPos", end.field="SNPPos", seqnames.field="SNPChr")
    summary_stats = lift_over(summary_stats, dir="37_to_38")
    write_tsv(summary_stats, paste0("data/", hits$variant[variant_ix] ,"/eqtl_", to_pull$ensembl_id[i], "_", to_pull$study[i], "_", to_pull$tissue[i], ".txt"))
    
  }
}

# get the imaging summary stats for the locus (this is the same data as the locus zoom plots)

gwas = read.table(pipe(paste0("awk 'NR==1 {print}; $2==", hits$chr[variant_ix], " && $3>", hits$pos_37[variant_ix]-win, " && $3<", hits$pos_37[variant_ix]+win, " {print}' /gpfs01/bhcbio/projects/UK_Biobank/20190102_UK_Biobank_Imaging/Results/GWAS/GWAS_diastolic_BOLT/Results/bolt_", hits$gwas_root[variant_ix], "_full.bgen.stats")), header=TRUE)
gwas = makeGRangesFromDataFrame(gwas, keep.extra.columns=TRUE, start.field="BP", end.field="BP")
gwas = lift_over(gwas, dir="37_to_38")
write_tsv(gwas, paste0("data/", hits$variant[variant_ix], "/gwas.txt"))

ggplot(gwas, aes(x=start, y=-log(as.numeric(P_BOLT_LMM),10))) + geom_point() + geom_vline(xintercept=hits$pos_38[variant_ix]) + theme_thesis(15) + ylab("-log10(P)") + xlab("")

qplot(-log(as.numeric(gwas$P_BOLT_LMM),10), -log(summary_stats$Pvalue,10)[match(gwas$start, summary_stats$start)]) + theme_thesis(15) + xlab("GWAS") + ylab("eQTL")

qtl_dat = vector("list", dim(to_pull)[1])
coloc_res = qtl_dat
names(qtl_dat) = apply(to_pull, 1, function(x) paste(x, collapse="_"))
names(coloc_res) = names(qtl_dat)

for(i in 1:length(qtl_dat)) {
  
  if(to_pull$source[i]=="none") next
  
  dat = read_tsv(paste0("data/", hits$variant[variant_ix] ,"/eqtl_", to_pull$ensembl_id[i], "_", to_pull$study[i], "_", to_pull$tissue[i], ".txt"))
  if(to_pull$source[i]=="api") {
    dat = dat %>% dplyr::rename(Pvalue="pvalue", start="position")
    common_variants = intersect(gwas$start, dat$start)
    coloc_input = data.frame(
      beta1 = as.numeric(gwas$BETA[match(common_variants, gwas$start)]),
      se1 = as.numeric(gwas$SE[match(common_variants, gwas$start)]),
      beta2 = as.numeric(dat$beta[match(common_variants, dat$start)]),
      se2 = as.numeric(dat$se[match(common_variants, dat$start)])
    )
    coloc_res[[i]] = coloc_wrapper(coloc_input)
  }
  
  # dat = dat %>% dplyr::select(Pvalue, start) %>% mutate(p_value = -log(Pvalue,10) / max(-log(Pvalue,10)))
  dat = dat %>% dplyr::select(Pvalue, start) %>% mutate(p_value = -log(Pvalue,10))
  qtl_dat[[i]] = dat
}

lapply(coloc_res, function(x) x$posterior)

coloc_plot = bind_rows(
  # gwas %>% dplyr::select(P_BOLT_LMM, start) %>% mutate(p_value = -log(P_BOLT_LMM,10) / max(-log(P_BOLT_LMM,10)), Group="GWAS"),
  gwas %>% dplyr::select(P_BOLT_LMM, start) %>% mutate(p_value = -log(P_BOLT_LMM,10), Group="GWAS"),
  bind_rows(qtl_dat, .id="Group")
)

scale_this <- function(x) {
  return(x / max(x))
}

coloc_plot = coloc_plot %>% group_by(Group) %>% mutate(p_value_scaled=scale_this(p_value))
p1 = ggplot(coloc_plot, aes(x=start, y=p_value, color=Group)) + geom_point(alpha=0.5, size=1) + theme_thesis(15) + ylab("-log10(P)") + xlab("") + geom_hline(yintercept = -log(5e-8, base=10), alpha=0.5, lty=2, color="grey") + geom_vline(xintercept=hits$pos_38[variant_ix], alpha=0.5, color="grey", lty=2) + theme(legend.position="None")
p2 = ggplot(coloc_plot, aes(x=start, y=p_value_scaled, color=Group)) + geom_point(alpha=0.5, size=1) + theme_thesis(15) + ylab("-log10(P) Scaled") + xlab("") + geom_vline(xintercept=hits$pos_38[variant_ix], alpha=0.5, color="grey", lty=2) + theme(legend.position="None")

# phewas plot?
# run phewas and pull signal for anything significant and add to colocalisation step?

snp_mat = getImputedV3GenotypesForVariants(hits$variant[variant_ix])
attr(snp_mat, "metadata")
snpStats::col.summary(snp_mat)
geno_data = getTibbleFromSnpMatrix(snp_mat)

# compute associations

phewas = bind_rows(
  getPheWASResults(geno_data %>% dplyr::select(SID, !!as.name(hits$variant[variant_ix])), slurm=F, cores=8),
  getQuantPheWASResults(geno_data %>% dplyr::select(SID, !!as.name(hits$variant[variant_ix])))
) %>% tbl_df()

phewas %>% arrange(p.value) %>% filter(p.value < 0.05) %>% View()
pheWASForest(phewas %>% filter(p.value < 0.05) %>% dplyr::select(-n) %>% arrange(p.value))

# pulse rate!

traits = c("sbp_adj","pulse_rate","dbp_adj")
gwas_dat = vector("list", length(traits))
coloc_res_gwas = gwas_dat
names(gwas_dat) = traits
names(coloc_res_gwas) = names(gwas_dat)

for(i in 1:length(gwas_dat)) {
  
  gwas_comp = loadGWAS(trait=traits[i], type="linear")
  gwas_comp = as.data.frame(gwas_comp) %>% dplyr::filter(CHR==hits$chr[variant_ix], (BP > hits$pos_37[variant_ix]-win & BP < hits$pos_37[variant_ix]+win)) %>% arrange(BP)
  gwas_comp = makeGRangesFromDataFrame(gwas_comp, keep.extra.columns=TRUE, start.field="BP", end.field="BP")
  gwas_comp = lift_over(gwas_comp, dir="37_to_38")
  
  # common_variants = intersect(gwas$start, gwas_comp$start)
  # coloc_input = data.frame(
  #   beta1 = as.numeric(gwas$BETA[match(common_variants, gwas$start)]),
  #   se1 = as.numeric(gwas$SE[match(common_variants, gwas$start)]),
  #   beta2 = as.numeric(gwas_comp$ESTIMATE[match(common_variants, gwas_comp$start)]),
  #   se2 = as.numeric(gwas_comp$SE[match(common_variants, gwas_comp$start)])
  # )
  # coloc_res_gwas[[i]] = coloc_wrapper(coloc_input)
  
  gwas_comp = gwas_comp %>% dplyr::select(P, start) %>% mutate(p_value_scaled = -log(P,10) / max(-log(P,10)), p_value = -log(P,10))
  gwas_dat[[i]] = gwas_comp
}

win = 1e6
load("~/links/bullseye/r_data/t_list.RData")
t_list_filtered = t_list %>% filter(chromosome_name==hits$chr[variant_ix], exon_chrom_start > hits$pos_38[variant_ix]-win, exon_chrom_start < hits$pos_38[variant_ix]+win)
pick_t = t_list_filtered %>% group_by(external_gene_name) %>% summarise(pick_t=ensembl_transcript_id[1]) %>% dplyr::select(pick_t) %>% unlist()
t_list_filtered = t_list_filtered %>% dplyr::filter(ensembl_transcript_id %in% pick_t)
t_list_filtered$strand = "*"
t_list_filtered$ensembl_transcript_id = t_list_filtered$external_gene_name
gene_track = makeGRangesFromDataFrame(t_list_filtered, keep.extra.columns=TRUE)
gene_track$model = "exon"

p3 = ggbio::autoplot(gene_track, aes(tgene_trackpe=model, group=ensembl_transcript_id)) + theme_thesis(10, angle_45=FALSE)
p3
dir.create(paste0("tex/images/", hits$variant[variant_ix]))
pdf(file=paste0("tex/images/", hits$variant[variant_ix], "/coloc_plot.pdf"))
ggbio::tracks(p1, p2, p3, heights=c(1,1,1))
dev.off()

coloc_plot_gwas = bind_rows(
  gwas %>% dplyr::select(P_BOLT_LMM, start) %>% mutate(p_value_scaled = -log(P_BOLT_LMM,10) / max(-log(P_BOLT_LMM,10)), Group="GWAS", p_value = -log(P_BOLT_LMM,10)),
  bind_rows(gwas_dat, .id="Group")
)
p4 = ggplot(coloc_plot_gwas, aes(x=start, y=p_value, color=Group)) + geom_point(alpha=0.5, size=1) + theme_thesis(15) + ylab("-log10(P)") + xlab("") + geom_hline(yintercept = -log(5e-8, base=10), alpha=0.5, lty=2, color="grey") + geom_vline(xintercept=hits$pos_38[variant_ix], alpha=0.5, color="grey", lty=2) + theme(legend.position="None")
p5 = ggplot(coloc_plot_gwas, aes(x=start, y=p_value_scaled, color=Group)) + geom_point(alpha=0.5, size=1) + theme_thesis(15) + ylab("-log10(P) Scaled") + xlab("") + geom_vline(xintercept=hits$pos_38[variant_ix], alpha=0.5, color="grey", lty=2) + theme(legend.position="None")
pdf(file=paste0("tex/images/", hits$variant[variant_ix], "/coloc_plot_gwas.pdf"))
ggbio::tracks(p4, p5, p3, heights=c(1,1,1))
dev.off()


# RADIAL - RS2234962 ------------------------------------------------------

# get closest gene from open targets genetics
variant_ix = 11

v2g = get_V2G_data(hits$id[variant_ix])
hits$closest_gene[variant_ix] = v2g$gene[which.min(unlist(lapply(v2g$distances, function(x) x$tissues))),]

# make a granges object for the variant

snp_pos = hits$pos_38[variant_ix]
win = 2e6
region_granges = GenomicRanges::GRanges(
  seqnames = hits$chr[variant_ix], 
  ranges = IRanges::IRanges(start=snp_pos-win, end=snp_pos+win), 
  strand = "*")
region_granges

# "Alasoo_2018","BLUEPRINT","BrainSeq","CEDAR","Fairfax_2012","Fairfax_2014","GENCORD","GEUVADIS","HipSci","Kasela_2017","Lepik_2017","Naranbhai_2015","Nedelec_2016","Quach_2016","ROSMAP","Schmiedel_2018","Schwartzentruber_2018","TwinsUK","van_de_Bunt_2015"
# "macrophage","monocyte","neutrophil","CD4+,T,cell","DLPFC","CD8+,T,cell","transverse,colon","platelet","rectum","B,cell","ileum","LCL","fibroblast","T,cell","iPSC","blood","Tfh,cell","Th17,cell","Th1,cell","Th2,cell","Treg,naive","Treg,memory","CD16+,monocyte","NK,cell","sensory,neuron","adipose","skin","pancreatic,islet"

# pull in eqtl
to_pull = data.frame(
  ensembl_id = c(rep("ENSG00000197771",6), "ENSG00000151929"),
  study = c("Fairfax_2012","CEDAR","Fairfax_2014","CEDAR","CEDAR","eQTLGen","eQTLGen"),
  tissue = c("B cell","monocyte","monocyte","CD4+ T cell","CD8+ T cell","Blood","Blood"),
  source = c("api", "api", "api", "api", "api", "eqtlgen", "eqtlgen")
)

for(i in 1:dim(to_pull)[1]) {
  
  # if(to_pull$source[i]=="api") {
  #   
  #   if(to_pull$study[i]=="GTEx_V8") {
  #     eqtl_df = dplyr::filter(imported_tabix_paths, study==to_pull$study[i], tissue_label==to_pull$tissue[i])  
  #   } else {
  #     eqtl_df = dplyr::filter(tabix_paths, study==to_pull$study[i], tissue_label==to_pull$tissue[i])  
  #   }
  #   column_names = colnames(readr::read_tsv(eqtl_df$ftp_path[1], n_max=1))
  #   summary_stats = import_eqtl_catalog(eqtl_df$ftp_path[1], region_granges, selected_gene_id=to_pull$ensembl_id[i], column_names)
  #   print(ggplot(summary_stats, aes(x=position, y=-log(pvalue,10))) + geom_point() + geom_vline(xintercept=hits$pos_38[variant_ix]) + theme_thesis(15) + ggtitle(paste(unlist(to_pull[i,]), collapse=", ")))
  #   write_tsv(summary_stats, paste0("data/", hits$variant[variant_ix] ,"/eqtl_", to_pull$ensembl_id[i], "_", to_pull$study[i], "_", to_pull$tissue[i], ".txt"))
  #   
  # }
  
  if(to_pull$source[i]=="eqtlgen") {
    
    summary_stats = eqtlgen %>% filter(Gene==to_pull$ensembl_id[i])
    summary_stats = makeGRangesFromDataFrame(summary_stats, keep.extra.columns=TRUE, start.field="SNPPos", end.field="SNPPos", seqnames.field="SNPChr")
    summary_stats = lift_over(summary_stats, dir="37_to_38")
    write_tsv(summary_stats, paste0("data/", hits$variant[variant_ix] ,"/eqtl_", to_pull$ensembl_id[i], "_", to_pull$study[i], "_", to_pull$tissue[i], ".txt"))
    
  }
}

# get the imaging summary stats for the locus (this is the same data as the locus zoom plots)

gwas = read.table(pipe(paste0("awk 'NR==1 {print}; $2==", hits$chr[variant_ix], " && $3>", hits$pos_37[variant_ix]-win, " && $3<", hits$pos_37[variant_ix]+win, " {print}' /gpfs01/bhcbio/projects/UK_Biobank/20190102_UK_Biobank_Imaging/Results/GWAS/GWAS_diastolic_BOLT/Results/bolt_", hits$gwas_root[variant_ix], "_full.bgen.stats")), header=TRUE)
gwas = makeGRangesFromDataFrame(gwas, keep.extra.columns=TRUE, start.field="BP", end.field="BP")
gwas = lift_over(gwas, dir="37_to_38")
write_tsv(gwas, paste0("data/", hits$variant[variant_ix], "/gwas.txt"))

ggplot(gwas, aes(x=start, y=-log(as.numeric(P_BOLT_LMM),10))) + geom_point() + geom_vline(xintercept=hits$pos_38[variant_ix]) + theme_thesis(15) + ylab("-log10(P)") + xlab("")

qtl_dat = vector("list", dim(to_pull)[1])
coloc_res = qtl_dat
names(qtl_dat) = apply(to_pull, 1, function(x) paste(x, collapse="_"))
names(coloc_res) = names(qtl_dat)

for(i in 1:length(qtl_dat)) {
  
  if(to_pull$source[i]=="none") next
  
  dat = read_tsv(paste0("data/", hits$variant[variant_ix] ,"/eqtl_", to_pull$ensembl_id[i], "_", to_pull$study[i], "_", to_pull$tissue[i], ".txt"))
  if(to_pull$source[i]=="api") {
    dat = dat %>% rename(pvalue="Pvalue", position="start")
    # common_variants = intersect(gwas$start, dat$start)
    # coloc_input = data.frame(
    #   beta1 = as.numeric(gwas$BETA[match(common_variants, gwas$start)]),
    #   se1 = as.numeric(gwas$SE[match(common_variants, gwas$start)]),
    #   beta2 = as.numeric(dat$beta[match(common_variants, dat$start)]),
    #   se2 = as.numeric(dat$se[match(common_variants, dat$start)])
    # )
    # coloc_res[[i]] = coloc_wrapper(coloc_input)
  }
  
  dat = dat %>% dplyr::select(Pvalue, start) %>% mutate(p_value = -log(Pvalue,10) / max(-log(Pvalue,10)))
  qtl_dat[[i]] = dat
}

lapply(coloc_res, function(x) x$posterior)

coloc_plot = bind_rows(
  gwas %>% dplyr::select(P_BOLT_LMM, start) %>% mutate(p_value = -log(P_BOLT_LMM,10), Group="GWAS"),
  bind_rows(qtl_dat, .id="Group")
)

scale_this <- function(x) {
  return(x / max(x))
}

coloc_plot = coloc_plot %>% group_by(Group) %>% mutate(p_value_scaled=scale_this(p_value))
p1 = ggplot(coloc_plot, aes(x=start, y=p_value, color=Group)) + geom_point(alpha=0.5, size=1) + theme_thesis(15) + ylab("-log10(P)") + xlab("") + geom_hline(yintercept = -log(5e-8, base=10), alpha=0.5, lty=2, color="grey") + geom_vline(xintercept=hits$pos_38[variant_ix], alpha=0.5, color="grey", lty=2) + theme(legend.position="None")
p2 = ggplot(coloc_plot, aes(x=start, y=p_value_scaled, color=Group)) + geom_point(alpha=0.5, size=1) + theme_thesis(15) + ylab("-log10(P) Scaled") + xlab("") + geom_vline(xintercept=hits$pos_38[variant_ix], alpha=0.5, color="grey", lty=2) + theme(legend.position="None")

# phewas plot?
# run phewas and pull signal for anything significant and add to colocalisation step?

snp_mat = getImputedV3GenotypesForVariants(hits$variant[variant_ix])
attr(snp_mat, "metadata")
snpStats::col.summary(snp_mat)
geno_data = getTibbleFromSnpMatrix(snp_mat)

# compute associations

phewas = bind_rows(
  getPheWASResults(geno_data %>% dplyr::select(SID, !!as.name(hits$variant[variant_ix])), slurm=F, cores=8),
  getQuantPheWASResults(geno_data %>% dplyr::select(SID, !!as.name(hits$variant[variant_ix])))
) %>% tbl_df()

phewas %>% arrange(p.value) %>% filter(p.value < 0.05) %>% View()
pheWASForest(phewas %>% filter(p.value < 0.05) %>% dplyr::select(-n) %>% arrange(p.value))

# pulse rate!

traits = c("sbp_adj","pulse_rate","dbp_adj")
gwas_dat = vector("list", length(traits))
coloc_res_gwas = gwas_dat
names(gwas_dat) = traits
names(coloc_res_gwas) = names(gwas_dat)

for(i in 1:length(gwas_dat)) {
  
  gwas_comp = loadGWAS(trait=traits[i], type="linear")
  gwas_comp = as.data.frame(gwas_comp) %>% dplyr::filter(CHR==hits$chr[variant_ix], (BP > hits$pos_37[variant_ix]-win & BP < hits$pos_37[variant_ix]+win)) %>% arrange(BP)
  gwas_comp = makeGRangesFromDataFrame(gwas_comp, keep.extra.columns=TRUE, start.field="BP", end.field="BP")
  gwas_comp = lift_over(gwas_comp, dir="37_to_38")
  
  common_variants = intersect(gwas$start, gwas_comp$start)
  coloc_input = data.frame(
    beta1 = as.numeric(gwas$BETA[match(common_variants, gwas$start)]),
    se1 = as.numeric(gwas$SE[match(common_variants, gwas$start)]),
    beta2 = as.numeric(gwas_comp$ESTIMATE[match(common_variants, gwas_comp$start)]),
    se2 = as.numeric(gwas_comp$SE[match(common_variants, gwas_comp$start)])
  )
  coloc_res_gwas[[i]] = coloc_wrapper(coloc_input)
  
  gwas_comp = gwas_comp %>% dplyr::select(P, start) %>% mutate(p_value = -log(P,10) / max(-log(P,10)))
  gwas_dat[[i]] = gwas_comp
}

win = 1e6
load("~/links/bullseye/r_data/t_list.RData")
t_list_filtered = t_list %>% filter(chromosome_name==hits$chr[variant_ix], exon_chrom_start > hits$pos_38[variant_ix]-win, exon_chrom_start < hits$pos_38[variant_ix]+win)
pick_t = t_list_filtered %>% group_by(external_gene_name) %>% summarise(pick_t=ensembl_transcript_id[1]) %>% dplyr::select(pick_t) %>% unlist()
t_list_filtered = t_list_filtered %>% dplyr::filter(ensembl_transcript_id %in% pick_t)
t_list_filtered$strand = "*"
t_list_filtered$ensembl_transcript_id = t_list_filtered$external_gene_name
gene_track = makeGRangesFromDataFrame(t_list_filtered, keep.extra.columns=TRUE)
gene_track$model = "exon"

p3 = ggbio::autoplot(gene_track, aes(tgene_trackpe=model, group=ensembl_transcript_id)) + theme_thesis(10, angle_45=FALSE) + geom_vline(xintercept=hits$pos_38[variant_ix], alpha=0.5, color="grey", lty=2) 
p3
dir.create(paste0("tex/images/", hits$variant[variant_ix]))
pdf(file=paste0("tex/images/", hits$variant[variant_ix], "/coloc_plot.pdf"))
ggbio::tracks(p1, p2, p3, heights=c(1,1,1))
dev.off()

coloc_plot_gwas = bind_rows(
  gwas %>% dplyr::select(P_BOLT_LMM, start) %>% mutate(p_value_scaled = -log(P_BOLT_LMM,10) / max(-log(P_BOLT_LMM,10)), Group="GWAS", p_value = -log(P_BOLT_LMM,10)),
  bind_rows(gwas_dat, .id="Group")
)
p4 = ggplot(coloc_plot_gwas, aes(x=start, y=p_value, color=Group)) + geom_point(alpha=0.5, size=1) + theme_thesis(15) + ylab("-log10(P)") + xlab("") + geom_hline(yintercept = -log(5e-8, base=10), alpha=0.5, lty=2, color="grey") + geom_vline(xintercept=hits$pos_38[variant_ix], alpha=0.5, color="grey", lty=2) + theme(legend.position="None")
p5 = ggplot(coloc_plot_gwas, aes(x=start, y=p_value_scaled, color=Group)) + geom_point(alpha=0.5, size=1) + theme_thesis(15) + ylab("-log10(P) Scaled") + xlab("") + geom_vline(xintercept=hits$pos_38[variant_ix], alpha=0.5, color="grey", lty=2) + theme(legend.position="None")
pdf(file=paste0("tex/images/", hits$variant[variant_ix], "/coloc_plot_gwas.pdf"))
ggbio::tracks(p4, p5, p3, heights=c(1,1,1))
dev.off()


# RADIAL - RS11170519 -----------------------------------------------------

# get closest gene from open targets genetics
variant_ix = 12

v2g = get_V2G_data(hits$id[variant_ix])
v2g = v2g[!unlist(lapply(v2g$qtls, is_empty)),] # only take anything with qtls
hits$closest_gene[variant_ix] = v2g$gene[which.min(unlist(lapply(v2g$distances, function(x) x$tissues))),]

# make a granges object for the variant

snp_pos = hits$pos_38[variant_ix]
win = 2e6
region_granges = GenomicRanges::GRanges(
  seqnames = hits$chr[variant_ix], 
  ranges = IRanges::IRanges(start=snp_pos-win, end=snp_pos+win), 
  strand = "*")
region_granges

# "Alasoo_2018","BLUEPRINT","BrainSeq","CEDAR","Fairfax_2012","Fairfax_2014","GENCORD","GEUVADIS","HipSci","Kasela_2017","Lepik_2017","Naranbhai_2015","Nedelec_2016","Quach_2016","ROSMAP","Schmiedel_2018","Schwartzentruber_2018","TwinsUK","van_de_Bunt_2015"
# "macrophage","monocyte","neutrophil","CD4+,T,cell","DLPFC","CD8+,T,cell","transverse,colon","platelet","rectum","B,cell","ileum","LCL","fibroblast","T,cell","iPSC","blood","Tfh,cell","Th17,cell","Th1,cell","Th2,cell","Treg,naive","Treg,memory","CD16+,monocyte","NK,cell","sensory,neuron","adipose","skin","pancreatic,islet"

# pull in eqtl
to_pull = data.frame(
  ensembl_id = v2g$gene$id,
  gene = mapping$Symbol[match(v2g$gene$id, mapping$ENSEMBL_ID)],
  study = c("GTEx_V8", "GTEx_V8", "eQTLGen", "GTEx_V8", "BLUEPRINT", "eQTLGen", "eQTLGen"),
  tissue = c("Cells - Cultured fibroblasts", "Esophagus - Mucosa", "Blood", "Testis", "monocyte", "Blood", "Blood"),
  source = c("api", "api", "eqtlgen", "api", "api", "eqtlgen", "eqtlgen")
)

for(i in 1:dim(to_pull)[1]) {
  
  if(to_pull$source[i]=="api") {
    
    if(to_pull$study[i]=="GTEx_V8") {
      eqtl_df = dplyr::filter(imported_tabix_paths, study==to_pull$study[i], tissue_label==to_pull$tissue[i])  
    } else {
      eqtl_df = dplyr::filter(tabix_paths, study==to_pull$study[i], tissue_label==to_pull$tissue[i])  
    }
    column_names = colnames(readr::read_tsv(eqtl_df$ftp_path[1], n_max=1))
    summary_stats = import_eqtl_catalog(eqtl_df$ftp_path[1], region_granges, selected_gene_id=to_pull$ensembl_id[i], column_names)
    print(ggplot(summary_stats, aes(x=position, y=-log(pvalue,10))) + geom_point() + geom_vline(xintercept=hits$pos_38[variant_ix]) + theme_thesis(15) + ggtitle(paste(unlist(to_pull[i,]), collapse=", ")))
    write_tsv(summary_stats, paste0("data/", hits$variant[variant_ix] ,"/eqtl_", to_pull$ensembl_id[i], "_", to_pull$study[i], "_", to_pull$tissue[i], ".txt"))
    
  }
  
  if(to_pull$source[i]=="eqtlgen") {
    
    summary_stats = eqtlgen %>% filter(Gene==to_pull$ensembl_id[i])
    summary_stats = makeGRangesFromDataFrame(summary_stats, keep.extra.columns=TRUE, start.field="SNPPos", end.field="SNPPos", seqnames.field="SNPChr")
    summary_stats = lift_over(summary_stats, dir="37_to_38")
    write_tsv(summary_stats, paste0("data/", hits$variant[variant_ix] ,"/eqtl_", to_pull$ensembl_id[i], "_", to_pull$study[i], "_", to_pull$tissue[i], ".txt"))
    
  }
}

# get the imaging summary stats for the locus (this is the same data as the locus zoom plots)

gwas = read.table(pipe(paste0("awk 'NR==1 {print}; $2==", hits$chr[variant_ix], " && $3>", hits$pos_37[variant_ix]-win, " && $3<", hits$pos_37[variant_ix]+win, " {print}' /gpfs01/bhcbio/projects/UK_Biobank/20190102_UK_Biobank_Imaging/Results/GWAS/GWAS_diastolic_BOLT/Results/bolt_", hits$gwas_root[variant_ix], "_full.bgen.stats")), header=TRUE)
gwas = makeGRangesFromDataFrame(gwas, keep.extra.columns=TRUE, start.field="BP", end.field="BP")
gwas = lift_over(gwas, dir="37_to_38")
write_tsv(gwas, paste0("data/", hits$variant[variant_ix], "/gwas.txt"))

ggplot(gwas, aes(x=start, y=-log(as.numeric(P_BOLT_LMM),10))) + geom_point() + geom_vline(xintercept=hits$pos_38[variant_ix]) + theme_thesis(15) + ylab("-log10(P)") + xlab("")

qtl_dat = vector("list", dim(to_pull)[1])
coloc_res = qtl_dat
names(qtl_dat) = apply(to_pull, 1, function(x) paste(x, collapse="_"))
names(coloc_res) = names(qtl_dat)

for(i in 1:length(qtl_dat)) {
  
  if(to_pull$source[i]=="none") next
  
  dat = read_tsv(paste0("data/", hits$variant[variant_ix] ,"/eqtl_", to_pull$ensembl_id[i], "_", to_pull$study[i], "_", to_pull$tissue[i], ".txt"))
  if(to_pull$source[i]=="api") {
    dat = dat %>% rename(pvalue="Pvalue", position="start")
    common_variants = intersect(gwas$start, dat$start)
    if(!to_pull$study[i] %in% c("BLUEPRINT")) { # blueprint does not have se
      coloc_input = data.frame(
        beta1 = as.numeric(gwas$BETA[match(common_variants, gwas$start)]),
        se1 = as.numeric(gwas$SE[match(common_variants, gwas$start)]),
        beta2 = as.numeric(dat$beta[match(common_variants, dat$start)]),
        se2 = as.numeric(dat$se[match(common_variants, dat$start)])
      )
      coloc_res[[i]] = coloc_wrapper(coloc_input)
    }
  }
  
  dat = dat %>% dplyr::select(Pvalue, start) %>% mutate(p_value = -log(Pvalue,10) / max(-log(Pvalue,10)))
  qtl_dat[[i]] = dat
  
}

lapply(coloc_res, function(x) x$posterior)

coloc_plot = bind_rows(
  gwas %>% dplyr::select(P_BOLT_LMM, start) %>% mutate(p_value = -log(P_BOLT_LMM,10) / max(-log(P_BOLT_LMM,10)), Group="GWAS"),
  bind_rows(qtl_dat, .id="Group")
)
ggplot(coloc_plot, aes(x=start, y=p_value, color=Group)) + geom_point(alpha=0.5, size=1) + theme_thesis(15) + ylab("-log10(P)") + xlab("")

# phewas plot?
# run phewas and pull signal for anything significant and add to colocalisation step?

snp_mat = getImputedV3GenotypesForVariants(hits$variant[variant_ix])
attr(snp_mat, "metadata")
snpStats::col.summary(snp_mat)
geno_data = getTibbleFromSnpMatrix(snp_mat)

# compute associations

phewas = bind_rows(
  getPheWASResults(geno_data %>% dplyr::select(SID, !!as.name(hits$variant[variant_ix])), slurm=F, cores=8),
  getQuantPheWASResults(geno_data %>% dplyr::select(SID, !!as.name(hits$variant[variant_ix])))
) %>% tbl_df()

phewas %>% arrange(p.value) %>% filter(p.value < 0.05) %>% View()
pheWASForest(phewas %>% filter(p.value < 0.05) %>% dplyr::select(-n) %>% arrange(p.value))

# pulse rate!
# add coloc for phenotype data as well

# Sys.setenv(R_CONFIG_ACTIVE="standard")
# gwas_2 = loadGWAS(trait="pulse_rate_adj", type="logistic") # not working


# RADIAL - RS369533272 ----------------------------------------------------

# get closest gene from open targets genetics
variant_ix = 13

v2g = get_V2G_data(hits$id[variant_ix])
hits$closest_gene[variant_ix] = v2g$gene[which.min(unlist(lapply(v2g$distances, function(x) x$tissues))),]

# make a granges object for the variant

snp_pos = hits$pos_38[variant_ix]
win = 2e6
region_granges = GenomicRanges::GRanges(
  seqnames = hits$chr[variant_ix], 
  ranges = IRanges::IRanges(start=snp_pos-win, end=snp_pos+win), 
  strand = "*")
region_granges

# "Alasoo_2018","BLUEPRINT","BrainSeq","CEDAR","Fairfax_2012","Fairfax_2014","GENCORD","GEUVADIS","HipSci","Kasela_2017","Lepik_2017","Naranbhai_2015","Nedelec_2016","Quach_2016","ROSMAP","Schmiedel_2018","Schwartzentruber_2018","TwinsUK","van_de_Bunt_2015"
# "macrophage","monocyte","neutrophil","CD4+,T,cell","DLPFC","CD8+,T,cell","transverse,colon","platelet","rectum","B,cell","ileum","LCL","fibroblast","T,cell","iPSC","blood","Tfh,cell","Th17,cell","Th1,cell","Th2,cell","Treg,naive","Treg,memory","CD16+,monocyte","NK,cell","sensory,neuron","adipose","skin","pancreatic,islet"

# pull in eqtl
to_pull = data.frame(
  ensembl_id = "ENSG00000134779",
  study = c("CEDAR","Fairfax_2012","GEUVADIS","Fairfax_2014","CEDAR","CEDAR"),
  tissue = c("B cell","B cell","LCL","monocyte","CD8+ T cell","CD8+ T cell"),
  source = c("api", "api", "api", "api", "api", "api")
)

for(i in 1:dim(to_pull)[1]) {
  
  if(to_pull$source[i]=="api") {

    if(to_pull$study[i]=="GTEx_V8") {
      eqtl_df = dplyr::filter(imported_tabix_paths, study==to_pull$study[i], tissue_label==to_pull$tissue[i])
    } else {
      eqtl_df = dplyr::filter(tabix_paths, study==to_pull$study[i], tissue_label==to_pull$tissue[i])
    }
    column_names = colnames(readr::read_tsv(eqtl_df$ftp_path[1], n_max=1))
    summary_stats = import_eqtl_catalog(eqtl_df$ftp_path[1], region_granges, selected_gene_id=to_pull$ensembl_id[i], column_names)
    print(ggplot(summary_stats, aes(x=position, y=-log(pvalue,10))) + geom_point() + geom_vline(xintercept=hits$pos_38[variant_ix]) + theme_thesis(15) + ggtitle(paste(unlist(to_pull[i,]), collapse=", ")))
    write_tsv(summary_stats, paste0("data/", hits$variant[variant_ix] ,"/eqtl_", to_pull$ensembl_id[i], "_", to_pull$study[i], "_", to_pull$tissue[i], ".txt"))

  }
  
  if(to_pull$source[i]=="eqtlgen") {
    
    summary_stats = eqtlgen %>% filter(Gene==to_pull$ensembl_id[i])
    summary_stats = makeGRangesFromDataFrame(summary_stats, keep.extra.columns=TRUE, start.field="SNPPos", end.field="SNPPos", seqnames.field="SNPChr")
    summary_stats = lift_over(summary_stats, dir="37_to_38")
    write_tsv(summary_stats, paste0("data/", hits$variant[variant_ix] ,"/eqtl_", to_pull$ensembl_id[i], "_", to_pull$study[i], "_", to_pull$tissue[i], ".txt"))
    
  }
}

# get the imaging summary stats for the locus (this is the same data as the locus zoom plots)

gwas = read.table(pipe(paste0("awk 'NR==1 {print}; $2==", hits$chr[variant_ix], " && $3>", hits$pos_37[variant_ix]-win, " && $3<", hits$pos_37[variant_ix]+win, " {print}' /gpfs01/bhcbio/projects/UK_Biobank/20190102_UK_Biobank_Imaging/Results/GWAS/GWAS_diastolic_BOLT/Results/bolt_", hits$gwas_root[variant_ix], "_full.bgen.stats")), header=TRUE)
gwas = makeGRangesFromDataFrame(gwas, keep.extra.columns=TRUE, start.field="BP", end.field="BP")
gwas = lift_over(gwas, dir="37_to_38")
write_tsv(gwas, paste0("data/", hits$variant[variant_ix], "/gwas.txt"))

ggplot(gwas, aes(x=start, y=-log(as.numeric(P_BOLT_LMM),10))) + geom_point() + geom_vline(xintercept=hits$pos_38[variant_ix]) + theme_thesis(15) + ylab("-log10(P)") + xlab("")

qtl_dat = vector("list", dim(to_pull)[1])
coloc_res = qtl_dat
names(qtl_dat) = apply(to_pull, 1, function(x) paste(x, collapse="_"))
names(coloc_res) = names(qtl_dat)

for(i in 1:length(qtl_dat)) {
  
  if(to_pull$source[i]=="none") next
  
  dat = read_tsv(paste0("data/", hits$variant[variant_ix] ,"/eqtl_", to_pull$ensembl_id[i], "_", to_pull$study[i], "_", to_pull$tissue[i], ".txt"))
  if(to_pull$source[i]=="api") {
    dat = dat %>% rename(pvalue="Pvalue", position="start")
    # common_variants = intersect(gwas$start, dat$start)
    # coloc_input = data.frame(
    #   beta1 = as.numeric(gwas$BETA[match(common_variants, gwas$start)]),
    #   se1 = as.numeric(gwas$SE[match(common_variants, gwas$start)]),
    #   beta2 = as.numeric(dat$beta[match(common_variants, dat$start)]),
    #   se2 = as.numeric(dat$se[match(common_variants, dat$start)])
    # )
    # coloc_res[[i]] = coloc_wrapper(coloc_input)
  }
  
  dat = dat %>% dplyr::select(Pvalue, start) %>% mutate(p_value = -log(Pvalue,10) / max(-log(Pvalue,10)))
  qtl_dat[[i]] = dat
}

lapply(coloc_res, function(x) x$posterior)

coloc_plot = bind_rows(
  gwas %>% dplyr::select(P_BOLT_LMM, start) %>% mutate(p_value = -log(P_BOLT_LMM,10) / max(-log(P_BOLT_LMM,10)), Group="GWAS"),
  bind_rows(qtl_dat, .id="Group")
)
ggplot(coloc_plot, aes(x=start, y=p_value, color=Group)) + geom_point(alpha=0.5, size=1) + theme_thesis(15) + ylab("-log10(P)") + xlab("")

# phewas plot?
# run phewas and pull signal for anything significant and add to colocalisation step?

snp_mat = getImputedV3GenotypesForVariants(hits$variant[variant_ix])
attr(snp_mat, "metadata")
snpStats::col.summary(snp_mat)
geno_data = getTibbleFromSnpMatrix(snp_mat)

# compute associations

phewas = bind_rows(
  getPheWASResults(geno_data %>% dplyr::select(SID, !!as.name(hits$variant[variant_ix])), slurm=F, cores=8),
  getQuantPheWASResults(geno_data %>% dplyr::select(SID, !!as.name(hits$variant[variant_ix])))
) %>% tbl_df()

phewas %>% arrange(p.value) %>% filter(p.value < 0.05) %>% View()
pheWASForest(phewas %>% filter(p.value < 0.05) %>% dplyr::select(-n) %>% arrange(p.value))

# pulse rate!

traits = c("sbp_adj","pulse_rate","dbp_adj")
gwas_dat = vector("list", length(traits))
coloc_res_gwas = gwas_dat
names(gwas_dat) = traits
names(coloc_res_gwas) = names(gwas_dat)

for(i in 1:length(gwas_dat)) {
  
  gwas_comp = loadGWAS(trait=traits[i], type="linear")
  gwas_comp = as.data.frame(gwas_comp) %>% dplyr::filter(CHR==hits$chr[variant_ix], (BP > hits$pos_37[variant_ix]-win & BP < hits$pos_37[variant_ix]+win)) %>% arrange(BP)
  gwas_comp = makeGRangesFromDataFrame(gwas_comp, keep.extra.columns=TRUE, start.field="BP", end.field="BP")
  gwas_comp = lift_over(gwas_comp, dir="37_to_38")
  
  common_variants = intersect(gwas$start, gwas_comp$start)
  coloc_input = data.frame(
    beta1 = as.numeric(gwas$BETA[match(common_variants, gwas$start)]),
    se1 = as.numeric(gwas$SE[match(common_variants, gwas$start)]),
    beta2 = as.numeric(gwas_comp$ESTIMATE[match(common_variants, gwas_comp$start)]),
    se2 = as.numeric(gwas_comp$SE[match(common_variants, gwas_comp$start)])
  )
  coloc_res_gwas[[i]] = coloc_wrapper(coloc_input)
  
  gwas_comp = gwas_comp %>% dplyr::select(P, start) %>% mutate(p_value = -log(P,10) / max(-log(P,10)))
  gwas_dat[[i]] = gwas_comp
}

coloc_plot_gwas = bind_rows(
  gwas %>% dplyr::select(P_BOLT_LMM, start) %>% mutate(p_value = -log(P_BOLT_LMM,10) / max(-log(P_BOLT_LMM,10)), Group="GWAS"),
  bind_rows(gwas_dat, .id="Group")
)
ggplot(coloc_plot_gwas, aes(x=start, y=p_value, color=Group)) + geom_point(alpha=0.5, size=1) + theme_thesis(15) + ylab("-log10(P)") + xlab("")


# LONGITUDINAL - RS2275950 ------------------------------------------------

# get closest gene from open targets genetics
variant_ix = 4

v2g = get_V2G_data(hits$id[variant_ix])
v2g = v2g[!unlist(lapply(v2g$qtls, is_empty)),] # only take anything with qtls
hits$closest_gene[variant_ix] = v2g$gene[which.min(unlist(lapply(v2g$distances, function(x) x$tissues))),]

# make a granges object for the variant

snp_pos = hits$pos_38[variant_ix]
win = 2e6
region_granges = GenomicRanges::GRanges(
  seqnames = hits$chr[variant_ix], 
  ranges = IRanges::IRanges(start=snp_pos-win, end=snp_pos+win), 
  strand = "*")
region_granges

# "Alasoo_2018","BLUEPRINT","BrainSeq","CEDAR","Fairfax_2012","Fairfax_2014","GENCORD","GEUVADIS","HipSci","Kasela_2017","Lepik_2017","Naranbhai_2015","Nedelec_2016","Quach_2016","ROSMAP","Schmiedel_2018","Schwartzentruber_2018","TwinsUK","van_de_Bunt_2015"
# "macrophage","monocyte","neutrophil","CD4+,T,cell","DLPFC","CD8+,T,cell","transverse,colon","platelet","rectum","B,cell","ileum","LCL","fibroblast","T,cell","iPSC","blood","Tfh,cell","Th17,cell","Th1,cell","Th2,cell","Treg,naive","Treg,memory","CD16+,monocyte","NK,cell","sensory,neuron","adipose","skin","pancreatic,islet"

# pull in eqtl
to_pull = data.frame(
  ensembl_id = v2g$gene$id,
  gene = mapping$Symbol[match(v2g$gene$id, mapping$ENSEMBL_ID)],
  study = c("eQTLGen","eQTLGen","eQTLGen","eQTLGen","Fairfax_2014","eQTLGen","GTEx_V8","eQTLGen","BLUEPRINT","eQTLGen","GENCORD"),
  tissue = c("Blood","Blood","Blood","Blood","monocyte","Blood","Skin - Sun Exposed (Lower leg)","Blood","monocyte","Blood","fibroblast"),
  source = c("eqtlgen", "eqtlgen", "eqtlgen", "eqtlgen", "api", "eqtlgen","api","eqtlgen","api","eqtlgen","api")
)

for(i in 1:dim(to_pull)[1]) {
  
  if(to_pull$source[i]=="api") {
    
    if(to_pull$study[i]=="GTEx_V8") {
      eqtl_df = dplyr::filter(imported_tabix_paths, study==to_pull$study[i], tissue_label==to_pull$tissue[i])
    } else {
      eqtl_df = dplyr::filter(tabix_paths, study==to_pull$study[i], tissue_label==to_pull$tissue[i])
    }
    column_names = colnames(readr::read_tsv(eqtl_df$ftp_path[1], n_max=1))
    summary_stats = import_eqtl_catalog(eqtl_df$ftp_path[1], region_granges, selected_gene_id=to_pull$ensembl_id[i], column_names)
    print(ggplot(summary_stats, aes(x=position, y=-log(pvalue,10))) + geom_point() + geom_vline(xintercept=hits$pos_38[variant_ix]) + theme_thesis(15) + ggtitle(paste(unlist(to_pull[i,]), collapse=", ")))
    write_tsv(summary_stats, paste0("data/", hits$variant[variant_ix] ,"/eqtl_", to_pull$ensembl_id[i], "_", to_pull$study[i], "_", to_pull$tissue[i], ".txt"))
    
  }
  
  if(to_pull$source[i]=="eqtlgen") {
    
    summary_stats = eqtlgen %>% filter(Gene==to_pull$ensembl_id[i])
    summary_stats = makeGRangesFromDataFrame(summary_stats, keep.extra.columns=TRUE, start.field="SNPPos", end.field="SNPPos", seqnames.field="SNPChr")
    summary_stats = lift_over(summary_stats, dir="37_to_38")
    write_tsv(summary_stats, paste0("data/", hits$variant[variant_ix] ,"/eqtl_", to_pull$ensembl_id[i], "_", to_pull$study[i], "_", to_pull$tissue[i], ".txt"))
    
  }
}

# get the imaging summary stats for the locus (this is the same data as the locus zoom plots)

gwas = read.table(pipe(paste0("awk 'NR==1 {print}; $2==", hits$chr[variant_ix], " && $3>", hits$pos_37[variant_ix]-win, " && $3<", hits$pos_37[variant_ix]+win, " {print}' /gpfs01/bhcbio/projects/UK_Biobank/20190102_UK_Biobank_Imaging/Results/GWAS/GWAS_diastolic_BOLT/Results/bolt_", hits$gwas_root[variant_ix], "_full.bgen.stats")), header=TRUE)
gwas = makeGRangesFromDataFrame(gwas, keep.extra.columns=TRUE, start.field="BP", end.field="BP")
gwas = lift_over(gwas, dir="37_to_38")
write_tsv(gwas, paste0("data/", hits$variant[variant_ix], "/gwas.txt"))

ggplot(gwas, aes(x=start, y=-log(as.numeric(P_BOLT_LMM),10))) + geom_point() + geom_vline(xintercept=hits$pos_38[variant_ix]) + theme_thesis(15) + ylab("-log10(P)") + xlab("")

qtl_dat = vector("list", dim(to_pull)[1])
coloc_res = qtl_dat
names(qtl_dat) = apply(to_pull, 1, function(x) paste(x, collapse="_"))
names(coloc_res) = names(qtl_dat)

for(i in 1:length(qtl_dat)) {
  
  if(to_pull$source[i]=="none") next
  
  dat = read_tsv(paste0("data/", hits$variant[variant_ix] ,"/eqtl_", to_pull$ensembl_id[i], "_", to_pull$study[i], "_", to_pull$tissue[i], ".txt"))
  if(nrow(dat)==0) next
  
  if(to_pull$source[i]=="api") {
    dat = dat %>% rename(pvalue="Pvalue", position="start")
    # common_variants = intersect(gwas$start, dat$start)
    # coloc_input = data.frame(
    #   beta1 = as.numeric(gwas$BETA[match(common_variants, gwas$start)]),
    #   se1 = as.numeric(gwas$SE[match(common_variants, gwas$start)]),
    #   beta2 = as.numeric(dat$beta[match(common_variants, dat$start)]),
    #   se2 = as.numeric(dat$se[match(common_variants, dat$start)])
    # )
    # coloc_res[[i]] = coloc_wrapper(coloc_input)
  }
  
  dat = dat %>% dplyr::select(Pvalue, start) %>% mutate(p_value = -log(Pvalue,10) / max(-log(Pvalue,10)))
  qtl_dat[[i]] = dat
}

lapply(coloc_res, function(x) x$posterior)

coloc_plot = bind_rows(
  gwas %>% dplyr::select(P_BOLT_LMM, start) %>% mutate(p_value = -log(P_BOLT_LMM,10) / max(-log(P_BOLT_LMM,10)), Group="GWAS"),
  bind_rows(qtl_dat, .id="Group")
)
ggplot(coloc_plot, aes(x=start, y=p_value, color=Group)) + geom_point(alpha=0.5, size=1) + theme_thesis(15) + ylab("-log10(P)") + xlab("")

# phewas plot?
# run phewas and pull signal for anything significant and add to colocalisation step?

snp_mat = getImputedV3GenotypesForVariants(hits$variant[variant_ix])
attr(snp_mat, "metadata")
snpStats::col.summary(snp_mat)
geno_data = getTibbleFromSnpMatrix(snp_mat)

# compute associations

phewas = bind_rows(
  getPheWASResults(geno_data %>% dplyr::select(SID, !!as.name(hits$variant[variant_ix])), slurm=F, cores=8),
  getQuantPheWASResults(geno_data %>% dplyr::select(SID, !!as.name(hits$variant[variant_ix])))
) %>% tbl_df()

phewas %>% arrange(p.value) %>% filter(p.value < 0.05) %>% View()
pheWASForest(phewas %>% filter(p.value < 0.05) %>% dplyr::select(-n) %>% arrange(p.value))

# pulse rate!

traits = c("sbp_adj","pulse_rate","dbp_adj")
gwas_dat = vector("list", length(traits))
coloc_res_gwas = gwas_dat
names(gwas_dat) = traits
names(coloc_res_gwas) = names(gwas_dat)

for(i in 1:length(gwas_dat)) {
  
  gwas_comp = loadGWAS(trait=traits[i], type="linear")
  gwas_comp = as.data.frame(gwas_comp) %>% dplyr::filter(CHR==hits$chr[variant_ix], (BP > hits$pos_37[variant_ix]-win & BP < hits$pos_37[variant_ix]+win)) %>% arrange(BP)
  gwas_comp = makeGRangesFromDataFrame(gwas_comp, keep.extra.columns=TRUE, start.field="BP", end.field="BP")
  gwas_comp = lift_over(gwas_comp, dir="37_to_38")
  
  common_variants = intersect(gwas$start, gwas_comp$start)
  coloc_input = data.frame(
    beta1 = as.numeric(gwas$BETA[match(common_variants, gwas$start)]),
    se1 = as.numeric(gwas$SE[match(common_variants, gwas$start)]),
    beta2 = as.numeric(gwas_comp$ESTIMATE[match(common_variants, gwas_comp$start)]),
    se2 = as.numeric(gwas_comp$SE[match(common_variants, gwas_comp$start)])
  )
  coloc_res_gwas[[i]] = coloc_wrapper(coloc_input)
  
  gwas_comp = gwas_comp %>% dplyr::select(P, start) %>% mutate(p_value = -log(P,10) / max(-log(P,10)), log_10=-log(P,10))
  gwas_dat[[i]] = gwas_comp
}

coloc_plot_gwas = bind_rows(
  gwas %>% dplyr::select(P_BOLT_LMM, start) %>% mutate(p_value = -log(P_BOLT_LMM,10) / max(-log(P_BOLT_LMM,10)), Group="GWAS", log_10=-log(P_BOLT_LMM,10)),
  bind_rows(gwas_dat, .id="Group")
)
ggplot(coloc_plot_gwas, aes(x=start, y=p_value, color=Group)) + geom_point(alpha=0.5, size=1) + theme_thesis(15) + ylab("-log10(P)") + xlab("")
ggplot(coloc_plot_gwas, aes(x=start, y=log_10, color=Group)) + geom_point(alpha=0.5, size=1) + theme_thesis(15) + ylab("-log10(P)") + xlab("") + geom_hline(yintercept = -log(5e-8, base=10))


# LONGITUDINAL - RS11970286 -----------------------------------------------

# get closest gene from open targets genetics
variant_ix = 5

v2g = get_V2G_data(hits$id[variant_ix])
v2g = v2g[!unlist(lapply(v2g$qtls, is_empty)),] # only take anything with qtls
hits$closest_gene[variant_ix] = v2g$gene[which.min(unlist(lapply(v2g$distances, function(x) x$tissues))),]

# make a granges object for the variant

snp_pos = hits$pos_38[variant_ix]
win = 2e6
region_granges = GenomicRanges::GRanges(
  seqnames = hits$chr[variant_ix], 
  ranges = IRanges::IRanges(start=snp_pos-win, end=snp_pos+win), 
  strand = "*")
region_granges

# "Alasoo_2018","BLUEPRINT","BrainSeq","CEDAR","Fairfax_2012","Fairfax_2014","GENCORD","GEUVADIS","HipSci","Kasela_2017","Lepik_2017","Naranbhai_2015","Nedelec_2016","Quach_2016","ROSMAP","Schmiedel_2018","Schwartzentruber_2018","TwinsUK","van_de_Bunt_2015"
# "macrophage","monocyte","neutrophil","CD4+,T,cell","DLPFC","CD8+,T,cell","transverse,colon","platelet","rectum","B,cell","ileum","LCL","fibroblast","T,cell","iPSC","blood","Tfh,cell","Th17,cell","Th1,cell","Th2,cell","Treg,naive","Treg,memory","CD16+,monocyte","NK,cell","sensory,neuron","adipose","skin","pancreatic,islet"

# pull in eqtl
to_pull = data.frame(
  ensembl_id = "ENSG00000111860",
  study = c("eQTLGen"),
  tissue = c("Blood"),
  source = c("eqtlgen")
)

for(i in 1:dim(to_pull)[1]) {
  
  if(to_pull$source[i]=="api") {
    
    if(to_pull$study[i]=="GTEx_V8") {
      eqtl_df = dplyr::filter(imported_tabix_paths, study==to_pull$study[i], tissue_label==to_pull$tissue[i])
    } else {
      eqtl_df = dplyr::filter(tabix_paths, study==to_pull$study[i], tissue_label==to_pull$tissue[i])
    }
    column_names = colnames(readr::read_tsv(eqtl_df$ftp_path[1], n_max=1))
    summary_stats = import_eqtl_catalog(eqtl_df$ftp_path[1], region_granges, selected_gene_id=to_pull$ensembl_id[i], column_names)
    print(ggplot(summary_stats, aes(x=position, y=-log(pvalue,10))) + geom_point() + geom_vline(xintercept=hits$pos_38[variant_ix]) + theme_thesis(15) + ggtitle(paste(unlist(to_pull[i,]), collapse=", ")))
    write_tsv(summary_stats, paste0("data/", hits$variant[variant_ix] ,"/eqtl_", to_pull$ensembl_id[i], "_", to_pull$study[i], "_", to_pull$tissue[i], ".txt"))
    
  }
  
  if(to_pull$source[i]=="eqtlgen") {
    
    summary_stats = eqtlgen %>% filter(Gene==to_pull$ensembl_id[i])
    summary_stats = makeGRangesFromDataFrame(summary_stats, keep.extra.columns=TRUE, start.field="SNPPos", end.field="SNPPos", seqnames.field="SNPChr")
    summary_stats = lift_over(summary_stats, dir="37_to_38")
    write_tsv(summary_stats, paste0("data/", hits$variant[variant_ix] ,"/eqtl_", to_pull$ensembl_id[i], "_", to_pull$study[i], "_", to_pull$tissue[i], ".txt"))
    
  }
}

# get the imaging summary stats for the locus (this is the same data as the locus zoom plots)

gwas = read.table(pipe(paste0("awk 'NR==1 {print}; $2==", hits$chr[variant_ix], " && $3>", hits$pos_37[variant_ix]-win, " && $3<", hits$pos_37[variant_ix]+win, " {print}' /gpfs01/bhcbio/projects/UK_Biobank/20190102_UK_Biobank_Imaging/Results/GWAS/GWAS_diastolic_BOLT/Results/bolt_", hits$gwas_root[variant_ix], "_full.bgen.stats")), header=TRUE)
gwas = makeGRangesFromDataFrame(gwas, keep.extra.columns=TRUE, start.field="BP", end.field="BP")
gwas = lift_over(gwas, dir="37_to_38")
write_tsv(gwas, paste0("data/", hits$variant[variant_ix], "/gwas.txt"))

ggplot(gwas, aes(x=start, y=-log(as.numeric(P_BOLT_LMM),10))) + geom_point() + geom_vline(xintercept=hits$pos_38[variant_ix]) + theme_thesis(15) + ylab("-log10(P)") + xlab("")

qtl_dat = vector("list", dim(to_pull)[1])
coloc_res = qtl_dat
names(qtl_dat) = apply(to_pull, 1, function(x) paste(x, collapse="_"))
names(coloc_res) = names(qtl_dat)

for(i in 1:length(qtl_dat)) {
  
  if(to_pull$source[i]=="none") next
  
  dat = read_tsv(paste0("data/", hits$variant[variant_ix] ,"/eqtl_", to_pull$ensembl_id[i], "_", to_pull$study[i], "_", to_pull$tissue[i], ".txt"))
  if(nrow(dat)==0) next
  
  if(to_pull$source[i]=="api") {
    dat = dat %>% rename(pvalue="Pvalue", position="start")
    # common_variants = intersect(gwas$start, dat$start)
    # coloc_input = data.frame(
    #   beta1 = as.numeric(gwas$BETA[match(common_variants, gwas$start)]),
    #   se1 = as.numeric(gwas$SE[match(common_variants, gwas$start)]),
    #   beta2 = as.numeric(dat$beta[match(common_variants, dat$start)]),
    #   se2 = as.numeric(dat$se[match(common_variants, dat$start)])
    # )
    # coloc_res[[i]] = coloc_wrapper(coloc_input)
  }
  
  dat = dat %>% dplyr::select(Pvalue, start) %>% mutate(p_value = -log(Pvalue,10) / max(-log(Pvalue,10)), log_10=-log(Pvalue,10))
  qtl_dat[[i]] = dat
}

lapply(coloc_res, function(x) x$posterior)

coloc_plot = bind_rows(
  gwas %>% dplyr::select(P_BOLT_LMM, start) %>% mutate(p_value = -log(P_BOLT_LMM,10) / max(-log(P_BOLT_LMM,10)), Group="GWAS", log_10=-log(P_BOLT_LMM,10)),
  bind_rows(qtl_dat, .id="Group")
)

ggplot(coloc_plot, aes(x=start, y=p_value, color=Group)) + geom_point(alpha=0.5, size=1) + theme_thesis(15) + ylab("-log10(P)") + xlab("")
ggplot(coloc_plot, aes(x=start, y=log_10, color=Group)) + geom_point(alpha=0.5, size=1) + theme_thesis(15) + ylab("-log10(P)") + xlab("") + geom_hline(yintercept = -log(5e-8, base=10))

# phewas plot?
# run phewas and pull signal for anything significant and add to colocalisation step?

snp_mat = getImputedV3GenotypesForVariants(hits$variant[variant_ix])
attr(snp_mat, "metadata")
snpStats::col.summary(snp_mat)
geno_data = getTibbleFromSnpMatrix(snp_mat)

# compute associations

phewas = bind_rows(
  getPheWASResults(geno_data %>% dplyr::select(SID, !!as.name(hits$variant[variant_ix])), slurm=F, cores=8),
  getQuantPheWASResults(geno_data %>% dplyr::select(SID, !!as.name(hits$variant[variant_ix])))
) %>% tbl_df()

phewas %>% arrange(p.value) %>% filter(p.value < 0.05) %>% View()
pheWASForest(phewas %>% filter(p.value < 0.05) %>% dplyr::select(-n) %>% arrange(p.value))

# pulse rate!

traits = c("sbp_adj","pulse_rate","dbp_adj")
gwas_dat = vector("list", length(traits))
coloc_res_gwas = gwas_dat
names(gwas_dat) = traits
names(coloc_res_gwas) = names(gwas_dat)

for(i in 1:length(gwas_dat)) {
  
  gwas_comp = loadGWAS(trait=traits[i], type="linear")
  gwas_comp = as.data.frame(gwas_comp) %>% dplyr::filter(CHR==hits$chr[variant_ix], (BP > hits$pos_37[variant_ix]-win & BP < hits$pos_37[variant_ix]+win)) %>% arrange(BP)
  gwas_comp = makeGRangesFromDataFrame(gwas_comp, keep.extra.columns=TRUE, start.field="BP", end.field="BP")
  gwas_comp = lift_over(gwas_comp, dir="37_to_38")
  
  common_variants = intersect(gwas$start, gwas_comp$start)
  coloc_input = data.frame(
    beta1 = as.numeric(gwas$BETA[match(common_variants, gwas$start)]),
    se1 = as.numeric(gwas$SE[match(common_variants, gwas$start)]),
    beta2 = as.numeric(gwas_comp$ESTIMATE[match(common_variants, gwas_comp$start)]),
    se2 = as.numeric(gwas_comp$SE[match(common_variants, gwas_comp$start)])
  )
  coloc_res_gwas[[i]] = coloc_wrapper(coloc_input)
  
  gwas_comp = gwas_comp %>% dplyr::select(P, start) %>% mutate(p_value = -log(P,10) / max(-log(P,10)), log_10=-log(P,10))
  gwas_dat[[i]] = gwas_comp
}

coloc_plot_gwas = bind_rows(
  gwas %>% dplyr::select(P_BOLT_LMM, start) %>% mutate(p_value = -log(P_BOLT_LMM,10) / max(-log(P_BOLT_LMM,10)), Group="GWAS", log_10=-log(P_BOLT_LMM,10)),
  bind_rows(gwas_dat, .id="Group")
)
ggplot(coloc_plot_gwas, aes(x=start, y=p_value, color=Group)) + geom_point(alpha=0.5, size=1) + theme_thesis(15) + ylab("-log10(P)") + xlab("")
ggplot(coloc_plot_gwas, aes(x=start, y=log_10, color=Group)) + geom_point(alpha=0.5, size=1) + theme_thesis(15) + ylab("-log10(P)") + xlab("") + geom_hline(yintercept = -log(5e-8, base=10))


# LONGITUDINAL - RS10261575 -----------------------------------------------

# get closest gene from open targets genetics
variant_ix = 6

v2g = get_V2G_data(hits$id[variant_ix])
v2g = v2g[!unlist(lapply(v2g$qtls, is_empty)),] # only take anything with qtls
hits$closest_gene[variant_ix] = v2g$gene[which.min(unlist(lapply(v2g$distances, function(x) x$tissues))),]

# make a granges object for the variant

snp_pos = hits$pos_38[variant_ix]
win = 2e6
region_granges = GenomicRanges::GRanges(
  seqnames = hits$chr[variant_ix], 
  ranges = IRanges::IRanges(start=snp_pos-win, end=snp_pos+win), 
  strand = "*")
region_granges

# "Alasoo_2018","BLUEPRINT","BrainSeq","CEDAR","Fairfax_2012","Fairfax_2014","GENCORD","GEUVADIS","HipSci","Kasela_2017","Lepik_2017","Naranbhai_2015","Nedelec_2016","Quach_2016","ROSMAP","Schmiedel_2018","Schwartzentruber_2018","TwinsUK","van_de_Bunt_2015"
# "macrophage","monocyte","neutrophil","CD4+,T,cell","DLPFC","CD8+,T,cell","transverse,colon","platelet","rectum","B,cell","ileum","LCL","fibroblast","T,cell","iPSC","blood","Tfh,cell","Th17,cell","Th1,cell","Th2,cell","Treg,naive","Treg,memory","CD16+,monocyte","NK,cell","sensory,neuron","adipose","skin","pancreatic,islet"

# pull in eqtl
to_pull = data.frame(
  ensembl_id = "ENSG00000111860",
  study = c("eQTLGen"),
  tissue = c("Blood"),
  source = c("eqtlgen")
)

for(i in 1:dim(to_pull)[1]) {
  
  if(to_pull$source[i]=="api") {
    
    if(to_pull$study[i]=="GTEx_V8") {
      eqtl_df = dplyr::filter(imported_tabix_paths, study==to_pull$study[i], tissue_label==to_pull$tissue[i])
    } else {
      eqtl_df = dplyr::filter(tabix_paths, study==to_pull$study[i], tissue_label==to_pull$tissue[i])
    }
    column_names = colnames(readr::read_tsv(eqtl_df$ftp_path[1], n_max=1))
    summary_stats = import_eqtl_catalog(eqtl_df$ftp_path[1], region_granges, selected_gene_id=to_pull$ensembl_id[i], column_names)
    print(ggplot(summary_stats, aes(x=position, y=-log(pvalue,10))) + geom_point() + geom_vline(xintercept=hits$pos_38[variant_ix]) + theme_thesis(15) + ggtitle(paste(unlist(to_pull[i,]), collapse=", ")))
    write_tsv(summary_stats, paste0("data/", hits$variant[variant_ix] ,"/eqtl_", to_pull$ensembl_id[i], "_", to_pull$study[i], "_", to_pull$tissue[i], ".txt"))
    
  }
  
  if(to_pull$source[i]=="eqtlgen") {
    
    summary_stats = eqtlgen %>% filter(Gene==to_pull$ensembl_id[i])
    summary_stats = makeGRangesFromDataFrame(summary_stats, keep.extra.columns=TRUE, start.field="SNPPos", end.field="SNPPos", seqnames.field="SNPChr")
    summary_stats = lift_over(summary_stats, dir="37_to_38")
    write_tsv(summary_stats, paste0("data/", hits$variant[variant_ix] ,"/eqtl_", to_pull$ensembl_id[i], "_", to_pull$study[i], "_", to_pull$tissue[i], ".txt"))
    
  }
}

# get the imaging summary stats for the locus (this is the same data as the locus zoom plots)

gwas = read.table(pipe(paste0("awk 'NR==1 {print}; $2==", hits$chr[variant_ix], " && $3>", hits$pos_37[variant_ix]-win, " && $3<", hits$pos_37[variant_ix]+win, " {print}' /gpfs01/bhcbio/projects/UK_Biobank/20190102_UK_Biobank_Imaging/Results/GWAS/GWAS_diastolic_BOLT/Results/bolt_", hits$gwas_root[variant_ix], "_full.bgen.stats")), header=TRUE)
gwas = makeGRangesFromDataFrame(gwas, keep.extra.columns=TRUE, start.field="BP", end.field="BP")
gwas = lift_over(gwas, dir="37_to_38")
write_tsv(gwas, paste0("data/", hits$variant[variant_ix], "/gwas.txt"))

ggplot(gwas, aes(x=start, y=-log(as.numeric(P_BOLT_LMM),10))) + geom_point() + geom_vline(xintercept=hits$pos_38[variant_ix]) + theme_thesis(15) + ylab("-log10(P)") + xlab("")

qtl_dat = vector("list", dim(to_pull)[1])
coloc_res = qtl_dat
names(qtl_dat) = apply(to_pull, 1, function(x) paste(x, collapse="_"))
names(coloc_res) = names(qtl_dat)

for(i in 1:length(qtl_dat)) {
  
  if(to_pull$source[i]=="none") next
  
  dat = read_tsv(paste0("data/", hits$variant[variant_ix] ,"/eqtl_", to_pull$ensembl_id[i], "_", to_pull$study[i], "_", to_pull$tissue[i], ".txt"))
  if(nrow(dat)==0) next
  
  if(to_pull$source[i]=="api") {
    dat = dat %>% rename(pvalue="Pvalue", position="start")
    # common_variants = intersect(gwas$start, dat$start)
    # coloc_input = data.frame(
    #   beta1 = as.numeric(gwas$BETA[match(common_variants, gwas$start)]),
    #   se1 = as.numeric(gwas$SE[match(common_variants, gwas$start)]),
    #   beta2 = as.numeric(dat$beta[match(common_variants, dat$start)]),
    #   se2 = as.numeric(dat$se[match(common_variants, dat$start)])
    # )
    # coloc_res[[i]] = coloc_wrapper(coloc_input)
  }
  
  dat = dat %>% dplyr::select(Pvalue, start) %>% mutate(p_value = -log(Pvalue,10) / max(-log(Pvalue,10)), log_10=-log(Pvalue,10))
  qtl_dat[[i]] = dat
}

lapply(coloc_res, function(x) x$posterior)

coloc_plot = bind_rows(
  gwas %>% dplyr::select(P_BOLT_LMM, start) %>% mutate(p_value = -log(P_BOLT_LMM,10) / max(-log(P_BOLT_LMM,10)), Group="GWAS", log_10=-log(P_BOLT_LMM,10)),
  bind_rows(qtl_dat, .id="Group")
)

ggplot(coloc_plot, aes(x=start, y=p_value, color=Group)) + geom_point(alpha=0.5, size=1) + theme_thesis(15) + ylab("-log10(P)") + xlab("")
ggplot(coloc_plot, aes(x=start, y=log_10, color=Group)) + geom_point(alpha=0.5, size=1) + theme_thesis(15) + ylab("-log10(P)") + xlab("") + geom_hline(yintercept = -log(5e-8, base=10))

# phewas plot?
# run phewas and pull signal for anything significant and add to colocalisation step?

snp_mat = getImputedV3GenotypesForVariants(hits$variant[variant_ix])
attr(snp_mat, "metadata")
snpStats::col.summary(snp_mat)
geno_data = getTibbleFromSnpMatrix(snp_mat)

# compute associations

phewas = bind_rows(
  getPheWASResults(geno_data %>% dplyr::select(SID, !!as.name(hits$variant[variant_ix])), slurm=F, cores=8),
  getQuantPheWASResults(geno_data %>% dplyr::select(SID, !!as.name(hits$variant[variant_ix])))
) %>% tbl_df()

phewas %>% arrange(p.value) %>% filter(p.value < 0.05) %>% View()
pheWASForest(phewas %>% filter(p.value < 0.05) %>% dplyr::select(-n) %>% arrange(p.value))

# pulse rate!

traits = c("sbp_adj","pulse_rate","dbp_adj")
gwas_dat = vector("list", length(traits))
coloc_res_gwas = gwas_dat
names(gwas_dat) = traits
names(coloc_res_gwas) = names(gwas_dat)

for(i in 1:length(gwas_dat)) {
  
  gwas_comp = loadGWAS(trait=traits[i], type="linear")
  gwas_comp = as.data.frame(gwas_comp) %>% dplyr::filter(CHR==hits$chr[variant_ix], (BP > hits$pos_37[variant_ix]-win & BP < hits$pos_37[variant_ix]+win)) %>% arrange(BP)
  gwas_comp = makeGRangesFromDataFrame(gwas_comp, keep.extra.columns=TRUE, start.field="BP", end.field="BP")
  gwas_comp = lift_over(gwas_comp, dir="37_to_38")
  
  common_variants = intersect(gwas$start, gwas_comp$start)
  coloc_input = data.frame(
    beta1 = as.numeric(gwas$BETA[match(common_variants, gwas$start)]),
    se1 = as.numeric(gwas$SE[match(common_variants, gwas$start)]),
    beta2 = as.numeric(gwas_comp$ESTIMATE[match(common_variants, gwas_comp$start)]),
    se2 = as.numeric(gwas_comp$SE[match(common_variants, gwas_comp$start)])
  )
  coloc_res_gwas[[i]] = coloc_wrapper(coloc_input)
  
  gwas_comp = gwas_comp %>% dplyr::select(P, start) %>% mutate(p_value = -log(P,10) / max(-log(P,10)), log_10=-log(P,10))
  gwas_dat[[i]] = gwas_comp
}

coloc_plot_gwas = bind_rows(
  gwas %>% dplyr::select(P_BOLT_LMM, start) %>% mutate(p_value = -log(P_BOLT_LMM,10) / max(-log(P_BOLT_LMM,10)), Group="GWAS", log_10=-log(P_BOLT_LMM,10)),
  bind_rows(gwas_dat, .id="Group")
)
ggplot(coloc_plot_gwas, aes(x=start, y=p_value, color=Group)) + geom_point(alpha=0.5, size=1) + theme_thesis(15) + ylab("-log10(P)") + xlab("")
ggplot(coloc_plot_gwas, aes(x=start, y=log_10, color=Group)) + geom_point(alpha=0.5, size=1) + theme_thesis(15) + ylab("-log10(P)") + xlab("") + geom_hline(yintercept = -log(5e-8, base=10))


# LONGITUDINAL - RS11535974 -----------------------------------------------



# LONGITUDINAL - RS499715 -------------------------------------------------



# LAV - RS59985551 --------------------------------------------------------

# get closest gene from open targets genetics
variant_ix = 1

v2g = get_V2G_data(hits$id[variant_ix])
hits$closest_gene[variant_ix] = v2g$gene[which.min(unlist(lapply(v2g$distances, function(x) x$tissues))),]

# make a granges object for the variant

snp_pos = hits$pos_38[1]
win = 2e6
region_granges = GenomicRanges::GRanges(
  seqnames = hits$chr[1], 
  ranges = IRanges::IRanges(start=snp_pos-win, end=snp_pos+win), 
  strand = "*")
region_granges

# pull in eqtl

snp_pos = hits$pos_38[variant_ix]
win = 2e6
region_granges = GenomicRanges::GRanges(
  seqnames = hits$chr[variant_ix], 
  ranges = IRanges::IRanges(start=snp_pos-win, end=snp_pos+win), 
  strand = "*")
region_granges

# "Alasoo_2018","BLUEPRINT","BrainSeq","CEDAR","Fairfax_2012","Fairfax_2014","GENCORD","GEUVADIS","HipSci","Kasela_2017","Lepik_2017","Naranbhai_2015","Nedelec_2016","Quach_2016","ROSMAP","Schmiedel_2018","Schwartzentruber_2018","TwinsUK","van_de_Bunt_2015"
# "macrophage","monocyte","neutrophil","CD4+,T,cell","DLPFC","CD8+,T,cell","transverse,colon","platelet","rectum","B,cell","ileum","LCL","fibroblast","T,cell","iPSC","blood","Tfh,cell","Th17,cell","Th1,cell","Th2,cell","Treg,naive","Treg,memory","CD16+,monocyte","NK,cell","sensory,neuron","adipose","skin","pancreatic,islet"

# pull in eqtl
to_pull = data.frame(
  ensembl_id = c("ENSG00000115380"), # candidate is efemp1
  study = c("GTEx_V8"),
  tissue = c("Thyroid")
)

for(i in 1:dim(to_pull)[1]) {
  if(to_pull$study[i]=="GTEx_V8") {
    eqtl_df = dplyr::filter(imported_tabix_paths, study==to_pull$study[i], tissue_label==to_pull$tissue[i])  
  } else {
    eqtl_df = dplyr::filter(tabix_paths, study==to_pull$study[i], tissue_label==to_pull$tissue[i])  
  }
  column_names = colnames(readr::read_tsv(eqtl_df$ftp_path[1], n_max=1))
  summary_stats = import_eqtl_catalog(eqtl_df$ftp_path[1], region_granges, selected_gene_id=to_pull$ensembl_id[i], column_names)
  print(ggplot(summary_stats, aes(x=position, y=-log(pvalue,10))) + geom_point() + geom_vline(xintercept=hits$pos_38[variant_ix]) + theme_thesis(15) + ggtitle(paste(unlist(to_pull[i,]), collapse=", ")))
  write_tsv(summary_stats, paste0("data/", hits$variant[variant_ix] ,"/eqtl_", to_pull$ensembl_id[i], "_", to_pull$study[i], "_", to_pull$tissue[i], ".txt"))
}
# get the imaging summary stats for the locus (this is the same data as the locus zoom plots)

gwas = read.table(pipe(paste0("awk 'NR==1 {print}; $2==", hits$chr[variant_ix], " && $3>", hits$pos_37[variant_ix]-win, " && $3<", hits$pos_37[variant_ix]+win, " {print}' /gpfs01/bhcbio/projects/UK_Biobank/20190102_UK_Biobank_Imaging/Results/GWAS/GWAS_diastolic_BOLT/Results/bolt_", hits$gwas_root[variant_ix], "_full.bgen.stats")), header=TRUE)
write_tsv(gwas, paste0("data/", hits$variant[variant_ix], "/gwas.txt"))

gwas = makeGRangesFromDataFrame(gwas, keep.extra.columns=TRUE, start.field="BP", end.field="BP")
gwas = lift_over(gwas, dir="37_to_38")

ggplot(gwas, aes(x=start, y=-log(as.numeric(P_BOLT_LMM),10))) + geom_point() + geom_vline(xintercept=hits$pos_38[variant_ix]) + theme_thesis(15) + ylab("-log10(P)") + xlab("")

summary_stats = read_tsv(paste0("data/", hits$variant[variant_ix] ,"/eqtl_", to_pull$ensembl_id[i], "_", to_pull$study[i], "_", to_pull$tissue[i], ".txt"))

qplot(-log(as.numeric(gwas$P_BOLT_LMM),10), -log(summary_stats$pvalue,10)[match(gwas$start, summary_stats$position)]) + theme_thesis(15) + xlab("GWAS") + ylab("eQTL")

coloc_plot = bind_rows(
  GWAS = gwas %>% dplyr::select(P_BOLT_LMM, start) %>% mutate(p_value = -log(P_BOLT_LMM,10) / max(-log(P_BOLT_LMM,10))),
  QTL = summary_stats %>% dplyr::select(pvalue, position) %>% rename(position="start") %>% mutate(p_value = -log(pvalue,10) / max(-log(pvalue,10))),
  .id = "Group"
)

ggplot(coloc_plot, aes(x=start, y=p_value, color=Group)) + geom_point(alpha=0.5, size=1) + theme_thesis(15) + ylab("-log10(P)") + xlab("")

common_variants = intersect(gwas$start, summary_stats$position)
coloc_input = data.frame(
  beta1 = as.numeric(gwas$BETA[match(common_variants, gwas$start)]),
  se1 = as.numeric(gwas$SE[match(common_variants, gwas$start)]),
  beta2 = as.numeric(summary_stats$beta[match(common_variants, summary_stats$position)]),
  se2 = as.numeric(summary_stats$se[match(common_variants, summary_stats$position)])
)

coloc_wrapper(coloc_input)

# phewas plot?
# run phewas and pull signal for anything significant and add to colocalisation step?

snp_mat = getImputedV3GenotypesForVariants(hits$variant[variant_ix])
attr(snp_mat, "metadata")
snpStats::col.summary(snp_mat)
geno_data = getTibbleFromSnpMatrix(snp_mat)

# compute associations

phewas = bind_rows(
  getPheWASResults(geno_data %>% dplyr::select(SID, !!as.name(hits$variant[1])), slurm=F, cores=8),
  getQuantPheWASResults(geno_data %>% dplyr::select(SID, !!as.name(hits$variant[1])))
) %>% tbl_df()

phewas %>% arrange(p.value) %>% filter(p.value < 0.05) %>% View()
pheWASForest(phewas %>% filter(p.value < 0.05) %>% dplyr::select(-n) %>% arrange(p.value))


