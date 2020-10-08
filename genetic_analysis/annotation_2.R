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

mapping = id_mapping()

scale_this <- function(x) {
  return(x / max(x))
}


# VARIANT DATA ------------------------------------------------------------

load("r_data/hits.RData")
variant_ix = 5
list.files(paste0("data/", hits$variant[variant_ix]))
snp_pos = hits$pos_38[variant_ix]
win = 2e6

v2g = get_V2G_data(hits$id[variant_ix])
v2g = v2g[!unlist(lapply(v2g$qtls, is_empty)),]

region_granges = GenomicRanges::GRanges(
  seqnames = hits$chr[variant_ix], 
  ranges = IRanges::IRanges(start=snp_pos-win, end=snp_pos+win), 
  strand = "*")
region_granges


# IMAGING REGION ----------------------------------------------------------

gwas = read_tsv(paste0("data/", hits$variant[variant_ix], "/gwas.txt"))


# EQTL DATA ---------------------------------------------------------------

source("eqtl_data.R")

qtl_dat = vector("list", dim(to_pull[[variant_ix]])[1])
names(qtl_dat) = apply(to_pull[[variant_ix]], 1, function(x) paste(x, collapse="_"))

for(i in 1:length(qtl_dat)) {
  if(to_pull[[variant_ix]]$source[i]=="none") next
  dat = read_tsv(paste0("data/", hits$variant[variant_ix] ,"/eqtl_", to_pull[[variant_ix]]$ensembl_id[i], "_", to_pull[[variant_ix]]$study[i], "_", to_pull[[variant_ix]]$tissue[i], ".txt"))
  if(to_pull[[variant_ix]]$source[i]=="api") {
    dat = dat %>% rename(pvalue="Pvalue", position="start")
  }
  dat = dat %>% dplyr::select(Pvalue, start) %>% mutate(p_value = -log(Pvalue,10) / max(-log(Pvalue,10)))
  qtl_dat[[i]] = dat
}

coloc_plot = bind_rows(
  gwas %>% dplyr::select(P_BOLT_LMM, start) %>% mutate(p_value = -log(P_BOLT_LMM,10), Group="GWAS"),
  bind_rows(qtl_dat, .id="Group")
)

coloc_plot = coloc_plot %>% group_by(Group) %>% mutate(p_value_scaled=scale_this(p_value))
p1 = ggplot(coloc_plot, aes(x=start, y=p_value, color=Group)) + geom_point(alpha=0.5, size=1) + theme_thesis(15) + ylab("-log10(P)") + xlab("") + geom_hline(yintercept = -log(5e-8, base=10), alpha=0.5, lty=2, color="grey") + geom_vline(xintercept=hits$pos_38[variant_ix], alpha=0.5, color="grey", lty=2) + theme(legend.position="None")
p2 = ggplot(coloc_plot, aes(x=start, y=p_value_scaled, color=Group)) + geom_point(alpha=0.5, size=1) + theme_thesis(15) + ylab("-log10(P) Scaled") + xlab("") + geom_vline(xintercept=hits$pos_38[variant_ix], alpha=0.5, color="grey", lty=2) + theme(legend.position="None")

save(coloc_plot, file=paste0("data/", hits$variant[variant_ix], "/coloc_plot.RData"))
load(file=paste0("data/", hits$variant[variant_ix], "/coloc_plot.RData"))


# TRAITS OF INTEREST ------------------------------------------------------

traits = c("sbp_adj","pulse_rate","dbp_adj")
gwas_dat = vector("list", length(traits))
names(gwas_dat) = traits

for(i in 1:length(gwas_dat)) {
  
  gwas_comp = loadGWAS(trait=traits[i], type="linear")
  gwas_comp = as.data.frame(gwas_comp) %>% dplyr::filter(CHR==hits$chr[variant_ix], (BP > hits$pos_37[variant_ix]-win & BP < hits$pos_37[variant_ix]+win)) %>% arrange(BP)
  gwas_comp = makeGRangesFromDataFrame(gwas_comp, keep.extra.columns=TRUE, start.field="BP", end.field="BP")
  gwas_comp = lift_over(gwas_comp, dir="37_to_38")
  gwas_comp = gwas_comp %>% dplyr::select(P, start) %>% mutate(p_value_scaled = -log(P,10) / max(-log(P,10)), p_value = -log(P,10))
  gwas_dat[[i]] = gwas_comp
}

coloc_plot_gwas = bind_rows(
  gwas %>% dplyr::select(P_BOLT_LMM, start) %>% mutate(p_value_scaled = -log(P_BOLT_LMM,10) / max(-log(P_BOLT_LMM,10)), Group="GWAS", p_value = -log(P_BOLT_LMM,10)),
  bind_rows(gwas_dat, .id="Group")
)

save(coloc_plot_gwas, file=paste0("data/", hits$variant[variant_ix], "/coloc_plot_gwas.RData"))
load(file=paste0("data/", hits$variant[variant_ix], "/coloc_plot_gwas.RData"))


# GENE ANNOTATION ---------------------------------------------------------

win = 1e6
load("~/OneDrive - Bayer/code/bullseye/r_data/t_list.RData")
t_list_filtered = t_list %>% filter(chromosome_name==hits$chr[variant_ix], exon_chrom_start > hits$pos_38[variant_ix]-win, exon_chrom_start < hits$pos_38[variant_ix]+win)
pick_t = t_list_filtered %>% group_by(external_gene_name) %>% summarise(pick_t=ensembl_transcript_id[1]) %>% dplyr::select(pick_t) %>% unlist()
t_list_filtered = t_list_filtered %>% dplyr::filter(ensembl_transcript_id %in% pick_t)
t_list_filtered$strand = "*"
t_list_filtered$ensembl_transcript_id = t_list_filtered$external_gene_name
gene_track = makeGRangesFromDataFrame(t_list_filtered, keep.extra.columns=TRUE)
gene_track$model = "exon"

p3 = ggbio::autoplot(gene_track, aes(tgene_trackpe=model, group=ensembl_transcript_id)) + theme_thesis(10, angle_45=FALSE) + geom_vline(xintercept=hits$pos_38[variant_ix], alpha=0.5, color="grey", lty=2) 
p3


# FINAL PLOTS -------------------------------------------------------------

p4 = ggplot(coloc_plot_gwas, aes(x=start, y=p_value, color=Group)) + geom_point(alpha=0.5, size=1) + theme_thesis(15) + ylab("-log10(P)") + xlab("") + geom_hline(yintercept = -log(5e-8, base=10), alpha=0.5, lty=2, color="grey") + geom_vline(xintercept=hits$pos_38[variant_ix], alpha=0.5, color="grey", lty=2) + theme(legend.position="None")
p5 = ggplot(coloc_plot_gwas, aes(x=start, y=p_value_scaled, color=Group)) + geom_point(alpha=0.5, size=1) + theme_thesis(15) + ylab("-log10(P) Scaled") + xlab("") + geom_vline(xintercept=hits$pos_38[variant_ix], alpha=0.5, color="grey", lty=2) + theme(legend.position="None")

pdf(file=paste0("tex/images/", hits$variant[variant_ix], "_coloc_plot_gwas.pdf"))
ggbio::tracks(p4, p5, p3, heights=c(1,1,1))
dev.off()

pdf(file=paste0("tex/images/", hits$variant[variant_ix], "_coloc_plot.pdf"))
ggbio::tracks(p1, p2, p3, heights=c(1,1,1))
dev.off()


