# 'gwas_res', list of dfs with 7 columns
# snp, chr, bp, p, estimate, se, n

locus_zoom_region_temp <- function(gwas_res, CHR, start, stop, rdir=".", sig_line1=7.3, sig_line2=5, clean=FALSE) {
  
  signiLine = paste0("\"", sig_line1, ",", sig_line2, "\"")
  dir.create(paste(rdir, "/locuszoom/", sep=""))
  
  gwas_res$P <- as.numeric(gwas_res$P)
  nulls = which(gwas_res$P == 0)
  if(length(nulls) > 0) {
    gwas_res$P[nulls] = .Machine$double.xmin
  }
  
  path = paste0(rdir, "/locuszoom/metal_file")
  # write.table(gwas_res, file=path, sep="\t", quote=FALSE, row.names=F, col.names=TRUE)
  
  for (i in 1:length(CHR)) {
  
    lzbash = paste(rdir, "/locuszoom/run_Locuszoom_", i, ".sh", sep="")
    job.out = paste(rdir, "/locuszoom/myjob-%J.out", sep="")
    
    # create bash file
    cat("#!/bin/bash", file=lzbash, sep="\n")
    lz.prefix = paste(rdir, "/locuszoom/", paste(CHR[i], start[i], stop[i], sep="_"), sep="")
 
    cmd = paste("module load R; module load PLINK2/1.9b_4.1-x86_64;", "/gpfs01/sw/eb20150106/software/locuszoom/1.3-foss-2014b/bin/locuszoom", "--metal", path, "--markercol SNP --pvalcol P --delim tab --chr", CHR[i], "--start", start[i], "--end", stop[i], "--source 1000G_March2012 --cache None", "--plotonly --pop EUR --build hg19 signifLine=", signiLine, "signifLineColor=\"red,blue\" --prefix", paste0("\"", lz.prefix, "\""), " &", sep=" ")
    cat(cmd, file=lzbash, sep="\n", append=TRUE)
    
  }
  
}


