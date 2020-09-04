# 'gwas_res', list of dfs with 7 columns
# snp, chr, bp, p, estimate, se, n

locus_zoom_region <- function(gwas_res, CHR, start, stop, rdir=".", sig_line1=7.3, sig_line2=5, clean=FALSE) {
  
  nodes = length(CHR)
  time = ceiling(length(CHR)/16) * 15 + 1
  signiLine = paste0("\"", sig_line1, ",", sig_line2, "\"")
  dir.create(paste(rdir, "/locuszoom/", sep=""))
  lzbash = paste(rdir, "/locuszoom/run_Locuszoom.sh", sep="")
  job.out = paste(rdir, "/locuszoom/myjob-%J.out", sep="")
  
  # create bash file
  cat("#!/bin/bash", file=lzbash, sep="\n")
  cat(paste("#SBATCH --partition=\"medium_priority\"", sep=" "), file=lzbash, sep="\n", append=TRUE)
  cat(paste("#SBATCH --time=\"", time, "\"", sep=""), file=lzbash, sep="\n", append=TRUE)
  cat(paste("#SBATCH -n", nodes, sep=" "), file=lzbash, sep="\n", append=TRUE)
  
  gwas_res$P <- as.numeric(gwas_res$P)
  nulls = which(gwas_res$P == 0)
  if(length(nulls) > 0) {
    gwas_res$P[nulls] = .Machine$double.xmin
  }
  
  gwas_res$P <- as.numeric(gwas_res$P)
  nulls = which(gwas_res$P == 0)
  if(length(nulls) > 0) {
    gwas_res$P[nulls] = .Machine$double.xmin
  }
  
  path = paste0(rdir, "/locuszoom/metal_file")
  write.table(gwas_res, file=path, sep="\t", quote=FALSE, row.names=F, col.names=TRUE)
  
  for (i in 1:length(CHR)) {
    
    lz.prefix = paste(rdir, "/locuszoom/", sep="")
    cmd = paste("srun --job-name", paste0(CHR[i], start[i], stop[i], sep="_"), "-o", paste0("'", job.out, "'"), "-N 1 -n 1 --exclusive", "bash -c 'module load R; module load PLINK2/1.9b_4.1-x86_64;", "/gpfs01/sw/eb20150106/software/locuszoom/1.3-foss-2014b/bin/locuszoom", "--metal", path, "--markercol SNP --pvalcol P --delim tab --chr", CHR[i], "--start", start[i], "--end", stop[i], "--source 1000G_March2012 --cache None", "--plotonly --pop EUR --build hg19 signifLine=", signiLine, "signifLineColor=\"red,blue\" --prefix", paste0("\"", lz.prefix, "\""), "' &", sep = " ")
    cat(cmd, file=lzbash, sep="\n", append=TRUE)
    
  }
  
  if(clean == TRUE) {
    cat("wait", file=lzbash, sep="\n", append=TRUE)
  }
  
  Sys.chmod(lzbash, mode="0777")
  system(paste0("\"", lzbash, "\""))
  
  if(clean == TRUE) {
    file.list = c(list.files(paste(rdir, "/locuszoom", sep=""), pattern="myjob", full.names=TRUE), list.files(paste(rdir, "/locuszoom", sep=""), pattern="metal", full.names=TRUE), list.files(paste(rdir, "/locuszoom", sep=""), pattern=".sh", full.names=TRUE))
    invisible(file.remove(file.list))
  }
  
  if(clean == FALSE) {
    print("Please clean up all files manually; specifically the metal file (the GWAS results)!")
  }
  
}


