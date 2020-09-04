# description of function: https://cnsgenomics.com/software/gcta/#COJO (including desirable structure of ma_file)

conditional_analysis <- function(ma_file, dir, thres=5*10^-8, chr) {
  
  # create a tempdir within dir that will be deleted again afterwards
  rn = sample(1:1000000,1)
  temp_dir <- paste0(dir, "/temp__", rn)
  dir.create(temp_dir)
  
  # write .ma file
  data.table::fwrite(ma_file, paste0(temp_dir, "/madata"), sep="\t")
  
  # export snp file
  snps_region = matrix(ma_file$SNP, ncol=1)
  write.table(snps_region, paste0(temp_dir, "/snpsinclude.txt"), col.names=F, row.names=F, quote=F)
  plink_file = filter_bgen_to_bed(chrom=chr, rs_file=paste0(temp_dir, "/snpsinclude.txt"), excludedSamples=NULL, outdir=temp_dir)
  
  # list of snps to condition on
  snps_cond = matrix(ma_file %>% dplyr::arrange(p) %>% dplyr::slice(1) %>% dplyr::pull(SNP), ncol=1)
  
  # first line of export results:
  res <- ma_file %>% dplyr::arrange(p) %>% dplyr::slice(1) %>% dplyr::select(SNP, p, b) %>% dplyr::mutate(pC=NA, bC=NA)
  o = 1
  
  while(o > 0) {
    
    o = 0
    write.table(snps_cond, paste0(temp_dir, "/snpscondition.txt"), col.names=F, row.names=F, quote=F)
    
    # analyse data
    cmd = paste0("/gpfs01/bhcbio/projects/UK_Biobank/20170616_UK_Biobank_Genotyping/Programs/gcta_1.93.1beta/gcta64 --bfile ", plink_file, " --cojo-file ", paste0(temp_dir, "/madata"), " --cojo-cond ", paste0(temp_dir, "/snpscondition.txt"), " --out ", paste0(temp_dir, "/conditional"))
    system(cmd)
    
    # re-import output
    out = data.table::fread(paste0(temp_dir, "/conditional.cma.cojo")) %>% dplyr::filter(!is.na(pC))
    if(any(out$pC < thres)) {
      res = rbind(res, out %>% dplyr::arrange(pC) %>% dplyr::slice(1) %>% dplyr::select(SNP, p, b, pC, bC))
      snps_cond = rbind(snps_cond, matrix(out %>% dplyr::arrange(pC) %>% dplyr::slice(1) %>% dplyr::pull(SNP), ncol=1))
      o = 1
    }
    
  }
  
  unlink(temp_dir, recursive=TRUE)
  return(res)
  
}

