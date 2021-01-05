filter_bgen_to_bed_temp <- function (chrom, rs_file, excludedSamples=NULL, outdir) {
  
  if(missing(excludedSamples)) {
    excludedSamples = getExcludedSamplesFromWhiteBritish()
  }
  
  bgen_file = file.path(config::get("imputation_dir"), paste0("_004_ukb_imp_chr", chrom, "_v3.bgen"))
  filtered_bgen = tempfile(paste0("filtered_", stringr::str_replace(base::basename(rs_file), "\\.", "_")), outdir, fileext=".bgen")
  
  cmd = config::get("bgenix_exe")
  cmd = paste0(cmd, " -g ", bgen_file, " -incl-rsids ", rs_file, " > ", filtered_bgen)
  system(cmd)
  print(cmd)
  
  plinkfile = tempfile("plink", outdir)
  cmd <- config::get("plink20_exe")
  cmd = paste0(cmd, " --out ", plinkfile, " ")
  cmd = paste0(cmd, "--bgen ", filtered_bgen, " ")
  cmd = paste0(cmd, "--sample ", config::get("sample_autosome_file"), " ")
  cmd = paste0(cmd, "--remove-fam ", config::get("withdrawn_sample_file"), " ")
  cmd = paste0(cmd, "--real-ref-alleles ")
  cmd = paste0(cmd, "--rm-dup exclude-all ")
  cmd = paste0(cmd, "--make-bed")
  system(cmd, ignore.stdout = TRUE)
  return(plinkfile)
  
}

