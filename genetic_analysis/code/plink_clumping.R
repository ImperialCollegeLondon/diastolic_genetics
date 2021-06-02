############  Do plink clumping ############
#'
#' @param plink_basename A bfile base name with genotype information used for LD calculations
#' @param score_file A GWAS result file with columns SNP and P that is used for clumping
#' @param outdir Where to store temporary files
#' @param clump_p1 The clump-p1 argument for plink, i.e. the maximum p-value for a lead SNP within a clump
#' @param clump_kb variants which are >= clump_kb distance or have LD r2<0.1 with the existing index variants are assigned to a new index variants. Details see https://www.cog-genomics.org/plink/1.9/postproc#clump
#' @param clump_r2 variants which are >= clump_kb distance or have LD r2<0.1 with the existing index variants are assigned to a new index variants. Details see https://www.cog-genomics.org/plink/1.9/postproc#clump
#'
#' @return A data table with lead SNPs and p-values for each clump.
#' @export
plink_clumping <- function(plink_basename, score_file, outdir, clump_p1, clump_kb, clump_r2) {
  plinkfile = tempfile("plink", outdir)

  cmd <- config::get("plink19_exe")
  cmd = paste0(cmd, ' --out "', plinkfile, '" ')
  cmd = paste0(cmd, ' --bfile "', plink_basename, '" ')

  cmd = paste0(cmd, '--remove-fam "', config::get("withdrawn_sample_file"), '" ')
  cmd = paste0(cmd, '--clump "', score_file, '" ')
  cmd = paste0(cmd, '--clump-snp-field ID ')
  cmd = paste0(cmd, '--clump-p1 ', clump_p1, ' ')
  cmd = paste0(cmd, '--clump-p2 0.001 ')
  cmd = paste0(cmd, '--clump-kb ', clump_kb, ' ')
  cmd = paste0(cmd, '--clump-r2 ', clump_r2, ' ')

  system(cmd, ignore.stdout = TRUE)

  # Read lead SNPs and p-values from clumps
  ret = data.table::fread((list.files(path = outdir, pattern = paste0(basename(plinkfile), ".clumped"), full.names = T))[1]) %>% dplyr::select(SNP, P)

  # remove temporary files
  purrr::map(list.files(path = dirname(plink_basename), pattern = paste0(basename(plink_basename), ".*"), full.names = T), file.remove)

  return(ret)
}

