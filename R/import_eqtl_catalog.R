
import_eqtl_catalog <- function(ftp_path, region, selected_gene_id, column_names, verbose=TRUE) {
  
  if(verbose){
    print(ftp_path)
  }
  
  # fetch summary statistics with Rsamtools
  summary_stats = scan_tabix_df(ftp_path, region, col_names=column_names)[[1]] %>% dplyr::filter(gene_id==selected_gene_id)
  
  # remove rsid duplicates and multi-allelic variant
  summary_stats = dplyr::select(summary_stats, -rsid) %>% 
    dplyr::distinct() %>% # mrsid duplicates
    dplyr::mutate(id = paste(chromosome, position, sep=":")) %>% 
    dplyr::group_by(id) %>% 
    dplyr::mutate(row_count = n()) %>% dplyr::ungroup() %>% 
    dplyr::filter(row_count == 1) # ultialllics
  
  return(summary_stats)
  
}

