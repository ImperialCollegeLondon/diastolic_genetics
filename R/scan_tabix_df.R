#' A general function to quickly import tabix indexed tab-separated files into data_frame
#'
#' @param tabix_file Path to tabix-indexed text file
#' @param param An instance of GRanges, RangedData, or RangesList
#' provide the sequence names and regions to be parsed. Passed onto Rsamtools::scanTabix()
#' @param ... Additional parameters to be passed on to readr::read_delim()
#'
#' @return List of data_frames, one for each entry in the param GRanges object.
#' @export

scan_tabix_df <- function(tabix_file, param, ...) {
  
  tabix_list = Rsamtools::scanTabix(tabix_file, param=param)
  df_list = lapply(tabix_list, function(x, ...) {
    if(length(x) > 0) {
      if(length(x) == 1) {
        
        # hack to make sure that it also works for data frames with only one row
        # adds an empty row and then removes it
        
        result = paste(paste(x, collapse = "\n"),"\n",sep = "")
        result = readr::read_delim(result, delim = "\t", ...)[1,]
      } else {
        result = paste(x, collapse="\n")
        result = readr::read_delim(result, delim="\t", ...)
      }
    } else {
      # return NULL if the nothing is returned from tabix file
      result = NULL
    }
    return(result)
  }, ...)
  
  return(df_list)
  
}

