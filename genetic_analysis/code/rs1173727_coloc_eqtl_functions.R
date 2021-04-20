#### Coloc functions ####

#' Colocalisation test 
#' 
#' Based on \url{https://github.com/tobyjohnson/gtx/blob/master/R/coloc.R}
#' A test for the support of overlapping GWAS signals. The result is the posterior probabilities 
#' for H0 (no signals), H1 (a GWAS signal for the first sample only), H2 (a GWAS signal for the 
#' second sample only), H1,2 (a GWAS signal in both samples but not overlapping), and H12 (a 
#' colocalisation in the 2 signals).
#'
#' @param coloc_input a data frame with beta1, se1, beta2, se2
#' @param priorsd1 Double. Default \code{1}
#' @param priorsd2 Double. Default \code{1} 
#' @param priorc1  Double. Default \code{1e-4}
#' @param priorc2  Double. Default \code{1e-4}
#' @param priorc12  Double. Default \code{1e-5}
#' @param join_type String. Default \code{'inner'}
#' @param summary_only Boolean. Default \code{FALSE}
#' 
#' @return List with six elements:
#'   \code{list(nv=nv, prior=prior, bf=bf, posterior=posterior, alpha12=alpha12, alpha21=alpha21}
#'    Contains posteriors for hypotheses H0, H1, H2, H1,2, H12 in element posterior.
coloc_wrapper <- function(coloc_input, priorsd1=1, priorsd2=1, 
                          priorc1=1e-4, priorc2=1e-4, priorc12=1e-5, 
                          join_type='inner', summary_only=FALSE) {
  
  # functions needed by coloc_compute
  
  norm1 <- function(x, log=FALSE) {
    if (all(is.na(x))) return(x)
    if (log) {
      x <- x - max(x, na.rm=TRUE)
      x <- exp(x)    
    } else {
      # this does not work if x contains nans or +infs
      stopifnot(all(x >= 0, na.rm=TRUE))
      x <- x / max(x, na.rm=TRUE)
    }
    return(x / sum(x, na.rm=TRUE))
  }
  
  abf.Wakefield <- function(beta, se, priorsd, log=FALSE) {
    if (log) {
      return(log(sqrt(se^2/(se^2 + priorsd^2))) + 
               (beta/se)^2/2 * priorsd^2/(se^2 + priorsd^2))
    } else {
      return(sqrt(se^2/(se^2+priorsd^2)) * exp((beta/se)^2/2 * priorsd^2/(se^2+priorsd^2)))
    }
  }
  
  # method
  stopifnot(is.data.frame(coloc_input))
  stopifnot(all(c('beta1', 'se1', 'beta2', 'se2') %in% names(coloc_input)))
  stopifnot(join_type %in% c('inner', 'left', 'right', 'outer'))
  
  coloc_input <- within(coloc_input, {
    lnabf1 <- abf.Wakefield(beta1, se1, priorsd1, log=TRUE)
    lnabf2 <- abf.Wakefield(beta2, se2, priorsd2, log=TRUE)
    inc1 <- if (join_type == 'inner' | join_type == 'left') !is.na(lnabf1) else TRUE
    inc2 <- if (join_type == 'inner' | join_type == 'right') !is.na(lnabf2) else TRUE
    inc <- inc1 & inc2
  })
  
  coloc_input <- coloc_input[which(coloc_input$inc),]
  coloc_input <- within(coloc_input, {
    inc1 <- NULL
    inc2 <- NULL
    inc <- NULL
  })
  
  abf1 <- norm1(c(0, coloc_input$lnabf1), log=TRUE)
  abf2 <- norm1(c(0, coloc_input$lnabf2), log=TRUE)
  abf1[is.na(abf1)] <- 0
  abf2[is.na(abf2)] <- 0
  coloc_input$abf1 <- abf1[-1]
  coloc_input$abf2 <- abf2[-1]
  attr(coloc_input, 'nullabf1') <- abf1[1]
  attr(coloc_input, 'nullabf2') <- abf2[1]
  
  nv <- nrow(coloc_input)
  
  # compute colocalization model probabilities
  prior <- norm1(c(1, priorc1*nv, priorc2*nv, priorc1*priorc2*nv*(nv-1), priorc12*nv))
  names(prior) <- c('H0', 'H1', 'H2', 'H1,2', 'H12')
  
  bf = if (nv > 0) c(
    abf1[1] * abf2[1], 
    sum(abf1[-1]) * abf2[1]/nv, 
    abf1[1] * sum(abf2[-1])/nv, 
    (sum(abf1[-1])*sum(abf2[-1]) - sum(abf1[-1]*abf2[-1]))/(nv*(nv - 1)), 
    sum(abf1[-1] * abf2[-1])/nv
  ) else rep(NA, 5)
  
  bf <- bf/max(bf)
  names(bf) <- names(prior)
  posterior <- norm1(prior*bf)
  names(posterior) <- names(prior)
  
  # compute model averaged effect size ratios
  mw <- coloc_input$abf1*coloc_input$abf2 # model weights
  alpha12 <- sum(mw*coloc_input$beta1/coloc_input$beta2)/sum(mw)
  alpha21 <- sum(mw*coloc_input$beta2/coloc_input$beta1)/sum(mw)
  
  return(list(nv=nv, prior=prior, bf=bf, posterior=posterior, alpha12=alpha12, alpha21=alpha21))
  
}


#' Colocalization worker function
#'
#' @param tab1 Data.table. Requires columns \code{SNP, BP, ESTIMATE, SE}
#' @param tab2 Data.table. Requires columns \code{SNP, BP, ESTIMATE, SE} 
#' @param pstart Integer. Start position
#' @param pend Integer. End position
#'
#' @return list with two elements \code{coloc, result}. The former is the 
#'   full result object returned by function \code{coloc_wrapper}, the latter
#'   the rounded posterior probabilities
make_coloc <- function(tab1, tab2, pstart, pend) {
  lsel <- tab1[BP >= pstart & BP <= pend]
  ukbb <- tab2[BP >= pstart & BP <= pend]
  coloc_input <- tibble::as_tibble(merge(lsel, ukbb, by="SNP"))
  coloc_input <- coloc_input %>% dplyr::rename(beta1 = ESTIMATE.x, se1 = SE.x, beta2 = ESTIMATE.y, se2 = SE.y) 
  #rownames(coloc_input) <- coloc_input$SNP
  coloc_input <- coloc_input %>% dplyr::select(beta1, se1, beta2, se2)
  coloc_input <- unique(coloc_input)
  cloc <- coloc_wrapper(coloc_input, priorsd1=1, priorsd2=1, priorc1=1e-4, priorc2=1e-4, priorc12=1e-5, join_type='inner', summary_only=FALSE)
  # return(list(coloc=cloc, result=format(cloc$posterior, scientific=FALSE, digits=5)))
  return(list(coloc=cloc, result=round(cloc$posterior, digits=4)))
}


#' Fetch all data to be put into coloclization analysis. Note: BAYER SPECIFIC
#' 
#' Note: This function is specific to the Bayer internal data storage. To run a similar analysis,
#' one has to adjust this function in a way that \code{coloc_partner} and \code{coloc_data} are 
#' fetched and saved into the returned object. That means, \code{}
#' 
#' @param gene_info Data frame with columns start and end and chr, holding genomic positions 
#'   of target region to be fetched.
#' @param flanking Integer. Retrieve data in region plus/minus flanking base pairs.
#' @param subset Integer vector. Can be used to subset available UKBB traits.
#'
#' @return Returns a list: 
#'   \code{list(coloc_data=list(trait1=data.table(SNP="rsxx", 
#'                                                CHR=5,
#'                                                BP=32828846,
#'                                                P=0.0085,
#'                                                ESTIMATE=1.42272e-02,
#'                                                SE=0.004268,
#'                                                N=25874),
#'                              trait2=data.table(...))), 
#'              coloc_partners=list(ref1=data.table(...), 
#'                                  ref2=data.table(...)),
#'              lst=lst,
#'              len=len,
#'              chr=chr)}.\cr
#'    \code{coloc_data} is a list with one data.table for each retrieved phenotype. Names of the list are 
#'      the names of the phenotypes. Each data.table requires columns \code{SNP,CHR,BP,P,ESTIMATE,SE,N}\cr
#'    \code{coloc_partners} List with the same format as \code{coloc_data}, holding the 'reference' or 
#'      partner summary stats for the colocalization test. Here, we compared all phenotypes to the 
#'      three blood pressure related traits \code{diastolic, systolic and mean arterial pressure} (dbp, sbp, map)\cr
#'    \code{lst,len,chr} are start, end and chromosome of genomic region.
coloc_get_my_data <- function(gene_info,
                              flanking=200000,
                              subset=NULL) {
  
  # retrieve start and end position of region by gene information
  lst <- as.numeric(gene_info$start) - flanking
  len <- as.numeric(gene_info$end) + flanking
  chr <- as.numeric(gene_info$chr)
  
  # coloc partner gwas stats
  # !! BAYER SPECIFIC FUNCTIONS !!
  ukbb_dbp <- UKBRlib::loadGWAS(trait = "dbp_adj", chromosome = chr, start = lst, end = len, type = "linear")
  ukbb_sbp <- UKBRlib::loadGWAS(trait = "sbp_adj", chromosome = chr, start = lst, end = len, type = "linear")
  ukbb_map <- UKBRlib::loadGWAS(trait = "map_adj", chromosome = chr, start = lst, end = len, type = "linear")
  
  # all traits available
  # !! BAYER SPECIFIC FUNCTIONS !!
  mytr <- setdiff(plink_get_traits(),"")
  if(!is.null(subset)) {
    mytr <- mytr[subset]
  }
  
  # which are quantitative
  # !! BAYER SPECIFIC FUNCTIONS !!
  myqtr <- UKBRlib::getQuantitativeTraitMetaData()$Trait
  
  # assemble data object
  coloc_data <- list()
  for(i in 1:length(mytr)) {
    mytrx <- mytr[i]
    # skip over these traits, are the 'coloc partners'
    if(mytrx %in% c("dbp", "sbp", "map", "dbp_adj", "sbp_adj", "map_adj")) {
      next
    }
    print(sprintf("%s/%s: %s", i, length(mytr), mytrx))
    mytype <- ifelse(mytrx %in% myqtr, "linear", "logistic")
    
    # !! BAYER SPECIFIC FUNCTIONS !!
    xx <- try(ukbb <- UKBRlib::loadGWAS(trait = mytrx, chromosome = chr, start = lst, end = len, type = mytype))
    if("try-error" %in% class(xx)) {
      next
    }
    
    coloc_data[[mytrx]] <- ukbb
  }
  return(list(coloc_data=coloc_data, 
              coloc_partners=list(ukbb_dbp=ukbb_dbp,
                                  ukbb_sbp=ukbb_sbp,
                                  ukbb_map=ukbb_map),
              lst=lst, len=len, chr=chr))
}

#' Run colocalization analyses given summary stats for target and reference traits
#'
#' @param input_data List generated as described in \code{coloc_get_bayer_data}. Format example:\cr
#'   \code{list(coloc_data=list(trait1=data.table(SNP="rsxx", 
#'                                                CHR=5,
#'                                                BP=32828846,
#'                                                P=0.0085,
#'                                                ESTIMATE=1.42272e-02,
#'                                                SE=0.004268,
#'                                                N=25874),
#'                              trait2=data.table(...))), 
#'              coloc_partners=list(ref1=data.table(...), 
#'                                  ref2=data.table(...)),
#'              lst=lst,
#'              len=len,
#'              chr=chr)}.\cr
#'    \code{coloc_data} is a list with one data.table for each retrieved phenotype. Names of the list are 
#'      the names of the phenotypes. Each data.table requires columns \code{}\cr
#'    \code{coloc_partners} List with the same format as \code{coloc_data}, holding the 'reference' or 
#'      partner summary stats for the colocalization test. Here, we compared all phenotypes to the 
#'      three blood pressure related traits \code{diastolic, systolic and mean arterial pressure} (dbp, sbp, map)\cr
#'    \code{lst,len,chr} are start, end and chromosome of genomic region.
#'
#' @return Tibble with columns \code{traits, H0, H1, H2, 'H1,2', H12}
coloc_loop <- function(input_data) {
  # i <- 1
  coloc_result <- list()
  for(i in 1:length(input_data$coloc_data)) {
    mytrx <- names(input_data$coloc_data)[i]
    ukbb <- input_data$coloc_data[[mytrx]]
    print(sprintf("%s/%s: %s", i, length(input_data$coloc_data), mytrx))
    
    # skip over these traits, are the 'coloc partners'
    if(mytrx %in% c("dbp", "sbp", "map", "dbp_adj", "sbp_adj", "map_adj")) {
      next
    }
    
    for(j in 1:length(input_data$coloc_partners)) {
      myprtx <- names(input_data$coloc_partners)[j]
      cref <- input_data$coloc_partners[[myprtx]]
      print(sprintf("... vs %s", myprtx))
      coloc_result[[sprintf("%s_%s",mytrx,myprtx)]] <- make_coloc(ukbb, cref, pstart=input_data$lst, pend=input_data$len)$result
    }
  }
  
  coloc_out <- tibble::as_tibble(cbind(traits=names(coloc_result), do.call("rbind", coloc_result)))
  return(coloc_out)
}






#### eQTL functions ####

#' Get eqtl data for gene from ensembl
#'
#' @param gene String. Gene symbol
#'
#' @return eqtl data
get_eqtl_ensembl <- function(gene) {
  # map to ensemble id
  stopifnot(require("org.Hs.eg.db")) # remember to install it if you don't have it already
  
  gene_ens <- mapIds(org.Hs.eg.db, keys = gene, keytype = "SYMBOL", column="ENSEMBL")
  
  # get all eqtl pvalues for gene
  server <- "https://rest.ensembl.org"
  ext <- sprintf("/eqtl/id/homo_sapiens/%s?statistic=p-value",gene_ens)
  r <- httr::GET(paste(server, ext, sep = ""), content_type("application/json"))
  httr::stop_for_status(r)
  
  # Convert to data table
  tabj <- jsonlite::fromJSON(jsonlite::toJSON(httr::content(r)), flatten = TRUE)
  tab <- data.table::as.data.table(do.call("cbind", lapply(tabj, unlist)))
  
  # convert some columns to numeric
  mycol <- c("value","seq_region_start", "seq_region_name", "seq_region_end", "minus_log10_p_value")
  for(mc in mycol) {
    tab[[mc]] <- as.numeric(tab[[mc]])
  }
  return(tab)
}

#' Make an eqtl plot for a gene with flanking regions
#'
#' @param tab eqtl data
#' @param gene_info gene info data frame
#' @param flank_left integer
#' @param flank_right integer
#' @param mark_region gene info data frame to be highlighted in plot
#' @param mark_snps rs ids to be highlighted in the plot
#' @param mark_snps_nudge_x x shift of marked labels
#' @param mark_snps_size Numeric. Size of snp marks
#' @param sig_snps_size Numeric. Size of snp marks
#' @param sig_snps_p_cut Positive numeric. -log10 pvalue cutoff for sig snps
#' @param debug boolean. If TRUE, enter debug browser
#'
#' @return ggplot2 object. 
plot_eqtl_flanking_region <- function(tab, gene_info, 
                                      flank_left=200000, # 200kb left/right of gene
                                      flank_right=flank_left,
                                      mark_region=NULL,
                                      mark_snps=NULL,
                                      mark_snps_nudge_x=100000, # x shift of marked labels
                                      mark_snps_size=4,
                                      sig_snps_size=2,
                                      sig_snps_p_cut=5,
                                      debug=FALSE) {
  if(debug) {
    browser()
  }
  require(ggplot2)
  require(grid)
  if(missing(gene_info)) {
    stop("Require information on the gene. Please add gene_info argument!")
  }
  stopifnot("data.table" %in% class(gene_info))
  stopifnot(all(c("start", "end", "external_name") %in% names(gene_info)))

  # get gene info
  gene_start <- as.numeric(gene_info$start)
  gene_end <- as.numeric(gene_info$end)
  gene_name <- gene_info$external_name
  gene_y <- 1
  
  # set gene position
  sstart <- gene_start - flank_left
  send <- gene_end + flank_right
  tab_reg <- tab[seq_region_start>sstart & seq_region_end<send]
  
  # plot
  p <- ggplot(data=tab_reg, aes(x=seq_region_start, y=minus_log10_p_value), group=tissue)
  p <- p + ggplot2::geom_point(color="lightgray")
  p <- p + ggplot2::geom_point(data = tab_reg[minus_log10_p_value>sig_snps_p_cut], aes(color=tissue))
  
  if(!is.null(mark_snps)) {
    p <- p + ggrepel::geom_text_repel(data = tab_reg[minus_log10_p_value>sig_snps_p_cut & snp %in% mark_snps], 
                                      aes(label=snp), 
                                      nudge_x           = mark_snps_nudge_x ,
                                      nudge_y = 1,
                                      direction         = "y",
                                      hjust             = 0,
                                      segment.size      = 0.2,
                                      color="red", size=mark_snps_size)
  } else {
    p <- p + ggrepel::geom_text_repel(data = tab_reg[minus_log10_p_value>sig_snps_p_cut], 
                                      aes(color=tissue, label=snp), size=sig_snps_size)
  }
  p <- p + ggtitle(label = sprintf("%s (arrows indicate gene position)", gene_name))
  if(!is.null(gene_info)) {
    # Add segment if given
    p <- p + geom_segment(aes(x = gene_start, y = gene_y, xend = gene_end, yend = gene_y), 
                          arrow=arrow(length = unit(0.2, "cm"), ends = "both"), 
                          color = "black")
    p <- p + ggplot2::geom_vline(xintercept = c(gene_start,gene_end), linetype="dotted", 
                                 color = "black", size=0.7)
    p <- p + ggplot2::annotate("text", label=gene_name, x = gene_start, y = 0.2, size=6, colour="black")
  }
  if(!is.null(mark_region)) {
    # Add segment for additional part if given
    ms <- as.numeric(mark_region$start)
    me <- as.numeric(mark_region$end)
    mn <- mark_region$external_name
    p <- p + geom_segment(aes(x = ms, y = gene_y, xend = me, yend = gene_y), 
                          arrow=arrow(length = unit(0.2, "cm"), ends = "both"), 
                          color = "blue")
    p <- p + ggplot2::geom_vline(xintercept = c(ms,me), linetype="dotted", 
                                 color = "blue", size=0.7)
    p <- p + ggplot2::annotate("text", label=mn, x = ms, y = 0.2, size=6, colour="blue")
  }
  print(p)
}

plot_eqtl_wrap <- function(gene, tab=NULL,
                           flank_left=200000, 
                           flank_right=flank_left,
                           mark_region=NULL,
                           mark_snps=NULL,
                           mark_snps_nudge_x=100000, # x shift of marked labels
                           mark_snps_size=4,
                           sig_snps_size=2,
                           sig_snps_p_cut=5) {
  if(is.null(tab)) {
    tab <- get_eqtl_ensembl(gene=mygene)
  }
  gene_info <- get_gene_ensembl(gene = gene)
  plot_eqtl_flanking_region(tab, flank_left=flank_left,
                            flank_right=flank_right,
                            gene_info=gene_info,
                            mark_region=mark_region,
                            mark_snps=mark_snps,
                            mark_snps_nudge_x=mark_snps_nudge_x, # x shift of marked labels
                            mark_snps_size=mark_snps_size,
                            sig_snps_size=sig_snps_size,
                            sig_snps_p_cut=sig_snps_p_cut)
}


#' Fetch gene info from ensembl
#'
#' @param gene String. Gene symbol
#'
#' @return Gene info data frame.
get_gene_ensembl <- function(gene) {
  stopifnot(require("org.Hs.eg.db")) # remember to install it if you don't have it already
  require(httr)
  require(jsonlite)
  require(xml2)

  # symbol mapping  
  gene_ens <- mapIds(org.Hs.eg.db, keys = gene, keytype = "SYMBOL", column="ENSEMBL")
  
  server <- "https://rest.ensembl.org"
  ext <- sprintf("/overlap/id/%s?feature=gene", gene_ens)
  r <- GET(paste(server, ext, sep = ""), content_type("application/json"))
  stop_for_status(r)
  # Convert to data table
  tabj <- jsonlite::fromJSON(jsonlite::toJSON(httr::content(r)), flatten = TRUE)
  tabx <- data.table::as.data.table(do.call("cbind", lapply(tabj, unlist)))
  # final selection:
  tabx <- tabx[biotype=="protein_coding" & external_name == gene]
  if(nrow(tabx)>1) {
    stop("Oops, something went wrong, more than one row in the gene information.")
  }
  return(tabx)
}
