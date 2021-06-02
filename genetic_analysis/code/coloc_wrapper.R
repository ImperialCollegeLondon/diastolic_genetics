############ colocalisation test ############
# coloc_input is a coloc_input frame with beta1, se1, beta2, se2
#' @title Colocalisation test based on https://github.com/tobyjohnson/gtx/blob/master/R/coloc.R
#' @description A test for the support of overlapping GWAS signals. The result is the posterior probabilities for H0 (no signals), H1 (a GWAS signal for the first sample only), H2 (a GWAS signal for the second sample only), H1,2 (a GWAS signal in both samples but not overlapping), and H12 (a colocalisation in the 2 signals)
#' @return List with posteriors for hypotheses H0, H1, H2, H1,2, H12
#' @param coloc_input a data frame with beta1, se1, beta2, se2
#' @export

coloc_wrapper <- function(coloc_input, priorsd1=1, priorsd2=1, priorc1=1e-4, priorc2=1e-4, priorc12=1e-5, join_type='inner', summary_only=FALSE) {

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
