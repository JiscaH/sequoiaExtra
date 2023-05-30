#' @title Find most-likely relationship based on output of pedigree_checker.f90
#' 
#' @description For each 'parent1' sum over all relationships with 'parent2', and 
#' likewise for 'parent2'. 
#' 
#' @param FileName  output file from pedigree_checker.f90, with or without 
#'  --noFS
#' @param Threshold Number between 0 and 1. If the summed probability of the
#'   most-likely relationship is below this threshold, 'TopRel' will be changed
#'   to 'XX' and the associated probability to <NA>.
#'   
#' @return a dataframe with
#' \item the 3 ID columns in the file 
#' \item{<parent1>_TopRel} highest probability relationship for parent1: PO,
#' (FS), GP (=any 2nd degree relationship), HA (=any 3rd degree relationship),
#' UU
#' \item{<parent1>_TopRel_prob} summed probability; e.g. if TopRel='PO' this 
#'  is prob_PO_PO + prob_PO_GP + prob_PO_HA + prob_PO_UU (+prob_PO_FS). 
#'  

pedChecker_topRel <- function(FileName, Threshold=0.5) 
{
  trio_probs <- read.table(FileName, header=TRUE, stringsAsFactors=FALSE,
                           na.strings=c('NA', '-9', '999.0000', 'NaN'))
  
  par_colnames <- colnames(trio_probs)[2:3]
  prob_colnames <- grep('^prob_', colnames(trio_probs), value=TRUE)
  if (any(grepl('FS', prob_colnames))) {
    RelNames <- c('PO', 'FS', 'GP', 'HA', 'UU')
  } else {
    RelNames <- c('PO', 'GP', 'HA', 'UU') 
  }
  
  probA <- plyr::aaply(as.matrix(trio_probs[,prob_colnames]), .margins=1,
                       .fun = function(V) matrix(V, length(RelNames), length(RelNames)))
  dimnames(probA) <- list(1:nrow(trio_probs), rel_id2 = RelNames, rel_id2 = RelNames)
  
  get_toprel <- function(M, d) {
    probs <- apply(M,d,sum)
    data.frame(TopRel = names(which.max(probs)),
               TopRel_prob = max(probs))
  }
  
  toprel <- list(par1 = plyr::adply(probA, 1, get_toprel, 1), 
                 par2 = plyr::adply(probA, 1, get_toprel, 2))
  
  for (x in 1:2) {
    # apply threshold
    TooLow <- toprel[[x]]$TopRel_prob < Threshold
    toprel[[x]]$TopRel[TooLow] <- 'XX'
    toprel[[x]]$TopRel_prob[TooLow] <- NA
  }
  
  # rename relationships
  # TODO? GP -> '2nd' etc, but than valid name not starting with number
  
  # turn into factor (fixes table order)
  if (any(c(toprel[[1]]$TopRel, toprel[[2]]$TopRel) == 'XX')) {
    for (x in 1:2) {
      toprel[[x]]$TopRel <- factor(toprel[[x]]$TopRel, levels = c(RelNames, 'XX'))
    }
  }
  
  # rename columns
  for (x in 1:2) {
    colnames(toprel[[x]]) <- paste0(par_colnames[x],'_',colnames(toprel[[x]]))
  }
  
  # combine into dataframe
  out <- cbind(trio_probs[,1:3],
               toprel[[1]][,-1],
               toprel[[2]][,-1])
  return( out )
}