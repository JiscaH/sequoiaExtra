#' @title Assignment confidence by offspring number
#'
#' @description Use output from \code{EstConf} to calculate the assignment 
#'  confidence split by dummy parent's offspring number (category \code{GD}) or
#'  focal dummy's offspring number (\code{DG} and \code{DD}). 
#'
#' @param Conf Output from \code{EstConf}
#'
#' @return an array with 4 dimensions:
#'   \item{1}{'conf' for proportion of assignments that is correct, and 'N' for 
#'   total number of individuals in that category with that number of offspring}
#'   \item{parent}{'dam' vs 'sire'}
#'   \item{cat}{GD, DG, or DD, where G=Genotyped, D=Dummy, first letter = focal 
#'   individual, second letter=assigned parent}
#'   \item{nOff}{Offspring number, of the parent (cat GD) or of the focal 
#'   individual (cat DG and DD). cat DD ignores parent offspring number for 
#'   simplicity (for now). ONLY GENOTYPED OFFSPRING ARE COUNTED.}
#'   
#' @details This function is a beta version with considerable scope for 
#'  improvement and more nuance. 
#'
#' @examples
#' conf_grif <- EstConf(Pedigree = Ped_griffin,
#'                      LifeHistData = LH_griffin,
#'                      args.sim = list(nSnp = 200, SnpError = 0.01, ParMis=c(0.2, 0.5)),
#'                      args.seq = list(Module = 'ped', Err=0.01),
#'                      nSim = 2)
#'                      
#' conf_grif_nOff <- Conf_by_nOff(conf_grif)
#' conf_grif_nOff['conf',,'GD',]
#' conf_grif_nOff['N',,'GD',]

Conf_by_nOff <- function(Conf) {
  
  nSim <- Conf$RunParams$EstConf$nSim
  
  tbl_nOff <- with(Conf$Pedigree.reference, c(table(dam), table(sire))) 
  mxOff <- max(tbl_nOff) +2
  
  ClassNames <- c("Match", "Mismatch", "P1only", "P2only", "_")
  counts_by_nOff <- array(dim = c(nSim, 2, 3, mxOff, 5),
                          dimnames = list(iter = seq_len(nSim),
                                          parent = c('dam', 'sire'),
                                          cat = c('GD',  # dam's/sire's no. offspring
                                                  'DG', 'DD'),  # id's no. offspring
                                          nOff = seq_len(mxOff),
                                          class = ClassNames))
  for (i in 1:nSim) {
    PedComp_i <- sequoia::PedCompare(Ped1 = Conf$Pedigree.inferred[[i]],
                            Ped2 = Conf$Pedigree.reference,
                            SNPd = Conf$SimSNPd[[i]],
                            Symmetrical=FALSE, Plot=FALSE)
    PedM <- PedComp_i$MergedPed
    # count genotyped offspring only (else mismatch w getAssignable)
    PedG <- PedM[PedM$id %in% Conf$SimSNPd[[i]], ]
    
    PedM$nOff.dam <- table(PedG$dam.1)[PedM$dam.1]
    PedM$nOff.sire <- table(PedG$sire.1)[PedM$sire.1]
    PedM$nOff.id <- c(table(PedG$dam.1), table(PedG$sire.1))[PedM$id.r]
    PedM$nOff.id[is.na(PedM$nOff.id)] <- 0
    
    for (p in c('dam', 'sire')) {
      counts_by_nOff[i,p,'GD', , ] <- table(factor(PedM[, paste0('id.',p,'.cat')] == 'GD', levels=c(TRUE,FALSE)),
                                            factor(PedM[, paste0('nOff.',p)], levels=1:mxOff),
                                            factor(PedM[, paste0(p,'.class')], levels=ClassNames))['TRUE',,]
      for (a in c('DG', 'DD')) {
        counts_by_nOff[i,p,a, , ] <- table(factor(PedM[, paste0('id.',p,'.cat')] == a, levels=c(TRUE,FALSE)),
                                           factor(PedM[, 'nOff.id'], levels=1:mxOff),
                                           factor(PedM[, paste0(p,'.class')], levels=ClassNames))['TRUE',,]
      }
    }
  }
  
  conf_by_nOff <- apply(counts_by_nOff, MARGIN=c('parent', 'cat', 'nOff'),
                        function(M) c(conf = sum(M[,'Match'])/sum(M), N = sum(M)))
  
  return(conf_by_nOff)
}
