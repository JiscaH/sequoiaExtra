#' @title compare output of pedigree_checker.f90 with known pedigree
#'
#' @description Useful to e.g. test accuracy for a specific pedigree, SNP
#'   number, etc. This function calls pedChecker_topRel()
#'
#' @param FileName  output file from pedigree_checker.f90, with or without
#'  --noFS
#' @param Pedigree reference pedigree with columns id-dam-sire. Passed to
#'   sequoia::GetRelM()
#' @param Threshold Number between 0 and 1, passed to pedChecker_topRel()
#'
#' @return an array with 3 dimensions: 
#'  \item parent1 vs parent2 (i.e. 2nd and 3rd column in FileName) 
#'  \item pedigree relationship: PO,FS,HS,GP,FA,HA,FC1,U 
#'  \item most probable relationship according to pedigree_checker.f90: PO, FS, 
#'  GP(=HS=FA), HA(=3rd degree), UU
#'  
#' @examples
#' Pedigree <- sequoia::Ped_griffin
#' N <- 5000
#' L <- 300

#' GenoM <- sequoia::SimGeno(Pedigree, nSnp=L, ParMis=0)

#' trios <- data.frame(id1 = sample(Pedigree$id, size=N, replace=TRUE),
#'                     id2 = sample(Pedigree$id, size=N, replace=TRUE),
#'                     id3 = sample(Pedigree$id, size=N, replace=TRUE))
#' # exclude any rows with same ID in more than 1 column
#' all_different <- apply(trios, 1, function(x) length(unique(x))==3)
#' trios <- trios[all_different, ]
#' 
#' write.table(GenoM, 'griffin/Geno.txt', row.names=TRUE, col.names=FALSE, quote=FALSE)
#' write.table(trios, 'griffin/trios.txt', row.names=FALSE, col.names=TRUE, quote=FALSE)
#' 
#' # in linux terminal / cygwin:
#' # ./../PedChecker --trios trios.txt
#' 
#' tbl <- pedChecker_compare('griffin/Pedigree_OUT.txt',
#'                           Pedigree=sequoia::Ped_griffin, Threshold=0.8)
#' 

pedChecker_compare <- function(FileName, Pedigree, Threshold=0.5) 
{
  
  # output from pedigree_checker.f90
  topRels <- pedChecker_topRel(FileName, Threshold)

  # Pedigree
  RelM <- sequoia::GetRelM(Pedigree, GenBack=2, patmat=FALSE, Return='Matrix')
  # combine GP+GO, MP+O, etc.
  RelRename <- c('MP' = 'PO',
                 'O' = 'PO',
                 'GO' = 'GP',
                 'FN' = 'FA',
                 'HN' = 'HA')
  tmp <- ifelse(RelM %in% names(RelRename), RelRename[RelM], RelM)
  RelM <- matrix(tmp, nrow(RelM), dimnames = dimnames(RelM))
  
  # Combine  
  par_cn <- colnames(topRels)[2:3]
  pedrel_cn <- paste0(par_cn,'_pedrel')
  toprel_cn <- paste0(par_cn,'_TopRel')
  for (x in 1:2) {
    topRels[,pedrel_cn[x]] <- RelM[ cbind(topRels[,1], topRels[,x+1]) ]
    topRels[,pedrel_cn[x]] <- factor(topRels[,pedrel_cn[x]], 
                                        levels = c("PO", "FS", "HS", 'GP', 'FA', 'HA', 'FC1', 'U'))
  }
 
  out <- list(table(topRels[,pedrel_cn[1]], topRels[,toprel_cn[1]]),
              table(topRels[,pedrel_cn[2]], topRels[,toprel_cn[2]]))
  
  out <- plyr::laply(out, function(x) x)
  dimnames(out)[[1]] <- par_cn
  
  return( out )
}