#' @title Find pedigree loop after 'an individual is its own ancestor' error
#'
#' @description When \code{sequoia} function \code{PedPolish} or
#'  \code{getGenerations} warns that an individual is its own ancestor, this 
#'  function identifies the pedigree loop(s). 
#'
#' @param Pedigree dataframe with columns id - dam - sire. 
#'
#' @return a list with one or more sub-pedigrees of the pedigree loops. The
#'   individual on the first row of (each/the) pedigree loop(s) is arbitrary. On
#'   each subsequent row is the id either parent of the id on the previous row,
#'   until a parent is the id on the first row and the loop is closed.
#'
#' @examples
#' Ped_grif_X <- Ped_griffin
#' Ped_grif_X$sire[Ped_grif_X$id == 'i048_2003_F'] <- 'i142_2008_M'
#' Find_ped_loops(Ped_grif_X)
#' 
 
Find_ped_loops <- function(Pedigree) {
  
  G <- suppressWarnings( sequoia::getGenerations(Pedigree, StopIfInvalid=FALSE) )
  if (!any(is.na(G)))  stop('Pedigree has no individual-is-own-ancestor loops')
  
  IsOwnAnc <- NA
  for (i in 1:nrow(Pedigree)) {
    if (!is.na(G[i])) {
      IsOwnAnc[i] <- FALSE
    } else {
      Anc_i <- sequoia:::GetAnc(i, Pedigree)
      if (any(unlist(Anc_i) == as.character(Pedigree$id)[i])) {
        IsOwnAnc[i] <- TRUE
      } else {
        IsOwnAnc[i] <- FALSE
      }
    }
  }
  
  Ped_loop_unsorted <- as.matrix(Pedigree[IsOwnAnc, 1:3])
  rownames(Ped_loop_unsorted) <- NULL
  
  Ped_loops <- list()
  for (l in 1:sum(IsOwnAnc))  {  # number of separate loops
    
    # sort: keep 1st row, 2nd row = its parents, etc.
    Ped_loop_sorted <- Ped_loop_unsorted  # same dimensions
    Ped_loop_sorted[] <- NA
    Ped_loop_sorted[1,] <- Ped_loop_unsorted[1,]
    Ped_loop_unsorted[1,] <- NA  # done. 
    
    y <- 1  
    for (x in 1:nrow(Ped_loop_unsorted)) {
      for (p in c('dam', 'sire')) {
        if (is.na(Ped_loop_sorted[x,p]))  next
        z <- match(Ped_loop_sorted[x,p], Ped_loop_unsorted[,'id'], nomatch=0)
        if (z > 0) {
          y <- y+1
          Ped_loop_sorted[y,] <- Ped_loop_unsorted[z,]
          Ped_loop_unsorted[z,] <- NA  # done
        }
      }
    }
    Ped_loops[[l]] <- as.data.frame(Ped_loop_sorted[!is.na(Ped_loop_sorted[,'id']), ])
    
    if (all(is.na(Ped_loop_unsorted[,'id']))) break
    Ped_loop_unsorted <- Ped_loop_unsorted[!is.na(Ped_loop_unsorted[,'id']), ]
  }
  
  return(Ped_loops)
  
}
