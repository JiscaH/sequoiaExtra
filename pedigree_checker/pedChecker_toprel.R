#' @title Find most-likely relationship based on output of pedigree_checker.f90
#'
#' @description For each 'parent1' sum over all relationships with 'parent2',
#'   and likewise for 'parent2'.
#'
#' @param FileName  output file from pedigree_checker.f90, with or without
#'  --noFS
#' @param Threshold Number between 0 and 1 (inclusive). If the (summed)
#'   probability of the most-likely relationship is below this threshold,
#'   'TopRel' will be changed to 'XX' and the associated probability to <NA>.
#'
#' @return a dataframe with
#' \item the 3 ID columns in the file
#' \item{pair_TopRel} relationship combination with highest probability. For
#'   non-genotyped parents, 'ZZ' is used instead of 'UU'.
#' \item{pair_TopRel_prob} associated probability
#' \item{<parent1>_TopRel} relationship with highest probability for parent1:
#'   PO, (FS), GP (=any 2nd degree relationship), HA (=any 3rd degree
#'   relationship), UU. For non-genotyped parents, 'ZZ' is used instead of 'UU'.
#' \item{<parent1>_TopRel_prob} summed probability; e.g. if TopRel='PO' this
#'  is prob_PO_PO + prob_PO_GP + prob_PO_HA + prob_PO_UU (+prob_PO_FS).
#' \item{<parent2>_TopRel}
#' \item{<parent2>_TopRel_prob}
#'
#'
#' @details
#' Abbreviations used: PO = parent-offspring, FS = full siblings,
#' GP = grand-parental, HA = half avuncular (niece/nephew -- aunt/uncle),
#' UU = unrelated.
#' Grand-parental relationships are indistinguishable from the other second
#' degree relationships, half-siblings and full avuncular. Similarly,
#' half-avuncular is indistinguishable from other third degree relationships
#'
#' In the output column names, <parent1> and <parent2> are the second and third
#' column name in the input file.



pedChecker_topRel <- function(FileName, Threshold=0.5)
{
  trio_probs <- read.table(FileName, header=TRUE, stringsAsFactors=FALSE,
                           na.strings=c('NA', '-9', '999.0000', 'NaN'))
  trio_IDs <- trio_probs[,1:3]

  prob_colnames <- grep('^prob_', colnames(trio_probs), value=TRUE)
  if (any(grepl('FS', prob_colnames))) {
    RelNames <- c('PO', 'FS', 'GP', 'HA', 'UU')
  } else {
    RelNames <- c('PO', 'GP', 'HA', 'UU')
  }
  trio_probs <- as.matrix(trio_probs[,prob_colnames])

  # most likely for each 'offspring'-'parent' pair separately
  probA <- plyr::aaply(trio_probs, .margins=1,
                       .fun = function(V) matrix(V, length(RelNames), length(RelNames)))
  dimnames(probA) <- list(1:nrow(trio_probs), rel_id2 = RelNames, rel_id2 = RelNames)

  get_toprel <- function(M, d) {
    probs <- apply(M, d,
                   function(v) ifelse(all(is.na(v)), NA, sum(v, na.rm=TRUE)))
    if (any(is.na(probs))) {  # parent not genotyped
      data.frame(TopRel = 'ZZ', TopRel_prob = NA)
    } else {
      data.frame(TopRel = names(which.max(probs)),
               TopRel_prob = max(probs, na.rm=TRUE))
    }
  }

  # most likely for (parent-)pair
  trio_probs <- cbind(trio_probs, 'ZZ_ZZ' = NA)
  prob_colnames <- c(gsub('prob_', '', prob_colnames), 'ZZ_ZZ')
  toprel_pair_w <- apply(trio_probs, 1,
                       function(x) ifelse(any(!is.na(x)), which.max(x), length(prob_colnames)))
  toprel_pair <- data.frame(TopRel = prob_colnames[toprel_pair_w],
                            TopRel_prob = trio_probs[cbind(1:nrow(trio_probs), toprel_pair_w)])

  # combine
  toprel <- list(pair = toprel_pair,
                 par1 = plyr::adply(probA, 1, get_toprel, 1),
                 par2 = plyr::adply(probA, 1, get_toprel, 2))

  # pair UU --> ZZ if parent not genotyped
  toprel$pair$TopRel <- ifelse(toprel$par1$TopRel=='ZZ',
                               gsub('^UU_', 'ZZ_', toprel$pair$TopRel),
                               toprel$pair$TopRel)
  toprel$pair$TopRel <- ifelse(toprel$par2$TopRel=='ZZ',
                               gsub('_UU$', '_ZZ', toprel$pair$TopRel),
                               toprel$pair$TopRel)


  for (x in names(toprel)) {
    # apply threshold
    TooLow <- toprel[[x]]$TopRel_prob < Threshold
    toprel[[x]]$TopRel[TooLow] <- 'XX'
    toprel[[x]]$TopRel_prob[TooLow] <- NA
  }


  # turn into factor (fixes table order)
  RelNames <- union(RelNames,
                    unique(c(toprel[[2]]$TopRel, toprel[[3]]$TopRel)))  # add XX and/or ZZ
  for (x in c('par1', 'par2')) {
    toprel[[x]]$TopRel <- factor(toprel[[x]]$TopRel, levels = RelNames)
  }
  pair_RelNames <- intersect(c(t(outer(c(RelNames,'ZZ'), c(RelNames,'ZZ'), paste, sep='_')),'XX'),
                             unique(toprel[[1]]$TopRel))
  toprel[['pair']]$TopRel <- factor(toprel[['pair']]$TopRel, levels=pair_RelNames)


  # rename columns
  par_colnames <- c(pair='pair', par1=colnames(trio_IDs)[2], par2=colnames(trio_IDs)[3])
  for (x in names(toprel)) {
    colnames(toprel[[x]]) <- paste0(par_colnames[x],'_',colnames(toprel[[x]]))
  }

  # combine into dataframe
  out <- cbind(trio_IDs,
               toprel[['pair']],
               toprel[['par1']][,-1],
               toprel[['par2']][,-1])
  return( out )
}
