#' @title Opposing homozygous (OH) loci between all individuals
#'
#' @description Count the number of opposing homozygous (OH) loci between all 
#'  pairs of individuals, i.e. one has 0 copies of the minor allele and the 
#'  other has 2. Parent-offspring pairs are never OH at any locus, except due to
#'  genotyping errors. 
#'
#' @param GenoM numeric matrix with genotype data: One row per individual, and 
#' one column per SNP, coded as 0, 1, 2. Missing value is irrelevant for this 
#' function. 
#'
#' @return a matrix with number of rows and columns equal to the number of 
#'  individuals in \code{GenoM}, with values between 0 and the number of SNPs.
#'
#' @examples
#' Geno_HS <- sequoia::SimGeno(Ped=Ped_HSg5, nSnp=200, ParMis=0.2, SnpError=5e-3)
#' OHM <- CalcOHM(Geno_HS)
#' hist(OHM, breaks=c(-1:max(OHM))+.5, col="grey")

CalcOHM <- function(GenoM) {
  CountOH <- function(GA, GB)  sum((GA==2 & GB==0) | (GA==0 & GB==2), na.rm=TRUE)
  
  OHM <- matrix(NA, nrow(GenoM), nrow(GenoM))
  for (i in 1:nrow(GenoM)) {
    if (i %% 50 == 0) cat(i, "\t")
    for (j in 1:nrow(GenoM)) {
      OHM[i,j] <- CountOH(GenoM[i,], GenoM[j,])
    }
  }
  return(OHM)
}





################################################################################
################################################################################
################################################################################

#' @title Shared non-missing loci between all individuals
#'
#' @description Count the number of loci that are non-missing for both 
#' individuals, for all pairs of individuals
#'
#' @param GenoM numeric matrix with genotype data: One row per individual, and 
#'  one column per SNP, coded as 0, 1, 2, and \code{MissingValue}.
#' @param MissingValue the number or letter denoting missing values. May be a 
#'  vector with several values. 
#'
#' @return a matrix with number of rows and columns equal to the number of 
#'  individuals in \code{GenoM}, with values between 0 and the number of SNPs.
#'
#' @examples
#' Geno_HS <- sequoia::SimGeno(Ped=Ped_HSg5, nSnp=200, CallRate=0.8)
#' NotMisM <- CalcNotMissing(Geno_HS)
#' hist(NotMisM, breaks=c(-1:max(NotMisM))+.5, col="grey")

CalcNotMissing <- function(GenoM, MissingValue=c(NA,-9)) {
  CountNotMis <- function(GA, GB)  sum(!GA %in% MissingValue & !GB %in% MissingValue) 
  
  NotMisM <- matrix(NA, nrow(GenoM), nrow(GenoM))
  for (i in 1:nrow(GenoM)) {
    if (i %% 50 == 0) cat(i, "\t")
    for (j in 1:nrow(GenoM)) {
      NotMisM[i,j] <- CountNotMis(GenoM[i,], GenoM[j,])
    }
  }
  return(NotMisM)
}
