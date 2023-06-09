#' @title Turn birth date into quarter number, starting Jan-March of Year1
#'
#' @description For species that breed multiple times per year and/or with a age
#' of first reproduction of less than a year, one cannot work with birth year in 
#' \code{sequoia}. The quarter number may be used instead. 
#'
#' @param BirthDate  either a named vector with dates (YY-MM-DD), or a matrix 
#'  with columns YY (year) and MM (month) and with IDs as rownames.
#' @param Year1  starting year
#'
#' @return a named vector with quarter numbers since Jan-March of Year1, starting counting at 1.
#'
#' @examples
#' BirthDate <- c('Yellow' = '2004-03-12', 'Green' = '2006-12-24', 'Red'='2005-01-08')
#' BirthDate2Quarter(BirthDate)  

BirthDate2Quarter <- function(BirthDate,
                                  Year1 = 2000) 
{
  if (inherits(BirthDate, 'matrix')) {
    IDs <- rownames(BirthDate)
    B.Year <- BirthDate[,1] - Year1
    B.Month <- BirthDate[,2]
  } else {
    BirthDate <- as.Date(BirthDate)
    IDs <- names(BirthDate)
    B.Year <- as.numeric(format(BirthDate, '%Y')) - Year1
    B.Month <- as.numeric(format(BirthDate, '%m'))
  }
  if (any(!is.na(B.Year) & B.Year < 0))  stop('Some BirthDates are before Year1')
  
  Quart.num <- 4*B.Year + ceiling(B.Month/3)  # round up
  
  names(Quart.num) <- IDs
  return(Quart.num)
}





################################################################################
################################################################################
################################################################################

#' @title Reverse \code{BirthDate2Quarter}
#'
#' @description Turn quarter numbers as generated by \code{BirthDate2Quarter}
#' back into calendar year + quarter (1-4). The same \code{Year1} must
#' be used!
#'
#' @param QuarterNumber a named vector with quarter numbers since Year1. 
#' @param Year1  starting year 
#'
#' @return a matrix with columns YY and QQ, with IDs in rownames. 
#' 

Quarter2BirthDate <- function(QuarterNumber,
                              Year1 = 2000) 
{
  B.Year <- floor( (QuarterNumber-1)/4) + Year1   
  B.Quarter <- (QuarterNumber-1) %% 4   
  
  BirthDate <- cbind(YY = B.Year,
                     Quarter = B.Quarter)
  rownames(BirthDate) <- names(QuarterNumber)
  
  return(BirthDate)
}
