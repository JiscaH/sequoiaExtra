# This function is adapted from the R script at https://yanglab.westlake.edu.cn/software/gcta/#MakingaGRM
# by Jisca Huisman, 2023-05-19. jisca.huisman@gmail.com .
# You are free to re-use, edit, and share.

ReadGRMBin <- function(prefix, AllN=FALSE, size=4, Return='list')
{
  BinFileName <- paste0(prefix, ".grm.bin")
  NFileName <- paste0(prefix, ".grm.N.bin")
  IDFileName <- paste0(prefix,".grm.id")
  id <- read.table(IDFileName)
  n <- dim(id)[1]
  BinFile <- file(BinFileName, "rb");
  grm <- readBin(BinFile, n=n*(n+1)/2, what=numeric(0), size=size)
  NFile <- file(NFileName, "rb");
  if(AllN) {
    N <- readBin(NFile, n=n*(n+1)/2, what=numeric(0), size=size)
  } else {
    N <- readBin(NFile, n=1, what=numeric(0), size=size)
  }

  if (Return=='list') {
    i <- sapply(1:n, function(i) sum(1:i))
    return(list(diag=grm[i], off=grm[-i], id=id, N=N))

  } else {
    grm.M <- matrix(NA, n, n, dimnames=list(id[,2],id[,2]))
    grm.M[upper.tri(grm.M, diag=TRUE)] <- grm
    if (Return=='matrix') {
      grm.M[lower.tri(grm.M)] <- t(grm.M)[lower.tri(grm.M)]
      return(grm.M)

    } else if (Return == 'dataframe') {
      grm.df <- plyr::adply(grm.M, .margins=c(1,2))
      colnames(grm.df) <- c('IID1', 'IID2', 'R.grm')
      return(grm.df[!is.na(grm.df$R.grm), ])
    }
  }
}
