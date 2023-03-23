#============================================================================
#============================================================================
#' @title Confidence Probabilities, with some experimental options
#'
#' @description Estimate confidence probabilities ('backward') and assignment
#'   error rates ('forward') per category (genotyped/dummy) by repeatedly
#'   simulating genotype data from a reference pedigree using
#'   \code{\link{SimGeno}}, reconstruction a pedigree from this using
#'   \code{\link{sequoia}}, and counting the number of mismatches using
#'   \code{\link{PedCompare}}.
#'
#' @details The confidence probability is taken as the number of correct
#'   (matching) assignments, divided by all assignments made in the
#'   \emph{observed} (inferred-from-simulated) pedigree. In contrast, the false
#'   negative & false positive assignment rates are proportions of the number of
#'   parents in the \emph{true} (reference) pedigree. Each rate is calculated
#'   separatedly for dams & sires, and separately for each category
#'   (\strong{G}enotyped/\strong{D}ummy(fiable)/\strong{X} (none)) of
#'   individual, parent and co-parent.
#'
#'  This function does not know which individuals in the actual \code{Pedigree}
#'  are genotyped, so the confidence probabilities need to be added to the
#'  \code{Pedigree} as shown in the example at the bottom.
#'
#'  A confidence of \eqn{1} means all assignments on simulated data were correct for
#'  that category-combination. It should be interpreted as (and perhaps modified
#'  to) \eqn{> 1 - 1/N}, where sample size \code{N} is given in the last column
#'  of the \code{ConfProb} and \code{PedErrors} dataframes in the output. The
#'  same applies for a false negative/positive rate of \eqn{0} (i.e. to be
#'  interpreted as \eqn{< 1/N}).
#'
#' @section Assumptions:
#'   Because the actual true pedigree is (typically) unknown, the provided
#'   reference pedigree is used as a stand-in and assumed to be the true
#'   pedigree, with unrelated founders. It is also assumed that the probability
#'   to be genotyped is equal for all parents; in each iteration, a new random
#'   set of parents (proportion set by \code{ParMis}) is mimicked to be
#'   non-genotyped. In addition, SNPs are assumed to segregate independently.
#'
#' @section Object size:
#'   The size in Kb of the returned list can become pretty big, as each of the
#'   inferred pedigrees is included. When running \code{EstConf} many times for
#'   a range of parameter values, it may be prudent to save the required summary
#'   statistics for each run rather than the full output.
#'
#'
#' @param Pedigree reference pedigree from which to simulate, dataframe with
#'   columns id-dam-sire. Additional columns are ignored.
#' @param LifeHistData dataframe with id, sex (1=female, 2=male, 3=unknown),
#' birth year, and optionally BY.min - BY.max - YearLast.
#' @param args.sim  list of arguments to pass to \code{\link{SimGeno}}, such as
#'   \code{nSnp} (number of SNPs), \code{SnpError} (genotyping error rate) and
#'   \code{ParMis} (proportion of non-genotyped parents). Set to \code{NULL} to
#'   use all default values.
#' @param args.seq  list of arguments to pass to \code{\link{sequoia}}, such as
#'   \code{Module} ('par' or 'ped'), \code{Err} (assumed genotyping error rate),
#'   and \code{Complex}. May include (part of) \code{SeqList}, a list of sequoia
#'   output (i.e. as a list-within-a-list). Set to \code{NULL} to use all
#'   default values.
#' @param nSim number of iterations of simulate - reconstruct - compare to
#'   perform, i.e. number of simulated datasets.
#' @param nCores number of computer cores to use. If \code{>1}, package
#'   \pkg{parallel} is used. Set to NULL to use all but one of the available
#'   cores, as detected by \code{parallel::detectCores()} (using all cores tends
#'   to freeze up your computer).
#' @param CallRate.sim  matrix with for each individual (rows) and each
#'   iteration (columns) the callrate to simulate. Allows fine-grained control
#'   to e.g. estimate confidence for individuals with low call rate.
#'   (Experimental)
#' @param groupX  a list of length \code{nSim} with vectors of IDs of
#'   individuals only genotyped for the first half of the SNPs, i.e. the rest
#'   is genotyped for an additional set of SNPs. (Experimental)
#' @param quiet suppress messages. \code{TRUE} runs \code{SimGeno} and
#'   \code{sequoia} quietly, \code{'very'} also suppresses other messages and
#'   the iteration counter when \code{nCores=1} (there is no iteration counter
#'   when \code{nCores>1}).
#'
#' @return A list, with elements:
#'   \item{ConfProb}{See below}
#'   \item{PedErrors}{See below}
#'   \item{Pedigree.reference}{the pedigree from which data was simulated}
#'   \item{LifeHistData}{}
#'   \item{Pedigree.inferred}{a list with for each iteration the inferred
#'     pedigree based on the simulated data}
#'   \item{SimSNPd}{a list with for each iteration the IDs of the individuals
#'     simulated to have been genotyped}
#'   \item{PedComp.fwd}{array with \code{Counts} from the 'forward'
#'     \code{PedCompare}, from which \code{PedErrors} is calculated}
#'   \item{RunParams}{a list with the call to \code{EstConf} as a semi-nested
#'   list (args.sim, args.seq, nSim, nCores), as well as the default parameter
#'   values for \code{SimGeno} and \code{sequoia}.}
#'   \item{RunTime}{\code{sequoia} runtime per simulation in seconds, as
#'     measured by \code{\link{system.time}()['elapsed']}.}
#'
#' Dataframe \code{ConfProb} has 7 columns:
#' \item{id.cat, dam.cat, sire.cat}{Category of the focal individual, dam, and
#'   sire, in the pedigree inferred based on the simulated data. Coded as
#'   G=genotyped, D=dummy, X=none}
#' \item{dam.conf}{Probability that the dam is correct, given the categories of
#'   the assigned dam and sire (ignoring whether or not the sire is correct)}
#' \item{sire.conf}{as \code{dam.conf}, for the sire}
#' \item{pair.conf}{Probability that both dam and sire are correct, given their
#'   categories}
#' \item{N}{Number of individuals per category-combination, across all
#'   \code{nSim} iterations}
#'
#' Array \code{PedErrors} has three dimensions:
#' \item{class}{\itemize{
#'   \item \code{FalseNeg}(atives): could have been assigned but was not
#' (individual + parent both genotyped or dummyfiable; P1only in
#' \code{PedCompare}).
#'   \item \code{FalsePos}(itives): no parent in reference pedigree, but
#' one was assigned based on the simulated data (P2only)
#'   \item \code{Mismatch}: different parents between the pedigrees
#'   }}
#' \item{cat}{Category of individual + parent, as a two-letter code where the
#'   first letter indicates the focal individual and the second the parent;
#'   G=Genotyped, D=Dummy, T=Total}
#' \item{parent}{dam or sire}
#'
#'
#' @seealso \code{\link{SimGeno}, \link{sequoia}, \link{PedCompare}}.
#'
#' @importFrom plyr adply
#'
#' @examples
#' NSIM <- 5
#' # 60% of individuals only genotyped for (same) half of SNPs
#' groupX <- list()
#' for (i in 1:NSIM) {
#'   groupX[[i]] <- sample(Ped_griffin$id, size = 0.6*nrow(Ped_griffin))
#' }
#' conf_x <- EstConf2(Pedigree = Ped_griffin,
#'                LifeHistData = LH_griffin,
#'              args.sim = list(nSnp = 2*96, CallRate=0.99, SnpError=1e-3, 
#'                            ParMis=0.2),
#'                args.seq = list(Module="par"),  
#'                groupX = groupX,
#'                nSim = NSIM,
#'                nCores=5)
#'
#' conf_xg <- Conf_by_group(conf_x)
#' dimnames(conf_xg)                

EstConf2 <- function(Pedigree = NULL,
                     LifeHistData = NULL,
                     args.sim = list(nSnp = 400, SnpError = 1e-3, ParMis=c(0.4, 0.4)),
                     args.seq = list(Module="ped", Err=1e-3, Tassign=0.5, CalcLLR = FALSE),
                     nSim = 10,
                     nCores = 1,
                     CallRate.sim = NULL,
                     groupX = NULL,
                     quiet=TRUE)
{
  
  # check input ----
  if (is.null(Pedigree))  stop("Please provide Pedigree")
  if (is.null(LifeHistData))  stop("Please provide LifeHistData")
  if (!is.null(args.sim) & !is.list(args.sim))  stop("args.sim should be a list or NULL")
  if (!is.null(args.seq) & !is.list(args.seq))  stop("args.seq should be a list or NULL")
  if (!sequoia:::is.wholenumber(nSim) || nSim<1 || length(nSim)>1)
    stop("nSim must be a single positive number")
  
  if (!quiet %in% c(TRUE, FALSE, "very"))  stop("'quiet' must be TRUE, FALSE, or 'very'")
  quiet.EC <- ifelse(quiet == "very", TRUE, FALSE)
  quiet <- ifelse(quiet %in% c("very", TRUE), TRUE, FALSE)
  if (!"quiet" %in% names(args.sim))  args.sim <- c(args.sim, list(quiet = quiet))
  if (!"quiet" %in% names(args.seq))  args.seq <- c(args.seq, list(quiet = quiet))
  
  if ("Err" %in% names(args.sim)) {
    args.sim[["SnpError"]] <- args.sim[["Err"]]
    args.sim[["Err"]] <- NULL    # common confusion, otherwise fuzy matching with 'ErrorFM'.
  }
  
  Ped.ref <- sequoia::PedPolish(Pedigree, KeepAllColumns=FALSE)
  if (any(substr(unlist(Ped.ref),1,6) %in% c("sim_F0", "sim_M0"))) {
    stop("Please don't use 'sim_F' or 'sim_M' in reference pedigree")
  }
  
  if (!is.null(CallRate.sim)) {
    if (!inherits(CallRate.sim, 'matrix'))
      stop("'CallRate.sim' must be a matrix; specify single callRate via args.sim")
    if ('CallRate' %in% names(args.sim))
      stop("Please provide 'CallRate' either via 'args.sim' or via 'CallRate.sim', not both")
    if (!all(Ped.ref$id %in% rownames(CallRate.sim)))
      stop("'CallRate.sim' must include a row for every individual in 'Pedigree'")
  }
  
  if (!is.null(groupX)) {
    if (!is.list(groupX) | length(groupX) < nSim)
      stop("'groupX' must be a list of length nSim")
  }
  
  
  if ("Module" %in% names(args.seq)) {
    ParSib <- ifelse(args.seq$Module == "ped", "sib", "par")
  } else if ("MaxSibIter" %in% names(args.seq)) {
    ParSib <- ifelse(args.seq$MaxSibIter > 0, "sib", "par")
  } else {
    ParSib <- "sib"   # default Module = "ped"
  }
  
  if (!quiet.EC) {
    if (ParSib == "par") {
      message("Simulating parentage assignment only ...")
    } else {
      message("Simulating full pedigree reconstruction ...")
    }
  }
  
  
  # no. cores ----
  if (!is.null(nCores) && (!sequoia:::is.wholenumber(nCores) || nCores<1))
    stop("nCores must be a positive number, or NULL")
  
  if (is.null(nCores) || nCores>1) {
    if (!requireNamespace("parallel", quietly = TRUE)) {
      if (interactive() & !quiet.EC) {
        message("Installing pkg 'parallel' to speed things up... ")
      }
      utils::install.packages("parallel")
    }
    maxCores <- parallel::detectCores()
    if (is.null(nCores)) {
      nCores <- maxCores -1
    } else if (nCores > maxCores) {
      nCores <- maxCores
      warning("Reducing 'nCores' to ", maxCores, ", as that's all you have",
              immediate.=TRUE)
    }
    if (nCores > nSim) {
      nCores <- nSim
    }
    if (!quiet.EC)  message("Using ", nCores, " out of ", maxCores, " cores")
  }
  
  
  utils::flush.console()
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # function to simulate genotypes & infer pedigree ----
  
  #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  SimInfer <- function(i, RefPedigree, args.sim,
                       CallRate.sim, groupX,
                       LifeHistData, args.seq,
                       quiet.EC, ParSib)
  {
    if (!quiet.EC)  cat("i=", i, "\t", format(Sys.time(), "%H:%M:%S"), "\n")
    # nothing printed by parallel::parSapply (?)
    OUT.i <- list()
    if (!is.null(CallRate.sim)) {
      args.sim <- c(args.sim, list(CallRate = CallRate.sim[,i]))
    }
    GM <- do.call(sequoia::SimGeno, c(list(Pedigree=RefPedigree), args.sim))
    OUT.i$SimSNPd <- rownames(GM)
    if (!is.null(groupX)) {
      GM[intersect(groupX[[i]], rownames(GM)), 1:floor(ncol(GM)/2)] <- -9
    }
    
    OUT.i$RunTime <- system.time(Seq.i <- do.call(sequoia::sequoia,
                                                  c(list(GenoM = GM,
                                                         LifeHistData = LifeHistData,
                                                         DummyPrefix = c("sim_F", "sim_M"),
                                                         Plot = FALSE),
                                                    args.seq) ))["elapsed"]
    rm(GM)
    gc()
    if (ParSib == "par") {
      OUT.i$Pedigree.inferred <- Seq.i[["PedigreePar"]]
    } else {
      OUT.i$Pedigree.inferred <- Seq.i[["Pedigree"]]
    }
    return( OUT.i )
  }
  #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  
  # call above function, for 1 or >1 cores ----
  if (nCores>1) {
    cl <- parallel::makeCluster(nCores)
    AllOUT <- parallel::parLapply(cl, X=seq.int(nSim), fun=SimInfer,
                                  RefPedigree = Ped.ref,
                                  args.sim, CallRate.sim, groupX,
                                  LifeHistData, args.seq,
                                  quiet.EC=TRUE, ParSib)
    #                                  chunk.size=1)
    parallel::stopCluster(cl)
  } else {
    AllOUT <- plyr::llply(seq.int(nSim), .fun=SimInfer,
                          RefPedigree = Ped.ref,
                          args.sim, CallRate.sim, groupX,
                          LifeHistData, args.seq,
                          quiet.EC, ParSib)
  }
  RunTime <- sapply(AllOUT, "[[", "RunTime")
  Pedigree.inferred <- sapply(AllOUT, "[", "Pedigree.inferred")
  SimSNPd <- sapply(AllOUT, "[", "SimSNPd")
  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # confidence probabilities ----
  nSimz <- ifelse(nSim>1, nSim,2)  # else problems w R auto-dropping dimension
  CatNames <- c("G", "D", "X")
  ClassNames <- c("Match", "Mismatch", "P1only", "P2only", "_")
  
  PC.rev.cd <- array(0, dim = c(nSimz, 3,3,3, 5,5),
                     dimnames = list(iter = seq_len(nSimz),
                                     id.cat = CatNames, dam.cat = CatNames, sire.cat = CatNames,
                                     dam.class = ClassNames, sire.class = ClassNames))
  for (i in 1:nSim) {
    PC.rev.cd[i,,,,,] <- sequoia::PedCompare(Ped1 = Pedigree.inferred[[i]],
                                             Ped2 = Ped.ref,
                                             SNPd = SimSNPd[[i]],
                                             Symmetrical=FALSE, Plot=FALSE)$Counts.detail
  }
  
  Ntot <- apply(PC.rev.cd, c('id.cat', 'dam.cat', 'sire.cat'), sum)
  OK <- list('G' = 'Match',
             'D' = 'Match',
             'X' = c('P2only', '_'))
  confA <- array(dim = c(3,3,3,3),
                 dimnames = c(list(paste0(c('dam', 'sire', 'pair'), '.conf')),
                              dimnames(PC.rev.cd)[2:4]))
  for (i in c('G','D','X')) {
    confA['dam.conf' ,,i,] <- apply(PC.rev.cd[,,i,,OK[[i]],], c('id.cat', 'sire.cat'), sum) / Ntot[,i,]
    confA['sire.conf',,,i] <- apply(PC.rev.cd[,,,i,,OK[[i]]], c('id.cat', 'dam.cat'), sum) / Ntot[,,i]
    for (j in c('G','D','X')) {
      confA['pair.conf',,i,j] <- apply(PC.rev.cd[,,i,j,OK[[i]],OK[[j]]], 'id.cat', sum) / Ntot[,i,j]
    }
  }
  
  confA[c("dam.conf" , "pair.conf"),,"X",] <- NA  # no dam
  confA[c("sire.conf", "pair.conf"),,,"X"] <- NA  # no sire
  
  Conf.df <- plyr::adply(confA, .margins=2:4)
  Conf.df <- merge(Conf.df,
                   plyr::adply(Ntot, .margins=3:1, function(x) data.frame(N=x)))
  Conf.df <- Conf.df[Conf.df$id.cat != 'X',]
  if (ParSib == "par") {
    Conf.df <- Conf.df[Conf.df$id.cat == 'G' & Conf.df$dam.cat %in% c('G','X') &
                         Conf.df$sire.cat %in% c('G','X'), ]
  }
  Conf.df <- Conf.df[order(Conf.df$id.cat, Conf.df$dam.cat, Conf.df$sire.cat), ]
  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # assignment errors (ignores co-parent) ----
  PedComp.fwd <-  array(0, dim=c(nSimz, 7,5,2),
                        dimnames = list(iter = seq_len(nSimz),
                                        cat = c("GG", "GD", "GT", "DG", "DD", "DT", "TT"),
                                        class = c("Total", "Match", "Mismatch", "P1only", "P2only"),
                                        parent = c("dam", "sire")))
  for (i in 1:nSim) {
    PedComp.fwd[i,,,] <- sequoia::PedCompare(Ped1 = Ped.ref,
                                             Ped2 = Pedigree.inferred[[i]],
                                             SNPd = SimSNPd[[i]],
                                             Symmetrical=FALSE, Plot=FALSE)$Counts
  }
  
  PedComp.tmp <- apply(PedComp.fwd, 2:4, sum)
  PedErrors <- sweep(PedComp.tmp[,c("P1only", "P2only","Mismatch"),], c(1,3),
                     PedComp.tmp[,"Total",], "/")
  PedErrors[c("GG", "GD", "DG", "DD"), "P2only", ] <- NA  # if parent in ref. pedigree, by def not P2only
  dimnames(PedErrors)[['class']] <- c("FalseNeg", "FalsePos", "Mismatch")
  
  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # out ----
  RunParams <- list(EstConf = sequoia:::namedlist(args.sim, args.seq, nSim, nCores),
                    SimGeno_default = formals(sequoia::SimGeno),
                    sequoia_default = formals(sequoia::sequoia),
                    sequoia_version = as.character(utils::packageVersion("sequoia")))
  
  return( list(ConfProb = Conf.df,
               PedErrors = PedErrors,
               Pedigree.reference = Ped.ref,
               LifeHistData = LifeHistData,
               Pedigree.inferred = Pedigree.inferred,
               SimSNPd = SimSNPd,
               PedComp.fwd = PedComp.fwd,
               RunParams = RunParams,
               CallRate.sim = CallRate.sim,
               groupX = groupX,
               RunTime = RunTime) )
}






################################################################################
################################################################################
################################################################################

#' @title Calculate confidence, split by group of id/dam/sire (from 
#'   \code{EstConf2})
#'
#' @description Post-processing of output from \code{EstConf2}, to split 
#'   confidence and assignment rate estimates by 'group', analogous to G/D/X.
#'
#' @param Conf_out output from \code{EstConf2}, includes list element 'groupX'.
#'
#' @return a large array
#'
#' @details 
#' Confidence:  proportion of *assigned* parents that are correct
#' Assignment rate : proportion of *true* parents that are correctly assigned
#' 

Conf_by_group <- function(Conf_out)
{
  nSim <- length(Conf_out$Pedigree.inferred)
  groupX <- Conf_out$groupX
  Ped.ref <- Conf_out$Pedigree.reference
  
  groupNames <- c('full','half','X')
  classNames <- c("Match", "Mismatch", "P1only", "P2only", "_")
  
  CalcPropMatch <- function(cnts) {
    A.out <- array(dim=c(2,3,3,3))
    for (v in 1:2) {  # id
      for (w in 1:3) {  # dam
        for (x in 1:3) {  # sire
          A.out[v,w,x,1] <- sum(cnts[v,w,x,'Match',])/sum(cnts[v,w,x,,])  # dam
          A.out[v,w,x,2] <- sum(cnts[v,w,x,,'Match'])/sum(cnts[v,w,x,,])  # sire
          if (w=='X' | x=='X')  next
          A.out[v,w,x,3] <- sum(cnts[v,w,x,'Match','Match'])/sum(cnts[v,w,x,,])  # pair
        }
      }
    }
    return(A.out)
  }
  
  ConfAR <- array(NA, dim=c(2,2,3,3,3,nSim),
                  dimnames = list(metric=c('conf', 'AR'),id.group = groupNames[1:2],
                                  dam.group = groupNames, sire.group = groupNames,
                                  c('dam.conf', 'sire.conf', 'pair.conf'), 1:nSim))
  
  # Confidence:  proportion of *assigned* parents that are correct
  for (i in 1:nSim) {
    Ped.i <- Conf_out$Pedigree.inferred[[i]][, 1:3]
    
    for (p in c('id', 'dam', 'sire')) {
      Ped.i[, paste0(p, '.group')] <- ifelse(Ped.i[,p] %in% groupX[[i]], 'half',
                                             ifelse(!is.na(Ped.i[,p]), 'full', 'X'))
      Ped.i[, paste0(p, '.group')] <- factor(Ped.i[, paste0(p, '.group')], levels = groupNames)
    }
    
    PedM <- sequoia::PedCompare(Ped.i, Ped.ref,
                                SNPd =Conf_out$SimSNPd[[i]],
                                Symmetrical=FALSE, Plot=FALSE)$MergedPed
    for(p in c('dam', 'sire')) {
      PedM[,paste0(p,'.class')] <- factor(PedM[,paste0(p,'.class')], levels = classNames)
    }
    
    PedM <- merge(PedM, Ped.i[, c('id', 'id.group', 'dam.group', 'sire.group')])
    
    cnts <- with(PedM, table(id.group, dam.group, sire.group, dam.class, sire.class))
    ConfAR['conf',,,,,i] <- CalcPropMatch(cnts)
  }
  
  
  # Assignment rate : proportion of *true* parents that are correctly assigned
  for (i in 1:nSim) {
    Ped.i <- Conf_out$Pedigree.inferred[[i]][, 1:3]
    
    for (p in c('id', 'dam', 'sire')) {
      Ped.ref[, paste0(p, '.group')] <- ifelse(Ped.ref[,p] %in% groupX[[i]], 'half',
                                               ifelse(!is.na(Ped.ref[,p]), 'full', 'X'))
      Ped.ref[, paste0(p, '.group')] <- factor(Ped.ref[, paste0(p, '.group')], levels = groupNames)
    }
    
    PedM <- sequoia::PedCompare(Ped.ref, Ped.i,
                                SNPd =Conf_out$SimSNPd[[i]],
                                Symmetrical=FALSE, Plot=FALSE)$MergedPed
    
    PedM <- merge(PedM,
                  Ped.ref[, c('id', 'id.group', 'dam.group', 'sire.group')])
    for(p in c('dam', 'sire')) {
      PedM[,paste0(p,'.class')] <- factor(PedM[,paste0(p,'.class')], levels = classNames)
    }
    
    cnts <- with(PedM, table(id.group, dam.group, sire.group, dam.class, sire.class))
    ConfAR['AR',,,,,i] <- CalcPropMatch(cnts)
  }
  
  return(ConfAR)
}