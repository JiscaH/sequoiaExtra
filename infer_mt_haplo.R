# Infer mitochondrial haplotype for as many individuals as possible based on
# a pedigree + known mt haplotype for some individuals.
# Infers for both ancestors, descendants, and other matrilineal relatives.

# IN:
# Pedigree: data.frame with columns id - dam - sire
# mtHaps: named vector

# OUT:
# mtHaps extended with haplotypes for additional individuals, if any were found.


infer_mt_haplotype <- function(Pedigree, mtHaps) {
  Ped <- sequoia::PedPolish(Pedigree, KeepAllColumns=FALSE)
  Ped$G <- sequoia::getGenerations(Ped)
  nG <- max(Ped$G)

  # go bottom -> top to derive haplotypes of ancestors,
  # then top -> bottom to get haplotypes of all their descendants
  gg <- c(nG:1, 1:nG)
  for (x in seq_along(gg)) {
    these <- with(Ped, G==gg[x] & !is.na(dam))
    if (x <= nG) {  # way up
      these <- these & Ped$id %in% names(mtHaps)
      mtHaps_new <- with(Ped, setNames( mtHaps[id[these]], dam[these] ) )
      mtHaps_new <- mtHaps_new[!duplicated(names(mtHaps_new))]
    } else {  # way down
      these <- these & Ped$dam %in% names(mtHaps)
      mtHaps_new <- with(Ped, setNames( mtHaps[dam[these]], id[these] ) )
    }
    if (all(is.na(mtHaps_new)))  next   # especially possible on way up

    # check for conflicts
    z <- intersect(names(mtHaps), names(mtHaps_new))
    if (any(mtHaps[z] != mtHaps_new[z])) {
      stop('Conflicting mt haplotype for individuals: ', names(which(mtHaps[z] != mtHaps_new[z])),
           ' Please check maternal pedigree links and/or mt haplotypes!')
    }

    mtHaps <- c(mtHaps, mtHaps_new)
    mtHaps <- mtHaps[!duplicated(names(mtHaps))]
  }

  return( mtHaps )
}



#===============================================================================
# simulate mtSame matrix for a pedigree, assuming all female founders have a
# unique mitochondrial haplotype

sim_mtSame <- function(Pedigree)
{
  library(magrittr)  # for %>%

  # check & format pedigree
  Pedigree <- sequoia::PedPolish(Pedigree)
  # get generation number for each individual
  Pedigree$G <- sequoia::getGenerations(Pedigree)

  # make up unique mt haplotype for each founder
  mtHaps_founders <- outer(LETTERS, LETTERS, FUN='paste0') %>%
    sample(size = sum(Pedigree$G==0)) %>%
    setNames(Pedigree$id[Pedigree$G==0])
  # infer mt hap for all descendants
  mtHaps_all <- infer_mt_haplotype(Pedigree, mtHaps_founders)
  # add to pedigree
  Pedigree$mt <- mtHaps_all[Pedigree$id]
  # create matrix indicating for each pair of individuals whether they have the
  # same (1) or different (0) haplotype
  mtSame <- outer(mtHaps_all, mtHaps_all, FUN = function(x,y) as.numeric(x==y))

  return( mtSame )
}
