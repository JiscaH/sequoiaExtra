# divide genotype data into N blocks of SNPs, do imputation in parallel, and combine results

parallel_impute <- function(geno_file = NA,      # path to genotype file, without file extension.
                                                 #   Currently only ped/map supported
                            pedigree_file = NA,  # file with id - parent1 - parent2
                            install_dir = NA,    # directory with imputator.f90, sqa_fileIO.fio, & sqa_general.f90
                            working_dir = NA,    # directory where the N temporary genotype subfiles will the created
                            n_cores = 1,         # number of cores = number of subfiles
                            compile = TRUE,
                            compilation_args = '-O3',  # optimisation level, or debugging flags
                            plink_cmd = 'plink',  # command to run plink
                            imputator_args = '--err 0.001 --T-impute 0',
                            output_file = 'geno_imputed',
                            quiet = FALSE)      # messages in R
{
  # if geno_file & pedigree_file are file names without path:
  # specify full path, in current working directory
  # (will be changing working directory a few times)
  current_dir <- getwd()
  if (dirname(geno_file) == '.')  geno_file <- file.path(current_dir, geno_file)
  if (dirname(pedigree_file) == '.')  pedigree_file <- file.path(current_dir, pedigree_file)
  if (dirname(output_file) == '.')  output_file <- file.path(current_dir, output_file)

  if (!dir.exists(working_dir))  dir.create(working_dir)

  files_needed <- c(paste0(geno_file,'.ped'), paste0(geno_file,'.map'), pedigree_file)
  for (f in files_needed) {
    if (!file.exists(f))  stop('could not find ', f)
  }


  # === compile imputator ====

  if (compile) {
    # check if imputator source files can be found
    setwd(install_dir)
    for (sf in c('imputator.f90', 'sqa_fileIO.f90', 'sqa_general.f90')) {
      if (!file.exists(sf)) stop("Could not find source file ", sf)
    }

    # compile
    system(glue::glue("gfortran {compilation_args} sqa_general.f90 sqa_fileIO.f90 imputator.f90 -o imputator"))
    # TODO? clean up .mod files
  }

  # === check number of cores ===

  if (!is.null(n_cores) && (!is.numeric(n_cores) || n_cores<1))
    stop("n_cores must be a positive number, or NULL")

  if (is.null(n_cores) || n_cores>1) {
    if (!requireNamespace("parallel", quietly = TRUE)) {
      if (interactive() & !quiet) {
        message("Installing pkg 'parallel' to speed things up... ")
      }
      utils::install.packages("parallel")
    }
    max_cores <- parallel::detectCores()
    if (is.null(n_cores)) {
      n_cores <- max_cores -1
    } else if (n_cores > max_cores) {
      n_cores <- max_cores
      warning("Reducing 'n_cores' to ", max_cores, ", as that's all you have",
              immediate.=TRUE)
    }
    if (!quiet)  message("Using ", n_cores, " out of ", max_cores, " cores")
  }


  # === split genotype file ====

  # use plink option --extract, requires text file with 1 column with SNP names
  # https://www.cog-genomics.org/plink/1.9/filter#snp
  # get SNP names from .map file, which has SNP name in 2nd column
  # (https://www.cog-genomics.org/plink/1.9/formats#map)

  # check if plink is working
  plink_test <- system(glue::glue("{plink_cmd} --version"))
  if (plink_test != 0)  stop(plink_cmd, " not found")

  setwd(working_dir)
  library(dplyr)

  # run ped cleaning on full dataset
  # TODO

  # divide SNPs into (nearly) equal batches
  map <- read.table(glue::glue("{geno_file}.map"))
  n_snps <- nrow(map)
  batch_size <- rep(floor(n_snps/n_cores), n_cores)
  n_left <- n_snps - sum(batch_size)
  batch_size[1:n_left] <- batch_size[1:n_left] +1
  snp_batches <- lapply(1:n_cores, function(x) seq.int(from = c(0,cumsum(batch_size))[x] +1,
                                                       length.out = batch_size[x]))

  # vector with batch names, for use in file names
  batch_names <- sapply(1:n_cores, function(x) glue::glue('batch_{formatC(x,width=3,flag=0)}'))

  for (x in 1:n_cores) {
    write.table(map[snp_batches[[x]], 2],
                file=glue::glue('snp_{batch_names[x]}.txt'),
                row.names=FALSE, col.names=FALSE, quote=FALSE)
  }

  # call plink to create subsets
  if (!quiet)  message('splitting genotype file ...')
  for (x in 1:n_cores) {
    plink_stuff <- system(glue::glue('{plink_cmd} --file {geno_file} --extract snp_{batch_names[x]}.txt ',
                   '--recode --out geno_{batch_names[x]}'), intern=TRUE)  # --recode makes .ped + .map
  }
  # clean up .nosex & .log files
  system('rm geno_batch_*.nosex geno_batch_*.log')


  # === impute each batch ===
  if (!quiet)  message('imputing ...')

  run_impute <- function(x, imp_args, i_dir, ped_file)
  {
    batch_name_x <- glue::glue('batch_{formatC(x,width=3,flag=0)}')
    # specify the command separately, otherwise it doesn't seem to recognize that the
    # input arguments are being used
    cmd <- glue::glue('{i_dir}/imputator {imp_args} --pedigree {ped_file} ',
                   '--geno geno_{batch_name_x} --edits-out {batch_name_x}.edits ',
                   '--informat PED --no-geno-out --quiet')
    system(cmd)
    return('OK')
  }

  cl <- parallel::makeCluster(n_cores)
  AllOUT <- parallel::parLapply(cl, X = seq.int(n_cores),
                                fun = run_impute,
                                imp_args = imputator_args,
                                i_dir = install_dir,
                                ped_file = pedigree_file)
  parallel::stopCluster(cl)


  # === combine imputation logs ===
  # NOTE: these files may possibly to get too big to be handled (comfortably) by R,
  # then they need to be combined using cat & headers somehow discarded
  # https://www.tutorialspoint.com/how-to-append-contents-of-multiple-files-into-one-file-on-linux
  if (!quiet)  message('combining edit lists ...')
  edits_L <- list()
  for (x in 1:n_cores) {
    edits_L[[x]] <- read.table(glue::glue("{batch_names[x]}.edits"), header=TRUE)
  }
  edits_all <- plyr::ldply(edits_L)
  # change SNP column indices to correspond to the original genotype file
  edits_all$snp_index <- match(edits_all$snp_name, map[,2])

  write.table(edits_all, file='combined.edits', row.names=FALSE, col.names=TRUE, quote=FALSE)


  # === apply edits to big genotype file  ===
  if (!quiet)  message('applying edits ...')
  system(glue::glue('{install_dir}/imputator {imputator_args} ',  # repeat, in case it includes outFormat
                   '--geno-in {geno_file} --edits-in combined.edits ',
                   '--informat PED --quiet --geno-out {output_file}'))


  # === clean up  ===
  system('rm snp_batch_*.txt geno_batch_*.ped geno_batch_*.map batch_*.edits')
  setwd(current_dir)

}
