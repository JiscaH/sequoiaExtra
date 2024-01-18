# Manual for `grm_tool.f90`

### A fortran program to summarise and filter GRMs

Jisca Huisman, 2024-01-14

## Description

This program can calculate various summary statistics from (potentially
huge) GRMs, and/or filter out pairs with R values above or below a
specified threshold. Additional functionality may be added in the
future.

## Preparations

### Requirements

- a Fortran compiler, e.g. gfortran
- gzip (or pigz), to decompress the .grm.gz file

### Compilation

Run `gfortran -O3 grm_tool.f90 -o grm_tool` (where you may substitute
gfortran by another compiler of your choice)

## Input files

- `--in` : input file pair with GRM, as returned by
  `plink --make-grm-gz` or `gcta --make-grm`. Extensions .grm.id and
  .grm.gz are added.
- `--only` : individual subset, only pairs where either individual is
  listed are considered. A text file with no header, and either two
  columns with IDs in the second column (same format as for
  `plink --keep`), or a single column with IDs (any additional columns
  are ignored if there are \>2 columns).
- `-only-among` : as `--only`, but only pairs where *both* individuals
  are listed are considered. grm

## Program options

### Summary statistics

Currently, the only command line options are to turn calculation off
(`--no-summary`) and to specify the outfile (`--summary-out`). Turning
the calculation off potentially allows filtering of GRMs with more
individuals, see section ‘Method’.

The summary output currently looks as follows:

              Group           Part        N_pairs         N_SNPs   minimum        maximum        mean           std_dev        count_<-0.5    count_<0.625   count_>0.875   count_>1.25 
              total       diagonal            200         495.07       0.796882       1.265416       0.984194       0.077524              0              0            189              1
              total        between          19900         490.21      -0.221723       0.756767      -0.004994       0.101081              0          19892              0              0
             subset       diagonal             11         495.18       1.001152       1.155578       1.077955       0.043230              0              0             11              0
             subset        between           2134         490.32      -0.199109       0.589298      -0.006551       0.085294              0           2134              0              0

where

- total: in the full GRM
- subset: among members of the `--only` list, and unless `--only-among`
  also between members of the subset and all other individuals

When the program is run without `--only`, the ‘subset’ rows are omitted,
but the columns are identical.

The count thresholds can be changed easily (search for ‘sumstats_counts’
in the source code); do not forget to also change the output column
header accordingly, via ‘sumstat_lbls’. Adding more summary statistics
is also quite straight forward; please contact me for assistance if you
are not familiar (enough) with Fortran to make the necessary edits
yourself.

Note that when reading this file into R, its default behaviour is to
change the ‘\<’, ‘\>’ and ‘-’ in the column names into a dot. To prevent
this, use `read.table(file, header=TRUE, check.names=FALSE)`.

### Filtering

A subset of the GRM can be written to a text file based on their R
values, e.g. to check any weird values (e.g. R \> 1.5 or R \< -1.0) or
to find off-diagonal pairs which may be 1st degree relatives but not
duplicates or inbred (e.g. 0.35 \< R \< 0.65).

NOTE: This is meant as a tool for data quality control, and I strongly
advise against using this as primary method to identify relatives. For
that, please use sequoia (<https://github.com/JiscaH/sequoia_notR>),
`find_pairs.f90` (in <https://github.com/JiscaH/sequoiaExtra> ), or
other likelihood-based methods.

Separately for the diagonal and off-diagonal (‘between’) elements, one
or two thresholds can be set:

- `diag-lower`, `betw-lower` : the lower threshold, pairs with values
  above this thresholds are written to the output file;
- `diag-upper`, `betw-upper` : the upper threshold.

One threshold may be left unspecified, and the `upper` threshold may be
higher or lower than the `lower` threshold, giving four different
possibilities as illustrated as A–D in the figure (grey = selected pairs
written to the output).

<img src="grm_tool_manual_files/figure-gfm/filter_illustration-1.png" width="60%" />

#### Infinity

NOTE: Filtering will currently ignore any values of `+Inf` or `-Inf`, as
the thresholds do not default to infinity but to ‘HUGE(0D0)’, which is
the largest ‘double precision’ value that Fortran can store. Whether
those are present in the data can be checked in the summary statistics
as

``` fortran
write(42,*)  'count_inf ',  COUNT(GRM >= HUGE(0D0))
```

If this is an issue, please let me know as it would be fairly
straightforward to change.

### Histogram

Counts per histogram bin are calculated with `--hist`, with separate
counts for the diagonal and off-diagonal (‘between’) elements written to
a 3-column text file (default hist_counts.txt, specify with
`--hist-out`). The first column in the output is the lower bound of each
bin. These are equidistant and default from -1.5 to +2.0 with a stepsize
of 0.05. These can be adjusted as optional arguments to `--hist`, in
which case all 3 must be provided, in the order first, last, step.

The output is (should be) identical to the output of the R command

``` r
  table(cut(GRM, breaks=seq(first, last, step)))
```

All bins are closed on the right and open on the left, any values \<=
first or \> last are discarded, with a warning. It might therefore be
advisable to edit the summary statistics to count the number of values
on and off the diagonal exceeding ‘last’. Note that the warning will be
given only after all the data has been processed, which may take more
than an hour for very large GRMs.

The output can e.g. be visualised in R as follows:

``` r
H <- read.table('hist_counts.txt', header=TRUE)
# remove superflous head & tail
min_R <- with(H, min(lower_bound[diagonal>0 | between>0]))
max_R <- with(H, max(lower_bound[diagonal>0 | between>0]))
H <- H[H$lower_bound >= min_R & H$lower_bound <= max_R, ]

# plot
plot(H$lower_bound, H$diagonal, type='h', xlim=c(-.4, 1.5),
    main='Diagonal', ylab='Count', xlab='R' )
# - or - 
barplot(H$between, space=0, names.arg=H$lower_bound,
        main='Off-diagonal', ylab='Count', xlab='R')
# - or -
HH <- list(breaks = c(H$lower_bound,
                      H$lower_bound[nrow(H)] +
                        (H$lower_bound[2]-H$lower_bound[1])),
           counts = H$diagonal,
           xname = 'R diagonal')
class(HH) <- 'histogram'
plot(HH, xlim=c(0, 1.3))
```

#### with `--only`

When the program is run with `--only`, the output file contains the
following columns:

- lower_bound: lower bound of the histogram bin
- total_diagonal: counts on the diagonal of the full GRM
- total_between: counts off-diagonal of the full GRM
- subset_diagonal: counts on the diagonal within the subset
- subset_between: counts off-diagonal, between members of the subset and
  unless `--only-among` also between members of the subset and all other
  individuals

Currently the summary statistics include only those individuals where
either or both individuals are on the `--only` list. A planned future
upgrade is to return two summaries (and two sets of histograms): one for
the full dataset, and one for the `--only` subset.

## Method

The compressed .grm.gz file is decompressed and processed as a
continuous stream of data via a so-called ‘named pipe’, as described at
<https://genomeek.wordpress.com/2012/05/07/tips-reading-compressed-file-with-fortran-and-named-pipe/>

The speed of the program is limited by the decompression speed. Whereas
compression of data can be done in parallel on multiple threads,
decompression cannot (according to <https://zlib.net/pigz/pigz.pdf> )

For the filtering, the data is taken from the stream, checked against
the specified thresholds, written to file (or not), and discarded. The
amount of computer memory required by the program does thus not increase
with the size of the GRM, although the output file may get very large
depending on the thresholds used.

For the summary statistics and `--hist`, after 5%, 10%, …, of the data
is read in, summary statistics are calculated and memory space reused
for the next 5% of data. After all data are processed, the
sum/min/max/mean is calculated over the 21 data chunks (20 equal sized +
left overs). For the standard deviation, the ‘extension to K groups’ at
<https://stats.stackexchange.com/questions/55999/is-it-possible-to-find-the-combined-standard-deviation?noredirect=1&lq=1>
is used.

The program currently uses gzip for decompression, but this can easily
be changed by changing `gzip` to e.g. `pigz` at the following line in
the source code:

``` fortran
call EXECUTE_COMMAND_LINE("(gzip -dc  "//trim(grmFile)//".grm.gz > grmpipe) &")
```

## Example data

A small mock GRM is available to try out this program and its various
setting. The files ‘griffin.grm.gz’ and ‘griffin.grm.id’ contain a GRM
generated from SNP data simulated from ‘Ped_griffin’ in the sequoia R
package. It includes 200 individuals, including 2 inbred pairs
(R\>0.75).
