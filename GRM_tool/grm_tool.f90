! read in a GRM file pair as returned by PLINK --make-grm-gz,
! and return various summary statistics, calculated on ln 286-297

! by: Jisca Huisman, Sept. 2023  jisca.huisman@gmail.com
! updated Jan. 2024: separate counts of R>1.5 etc. on diagonal vs off-diagonal

! The named pipe approach is derived from
! https://genomeek.wordpress.com/2012/05/07/tips-reading-compressed-file-with-fortran-and-named-pipe/
 
 
! .grm.gz data format:
! 1. index of first sample
! 2. index of second sample
! 3. Number of variants where neither sample has a missing call
! 4. Relationship value

! .grm.id data format:
! 1. FID (family ID, not used)
! 2. IID


! compile: gfortran -O3 grm_tool.f90 -o grmtool
! debug: gfortran -Wall -pedantic -fbounds-check -g -Og grm_tool.f90 -o grmtool

!===============================================================================
module Global_vars
  implicit none
  
  integer, parameter :: nchar_filename = 2000, nchar_ID = 40
  integer,parameter :: ik10 = selected_int_kind(10)
  integer(kind=ik10) :: nrows_grm
  integer, allocatable :: nSnp(:), indx(:,:)
  character(len=nchar_ID), allocatable :: ID(:)
  logical :: DoSummary, DoFilter(2), OnlyAmong, quiet
  logical, allocatable :: skip(:), IsDiagonal(:), InSubset(:)
  double precision :: lowr_d, upr_d, lowr_b, upr_b
  double precision, allocatable :: GRM(:)
  
end module Global_vars

!===============================================================================

module Fun
  use Global_vars
  implicit none

  contains
  
  ! option documentation
  subroutine print_help()
    print '(a, /)', ' ~~ Calculate summary statistics and/or filter pairs from (very) large GRMs ~~'
    print '(a, /)', 'command-line options:'
    print '(a)',    '  --help         print usage information and exit'
    print '(a)',    '  --in <grmfile>    input file with GRM; extensions .grm.id and .grm.gz are added' 
    print '(a)',    '  --out_prefix <file>  prefix for output files; defaults to <grmfile> '    
    print '(a)',    '  --summary-out <file>  output file for summary statistics, default: ', &
                    '                    <grmfile>_summary_stats.txt ' 
   print '(a)',    '   --no-summary   do not calculate various summary statistics, e.g. when GRM is ', &
                    '                    too large to fit in memory. See/adjust the calculated ', &
                    '                    statistics on line 351-364 of the source code'
    print '(a)',    '  --hist <f,l,s> counts per histogram bin, separated into diagonal and between ',&
                    '                    (off-diagonal). First, Last, and Step are optional arguments',&
                    '                    and default to -1.5, +2.0, and 0.05'
    print '(a)',    '  --hist-out <file>  output file for histogram counts, defaults to <grmfile>_hist_counts.txt'
    print '(a)',    '  --diag-lower   lower bound of R value of exported individuals (on diagonal)'
    print '(a)',    '  --diag-upper   upper bound of R value of exported individuals (on diagonal)'
    print '(a)',    '  --betw-lower   lower bound of R value between exported pairs (off-diagonal)'
    print '(a)',    '  --betw-upper   upper bound of R value between exported pairs (off-diagonal)'
    print '(a)',    '  --filter-out <file>  output file with pairs under lower / above upper threshold,'&
                    '                    defaults to <grmfile>_filter_output.txt'
    print '(a)',    '  --only <file>  only consider pairs with one or both individuals listed,',&
                    '                    IDs in first/single column or columns FID (ignored) + IID'
    print '(a)',    '  --only-among <file>  only consider pairs with both individuals listed'
    print '(a)',    '  --quiet        hide messages'
    print '(a)',    ''
  end subroutine print_help
    
     
  ! determine number of rows in a file
  integer function FileNumRow(FileName)
    character(len=*), intent(IN) :: FileName
    integer :: nrow, ios  
    character(len=42) :: dumC

    nrow = 0
    open(unit=102, file=trim(FileName), status="old")
    do 
      read(102, *, IOSTAT=ios) dumC
      if (ios < 0) then
        exit  ! EOF
      else
        nrow = nrow +1  
      end if
    enddo
    close(102)
    FileNumRow = nrow

  end function FileNumRow
  
  integer function FileNumCol(FileName)
    implicit none

    character(len=*), intent(IN) :: FileName
    integer :: j, strLen, numcol
    character(:), allocatable :: line

    allocate(character(200000) :: line)

    open(unit=102, file=trim(FileName), status="old")
    read(102, '(a)' ) line
    close(102) 

    strLen = len_trim(line)
    if (strLen >= 200000)  print *, 'WARNING: EXCEEDING MAXIMUM NUMBER OF SNPs!'
    if (strLen  == 0) then
      FileNumCol = 0
      return
    endif

    numcol = 0   ! first column (no space 'after')  achar(9) = \t
    do j=1, strLen-1
      if (j==1 .and. line(j:j) /= ' ' .and. line(j:j) /= achar(9)) then
        numcol = numcol +1
      endif
      if (line(j:j) == ' ' .or. line(j:j) == achar(9)) then
        if (line((j+1):(j+1)) /= ' ' .and. line((j+1):(j+1)) /= achar(9)) then
          numcol = numcol +1    ! n(ew column starts at j+1
        endif
      endif
    enddo
    FileNumCol = numcol  

    deallocate(line)

  end function FileNumCol
  
  
  ! round number x to d significant digits
  integer function roundit(x,d)
    double precision, intent(IN) :: x
    integer, intent(IN) :: d
    double precision :: z
    integer :: i
    
    z = x
    i = 0
    do
      i = i+1
      z = z/10.0
      if (z < 10**d)  exit
    enddo
  
    roundit = NINT(z) * 10**i
      
  end function roundit
  
  
  ! calculate standard deviation
  double precision function SD(x, mask)
    double precision, intent(IN) :: x(:)
    logical, intent(IN) :: mask(:)
    integer :: n
    double precision :: mean
    
    n = count(mask)   ! size(x)
    mean = sum(x, mask=mask)/n
    SD = sqrt(sum((x - mean)**2, mask=mask)/n)

  end function SD
  
  
  subroutine timestamp(spaceIN)
    implicit none
    
    logical, intent(IN), OPTIONAL :: spaceIN
    logical :: space
    integer :: date_time_values(8)
    
    if(present(spaceIN))then
      space=spaceIN
    else
      space=.FALSE.
    endif
    ! NOTE: print *, & write(*,*)  have initital space, write(*,'(...)')  does not
    
    call date_and_time(VALUES=date_time_values)
    write(*,'(i2.2,":",i2.2,":",i2.2, 1X)', advance='no') date_time_values(5:7)
    if (space) write(*, '(1X)', advance='no')
     
  end subroutine timestamp
  
  subroutine printt(text)
    character(len=*), intent(IN) :: text
  
    call timestamp()
    print *, text
  
  end subroutine printt

end module Fun

!===============================================================================

module histogram
  use Global_vars, ONLY: ik10, GRM
  implicit none

  contains
  ! counts per bin (histogram)
  function grm_hist(MASK, breaks)
    double precision, intent(IN) :: breaks(:)
    logical, intent(IN) :: MASK(:)
    integer(kind=ik10), allocatable ::  grm_hist(:)
    integer :: nBins, b
    
    nBins = SIZE(breaks) -1
    
    if (any(GRM <= breaks(1))) print *, 'hist() WARNING: some data <= first'
    if (any(GRM > breaks(nBins)))   print *, 'hist() WARNING: some data > last'
    
    allocate(grm_hist(nBins))
    grm_hist = 0  
    do b=1, nBins
      grm_hist(b) = COUNT(GRM > breaks(b) .and. GRM <= breaks(b+1) .and. MASK) 
    enddo

 !   deallocate(hist)  DO NOT DEALLOCATE
    
  end function grm_hist
  
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  function breaks(first, last, step)
    double precision, intent(IN) :: first, last, step 
    integer :: nBins, b
    double precision, allocatable :: breaks(:)
    
    nBins = NINT((last - first)/step)
    allocate(breaks(nBins+1))
    breaks(1) = first
    do b=2, nBins+1
      breaks(b) = breaks(b-1) + step
    enddo
  
  end function breaks

end module histogram

!===============================================================================

program main
  use Fun
  use histogram
  implicit none
   
  integer :: x, i,j, nArg, nInd, nBins
  double precision :: hist_opts(3)
  character(len=32) :: arg, argOption
  character(len=2) :: chk
  character(len=15) :: sumstat_lbls(8)
  character(len=nchar_filename) :: grmFile, filterFile, summaryFile, onlyFile, histFile, outPrefix
  logical :: FileExists, DoHist
  logical, allocatable :: hist_mask(:), sum_mask(:)
  double precision, allocatable :: hist_brks(:)
  integer(kind=ik10), allocatable :: hist_counts(:,:)
  
  write(*,*) ''
  write(*,*) REPEAT('~',24)
  write(*,*) '   >>>  GRM TOOL  <<<'
  write(*,*) REPEAT('~',24)
  write(*,*) '||  Jisca Huisman  |  https://github.com/JiscaH  |  GNU GPLv3  ||'
  write(*,*) ''
  
  ! set default values
  grmFile = 'nofile'
  outPrefix = grmFile
  onlyFile = 'nofile'
  summaryFile = 'default'
  filterFile = 'default'
  histFile = 'default'
  
  OnlyAmong = .FALSE.
  DoSummary = .TRUE.
  DoFilter = .FALSE.
  DoHist = .FALSE.
  quiet = .FALSE.
  
  lowr_d = -HUGE(0D0)
  upr_d  = HUGE(0D0)
  lowr_b = -HUGE(0D0)
  upr_b  = HUGE(0D0)
  
  hist_opts = (/-1.5d0, 2.0d0, 0.05d0/)  ! first, last, step

  
  ! read command line arguments:
  nArg = command_argument_count()
  i = 0
  do x = 1, nArg
    i = i+1
    if (i > nArg)  exit
    call get_command_argument(i, arg)
    
    select case (arg)
      case ('--help')
        call print_help()
        stop
        
      case ('--in')
        i = i+1
        call get_command_argument(i, grmFile)
        
      case ('--out-prefix')
        i = i+1
        call get_command_argument(i, outPrefix)
        
      case ('--summary-out')
        i = i+1
        call get_command_argument(i, summaryFile)
        
      case ('--no-summary')
        DoSummary = .FALSE.
        
      case ('--hist')
        DoHist = .TRUE.
        call get_command_argument(i+1, argOption)
        read(argOption, '(a2)') chk
        if (chk /= '--' .and. argOption/='') then  ! optional arguments to --hist
          do j=1,3
            i = i+1
            call get_command_argument(i, argOption)
            if (argOption == '' .or. argOption=='--')  then
              print *, '--hist requires either 0 or 3 arguments: first, last, step'
              stop
            endif
            read(argOption, *) hist_opts(j)
          enddo
        endif       
        
      case ('--hist-out')
        i = i+1
        call get_command_argument(i, histFile)
        
      ! case ('--lower')
        ! i = i+1
        ! call get_command_argument(i, argOption)
        ! read(argOption, *)  lowr_d
        ! lowr_b = lowr_d
        
      ! case ('--upper')
        ! i = i+1
        ! call get_command_argument(i, argOption)
        ! read(argOption, *)  upr_d
        ! upr_b = upr_d
        
      case ('--diag-lower')
        i = i+1
        call get_command_argument(i, argOption)
        read(argOption, *)  lowr_d
        DoFilter(1) = .TRUE.
        
      case ('--diag-upper')
        i = i+1
        call get_command_argument(i, argOption)
        read(argOption, *)  upr_d
        DoFilter(1) = .TRUE.
      
      case ('--betw-lower')
        i = i+1
        call get_command_argument(i, argOption)
        read(argOption, *)  lowr_b
        DoFilter(2) = .TRUE.
        
      case ('--betw-upper')
        i = i+1
        call get_command_argument(i, argOption)
        read(argOption, *)  upr_b
        DoFilter(2) = .TRUE.
        
      case ('--filter-out')
        i = i+1
         call get_command_argument(i, filterFile)
        
      case ('--only ')
        i = i+1
        call get_command_argument(i, onlyFile)
      
      case ('--only-among')
        i = i+1
        call get_command_argument(i, onlyFile)
        OnlyAmong = .TRUE.
        
      case ('--quiet')
        quiet = .TRUE.
        
      case default
        print '(2a, /)', 'Unrecognised command-line option: ', arg
        call print_help()
        stop

    end select
  enddo        
  
  if (trim(outPrefix) == 'nofile') then
    outPrefix = grmFile
  endif
  if (trim(summaryFile)=='default')  summaryFile = trim(outPrefix)//'_summary_stats.txt'
  if (trim(filterFile)=='default')  filterFile = trim(outPrefix)//'_filter_output.txt'
  if (trim(histFile)=='default')  histFile = trim(outPrefix)//'_hist_counts.txt'
  
  ! input check ------
  if (grmFile == 'nofile') then
    write (*,*)  "Please specify an input file, without file extensions"
    write (*,*)
    call print_help()
    stop
  else
    inquire(file=trim(grmFile)//'.grm.id', exist = FileExists)
    if (.not. FileExists) then
      write(*,*)  "Input file ", trim(grmFile), ".grm.id not found"
      stop
    endif
    inquire(file=trim(grmFile)//'.grm.gz', exist = FileExists)
    if (.not. FileExists) then
      write(*,*)  "Input file ", trim(grmFile), ".grm.gz not found"
      stop
    endif
  endif

  if (DoHist .and. .not. DoSummary) then
    print *, '--hist cannot be combined with --no-summary'
    stop
  endif
  
  if (DoHist) then
    if (hist_opts(2) < hist_opts(1)) then
      print *, '--hist: last must be larger than first'
      stop
    endif
    if (hist_opts(3) <= 0d0 .or.  hist_opts(3) > (hist_opts(2) - hist_opts(1))) then
      print *, '--hist: step must be >0 and smaller than last - first'
      stop
    endif
  endif
  
  if (ANY(DoFilter) .and. .not. quiet) then  ! diagonal
    call timestamp()
    print *, 'Using the following filtering thresholds:'
    if (lowr_d > -HUGE(0D0))  write(*, "(3x, 'Diagonal: R >', f7.4)") lowr_d
    if ( upr_d <  HUGE(0D0))  write(*, "(3x, 'Diagonal: R <', f7.4)") upr_d
    if (lowr_b > -HUGE(0D0))  write(*, "(3x, 'Off-diagonal: R >', f7.4)") lowr_b
    if ( upr_b <  HUGE(0D0))  write(*, "(3x, 'Off-diagonal: R <', f7.4)") upr_b
    print *, ''
  endif
  
  if (.not. DoSummary .and. .not. any(DoFilter)) then
    print *, '--no-summary and no filter set: nothing to be done. See --help'
    stop
  endif
  
   ! if (DoSummary .and. .not. quiet) then
    ! inquire(file=trim(summaryFile), exist = FileExists)
    ! if (FileExists)  print *, 'WARNING: '//trim(summaryFile)//' will be overwritten' 
  ! endif
  ! if (DoFilter .and. .not. quiet) then
    ! inquire(file=trim(filterFile), exist = FileExists)
    ! if (FileExists)  print *, 'WARNING: '//trim(filterFile)//' will be overwritten' 
  ! endif
  ! if (DoHist .and. .not. quiet) then
    ! inquire(file=trim(histFile), exist = FileExists)
    ! if (FileExists)  print *, 'WARNING: '//trim(histFile)//' will be overwritten' 
  ! endif

  
  ! read in .grm.id -----
  nInd = FileNumRow( trim(grmFile)//'.grm.id' )
  allocate(ID(nInd))
  
  open(12, file = trim(grmFile)//'.grm.id')
    do i=1, nInd
      read(12, *)  x, ID(i)
    enddo
  close(12)
  
  nrows_grm = (int(nInd, kind=ik10) * (nInd-1)/2 + nInd)
  
  if (.not. quiet) then
    call timestamp()
    print *, 'Read in IDs, N individuals =', nInd
    call timestamp()
    print *, '--> # rows grm =' ,  nrows_grm
  endif
  
  if (nrows_grm < 0) then
    print *, '# rows exceeds current storage mode!'
    print *, 'Increase selected_int_kind() on line 103'
    print *, ' or contact jisca.huisman@gmail.com'
    stop
  endif

  if (DoSummary) then
    allocate(GRM(nrows_grm))
    GRM = -999D0
    allocate(indx(2, nrows_grm))
    indx = 0
    allocate(nSnp(nrows_grm))
    nSnp = 0
    allocate(IsDiagonal(nrows_grm))
    IsDiagonal = .FALSE.  
    allocate(InSubset(nrows_grm))
    InSubset = .TRUE.
  endif
    
  
  ! read in --only list
  allocate(skip(nInd))
  if (trim(onlyFile)/= "nofile") then
    if (.not. quiet) then
      call printt("Reading individuals in --only file "//trim(onlyFile)//" ... ")
    endif
    call ReadOnlyList(onlyFile)
  else
    skip = .FALSE.
  endif

  ! read in GRM & filter
  call ReadFilterGRM(grmFile, filterFile)
 
  
  ! calculate & write summary statistics
  if (DoSummary) then
    sumstat_lbls(1) = 'minimum'
    sumstat_lbls(2) = 'maximum'
    sumstat_lbls(3) = 'mean'
    sumstat_lbls(4) = 'std_dev'
    sumstat_lbls(5) = 'count_<-0.5'
    sumstat_lbls(6) = 'count_<0.625'
    sumstat_lbls(7) = 'count_>0.875'
    sumstat_lbls(8) = 'count_>1.25'   ! TODO: check if labels read in OK by R  

    allocate(sum_mask(nrows_grm))
  
    open(42, file=trim(summaryfile), action='write')
!     write(42,*)  'grmFile   ',  trim(grmFile)
!      write(42,*)  'OnlyFile  ',  trim(onlyFile)
      write(42,'(4a15, 20(3x,a12))') 'Group', 'Part', 'N_pairs', 'N_SNPs', sumstat_lbls
      sum_mask = IsDiagonal
      write(42, '(2a15, i15, f15.2, 4f15.6, 4i15)') 'total', 'diagonal', &
        COUNT(sum_mask), mean_SNPs(sum_mask), sumstats(sum_mask), sumstats_counts(sum_mask)
      sum_mask = .not. IsDiagonal
      write(42, '(2a15, i15, f15.2, 4f15.6, 4i15)') 'total', 'between', &
        COUNT(sum_mask), mean_SNPs(sum_mask), sumstats(sum_mask), sumstats_counts(sum_mask)
      if (any(skip)) then 
        sum_mask = IsDiagonal .and. InSubset
        write(42, '(2a15, i15, f15.2, 4f15.6, 4i15)') 'subset', 'diagonal', &
          COUNT(sum_mask), mean_SNPs(sum_mask), sumstats(sum_mask), sumstats_counts(sum_mask)
        sum_mask = .not. IsDiagonal .and. InSubset
        write(42, '(2a15, i15, f15.2, 4f15.6, 4i15)') 'subset', 'between', &
          COUNT(sum_mask), mean_SNPs(sum_mask), sumstats(sum_mask), sumstats_counts(sum_mask) 
      endif
    close(42)
    if (.not. quiet) then
      call printt("summary statistics written to "//trim(summaryfile))
    endif
    deallocate(sum_mask)
  endif
  
  
  if (DoHist) then   
    hist_brks = breaks(first=hist_opts(1), last=hist_opts(2), step=hist_opts(3))
    nBins = size(hist_brks) -1
    if (any(skip)) then
      allocate(hist_counts(nBins,4))
    else
      allocate(hist_counts(nBins,2))
    endif
    allocate(hist_mask(nrows_GRM))
    
    hist_mask = IsDiagonal
    hist_counts(:,1) = grm_hist(MASK=hist_mask, breaks=hist_brks)
    hist_mask = .not. IsDiagonal
    hist_counts(:,2) = grm_hist(MASK=hist_mask, breaks=hist_brks)
    if (any(skip)) then
      hist_mask = IsDiagonal .and. InSubset
      hist_counts(:,3) = grm_hist(MASK=hist_mask, breaks=hist_brks)
      hist_mask = .not. IsDiagonal .and. InSubset
      hist_counts(:,4) = grm_hist(MASK=hist_mask, breaks=hist_brks)
    endif
    
    open(101, file=trim(histFile), action='write')
      if (any(skip)) then
        write(101, '(a12, 4a25)')  'lower_bound', 'total_diagonal', 'total_between', &
          'subset_diagonal', 'subset_between'
      else
        write(101, '(a12, 2a25)')  'lower_bound', 'diagonal', 'between'
      endif
      do x=1, nBins
        write(101, '(f12.4, 4i25)')  hist_brks(x), hist_counts(x,:)
      enddo
    close(101)

    if (.not. quiet)  call printt("histogram counts written to "//trim(histFile))   
  endif
  
 
  ! clean up
  if (allocated(ID))    deallocate(ID)  
  if (allocated(skip))  deallocate(skip)
  if (allocated(GRM))   deallocate(GRM)
  if (allocated(IsDiagonal))  deallocate(IsDiagonal)
  if (allocated(hist_brks))   deallocate(hist_brks)
  if (allocated(hist_mask))   deallocate(hist_mask)
  if (allocated(nSnp))   deallocate(nSnp)
  if (allocated(InSubset))  deallocate(InSubset)
  
  if (.not. quiet)  call printt('Done.')
  
  !=========================
  contains
    function mean_SNPs(MASK)
      logical, intent(IN) :: MASK(:)
      double precision :: mean_SNPs
      
      mean_SNPs = SUM(nSnp/1000D0, MASK=MASK) / COUNT(MASK=MASK) * 1000D0   ! got rounded to INF..     
    end function mean_SNPs
  
    function sumstats(MASK)
      logical, intent(IN) :: MASK(:)
      double precision :: sumstats(4)
      
      sumstats(1) = MINVAL(GRM, MASK=MASK)
      sumstats(2) = MAXVAL(GRM, MASK=MASK)
      sumstats(3) = SUM(GRM, MASK=MASK)/COUNT(MASK=MASK)  ! mean
      sumstats(4) = SD(GRM, MASK=MASK)   ! standard deviation; function defined in module Fun
    end function sumstats

    function sumstats_counts(MASK)
      logical, intent(IN) :: MASK(:)
      integer :: sumstats_counts(4)

      sumstats_counts(1) = COUNT(GRM < -0.5  .and. MASK)
      sumstats_counts(2) = COUNT(GRM < 0.625 .and. MASK)
      sumstats_counts(3) = COUNT(GRM > 0.875 .and. MASK)
      sumstats_counts(4) = COUNT(GRM > 1.25  .and. MASK)
    end function sumstats_counts

        
end program main

!========================================================================

subroutine ReadFilterGRM(grmFile, filterFile)
  use Global_vars
  use Fun, ONLY: timestamp, roundit
  implicit none
  
  character(len=*), intent(IN) :: grmFile, filterFile
  integer :: n, x, p, i, j, z, ios
  integer(kind=ik10) :: print_chunk, timing_y, y
  logical :: WritePair
  double precision :: r, CurrentTime(2)
  
  
    ! create & open named pipe with data from .grm.gz
  ! NOTE: EXECUTE_COMMAND_LINE() is fortran 2008 standard, 
  ! and possibly not supported by ifort. 
  ! SYSTEM() is gnu extension and possibly supported by both gfortran & ifort
  
  ! create a named pipe
  ! see https://www.linuxjournal.com/article/2156
  call EXECUTE_COMMAND_LINE("rm -f grmpipe ; mkfifo grmpipe")

  ! decompression instruction, this forms the flow into the pipe
  ! between brackets: run in separate subshell
  ! &: put the process in background
  !call EXECUTE_COMMAND_LINE("(pigz -dc  "//trim(grmFile)//".grm.gz > grmpipe) &")
  call EXECUTE_COMMAND_LINE("(gzip -dc  "//trim(grmFile)//".grm.gz > grmpipe) &")

  ! open a read (outflow) connection to the pipe
  open(11, file="grmpipe", action='read')  
  
  if (ANY(DoFilter)) then
    open(42, file=trim(filterFile), action='write')  
    write(42, '(2a10, 2X, 2a40, a10, a15)') 'index1', 'index2', 'ID1', 'ID2', 'nSNP', 'R'
  endif  
    
    n = 0
    x = 0
    print_chunk = roundit(nrows_grm/20D0,2)  ! print at approx every 5% progress
    p = 0
    timing_y = roundit(nrows_grm/50D0,1)   ! round at which to estimate & print total runtime
    call cpu_time(CurrentTime(1))

    do y = 1, nrows_grm
      
      if (.not. quiet) then
        if (MOD(y, print_chunk)==0) then
          p = p +1
          call timestamp()
          print *, y, '  ', p*5, '%'
        else if (y == timing_y) then
          call cpu_time(CurrentTime(2))
          print *, ''
          call timestamp(.TRUE.)
          write(*,'("Estimated total runtime (min): ", f7.1)')  (CurrentTime(2) - CurrentTime(1))*50/60
          print *, ''
        endif
      endif
      
      read(11, *, iostat=ios) i,j,z,r  
      if (ios/=0) exit   ! stop if end of file / incorrect entry
      if (skip(i) .and. skip(j)) then
        if (DoSummary) then
          InSubset(y) = .FALSE.
        else
          cycle 
        endif
      endif
      if (OnlyAmong) then
        if (skip(i) .or. skip(j)) then
          if (DoSummary) then
            InSubset(y) = .FALSE.
          else
            cycle
          endif
        endif
      endif
      x = x+1
      if (DoSummary) then
        indx(:,y) = (/i,j/)
        nSnp(y) = z
        GRM(y) = r
        if (i==j)  Isdiagonal(y) = .TRUE.
      endif
      
      if (.not. ANY(DoFilter)) cycle
      ! write to outfile entries that meet criteria
      WritePair = .FALSE.
      if (i==j .and. DoFilter(1)) then  ! diagonal
        if (upr_d > lowr_d) then
          if (r > lowr_d .and. r < upr_d)  WritePair = .TRUE.
        else
          if (r > lowr_d .or. r < upr_d)   WritePair = .TRUE.
        endif
      else if (i/=j .and. DoFilter(2)) then
        if (upr_b > lowr_b) then
          if (r > lowr_b .and. r < upr_b)   WritePair = .TRUE.
        else
          if (r > lowr_b .or. r < upr_b)   WritePair = .TRUE.
        endif
      endif      
     if (WritePair) then
       n = n+1
       write(42, '(2i10, 2X, 2a40, i10, e15.6)')  i, j, ID(i), ID(j), z, r
     endif 
  
    end do
  
  if (ANY(DoFilter))  close(42)                                          
  close(11)
  
   if (ANY(DoFilter) .and. .not. quiet) then
    write(*,*)  ""
    call timestamp()
    write(*,*)  "Found ", n, " pairs matching the criteria, written to "//trim(filterFile)
    write(*,*)  ""
  endif                                                                                    
  
  ! Warning if number of entries read from .grm.gz does not match number of 
  ! IDs read from .grm.id
  if (y < nrows_grm) then
    write(*,*)  ""
    write(*,*)  "   WARNING !!!"
    write(*,*)  "Number of entries read from "//trim(grmFile)//".grm.gz ", &
      "( ", y, " ) does not match number of individuals in ", &
      trim(grmFile)//".grm.id ( ", size(ID), " )"
    write(*,*)  ""
  endif
  
  ! remove named pipe
  CALL EXECUTE_COMMAND_LINE("rm -f grmpipe")
  
end subroutine ReadFilterGRM

!===============================================================================

subroutine ReadOnlyList(FileName)
  use Global_vars, ONLY: Id, skip
  use Fun
  implicit none

  character(len=nchar_filename), intent(IN) :: FileName
  integer :: x, i, nrows, IOerr, ncol
  character(len=nchar_ID) :: tmpC, tmpX

  nrows = FileNumRow(trim(FileName))  ! no header
  ncol  = FileNumCol(trim(FileName))
  
  if (.not. quiet) then
    if (ncol==2) then
      call printt("--only file in 2-column format, assuming IDs in column 2 ...")
    else
      call printt("--only file in 1-column or multi-column format, assuming IDs in column 1 ...")
    endif
  endif

  skip = .TRUE.

  ! single column (ignore all other columns)
!  call printt("Reading individuals in --only file "//trim(FileName)//" ... ")
  open(unit=103, file=trim(FileName), status="old")
    do x=1, nrows
      if (ncol==2) then   ! PLINK format: FID + IID column
        read(103, *,IOSTAT=IOerr)  tmpX, tmpC
      else
        read(103, *,IOSTAT=IOerr)  tmpC
      endif
      if (IOerr > 0) then
        print *, "Wrong input in file "//trim(FileName)//" on line ", x
        stop
      else if (IOerr < 0) then
        exit   ! EOF
      else
        do i=1, size(Id)
          if (Id(i) == tmpC) skip(i) = .FALSE.
        enddo
      endif
    enddo
  close(103)

  if (.not. quiet) then
    call timestamp()
    print *, 'read ', COUNT(.not. skip) ,' unique individuals present in GRM from --only file'
  endif

end subroutine ReadOnlyList

!===============================================================================

! end.