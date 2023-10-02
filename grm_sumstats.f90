! read in a GRM file pair as returned by PLINK --make-grm-gz,
! and return rows for which R < lower bound or R > upper bound.
! useful to identify and investigate outlier values, or 
! to select close relatives.

! by: Jisca Huisman, Sept. 2023  jisca.huisman@gmail.com

! The named pipe approach is derived from
! https://genomeek.wordpress.com/2012/05/07/tips-reading-compressed-file-with-fortran-and-named-pipe/
 


!===============================================================================
module Fun
  implicit none

  integer, parameter :: nchar_filename = 2000, nchar_ID = 40
  integer :: nInd
  character(len=nchar_ID), allocatable :: ID(:)
  logical, allocatable :: skip(:)

  contains
  
  ! option documentation
  subroutine print_help()
    print '(a, /)', 'command-line options:'
    print '(a)',    '  --help        print usage information and exit'
    print '(a)',    '  --in          input file with grm; extensions .grm.id and .grm.gz are added'
    print '(a)',    '  --out         output file'
    print '(a)',    '  --only        only consider pairs with one or both individuals listed'
    print '(a)',    '  --quiet       hide counter while running and message when done.'
    print '(a)',    ''
  end subroutine print_help
    
     
  ! determine number of rows in a file
  integer function FileNumRow(FileName)
    implicit none

    character(len=*), intent(IN) :: FileName
    integer :: nrow, ios  
    character(len=5000) :: dumC

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
  
  
  ! round number to 1 significant digit
  integer function roundit(x)
    implicit none
    
    double precision, intent(IN) :: x
    double precision :: z
    integer :: i
    
    z = x
    i = 0
    do
      i = i+1
      z = z/10.0
      if (z < 10)  exit
    enddo
  
    roundit = NINT(z) * 10**i
      
  end function roundit
  
  
  ! calculate standard deviation
  double precision function SD(x, mask)
    implicit none
    
    double precision, intent(IN) :: x(:)
    logical, intent(IN) :: mask(:)
    integer :: n
    double precision :: mean
    
    n = count(mask)   ! size(x)
    mean = sum(x, mask=mask)/n
    SD = sqrt(sum((x - mean)**2, mask=mask)/n)

  end function SD

end module Fun

!===============================================================================

program grm_sumstats
  use Fun
  implicit none
   
  integer :: x, i,j, z, n, ios, nArg, nrows_grm, print_chunk, p, timing_x
  double precision :: r, CurrentTime(2)
  character(len=32) :: arg, argOption
  character(len=nchar_filename) :: infile, outfile, OnlyListFileName
  logical :: FileOK, quiet
  double precision, allocatable, dimension(:) :: GRM
  integer, allocatable :: indx(:,:), nSnp(:)
  logical, allocatable, dimension(:) :: IsDiagonal
  
  ! set default values
  infile = 'nofile'
  outfile = 'grm_summary_stats.txt'
  OnlyListFileName = 'nofile'
  quiet = .FALSE.

  
  ! read command line arguments:
  nArg = command_argument_count()
  i = 0
  do x = 1, nArg
    i = i+1
    if (i > nArg)  exit
    call get_command_argument(i, arg)
    
    select case (arg)
      case ('-h', '--help')
        call print_help()
        stop
        
      case ('--in')
        i = i+1
        call get_command_argument(i, infile)
        
      case ('--out')
        i = i+1
        call get_command_argument(i, outfile)
        
      case ('--only ')
        i = i+1
        call get_command_argument(i, OnlyListFileName)
        
      case ('--quiet')
        quiet = .TRUE.
        
      case default
        print '(2a, /)', 'Unrecognised command-line option: ', arg
        call print_help()
        stop

    end select
  enddo        
  
  
  if (infile == 'nofile') then
    write (*,*)  "Please specify an input file, without file extensions"
    write (*,*)
    call print_help()
    stop
  else
    inquire(file=trim(infile)//'.grm.id', exist = FileOK)
    if (.not. FileOK) then
      write(*,*)  "Input file ", trim(infile), ".grm.id not found"
      stop
    endif
    inquire(file=trim(infile)//'.grm.gz', exist = FileOK)
    if (.not. FileOK) then
      write(*,*)  "Input file ", trim(infile), ".grm.gz not found"
      stop
    endif
  endif
  
  ! .grm.gz data format:
  ! 1. index of first sample
  ! 2. index of second sample
  ! 3. Number of variants where neither sample has a missing call
  ! 4. Relationship value
  
  ! .grm.id data format:
  ! 1. FID (family ID, not used)
  ! 2. IID
  
  
  ! read in .grm.id
  nInd = FileNumRow( trim(infile)//'.grm.id' )
  allocate(ID(nInd))
  
  open(12, file = trim(infile)//'.grm.id')
    do i=1, nInd
      read(12, *)  x, ID(i)
    enddo
  close(12)
  
  
  nrows_grm = (nInd * (nInd-1)/2 + nInd)
  allocate(GRM(nrows_grm))
  GRM = -999D0
  allocate(indx(2, nrows_grm))
  indx = 0
  allocate(nSnp(nrows_grm))
  nSnp = 0
  allocate(IsDiagonal(nrows_grm))
  IsDiagonal = .FALSE.
    
  
  ! read in --only list
  allocate(skip(nInd))
  if (trim(OnlyListFileName)/= "nofile") then
    if (.not. quiet)  print *, "Reading individuals in --only file "//trim(OnlyListFileName)//" ... "
    call ReadOnlyList(OnlyListFileName)
  else
    skip = .FALSE.
  endif


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
  !call EXECUTE_COMMAND_LINE("(pigz -dc  "//trim(infile)//".grm.gz > grmpipe) &")
  call EXECUTE_COMMAND_LINE("(gzip -dc  "//trim(infile)//".grm.gz > grmpipe) &")

  ! open a read (outflow) connection to the pipe
  open(11, file="grmpipe", action='read')  
    
    n = 0
    print_chunk = roundit(nrows_grm/20D0)  ! print at approx every 5% progress
    p = 0
    timing_x = roundit(nrows_grm/50D0)   ! round at which to estimate & print total runtime
    call cpu_time(CurrentTime(1))
    
    do x = 1, nrows_grm
      
      if (.not. quiet) then
        if (MOD(x, print_chunk)==0) then
          p = p +1
          print *, x, '  ', p*5, '%'
        else if (x == timing_x) then
          call cpu_time(CurrentTime(2))
          print *, ''
          write(*,'("Estimated total runtime (min): ", f7.1)')  (CurrentTime(2) - CurrentTime(1))*50/60
          print *, ''
        endif
      endif
      
      read(11, *, iostat=ios) i,j,z,r  
      if (ios/=0) exit   ! stop if end of file / incorrect entry
      if (skip(i) .and. skip(j))  cycle 
      indx(:,x) = (/i,j/)
      nSnp(x) = z
      GRM(x) = r
      if (i==j)  Isdiagonal(x) = .TRUE.
  
    end do
  
  close(11)
  
  
  ! Warning if number of entries read from .grm.gz does not match number of 
  ! IDs read from .grm.id
  if (x < (nInd * (nInd-1)/2 + nInd)) then
    write(*,*)  ""
    write(*,*)  "   WARNING !!!"
    write(*,*)  "Number of entries read from ", trim(infile),".grm.gz ", &
      "( ", x, " ) does not match number of individuals in ", &
      trim(infile),".grm.id ( ", nInd, " )"
    write(*,*)  ""
  endif
  
  ! remove named pipe
  CALL EXECUTE_COMMAND_LINE("rm -f grmpipe")
  
  
  ! calculate & write summary statistics
  open(42, file=trim(outfile), action='write')
    write(42,*)  'InFile    ',  trim(infile)
    write(42,*)  'OnlyFile  ',  trim(OnlyListFileName)
    write(42,*)  'N         ',  nInd
    write(42,*)  'mean_SNPs ',  SUM(nSnp/1000D0) / nrows_grm * 1000D0   ! got rounded to INF..
    write(42,*)  'min_R     ',  MINVAL(GRM)
    write(42,*)  'max_R     ',  MAXVAL(GRM)
    write(42,*)  'count_1.5 ',  COUNT(GRM > 1.5)
    write(42,*)  'count_5   ',  COUNT(GRM > 5)
    write(42,*)  'mean_diag ',  SUM(GRM, MASK=IsDiagonal)/COUNT(IsDiagonal)
    write(42,*)  'sd_diag   ',  SD(GRM, MASK=IsDiagonal)   ! function defined in module Fun
    write(42,*)  'mean_betw ',  SUM(GRM, MASK=(.not. IsDiagonal))/COUNT(.not. IsDiagonal)
    write(42,*)  'sd_betw   ',  SD(GRM, MASK=(.not. IsDiagonal)) 
  close(42)
  
  
  
  ! clean up
  if (allocated(ID))  deallocate(ID)  
  if (allocated(skip))  deallocate(skip)
  if (allocated(GRM))  deallocate(GRM)
  if (allocated(IsDiagonal))  deallocate(IsDiagonal)
        
end program grm_sumstats


!===============================================================================

subroutine ReadOnlyList(FileName)
  use Fun
  implicit none

  character(len=nchar_filename), intent(IN) :: FileName
  integer :: i, ios, x
  character(len=nchar_ID) :: tmpC

  skip = .TRUE.

  ! single column (ignore all other columns)
  open(unit=103, file=trim(FileName), status="old")
    x=0
    do 
      x = x+1
      read(103, *,IOSTAT=ios)  tmpC
      if (ios > 0) then
        print *, "Unexpected input in file "//trim(FileName)//" on line ", x
        stop
      else if (ios < 0) then
        exit   ! EOF
      else
        do i=1, nInd
          if (Id(i) == tmpC) skip(i) = .FALSE.
        enddo
      endif
    enddo
  close(103)

end subroutine ReadOnlyList

!===============================================================================