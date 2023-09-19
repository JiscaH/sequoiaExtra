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
    print '(a)',    '  --help       print usage information and exit'
    print '(a)',    '  --in         input file with grm; extensions .grm.id and .grm.gz are added'
    print '(a)',    '  --out        output file'
    print '(a)',    '  --only       only consider pairs with one or both individuals listed'
    print '(a)',    '  --lower      export pairs with R value below this threshold'
    print '(a)',    '  --upper      export pairs with R value above this threshold'
    print '(a)',    '  --quiet      hide counter while running and message when done.'
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
  

end module Fun

!===============================================================================

program filter_grm
  use Fun
  implicit none
   
  integer :: x, i,j, z, n, ios, nArg, nrows_grm, print_chunk, p
  double precision :: r, lowr, upr
  character(len=32) :: arg, argOption
  character(len=nchar_filename) :: infile, outfile, OnlyListFileName
  logical :: FileOK, quiet
  
  ! set default values
  lowr = -HUGE(0D0)
  upr = HUGE(0D0)
  infile = 'nofile'
  outfile = 'grm_filter_output.txt'
  OnlyListFileName = 'nofile'
  quiet = .FALSE.
  
  ! read command line arguments: --lower, --upper, --in, --out, --only
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
        
      case ('--lower')
        i = i+1
        call get_command_argument(i, argOption)
        read(argOption, *)  lowr
        
      case ('--upper')
        i = i+1
        call get_command_argument(i, argOption)
        read(argOption, *)  upr
        
      case ('--only')
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
  call EXECUTE_COMMAND_LINE("(gzip -dc  "//trim(infile)//".grm.gz > grmpipe) &")
  
  ! open a read (outflow) connection to the pipe
  open(11, file="grmpipe", action='read')  
   
  ! open a textfile to write selected subset to
  open(42, file=trim(outfile), action='write')  
    write(42, '(2a10, 2X, 2a40, a10, a15)') 'index1', 'index2', 'ID1', 'ID2', 'nSNP', 'R'

    nrows_grm = (nInd * (nInd-1)/2 + nInd)
    n = 0
    print_chunk = roundit(nrows_grm/20D0)  ! print at approx every 5% progress
    p = 0
    
    do x = 1, nrows_grm
      if (.not. quiet .and. MOD(x, print_chunk)==0) then
        p = p +1
        print *, x, '  ', p*5, '%'
      endif
      read(11, *, iostat=ios) i,j,z,r  
      if (ios/=0) exit   ! stop if end of file / incorrect entry
      if (skip(i) .and. skip(j))  cycle
      ! write to outfile entries that meet criteria
     if (r < lowr .or. r > upr) then
       n = n+1
       write(42, '(2i10, 2X, 2a40, i10, e15.6)')  i, j, ID(i), ID(j), z, r
     endif   
    end do
  
  close(42)
  close(11)
  
  if (.not. quiet) then
    write(*,*)  ""
    write(*,*)  "Found ", n, " pairs matching the criteria, written to ", trim(outfile)
    write(*,*)  ""
  endif
  
  
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

  if (allocated(ID))  deallocate(ID)  
  if (allocated(skip))  deallocate(skip)
        
end program filter_grm


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