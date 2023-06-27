! Find likely parent-offspring pairs based on OH count
! While sequoia may run out of memory space for very large datasets,
! this minimalistic program is much less likely too. 
! 
! Input:
! - genotype file
! Command line arguments:
! --geno : genotype filename
! --dup : search for duplicates
! --po  : search for parent-offspring pairs
! --max_dup : maximum SNPs at which duplicates differ
! --max_oh : maximum OH count for likely parent-offspring pairs
! --only : only find potential parents or offspring for those individuals 
! --out : output filename  (_PO.txt & _DUP.txt)
! Output:
! - text file with ID1 - ID2 - OH - SnpdBoth  for pairs with OH <= maxOH
!
! to compile:
! gfortran -std=f95 -fall-intrinsics -O3 find_PO_dups.f90 -o findPO
! to debug:
! gfortran -std=f95 -fall-intrinsics -Wall -pedantic -fbounds-check -g -Og find_PO_dups.f90 -o findPO

!===============================================================================
!===============================================================================

module Global_variables
  implicit none

  integer :: nInd, nSnp, maxOH, maxDUP
  integer, parameter :: nchar_filename = 2000
  logical :: quiet
  character(len=20), allocatable, dimension(:) :: Id
  logical, allocatable, dimension(:) :: skip
  integer, allocatable, dimension(:,:) :: Genos

end module Global_variables

!===============================================================================

module FileDim
  implicit none

  contains
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    integer function FileNumCol(FileName)
      implicit none

      character(len=*), intent(IN) :: FileName
      integer :: j, strLen, numcol
      character(len=5000) :: line

      open(unit=102, file=trim(FileName), status="old")
      read(102, '(a)' ) line
      close(102) 

      strLen = len_trim(line)
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

    end function FileNumCol

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    integer function FileNumRow(FileName)
      implicit none

      character(len=*), intent(IN) :: FileName
      integer :: nrow, i, maxRow, IOerr
      character(len=5000) :: dumC

      maxRow = 5000000  ! fail safe
      nrow = 0
      open(unit=102, file=trim(FileName), status="old")
      do i=1, maxRow
        read(102,*,IOSTAT=IOerr) dumC
        if (IOerr < 0) then
          exit  ! EOF
        else
          nrow = nrow +1  
        end if
      enddo
      close(102)
      FileNumRow = nrow

    end function FileNumRow
  
end module FileDim

!===============================================================================
!===============================================================================

program main
  use Global_variables
  implicit none
  
  integer :: nArg, i, x
  character(len=32) :: arg, argOption
  character(len=nchar_filename) :: GenoFileName, OnlyListFileName, OutFileName
  logical :: FileOK, DoDup, DoPO
  
  ! defaults  ~~~~
  GenoFileName = 'Geno.txt'
  OnlyListFileName = 'NoFile'
  OutFileName = 'Pairs_maybe'
  maxOH = 0
  maxDUP = 0
  DoPO = .FALSE.
  DoDup = .FALSE.
  quiet = .FALSE.

  ! command line argumens  ~~~~
  nArg = command_argument_count()
  if (nArg > 0) then
    i = 0
    do x = 1, nArg
      i = i+1
      if (i > nArg)  exit
      call get_command_argument(i, arg)
      
      select case (arg)
        case ('-h', '--help')
          call print_help()
          stop       
           
        case ('--geno')  
          i = i+1
          call get_command_argument(i, GenoFileName)
         
        case('--po', '--PO')
          DoPO = .TRUE.

        case('--dup', '--DUP')
          DoDup = .TRUE.
          
        case ('--max_diff', '--max_dup', '--maxDUP')
          i = i+1
          call get_command_argument(i, argOption)
          read(argOption, *)  maxDUP
          
        case ('--max_oh', '--maxOH')
          i = i+1
          call get_command_argument(i, argOption)
          read(argOption, *)  maxOH        
         
        case ('--only')
          i = i+1
          call get_command_argument(i, OnlyListFileName)
          
        case ('-o', '--out')
          i = i+1
          call get_command_argument(i, OutFileName)

        case ('--quiet')
          quiet = .TRUE.
          
        case default
            print '(2a, /)', 'Unrecognised command-line option: ', arg
            call print_help()
            stop

      end select
    end do
  endif
  
  
  ! input check  ~~~~
  if (.not. DoDup .and. .not. DoPO) then
    print *, 'Please select --dup and/or --po'
    print *, ''
    call print_help()
    stop
  endif
  
  if (maxOH < 0) then
    print *, 'maxOH must be positive!'
    stop
  endif
  inquire(file=trim(GenoFileName), exist = FileOK)
  if (.not. FileOK) then
    print *, 'File ', trim(GenoFileName), ' not found!'
    stop
  endif
  if (OnlyListFileName /= 'NoFile') then
    inquire(file=trim(OnlyListFileName), exist = FileOK)
    if (.not. FileOK) then
      print *, 'File ', trim(OnlyListFileName), ' not found!'
      stop
    endif
  endif

  
  ! run  ~~~~
  if (.not. quiet)  print *, 'Reading genotype data in '//trim(GenoFileName)//' ...'
  call ReadGeno(GenoFileName)
 
  allocate(skip(nInd))
  skip = .FALSE.
  if (OnlyListFileName /= 'NoFile') then
    if (.not. quiet)  print *, "Reading individuals in --only file "//trim(OnlyListFileName)//" ... " 
    call ReadOnlyList(OnlyListFileName)
  endif
  
  if (DoDup) then
    if (.not. quiet)  print *, 'Searching for duplicates ...'
    call find_dups(OutFileName)
  endif
  
  if (DoPO) then
    if (.not. quiet)  print *, 'Searching for parent-offspring pairs ...'
    call find_PO(OutFileName)
  endif  
  
  
  ! free up allocated variables  ~~~~
  deallocate(Id)
  deallocate(skip)
  deallocate(Genos)
  

  ! help  ~~~~
  contains
    subroutine print_help()
        print '(a)',  'Find potential parent-offspring pairs based on OH count'
        print '(a, /)', 'command-line options:'
        print '(a)',    '  -h, --help          print usage information and exit'
        print '(a)',    '  --geno <filename>   file with genotype data. Default: Geno.txt'
        print '(a)',    '  --dup               search for duplicates'   
        print '(a)',    '  --po                search for potential parent-offspring pairs' 
        print '(a)',    '  --max_oh <n>        maximum OH count for potential PO pairs. Default: 1', &
                        '                        see ?sequoia::MaxMismatch in R' 
        print '(a)',    '  --max_diff <n>      maximum number of differences between duplicate samples'    
        print '(a)',    '  --max_dup <n>       maximum number of differences between duplicate samples'
        print '(a)',    '  --only <filename>   only calculate OH when one or both are in this subset'
        print '(a)',    '  --out <filename>    output file name. Default: Pairs_maybe_PO.txt' 
!        print '(a)',    '  --quiet             suppress all messages'
    end subroutine print_help

end program main


!===============================================================================
!===============================================================================

subroutine find_dups(FileName_part)
  use Global_variables
  implicit none
  
  character(len=nchar_filename), intent(IN) :: FileName_part 
  character(len=nchar_filename+7) :: FileName
  integer :: i,j, Lboth, Diff_ij, nPairs, ID_len, l
  character(len=200) :: HeaderFMT, DataFMT
  
  FileName = trim(FileName_part)//'_DUP.txt'
  
  ID_len = 20
  write(HeaderFMT, '( "(2(a8, 4X, a", I0, ", 4X), 2a10)" )')  ID_len
  write(DataFMT,   '( "(2(i8, 4X, a", I0, ", 4X), 2i10)" )')  ID_len  
  
  nPairs = 0
  open(unit=201, file=trim(FileName), status='unknown')
    write (201, HeaderFMT) 'Row1', 'ID1', 'Row2', 'ID2', 'nDiffer', 'SnpdBoth'
    do i=1, nInd-1
      if (.not. quiet .and. mod(i, 2000)==0)  print *, i
      do j=i+1, nInd
        if (skip(i) .and. skip(j))  cycle
        Lboth = 0
        Diff_ij = 0
        do l=1, nSnp
          if (Genos(l,i)==-1 .or. Genos(l,j)==-1)  cycle
          Lboth = Lboth +1
          if (Genos(l,i) /= Genos(l,j)) then
            Diff_ij = Diff_ij +1
            if (Diff_ij > maxDup)  exit
          endif
        enddo  
        if (Diff_ij > maxDup)  cycle
        if (dble(Diff_ij)/Lboth  > 2.0*dble(MaxDUP)/nSnp)  cycle    
        ! if arrived here, i+j are potential duplicate samples from same individual
        write(201, DataFMT) i, Id(i), j, Id(j), Diff_ij, Lboth
        nPairs = nPairs +1         
      enddo
    enddo
  
  close(201)
  
  write(*,'("Found ", i8, " duplicate pairs using maxDUP=", i4)') nPairs, maxDUP 

end subroutine find_dups


!===============================================================================

subroutine find_PO(FileName_part)
  use Global_variables
  implicit none
  
  character(len=nchar_filename), intent(IN) :: FileName_part
  character(len=nchar_filename+7) :: FileName  
  integer :: i,j, Lboth, OH_ij, nPairs, ID_len, l
  character(len=200) :: HeaderFMT, DataFMT
  
  FileName = trim(FileName_part)//'_PO.txt'
  
  ID_len = 20
  write(HeaderFMT, '( "(2(a8, 4X, a", I0, ", 4X), 2a10)" )')  ID_len
  write(DataFMT,   '( "(2(i8, 4X, a", I0, ", 4X), 2i10)" )')  ID_len
  
  nPairs = 0
  open(unit=201, file=trim(FileName), status='unknown')
    write(201, HeaderFMT) 'Row1', 'ID1', 'Row2', 'ID2', 'OH', 'SnpdBoth'
  
    do i=1, nInd-1
      if (.not. quiet .and. Mod(i, 5000)==0)  print *, i
      do j=i+1, nInd
        if (skip(i) .and. skip(j))  cycle
        Lboth = 0
        OH_ij = 0
        do l=1,nSnp
          if (Genos(l,i)==-1 .or. Genos(l,j)==-1)  cycle
          Lboth = Lboth +1
          if ((Genos(l,i)==0 .and.Genos(l,j)==2) .or. &
           (Genos(l,i)==2 .and. Genos(l,j)==0)) then
            OH_ij = OH_ij +1
            if (OH_ij > maxOH) exit
          endif  
        enddo
        if (OH_ij > maxOH)  cycle
        if (Lboth < nSnp/3.0)   cycle   ! >2/3th of markers missing for one or both 
        if (dble(OH_ij)/Lboth  > 2.0*dble(MaxOH)/nSnp)  cycle        
        ! if arrived here, i+j are potential parent-offspring pair
        write(201, DataFMT) i, Id(i), j, Id(j), OH_ij, Lboth
        nPairs = nPairs +1
          
      enddo
    enddo
  
  close(201)
  
  write(*,'("Found ", i8, " parent-offspring pairs using maxOH=", i4)') nPairs, maxOH 

end subroutine find_PO

!===============================================================================
!===============================================================================

subroutine ReadGeno(GenoFileName)
  use Global_variables
  use FileDim
  implicit none

  character(len=nchar_filename), intent(IN) :: GenoFileName
  integer :: i, l
  integer, allocatable, dimension(:) :: GenosV

  nSnp = FileNumCol(trim(GenoFileName)) -1  ! column 1 = IDs
  nInd = FileNumRow(trim(GenoFileName))  

  allocate(GenosV(nSnp))
  allocate(Genos(nSnp, nInd))   ! transpose: faster
  Genos = -1
  allocate(Id(nInd))
  Id = "NA"

  open (unit=101,file=trim(GenoFileName),status="old")
  do i=1,nInd
    read (101,*)  Id(i), GenosV
    do l=1,nSnp
      if (GenosV(l)>=0) then
        Genos(l,i) = GenosV(l)  
      endif
    enddo
  enddo
  close (101)
  
  deallocate(GenosV)
  
end subroutine ReadGeno

!===============================================================================

subroutine ReadOnlyList(FileName)
  use Global_variables
  use FileDim
  implicit none

  character(len=nchar_filename), intent(IN) :: FileName
  integer :: x, i, nrows, IOerr
  character(len=20) :: tmpC

  nrows = FileNumRow(trim(FileName))  ! no header

  skip = .TRUE.

  ! single column (ignore all other columns)
  open(unit=103, file=trim(FileName), status="old")
    do x=1, nrows
      read(103, *,IOSTAT=IOerr)  tmpC
      if (IOerr > 0) then
        print *, "Wrong input in file "//trim(FileName)//" on line ", i
        stop
      else if (IOerr < 0) then
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
! end. 