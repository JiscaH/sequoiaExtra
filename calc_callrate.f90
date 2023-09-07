! input: 
! - filename of genotype file in sequoia format
! - filename of text file with IDs (optional)
!
! output:
! - text file with 3 columns: rownum in geno file, id, callrate

!===============================================================================

module CR
  implicit none
  
  integer :: nInd, nSnp, ID_len
  integer, parameter :: nchar_filename = 2000, nchar_ID = 40
  integer,allocatable,dimension(:,:) :: Genos
  logical, allocatable, dimension(:) :: skip
  character(len=nchar_ID), allocatable, dimension(:) :: Id
  
   
  contains
  
  integer function FileNumCol(FileName)
    implicit none

    character(len=*), intent(IN) :: FileName
    integer :: j, strLen, numcol
    character(len=5000) :: line

    open(unit=102, file=trim(FileName), status="old")
    read(102, '(a)' ) line   ! skip header (more likely to be >5k chars)                                                                    
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
          numcol = numcol +1    ! new column starts at j+1
        endif
      endif
    enddo
    FileNumCol = numcol
  end function FileNumCol

  
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
   
end module CR

!===============================================================================

program CallRate
  use CR
  implicit none
  
  integer :: i, x, nArg
  character(len=32) :: arg
  character(len=nchar_filename) :: GenoFileName, OnlyListFileName, OutFileName
    
  GenoFileName = 'Geno.txt'
  OnlyListFileName = 'NoFile'
  OutFileName = 'indiv_callrates.txt'
  
  ! read command line arguments
  nArg = command_argument_count()
  
  if (nArg > 0) then
    i = 0
    do x = 1, nArg
      i = i+1
      if (i > nArg)  exit
      call get_command_argument(i, arg)
      
      select case (arg)
      case ('--geno')  
        i = i+1
        call get_command_argument(i, GenoFileName)
        
      case ('--only')
        i = i+1
        call get_command_argument(i, OnlyListFileName)
          
      case ('--out')
        i = i+1
        call get_command_argument(i, OutFileName)    
          
      case default
        print '(2a, /)', 'Unrecognised command-line option: ', arg
        stop

      end select
      
    enddo
  endif
  
  call ReadGeno(GenoFileName)
  print *, 'nInd: ', nInd, '  nSnp: ', nSnp
  
  allocate(skip(nInd))  
  if (OnlyListFileName /= 'NoFile') then
    call ReadOnlyList(OnlyListFileName)
  else 
    skip = .FALSE.
  endif
   
  call calc_CR(OutFileName)
  
  call deallocall()
  
  print *, 'done'
  
end program CallRate

!===============================================================================

subroutine calc_CR(OutFileName)
  use CR
  implicit none
  
   character(len=nchar_filename), intent(IN) :: OutFileName
  integer :: i, nScored(nInd)
  double precision :: CRi(nInd)
  
  nScored = -1
  CRi = -1D0
  do i = 1, nInd
    if (skip(i))  cycle
    nScored(i) = COUNT(Genos(:,i) >= 0)
  enddo
  CRi = nScored / dble(nSnp)

  ! write to file
  open(unit=101, file = trim(OutFileName), status='unknown')
    write(101, '(a6, 2X, a20, a9)')  'row', 'id', 'callrate'
    do i=1, nInd
      if (skip(i))  cycle
      write(101, '(i6, 2X, a20, f8.3)')  i, ID(i), CRi(i)
    enddo
  close(101)

end subroutine calc_CR

!===============================================================================

subroutine ReadGeno(GenoFileName)
  use CR
  implicit none

  character(len=nchar_filename), intent(IN) :: GenoFileName
  integer :: i, l
  integer, allocatable, dimension(:) :: GenosV
  character(len=3) :: maxchar_ID
  character(len=nchar_ID) :: IDx

  nInd = FileNumRow(trim(GenoFileName))
  nSnp = FileNumCol(trim(GenoFileName)) -1  ! column 1 = IDs    

  allocate(GenosV(nSnp))
  allocate(Genos(nSnp, nInd))   ! transpose: faster
  Genos = -1
  allocate(Id(0:nInd))
  Id = "NA"

  ID_len = 3

  open (unit=101,file=trim(GenoFileName),status="old")
    do i=1,nInd
      read (101,*)  IDx, GenosV
      if (ANY(Id == IDx)) then
        print *, "ERROR! IDs in genotype file must be unique"
        stop
      endif
      Id(i) = IDx
      if (LEN_TRIM(Id(i)) > ID_len)  ID_len = LEN_TRIM(Id(i))
      do l=1,nSnp
        if (GenosV(l)>=0 .and. GenosV(l)<=2) then
          Genos(l,i) = GenosV(l)  
        endif
      enddo
    enddo
  close (101)
  deallocate(GenosV)

  if (ID_len > nchar_ID) then
    write(maxchar_ID, '(i3)')  nchar_ID
    print *, "Max length for IDs is "//maxchar_ID//" characters", .FALSE.
    stop
  endif

end subroutine ReadGeno


!===============================================================================

subroutine ReadOnlyList(FileName)
use CR
implicit none

character(len=nchar_filename), intent(IN) :: FileName
integer :: x, i, nrows, IOerr
character(len=nchar_ID) :: tmpC

nrows = FileNumRow(trim(FileName))  ! no header

skip = .TRUE.

! single column (ignore all other columns)
print *, "Reading individuals in --only file "//trim(FileName)//" ... "
open(unit=103, file=trim(FileName), status="old")
  do x=1, nrows
    read(103, *,IOSTAT=IOerr)  tmpC
    if (IOerr > 0) then
      print *, "Wrong input in file "//trim(FileName)//" on line ", x
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

subroutine deallocall
  use CR
  implicit none

  if (allocated(skip)) deallocate(skip)
  if (allocated(Genos)) deallocate(Genos)
  if (allocated(Id)) deallocate(Id)

end subroutine deallocall

!===============================================================================