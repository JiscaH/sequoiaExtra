! Estimate rate of Transmission Ratio Distortion for each SNP
! among offspring-dam-sire trios, with all three genotyped
! Use chi-square test (?)

! note: TR may appear distorted due to genotyping errors


! ##############################################################################

module global
  implicit none

  integer :: nInd, nSnp, ID_len
  integer, parameter :: nchar_filename = 2000, nchar_ID = 40
  integer,allocatable,dimension(:,:) :: Genos, Parent
  logical :: quiet
  double precision :: Er, E(0:2,0:2,0:2)
  character(len=nchar_ID), allocatable, dimension(:) :: Id 
  
    !=========================
  contains
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! count number of columns in text file
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    integer function FileNumCol(FileName)
      implicit none

      character(len=*), intent(IN) :: FileName
      integer :: j, strLen, numcol
      character(len=500000) :: line

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
            numcol = numcol +1    ! new column starts at j+1
          endif
        endif
      enddo
      FileNumCol = numcol

    end function FileNumCol
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! count number of rows in text file
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
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

end module global


! ##############################################################################
! ##  Main program   ##
! ##############################################################################

program TRD
  use Global
  implicit none

  ! input
  integer :: x, i, nArg
  character(len=32) :: arg, argOption
  character(len=nchar_filename) :: PedFileName, GenoFileName, OutFileName
  logical :: FileOK
  ! output

  ! set default values
  PedFileName = 'Pedigree.txt'
  GenoFileName = 'Geno.txt'
  quiet = .FALSE.
  Er = 1e-4
  
  
  ! read arguments from command line
  nArg = command_argument_count()

  if (nArg == 0) then
    if (.not. quiet)  print *, 'Please specify input files!'
    call print_help()
    stop
    
  else 

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
          
        case ('--pedigreeIN')    
          i = i+1
          call get_command_argument(i, PedFileName)    

        case ('--err')
          i = i+1
          call get_command_argument(i, argOption)
          read(argOption, *)  Er   ! TODO: length 3 vector
            
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
    enddo

  endif

  !=========================
  if (Er <= 0D0 .or. Er > 0.5) then
    write(*,*) "Er must be a value between 0 and 0.5"
    stop
  endif
  
  ! check if input files exist
  inquire(file=trim(PedFileName), exist = FileOK)
  if (.not. FileOK) then
    write(*,*)  "Input file ", trim(PedFileName), " not found"
    stop
  endif

  inquire(file=trim(GenoFileName), exist = FileOK)
  if (.not. FileOK) then
    write(*,*)  "Genotype file ", trim(GenoFileName), " not found"
    stop
  endif
  
  !=========================
  ! Do stuff
  
  if (.not. quiet)  print *, "Reading genotype data in "//trim(GenoFileName)//" ... "
  call ReadGeno(GenoFileName) 
  if (.not. quiet)  print *, "Read ", nInd, " Individuals and ", nSnp, " SNPs"
  
  if (.not. quiet)  print *, "Reading pedigree in "//trim(PedFileName)//" ... "
  call ReadPedFile(PedFileName)
  if (.not. quiet)  print *, "Read ", COUNT(Parent(1,:)>0 .and. Parent(2,:)>0), " id-dam-sire trios "
  
  if (.not. quiet)  print *, "Counting genotype combinations ... "
  call prep_E()
  call Calc_TRD(OutFileName)
  if (.not. quiet)  print *, "Output written to ", trim(OutFileName)
  
  call deallocall()
  
  
  
  !=========================
  contains
    subroutine print_help()
        print *, ''
        print '(a)',    '<<< Count transmission ratios >>>'
        print *, ''
        print '(a, /)', 'command-line options:'
        print '(a)',    '  -h, --help          print usage information and exit'
        print '(a)',    '  --pedigreeIN <filename>  input file with id, dam, sire'
        print '(a)',    '  --geno <filename>   input file with genotype data. Default: Geno.txt'  
        print '(a)',    '  --err               presumed genotyping error rate; default: 0.0001'        
        print '(a)',    '  --out <filename>    output file '
        print '(a)',    '  --quiet             suppress all messages'
    end subroutine print_help 
    
end program TRD

    
! ##############################################################################
! ##   read in data   ##
! ##############################################################################


subroutine ReadGeno(GenoFileName)
  use Global
  implicit none

  character(len=nchar_filename), intent(IN) :: GenoFileName
  integer :: i, l
  integer, allocatable, dimension(:) :: GenosV
  character(len=3) :: maxchar_ID
  character(len=nchar_ID) :: IDx
    
  nSnp = FileNumCol(trim(GenoFileName)) -1  ! column 1 = IDs
  nInd = FileNumRow(trim(GenoFileName))   

  allocate(GenosV(nSnp))
  allocate(Genos(nSnp, nInd))   ! transpose: faster
  Genos = -1
  allocate(Id(nInd))
  Id = "NA"

  ID_len = 5

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
    print *, "ERROR! Max length for IDs is "//maxchar_ID//" characters"
    stop
  endif

end subroutine ReadGeno

!===============================================================================

subroutine ReadPedFile(FileName)
  use Global
  implicit none

  character(len=*), intent(IN) :: FileName
  integer :: i, j, k, x, IOerr, nIndP
  character(len=nchar_ID) :: tmpC(3)
  character(len=nchar_ID), allocatable, dimension(:,:) :: NamePed

  nIndP = FileNumRow(trim(FileName)) -1  ! 1st row = header
  allocate(NamePed(3, nIndP))
  NamePed = 'NA'

  open(unit=103, file=trim(FileName), status="old")
  read(103,*)    ! header              
  do i=1,nIndP
    read(103, *,IOSTAT=IOerr)  tmpC
    if (IOerr > 0) then
      print *, "Wrong input in file "//trim(FileName)//" on line ", i
      stop
    else if (IOerr < 0) then
      exit   ! EOF
    else
      NamePed(:,i) = tmpC 
    end if
  enddo
  close(103)

  ! Pedigree names to row numbers in genotype file
  allocate(Parent(2, nInd))
  Parent = 0
  do i = 1, nIndP
    j = 0
    do x= 1,nInd
      if (NamePed(1,i) == Id(x)) then
        j = x
        exit
      endif
    enddo
    if (j==0)  cycle
  
    do k = 1,2
      do x=1, nInd
        if (NamePed(k+1,i) == Id(x)) then
          Parent(k,j) = x
          exit
        endif
      enddo
    enddo
  enddo
  
  if (allocated(NamePed)) deallocate(NamePed)
 
end subroutine ReadPedFile


! ##############################################################################
! ##   Count observed genotype combo's   ##
! ##############################################################################


subroutine Calc_TRD(FileName)
  use Global
  implicit none
  
  character(len=*), intent(IN) :: FileName
  double precision :: Chi2(nSnp), O_prop(0:2,0:2,0:2)
  integer :: i,j,h,l, hdr_A(0:2,0:2,0:2)
  integer :: obs_count(-1:2, -1:2, -1:2, nSnp), nObs(nSnp), O(-1:2,-1:2,-1:2)
  character(len=4) :: header(3*3*3)

  
  ! observed count
  obs_count = 0
  nObs = 0
  do l=1, nSnp
    if (MOD(l, 2000)==0)  print *, l 
    O = 0    
    do i = 1, nInd
      j = Parent(1,i)
      h = Parent(2,i)
      if (j==0 .or. h==0)  cycle
      O(Genos(l,i), Genos(l,j), Genos(l,h)) =  O(Genos(l,i), Genos(l,j), Genos(l,h)) +1
    enddo
    obs_count(:,:,:, l) =  O
    nObs(l) = SUM(O(0:2, 0:2, 0:2))  ! count of id-dam-sire all three genotyped
    ! Calc Chi-square
    O_prop = O(0:2, 0:2, 0:2) / dble(nObs(l))
    Chi2(l) = SUM( (O_prop - E)**2 / E )   ! division by zero ...
!    Chi2(l) = SUM( (O_prop - E)**2 )
  enddo
  
  ! Write array to a text file; do further analysis in R (incl P values)
  
  ! make header
  do h=0,2
    do j=0,2
      do i=0,2
        hdr_A(i,j,h) = 100*i + 10*j + h
      enddo
    enddo
  enddo  
  write(header, '("G",i3.3)')  RESHAPE(hdr_A, (/3*3*3/)) 
   
  open(unit=42, file=trim(FileName), status='unknown')
    write(42, '(2a6, a10, 50a6)')  'SNP', 'Nobs', 'Chi2', header
    do l=1, nSnp
      write (42,'(2i6, f10.3, 50i6)')  l, nObs(l), Chi2(l), RESHAPE(obs_count(0:2, 0:2, 0:2, l), (/3*3*3/)) 
    enddo
  close(42)
  
end subroutine Calc_TRD


! ##############################################################################

subroutine prep_E
  use Global
  implicit none
  ! calc expected proportions of observed genotype combinations, with genotyping errors,
  ! (else division by zero issues), conditional on the trio being id-dam-sire
  
  double precision :: AKA2P(0:2,0:2,0:2), OcA(-1:2,0:2), q, AHWE(0:2), &
    AcO(0:2, -1:2), Tmp(0:2, 0:2, 0:2)
  integer :: i,j,h, x,y,z
   
  ! inheritance conditional on both parents  (actual genotypes)
  AKA2P(0,0,:) = dble((/ 1.0, 0.5, 0.0 /))
  AKA2P(0,1,:) = dble((/ 0.5, 0.25, 0.0 /))
  AKA2P(0,2,:) = dble((/ 0.0, 0.0, 0.0 /))

  AKA2P(1,0,:) = dble((/ 0.0, 0.5, 1.0 /))
  AKA2P(1,1,:) = dble((/ 0.5, 0.5, 0.5 /))
  AKA2P(1,2,:) = dble((/ 1.0, 0.5, 0.0 /))

  AKA2P(2,0,:) = dble((/ 0.0, 0.0, 0.0 /))
  AKA2P(2,1,:) = dble((/ 0.0, 0.25, 0.5 /))
  AKA2P(2,2,:) = dble((/ 0.0, 0.5, 1.0 /))
  
  ! Prob. observed (rows) conditional on actual (columns) (ErrFlavour' = 2.0)
  OcA(-1,:) = 1.0D0      ! missing 
  OcA(0:2, 0) = (/ (1-Er/2)**2, Er*(1-Er/2), (Er/2)**2 /)   ! act=0
  OcA(0:2, 1) = (/ Er/2, 1-Er, Er/2 /)                      ! act=1
  OcA(0:2, 2) = (/ (Er/2)**2, Er*(1-Er/2),  (1-Er/2)**2 /)  ! act=2
 
  ! For parents, calc P(act|obs) = P(obs|act) * P(act) / P(obs)
  ! P(act) = under HWE, P(obs) = scaling factor
  ! for simplicity, assume all SNPs have allele frequency 0.5
  ! (TODO later: read in allele frequencies & see if it matters)
  q = 0.5
  AHWE = (/ (1-q)**2, 2*q*(1-q) , q**2 /)
 
  AcO = 0D0
  do y=-1,2  ! observed genotype
    do x=0,2  ! actual genotype    
      AcO(x,y) = OcA(y,x) * AHWE(x)
    enddo
    AcO(:,y) = AcO(:,y) / SUM(AcO(:,y))   ! scale to sum to 1
  enddo   
 
  E = 0D0
  do h=0,2  ! sire
    do j=0,2  ! dam
      do i=0,2  ! offspring obs
        do z=0,2  ! sire act
          do y=0,2  ! dam act
            do x=0,2  ! offspring act
              Tmp(x,y,z) = OcA(i,x) * AKA2P(x,y,z) * AcO(y,j) * AcO(z,h)
            enddo
          enddo
        enddo
        E(i,j,h) = SUM(Tmp)
      enddo
    enddo
  enddo

end subroutine prep_E

! ##############################################################################

subroutine deallocall
  use Global
  implicit none

  if (allocated(Genos)) deallocate(Genos)
  if (allocated(Parent)) deallocate(Parent) 
  if (allocated(Id)) deallocate(Id)

end subroutine deallocall