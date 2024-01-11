! check if the pedigree is compatible with the genetic data
! pedigree: animal_id, dam_id, sire_id, mat_grandsire_id
! genotype: 1 row per individual, 1st column is ID, then 1 column per SNP; 0/1/2 copies of minor allele, -9=missing
! other parameters: runtime options (a.o. genotyping error)

! This is an extremely minimalistic version of sequoia, 
! useful for finding errors even in very large pedigrees

! compile: gfortran -std=f95 -fall-intrinsics -O3 pedigree_checker.f90 -o PedChecker
! gfortran -std=f95 -fall-intrinsics -Wall -pedantic -fbounds-check -g -Og pedigree_checker.f90 -o PedChecker

! ##############################################################################
! ##  Global variables  ##
! ##############################################################################

module global
  implicit none

  integer :: nInd, nSnp, nTrios, ID_len, nRel 
  integer,allocatable,dimension(:,:) :: Geno, Trios
  logical :: quiet
  double precision, parameter :: Missing = 999D0
  double precision, allocatable, dimension(:,:,:,:,:,:) :: LL_trio
  
  contains
    subroutine timestamp(spaceIN)
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

end module global

!===============================================================================

module FileIO
  use Global
  implicit none
  
  integer, parameter :: nchar_filename = 2000, nchar_ID = 40
  character(len=nchar_ID), allocatable, dimension(:) :: Id 
  character(len=nchar_ID), allocatable, dimension(:,:) :: TrioNames
  character(len=7) :: Ped_header(3)

  contains
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    integer function FileNumCol(FileName)
      character(len=*), intent(IN) :: FileName
      integer :: j, strLen, numcol
      character(:), allocatable :: line
      
      allocate(character(500000) :: line)

      open(unit=102, file=trim(FileName), status="old")
      read(102, '(a)' ) line
      close(102) 

      strLen = len_trim(line)
      if (strLen >= 500000)  print *, 'WARNING: EXCEEDING MAXIMUM NUMBER OF SNPs!'
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

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    integer function FileNumRow(FileName)
      character(len=*), intent(IN) :: FileName
      integer :: nrow, i, maxRow, IOerr
      character(len=42) :: dumC

      maxRow = 1e6  ! fail safe
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
  
end module FileIO

!===============================================================================

module Probs
  implicit none
  private :: getAF

  double precision :: OcA(-1:2,3), OKA2P(-1:2,3,3), AKA2P(3,3,3)
  double precision, allocatable, dimension(:,:) :: AHWE, OHWE
  double precision, allocatable, dimension(:,:,:) :: AKAP, OKAP
  
  contains
    subroutine PrecalcProbs(Er, AF_FileName)
    use Global, ONLY: nSnp

    double precision, intent(IN) :: Er
    character(len=*), intent(IN) :: AF_FileName
    integer :: h,i,j,k,l
    double precision ::  AF(nSnp), Tmp(3)

    ! allele frequencies
    AF = getAF(AF_FileName)
    

    !###################
    ! probabilities actual genotypes under HWE
    allocate(AHWE(3,nSnp))
    do l=1,nSnp
      AHWE(1,l)=(1 - AF(l))**2 
      AHWE(2,l)=2*AF(l)*(1-AF(l)) 
      AHWE(3,l)=AF(l)**2 
    enddo


    ! Prob. observed (rows) conditional on actual (columns) 
    ! 'ErrFlavour' = 2.0
    OcA(-1,:) = 1.0D0      ! missing 
    OcA(0:2, 1) = (/ (1-Er/2)**2, Er*(1-Er/2), (Er/2)**2 /)   ! act=0
    OcA(0:2, 2) = (/ Er/2, 1-Er, Er/2 /)                      ! act=1
    OcA(0:2, 3) = (/ (Er/2)**2, Er*(1-Er/2),  (1-Er/2)**2 /)  ! act=2


    ! probabilities observed genotypes under HWE  + genotyping error pattern
    allocate(OHWE(-1:2,nSnp))
    do l=1,nSnp
      do i=-1,2    ! obs
        OHWE(i,l) = SUM( OcA(i,:) * AHWE(:, l) )
      enddo
    enddo


    ! ########################
    ! inheritance conditional on 1 parent
    allocate(AKAP(3,3,nSnp))
    allocate(OKAP(-1:2,3,nSnp))

    do l=1,nSnp
      AKAP(1, :, l) = (/ 1-AF(l), (1-AF(l))/2, 0.0D0 /)
      AKAP(2, :, l) = (/ AF(l), 0.5D0, 1-AF(l) /)
      AKAP(3, :, l) = (/ 0D0, AF(l)/2, AF(l) /)
    enddo

    do l=1,nSnp
      do i=-1,2  ! obs offspring
        do j=1,3    ! act parent
          Tmp=0D0
          do k=1,3    ! act offspring
            Tmp(k) = OcA(i,k) * AKAP(k,j,l)
          enddo
          OKAP(i,j,l) = SUM(Tmp)
        enddo
      enddo
    enddo


    ! #########################
    ! inheritance conditional on both parents

    AKA2P(1,1,:) = dble((/ 1.0, 0.5, 0.0 /))
    AKA2P(1,2,:) = dble((/ 0.5, 0.25, 0.0 /))
    AKA2P(1,3,:) = dble((/ 0.0, 0.0, 0.0 /))

    AKA2P(2,1,:) = dble((/ 0.0, 0.5, 1.0 /))
    AKA2P(2,2,:) = dble((/ 0.5, 0.5, 0.5 /))
    AKA2P(2,3,:) = dble((/ 1.0, 0.5, 0.0 /))

    AKA2P(3,1,:) = dble((/ 0.0, 0.0, 0.0 /))
    AKA2P(3,2,:) = dble((/ 0.0, 0.25, 0.5 /))
    AKA2P(3,3,:) = dble((/ 0.0, 0.5, 1.0 /))

    do i=-1,2  ! obs offspring
      do j=1,3    ! act parent 1
        do h=1,3    !act parent 2
          Tmp=0D0
          do k=1,3    ! act offspring
            Tmp(k) = OcA(i,k) * AKA2P(k,j,h) 
          enddo
          OKA2P(i,j,h) = SUM(Tmp)
        enddo
      enddo
    enddo

  end subroutine PrecalcProbs
  
  !-----------------------------------
  
  function getAF(FileName)
    use Global, ONLY: nSnp, Geno
    use FileIO
    
    character(len=*), intent(IN) :: FileName
    double precision :: getAF(nSnp)
    integer :: l, nCol, nRow, AFcol, k, IOerr
    character(len=50), allocatable :: header(:), tmpC(:)
    double precision :: tmpD, AF(nSnp)
    
    AF = 1D0
    
    if (FileName == 'NoFile') then
    
      do l=1,nSnp
        if (ANY(Geno(l,:)/=-1)) then
          AF(l)=float(SUM(Geno(l,:), MASK=Geno(l,:)/=-1))/(COUNT(Geno(l,:)/=-1)*2)
        endif
      enddo
    
    else
      if (.not. quiet)  print *, "Reading allele frequencies in "//trim(FileName)//" ... "
    
      nCol = FileNumCol(trim(FileName))
      nRow = FileNumRow(trim(FileName))
      if ((nCol==1 .and. nRow /= nSnp) .or. (nCol>1 .and. nRow /= (nSnp+1))) then
        print *, "MAF file "//trim(FileName)//" has different number of SNPs than genotype file!"
        stop
      endif
      allocate(header(nCol))
      header = 'NA'
      AFcol = 0
      
      open(unit=103, file=trim(FileName), status="old")     
        if (nCol == 1) then
          AFcol = 1
        else
          read(103,*)  header
          do k=1, nCol
            if (header(k) == 'MAF' .or. header(k)=='AF' .or. header(k)=='Frequency') then
              AFcol = k
            endif
          enddo
        endif
        if (AFcol > 1)  allocate(tmpC(AFcol -1))
        
        do l=1, nSnp
          if (AFcol == 1) then
            read(103, *,IOSTAT=IOerr)  tmpD
          else
            read(103, *,IOSTAT=IOerr)  tmpC, tmpD
          endif
          if (IOerr > 0) then
            print *, "Wrong input in file "//trim(FileName)//" on line ", l
            stop
          else if (IOerr < 0) then  
            exit  ! EOF
          else
            AF(l) = tmpD
          end if
        enddo
      close(103)
    endif
    
    getAF = AF

  end function getAF

end module Probs


! ##############################################################################
! ##  Main program   ##
! ##############################################################################

program pedigree_checker
  use Global
  use FileIO
  use Probs, ONLY: PrecalcProbs
  implicit none

  ! input
  integer :: x, i, nArg
  double precision :: Er
  character(len=32) :: arg, argOption
  character(len=nchar_filename) :: PedFileName, GenoFileName, OutFileName, AF_FileName
  logical :: CalcProbs, FileOK
  ! output
  integer, allocatable, dimension(:,:) :: OppHom
  double precision, allocatable, dimension(:,:,:) :: LL_out_array
  
  write(*,*) ''
  write(*,*) REPEAT('~',32)
  write(*,*) '   >>>  PEDIGREE CHECKER  <<<'
  write(*,*) REPEAT('~',32)
  write(*,*) '||  Jisca Huisman  |  https://github.com/JiscaH  |  GNU GPLv3  ||'
  write(*,*) ''
  
    
  ! set default values
  PedFileName = 'Pedigree.txt'
  GenoFileName = 'Geno.txt'
  OutFileName = 'Pedigree_OUT.txt'
  AF_FileName = 'NoFile'
  Er = 0.005  
  CalcProbs = .TRUE.   ! transform log10-likelihoods into probabilities
  nRel = 5  ! number of relationships: PO, GP/FA/HS, HA/3rd, UU, (FS)
  quiet = .FALSE.


  ! read arguments from command line
  nArg = command_argument_count()

  if (nArg == 0) then
    if (.not. quiet)  call printt('Using default values, see --help')
    
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
          
        case ('--pedigreeIN', '--trios')    
          i = i+1
          call get_command_argument(i, PedFileName)        
            
        case ('-o', '--out')
          i = i+1
          call get_command_argument(i, OutFileName)
          
        case ('--err')
          i = i+1
          call get_command_argument(i, argOption)
          read(argOption, *)  Er   ! TODO: 1 or 3 
          
        case ('--noFS')
          nRel = 4
          
        case ('--LLR')
          CalcProbs = .FALSE.
          
        case ('--af', '--maf', '--freq')
          i = i+1
          call get_command_argument(i, AF_FileName)         
          
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


  ! read in data & general prep
  ! printt: timestamp() + print *,
  if (.not. quiet)  call printt("Reading genotype data in "//trim(GenoFileName)//" ... ")
  call ReadGeno(GenoFileName)
  
  call PrecalcProbs(Er, AF_FileName) 
  
  if (.not. quiet)  call printt("Reading trios in "//trim(PedFileName)//" ... ")
  call ReadPedFile(PedFileName)
  if (.not. quiet) then
    call timestamp()
    print *, "Read ", nTrios, " trios. "
  endif

  ! count number of opposing homozygous loci for each id-dam and id-sire pair, and id-dam-sire trio
  if (.not. quiet)  call printt("Counting number of opposing homozygous loci ... ")
  allocate(OppHom(3, nTrios))
  call CalcOppHom(OppHom)

  ! calculate log-likelihoods
  if (.not. quiet)  call printt("Calculating log-likelihoods ... ")
  allocate(LL_out_array(nRel,nRel,nTrios))
  call Precalc_trioLLs()
  call CalcLL(LL_out_array)

  ! write to file
  if (.not. quiet)  call printt("Writing output to file ... ")
  call writeped(LL_out_array, OppHom, CalcProbs, OutFileName)
  call DeAllocAll()
  if (.not. quiet)  call printt("")
  if (.not. quiet)  call printt("Done.")

  
  !=========================
  contains
    subroutine print_help()
        print '(a)',    'Calculate relationship probabilities of candidate parents'
        print '(a, /)', 'command-line options:'
        print '(a)',    '  -h, --help          print usage information and exit'
        print '(a)',    '  --trios <filename>  input file with e.g. pedigree: animal_id, dam_id,',&
                        '                       sire_id, or animal_id, matGS_id, sire_id;',&
                        '                       each animal_id may occur multiple times.',&
                        '                       Default: Pedigree.txt'
        print '(a)',    '  --pedigreeIN <filename>    same as --trios'             
        print '(a)',    '  --geno <filename>   input file with genotype data. Default: Geno.txt'       
        print '(a)',    '  --err               presumed genotyping error rate; default: 0.005'
        print '(a)',    '  --LLR               do not transform LLRs to probabilities'
        print '(a)',    '  --noFS              assume pairs cannot be full siblings'
        print '(a)',    '  --out <filename>    output file with pedigree + OH counts + Likelihoods;',&
                        '                       default: Pedigree_OUT.txt'
        print '(a)',    '  --af <filename>     optional input file with allele frequencies. Either',&
                        '                        1 column and no header, or multiple columns with column',&
                        '                        MAF, AF, or Frequency. E.g. output from plink --freq.'                       
        print '(a)',    '  --quiet             suppress all messages'
    end subroutine print_help 

end program pedigree_checker



! ##############################################################################
! ##   read in data   ##
! ##############################################################################

subroutine ReadGeno(GenoFileName)
  use Global
  use FileIO
  implicit none

  character(len=nchar_filename), intent(IN) :: GenoFileName
  integer :: i, l
  integer, allocatable, dimension(:) :: GenoV
  character(len=3) :: maxchar_ID
  character(len=nchar_ID) :: IDx
  
  
  nSnp = FileNumCol(trim(GenoFileName)) -1  ! column 1 = IDs
  nInd = FileNumRow(trim(GenoFileName))   

  allocate(GenoV(nSnp))
  allocate(Geno(nSnp, nInd))   ! transpose: faster
  Geno = -1
  allocate(Id(nInd))
  Id = "NA"

  ID_len = 5

  open (unit=101,file=trim(GenoFileName),status="old")
    do i=1,nInd
      read (101,*)  IDx, GenoV
      if (ANY(Id == IDx)) then
        print *, "ERROR! IDs in genotype file must be unique"
        stop
      endif
      Id(i) = IDx
      if (LEN_TRIM(Id(i)) > ID_len)  ID_len = LEN_TRIM(Id(i))
      do l=1,nSnp
        if (GenoV(l)>=0 .and. GenoV(l)<=2) then
          Geno(l,i) = GenoV(l)  
        endif
      enddo 
    enddo
  close (101)
  
  deallocate(GenoV)

  if (ID_len > nchar_ID) then
    write(maxchar_ID, '(i3)')  nchar_ID
    print *, "ERROR! Max length for IDs is "//maxchar_ID//" characters"
    stop
  endif

end subroutine ReadGeno

!===============================================================================

subroutine ReadPedFile(FileName)
  use Global
  use FileIO
  implicit none

  character(len=*), intent(IN) :: FileName
  integer :: i, j, k, IOerr
  character(len=nchar_ID) :: tmpC(3)

  nTrios = FileNumRow(trim(FileName)) -1  ! 1st row = header
  
  allocate(TrioNames(3,nTrios))   ! typically offspring - dam - sire, but not necessarily
  TrioNames = "NA"

  open(unit=103, file=trim(FileName), status="old")
  read(103,*)  Ped_header                         
  do i=1,nTrios
    read(103, *,IOSTAT=IOerr)  tmpC
    if (IOerr > 0) then
      print *, "Wrong input in file "//trim(FileName)//" on line ", i
      stop
    else if (IOerr < 0) then
      exit   ! EOF
    else
      TrioNames(:,i) = tmpC 
    end if
  enddo
  close(103)


  ! Trios names to row numbers in genotype file
  allocate(Trios(3,nTrios))  ! row numbers in genotype data
  Trios = 0
  do i = 1, nTrios
    do k = 1,3
      if (TrioNames(k,i)=='NA')  cycle
      do j = 1, nInd
        if (TrioNames(k,i) == Id(j)) then
          Trios(k,i) = j
          exit  ! break from inner loop
        endif
      enddo
    enddo
  enddo
  
  
end subroutine ReadPedFile



! ##############################################################################
! ##   General prep   ##
! ##############################################################################



!===============================================================================



! ##############################################################################
! ##   Count OH   ##
! ##############################################################################

subroutine CalcOppHom(OppHomDF)
  use Global
  implicit none

  integer, intent(OUT) :: OppHomDF(3, nTrios)
  integer :: i, k

  OppHomDF = -9
  do i=1,nTrios
    if (Trios(1,i)==0)  cycle
    do k=1,2
      if (Trios(k+1,i) == 0)  cycle
      call calcOH(Trios(1,i), Trios(k+1,i), OppHomDF(k,i))
    enddo
    if (Trios(2,i)/=0 .and. Trios(3,i)/=0) then
      call CalcTrioErr(Trios(1,i), Trios(:,i), OppHomDF(3,i))
    endif
  enddo

end subroutine CalcOppHom

!===============================================================================

subroutine CalcOH(A,B,OH)
use Global
implicit none

integer, intent(IN) :: A, B
integer, intent(OUT) :: OH
integer :: l

OH = 0
do l=1,nSnp
  if ((Geno(l,A)==0 .and.Geno(l,B)==2) .or. (Geno(l,A)==2 .and. Geno(l,B)==0)) then
    OH = OH+1
  endif                       
enddo

end subroutine CalcOH

!===============================================================================

subroutine CalcTrioErr(A,Par, ME)   ! Mendelian errors in offspring-parent-parent trios
use Global
implicit none

integer, intent(IN) :: A, Par(2)
integer, intent(OUT) :: ME
integer :: l, k, Ecnt(3,3,3)   ! offspring - dam - sire

Ecnt(:,1,1) = (/ 0, 1, 2 /)
Ecnt(:,1,2) = (/ 0, 0, 1 /)
Ecnt(:,1,3) = (/ 1, 0, 1 /)

Ecnt(:,2,1) = (/ 0, 0, 1 /)
Ecnt(:,2,2) = (/ 0, 0, 0 /)
Ecnt(:,2,3) = (/ 1, 0, 0 /)

Ecnt(:,3,1) = (/ 1, 0, 1 /)
Ecnt(:,3,2) = (/ 1, 0, 0 /)
Ecnt(:,3,3) = (/ 2, 1, 0 /)

ME = 0
do l=1,nSnp
  if (Geno(l,A)==-1 .or. ALL(Geno(l, Par)==-1)) then
    cycle
  else if (ANY(Geno(l, Par)==-1)) then
    do k=1,2
      if (Geno(l, Par(k))==-1) then
        if (((Geno(l,A)==0).and.(Geno(l,Par(3-k))==2)) .or. &
         ((Geno(l,A)==2).and.(Geno(l,Par(3-k))==0))) then
          ME = ME +1
          cycle
        endif
      endif
    enddo
  else
    ME = ME + Ecnt(Geno(l,A)+1, Geno(l, Par(1))+1, Geno(l, Par(2))+1)
  endif
enddo

end subroutine CalcTrioErr

! ##############################################################################
! ##   Calculate likelihoods   ##
! ##############################################################################

subroutine Precalc_trioLLs
  use Global
  use Probs
  implicit none

  ! calculate likelihoods for PO+PO, PO+FA, PO+HA, PO+U, etc
  ! store in look-up table: G_off (0/1/2/-9) + G_dam (P0/P1/P2) + G_sire (P0/P1/P2)

  ! when not conditioning on parents of both individuals, LL(FA)=LL(GP)=LL(HS),
  ! and similarly all 3rd degree relationships have same LL (log likelihood)

  ! assume that the two parents with genotypes y & z are unrelated. 

  integer :: l, x, y, z, w, v, u
  double precision :: Tmp(3), Tmp2(3,3), Tmp3(3,3,3)

  allocate(LL_trio(-1:2,3,3, nSnp, nRel,nRel))  ! G_Off, G_dam, G_sire, snp, rel_1, rel_2 (PO/GP/HA/U)
  LL_trio = Missing


  do l = 1, nSnp
    do x=-1,2  ! observed genotype offspring
      do y=1,3    ! actual genotype parent 1
        do z=1,3    ! actual genotype parent 2
          ! PO + PO
          LL_trio(x,y,z,l, 1,1) = OKA2P(x,y,z)
          
          ! PO + GP = PO + FA
          LL_trio(x,y,z,l, 1,2) = SUM( OKA2P(x,y,:) * AKAP(:,z,l)  )
          
          ! PO + GGP = PO + HA
          do w=1,3
            Tmp(w) = OKA2P(x,y,w) * SUM( AKAP(w,:,l) * AKAP(:,z,l) ) 
          enddo
          LL_trio(x,y,z,l, 1,3) = SUM( Tmp )
          
          ! PO + U
          LL_trio(x,y,z,l, 1,4) = OKAP(x,y,l)
          
          ! GP + GP
          do w=1,3
            Tmp(w) = SUM(OKA2P(x,w,:) * AKAP(w,y,l) * AKAP(:,z,l) ) 
          enddo
          LL_trio(x,y,z,l, 2,2) = SUM( Tmp )
          
          ! GP + GGP
          do w=1,3
            do v=1,3
              Tmp2(v,w) = OKA2P(x,w,v) * AKAP(w,y,l) * SUM( AKAP(v,:,l) * AKAP(:,z,l) )
            enddo
          enddo
          LL_trio(x,y,z,l, 2,3) = SUM( Tmp2 )
          
          ! GP + U
          LL_trio(x,y,z,l, 2,4) = SUM( OKAP(x,:,l) * AKAP(:,y,l) )
          
          ! GGP + GGP
          do w=1,3
            do v=1,3
              Tmp2(v,w) = OKA2P(x,w,v) * SUM( AKAP(w,:,l) * AKAP(:,y,l) ) * SUM( AKAP(v,:,l) * AKAP(:,z,l) )
            enddo
          enddo
          LL_trio(x,y,z,l, 3,3) = SUM( Tmp2 )      
          
          ! GGP + U
          do w=1,3
            Tmp(w) = OKAP(x,w,l) * SUM( AKAP(w,:,l) * AKAP(:,y,l) )
          enddo
          LL_trio(x,y,z,l, 3,4) = SUM( Tmp )
        
          ! U + U
          LL_trio(x,y,z,l, 4,4) = OHWE(x,l)
          
          if (nRel == 4)  cycle  ! no FS
          !~~~~~~~~~~~~~~~~~~~~
          
          ! PO + FS
          do w=1,3
            Tmp(w) = OKA2P(x,y,w) * AKA2P(z,y,w) * AHWE(w,l) / AHWE(z,l)
          enddo
          LL_trio(x,y,z,l, 1,5) = SUM( Tmp )
          
          ! GP + FS
          do w=1,3
            do v=1,3
              Tmp2(v,w) = OKA2P(x,w,v) * AKA2P(z,w,v) * AKAP(w,y,l) * AHWE(v,l) / AHWE(z,l)
            enddo
          enddo
          LL_trio(x,y,z,l, 2,5) = SUM( Tmp2 )
          
          ! GGP + FS
          do w=1,3
            do v=1,3
              do u=1,3
                Tmp3(u,v,w) = OKA2P(x,w,v) * AKA2P(z,w,v) * AKAP(w,u,l) * AKAP(u,y,l) &
                 * AHWE(v,l) / AHWE(z,l)
              enddo
            enddo
          enddo
          LL_trio(x,y,z,l, 3,5) = SUM( Tmp3 )
          
          ! U + FS
          do w=1,3
            do v=1,3
              Tmp2(v,w) = OKA2P(x,w,v) * AKA2P(z,w,v) * AHWE(w,l) * AHWE(v,l) / AHWE(z,l)
            enddo
          enddo
          LL_trio(x,y,z,l, 4,5) = SUM( Tmp2 )
          
          ! FS + FS
          do w=1,3
            do v=1,3
              Tmp2(v,w) = OKA2P(x,w,v) * AKA2P(z,w,v) * AKA2P(y,w,v) * AHWE(w,l) * AHWE(v,l) &
               / (AHWE(z,l) * AHWE(y,l))
            enddo
          enddo
          LL_trio(x,y,z,l, 5,5) = SUM( Tmp2 )

        enddo
      enddo
    enddo
    
    ! PO + GP --> GP + PO etc.  
    do w=1,nRel-1
      do v=w+1,nRel
        do z=1,3
          do y=1,3
            LL_trio(:,y,z,l, v,w) = LL_trio(:,z,y,l, w,v)
          enddo
        enddo
      enddo
    enddo
  
  enddo
 

end subroutine Precalc_trioLLs

!===============================================================================

subroutine CalcLL(LL_array)
  use Global
  use Probs, ONLY: AHWE, OcA
  implicit none

  double precision, intent(OUT) :: LL_array(nRel,nRel,nTrios)
  integer :: i, j, rel_d, rel_s, l, y,z
  double precision :: PrL(nSnp), PGeno(3, nSnp, 0:nInd), Tmp1(3), Tmp2(3,3)


  ! PGeno: 
  ! probability that parent i has observed genotype o at locus l = 
  ! probability that parent has o given actual genotype a, times
  ! probability to have a as a random draw from population
  ! not conditional on its parents (when missing), else wrong grandparent can confuse signal
  
  do i=0, nInd
    do l = 1, nSnp
      if (i==0) then
        PGeno(:,l,0) = AHWE(:,l)   ! when no parent
      else
        Tmp1 = OcA(Geno(l,i), :) * AHWE(:, l)     
        PGeno(:, l, i) = Tmp1 / SUM(Tmp1)
      endif
    enddo
  enddo
  
  
  LL_array = Missing

  do j=1, nTrios
    if (.not. quiet .and. MOD(j,2000)==0)  print *, j
    if (Trios(1,j)==0)  cycle  ! offspring not genotyped
    do rel_s = 1,nRel
      if (Trios(3,j)==0 .and. rel_s/=4)  cycle  ! no sire in pedigree / not genotyped
      
      do rel_d = 1,nRel
        if (Trios(2,j)==0 .and. rel_d/=4)  cycle  ! no dam in pedigree / not genotyped
       
        
        PrL = 0D0
        do l=1, nSnp
          Tmp2 = 0D0
          do y=1,3  ! dam genotype
            do z=1,3  ! sire genotype
              Tmp2(y,z) = LL_trio(Geno(l,Trios(1,j)), y,z, l, rel_d, rel_s) * &
                PGeno(y,l, Trios(2,j)) * PGeno(z,l, Trios(3,j))
            enddo
          enddo    
          
          PrL(l) = LOG10( SUM( Tmp2 ))
          
        enddo
      
        LL_array(rel_d, rel_s, j) = SUM(PrL)
      
      enddo
    enddo
  enddo


end subroutine CalcLL


! ##############################################################################
! ##  Write output   ##
! ##############################################################################


subroutine writeped(LL_array, OppHom, CalcProbs, FileName)
  use Global
  use FileIO
  implicit none

  double precision, intent(IN) :: LL_array(nRel,nRel, nTrios)
  integer, intent(IN) :: OppHom(3, nTrios)
  logical, intent(IN) :: CalcProbs
  character(len=*), intent(IN) :: FileName
  integer :: i, rel_d, rel_s, rel_order(5), x
  double precision :: out_array(nRel,nRel, nTrios), tmp_array(nRel,nRel, nTrios)
  character(len=200) :: HeaderFMT, DataFMT
  character(len=2) :: relnames(nRel)
  character(len=10) :: data_header(nRel,nRel)


  out_array = Missing
  tmp_array = Missing
  do i=1,nTrios
    if (trios(1,i)==0)  cycle
    WHERE (LL_array(:,:,i)/=Missing)  tmp_array(:,:,i) = LL_array(:,:,i) - LL_array(4,4,i)
    if (CalcProbs) then
      do rel_s = 1,nRel
        if (Trios(3,i)==0 .and. rel_s/=4)  cycle 
        do rel_d = 1,nRel
          if (Trios(2,i)==0 .and. rel_d/=4)  cycle                
          out_array(rel_d, rel_s, i) = 10**tmp_array(rel_d, rel_s,i) / &
              SUM(10**tmp_array(:,:,i), MASK=LL_array(:, :, i)/=Missing)
        enddo
      enddo
    else
      out_array(:,:,i) = tmp_array(:,:,i)
    endif
  enddo
  
  ! swap column order if nRel=5: PO/GP/HA/UU/FS --> PO/FS/GP/HA/UU
  tmp_array = Missing
  if (nRel==5) then
    rel_order = (/1,5,2,3,4/)
    do x=1,5
      tmp_array(x,:,:) = out_array(rel_order(x),:,:)
    enddo
    do x=1,5
      out_array(:,x,:) = tmp_array(:,rel_order(x),:)
    enddo
  endif

  write(HeaderFMT, '( "(3(a", I0, ", 4X), 3a11, 100a12)" )')  ID_len
  if (CalcProbs) then
    write(DataFMT, '( "(3(a", I0, ", 4X), 3i11, 100f12.4)" )')  ID_len  ! round to 4 decimals
  else
    write(DataFMT, '( "(3(a", I0, ", 4X), 3i11, 100f12.2)" )')  ID_len  
  endif
  
  if (nRel==4) then
    relnames = (/'PO','GP','HA','UU'/)
  else
    relnames = (/'PO','FS','GP','HA','UU'/)
  endif
  do rel_s = 1,nRel
    do rel_d = 1,nRel
      if (CalcProbs) then
        data_header(rel_d, rel_s) = 'prob_'//relnames(rel_d)//'_'//relnames(rel_s)
      else
        data_header(rel_d, rel_s) = 'LLR_'//relnames(rel_d)//'_'//relnames(rel_s)
      endif
    enddo
  enddo

  open (unit=201,file=trim(FileName), status="unknown") 
    write (201, HeaderFMT)  Ped_header, 'OH_'//trim(Ped_header(2)), 'OH_'//trim(Ped_header(3)), 'ME_pair', & 
     RESHAPE(data_header, (/nRel*nRel/))
    do i=1,nTrios
      write (201,DataFMT) TrioNames(:,i), OppHom(:,i), RESHAPE(out_array(:,:,i), (/nRel*nRel/))
    enddo     
  close (201)

end subroutine writeped

!===============================================================================

subroutine deallocall
  use Global
  use FileIO
  use Probs
  implicit none

  if (allocated(Geno)) deallocate(Geno)
  if (allocated(Trios)) deallocate(Trios) 

  if (allocated(AHWE)) deallocate(AHWE)
  if (allocated(OHWE)) deallocate(OHWE) 
  if (allocated(AKAP)) deallocate(AKAP)
  if (allocated(OKAP)) deallocate(OKAP)

  if (allocated(LL_trio)) deallocate(LL_trio)
  if (allocated(Id)) deallocate(Id)
  if (allocated(TrioNames)) deallocate(TrioNames)

end subroutine deallocall