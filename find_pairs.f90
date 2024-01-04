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
! --min_prob : optional, select also based on LLR-based probability to be PO/dup
! --err : genotyping error rate; only relevant in combination with min_prob
! --af : optional, only relevant in combination with min_prob
! --out : output filename  (_PO.txt & _DUP.txt)
! Output:
! - text file with ID1 - ID2 - OH - SnpdBoth  for pairs with OH <= maxOH
!
! to compile:
! gfortran -O3 find_PO_dups.f90 -o findPO
! to debug:
! gfortran -Wall -pedantic -fbounds-check -g -Og find_PO_dups.f90 -o findPO

! TODO: timestamps

!===============================================================================
!===============================================================================

module Global_variables
  implicit none

  integer :: nInd, nSnp, maxOH, maxDIF
  double precision :: minProb
  integer, parameter :: nchar_filename = 2000, nRel=6
  logical :: quiet, DoProbs
  character(len=20), allocatable, dimension(:) :: Id
  logical, allocatable, dimension(:) :: skip
  integer(kind=1), allocatable, dimension(:,:) :: Geno
  integer, parameter :: ik8 = selected_int_kind(8)
  character(len=2) :: rel_suffix(0:(nRel-1))
  character(len=32) :: rel_lbls(0:(nRel-1))

end module Global_variables

!===============================================================================

module Calc
  use Global_variables
  implicit none
  
  integer :: IsBothScored(-1:2,-1:2), IsOppHom(-1:2,-1:2), IsDifferent(-1:2, -1:2)   

  contains 
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine init_lookups
  
    IsBothScored = 1
    IsBothScored(-1,:) = 0
    IsBothScored(:,-1) = 0
    
    IsDifferent = 1
    IsDifferent(-1,:) = 0
    IsDifferent(:,-1) = 0
    IsDifferent(0,0) = 0
    IsDifferent(1,1) = 0
    IsDifferent(2,2) = 0

    IsOppHom = 0
    IsOppHom(0,2) = 1
    IsOppHom(2,0) = 1 
  
  end subroutine init_lookups
  
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  pure function nBothSNPd(A,B)
    integer, intent(IN) :: A, B
    integer :: nBothSNPd
    integer :: l
    
    nBothSNPd = 0
    do l=1,nSnp
      nBothSNPd = nBothSNPd + IsBothScored(Geno(l,A), Geno(l,B))
    enddo

  end function nBothSNPd  
  
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  pure function nDiffer(A,B)
    integer, intent(IN) :: A, B
    integer :: nDiffer
    integer :: l
    
    nDiffer = 0
    do l=1,nSnp
      nDiffer = nDiffer + IsDifferent(Geno(l,A), Geno(l,B))
      if (nDiffer > maxDIF)  exit
    enddo

  end function nDiffer
  
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  pure function nOH(A,B)
    integer, intent(IN) :: A, B
    integer :: nOH
    integer :: l
  
    nOH = 0
    do l=1,nSnp
      nOH = nOH + IsOppHom(Geno(l,A), Geno(l,B))
      if (nOH > maxOH) exit               
    enddo

  end function nOH 

end module Calc

!===============================================================================

module FileDim
  implicit none

  contains
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    integer function FileNumCol(FileName)
      implicit none

      character(len=*), intent(IN) :: FileName
      integer :: j, strLen, numcol
      character(len=:), allocatable :: line
      
      allocate(character(len=500000) :: line)

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
      
      deallocate(line)

    end function FileNumCol

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    integer function FileNumRow(FileName)
      implicit none

      character(len=*), intent(IN) :: FileName
      integer :: nrow, i, maxRow, IOerr 
      character(len=:), allocatable :: dumC  
      
      allocate(character(len=500000) :: dumC)

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
      
      deallocate(dumC)

    end function FileNumRow
  
end module FileDim


!===============================================================================
! Probabilities/likelihoods, only used when --min_prob is invoked
!===============================================================================

Module Probs_module
  use Global_variables
  implicit none
  private :: getAF
  
  double precision :: OcA(-1:2,3), OKA2P(-1:2,3,3), AKA2P(3,3,3)
  double precision, allocatable, dimension(:,:) :: AHWE, OHWE
  double precision, allocatable, dimension(:,:,:) :: AKAP, OKAP
  double precision, allocatable, dimension(:,:,:,:) :: LL_pair 
  
  !======================
  contains
  
  function CalcProb(i,j)
    use Global_variables
    
    integer, intent(IN) :: i, j
    double precision :: CalcProb(nRel)
    integer :: r, l
    double precision :: PrL(nSnp), LL(nRel), tmp(nRel), Probs(nRel)
    
    LL = 0D0
    do r=1, nRel
      PrL = 0D0
      do l=1, nSnp
        PrL(l) = LOG10( LL_Pair(Geno(l,i), Geno(l,j), l, r) )
      enddo
      LL(r) = SUM(PrL)
    enddo    
        
    ! likelihoods --> probabilities to be dup/PO/FS/...   
    tmp = LL - MAXVAL(LL)  ! scale to minimise problems with rounding to 0
    Probs = 10**tmp
    CalcProb = Probs / SUM(Probs)   ! scale to sum to 1 

    ! system_maxval = log10(HUGE(0D0))     use to avoid rounding issues ?

  end function CalcProb

  !=============================================================================

  subroutine Precalc_pairLLs
    ! calculate likelihoods for: self, PO, FS, HS/GP/FA, HA/3rd, U
    ! store in look-up table: G_IID1 (0/1/2/-9) + G_IID2 (0/1/2/-9)

    ! when not conditioning on parents of both individuals, LL(FA)=LL(GP)=LL(HS),
    ! and similarly all 3rd degree relationships have same LL (log likelihood)

    ! assume no inbreeding, and otherwise unrelated. 

    integer :: l, x, y, z, w, v
    double precision :: Tmp(3), Tmp2(3,3), Tmp3(3,3,3)

    allocate(LL_Pair(-1:2,-1:2, nSnp, nRel))  ! G_IID1, G_IID2, snp, rel (S/PO/FS/GP/HA/U)
    LL_Pair = 0D0

    do l = 1, nSnp
      do x=-1,2  ! observed genotype IID1
        do y=-1,2    ! observed genotype IID2
   
          ! Self
          do z=1,3  ! actual genotype IID2
            Tmp(z) = OcA(x,z) * OcA(y,z) * AHWE(z,l)
          enddo
          LL_Pair(x,y,l, 1) = SUM(Tmp)
          
          ! PO
          do z=1,3  ! actual genotype IID2
            Tmp(z) = OKAP(x,z,l) * OcA(y,z) * AHWE(z,l)
          enddo
          LL_Pair(x,y,l, 2) = SUM(Tmp)
          
          ! FS
          do v=1,3  ! actual genotype shared dam
            do w=1,3  ! actual genotype shared sire
              Tmp2(v,w) = OKA2P(x,v,w) * OKA2P(y,v,w) * AHWE(v,l) * AHWE(w,l)
            enddo
          enddo
          LL_Pair(x,y,l, 3) = SUM(Tmp2)
          
          ! HS/GP/FA
          do z=1,3  ! actual genotype IID2
            do w=1,3  ! parent of IID1 = offspring of IID2
              Tmp2(z,w) = OKAP(x,w,l) * AKAP(w,z,l) * OcA(y,z) * AHWE(z,l)
            enddo
          enddo
          LL_Pair(x,y,l, 4)  = SUM(Tmp2)
                   
          ! GGP/HA/3rd
          do z=1,3  ! actual genotype IID2
            do w=1,3  ! parent of IID1
              do v=1,3  ! grandparent of IID1 = offspring of IID2
                Tmp3(z,w,v) = OKAP(x,w,l) * AKAP(w,v,l) * AKAP(v,z,l) * OcA(y,z) * AHWE(z,l)
              enddo
            enddo
          enddo
          LL_Pair(x,y,l, 5)  = SUM(Tmp3)
          
          ! U
          LL_Pair(x,y,l, 6) = OHWE(x,l) * OHWE(y,l)
          
        enddo
      enddo
    enddo

  end subroutine Precalc_pairLLs

  !===============================================================================

  subroutine PrecalcProbs(Er, AF_FileName)
    use Global_variables
    implicit none

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
  
  !=============================================================================
  function getAF(FileName)
    use FileDim    
    character(len=*), intent(IN) :: FileName
    double precision :: getAF(nSnp), AF(nSnp)
    integer :: l, nCol, nRow, AFcol, k, IOerr, Acnt
    character(len=50), allocatable :: header(:), tmpC(:)
    double precision :: tmpD
    
    AF = 1D0
    
    if (FileName == 'NoFile') then
    
      do l=1,nSnp
        if (ANY(Geno(l,:)/=-1)) then
          Acnt = SUM(int(Geno(l,:),kind=ik8), MASK=Geno(l,:)/=-1)
          AF(l)=dble(Acnt)/(COUNT(Geno(l,:)/=-1)*2)
        else
          AF(l) = 1D0
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
  
end module Probs_module

!===============================================================================
!===============================================================================

program main
  use Global_variables
  use Probs_module, ONLY: Precalc_pairLLs, PrecalcProbs
  implicit none
  
  integer :: nArg, i, x, r
  character(len=32) :: arg, argOption
  character(len=nchar_filename) :: GenoFileName, OnlyListFileName, OutFileName, AF_FileName
  logical :: FileOK, focal_rels(0:(nRel-1))  ! , DoDup, DoPO
  double precision :: Er
  
  ! defaults  ~~~~
  GenoFileName = 'Geno.txt'
  OnlyListFileName = 'NoFile'
  AF_FileName = 'NoFile'
  OutFileName = 'Pairs_maybe'
!  DoPO = .FALSE.
!  DoDup = .FALSE.
  DoProbs = .FALSE.
  quiet = .FALSE.
  maxOH = 0
  maxDIF = 0
  minProb = 0.0D0
  Er = 0.0D0
  focal_rels = .FALSE.
  
  rel_suffix = (/'S ', 'PO', 'FS', 'GP', 'HA', 'UU'/)
  rel_lbls(0) = 'duplicate'
  rel_lbls(1) = 'parent-offspring'
  rel_lbls(2) = 'full sibling'
  rel_lbls(3) = '2nd degree relatives'
  rel_lbls(4) = '3rd degree relatives'
  rel_lbls(5) = 'unrelated' 

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
          
        case('--dup')
          focal_rels(0) = .TRUE.
         
        case('--po')
          focal_rels(1) = .TRUE.
          
        case ('--focal')
          i = i+1
          call get_command_argument(i, argOption)
          if (any(rel_suffix == argOption(1:2))) then
            do r=0,5
              if (rel_suffix(r) == argOption(1:2)) then
                focal_rels(r) = .TRUE.
              endif
            enddo
          else
            print *, 'Invalid option for --focal!'
            write(*, '("Valid are: ", 6a4)') rel_suffix
            stop
          endif
          
        case ('--max_dif')
          i = i+1
          call get_command_argument(i, argOption)
          read(argOption, *)  maxDIF
          
        case ('--max_oh', '--max_OH')
          i = i+1
          call get_command_argument(i, argOption)
          read(argOption, *)  maxOH    
        
        case ('--min_prob', '--minProb')
          i = i+1
          call get_command_argument(i, argOption)
          read(argOption, *)  minProb
          DoProbs = .TRUE.
          
        case ('--err')
          i = i+1
          call get_command_argument(i, argOption)
          read(argOption, *)  Er   ! TODO: length 3 
          
        case ('--af', '--maf', '--freq')
          i = i+1
          call get_command_argument(i, AF_FileName)
         
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
  if (.not. any(focal_rels)) then
    print *, 'Please select a focal relationship, via --dup, --po and/or --focal'
    print *, 'See --help for help'
    stop
  endif
  
  if (maxOH < 0) then
    print *, 'max_oh must be positive!'
    stop
  endif
  if (maxDIF < 0) then
    print *, 'max_dup must be positive!'
    stop
  endif
  if (any(focal_rels(2:5)) .and. .not. DoProbs) then
    print *, 'with --focal different from S or PO, --min_prob and --err must be specified'
    stop
  endif
  if (DoProbs) then
    if (minProb < 0.0 .or. minProb > 1.0) then
      print *, 'min_prob must be between 0 and 1!'
      stop
    else if (Er <= 0.0 .or. Er > 0.5) then
      print *, 'when using --min_prob, please provide a genotyping error rate --err >0 and <0.5'
      stop
    endif
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
!  if (.not. quiet)  write(*, '("Selected forcal relationships: ", 6a4)') rel_suffix  TODO
  if (.not. quiet)  print *, 'Reading genotype data in '//trim(GenoFileName)//' ...'
  call ReadGeno(GenoFileName)
  if (.not. quiet)  print *, "Read ", nInd, "Individuals and ", nSnp, " SNPs."
 
  allocate(skip(nInd))
  skip = .FALSE.
  if (OnlyListFileName /= 'NoFile') then
    if (.not. quiet)  print *, "Reading individuals in --only file "//trim(OnlyListFileName)//" ... " 
    call ReadOnlyList(OnlyListFileName)
  endif
  
  if (DoProbs) then
    if (.not. quiet)  print *, "Pre-calculating log-likelihoods ... "
    call PrecalcProbs(Er, AF_FileName) 
    call Precalc_pairLLs() 
  endif
  
  do r=0,5
    if (.not. focal_rels(r))  cycle
    if (.not. quiet)  print *, 'Searching for '//trim(rel_lbls(r))//' pairs ...'
    call find_pairs(OutFileName, r)
  enddo
  
  ! free up allocated variables  ~~~~
  call deallocall()
  

  ! help  ~~~~
  ! TODO update
  contains
    subroutine print_help()
        print '(a)',  'Find potential parent-offspring pairs based on OH count'
        print '(a, /)', 'command-line options:'
        print '(a)',    '  -h, --help          print usage information and exit'
        print '(a)',    '  --geno <filename>   file with genotype data. Default: Geno.txt'
        print '(a)',    '  --dup               search for duplicates'   
        print '(a)',    '  --po                search for potential parent-offspring pairs' 
        print '(a)',    '  --focal <REL>       search for pairs with relationship REL, where ', &
                        '                        REL is S, PO, FS, GP, HA, or UU'
        print '(a)',    '  --max_dup <n>       maximum number of differences between duplicate samples' 
        print '(a)',    '  --max_oh <n>        maximum OH count for potential PO pairs. Default: 1', &
                        '                        see ?sequoia::MaxMismatch in R'  
        print '(a)',    '  --min_prob <p>      optional, select also based on LLR-based probability'  
        print '(a)',    '  --err <value>       presumed genotyping error rate; only relevant in ',&
                        '                        combination with --min_prob'
        print '(a)',    '  --af <filename>     optional input file with allele frequencies; only relevant',&
                        '                        in combination with --min_prob. Either 1 column and no header,',&
                        '                        or multiple columns with a column MAF, AF, or Frequency',&
                        '                        E.g. output from plink --freq.'                      
        print '(a)',    '  --only <filename>   only calculate OH when one or both are in this subset'
        print '(a)',    '  --out <filename>    output file name. Default: Pairs_maybe_PO.txt' 
        print '(a)',    '  --quiet             suppress all messages'
    end subroutine print_help

end program main


!===============================================================================
!===============================================================================

subroutine find_pairs(FileName_part, focal)
  use Global_variables
  use Probs_module, ONLY: CalcProb
  use Calc
  implicit none
  
  character(len=nchar_filename), intent(IN) :: FileName_part 
  integer, intent(IN) :: focal  ! 0=self (duplicates), 1=PO, 2=FS, 3=GP, 4=HA, 5=U
  character(len=nchar_filename+7) :: FileName
  integer :: i,j, Lboth, nPairs, cnt_ij, max_cnt
  double precision :: probs_ij(nRel)
  character(len=:), allocatable :: summary_lbl
  character(len=:), allocatable :: nx
  
  if (focal==0) then
    FileName = trim(FileName_part)//'_DUP.txt'
  else if (focal > 0 .and. focal <= 5) then
    FileName = trim(FileName_part)//'_'//rel_suffix(focal)//'.txt'
  else
    print *, "invalid focal relationship!"
    stop
  endif   
  
  if (focal==0) then
    max_cnt = maxDIF
    nx = 'nDiff'
  else if (focal==1) then
    max_cnt = maxOH
    nx = 'OH'
  else
    max_cnt = 0  ! filter on probability only
    nx = 'xx'
  endif
  
  call init_lookups()
  
  nPairs = 0
  cnt_ij = 0
  open(unit=201, file=trim(FileName), status='unknown')
    if (DoProbs) then
      write(201, '(2a7, 2a20, 5X, 2a10, 6(5X, a7))') 'row1', 'row2', 'ID1', 'ID2', nx, 'SnpdBoth', &
       'prob_S ', 'prob_PO', 'prob_FS', 'prob_GP', 'prob_HA', 'prob_UU'
    else
      write(201, '(2a7, 2a20, 5X, 2a10)') 'row1', 'row2', 'ID1', 'ID2', nx, 'SnpdBoth'
    endif
  
    do i=1, nInd-1
      if (.not. quiet .and. Modulo(i, 5000)==0)  print *, i
      do j=i+1, nInd
        if (skip(i) .and. skip(j))  cycle
        if (focal==0) then
          cnt_ij = nDiffer(i,j)
        else if (focal==1) then
          cnt_ij = nOH(i,j)
        endif
        if (cnt_ij > max_cnt)  cycle
        if (DoProbs) then
          probs_ij = CalcProb(i,j)  ! probs: self, PO, FS, HS/GP/FA, HA/3rd, U
          if (probs_ij(focal+1) < minProb)  cycle
        endif  
        
        ! if arrived here, i+j are potential duplicate samples from same individual
        Lboth = nBothSNPd(i,j)  
        if (DoProbs) then
          write(201, '(2i7, 5X, 2a20, 2i10, 6f12.4)') i,j, Id(i), Id(j), cnt_ij, Lboth, probs_ij
        else
          write(201, '(2i7, 5X, 2a20, 2i10)') i,j, Id(i), Id(j), cnt_ij, Lboth
        endif
        nPairs = nPairs +1         
      enddo
    enddo
  
  close(201)
  
  summary_lbl = '(" Found ", i7, " '//trim(rel_lbls(focal))//' pairs using '
  if (focal == 0) then
    if (DoProbs)       summary_lbl = summary_lbl//'maxDIF=", i4, " and '
    if (.not. DoProbs) summary_lbl = summary_lbl//'maxDIF=", i4)'
  else if (focal == 1) then
    if (DoProbs)       summary_lbl = summary_lbl//'maxOH=", i4, " and '
    if (.not. DoProbs) summary_lbl = summary_lbl//'maxOH=", i4)'
  endif
  if (DoProbs)     summary_lbl = summary_lbl//'min_Prob=", f7.3)'
    
  if (focal==0 .and. .not. DoProbs)  write(*, summary_lbl)   nPairs, maxDIF
  if (focal==0 .and. DoProbs)        write(*, summary_lbl)   nPairs, maxDIF, minProb
  if (focal==1 .and. .not. DoProbs)  write(*, summary_lbl)   nPairs, maxOH
  if (focal==1 .and. DoProbs)        write(*, summary_lbl)   nPairs, maxOH, minProb
  if (focal>1)                       write(*, summary_lbl)   nPairs, minProb
  
end subroutine find_pairs


!===============================================================================
!===============================================================================

subroutine ReadGeno(GenoFileName)
  use Global_variables
  use FileDim
  implicit none

  character(len=nchar_filename), intent(IN) :: GenoFileName
  integer :: i, l
  integer(kind=1), allocatable, dimension(:) :: GenoV

  nSnp = FileNumCol(trim(GenoFileName)) -1  ! column 1 = IDs
  nInd = FileNumRow(trim(GenoFileName))  

  allocate(GenoV(nSnp))
  allocate(Geno(nSnp, nInd))  
  Geno = -1
  allocate(Id(nInd))
  Id = "NA"

  open (unit=101,file=trim(GenoFileName),status="old")
  do i=1,nInd
    if (nInd > 5000) then
      if (.not. quiet .and. MODULO(i,2000)==0)  write(*,'(i10, 2x)', advance='no') i
    endif
    read (101,*)  Id(i), GenoV
    do l=1,nSnp
      if (GenoV(l)>=0) then
        Geno(l,i) = GenoV(l)  
      endif
    enddo
  enddo
  close (101)
  
  deallocate(GenoV)
  if (nInd > 5000) write(*,*) ''  ! advance to next line
  
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

subroutine deallocall
  use Global_variables
  use Probs_module

  if (allocated(ID)) deallocate(ID)
  if (allocated(skip)) deallocate(skip)
  if (allocated(Geno)) deallocate(Geno)
  if (allocated(AHWE)) deallocate(AHWE)
  if (allocated(OHWE)) deallocate(OHWE)
  if (allocated(AKAP)) deallocate(AKAP)
  if (allocated(OKAP)) deallocate(OKAP)
  if (allocated(LL_pair)) deallocate(LL_pair)

end subroutine deallocall

!===============================================================================
! end. 