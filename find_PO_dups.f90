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
! gfortran -std=f95 -fall-intrinsics -O3 find_PO_dups.f90 -o findPO
! to debug:
! gfortran -std=f95 -fall-intrinsics -Wall -pedantic -fbounds-check -g -Og find_PO_dups.f90 -o findPO

!===============================================================================
!===============================================================================

module Global_variables
  implicit none

  integer :: nInd, nSnp, maxOH, maxDUP
  double precision :: minProb
  integer, parameter :: nchar_filename = 2000, nRel=6
  logical :: quiet, DoProbs
  character(len=20), allocatable, dimension(:) :: Id
  logical, allocatable, dimension(:) :: skip
  integer, allocatable, dimension(:,:) :: Genos

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
      nBothSNPd = nBothSNPd + IsBothScored(Genos(l,A), Genos(l,B))
    enddo

  end function nBothSNPd  
  
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  pure function nDiffer(A,B)
    integer, intent(IN) :: A, B
    integer :: nDiffer
    integer :: l
    
    nDiffer = 0
    do l=1,nSnp
      nDiffer = nDiffer + IsDifferent(Genos(l,A), Genos(l,B))
      if (nDiffer > maxDUP)  exit
    enddo

  end function nDiffer
  
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  pure function nOH(A,B)
    integer, intent(IN) :: A, B
    integer :: nOH
    integer :: l
  
    nOH = 0
    do l=1,nSnp
      nOH = nOH + IsOppHom(Genos(l,A), Genos(l,B))
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
! Probabilities/likelihoods, only used when --min_prob is invoked
!===============================================================================

Module Probs_module
  implicit none
  
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
        PrL(l) = LOG10( LL_Pair(Genos(l,i), Genos(l,j), l, r) )
      enddo
      LL(r) = SUM(PrL)
    enddo    
        
    ! likelihoods --> probabilities to be dup/PO/FS/...
    Probs = 0.0D0    
    tmp = LL - MAXVAL(LL)  ! scale to minimise problems with rounding to 0
    Probs = 10**tmp
    CalcProb = Probs / SUM(Probs)   ! scale to sum to 1 
    
    ! system_maxval = log10(HUGE(0D0))     use to avoid rounding issues ?

  end function CalcProb

  !===============================================================================

  subroutine Precalc_pairLLs
    use Global_variables
    implicit none

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

  subroutine PrecalcProbs(Er)
    use Global_variables
    implicit none

    double precision, intent(IN) :: Er
    integer :: h,i,j,k,l
    double precision ::  AF(nSnp), Tmp(3)

    ! allele frequencies
    do l=1,nSnp
      if (ANY(Genos(l,:)/=-1)) then
        AF(l)=float(SUM(Genos(l,:), MASK=Genos(l,:)/=-1))/(COUNT(Genos(l,:)/=-1)*2)
      else
        AF(l) = 1D0
      endif
    enddo
    ! TODO: allele frequencies from file


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
  
end module Probs_module

!===============================================================================
!===============================================================================

program main
  use Global_variables
  use Probs_module
  implicit none
  
  integer :: nArg, i, x
  character(len=32) :: arg, argOption
  character(len=nchar_filename) :: GenoFileName, OnlyListFileName, OutFileName, AF_FileName
  logical :: FileOK, DoDup, DoPO
  double precision :: Er
  
  ! defaults  ~~~~
  GenoFileName = 'Geno.txt'
  OnlyListFileName = 'NoFile'
  AF_FileName = 'NoFile'
  OutFileName = 'Pairs_maybe'
  DoPO = .FALSE.
  DoDup = .FALSE.
  DoProbs = .FALSE.
  quiet = .FALSE.
  maxOH = 0
  maxDUP = 0
  minProb = 0.0D0
  Er = 0.0D0

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
          
        case ('--max_dup', '--maxDUP', '--max_DUP')
          i = i+1
          call get_command_argument(i, argOption)
          read(argOption, *)  maxDUP
          
        case ('--max_oh', '--maxOH', '--max_OH')
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
  if (.not. DoDup .and. .not. DoPO) then
    print *, 'Please select --dup and/or --po'
    print *, ''
    call print_help()
    stop
  endif
  
  if (maxOH < 0) then
    print *, 'max_oh must be positive!'
    stop
  endif
  if (maxDup < 0) then
    print *, 'max_dup must be positive!'
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
    call PrecalcProbs(Er) 
    call Precalc_pairLLs() 
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

subroutine find_dups(FileName_part)
  use Global_variables
  use Probs_module
  use Calc
  implicit none
  
  character(len=nchar_filename), intent(IN) :: FileName_part 
  character(len=nchar_filename+7) :: FileName
  integer :: i,j, Lboth, Diff_ij, nPairs
  double precision :: probs_ij(nRel)
  
  FileName = trim(FileName_part)//'_DUP.txt'
  
  call init_lookups()
  
  nPairs = 0
  open(unit=201, file=trim(FileName), status='unknown')
    if (DoProbs) then
      write(201, '(2a7, 5X, 2a20, 2a10, 6(5X, a7))') 'row1', 'row2', 'ID1', 'ID2', 'nDiff', 'SnpdBoth', &
       'prob_S ', 'prob_PO', 'prob_FS', 'prob_GP', 'prob_HA', 'prob_UU'
    else
      write(201, '(2a7, 5X, 2a20, 2a10)') 'row1', 'row2', 'ID1', 'ID2', 'nDiff', 'SnpdBoth'
    endif
  
    do i=1, nInd-1
      if (.not. quiet .and. Modulo(i, 5000)==0)  print *, i
      do j=i+1, nInd
        if (skip(i) .and. skip(j))  cycle
        Diff_ij = nDiffer(i,j)
        if (Diff_ij > maxDUP)  cycle
        if (DoProbs) then
          probs_ij = CalcProb(i,j)  ! probs: self, PO, FS, HS/GP/FA, HA/3rd, U
          if (probs_ij(1) < minProb)  cycle
        endif  
        
        ! if arrived here, i+j are potential duplicate samples from same individual
        Lboth = nBothSNPd(i,j)  
        if (DoProbs) then
          write(201, '(2i7, 5X, 2a20, 2i7, 6f12.4)') i,j, Id(i), Id(j), Diff_ij, Lboth, probs_ij
        else
          write(201, '(2i7, 5X, 2a20, 2i7)') i,j, Id(i), Id(j), Diff_ij, Lboth
        endif
        nPairs = nPairs +1
          
      enddo
    enddo
  
  close(201)
  
  if (DoProbs) then
    write(*,'("Found ", i7, " duplicate pairs using maxDUP=", i4, " and min_Prob=", f7.3)') &
      nPairs, maxDUP, minProb
  else
    write(*,'("Found ", i7, " duplicate pairs using maxDUP=", i4)') nPairs, maxDUP 
  endif
  
end subroutine find_dups


!===============================================================================

subroutine find_PO(FileName_part)
  use Global_variables
  use Probs_module
  use Calc
  implicit none
  
  character(len=nchar_filename), intent(IN) :: FileName_part
  character(len=nchar_filename+7) :: FileName  
  integer :: i,j, Lboth, OH_ij, nPairs
  double precision :: probs_ij(nRel)
  
  FileName = trim(FileName_part)//'_PO.txt'
  
  call init_lookups()
  
  nPairs = 0
  open(unit=201, file=trim(FileName), status='unknown')
    if (DoProbs) then
      write(201, '(2a7, 5X, 2a20, 2a10, 6(5X, a7))') 'row1', 'row2', 'ID1', 'ID2', 'OH', 'SnpdBoth', &
       'prob_S ', 'prob_PO', 'prob_FS', 'prob_GP', 'prob_HA', 'prob_UU'
    else
      write(201, '(2a7, 5X, 2a20, 2a10)') 'row1', 'row2', 'ID1', 'ID2', 'OH', 'SnpdBoth'
    endif
  
    do i=1, nInd-1
      if (.not. quiet .and. Modulo(i, 5000)==0)  print *, i
      do j=i+1, nInd
        if (skip(i) .and. skip(j))  cycle
        OH_ij = nOH(i,j)
        if (OH_ij > maxOH)  cycle        
        if (DoProbs) then
          probs_ij = CalcProb(i,j)  ! probs: self, PO, FS, HS/GP/FA, HA/3rd, U
          if (probs_ij(2) < minProb)  cycle
        endif  
        
        ! if arrived here, i+j are potential parent-offspring pair
        Lboth = nBothSNPd(i,j) 
        if (DoProbs) then
          write(201, '(2i7, 5X, 2a20, 2i7, 6f12.4)') i, j, Id(i), Id(j), OH_ij, Lboth, probs_ij
        else
          write(201, '(2i7, 5X, 2a20, 2i7)') i, j, Id(i), Id(j), OH_ij, Lboth
        endif
        nPairs = nPairs +1
          
      enddo
    enddo
  
  close(201)
  
  if (DoProbs) then
    write(*,'("Found ", i7, " parent-offspring pairs using maxOH=", i4, " and min_Prob=", f7.3)') &
      nPairs, maxOH, minProb
  else
    write(*,'("Found ", i7, " parent-offspring pairs using maxOH=", i4)') nPairs, maxOH 
  endif

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