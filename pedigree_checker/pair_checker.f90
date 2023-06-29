! as pedigree_checker, but only for pairs, not trios
! adds 'self'

! input: either
! - file with 2 columns with IDs
! - file with multiple columns, incl 'IID1' & 'IID2', e.g. output from --dup

! compile: gfortran -std=f95 -fall-intrinsics -O3 pair_checker.f90 -o PairChecker

! ##############################################################################
! ##  Global variables  ##
! ##############################################################################

module global
  implicit none

  integer :: nInd, nSnp, nPairs, ID_len
  integer, parameter :: nchar_filename = 2000, nchar_ID = 40, nRel=6
  integer,allocatable,dimension(:) :: nBoth, nDiff, nOH
  integer,allocatable,dimension(:,:) :: Genos, Pairs
  logical :: quiet
  double precision :: OcA(-1:2,3), OKA2P(-1:2,3,3), AKA2P(3,3,3)
  double precision, parameter :: Missing = 999D0
  double precision, allocatable, dimension(:,:) :: AHWE, OHWE
  double precision, allocatable, dimension(:,:,:) :: AKAP, OKAP
  double precision, allocatable, dimension(:,:,:,:) :: LL_pair
  character(len=nchar_ID), allocatable, dimension(:) :: Id 
  character(len=nchar_ID), allocatable, dimension(:,:) :: PairNames
  
    !=========================
  contains
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! count number of columns in text file
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

program pair_checker
  use Global
  implicit none

  ! input
  integer :: x, i, nArg
  double precision :: Er
  character(len=32) :: arg, argOption
  character(len=nchar_filename) :: PairFileName, GenoFileName, OutFileName
  logical :: CalcProbs, FileOK
  ! output
  double precision, allocatable, dimension(:,:) :: LL_out_array 
    
    
  ! set default values
  PairFileName = 'Pairs.txt'
  GenoFileName = 'Geno.txt'
  OutFileName = 'Pairs_OUT.txt'
  Er = 0.005  
  CalcProbs = .TRUE.   ! transform log10-likelihoods into probabilities
  quiet = .FALSE.


  ! read arguments from command line
  nArg = command_argument_count()

  if (nArg == 0) then
    if (.not. quiet)  print *, 'Using default values, see --help'
    
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
          
        case ('--pairs')    
          i = i+1
          call get_command_argument(i, PairFileName)        
            
        case ('-o', '--out')
          i = i+1
          call get_command_argument(i, OutFileName)
          
        case ('--err')
          i = i+1
          call get_command_argument(i, argOption)
          read(argOption, *)  Er   ! TODO: 1 or 3 
          
        case ('--LLR')
          CalcProbs = .FALSE.
          
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
  inquire(file=trim(PairFileName), exist = FileOK)
  if (.not. FileOK) then
    write(*,*)  "Input file ", trim(PairFileName), " not found"
    stop
  endif

  inquire(file=trim(GenoFileName), exist = FileOK)
  if (.not. FileOK) then
    write(*,*)  "Genotype file ", trim(GenoFileName), " not found"
    stop
  endif

  ! read in data
  if (.not. quiet)  print *, "Reading genotype data in "//trim(GenoFileName)//" ... "
  call ReadGeno(GenoFileName) 
  
  if (.not. quiet)  print *, "Reading pairs in "//trim(PairFileName)//" ... "
  call ReadPairs(PairFileName)
  if (.not. quiet)  print *, "Read ", nPairs, " pairs. "

  ! count number of opposing homozygous loci for each pair
  if (.not. quiet)  print *, "Counting number of opposing homozygous loci ... "
  allocate(nBoth(nPairs))
  allocate(nDiff(nPairs))
  allocate(nOH(nPairs))
  call CountDifs()

  ! calculate log-likelihoods
  if (.not. quiet)  print *, "Calculating log-likelihoods ... "
  call PrecalcProbs(Er) 
  allocate(LL_out_array(nRel,nPairs))
  call Precalc_pairLLs()
  call CalcLL(LL_out_array)

  ! write to file
  if (.not. quiet)  print *, "Writing output to file ... "
  call writePairs(LL_out_array, CalcProbs, OutFileName)
  call DeAllocAll()
  if (.not. quiet)  print *, ""
  if (.not. quiet)  print *, "Done."

  
  !=========================
  contains
    subroutine print_help()
        print '(a)',    'Calculate relationship probabilities of candidate parents'
        print '(a, /)', 'command-line options:'
        print '(a)',    '  -h, --help          print usage information and exit'
        print '(a)',    '  --pairs <filename>  input file with pairs: either 2 columns,',&
                        '                        or multiple columns including IID1 & IID2',& 
                        '                        Default: Pairs.txt'             
        print '(a)',    '  --geno <filename>   input file with genotype data. Default: Geno.txt'       
        print '(a)',    '  --err               presumed genotyping error rate; default: 0.005'
        print '(a)',    '  --LLR               do not transform LLRs to probabilities'
        print '(a)',    '  --noFS              assume pairs cannot be full siblings'
        print '(a)',    '  --out <filename>    output file with pairs + OH counts + Likelihoods;',&
                        '                       default: Pairs_OUT.txt'
        print '(a)',    '  --quiet             suppress all messages'
    end subroutine print_help 

end program pair_checker


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

subroutine ReadPairs(FileName)
  use Global
  implicit none

  character(len=*), intent(IN) :: FileName
  integer :: i, j, k, IOerr, nCol, IDcol(2)
  character(len=nchar_ID), allocatable :: tmpC(:)
  character(len=50), allocatable :: header(:)

  nPairs = FileNumRow(trim(FileName)) -1  ! 1st row = header 
  allocate(PairNames(2,nPairs))   ! typically offspring - dam - sire, but not necessarily
  PairNames = "NA"

  nCol = FileNumCol(trim(FileName))
  allocate(tmpC(nCol))
  allocate(header(nCol))
  header = 'NA'
  IDcol = 0
  
  open(unit=103, file=trim(FileName), status="old")
    read(103,*)  header          
    if (nCol == 2) then
      IDcol = (/1,2/)
    else
      do k=1, nCol
        if (header(k)=='IID1' .or. header(k)=='ID1')  IDcol(1) = k
        if (header(k)=='IID2' .or. header(k)=='ID2')  IDcol(2) = k  
      enddo
    endif   
    
    do i=1,nPairs
      read(103, *,IOSTAT=IOerr)  tmpC
      if (IOerr > 0) then
        print *, "Wrong input in file "//trim(FileName)//" on line ", i
        stop
      else if (IOerr < 0) then
        exit   ! EOF
      else
        do k=1,2
          PairNames(k,i) = tmpC( IDcol(k) )
        enddo
      end if
    enddo
  close(103)

  ! Pairs names to row numbers in genotype file
  allocate(Pairs(2,nPairs))  ! row numbers in genotype data
  Pairs = 0
  do i = 1, nPairs
    do k = 1,2
      if (PairNames(k,i)=='NA')  cycle
      do j = 1, nInd
        if (PairNames(k,i) == Id(j)) then
          Pairs(k,i) = j
          exit  ! break from inner loop
        endif
      enddo
    enddo
  enddo
  
end subroutine ReadPairs


! ##############################################################################
! ##   Count scored both + diffs + OH   ##
! ##############################################################################

subroutine CountDifs
  use Global
  implicit none

  integer :: i,l, A,B, nBoth_i, nDiff_i, nOH_i

  nBoth = -9
  nDiff = -9
  nOH = -9
  do i=1,nPairs
    if (Pairs(1,i)==0 .or. Pairs(2,i)==0) then
      nBoth(i) = 0
      cycle
    endif
    A = Pairs(1,i)
    B = Pairs(2,i)
    nBoth_i = 0
    nDiff_i = 0
    nOH_i = 0
    do l=1,nSnp
      if (Genos(l,A)==-1 .or. Genos(l,B)==-1)  cycle
      nBoth_i = nBoth_i +1
      if (Genos(l,A) == Genos(l,B))  cycle
      nDiff_i = nDiff_i +1      
      if ((Genos(l,A)==0 .and.Genos(l,B)==2) .or. (Genos(l,A)==2 .and. Genos(l,B)==0)) then
        nOH_i = nOH_i+1
      endif                       
    enddo
    nBoth(i) = nBoth_i
    nDiff(i) = nDiff_i
    nOH(i) = nOH_i
  enddo

end subroutine CountDifs


! ##############################################################################
! ##   Calculate likelihoods   ##
! ##############################################################################

subroutine PrecalcProbs(Er)
  use Global
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

!===============================================================================

subroutine Precalc_pairLLs
  use Global
  implicit none

  ! calculate likelihoods for: self, PO, FS, HS/GP/FA, HA/3rd, U
  ! store in look-up table: G_IID1 (0/1/2/-9) + G_IID2 (0/1/2/-9)

  ! when not conditioning on parents of both individuals, LL(FA)=LL(GP)=LL(HS),
  ! and similarly all 3rd degree relationships have same LL (log likelihood)

  ! assume no inbreeding, and otherwise unrelated. 

  integer :: l, x, y, z, w, v
  double precision :: Tmp(3), Tmp2(3,3), Tmp3(3,3,3)

  allocate(LL_Pair(-1:2,-1:2, nSnp, nRel))  ! G_IID1, G_IID2, snp, rel (S/PO/FS/GP/HA/U)
  LL_Pair = Missing

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

subroutine CalcLL(LL_array)
  use Global
  implicit none

  double precision, intent(OUT) :: LL_array(nRel,nPairs)
  integer :: i, r, l
  double precision :: PrL(nSnp)
  
  LL_array = Missing

  do i=1, nPairs
    if (.not. quiet .and. MOD(i,2000)==0)  print *, i
    if (Pairs(1,i)==0 .or. Pairs(2,i)==0)  cycle   ! either one not genotyped
    do r=1, nRel
      PrL = 0D0
      do l=1, nSnp
        PrL(l) = LOG10( LL_Pair(Genos(l, Pairs(1,i)), Genos(l, Pairs(2,i)), l, r) )
      enddo
      LL_array(r, i) = SUM(PrL)
    
    enddo
  enddo

end subroutine CalcLL

! ##############################################################################
! ##  Write output   ##
! ##############################################################################


subroutine writePairs(LL_array, CalcProbs, FileName)
  use Global
  implicit none

  double precision, intent(IN) :: LL_array(nRel, nPairs)
  logical, intent(IN) :: CalcProbs
  character(len=*), intent(IN) :: FileName
  integer :: i, r
  double precision :: out_array(nRel, nPairs), tmp_array(nRel,nPairs)
  character(len=200) :: HeaderFMT, DataFMT
  character(len=2) :: relnames(nRel)
  character(len=10) :: data_header(nRel)


  out_array = Missing
  tmp_array = 0D0
  do i=1,nPairs
    if (pairs(1,i)==0 .or. Pairs(2,i)==0)  cycle
    if (CalcProbs) then
      ! scale to prevent problems with rounding to 0
      WHERE (LL_array(:,i)/=Missing)  tmp_array(:,i) = LL_array(:,i) - maxval(LL_array(:,i))
      do r=1, nRel
        if (LL_array(r,i) == Missing)  cycle          
        out_array(r,i) = 10**tmp_array(r,i) / SUM(10**tmp_array(:,i), MASK=LL_array(:,i)/=Missing)
      enddo 
      ! make sure sums to 1  
      WHERE (out_array(:,i)/=Missing)  out_array(:,i) = out_array(:,i) / SUM(out_array(:,i))
    else
      ! scale all LLs by LL(U/U)
      WHERE (LL_array(:,i)/=Missing) out_array(:,i) = LL_array(:,i) - LL_array(nRel,i)
    endif
  enddo
 

  write(HeaderFMT, '( "(2(a", I0, ", 4X), 3a11, 6(5X, a7))" )')  ID_len
  if (CalcProbs) then
    write(DataFMT, '( "(2(a", I0, ", 4X), 3i11, 6f12.4)" )')  ID_len  ! round to 4 decimals
  else
    write(DataFMT, '( "(2(a", I0, ", 4X), 3i11, 6f12.2)" )')  ID_len  
  endif
  
  relnames = (/' S', 'PO', 'FS', 'GP','HA','UU'/)
  do r = 1,nRel
    if (CalcProbs) then
      data_header(r) = 'prob_'//relnames(r)
    else
      data_header(r) = 'LLR_'//relnames(r)
    endif
  enddo

  open (unit=201,file=trim(FileName), status="unknown") 
    write (201, HeaderFMT) 'ID1', 'ID2', 'nBoth', 'nDiff', 'nOH', data_header
    do i=1,nPairs
      write (201,DataFMT) PairNames(:,i), nBoth(i), nDiff(i), nOH(i), out_array(:,i)
    enddo     
  close (201)

end subroutine writePairs

!===============================================================================

subroutine deallocall
  use Global
  implicit none

  if (allocated(Genos)) deallocate(Genos)
  if (allocated(Pairs)) deallocate(Pairs) 

  if (allocated(AHWE)) deallocate(AHWE)
  if (allocated(OHWE)) deallocate(OHWE) 
  if (allocated(AKAP)) deallocate(AKAP)
  if (allocated(OKAP)) deallocate(OKAP)

  if (allocated(LL_pair)) deallocate(LL_pair)
  if (allocated(Id)) deallocate(Id)
  if (allocated(PairNames)) deallocate(PairNames)

end subroutine deallocall