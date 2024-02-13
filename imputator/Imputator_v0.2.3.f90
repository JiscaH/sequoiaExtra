! clean & impute SNP genotype data, using a pedigree
! each SNP is done separately, no use is made of LD
!
! Jisca Huisman, jisca.huisman@gmail.com
! 2024-02-12
!
! This code is available under GNU General Public License v3
!
! Compilation: 
! gfortran -O3 sqa_general.f90 sqa_fileIO.f90 Imputator_v0.2.3.f90 -o imputator
!
!===============================================================================

module global_variables 
  use sqa_general, ONLY: ishort, nchar_ID
  implicit none

  integer :: nIndG, nInd_max, nIndT, nSnp, nMatings
  integer(kind=ishort), allocatable, target :: Geno(:,:)
  character(len=nchar_ID), allocatable :: SNP_names(:)
  double precision :: Er
  double precision, allocatable :: AF(:)
  logical :: quiet
  
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! component type of pedigree vector
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  type individual
    integer :: index
    character(len=nchar_ID) :: ID='NA'
    integer :: sex=3  ! 1=female, 2=male, 3=unknown, 4=hermaphrodite
    integer, dimension(2) :: parent = 0   ! dam, sire
    integer :: nMates = 0
    integer, allocatable :: matings_m(:)  ! index in matings array of mates & offspring
 !   integer :: parent_m = 0   !  index in matings array of parents & full siblings
    type(matepair), pointer :: parentpair => null()
    ! TODO: birth year 
  end type individual
  
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! component type of matings vector
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  type matepair
    integer :: index
    integer :: nOffspring = 0, parent_q(2) = 0   ! parent_i = (/dam_ped_index, sire_ped_index/)
    integer, allocatable :: offspring_q(:)
  end type matepair

end module global_variables

!===============================================================================

module pedigree_fun
  use sqa_general, ONLY: chunk_size_large, nchar_ID
  use global_variables, ONLY: nIndG, nIndT, nInd_max, nMatings, individual, matepair
  implicit none
   
  type(individual), allocatable, target :: pedigree(:)
  type(individual), target, private :: indiv0  
  type(matepair), allocatable, target :: matings(:)
  
  private :: ped_add, grow_pedigree
  
contains
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! initialise pedigree with individuals read from genotype file
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine init_pedigree(IdV)  
    character(len=nchar_ID), intent(IN) :: IdV(nIndG)
    integer :: i
    
    indiv0 = individual(index=0)  
    nInd_max = nIndG
    nIndT = nIndG
    allocate(pedigree(0:nInd_max))      
    do i=1,nIndG
      call ped_add(individual(index=i, ID=IdV(i), sex=3))
    enddo
  end subroutine init_pedigree
  
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! add individuals to pedigree; only at initalisation. no parents.
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine ped_add(ind_new)
    type(individual), intent(IN) :: ind_new  
    integer :: i
    
    i = ind_new%index 
    
    if (i > nInd_max) then
      call grow_pedigree()  ! increase pedigree with another chunk
    endif
    if (i == nIndT +1) then
      nIndT = nIndT +1
    else if (i > nIndT) then
      stop 'ped_add: invalid index'
    endif
        
    pedigree(i) = ind_new
    
  end subroutine ped_add
  
  
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! increase pedigree size by chunk_size_large
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine grow_pedigree()
    type(individual), allocatable :: tmp(:)
    
    allocate(tmp(0:(nInd_max + chunk_size_large)))
    tmp(0:nInd_max) = pedigree
    call move_alloc(tmp, pedigree)  ! from, to
    nInd_max = nInd_max + chunk_size_large 
    ! TODO: also grow other arrays
  end subroutine grow_pedigree
  
  
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! for individual with index iq, get index of dam (p=1) or sire (p=2)
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  pure function get_par(iq, p) 
    integer :: get_par
    integer, intent(IN) :: iq  ! focal (offspring) individual index
    integer, intent(IN) :: p  ! parent sex: 1=dam, 2=sire
    
    get_par = pedigree(iq)%parent(p)
  
  end function get_par
  
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! for individual with index i, update dam (p=1) or sire (p=2) in pedigree
  ! to individual with index par_q
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine set_par(iq, par_q, p)
    integer, intent(IN) :: iq, par_q
    type(individual), pointer :: par
    integer, intent(IN), optional :: p
    integer :: pp

    if (par_q == 0) then
      if (.not. present(p)) then
        stop 'set_par: if par_q is 0, p must be provided'
      else
        pp = p
      endif
    else 
      par => pedigree(par_q)
      if (present(p)) then
        if (any((/1,2/) == par%sex) .and. p/=par%sex) then
          stop 'set_par: parent sex does not match p'
        else
          pp = p
          if (.not. any((/1,2/) == par%sex)) then
            par%sex = p
          endif
        endif
      else
        if (.not. any((/1,2/) == par%sex)) then
          stop 'set_par: parent sex must be 1 or 2, or p provided'
        else
          pp = par%sex
        endif
      endif
    endif
    
    if (pp==1 .or. pp==2) then
      Pedigree(iq)%parent(pp) = par_q
    else
      stop 'set_par: p must be 1 (dam) or 2 (sire)'
    endif
    
    ! TODO: check: does matepair pointer get automatically updated?
    
  end subroutine set_par


  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! read pedigree from file & update pedigree array
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine read_pedigree(FileName) 
    use sqa_fileIO, ONLY: checkFile, IOstat_handler, FileNumRow
    character(len=*), intent(IN) :: FileName
    integer :: i, k, IOerr, nIndP, focal_index, parent_index
    integer, parameter :: no_index=-999
    character(len=nchar_ID) :: focal_id, parent_id(2)
    
    call CheckFile(FileName)
    
    nIndP = FileNumRow(trim(FileName)) -1  ! 1st row = header

    open(unit=103, file=trim(FileName), status="old")
      read(103,*)   ! header                              
      do i=1,nIndP
        read(103, *,IOSTAT=IOerr)  focal_id, parent_id
        if (IOerr /= 0)  call IOstat_handler(IOerr, i, FileName)
        focal_index = get_index(focal_id)
        if (focal_index == no_index) then
          call ped_add(individual(index=nIndT+1, ID=focal_id))
          focal_index = nIndT   ! nIndT updated by ped_add
        endif
        do k=1,2
          parent_index = get_index(parent_id(k))
          if (parent_index == 0)  cycle
          if (parent_index == no_index) then
            call ped_add(individual(index=nIndT+1, ID=parent_id(k)))
            parent_index = nIndT
          endif
          call set_par(focal_index, parent_index, k)
        enddo  
      enddo
    close(103)
    
    
  contains
    ! get pedigree index (=genotype row number, if any) for a specific ID
    pure function get_index(this_ID)
      integer :: get_index
      character(len=nchar_id), intent(IN) :: this_ID
      integer :: i
      
      if (this_ID == 'NA' .or. this_ID=='0') then
        get_index = 0
        return
      endif
      
      get_index = no_index
      do i=1,nIndT
        if (Pedigree(i)%ID == this_ID) then
          get_index = i
          return
        endif
      enddo   
    end function get_index

  end subroutine read_pedigree


  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! initialise vector with all unique dam-sire pairings & resulting offspring 
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine init_matings()
    use sqa_general, ONLY: chunk_size_small
    integer :: i,j,m
    type(matepair), allocatable, target :: tmp(:)
    logical :: pair_exists
    
    allocate(tmp(nIndT))   ! no. pairs < nIndT
         
    m=0
    do i=1,nIndT
      ! do not include founders as offspring of a unique pairing
      if (all(pedigree(i)%parent == 0))  cycle
      if (m==0 .or. any(pedigree(i)%parent == 0)) then  
        ! either parent is NA: always a unique pair
        m = m+1
        tmp(m) = matepair(index=m, parent_q = pedigree(i)%parent, nOffspring=1)
        allocate(tmp(m)%offspring_q(chunk_size_small))
        tmp(m)%offspring_q = 0
        tmp(m)%offspring_q(1) = i
!        pedigree(i)%parent_m = m
      else 
      ! both parents are known: check if pair already exists in list
        pair_exists = .FALSE.
        do j=1,m
          if (all(pedigree(i)%parent == tmp(j)%parent_q)) then
            pair_exists = .TRUE.
            exit
          endif
        enddo
        if (.not. pair_exists) then
          m = m+1
          j = m
          tmp(m) = matepair(index=m, parent_q = pedigree(i)%parent)
          allocate(tmp(m)%offspring_q(chunk_size_small))
          tmp(m)%offspring_q = 0
        endif
        if (pair_exists) then
          if (modulo(tmp(j)%nOffspring, chunk_size_small)==0) then  ! increase size
            call grow_alloc_V(tmp(j)%offspring_q, chunk_size_small)
          endif
        endif       
        tmp(j)%nOffspring = tmp(j)%nOffspring +1
        tmp(j)%offspring_q( tmp(j)%nOffspring ) = i
!        pedigree(i)%parent_m = j
      endif
    enddo
    
    nMatings = m
    allocate(matings(0:nMatings))   
    matings(0) = matepair(index=0)
    matings(1:nMatings) = tmp(1:nMatings)
    deallocate(tmp)
    
    call matings_to_pedigree()
    
  contains
    !~~~~~~~~~~~~~~~~~~~~~~ 
    ! write matings to pedigree individuals
    subroutine matings_to_pedigree()
      integer :: m,k,j,n
    
      ! parents: write nMates & matings_m
      do m=1, nMatings
        do k=1,2
          if (matings(m)%parent_q(k) == 0) cycle
          j = matings(m)%parent_q(k)
          n = pedigree(j)%nMates +1
          pedigree(j)%nMates = n
          if (n==1) then
            allocate(pedigree(j)%matings_m(chunk_size_small))
            pedigree(j)%matings_m = 0
          endif
          if (modulo(n, chunk_size_small)==0) then
            ! increase size of pedigree(j)%matings_m
            call grow_alloc_V(pedigree(j)%matings_m, chunk_size_small)
          endif
          pedigree(j)%matings_m(n) = m
        enddo
      enddo
      
      ! offspring: point to parentpair
      do m=1,nMatings
        do j=1, matings(m)%nOffspring
          pedigree(matings(m)%offspring_q(j))%parentpair => matings(m)
        enddo
      enddo 
    
    end subroutine matings_to_pedigree

  end subroutine init_matings
  
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine grow_alloc_V(V, chunk_size)
    integer, allocatable, intent(INOUT) :: V(:)    ! pointer? class(*)?
    integer, intent(IN) :: chunk_size
    integer, allocatable :: tmp(:)
    
    allocate(tmp(LBOUND(V,1) : (UBOUND(V,1) + chunk_size)))
    tmp(LBOUND(V,1) : UBOUND(V,1)) = V
    call move_alloc(tmp, V)   ! from, to 
  end subroutine grow_alloc_V

end module pedigree_fun



!===============================================================================
!===============================================================================

module check_pedigree
  use sqa_general, ONLY: OcA, AHWE, OHWE, OKA2P, AKAP, OKAP
  use global_variables, ONLY: nSnp, Er, AF, Geno
  implicit none
  private
  public:: clean_pedigree
  
  integer, parameter :: nRel=6, PO_rel = 2
  double precision, allocatable, dimension(:,:,:,:) :: LL_pair

contains
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! Calculate likelihoods for various relationships & store in look-up table
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine init_pairLL
    ! calculate likelihoods for: self, PO, FS, HS/GP/FA, HA/3rd, U
    ! store in look-up table: G_IID1 (-1/0/1/2) + G_IID2 (-1/0/1/2)

    ! when not conditioning on parents of both individuals, LL(FA)=LL(GP)=LL(HS),
    ! and similarly all 3rd degree relationships have same LL (log likelihood)

    ! assume no inbreeding, and otherwise unrelated.    <--- TODO: CHECK/ADD INBRED

    integer :: l, x, y, z, w, v
    double precision :: Tmp(0:2), Tmp2(0:2,0:2), Tmp3(0:2,0:2,0:2)
       

    allocate(LL_Pair(-1:2,-1:2, nSnp, nRel))  ! G_IID1, G_IID2, snp, rel (S/PO/FS/GP/HA/U)
    LL_Pair = 0D0

    do l = 1, nSnp
      do x=-1,2  ! observed genotype IID1
        do y=-1,2    ! observed genotype IID2
   
          ! Self
          forall (z=0:2)  Tmp(z) = OcA(z,x) * OcA(z,y) * AHWE(z,l)
          LL_Pair(x,y,l, 1) = SUM(Tmp)
          
          ! PO
          forall (z=0:2)  Tmp(z) = OKAP(x,z,l) * OcA(z,y) * AHWE(z,l)
          LL_Pair(x,y,l, 2) = SUM(Tmp)
          
          ! FS
          forall (v=0:2, w=0:2)  Tmp2(v,w) = OKA2P(x,v,w) * OKA2P(y,v,w) * AHWE(v,l) * AHWE(w,l)
          LL_Pair(x,y,l, 3) = SUM(Tmp2)
          
          ! HS/GP/FA
          forall (z=0:2, w=0:2)  Tmp2(z,w) = OKAP(x,w,l) * AKAP(w,z,l) * OcA(z,y) * AHWE(z,l)
          LL_Pair(x,y,l, 4)  = SUM(Tmp2)
                   
          ! GGP/HA/3rd
          forall (z=0:2, w=0:2, v=0:2)  Tmp3(z,w,v) = OKAP(x,w,l) * AKAP(w,v,l) * &
            AKAP(v,z,l) * OcA(z,y) * AHWE(z,l)
          LL_Pair(x,y,l, 5)  = SUM(Tmp3)
          
          ! U
          LL_Pair(x,y,l, 6) = OHWE(x,l) * OHWE(y,l)
          
        enddo
      enddo 
    enddo
   
  end subroutine init_pairLL
  
  
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! calculate relatedness probabilities for individuals i + j
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  pure function R_probs(i,j)
    double precision :: R_probs(nRel)
    integer, intent(IN) :: i,j   ! pedigree indices of the two individuals
    integer :: l, r
    double precision :: LL_tmp(nSnp), R_LL(nRel), R_tmp(nRel)
    
    R_LL = 0D0
    do r=1,nRel
      LL_tmp = 0D0    
      do l=1,nSnp
        LL_tmp(l) = LOG(LL_Pair(Geno(i,l), Geno(j,l), l, r))
      enddo
      R_LL(r) = SUM(LL_tmp)
    enddo
  
    ! scale to probabilities summing to 1
    R_tmp = R_LL - MAXVAL(R_LL)   ! minimises problems with rounding to 0
    R_tmp = EXP(R_tmp)
    R_probs =  R_tmp / SUM(R_tmp) 
  end function R_probs  
  
  
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! find & remove low-probability parent-offspring links in the pedigree
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine clean_pedigree(LogFileName, Threshold)
    ! NOTE: only checking genotyped-genotyped parent-offspring pairs
    use global_variables, ONLY : nIndG
    use pedigree_fun, ONLY: pedigree, get_par
    
    character(len=*), intent(IN) :: LogFileName
    real, intent(IN) :: Threshold
    integer :: i, k, j
    double precision :: probs_ij(nRel)
    
    call init_pairLL()
    
    ! write to log which links have been removed, with R_probs
    open(42, file=trim(LogFileName), status='unknown', action='write')
      write(42, '(2a40, a12, 6(5X, a7))') 'id', 'parent_id', 'parent_sex', &
        'prob_S ', 'prob_PO', 'prob_FS', 'prob_GP', 'prob_HA', 'prob_UU'
      do i=1,nIndG  
        do k=1,2   ! dam,sire
          j = get_par(i, k)
          if (j==0 .or. j > nIndG)  cycle   ! parent unknown / not genotyped
          
          probs_ij = R_probs(i,j)
          
          if (probs_ij(PO_rel) < Threshold) then 
            write(42, '(2a40, i12, 6(5X, f7.4))')  Pedigree(i)%ID, Pedigree(j)%ID, k, probs_ij
            pedigree(i)%parent(k) = 0
            ! NOTE: this is done before matingpairs are initialised
          endif
        enddo
      enddo
    close(42)
    
    if (allocated(LL_pair)) deallocate(LL_pair)

  end subroutine clean_pedigree

end module check_pedigree


!===============================================================================
!===============================================================================

module Generations
  use global_variables, ONLY: nIndT, nMatings, quiet
  use pedigree_fun, ONLY: pedigree, Matings
  implicit none
  private

  integer, parameter :: max_gen = 999
  integer, allocatable, dimension(:) :: Gen_down, Gen_up, Parent(:,:)
  integer, allocatable, dimension(:), public ::Gen_rank_down, Gen_rank_up
  
  public :: calc_Gens

contains  
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! rank individuals & matings by generation number, down resp. upwards
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine calc_Gens()
    ! TODO: use functions instead of subroutines?
    use Sort, ONLY: getRank
  
    call mk_parent_matrix()
    
    call calc_gens_down()
     ! from generation numbers per individual to order in which to do peeling
    allocate(Gen_rank_down(1:nIndT))
    Gen_rank_down = getRank(Gen_down(1:nIndT))
    if (.not. quiet)  print *, 'max gen down: ', maxval(gen_down)
    deallocate(Gen_down)
    
    call calc_gens_up()   
    allocate(Gen_rank_up(1:nMatings))
    Gen_rank_up = getRank(Gen_up)
    if (.not. quiet)  print *, 'max gen up: ', maxval(gen_up)
    deallocate(Gen_up)
    
    deallocate(Parent)
    
  end subroutine calc_Gens
  
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! make array with parent indices
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine mk_parent_matrix()
    use pedigree_fun, ONLY: get_par
    integer :: i,p
  
    allocate(Parent(2,nIndT))
    forall (i=1:nIndT, p=1:2)  Parent(p,i) = get_par(i,p)
  
  end subroutine mk_parent_matrix
  
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! generation numbers down: from earliest ancestor to most recent descendant
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
  subroutine calc_gens_down
    integer :: i,g
    
    allocate(Gen_down(0:nIndT))
    Gen_down = max_gen
    Gen_down(0) = 0
    
    do i=1,nIndT
      if (ALL(Parent(:,i)==0))  Gen_down(i) = 0  ! founder
    enddo

    do g = 0,max_gen
      do i=1, nIndT
        if (Gen_down(i) < max_gen)  cycle
        if (ALL(Gen_down(Parent(:,i)) <= g)) then   ! including Parent(i,m)==0
          Gen_down(i) = g+1
        endif    
      enddo
      if (ALL(Gen_down < max_gen))  exit   ! all done
    enddo

  end subroutine calc_gens_down
  
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! generation numbers up: from most recent descendant to earliest ancestor
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! no offspring --> lp_post not used --> generation number irrelevant
  subroutine calc_gens_up()
    integer :: i,g, max_gen_off, m, no, k, mk 
    integer, allocatable, dimension(:) :: Gen_up_i

    allocate(Gen_up_i(0:nIndT))  ! use temporary per-individual array to make things easier
    Gen_up_i(0) = 0
    Gen_up_i = max_gen

    do i=1,nIndT
      if (.not. any(Parent==i))  Gen_up_i(i) = 0  ! individual without offspring
    enddo
    
    ! for each mating, check if all offspring have a generation number yet. 
    ! If so, Gen_up(i) = offspring max +1
    do g=0,max_gen
      do i=1,nIndT
        if (Gen_up_i(i) < max_gen)  cycle   ! already done/no offspring
        max_gen_off = -1
        do k=1, Pedigree(i)%nMates
          mk = Pedigree(i)%matings_m(k)
          no = Matings(mk)%nOffspring
          max_gen_off = MAX(max_gen_off, MAXVAL(Gen_up_i( Matings(mk)%offspring_q(1:no) )) ) 
          if (max_gen_off == max_gen)  exit  ! one or more offspring with not-yet-known Gen (> g)
        enddo
        if (max_gen_off <= g) then
          Gen_up_i(i) = g +1
        endif
      enddo
      if (ALL(Gen_up_i < max_gen)) exit
    enddo
    
    ! generation number for each mating
    allocate(Gen_up(1:nMatings))    
    Gen_up = max_gen
    do m=1, nMatings
      no = Matings(m)%nOffspring
      Gen_up(m) = MAXVAL(Gen_up_i( Matings(m)%offspring_q(1:no) ))       
    enddo
    
    deallocate(Gen_up_i)
  end subroutine calc_gens_up

end module Generations

!===============================================================================
!===============================================================================
module impute_fun
  ! TODO: divide into submodules?
  use global_variables, ONLY: nSnp, Geno, SNP_names, quiet
  use Generations, ONLY: calc_Gens, Gen_rank_down, Gen_rank_up
  ! for prob_ant_post
  use sqa_general, ONLY: ishort, logSumExp, logScale, OcA, AKA2P, AHWE
  use pedigree_fun
  
  implicit none
  private

  logical, public :: do_snpclean, do_impute, do_geno_out, do_probs_out, do_impute_all, do_quick
  real, public :: Threshold_snpclean, Threshold_impute
  double precision, public :: tol
  character(len=3), public :: imp_default
  double precision, allocatable :: Gprob(:,:), Gprob_prev(:,:)
  integer :: N
  ! for prob_ant_post
  integer(kind=ishort), allocatable :: Gl(:)   ! genotypes at SNP l
  double precision, allocatable :: lp_ant(:,:), lp_post(:,:,:)  ! log-probabilities
  double precision :: lOcA(0:2, -1:2), lAKA2P(0:2,0:2,0:2)
  
  public :: clean_n_impute, apply_edits

contains
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! snp cleaning & imputation across all SNPs
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine clean_n_impute(LogFile)

    character(len=*), intent(IN) :: LogFile
    double precision, parameter :: doubt_threshold = 0.49  ! if two genotypes have prob > doubt_threshold, impute as --when-in-doubt
    integer :: l,x, unit_log
    character(len=20) :: nIndT_char
    character(len=500), allocatable :: fmt_gprob_header, fmt_gprob_data
    
    ! allocate once & re-use for each SNP
    allocate(Gprob(0:2,1:nIndT)) 
    allocate(Gprob_prev(0:2,1:nIndT))    
    allocate(Gl(0:nIndT))
    allocate(lp_ant(0:2, 0:nIndT))
    allocate(lp_post(0:2,2, 0:nMatings))
    
    if (do_probs_out) then
      write(nIndT_char, *) nIndT
      fmt_gprob_header = '(a10, a3, 2x, '//trim(nIndT_char)//'a20)'
      fmt_gprob_data   = '(a10, i3, 2x, '//trim(nIndT_char)//'f8.5)'
      open(unit=3, file='geno_probs.txt')
      write(3, fmt_gprob_header)  'SNP', 'G', Pedigree(1:nIndT)%ID  
    endif
    
    if (do_snpclean .or. do_impute) then
      unit_log = 4
      open(unit=unit_log, file=trim(LogFile))
      write(unit_log, '(5(a9,2x), 29x,2(a5,2x),3a9)') 'snp_index', 'snp_name',  &
        'threshold', 'id_index', 'id_name', 'g_in', 'g_out', 'prob_0', 'prob_1', 'prob_2'
    endif
    
    if (do_impute_all) then
      N = nIndT
    else
      N = nIndG
    endif
    
    if (.not. quiet)  print *, 'getting generation numbers ... '
    call calc_Gens()
    
    if (.not. quiet) print *, 'cleaning and/or imputing genotype data... '
    do l=1, nSnp 
      if (.not. quiet .and. mod(l,50)==0) write(*,'(i5, 2x)', advance='no') l
      ! only first iteration, if looping over several thresholds:
      Gl = -1
      Gl(0:nIndG) = Geno(:,l) 
      Gprob = 0D0
      
      call init_gprobs(l)
      
      if (do_quick) then
        call quick_set_ant()
      else
        call peeler(tol)
      endif
      
      if (do_probs_out) then
        do x=0,2
          write(3, fmt_gprob_data) SNP_names(l), x, exp(Gprob(x,:))
        enddo
      endif
     
      ! remove genotyping errors
      if (do_snpclean) then
        call clean_snp(l, Threshold_snpclean)
        call peeler(tol)
      endif 
      
      ! impute
      if (do_impute)  call impute_snp(l, Threshold_impute)
      
    enddo
    if (do_probs_out)  close(3)
    if (do_snpclean .or. do_impute)  close(unit_log) 
    if (.not. quiet)  print *, ''   ! line break
    
    deallocate(Gprob_prev) 
    deallocate(Gprob)
    if (allocated(Gl)) deallocate(Gl)
    if (allocated(lp_ant)) deallocate(lp_ant)
    if (allocated(lp_post)) deallocate(lp_post)
  
  contains
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! initialise arrays
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    subroutine init_gprobs(l)
      integer, intent(IN) :: l  ! SNP number
      integer :: i,m
      
      forall (i=0:nIndT)      lp_ant(:,i)    = log(AHWE(:,l))  
      forall (m=0:nMatings)   lp_post(:,:,m) = 0D0  
      lOcA = log(OcA)
      lAKA2P = log(AKA2P) 
    
    end subroutine init_gprobs
    
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! set genotypes with Gprob <= Threshold to missing
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    subroutine clean_snp(l, Threshold)
      integer, intent(IN) :: l
      real, intent(IN) :: Threshold
      integer :: i
      integer(kind=ishort) :: g_new
      
      do i=1,nIndG
        if (Gl(i) == -1)  cycle
        if (Gprob(Gl(i),i) > log(Threshold))  cycle
 !       Gprob_excl = calc_g_prob(i, excl_self=.TRUE.)  
        ! DO NOT EXLUDE SELF: else doesn't resolve whether parent or offspring
        ! is more likely to have error when they conflict. 
        g_new = -1
        write(unit_log, '(i9,2x,a9,2x,f9.5,2x,i9,2x,a40,2(i5,2x),3f9.5)') l, SNP_names(l), Threshold, &
         i, Pedigree(i)%ID, Gl(i), g_new, exp(Gprob(:,i))
        Gl(i) = g_new
      enddo
    
    end subroutine clean_snp
  
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! impute missing genotypes with Gprob > Threshold
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    subroutine impute_snp(l, Threshold)
      integer, intent(IN) :: l
      real, intent(IN) :: Threshold 
      integer :: i
      integer(kind=ishort) :: g_new, g_common
      
      g_common = -1
      if (imp_default == 'com') then  ! determine most common genotype
        g_common = MAXLOC( G_freq(l), DIM=1, KIND=ishort) -1_ishort
      endif
      
      do i=1,N
        if (Gl(i) /= -1)  cycle
        if (ALL(Gprob(:,i) < log(Threshold)))  cycle
        if (COUNT(Gprob(:,i) >= log(doubt_threshold)) > 1) then   ! implies threshold <= doubt_threshold 
          select case (imp_default)
            case ('het')  
              g_new = 1
            case ('hom')
              if (Gprob(0,i) > Gprob(2,i)) then
                g_new = 0
              else
                g_new = 2
              endif
            case ('com')
              g_new = g_common
          end select
        else
          g_new = MAXLOC(Gprob(:,i), DIM=1, KIND=ishort) -1_ishort
        endif
        write(unit_log, '(i9,2x,a9,2x,f9.5,2x,i9,2x,a40,2(i5,2x),3f9.5)') l, SNP_names(l), Threshold, &
         i, Pedigree(i)%ID, Gl(i), g_new, exp(Gprob(:,i))
        Gl(i) = g_new
      enddo 
    
    end subroutine impute_snp   
      
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! genotype frequencies at SNP l
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    pure function G_freq(l)
      integer, intent(IN) :: l
      double precision :: G_freq(0:2)
      integer :: x

      do x=0,2
        G_freq(x) = COUNT(Geno(1:nIndG,l)==x)/dble(nIndG)
      enddo

    end function G_freq
  
  end subroutine clean_n_impute
  
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! iteratively calculate anterior & posterior probabilities until genotype 
  ! probabilities converge
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine peeler(tol)
    double precision, intent(IN) :: tol  ! tolerance level to declare convergence
    integer, parameter ::  r_max = 50  ! max number of iterations (fail safe)
    integer :: r,x,i,m
       
    do r=1,r_max
      Gprob_prev = Gprob
      
      ! peel down
      do x=1,nIndT
        i = Gen_rank_down(x)
        lp_ant(:,i) = calc_p_ant(i)
      enddo  
      
      ! peel up
      do x=1,nMatings
        m = Gen_rank_up(x)
        lp_post(:,:,m) = calc_p_post(m)
      enddo 
      
      ! calc genotype probabilities
      do i=1,nIndT
        Gprob(:,i) = calc_g_prob(i)
      enddo
 
      ! check for convergence
!      write(*,'(2i4, "  total abs diff: ", f15.5)') l, r, SUM(abs(Gprob - Gprob_prev))
      if (all(abs(Gprob - Gprob_prev) < tol))  exit
    enddo
    
  end subroutine peeler
  
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! set anterior probabilities based on parental genotypes only
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine quick_set_ant()
    integer :: x,i
    
    do x=1,nIndT
      i = Gen_rank_down(x)
      lp_ant(:,i) = quick_p_ant(i)
    enddo
    
    do i=1,nIndT
      Gprob(:,i) = calc_g_prob(i)
    enddo
       
       
  contains
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    function quick_p_ant(i)
      double precision :: quick_p_ant(0:2)
      integer, intent(IN) :: i
      integer :: p, par_q(2),x,y,u
      double precision :: part_par(0:2,2), tmpA(0:2,0:2), tmpB(0:2,0:2)
      
      forall (p=1:2)  par_q(p) = get_par(i,p) 
      forall (p=1:2)  part_par(:,p) = lp_ant(:,par_q(p)) + lOcA(:,Gl(par_q(p)))
      ! combine the parts. y=sire genotype, x=dam genotype
      do y=0,2  ! sire genotype
        forall (x=0:2)  tmpA(:,x) = part_par(x,1) + part_par(y,2) + lAKA2P(:,x,y)
        forall (u=0:2)  tmpB(y,u) = logSumExp(tmpA(u,:))
      enddo
    forall (u=0:2)  quick_p_ant(u) = logSumExp(tmpB(:,u))
    
    ! scale to sum to 1 on non-log scale 
    quick_p_ant = logScale(quick_p_ant)

    end function quick_p_ant

  end subroutine quick_set_ant  
  
  ! TODO? quick_p_post
  
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! apply edits in log file to Geno & write to file
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine apply_edits(editFile, GenoOutFile, GenoOutFormat)
    use sqa_fileIO
    
    character(len=*), intent(IN) :: editFile, GenoOutFile, GenoOutFormat
    integer :: i, l, x, N_edits
    integer(kind=ishort) :: g_in, g_out
    character(len=40) :: dumC(2)
    real :: dumR
    
    N_edits = FileNumRow(trim(editFile)) -1  ! 1st row = header
    
    open (5, file=trim(editFile), action='read')
      read (5,*)  ! header
      do x=1, N_edits
        read (5,*) l, dumC(1), dumR, i, dumC(2), g_in, g_out 
        if (geno(i,l) == g_in) then
          geno(i,l) = g_out
        else
          print *, 'g_in in editfile does not match genotype matrix!'
          print *, 'at line ', x, ': ', dumC, ' in: ', g_in, ' out: ', g_out
        endif
      enddo   
    close(5)  

    ! write to file
    call writeGeno(Geno=transpose(geno), nInd=SIZE(geno,DIM=1), nSnp=SIZE(geno,DIM=2),&
      ID=Pedigree(1:nIndG)%ID, FileName=GenoOutFile, FileFormat=GenoOutFormat)
  
  end subroutine apply_edits


!submodule (impute_fun) prob_ant_pos  !Gprob_mod   
!  use sqa_general, ONLY: ishort, logSumExp, logScale, OcA, AKA2P, AHWE
!  use pedigree_fun
!  implicit none
!contains
   
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! probability that i has genotype u, given all the data
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  function calc_g_prob(i)
    double precision :: calc_g_prob(0:2)
    integer, intent(IN) :: i
       
    calc_g_prob = lp_ant(:,i) + lOcA(:,Gl(i)) + mates_post(i)
    calc_g_prob = logScale(calc_g_prob)  ! scale to sum to 1
    
!    call chk_logprob_OK(calc_g_prob, lbl='calc_g_prob', lbl_i = i)

  end function calc_g_prob
  
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! anterior genotype probability (based on ancestors)
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  function calc_p_ant(i)
    double precision :: calc_p_ant(0:2)
    integer, intent(IN) :: i
    double precision :: part_par(0:2,2), part_sibs(0:2,0:2), tmpA(0:2,0:2), tmpB(0:2,0:2)
    integer :: p, x,y, par_q(2), u

    forall (p=1:2)  par_q(p) = get_par(i,p) 
    do p=1,2
      part_par(:,p) = lp_ant(:,par_q(p)) + lOcA(:,Gl(par_q(p))) + &
        mates_post(par_q(p), q_not=par_q(3-p))
    enddo

    ! full siblings 
    part_sibs = inh_offspr(pedigree(i)%parentpair, q_not=i)   ! TODO: pointer? chk speed & mem use

    ! combine the parts. y=sire genotype, x=dam genotype
    do y=0,2  ! sire genotype
      forall (x=0:2)  tmpA(:,x) = part_par(x,1) + part_par(y,2) + lAKA2P(:,x,y) + part_sibs(x,y)
      forall (u=0:2)  tmpB(y,u) = logSumExp(tmpA(u,:))
    enddo
    forall (u=0:2)  calc_p_ant(u) = logSumExp(tmpB(:,u))
    
    ! scale to sum to 1 on non-log scale 
    calc_p_ant = logScale(calc_p_ant)
    
 !   call chk_logprob_OK(calc_p_ant, lbl='calc_p_ant', lbl_i = i)
     
  end function calc_p_ant
  
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! posterior genotype probability (based on descendants) for mating pair m
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  function calc_p_post(m)
    double precision :: calc_p_post(0:2, 2)  ! genotypes; dam-sire  
    integer, intent(IN) :: m
    type(matepair), pointer :: mp
    double precision :: part_par(0:2,2), part_off(0:2,0:2)
    integer :: p, par_q(2), x,y
    
    mp => Matings(m)
    if (associated(mp)) then
      par_q = mp%parent_q
    else
      par_q = 0
    endif
    do p=1,2  
      part_par(:,p) = lp_ant(:,par_q(p)) + lOcA(:,Gl(par_q(p))) + &
        mates_post(par_q(p), q_not=par_q(3-p))
    enddo
    
    ! offspring
    part_off = inh_offspr(mp)
            
    ! combine parts
    ! NOTE: dam prob sums over possible sire genotypes and vv
    forall (x=0:2)  calc_p_post(x,1) = logSumExp(part_par(:,2) + part_off(x,:))
    forall (y=0:2)  calc_p_post(y,2) = logSumExp(part_par(:,1) + part_off(:,y))
        
    ! scale to sum to 1
    forall (p=1:2)  calc_p_post(:,p) = logScale(calc_p_post(:,p))
    
 !   call chk_logprob_OK(RESHAPE(calc_p_post, (/3*2/)), lbl='calc_p_post', lbl_i = m)
    
  end function calc_p_post
  
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! calc g_i(obs|act) * prod_k(ppost_ik) over mates k
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  function mates_post(i, q_not)
    double precision :: mates_post(0:2)
    integer, intent(IN) :: i  ! focal individual
    integer, intent(IN), optional :: q_not  ! mate to exclude from calculation
    integer :: k,km
    double precision, allocatable :: post_per_mate(:,:)
    type(individual), pointer :: ip
    
    ip => pedigree(i)
    
    allocate(post_per_mate(0:2, ip%nMates))
    post_per_mate = 0D0
    do k=1,ip%nMates    ! if nMates=0, nothing is done.
      km = ip%matings_m(k)    ! index in matings array (type matepair)
      if (present(q_not)) then
        if (any(Matings(km)%parent_q == q_not))  cycle
      endif
      post_per_mate(:,k) = lp_post(:, ip%sex, km)
    enddo    
    mates_post = SUM(post_per_mate, DIM=2)

  end function mates_post
  
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! calc prod_o(inh_o * mates_post)
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  function inh_offspr(mp, q_not)
    double precision :: inh_offspr(0:2,0:2)
    type(matepair), pointer, intent(IN) :: mp
!    integer, intent(IN) :: m
    integer, intent(IN), optional :: q_not  ! offspring to exclude from calculation
    integer :: o,oq,x,y
    double precision, allocatable :: per_off(:,:,:)   ! TODO? allocate once & reuse
    double precision :: off_tmp(0:2)
!    type(matepair), pointer :: mp
    
!    mp => Matings(m)
    
    if (.not. associated(mp)) then  ! no mating (both parents unknown)
!    if (m==0) then
      inh_offspr = 0D0
      return
    endif
    
    allocate(per_off(0:2, 0:2, mp%nOffspring))
    per_off = 0D0
    do o=1,mp%nOffspring
      oq = mp%Offspring_q(o)
      if (present(q_not)) then
        if (oq == q_not)  cycle
      endif
      off_tmp = lOcA(:,Gl(oq)) + mates_post(oq)
      ! sum over possible offspring genotypes; y=sire genotype, x=dam genotype
      forall (y=0:2, x=0:2)  per_off(x,y,o) = logSumExp(off_tmp + lAKA2P(:,x,y))
    enddo
    
    ! multiply across offspring
    inh_offspr = SUM(per_off, DIM=3)
  end function inh_offspr
  
  
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! check if probabilities are valid (debugging tool)
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine chk_logprob_OK(logprob, lbl, lbl_i)
    double precision, intent(IN), dimension(:) :: logprob
    character(len=*), intent(IN) :: lbl
    integer, intent(IN) :: lbl_i
    
    if (any(logprob /= logprob) .or. any(logprob > 2*EPSILON(0d0)) .or. all(logprob < -HUGE(0d0))) then
      print *, ''
      print *, trim(lbl), ' not a valid probability! ', lbl_i
      write(*,'(10f9.4)') logprob
      print *, ''
      stop
    endif
  
  end subroutine chk_logprob_OK
  
end module impute_fun


!===============================================================================
!===============================================================================

program main
  use sqa_fileIO, ONLY: ReadGeno
  use sqa_general, ONLY: PreCalcProbs
  use global_variables
  use pedigree_fun, ONLY: init_pedigree, read_pedigree, init_matings
  use check_pedigree
  use impute_fun
  implicit none
  
 ! integer :: i
  character(len=200) :: GenoFile, PedigreeFile, AFFile, GenoOutFile, LogFile
  character(len=3) :: GenoInFormat, GenoOutFormat
  character(len=nchar_ID), allocatable :: IdV(:)
  real :: Threshold_pedclean
  logical :: do_pedclean
    
  call read_args()   ! set defaults & read command line arguments
  
  if (.not. quiet)  print *, 'reading genotypes ...'
  call ReadGeno(Geno=Geno, ID=IdV, SNP_names=SNP_names, FileName=GenoFile, &
    FileFormat = GenoInFormat, transp=.FALSE.)
  nIndG = SIZE(Geno, DIM=1) -1 ! dim 0:nIndG
  nSnp  = SIZE(Geno, DIM=2)
  
  if (.not. quiet)  print *, 'initializing pedigree ...'
  call init_pedigree(IdV)
 
  if (.not. quiet)  call print_sumstats('IN') 
 
  if (.not. quiet)  print *, 'reading pedigree file ... '
  call read_pedigree(PedigreeFile)
  
  allocate(AF(nSnp))
  AF = getAF(AFFile) 
  call PreCalcProbs(nSnp, Er, AF)
  
  if (do_pedclean) then
    if (.not. quiet)  print *, 'cleaning pedigree ... '   
    call clean_pedigree('log_ped_clean.txt', Threshold_pedclean)
  else
    if (.not. quiet)  print *, 'skipping pedigree cleaning ...'
  endif
  
  call init_matings()
  
  call clean_n_impute(LogFile)
  
  if (.not. quiet)  print *, 'writing new genotypes to file ...'
  if (do_geno_out)  call apply_edits(LogFile, GenoOutFile, GenoOutFormat)
  
  if (.not. quiet)  print *, 'done.'
  
  if (.not. quiet .and. do_geno_out)  call print_sumstats('OUT') 
  
  call deallocall()
  deallocate(IdV)


contains
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !  read in command line arguments
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine read_args()
    use sqa_fileIO, only: valid_formats
    integer :: nArg, i, x
    character(len=32) :: arg, argOption 
    
    ! set defaults
    GenoFile = 'Geno'
    GenoInFormat = 'SEQ'
    PedigreeFile = 'Pedigree.txt'  ! 'ped_griffin.txt'
    AFFile = 'NoFile'
    GenoOutFile = 'Geno_imputed'
    GenoOutFormat = 'SEQ'
    LogFile = 'Imputation_log.txt'
    Er = 0.001
    do_pedclean = .TRUE.
    Threshold_pedclean = 0.05
    do_snpclean = .TRUE.
    Threshold_snpclean = 0.01
    do_impute = .TRUE.
    do_impute_all = .FALSE.
    Threshold_impute = 0.9
    imp_default = 'het'
    tol = 0.0001d0
    do_geno_out = .TRUE.
    do_probs_out = .FALSE.
    do_quick = .FALSE.
    quiet = .FALSE.
  
    nArg = command_argument_count()
    if (nArg > 0) then
      i = 0
      do x = 1, nArg
        i = i+1
        if (i > nArg)  exit
        call get_command_argument(i, arg)
        
        select case (arg)
          case ('--help')
            call print_help()
            stop       
             
          case ('--geno', '--geno-in')  
            i = i+1
            call get_command_argument(i, GenoFile)
            
          case ('--informat', '--inFormat')
            i = i+1
            call get_command_argument(i, GenoInFormat)
            if (.not. any(valid_formats == GenoInFormat)) then
              print *, 'ERROR: informat must be one of: ', valid_formats
              stop
            endif
            
          case ('--pedigree')  
            i = i+1
            call get_command_argument(i, PedigreeFile)
            
          case ('--err')
            i = i+1
            call get_command_argument(i, argOption)
            read(argOption, *)  Er   ! TODO: length 3
            if (Er <= 0.0 .or. Er > 0.5)  stop 'please provide a genotyping error rate --err >0 and <0.5'            
            
          case ('--af', '--maf', '--freq')
            i = i+1
            call get_command_argument(i, AFFile)
            
          case ('--T-pedclean')
            i = i+1
            call get_command_argument(i, argOption)
            read(argOption, *)  Threshold_pedclean
            if (Threshold_Pedclean < 0.0 .or. Threshold_pedclean > 1.0) then
              stop '--T-pedclean must be between 0 and 1, inclusive'
            endif
            
          case ('--no-pedclean')
            do_pedclean = .FALSE.
                                 
          case ('--T-snpclean')
            i = i+1
            call get_command_argument(i, argOption)
            read(argOption, *)  Threshold_snpclean
            if (Threshold_snpclean < 0.0 .or. Threshold_snpclean > 1.0) then
              stop '--T-snpclean must be between 0 and 1, inclusive'
            endif
            
          case ('--no-snpclean')
            do_snpclean = .FALSE.
            
          case ('--T-impute')
            i = i+1
            call get_command_argument(i, argOption)
            read(argOption, *)  Threshold_impute
            if (Threshold_impute < 0.0 .or. Threshold_impute > 1.0) then
              stop '--T-impute must be between 0 and 1, inclusive'
            endif
          
          case ('--no-impute')
            do_impute = .FALSE.   
            
           case ('--impute-all')
            do_impute_all = .TRUE.
            
          case ('--quick')
            do_quick = .TRUE.

          case ('--when-in-doubt')
            i = i+1
            call get_command_argument(i, imp_default)
            if (.not. any((/'het','hom','com'/) == imp_default)) then
              stop '--when-in-doubt must be one of "het", "hom", "com"'
            endif
            
          case ('--tol')
            i = i+1
            call get_command_argument(i, argOption)
            read(argOption, *)  tol          
                     
          case ('--out', '--geno-out')
            i = i+1
            call get_command_argument(i, GenoOutFile)
            
          case ('--outformat', '--outFormat')
            i = i+1
            call get_command_argument(i, GenoOutFormat)
            if (.not. any(valid_formats == GenoOutFormat)) then
              print *, 'ERROR: outFormat must be one of: ', valid_formats
              stop
            endif
            
          case ('--no-geno-out')
            do_geno_out = .FALSE.
            
          case ('--probs-out')
            do_probs_out = .TRUE.
         
          case ('--log')
            i = i+1
            call get_command_argument(i, LogFile)
            
!         case ('--log-in')
            ! TODO

          case ('--quiet')
            quiet = .TRUE.
            
          case default
            print '(2a, /)', 'Unrecognised command-line option: ', arg
            call print_help()
            stop

        end select
      end do
    endif

    if (do_quick) then
 !     do_pedclean = .FALSE.
      do_snpclean = .FALSE.
    endif
 
  end subroutine read_args

  
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine print_help()
    print '(a)',  'Impute genotypes and fix genotyping errors'
    print '(a, /)', 'command-line options:'
    print '(a)',    '  --help               print usage information and exit'
    print '(a)',    '  --geno <filename>    file with genotype data, file extension will be added based on',&
                    '                         --informat. Default: Geno'
    print '(a)',    '  --informat <x>       SEQ: no header, 0/1/2/-9, IDs in column 1; .txt (default)', &  
                    '                       PED: no header, 11/12/22/00, IDs in column 2 of 6 ', &
                    '                        non-SNP columns; .ped (+.map)', &
                    '                       RAW: header, 0/1/2/NA, IDs in column 2 of 6 non-SNP columns; .raw', &
                    '                       LMT: no header, 0/1/2 without spacing, IDs in separate file;', &
                    '                        .geno + .id'                 
    print '(a)',    '  --pedigree <filename>  file with pedigree, with columns id-parent1-parent2.',&
                    '                          Default: Pedigree.txt'
    print '(a)',    '  --err <value>        presumed genotyping error rate',&
                    '                        (Will be estimated from data in future version)'  
    print '(a)',    '  --af <filename>      optional input file with allele frequencies; only relevant',&
                    '                        in combination with --min_prob. Either 1 column and no header,',&
                    '                        or multiple columns with a column MAF, AF, or Frequency',&
                    '                        E.g. output from plink --freq.'                     
    print '(a)',    '  --T-pedclean <value>  Threshold: minimum prob_PO when cleaning the pedigree',&
                    '                         Default: 0.05'
    print '(a)',    '  --no-pedclean        no pedigree cleaning (removal of erroneous assignments)'
    print '(a)',    '  --T-snpclean <value> Threshold: maximum probability of observed genotype to declare',&
                    '                         a genotyping error and replace/set to missing. Default: 0.001'
    print '(a)',    '  --no-snpclean        no removal of probable genotyping errors'                
    print '(a)',    '  --T-impute <value>   threshold: minimum probability of estimated genotype for',&
                    '                        Imputation (Future version will use set of decreasing tresholds)'
    print '(a)',    '  --no-impute          no imputation' 
    print '(a)',    '  --impute-all         impute all individuals in the pedigree (default: only those in genotype file)'
    print '(a)',    '  --quick              base imputation on parent genotypes only, no iterative peeling'  
    print '(a)',    '  --when-in-doubt <x>  inference when 2 genotypes are equally likely (e.g. offspring of',&
                    '                         hom x het): "het" (default), "hom", or "com" (most common',&
                    '                         genotype at that SNP)'  
! TODO: doubt threshold
    print '(a)',    '  --tol <value>        tolerance to determine convergence of genotype probabilities.',&
                    '                        default: 0.0001'                    
    print '(a)',    '  --out <filename>     output file name, extension will be added. Default: Geno_imputed'
    print '(a)',    '  --outformat <x>      same options as for --informat. Default: SEQ'     
    print '(a)',    '  --no-geno-out        no genotype output file with edits (only log, which will be',& 
                    '                         possible to read in to apply the edits'
    print '(a)',    '  --probs-out          write all genotype probabilities to a text file, with 3 rows',&
                    '                         per SNP and 1 column per ID'
    print '(a)',    '  --log <filename>     log file name. Default: Imputation_log.txt'
    print '(a)',    '  --quiet              suppress all messages'
  end subroutine print_help
  
  
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine print_sumstats(lbl)
    use Global_variables, ONLY : nSnp, nIndG
!    use SNP_probs, ONLY : calcAF
    implicit none

    character(len=*), intent(IN), optional :: lbl
!    double precision :: AF(nSnp)
    real :: prop_missing(nSnp)
    
!    AF = calcAF()
    prop_missing = COUNT(Geno(1:nIndG,:)==-1, DIM=1)/float(nIndG)
     
    write(*,*) "========================="
    if (present(lbl))  print *, lbl
    write(*,*) "========================="
    write(*,'("Genotype matrix: ", i6, " SNPs X ", i6, " individuals")')  nSNP, nIndG
    write(*,*) "Genotype frequencies: "
    write(*, '(4a7)') '0','1','2','?'
    write(*,'(4f7.3)') G_freq_all()
    write(*,'("Missingness per SNP: min ", f7.3, "  max ", f7.3 )')  MINVAL(prop_missing), MAXVAL(prop_missing)
    write(*,*) "========================="

  end subroutine print_sumstats 
  
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! genotype frequencies across all SNPs
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  pure function G_freq_all()
    use Global_variables, ONLY : nSnp, nIndG, Geno
    double precision :: G_freq_all(0:3), G_freq_tmp(-1:2)
    integer :: x

    forall (x=-1:2)  G_freq_tmp(x) = COUNT(Geno(1:nIndG,:)==x)/dble(nIndG*nSnp)
    ! make missing count last item in vector
    G_freq_all(0:2) = G_freq_tmp(0:2)
    G_freq_all(3)   = G_freq_tmp(-1)

  end function G_freq_all
  
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! calculate allele frequencies (SNPs in D2)
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  pure function calcAF()
    use sqa_general, ONLY: ishort, ilong
    
    double precision :: calcAF(nSnp)
    integer :: l
    
    calcAF = 1D0
    do l=1,nSnp
      if (ALL(Geno(:,l)==-1)) cycle
      calcAF(l) = dble(SUM(int(Geno(1:,l),kind=ilong), MASK=Geno(1:,l)/=-1))/(COUNT(Geno(1:,l)/=-1)*2)
    enddo
  end function calcAF

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! get allele frequencies, from file or calculated 
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  function getAF(FileName)
    use sqa_fileIO, ONLY: readAF
    
    character(len=*), intent(IN), optional :: FileName
    double precision :: getAF(nSnp)
    double precision, allocatable :: AF_tmp(:)
    logical :: FromFile
    
    AF = 1D0
    if (.not. present(FileName)) then
      FromFile = .FALSE.
    else if (FileName == 'NoFile') then
      FromFile = .FALSE.
    else 
      FromFile = .TRUE.
    endif
    
    if (.not. FromFile) then  
      getAF = calcAF()
    else
      if (.not. quiet) print *, "Reading allele frequencies in "//trim(FileName)//" ..."   
      AF_tmp = readAF(trim(FileName))
      if (SIZE(AF_tmp) /= nSnp) then
        stop "MAF file "//trim(FileName)//" has different number of SNPs than genotype file!"
      else
        getAF = AF_tmp
        deallocate(AF_tmp)
      endif
    endif

  end function getAF
  
end program main


!===============================================================================   

subroutine deallocall
  use sqa_general
  use global_variables
  use pedigree_fun

  ! global_varirables
  if (allocated(Geno))  deallocate(Geno)
  if (allocated(SNP_names))  deallocate(SNP_names) 
  if (allocated(AF)) deallocate(AF)  
  ! pedigree_fun
  if (allocated(pedigree))  deallocate(pedigree)
  if (allocated(matings))  deallocate(matings)
  ! sqa_general
  if (allocated(AHWE)) deallocate(AHWE)
  if (allocated(OHWE)) deallocate(OHWE)
  if (allocated(AKAP)) deallocate(AKAP)
  if (allocated(OKAP)) deallocate(OKAP)

end subroutine deallocall
  
!===============================================================================
! end. 