! clean & impute SNP genotype data, using a pedigree
! each SNP is done separately, no use is made of LD
!
! Jisca Huisman, jisca.huisman@gmail.com
! 2024-02-12
!
! This code is available under GNU General Public License v3
!
! Compilation: 
! gfortran -O3 sqa_general.f90 sqa_fileIO.f90 Imputator.f90 -o imputator
!
!===============================================================================

module global_variables 
  use sqa_general, ONLY: ishort, nchar_ID
  implicit none

  character(len=30), parameter :: version = "Imputator v0.3.2 (28 Feb 2024)"
  integer, parameter :: chunk_size_large = 100, chunk_size_small=10
  integer :: nIndG, nInd_max, nIndT, nSnp, nMatings, nMat_max
  integer(kind=ishort), allocatable :: Geno(:,:)
  integer, allocatable :: indiv2mating(:)
  character(len=nchar_ID), allocatable :: IdV(:), SNP_names(:)
  double precision, allocatable :: AF(:)
  logical :: quiet
  
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! component type of population vector
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  type individual
    integer :: index
    character(len=nchar_ID) :: ID='NA'
    integer :: sex=3  ! 1=female, 2=male, 3=unknown, 4=hermaphrodite
!    integer :: parent_m = 0  !  matingnode with its parents
    integer, dimension(:), pointer :: offspring_m => null()  ! matingnodes with its offspring
    ! TODO: birth year 
  end type individual
  
  type matingnode
    integer :: index
    integer :: parent(2) = 0
    integer, dimension(:), pointer :: offspring => null()  ! pop individual indices of offspring
  end type matingnode
  
end module global_variables 

!===============================================================================
!===============================================================================

module pedigree_fun
  ! TODO: separate source file
  use sqa_general, ONLY: nchar_ID
  use global_variables, ONLY: nIndG, nIndT, nInd_max, nMatings, individual, &
    chunk_size_small, chunk_size_large, nMat_max, matingnode, indiv2mating  ! , mating2parent
  implicit none
  
  type(individual), allocatable, target :: pop(:)  
  type(matingnode), allocatable, target :: pedigree(:)
  integer, allocatable, private :: ped_array(:,:)
  
  ! TODO: set private/public
  
contains
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! get all locations of value in vector  (PACK is slow?)
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  pure function get_locs(vector, lowerb, upperb, value)
    integer, intent(IN) :: lowerb, upperb, value
    integer, intent(IN) :: vector(lowerb:upperb)
    integer, allocatable :: get_locs(:), v_indx(:)
    integer :: x
    
    allocate(v_indx(lowerb:upperb))
    v_indx = (/ (x, x=lowerb, upperb) /)
    get_locs = PACK(v_indx, MASK = vector == value)
  
  end function get_locs

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! initialise population with individuals read from genotype file
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine init_pop(IdV)  
    character(len=nchar_ID), intent(IN) :: IdV(nIndG)
    integer :: i
    
    allocate(pop(0:nIndG)) 
    pop(0) = individual(index=0, ID='NA')
    allocate(pop(0)%offspring_m(0))
    do i=1,nIndG
      pop(i) = individual(index=i, ID=IdV(i))
    enddo   
    nInd_max = nIndG
    nIndT = nIndG      
    allocate(ped_array(2,1:nInd_max))
    ped_array = 0
  end subroutine init_pop

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! read pedigree from file & init parent array
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine read_pedigree(FileName) 
    ! TODO: speed up
    use sqa_fileIO, ONLY: checkFile, IOstat_handler, FileNumRow
    
    character(len=*), intent(IN) :: FileName
    integer :: i,x,k, ios, nIndP, par_i, r, nInd_R(3)
    integer, parameter :: no_index=-999
    character(len=nchar_ID), allocatable :: names_pedigree(:,:)
    
                                           
    call CheckFile(FileName)    
    nIndP = FileNumRow(trim(FileName)) -1  ! 1st row = header
    allocate(names_pedigree(3,nIndP))
    
    open(unit=103, file=trim(FileName), status="old")
      read(103,*)   ! header                              
      do x=1,nIndP
        read(103, *,IOSTAT=ios)  names_pedigree(:,x)
        if (ios /= 0)  call IOstat_handler(ios, x, FileName)
      enddo
    close(103)
    
    ! add parents of genotyped individuals only + grandparents (= parents of those added in R1)
    nInd_R = nIndG
    do r=1,2
      do x=1,nIndP
        i = get_index(names_pedigree(1,x))
        if (i <= 0 .or. i>nInd_R(r))  cycle  ! not genotyped (R1) / not parent-of-genotyped (R2)
        do k=1,2
          par_i = get_index(names_pedigree(1+k,x))
          if (par_i == 0)  cycle
          if (par_i == no_index) then  
            call pop_add(individual(index=nIndT+1, ID=names_pedigree(1+k,x), sex=k))
            par_i = nIndT
          endif
          ped_array(k,i) = par_i
          call set_sex(par_i, k) 
        enddo
      enddo
      nInd_R(r+1) = nIndT   ! nInd + individuals added in round 1 
    enddo
   
    ! print *, 'nInd_R: ', nInd_R
    ! print *, '# parents of genotyped indivs: ', nIndG, COUNT(ped_array(1,1:nIndG)/=0), COUNT(ped_array(2,1:nIndG)/=0)
    ! print *, '# parents of non-genotyped indivs: ', nInd_R(2)-nIndG, &
      ! COUNT(ped_array(1,(nIndG+1):nInd_R(2))/=0), COUNT(ped_array(2,(nIndG+1):nInd_R(2))/=0)
    ! print *, '# parents of those: ', nInd_R(3)-nInd_R(2), &
      ! COUNT(ped_array(1,(nInd_R(2)+1):nInd_R(3))/=0), COUNT(ped_array(2,(nInd_R(2)+1):nInd_R(3))/=0)
   
    deallocate(names_pedigree)
  
                       
  contains
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! get pedigree index (=genotype row number, if any) for a specific ID
    pure function get_index(this_ID)
      integer :: get_index
      character(len=nchar_id), intent(IN) :: this_ID
      integer :: i
      
      if (this_ID == 'NA' .or. this_ID=='0') then
        get_index = 0
      else 
        get_index = no_index
        if (any(pop%ID == this_id)) then    
          do i=1,nIndT
            if (pop(i)%ID == this_ID) then
              get_index = i
              return
            endif
          enddo
        endif
      endif
      
    end function get_index
    
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! add new individual to the population
    subroutine pop_add(ind_new)
      type(individual), intent(IN) :: ind_new  
      integer :: i
      
      i = ind_new%index
      if (i >= nInd_max) then
        call grow_pop()
      endif
      
      if (i == nIndT +1) then
        nIndT = nIndT +1
      else if (i > nIndT) then
        stop 'pop_add: invalid index'
      endif
      
      pop(i) = ind_new
      
    end subroutine pop_add
   
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! increase population size by chunk_size_large
    subroutine grow_pop()
      type(individual), allocatable :: tmp_pop(:)
      integer, allocatable :: tmp_ped_array(:,:)

      allocate(tmp_pop(0:(nInd_max + chunk_size_large)))
      tmp_pop(0:nInd_max) = pop
      call move_alloc(tmp_pop, pop)  ! from, to   
      
      allocate(tmp_ped_array(2,1:(nInd_max + chunk_size_large)))
      tmp_ped_array = 0
      tmp_ped_array(:,1:nInd_max) = ped_array
      call move_alloc(tmp_ped_array, ped_array)
      
      nInd_max = nInd_max + chunk_size_large

    end subroutine grow_pop
    
  end subroutine read_pedigree
  
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! initialise pedigree & indiv2mating arrays
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine init_pedigree()
    integer, allocatable :: m2p(:,:)  ! temp stand-in for mating2parent
    integer :: i,m, par(2),p, n_off
    logical :: pair_exists
    integer, allocatable :: i_indx(:), m_indx(:) 
    
    allocate(m2p(2,0:nIndT))  ! max. no unique pairs < nIndT
    m2p = 0
    allocate(indiv2mating(0:nIndT))  
    indiv2mating = 0 
    
    if (.not. allocated(ped_array)) then  ! no pedigree read from file
      allocate(ped_array(2,1:nIndG))
      ped_array = 0
    endif
    
    nMatings = 0  
    do i=1, nIndG
      par = ped_array(:,i)
      if (all(par==0))  cycle  ! mating=0
      if (any(par==0)) then    ! always unique mating
        nMatings = nMatings +1
        m2p(:,nMatings) = par
        m = nMatings
      else
        pair_exists = .FALSE.
        do m=1,nMatings
          if (all(par == m2p(:,m))) then
            pair_exists = .TRUE.
            exit
          endif
        enddo
        if (.not. pair_exists) then
          nMatings = nMatings +1
          m2p(:,nMatings) = par
          m = nMatings
        endif
      endif
!        pop(i)%parent_m = m
      indiv2mating(i) = m
    enddo
  
    ! create pedigree with matingnodes: 2 parents + >=1 offspring
    nMat_max = nMatings
    allocate(i_indx(0:nIndT))
    i_indx = (/ (i,i=0,nIndT) /)
    allocate(pedigree(0:nMat_max))  ! allocated by init_pop  (also needed when .not. do_pedigree)   
    pedigree(0) = matingnode(index=0)
    allocate(pedigree(0)%offspring(0))
    do m=1,nMatings
      pedigree(m) = matingnode(index=m, parent = m2p(:,m))
      n_off = COUNT(indiv2mating == m)
      allocate(pedigree(m)%offspring(n_off))
      pedigree(m)%offspring = PACK(i_indx, MASK = indiv2mating==m)
    enddo
    deallocate(m2p)

    ! add matingnode indices to pop individual nodes, for quicker look-up
    allocate(m_indx(0:nMatings))
    m_indx = (/ (i,i=0,nMatings) /)
    do i=1,nIndT
      p = pop(i)%sex
      if (p==3) then
        n_Off = 0
      else
        n_off = COUNT(pedigree%parent(p) == i)
      endif
      allocate(pop(i)%offspring_m(n_off))
      if (n_off>0)  pop(i)%offspring_m = PACK(m_indx, MASK = pedigree%parent(p)==i)
    enddo
    
    deallocate(ped_array)
  end subroutine init_pedigree 
  
  
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! check if sex is consistent across pedigree; set sex if currently unknown
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine set_sex(fcl_i, p)  
    integer, intent(IN) :: fcl_i, p  ! p=supposed sex
    integer :: cur_sex
    
    if (.not. any((/1,2/) == p))  stop 'set_sex: p must be 1 (dam) or 2 (sire)'   
    cur_sex = pop(fcl_i)%sex
    if (any((/1,2/) == cur_sex) .and. p/=cur_sex)  stop 'set_par: parent sex does not match p'
    if (cur_sex==3)  pop(fcl_i)%sex = p    
  end subroutine set_sex   

  
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! for individual with index fcl_i, update dam (p=1) or sire (p=2) in pedigree
  ! to individual with index par_i
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine set_par(fcl_i, par_i, p)
    integer, intent(IN) :: fcl_i, par_i
    integer, intent(IN), optional :: p
    integer :: pp, par_sex, m_old, m_new, nOff,j, curpar(2)
    logical :: pair_exists

    if (par_i == 0) then
      if (.not. present(p))  stop 'set_par: if par_i is 0, p must be provided'
      pp = p
    else 
      par_sex = pop(par_i)%sex 
      ! check if sex matches
      if (present(p)) then
        if (any((/1,2/) == par_sex) .and. p/=par_sex)  stop 'set_par: parent sex does not match p'
        pp = p
        if (par_sex==3)  pop(par_i)%sex = p
      else
        if (.not. any((/1,2/) == par_sex))  stop 'set_par: parent sex must be 1 or 2, or p provided'
        pp = par_sex
      endif
    endif
    if (.not. any((/1,2/) == pp))  stop 'set_par: p must be 1 (dam) or 2 (sire)'
       
!    m_old = pop(fcl_i)%parent_m
    m_old = indiv2mating(fcl_i)
    curpar = pedigree(m_old)%parent
    m_new = m_old
    if (m_old == 0) then
      if (par_i/=0) then
        nMatings = nMatings +1
        if (nMatings >= nMat_max)  call grow_pedigree()
        m_new = nMatings   
      endif  ! else nothing changes
    else if (curpar(3-pp)/=0 .and. curpar(pp)==0) then  ! there is a co-parent
      m_new = m_old
    else if (curpar(3-pp)==0 .and. curpar(pp)/=0) then
      if (par_i==0) then
        m_new = 0
 !       deallocate(pedigree(m_old)%offspring)
        pedigree(m_old)%offspring => null()
        allocate(pedigree(m_old)%offspring(0))
 !       call remove_mating(m_old)  ! remove m_old & shift all subsequent ones in pedigree  TODO
      endif
    else ! all(mating2parent(:, m)/=0)
      nOff = SIZE(pedigree(m_old)%offspring)  ! COUNT(pop%parent_m == m_old)
      ! if nOff==1 no full siblings: change parents of mating
      ! with full siblings: move individual to a new/different mating
      if (nOff > 1) then  
        pair_exists = .FALSE.
        if (par_i/=0) then
          do j=1,nMatings
            if (pedigree(j)%parent(pp)==par_i .and. pedigree(j)%parent(3-pp)==curpar(3-pp)) then
              pair_exists = .TRUE.
              exit
            endif
          enddo
          if (pair_exists) m_new = j
        endif
        if (.not. pair_exists) then  ! make a new mating pair
          nMatings = nMatings +1
          if (nMatings >= nMat_max)  call grow_pedigree()
          m_new = nMatings
          pedigree(nMatings)%parent(3-pp) = curpar(3-pp)  ! bring old co-parent along
        endif
      endif
    endif
    
    if (m_new/=0)  pedigree(m_new)%parent(pp) = par_i
!    pop(fcl_i)%parent_m = m_new 
    indiv2mating(fcl_i) = m_new
    if (m_new == nMatings) then
      allocate(pedigree(m_new)%offspring(1))
      pedigree(m_new)%offspring = fcl_i
 !     call add_mating(par_i, m_new)
     if (par_i/=0) then
        pop(par_i)%offspring_m = get_locs(vector=pedigree%parent(pp), lowerb=0, &
            upperb=nMatings, value=par_i)
      endif
    else
      if (m_new/=0) then 
        ! call add_offspring(m_new, fcl_i)    ! TODO: speed test
        pedigree(m_new)%offspring = get_locs(vector=indiv2mating, lowerb=0, upperb=nIndT, value=m_new)
      endif     
    endif
    
    ! THIS IS A PAIN, KEEPING TRACK OF THINGS IN TWO DIFFERENT ARRAYS
    ! USE POINTERS; SHOULD GET UPDATED AUTOMAGICALLY? OR SOME OTHER WAY THAT IS 
    ! FAST & CLEAR & LOW MEMORY
     
  contains
    ! !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! ! add individual fcl_i as offspring to matingnode m
    ! subroutine add_offspring(m, fcl_i)
      ! integer, intent(IN) :: m, fcl_i
      ! integer :: nOff
      ! integer, allocatable :: off_tmp(:)
      
      ! nOff = SIZE(pedigree(m)%offspring)
 ! !     pedigree(m)%offspring(nOff +1) = fcl_i
      ! allocate(off_tmp(nOff+1))
      ! off_tmp(1:nOff) = pedigree(m)%offspring
      ! off_tmp(nOff+1) = fcl_i
! !      move_alloc(off_tmp, pedigree(m)%offspring)  ! Error: Unclassifiable statement
      ! deallocate(pedigree(m)%offspring)
      ! allocate(pedigree(m)%offspring(nOff+1))
      ! pedigree(m)%offspring = off_tmp
      ! deallocate(off_tmp)
    ! end subroutine add_offspring
    
    
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! increase size of pedigree by another chunk
    subroutine grow_pedigree()
      type(matingnode), allocatable :: tmp(:)

      allocate(tmp(0:(nMat_max + chunk_size_small)))
      tmp(0:nMat_max) = pedigree
      call move_alloc(tmp, pedigree)  ! from, to
      nMat_max = nMat_max + chunk_size_small     
    end subroutine grow_pedigree

  end subroutine set_par
  
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! for individual with index fcl_i, get dam + sire
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    pure function get_par(fcl_i) result(par)
      integer, intent(IN) :: fcl_i
      integer :: par(2), m
      
      m = indiv2mating(fcl_i)
      par = pedigree(m)%parent

    end function get_par
  
  ! no other getters: big decrease in runtime
  
    
  
  ! !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! ! keep only genotyped individuals + their parents (+ grandparents?)
  ! !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! subroutine prune_pedigree()
    ! logical, allocatable :: keep_i(:), keep_m(:)
    ! integer :: i,m,j
    ! type(individual), allocatable :: tmp_pop(:)
      
    ! allocate(keep_m(0:nMatings))
    ! keep_m = .FALSE.
    ! keep_m(0) = .TRUE.
    ! do m=1,nMatings
      ! if (any(pedigree(m)%offspring <= nIndG))  keep_m(m) = .TRUE.
    ! enddo
    
 ! !   print *, 'nMatings: ', nMatings, ' keep: ', COUNT(keep_m)
    
    ! ! keep if non-genotyped individual has at least 1 non-genotyped offspring
    ! allocate(keep_i((nIndG+1):nIndT))
    ! keep_i = .FALSE.
    ! do i=nIndG+1, nIndT
      ! if (any(keep_m .and. (pedigree%parent(1) == i .or. pedigree%parent(2) == i))) then
        ! keep_i(i) = .TRUE.
      ! endif    
    ! enddo   
    
 ! !   print *, 'n non-genotyped: ', nIndT-nIndG, ' keep: ', COUNT(keep_i)
    
    ! allocate(tmp_pop(0:(nIndG+COUNT(keep_i))))
  
  
    
  
  ! end subroutine prune_pedigree
  
end module pedigree_fun


!===============================================================================
!===============================================================================

module check_pedigree
  use sqa_general, ONLY: OcA, AHWE, OHWE, OKA2P, AKAP, OKAP
  use global_variables, ONLY: nSnp, AF, Geno
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
  subroutine clean_pedigree(EditsFileName, Threshold)
    ! NOTE: only checking genotyped-genotyped parent-offspring pairs
    use global_variables, ONLY : nIndG
    use pedigree_fun, ONLY: pop, get_par, set_par
    
    character(len=*), intent(IN) :: EditsFileName
    real, intent(IN) :: Threshold
    integer :: i, k, par(2)
    double precision :: probs_ip(nRel)
    
    call init_pairLL()
    
    ! write to log which links have been removed, with R_probs
    open(42, file=trim(EditsFileName), status='unknown', action='write')
      write(42, '(2a40, a12, 6(5X, a7))') 'id', 'parent_id', 'parent_sex', &
        'prob_S ', 'prob_PO', 'prob_FS', 'prob_GP', 'prob_HA', 'prob_UU'
      do i=1,nIndG  
        par = get_par(i)
        do k=1,2   ! dam,sire
          if (par(k)==0 .or. par(k) > nIndG)  cycle   ! parent unknown / not genotyped        
          probs_ip = R_probs(i,par(k))
          if (probs_ip(PO_rel) < Threshold) then 
            write(42, '(2a40, i12, 6(5X, f7.4))')  pop(i)%ID, pop(par(k))%ID, k, probs_ip
            call set_par(fcl_i=i, par_i=0, p=k)
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
  use pedigree_fun, ONLY: pop, pedigree
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
    integer :: i
  
    allocate(Parent(2,1:nIndT))
    forall (i=1:nIndT)  Parent(:,i) = get_par(i)

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
      if (all(parent(:,i)==0))  Gen_down(i) = 0  ! founder  
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
!    use pedigree_fun, ONLY: pop, pedigree 
    integer :: i,g, max_gen_off, m, k
    integer, allocatable, dimension(:) :: Gen_up_i
    integer, pointer :: offspring(:), mates(:)

    allocate(Gen_up_i(0:nIndT))  ! use temporary per-individual array to make things easier
    Gen_up_i(0) = 0
    Gen_up_i = max_gen

    do i=1,nIndT
      if (.not. any(parent==i))  Gen_up_i(i) = 0  ! individual without offspring
    enddo
    
    ! for each individual, check if all offspring have a generation number yet. 
    ! If so, Gen_up(i) = offspring max +1
    do g=0,max_gen
      do i=1, nIndT
        if (Gen_up_i(i) < max_gen)  cycle   ! already done/no offspring
        max_gen_off = -1
        mates => pop(i)%offspring_m 
        do k=1,SIZE(mates)         
          offspring => pedigree(mates(k))%offspring 
          max_gen_off = MAX(max_gen_off, MAXVAL(Gen_up_i( offspring )) )
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
      offspring => pedigree(m)%offspring 
      Gen_up(m) = MAXVAL(Gen_up_i( offspring ))       
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

  logical, public :: do_snpclean, do_impute, do_geno_out, do_probs_out, &
    do_impute_all, do_quick
  real, public :: Threshold_snpclean, Threshold_impute
  double precision, public :: tol
  character(len=3), public :: imp_default
  double precision, allocatable :: Gprob(:,:), Gprob_prev(:,:)
  integer :: N
  ! for prob_ant_post
  integer(kind=ishort), allocatable :: Gl(:)   ! genotypes at SNP l
  double precision, allocatable :: lp_ant(:,:), lp_post(:,:,:)  ! log-probabilities
  double precision :: lOcA(0:2, -1:2), lAKA2P(0:2,0:2,0:2), lAHWE(0:2)
  
public :: clean_n_impute, apply_edits

contains
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! snp cleaning & imputation across all SNPs
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine clean_n_impute(EditsFile)  ! , do_pedigree
    use sqa_general, ONLY: printt

    character(len=*), intent(IN) :: EditsFile
 !   logical, intent(IN) :: do_pedigree
    double precision, parameter :: doubt_threshold = 0.49 !log(1d0/3) !0.49  ! if two genotypes have prob > doubt_threshold, impute as --when-in-doubt
    integer :: l,x, unit_log
    character(len=20) :: nIndT_char
    character(len=500), allocatable :: fmt_gprob_header, fmt_gprob_data
!    integer :: i, max_mates_per_indiv, n_mates_i
       
    ! allocate once & re-use for each SNP
    allocate(Gprob(0:2,1:nIndT)) 
    allocate(Gprob_prev(0:2,1:nIndT))    
    allocate(Gl(0:nIndT))
    allocate(lp_ant(0:2, 0:nIndT))
    allocate(lp_post(0:2,2, 0:nMatings))

    ! max_mates_per_indiv = 0
    ! do i=1,nIndT
      ! n_mates_i = SIZE(pop(i)%offspring_m)
      ! if (n_mates_i >  max_mates_per_indiv)   max_mates_per_indiv = n_mates_i
    ! enddo   
    ! allocate(post_per_mate(0:2, max_mates_per_indiv))
    ! allocate(im_indx(0:2, max_mates_per_indiv))
    ! forall (x=0:2)  im_indx(x,:) = (/ (i, i=1,max_mates_per_indiv) /)

    
    if (do_probs_out) then
      write(nIndT_char, *) nIndT
      fmt_gprob_header = '(a10, a3, 2x, '//trim(nIndT_char)//'a20)'
      fmt_gprob_data   = '(a10, i3, 2x, '//trim(nIndT_char)//'f8.5)'
      open(unit=3, file='geno_probs.txt')
      write(3, fmt_gprob_header)  'SNP', 'G', pop(1:nIndT)%ID  
      deallocate(fmt_gprob_header)
    endif
    
    if (do_snpclean .or. do_impute) then
      unit_log = 4
      open(unit=unit_log, file=trim(EditsFile))
      write(unit_log, '(5(a9,2x), 29x,2(a5,2x),3(3a9,2x))') 'snp_index', 'snp_name',  &
        'threshold', 'id_index', 'id_name', 'g_in', 'g_out', 'prob_0', 'prob_1', 'prob_2', &
        'p_ant_0', 'p_ant_1', 'p_ant_2', 'p_post_0', 'p_post_1', 'p_post_2'
    endif
    
    if (do_impute_all) then
      N = nIndT
    else
      N = nIndG
    endif
    
    if (.not. quiet)  call printt('getting generation numbers ... ')
    call calc_Gens()
    
    if (.not. quiet) call printt('cleaning and/or imputing genotype data... ')
    
    do l=1, nSnp 
      if (.not. quiet .and. mod(l,50)==0) call printt(text='l=', int=l)  ! write(*,'(i5, 2x)', advance='no') l
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
    
    deallocate(Gprob_prev) 
    deallocate(Gprob)
    if (allocated(Gl)) deallocate(Gl)
    if (allocated(lp_ant)) deallocate(lp_ant)
    if (allocated(lp_post)) deallocate(lp_post)
    if(allocated(fmt_gprob_data)) deallocate(fmt_gprob_data)
    deallocate(Gen_rank_down)
    deallocate(Gen_rank_up)
  
  contains
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! initialise arrays
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    subroutine init_gprobs(l)
      integer, intent(IN) :: l  ! SNP number
      integer :: i,m,x, par(2)
      double precision :: lOHWE(0:2,-1:2)
      
      lOcA = log(OcA)
      lAKA2P = log(AKA2P)
      lAHWE = log(AHWE(:,l))
      
      forall (x=0:2)  lOHWE(x,:) = lAHWE(x) + lOcA(x,:)
      forall (x=-1:2) lOHWE(:,x) = logScale(lOHWE(:,x))      
      forall (i=0:nIndT) lp_ant(:,i)   = lOHWE(:,Gl(i))     
      do m=0, nMatings
        par = pedigree(m)%parent
        forall (i=1:2)  lp_post(:,i,m) = lOHWE(:,Gl(par(i)))
      enddo  
      
    end subroutine init_gprobs
    
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! set genotypes with Gprob <= Threshold to missing
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    subroutine clean_snp(l, Threshold)
      integer, intent(IN) :: l
      real, intent(IN) :: Threshold
      integer :: i!,x
      integer(kind=ishort) :: g_new
      double precision :: p_not_err, i_ap(2,0:2), lp_err!, lp_not_err, ref_unif
        
      p_not_err = (OcA(0,0) + OcA(1,1) + OcA(2,2))/3d0
      lp_err = log(1-p_not_err)
 !     lp_not_err = log(p_not_err)
   !     forall (x=0:2)  ref_lp(x) = log(1-OcA(x,x))
 !     ref_unif = log(1d0/3)   ! too arbitrary?   ! = new doubt_threshold
 !     write(*, '("ref_lp:  ", 3f12.5)')  lp_err, lp_not_err, ref_unif
      ! print *, ''

      do i=1,nIndG
        ! TODO: loop over matings instead?
        if (Gl(i) == -1)  cycle
        i_ap(1,:) = lp_ant(:,i)
        i_ap(2,:) = logScale(mates_post(i))
        
        ! if (i==79) then
          ! write(*,'(i5, 2x, a12, 2x, i5, ";", i5, ";", 10i5)') i, pop(i)%ID, gl(i), &
            ! indiv2mating(i), pop(i)%offspring_m
          ! write(*, '("ant:  ", 3f10.4)') lp_ant(:,i)
          ! write(*, '("post: ", 3f10.4)') i_ap(2,:)
          ! write(*, '("gpr:  ", 3f10.4)') Gprob(:,i)
          ! print *, ''
          ! write(*, '("ant:  ", 3f9.5)') exp(lp_ant(:,i))
          ! write(*, '("post: ", 3f9.5)') exp(i_ap(2,:))
          ! write(*, '("gpr:  ", 3f9.5)') exp(Gprob(:,i))
          ! print *, ''
        ! endif        
        
        if (Gprob(Gl(i),i) < lp_err) then  !log(Threshold))  then
   !       Gprob_excl = calc_g_prob(i, excl_self=.TRUE.)  
          ! DO NOT EXLUDE SELF: else doesn't resolve whether parent or offspring
          ! is more likely to have error when they conflict. 
    !      if (any(Gprob(:,i) >= lp_not_err)) then
    !        g_new = MAXLOC(Gprob(:,i), DIM=1, KIND=ishort) -1_ishort
    !      else
            g_new = -1
     !     endif
          write(unit_log, '(i9,2x,a9,2x,f9.5,2x,i9,2x,a40,2(i5,2x),3(3f9.5,2x))') l, SNP_names(l), &
           Threshold, i, pop(i)%ID, Gl(i), g_new, exp(Gprob(:,i)), exp(lp_ant(:,i)), exp(i_ap(2,:))
          Gl(i) = g_new
  !      else if (any(i_ap(:,Gl(i)) < lp_err) .and. all(i_ap(:,Gl(i)) < ref_unif)) then
  !        g_new = -1  !Gl(i)
  !        write(unit_log, '(i9,2x,a9,2x,f9.5,2x,i9,2x,a40,2(i5,2x),3(3f9.5,2x))') l, SNP_names(l), &
  !         Threshold, i, pop(i)%ID, Gl(i), g_new, exp(Gprob(:,i)), exp(lp_ant(:,i)), exp(i_ap(2,:))  
  !        Gl(i) = g_new           
        endif       
      enddo
    
    end subroutine clean_snp
  
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! get posterior prob for individual i
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    pure function get_post(i)
      double precision :: get_post(3)
      integer, intent(IN) :: i
      integer :: isex
      
      isex = pop(i)%sex
      if (isex==3) then
        get_post = lp_post(:,1,0)
      else
        get_post = lp_post(:,isex,indiv2mating(i))
      endif
    end function get_post 
  
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
        if (COUNT(Gprob(:,i) >= doubt_threshold) > 1) then    
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
         write(unit_log, '(i9,2x,a9,2x,f9.5,2x,i9,2x,a40,2(i5,2x),3(3f9.5,2x))') l, SNP_names(l), &
         Threshold, i, pop(i)%ID, Gl(i), g_new, exp(Gprob(:,i)), exp(lp_ant(:,i)), exp(get_post(i))
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
    integer, parameter ::  r_max = 20  ! max number of iterations (fail safe)
    integer :: r,x,i,m
       
    do r=1,r_max
      Gprob_prev = Gprob
      
!       call timestamp2()
!      print *, 'go'                
      ! peel down
      do x=1,nIndT
        i = Gen_rank_down(x)
        lp_ant(:,i) = calc_p_ant(i)
      enddo  
      
!      call timestamp2()
 !    print *, r, 'lp_ant :', SUM(lp_ant)                 
      ! peel up
      do x=1,nMatings
        m = Gen_rank_up(x)
        lp_post(:,:,m) = calc_p_post(m)
      enddo 
      
!       call timestamp2()
  !    print *, r, 'lp_post :', SUM(lp_post)                                    
      ! calc genotype probabilities
      do i=1,nIndT
        Gprob(:,i) = calc_g_prob(i)
      enddo
 
      ! check for convergence
  !    write(*,'(i4, "  total abs diff: ", f15.5)') r, SUM(abs(Gprob - Gprob_prev))
 !     print *, ''
      if (all(abs(Gprob - Gprob_prev) < tol))  exit
    enddo
    
   ! contains
    ! subroutine timestamp2()
      ! character(len=10) :: time
    
      ! call date_and_time(TIME=time)
      ! write(*, '(a10,1x)', advance='no')  time
       
    ! end subroutine timestamp2
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
      integer :: p, par_q(2),x,y,u, m
      double precision :: part_par(0:2,2), tmpA(0:2,0:2), tmpB(0:2,0:2)
      
      m = indiv2mating(i)
      ! if (m==0) then
        ! calc_p_ant = lAOHWE(:,Gl(i))  ! internal compiler error
        ! return
      ! endif
      par_q = pedigree(m)%parent
      forall (p=1:2)  part_par(:,p) = lp_ant(:,par_q(p)) + lOcA(:,Gl(par_q(p)))
      ! combine the parts. y=sire genotype, x=dam genotype
      do y=0,2 
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
    use global_variables, ONLY: IdV
    
    character(len=*), intent(IN) :: editFile, GenoOutFile, GenoOutFormat
    integer :: i, l, x, N_edits
    integer(kind=ishort) :: g_in, g_out
    character(len=40) :: dumC(2)
    real :: dumR
    
    call CheckFile(trim(editFile))    
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
          stop  
        endif
      enddo   
    close(5)  

    ! write to file  (excl indiv 0)
    call write_geno(Geno=transpose(geno(1:,:)), nInd=SIZE(geno,DIM=1)-1, nSnp=SIZE(geno,DIM=2),&
      ID=IdV, FileName=GenoOutFile, FileFormat=GenoOutFormat)
  
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
       
    calc_g_prob = lp_ant(:,i) + lOcA(:,Gl(i)) + logScale(mates_post(i))
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
    integer :: p, x,y, par_q(2), u, m

    m = indiv2mating(i)
    if (m==0) then
      calc_p_ant = lAHWE
      return
    endif
    par_q = pedigree(m)%parent  
    do p=1,2
      part_par(:,p) = lp_ant(:,par_q(p)) + lOcA(:,Gl(par_q(p))) + &
        mates_post(par_q(p), q_not=par_q(3-p))
    enddo

    ! full siblings   
    part_sibs = inh_offspr(m, q_not=i)

    ! combine the parts. y=sire genotype, x=dam genotype
    do y=0,2  
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
    double precision :: part_par(0:2,2), part_off(0:2,0:2)
    integer :: p, x,y, par_q(2)
    
    par_q = pedigree(m)%parent  ! pointer?
    do p=1,2  
      part_par(:,p) = lp_ant(:,par_q(p)) + lOcA(:,Gl(par_q(p))) + &
        mates_post(par_q(p), q_not=par_q(3-p))
    enddo
    
    ! offspring
    part_off = inh_offspr(m)
            
    ! combine parts; x=dam genotype, y=sire genotype   
    ! NOTE: dam prob sums over possible sire genotypes and vv  (marginal probs)
    forall (x=0:2)  calc_p_post(x,1) = logSumExp(part_par(:,2) + part_off(x,:))
    forall (y=0:2)  calc_p_post(y,2) = logSumExp(part_par(:,1) + part_off(:,y))
        
    ! scale to sum to 1
    forall (p=1:2)  calc_p_post(:,p) = logScale(calc_p_post(:,p))
    
!    call chk_logprob_OK(RESHAPE(calc_p_post, (/3*2/)), lbl='calc_p_post', lbl_i = m)
    
  end function calc_p_post
  
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! calc g_i(obs|act) * prod_k(ppost_ik) over mates k
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  function mates_post(i, q_not)
    double precision :: mates_post(0:2)
    integer, intent(IN) :: i  ! focal individual
    integer, intent(IN), optional :: q_not  ! mate to exclude from calculation
    integer :: k, nmates, psex
    double precision, allocatable :: post_per_mate(:,:)
    integer, pointer :: mates(:)
    
    mates => pop(i)%offspring_m 
    nmates = SIZE(mates)
    psex = pop(i)%sex  
    
    if (nmates < 2) then
      if (nmates==1 .and. present(q_not)) then
        if (pedigree(mates(1))%parent(3-psex) == q_not)  nmates = 0
      endif
      if (nmates==0) then
        mates_post = 0D0
      else if (nmates==1) then
        mates_post = lp_post(:, psex, mates(1))
      endif
      return
    endif
  
    allocate(post_per_mate(0:2, nmates))     
    post_per_mate = 0D0                          
    do k=1,nmates
      if (present(q_not)) then 
        if (pedigree(mates(k))%parent(3-psex) == q_not)  cycle
      endif
      post_per_mate(:,k) = lp_post(:, psex, mates(k))
    enddo    
    mates_post = SUM(post_per_mate, DIM=2)

  end function mates_post
  
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! calc prod_o(inh_o * mates_post)
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  function inh_offspr(m, q_not)
    double precision :: inh_offspr(0:2,0:2)
    integer, intent(IN) :: m
    integer, intent(IN), optional :: q_not  ! offspring to exclude from calculation
    integer :: o,x,y, noff
    double precision, allocatable :: per_off(:,:,:)   ! TODO? allocate once & reuse
    double precision :: off_tmp(0:2)
    integer, pointer :: offspring(:)
    
    offspring => pedigree(m)%offspring
    noff = SIZE(offspring)
    if (noff == 1 .and. present(q_not)) then
      if (offspring(1) == q_not)  noff = 0
    endif
    if (noff == 0) then
      inh_offspr = 0D0
     return
   endif
   
    allocate(per_off(0:2, 0:2, noff))
    per_off = 0D0
    do o=1,noff
      if (present(q_not)) then
        if (offspring(o) == q_not)  cycle
      endif
      off_tmp = lOcA(:,Gl(offspring(o))) + mates_post(offspring(o))
      ! sum over possible offspring genotypes; y=sire genotype, x=dam genotype
      forall (y=0:2, x=0:2)  per_off(x,y,o) = logSumExp(off_tmp + lAKA2P(:,x,y))
    enddo
    
    ! multiply across offspring
    inh_offspr = SUM(per_off, DIM=3)   
  end function inh_offspr
  
  
  ! !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! ! check if probabilities are valid (debugging tool)
  ! !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! subroutine chk_logprob_OK(logprob, lbl, lbl_i)
    ! double precision, intent(IN), dimension(:) :: logprob
    ! character(len=*), intent(IN) :: lbl
    ! integer, intent(IN) :: lbl_i
    
    ! if (any(logprob /= logprob) .or. any(logprob > sqrt(EPSILON(0d0))) .or. all(logprob < -HUGE(0d0))) then
      ! print *, ''
      ! print *, trim(lbl), ' not a valid probability! ', lbl_i
      ! write(*,'(10f9.4)') logprob
      ! print *, ''
      ! stop
    ! endif
  
  ! end subroutine chk_logprob_OK
  
end module impute_fun


!===============================================================================
!===============================================================================

program main
  use sqa_fileIO, ONLY: read_geno
  use sqa_general, ONLY: PreCalcProbs, printt
  use global_variables
  use pedigree_fun, ONLY: init_pop, read_pedigree, init_pedigree
  use check_pedigree
  use impute_fun
  implicit none
  
 ! integer :: i
  character(len=200) :: GenoFile, PedigreeFile, AFFile, GenoOutFile, EditsFile
  character(len=3) :: GenoInFormat, GenoOutFormat
  real :: Threshold_pedclean
  logical :: do_pedclean, do_pedigree
  double precision :: Er, ErV(3)
!  integer :: i,m
!  integer, allocatable :: mates(:)
    
  call read_args()   ! set defaults & read command line arguments
  
  if (.not. quiet)  call printt('reading genotypes ...')
  call read_geno(Geno=Geno, ID=IdV, SNP_names=SNP_names, FileName=GenoFile, &
    FileFormat = GenoInFormat, transp=.FALSE.)
  nIndG = SIZE(Geno, DIM=1) -1 ! dim 0:nIndG
  nSnp  = SIZE(Geno, DIM=2)
  
  if (do_pedclean .or. do_snpclean .or. do_impute) then
    if (.not. quiet)  call printt('initializing population ...')
    call init_pop(IdV)
  endif
 
  if (.not. quiet)  call print_sumstats('IN') 
 
  if (do_pedigree) then
    if (.not. quiet)  call printt('reading pedigree file ...')
    call read_pedigree(PedigreeFile)
    if (.not. quiet) print *, 'Total # individuals: ', nIndT
  endif
  call init_pedigree()  
  ! if (do_pedigree) then
    ! if (.not. quiet) call printt('pruning pedigree ...')
    
    ! call prune_pedigree()
    ! if (.not. quiet) print *, 'Total # individuals post-prune: ', nIndT
  ! endif
  
  if (do_pedclean .or. do_snpclean .or. do_impute) then  ! not if do_read_edits
    allocate(AF(nSnp))
    AF = getAF(AFFile) 
    
    if (erV(1) > TINY(0D0)) then
      call PreCalcProbs(nSnp, AF, ErV=ErV)
    else
      call PreCalcProbs(nSnp, AF, Er=Er)  
    endif
    
  endif
  
  if (do_pedclean) then
    if (.not. quiet)  call printt('cleaning pedigree ... ')   
    call clean_pedigree('pedclean_edits.txt', Threshold_pedclean)
  else
    if (.not. quiet)  call printt('skipping pedigree cleaning ...')
  endif
  
  if (do_pedclean .or. do_snpclean .or. do_impute) then 
    call clean_n_impute(EditsFile)  ! , do_pedigree
  endif

  if (do_geno_out) then
    if (.not. quiet)  call printt('writing new genotypes to file ...')
    call apply_edits(EditsFile, GenoOutFile, GenoOutFormat)
  endif
  
  if (.not. quiet)  call printt('done.')
  
  if (.not. quiet .and. do_geno_out)  call print_sumstats('OUT') 
  
  call deallocall()


contains
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !  read in command line arguments
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine read_args()
    use sqa_fileIO, only: valid_formats
    integer :: nArg, i, x, z
    character(len=32) :: arg, argOption
    logical :: only_read_edits    
    
    ! set defaults
    GenoFile = 'Geno'
    GenoInFormat = 'SEQ'
    do_pedigree = .TRUE.
    PedigreeFile = 'Pedigree.txt'  ! 'ped_griffin.txt'
    AFFile = 'NoFile'
    GenoOutFile = 'geno_imputed'
    GenoOutFormat = 'SEQ'
    EditsFile = 'imputation_edits.txt'
    Er = 0.001
    ErV = (/ 0.0d0, 0.0d0, 0.0d0 /)  ! hom|other hom, het|hom, and hom|het
    do_pedclean = .TRUE.
    Threshold_pedclean = 0.01
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
    only_read_edits   = .FALSE.
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

          case ('--version')
            print *, version
             
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
            
          case ('--no-pedigree')
            do_pedigree = .FALSE.
            
          case ('--err')
            i = i+1
            call get_command_argument(i, argOption)
            read(argOption, *)  Er   
            if (Er <= 0.0 .or. Er > 0.5)  stop 'please provide a genotyping error rate --err >0 and <=0.5'        

          case ('--errV')
            do z=1,3
              i = i+1
              call get_command_argument(i, argOption)
              read(argOption, *)  ErV(z)
              if (ErV(z) <= 0.0 .or. ErV(z) > 0.5)  stop 'Genotyping error rates must be >0 and <=0.5'
            enddo           
            
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
         
          case ('--edits', '--edits-out')
            i = i+1
            call get_command_argument(i, EditsFile)
            
          case ('--edits-in')
            i = i+1
            call get_command_argument(i, EditsFile)
            only_read_edits   = .TRUE.

          case ('--quiet')
            quiet = .TRUE.
            
          case default
            print '(2a, /)', 'Unrecognised command-line option: ', arg
            print *, 'see --help for valid options'
!            call print_help()
            stop

        end select
      end do
    endif

    if (do_quick) then
 !     do_pedclean = .FALSE.
      do_snpclean = .FALSE.
    endif
    
    if (only_read_edits) then
      do_pedigree = .FALSE.
      do_impute = .FALSE.
      do_geno_out = .TRUE.
    endif
    
    if (.not. do_pedigree) then
      do_pedclean = .FALSE.
      do_snpclean = .FALSE.   ! relies on pedigree to detect genotyping errors
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
    print '(a)',    '  --no-pedigree        impute without pedigree, use `--when-in-doubt` everywhere'  
    print '(a)',    '  --err <value>        presumed genotyping error rate',&
                    '                        (Will be estimated from data in future version)'  
     print '(a)',    '  --errV <3 values>   alternative to --err: P(observed|actual) for',&
                    '                        hom|other hom, het|hom, and hom|het'                
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
    print '(a)',    '  --no-geno-out        no genotype output file with edits (only imputation_edits.txt), which',& 
                    '                          can be read in to apply the edits'
    print '(a)',    '  --probs-out          write all genotype probabilities to a text file, with 3 rows',&
                    '                         per SNP and 1 column per ID'
    print '(a)',    '  --edits-out <filename>  file name for list with edits. Default: imputation_edits.txt'
    print '(a)',    '  --edits-in <filename>   edits file as input'
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
  if (allocated(indiv2mating))  deallocate(indiv2mating)
  if (allocated(IdV))  deallocate(IdV) 
  if (allocated(SNP_names))  deallocate(SNP_names) 
  if (allocated(AF)) deallocate(AF)  
   if (allocated(pedigree))  deallocate(pedigree)
  ! pedigree_fun
  if (allocated(pop))  deallocate(pop)
  if (allocated(pedigree))  deallocate(pedigree)
  ! sqa_general
  if (allocated(AHWE)) deallocate(AHWE)
  if (allocated(OHWE)) deallocate(OHWE)
  if (allocated(AKAP)) deallocate(AKAP)
  if (allocated(OKAP)) deallocate(OKAP)
  
 

end subroutine deallocall
  
!===============================================================================
! end. 