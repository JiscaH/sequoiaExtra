! general parameters and subroutines 
!
! Jisca Huisman, jisca.huisman@gmail.com
!
! This code is available under GNU General Public License v3
!
!===============================================================================
module sqa_general
  implicit none
  
  integer, parameter :: ishort = selected_int_kind(1), ilong=selected_int_kind(8) 
  integer, parameter :: nchar_ID=40
  ! inheritance rules/probabilities; O=observed, A=actual
  double precision :: OcA(0:2,-1:2), OKA2P(-1:2,0:2,0:2), AKA2P(0:2,0:2,0:2)  
  double precision, allocatable :: AHWE(:,:), OHWE(:,:), AKAP(:,:,:), OKAP(:,:,:) 
  
  interface mk_OcA
    module procedure mk_OcA_s, mk_OcA_v 
  end interface
  
contains
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! LogSumExp
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! LSE(logx1, logx2) = log(x1 + x2)
  pure function logSumExp(logX)
    double precision, intent(IN) :: logX(0:2)
    double precision :: logSumExp, max_logX, d(0:2), s
  
    ! https://en.wikipedia.org/wiki/LogSumExp
    max_logX = MAXVAL(logX)
    if (max_logX < -HUGE(0D0)) then   ! all values are 0 on non-log scale
      logSumExp = max_logX
      return
    endif
    d = logX - max_logX
    s = SUM( EXP(d) )
    if (s > HUGE(0D0)) then
      logSumExp = max_logX
    else 
      logSumExp = max_logX + LOG(s)
    endif
  
  end function logSumExp
  
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! scale logx to sum to unity on non-logscale
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
  pure function logScale(logX)   
    double precision, intent(IN) :: logX(0:2)
    double precision :: logScale(0:2), logS
    
    logS = logSumExp(logX)
    if (logS < -HUGE(0D0)) then
      logScale = logX
    else
      logScale = LOG( EXP(logX) / EXP(logS) )
    endif
    ! use unscaled tiny value to allow convergence check
    WHERE(logScale < -HUGE(0D0))  logScale = logX

  end function logScale 
  
  
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! initialise arrays with Mendelian inheritance rules and probabilities
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine PreCalcProbs(nSnp, AF, Er, ErV)
    integer, intent(IN) :: nSnp
    double precision, intent(IN) :: AF(nSnp)
    double precision, intent(IN), optional :: Er
    double precision, intent(IN), optional :: ErV(3)   ! either one must be specified
    integer :: l,i,j

    ! Probability observed genotype (dim2) conditional on actual genotype (dim1)
    if (present(ErV)) then
      OcA = mk_OcA(ErV)
    else if (present(Er)) then
      OcA = mk_OcA(Er)
    else
      stop 'precalcProbs needs either Er or ErV'
    endif 

    ! probabilities actual genotypes under HWE
    allocate(AHWE(0:2,nSnp))
    do l=1,nSnp
      AHWE(0,l) = (1 - AF(l))**2 
      AHWE(1,l) = 2*AF(l)*(1-AF(l)) 
      AHWE(2,l) = AF(l)**2 
    enddo
    
    ! inheritance conditional on both parents  (actual genotypes)
    AKA2P(0,0,:) = dble((/ 1.0, 0.5, 0.0 /))
    AKA2P(0,1,:) = dble((/ 0.5, 0.25,0.0 /))
    AKA2P(0,2,:) = dble((/ 0.0, 0.0, 0.0 /))

    AKA2P(1,0,:) = dble((/ 0.0, 0.5, 1.0 /))
    AKA2P(1,1,:) = dble((/ 0.5, 0.5, 0.5 /))
    AKA2P(1,2,:) = dble((/ 1.0, 0.5, 0.0 /))

    AKA2P(2,0,:) = dble((/ 0.0, 0.0, 0.0 /))
    AKA2P(2,1,:) = dble((/ 0.0, 0.25,0.5 /))
    AKA2P(2,2,:) = dble((/ 0.0, 0.5, 1.0 /))
    
    ! probabilities observed genotypes under HWE  + genotyping error pattern
    allocate(OHWE(-1:2,nSnp))
    forall (l=1:nSnp, i=-1:2)  OHWE(i,l) = SUM( OcA(:,i) * AHWE(:, l) ) 
    
    ! inheritance conditional on 1 parent (Actual & Observed refer to genotype)
    allocate(AKAP(0:2,0:2,nSnp))  ! Actual Kid Actual Parent
    allocate(OKAP(-1:2,0:2,nSnp)) ! Observed Kid Actual Parent
    
    do l=1,nSnp
      AKAP(0, :, l) = (/ 1-AF(l), (1-AF(l))/2, 0.0D0 /)
      AKAP(1, :, l) = (/ AF(l), 0.5D0, 1-AF(l) /)
      AKAP(2, :, l) = (/ 0D0, AF(l)/2, AF(l) /)
    enddo   
   
   forall (l=1:nSnp, i=-1:2, j=0:2)  OKAP(i,j,l) = SUM(OcA(:,i) * AKAP(:,j,l)) 

  end subroutine PreCalcProbs
  
  
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! Observed Conditional on Actual matrix, from Er scalar or ErV vector
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  function mk_OcA_s(Er)  result(OcA)
    double precision :: OcA(0:2,-1:2)
    double precision, intent(IN) :: Er
    
    ! Probability observed genotype (dim2) conditional on actual genotype (dim1)
    ! (ErrFlavour' = 2.9)
    OcA(:,-1) = 1.0D0      ! missing 
    OcA(0, 0:2) = (/ 1-Er     , Er-(Er/2)**2, (Er/2)**2 /)  ! act=0
    OcA(1, 0:2) = (/ Er/2     , 1-Er        , Er/2 /)       ! act=1
    OcA(2, 0:2) = (/ (Er/2)**2, Er-(Er/2)**2, 1-Er /)       ! act=2   
  end function mk_OcA_s
  
  function mk_OcA_v(Er)  result(OcA)
    double precision :: OcA(0:2,-1:2)
    double precision, intent(IN) :: Er(3)  ! hom|other hom, het|hom, and hom|het
    
    ! Probability observed genotype (dim2) conditional on actual genotype (dim1)
    ! (ErrFlavour' = 2.9)
    OcA(:,-1) = 1.0D0      ! missing 
    OcA(0, 0:2) = (/ 1-Er(1)-Er(2), Er(2)    , Er(1) /)          ! act=0
    OcA(1, 0:2) = (/ Er(3)        , 1-2*Er(3), Er(3) /)          ! act=1
    OcA(2, 0:2) = (/ Er(1)        , Er(2)    , 1-Er(1)-Er(2) /)  ! act=2   
  end function mk_OcA_v 
  
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! print timestamp, without carriage return
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine timestamp(add_blank)
    logical, intent(IN), OPTIONAL :: add_blank
    logical :: do_blank
    integer :: date_time_values(8)
    
    if(present(add_blank))then
      do_blank=add_blank
    else
      do_blank=.FALSE.
    endif
    ! NOTE: print *, & write(*,*) have leading blank, write(*,'(...)') does not
    
    call date_and_time(VALUES=date_time_values)
    write(*,'(i2.2,":",i2.2,":",i2.2, 1X)', advance='no') date_time_values(5:7)
    if (do_blank) write(*, '(1X)', advance='no')
     
  end subroutine timestamp
    
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! print text, preceded by timestamp
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    
    subroutine printt(text, int)
      character(len=*), intent(IN), optional :: text
      integer, intent(IN), optional :: int   ! e.g. counter in a loop
      character(len=10) :: intchar
         
      call timestamp()
      if (present(text) .and. present(int)) then
        write(intchar,'(i10)') int
        print *, text//adjustl(trim(intchar))
      else if (present(text)) then
        print *, text
      else if (present(int)) then
        print *, int
      endif
    
    end subroutine printt

end module sqa_general


!===============================================================================
!===============================================================================

module Sort
  implicit none
  private :: Partition

 contains
 !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  recursive subroutine QsortC(A, Rank)
    double precision, intent(in out), dimension(:) :: A
    integer, intent(in out), dimension(:) :: Rank
    integer :: iq

    if(size(A) > 1) then
     call Partition(A, iq, Rank)
     call QsortC(A(:iq-1), Rank(:iq-1))
     call QsortC(A(iq:), Rank(iq:))
    endif
  end subroutine QsortC

  subroutine Partition(A, marker, Rank)
    double precision, intent(in out), dimension(:) :: A
    integer, intent(in out), dimension(:) :: Rank
    integer, intent(out) :: marker
    integer :: i, j, TmpI
    double precision :: temp
    double precision :: x      ! pivot point
    x = A(1)
    i= 0
    j= size(A) + 1
    do
     j = j-1
     do
      if (j < 1) exit
      if (A(j) <= x) exit
      j = j-1
     end do
     i = i+1
     do
      if (i >= size(A)) exit
      if (A(i) >= x) exit
      i = i+1
     end do
     if (i < j) then
      ! exchange A(i) and A(j)
      temp = A(i)
      A(i) = A(j)
      A(j) = temp 
      
      TmpI = Rank(i) 
      Rank(i) = Rank(j)
      Rank(j) = TmpI
     elseif (i == j) then
      marker = i+1
      return
     else
      marker = i
      return
     endif
    end do

  end subroutine Partition
  
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  function getRank(V)   ! get ranking based on integer values
    integer, intent(IN) :: V(:)
    integer :: getRank(SIZE(V)) 
    double precision :: V_dbl(SIZE(V))
    integer :: i, V_rank(SIZE(V))
 
    V_dbl = dble(V)
    V_rank = (/ (i, i=1, SIZE(V), 1) /)
    call QsortC(V_dbl, V_rank)
    getRank = V_rank
  end function getRank
  
end module sort

!===============================================================================