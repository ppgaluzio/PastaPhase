!
! =============================================================
!
!                        MODULE POLYEOS
!
! =============================================================
!
!
!
! --------------------------NOTES------------------------------
!
!  Calculate the polytropic equation of state
!
!
! --------------------Public subroutines-----------------------
!
!
! -------------------------------------------------------------
! Code by Jerome Margueron
! INT Seattle, USA and IPN Lyon, France
! Date: 22-01-2018
! Latest revision - 22-01-2018
! -------------------------------------------------------------

module polyEosT0

  use acc; use cst;
  use polyEosT0Type;

  implicit none
  
  public

contains



! ======================================================
! ======================================================
! Read inputs as arguments
!
! args1 is the value of the symmetry energy (integer)
! args2 can be: BSTph, BSTmicro, EFT1 (Ingo), EFT2 (Drischler)
!
! ======================================================
! ======================================================
!
subroutine polyeos_read_inputs(coef_poly)

!  use acc; use cst;

  implicit none
  !
  ! ------------------------------
  ! Output variables
  ! ------------------------------
  !
  TYPE (polyEos_coef), intent(out) :: coef_poly
  !
  ! ------------------------------
  ! Local variables
  ! ------------------------------
  !
  integer, dimension(polyEos_neparam)   :: eparam
  !
  integer :: num_args, ix, i
  character(len=10), dimension(:), allocatable :: args
  !
  ! ======================================================
  ! Start routine
  ! ======================================================
  !
  ! pass the values of the empirical parameters (integer numbers)
  !
  num_args = command_argument_count()
  if (num_args.eq.0) then
     if (polyEos_lverb) write(*,*)"  read args from input file"
     open(unit=10, file="polyEos.in", status="old")
     read(10,*)eparam(:)
     close(10)
  else if (num_args.eq.12) then
     if (polyEos_lverb) write(*,*)"  read args from prompt"
     allocate(args(num_args))
     do ix = 1, num_args
        call get_command_argument(ix,args(ix))
     end do
     if (polyEos_lverb) write(*,*)"  args():",args(:)
     do i = 1, polyEos_neparam
        read(args(i),'(i8)')eparam(i)
     enddo
     deallocate(args)
  else
     stop "incorrect number of input parameters"
  endif
  if (polyEos_lverb) write(*,*)"  Read arguments:"
  if (polyEos_lverb) write(*,'("   param():",12i8)')eparam(:)
  !
  ! Define the coefficients of the polytropic EOS
  !
  coef_poly%logp1     = 1.d-3*eparam(1)
  coef_poly%gam(1)    = 1.d-3*eparam(2)
  coef_poly%logrho(1) = 1.d-1*eparam(3)
  coef_poly%gam(2)    = 1.d-3*eparam(4)
  coef_poly%logrho(2) = 1.d-1*eparam(5)
  coef_poly%gam(3)    = 1.d-3*eparam(6)

  coef_poly%p(1) = 10.d0**coef_poly%logp1
  coef_poly%trho(1) = 10.d0**coef_poly%logrho(1)
  coef_poly%trho(2) = 10.d0**coef_poly%logrho(2)
  
  coef_poly%k(1) = coef_poly%p(1) / coef_poly%trho(1)**coef_poly%gam(1)
  coef_poly%k(2) = coef_poly%k(1) * coef_poly%trho(1)**( coef_poly%gam(1) - coef_poly%gam(2) )
  coef_poly%k(3) = coef_poly%k(2) * coef_poly%trho(2)**( coef_poly%gam(2) - coef_poly%gam(3) )

  coef_poly%p(2) = coef_poly%k(2) * coef_poly%trho(2)**coef_poly%gam(2)
  
  coef_poly%a(1) = 0.d0
  coef_poly%trho(1) = coef_poly%trho(1) + coef_poly%p(1) / ( coef_poly%gam(1) - 1.0 )
  coef_poly%a(2) = coef_poly%trho(1) / coef_poly%trho(1) - 1.0 - &
       &coef_poly%k(2) / ( coef_poly%gam(2) - 1.0 ) * coef_poly%trho(1)**(coef_poly%gam(2)-1.0) 
  coef_poly%trho(2) = ( 1.0 + coef_poly%a(2) ) * coef_poly%trho(2) + coef_poly%p(2) / ( coef_poly%gam(2) - 1.0 )
  coef_poly%a(3) = coef_poly%trho(2) / coef_poly%trho(2) - 1.0 - &
       &coef_poly%k(3) / ( coef_poly%gam(3) - 1.0 ) * coef_poly%trho(2)**(coef_poly%gam(3)-1.0) 
  !
  coef_poly%nsat = 0.16
  !
  write(*,'(a8,3d12.3)')"p(i)",coef_poly%p(:)
  write(*,'(a8,3d12.3)')"rho(i)",coef_poly%trho(:)
  write(*,'(a8,3d12.3)')"k(i)",coef_poly%k(:)
  write(*,'(a8,3d12.3)')"gam(i)",coef_poly%gam(:)
  write(*,'(a8,3d12.3)')"a(i)",coef_poly%a(:)
  write(*,'(a8,3d12.3)')"rho(i)",coef_poly%trho(:)
  !
  ! ======================================================
  ! End routine
  ! ======================================================
  !
end subroutine polyeos_read_inputs








subroutine polyeos_compute_T0(coef_poly)
  
  TYPE (polyEos_coef), intent(in) :: coef_poly
  !
  ! local variables
  !
  TYPE (polyEos_eos)     :: eos
  !
  ! define the max density where the physical filter is applied
  integer        :: nden
  real (kind=pr) :: den_step
  !
  real (kind=pr) :: den, rho
  integer :: iden
  real (kind=pr), dimension(:,:), allocatable :: aeos
  !
  ! parameters
  !
  den_step = CST_nsat / 5._pr
  nden = 20
  allocate(aeos(nden,5))
  !
  do iden = 1, nden
     den = iden * den_step
     rho = den * CST_mnuc2 * CONV_MeV_fm3_to_g_cm3
     !
     ! calculate beta equilibrium with n, p, e, muon
     !
     call polyeos_T0(rho,coef_poly,eos)
     !
     aeos(iden,1) = eos%p
     aeos(iden,2) = eos%rho
     aeos(iden,3) = eos%enthalpy
     aeos(iden,4) = eos%cs2
     aeos(iden,5) = eos%e2a

     write(*,'(i4,11f14.4)')iden,den,rho,aeos(iden,1),aeos(iden,2),aeos(iden,3),aeos(iden,4),aeos(iden,5)

  enddo

end subroutine polyeos_compute_T0


  




! -------------------------------------------------------------
! This routine calculate the meta equation of state for given
! values of the densities given in eosDen
! Results in o_eosb0
! -------------------------------------------------------------

  subroutine polyeos_T0(i_rho,i_coef,o_eos)

    implicit none

    ! input variables
    real (kind=pr),      intent(in)       :: i_rho
    type (polyEos_coef), intent(in)       :: i_coef
    ! output variables
    type (polyEos_eos), intent(out)   :: o_eos
    !
    ! local variables
    !
    integer :: idx

    if (i_rho.le.i_coef%trho(1)) then
       idx = 1
    else if (i_rho.le.i_coef%trho(2)) then
       idx = 2
    else
       idx = 3
    endif

    o_eos%p = i_coef%k(idx) * i_rho**i_coef%gam(idx)
    o_eos%rho = ( 1.0 + i_coef%a(idx) ) * i_rho + o_eos%p / ( i_coef%gam(idx) - 1.0 )
    o_eos%enthalpy = 1.0 + i_coef%a(idx) + i_coef%gam(idx) / i_coef%trho(idx) * o_eos%p / (i_coef%gam(idx) - 1.0 )
    o_eos%cs2 = i_coef%gam(idx) * o_eos%p / (i_rho * o_eos%enthalpy )
    o_eos%e2a = i_coef%a(idx) + o_eos%p / i_coef%trho(idx) / ( i_coef%gam(idx) - 1.0 )
    
  end subroutine polyeos_T0
     
     

end module polyEosT0
