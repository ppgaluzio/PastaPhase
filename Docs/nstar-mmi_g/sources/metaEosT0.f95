!
! =============================================================
!
!                        MODULE EEOS !NEW!
!
! =============================================================
!
!
!
! --------------------------NOTES------------------------------
!
!  Calculate the meta equation of state
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

module metaEosT0

  use acc; use cst;
  use metaEosT0Type;
!  use fgas;

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
subroutine metaeos_read_inputs(coef)

!  use acc; use cst;

  implicit none
  !
  ! ------------------------------
  ! Output variables
  ! ------------------------------
  !
  TYPE (metaEos_coef), intent(out) :: coef
  !
  ! ------------------------------
  ! Local variables
  ! ------------------------------
  !
  integer, dimension(metaEos_neparam)   :: eparam
  TYPE (metaEos_empirical)      :: emp
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
     if (metaEos_lverb) write(*,*)"  read args from input file"
     open(unit=10, file="metaEos.in", status="old")
     read(10,*)eparam(:)
     close(10)
  else if (num_args.eq.metaEos_neparam) then
     if (metaEos_lverb) write(*,*)"  read args from prompt"
     allocate(args(num_args))
     do ix = 1, num_args
        call get_command_argument(ix,args(ix))
     end do
     if (metaEos_lverb) write(*,*)"  args():",args(:)
     do i = 1, metaEos_neparam
        read(args(i),'(i8)')eparam(i)
     enddo
     deallocate(args)
  else
     stop "incorrect number of input parameters"
  endif
  if (metaEos_lverb) write(*,*)"  Read arguments:"
  if (metaEos_lverb) write(*,'("   param():",15i8)')eparam(1:5)
  if (metaEos_lverb) write(*,'("   param():",15i8)')eparam(6:10)
  if (metaEos_lverb) write(*,'("   param():",15i8)')eparam(11:12)
  if (metaEos_lverb) write(*,'("   param():",15i8)')eparam(13:17)
  if (metaEos_lverb) write(*,'("   param():",15i8)')eparam(18:19)
  if (metaEos_lverb) write(*,'("   param():",15i8)')eparam(20:21)
  if (metaEos_lverb) write(*,'("   param():",15i8)')eparam(22:25)
  !
  ! Define the empirical parameters (real numbers)
  !
  emp%eparam(:) = eparam(:)
  !
  emp%param(0,0) =-0.01   * eparam(1)
  emp%nsat       = 0.0001 * eparam(2)
  emp%param(0,2) = 0.1    * eparam(3)
  emp%param(0,3) = 0.1    * eparam(4)
  emp%param(0,4) = 1.0    * eparam(5)
  emp%param(1,0) = 0.01   * eparam(6)
  emp%param(1,1) = 0.01   * eparam(7)
  emp%param(1,2) = 0.1    * eparam(8)
  emp%param(1,3) = 1.0    * eparam(9)
  emp%param(1,4) = 1.0    * eparam(10)
  emp%ms         = 0.0001 * eparam(11)
  emp%dm         = 0.0001 * eparam(12)
  ! non-quadratic terms:
  emp%param(2,0) = 0.01   * eparam(13)
  emp%param(2,1) = 0.01   * eparam(14)
  emp%param(2,2) = 0.1    * eparam(15)
  emp%param(2,3) = 1.0    * eparam(16)
  emp%param(2,4) = 1.0    * eparam(17)
  ! b coefficient:
  emp%bsat         = 0.01   * eparam(18) 
  emp%bsym         = 0.01   * eparam(19) 
  ! quark matter parameters
  emp%LQCD         = eparam(20) 
  emp%kQCD         = 0.1 * eparam(21) 
  !
  if (metaEos_lverb) write(*,*)"  Read arguments:"
  if (metaEos_lverb) write(*,'("   param():",15f12.2)')emp%param(0,0:4)
  if (metaEos_lverb) write(*,'("   param():",15f12.2)')emp%param(1,0:4)
  if (metaEos_lverb) write(*,'("   param():",15f12.2)')emp%ms,emp%dm
  if (metaEos_lverb) write(*,'("   param():",15f12.2)')emp%bsat,emp%bsym
  if (metaEos_lverb) write(*,'("   param():",15f12.2)')emp%LQCD,emp%kQCD
  !
  ! With or without kinetic energy term
  !
  emp%ikin = eparam(22)
  !
  ! EOS with or without muons
  !
  emp%imuon = eparam(23)
  !
  ! Non-relativistic (0) or relativistic (1) kinetic energy for the baryons
  !
  emp%irel = eparam(24)
  !
  ! Quak matter model
  !
  emp%iquark = eparam(25)
  ! Fixes the constant for the low density correction term
  !
  !emp%xb = 10.0*log(2.0)
  !
  ! Set up the coefficients of the metamodeling
  !
  call metaeos_setup_coeff(emp,coef)
  !
  ! ======================================================
  ! End routine
  ! ======================================================
  !
end subroutine metaeos_read_inputs





! -------------------------------------------------------------
! This routine calculate the meta-model coefficients from the
! empirical parameter defined in metaeos_emp -> i_emp
! -------------------------------------------------------------

  subroutine metaeos_setup_coeff(i_emp,o_coef)

!    use acc; use CST;
!    use metaEosType;

    implicit none

    ! input variables
    type (metaEos_empirical), intent(in) :: i_emp
    ! output variables
    type (metaEos_coef) :: o_coef
    ! local variables
    integer :: n
    real (kind=pr) :: c0tau, c1tau, c1, c2
    real (kind=pr) :: ccs, ccv, coefc
    real (kind=pr) :: gSM, gNM, dNM, dSM
    real (kind=pr) :: mspSM, msnSM, msNM, kFSM, kFNM, zNM, znSM, zpSM ! some constants in neutrons (NM) and symetric matter (SM)
    real (kind=pr) :: d, m, ms, z, kF, kappa, f, epsk ! useful variables
    real (kind=pr) :: d1n_ms, d2n_ms, d3n_ms, d4n_ms ! derivatives of ms with respect to density at constant asymetry 
    real (kind=pr) :: d1n_z, d2n_z, d3n_z, d4n_z ! derivatives of z with respect to density at constant asymetry 
    real (kind=pr) :: d1z_f, d2z_f, d3z_f, d4z_f ! derivatives of f with respect to z at constant effective mass
    real (kind=pr) :: d1n_f, d2n_f, d3n_f, d4n_f !derivatives of f with respect to density at constant asymetry
    real (kind=pr) :: d1n_epsk, d2n_epsk, d3n_epsk, d4n_epsk !derivatives of kinetic energy volumic density with respect to density at constant asymetry
    !
    ! Define the cut-off in order derivative (usually = 4)
    !
    o_coef%nmax = 4 ! cut-off in the density expansion
    !
    ! factorial
    !
    o_coef%facto(0) = 1
    do n = 1, 5
       o_coef%facto(n) = n * o_coef%facto(n-1)
    enddo
    !
    !
    !
    o_coef%bsat   = i_emp%bsat
    o_coef%bsym   = i_emp%bsym
    o_coef%nsat = i_emp%nsat
    o_coef%ms   = i_emp%ms
    o_coef%dm   = i_emp%dm
    !
    !
    !
    o_coef%kappas = 1.0 / i_emp%ms - 1.0
    o_coef%kappasym = 0._pr
    if (i_emp%dm.ne.0._pr) o_coef%kappasym = ( 1.0 - sqrt( 1.0 + (i_emp%dm**2) * ( 1.0 + o_coef%kappas )**2 ) ) / i_emp%dm
    o_coef%kappav = o_coef%kappas - o_coef%kappasym
    o_coef%kappaNM = o_coef%kappas + o_coef%kappasym
    !
    o_coef%mb   = o_coef%kappas
    o_coef%db   = o_coef%kappasym

!    write(*,*)"kappa s",kappas
!    write(*,*)"kappa v",kappav

    c0tau = CST_htm / 2. * o_coef%kappas / i_emp%nsat
    c1tau = CST_htm / 2. * ( o_coef%kappas - o_coef%kappav ) / i_emp%nsat
!    write(*,*)"mean c0tau ",c0tau
!    write(*,*)"mean c1tau ",c1tau
    c1 = o_coef%kappas / i_emp%nsat
    c2 = ( o_coef%kappas - o_coef%kappav ) / i_emp%nsat

    !
    ! non-relativistic case
    !
    if (i_emp%irel.eq.0) then
       !
       if (i_emp%ikin.eq.0) then
          if (metaEos_lverb) write(*,*)"Non-relativistic Meta-EoS without kinetic energy term"
          o_coef%TFGsat = 0._pr
          o_coef%TFGSsat = 0._pr
       elseif (i_emp%ikin.eq.1) then
          if (metaEos_lverb) write(*,*)"Non-relativistic Meta-EoS with kinetic energy term"
          o_coef%TFGsat = 3. / 10. * CST_htm * ( 3.0 * CST_pi2 / 2.0 * i_emp%nsat )**(2.0/3.0)
          o_coef%TFGSsat = o_coef%TFGsat * ( 1.0 + o_coef%mb )
       else
          stop "i_emp%ikin badly defined (different from 0 or 1) in metaeos_setup_coeff()"
       endif
       !
       ccs = o_coef%kappas
       o_coef%VCOEF(0,0) = i_emp%param(0,0) -     o_coef%TFGsat * ( 1 + ccs)
       o_coef%VCOEF(0,1) =                  -     o_coef%TFGsat * ( 2 + 5*ccs)
       o_coef%VCOEF(0,2) = i_emp%param(0,2) - 2 * o_coef%TFGsat * (-1 + 5*ccs)
       o_coef%VCOEF(0,3) = i_emp%param(0,3) - 2 * o_coef%TFGsat * ( 4 - 5*ccs)
       o_coef%VCOEF(0,4) = i_emp%param(0,4) - 8 * o_coef%TFGsat * (-7 + 5*ccs)
       !
       ccv = o_coef%kappaNM
       o_coef%VCOEF(1,0) = i_emp%param(1,0) -     o_coef%TFGsat * ( 2**CST_p23*( 1+  ccv) - ( 1+  ccs) ) - i_emp%param(2,0)
       o_coef%VCOEF(1,1) = i_emp%param(1,1) -     o_coef%TFGsat * ( 2**CST_p23*( 2+5*ccv) - ( 2+5*ccs) ) - i_emp%param(2,1)
       o_coef%VCOEF(1,2) = i_emp%param(1,2) - 2 * o_coef%TFGsat * ( 2**CST_p23*(-1+5*ccv) - (-1+5*ccs) ) - i_emp%param(2,2)
       o_coef%VCOEF(1,3) = i_emp%param(1,3) - 2 * o_coef%TFGsat * ( 2**CST_p23*( 4-5*ccv) - ( 4-5*ccs) ) - i_emp%param(2,3)
       o_coef%VCOEF(1,4) = i_emp%param(1,4) - 8 * o_coef%TFGsat * ( 2**CST_p23*(-7+5*ccv) - (-7+5*ccs) ) - i_emp%param(2,4)
       !
       !    write(*,*)"Check:", o_coef%kappas,o_coef%kappasym
       !    write(*,*)"Check:",ccs,ccv,i_emp%param(1,0),i_emp%param(2,0)
       !    write(*,*)'check:',o_coef%kappas,o_coef%kappaNM,o_coef%TFGsat,2**CST_p23
       !    write(*,*)'check:',o_coef%TFGsat * ( 1.0_pr +  ccv) !( 2**CST_p23*( 1+  ccs) )! - ( 1+  ccv) )
       !
       o_coef%VCOEF(2,0) = i_emp%param(2,0)
       o_coef%VCOEF(2,1) = i_emp%param(2,1)
       o_coef%VCOEF(2,2) = i_emp%param(2,2)
       o_coef%VCOEF(2,3) = i_emp%param(2,3)
       o_coef%VCOEF(2,4) = i_emp%param(2,4)
       !
       ! Relativistic case
       !
    elseif (i_emp%irel.eq.1) then
       !
       if (i_emp%ikin.eq.0) then
          if (metaEos_lverb) write(*,*)"Relativistic Meta-EoS without kinetic energy term"
          coefC = 0._pr
          o_coef%TFGsat = 0._pr
          o_coef%TFGSsat = 0._pr
       elseif (i_emp%ikin.eq.1) then
          if (metaEos_lverb) write(*,*)"Relativistic Meta-EoS with kinetic energy term"
          coefC = 1.0_pr / (8.0_pr * CST_pi2 * CST_hbc**3)
          o_coef%TFGsat = 3. / 10. * CST_htm * ( 3.0 * CST_pi2 / 2.0 * i_emp%nsat )**(2.0/3.0)
          o_coef%TFGSsat = o_coef%TFGsat * ( 1.0 + o_coef%mb )
       else
          stop "i_emp%ikin badly defined (different from 0 or 1) in metaeos_setup_coeff()"
       endif
!
       !definitions of msSM, msNM, kFSM, kFNM,zNM,zSM 
!       dSM = 4._pr
!       dNM = 2._pr
       gSM = ( 1.0_pr + o_coef%mb)
       gNM = ( 1.0_pr + o_coef%mb + o_coef%db )
       msNM = CST_mnc2/gNM
       msnSM = CST_mnc2/gSM
       mspSM = CST_mpc2/gSM
       kFSM = (CST_pi2*3._pr*CST_nsat/2.0_pr)**CST_p13
       kFNM = (CST_pi2*3._pr*CST_nsat)**CST_p13
       znSM = CST_hbc*kFSM/msnSM
       zpSM = CST_hbc*kFSM/mspSM
       zNM = CST_hbc*kFNM/msNM
!       write(*,*)"essai",zNM,kFNM,msNM,CST_mnc2,gNM,o_coef%mb,o_coef%db
!
       !allocation of some useful variables for symmetric matter
       ! for neutrons
!       d = dSM
       ms = msnSM
       m = CST_mnc2
!       ms = m
       z = znSM
       kF = kFSM
       kappa = o_coef%mb
       f = z*(1._pr+2._pr*z**2)*dsqrt(1._pr+z**2) - dlog(z+dsqrt(1._pr+z**2))
       epsk = coefC * (ms)**4 * f - CST_nsat/2._pr * m
!       if (i_emp%ikin.eq.0) epsk = 0._pr
!
       !derivatives of effective mass with respect to density at constant asymetry
       d1n_ms = -m*(kappa/CST_nsat)*(1.0_pr/(1.0_pr + kappa))**2
       d2n_ms = m*2.0_pr*(kappa/CST_nsat)**2*(1.0_pr/(1.0_pr + kappa))**3
       d3n_ms = -m*3.0_pr*2.0_pr*(kappa/CST_nsat)**3*(1.0_pr/(1.0_pr + kappa))**4
       d4n_ms = m*4.0_pr*3.0_pr*2.0_pr*(kappa/CST_nsat)**4*(1.0_pr/(1.0_pr + kappa))**5
!
       !derivatives of z with respect to density at constant asymetry
       d1n_z = z/(3.0_pr*CST_nsat) + (CST_hbc/m)*(kappa/CST_nsat)*kF
       d2n_z = -z/(3.0_pr*CST_nsat**2) + d1n_z/(3.0_pr*CST_nsat) + (CST_hbc/(3.0_pr*m))*(kappa/CST_nsat**2)*kF
       d3n_z = 2.0_pr*z/(3.0_pr*CST_nsat**3) - 2.0_pr*d1n_z/(3.0_pr*CST_nsat**2) + d2n_z/(3.0_pr*CST_nsat) &
            &- (2.0_pr*CST_hbc/(9.0_pr*m))*(kappa/CST_nsat**3)*kF
       d4n_z = -2.0_pr*z/(CST_nsat**4) + 2.0_pr*d1n_z/(CST_nsat**3) - d2n_z/(CST_nsat**2) + d3n_z/(3.0_pr*CST_nsat) &
            &+ (10.0_pr*CST_hbc/(27.0_pr*m))*(kappa/CST_nsat**4)*kF
!
       !derivatives of f with respect to z at constant effective mass
       d1z_f = 8._pr*z**2*(dsqrt(1._pr+z**2))
       d2z_f = 8._pr*z*((2._pr+3._pr*z**2)/dsqrt(1._pr+z**2))
       d3z_f = 8._pr*((2._pr+9._pr*z**2+6._pr*z**4)/(1._pr+z**2)**(3._pr/2._pr))
       d4z_f = 24._pr*z*(2._pr*z**4+5._pr*z**2+4._pr)/(1._pr+z**2)**(5._pr/2._pr)
!
       !derivatives of f with respect to density at constant asymetry
       d1n_f = d1z_f*d1n_z 
       d2n_f = d2z_f*d1n_z**2 + d1z_f*d2n_z 
       d3n_f = d3z_f*d1n_z**3 + 3._pr*d2z_f*d1n_z*d2n_z + d1z_f*d3n_z
       d4n_f = d4z_f*d1n_z**4 + 6._pr*d3z_f*d1n_z**2*d2n_z + 3._pr*d2z_f*d2n_z**2 + 4._pr*d2z_f*d1n_z*d3n_z + d1z_f*d4n_z 
!
       !derivatives of kinetic energy volumique density with respect to density at constant asymetry
       d1n_epsk = 4._pr*coefC*ms**3*d1n_ms*f + coefC*ms**4*d1n_f - m / 2._pr
       d2n_epsk = 12._pr*coefC*ms**2*d1n_ms**2*f + 4._pr*coefC*ms**3*(d2n_ms*f + 2._pr*d1n_ms*d1n_f) + coefC*ms**4*d2n_f
       d3n_epsk = 24._pr*coefC*ms*d1n_ms**3*f + 36._pr*coefC*ms**2*d1n_ms*(d2n_ms*f + d1n_ms*d1n_f) &
            &+ 4._pr*coefC*ms**3*(d3n_ms*f + 3._pr*d2n_ms*d1n_f + 3._pr*d1n_ms*d2n_f) + coefC*ms**4*d3n_f
       d4n_epsk = 24._pr*coefC*d1n_ms**4*f + 24._pr*coefC*ms*(6._pr*d1n_ms**2*d2n_ms*f + 4._pr*d1n_ms**3*d1n_f) &
            &+ 12._pr*coefC*ms**2*(3._pr*d2n_ms**2*f + 4._pr*d1n_ms*d3n_ms*f + 12._pr*d1n_ms*d2n_ms*d1n_f + 6._pr*d1n_ms**2*d2n_f) &
            &+ 4._pr*coefC*ms**3*(d4n_ms*f + 4._pr*d3n_ms*d1n_f + 6._pr*d2n_ms*d2n_f + 4._pr*d1n_ms*d3n_f) + coefC*ms**4*d4n_f
       !

       !isoscalar coefficients
       o_coef%VCOEF(0,0) = i_emp%param(0,0) - epsk/CST_nsat
       o_coef%VCOEF(0,1) =                  + 3._pr*(epsk/CST_nsat - d1n_epsk)
       o_coef%VCOEF(0,2) = i_emp%param(0,2) - 9._pr*(2._pr*epsk/CST_nsat - 2._pr*d1n_epsk + CST_nsat*d2n_epsk)
      
       o_coef%VCOEF(0,3) = i_emp%param(0,3) + 27._pr*(6._pr*epsk/CST_nsat - 6._pr*d1n_epsk + 3._pr*CST_nsat*d2n_epsk &
            &- CST_nsat**2*d3n_epsk)
       o_coef%VCOEF(0,4) = i_emp%param(0,4) - 81._pr*(24._pr*epsk/CST_nsat - 24._pr*d1n_epsk + 12._pr*CST_nsat*d2n_epsk &
            &- 4._pr*CST_nsat**2*d3n_epsk + CST_nsat**3*d4n_epsk)
!
!
!
!
!
!
       ! for protons
!       d = dSM
       ms = mspSM
       m = CST_mpc2
!       ms = m
       z = zpSM
       kF = kFSM
       kappa = o_coef%mb
       f = z*(1._pr+2._pr*z**2)*dsqrt(1._pr+z**2) - dlog(z+dsqrt(1._pr+z**2))
       epsk = coefC * (ms)**4 * f - CST_nsat/2._pr * m
!       if (i_emp%ikin.eq.0) epsk = 0._pr
!
       !derivatives of effective mass with respect to density at constant asymetry
	   d1n_ms = -m*(kappa/CST_nsat)*(1.0_pr/(1.0_pr + kappa))**2
       d2n_ms = m*2.0_pr*(kappa/CST_nsat)**2*(1.0_pr/(1.0_pr + kappa))**3
       d3n_ms = -m*3.0_pr*2.0_pr*(kappa/CST_nsat)**3*(1.0_pr/(1.0_pr + kappa))**4
       d4n_ms = m*4.0_pr*3.0_pr*2.0_pr*(kappa/CST_nsat)**4*(1.0_pr/(1.0_pr + kappa))**5
!
       !derivatives of z with respect to density at constant asymetry
	   d1n_z = z/(3.0_pr*CST_nsat) + (CST_hbc/m)*(kappa/CST_nsat)*kF
       d2n_z = -z/(3.0_pr*CST_nsat**2) + d1n_z/(3.0_pr*CST_nsat) + (CST_hbc/(3.0_pr*m))*(kappa/CST_nsat**2)*kF
       d3n_z = 2.0_pr*z/(3.0_pr*CST_nsat**3) - 2.0_pr*d1n_z/(3.0_pr*CST_nsat**2) + d2n_z/(3.0_pr*CST_nsat) &
            &- (2.0_pr*CST_hbc/(9.0_pr*m))*(kappa/CST_nsat**3)*kF
       d4n_z = -2.0_pr*z/(CST_nsat**4) + 2.0_pr*d1n_z/(CST_nsat**3) - d2n_z/(CST_nsat**2) + d3n_z/(3.0_pr*CST_nsat) &
            &+ (10.0_pr*CST_hbc/(27.0_pr*m))*(kappa/CST_nsat**4)*kF
!
       !derivatives of f with respect to z at constant effective mass
	   d1z_f = 8._pr*z**2*(dsqrt(1._pr+z**2))
       d2z_f = 8._pr*z*((2._pr+3._pr*z**2)/dsqrt(1._pr+z**2))
       d3z_f = 8._pr*((2._pr+9._pr*z**2+6._pr*z**4)/(1._pr+z**2)**(3._pr/2._pr))
       d4z_f = 24._pr*z*(2._pr*z**4+5._pr*z**2+4._pr)/(1._pr+z**2)**(5._pr/2._pr)
!
       !derivatives of f with respect to density at constant asymetry
	   d1n_f = d1z_f*d1n_z 
       d2n_f = d2z_f*d1n_z**2 + d1z_f*d2n_z 
       d3n_f = d3z_f*d1n_z**3 + 3._pr*d2z_f*d1n_z*d2n_z + d1z_f*d3n_z
       d4n_f = d4z_f*d1n_z**4 + 6._pr*d3z_f*d1n_z**2*d2n_z + 3._pr*d2z_f*d2n_z**2 + 4._pr*d2z_f*d1n_z*d3n_z + d1z_f*d4n_z 
!
       !derivatives of kinetic energy volumique density with respect to density at constant asymetry
       d1n_epsk = 4._pr*coefC*ms**3*d1n_ms*f + coefC*ms**4*d1n_f - m / 2._pr
       d2n_epsk = 12._pr*coefC*ms**2*d1n_ms**2*f + 4._pr*coefC*ms**3*(d2n_ms*f + 2._pr*d1n_ms*d1n_f) + coefC*ms**4*d2n_f
       d3n_epsk = 24._pr*coefC*ms*d1n_ms**3*f + 36._pr*coefC*ms**2*d1n_ms*(d2n_ms*f + d1n_ms*d1n_f) &
            &+ 4._pr*coefC*ms**3*(d3n_ms*f + 3._pr*d2n_ms*d1n_f + 3._pr*d1n_ms*d2n_f) + coefC*ms**4*d3n_f
       d4n_epsk = 24._pr*coefC*d1n_ms**4*f + 24._pr*coefC*ms*(6._pr*d1n_ms**2*d2n_ms*f + 4._pr*d1n_ms**3*d1n_f) &
            &+ 12._pr*coefC*ms**2*(3._pr*d2n_ms**2*f + 4._pr*d1n_ms*d3n_ms*f + 12._pr*d1n_ms*d2n_ms*d1n_f + 6._pr*d1n_ms**2*d2n_f) &
            &+ 4._pr*coefC*ms**3*(d4n_ms*f + 4._pr*d3n_ms*d1n_f + 6._pr*d2n_ms*d2n_f + 4._pr*d1n_ms*d3n_f) + coefC*ms**4*d4n_f
       !
       !write(*,*)"TEST (nsat):",d1n_epsk-epsk/CST_nsat
       !isoscalar coefficients
       o_coef%VCOEF(0,0) = o_coef%VCOEF(0,0) - epsk/CST_nsat
       o_coef%VCOEF(0,1) = o_coef%VCOEF(0,1) + 3._pr*(epsk/CST_nsat - d1n_epsk)
     if (metaEos_lverb) write(*,*)"param K prot (IS):",9._pr*(2._pr*epsk/CST_nsat - 2._pr*d1n_epsk + CST_nsat*d2n_epsk),&
            & 2._pr*epsk/CST_nsat,- 2._pr*d1n_epsk,CST_nsat*d2n_epsk
       o_coef%VCOEF(0,2) = o_coef%VCOEF(0,2) - 9._pr*(2._pr*epsk/CST_nsat - 2._pr*d1n_epsk + CST_nsat*d2n_epsk)
       o_coef%VCOEF(0,3) = o_coef%VCOEF(0,3) + 27._pr*(6._pr*epsk/CST_nsat - 6._pr*d1n_epsk + 3._pr*CST_nsat*d2n_epsk &
            &- CST_nsat**2*d3n_epsk)
       o_coef%VCOEF(0,4) = o_coef%VCOEF(0,4) - 81._pr*(24._pr*epsk/CST_nsat - 24._pr*d1n_epsk + 12._pr*CST_nsat*d2n_epsk &
            &- 4._pr*CST_nsat**2*d3n_epsk + CST_nsat**3*d4n_epsk)
!
     if (metaEos_lverb) write(*,*)"VCOEF (IS)"
     if (metaEos_lverb) write(*,*)o_coef%VCOEF(0,0),o_coef%VCOEF(0,1),o_coef%VCOEF(0,2),o_coef%VCOEF(0,3),o_coef%VCOEF(0,4),epsk
!       write(*,*)
!
!
       !reallocation of variables for neutron matter
       ! 
       ! for neutrons
!       d = dNM
       ms = msNM
       m = CST_mnc2
!       ms = m
       z = zNM
       kF = kFNM
       kappa = o_coef%mb + o_coef%db 
       f = z*(1._pr+2._pr*z**2)*dsqrt(1._pr+z**2) - dlog(z+dsqrt(1._pr+z**2))
       epsk = coefC * (ms)**4 * f - CST_nsat * m
       !write(*,*)"here",z,kF,f,epsk
!       if (i_emp%ikin.eq.0) epsk = 0._pr
!
       !derivatives of effective mass with respect to density at constant asymetry
       d1n_ms = -m*(kappa/CST_nsat)*(1.0_pr/(1.0_pr + kappa))**2
       d2n_ms = m*2.0_pr*(kappa/CST_nsat)**2*(1.0_pr/(1.0_pr + kappa))**3
       d3n_ms = -m*3.0_pr*2.0_pr*(kappa/CST_nsat)**3*(1.0_pr/(1.0_pr + kappa))**4
       d4n_ms = m*4.0_pr*3.0_pr*2.0_pr*(kappa/CST_nsat)**4*(1.0_pr/(1.0_pr + kappa))**5
!
       !derivatives of z with respect to density at constant asymetry
       d1n_z = z/(3.0_pr*CST_nsat) + (CST_hbc/m)*(kappa/CST_nsat)*kF
       d2n_z = -z/(3.0_pr*CST_nsat**2) + d1n_z/(3.0_pr*CST_nsat) + (CST_hbc/(3.0_pr*m))*(kappa/CST_nsat**2)*kF
       d3n_z = 2.0_pr*z/(3.0_pr*CST_nsat**3) - 2.0_pr*d1n_z/(3.0_pr*CST_nsat**2) + d2n_z/(3.0_pr*CST_nsat) &
            &- (2.0_pr*CST_hbc/(9.0_pr*m))*(kappa/CST_nsat**3)*kF
       d4n_z = -2.0_pr*z/(CST_nsat**4) + 2.0_pr*d1n_z/(CST_nsat**3) - d2n_z/(CST_nsat**2) + d3n_z/(3.0_pr*CST_nsat) &
            &+ (10.0_pr*CST_hbc/(27.0_pr*m))*(kappa/CST_nsat**4)*kF
!
       !derivatives of f with respect to z at constant effective mass
       d1z_f = 8._pr*z**2*(dsqrt(1._pr+z**2))
       d2z_f = 8._pr*z*((2._pr+3._pr*z**2)/dsqrt(1._pr+z**2))
       d3z_f = 8._pr*((2._pr+9._pr*z**2+6._pr*z**4)/(1._pr+z**2)**(3._pr/2._pr))
       d4z_f = 24._pr*z*(2._pr*z**4+5._pr*z**2+4._pr)/(1._pr+z**2)**(5._pr/2._pr)
!
       !derivatives of f with respect to density at constant asymetry
       d1n_f = d1z_f*d1n_z 
       d2n_f = d2z_f*d1n_z**2 + d1z_f*d2n_z 
       d3n_f = d3z_f*d1n_z**3 + 3._pr*d2z_f*d1n_z*d2n_z + d1z_f*d3n_z
       d4n_f = d4z_f*d1n_z**4 + 6._pr*d3z_f*d1n_z**2*d2n_z + 3._pr*d2z_f*d2n_z**2 + 4._pr*d2z_f*d1n_z*d3n_z + d1z_f*d4n_z 
!
       !derivatives of kinetic energy volumique density with respect to density at constant asymetry
       d1n_epsk = 4._pr*coefC*ms**3*d1n_ms*f + coefC*ms**4*d1n_f - m
       d2n_epsk = 12._pr*coefC*ms**2*d1n_ms**2*f + 4._pr*coefC*ms**3*(d2n_ms*f + 2._pr*d1n_ms*d1n_f) + coefC*ms**4*d2n_f
       d3n_epsk = 24._pr*coefC*ms*d1n_ms**3*f + 36._pr*coefC*ms**2*d1n_ms*(d2n_ms*f + d1n_ms*d1n_f) &
            &+ 4._pr*coefC*ms**3*(d3n_ms*f + 3._pr*d2n_ms*d1n_f + 3._pr*d1n_ms*d2n_f) + coefC*ms**4*d3n_f
       d4n_epsk = 24._pr*coefC*d1n_ms**4*f + 24._pr*coefC*ms*(6._pr*d1n_ms**2*d2n_ms*f + 4._pr*d1n_ms**3*d1n_f) &
            &+ 12._pr*coefC*ms**2*(3._pr*d2n_ms**2*f + 4._pr*d1n_ms*d3n_ms*f + 12._pr*d1n_ms*d2n_ms*d1n_f + 6._pr*d1n_ms**2*d2n_f) &
            &+ 4._pr*coefC*ms**3*(d4n_ms*f + 4._pr*d3n_ms*d1n_f + 6._pr*d2n_ms*d2n_f + 4._pr*d1n_ms*d3n_f) + coefC*ms**4*d4n_f
       !
       !isovector coefficients
       o_coef%VCOEF(1,0) = i_emp%param(0,0) + i_emp%param(1,0) - epsk/CST_nsat - o_coef%VCOEF(0,0) - i_emp%param(2,0)
       o_coef%VCOEF(1,1) = i_emp%param(1,1) + 3._pr*(epsk/CST_nsat - d1n_epsk) - o_coef%VCOEF(0,1) - i_emp%param(2,1)
       o_coef%VCOEF(1,2) = i_emp%param(0,2) + i_emp%param(1,2) - 9._pr*(2._pr*epsk/CST_nsat - 2._pr*d1n_epsk &
            &+ CST_nsat*d2n_epsk) - o_coef%VCOEF(0,2) - i_emp%param(2,2)
       o_coef%VCOEF(1,3) = i_emp%param(0,3) + i_emp%param(1,3) + 27._pr*(6._pr*epsk/CST_nsat - 6._pr*d1n_epsk & 
            &+ 3._pr*CST_nsat*d2n_epsk - CST_nsat**2*d3n_epsk) - o_coef%VCOEF(0,3) - i_emp%param(2,3)
       o_coef%VCOEF(1,4) = i_emp%param(0,4) + i_emp%param(1,4) - 81._pr*(24._pr*epsk/CST_nsat - 24._pr*d1n_epsk &
            &+ 12._pr*CST_nsat*d2n_epsk - 4._pr*CST_nsat**2*d3n_epsk + CST_nsat**3*d4n_epsk) - o_coef%VCOEF(0,4) - i_emp%param(2,4)

     if (metaEos_lverb) write(*,*)"VCOEF (IV)"
     if (metaEos_lverb) write(*,*)o_coef%VCOEF(1,0),o_coef%VCOEF(1,1),o_coef%VCOEF(1,2),o_coef%VCOEF(1,3),o_coef%VCOEF(1,4)
     if (metaEos_lverb) write(*,*)" test epsk"
     if (metaEos_lverb) write(*,*)epsk,f,d1n_f,d2n_f,d3n_f,d4n_f,z
!
	   ! non-quadratic contribution
       o_coef%VCOEF(2,0) = i_emp%param(2,0)
       o_coef%VCOEF(2,1) = i_emp%param(2,1)
       o_coef%VCOEF(2,2) = i_emp%param(2,2)
       o_coef%VCOEF(2,3) = i_emp%param(2,3)
       o_coef%VCOEF(2,4) = i_emp%param(2,4)

    endif

    ! quark model coefficients
    o_coef%LQCD = i_emp%LQCD
    o_coef%kQCD = i_emp%kQCD
    !
    o_coef%eparam(:) = i_emp%eparam(:)
    o_coef%param(:,:) = i_emp%param(:,:)
    !
    o_coef%ikin = i_emp%ikin
    o_coef%imuon = i_emp%imuon
    o_coef%irel = i_emp%irel
    o_coef%iquark = i_emp%iquark

!    o_coef%VCOEF(:,:)=0.0_pr
   

    if (i_emp%iquark.eq.0) then
        if (metaEos_lverb) write(*,*)"no quark model"
    elseif (i_emp%iquark.eq.1) then
        if (metaEos_lverb) write(*,*)"quarkyonic model"
    endif
    
  end subroutine metaeos_setup_coeff





  

subroutine metaeos_compute_T0_MCMC(coef,withq)

! use acc; use cst;
!  use metaeosT0; use metaEosType;
  
  TYPE (metaEos_coef), intent(in) :: coef
  logical,             intent(in) :: withq 
  !
  ! local variables
  !
  TYPE (metaEos_densities)     :: eosDen
  TYPE (metaeos_Baryons)       :: eosb0
  TYPE (metaEos_Q1)            :: eosq0
  !
  ! define the max density where the physical filter is applied
  integer        :: nden
  real (kind=pr) :: den_step
  !
  real (kind=pr) :: den, xd, xmuon
  integer :: iden
  !
  ! parameters
  !
  nden = 10
  !
  xmuon = 0._pr
  !
  if (coef%iquark.eq.0) then
     !
     den_step = 0.02
     !
     ! SM
     !
     xd = 0.0
     if (metaEos_lverb) write(*,*)"Symmetric matter"
     if (metaEos_lverb) write(*,'(a4,12a14)')"iden","density","E/A","P","P_test","cs^2","Esym","Esym2","K","mun","mup","u_n"
     !
     do iden = 1, nden
        !
        den = iden * den_step
        !
        call metaeos_T0_set_densities(den,xd,xmuon,coef,withq,eosDen)
        call metaeos_T0_baryons(coef,eosDen,.false.,withq,eosb0,eosq0)
        write(*,*)eosDen%den_nuc,eosb0%e2a_b
        write(*,*)eosDen%den_nuc,eosb0%p_b
        !
     enddo
     !
     ! NM
     !
     xd = 1.0 ! 0.99
     !
     if (metaEos_lverb) write(*,*)"Neutron matter"
     !
     do iden = 1, nden
        den = iden * den_step
        call metaeos_T0_set_densities(den,xd,xmuon,coef,withq,eosDen)
        call metaeos_T0_baryons(coef,eosDen,.false.,withq,eosb0,eosq0)
        write(*,*)eosDen%den_nuc,eosb0%e2a_b
        write(*,*)eosDen%den_nuc,eosb0%p_b
        !
     enddo
     !
  endif
  !
end subroutine metaeos_compute_T0_MCMC






subroutine metaeos_compute_T0_SM_NM(coef,withq)

! use acc; use cst;
!  use metaeosT0; use metaEosType;
  
  TYPE (metaEos_coef), intent(in) :: coef
  logical,             intent(in) :: withq 
!  integer, dimension(metaEos_neparam), intent(in)           :: eparam
  !
  ! local variables
  !
  TYPE (metaEos_densities)     :: eosDen
  TYPE (metaeos_Baryons)       :: eosb0
  TYPE (metaEos_Q1)            :: eosq0
  !
  ! define the max density where the physical filter is applied
  integer        :: nden
  real (kind=pr) :: den_step
  !
  real (kind=pr) :: kappas, kappav, kappasym
  real (kind=pr) :: den, xd, denn, denp, xmuon, xt, xx
  integer :: iden
  real (kind=pr), dimension(:,:,:), allocatable :: aeos, aeos2
  real (kind=pr) :: e00, p00, de, dp, cs, cs2
  real (kind=pr) :: ptest
  !
  real (kind=pr) :: xpre0, xrho0, cs_test, p_test
  !
  character (len=1) :: cxd
  !
  ! parameters
  !
  nden = 120
  !
  allocate(aeos(nden,0:1,12),aeos2(nden,0:1,12))
  !
  xmuon = 0._pr
  !
  e00 = 0._pr
  p00 = 0._pr
  !
  if (coef%iquark.eq.0) then
     !
     den_step = CST_nsat / 10._pr
     !
!     write(*,*)"no quark model"
     !
     ! SM
     !
     xd = 0.0
     if (metaEos_lverb) write(*,*)"Symmetric matter"
     if (metaEos_lverb) write(*,'(a4,12a14)')"iden","density","E/A","P","P_test","cs^2","Esym","Esym2","K","mun","mup","u_n"
     open(81,file="metaEos-res/nm-sm.dat",status="unknown")
     write(81,*)"# symmetric nucleonic matter with deltaN=",xd
     write(81,'(20a14)')"n_N","kFN","n_n","n_p","rho_B","E/A","Esym","Esym2","mu_n","P_B","K_B","(vs/c)^2"
     !
     do iden = 1, nden
        !
        den = iden * den_step
        !
        call metaeos_T0_set_densities(den,xd,xmuon,coef,withq,eosDen)
        call metaeos_T0_baryons(coef,eosDen,.false.,withq,eosb0,eosq0)
        aeos(iden,0,1) = eosb0%e2a_b
        aeos(iden,0,2) = eosb0%p_b
        aeos(iden,0,3) = eosb0%cs_b
        aeos(iden,0,4) = eosb0%esym_b
        aeos(iden,0,5) = eosb0%esym2_b
        aeos(iden,0,6) = eosb0%K_b
        aeos(iden,0,7) = eosb0%mu_n
        aeos(iden,0,8) = eosb0%mu_p
        aeos(iden,0,9) = eosb0%ms_n
        aeos(iden,0,10) = eosb0%ms_p
        aeos(iden,0,11)= eosb0%umf_n
        aeos(iden,0,12)= eosb0%umf_p
        ptest = -den*eosb0%e2a_b + eosDen%den_n * eosb0%mu_n + eosDen%den_p * eosb0%mu_p
        if (metaEos_lverb) write(*,'(i4,13f14.4)')iden,den,aeos(iden,0,1),aeos(iden,0,2),ptest,&
             !          & -den*eosb0%e2a_b,eosDen%den_n * eosb0%mu_n + eosDen%den_p * eosb0%mu_p,&
             & aeos(iden,0,3),aeos(iden,0,4),&
             & aeos(iden,0,5),aeos(iden,0,6),aeos(iden,0,7),aeos(iden,0,8),aeos(iden,0,11)
        write(81,'(20f14.3)')eosDen%den_nuc,eosDen%kF_nuc,eosDen%den_n,eosDen%den_p,&
             & eosb0%rho_b,eosb0%e2a_b,eosb0%esym_b,eosb0%esym2_b,eosb0%mu_n,eosb0%mu_np,eosb0%p_b,eosb0%K_b,eosb0%cs_b
        !
  enddo
  close(81)
     !
     ! NM
     !
     xd = 1.0 ! 0.99
     !
     if (metaEos_lverb) write(*,*)"Neutron matter"
     if (metaEos_lverb) write(*,'(a4,12a14)')"iden","density","E/A","P","P_test","cs^2","Esym","Esym2","K","mun","mup","u_n"
     open(81,file="metaEos-res/nm-nm.dat",status="unknown")
     write(81,*)"# neutron nucleonic matter with deltaN=",xd
     write(81,'(20a14)')"n_N","kFN","n_n","n_p","rho_B","E/A","Esym","Esym2","mu_n","P_B","K_B","(vs/c)^2"
     !
     do iden = 1, nden
        den = iden * den_step
        call metaeos_T0_set_densities(den,xd,xmuon,coef,withq,eosDen)
        call metaeos_T0_baryons(coef,eosDen,.false.,withq,eosb0,eosq0)
        aeos(iden,1,1) = eosb0%e2a_b
        aeos(iden,1,2) = eosb0%p_b
        aeos(iden,1,3) = eosb0%cs_b
        aeos(iden,1,4) = eosb0%esym_b
        aeos(iden,1,5) = eosb0%esym2_b
        !     aeos(iden,1,6) = eosb0%K_b-18/den*eosb0%p_b
        aeos(iden,1,6) = eosb0%K_b
        aeos(iden,1,7) = eosb0%mu_n
        aeos(iden,1,8) = eosb0%mu_p
        aeos(iden,1,9) = eosb0%ms_n
        aeos(iden,1,10) = eosb0%ms_p
        aeos(iden,1,11)= eosb0%umf_n
        aeos(iden,1,12)= eosb0%umf_p
        ptest = -den*eosb0%e2a_b + eosDen%den_n * eosb0%mu_n + eosDen%den_p * eosb0%mu_p
        if (metaEos_lverb) write(*,'(i4,11f14.4)')iden,den,aeos(iden,1,1),aeos(iden,1,2),ptest,&
             & aeos(iden,1,3),aeos(iden,1,4),&
             & aeos(iden,1,5),aeos(iden,1,6),aeos(iden,1,7),aeos(iden,1,8),aeos(iden,1,11)
        write(81,'(20f14.3)')eosDen%den_nuc,eosDen%kF_nuc,eosDen%den_n,eosDen%den_p,&
             & eosb0%rho_b,eosb0%e2a_b,eosb0%esym_b,eosb0%esym2_b,eosb0%mu_n,eosb0%mu_np,eosb0%p_b,eosb0%K_b,eosb0%cs_b
     enddo

  elseif (coef%iquark.ge.1) then
     !     write(*,*)"quarkyonic model"
     den_step = 3._pr*CST_nsat / 5._pr
     !
     xd = 0.99 ! 0.0, 0.5, 0.99
     write(cxd,'(i1)')int(xd*10)
     if (metaEos_lverb) write(*,*)"file:",cxd
     !
     if (metaEos_lverb) write(*,*)"Asymmetric matter with deltaN=",xd
     if (coef%iquark.eq.2) then
        open(80,file="metaEos-res/qm-qnm.dat",status="unknown")
        write(80,*)"# Asymmetric matter with deltaN=",xd
        write(80,'(20a14)')"n_B","n_N","n_Q","Delta","kFN","rho_B","E/A","mu_B","P_B","K_B","(vs/c)^2"
     elseif (coef%iquark.eq.5) then
        open(80,file="metaEos-res/qm-xd"//cxd//".dat",status="unknown")
        write(80,*)"# Asymmetric matter with deltaN=",xd
        write(80,'(20a14)')"n_B","n_N","n_Q","Delta","kFN","n_n","n_p","n_u","n_d","rho_B","E/A","mu_n","mu_np","P_B","K_B","(vs/c)^2"
     endif
        !     write(80,'(20a14)')"n_noq","n_N","n_Q","n_B","rho_noq","rho_B","E/A_nq","E/A","mu_n","mu_n_QP","mu_u","mu_u_qp","mu_B","pn","pn_QP","P_Q","P_B"
     if (metaEos_lverb) write(*,'(a4,20a14)')"iden","density","DeltaQCD","nB","nN","nQ","rhoB c^2","E/A","P","K","cs","mu_B","cs_test"
     !
     xrho0 = 0._pr
     xpre0 = 0._pr
     !
     do iden = 1, nden
        den = iden * den_step
        call metaeos_T0_set_densities(den,xd,xmuon,coef,withq,eosDen)
        call metaeos_T0_baryons(coef,eosDen,.false.,withq,eosb0,eosq0)
        cs_test=(eosq0%p_B-xpre0)/(eosq0%rho_B-xrho0)
        aeos(iden,1,1) = eosq0%rho_B
        aeos(iden,1,2) = eosq0%e2a
        aeos(iden,1,3) = eosq0%p_B
        aeos(iden,1,4) = eosq0%K_B
        aeos(iden,1,5) = eosq0%cs
!        aeos(iden,1,4) = eosb0%esym_b
!        aeos(iden,1,5) = eosb0%esym2_b
!        aeos(iden,1,6) = eosb0%K_b
!        aeos(iden,1,7) = eosb0%mu_n
!        aeos(iden,1,8) = eosb0%mu_p
!        aeos(iden,1,9) = eosb0%ms_n
!        aeos(iden,1,10) = eosb0%ms_p
!        aeos(iden,1,11)= eosb0%umf_n
        !        aeos(iden,1,12)= eosb0%umf_p
        !ptest = -eosq0%rho_B + eosDen%den_n * eosb0%mu_n + eosDen%den_p * eosb0%mu_p
     if (metaEos_lverb) write(*,'(i4,20f14.4)')iden,den,eosDen%q1_delta,eosDen%den_b,eosDen%den_nuc,&
             & eosDen%q1_nq,&
             & aeos(iden,1,1),aeos(iden,1,2),aeos(iden,1,3),aeos(iden,1,4),aeos(iden,1,5),eosq0%mu_B,cs_test

!        write(69,'(20(f14.2,a3))')den,"&",eosDen%q1_delta,"&",eosDen%den_n,"&",eosDen%den_p,"&",eosDen%den_nuc,&
!             &"&",eosDen%xd,"&",1000*eosDen%q1_nu,"&",1000*eosDen%q1_nd,"&",1000*eosDen%q1_nq,"&",eosDen%q1_xdq,"&",&
!             & eosDen%den_b,"&",eosDen%xd_b,"&",aeos(iden,1,1),"&", aeos(iden,1,2),"&", aeos(iden,1,3),"\\"

        if (coef%iquark.eq.2) then
           write(80,'(20f14.3)')eosDen%den_b,eosDen%den_nuc,100*eosDen%q1_nq,eosDen%q1_delta,eosDen%kF_nuc,eosq0%rho_B,&
                &eosq0%e2a,eosq0%mu_B,eosq0%p_B,eosq0%K_B,eosq0%cs
        elseif (coef%iquark.eq.5) then
           write(80,'(20f14.3)')eosDen%den_b,eosDen%den_nuc,100*eosDen%q1_nq,eosDen%q1_delta,eosDen%kF_nuc,&
                eosDen%den_n,eosDen%den_p,eosDen%q1_nu,eosDen%q1_nd,&
                &eosq0%rho_B,eosq0%e2a,eosq0%mun,eosq0%munp,eosq0%p_B,eosq0%K_B,eosq0%cs
        endif
        !
        xpre0 = eosq0%p_B
        xrho0 = eosq0%rho_B
        !
     enddo

     close(80)

  endif
  
     
end subroutine metaeos_compute_T0_SM_NM





subroutine metaeos_compute_T0_beta(coef,withq)

!  use acc; use cst;
!  use metaeosT0; use metaEosType;
  
!  integer, dimension(metaEos_neparam), intent(in)           :: eparam
  TYPE (metaEos_coef), intent(in) :: coef
  logical,             intent(in) :: withq 
  !
  ! local variables
  !
!  TYPE (metaEos_empirical)     :: emp
!  TYPE (metaEos_coef)          :: coef
  TYPE (metaEos_densities)     :: eosDen
  TYPE (metaEos_Baryons)       :: eosb0
  TYPE (metaEos_Q1)            :: eosq0
  TYPE (metaEos_Leptons)       :: eosl0
  TYPE (metaEos_Photons)       :: eosp0
  TYPE (metaEos_All)           :: eosa0
  !
  ! define the max density where the physical filter is applied
  integer        :: nden, nden_max, ibeta
  real (kind=pr) :: den_step
  !
  real (kind=pr) :: kappas, kappav, kappasym
  real (kind=pr) :: den, xd, denn, denp, xmuon, xt, xx
  integer :: iden
  real (kind=pr), dimension(:,:,:), allocatable :: aeos
  real (kind=pr), dimension(:,:), allocatable   :: aden
  !  real (kind=pr) :: e00, p00, de, dp, cs
  ! check of the sound velocity
  real (kind=pr) :: cs2, deps, dp
  real (kind=pr) :: ptest

  !
  ! parameters
  !
  nden = 120
  !
  allocate(aeos(nden,0:3,12))
  allocate(aden(nden,7))

  !
  xd = 0.9_pr
  xmuon = 0._pr
  !
  if (metaEos_lverb) write(*,*)"Beta-equilibrium with muons"
  if (metaEos_lverb) write(*,*)"Baryon quantities (no leptons)"
  if (coef%iquark.eq.0.and.metaEos_lverb) write(*,'(a4,14a14)')"iden","nB","delta","xmuon","E/A","P","P_test","cs^2","Esym","Esym2","K","mu_n","m^*_n","u_n"
  if (coef%iquark.ge.1.and.metaEos_lverb) write(*,'(a4,14a14)')"iden","nB","nN","delta","xmuon","E/A","P","cs^2","K","mu_n","munp","mu_e","mu_muon"
  !
  if (coef%iquark.eq.0) then
     den_step = CST_nsat / 10._pr
     open(80,file="metaEos-res/nm-beta.dat",status="unknown")
     write(80,*)"# Beta-eq. matter without neutrinos"
     write(80,'(20a14)')"n_N","kFN","n_n","n_p","rho_B","E/A","mu_n","mu_np","P_B","K_B","(vs/c)^2"
  elseif (coef%iquark.eq.5) then
     den_step = 3._pr * CST_nsat / 5._pr
     open(80,file="metaEos-res/qm-qbm.dat",status="unknown")
     write(80,*)"# Beta-eq. matter without neutrinos"
     write(80,'(20a14)')"n_B","n_N","n_Q","Delta","kFN","n_n","n_p","n_u","n_d","rho_B","E/A","mu_n","mu_np","P_B","K_B","(vs/c)^2"
  endif

  !e00 = 0._pr
  !p00 = 0._pr
  ! Calculate beta-equilibrium
  nden_max = nden
  do iden = 1, nden_max
     den = iden * den_step
     !     den = 0.44 + iden*den_step
!     write(*,*)"beta-eq:",den
     !
     ! Include or not muons
     !
     if (coef%imuon.eq.0) then
        ! calculate beta equilibrium with n, p, e only
!        write(*,*)"beta npe"
        call metaeos_T0_beta_npe(den,xd,coef,withq,eosDen,ibeta)
     else if (coef%imuon.eq.1) then
        ! calculate beta equilibrium with n, p, e, muon
!        write(*,*)"beta npemuon"
        call metaeos_T0_beta_npemuon(den,xd,xmuon,coef,withq,eosDen,ibeta)
     else
        stop "coef%imuon badly defined (different from 0 or 1)"
     endif
!     write(*,*)den,eosDen%xd,eosDen%xmuon
     if (ibeta.eq.0) then
        nden_max=iden-1
        write(*,*)"den_max",nden_max
        exit
     endif
     xd = eosDen%xd
     xmuon = eosDen%xmuon
     !
     aden(iden,1) = den
     aden(iden,2) = eosDen%xd
     aden(iden,3) = eosDen%xmuon
     aden(iden,4) = eosDen%den_n
     aden(iden,5) = eosDen%den_p
     aden(iden,6) = eosDen%den_e
     aden(iden,7) = eosDen%den_muon
     !
     ! Calculate Baryon properties
     !
     call metaeos_T0_baryons(coef,eosDen,.false.,withq,eosb0,eosq0)
     !
     ! Calculate Lepton properties
     !
     call metaeos_T0_leptons(coef,eosDen,.false.,eosl0)
     !
     ! Calculate Photon properties
     !
     call metaeos_T0_Photons(eosDen,.false.,eosp0)
     !
     ! Calculate global properties
     !
     call metaeos_T0_all(coef,eosDen,.false.,eosb0,eosq0,eosl0,eosp0,eosa0)         
     !
     ! Baryon contribution
     if (coef%iquark.eq.0) then
        aeos(iden,0,1) = eosb0%e2a_b
        aeos(iden,0,2) = eosb0%p_b
        aeos(iden,0,3) = eosb0%cs_b
        aeos(iden,0,4) = eosb0%esym_b
        aeos(iden,0,5) = eosb0%esym2_b
        aeos(iden,0,6) = eosb0%K_b
        aeos(iden,0,7) = eosb0%mu_n
        aeos(iden,0,8) = eosb0%mu_p
        aeos(iden,0,9) = eosb0%ms_n
        aeos(iden,0,10) = eosb0%ms_p
        aeos(iden,0,11)= eosb0%umf_n
        aeos(iden,0,12)= eosb0%umf_p
        ptest = -den*eosb0%e2a_b + eosDen%den_n * eosb0%mu_n + eosDen%den_p * eosb0%mu_p
     elseif (coef%iquark.ge.1) then
        aeos(iden,0,1) = eosq0%e2a
        aeos(iden,0,2) = eosq0%p_b
        aeos(iden,0,3) = eosq0%cs
        aeos(iden,0,4) = eosb0%esym_b
        aeos(iden,0,5) = eosb0%esym2_b
        aeos(iden,0,6) = eosq0%K_b
        aeos(iden,0,7) = eosq0%mun
        aeos(iden,0,8) = eosq0%mup
!        aeos(iden,0,9) = eosb0%ms_n
!        aeos(iden,0,10) = eosb0%ms_p
!        aeos(iden,0,11)= eosb0%umf_n
!        aeos(iden,0,12)= eosb0%umf_p
     endif
     if (coef%iquark.eq.0.and.metaEos_lverb) write(*,'(i4,14f14.4)')iden,eosDen%den_b,eosDen%xd,eosDen%xmuon,&
          & aeos(iden,0,1),aeos(iden,0,2),ptest,&
          & aeos(iden,0,3),aeos(iden,0,4),&
          & aeos(iden,0,5),aeos(iden,0,6),aeos(iden,0,7),aeos(iden,0,9),aeos(iden,0,11)
     if (coef%iquark.ge.1.and.metaEos_lverb) write(*,'(i4,18f14.4)')iden,eosDen%den_b,eosDen%den_nuc,eosDen%xd,&
          & eosDen%xmuon*eosDen%den_nuc/eosDen%den_b,aeos(iden,0,1),aeos(iden,0,2),&
          & aeos(iden,0,3),aeos(iden,0,6),aeos(iden,0,7),eosq0%munp,eosl0%mu_e,eosl0%mu_muon
     ! electron contribution
     aeos(iden,1,1) = eosl0%mu_e
     aeos(iden,1,2) = eosl0%e2a_e
     aeos(iden,1,3) = eosl0%e2a_e_UR
     aeos(iden,1,4) = eosl0%p_e
     aeos(iden,1,5) = eosl0%p_e_UR
     aeos(iden,1,6) = eosl0%K_e
     aeos(iden,1,7) = eosl0%K_e_UR

     ! muon contribution
     aeos(iden,2,1) = eosl0%mu_muon
     aeos(iden,2,2) = eosl0%e2a_muon
     aeos(iden,2,3) = eosl0%e2a_muon_UR
     aeos(iden,2,4) = eosl0%p_muon
     aeos(iden,2,5) = eosl0%p_muon_UR
     aeos(iden,2,6) = eosl0%K_muon
     aeos(iden,2,7) = eosl0%K_muon_UR

     ! total contribution
     aeos(iden,3,1) = eosa0%e2a
     aeos(iden,3,2) = eosa0%p
     aeos(iden,3,3) = eosa0%cs
     aeos(iden,3,4) = eosa0%gamma
     aeos(iden,3,5) = eosa0%rho
     aeos(iden,3,6) = eosa0%enthalpy
     aeos(iden,3,7) = eosl0%e2a_e
     aeos(iden,3,8) = eosl0%e2a_e_UR
     aeos(iden,3,9) = eosl0%nu_e

     if (coef%iquark.eq.0) then
        write(80,'(20f14.3)')eosDen%den_nuc,eosDen%kF_nuc,eosDen%den_n,eosDen%den_p,&
             &eosa0%rho,eosa0%e2a,eosb0%mu_n,eosb0%mu_np,eosa0%p,eosa0%K,eosa0%cs
     elseif (coef%iquark.eq.5) then
        write(80,'(20f14.3)')eosDen%den_b,eosDen%den_nuc,100*eosDen%q1_nq,eosDen%q1_delta,eosDen%kF_nuc,&
             eosDen%den_n,eosDen%den_p,eosDen%q1_nu,eosDen%q1_nd,&
             &eosa0%rho,eosa0%e2a,eosq0%mun,eosq0%munp,eosa0%p,eosa0%K,eosa0%cs
     endif

  enddo

  if (metaEos_lverb) write(*,*)"Beta-equilibrium: Electron contribution"
  if (metaEos_lverb) write(*,'(a4,14a14)')"iden","density","delta","mu","E/A","E/A UR","P","P UR","K","K UR"
  do iden = 1, nden_max
     if (metaEos_lverb) write(*,'(i4,10f14.4)')iden,aden(iden,1),aden(iden,2),aeos(iden,1,1),aeos(iden,1,2),&
          & aeos(iden,1,3),aeos(iden,1,4),aeos(iden,1,5),aeos(iden,1,6),aeos(iden,1,7)
  enddo

  if (metaEos_lverb) write(*,*)"Beta-equilibrium: Muon contribution"
  if (metaEos_lverb) write(*,'(a4,14a14)')"iden","density","x_muon","mu","E/A","E/A UR","P","P UR","K","K UR"
  do iden = 1, nden_max
     if (metaEos_lverb) write(*,'(i4,10f14.4)')iden,aden(iden,1),aden(iden,3),aeos(iden,2,1),aeos(iden,2,2),&
          & aeos(iden,2,3),aeos(iden,2,4),aeos(iden,2,5),aeos(iden,2,6),aeos(iden,2,7)
  enddo
  
  if (metaEos_lverb) write(*,*)"Total quantities (Baryons+Leptons)"
  if (metaEos_lverb) write(*,'(a4,14a14)')"iden","density","E/A","P","cs^2","cs^2(num.)","E/A(e)","E/A(e) UR"
  do iden = 1, nden_max
     den = iden * den_step
     dp = 0.5 * ( aeos(iden+1,3,2) - aeos(iden-1,3,2) )
!     deps = 0.5 * ( (den + den_step) * aeos(iden+1,1,1) - (den - den_step) * aeos(iden-1,1,1) )
     deps = 0.5 * ( aeos(iden+1,3,5) - aeos(iden-1,3,5) )
     cs2 = dp / deps
     if (metaEos_lverb) write(*,'(i4,10f14.4)')iden,den,aeos(iden,3,1),aeos(iden,3,2),aeos(iden,3,3),cs2,&
          & aeos(iden,3,7),aeos(iden,3,8)
  enddo

end subroutine metaeos_compute_T0_beta


  

! setup initialisation of all variables to 0
! define o_den0
! define o_denT
! define o_eosb0
! define o_eosbT
! define o_eosl0
! define o_eoslT
! define o_eostot

  subroutine metaeos_T0_set_densities(i_denb,i_xd,i_xmuon,i_coef,i_withq,o_eosDen)

!    use acc; use CST;
!    use metaEosType;

    implicit none

    ! input variables
    real (kind=pr), intent(in)            :: i_denb, i_xd, i_xmuon
    type (metaEos_coef), intent(in)       :: i_coef
    logical, intent(in)                   :: i_withq
    ! output variables
    type (metaEos_densities), intent(out) :: o_eosDen
    ! local variables
    real (kind=pr)                        :: epsn, epsp
    real (kind=pr)                        :: kf_u, kf_d
    !
    o_eosDen%den_noq = i_denb
    !
    ! Nucleonic model (no quark)
    if (i_coef%iquark.eq.0) then
       ! densities
       o_eosDen%den_nuc    = i_denb
       o_eosDen%kf_nuc    = ( 1.5_pr * CST_pi2 * o_eosDen%den_nuc )**CST_p13
       ! nucleon fractions
       o_eosDen%xd = i_xd
       o_eosDen%xp = ( 1._pr - i_xd ) / 2._pr
       o_eosDen%xn = 1._pr - o_eosDen%xp
       ! nucleon densities
       o_eosDen%den_p    = o_eosDen%xp * o_eosDen%den_nuc
       o_eosDen%den_n    = o_eosDen%xn * o_eosDen%den_nuc
       ! Fermi momenta (at zero T of course!)
       o_eosDen%kf_n    = ( 3.0_pr * CST_pi2 * o_eosDen%den_n )**CST_p13
       o_eosDen%kf_p    = ( 3.0_pr * CST_pi2 * o_eosDen%den_p )**CST_p13
       !
       ! quark densities
       o_eosDen%q1_kf_qd = 0._pr
       o_eosDen%q1_kf_qu = 0._pr
       o_eosDen%q1_nd = 0._pr
       o_eosDen%q1_nu = 0._pr
       o_eosDen%q1_nq = 0._pr
       o_eosDen%q1_xdq = 0._pr
       !
       ! total Baryon densities
       o_eosDen%xd_b = o_eosDen%xd
       o_eosDen%den_b = o_eosDen%den_nuc
       !
    ! Quarkyonic model in NM (PRL 122, 122701)
    ! Symmetric matter
    elseif (i_coef%iquark.eq.1) then
       !
       o_eosDen%xd = 0.0_pr
       o_eosDen%xp = 0.5_pr
       o_eosDen%xn = 0.5_pr
       !
       ! Set kFn and kFp (for given i_denb and i_xd)
       !
       o_eosDen%kf_nuc    = ( 3._pr / 2._pr * CST_pi2 * i_denb )**CST_p13
       o_eosDen%kf_n    = o_eosDen%kf_nuc
       o_eosDen%kf_p    = o_eosDen%kf_nuc
       !
       ! gap delta:
       o_eosDen%q1_delta = i_coef%LQCD**3 / ( CST_hbc**3 * o_eosDen%kf_nuc**2 ) + &
            & i_coef%kQCD * i_coef%LQCD / ( CST_hbc * CST_Nc**2 )
       ! nucleon density
       o_eosDen%den_nuc = 2._pr / ( 3._pr * CST_pi2) * o_eosDen%kf_nuc**3
       if ((o_eosDen%kf_nuc - o_eosDen%q1_delta).gt.0._pr.and.i_withq) then
          o_eosDen%den_nuc = o_eosDen%den_nuc - 2._pr / ( 3._pr * CST_pi2) * ( o_eosDen%kf_nuc - o_eosDen%q1_delta )**3
       endif
       o_eosDen%den_p    = o_eosDen%xp * o_eosDen%den_nuc
       o_eosDen%den_n    = o_eosDen%xn * o_eosDen%den_nuc
       ! quark density
       ! u and d quarks
       o_eosDen%q1_kf_qd = 0._pr
       o_eosDen%q1_kf_qu = 0._pr
       if ((o_eosDen%kf_nuc - o_eosDen%q1_delta).gt.0._pr.and.i_withq) then
          o_eosDen%q1_kf_qd = (o_eosDen%kf_nuc - o_eosDen%q1_delta) / CST_Nc
          o_eosDen%q1_kf_qu = o_eosDen%q1_kf_qd
       endif
       o_eosDen%q1_nd = 1._pr / ( 3._pr * CST_pi2) * o_eosDen%q1_kf_qd**3
       o_eosDen%q1_nu = 1._pr / ( 3._pr * CST_pi2) * o_eosDen%q1_kf_qu**3
       ! u and d quarks
       o_eosDen%q1_nq = o_eosDen%q1_nu + o_eosDen%q1_nd
       o_eosDen%q1_xdq = 0._pr
       if (o_eosDen%q1_nq.gt.0._pr) o_eosDen%q1_xdq = ( o_eosDen%q1_nd - o_eosDen%q1_nu ) / o_eosDen%q1_nq
       !
       ! Total baryon number density
       !
       o_eosDen%den_b = o_eosDen%den_nuc + o_eosDen%q1_nq
       o_eosDen%xd_b = ( o_eosDen%den_nuc * o_eosDen%xd + CST_Nc * o_eosDen%q1_nq * o_eosDen%q1_xdq ) / o_eosDen%den_b
       !
    ! Quarkyonic model in NM (PRL 122, 122701)
    ! Neutron matter
    elseif (i_coef%iquark.eq.2) then
       !
       o_eosDen%xd = 1.0_pr
       o_eosDen%xp = 0._pr
       o_eosDen%xn = 1._pr
       !
       ! Set kFn and kFp (for given i_denb and i_xd)
       !
       o_eosDen%kf_nuc    = ( 3._pr * CST_pi2 * i_denb )**CST_p13
       o_eosDen%kf_n    = o_eosDen%kf_nuc
       o_eosDen%kf_p    = 0._pr
       !
       ! gap delta:
       o_eosDen%q1_delta = i_coef%LQCD**3 / ( CST_hbc**3 * o_eosDen%kf_nuc**2 ) + &
            & i_coef%kQCD * i_coef%LQCD / ( CST_hbc * CST_Nc**2 )
       ! nucleon density
       o_eosDen%den_nuc = 1._pr / ( 3._pr * CST_pi2) * o_eosDen%kf_nuc**3
       if ((o_eosDen%kf_nuc - o_eosDen%q1_delta).gt.0._pr.and.i_withq) then
          o_eosDen%den_nuc = o_eosDen%den_nuc - 1._pr / ( 3._pr * CST_pi2) * ( o_eosDen%kf_nuc - o_eosDen%q1_delta )**3
       endif
       o_eosDen%den_p    = o_eosDen%xp * o_eosDen%den_nuc
       o_eosDen%den_n    = o_eosDen%xn * o_eosDen%den_nuc
       ! Fermi momenta (at zero T of course!)
       o_eosDen%kf_n    = o_eosDen%kf_nuc
       o_eosDen%kf_p    = 0._pr
       ! quark density
       ! d quark
       o_eosDen%q1_kf_qd = 0._pr
       if ((o_eosDen%kf_n - o_eosDen%q1_delta).gt.0._pr.and.i_withq) then
          o_eosDen%q1_kf_qd = (o_eosDen%kf_n - o_eosDen%q1_delta) / CST_Nc
       endif
       o_eosDen%q1_nd = 1._pr / ( 3._pr * CST_pi2) * o_eosDen%q1_kf_qd**3
       ! u quark
       o_eosDen%q1_kf_qu = o_eosDen%q1_kf_qd / 2**CST_p13
       o_eosDen%q1_nu = 1._pr / ( 3._pr * CST_pi2) * o_eosDen%q1_kf_qu**3
       ! u and d quarks
       o_eosDen%q1_nq = o_eosDen%q1_nu + o_eosDen%q1_nd
       o_eosDen%q1_xdq = 0._pr
       if (o_eosDen%q1_nq.gt.0._pr) o_eosDen%q1_xdq = ( o_eosDen%q1_nd - o_eosDen%q1_nu ) / o_eosDen%q1_nq
       !
       ! Total baryon number density
       !
       o_eosDen%den_b = o_eosDen%den_nuc + o_eosDen%q1_nq
       o_eosDen%xd_b = ( o_eosDen%den_nuc * o_eosDen%xd + CST_Nc * o_eosDen%q1_nq * o_eosDen%q1_xdq ) / o_eosDen%den_b
       !
    ! Isovector Quarkyonic model in asymmetric matter: non-conservation of isospin number
    elseif (i_coef%iquark.eq.3) then
       !
       ! Set kFn and kFp (for given i_denb and i_xd)
       !
       o_eosDen%kf_nuc    = ( 1.5_pr * CST_pi2 * i_denb )**CST_p13
       o_eosDen%kf_n    = ( 1._pr + i_xd )**CST_p13 * o_eosDen%kf_nuc
       o_eosDen%kf_p    = ( 1._pr - i_xd )**CST_p13 * o_eosDen%kf_nuc
       !
       ! Calculate Quarkyonic model
       !
       ! gap delta:
       o_eosDen%q1_delta = i_coef%LQCD**3 / CST_hbc**3 / ( o_eosDen%kf_n**3 + o_eosDen%kf_p**3 )**CST_p23 + &
            & i_coef%kQCD * i_coef%LQCD / ( CST_hbc * CST_Nc**2 )
       !
       ! Nucleon density
       !
       ! neutron
       o_eosDen%den_n = o_eosDen%kf_n**3 / ( 3._pr * CST_pi2)
       if ((o_eosDen%kf_n - o_eosDen%q1_delta).gt.0._pr.and.i_withq) then
          o_eosDen%den_n = o_eosDen%den_n - ( o_eosDen%kf_n - o_eosDen%q1_delta )**3 / ( 3._pr * CST_pi2)
       endif
       ! proton
       o_eosDen%den_p = o_eosDen%kf_p**3 / ( 3._pr * CST_pi2)
       if ((o_eosDen%kf_p - o_eosDen%q1_delta).gt.0._pr.and.i_withq) then
          o_eosDen%den_p = o_eosDen%den_p - ( o_eosDen%kf_p - o_eosDen%q1_delta )**3 / ( 3._pr * CST_pi2)
       endif
       o_eosDen%den_nuc = o_eosDen%den_n + o_eosDen%den_p
       o_eosDen%xn = o_eosDen%den_n / o_eosDen%den_nuc
       o_eosDen%xp = o_eosDen%den_p / o_eosDen%den_nuc
       o_eosDen%xd = 1._pr - 2._pr * o_eosDen%xp
       !
       ! Quark densities
       !
       epsn = CST_mnuc2
       if ((o_eosDen%kf_n - o_eosDen%q1_delta).gt.0._pr.and.i_withq) then
          epsn = dsqrt( CST_mnuc2**2 + CST_hbc**2 * ( o_eosDen%kf_n - o_eosDen%q1_delta)**2)
       endif
       epsp = CST_mnuc2
       if ((o_eosDen%kf_p - o_eosDen%q1_delta).gt.0._pr.and.i_withq) then
          epsp = dsqrt( CST_mnuc2**2 + CST_hbc**2 * (o_eosDen%kf_p - o_eosDen%q1_delta)**2)
       endif
       o_eosDen%q1_kf_qu = 0._pr
       if ((dabs(2*epsp-epsn)-CST_Nc*CST_mqc2).gt.0._pr.and.i_withq) then
          o_eosDen%q1_kf_qu = dsqrt( ( 2 * epsp - epsn )**2 - ( CST_Nc * CST_mqc2 )**2 ) / ( CST_Nc * CST_hbc )
       endif
       o_eosDen%q1_kf_qd = 0._pr
       if ((dabs(2*epsn-epsp)-CST_Nc*CST_mqc2).gt.0._pr.and.i_withq) then
          o_eosDen%q1_kf_qd = dsqrt( ( 2 * epsn - epsp )**2 - ( CST_Nc * CST_mqc2 )**2 ) / ( CST_Nc * CST_hbc )
       endif
!       write(*,*)i_denb,epsn,epsp,2*epsp-epsn,2*epsn-epsp
       !
       o_eosDen%q1_nu = 1._pr / ( 3._pr * CST_pi2) * o_eosDen%q1_kf_qu**3
       o_eosDen%q1_nd = 1._pr / ( 3._pr * CST_pi2) * o_eosDen%q1_kf_qd**3
       o_eosDen%q1_nq = o_eosDen%q1_nu + o_eosDen%q1_nd
       o_eosDen%q1_xdq = 0._pr
       if (o_eosDen%q1_nq.gt.0._pr) o_eosDen%q1_xdq = ( o_eosDen%q1_nd - o_eosDen%q1_nu ) / o_eosDen%q1_nq
!       write(*,*)i_denb,CST_hbc*o_eosDen%kf_n,o_eosDen%den_n,o_eosDen%den_b,o_eosDen%q1_kf_qu,o_eosDen%q1_nu
       !
       ! Total baryon number density
       !
       o_eosDen%den_b = o_eosDen%den_nuc + o_eosDen%q1_nq
       o_eosDen%xd_b = ( o_eosDen%den_nuc * o_eosDen%xd + CST_Nc * o_eosDen%q1_nq * o_eosDen%q1_xdq ) / o_eosDen%den_b
       !
    ! Isovector Quarkyonic model in asymmetric matter: conservation of isospin number
    elseif (i_coef%iquark.eq.4) then
       !
       ! Set kFn and kFp (for given i_denb and i_xd)
       !
       o_eosDen%kf_nuc    = ( 1.5_pr * CST_pi2 * i_denb )**CST_p13
       o_eosDen%kf_n    = ( 1._pr + i_xd )**CST_p13 * o_eosDen%kf_nuc
       o_eosDen%kf_p    = ( 1._pr - i_xd )**CST_p13 * o_eosDen%kf_nuc
       !
       ! Calculate Quarkyonic model
       !
       ! gap delta:
       o_eosDen%q1_delta = i_coef%LQCD**3 / CST_hbc**3 / ( o_eosDen%kf_n**3 + o_eosDen%kf_p**3 )**CST_p23 + &
            & i_coef%kQCD * i_coef%LQCD / ( CST_hbc * CST_Nc**2 )
       !
       ! Nucleon density
       !
       ! neutron
       o_eosDen%den_n = o_eosDen%kf_n**3 / ( 3._pr * CST_pi2)
       if ((o_eosDen%kf_n - o_eosDen%q1_delta).gt.0._pr.and.i_withq) then
          o_eosDen%den_n = o_eosDen%den_n - ( o_eosDen%kf_n - o_eosDen%q1_delta )**3 / ( 3._pr * CST_pi2)
       endif
       ! proton
       o_eosDen%den_p = o_eosDen%kf_p**3 / ( 3._pr * CST_pi2)
       if ((o_eosDen%kf_p - o_eosDen%q1_delta).gt.0._pr.and.i_withq) then
          o_eosDen%den_p = o_eosDen%den_p - ( o_eosDen%kf_p - o_eosDen%q1_delta )**3 / ( 3._pr * CST_pi2)
       endif
       o_eosDen%den_nuc = o_eosDen%den_n + o_eosDen%den_p
       o_eosDen%xn = o_eosDen%den_n / o_eosDen%den_nuc
       o_eosDen%xp = o_eosDen%den_p / o_eosDen%den_nuc
       o_eosDen%xd = 1._pr - 2._pr * o_eosDen%xp
       !
       ! Quark densities
       !
       epsn = CST_mnuc2
       if ((o_eosDen%kf_n - o_eosDen%q1_delta).gt.0._pr.and.i_withq) then
          epsn = dsqrt( CST_mnuc2**2 + CST_hbc**2 * (o_eosDen%kf_n - o_eosDen%q1_delta)**2)
       endif
       epsp = CST_mnuc2
       if ((o_eosDen%kf_p - o_eosDen%q1_delta).gt.0._pr.and.i_withq) then
          epsp = dsqrt( CST_mnuc2**2 + CST_hbc**2 * (o_eosDen%kf_p - o_eosDen%q1_delta)**2)
       endif
       o_eosDen%q1_kf_qu = 0._pr
       if ((dabs(2*epsp-epsn)-CST_Nc*CST_mqc2).gt.0._pr.and.i_withq) then
          o_eosDen%q1_kf_qu = dsqrt( ( 2 * epsp - epsn )**2 - ( CST_Nc * CST_mqc2 )**2 ) / ( CST_Nc * CST_hbc )
       endif
       o_eosDen%q1_kf_qd = 0._pr
       if ((dabs(2*epsn-epsp)-CST_Nc*CST_mqc2).gt.0._pr.and.i_withq) then
          o_eosDen%q1_kf_qd = dsqrt( ( 2 * epsn - epsp )**2 - ( CST_Nc * CST_mqc2 )**2 ) / ( CST_Nc * CST_hbc )
       endif
       ! impose the conservation of isospin number
       if (o_eosDen%q1_kf_qu.gt.o_eosDen%q1_kf_qd) then
          o_eosDen%q1_kf_qd = ( ( CST_Nc + o_eosDen%xd ) / ( CST_Nc - o_eosDen%xd ) )**CST_p13 * o_eosDen%q1_kf_qu
       else
          o_eosDen%q1_kf_qu = ( ( CST_Nc - o_eosDen%xd ) / ( CST_Nc + o_eosDen%xd ) )**CST_p13 * o_eosDen%q1_kf_qd
       endif
!       write(*,*)i_denb,epsn,epsp,2*epsp-epsn,2*epsn-epsp
       !
       o_eosDen%q1_nu = 1._pr / ( 3._pr * CST_pi2) * o_eosDen%q1_kf_qu**3
       o_eosDen%q1_nd = 1._pr / ( 3._pr * CST_pi2) * o_eosDen%q1_kf_qd**3
       o_eosDen%q1_nq = o_eosDen%q1_nu + o_eosDen%q1_nd
       o_eosDen%q1_xdq = 0._pr
       if (o_eosDen%q1_nq.gt.0._pr) o_eosDen%q1_xdq = ( o_eosDen%q1_nd - o_eosDen%q1_nu ) / o_eosDen%q1_nq
!       write(*,*)i_denb,CST_hbc*o_eosDen%kf_n,o_eosDen%den_n,o_eosDen%den_b,o_eosDen%q1_kf_qu,o_eosDen%q1_nu
       !
       ! Total baryon number density
       !
       o_eosDen%den_b = o_eosDen%den_nuc + o_eosDen%q1_nq
       o_eosDen%xd_b = ( o_eosDen%den_nuc * o_eosDen%xd + CST_Nc * o_eosDen%q1_nq * o_eosDen%q1_xdq ) / o_eosDen%den_b
       !
    ! Isoscalar Quarkyonic model in asymmetric matter
    elseif (i_coef%iquark.eq.5) then
       !
       ! Set kFb (for given i_denb and i_xd)
       !
       o_eosDen%kf_nuc    = ( 3._pr / 2._pr * CST_pi2 * i_denb )**CST_p13
       o_eosDen%xd = i_xd
       o_eosDen%xp = ( 1._pr - i_xd ) / 2._pr
       o_eosDen%xn = 1._pr - o_eosDen%xp
       o_eosDen%q1_xdo = 0._pr
       o_eosDen%q1_xup = 0._pr
       !
       o_eosDen%den_nuc   = i_denb
       o_eosDen%den_p   = o_eosDen%xp * o_eosDen%den_nuc
       o_eosDen%den_n   = o_eosDen%xn * o_eosDen%den_nuc
       !
       ! Fermi momenta (at zero T of course!)
       o_eosDen%kf_n    = ( 3.0_pr * CST_pi2 * o_eosDen%den_n )**CST_p13
       o_eosDen%kf_p    = ( 3.0_pr * CST_pi2 * o_eosDen%den_p )**CST_p13
       !
       ! gap delta:
       !
       o_eosDen%q1_delta = i_coef%LQCD**3 / ( CST_hbc**3 * o_eosDen%kf_nuc**2 ) + &
            & i_coef%kQCD * i_coef%LQCD / ( CST_hbc * CST_Nc**2)

!       o_eosDen%q1_delta = i_coef%LQCD**3 / ( CST_hbc**3 * 2.5**2 ) + &
!            & i_coef%kQCD * i_coef%LQCD / ( CST_hbc * CST_Nc**2)
       !
       o_eosDen%q1_kf_q = 0._pr ;
       o_eosDen%q1_nq   = 0._pr ;
       o_eosDen%q1_xdo  = 0._pr ;
       o_eosDen%q1_xup = 0._pr ;
       o_eosDen%q1_nd   = 0._pr ;
       o_eosDen%q1_nu   = 0._pr ;
       o_eosDen%q1_xdq  = 0._pr ;
       o_eosDen%q1_kf_qd = 0._pr;
       o_eosDen%q1_kf_qu = 0._pr;
       o_eosDen%kf_n_min = 0._pr;
       o_eosDen%kf_p_min = 0._pr;
       !
       ! case quark are present
       !
       if ((o_eosDen%kf_nuc - o_eosDen%q1_delta).gt.0._pr.and.i_withq) then
          ! quark densities
          o_eosDen%q1_kf_q = (o_eosDen%kf_nuc - o_eosDen%q1_delta) / CST_Nc
          o_eosDen%q1_nq = 2._pr / ( 3._pr * CST_pi2 ) * o_eosDen%q1_kf_q**3
          o_eosDen%q1_xdo = ( 1._pr + i_xd/CST_Nc ) / 2._pr
          o_eosDen%q1_xup = ( 1._pr - i_xd/CST_Nc ) / 2._pr
          o_eosDen%q1_nd = o_eosDen%q1_xdo * o_eosDen%q1_nq
          o_eosDen%q1_nu = o_eosDen%q1_xup * o_eosDen%q1_nq
          o_eosDen%q1_xdq = i_xd/CST_Nc
          o_eosDen%q1_kf_qd = ( 3.0_pr * CST_pi2 * o_eosDen%q1_nd )**CST_p13
          o_eosDen%q1_kf_qu = ( 3.0_pr * CST_pi2 * o_eosDen%q1_nu )**CST_p13
          !write(*,*)"!",o_eosDen%kf_nuc,o_eosDen%q1_delta,o_eosDen%q1_nq,o_eosDen%q1_kf_q
          ! correction to the nucleon densities
          o_eosDen%den_nuc = 2._pr/(3.*CST_pi2) * ( o_eosDen%kf_nuc**3 - (o_eosDen%kf_nuc-o_eosDen%q1_delta)**3 )
          o_eosDen%den_p   = o_eosDen%xp * o_eosDen%den_nuc
          o_eosDen%den_n   = o_eosDen%xn * o_eosDen%den_nuc
          o_eosDen%kf_n = ( 2*o_eosDen%xn )**CST_p13 * o_eosDen%kf_nuc
          o_eosDen%kf_p = ( 2*o_eosDen%xp )**CST_p13 * o_eosDen%kf_nuc
          o_eosDen%kf_n_min = (o_eosDen%kf_nuc - o_eosDen%q1_delta) * (2*o_eosDen%xn)**CST_p13
          o_eosDen%kf_p_min = (o_eosDen%kf_nuc - o_eosDen%q1_delta) * (2*o_eosDen%xp)**CST_p13
       endif
       !
       ! Total baryon number density
       !
       o_eosDen%den_b = o_eosDen%den_nuc + o_eosDen%q1_nq
       o_eosDen%xd_b = ( o_eosDen%den_nuc * o_eosDen%xd + CST_Nc * o_eosDen%q1_nq * o_eosDen%q1_xdq ) / o_eosDen%den_b
       !
    endif
    !
    o_eosDen%xx       = ( o_eosDen%den_nuc - i_coef%nsat ) / ( 3.0 * i_coef%nsat )
    ! lepton densities
    !
    ! electroneutrality condition with muons
    o_eosDen%xe = o_eosDen%xp - i_xmuon
    ! electroneutrality condition with muons and quarks
    if (i_coef%iquark.eq.5) o_eosDen%xe = o_eosDen%xp - i_xmuon + ( 2. * o_eosDen%q1_nu - o_eosDen%q1_nd ) / o_eosDen%den_nuc
    !
    o_eosDen%xmuon = i_xmuon
    o_eosDen%den_e    = o_eosDen%xe    * o_eosDen%den_nuc
    o_eosDen%den_muon = o_eosDen%xmuon * o_eosDen%den_nuc
    o_eosDen%kf_e    = ( 3.0_pr * CST_pi2 * o_eosDen%den_e )**CST_p13
    o_eosDen%kf_muon = ( 3.0_pr * CST_pi2 * o_eosDen%den_muon )**CST_p13

    !
  end subroutine metaeos_T0_set_densities









! -------------------------------------------------------------
! This routine calculate the meta equation of state for given
! values of the densities given in eosDen
! Results in o_eosb0
! -------------------------------------------------------------

  subroutine metaeos_T0_baryons(i_coef,i_eosDen,i_opt_beta,i_withq,o_eosb0,o_eosq0)

    implicit none

    ! input variables
    type (metaEos_coef),      intent(in)  :: i_coef
    type (metaEos_densities), intent(in)  :: i_eosDen
    logical,                  intent(in)  :: i_opt_beta
    logical,                  intent(in)  :: i_withq
    ! output variables
    type (metaEos_Baryons),   intent(out) :: o_eosb0
    type (metaEos_Q1),        intent(out) :: o_eosq0

    if (i_coef%iquark.eq.0) then
       if (i_coef%irel.eq.0) call metaeos_T0_baryons_NR(i_coef,i_eosDen,i_opt_beta,o_eosb0)
       if (i_coef%irel.eq.1) call metaeos_T0_baryons_R(i_coef,i_eosDen,i_opt_beta,o_eosb0)
    elseif (i_coef%iquark.ge.1) then
       call metaeos_T0_baryons_R(i_coef,i_eosDen,i_opt_beta,o_eosb0)
       call metaeos_T0_q1(i_coef,i_eosDen,i_opt_beta,i_withq,o_eosb0,o_eosq0)
    endif
       
  end subroutine metaeos_T0_baryons






  
  subroutine metaeos_T0_q1(i_coef,i_eosDen,i_opt_beta,i_withq,i_eosb0,o_eosq0)

    implicit none

    ! input variables
    type (metaEos_coef), intent(in)       :: i_coef
    type (metaEos_densities), intent(in)  :: i_eosDen
    logical,                  intent(in)  :: i_opt_beta
    logical,                  intent(in)  :: i_withq
    type (metaEos_Baryons),   intent(in)  :: i_eosb0
    ! output variables
    type (metaEos_Q1), intent(out)   :: o_eosq0
    ! local variables
    real (kind=pr)                   :: coefC
    real (kind=pr)                   :: zN1, zN2, zP1, zP2, zQd, zQu, zQ
    real (kind=pr)                   :: fN1, fN2, fP1, fP2, fQd, fQu
    real (kind=pr)                   :: d1kFN_rhokin, d1nB_kFN, d1kFN_nN, d1kFN_nQ, d1kFN_Delta, d1kFN_pot
    real (kind=pr)                   :: d2kFN_rhokin, d2nB_kFN, d2kFN_nN, d2kFN_nQ, d2kFN_Delta, d2kFN_pot
    real (kind=pr)                   :: C_2, C_4
    real (kind=pr)                   :: d1z_fQu, d1z_fQd, d1z_fN1, d1z_fN2, d1z_fP1, d1z_fP2
    real (kind=pr)                   :: d2z_fQu, d2z_fQd, d2z_fN1, d2z_fN2, d2z_fP1, d2z_fP2
    real (kind=pr)                   :: fN, fQ
    real (kind=pr)                   :: xat, xbt
    real (kind=pr)                   :: d1kFN_rhoN, d1nN_rhoN, d1nN_kFN, d1deltaN_rhoN, d1nN_pot, d1deltaN_pot
    real (kind=pr)                   :: d1kFN_kFQ, d1kFQ_rhoQ, d1nQ_rhoQ, d1nQ_kFQ, d1deltaQ_rhoQ
    real (kind=pr)                   :: d1nB_rhoN, d1nB_pot, d1nB_rhoQ, d1kFN_rhoQ
    real (kind=pr)                   :: d2kFN_rhoN, d2kFQ_rhoQ, d2kFN_kFQ, d2kFN_rhoQ
    real (kind=pr)                   :: d2kFN_rhoB, d2nB_rhoB
    real (kind=pr) :: d1kFN_nN1, d1kFN_nN2
    real (kind=pr) :: xtest

    !
    ! Coeficcient for the relativitic Fermi gas Energy
    coefC = 1.0_pr / (8.0_pr * CST_pi2 * CST_hbc**3)
    if (i_coef%ikin.eq.0) then
       coefC = 0.0_pr
    endif
    C_2 = coefC * CST_hbc
    C_4 = 2._pr * C_2
    !
    ! energy-density:
    !
    ! Symmetric matter
    if (i_coef%iquark.eq.1) then
       ! baryon contribution to the energy density
       zN1 = CST_hbc * i_eosDen%kf_n / CST_mnuc2
       fN1 = 2*coefC * (CST_mnuc2)**4 * ( zN1*(1._pr+2._pr*zN1**2)*dsqrt(1._pr+zN1**2) - dlog(zN1+dsqrt(1._pr+zN1**2)) )
       zN2 = 0._pr
       fN2 = 0._pr
       if (i_eosDen%q1_kf_qd.gt.0._pr.and.i_withq) then
!          zb2 = CST_hbc * CST_Nc * i_eosDen%q1_kf_q / CST_mnuc2
          zN2 = CST_hbc * ( i_eosDen%kf_n - i_eosDen%q1_delta ) / CST_mnuc2
          fN2 = 2*coefC * (CST_mnuc2)**4 * ( zN2*(1._pr+2._pr*zN2**2)*dsqrt(1._pr+zN2**2) - dlog(zN2+dsqrt(1._pr+zN2**2)) )
       endif
       fN = fN1 - fN2
       ! quark contribution to the energy density
       fQ = 0._pr
       if (i_eosDen%q1_kf_qd.gt.0._pr.and.i_withq) then
          zQ = CST_hbc * i_eosDen%q1_kf_qd / CST_mqc2
          fQ = CST_Nc * 2*coefC * (CST_mqc2)**4 * ( zQ*(1._pr+2._pr*zQ**2)*dsqrt(1._pr+zQ**2) - dlog(zQ+dsqrt(1._pr+zQ**2)) )
       endif
       !write(*,*)"test",fN1
       !    write(*,*)"test",i_eosDen%den_nuc,fq,fb1,fb2
       o_eosq0%rho_kin_B = fN + fQ
       !
       ! Potential energy-density in symmetric matter
       !
       !write(*,*)"test",i_coef%param(0,0),i_coef%param(0,2),i_eosDen%xx
       o_eosq0%rho_pot = ( i_coef%VCOEF(0,0) + i_coef%VCOEF(0,1) * i_eosDen%xx &
            & + 0.5 * i_coef%VCOEF(0,2) * i_eosDen%xx**2 ) * i_eosDen%den_nuc
!       write(*,*)"test rho_pot",i_eosDen%den_nuc,o_eosq0%rho_pot,i_eosb0%epot2v_b
       o_eosq0%rho_pot = i_eosb0%epot2v_b
       !
       o_eosq0%rho_B = o_eosq0%rho_kin_B + o_eosq0%rho_pot
       !write(*,*)"test",o_eosq0%rho_kin,o_eosq0%rho_pot,i_eosDen%den_b * CST_mnuc2
       !
       ! energy-density and energy per particle (without the rest mass)
       !
       o_eosq0%e2v = o_eosq0%rho_B - i_eosDen%den_b * CST_mnuc2
       o_eosq0%e2a = o_eosq0%e2v / i_eosDen%den_b
       !
       ! chemical potential
       !
       d1kFN_nN = 2 * i_eosDen%kf_n**2/CST_pi2
       !d1kFN_pot = 0._pr
       !d1kFN_pot = (2._pr*xat * (i_eosDen%den_n/0.16) + 3._pr*xbt * (i_eosDen%den_n/0.16)**2)*d1kFN_nN
       d1z_fN1 = 8._pr*zN1**2*dsqrt(1._pr+zN1**2)
       d1kFN_rhokin = C_4*CST_mnuc2**3*d1z_fN1
       d1nB_kFN = 1/d1kFN_nN 
       if (i_eosDen%q1_kf_qd.gt.0._pr.and.i_withq) then
          d1kFN_Delta = -2._pr*i_coef%LQCD**3 / CST_hbc**3 / i_eosDen%kf_n**3 
          d1kFN_nN =  2 * ( i_eosDen%kf_n**2 - (i_eosDen%kf_n - i_eosDen%q1_delta)**2*(1-d1kFN_Delta))/CST_pi2
          !d1kFN_pot = 0._pr
          !d1kFN_pot = (2._pr*xat * (i_eosDen%den_n/0.16) + 3._pr*xbt * (i_eosDen%den_n/0.16)**2)*d1kFN_nN
          d1kFN_nQ = (2._pr*i_eosDen%q1_kf_qd**2*(1-d1kFN_Delta))/CST_pi2/CST_Nc
          d1z_fQd = 8._pr*zQd**2*dsqrt(1._pr+zQd**2)
          d1z_fN2 = 8._pr*zN2**2*dsqrt(1._pr+zN2**2)
          d1kFN_rhokin = d1kFN_rhokin + C_4*(-CST_mnuc2**3*d1z_fN2 + CST_mqc2**3*d1z_fQd)*(1-d1kFN_Delta)
          d1nB_kFN = 1/(d1kFN_nN + d1kFN_nQ)
       endif
       o_eosq0%mu_kin_B = d1kFN_rhokin*d1nB_kFN
       !
       ! pressure
       !
       o_eosq0%p_kin_B =  i_eosDen%den_b * o_eosq0%mu_kin_B - o_eosq0%rho_kin_B
       o_eosq0%p_B = o_eosq0%p_kin_B ! + o_eosq0%p_pot
       !
       ! incompressibility
       !
       o_eosq0%k_kin = 0._pr
       !
       ! sound velocity
       !
       o_eosq0%cs = 0._pr
       !
    ! neutron matter
    elseif (i_coef%iquark.eq.2) then
       ! baryon contribution to the energy density
       zN1 = CST_hbc * i_eosDen%kf_n / CST_mnuc2
       fN1 = coefC * (CST_mnuc2)**4 * ( zN1*(1._pr+2._pr*zN1**2)*dsqrt(1._pr+zN1**2) - dlog(zN1+dsqrt(1._pr+zN1**2)) )
       zN2 = 0._pr
       fN2 = 0._pr
       if (i_eosDen%q1_kf_qd.gt.0._pr.and.i_withq) then
!          zb2 = CST_hbc * CST_Nc * i_eosDen%q1_kf_q / CST_mnuc2
          zN2 = CST_hbc * ( i_eosDen%kf_n - i_eosDen%q1_delta ) / CST_mnuc2
          fN2 = coefC * (CST_mnuc2)**4 * ( zN2*(1._pr+2._pr*zN2**2)*dsqrt(1._pr+zN2**2) - dlog(zN2+dsqrt(1._pr+zN2**2)) )
       endif
       fN = fN1 - fN2
       o_eosq0%rho_kin_N = fN
       ! quark contribution to the energy density
       fQd = 0._pr
       fQu = 0._pr
       if (i_eosDen%q1_kf_qd.gt.0._pr.and.i_withq) then
          zQd = CST_hbc * i_eosDen%q1_kf_qd / CST_mqc2
          fQd = CST_Nc * coefC * (CST_mqc2)**4 * ( zQd*(1._pr+2._pr*zQd**2)*dsqrt(1._pr+zQd**2) - dlog(zQd+dsqrt(1._pr+zQd**2)) )
          zQu = CST_hbc * i_eosDen%q1_kf_qu / CST_mqc2
          fQu = CST_Nc * coefC * (CST_mqc2)**4 * ( zQu*(1._pr+2._pr*zQu**2)*dsqrt(1._pr+zQu**2) - dlog(zQu+dsqrt(1._pr+zQu**2)) )
       endif
       fQ = fQu + fQd
       o_eosq0%rho_Q = fQ
       !write(*,*)"test",i_eosDen%den_nuc,i_eosDen%q1_kf_qd,fN1,fN2
       !    write(*,*)"test",i_eosDen%den_nuc,fq,fb1,fb2
       o_eosq0%rho_kin_B = o_eosq0%rho_kin_N + o_eosq0%rho_Q
       !
       ! Potential energy-density
       !
       xat = -28.8 ! MeV
       xbt = 10.0 ! MeV
       o_eosq0%rho_pot = xat * i_eosDen%den_n * (i_eosDen%den_n/0.16) + xbt * i_eosDen%den_n * (i_eosDen%den_n/0.16)**2
       !
       o_eosq0%rho_B = o_eosq0%rho_kin_B + o_eosq0%rho_pot
       !write(*,*)"test", o_eosq0%rho_kin, o_eosq0%rho_pot
       !
       ! energy-density and energy per particle (without the rest mass)
       !
       o_eosq0%e2v = o_eosq0%rho_B - i_eosDen%den_b * CST_mnuc2
       o_eosq0%e2a = o_eosq0%e2v / i_eosDen%den_b
       !
       ! chemical potential
       !
       ! Kinetic term
       C_2 = coefC*CST_hbc
       d1kFN_nN = i_eosDen%kf_n**2/CST_pi2
       d1z_fN1 = 8._pr*zN1**2*dsqrt(1._pr+zN1**2)
       d1kFN_rhokin = C_2*CST_mnuc2**3*d1z_fN1
       d1nB_kFN = 1._pr/d1kFN_nN
       ! noquark
       o_eosq0%mu_B_noq = d1kFN_rhokin * d1nB_kFN
       if (i_eosDen%q1_kf_qd.gt.0._pr.and.i_withq) then
          d1kFN_Delta = -2._pr*i_coef%LQCD**3 / CST_hbc**3 / i_eosDen%kf_n**3 
          d1kFN_nN =  (i_eosDen%kf_n**2 - (i_eosDen%kf_n - i_eosDen%q1_delta)**2*(1-d1kFN_Delta))/CST_pi2
          d1kFN_nQ = (1.5_pr*i_eosDen%q1_kf_qd**2*(1-d1kFN_Delta))/CST_pi2/CST_Nc
          d1z_fQu = 8._pr*zQu**2*dsqrt(1._pr+zQu**2)
          d1z_fQd = 8._pr*zQd**2*dsqrt(1._pr+zQd**2)
          d1z_fN2 = 8._pr*zN2**2*dsqrt(1._pr+zN2**2)
          d1kFN_rhokin = d1kFN_rhokin + C_2*(-CST_mnuc2**3*d1z_fN2 + CST_mqc2**3*(d1z_fQu/2._pr**CST_p13 + d1z_fQd))*(1-d1kFN_Delta)
          d1nB_kFN = 1._pr/(d1kFN_nN + d1kFN_nQ)
       endif
       o_eosq0%mu_kin_B = d1kFN_rhokin * d1nB_kFN
       ! Potential term
       d1kFN_pot = (2._pr*xat * (i_eosDen%den_n/0.16) + 3._pr*xbt * (i_eosDen%den_n/0.16)**2)*d1kFN_nN
       o_eosq0%mu_pot_B = d1kFN_pot * d1nB_kFN
       ! Total
       o_eosq0%mu_B = o_eosq0%mu_kin_B + o_eosq0%mu_pot_B
       ! noquark
!       o_eosq0%mu_B_noq = o_eosq0%mu_B_noq + d1kFN_pot * CST_pi2 / i_eosDen%kf_n**2
       !
       ! pressure
       !
       o_eosq0%p_kin_B =  i_eosDen%den_b * o_eosq0%mu_kin_B - o_eosq0%rho_kin_B
       o_eosq0%p_pot_B =  i_eosDen%den_b * o_eosq0%mu_pot_B - o_eosq0%rho_pot
       o_eosq0%p_B = o_eosq0%p_kin_B + o_eosq0%p_pot_B
!       write(*,*)"test",o_eosq0%e2v,o_eosq0%p_B
       !
       ! Incompressibility
       !
       ! Kinetic term
!       d2z_fN1 = 8._pr*zN1*(2._pr+3._pr*zN1**2)/dsqrt(1._pr+zN1**2)
!       d2kFN_rhoN = coefC*CST_hbc**2*CST_mnuc2**2 * d2z_fN1 
!       d2kFN_rhoQ = 0._pr
!       d2nB_kFN = -2._pr/CST_pi2*d1nB_kFN**3*i_eosDen%kf_n
!       d2kFN_nN = 2._pr/CST_pi2*i_eosDen%kf_n
!       if ((i_eosDen%kf_nuc - i_eosDen%q1_delta).gt.0._pr.and.i_withq) then
!          d2z_fN2 = 8._pr*zN2*(2._pr+3._pr*zN2**2)/dsqrt(1._pr+zN2**2)
!          d2z_fQd = 8._pr*zQd*(2._pr+3._pr*zQd**2)/dsqrt(1._pr+zQd**2)
!          d2z_fQu = 8._pr*zQu*(2._pr+3._pr*zQu**2)/dsqrt(1._pr+zQu**2)
!          d2kFN_Delta = 6._pr*i_coef%LQCD**3 / CST_hbc**3 / i_eosDen%kf_n**4
!          d2kFN_rhoN = d2kFN_rhoN &
!                   & - coefC*CST_hbc**2*CST_mnuc2**2 *(1._pr-d1kFN_Delta)**2* d2z_fN2 &
!                   & + coefC*CST_hbc*CST_mnuc2**3 *d2kFN_Delta*  d1z_fN2 
!          d2kFQ_rhoQ = CST_Nc*coefC*CST_hbc**2*CST_mqc2**2 * ( d2z_fQd + d2z_fQu/2**CST_p13 )
!          d2kFN_kFQ = -1._pr/CST_Nc*d2kFN_Delta
!          d2kFN_rhoQ = d2kFQ_rhoQ*d1kFN_kFQ**2 + d1kFQ_rhoQ*d2kFN_kFQ
!          d2kFN_nN = 2._pr/CST_pi2* ( 2._pr*i_eosDen%kf_nuc + d2kFN_Delta*(i_eosDen%kf_nuc - i_eosDen%q1_delta)**2 &
!                 & - 2._pr*(1._pr-d1kFN_Delta)**2*(i_eosDen%kf_nuc - i_eosDen%q1_delta)  )
!          d2nB_kFN = -2._pr/CST_pi2*d1nB_kFN**3* ( 2._pr*i_eosDen%kf_nuc + (1._pr - 1._pr/CST_Nc**3)*( d2kFN_Delta*(i_eosDen%kf_nuc - i_eosDen%q1_delta)**2 &
!                 & - 2._pr*(1._pr-d1kFN_Delta)**2*(i_eosDen%kf_nuc - i_eosDen%q1_delta) ) )         
!       endif
!       d2kFN_rhoB = d2kFN_rhoN + d2kFN_rhoQ
!       d2nB_rhoB = d2kFN_rhoB*d1nB_kFN**2 + (d1kFN_rhoN+d1kFN_rhoQ)*d2nB_kFN
!       o_eosq0%K_kin = 9._pr*i_eosDen%den_b*d2nB_rhoB
       d2kFN_nN = 2._pr * i_eosDen%kf_n / CST_pi2
       d2z_fN1 = 8._pr * zN1 * (2._pr+3._pr*zN1**2) / dsqrt(1._pr+zN1**2)
       d2kFN_rhokin = C_2 * CST_hbc * CST_mnuc2**2 * d2z_fN1
       d2nB_kFN = - d1nB_kFN**3 * d2kFN_nN
       if (i_eosDen%q1_kf_qd.gt.0._pr.and.i_withq) then
          d2kFN_Delta = 6._pr * i_coef%LQCD**3 / CST_hbc**3 / i_eosDen%kf_n**4
          d2kFN_nN = d2kFN_nN - 2._pr * ( i_eosDen%kf_n - i_eosDen%q1_delta) * ( 1._pr - d1kFN_Delta)**2 / CST_pi2 &
                 & + d2kFN_Delta / CST_pi2 * ( i_eosDen%kf_n - i_eosDen%q1_delta )**2
          d2kFN_nQ = 3._pr * i_eosDen%q1_kf_qd * ( 1._pr - d1kFN_Delta )**2 / CST_pi2 / CST_Nc**2 &
                 & - 1.5_pr * i_eosDen%q1_kf_qd**2 * d2kFN_Delta / CST_pi2 / CST_Nc
          d2z_fQu = 8._pr * zQu * (2._pr+3._pr*zQu**2) / dsqrt(1._pr+zQu**2)
          d2z_fQd = 8._pr * zQd * (2._pr+3._pr*zQd**2) / dsqrt(1._pr+zQd**2)
          d2z_fN2 = 8._pr * zN2 * (2._pr+3._pr*zN2**2) / dsqrt(1._pr+zN2**2)
          d2kFN_rhokin = d2kFN_rhokin + C_2 * &
               & ( - CST_hbc * CST_mnuc2**2 * d2z_fN2 * (1._pr - d1kFN_Delta)**2 &
               & + CST_mnuc2**3 * d1z_fN2 * d2kFN_Delta &
               & + CST_hbc * CST_mqc2**2 * ( d2z_fQu / 2._pr**CST_p23 + d2z_fQd ) * (1 - d1kFN_Delta)**2 / CST_Nc &
               & - CST_mqc2**3 * ( d1z_fQu / 2._pr**CST_p13 + d1z_fQd ) * d2kFN_Delta )
          d2nB_kFN = - d1nB_kFN**3 * ( d2kFN_nN + d2kFN_nQ )
       endif      
       !write(*,*)
       o_eosq0%K_kin = 9._pr * i_eosDen%den_b * ( d2kFN_rhokin * d1nB_kFN**2 + d1kFN_rhokin * d2nB_kFN )
       ! Potential term
       d2kFN_pot = ( 2._pr*xat /0.16 + 6._pr*xbt * i_eosDen%den_n/(0.16)**2 ) * d1kFN_nN**2 + &
            & d1kFN_pot / d1kFN_nN * d2kFN_nN
       o_eosq0%K_pot = 9._pr * i_eosDen%den_b * ( d2kFN_pot * d1nB_kFN**2 + d1kFN_pot * d2nB_kFN )
       ! Total
       o_eosq0%K_B = o_eosq0%K_kin + o_eosq0%K_pot
       !
       ! Enthalpy
       !
       o_eosq0%h_B = ( o_eosq0%rho_B + o_eosq0%p_B ) / i_eosDen%den_b 
       !
       ! sound velocity
       !
       o_eosq0%cs = o_eosq0%K_B / ( 9._pr * o_eosq0%h_B )
       !
    ! ---------------------
    ! Asymmetric matter
    ! ---------------------
    elseif (i_coef%iquark.ge.3) then
       !
       ! ---------------------
       ! Energies
       ! ---------------------
       !
       !
       ! neutron contribution to the energy density
       zN1 = CST_hbc * i_eosDen%kf_n / CST_mnuc2
       fN1 = coefC * (CST_mnuc2)**4 * ( zN1*(1._pr+2._pr*zN1**2)*dsqrt(1._pr+zN1**2) - dlog(zN1+dsqrt(1._pr+zN1**2)) )
       zN2 = 0._pr
       fN2 = 0._pr
       if ((i_eosDen%kf_nuc - i_eosDen%q1_delta).gt.0._pr.and.i_withq) then
!          zb2 = CST_hbc * CST_Nc * i_eosDen%q1_kf_q / CST_mnuc2
          zN2 = CST_hbc * ( i_eosDen%kf_n - i_eosDen%q1_delta*(1._pr + i_eosDen%xd)**CST_p13 ) / CST_mnuc2
          fN2 = coefC * (CST_mnuc2)**4 * ( zN2*(1._pr+2._pr*zN2**2)*dsqrt(1._pr+zN2**2) - dlog(zN2+dsqrt(1._pr+zN2**2)) )
       endif
       ! proton contribution to the energy density
       zP1 = CST_hbc * i_eosDen%kf_p / CST_mnuc2
       fP1 = coefC * (CST_mnuc2)**4 * ( zP1*(1._pr+2._pr*zP1**2)*dsqrt(1._pr+zP1**2) - dlog(zP1+dsqrt(1._pr+zP1**2)) )
       zP2 = 0._pr
       fP2 = 0._pr
       if ((i_eosDen%kf_nuc - i_eosDen%q1_delta).gt.0._pr.and.i_withq) then
          zP2 = CST_hbc * ( i_eosDen%kf_p - i_eosDen%q1_delta*(1._pr - i_eosDen%xd)**CST_p13 ) / CST_mnuc2
          fP2 = coefC * (CST_mnuc2)**4 * ( zP2*(1._pr+2._pr*zP2**2)*dsqrt(1._pr+zP2**2) - dlog(zP2+dsqrt(1._pr+zP2**2)) )
       endif
       fN = fN1 - fN2 + fP1 - fP2
       o_eosq0%rho_kin_N = fN
       o_eosq0%rho_kin_noq = fN1 + fP1
       ! quark contribution to the energy density
       fQd = 0._pr
       fQu = 0._pr
       if (i_eosDen%q1_kf_qd.gt.0._pr.and.i_withq) then
          zQd = CST_hbc * i_eosDen%q1_kf_qd / CST_mqc2
          fQd = CST_Nc * coefC * (CST_mqc2)**4 * ( zQd*(1._pr+2._pr*zQd**2)*dsqrt(1._pr+zQd**2) - dlog(zQd+dsqrt(1._pr+zQd**2)) )
          zQu = CST_hbc * i_eosDen%q1_kf_qu / CST_mqc2
          fQu = CST_Nc * coefC * (CST_mqc2)**4 * ( zQu*(1._pr+2._pr*zQu**2)*dsqrt(1._pr+zQu**2) - dlog(zQu+dsqrt(1._pr+zQu**2)) )
       endif
       fQ = fQu + fQd
       o_eosq0%rho_Q = fQ
!       write(*,*)"testQ",i_eosDen%q1_kf_qd,fQu,fQd,fQ
       !    write(*,*)"test",i_eosDen%den_nuc,fq,fb1,fb2
       o_eosq0%rho_kin_B = o_eosq0%rho_kin_N + o_eosq0%rho_Q
       !
       ! Potential energy-density
       !
!       o_eosq0%rho_pot = (i_eosDen%xd**2 * ( i_coef%VCOEF(1,0) + i_coef%VCOEF(1,1)* i_eosDen%xx &
!            & + 0.5 *i_coef%VCOEF(1,2)* i_eosDen%xx**2) + &
!            & ( i_coef%VCOEF(0,0) + i_coef%VCOEF(0,1) * i_eosDen%xx & 
!            & + 0.5 * i_coef%VCOEF(0,2) * i_eosDen%xx**2 )) * i_eosDen%den_nuc
!       write(*,*)"test rho_pot",i_eosDen%den_nuc,o_eosq0%rho_pot,i_eosb0%epot2v_b
       o_eosq0%rho_pot = i_eosb0%epot2v_b
       !
       !
       o_eosq0%rho_B = o_eosq0%rho_kin_B + o_eosq0%rho_pot
       o_eosq0%rho_N = o_eosq0%rho_kin_N + o_eosq0%rho_pot
       o_eosq0%rho_noq = o_eosq0%rho_kin_noq + o_eosq0%rho_pot
       !write(*,*)"test",i_eosDen%den_nuc,o_eosq0%rho_kin/i_eosDen%den_b,o_eosq0%rho_pot/i_eosDen%den_b,i_coef%VCOEF(0,0)
       !
       ! energy-density and energy per particle (without the rest mass)
       !
       o_eosq0%e2v = o_eosq0%rho_B - i_eosDen%den_b * CST_mnuc2
       o_eosq0%e2a = o_eosq0%e2v / i_eosDen%den_b
       ! no quarks
!       o_eosq0%e2v_noq = o_eosq0%rho_noq - i_eosDen%den_noq * CST_mnuc2
!       o_eosq0%e2a_noq = o_eosq0%e2v_noq / i_eosDen%den_noq

       !
       ! ---------------------
       ! chemical potentials
       ! ---------------------
       !
       ! TD kinetic part for nucleons 
       !
       d1z_fN1 = 8._pr*zN1**2*dsqrt(1._pr+zN1**2)
       d1z_fP1 = 8._pr*zP1**2*dsqrt(1._pr+zP1**2)
       d1kFN_rhoN = coefC*CST_hbc*CST_mnuc2**3 * ( (1._pr+i_eosDen%xd)**CST_p13 * d1z_fN1 + (1._pr-i_eosDen%xd)**CST_p13 * d1z_fP1 ) 
       d1nN_kFN = CST_pi2/2._pr / i_eosDen%kf_nuc**2
       if ( i_eosDen%xd.eq.1.) then
          d1deltaN_rhoN = coefC*CST_hbc*CST_mnuc2**3 * i_eosDen%kf_nuc/3._pr / (1._pr+i_eosDen%xd)**CST_p23 * d1z_fN1
       elseif ( i_eosDen%xd.eq.-1.) then
          d1deltaN_rhoN = -coefC*CST_hbc*CST_mnuc2**3 * i_eosDen%kf_nuc/3._pr / (1._pr-i_eosDen%xd)**CST_p23 * d1z_fP1 
       else
          d1deltaN_rhoN = coefC*CST_hbc*CST_mnuc2**3 * i_eosDen%kf_nuc/3._pr * ( d1z_fN1 / (1._pr+i_eosDen%xd)**CST_p23 &
                      & - d1z_fP1 / (1._pr-i_eosDen%xd)**CST_p23 )
       endif
       if ((i_eosDen%kf_nuc - i_eosDen%q1_delta).gt.0._pr.and.i_withq) then
          d1kFN_Delta = -2._pr*i_coef%LQCD**3 / CST_hbc**3 / i_eosDen%kf_nuc**3 
          d1z_fN2 = 8._pr*zN2**2*dsqrt(1._pr+zN2**2)  
          d1z_fP2 = 8._pr*zP2**2*dsqrt(1._pr+zP2**2) 
          d1kFN_rhoN = d1kFN_rhoN - coefC*CST_hbc*CST_mnuc2**3 * (1._pr-d1kFN_Delta) * ( (1._pr+i_eosDen%xd)**CST_p13 * d1z_fN2 &
                   & + (1._pr-i_eosDen%xd)**CST_p13 * d1z_fP2 )
          d1nN_kFN = CST_pi2/2._pr / (i_eosDen%kf_nuc**2 - (1._pr-d1kFN_Delta) * (i_eosDen%kf_nuc - i_eosDen%q1_delta)**2 )
          if ( i_eosDen%xd.eq.1.) then
             d1deltaN_rhoN = d1deltaN_rhoN &
               & - coefC*CST_hbc*CST_mnuc2**3 * (i_eosDen%kf_nuc-i_eosDen%q1_delta)/3._pr / (1._pr+i_eosDen%xd)**CST_p23 * d1z_fN1
          elseif ( i_eosDen%xd.eq.-1.) then
             d1deltaN_rhoN = d1deltaN_rhoN &
               & + coefC*CST_hbc*CST_mnuc2**3 * (i_eosDen%kf_nuc-i_eosDen%q1_delta)/3._pr / (1._pr-i_eosDen%xd)**CST_p23 * d1z_fP1 
          else
             d1deltaN_rhoN = d1deltaN_rhoN - coefC*CST_hbc*CST_mnuc2**3 * (i_eosDen%kf_nuc - i_eosDen%q1_delta)/3._pr &
                         & * ( d1z_fN1 / (1._pr+i_eosDen%xd)**CST_p23 - d1z_fP1 / (1._pr-i_eosDen%xd)**CST_p23 )
          endif               
       endif
       d1nN_rhoN = d1kFN_rhoN*d1nN_kFN
       !
       o_eosq0%mun_kin = d1nN_rhoN + (1._pr-i_eosDen%xd)/i_eosDen%den_nuc * d1deltaN_rhoN
       o_eosq0%mup_kin = d1nN_rhoN - (1._pr+i_eosDen%xd)/i_eosDen%den_nuc * d1deltaN_rhoN
       o_eosq0%munp_kin = 2._pr * d1deltaN_rhoN / i_eosDen%den_nuc
       !
       ! TD potential part for nucleons
       !
!       d1nN_pot = o_eosq0%rho_pot/i_eosDen%den_nuc &
!              & + i_eosDen%den_nuc/3._pr/CST_nsat * ( i_eosDen%xd**2 * ( i_coef%VCOEF(1,1) + i_coef%VCOEF(1,2)*i_eosDen%xx ) &
!              & + i_coef%VCOEF(0,1) + i_coef%VCOEF(0,2)*i_eosDen%xx )
!       d1deltaN_pot = 2._pr*i_eosDen%xd*i_eosDen%den_nuc * ( i_coef%VCOEF(1,0) + i_coef%VCOEF(1,1)*i_eosDen%xx &
!            & + 0.5 *i_coef%VCOEF(1,2)* i_eosDen%xx**2 )
       !
!       write(*,*)"test mu",i_eosDen%den_nuc,d1nN_pot,i_eosb0%mu_pot_is,d1deltaN_pot,0.5_pr * i_eosDen%den_nuc * i_eosb0%mu_pot_np
       d1nN_pot = i_eosb0%mu_pot_is
       d1deltaN_pot = 0.5_pr * i_eosDen%den_nuc * i_eosb0%mu_pot_np
       !
       o_eosq0%mun_pot = d1nN_pot + (1._pr-i_eosDen%xd)/i_eosDen%den_nuc * d1deltaN_pot
       o_eosq0%mup_pot = d1nN_pot - (1._pr+i_eosDen%xd)/i_eosDen%den_nuc * d1deltaN_pot

       o_eosq0%mun_pot = i_eosb0%mu_pot_n
       o_eosq0%mup_pot = i_eosb0%mu_pot_p
       o_eosq0%munp_pot = i_eosb0%mu_pot_np

       !
       ! TD Total for nucleons
       !
       o_eosq0%mun = o_eosq0%mun_kin + o_eosq0%mun_pot
       o_eosq0%mup = o_eosq0%mup_kin + o_eosq0%mup_pot 
       o_eosq0%munp = o_eosq0%munp_kin + o_eosq0%munp_pot 
       !write(*,*)"test",d1deltaN_rhoN,d1nN_rhoN,o_eosq0%mun,o_eosq0%mup,d1nN_pot,d1deltaN_pot
       !
       o_eosq0%mun_qp = dsqrt( CST_mnuc2**2 + (CST_hbc*i_eosDen%kf_n)**2) + o_eosq0%mun_pot
       o_eosq0%mup_qp = dsqrt( CST_mnuc2**2 + (CST_hbc*i_eosDen%kf_p)**2) + o_eosq0%mup_pot
       !
       ! chemical potential for quarks
       !
       d1kFQ_rhoQ = 0._pr 
       d1nQ_kFQ = 0._pr
       d1deltaQ_rhoQ = 0._pr
       if ((i_eosDen%kf_nuc - i_eosDen%q1_delta).gt.0._pr.and.i_withq) then
          d1z_fQu = 8._pr*zQu**2*dsqrt(1._pr+zQu**2)
          d1z_fQd = 8._pr*zQd**2*dsqrt(1._pr+zQd**2)
          d1kFQ_rhoQ = CST_Nc*coefC*CST_hbc*CST_mqc2**3 * ( (1._pr+i_eosDen%xd/CST_Nc)**CST_p13 * d1z_fQd &
                   & + (1._pr-i_eosDen%xd/CST_Nc)**CST_p13 * d1z_fQu )
          d1nQ_kFQ = CST_pi2/2._pr / i_eosDen%q1_kf_q**2
          d1deltaQ_rhoQ = coefC*CST_hbc*CST_mqc2**3 * i_eosDen%q1_kf_q/3._pr * ( d1z_fQd / (1._pr+i_eosDen%xd/CST_Nc)**CST_p23 &
                      & - d1z_fQu / (1._pr-i_eosDen%xd/CST_Nc)**CST_p23 )
       endif
       d1nQ_rhoQ = d1kFQ_rhoQ * d1nQ_kFQ
       !
       o_eosq0%mud = 0._pr
       o_eosq0%muu = 0._pr
       o_eosq0%mud_qp = 0._pr
       o_eosq0%muu_qp = 0._pr
       if ((i_eosDen%kf_nuc - i_eosDen%q1_delta).gt.0._pr.and.i_withq) then
          o_eosq0%mud = d1nQ_rhoQ + (1._pr-i_eosDen%xd/CST_Nc)/i_eosDen%q1_nq * d1deltaQ_rhoQ
          o_eosq0%muu = d1nQ_rhoQ - (1._pr+i_eosDen%xd/CST_Nc)/i_eosDen%q1_nq * d1deltaQ_rhoQ
       o_eosq0%muu_qp = CST_Nc * dsqrt( CST_mqc2**2 + (CST_hbc*i_eosDen%q1_kf_qu)**2)
       o_eosq0%mud_qp = CST_Nc * dsqrt( CST_mqc2**2 + (CST_hbc*i_eosDen%q1_kf_qd)**2)
       endif         
       !
       ! TD Total baryonic chemical potential (as in SM and NM) 
       ! TD kinetic for nucleon part
       !
       d1z_fN1 = 8._pr*zN1**2*dsqrt(1._pr+zN1**2)
       d1z_fP1 = 8._pr*zP1**2*dsqrt(1._pr+zP1**2)
       d1kFN_rhoN = coefC*CST_hbc*CST_mnuc2**3 * ( (1+i_eosDen%xd)**CST_p13 * d1z_fN1 + (1._pr-i_eosDen%xd)**CST_p13 * d1z_fP1 ) 
       d1nN_kFN = CST_pi2/2._pr / i_eosDen%kf_nuc**2
       d1nB_kFN = d1nN_kFN
       if ((i_eosDen%kf_nuc - i_eosDen%q1_delta).gt.0._pr.and.i_withq) then
          d1kFN_Delta = -2._pr * i_coef%LQCD**3 / CST_hbc**3 / i_eosDen%kf_nuc**3 
          d1z_fN2 = 8._pr*zN2**2*dsqrt(1._pr+zN2**2)  
          d1z_fP2 = 8._pr*zP2**2*dsqrt(1._pr+zP2**2) 
          d1kFN_rhoN = d1kFN_rhoN - coefC*CST_hbc*CST_mnuc2**3 * (1._pr-d1kFN_Delta) * ( (1._pr+i_eosDen%xd)**CST_p13 * d1z_fN2 &
                   & + (1._pr-i_eosDen%xd)**CST_p13 * d1z_fP2 )
          d1nN_kFN = CST_pi2/2._pr / (i_eosDen%kf_nuc**2 &
                 & - (1._pr-d1kFN_Delta) * (i_eosDen%kf_nuc - i_eosDen%q1_delta)**2 )
!          xtest = 1._pr / d1nN_kFN - 2/CST_pi2 * i_eosDen%kf_nuc**2
!          d1kFN_nN1 = 2._pr/CST_pi2 * i_eosDen%kf_nuc**2
!          d1kFN_nN2 = 2._pr/CST_pi2 * (i_eosDen%kf_nuc-i_eosDen%q1_delta)**2 * (1._pr-d1kFN_Delta)
!          d1kFN_nN = d1kFN_nN1 - d1kFN_nN2
          d1nN_kFN = CST_pi2/2._pr / (i_eosDen%kf_nuc**2 &
                 & - (1._pr-d1kFN_Delta) * (i_eosDen%kf_nuc - i_eosDen%q1_delta)**2 )
          d1nB_kFN = CST_pi2/2._pr / (i_eosDen%kf_nuc**2 &
                 & - (1._pr - 1._pr/CST_Nc**3) * (1._pr-d1kFN_Delta) * (i_eosDen%kf_nuc - i_eosDen%q1_delta)**2 )               
!          write(*,*)"analyt.:",i_eosDen%den_B,1.0/d1kFN_nN1,1.0/d1kFN_nN2,1.0/(d1kFN_nN1-d1kFN_nN2),d1nN_kFN
       endif
       d1nB_rhoN = d1nB_kFN*d1kFN_rhoN
       !
!       o_eosq0%mu_kin_N = d1nB_rhoN
       !
       ! TD potential for nucleon part
       !
!       d1nN_pot = o_eosq0%rho_pot/i_eosDen%den_nuc &
!              & + 1._pr/3._pr/CST_nsat * ( i_eosDen%xd**2 * ( i_coef%VCOEF(1,1) + i_coef%VCOEF(1,2)*i_eosDen%xx ) &
!              & + i_coef%VCOEF(0,1) + i_coef%VCOEF(0,2)*i_eosDen%xx )
       d1nB_pot = d1nN_pot / d1nN_kFN * d1nB_kFN
       !
       o_eosq0%mu_pot_B = d1nB_pot
       !
       !o_eosq0%mu_N = o_eosq0%mu_kin_N + o_eosq0%mu_pot_N
       !
       !o_eosq0%mu_kin_N_qp = dsqrt( CST_mnuc2**2 + (CST_hbc*i_eosDen%kf_nuc)**2 )
       !o_eosq0%mu_N_qp = o_eosq0%mu_kin_N_qp + o_eosq0%mu_pot_N
       !
       ! quark part
       !
       d1kFN_rhoQ = 0._pr
       if ((i_eosDen%kf_nuc - i_eosDen%q1_delta).gt.0._pr.and.i_withq) then
          d1z_fQd = 8._pr*zQd**2*dsqrt(1._pr+zQd**2)
          d1z_fQu = 8._pr*zQu**2*dsqrt(1._pr+zQu**2) 
          d1kFQ_rhoQ = CST_Nc*coefC*CST_hbc*CST_mqc2**3 * ( (1+i_eosDen%xd/CST_Nc)**CST_p13 * d1z_fQd &
                   & + (1._pr-i_eosDen%xd/CST_Nc)**CST_p13 * d1z_fQu )   
          d1kFN_kFQ = (1._pr-d1kFN_Delta)/CST_Nc
          d1kFN_rhoQ = d1kFQ_rhoQ*d1kFN_kFQ
       endif
       d1nB_rhoQ = d1nB_kFN*d1kFN_rhoQ
       !
!       o_eosq0%mu_Q = d1nB_rhoQ
       !
       ! Total 
       !
       o_eosq0%mu_kin_B = d1nB_rhoN + d1nB_rhoQ
       !       o_eosq0%mu_kin_B = o_eosq0%mu_kin_N + o_eosq0%mu_Q
!       o_eosq0%mu_pot_B = o_eosq0%mu_pot_N
       o_eosq0%mu_B = o_eosq0%mu_kin_B + o_eosq0%mu_pot_B
       !
       ! ---------------------
       ! Pressure
       ! ---------------------
       !
       ! partial pressures
       !
       o_eosq0%pn_kin = i_eosDen%den_n * o_eosq0%mun_kin + i_eosDen%den_p * o_eosq0%mup_kin - o_eosq0%rho_kin_N
       o_eosq0%pn_pot = i_eosDen%den_n * o_eosq0%mun_pot + i_eosDen%den_p * o_eosq0%mup_pot - o_eosq0%rho_pot
       o_eosq0%pn = o_eosq0%pn_kin + o_eosq0%pn_pot
       !write(*,*)"test mu",o_eosq0%mu_N,o_eosq0%mun,o_eosq0%mup
!       o_eosq0%pn =  i_eosDen%den_nuc * o_eosq0%mu_N - o_eosq0%rho_N
       o_eosq0%pn_qp =  i_eosDen%den_n * o_eosq0%mun_qp + i_eosDen%den_p * o_eosq0%mup_qp - o_eosq0%rho_N
!       write(*,*)"test pression",o_eosq0%pn_qp,i_eosDen%den_nuc * o_eosq0%mu_N - o_eosq0%rho_N
       o_eosq0%pq = i_eosDen%q1_nu * o_eosq0%muu + i_eosDen%q1_nd * o_eosq0%mud - o_eosq0%rho_Q
       !
       ! total pressure
       !
       o_eosq0%p_kin_B = i_eosDen%den_B * o_eosq0%mu_kin_B - o_eosq0%rho_kin_B
       o_eosq0%p_pot_B = i_eosDen%den_B * o_eosq0%mu_pot_B - o_eosq0%rho_pot
       o_eosq0%p_B = o_eosq0%p_kin_B + o_eosq0%p_pot_B
!       o_eosq0%p_B_qp =  i_eosDen%den_B * ( o_eosq0%mu_N_qp + o_eosq0%mu_Q ) - o_eosq0%rho_B
!       write(*,*)"test"!,o_eosq0%mu_B,o_eosq0%mun,o_eosq0%mup,o_eosq0%muu,o_eosq0%mud
!       write(*,*)"test"," P_N = ",o_eosq0%pn,"P_Q = ",o_eosq0%pq,"P_N + P_Q = ",o_eosq0%pn+o_eosq0%pq,"P_B = ",o_eosq0%p_B
!       write(*,*)"test"!,o_eosq0%pn_pot,o_eosq0%p_pot,o_eosq0%mun_pot,o_eosq0%mup_pot,o_eosq0%mu_pot,o_eosq0%rho_pot
       !write(*,*)"test",i_eosDen%den_n * o_eosq0%mun_pot + i_eosDen%den_p * o_eosq0%mup_pot,i_eosDen%den_B * o_eosq0%mu_pot
       !---------------------------------------------------------------------------------------------------------------------------------------------
       ! incompressibility
       !
	   !
       d2z_fN1 = 8._pr*zN1*(2._pr+3._pr*zN1**2)/dsqrt(1._pr+zN1**2)
       d2z_fP1 = 8._pr*zP1*(2._pr+3._pr*zP1**2)/dsqrt(1._pr+zP1**2)
       d2kFN_rhoN = coefC*CST_hbc**2*CST_mnuc2**2 * ( (1+i_eosDen%xd)**CST_p23 * d2z_fN1 + (1._pr-i_eosDen%xd)**CST_p23 * d2z_fP1 ) 
       d2kFN_rhoQ = 0._pr
       d2nB_kFN = -4._pr/CST_pi2*d1nB_kFN**3*i_eosDen%kf_nuc
       d2kFN_nN = 4._pr/CST_pi2*i_eosDen%kf_nuc
       if ((i_eosDen%kf_nuc - i_eosDen%q1_delta).gt.0._pr.and.i_withq) then
          d2z_fN2 = 8._pr*zN2*(2._pr+3._pr*zN2**2)/dsqrt(1._pr+zN2**2)
          d2z_fP2 = 8._pr*zP2*(2._pr+3._pr*zP2**2)/dsqrt(1._pr+zP2**2)
          d2z_fQd = 8._pr*zQd*(2._pr+3._pr*zQd**2)/dsqrt(1._pr+zQd**2)
          d2z_fQu = 8._pr*zQu*(2._pr+3._pr*zQu**2)/dsqrt(1._pr+zQu**2)
          d2kFN_Delta = 6._pr*i_coef%LQCD**3 / CST_hbc**3 / i_eosDen%kf_nuc**4
          d2kFN_rhoN = d2kFN_rhoN &
                   & - coefC*CST_hbc**2*CST_mnuc2**2 *(1._pr-d1kFN_Delta)**2* ( (1+i_eosDen%xd)**CST_p23 * d2z_fN2 + (1._pr-i_eosDen%xd)**CST_p23 * d2z_fP2 ) &
                   & + coefC*CST_hbc*CST_mnuc2**3 *d2kFN_Delta* ( (1+i_eosDen%xd)**CST_p13 * d1z_fN2 + (1._pr-i_eosDen%xd)**CST_p13 * d1z_fP2 )
          d2kFQ_rhoQ = CST_Nc*coefC*CST_hbc**2*CST_mqc2**2 * ( (1+i_eosDen%xd/CST_Nc)**CST_p23 * d2z_fQd + (1._pr-i_eosDen%xd/CST_Nc)**CST_p23 * d2z_fQu )
          d2kFN_kFQ = -1._pr/CST_Nc*d2kFN_Delta
          d2kFN_rhoQ = d2kFQ_rhoQ*d1kFN_kFQ**2 + d1kFQ_rhoQ*d2kFN_kFQ
          d2kFN_nN = 2._pr/CST_pi2* ( 2._pr*i_eosDen%kf_nuc + d2kFN_Delta*(i_eosDen%kf_nuc - i_eosDen%q1_delta)**2 &
                 & - 2._pr*(1._pr-d1kFN_Delta)**2*(i_eosDen%kf_nuc - i_eosDen%q1_delta)  )
          d2nB_kFN = -2._pr/CST_pi2*d1nB_kFN**3* ( 2._pr*i_eosDen%kf_nuc + (1._pr - 1._pr/CST_Nc**3)*( d2kFN_Delta*(i_eosDen%kf_nuc - i_eosDen%q1_delta)**2 &
                 & - 2._pr*(1._pr-d1kFN_Delta)**2*(i_eosDen%kf_nuc - i_eosDen%q1_delta) ) )         
       endif
       d2kFN_rhoB = d2kFN_rhoN + d2kFN_rhoQ
       d2nB_rhoB = d2kFN_rhoB*d1nB_kFN**2 + (d1kFN_rhoN+d1kFN_rhoQ)*d2nB_kFN
       o_eosq0%k_kin = 9._pr*i_eosDen%den_b*d2nB_rhoB
       !
       d2kFN_pot = i_eosb0%K_pot_b/9._pr/i_eosDen%den_nuc / d1nN_kFN**2 + d1nN_pot * d2kFN_nN
       o_eosq0%k_pot = 9._pr*i_eosDen%den_b* ( d2kFN_pot * d1nB_kFN**2 + d1nN_pot/d1nN_kFN * d2nB_kFN )
       o_eosq0%K_B = o_eosq0%k_kin + o_eosq0%k_pot
       !
       ! Enthalpy
       !
       o_eosq0%h_B = ( o_eosq0%rho_B + o_eosq0%p_B ) / i_eosDen%den_b 
       !
       ! sound velocity
       !
       o_eosq0%cs = o_eosq0%K_B / ( 9._pr * o_eosq0%h_B )
!       write(*,*)"test",o_eosq0%K_pot,o_eosq0%K_kin,o_eosq0%cs
      ! write(*,*)"test",coefC*CST_hbc**2*CST_mnuc2**2 * ( d2z_fN1 + d2z_fP1 ) 
      ! write(*,*)"test",-coefC*CST_hbc**2*CST_mnuc2**2 *(1._pr-d1kFN_Delta)**2* ( d2z_fN2 + d2z_fP2 )
      ! write(*,*)"test",coefC*CST_hbc*CST_mnuc2**3 *d2kFN_Delta* ( (1+i_eosDen%xd)**CST_p13 * d1z_fN2 + (1._pr-i_eosDen%xd)**CST_p13 * d1z_fP2 ) 
       !---------------------------------------------------------------------------------------------------------------------------------------------
    endif
    !
  end subroutine metaeos_T0_q1

  
  subroutine metaeos_T0_baryons_NR(i_coef,i_eosDen,i_opt_beta,o_eosb0)

!    use acc; use cst;
!    use metaEosType;
    
    implicit none

    ! input variables
    type (metaEos_coef), intent(in)       :: i_coef
    type (metaEos_densities), intent(in)  :: i_eosDen
    logical, intent(in)                   :: i_opt_beta
    ! output variables
    type (metaEos_Baryons), intent(out)   :: o_eosb0
    !
    ! local variables
    !
!    logical                   :: i_opt_beta2

    integer                          :: n
!    real (kind=pr), dimension(0:5)  :: fac
    real (kind=pr)                   :: tmp, sum
    ! kinetic energy factors
    real (kind=pr)                   :: ffg, d1_ffg, d2_ffg
    real (kind=pr)                   :: fsky, d1_fsky, d2_fsky
    real (kind=pr)                   :: ffg_sym, fsky_sym
    ! first order derivatives
    real (kind=pr)                   :: d1_t2aFG_b_den
    real (kind=pr)                   :: d1_t2a_b_den, d1_t2a_b_xe
    ! second order derivatives
    real (kind=pr)                   :: d2_t2aFG_b_den
    real (kind=pr)                   :: d2_t2a_b_den, d2_t2a_b_den_xe, d2_t2a_b_xe
    real (kind=pr), dimension(0:4)   :: d2_e2a_b_den, d2_e2a_b_den_xe, d2_e2a_b_xe
    ! third order derivatives
    real (kind=pr)                   :: d3_t2a_b_den
    real (kind=pr)                   :: chi, kappa, enthalpy
    real (kind=pr)                   :: coefn, coefp
    ! parameters
    real (kind=pr), parameter        :: p13=1._pr/3._pr,  p23=2./3.,  p43=4./3.,  p53=5./3.
    real (kind=pr), parameter        :: m13=-1./3., m23=-2./3., m43=-4./3., m53=-5./3., m73 = -7./3.
    ! vcoef: interaction coefficient
    real (kind=pr), dimension(0:4)   :: vcoef
    ! 1+3x
    real (kind=pr) :: up3x, up3xp23
   
    ! *******************
    ! *******************
    ! Baryons at T = 0
    ! *******************
    ! *******************
!    write(66,*)"Baryons"
!    i_opt_beta2=.true.

   
    ! For Fermi gas Energy
    if (.not.i_opt_beta) then
       ffg   = ( 1.0 + i_eosDen%xd )**p53 + ( 1.0 - i_eosDen%xd )**p53
    endif
    d1_ffg  = 5.0/3.0 * ( ( 1.0 + i_eosDen%xd )**p23 - ( 1.0 - i_eosDen%xd )**p23 )
    if (.not.i_opt_beta) then
       d2_ffg  = 10.0/9.0 * ( ( 1.0 + i_eosDen%xd )**m13 + (1.0 - i_eosDen%xd )**m13 )

       ! For kinetic energy with effective mass effect
       fsky  = ( i_coef%mb + i_coef%db * i_eosDen%xd ) * ( 1.0 + i_eosDen%xd )**p53 + &
            & ( i_coef%mb - i_coef%db * i_eosDen%xd ) * ( 1.0 - i_eosDen%xd )**p53
    endif
    d1_fsky = 5.0/3.0 * ( ( i_coef%mb + i_coef%db * i_eosDen%xd ) * ( 1.0 + i_eosDen%xd )**p23 - &
         ( i_coef%mb - i_coef%db * i_eosDen%xd ) * ( 1.0 - i_eosDen%xd )**p23 ) + & 
         i_coef%db * ( ( 1.0 + i_eosDen%xd )**p53 - (1.0 - i_eosDen%xd )**p53 )
    if (.not.i_opt_beta) then
       d2_fsky = 10.0/9.0 * ( ( i_coef%mb + i_coef%db * i_eosDen%xd ) * ( 1.0 + i_eosDen%xd )**m13 + &
            ( i_coef%mb - i_coef%db * i_eosDen%xd ) * ( 1.0 - i_eosDen%xd )**m13 ) + & 
            10.0/3.0 * i_coef%db * ( ( 1.0 + i_eosDen%xd )**p23 + ( 1.0 - i_eosDen%xd )**p23 )
       ffg_sym = 20/9._pr
       fsky_sym = 20/9._pr * ( i_coef%mb + 3 * i_coef%db )
       !
       ! kinetic densities
       if (i_coef%ikin.eq.0) then
          o_eosb0%tau_n   = 0._pr
          o_eosb0%tau_p   = 0._pr
       elseif (i_coef%ikin.eq.1) then
          o_eosb0%tau_n   = 3._pr / 5._pr *  i_eosDen%den_n * i_eosDen%kf_n**2
          o_eosb0%tau_p   = 3._pr / 5._pr *  i_eosDen%den_p * i_eosDen%kf_p**2
       else
          stop "i_coef%ikin ill-defined in metaeos_T0_baryons()"
       endif
       o_eosb0%tau_b   = o_eosb0%tau_n + o_eosb0%tau_p
    endif
    !
    ! -------------------
    ! Density correction u_alpha(x,delta)
    ! -------------------
    !
    do n = 0, i_coef%nmax
       o_eosb0%udc(n) = ( 1._pr - ( - 3._pr * i_eosDen%xx )**(i_coef%nmax+1-n) * &
            & dexp( - ( i_coef%bsat + i_coef%bsym * i_eosDen%xd**2 ) * i_eosDen%den_nuc / i_coef%nsat ) ) 
       o_eosb0%udcsat(n) = ( 1._pr - ( - 3._pr * i_eosDen%xx )**(i_coef%nmax+1-n) * &
            & dexp( - i_coef%bsat * i_eosDen%den_nuc / i_coef%nsat ) ) 
       o_eosb0%udcNM(n) = ( 1._pr - ( - 3._pr * i_eosDen%xx )**(i_coef%nmax+1-n) * &
            & dexp( - (i_coef%bsat+i_coef%bsym) * i_eosDen%den_nuc / i_coef%nsat ) ) 
       ! first derivative:
       o_eosb0%d1_udc(n) = 0.0_pr
       o_eosb0%d2_udc(n) = 0.0_pr
       if (i_eosDen%xx.ne.0.0_pr) then
          o_eosb0%d1_udc(n) = ( i_coef%nmax + 1 - n - &
               & 3 * ( i_coef%bsat + i_coef%bsym * i_eosDen%xd**2 ) * i_eosDen%xx ) * &
               & ( o_eosb0%udc(n) - 1.0_pr ) / i_eosDen%xx
          o_eosb0%d2_udc(n) = ( - ( i_coef%nmax + 1 - n ) * ( i_coef%nmax - n ) + &
               & 6 * ( i_coef%bsat + i_coef%bsym * i_eosDen%xd**2 ) * i_eosDen%xx * ( i_coef%nmax + 1 - n ) - &
               & 9 * ( i_coef%bsat + i_coef%bsym * i_eosDen%xd**2 )**2 * i_eosDen%xx**2 ) * &
               & ( 1.0_pr - o_eosb0%udc(n) ) / i_eosDen%xx**2
       endif
    enddo
    !
    ! -------------------
    ! full interaction parameter v(x,delta)
    ! -------------------
    !
    do n = 0, i_coef%nmax
       vcoef(n) = i_coef%VCOEF(0,n) + i_coef%VCOEF(1,n) * i_eosDen%xd**2 + i_coef%VCOEF(2,n) * i_eosDen%xd**4
    enddo
    !
    ! -------------------
    ! 1+3x, (1+3x)^2/3
    ! -------------------
    up3x = ( 1.0_pr + 3 * i_eosDen%xx )
    up3xp23 = ( 1.0_pr + 3 * i_eosDen%xx )**CST_p23
!       up3x = 1.0_pr
    !       up3xp23 = 1.0_pr
    if (.not.i_opt_beta) then
       !
       ! -------------------
       ! energy per particle
       ! -------------------
       !
       ! kinetic energy : Fermi Gas
       o_eosb0%t2aFG_b = 0.5_pr * i_coef%TFGsat * up3xp23 * ffg
       ! kinetic energy : total
       o_eosb0%t2a_b = 0.5_pr * i_coef%TFGsat * up3xp23 * ( ffg + ( 1.0 + 3 * i_eosDen%xx ) * fsky )
       ! energy per unit volum
       o_eosb0%t2v_b = 0.5_pr * CST_htm * ( o_eosb0%tau_b + i_coef%mb / i_coef%nsat * i_eosDen%den_nuc * o_eosb0%tau_b + &
            & i_coef%db / i_coef%nsat *  i_eosDen%xd  * i_eosDen%den_nuc * ( o_eosb0%tau_n - o_eosb0%tau_p ) ) 
       ! potential energy and total energy
       do n = 0, i_coef%nmax
          o_eosb0%a_v2a_order_n(n) = vcoef(n) * (i_eosDen%xx**n) / i_coef%facto(n) * o_eosb0%udc(n)
          if (n.eq.0) then
             o_eosb0%a_v2a(n) = o_eosb0%a_v2a_order_n(n)
          else
             o_eosb0%a_v2a(n) = o_eosb0%a_v2a(n-1) + o_eosb0%a_v2a_order_n(n)
          endif
          o_eosb0%a_e2a_b(n) = o_eosb0%t2a_b + o_eosb0%a_v2a(n)
          o_eosb0%a_e2v_b(n) = o_eosb0%t2v_b + i_eosDen%den_nuc * o_eosb0%a_v2a(n)
       enddo
       o_eosb0%epot2a_b = o_eosb0%a_v2a(i_coef%nmax)
       o_eosb0%epot2v_b = i_eosDen%den_nuc * o_eosb0%a_v2a(i_coef%nmax)
       o_eosb0%e2a_b = o_eosb0%a_e2a_b(i_coef%nmax)
       o_eosb0%e2v_b = o_eosb0%e2a_b * i_eosDen%den_nuc
!       o_eosb0%eps_b = CST_mnuc2 * i_eosDen%den_nuc + o_eosb0%e2v_b
       o_eosb0%rho_b = i_eosDen%den_p * CST_mpc2 + i_eosDen%den_n * CST_mnc2 + o_eosb0%e2v_b
       !
       ! -------------------
       ! symmetry energy per particle
       ! -------------------
       !
       ! FG contribution
       o_eosb0%tsymFG_b = i_coef%TFGsat * up3xp23 * ( 2**CST_p23 - 1.0_pr )
!       o_eosb0%tsym2FG_b = 0.25_pr * i_coef%TFGsat * ( 1.0 + 3 * i_eosDen%xx )**p23 * ffg_sym
       ! Kinetic contribution
!       o_eosb0%tsym_b = i_coef%TFGsat * ( 1.0_pr + i_coef%kappaNM )!* ( 2**CST_p23 * ( 1.0_pr + up3x*i_coef%kappaNM ) )!&
!       o_eosb0%tsym_b = i_coef%TFGsat * up3xp23 * ( 1.0_pr + up3x*i_coef%kappaNM )!* ( 2**CST_p23 * ( 1.0_pr + up3x*i_coef%kappaNM ) )!&
!         &            - ( 1.0_pr + up3x*i_coef%kappas ) )
       o_eosb0%tsym_b = i_coef%TFGsat * up3xp23 * ( 2**CST_p23 * ( 1.0_pr + up3x*i_coef%kappaNM ) - &
         &            ( 1.0_pr + up3x*i_coef%kappas ) )
!       write(*,*)'check2:',i_coef%kappas,i_coef%kappaNM,i_coef%TFGsat,2**CST_p23,up3xp23 
!       write(*,*)'check2:',i_eosDen%den_nuc,o_eosb0%tsym_b
!       write(*,*)'check2:', o_eosb0%udcsat(0), o_eosb0%udcNM(0)
       !
       ! -------------------
       ! symmetry energy S_2 per particle
       ! -------------------
       !
       ! FG contribution
       o_eosb0%tsym2FG_b = 5.0_pr/9.0_pr * i_coef%TFGsat * up3xp23 
!       o_eosb0%tsym2FG_b = 0.25_pr * i_coef%TFGsat * ( 1.0 + 3 * i_eosDen%xx )**p23 * ffg_sym
       ! Kinetic contribution
       o_eosb0%tsym2_b = 5.0_pr/9.0_pr * i_coef%TFGsat * up3xp23 * ( 1.0_pr &
            & + up3x * (i_coef%mb + 3 * i_coef%db) )
!       o_eosb0%tsym2_b = 0.25_pr * i_coef%TFGsat * ( 1.0 + 3 * i_eosDen%xx )**p23 * &
!            & ( ffg_sym + ( 1.0 + 3 * i_eosDen%xx ) * fsky_sym )
    endif
    ! potential contribution
    do n = 0, i_coef%nmax
       o_eosb0%a_vsym_order_n(n) = (i_eosDen%xx**n) / i_coef%facto(n) * ( &
            ( i_coef%VCOEF(0,n) + i_coef%VCOEF(1,n) + i_coef%VCOEF(2,n) ) * o_eosb0%udcNM(n) - &
            ( i_coef%VCOEF(0,n) ) * o_eosb0%udcsat(n) )
    enddo
    o_eosb0%a_vsym(0) = o_eosb0%a_vsym_order_n(0)
    do n = 1, i_coef%nmax
       o_eosb0%a_vsym(n) = o_eosb0%a_vsym(n-1) + o_eosb0%a_vsym_order_n(n)
    enddo
    ! potential contribution S_2
    do n = 0, i_coef%nmax
       o_eosb0%a_vsym2_order_n(n) = (i_eosDen%xx**n) / i_coef%facto(n) * ( i_coef%VCOEF(1,n) * o_eosb0%udcsat(n) + &
            i_coef%VCOEF(0,n) * i_coef%bsym * up3x * ( 1.0_pr - o_eosb0%udcsat(n) ) )
    enddo
    o_eosb0%a_vsym2(0) = o_eosb0%a_vsym2_order_n(0)
    do n = 1, i_coef%nmax
       o_eosb0%a_vsym2(n) = o_eosb0%a_vsym2(n-1) + o_eosb0%a_vsym2_order_n(n)
    enddo
    if (.not.i_opt_beta) then
       ! total symmetry energy: kinetic + potential
       do n = 0, i_coef%nmax
          o_eosb0%a_esym_b(n) = o_eosb0%tsym_b + o_eosb0%a_vsym(n)
       enddo
       o_eosb0%esym_b = o_eosb0%a_esym_b(i_coef%nmax)
       o_eosb0%esympot_b = o_eosb0%a_vsym(i_coef%nmax)
       ! total symmetry energy S_2: kinetic + potential
       do n = 0, i_coef%nmax
          o_eosb0%a_esym2_b(n) = o_eosb0%tsym2_b + o_eosb0%a_vsym2(n)
       enddo
       o_eosb0%esym2_b = o_eosb0%a_esym2_b(i_coef%nmax)
       o_eosb0%esym2pot_b = o_eosb0%a_vsym2(i_coef%nmax)
       !
       ! -------------------
       ! Entropy
       ! -------------------
       !
       o_eosb0%s2a_b = 0._pr
       o_eosb0%s2v_b = 0._pr
    endif
    !
    ! -------------------
    ! Chemical potentials
    ! -------------------
    ! mu_np
    o_eosb0%mut_np = i_coef%TFGsat * up3xp23 * ( d1_ffg +  up3x * d1_fsky )
    o_eosb0%mu_pot_np = 0.0_pr
    do n = 0, i_coef%nmax
       o_eosb0%mu_pot_np = o_eosb0%mu_pot_np + 2 * (i_eosDen%xx**n) / i_coef%facto(n) * ( &
            &   ( 2 * i_coef%VCOEF(1,n) * i_eosDen%xd + i_coef%VCOEF(2,n) * 4 * i_eosDen%xd**3 ) * o_eosb0%udc(n) + &
            &   2 * i_eosDen%xd * i_coef%bsym * vcoef(n) * up3x * ( 1 - o_eosb0%udc(n) ) )
    enddo
    o_eosb0%mu_np = o_eosb0%mut_np + o_eosb0%mu_pot_np
    if (.not.i_opt_beta) then
       ! mu_n and mu_p
       !
       o_eosb0%mut_n = i_coef%TFGsat / 6.0_pr * up3xp23 * ( 5 * ffg + 8 *  ( 1.0 + 3 * i_eosDen%xx ) * fsky ) &
            & + i_eosDen%den_p / i_eosDen%den_nuc * o_eosb0%mut_np
       o_eosb0%mut_p = i_coef%TFGsat / 6.0_pr * up3xp23 * ( 5 * ffg + 8 *  ( 1.0 + 3 * i_eosDen%xx ) * fsky ) &
            & - i_eosDen%den_n / i_eosDen%den_nuc * o_eosb0%mut_np
       !
       o_eosb0%mu_pot_n =   i_eosDen%den_p / i_eosDen%den_nuc * o_eosb0%mu_pot_np
       o_eosb0%mu_pot_p = - i_eosDen%den_n / i_eosDen%den_nuc * o_eosb0%mu_pot_np
       do n = 0, i_coef%nmax
          o_eosb0%mu_pot_n = o_eosb0%mu_pot_n + (i_eosDen%xx**n) / i_coef%facto(n) * vcoef(n) * &
               & ( o_eosb0%udc(n) +  up3x / 3.0_pr * o_eosb0%d1_udc(n) )
          o_eosb0%mu_pot_p = o_eosb0%mu_pot_p + (i_eosDen%xx**n) / i_coef%facto(n) * vcoef(n) * &
               & ( o_eosb0%udc(n) +  up3x / 3.0_pr * o_eosb0%d1_udc(n) )
       enddo
       if (i_coef%nmax.ge.1) then
          do n = 1, i_coef%nmax
             o_eosb0%mu_pot_n = o_eosb0%mu_pot_n + (i_eosDen%xx**(n-1)) / i_coef%facto(n-1) * vcoef(n) *&
                  & o_eosb0%udc(n) * up3x / 3.0_pr
             o_eosb0%mu_pot_p = o_eosb0%mu_pot_p + (i_eosDen%xx**(n-1)) / i_coef%facto(n-1) * vcoef(n) *&
                  & o_eosb0%udc(n) * up3x / 3.0_pr
          enddo
       endif
       !
       o_eosb0%mu_n = o_eosb0%mut_n + o_eosb0%mu_pot_n
       o_eosb0%mu_p = o_eosb0%mut_p + o_eosb0%mu_pot_p
       !
       ! -------------------
       ! Mean-Field
       ! -------------------
       !
       o_eosb0%umf_kin_n = 0.5_pr * CST_htm / i_coef%nsat * ( ( i_coef%mb - i_coef%db )  * o_eosb0%tau_b &
            & + 2._pr * i_coef%db * o_eosb0%tau_n )
       o_eosb0%umf_kin_p = 0.5_pr * CST_htm / i_coef%nsat * ( ( i_coef%mb - i_coef%db )  * o_eosb0%tau_b &
            & + 2._pr * i_coef%db * o_eosb0%tau_p )
       !
       o_eosb0%umf_n = o_eosb0%umf_kin_n + o_eosb0%mu_pot_n
       o_eosb0%umf_p = o_eosb0%umf_kin_p + o_eosb0%mu_pot_p
       !
       ! -------------------
       ! Effective mass
       ! -------------------
       !
       o_eosb0%ms_n = 1._pr / ( 1._pr + i_eosDen%den_nuc / i_coef%nsat * ( i_coef%mb + i_coef%db *  i_eosDen%xd  ) )
       o_eosb0%ms_p = 1._pr / ( 1._pr + i_eosDen%den_nuc / i_coef%nsat * ( i_coef%mb - i_coef%db *  i_eosDen%xd  ) )
       !
       ! -------------------
       ! Effective chemical potentials
       ! -------------------
       !
       o_eosb0%nu_n = 0.5_pr * CST_htm / o_eosb0%ms_n * i_eosDen%kf_n**2
       o_eosb0%nu_p = 0.5_pr * CST_htm / o_eosb0%ms_p * i_eosDen%kf_p**2
       o_eosb0%nu_b = 0.5_pr * CST_htm * i_eosDen%kf_nuc**2 ! baryon effective chemical potential without effective mass
       !
       ! -------------------
       ! First order derivatives
       ! -------------------
       !
       ! Kinetic energy : Fermi Gas
!       d1_t2aFG_b_den =  1._pr / 6._pr * i_coef%TFGsat / i_coef%nsat * ( 1.0 + 3 * i_eosDen%xx )**m13 * 2 * ffg
       ! Kinetic energy : total
!       d1_t2a_b_den = 1._pr / 6._pr * i_coef%TFGsat / i_coef%nsat * ( 1.0 + 3 * i_eosDen%xx )**m13 * &
!            & ( 2.0 * ffg + 5.0 * ( 1.0 + 3 * i_eosDen%xx ) * fsky )
    endif
!    d1_t2a_b_xe = - i_coef%TFGsat * ( 1.0 + 3 * i_eosDen%xx )**p23 * ( d1_ffg + ( 1.0 + 3 * i_eosDen%xx ) * d1_fsky )
!    if (.not.i_opt_beta) then
       ! Binding energy
!       o_eosb0%a_pv(:) = 0
!       o_eosb0%a_d1_e2a_b_den(0) =  d1_t2a_b_den
!       do n = 1, i_coef%nmax
!          tmp = p13 / i_coef%nsat * ( i_coef%VCOEF(0,n) + i_coef%VCOEF(1,n) * i_eosDen%xd**2 ) * &
!               & i_eosDen%xx**(n-1) / i_coef%facto(n-1)
!          o_eosb0%a_pv_order_n(n) = i_eosDen%den_nuc**2 * tmp
!          o_eosb0%a_pv(n) = o_eosb0%a_pv(n-1) + tmp
!          o_eosb0%a_d1_e2a_b_den(n) = o_eosb0%a_d1_e2a_b_den(n-1) + tmp
!       enddo
!    endif
!    do n = 0, i_coef%nmax
!       o_eosb0%a_d1_e2a_b_xe(n) = d1_t2a_b_xe - 4 * i_eosDen%xd * o_eosb0%a_vsym(n)
!    enddo
    if (.not.i_opt_beta) then
       !
       ! -------------------
       ! Pressure
       ! -------------------
       !
       o_eosb0%PFG_b = i_coef%nsat * i_coef%TFGsat / 3.0_pr * up3x * up3xp23 * ffg 
       o_eosb0%Pt_b = i_coef%nsat * i_coef%TFGsat / 3.0_pr * up3x * up3xp23 * &
            & ( ffg + 2.5_pr * up3x * fsky )

       o_eosb0%P_pot_b = vcoef(0) * o_eosb0%d1_udc(0)
       do n = 1, i_coef%nmax
          o_eosb0%P_pot_b = o_eosb0%P_pot_b + (i_eosDen%xx**(n-1)) / i_coef%facto(n) * &
               & vcoef(n) * ( n * o_eosb0%udc(n) + i_eosDen%xx * o_eosb0%d1_udc(n) )
       enddo
       o_eosb0%P_pot_b = i_coef%nsat / 3.0_pr * up3x**2 * o_eosb0%P_pot_b
       o_eosb0%p_b = o_eosb0%Pt_b + o_eosb0%P_pot_b
    endif
    !
    if (.not.i_opt_beta) then
       !
       ! -------------------
       ! Second order derivatives
       ! -------------------
       !
       ! Derivative of the Fermi Gas energy
!       d2_t2aFG_b_den =  - 1._pr / 9._pr * i_coef%TFGsat / (i_coef%nsat**2) * ( 1.0 + 3 * i_eosDen%xx )**m43 * ffg
       ! Derivative of the Kinetic energy
!       d2_t2a_b_den = 1._pr / 9._pr * i_coef%TFGsat / (i_coef%nsat**2) * ( 1.0 + 3 * i_eosDen%xx )**m43 * &
!            & ( - ffg + 5 * ( 1.0 + 3 * i_eosDen%xx ) * fsky )
!       d2_t2a_b_den_xe = m13 * i_coef%TFGsat / i_coef%nsat * ( 1.0 + 3 * i_eosDen%xx )**m13 * &
!            & ( 2.0 * d1_ffg + 5 * ( 1.0 + 3 * i_eosDen%xx ) * d1_fsky )
!       d2_t2a_b_xe = 2 * i_coef%TFGsat * ( 1.0 + 3 * i_eosDen%xx )**p23 * ( d2_ffg + ( 1.0 + 3 * i_eosDen%xx ) * d2_fsky )
       ! 2nd derivatives / den and potential incompressibility
!       sum = 0
!       o_eosb0%a_Kv(:) = 0
!       do n = 2, i_coef%nmax
!          tmp = ( i_coef%VCOEF(0,n) + i_coef%VCOEF(1,n) * i_eosDen%xd**2 ) * i_eosDen%xx**(n-2) / i_coef%facto(n-2)
!          sum = sum + tmp
!          d2_e2a_b_den(n) = d2_t2a_b_den + sum / ( 9.0 * i_coef%nsat**2 )
!          o_eosb0%a_Kv_order_n(n) = 9.0 * ( 1.0 + 3 * i_eosDen%xx )**2 * tmp + 18.0 / i_eosDen%den_nuc * o_eosb0%a_pv_order_n(n)
!          o_eosb0%a_Kv(n) = o_eosb0%a_Kv(n-1) + o_eosb0%a_Kv_order_n(n)
!       enddo
       ! 2nd derivative / den / xe
!       sum = 0
!       do n = 1, i_coef%nmax
!          sum = sum - 4.0 * i_eosDen%xd  * i_coef%VCOEF(1,n) * i_eosDen%xx**(n-1) / i_coef%facto(n-1)
!          d2_e2a_b_den_xe(n) = d2_t2a_b_den_xe + sum / ( 3 * i_coef%nsat )
!       enddo
       ! 2nd derivative / xe
!       do n = 0, i_coef%nmax
!          d2_e2a_b_xe(n) = d2_t2a_b_xe + 8.0 * o_eosb0%a_vsym(n)
!       enddo
       ! d1_mu_np_xe
!       do n = 1, i_coef%nmax
!          o_eosb0%a_d1_mu_np_xe(n) = - d2_e2a_b_xe(n)
!       enddo
       !
       ! -------------------
       ! Incompressibility
       ! -------------------
       !
       o_eosb0%KFG_b = - i_coef%TFGsat * up3xp23 * ffg + 18 * o_eosb0%PFG_b / i_eosDen%den_nuc
       o_eosb0%Kt_b = i_coef%TFGsat * up3xp23 * ( -ffg + 5 * ( 1.0 + 3 * i_eosDen%xx ) * fsky ) + &
            & 18 * o_eosb0%Pt_b / i_eosDen%den_nuc
       !
       o_eosb0%K_pot_b = 0.0_pr
       do n = 0, i_coef%nmax
          o_eosb0%K_pot_b = o_eosb0%K_pot_b + i_eosDen%xx**n / i_coef%facto(n) * &
               & vcoef(n) * o_eosb0%d2_udc(n)
       enddo
       if (i_coef%nmax.ge.1) then
          do n = 1, i_coef%nmax
             o_eosb0%K_pot_b = o_eosb0%K_pot_b + i_eosDen%xx**(n-1) / i_coef%facto(n-1) * &
                  & vcoef(n) * 2 * o_eosb0%d1_udc(n)
          enddo
       endif
       if (i_coef%nmax.ge.2) then
          do n = 2, i_coef%nmax
             o_eosb0%K_pot_b = o_eosb0%K_pot_b + i_eosDen%xx**(n-2) / i_coef%facto(n-2) * &
                  & vcoef(n) * o_eosb0%udc(n)
          enddo
       endif
       o_eosb0%K_pot_b = up3x**2 * o_eosb0%K_pot_b + 18 * o_eosb0%P_pot_b / i_eosDen%den_nuc
       !
       o_eosb0%K_b = o_eosb0%Kt_b +  o_eosb0%K_pot_b
       !       
       !    write(66,*)o_eosb0%k_b(4),d2e2a_b0_den2(4),o_eosb0%d1e2a_b_den(4),i_eosDen%den_nuc
       !
       ! -------------------
       ! sound velocity square (only for the baryon component)
       ! -------------------
       !
       o_eosb0%cs_b = o_eosb0%K_b / 9._pr / &
               &( CST_mnuc2 + o_eosb0%e2a_b + o_eosb0%p_b / i_eosDen%den_nuc )
       !
       ! -------------------
       ! Pressure derivatives
       ! -------------------
       !
!       do n = 2, i_coef%nmax
!          o_eosb0%a_d1_P_b_den(n) = 2.0 * i_eosDen%den_nuc * o_eosb0%a_d1_e2a_b_den(n) + i_eosDen%den_nuc**2 * d2_e2a_b_den(n)
!       enddo
!       do n = 1, i_coef%nmax
!          o_eosb0%a_d1_P_b_xe(n) =  i_eosDen%den_nuc**2 * d2_e2a_b_den_xe(n)
!       enddo
       !
       ! -------------------
       ! Third order derivative
       ! -------------------
       !
       !write(*,*)i_coef%nsat,1.0 + 3 * i_eosDen%xx,i_eosDen%den_nuc
!       d3_t2a_b_den =  i_coef%TFGsat / ( 27.0 * i_coef%nsat**3 ) * ( 1.0 + 3 * i_eosDen%xx )**m73 * &
!            & ( 4.0 * ffg - 5.0 * ( 1.0 + 3 * i_eosDen%xx ) * fsky )
!       o_eosb0%a_d3_e2a_b_den(:) = 0
!       o_eosb0%a_d3_e2a_b_den(2) = d3_t2a_b_den
!       sum = 0
!       do n = 3, i_coef%nmax
!          tmp = (i_coef%VCOEF(0,n) + i_coef%VCOEF(1,n) * i_eosDen%xd **2) * i_eosDen%xx**(n-3) / i_coef%facto(n-3)
!          sum = sum + tmp
!          o_eosb0%a_d3_e2a_b_den(n) = d3_t2a_b_den + sum / ( 27.0 * i_coef%nsat**3 )
!       enddo
       !
       ! -------------------
       ! Contribution of the rearrangement
       ! terms to the densities in den
       ! -------------------
       !
       !    if (i_rea.eq.1) then
       
       !    endif
       !

    endif
       
  end subroutine metaeos_T0_baryons_NR
     


  
  subroutine metaeos_T0_baryons_R(i_coef,i_eosDen,i_opt_beta,o_eosb0)

!    use acc; use cst;
!    use metaEosType;
    
    implicit none

    ! input variables
    type (metaEos_coef), intent(in)       :: i_coef
    type (metaEos_densities), intent(in)  :: i_eosDen
    logical, intent(in)                   :: i_opt_beta
    ! output variables
    type (metaEos_Baryons), intent(out)   :: o_eosb0
    !
    ! local variables
    !
!    logical                   :: i_opt_beta2

    integer                          :: n
!    real (kind=pr), dimension(0:5)  :: fac
    real (kind=pr)                   :: tmp, sum
    ! NEW
    ! kinetic energy factors
    real (kind=pr)                   :: coefC
    real (kind=pr)                   :: msnc2, mspc2, msNM, msSM
    real (kind=pr)                   :: gSM, gNM
    real (kind=pr)                   :: x, zn, zp, zsn, zsp, znSM, zpSM, zNM, zsnSM, zspSM, zsNM
    real (kind=pr)                   :: fn, fp, fsn, fsp, fnSM, fpSM, fNM, fsnSM, fspSM, fsNM
	real (kind=pr)					 :: T3,m,ms,ms_p,ms_n,z,z_p,z_n,kF,kF_p,f,kF_n,K,Kt_n,Kt_p,T3_p,T3_n
	real (kind=pr)					 :: rho_kin_n, rho_kin_p
	real (kind=pr)					 :: epsk, e2v_kin_n, e2v_kin_p,mterm,nterm
    ! OLD
    ! kinetic energy factors
    real (kind=pr)                   :: ffg, d1_ffg, d2_ffg
    real (kind=pr)                   :: fsky, d1_fsky, d2_fsky
    real (kind=pr)                   :: ffg_sym, fsky_sym
    ! first order derivatives
    real (kind=pr)                   :: d1_t2aFG_b_den
    real (kind=pr)                   :: d1_t2a_b_den, d1_t2a_b_xe
	real (kind=pr)					 :: d1n_ms_n,d1n_z_n,d1z_f_n,d1n_f_n,d1n_epsk_n
	real (kind=pr)					 :: d1n_ms_p,d1n_z_p,d1z_f_p,d1n_f_p,d1n_epsk_p
	real (kind=pr)					 :: d1d_f_n,d1d_ms_n,d1d_z_n,d1d_epsk_n,d1d_kF_n
	real (kind=pr)					 :: d1d_f_p,d1d_ms_p,d1d_z_p,d1d_epsk_p,d1d_kF_p
    ! second order derivatives
    real (kind=pr)                   :: d2_t2aFG_b_den
    real (kind=pr)                   :: d2_t2a_b_den, d2_t2a_b_den_xe, d2_t2a_b_xe
    real (kind=pr), dimension(0:4)   :: d2_e2a_b_den, d2_e2a_b_den_xe, d2_e2a_b_xe
   	real (kind=pr)					 :: d2n_ms_n,d2n_z_n,d2z_f_n,d2n_f_n,d2n_epsk_n
	real (kind=pr)					 :: d2n_ms_p,d2n_z_p,d2z_f_p,d2n_f_p,d2n_epsk_p
    ! third order derivatives
    real (kind=pr)                   :: d3_t2a_b_den
    real (kind=pr)                   :: chi, kappas, enthalpy
    real (kind=pr)                   :: coefn, coefp
    ! parameters
    real (kind=pr), parameter        :: p13=1._pr/3._pr,  p23=2./3.,  p43=4./3.,  p53=5./3.
    real (kind=pr), parameter        :: m13=-1./3., m23=-2./3., m43=-4./3., m53=-5./3., m73 = -7./3.
    ! vcoef: interaction coefficient
    real (kind=pr), dimension(0:4)   :: vcoef
    ! 1+3x
    real (kind=pr) :: up3x, up3xp23
   
    ! *******************
    ! *******************
    ! Baryons at T = 0
    ! *******************
    ! *******************
!    write(66,*)"Baryons"
!    i_opt_beta2=.true.
    ffg=0.0; d1_ffg=0.0; d2_ffg=0.0; fsky=0.0; d1_fsky=0.0; d2_fsky=0.0; ffg_sym=0.0; fsky_sym=0.0;
!

    ! For Relativitic Fermi gas Energy
    coefC = 1.0_pr / (8.0_pr * CST_pi2 * CST_hbc**3)
    if (i_coef%ikin.eq.0) then
       coefC = 0.0_pr
    endif

    !definitions of ms_p, ms_n, kF_p, kF_n,z_p,z_n 
    ms_p = CST_mpc2 / ( 1.0_pr + (i_coef%mb - i_coef%db*i_eosDen%xd)/CST_nsat*i_eosDen%den_nuc)
    ms_n = CST_mnc2 / ( 1.0_pr + (i_coef%mb + i_coef%db*i_eosDen%xd)/CST_nsat*i_eosDen%den_nuc)
    kF_p = ((CST_pi2*3._pr*i_eosDen%den_nuc)/2._pr*(1.0_pr - i_eosDen%xd))**CST_p13
    kF_n = ((CST_pi2*3._pr*i_eosDen%den_nuc)/2._pr*(1.0_pr + i_eosDen%xd))**CST_p13
    z_p = CST_hbc*kF_p/(ms_p)
    fsp = z_p*(1._pr+2._pr*z_p**2)*dsqrt(1._pr+z_p**2) - dlog(z_p+dsqrt(1._pr+z_p**2))
    z_n = CST_hbc*kF_n/(ms_n)
    rho_kin_p = coefC*ms_p**4*fsp
    e2v_kin_p = rho_kin_p - i_eosDen%den_nuc*(1.0_pr - i_eosDen%xd)/2._pr*CST_mpc2 
    fsn = z_n*(1._pr+2._pr*z_n**2)*dsqrt(1._pr+z_n**2) - dlog(z_n+dsqrt(1._pr+z_n**2))
    rho_kin_n = coefC*ms_n**4*fsn
    e2v_kin_n = rho_kin_n - i_eosDen%den_nuc*(1.0_pr + i_eosDen%xd)/2._pr*CST_mnc2
    T3_p = -1._pr
    T3_n = 1._pr
!

	   !allocation of some useful variables for symetric matter
	   !
       ! -------------------
       ! Proton 
       ! -------------------
	   !
    z = z_p
    K = i_coef%mb - i_coef%db*i_eosDen%xd
    ms = ms_p
!    if (i_eosDen%xd.ge.1) then
!       ms=0._pr
!    endif
    m = CST_mpc2
    kF=kF_p
    T3 = T3_p
    f=fsp
    epsk = e2v_kin_p
    mterm = m*(1._pr - i_eosDen%xd)/2._pr
    nterm = -i_eosDen%den_nuc/2._pr*CST_mpc2
!
    !derivatives of effective mass with respect to density at constant asymmetry
    d1n_ms_p = -m*(K/CST_nsat)*(ms/m)**2
    d2n_ms_p = m*2.0_pr*(K/CST_nsat)**2*(ms/m)**3
!
    !derivatives of z with respect to density at constant asymetry
    d1n_z_p = z/(3.0_pr*i_eosDen%den_nuc) + (CST_hbc/m)*(K)/CST_nsat*kF
    d2n_z_p = -z/(3.0_pr*i_eosDen%den_nuc**2) + d1n_z_p/(3.0_pr*i_eosDen%den_nuc) &
         &+ (CST_hbc/(3.0_pr*m))*(K)/CST_nsat*kF/i_eosDen%den_nuc    
!
    !derivatives of f with respect to z 
    d1z_f_p = 8._pr*z**2*(dsqrt(1._pr+z**2))
    d2z_f_p = 8._pr*z*((2._pr+3._pr*z**2)/dsqrt(1._pr+z**2))
!
    !derivatives of f with respect to density at constant asymetry
    d1n_f_p = d1z_f_p*d1n_z_p 
    d2n_f_p = d2z_f_p*d1n_z_p**2 + d1z_f_p*d2n_z_p
    !
    !derivatives of kinetic energy volumique density with respect to density at constant asymetry
    d1n_epsk_p = 4._pr*coefC*ms**3*d1n_ms_p*f + coefC*ms**4*d1n_f_p - mterm
    d2n_epsk_p = 12._pr*coefC*ms**2*d1n_ms_p**2*f + 4._pr*coefC*ms**3*(d2n_ms_p*f+2._pr*d1n_ms_p*d1n_f_p)+coefC*ms**4*d2n_f_p
!
    !derivative of kF with respect to constant constant asymetry at constant density 
    d1d_kF_p=T3*kF/(3._pr*(1._pr+T3*i_eosDen%xd))
    if (i_eosDen%xd.ge.1.0_pr) then
       d1d_kF_p=0._pr
    endif
!
    !derivative of ms with respect to constant asymetry at constant density
    d1d_ms_p=-T3*i_coef%db*i_eosDen%den_nuc/CST_nsat*m*(ms/m)**2
!
    !derivative of z with respect of constant asymetry at constant density 
    d1d_z_p=CST_hbc/(ms)*d1d_kF_p-CST_hbc*kF/(ms)**2*d1d_ms_p
!    if (i_eosDen%xd.ge.1) then
!       d1d_z_p=0._pr
!    endif
!
    !derivative of f fonction with respect to constant asymetry at constant density
    d1d_f_p=d1d_z_p*d1z_f_p
    !	   
    !derivatives of kinetic energy per unit of volume density with respect to asymetry at constant density 
    d1d_epsk_p=4._pr*coefC*(ms)**3*d1d_ms_p*f+coefC*(ms)**4*d1d_f_p - nterm
!
!
!
!
!
!
!
    ! -------------------
    ! Neutron 
    ! -------------------
    !
    z = z_n
    K = i_coef%mb + i_coef%db*i_eosDen%xd
    ms = ms_n
    m = CST_mnc2
    kF=kF_n
    T3 = T3_n
    f=fsn
    epsk = e2v_kin_n
    mterm = m*(1._pr + i_eosDen%xd)/2._pr
    nterm = i_eosDen%den_nuc/2._pr*CST_mnc2
!
    !derivatives of effective mass with respect to density at constant asymmetry
    d1n_ms_n = -m*(K/CST_nsat)*(ms/m)**2
    d2n_ms_n = m*2.0_pr*(K/CST_nsat)**2*(ms/m)**3   
!
    !derivatives of z with respect to density at constant asymetry
    d1n_z_n = z/(3.0_pr*i_eosDen%den_nuc) + (CST_hbc/m)*(K)/CST_nsat*kF
    d2n_z_n = -z/(3.0_pr*i_eosDen%den_nuc**2) + d1n_z_n/(3.0_pr*i_eosDen%den_nuc) &
         &+ (CST_hbc/(3.0_pr*m))*(K)/CST_nsat*kF/i_eosDen%den_nuc    
!
    !derivatives of f with respect to z 
    d1z_f_n = 8._pr*z**2*(dsqrt(1._pr+z**2))
    d2z_f_n = 8._pr*z*((2._pr+3._pr*z**2)/dsqrt(1._pr+z**2))
!
    !derivatives of f with respect to density at constant asymetry
    d1n_f_n = d1z_f_n*d1n_z_n 
    d2n_f_n = d2z_f_n*d1n_z_n**2 + d1z_f_n*d2n_z_n
    !
    !derivatives of kinetic energy volumique density with respect to density at constant asymetry
    d1n_epsk_n = 4._pr*coefC*ms**3*d1n_ms_n*f + coefC*ms**4*d1n_f_n - mterm
    d2n_epsk_n = 12._pr*coefC*ms**2*d1n_ms_n**2*f + 4._pr*coefC*ms**3*(d2n_ms_n*f+2._pr*d1n_ms_n*d1n_f_n)+coefC*ms**4*d2n_f_n
!
    !derivative of kF with respect to constant constant asymetry at constant density 
    d1d_kF_n=T3*kF/(3._pr*(1+T3*i_eosDen%xd))
    if (i_eosDen%xd.le.-1.0_pr) then
       d1d_kF_n = 0._pr
    endif
!
    !derivative of ms with respect to constant asymetry at constant density
    d1d_ms_n=-T3*i_coef%db*i_eosDen%den_nuc/CST_nsat*m*(ms/m)**2
!
    !derivative of z with respect of constant asymetry at constant density 
    d1d_z_n=CST_hbc/(ms)*d1d_kF_n-CST_hbc*kF/(ms)**2*d1d_ms_n
!
    !derivative of f fonction with respect to constant asymetry at constant density
    d1d_f_n=d1d_z_n*d1z_f_n
    !	   
    !derivatives of kinetic energy per unit of volume density with respect to asymetry at constant density 
    d1d_epsk_n=4._pr*coefC*(ms)**3*d1d_ms_n*f+coefC*(ms)**4*d1d_f_n - nterm
!
!       write(*,*)" test "
!       write(*,*)
    !
    ! -------------------
    ! Density correction u_alpha(x,delta)
    ! -------------------
    !
    do n = 0, i_coef%nmax
       o_eosb0%udc(n) = ( 1._pr - ( - 3._pr * i_eosDen%xx )**(i_coef%nmax+1-n) * &
            & dexp( - ( i_coef%bsat + i_coef%bsym * i_eosDen%xd**2 ) * i_eosDen%den_nuc / i_coef%nsat ) ) 
       o_eosb0%udcsat(n) = ( 1._pr - ( - 3._pr * i_eosDen%xx )**(i_coef%nmax+1-n) * &
            & dexp( - i_coef%bsat * i_eosDen%den_nuc / i_coef%nsat ) ) 
       o_eosb0%udcNM(n) = ( 1._pr - ( - 3._pr * i_eosDen%xx )**(i_coef%nmax+1-n) * &
            & dexp( - (i_coef%bsat+i_coef%bsym) * i_eosDen%den_nuc / i_coef%nsat ) ) 
       ! first derivative:
       o_eosb0%d1_udc(n) = 0.0_pr
       o_eosb0%d2_udc(n) = 0.0_pr
       if (i_eosDen%xx.ne.0.0_pr) then
          o_eosb0%d1_udc(n) = ( i_coef%nmax + 1 - n - &
               & 3 * ( i_coef%bsat + i_coef%bsym * i_eosDen%xd**2 ) * i_eosDen%xx ) * &
               & ( o_eosb0%udc(n) - 1.0_pr ) / i_eosDen%xx
          o_eosb0%d2_udc(n) = ( - ( i_coef%nmax + 1 - n ) * ( i_coef%nmax - n ) + &
               & 6 * ( i_coef%bsat + i_coef%bsym * i_eosDen%xd**2 ) * i_eosDen%xx * ( i_coef%nmax + 1 - n ) - &
               & 9 * ( i_coef%bsat + i_coef%bsym * i_eosDen%xd**2 )**2 * i_eosDen%xx**2 ) * &
               & ( 1.0_pr - o_eosb0%udc(n) ) / i_eosDen%xx**2
       endif
    enddo
    !
    ! -------------------
    ! full interaction parameter v(x,delta)
    ! -------------------
    !
    do n = 0, i_coef%nmax
       vcoef(n) = i_coef%VCOEF(0,n) + i_coef%VCOEF(1,n) * i_eosDen%xd**2 + i_coef%VCOEF(2,n) * i_eosDen%xd**4
    enddo
    !
    ! -------------------
    ! 1+3x, (1+3x)^2/3
    ! -------------------
    up3x = ( 1.0_pr + 3 * i_eosDen%xx )
    up3xp23 = ( 1.0_pr + 3 * i_eosDen%xx )**CST_p23
!       up3x = 1.0_pr
    !       up3xp23 = 1.0_pr
    !
    ! -------------------
    ! Effective mass
    ! -------------------
    !
    msnc2 = CST_mnc2 / ( 1.0_pr + (i_coef%mb + i_coef%db * i_eosDen%xd ) * up3x )
    mspc2 = CST_mpc2 / ( 1.0_pr + (i_coef%mb - i_coef%db * i_eosDen%xd ) * up3x )
    !
    ! -------------------
    ! function f: free Fermi gas
    ! -------------------
    !
    zn = CST_hbc * i_eosDen%kf_n / CST_mnc2
    zp = CST_hbc * i_eosDen%kf_p / CST_mpc2
    !
    fn = zn * ( 1.0_pr + 2.0_pr*zn**2 ) * dsqrt( 1.0_pr + zn**2 ) - dlog ( zn + dsqrt( 1.0_pr + zn**2 )) 
    fp = zp * ( 1.0_pr + 2.0_pr*zp**2 ) * dsqrt( 1.0_pr + zp**2 ) - dlog ( zp + dsqrt( 1.0_pr + zp**2 )) 
    !
    ! -------------------
    ! function fs: including effective mass
    ! -------------------
    !fsn and fsp have already be defined
    !
    if (.not.i_opt_beta) then
       !
       ! -------------------
       ! energy per particle
       ! -------------------
       !
       ! kinetic energy : Fermi Gas
!       o_eosb0%t2aFG_b = coefC * ( CST_mnc2**4 * fn + CST_mpc2**4 * fp ) / i_eosDen%den_nuc
       ! kinetic energy : total
!       o_eosb0%t2a_b = coefC * ( ms_n**4 * fsn + ms_p**4 * fsp ) / i_eosDen%den_nuc - CST_mnc2*(1._pr + i_eosDen%xd)/2._pr &
!                        & - CST_mpc2*(1._pr - i_eosDen%xd)/2._pr
       ! kinetic energy per unit volum : Fermi Gas
       !o_eosb0%t2vFG_b = o_eosb0%t2aFG_b * i_eosDen%den_nuc
       !
       ! kinetic energy per unit volum : total
!       o_eosb0%t2v_b = o_eosb0%t2a_b * i_eosDen%den_nuc
       o_eosb0%t2v_b = e2v_kin_n + e2v_kin_p
       ! kinetic energy : total
       o_eosb0%t2a_b = o_eosb0%t2v_b / i_eosDen%den_nuc

       ! potential energy and total energy
       do n = 0, i_coef%nmax
          o_eosb0%a_v2a_order_n(n) = vcoef(n) * (i_eosDen%xx**n) / i_coef%facto(n) * o_eosb0%udc(n)
          if (n.eq.0) then
             o_eosb0%a_v2a(n) = o_eosb0%a_v2a_order_n(n)
          else
             o_eosb0%a_v2a(n) = o_eosb0%a_v2a(n-1) + o_eosb0%a_v2a_order_n(n)
          endif
          o_eosb0%a_e2a_b(n) = o_eosb0%t2a_b + o_eosb0%a_v2a(n)
          o_eosb0%a_e2v_b(n) = o_eosb0%t2v_b + i_eosDen%den_nuc * o_eosb0%a_v2a(n)
       enddo
       o_eosb0%epot2a_b = o_eosb0%a_v2a(i_coef%nmax)
       o_eosb0%epot2v_b = i_eosDen%den_nuc * o_eosb0%a_v2a(i_coef%nmax)
       o_eosb0%e2a_b = o_eosb0%a_e2a_b(i_coef%nmax)
       o_eosb0%e2v_b = o_eosb0%e2a_b * i_eosDen%den_nuc
       o_eosb0%rho_b = i_eosDen%den_p * CST_mpc2 + i_eosDen%den_n * CST_mnc2 + o_eosb0%e2v_b
       !
       ! -------------------
       ! symmetry energy per particle
       ! -------------------
       !
       ! FG contribution
       zNM = CST_hbc / CST_mnc2 * i_eosDen%kf_nuc * 2**CST_p13 
       znSM = CST_hbc / CST_mnc2 * i_eosDen%kf_nuc 
       zpSM = CST_hbc / CST_mpc2 * i_eosDen%kf_nuc 
       x=zNM
       fNM = x*(1._pr+2._pr*x**2)*dsqrt(1._pr+x**2) - dlog(x+dsqrt(1._pr+x**2))
       x=znSM
       fnSM = x*(1._pr+2._pr*x**2)*dsqrt(1._pr+x**2) - dlog(x+dsqrt(1._pr+x**2))
       x=zpSM
       fpSM = x*(1._pr+2._pr*x**2)*dsqrt(1._pr+x**2) - dlog(x+dsqrt(1._pr+x**2))
       o_eosb0%tsymFG_b = coefC * CST_mnc2**4 * fNM - coefC * ( CST_mnc2**4 * fnSM + CST_mpc2**4 * fpSM ) &
            & - (  CST_mnc2 -CST_mpc2 ) / 2.0_pr
!       o_eosb0%tsymFG_b = i_coef%TFGsat * up3xp23 * ( 2**CST_p23 - 1.0_pr )
!       o_eosb0%tsym2FG_b = 0.25_pr * i_coef%TFGsat * ( 1.0 + 3 * i_eosDen%xx )**p23 * ffg_sym
       ! Kinetic contribution
!       o_eosb0%tsym_b = i_coef%TFGsat * ( 1.0_pr + i_coef%kappaNM )!* ( 2**CST_p23 * ( 1.0_pr + up3x*i_coef%kappaNM ) )!&
!       o_eosb0%tsym_b = i_coef%TFGsat * up3xp23 * ( 1.0_pr + up3x*i_coef%kappaNM )!* ( 2**CST_p23 * ( 1.0_pr + up3x*i_coef%kappaNM ) )!&
!         &            - ( 1.0_pr + up3x*i_coef%kappas ) )
       gSM = 1.0_pr + (i_coef%mb)*up3x
       gNM = 1.0_pr + (i_coef%mb + i_coef%db )*up3x
       zsNM = CST_hbc * i_eosDen%kf_nuc / CST_mnc2 * 2.0_pr**CST_p13 * gNM 
       zsnSM = CST_hbc * i_eosDen%kf_nuc / CST_mnc2 * gSM 
       zspSM = CST_hbc * i_eosDen%kf_nuc / CST_mpc2 * gSM 
       x=zsNM
       fsNM = x*(1._pr+2._pr*x**2)*dsqrt(1._pr+x**2) - dlog(x+dsqrt(1._pr+x**2))
       x=zsnSM
       fsnSM = x*(1._pr+2._pr*x**2)*dsqrt(1._pr+x**2) - dlog(x+dsqrt(1._pr+x**2))
       x=zspSM
       fspSM = x*(1._pr+2._pr*x**2)*dsqrt(1._pr+x**2) - dlog(x+dsqrt(1._pr+x**2))
       o_eosb0%tsym_b = coefC * (CST_mnc2/gNM)**4 * fsNM / i_eosDen%den_nuc &
            & - coefC * ( (CST_mnc2/gSM)**4 * fsnSM + (CST_mpc2/gSM)**4 * fspSM ) / i_eosDen%den_nuc &
            & - (  CST_mnc2 -CST_mpc2 ) / 2.0_pr
    !*i_eosDen%den_nuc - 2*coefC*msSM**4*fsSM/i_eosDen%den_nuc &
    !      &+ CST_mnuc2*i_eosDen%den_nuc
!       o_eosb0%tsym_b = i_coef%TFGsat * up3xp23 * ( 2**CST_p23 * ( 1.0_pr + up3x*i_coef%kappaNM ) - &
!         &            ( 1.0_pr + up3x*i_coef%kappas ) )
!       write(*,*)'check2:',i_coef%kappas,i_coef%kappaNM,i_coef%TFGsat,2**CST_p23,up3xp23 
!       write(*,*)'check2:',i_eosDen%den_nuc,o_eosb0%tsym_b
!       write(*,*)'check2:', o_eosb0%udcsat(0), o_eosb0%udcNM(0)
       !
       ! -------------------
       ! symmetry energy S_2 per particle
       ! -------------------
       !
       ! FG contribution
       o_eosb0%tsym2FG_b = 5.0_pr/9.0_pr * i_coef%TFGsat * up3xp23 
!       o_eosb0%tsym2FG_b = 0.25_pr * i_coef%TFGsat * ( 1.0 + 3 * i_eosDen%xx )**p23 * ffg_sym
       ! Kinetic contribution
       o_eosb0%tsym2_b = 5.0_pr/9.0_pr * i_coef%TFGsat * up3xp23 * ( 1.0_pr &
            & + up3x * (i_coef%mb + 3 * i_coef%db) )
!       o_eosb0%tsym2_b = 0.25_pr * i_coef%TFGsat * ( 1.0 + 3 * i_eosDen%xx )**p23 * &
!            & ( ffg_sym + ( 1.0 + 3 * i_eosDen%xx ) * fsky_sym )
    endif
    ! potential contribution
    do n = 0, i_coef%nmax
       o_eosb0%a_vsym_order_n(n) = (i_eosDen%xx**n) / i_coef%facto(n) * ( &
            ( i_coef%VCOEF(0,n) + i_coef%VCOEF(1,n) + i_coef%VCOEF(2,n) ) * o_eosb0%udcNM(n) - &
            ( i_coef%VCOEF(0,n) ) * o_eosb0%udcsat(n) )
    enddo
    o_eosb0%a_vsym(0) = o_eosb0%a_vsym_order_n(0)
    do n = 1, i_coef%nmax
       o_eosb0%a_vsym(n) = o_eosb0%a_vsym(n-1) + o_eosb0%a_vsym_order_n(n)
    enddo
    ! potential contribution S_2
    do n = 0, i_coef%nmax
       o_eosb0%a_vsym2_order_n(n) = (i_eosDen%xx**n) / i_coef%facto(n) * ( i_coef%VCOEF(1,n) * o_eosb0%udcsat(n) + &
            i_coef%VCOEF(0,n) * i_coef%bsym * up3x * ( 1.0_pr - o_eosb0%udcsat(n) ) )
    enddo
    o_eosb0%a_vsym2(0) = o_eosb0%a_vsym2_order_n(0)
    do n = 1, i_coef%nmax
       o_eosb0%a_vsym2(n) = o_eosb0%a_vsym2(n-1) + o_eosb0%a_vsym2_order_n(n)
    enddo
    if (.not.i_opt_beta) then
       ! total symmetry energy: kinetic + potential
       do n = 0, i_coef%nmax
          o_eosb0%a_esym_b(n) = o_eosb0%tsym_b + o_eosb0%a_vsym(n)
       enddo
       o_eosb0%esym_b = o_eosb0%a_esym_b(i_coef%nmax)
       o_eosb0%esympot_b = o_eosb0%a_vsym(i_coef%nmax)
       ! total symmetry energy S_2: kinetic + potential
       do n = 0, i_coef%nmax
          o_eosb0%a_esym2_b(n) = o_eosb0%tsym2_b + o_eosb0%a_vsym2(n)
       enddo
       o_eosb0%esym2_b = o_eosb0%a_esym2_b(i_coef%nmax)
       o_eosb0%esym2pot_b = o_eosb0%a_vsym2(i_coef%nmax)
       !
       ! -------------------
       ! Entropy
       ! -------------------
       !
       o_eosb0%s2a_b = 0._pr
       o_eosb0%s2v_b = 0._pr
    endif
    !
    ! -------------------
    ! Chemical potentials
    ! -------------------
    ! mu_np
    o_eosb0%mut_np = 2._pr/i_eosDen%den_nuc * (d1d_epsk_n + d1d_epsk_p)
    o_eosb0%mu_pot_np = 0.0_pr
    do n = 0, i_coef%nmax
       o_eosb0%mu_pot_np = o_eosb0%mu_pot_np + 2 * (i_eosDen%xx**n) / i_coef%facto(n) * ( &
            &   ( 2 * i_coef%VCOEF(1,n) * i_eosDen%xd + i_coef%VCOEF(2,n) * 4 * i_eosDen%xd**3 ) * o_eosb0%udc(n) + &
            &   2 * i_eosDen%xd * i_coef%bsym * vcoef(n) * up3x * ( 1 - o_eosb0%udc(n) ) )
    enddo
    o_eosb0%mu_np = o_eosb0%mut_np + o_eosb0%mu_pot_np
    if (.not.i_opt_beta) then
       ! mu_n and mu_p
       !
!       write(*,*)'test1 (mu):',epsk_kin_n,epsk_kin_p,d1n_epsk_n,d1n_epsk_p
       o_eosb0%mut_n= d1n_epsk_n + d1n_epsk_p + T3_n*(1._pr + T3_p*i_eosDen%xd)/i_eosDen%den_nuc *(d1d_epsk_n + d1d_epsk_p)
       o_eosb0%mut_p= d1n_epsk_p + d1n_epsk_n + T3_p*(1._pr + T3_n*i_eosDen%xd)/i_eosDen%den_nuc *(d1d_epsk_p + d1d_epsk_n) 
       !
       o_eosb0%mu_pot_is = 0._pr
       do n = 0, i_coef%nmax
          o_eosb0%mu_pot_is = o_eosb0%mu_pot_is + (i_eosDen%xx**n) / i_coef%facto(n) * vcoef(n) * &
               & ( o_eosb0%udc(n) +  up3x / 3.0_pr * o_eosb0%d1_udc(n) )
       enddo
       if (i_coef%nmax.ge.1) then
          do n = 1, i_coef%nmax
             o_eosb0%mu_pot_is = o_eosb0%mu_pot_is + (i_eosDen%xx**(n-1)) / i_coef%facto(n-1) * vcoef(n) *&
                  & o_eosb0%udc(n) * up3x / 3.0_pr
          enddo
       endif
       !
       o_eosb0%mu_pot_n = o_eosb0%mu_pot_is + i_eosDen%den_p / i_eosDen%den_nuc * o_eosb0%mu_pot_np
       o_eosb0%mu_pot_p = o_eosb0%mu_pot_is - i_eosDen%den_n / i_eosDen%den_nuc * o_eosb0%mu_pot_np
       !
!       do n = 0, i_coef%nmax
!          o_eosb0%mu_pot_n = o_eosb0%mu_pot_n + (i_eosDen%xx**n) / i_coef%facto(n) * vcoef(n) * &
!               & ( o_eosb0%udc(n) +  up3x / 3.0_pr * o_eosb0%d1_udc(n) )
!          o_eosb0%mu_pot_p = o_eosb0%mu_pot_p + (i_eosDen%xx**n) / i_coef%facto(n) * vcoef(n) * &
!               & ( o_eosb0%udc(n) +  up3x / 3.0_pr * o_eosb0%d1_udc(n) )
!       enddo
!       if (i_coef%nmax.ge.1) then
!          do n = 1, i_coef%nmax
!             o_eosb0%mu_pot_n = o_eosb0%mu_pot_n + (i_eosDen%xx**(n-1)) / i_coef%facto(n-1) * vcoef(n) *&
!                  & o_eosb0%udc(n) * up3x / 3.0_pr
!             o_eosb0%mu_pot_p = o_eosb0%mu_pot_p + (i_eosDen%xx**(n-1)) / i_coef%facto(n-1) * vcoef(n) *&
!                  & o_eosb0%udc(n) * up3x / 3.0_pr
!          enddo
!       endif
       !
       o_eosb0%mu_n = o_eosb0%mut_n + o_eosb0%mu_pot_n
       o_eosb0%mu_p = o_eosb0%mut_p + o_eosb0%mu_pot_p
       !
       ! -------------------
       ! Effective mass
       ! -------------------
       !
       o_eosb0%ms_n = ms_n
       o_eosb0%ms_p = ms_p
       !
       ! -------------------
       ! Mean-Field
       ! -------------------
       !
       o_eosb0%umf_kin_n = o_eosb0%mut_n-(CST_hbc*kF_n)**2/(2*ms_n)
       o_eosb0%umf_kin_p = o_eosb0%mut_p-(CST_hbc*kF_p)**2/(2*ms_p)
       !
       o_eosb0%umf_n = o_eosb0%umf_kin_n + o_eosb0%mu_pot_n
       o_eosb0%umf_p = o_eosb0%umf_kin_p + o_eosb0%mu_pot_p
       !
       ! -------------------
       ! Effective chemical potentials
       ! -------------------
       !
       o_eosb0%nu_n = 0.5_pr * CST_htm / o_eosb0%ms_n * i_eosDen%kf_n**2
       o_eosb0%nu_p = 0.5_pr * CST_htm / o_eosb0%ms_p * i_eosDen%kf_p**2
       o_eosb0%nu_b = 0.5_pr * CST_htm * i_eosDen%kf_nuc**2 ! baryon effective chemical potential without effective mass
       !
       ! -------------------
       ! First order derivatives
       ! -------------------
       !
       ! Kinetic energy : Fermi Gas
!       d1_t2aFG_b_den =  1._pr / 6._pr * i_coef%TFGsat / i_coef%nsat * ( 1.0 + 3 * i_eosDen%xx )**m13 * 2 * ffg
       ! Kinetic energy : total
!       d1_t2a_b_den = 1._pr / 6._pr * i_coef%TFGsat / i_coef%nsat * ( 1.0 + 3 * i_eosDen%xx )**m13 * &
!            & ( 2.0 * ffg + 5.0 * ( 1.0 + 3 * i_eosDen%xx ) * fsky )
    endif
!    d1_t2a_b_xe = - i_coef%TFGsat * ( 1.0 + 3 * i_eosDen%xx )**p23 * ( d1_ffg + ( 1.0 + 3 * i_eosDen%xx ) * d1_fsky )
!    if (.not.i_opt_beta) then
       ! Binding energy
!       o_eosb0%a_pv(:) = 0
!       o_eosb0%a_d1_e2a_b_den(0) =  d1_t2a_b_den
!       do n = 1, i_coef%nmax
!          tmp = p13 / i_coef%nsat * ( i_coef%VCOEF(0,n) + i_coef%VCOEF(1,n) * i_eosDen%xd**2 ) * &
!               & i_eosDen%xx**(n-1) / i_coef%facto(n-1)
!          o_eosb0%a_pv_order_n(n) = i_eosDen%den_nuc**2 * tmp
!          o_eosb0%a_pv(n) = o_eosb0%a_pv(n-1) + tmp
!          o_eosb0%a_d1_e2a_b_den(n) = o_eosb0%a_d1_e2a_b_den(n-1) + tmp
!       enddo
!    endif
!    do n = 0, i_coef%nmax
!       o_eosb0%a_d1_e2a_b_xe(n) = d1_t2a_b_xe - 4 * i_eosDen%xd * o_eosb0%a_vsym(n)
!    enddo
    if (.not.i_opt_beta) then
       !
       ! -------------------
       ! Pressure
       ! -------------------
       !
       !o_eosb0%PFG_b = i_coef%nsat * i_coef%TFGsat / 3.0_pr * up3x * up3xp23 * ffg
!       write(*,*)'test2 (p):',epsk_kin_n/i_eosDen%den_nuc,epsk_kin_p/i_eosDen%den_nuc,d1n_epsk_n,d1n_epsk_p
       o_eosb0%Pt_b = - e2v_kin_n + i_eosDen%den_nuc * d1n_epsk_n - e2v_kin_p + i_eosDen%den_nuc * d1n_epsk_p
!       write(*,*)'pression (kin):',o_eosb0%Pt_b
       o_eosb0%P_pot_b = vcoef(0) * o_eosb0%d1_udc(0)
       do n = 1, i_coef%nmax
          o_eosb0%P_pot_b = o_eosb0%P_pot_b + (i_eosDen%xx**(n-1)) / i_coef%facto(n) * &
               & vcoef(n) * ( n * o_eosb0%udc(n) + i_eosDen%xx * o_eosb0%d1_udc(n) )
       enddo
       o_eosb0%P_pot_b = i_coef%nsat / 3.0_pr * up3x**2 * o_eosb0%P_pot_b
!       write(*,*)'pression (pot):',o_eosb0%P_pot_b
       o_eosb0%p_b = o_eosb0%Pt_b + o_eosb0%P_pot_b
!       write(*,*)'pression (tot):',o_eosb0%p_b
    endif
    !
    if (.not.i_opt_beta) then
       !
       ! -------------------
       ! Second order derivatives
       ! -------------------
       !
       ! Derivative of the Fermi Gas energy
!       d2_t2aFG_b_den =  - 1._pr / 9._pr * i_coef%TFGsat / (i_coef%nsat**2) * ( 1.0 + 3 * i_eosDen%xx )**m43 * ffg
       ! Derivative of the Kinetic energy
!       d2_t2a_b_den = 1._pr / 9._pr * i_coef%TFGsat / (i_coef%nsat**2) * ( 1.0 + 3 * i_eosDen%xx )**m43 * &
!            & ( - ffg + 5 * ( 1.0 + 3 * i_eosDen%xx ) * fsky )
!       d2_t2a_b_den_xe = m13 * i_coef%TFGsat / i_coef%nsat * ( 1.0 + 3 * i_eosDen%xx )**m13 * &
!            & ( 2.0 * d1_ffg + 5 * ( 1.0 + 3 * i_eosDen%xx ) * d1_fsky )
!       d2_t2a_b_xe = 2 * i_coef%TFGsat * ( 1.0 + 3 * i_eosDen%xx )**p23 * ( d2_ffg + ( 1.0 + 3 * i_eosDen%xx ) * d2_fsky )
       ! 2nd derivatives / den and potential incompressibility
!       sum = 0
!       o_eosb0%a_Kv(:) = 0
!       do n = 2, i_coef%nmax
!          tmp = ( i_coef%VCOEF(0,n) + i_coef%VCOEF(1,n) * i_eosDen%xd**2 ) * i_eosDen%xx**(n-2) / i_coef%facto(n-2)
!          sum = sum + tmp
!          d2_e2a_b_den(n) = d2_t2a_b_den + sum / ( 9.0 * i_coef%nsat**2 )
!          o_eosb0%a_Kv_order_n(n) = 9.0 * ( 1.0 + 3 * i_eosDen%xx )**2 * tmp + 18.0 / i_eosDen%den_nuc * o_eosb0%a_pv_order_n(n)
!          o_eosb0%a_Kv(n) = o_eosb0%a_Kv(n-1) + o_eosb0%a_Kv_order_n(n)
!       enddo
       ! 2nd derivative / den / xe
!       sum = 0
!       do n = 1, i_coef%nmax
!          sum = sum - 4.0 * i_eosDen%xd  * i_coef%VCOEF(1,n) * i_eosDen%xx**(n-1) / i_coef%facto(n-1)
!          d2_e2a_b_den_xe(n) = d2_t2a_b_den_xe + sum / ( 3 * i_coef%nsat )
!       enddo
       ! 2nd derivative / xe
!       do n = 0, i_coef%nmax
!          d2_e2a_b_xe(n) = d2_t2a_b_xe + 8.0 * o_eosb0%a_vsym(n)
!       enddo
       ! d1_mu_np_xe
!       do n = 1, i_coef%nmax
!          o_eosb0%a_d1_mu_np_xe(n) = - d2_e2a_b_xe(n)
!       enddo
       !
       ! -------------------
       ! Incompressibility
       ! -------------------
       !
       !       o_eosb0%KFG_b = - i_coef%TFGsat * up3xp23 * ffg + 18 * o_eosb0%PFG_b / i_eosDen%den_nuc
!       write(*,*)"test K:",d2n_epsk_n,d2n_epsk_p

       o_eosb0%Kt_b = 9._pr * i_eosDen%den_nuc * d2n_epsk_n + 9._pr * i_eosDen%den_nuc * d2n_epsk_p
!       write(*,*)"test K:",d2n_epsk_n,d2n_epsk_p,o_eosb0%Kt_b
       !
       o_eosb0%K_pot_b = 0.0_pr
       do n = 0, i_coef%nmax
          o_eosb0%K_pot_b = o_eosb0%K_pot_b + i_eosDen%xx**n / i_coef%facto(n) * &
               & vcoef(n) * o_eosb0%d2_udc(n)
       enddo
       if (i_coef%nmax.ge.1) then
          do n = 1, i_coef%nmax
             o_eosb0%K_pot_b = o_eosb0%K_pot_b + i_eosDen%xx**(n-1) / i_coef%facto(n-1) * &
                  & vcoef(n) * 2 * o_eosb0%d1_udc(n)
          enddo
       endif
       if (i_coef%nmax.ge.2) then
          do n = 2, i_coef%nmax
             o_eosb0%K_pot_b = o_eosb0%K_pot_b + i_eosDen%xx**(n-2) / i_coef%facto(n-2) * &
                  & vcoef(n) * o_eosb0%udc(n)
          enddo
       endif
       o_eosb0%K_pot_b = up3x**2 * o_eosb0%K_pot_b + 18 * o_eosb0%P_pot_b / i_eosDen%den_nuc
       !
       o_eosb0%K_b = o_eosb0%Kt_b +  o_eosb0%K_pot_b
       !       
       !    write(66,*)o_eosb0%k_b(4),d2e2a_b0_den2(4),o_eosb0%d1e2a_b_den(4),i_eosDen%den_nuc
       !
       ! -------------------
       ! sound velocity square (only for the baryon component)
       ! -------------------
       !
       o_eosb0%cs_b = o_eosb0%K_b / 9._pr / &
               &( CST_mnuc2 + o_eosb0%e2a_b + o_eosb0%p_b / i_eosDen%den_nuc )
       !
       ! -------------------
       ! Pressure derivatives
       ! -------------------
       !
!       do n = 2, i_coef%nmax
!          o_eosb0%a_d1_P_b_den(n) = 2.0 * i_eosDen%den_nuc * o_eosb0%a_d1_e2a_b_den(n) + i_eosDen%den_nuc**2 * d2_e2a_b_den(n)
!       enddo
!       do n = 1, i_coef%nmax
!          o_eosb0%a_d1_P_b_xe(n) =  i_eosDen%den_nuc**2 * d2_e2a_b_den_xe(n)
!       enddo
       !
       ! -------------------
       ! Third order derivative
       ! -------------------
       !
       !write(*,*)i_coef%nsat,1.0 + 3 * i_eosDen%xx,i_eosDen%den_nuc
!       d3_t2a_b_den =  i_coef%TFGsat / ( 27.0 * i_coef%nsat**3 ) * ( 1.0 + 3 * i_eosDen%xx )**m73 * &
!            & ( 4.0 * ffg - 5.0 * ( 1.0 + 3 * i_eosDen%xx ) * fsky )
!       o_eosb0%a_d3_e2a_b_den(:) = 0
!       o_eosb0%a_d3_e2a_b_den(2) = d3_t2a_b_den
!       sum = 0
!       do n = 3, i_coef%nmax
!          tmp = (i_coef%VCOEF(0,n) + i_coef%VCOEF(1,n) * i_eosDen%xd **2) * i_eosDen%xx**(n-3) / i_coef%facto(n-3)
!          sum = sum + tmp
!          o_eosb0%a_d3_e2a_b_den(n) = d3_t2a_b_den + sum / ( 27.0 * i_coef%nsat**3 )
!       enddo
       !
       ! -------------------
       ! Contribution of the rearrangement
       ! terms to the densities in den
       ! -------------------
       !
       !    if (i_rea.eq.1) then
       
       !    endif
       !

    endif
       
  end subroutine metaeos_T0_baryons_R
     
     

  
  

  
  subroutine metaeos_T0_leptons(i_coef,i_eosDen,i_opt_beta,o_eosl0)

!    use acc; use cst;
!    use metaEosType;

    implicit none

    ! input variables
    type (metaEos_coef),      intent(in)  :: i_coef
    type (metaEos_densities), intent(in)  :: i_eosDen
    logical,                  intent(in)  :: i_opt_beta
    ! output variables
    type (metaEos_Leptons), intent(out)   :: o_eosl0
    !
    ! local variables
    !
    real (kind=pr)            :: xvar, xvar2, xfac
    !
    real(kind=pr) :: fen,fs,fp
    real(kind=pr) :: frhoA
    ! parameters
    real (kind=pr), parameter :: p13=1._pr/3._pr,  p23=2./3.,  p43=4./3.,  p53=5./3.
    real (kind=pr), parameter :: m13=-1./3., m23=-2./3., m43=-4./3., m53=-5./3., m73 = -7./3.

    ! *******************
    ! Electrons at T = 0
    ! *******************
    !
    ! Electronic chemical potential at T=0 (MeV)
    !
    o_eosl0%mu_e = dsqrt ( CST_mec2**2 + ( CST_hbc * i_eosDen%kf_e )**2 )
    o_eosl0%nu_e = CST_hbc * i_eosDen%kf_e
!    o_eosl0%nu_e = o_eosl0%mu_e - CST_mec2
    if (.not.i_opt_beta) then
       ! temporary variables
       xvar =  o_eosl0%nu_e / CST_mec2
       xvar2 = dsqrt ( 1._pr + xvar**2 )
       xfac = CST_mec2**4 / ( 8._pr * CST_pi2 * CST_hbc**3 ) 
       !
       ! energy
       !
       o_eosl0%rho_e = xfac * ( xvar * ( 1._pr + 2 * xvar**2 ) * xvar2 - asinh(xvar) )
       o_eosl0%e2v_e = o_eosl0%rho_e - i_eosDen%den_e * CST_mec2
       if (i_coef%iquark.eq.0.and.i_eosDen%den_nuc.ne.0._pr) then
          o_eosl0%e2a_e =  o_eosl0%e2v_e / i_eosDen%den_nuc 
       elseif (i_coef%iquark.ge.1.and.i_eosDen%den_b.ne.0._pr) then
          o_eosl0%e2a_e =  o_eosl0%e2v_e / i_eosDen%den_b 
       else
          o_eosl0%e2a_e = 0._pr
       endif
       !
       ! Entropy
       !
       o_eosl0%s2a_e = 0.0_pr
       !
       ! Pressure
       !
       o_eosl0%p_e = -o_eosl0%e2v_e + 8._pr * xfac / 3._pr * xvar**3 * xvar2 
       !
       ! Incompressibility: Ke=9dP_e/dn_b at xe fixed
       !
       if (i_coef%iquark.eq.0.and.i_eosDen%den_nuc.ne.0._pr) then
          o_eosl0%K_e = 8._pr * xfac / i_eosDen%den_nuc * &
               & xvar**3 * ( 4._pr * xvar**2 + 3._pr) / xvar2 &
               & - 9._pr / i_eosDen%den_nuc * ( o_eosl0%e2v_e + o_eosl0%p_e )
       elseif (i_coef%iquark.ge.1.and.i_eosDen%den_b.ne.0._pr) then
          o_eosl0%K_e = 8._pr * xfac / i_eosDen%den_b * &
               & xvar**3 * ( 4._pr * xvar**2 + 3._pr) / xvar2 &
               & - 9._pr / i_eosDen%den_b * ( o_eosl0%e2v_e + o_eosl0%p_e )
       else
          o_eosl0%K_e = 0._pr
       endif
       !
       ! UR limit
       !
       ! Energy per baryon (UR)
       !
       if (i_coef%iquark.eq.0.and.i_eosDen%den_nuc.ne.0._pr) then
          o_eosl0%e2a_e_UR = 3._pr/4._pr *  o_eosl0%nu_e * i_eosDen%xe
       elseif (i_coef%iquark.ge.1.and.i_eosDen%den_b.ne.0._pr) then
          o_eosl0%e2a_e_UR = 3._pr/4._pr *  o_eosl0%nu_e * i_eosDen%den_e / i_eosDen%den_b
       else
          o_eosl0%e2a_e_UR = 0._pr
       endif
       !
       ! Pressure (UR)
       !
!       o_eosl0%P_e_UR = 1._pr/4._pr * CST_hbc* i_eosDen%kf_e * i_eosDen%den_e
       o_eosl0%P_e_UR = o_eosl0%nu_e * i_eosDen%den_e / 4._pr
       !
       ! Incompressibility (UR)
       !
       if (i_coef%iquark.eq.0.and.i_eosDen%den_nuc.ne.0._pr) then
          o_eosl0%K_e_UR = 3._pr * i_eosDen%xe * o_eosl0%nu_e
       elseif (i_coef%iquark.ge.1.and.i_eosDen%den_b.ne.0._pr) then
          o_eosl0%K_e_UR = 3._pr * i_eosDen%den_e / i_eosDen%den_b * o_eosl0%nu_e
       else
          o_eosl0%K_e_UR = 0._pr
       endif
       !
       ! Derivatives d / d n_b
       !
       if (i_coef%iquark.eq.0.and.i_eosDen%den_nuc.ne.0._pr) then
          o_eosl0%d1_e2a_e_den = 1._pr/4._pr * i_eosDen%den_e / i_eosDen%den_nuc**2 * o_eosl0%mu_e
          o_eosl0%d1_mu_e_den  = p13 * o_eosl0%mu_e / i_eosDen%den_nuc
          o_eosl0%d1_p_e_den   = p13 * i_eosDen%den_e / i_eosDen%den_nuc * o_eosl0%mu_e
       elseif (i_coef%iquark.ge.1.and.i_eosDen%den_b.ne.0._pr) then
          o_eosl0%d1_e2a_e_den = 1._pr/4._pr * i_eosDen%den_e / i_eosDen%den_b**2 * o_eosl0%mu_e
          o_eosl0%d1_mu_e_den  = p13 * o_eosl0%mu_e / i_eosDen%den_b
          o_eosl0%d1_p_e_den   = p13 * i_eosDen%den_e / i_eosDen%den_b * o_eosl0%mu_e
       else
          o_eosl0%d1_e2a_e_den = 0._pr
          o_eosl0%d1_mu_e_den  = 0._pr
          o_eosl0%d1_p_e_den   = 0._pr
       endif
       !
       ! Derivatives d / d xe
       !
       if (i_coef%iquark.eq.0.and.i_eosDen%xe.ne.0._pr.and.i_eosDen%den_nuc.ne.0._pr) then
          o_eosl0%d1_e2a_e_xe = o_eosl0%mu_e
          o_eosl0%d1_mu_e_xe  = p13 * o_eosl0%mu_e / i_eosDen%xe
          o_eosl0%d1_p_e_xe   = p13 * o_eosl0%mu_e * i_eosDen%den_nuc 
       elseif (i_coef%iquark.eq.1.and.i_eosDen%xe.ne.0._pr.and.i_eosDen%den_b.ne.0._pr) then
          o_eosl0%d1_e2a_e_xe = o_eosl0%mu_e
          o_eosl0%d1_mu_e_xe  = p13 * o_eosl0%mu_e / i_eosDen%xe
          o_eosl0%d1_p_e_xe   = p13 * o_eosl0%mu_e * i_eosDen%den_b
       else
          o_eosl0%d1_e2a_e_xe = 0._pr
          o_eosl0%d1_mu_e_xe  = 0._pr
          o_eosl0%d1_p_e_xe   = 0._pr
       endif
    endif
    !
    !  
    ! *******************
    ! Muons at T = 0
    ! *******************
    !
    !
    ! Muon chemical potential at T=0 (MeV)
    o_eosl0%mu_muon = dsqrt ( CST_mmuonc2**2 + ( CST_hbc * i_eosDen%kf_muon )**2 )
    o_eosl0%nu_muon = CST_hbc * i_eosDen%kf_muon
!    o_eosl0%nu_muon = o_eosl0%mu_muon - CST_mmuonc2
    if (.not.i_opt_beta) then
       ! temporary variables
       xvar = o_eosl0%nu_muon / CST_mmuonc2
       xvar2 = dsqrt ( 1._pr + xvar**2 )
       xfac = CST_mmuonc2**4 / ( 8._pr * CST_pi2 * CST_hbc**3 )
       !
       ! energy 
       !
!       o_eosl0%e2v_muon = xfac * ( xvar * ( 1._pr + 2 * xvar**2 ) * xvar2 - asinh(xvar) )
       o_eosl0%rho_muon = xfac * ( xvar * ( 1._pr + 2 * xvar**2 ) * xvar2 - asinh(xvar) )
       o_eosl0%e2v_muon = o_eosl0%rho_muon - i_eosDen%den_muon * CST_mmuonc2
       if (i_eosDen%den_nuc.ne.0._pr) then
          o_eosl0%e2a_muon =  o_eosl0%e2v_muon / i_eosDen%den_nuc
       elseif (i_coef%iquark.ge.1.and.i_eosDen%den_b.ne.0._pr) then
		  o_eosl0%e2a_muon =  o_eosl0%e2v_muon / i_eosDen%den_b 
       else
          o_eosl0%e2a_muon = 0._pr
       endif
       !
       ! Entropy
       !
!       write(*,*)"P(R)",o_eosl0%p_muon,"P(UR)",o_eosl0%P_muon_UR
       o_eosl0%s2a_muon = 0.0_pr
       !
       ! Pressure
       !
       o_eosl0%p_muon = -o_eosl0%e2v_muon + 8._pr * xfac / 3._pr * xvar**3 * xvar2 
       !
       ! Incompressibility: Kmuon=9dP_muon/dn_b at xmuon fixed
       !
       if (i_coef%iquark.eq.0.and.i_eosDen%den_nuc.ne.0._pr) then
          o_eosl0%K_muon = 8._pr * xfac / i_eosDen%den_nuc * &
               & xvar**3 * ( 4._pr * xvar**2 + 3._pr) / xvar2 &
               & - 9._pr / i_eosDen%den_nuc * ( o_eosl0%e2v_muon + o_eosl0%p_muon )
       elseif (i_coef%iquark.ge.1.and.i_eosDen%den_b.ne.0._pr) then
          o_eosl0%K_muon = 8._pr * xfac / i_eosDen%den_b * &
               & xvar**3 * ( 4._pr * xvar**2 + 3._pr) / xvar2 &
               & - 9._pr / i_eosDen%den_b * ( o_eosl0%e2v_muon + o_eosl0%p_muon )
       else
          o_eosl0%K_muon = 0._pr
       endif
       !
       ! UR limit
       !
       ! Energy per baryon (UR)
       !
       if (i_eosDen%den_nuc.ne.0._pr) then
          o_eosl0%e2a_muon_UR = 3._pr/4._pr * o_eosl0%nu_muon * i_eosDen%den_muon / i_eosDen%den_nuc
       elseif (i_coef%iquark.ge.1.and.i_eosDen%den_b.ne.0._pr) then
          o_eosl0%e2a_muon_UR = 3._pr/4._pr * o_eosl0%nu_muon * i_eosDen%den_muon / i_eosDen%den_b
       else
          o_eosl0%e2a_muon_UR = 0._pr
       endif
       !
       ! Pressure (UR)
       !
       o_eosl0%P_muon_UR = o_eosl0%nu_muon * i_eosDen%den_muon / 4._pr 
       !
       ! Incompressibility (UR)
       !
       if (i_coef%iquark.eq.0.and.i_eosDen%den_nuc.ne.0._pr) then
          o_eosl0%K_muon_UR = 3._pr * i_eosDen%xmuon * o_eosl0%nu_muon
       elseif (i_coef%iquark.ge.1.and.i_eosDen%den_b.ne.0._pr) then
          o_eosl0%K_muon_UR = 3._pr * o_eosl0%nu_muon * i_eosDen%den_muon / i_eosDen%den_b
       else
          o_eosl0%K_muon_UR = 0._pr
       endif
       !
       ! Derivatives d / d n_b
       !
       if (i_eosDen%xmuon.ne.0._pr.and.i_eosDen%den_nuc.ne.0._pr) then
          o_eosl0%d1_e2a_muon_den = 1._pr/4._pr * i_eosDen%den_muon / i_eosDen%den_nuc**2 * o_eosl0%mu_muon
          o_eosl0%d1_mu_muon_den  = p13 * o_eosl0%mu_muon / i_eosDen%den_nuc
          o_eosl0%d1_p_muon_den   = p13 * i_eosDen%den_muon / i_eosDen%den_nuc * o_eosl0%mu_muon
       elseif (i_coef%iquark.ge.1.and.i_eosDen%den_b.ne.0._pr) then
          o_eosl0%d1_e2a_muon_den = 1._pr/4._pr * i_eosDen%den_muon / i_eosDen%den_b**2 * o_eosl0%mu_muon
          o_eosl0%d1_mu_muon_den  = p13 * o_eosl0%mu_muon / i_eosDen%den_b
          o_eosl0%d1_p_muon_den   = p13 * i_eosDen%den_muon / i_eosDen%den_b * o_eosl0%mu_muon
       else
          o_eosl0%d1_e2a_muon_den = 0._pr
          o_eosl0%d1_mu_muon_den  = 0._pr
          o_eosl0%d1_p_muon_den   = 0._pr
       endif
       !
       ! Derivatives d / d x_muon
       !
       if (i_eosDen%xmuon.ne.0._pr.and.i_eosDen%den_nuc.ne.0._pr) then
          o_eosl0%d1_e2a_muon_xmuon = o_eosl0%mu_muon
          o_eosl0%d1_mu_muon_xmuon  = p13 * o_eosl0%mu_muon / i_eosDen%xmuon
          o_eosl0%d1_p_muon_xmuon   = p13 * o_eosl0%mu_muon * i_eosDen%den_nuc
       elseif (i_coef%iquark.eq.1.and.i_eosDen%xe.ne.0._pr.and.i_eosDen%den_b.ne.0._pr) then
          o_eosl0%d1_e2a_muon_xmuon = o_eosl0%mu_muon
          o_eosl0%d1_mu_muon_xmuon  = p13 * o_eosl0%mu_muon / i_eosDen%xmuon
          o_eosl0%d1_p_muon_xmuon   = p13 * o_eosl0%mu_muon * i_eosDen%den_b
       else
          o_eosl0%d1_e2a_muon_xmuon = 0._pr
          o_eosl0%d1_mu_muon_xmuon  = 0._pr
          o_eosl0%d1_p_muon_xmuon   = 0._pr
       endif
    endif
    !
    ! -------------------
    ! neutrinos contribution
    ! -------------------
    !
    ! Chemical potential
    o_eosl0%mu_nu = 0.0_pr    
    if (.not.i_opt_beta) then
       ! energy 
       o_eosl0%e2a_nu = 0.0_pr 
       ! Pressure
       o_eosl0%p_nu = 0.0_pr
       ! Entropy
       o_eosl0%s2a_nu = 0.0_pr
       o_eosl0%K_nu = 0.0_pr
       ! other derivatives
       o_eosl0%d1_e2a_nu_den = 0.0_pr
       o_eosl0%d1_e2a_nu_xe = 0.0_pr
       o_eosl0%d1_e2a_nu_T = 0.0_pr
       o_eosl0%d1_p_nu_den = 0.0_pr
       o_eosl0%d1_p_nu_xe = 0.0_pr
       o_eosl0%d1_p_nu_T = 0.0_pr
    endif

  end subroutine metaeos_T0_leptons


     

  
  subroutine metaeos_T0_photons(i_eosDen,i_opt_beta,o_eosp0)

!    use acc; use cst;
!    use metaEosType;

    implicit none

    ! input variables
    type (metaEos_densities), intent(in)  :: i_eosDen
    logical, intent(in)                   :: i_opt_beta
    ! output variables
    type (metaEos_Photons), intent(out)   :: o_eosp0
    !
    ! local variables
    !
    ! parameters
    real (kind=pr), parameter :: p13=1._pr/3._pr,  p23=2./3.,  p43=4./3.,  p53=5./3.
    real (kind=pr), parameter :: m13=-1./3., m23=-2./3., m43=-4./3., m53=-5./3., m73 = -7./3.
    !
    ! -------------------
    ! photons
    ! -------------------
    ! Chemical potential
    o_eosp0%mu_phot = 0.0_pr
    ! energy 
    o_eosp0%e2a_phot = 0.0_pr 
    ! Pressure
    o_eosp0%p_phot = 0.0_pr
    ! Entropy
    o_eosp0%s2a_phot = 0.0_pr
    ! other derivatives
    o_eosp0%d1_e2a_phot_den = 0.0_pr
    o_eosp0%d1_e2a_phot_T = 0.0_pr
    o_eosp0%d1_p_phot_den = 0.0_pr
    o_eosp0%d1_p_phot_T = 0.0_pr
    
  end subroutine metaeos_T0_photons


  subroutine metaeos_T0_all(i_coef,i_eosDen,i_opt_beta,i_eosb,i_eosq,i_eosl,i_eosp,o_eosa)

!    use acc; use cst;
!    use metaEosType;

    implicit none

    ! input variables
    type (metaEos_coef), intent(in)       :: i_coef
    type (metaEos_densities), intent(in)  :: i_eosDen
    logical, intent(in)                   :: i_opt_beta
    type (metaEos_Baryons), intent(in)    :: i_eosb
    type (metaEos_Q1),      intent(in)    :: i_eosq
    type (metaEos_Leptons), intent(in)    :: i_eosl
    type (metaEos_Photons), intent(in)    :: i_eosp
    ! output variables
    type (metaEos_All), intent(out)       :: o_eosa
    !
    ! local variables
    !
    real (kind=pr) :: i_xt = 0._pr
    real (kind=pr) :: chi, kappa
    ! parameters
    real (kind=pr), parameter :: p13=1._pr/3._pr,  p23=2./3.,  p43=4./3.,  p53=5./3.
    real (kind=pr), parameter :: m13=-1./3., m23=-2./3., m43=-4./3., m53=-5./3., m73 = -7./3.


    if (i_coef%iquark.eq.0) then

       o_eosa%den_b = i_eosDen%den_nuc
       
       ! Total energy-density including the rest mass energy (in MeV):
       o_eosa%rho =  i_eosb%rho_b + i_eosl%rho_e + i_eosl%rho_muon !+ i_eosl%rho_nu + i_eosp%rho_phot

       ! Total energy-density without rest mass
       o_eosa%e2v = o_eosa%rho - i_eosDen%den_e * CST_mec2 - i_eosDen%den_p * CST_mpc2 - &
            & i_eosDen%den_n * CST_mnc2 - i_eosDen%den_muon * CST_mmuonc2
       ! Total energy per baryon (MeV)
       o_eosa%e2a = o_eosa%e2v / i_eosDen%den_nuc
       
       ! Total energy per baryon (MeV)
!       o_eosa%e2a = i_eosb%e2a_b + i_eosl%e2a_e + i_eosl%e2a_muon !+&
!         i_eosl%e2a_nu + i_eosp%e2a_phot
!       o_eosa%e2v = o_eosa%e2a * i_eosDen%den_nuc

       ! Total energy including the rest mass energy (in MeV):
!       o_eosa%eps2a =  o_eosa%e2a + &
!            & i_eosDen%xe * CST_mec2 + i_eosDen%xp * CST_mpc2 + i_eosDen%xn * CST_mnc2 + i_eosDen%xmuon * CST_mmuonc2
!    o_eosa%eps2a = i_eosDen%xn * CST_mnc2 + i_eosDen%xp * CST_mpc2 + i_eosDen%xe * CST_mec2 + o_eosa%e2a
!       o_eosa%rho = o_eosa%eps2a * i_eosDen%den_nuc

       ! Total pressure (MeV.fm-3)
       o_eosa%p = i_eosb%p_b + i_eosl%p_e + i_eosl%p_muon !+ i_eosl%p_nu + i_eosp%p_phot

       ! Total entropy
       o_eosa%s2a = i_eosb%s2a_b + i_eosl%s2a_e + i_eosl%s2a_muon !+ i_eosl%s2a_nu + i_eosp%s2a_phot

       ! Toal enthalpy
       o_eosa%enthalpy =  ( o_eosa%rho + o_eosa%p ) / i_eosDen%den_nuc

       ! Total incompressibility
       o_eosa%K = i_eosb%K_b + i_eosl%K_e + i_eosl%K_muon + i_eosl%K_nu
    
    ! derivative: d e2a_tot / d den
!    o_eosa%d1_e2a_den = i_eosb%a_d1_e2a_b_den(i_coef%nmax) + i_eosl%d1_e2a_e_den !+ i_eosl%d1_e2a_muon_den + &
!         & i_eosl%d1_e2a_nu_den + i_eosp%d1_e2a_phot_den

    ! derivative: d e2a_tot / d xe
!    o_eosa%d1_e2a_xe = i_eosb%a_d1_e2a_b_xe(i_coef%nmax) + i_eosl%d1_e2a_e_xe !+ i_eosl%d1_e2a_muon_xe + &
!         & i_eosl%d1_e2a_nu_xe + i_eosp%d1_e2a_phot_xe

    ! derivative: d e2a_tot / d T
!    o_eosa%d1_e2a_T = i_eosb%d1_e2a_b_T + i_eosl%d1_e2a_e_T !+ i_eosl%d1_e2a_muon_T + &
!         & i_eosl%d1_e2a_nu_T + i_eosp%d1_e2a_phot_T

    ! derivatives: d P_tot / d den
!    o_eosa%d1_p_den = i_eosb%a_d1_p_b_den(i_coef%nmax) + i_eosl%d1_p_e_den !+&
!         i_eosl%d1_p_muon_den + i_eosl%d1_p_nu_den + i_eosp%d1_p_phot_den
 
    ! derivatives: d P_tot / d xe
!    o_eosa%d1_p_xe = i_eosb%a_d1_p_b_xe(i_coef%nmax) + i_eosl%d1_p_e_xe !+&
!         i_eosl%d1_p_muon_xe + i_eosl%d1_p_nu_xe + i_eosp%d1_p_phot_xe

    ! derivative: d P_tot / d T
!    o_eosa%d1_p_T = i_eosb%d1p_b_T + i_eosl%d1p_e_T !+ &
!         & + i_eosl%d1_p_muon_T + i_eosl%d1_p_nu_T + i_eosp%d1_p_phot_T

    ! Gamma ! to check
!    o_eosa%gamma = i_eosDen%den_nuc / o_eosa%p * o_eosa%d1_p_den
!    if (i_xt .gt. 0._pr) o_eosa%gamma = o_eosa%gamma + i_xt * o_eosa%d1_p_T**2 &
!         & / ( i_eosDen%den_nuc * o_eosa%p * o_eosa%d1_e2a_T )

       ! Sound velocity
       if (i_eosDen%den_nuc.ne.0._pr) then
          chi = o_eosa%K / 9.0_pr
          kappa = 0._pr
!       chi  = i_eosDen%den_nuc * o_eosa%d1_p_den !- i_eosDen%xe * o_eosa%d1_p_xe
!       o_eosa%cs = chi / o_eosa%enthalpy
       ! at finite T
!       if (i_xt .gt. 0._pr) then
!          chi = chi + o_eosa%d1_p_T / o_eosa%d1_e2a_T * ( - o_eosa%d1_e2a_den )!+(1-i_xd)/(2*i_eosDen%den_nuc)*o_eosa%d1_e2a_xe)
!          kappa = o_eosa%d1_p_T / o_eosa%d1_e2a_T
          o_eosa%cs = ( chi + o_eosa%p / i_eosDen%den_nuc**2 * kappa ) / o_eosa%enthalpy
!       endif
       else
          o_eosa%cs = 0._pr
       endif

       elseif (i_coef%iquark.ge.1) then

          o_eosa%den_b = i_eosDen%den_b

          ! Total energy-density including the rest mass energy (in MeV):
          o_eosa%rho =  i_eosq%rho_B + i_eosl%rho_e + i_eosl%rho_muon !+ i_eosl%rho_nu + i_eosp%rho_phot

          ! Total energy-density without rest mass
          o_eosa%e2v = o_eosa%rho - i_eosDen%den_e * CST_mec2 - i_eosDen%den_b * CST_mnuc2 - &
               & i_eosDen%den_muon * CST_mmuonc2
          ! Total energy per baryon (MeV)
          o_eosa%e2a = o_eosa%e2v / i_eosDen%den_b
      
          ! Total pressure (MeV.fm-3)
          o_eosa%p = i_eosq%p_B + i_eosl%p_e + i_eosl%p_muon !+ i_eosl%p_nu + i_eosp%p_phot
          !write(*,*)i_eosb%p_B

          ! Total entropy
          !o_eosa%s2a = i_eosb%s2a_b + i_eosl%s2a_e + i_eosl%s2a_muon !+ i_eosl%s2a_nu + i_eosp%s2a_phot

          ! Toal enthalpy
          o_eosa%enthalpy =  ( o_eosa%rho + o_eosa%p ) / i_eosDen%den_b

          ! Total incompressibility
          o_eosa%K = i_eosq%K_b + i_eosl%K_e + i_eosl%K_muon + i_eosl%K_nu
          if (i_eosDen%den_nuc.ne.0._pr) then
             chi = o_eosa%K / 9.0_pr
             kappa = 0._pr
             o_eosa%cs = ( chi + o_eosa%p / i_eosDen%den_b**2 * kappa ) / o_eosa%enthalpy
          else
             o_eosa%cs = 0._pr
          endif
       
       endif
       
  end subroutine metaeos_T0_all
  


  

  ! Find the beta equilibrium with electrons and muons
  ! using dichotomy method
  ! o_xd has to be guessed initially
  ! Define an intervalle xd1 < xd2
  ! Define xdm = (xd1 + xd2)/2
  ! decide to take [xd1,xdm] or [xdm,xd2]

  ! note: at T=0, we have mu_mu**2 = m_mu**2 + nu_mu**2 = mu_e**2
  ! then nu_mu**2 = mu_e**2 - m_mu**2
  ! and nu_mu = hbc * kf_mu
  ! and kf_mu**3 = 3pi2 rho_b x_mu
  ! then x_mu = (mu_e**2 - m_mu**2)**3/2 / (3pi2 rho_b)

  ! The value of xd is obtained from the chemical potential equation
  ! The value of xmuon is obtained from mu_e=mu_muon
  
  subroutine metaeos_T0_beta_npemuon(i_den,i_xd,i_xmuon,i_coef,i_withq,o_eosDen,ibeta)

! The beta equilibrium at finite T is not implemented
! Fermi momentum of e at finite T shall be calculated

!    use acc; use CST;
!    use metaEosType;

    implicit none

    ! input variables
    real (kind=pr),      intent(in)  :: i_den, i_xd, i_xmuon
    type (metaEos_coef), intent(in)  :: i_coef
    logical,             intent(in)  :: i_withq
    ! output variables
    type (metaEos_Densities), intent(out)  :: o_eosDen
    integer,                  intent(out)  :: ibeta
    ! local variables
    type (metaEos_Baryons)   :: eosb
    type (metaEos_Q1)        :: eosq
    type (metaEos_Leptons)   :: eosl
    type (metaEos_Densities) :: eosDen
    integer             :: i
    integer             :: imax = 100
    real (kind=pr)      :: i_xt
    real (kind=pr)      :: xd
    real (kind=pr)      :: xmuon, xmuon0, dxmuon
    real (kind=pr)      :: xmuon1, xmuon2, xmuonm
    real (kind=pr)      :: xd1, xd2, xdm, dxd, xtt
    real (kind=pr)      :: f1, f2, fm
    real (kind=pr)      :: g1, g2, gm
    real (kind=pr)      :: den_muon, mu_muon
!    write(*,*)"init",i_xd,i_xmuon
    ! T=0
    i_xt = 0._pr
    ibeta = 1
    !
    ! Solution with xd1 min
    !
    xd1 = max(i_xd-0.5_pr,-0.99999)
    dxmuon = 1._pr
    xmuon0 = 1._pr
    xmuon = i_xmuon ! 0.1 * xd1
    ! Search xmuon for xd1 fixed imposing mu_muon=mu_e
    do while (dabs(dxmuon).gt.1.d-12)
       call metaeos_T0_set_densities(i_den,xd1,xmuon,i_coef,i_withq,eosDen)
       call metaeos_T0_leptons(i_coef,eosDen,.true.,eosl)
      ! condition at T=0
       if (i_xt.eq.0.0.and.eosl%mu_e.gt.CST_mmuonc2) then 
          xmuon = ( ( eosl%mu_e**2 - CST_mmuonc2**2 ) / CST_hbc**2 )**1.5 / &
               & ( 3._pr * CST_pi2 * eosDen%den_nuc ) 
          ! condition at T!=0
       else if (i_xt.ne.0.0.and.eosl%mu_e.gt.(CST_mmuonc2-5*i_xt)) then
          den_muon = xmuon * eosDen%den_nuc
          mu_muon = 0.5*(eosl%mu_e + eosl%mu_muon)
          !den_muon = fgas_R_rho(i_xt,den_muon,mu_muon,CST_mmuonc2)
          xmuon = den_muon / eosDen%den_nuc
       else
          xmuon = 0._pr
       endif
!       write(*,*)"sol1:",i_den,eosDen%den_nuc,xmuon,eosl%mu_e,CST_mmuonc2
!       call metaeos_T0_set_densities(i_den,xd1,xmuon,i_coef,eosDen)
!       call metaeos_T0_leptons(i_coef,eosDen,.true.,eosl)
       dxmuon = xmuon - xmuon0
       xmuon0 = xmuon
    enddo
    xmuon1 = xmuon
    call metaeos_T0_baryons(i_coef,eosDen,.true.,i_withq,eosb,eosq)
    if (i_coef%iquark.eq.0) f1  = eosl%mu_e - ( CST_mnc2 - CST_mpc2 + eosb%mu_np)
    if (i_coef%iquark.ge.1) f1  = eosl%mu_e - ( CST_mnc2 - CST_mpc2 + eosq%munp)
!    write(*,*)"sol1:",xd1,xmuon1,f1
    !
    ! Solution with xd2 max
    !
    xd2 = min(i_xd+0.5_pr,0.99999)
    dxmuon = 1._pr
    xmuon0 = 1._pr
    xmuon = i_xmuon ! 0.1 * xd2
    ! Search xmuon for xd2 fixed imposing mu_muon=mu_e
    do while (dabs(dxmuon).gt.1.d-12)
       call metaeos_T0_set_densities(i_den,xd2,xmuon,i_coef,i_withq,eosDen)
       call metaeos_T0_leptons(i_coef,eosDen,.true.,eosl)
       ! condition at T=0
       if (i_xt.eq.0.0.and.eosl%mu_e.gt.CST_mmuonc2) then 
          xmuon = ( ( eosl%mu_e**2 - CST_mmuonc2**2 ) / CST_hbc**2 )**1.5 / &
               & ( 3._pr * CST_pi2 * eosDen%den_nuc ) 
          ! condition at T!=0
       else if (i_xt.ne.0.0.and.eosl%mu_e.gt.(CST_mmuonc2-5*i_xt)) then
          den_muon = xmuon * i_den
          mu_muon = 0.5 * ( eosl%mu_e + eosl%mu_muon )
          !den_muon = fgas_R_rho(i_xt,den_muon,mu_muon,CST_mmuonc2)
          xmuon = den_muon / i_den
       else
          xmuon = 0._pr
       endif
!       call eos_T0_set_densities(i_den,xd2,xmuon,i_coef,eosDen)
!       call eos_T0_leptons(i_coef,eosDen,.true.,eosl)
       dxmuon = xmuon - xmuon0
       xmuon0 = xmuon
    enddo
    xmuon2 = xmuon
    call metaeos_T0_baryons(i_coef,eosDen,.false.,i_withq,eosb,eosq)
!    f2  = eosl%mu_e - ( CST_mnc2 - CST_mpc2 + eosb%mu_np)
    if (i_coef%iquark.eq.0) f2  = eosl%mu_e - ( CST_mnc2 - CST_mpc2 + eosb%mu_np)
    if (i_coef%iquark.ge.1) f2  = eosl%mu_e - ( CST_mnc2 - CST_mpc2 + eosq%munp)
!    write(*,*)"sol2:",xd2,xmuon2,f2
!    write(*,*)"esym ",eosb%esym_b
!    write(*,*)"test:",eosDen%den_nuc,f1,f2,xd1,xd2,eosb%mu_np,eosq%munp
!    if (i_coef%iquark.eq.0.and.eosb%mu_np.lt.0._pr) then
!       eosDen%xd = 1._pr
!       eosDen%xmuon = 0._pr
!    elseif (i_coef%iquark.eq.1.and.eosq%munp.lt.0._pr) then
!       eosDen%xd = 1._pr
!       eosDen%xmuon = 0._pr
!    write(*,*)"eq-b:",f1,f2
    if (f1*f2.gt.0._pr) then
       eosDen%xd = 1._pr
       eosDen%xmuon = 0._pr
       ibeta = 0
    else
       !
       ! Calculate xd by dichotomy:
       ! first define the boundaries xd1 and xd2
       ! (done before), then fix the mean value xdm
       ! Finally decide where is the zero in the
       ! smaller intervalle [xd1,xdm] or [xdm,xd2]
       !
       dxd = 1.0
       i = 1
       xmuon = i_xmuon !0.5 * ( xmuon1 + xmuon2 )
       !
       do while (dabs(f2).gt.1.d-8)
!       do while (dabs(f2).gt.1.d-8.and.dxd.gt.0.1)
          !
          ! Calculate xmuon imposing mu_muon=mu_e
          ! for a fixed value of xd = xdm.
          ! This relation is solved self-consistently
          ! since we first assume a value for xmuon
          ! then calculate mu_e, and deduce den_muon
          ! from den_muon, we deduce the new value for xmuon
          !
          xdm = 0.5 * ( xd1 + xd2 )
          dxmuon = 1.0
          xmuon0 = xmuon
          do while (dabs(dxmuon).gt.1.d-12)
             call metaeos_T0_set_densities(i_den,xdm,xmuon,i_coef,i_withq,eosDen)
             call metaeos_T0_leptons(i_coef,eosDen,.true.,eosl)
             ! condition at T=0
             if (i_xt.eq.0.0.and.eosl%mu_e.gt.CST_mmuonc2) then
!                write(*,*)"T=0",xmuon
                xmuon = ( ( eosl%mu_e**2 - CST_mmuonc2**2 ) / CST_hbc**2 )**1.5 / &
                     & ( 3._pr * CST_pi2 * eosDen%den_nuc ) 
                ! condition at T!=0
             else if (i_xt.ne.0.0.and.eosl%mu_e.gt.(CST_mmuonc2-5*i_xt)) then
                den_muon = xmuon * i_den
                mu_muon = 0.5 * ( eosl%mu_e + eosl%mu_muon )
                !den_muon = fgas_R_rho(i_xt,den_muon,mu_muon,CST_mmuonc2)
                xmuon = den_muon / i_den
             else
                xmuon = 0._pr
             endif
!             call metaeos_T0_set_densities(i_den,xdm,xmuon,i_coef,i_withq,eosDen)
             !          call eos_T0_set_densities(i_den,xdm,xmuon,i_coef,eosDen)
             !          call eos_T0_leptons(i_coef,eosDen,.true.,eosl)
             dxmuon = xmuon - xmuon0
             xmuon0 = xmuon
          enddo
          call metaeos_T0_baryons(i_coef,eosDen,.true.,i_withq,eosb,eosq)
!          fm  = eosl%mu_e - ( CST_mnc2 - CST_mpc2 + eosb%mu_np)
          if (i_coef%iquark.eq.0) fm  = eosl%mu_e - ( CST_mnc2 - CST_mpc2 + eosb%mu_np)
          if (i_coef%iquark.ge.1) fm  = eosl%mu_e - ( CST_mnc2 - CST_mpc2 + eosq%munp)
!          write(*,*)"    ",eosDen%den_nuc,i,fm,xdm,eosb%mu_np,eosq%munp,eosl%mu_e,eosl%mu_muon
          
          if ( f1 * fm .le. 0.0) then
             xd2 = xdm
             f2 = fm
          else if ( f2 * fm .le. 0.0) then
             xd1 = xdm
             f1 = fm
          else
             write(*,*)"!! issues min: ",eosDen%den_nuc
             write(*,*)"!! issues min: ",f1,fm,f2
             write(*,*)"!! issues min: ",xd1,xdm,xd2,dxd
             stop "!! issues with min"
          endif

          dxd = ( xd2 - xd1 )**2 + ( eosl%mu_e - eosl%mu_muon )**2
          i = i + 1
          if (i.gt.imax) then
             write(*,*)"!! issues ",i,dxd,f1,fm,f2
             stop "!! issues with the number of iterations (xd)"
          endif
      
       enddo

    endif

    o_eosDen = eosDen
    
    ! check: nu_e = nu_muon when muon are present
    if ( eosl%mu_e.gt.CST_mmuonc2 .and.dabs(eosl%mu_e-eosl%mu_muon).gt.1.0) then
       write(*,*)"!! issues nu: ",i,eosl%mu_e,eosl%mu_muon,i_den,dxd
       stop "!! issues with mu_e and mu_muon"
    endif

  end subroutine metaeos_T0_beta_npemuon




  
  subroutine metaeos_T0_beta_npe(i_den,i_xd,i_coef,i_withq,o_eosDen,ibeta)

! The beta equilibrium at finite T is not implemented
! Fermi momentum of e at finite T shall be calculated

!    use acc; use CST;
!    use metaEosType;

    implicit none

    ! input variables
    real (kind=pr),      intent(in)  :: i_den, i_xd
    type (metaEos_coef), intent(in)  :: i_coef
    logical,             intent(in)  :: i_withq
    ! output variables
    type (metaEos_Densities), intent(out)  :: o_eosDen
    integer,                  intent(out)  :: ibeta
    ! local variables
    type (metaEos_Baryons)   :: eosb
    type (metaEos_Q1)        :: eosq
    type (metaEos_Leptons)   :: eosl
    type (metaEos_Densities) :: eosDen
    integer             :: i, id
    integer             :: imax = 100
    real (kind=pr)      :: xd
    real (kind=pr)      :: xmuon, xmuon0, dxmuon
    real (kind=pr)      :: xmuon1, xmuon2, xmuonm
    real (kind=pr)      :: xd1, xd2, xdm, dxd, xtt
    real (kind=pr)      :: f1, f2, fm
    real (kind=pr)      :: g1, g2, gm
    real (kind=pr)      :: den_muon, mu_muon

    ! no muons
    xmuon = 0._pr
    ibeta = 1
    !
    ! represent f(xd):
    !
!    do id=1,10
!       xd=-1+id*0.2
!       call metaeos_T0_set_densities(i_den,xd,xmuon,i_coef,eosDen)
!       call metaeos_T0_leptons(i_coef,eosDen,.true.,eosl)
!       call metaeos_T0_baryons(i_coef,eosDen,.true.,eosb)
!       f1  = eosl%mu_e - ( CST_mnc2 - CST_mpc2 + eosb%mu_np)
!       write(*,*)i_den,xd,eosl%mu_e,eosb%mu_np,f1
!    enddo
!    stop
!    xd1 = max(i_xd-0.5_pr,-0.99999)
    xd1 = -0.99999999
!    xd1 = -0.999
!    write(*,*)"xd1",xd1
!    call eos_sim_xd(i_xx,xd1,i_xt,xmuon,i_param,o_eos)
    call metaeos_T0_set_densities(i_den,xd1,xmuon,i_coef,i_withq,eosDen)
    call metaeos_T0_leptons(i_coef,eosDen,.true.,eosl)
    call metaeos_T0_baryons(i_coef,eosDen,.true.,i_withq,eosb,eosq)
!    f1  = eosl%mu_e - ( CST_mnc2 - CST_mpc2 + eosb%mu_np)
    if (i_coef%iquark.eq.0) f1  = eosl%mu_e - ( CST_mnc2 - CST_mpc2 + eosb%mu_np)
    if (i_coef%iquark.ge.1) f1  = eosl%mu_e - ( CST_mnc2 - CST_mpc2 + eosq%munp)
!    write(*,*)'!!! f1 ',xd1,xmuon,f1,eosl%mu_e,eosb%mu_np

    !    xd2 = min(i_xd+0.5_pr,0.99999)
    xd2 = 0.99999999
!    write(*,*)"xd2",xd2
!    call eos_sim_xd(i_xx,xd2,i_xt,xmuon,i_param,o_eos)
    call metaeos_T0_set_densities(i_den,xd2,xmuon,i_coef,i_withq,eosDen)
    call metaeos_T0_leptons(i_coef,eosDen,.true.,eosl)
    call metaeos_T0_baryons(i_coef,eosDen,.false.,i_withq,eosb,eosq)
!    f2  = eosl%mu_e - ( CST_mnc2 - CST_mpc2 + eosb%mu_np)
    if (i_coef%iquark.eq.0) f2  = eosl%mu_e - ( CST_mnc2 - CST_mpc2 + eosb%mu_np)
    if (i_coef%iquark.ge.1) f2  = eosl%mu_e - ( CST_mnc2 - CST_mpc2 + eosq%munp)
!    write(*,*)'!!! f2 ',xd2,xmuon,f2,eosl%mu_e,eosb%mu_np


!    write(*,*)"Esym ",eosb%esym_b
!    if (eosb%mu_np.le.0._pr) then
       !
!       eosDen%xd = 1._pr
!       eosDen%xmuon = 0._pr
       !
    if (f1*f2.gt.0._pr) then
       eosDen%xd = 1._pr
       eosDen%xmuon = 0._pr
       ibeta = 0
    else
       !
       dxd = 1.0
       i = 1
       !
       ! Calculate xd by dichotomy:
       ! first define the boundaries xd1 and xd2
       ! (done above), then fix the mean value xdm
       ! Finally decide where is the zero in the
       ! smaller intervalle xd1-xdm or xdm-xd2
       !
       do while (dabs(f2).gt.1.d-8)
          !
          ! Calculate xmuon imposing mu_muon=mu_e
          ! for a fixed value of xd = xdm.
          ! This relation is solved self-consistently
          ! since we first assume a value for xmuon
          ! then calculate mu_e, and deduce den_muon
          ! from den_muon, we deduce the new value for xmuon
          !
          xdm = 0.5 * ( xd1 + xd2 )
          dxmuon = 1.0
          xmuon0 = 0._pr
          !       call eos_sim_xd(i_xx,xdm,i_xt,xmuon,i_param,o_eos)
          call metaeos_T0_set_densities(i_den,xdm,xmuon,i_coef,i_withq,eosDen)
          call metaeos_T0_leptons(i_coef,eosDen,.true.,eosl)
          call metaeos_T0_baryons(i_coef,eosDen,.true.,i_withq,eosb,eosq)
!          fm  = eosl%mu_e - ( CST_mnc2 - CST_mpc2 + eosb%mu_np)
          if (i_coef%iquark.eq.0) fm  = eosl%mu_e - ( CST_mnc2 - CST_mpc2 + eosb%mu_np)
          if (i_coef%iquark.ge.1) fm  = eosl%mu_e - ( CST_mnc2 - CST_mpc2 + eosq%munp)
!          write(*,*)'!!! fm ',xdm,xmuon,fm,eosl%mu_e,eosb%mu_np,dxd
          
          if ( f1 * fm .le. 0.0) then
             xd2 = xdm
             f2 = fm
          else if ( f2 * fm .le. 0.0) then
             xd1 = xdm
             f1 = fm
          else
             write(*,*)"!! issues min: no zero in the two segments"
             write(*,*)"!! issues min: den ",i_den
             write(*,*)"!! issues min: f1,fm,f2 ",f1,fm,f2
             write(*,*)"!! issues min: xd1,xdm,xd2,dxd ",xd1,xdm,xd2,dxd
             stop "!! issues with min"
          endif
          
          dxd = ( xd2 - xd1 )**2 + ( eosl%mu_e - eosl%mu_muon )**2
          i = i + 1
          !       write(*,*)xdm, dxd, f2
          if (i.gt.imax) then
             write(*,*)"!! issues ",i,dxd,f2,xdm,fm,eosb%mu_np,eosl%mu_e
             stop "!! issues with the number of iterations (xd)"
          endif
      
       enddo

    endif

    o_eosDen = eosDen

  end subroutine metaeos_T0_beta_npe








end module metaEosT0
