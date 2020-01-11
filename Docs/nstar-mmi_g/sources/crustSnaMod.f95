! Fortran95 line length stop there --->                                                                                          132
!
! =============================================================
!
!                        MODULE  crustSnaMod
!
! =============================================================
!
!
! --------------------------NOTES------------------------------
!
!  Build the EOS in the crust.
!  for the TOV solver.
!             nb   baryon number density in fm-3
!             rho  energy-density in g/cm3
!             pre  pressure in dynes/cm2
!
! -------------------------------------------------------------
! Code by Jerome Margueron
! Date: 15-09-2009
! Latest revision - 10-05-2011
! -------------------------------------------------------------

MODULE crustSnaMod

  use acc; use cst;
  use crustSnaType;
  use metaEosT0Type

  TYPE (crustSna_inputs), private :: fcn_crustSna_inp
  real (kind=pr),         private :: fcn_nb
  !
  TYPE (crustSna_crust), private  :: fcn_crust
  integer , private               :: iloop

  TYPE (metaEos_Coef)             :: fcn_coef_meta


!  IMPLICIT none

contains


! Instead of reading an input file, the data are hard coded
! This fasten the code
!
subroutine crustSna_read_inputs(crustSna_inp)

  use metaEosT0Type;

  IMPLICIT None
  !
  ! =====================
  ! DEFINE VARIABLES
  ! =====================
  !
  ! ---------------------
  ! OUTPUT VARIABLES
  ! ---------------------
  !
  TYPE (crustSna_inputs), intent(out) :: crustSna_inp
  !
  ! ---------------------
  ! LOCAL VARIABLES
  ! ---------------------
  !
  ! .....................
  ! get argument
  ! .....................
  !
  integer :: num_args
  character(len=20), dimension(:), allocatable :: args
  !
  ! .....................
  ! read input file
  ! .....................
  !
  character (len=20) :: model
  real (kind=pr)     :: nb, temp
  real (kind=pr)     :: Ye_init, ncl_init, I_init
  integer            :: dobeta, doloopnb
  real (kind=pr)     :: den_crust_min, den_crust_max
  real (kind=pr)     :: ac, apair
  real (kind=pr)     :: sigSurf, bSurf, pSurf
  namelist /crustPar/model, nb, temp, dobeta, Ye_init, ncl_init, I_init, doloopnb, den_crust_min, den_crust_max, ac, sigSurf, bSurf, pSurf
  integer            :: ix
  character (len=20) :: ifile_crustSna
  !
  ! =====================
  ! START ROUTINE
  ! =====================
  !
  ! Read arguments of the code
  !
  num_args = command_argument_count()
  !
  ! Select option according to the arguments
  !   If no argument, then take default values
  !   If argument, then read the parameter file
  !
  if (num_args.eq.0.or.num_args.eq.metaEos_neparam) then
!  if (crust_lfast) then
     !
     !  crustSna_inp%model     = "sly_nd.dat" ! choose between "sly.dat" and "fps.dat" for the crust EoS
     !  crustSna_inp%model     = "sly_lorene_sd.dat" ! choose between "sly.dat" and "fps.dat" for the crust EoS
!     crustSna_inp%model     = "fps.dat" ! choose between "sly.dat" and "fps.dat" for the crust EoS
!     crustSna_inp%model     = "fps.dat" ! choose between "sly.dat" and "fps.dat" for the crust EoS
!     crustSna_inp%model     = "sly.dat" ! choose between "sly.dat" and "fps.dat" for the crust EoS
     crustSna_inp%model  = "meos"    ! choose between "same" (=sly or fps) or "meos" for meta-EOS
     crustSna_inp%nb        = 0.001      ! baryon density in fm-3
     crustSna_inp%temp      = 0.d0      ! temperature (MeV)
     crustSna_inp%dobeta    =  0
     crustSna_inp%Ye_init    =  0.2
     crustSna_inp%ncl_init    =  0.155
     crustSna_inp%I_init    =  0.4
     crustSna_inp%doloopnb    =  0
     crustSna_inp%den_crust_min  = 1.d-4      ! min density for the empirical eos in unit of den0
     crustSna_inp%den_crust_max  = 0.01     ! max density for the empirical eos in unit of den0
     crustSna_inp%ac = 0.7
     crustSna_inp%sigSurf = 1.1
     crustSna_inp%bsurf = 10.0
     crustSna_inp%pSurf = 3.0
     crustSna_inp%apair = 12.0

  else if (num_args.eq.1.or.num_args.eq.2) then
     !
     allocate(args(num_args))
     do ix = 1, num_args
        call get_command_argument(ix,args(ix))
     end do
     !
     ! args contains the name of the input files:
     !    - args(1): crustSna.in
     !    - args(2): tov.in
     !
     read(args(1),'(a20)')ifile_crustSna
     !
     write(*,*)" crustSna input file :",ifile_crustSna
     open(unit=10,file=trim(adjustl(ifile_crustSna)),status="unknown")
     read( unit = 10, nml = crustPar )

     fcn_crustSna_inp%model    = model
     crustSna_inp%nb       = nb
     crustSna_inp%temp      = temp
     crustSna_inp%dobeta    =  dobeta
     crustSna_inp%Ye_init    =  Ye_init
     crustSna_inp%ncl_init    =  ncl_init
     crustSna_inp%I_init    =  I_init
     crustSna_inp%doloopnb    =  doloopnb
     crustSna_inp%den_crust_min  = den_crust_min
     crustSna_inp%den_crust_max  = den_crust_max
     crustSna_inp%ac = ac
     crustSna_inp%sigSurf = sigSurf
     crustSna_inp%bSurf = bSurf
     crustSna_inp%pSurf = pSurf
     crustSna_inp%apair = apair
     !
     close(10)
     !
  else
     !
     write(*,*)"setup_input: issue with arguments of the code"
     write(*,*)"Two options:"
     write(*,*)"   - No argument : take default arguments defined in setup_input"
     write(*,*)"   - 1 or 2 argument(s) : name of the input file(s) (less than 20 characters)"
     write(*,*)"Exemples:"
     write(*,*)"   - ./crustSna.e "
     write(*,*)"   - ./crustSna.e crustSna.in"
     write(*,*)"   - ./crustSna.e crustSna.in tov.in"
     write(*,*)"Code stopped"
     stop
     !
  end if

  if (crustSna_lverb) write(*,'(T4,a26,a20)') "crustSna_inp%model          = ",crustSna_inp%model
  if (crustSna_lverb) write(*,'(T4,a26,i2)')  "crustSna_inp%dobeta         = ",crustSna_inp%dobeta
  if (crustSna_lverb) write(*,'(T4,a26,f7.2)')"crustSna_inp%nb             = ",crustSna_inp%nb
  if (crustSna_lverb) write(*,'(T4,a26,f7.2)')"crustSna_inp%temp           = ",crustSna_inp%temp
  if (crustSna_lverb) write(*,'(T4,a26,f7.2)')"crustSna_inp%Ye_init        = ",crustSna_inp%Ye_init
  if (crustSna_lverb) write(*,'(T4,a26,f7.2)')"crustSna_inp%ncl_init       = ",crustSna_inp%ncl_init
  if (crustSna_lverb) write(*,'(T4,a26,f7.2)')"crustSna_inp%I_init         = ",crustSna_inp%I_init
  if (crustSna_lverb) write(*,'(T4,a26,f7.2)')"crustSna_inp%den_crust_min  = ",crustSna_inp%den_crust_min
  if (crustSna_lverb) write(*,'(T4,a26,f7.2)')"crustSna_inp%den_crust_max  = ",crustSna_inp%den_crust_max


  fcn_crustSna_inp = crustSna_inp
  !
  ! =====================
  ! END ROUTINE
  ! =====================
  !
end subroutine crustSna_read_inputs










! There are two possibilities for the crust EOS which is matched to
! Either the sly EOS or the FPS EOS
! The choice is determined by the value of inp%eos_crust

subroutine crustSna_compute_eos(coef_meta,coef_poly,crustSna_inp,withq,crust_eos)

  use metaEosT0Type;
  use polyEosT0Type;


  IMPLICIT None
  !
  ! =====================
  ! DEFINE VARIABLES
  ! =====================
  !
  ! ---------------------
  ! INPUT VARIABLES
  ! ---------------------
  !
  TYPE (metaEos_coef), intent(in)   :: coef_meta
  TYPE (polyEos_coef), intent(in)   :: coef_poly
  TYPE (crustSna_inputs), intent(in)   :: crustSna_inp
  logical,             intent(in)   :: withq
  !

  ! ---------------------
  ! OUTPUT VARIABLES
  ! ---------------------
  !
  TYPE (crustSna_eos_outputs),    intent(out)  :: crust_eos
  !
  !
  ! =====================
  ! START ROUTINE
  ! =====================
  !
  select case (trim(crustSna_inp%model))
     !
  case ("meos")
     !
     ! Set up the crust EOS from the SNA model
      !
    fcn_coef_meta = coef_meta
     call crustSna_compute_eos_sna(crustSna_inp,crust_eos)
     !
  case default
     !
     ! Or read the crust EoS from file
     !
     call crustSna_read_file(crustSna_inp,crust_eos)
     !
  end select
  !
  ! =====================
  ! END ROUTINE
  ! =====================
  !
end subroutine crustSna_compute_eos

!




! The array crustEos is defined as:
!   -crust%EosLog(1,i): baryon density (fm-3)
!   -crust%EosLog(2,i): total energy density (g cm-3)
!   -crust%EosLog(3,i): total pressure (dyn cm-2)
!   -crust%EosLog(4,i): square of the total sound velocity (in c2)



subroutine crustSna_compute_eos_sna(crustSna_inp,crust_eos)

  use spl;

  IMPLICIT None
  !
  ! =====================
  ! DEFINE VARIABLES
  ! =====================
  !
  ! ---------------------
  ! INPUT VARIABLES
  ! ---------------------
  !
  TYPE (crustSna_inputs), intent(in) :: crustSna_inp
  !
  ! ---------------------
  ! OUTPUT VARIABLES
  ! ---------------------
  !
  TYPE (crustSna_eos_outputs), intent(out) :: crust_eos
  !
  ! ---------------------
  ! LOCAL VARIABLES
  ! ---------------------
  !
!  TYPE (crustSna_crust) :: crust
  !
  integer        :: iden
!  real (kind=pr) :: den
  !
! loop variables
  !
  integer        :: NPOINT,ii
  real (kind=pr) :: den,DNB
  ! =====================
  ! START ROUTINE
  ! =====================
  !
  if (crustSna_lverb) write(*,*)" Crust EOS:"
  !
  ! compute the crust EOS here !!
  ! make a loop over the density or calculate at one 1 density
  if (fcn_crustSna_inp%doloopnb.eq.1) then
     den = fcn_crustSna_inp%den_crust_min
  else
     den = fcn_crustSna_inp%nb
  endif
  !
  iden = 1
  !
  if ((den.le.fcn_crustSna_inp%den_crust_max) .or. (fcn_crustSna_inp%doloopnb.eq.0)) then
   fcn_nb = den
     call crustSna_compute_sna(fcn_crustSna_inp,fcn_crust)
     ! store results here
     iden = iden + 1
      write(*,*)'loop',iden,den,fcn_crustSna_inp%doloopnb
   end if

  !
  ! =====================
  ! END ROUTINE
  ! =====================
  !

end subroutine crustSna_compute_eos_sna




subroutine crustSna_compute_sna(crustSna_inp,crust)

  use spl;  USE NLE; !USE crust_fcn; !USE NLE_FCN;
  USE metaEosT0;  USE metaEosT0Type;
  use crustSnaType

  IMPLICIT None
  !
  ! =====================
  ! DEFINE VARIABLES
  ! =====================
  !
  ! ---------------------
  ! INPUT VARIABLES
  ! ---------------------
  !
  TYPE (crustSna_inputs), intent(in) :: crustSna_inp
  !
  ! ---------------------
  ! OUTPUT VARIABLES
  ! ---------------------
  !
 TYPE (crustSna_crust), intent(out) :: crust
  !
  ! ---------------------
  ! LOCAL VARIABLES
  ! ---------------------
  integer                                    :: n, ifail
  real (kind=pr), dimension(:), allocatable  :: x, y, stp,fvec, v
  real (kind=pr)                             :: eps, nor,x1,x2,x3
  LOGICAL :: check
  ! local variable
  integer  :: i
!  real (kind=pr)                             ::

 !
  ! =====================
  ! START ROUTINE
  ! =====================
 ! if Ye fixed then
!===================================================================
  ! Solve 2 non-linear equations with nle.f90 solver

  ! set the dimension:
  n = 2
  ! set the accuracy:
  eps = 1.e-10
  ! allocate the vectors:
  allocate(x(n)); allocate(y(n)); allocate(stp(n));

!********************************************************************
!         INITIAL GUESS
!********************************************************************
  x(1) = fcn_crustSna_inp%ncl_init                  !n_cl
  x(2) = fcn_crustSna_inp%I_init                    !I_cl
 ! x(3) = fcn_crustSna_inp%Ye_init
  ! initialisation of the step and first guess of the solution:
  do i=1,n
     stp(i)=0.1_pr * x(i)
  enddo

 write(*,'(1x,a13,t15,3(f14.10,2x))')' Init. guess:',x(:)
!********************************************************************
!         CALL Non linear equations SOLVER
!********************************************************************

   iloop = 1
!  call nles(fcn4,x,stp,eps,n,nor,ifail)
   call nles(crust_sna,x,stp,eps,n,nor,ifail)
  fcn_crust%ncl = x(1)
  fcn_crust%xi  = x(2)

  ! write the solution:
!   write(*,'(1x,a13,t15,3(f14.10,2x))')'a,ng,ne,nb', fcn_crust%xa,fcn_crust%ng,fcn_crust%ne,fcn_crustSna_inp%nb
  if (ifail.eq.0) write(*,*)' Solution nle:',x(:)
  write(*,'(1x,a10,t15,d14.6)')' Accuracy:',nor
  ! deallocate vectors
  deallocate(x,y,stp)
!========================================================================


  if (crustSna_lverb) write(*,*)" Crust:"
  !
  ! compute the crust EOS here !!
  !
  !
  ! =====================
  ! END ROUTINE
  ! =====================
  !

end subroutine crustSna_compute_sna




! =============================================================
! 	             SUBROUTINE crust_sna
!
!    to be linked with nle.f90
!
!    solve system of 2 non linear eqs if Ye fixed
!    solve system of 3 non linear eqs if beta equilibrium
! =============================================================
! --------------------------NOTES------------------------------
!  The soubroutine 'crust_sna' takes nb from crustSna_compute_eos_sna suroutine
!  and the initial guess to ncl, ne and xi
!
!  Construct           ng = ng (den,ncl,ne,xi)
!  and mass number     xa = xa (ncl,ne,xi)
!
!  Defines functions to be minimized in terms of the
!  chemical potentials and pressures
!
!  Find the root of the functions using nle.f90 solver and
!  gives as output the variables (A,xi,ne,n_cl,ng)
!
!********************************************************************
 subroutine crust_sna(x,y,n)

  use acc; USE CST; USE crustSnaType;
  USE metaEosT0;  USE metaEosT0Type;

  implicit none

  ! ---------------------
  ! INPUT VARIABLES
  ! ---------------------

  !
! ---------------------
  ! LOCAL VARIABLES
  ! ---------------------

  real (kind=pr), dimension(n), intent(in)  :: x
  integer, intent(in)                       :: n
  real (kind=pr), dimension(n), intent(out) :: y

  TYPE (metaEos_Densities)     :: eosDen_el, eosDen_cl, eosDen_gas
  TYPE (metaEos_Baryons)       :: eosb0_cl, eosb0_gas
  TYPE (metaEos_Leptons)       :: eosl0


!  real (kind=pr) :: ncl, xi
  real (kind=pr) :: xg, xd
  real (kind=pr) :: Z_cl, u, f_u
  real (kind=pr) :: r_cl, Yp, sigma, d_sigma_d_I
  real (kind=pr) :: m_p, m_n, pi
  real (kind=pr) :: E_cl, E_surf, E_coul, ener_unif_cl
  real (kind=pr) :: d_E_unif_d_I, d_E_cl_d_I
  real (kind=pr) :: o3ncl23, denom

! declare pressures and chemical potentials
   real (kind=pr) :: mu_el, mu_el_unif, P_gas, P_cl, P_unif_cl

!     parameters from module 'cst'
	m_p  = CST_mpc2
	m_n  = CST_mnc2
	pi   = CST_pi

! reads initial guess
   fcn_crust%ncl = x(1)
   fcn_crust%xi  = x(2)
!  fcn_crust%Ye = x(3)
    fcn_crust%Ye = fcn_crustSna_inp%Ye_init

!   write(*,'(a9,10f10.4)')'xxx(:)',x(:)
!   write(*,*)'yyy(:)',y(:)
!   if(iloop.eq.800) stop

!  baryon and electron densities fixed for now
!   fcn_nb = fcn_crustSna_inp%nb
   fcn_crust%ne = fcn_nb * fcn_crust%Ye
!********************************************************************
!     defines gas density: ng = ng(nb,ne,xi,ncl)
!********************************************************************
 xg = 2._pr * fcn_crust%ne / (1._pr - fcn_crust%xi)

 fcn_crust%ng = ( fcn_nb - xg)/(1._pr - xg/fcn_crust%ncl)

!********************************************************************
!     defines xa : A = A(ne,ncl,I)
!********************************************************************
   u   = 2._pr * fcn_crust%ne/(fcn_crust%ncl * (1._pr - fcn_crust%xi))

   f_u = 1._pr - 1.5_pr * u**(1._pr/3._pr) + 0.5_pr * u

   Yp  = (1._pr - fcn_crust%xi)/2._pr

   sigma = fcn_crustSna_inp%sigSurf * (2._pr**(fcn_crustSna_inp%pSurf+1._pr) + fcn_crustSna_inp%bSurf)/(Yp**(- fcn_crustSna_inp%pSurf) + fcn_crustSna_inp%bSurf + (1._pr - Yp)**(- fcn_crustSna_inp%pSurf))


    o3ncl23     = (3._pr /fcn_crust%ncl)**(2._pr/3._pr)

    denom      =  (fcn_crustSna_inp%ac * f_u * (1._pr - fcn_crust%xi)**2 )

   fcn_crust%xa = ( 2._pr * sigma * (4._pr * pi)**(1._pr/3._pr))*o3ncl23/ denom



!********************************************************************
!     defines functions to be solved
!********************************************************************
!11111111111111111111111111111111111111111111111111111111111111111111
! first function > pressures equilibrium: P_cl = P_gas

! we take P_unif_cl and P_unif_g from the meta model code

 Z_cl = fcn_crust%xa * (1._pr - fcn_crust%xi)/2._pr

! cluster radius
 r_cl = ((3._pr * fcn_crust%xa)/(4._pr * pi * fcn_crust%ncl))**(1._pr/3._pr)

! defines surface energy (total)
 E_surf = 4._pr * pi * sigma * r_cl**2



! defines cluster pressure end energy density
! CALL META MODEL
! set densities : i_denb = ncl and xi=i_xd
! set output crust type variables ? o_eosDen = crustSna_crust ?

 call metaeos_T0_set_densities(fcn_crust%ncl,fcn_crust%xi,0._pr,fcn_coef_meta,.false.,eosDen_cl)
   ! if NR
 call metaeos_T0_baryons_NR(fcn_coef_meta,eosDen_cl,.false.,eosb0_cl)

 P_unif_cl = eosb0_cl%p_b

  ener_unif_cl = eosb0_cl%e2v_b

!  ener_unif_cl = eosb0_cl%rho_b

! write(*,'(a6,10f10.4)')'debug3',eosb0_cl%e2v_b,eosb0_cl%rho_b,fcn_crust%ncl,fcn_crust%xi

! derivative of e_unif with respect to xi is related to the mu_np

 d_E_unif_d_I = eosb0_cl%mu_np/2._pr

! defines cluster pressure

 P_cl = P_unif_cl + 0.5_pr * fcn_crustSna_inp%ac * Z_cl**2 * fcn_crust%xa**(- 4._pr/3._pr) * fcn_crust%ncl * u  * (u**(-2._pr/3._pr) - 1._pr) - 2._pr * E_surf * fcn_crust%ncl /(3._pr * fcn_crust%xa)


! write(*,'(a6,10f10.4)')'d',P_unif_cl , 0.5_pr * fcn_crustSna_inp%ac * Z_cl**2 * fcn_crust%xa**(- 4._pr/3._pr) * u * fcn_crust%ncl * (u**(-2._pr/3._pr) - 1._pr), 2._pr * E_surf * fcn_crust%ncl /(3._pr * fcn_crust%xa)

! defines neutron gas pressure and chemical potential from the meta model code
! CALL META MODEL
! set i_denb = ng and xi=0
 call metaeos_T0_set_densities(fcn_crust%ng,0._pr,0._pr,fcn_coef_meta,.false.,eosDen_gas)
! if NR
 call metaeos_T0_baryons_NR(fcn_coef_meta,eosDen_gas,.false.,eosb0_gas)


 P_gas = eosb0_gas%p_b

 fcn_crust%mu_g_n = eosb0_gas%mu_n


!2222222222222222222222222222222222222222222222222222222222222222222222222
! second function : chemical potential equilibrium: mu_n^cl = mu_n^g(unif) - P_gas/ncl


! derivative of sigma with respect to xi

 d_sigma_d_I = 0.5_pr * fcn_crustSna_inp%pSurf * sigma * ((1._pr - Yp)**(- fcn_crustSna_inp%pSurf-1._pr) - Yp**(- fcn_crustSna_inp%pSurf - 1._pr) )/(Yp**(- fcn_crustSna_inp%pSurf)+ fcn_crustSna_inp%bSurf + (1._pr - Yp)**(- fcn_crustSna_inp%pSurf))


! defines total coulomb energy: E_coul

 E_coul = fcn_crustSna_inp%ac * f_u * (Z_cl**2) /fcn_crust%xa**(1._pr/3._pr)

! defines total cluster energy: E_cl

 E_cl = ener_unif_cl * fcn_crust%xa /fcn_crust%ncl + E_coul + E_surf

! derivative of E_cl with respect to xi

  d_E_cl_d_I = d_E_unif_d_I * fcn_crust%xa &
 & + 0.125_pr * fcn_crustSna_inp%ac * fcn_crust%xa**(5._pr/3._pr) * (1._pr - fcn_crust%xi) * (u-u**(1._pr/3._pr) - 4._pr * f_u ) &
 & + 4._pr * pi * d_sigma_d_I * r_cl**2
! defines chemical potential of the neutrons in the cluster

  fcn_crust%mu_cl_n = d_E_cl_d_I * (1._pr - fcn_crust%xi)/fcn_crust%xa + E_cl/fcn_crust%xa

! call meta model to take electron chemical potential in uniform matter (i_opt_beta ?)
 !xd = 1._pr - 2._pr * fcn_crustSna_inp%Ye_init
  xd = 1._pr - 2._pr * fcn_crust%Ye

 call metaeos_T0_set_densities(fcn_nb,xd,0._pr,fcn_coef_meta,.false.,eosDen_el)

 call metaeos_T0_leptons(fcn_coef_meta,eosDen_el,.false.,eosl0)

 mu_el_unif = eosl0%mu_e


! write(*,'(a16,10f10.4)')'debug2-from meta',eosb0_cl%e2v_b,eosb0_cl%p_b,eosb0_cl%mu_np,eosb0_gas%p_b,eosl0%mu_e


! defines chemical potential of the electrons

  mu_el = mu_el_unif + (1._pr - u**(- 2._pr/3._pr)) * 0.5_pr * fcn_crust%ne * fcn_crustSna_inp%ac * fcn_crust%xa**(2._pr/3._pr)/fcn_crust%ncl

! defines chemical potential of the protons in the cluster

  fcn_crust%mu_cl_p = E_cl/fcn_crust%xa - (1._pr + fcn_crust%xi) * d_E_cl_d_I /fcn_crust%xa


!   write(*,*)'debug3',P_gas*CONV_MeV_fm3_to_dyn_cm2,P_cl*CONV_MeV_fm3_to_dyn_cm2,fcn_crust%mu_cl_n, (fcn_crust%mu_g_n + P_gas/fcn_crust%ncl),y(:)
!
!  Solve the system of non linear equations
!  Makes y(n) = 0
  y(1) = P_gas  - P_cl
  y(2) = fcn_crust%mu_cl_n - fcn_crust%mu_g_n + P_gas/fcn_crust%ncl
 ! y(3) = fcn_crust%mu_cl_p + mu_el +(CST_mpc2-CST_mnc2) - fcn_crust%mu_cl_n

! write(*,*)'debug3',-fcn_crust%mu_g_n,fcn_crust%mu_cl_n+P_gas/fcn_crust%ncl,y(:),x(:)
! write(*,*)'debug3', P_gas,P_cl,y(:),x(:)


  fcn_crust%press  = P_cl

!   x(1) = fcn_crust%ncl
!   x(2) = fcn_crust%xi

    iloop = iloop + 1
end subroutine crust_sna




subroutine crustSna_read_file(crustSna_inp,crust_eos)

  use spl;

  IMPLICIT None
  !
  ! =====================
  ! DEFINE VARIABLES
  ! =====================
  !
  ! ---------------------
  ! INPUT VARIABLES
  ! ---------------------
  !
  TYPE (crustSna_inputs), intent(in) :: crustSna_inp
  !
  ! ---------------------
  ! OUTPUT VARIABLES
  ! ---------------------
  !
  TYPE (crustSna_eos_outputs), intent(out) :: crust_eos
  !
  ! ---------------------
  ! LOCAL VARIABLES
  ! ---------------------
  !
  integer        :: ios
  real (kind=pr), dimension(:), allocatable :: xcrust_vs2, xcrust_vs2s
  real (kind=pr), dimension(:), allocatable :: res_data
  integer        :: iden, k
  integer        :: idum
  real (kind=pr) :: xden, rho, pre

  ! local variables
  !
  real (kind=pr) :: xy1, xh1, xy2, xh2, xdif
  real (kind=pr) :: xm, ym, xslope, xdeno
  real (kind=pr), dimension(:), allocatable :: help,au
  !
  ! =====================
  ! START ROUTINE
  ! =====================
  !
  if (crustSna_lverb) write(*,*)" Crust:"
  !
     !
     !  READ THE CRUST EOS FROM THE INPUT FILE: crusteos/inp%model
     !  xn  is the baryonic density (in fm-3)
     !  rho is the energy-density including the rest mass energy (in g/cm^3)
     !  pre is the pressure (in dyn/cm2)
     !
!     write(*,*)"  Read crust file : ","../nstar-EosCrust/"//trim(crustSna_inp%model)//"."
!     OPEN(UNIT=100,FILE="../nstar-EosCrust/"//trim(crustSna_inp%model))
     if (crustSna_lverb) write(*,*)"  Read crust file : ","nsEos-EosCrust/"//trim(crustSna_inp%model)//"."
     OPEN(UNIT=100,FILE="nsEos-EosCrust/"//trim(crustSna_inp%model))
     read(100,*)crust_eos%ndata ! number of points
     read(100,*)crust_eos%rho_cc ! core-crust energy-density in g/cm3
     read(100,*)crust_eos%rho_drip ! drip energy-density in g/cm3 (neutron drip density)
     !
     !
     do iden = 1, crust_eos%ndata
        READ(100,*,iostat=ios)idum,xden,rho,pre
        if (ios.ne.0) stop "tovEOS:setup_crust: exit crust file before the end"
        crust_eos%EosLog(1,iden)  = log10(xden)
        crust_eos%EosLog(2,iden) = log10(rho)
        crust_eos%EosLog(3,iden) = log10(pre)
        if (rho.le.crust_eos%rho_cc) crust_eos%i_cc = iden
     enddo
     close(100)
     !
     ! The maximum density is defined as the highest density in the EoS profile, or 12*nsat.
     !  nbmaxcu=min(10**xcunb(ndimcu),12.0*CST_nsat)
     !
     crust_eos%denb_max = 10**crust_eos%EosLog(1,crust_eos%ndata)
     if (crustSna_lverb) write(*,*)'  crust_eos%ndata = ',crust_eos%ndata
     if (crustSna_lverb) write(*,*)'  den_init = ',10**crust_eos%EosLog(1,1),'fm-3',10**crust_eos%EosLog(2,1),'g/cm3'
     if (crustSna_lverb) write(*,'(T4,a32,d10.3,a6)')"Neutron drip density          = ",crust_eos%rho_drip,"g/cm3"
     if (crustSna_lverb) write(*,'(T4,a32,d10.3,a6)')"Core-crust transition density = ",crust_eos%rho_cc,"g/cm3"
     if (crustSna_lverb) write(*,*)'  crust%i_cc = ',crust_eos%i_cc
     if (crustSna_lverb) write(*,'(a31,2(d14.4,a6))')'  transition energy-density : ',&
          &10**crust_eos%EosLog(2,crust_eos%i_cc),'g/cm3',10**crust_eos%EosLog(2,crust_eos%i_cc)/CST_rsat,'nsat'
     if (crustSna_lverb) write(*,'(a31,2(d14.4,a6))')'  maximum energy-density    : ',&
          &10**crust_eos%EosLog(2,crust_eos%ndata),'g/cm3',10**crust_eos%EosLog(2,crust_eos%ndata)/CST_rsat,'nsat'

     ! calculate the sound velocity using spline routine: vs2 = dP / drho
     allocate(help(crust_eos%ndata)); allocate(au(crust_eos%ndata)); allocate(res_data(crust_eos%ndata))
!     call spls3(xcurho,xcupre,ndimcu,xcurho,xcuvs2s,ndimcu,help,au,2,0)
     call spls3(crust_eos%EosLog(2,:),crust_eos%EosLog(3,:),crust_eos%ndata,crust_eos%EosLog(2,:),res_data(:),crust_eos%ndata,help,au,2,0)
     crust_eos%EosLog(4,:) = res_data(:) * 10.**crust_eos%EosLog(3,:) / ( 10.**crust_eos%EosLog(2,:) * CST_clight**2 )
     deallocate(help,au,res_data)

     allocate(xcrust_vs2(crust_eos%ndata))
     allocate(xcrust_vs2s(crust_eos%ndata))

     ! calculate the sound velocity by centered difference interpolation: vs2 = dP / drho
     xcrust_vs2(1)=0._pr
     DO iden = 2, crust_eos%ndata - 1
        xy1 =   10.**crust_eos%EosLog(3,iden+1) - 10.**crust_eos%EosLog(3,iden)
        xh1 = ( 10.**crust_eos%EosLog(2,iden+1) - 10.**crust_eos%EosLog(2,iden) ) * CST_clight**2
        xy2 =   10.**crust_eos%EosLog(3,iden)   - 10.**crust_eos%EosLog(3,iden-1)
        xh2 = ( 10.**crust_eos%EosLog(2,iden)   - 10.**crust_eos%EosLog(2,iden-1) ) * CST_clight**2
        xdif = ( xy1 / xh1 - xy2 / xh2 ) / ( xy1 / xh1 + xy2 / xh2 )
        ! check if there is not discontinuity in the derivatives
        if (dabs(xdif).le.0.5.or.iden.lt.12) then
           xcrust_vs2(iden) = 0.5 * xy1 / xh1 + 0.5 * xy2 / xh2
        else
           ! if there is a discontinuity
           ! use a simple linear extrapolation
           !xcuvs2(i) = xcuvs2(i-1) + ( xcuvs2(i-1) - xcuvs2(i-2) ) * &
           !     & ( 10.**xcurho(i) - 10.**xcurho(i-1) ) / &
           !     & ( 10.**xcurho(i-1) - 10.**xcurho(i-2) )
           ! use linear regression to extrapolate
           xm = 0.0
           ym = 0.0
           do k = iden - 4, iden - 1, 1
              ym = ym + xcrust_vs2(k)
              xm = xm + 10.**crust_eos%EosLog(2,k)
           enddo
           ym = ym / 4.0
           xm = xm / 4.0
           xslope = 0.0
           xdeno = 0.0
           do k = iden - 4, iden - 1, 1
              xslope = xslope + ( 10.**crust_eos%EosLog(2,k) - xm ) * (  xcrust_vs2(k) - ym )
              xdeno = xdeno + ( 10.**crust_eos%EosLog(2,k) - xm )**2
           enddo
           xslope = xslope / xdeno
           xcrust_vs2(iden) = ym + xslope * ( 10.**crust_eos%EosLog(2,iden) - xm )
           !xcrust_vs2(i) = 0.5 * xy1 / xh1 + 0.5 * xy2 / xh2
        endif
        !write(68,*)10.**crust%EosLog(2,iden),xy1/xh1,xy2/xh2,xcrust_vs2(iden),crust%EosLog(4,iden),xdif
     ENDDO
     xcrust_vs2(crust_eos%ndata) = 2. * xcrust_vs2(crust_eos%ndata-1) - xcrust_vs2(crust_eos%ndata-2)

     ! smooth the spline over 3 points to soften the oscillations
     xcrust_vs2s(1) = 0._pr
     DO iden = 2, crust_eos%ndata - 1
        !xcuvs2m(i) = ( xcuvs2s(i-1) + xcuvs2s(i) + xcuvs2s(i+1) ) / 3.0
        ym = crust_eos%EosLog(4,iden)
        xm = 1.0
        do k = iden - 2, iden + 2
           if (k.eq.iden) cycle
           ym = ym + crust_eos%EosLog(4,k) / dabs( crust_eos%EosLog(2,k) - crust_eos%EosLog(2,iden) )**2
           xm = xm + 1.0 / dabs( crust_eos%EosLog(2,k) - crust_eos%EosLog(2,iden) )**2
        enddo
        xcrust_vs2s(iden) = ym / xm
     ENDDO
     xcrust_vs2s(crust_eos%ndata) = 2.*xcrust_vs2s(crust_eos%ndata-1) - xcrust_vs2s(crust_eos%ndata-2)

     OPEN(UNIT=67,FILE="crustSna-res/crust-"//trim(crustSna_inp%model))
     write(67,*)" den (fm-3), rho (g/cm3), pression (Dyn.cm2), vs/c**2"
     DO iden = 1, crust_eos%ndata
        write(67,'(6d15.5)')10.**crust_eos%EosLog(1,iden),10.**crust_eos%EosLog(2,iden),10.**crust_eos%EosLog(3,iden),&
             &crust_eos%EosLog(4,iden),xcrust_vs2(iden),xcrust_vs2s(iden)
     ENDDO
     CLOSE(67)
     !
     ! write the include file for the fast code
     !
     OPEN(UNIT=100,FILE="crustSna-res/inc-"//trim(crustSna_inp%model))
     write(100,'("crust_eos%ndata    = ",i5)')crust_eos%ndata
     write(100,'("crust_eos%rho_cc   = ",d10.4)')crust_eos%rho_cc ! core-crust energy-density in g.cm3 neutron drip density
     write(100,'("crust_eos%i_cc     = ",i5)')crust_eos%i_cc
     write(100,'("crust_eos%rho_drip = ",d10.4)')crust_eos%rho_drip ! drip energy-density in g.cm3 neutron drip density
     do iden = 1, crust_eos%ndata
        write(100,'("crust_eos%EosLog(1,",i3,")=",f10.4)')iden,crust_eos%EosLog(1,iden)
        write(100,'("crust_eos%EosLog(2,",i3,")=",f10.4)')iden,crust_eos%EosLog(2,iden)
        write(100,'("crust_eos%EosLog(3,",i3,")=",f10.4)')iden,crust_eos%EosLog(3,iden)
        write(100,'("crust_eos%EosLog(4,",i3,")=",d10.4)')iden,crust_eos%EosLog(4,iden)
     enddo
     close(100)

     deallocate(xcrust_vs2,xcrust_vs2s)
  !
  ! =====================
  ! END ROUTINE
  ! =====================
  !

end subroutine crustSna_read_file




subroutine crustSna_write_eos(crust_eos)

  use spl;

  IMPLICIT None
  !
  ! =====================
  ! DEFINE VARIABLES
  ! =====================
  ! ---------------------
  ! INPUT VARIABLES
  ! ---------------------
  !
  TYPE (crustSna_eos_outputs), intent(in)  :: crust_eos
  !
  ! ---------------------
  ! LOCAL VARIABLES
  ! ---------------------
  !
  integer :: iden, ii
  integer :: stat_alo
  character :: err_alo
  !
  real (kind=pr) :: den
  real (kind=pr), dimension(1) :: denLog, res_den, res_rho, res_pre, res_cs2
  real (kind=pr), dimension(:), allocatable :: help,au
  !
  if (crustSna_lverb) write(*,*)" Write file:"
  !
  ! =====================
  ! START ROUTINE
  ! =====================
  !
  if (crustSna_lverb) write(*,*) "  Write: crustSna-res/eos.out"

  open(unit=20,file="crustSna-res/eos.out",status="unknown")
  write(20,*)"# Neutron drip density          = ",crust_eos%rho_drip,"g/cm3"
  write(20,*)"# Core-crust transition density = ",crust_eos%rho_cc,"g/cm3"
  write(20,*)"# Max value of baryon-density   = ",crust_eos%denb_max,"fm-3"
  write(20,*)"# Max value of baryon-density   = ",crust_eos%rho_max,"g/cm3"
  write(20,*)"# Number of data                = ",crust_eos%ndata
  write(20,*)"# i, den (fm-3), rho (g/cm3), pression (Dyn/cm2), vs/c**2"

  do iden = 1, crust_eos%ndata
     write(20,'(i6,4d15.5)') &
          iden, &
          10**crust_eos%eosLog(1,iden), &
          10**crust_eos%eosLog(2,iden), &
          10**crust_eos%eosLog(3,iden), &
          crust_eos%eosLog(4,iden)
  enddo
  close(20)

  if (crustSna_lverb) write(*,*) "  Write: crustSna-res/eos-lin.out"

  open(unit=20,file="crustSna-res/eos-lin.out",status="unknown")
  write(20,*)"# Neutron drip density          = ",crust_eos%rho_drip,"g/cm3"
  write(20,*)"# Core-crust transition density = ",crust_eos%rho_cc,"g/cm3"
  write(20,*)"# Max value of baryon-density   = ",crust_eos%denb_max,"fm-3"
  write(20,*)"# Number of data                = ",crust_eos%ndata
  write(20,*)"# i, den (fm-3), rho (g/cm3), pression (Dyn/cm2), vs/c**2"

  print*, "BEFORE ALLOCATE"

  allocate(help(crust_eos%ndata))
  allocate(au(crust_eos%ndata))

  print*, "NDATA = ", crust_eos%ndata

  den = 0.01
  do
     !     denLog(1) = log10(939*den*1.78266181e12) ! conversion from MeV.fm-3 to g cm-3

     denLog(1) = log10(den)
     call spls3(crust_eos%EosLog(1,:), crust_eos%EosLog(1,:), &
          crust_eos%ndata, denLog, res_den, 1, help, au, 1, 0)

     call spls3(crust_eos%EosLog(1,:), crust_eos%EosLog(2,:), &
          crust_eos%ndata, denLog, res_rho, 1, help, au, 1, 0)

     call spls3(crust_eos%EosLog(1,:),crust_eos%EosLog(3,:), &
          crust_eos%ndata,denLog,res_pre,1,help,au,1,0)

     call spls3(crust_eos%EosLog(1,:), crust_eos%EosLog(4,:), &
          crust_eos%ndata, denLog, res_cs2, 1, help, au, 1, 0)

     write(20,'(4d15.5)') den, 10**res_rho(1), 10**res_pre(1), res_cs2(1)

     den = den + 0.01

     !     if (crustSna_lverb) write(*,*)den,crust%denb_max,denLog(1),res_den(1),10**res_rho(1)
     if (den.gt.crust_eos%denb_max) exit
  enddo
  close(20)

  print*, "BEFORE DEALLOCATE"

  if(allocated(help)) then
     deallocate(help)
  endif

  print*, "DEALLOCATED FIRST OK"

  if(allocated(au)) then
     print*, "Allocated au"
     deallocate(au, stat=stat_alo, errmsg=err_alo)
     print*, stat_alo, err_alo
  endif

  print*, "AFTER DEALLOCATE"

  if (crustSna_lverb) write(*,*)"  Write: crustSna-res/eos-log.out"
  open(unit=20,file="crustSna-res/eos-log.out",status="unknown")
  write(20,*)"# Neutron drip density          = ",crust_eos%rho_drip,"g/cm3"
  write(20,*)"# Core-crust transition density = ",crust_eos%rho_cc,"g/cm3"
  write(20,*)"# Max value of baryon-density   = ",crust_eos%denb_max,"fm-3"
  write(20,*)"# Number of data                = ",crust_eos%ndata
  write(20,*)"# i, den (fm-3), rho (g/cm3), pression (Dyn/cm2), vs/c**2"
  !
  allocate(help(crust_eos%ndata)); allocate(au(crust_eos%ndata));
!  denLog(1) = -2.0*939*1.78266181e12 ! conversion from fm-3 to g cm-3
  denLog(1) = -12
  do
     den = 10**denLog(1)
     if (den.gt.crust_eos%denb_max) exit
     call spls3(crust_eos%EosLog(1,:),crust_eos%EosLog(2,:),crust_eos%ndata,denLog(1),res_rho(1),1,help,au,1,0)
     call spls3(crust_eos%EosLog(1,:),crust_eos%EosLog(3,:),crust_eos%ndata,denLog(1),res_pre(1),1,help,au,1,0)
     call spls3(crust_eos%EosLog(1,:),crust_eos%EosLog(4,:),crust_eos%ndata,denLog(1),res_cs2(1),1,help,au,1,0)
     write(20,'(4d15.5)')den,10**res_rho(1),10**res_pre(1),res_cs2(1)
!     if (crustSna_lverb) write(*,*)den,crust%denb_max,10**res_den(1),10**res_rho(1)
     denLog(1) = denLog(1) + 0.1
  enddo
  close(20)
  deallocate(help,au)
  !
  ! =====================
  ! END ROUTINE
  ! =====================
  !
end subroutine crustSna_write_eos




end MODULE crustSnaMod
