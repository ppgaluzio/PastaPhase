! =============================================================
! 	              MODULE fcn_crust
!
!    contains soubroutine to be linked with nle.f90
!
!    solve system of 2 non linear eqs if Ye fixed
!    solve system of 3 non linear eqs if beta equilibrium
! =============================================================
! --------------------------NOTES------------------------------
!  The soubroutine 'crust' takes n_b from the main program 
!  and the initial guess to ncl, ne and xi  
! 
!  Construct           ng = ng (den,ncl,ne,xi)   
!  and mass number     xa = xa (ne,ncl,xi) 
!
!  Defines functions to be minimized in terms of the
!  chemical potentials and pressures
!
!  Find the root of the functions using nle.f90 solver and
!  gives as output the variables (A,xi,ne,n_cl,ng)
!
! ----------------------private subroutines--------------------
!
! -------------------------------------------------------------
module crust_fcn


  implicit none

!  public

contains


!********************************************************************
! if Ye is fixed then = system of 2 functions to be solved
!********************************************************************

subroutine crust_sna(x,y,n)!fcn4(x,y,n)

  use acc; USE CST; USE crustSnaType;
  USE metaEosT0;  USE metaEosT0Type;

  implicit none

  ! ---------------------
  ! INPUT VARIABLES
  ! ---------------------
  !
  TYPE (crustSna_inputs), intent(in) :: crustSna_inp
  real (kind=pr),         intent(in) :: den
  !
  ! ---------------------
  ! OUTPUT VARIABLES
  ! ---------------------
  !
  TYPE (crustSna_crust), intent(out) :: crust

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
  real (kind=pr) :: r_cl, Y_p, sigma, d_sigma_d_I  
  real (kind=pr) :: m_p, m_n, pi  
  real (kind=pr) :: E_surf,  E_cl, E_coul
      
! declare pressures and chemical potentials
   real (kind=pr) :: mu_el, mu_el_unif, P_gas, P_cl, P_unif_cl    

!     parameters from module 'cst'

	m_p  = CST_mpc2
	m_n  = CST_mnc2
	pi   = CST_pi

! reads initial guess
   crust%ncl = x(1)
   crust%xi  = x(2)

!  baryon and electron densities fixed for now
   den = crustSna_inp%nb                      
   crust%ne = den * crustSna_inp%Ye_init
!********************************************************************
!     defines gas density: ng = ng(nb,ne,xi,ncl)   
!********************************************************************
 xg = 2._pr * crust%ne * (1._pr - crust%xi)
   
 crust%ng = ( den - xg)/(1._pr - xg/crust%ncl) 

!********************************************************************
!     defines xa : A = A(ne,ncl,I)    
!********************************************************************
   u   = 2._pr * crust%ne/(crust%ncl * (1._pr - crust%xi))

   f_u = 1._pr - 1.5_pr * u**(1._pr/3._pr) + 0.5_pr * u

   Yp  = (1._pr-crust%xi)/2._pr
	
   sigma = crustSna_inp%sigSurf * (2._pr**(crustSna_inp%pSurf+1._pr)+ crustSna_inp%bSurf)/(Yp**(-crustSna_inp%pSurf)+ crustSna_inp%bSurf + (1._pr - Yp)**(-crustSna_inp%pSurf))   

   crust%xa = ( 2._pr * sigma * (4._pr * pi)**(1._pr/3._pr) * (3._pr /crust%ncl)**(2._pr/3._pr))/(crustSna_inp%ac * f_u * (1._pr - crust%xi)**2 )

!********************************************************************
!     defines functions to be solved   
!********************************************************************
!11111111111111111111111111111111111111111111111111111111111111111111
! first function > pressures equilibrium: P_cl = P_gas

! we take P_unif_cl and P_unif_g from the meta model code

 Z_cl = crust%xa * (1._pr - crust%xi)/2._pr

! cluster radius
 r_cl = ((3._pr * crust%xa)/(4._pr * pi * crust%ncl))**(1._pr/3._pr)

! defines surface energy (total)
 E_surf = 4._pr * pi * sigma * r_cl**2

! defines cluster pressure
! CALL META MODEL
! set densities : i_denb = ncl and xi=i_xd
! set output crust type variables ? o_eosDen = crustSna_crust ?
 call metaeos_T0_set_densities(crust%ncl,crust%xi,0.,i_coef,i_withq,eosDen_cl)
! if NR
 call metaeos_T0_baryons_NR(i_coef,eosDen_cl,.false.,eosb0_cl)

 P_unif_cl = eosb0_cl%p_b

! defines cluster pressure
 
 P_cl = P_unif_cl + 0.5_pr * Z_cl**2 * crust%xa**(- 4._pr/3._pr) * u * crust%ncl * (u**(-2._pr/3._pr) - 1._pr) + 2._pr * E_surf * crust%ncl /(3._pr * crust%xa * r_cl)

! defines neutron gas pressure and chemical potential from the meta model code
! CALL META MODEL
! set i_denb = ng and xi=0
! set output crust type variables ? o_eosDen = crustSna_crust ?
 call metaeos_T0_set_densities(crust%ng,0._pr,0.,i_coef,i_withq,eosDen_gas)
! if NR
 call metaeos_T0_baryons_NR(i_coef,eosDen_gas,.false.,eosb0_gas)

 P_gas = eosb0_gas%p_b 

 crust%mu_g_n = eosb0_gas%mu_n

!22222222222222222222222222222222222222222222222222222222222222222222
! second function > chemical potential equilibrium: mu_n^cl = mu_n^g(unif) - P_gas/ncl

! we take the derivative of E_unif with respect to xi 
 
  d_E_unif_d_I = 

! derivative of sigma with respect to xi  

 d_sigma_d_I = 0.5_pr * crustSna_inp%pSurf * sigma * ((1._pr - Yp)**(-crustSna_inp%pSurf-1._pr) - Yp**(-crustSna_inp%pSurf-1._pr) )/(Yp**(-crustSna_inp%pSurf)+ crustSna_inp%bSurf + (1._pr - Yp)**(-crustSna_inp%pSurf))

! defines E_cl


 E_cl = 

! derivative of E_cl with respect to xi  
 
  d_E_cl_d_I = d_E_unif_d_I*crust%xa/crust%ncl - crustSna_inp%ac * f_u * (1._pr - crust%xi) * crust%xa**(5._pr/3._pr) + crust%ne * (1._pr - u**(- 2._pr/3._pr))/crust%ncl + 4._pr * pi * d_sigma_d_I * r_cl**2
! defines chemical potential of the neutrons in the cluster

  crust%mu_cl_n = d_E_cl_d_I * (1._pr - crust%xi)/crust%xa + E_cl/crust%xa

! call meta model to take electron chemical potential in uniform matter (i_opt_beta ?)
 xd = 1._pr - 2._pr * crustSna_inp%Ye_init
 call metaeos_T0_set_densities(den,xd,0.,i_coef,i_withq,eosDen_el)

 call metaeos_T0_leptons(i_coef,eosDen_el,i_opt_beta,eosl0)
 mu_el_unif = eosl0%mu_e

! defines chemical potential of the electrons

  mu_el = mu_el_unif + (1._pr - u**(- 2._pr/3._pr)) * 0.5_pr * crust%ne * crustSna_inp%ac * crust%xa**(2._pr/3._pr)/crust%ncl
  
! defines chemical potential of the protons in the cluster

  crust%mu_cl_n = - 2._pr * d_E_cl_d_I /crust%xa 
 
!
!  Solve the system of non linear equations
!  Makes y(n) = 0

  y(1) = P_gas - P_cl
  y(2) = crust%mu_cl_n - crust%mu_g_n + P_gas/crust%ncl
  
  crust%press  = P_cl

!   x(1) = ncl  
!   x(2) = xi   
   

end subroutine crust_sna


end module crust_fcn
