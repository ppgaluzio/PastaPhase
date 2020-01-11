!
! =============================================================
!
!                      MODULE eeostype
!
! =============================================================
!
! --------------------------NOTES------------------------------
!
!  Defines the types used in metaeos modules and modules
!  
!
! ----------------------Public types---------------------------
!
!
! -------------------------------------------------------------
! Code by Jerome Margueron
! INT Seattle, USA and IPN Lyon, France
! Date: 22-01-2018
! Latest revision - 22-01-2018
! -------------------------------------------------------------
!
!
module metaEosT0Type
  
  use ACC; use CST;

  implicit none

  !
  ! Parameters
  !
  ! Number of parameters to read: 13 empirical parameters + ikin + imuon
  ! plus 5 variables for non-quadraticities
  ! and b is changed in to bsat and bsym
  integer, parameter :: metaEos_neparam = 25

  ! verbose mode (=1), noverbose mode (=0)
!  logical, parameter :: metaEos_lverb = .true.
  logical, parameter :: metaEos_lverb = .false.

  
  TYPE metaEos_Empirical
     integer, dimension(metaEos_neparam) :: eparam
     real (kind=pr), dimension(0:2,0:4)  :: param
     real (kind=pr)                      :: nsat, ms, dm
     real (kind=pr)                      :: bsat, bsym
     real (kind=pr)                      :: LQCD, kQCD
     integer                             :: ikin, imuon, irel, iquark
  END type metaEos_Empirical
    
  TYPE metaEos_Coef
!     real (kind=pr), dimension(0:1,0:4) :: VSIM !, VASIM
     integer, dimension(metaEos_neparam) :: eparam
     real (kind=pr), dimension(0:2,0:4)  :: param
     real (kind=pr), dimension(0:2,0:4)  :: VCOEF
     real (kind=pr)                      :: TFGsat, TFGSsat
     real (kind=pr)                      :: nsat, ms, dm, mb, db
     real (kind=pr)                      :: kappas, kappasym, kappav, kappaNM
     real (kind=pr)                      :: bsat, bsym
     real (kind=pr)                      :: LQCD, kQCD
     integer                             :: ikin, imuon, irel, iquark
     integer                             :: nmax ! cut-off in the density expansion
     integer, dimension(0:5)             :: facto ! factorial
  END type metaEos_Coef

    ! ---------------------------------------
    ! densities and related global quantities
    ! ---------------------------------------
  TYPE metaEos_Densities
     ! nucleonic model
     real (kind=pr)                  :: xx, xd, xn, xp, xe, xmuon
     real (kind=pr)                  :: den_noq, den_nuc, den_n, den_p, den_e, den_muon
     real (kind=pr)                  :: kf_nuc, kf_n, kf_p, kf_e, kf_muon
     real (kind=pr)                  :: kf_n_min, kf_p_min
     ! quark model
     real (kind=pr)                  :: q1_xup, q1_xdo
     real (kind=pr)                  :: q1_delta, q1_kf_q, q1_kf_qu, q1_kf_qd
     real (kind=pr)                  :: q1_nu, q1_nd, q1_nq, q1_xdq
     ! Baryon number
     real (kind=pr)                  :: den_b, xd_b
  END type metaEos_Densities

    ! ---------------------------------------
    ! quarkyonic at T=0
    ! ---------------------------------------
    TYPE metaEos_Q1
       ! kinetic energy
       real (kind=pr)                  :: rho_kin_B, rho_pot, rho_B ! Kinetic, potential, total energy with effective mass effect
       real (kind=pr)                  :: rho_kin_N, rho_N, rho_Q ! Kinetic, potential, total energy with effective mass effect
       real (kind=pr)                  :: rho_kin_noq, rho_noq, e2v_noq, e2a_noq ! Kinetic, potential, total energy with effective mass effect
       real (kind=pr)                  :: e2v, e2a ! energy density and energy per particle
       real (kind=pr)                  :: mun_qp, mup_qp, mu_B_noq  ! total chemical potential with it's contributions
       real (kind=pr)                  :: mu_kin_B, mu_pot_B, mu_B  ! total chemical potential with it's contributions
       real (kind=pr)                  :: mun_kin, mun_pot, mun, mup_kin, mup_pot, mup, munp_pot, munp_kin, munp ! pseudo chemical potentiels of quarks and nucleons
       real (kind=pr)                  :: muu, mud, muu_qp, mud_qp ! pseudo chemical potentiels of quarks and nucleons
       real (kind=pr)                  :: p_kin_B, p_pot_B, p_B, p_B_qp ! pressure
       real (kind=pr)                  :: pn_kin, pn_pot, pn, pn_qp, pq ! partial pressures of quarks and nucleons
       real (kind=pr)                  :: k_kin, k_pot, k_B ! incompressibility
       real (kind=pr)                  :: h_B ! enthalpy
       real (kind=pr)                  :: cs ! sound velocity
    END type metaEos_Q1
  
    ! ---------------------------------------
    ! baryons at T=0 or finite T
    ! ---------------------------------------
    TYPE metaEos_Baryons
       ! kinetic energy
       real (kind=pr)                  :: tau_b, tau_n, tau_p ! kinetic density
       real (kind=pr)                  :: t2aFG_b, t2vFG_b ! Fermi Gas kinetic energy
       real (kind=pr)                  :: t2a_b, t2v_b ! Kinetic energy with effective mass effect
       ! density correction for the potential energy
       real (kind=pr), dimension(0:4)  :: udc, udcsat, udcNM ! u for the choosen delta, u for delta=0, u for delta=1
       real (kind=pr), dimension(0:4)  :: d1_udc, d2_udc
       ! potential energy
       real (kind=pr), dimension(0:4)  :: a_v2a_order_n, a_v2a ! potential energy: order n, cumulative sum
       real (kind=pr), dimension(0:4)  :: a_e2a_b, a_e2v_b ! total energy (kinetic+potential)
       real (kind=pr)                  :: epot2a_b, epot2v_b ! potential energy for n=4
       ! total baryon energy
       real (kind=pr)                  :: e2a_b, e2v_b ! total energy (kinetic+potential) for n=4
       real (kind=pr)                  :: rho_b ! total energy including the rest mass
       ! symmetry energy 
       real (kind=pr)                  :: tsymFG_b ! Fermi Gas Symmetry energy
       real (kind=pr)                  :: tsym_b ! Kinetic symmetry energy with effective mass effect
       real (kind=pr), dimension(0:4)  :: a_vsym_order_n, a_vsym ! potential symmetry energy: order n, cumulative sum
       real (kind=pr), dimension(0:4)  :: a_esym_b ! total symmetry energy (kinetic+potential)
       real (kind=pr)                  :: esym_b ! total symmetry energy (kinetic+potential) for n=4
       real (kind=pr)                  :: esympot_b ! total potential symmetry energy for n=4
       ! symmetry energy S_2
       real (kind=pr)                  :: tsym2FG_b ! Fermi Gas Symmetry energy
       real (kind=pr)                  :: tsym2_b ! Kinetic symmetry energy with effective mass effect
       real (kind=pr), dimension(0:4)  :: a_vsym2_order_n, a_vsym2 ! potential symmetry energy: order n, cumulative sum
       real (kind=pr), dimension(0:4)  :: a_esym2_b ! total symmetry energy (kinetic+potential)
       real (kind=pr)                  :: esym2_b ! total symmetry energy (kinetic+potential) for n=4
       real (kind=pr)                  :: esym2pot_b ! total potential symmetry energy for n=4
       ! entropy
       real (kind=pr)                  :: s2a_b, s2v_b
       ! Chemical potentials
!       real (kind=pr), dimension(0:4)  :: a_mu_n, a_mu_p, a_mu_np
       real (kind=pr)                  :: mut_n, mut_p, mut_np
       real (kind=pr)                  :: mu_pot_n, mu_pot_p, mu_pot_np, mu_pot_is
       real (kind=pr)                  :: mu_n, mu_p, mu_np
       ! Mean-field
       real (kind=pr)                  :: umf_kin_n, umf_kin_p ! n and p mean field for the kinetic energy contribution
       real (kind=pr)                  :: umf_pot_n, umf_pot_p ! n and p mean field for the potential energy contribution
       real (kind=pr)                  :: umf_n, umf_p ! n and p mean field for n=4
       ! effective mass
       real (kind=pr)                  :: ms_n, ms_p
       ! effective chemical potentials
       real (kind=pr)                  :: nu_n, nu_p ! n and p effective chemical potentials
       real (kind=pr)                  :: nu_b ! baryon effective chemicla potential without effective mass
       ! 1st order devivatives
!       real (kind=pr), dimension(0:4)  :: a_d1_e2a_b_den ! d E/A / d n_b
!       real (kind=pr), dimension(0:4)  :: a_d1_e2a_b_xe ! d E/A / d xe
       ! pressure
       real (kind=pr)                  :: PFG_b ! Fermi Gas pressure
       real (kind=pr)                  :: Pt_b ! kinetic pressure with effective mass effect
       real (kind=pr)                  :: P_pot_b ! kinetic pressure with effective mass effect
!       real (kind=pr), dimension(0:4)  :: a_pv_order_n, a_pv ! potential pressure: order n, cumulative sum
!       real (kind=pr), dimension(0:4)  :: a_p_b ! total pressure
       real (kind=pr)                  :: p_b ! total pressure for n=4
       ! 2nd order devivatives
       ! incompressibility
       real (kind=pr)                  :: KFG_b ! Fermi Gas incompressibility
       real (kind=pr)                  :: Kt_b ! kinetic incompressibility with effective mass effect
       real (kind=pr)                  :: K_pot_b ! kinetic incompressibility with effective mass effect
!       real (kind=pr), dimension(0:4)  :: a_Kv_order_n, a_Kv ! potential incompresibility: order n, cumulative sum
!       real (kind=pr), dimension(0:4)  :: a_K_b
       real (kind=pr)                  :: K_b ! total incompressibility for n=4
       ! Baryon sound velocity square
!       real (kind=pr), dimension(0:4)  :: a_cs_b ! sound velocity per order
       real (kind=pr)                  :: cs_b ! sound velocity for n=4
       ! pressure derivatives
!       real (kind=pr), dimension(0:4)  :: a_d1_P_b_den ! d P / d n_b
!       real (kind=pr), dimension(0:4)  :: a_d1_P_b_xe ! d P / d x_e
       ! chemical potential derivatives
!       real (kind=pr), dimension(0:4)  :: a_d1_mu_np_xe
       ! 3rd order devivatives
!       real (kind=pr), dimension(0:4)  :: a_d3_e2a_b_den ! d^3 E/A / d n_b^3

    END type metaEos_Baryons
       
    ! ---------------------------------------
    ! electrons, muons and neutrons at T=0 and finite T
    ! ---------------------------------------
    TYPE metaEos_Leptons
       ! ---------------------------------------
       ! electrons
       ! ---------------------------------------
       real (kind=pr)                  :: mu_e ! chemical potential (with rest mass)
       real (kind=pr)                  :: nu_e ! chemical potential without rest mass (mu_e - mec2)
       real (kind=pr)                  :: nu_e_UR ! UR chemical potential
       real (kind=pr)                  :: rho_e ! energy-density including rest mass
       real (kind=pr)                  :: e2a_e, e2v_e ! energy
       real (kind=pr)                  :: p_e ! pressure
       real (kind=pr)                  :: e2a_e_UR, p_e_UR ! UR limit for energy and pressure
       real (kind=pr)                  :: s2a_e
       real (kind=pr)                  :: K_e, K_e_UR
       real (kind=pr)                  :: d1_e2a_e_den ! d e/a / n_b
       real (kind=pr)                  :: d1_mu_e_den ! d mue / d n_b
       real (kind=pr)                  :: d1_p_e_den ! d P / d n_b
       real (kind=pr)                  :: d1_e2a_e_xe ! d e/a / d xe
       real (kind=pr)                  :: d1_mu_e_xe ! d mue / d xe
       real (kind=pr)                  :: d1_p_e_xe ! d P / d xe
       ! ---------------------------------------
       ! muons
       ! ---------------------------------------
       real (kind=pr)                  :: mu_muon ! chemical potential (with rest mass)
       real (kind=pr)                  :: nu_muon ! chemical potential without rest mass (mu_muon - mmuonc2)
       real (kind=pr)                  :: rho_muon ! energy-density including rest mass
       real (kind=pr)                  :: e2a_muon, e2v_muon ! energy
       real (kind=pr)                  :: p_muon ! pressure
       real (kind=pr)                  :: nu_muon_UR ! UR chemical potential
       real (kind=pr)                  :: e2a_muon_UR, p_muon_UR ! UR limit for energy and pressure
       real (kind=pr)                  :: s2a_muon
       real (kind=pr)                  :: K_muon, K_muon_UR
       real (kind=pr)                  :: d1_e2a_muon_den ! d e/a / d n_b
       real (kind=pr)                  :: d1_mu_muon_den ! d mu / d n_b
       real (kind=pr)                  :: d1_p_muon_den ! d P / d n_b
       real (kind=pr)                  :: d1_e2a_muon_xmuon ! d e/a / d x_muon
       real (kind=pr)                  :: d1_mu_muon_xmuon ! d mu / d x_muon
       real (kind=pr)                  :: d1_p_muon_xmuon ! d P_muon / d x_muon
       ! ---------------------------------------
       ! neutrinos
       ! ---------------------------------------
       real (kind=pr)                  :: e2a_nu
       real (kind=pr)                  :: p_nu
       real (kind=pr)                  :: mu_nu
       real (kind=pr)                  :: s2a_nu
       real (kind=pr)                  :: K_nu
       real (kind=pr)                  :: d1_e2a_nu_den, d1_e2a_nu_xe, d1_e2a_nu_T
       real (kind=pr)                  :: d1_p_nu_den, d1_p_nu_xe, d1_p_nu_T
    END type metaEos_Leptons

    ! ---------------------------------------
    ! photons at T=0
    ! ---------------------------------------
    type metaEOS_Photons
       real (kind=pr)                  :: mu_phot
       real (kind=pr)                  :: e2a_phot
       real (kind=pr)                  :: s2a_phot
       real (kind=pr)                  :: p_phot
       real (kind=pr)                  :: d1_e2a_phot_den, d1_e2a_phot_xe, d1_e2a_phot_T
       real (kind=pr)                  :: d1_p_phot_den, d1_p_phot_xe, d1_p_phot_T
    END type metaEOS_Photons  

    ! ---------------------------------------
    ! total
    ! ---------------------------------------
    type metaEos_All
       real (kind=pr)                  :: den_b
       real (kind=pr)                  :: e2a, e2v ! total energy per baryon, energy density
       real (kind=pr)                  :: rho ! total energy density including the rest mass
       real (kind=pr)                  :: p ! total pressure
       real (kind=pr)                  :: s2a ! total entropy
       real (kind=pr)                  :: enthalpy ! total enthalpy
       real (kind=pr)                  :: K
!       real (kind=pr)                  :: d1_e2a_den ! d e/a / d n_b
!       real (kind=pr)                  :: d1_e2a_xe ! d e/a / d xe
!       real (kind=pr)                  :: d1_e2a_T ! d e/a / d T
!       real (kind=pr)                  :: d1_p_den ! d P / d n_b
!       real (kind=pr)                  :: d1_p_xe ! d P / d xe
!       real (kind=pr)                  :: d1_p_T ! d P / d T
       real (kind=pr)                  :: gamma ! Gamma factor for a polytropic EOS (Gama = d ln P / d ln eps)
       real (kind=pr)                  :: cs ! sound velocity square
       ! ---------------------------------------
       ! adiabatic index along a given path
       ! ---------------------------------------
       real (kind=pr)                  :: gam_b0, gam_bT, gam_tot
    END type metaEos_All

  end module metaEosT0Type
