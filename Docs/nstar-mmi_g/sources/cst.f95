!
! =============================================================
!
!                        MODULE CST
!
! =============================================================
!
!
!
! --------------------------NOTES------------------------------
!
!  Defines usefull constants like pi, qp, hbc, htm, ...
!
!
! --------------------Public subroutines-----------------------
!
! CST_print()
!     print on the standard outputs the values of the constants
!     which are defined in this module.
!
! -------------------------------------------------------------
! Code by Jerome Margueron
! Date: 15-09-2009
! Latest revision - 10-05-2011
! -------------------------------------------------------------

module CST

  use ACC
  
  implicit none

  public

!  integer, parameter :: pr = selected_real_kind( p = 12 )

  real (kind=pr), parameter :: CST_pi     =    3.14159265358979_pr,&
                               CST_pi2    = CST_pi*CST_pi, &
                               CST_qp     = 4._pr* CST_pi, &
                               CST_hbc    =  197.32705_pr,  &  ! in MeV.fm
                               CST_mec2   =    0.510998928_pr,  &  ! electron mass in MeV
                               CST_mmuonc2=    105.6583715_pr,  &  ! muon mass in MeV
                               CST_mpc2   =  938.272013_pr,  &  ! proton mass in MeV
                               CST_mnc2   =  939.565346_pr, &  ! neutron mass in MeV
                               CST_mnuc2  =  0.5*(CST_mnc2+CST_mpc2), & ! nucleonic mass
                               CST_mlc2   = 1115.683_pr, &  ! lambda mass in MeV
!CST_mnc2 = 1.008665*931.494, CST_mpc2 = 1.007275*931.494, &
                               CST_mdeutc2= 1875.61282_pr, &  ! deuteron mass in MeV
                               CST_htm    = CST_hbc**2/CST_mnuc2,&
                               CST_htmn    = CST_hbc**2/CST_mpc2,&
                               CST_htmp    = CST_hbc**2/CST_mnc2,&
                               CST_html    = CST_hbc**2/CST_mlc2,&
                               CST_clight = 2.99792458d10, &  ! cm/s
                               CST_clight_mksa = 2.99792458d8, &  ! m/s
                               CST_Nc     = 3._pr, & ! number of color in QCD
                               CST_mqc2   = CST_mnuc2 / CST_Nc, & ! number of color in QCD
                               CST_p13    = 1._pr/3._pr,   &
                               CST_p23    = 2._pr/3._pr,   &
                               CST_p43    = 4._pr/3._pr,   &
                               CST_p53    = 5._pr/3._pr,   &
                               CST_G      = 6.67259d-8,    &  ! Grav. const in cm^3/g/s^2
                               CST_Msol   = 1.99d33,       &  ! solar mass in g
                               CST_Rsol   = 6.96d10,       &  ! solar radius in cm
                               CST_rshsol = 2._pr*CST_G*CST_Msol/(CST_clight**2),& ! Schw. radius of the sun in cm
                               CST_rsat   = 2.66d14,       &  ! saturation density (g/cm3)
                               CST_nsat   = 0.16_pr,       &  ! saturation density (fm-3)
                               CST_mfe    = 1.66d-24,      &  ! mass of 56Fe atom (in g)
                               CST_kb     = 8.617d-11,     &  ! Boltz. cst. (MeV.K^-1)
                               CST_ev2J   = 1.60217653d-19,&  ! 1 eV in Joule
                               CST_erg2J  = 1.d-7,         &  ! 1 erg in Joule
                               CST_ev2kg  = 1.78266181d-36,&  ! 1 eV/c2 in kg
                               CST_alpha  = 1.0/137.03599911,&! fine structure constant alpha=e^2/(4pi epsilon0 hbc) 
                               CST_alphahbc = CST_alpha * CST_hbc,& ! alpha in MeV.fm
                               CST_GF     = 1.16639d-11      ! Fermi const in MeV**{-2} hbc**3

  real (kind=pr), parameter :: CONV_MeV_fm3_to_g_cm3 = CST_ev2kg*1.d48,& !conversion from MeV.fm-3 -> g.cm-3
                               CONV_MeV_fm3_to_dyn_cm2 = CST_ev2kg*CST_clight_mksa**2*1.d52 ! conversion from MeV.fm-3 -> dyn.cm-2

contains

  !
  ! print constants
  !
  subroutine CST_print()

    implicit none

    REAL (KIND = PR), DIMENSION(5) :: rnd

    call random_number(rnd)
    write(*,'(1x,a23,5(f7.3,2x))')"Test random generator: ",rnd

    write(*,*)CST_pi,4.0_pr*atan(1.0_pr)
    write(*,*)CST_qp
    write(*,*)CST_hbc
    write(*,*)CST_mpc2
    write(*,*)CST_htm
    write(*,*)CST_p13
    write(*,*)CST_p43

  end subroutine CST_print



  SUBROUTINE init_random_seed()

    INTEGER :: i, n, clock
    INTEGER, DIMENSION(:), ALLOCATABLE :: seed
          
    CALL RANDOM_SEED(size = n)
    ALLOCATE(seed(n))
          
    CALL SYSTEM_CLOCK(COUNT=clock)
          
    seed = clock + 37 * (/ (i - 1, i = 1, n) /)
    CALL RANDOM_SEED(PUT = seed)
          
    DEALLOCATE(seed)

  END SUBROUTINE init_random_seed


  
end module CST

