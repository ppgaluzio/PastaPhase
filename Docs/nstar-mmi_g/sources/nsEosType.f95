! Fortran95 line length stop there --->                                                                                          132
! =============================================================
!
!                      MODULE tovtype
!
! =============================================================

! --------------------------NOTES------------------------------
!
!  Defines the types used in rsn modules and modules called by 
!  rsn.
!
! ----------------------Public types---------------------------
!
!
! -------------------------------------------------------------
! Code by Jerome Margueron
! Date: 15-09-2009
! Latest revision - 10-05-2011
! -------------------------------------------------------------


module nsEosType

  use acc;

  implicit none

  ! ===================================================
  ! PARAMETERS
  ! ===================================================
 
  ! Fixes the number of point density for the crust EOS
  integer, parameter   :: nsEos_ncrust = 500
  ! Fixes the number of point density for the core EOS
  integer, parameter   :: nsEos_ncore = 1000
  ! Fixes the number of point density for the crust+core EOS
  integer, parameter   :: nsEos_neos = nsEos_ncrust + nsEos_ncore
  ! Number of empirical parameters read
!  integer, parameter   :: neparam

  ! verbose mode (=1), noverbose mode (=0)
!  logical, parameter :: nsEos_lverb = .true.
  logical, parameter :: nsEos_lverb = .false.
  ! fast mode (=1)
  logical, parameter :: nsEos_lfast = .false.
!  logical, parameter :: nsEos_lfast = .true.


  ! ===================================================
  ! TYPED VARIABLES
  ! ===================================================
  
  TYPE nsEos_inputs
    character (len=20) :: crust, core
    real (kind=pr)     :: temp
    integer            :: dobeta
    real (kind=pr)     :: iden_core_min, iden_core_max !, iden_core_step
 END TYPE nsEos_inputs
  
  TYPE nsEos_crust
     ! number of points in the xxx.eos file
     integer                                :: ndata
     ! crust-core energy-density
     real (kind=pr)                         :: rho_cc
     integer                                :: i_cc
     ! drip energy-density
     real (kind=pr)                         :: rho_drip
     ! maximal density (in fm-3)
     real (kind=pr)                         :: denb_max
     ! array with the EOS
     !   -crust%EosLog(1,i): baryon density (fm-3)
     !   -crust%EosLog(2,i): total energy density (g cm-3)
     !   -crust%EosLog(3,i): total pressure (dyn cm-2)
     !   -crust%EosLog(4,i): square of the total sound velocity (in c2)
     real (kind=pr), dimension(4,nsEos_ncrust)    :: EosLog
  END type nsEos_crust

  TYPE nsEos_core
     ! number of calculated points
     integer                                :: ndata
     integer                                :: ndata_rho
     real (kind=pr)                         :: denb_max, rho_max
     ! array with the EOS
     !   -crust%EosLog(1,i): baryon density (fm-3)
     !   -crust%EosLog(2,i): total energy density (g cm-3)
     !   -crust%EosLog(3,i): total pressure (dyn cm-2)
     !   -crust%EosLog(4,i): square of the total sound velocity (in c2)
     real (kind=pr), dimension(4,nsEos_ncore)    :: EosLog
     logical                                :: lstab, lcaus
  END type nsEos_core

  TYPE nsEos_eos
     ! number of calculated points
     integer                                :: ndata
     integer                                :: nmax
     integer                                :: ndata_rho
     real (kind=pr)                         :: denb_max, rho_max
     ! crust-core energy-density
     real (kind=pr)                         :: rho_cc
     integer                                :: i_cc
     ! drip energy-density
     real (kind=pr)                         :: rho_drip
     ! array with the EOS
     !   -crust%EosLog(1,i): baryon density (fm-3)
     !   -crust%EosLog(2,i): total energy density (g cm-3)
     !   -crust%EosLog(3,i): total pressure (dyn cm-2)
     !   -crust%EosLog(4,i): square of the total sound velocity (in c2)
     real (kind=pr), dimension(4,nsEos_neos)    :: EosLog
  END type nsEos_eos
  

end module nsEosType
