! Fortran95 line length stop there --->                                                                                          132
! =============================================================
!
!                      MODULE crustSnatype
!
! =============================================================

! --------------------------NOTES------------------------------
!
!
! ----------------------Public types---------------------------
!
!
! -------------------------------------------------------------
! Code by Jerome Margueron
! Date: 15-09-2009
! Latest revision - 10-05-2011
! -------------------------------------------------------------


module crustSnaType

  use acc;

  implicit none

  ! ===================================================
  ! PARAMETERS
  ! ===================================================
 
  ! Fixes the number of point density for the crust EOS
  integer, parameter   :: crustSna_ncrust = 500

  ! verbose mode (=1), noverbose mode (=0)
  logical, parameter :: crustSna_lverb = .true.
!  logical, parameter :: crustSna_lverb = .false.


  ! ===================================================
  ! TYPED VARIABLES
  ! ===================================================

  TYPE crustSna_inputs
     character (len=20) :: model
     real (kind=pr)     :: nb, temp
     real (kind=pr)     :: Ye_init, ncl_init, I_init
     integer            :: dobeta, doloopnb
     real (kind=pr)     :: den_crust_min, den_crust_max
     real (kind=pr)     :: ac, apair
     real (kind=pr)     :: sigSurf, bSurf, pSurf
 END TYPE crustSna_inputs
  
  TYPE crustSna_eos_outputs
     ! number of points in the xxx.eos file
     integer                                :: ndata
     ! crust-core energy-density
     real (kind=pr)                         :: rho_cc
     integer                                :: i_cc
     ! drip energy-density
     real (kind=pr)                         :: rho_drip
     ! maximal density (in fm-3)
     real (kind=pr)                         :: denb_max
     ! maximal density (in fm-3)
     real (kind=pr)                         :: rho_max
     ! array with the EOS
     !   -crust%EosLog(1,i): baryon density (fm-3)
     !   -crust%EosLog(2,i): total energy density (g cm-3)
     !   -crust%EosLog(3,i): total pressure (dyn cm-2)
     !   -crust%EosLog(4,i): square of the total sound velocity (in c2)
     real (kind=pr), dimension(4,crustSna_ncrust)    :: EosLog
  END type crustSna_eos_outputs

  TYPE crustSna_crust
     real (kind=pr) :: xa, xi, ncl, Ye, ng
     real (kind=pr) :: ne
     real (kind=pr) :: mu_cl_n, mu_cl_p
     real (kind=pr) :: mu_g_n
     real (kind=pr) :: press
  END type crustSna_crust

  

end module crustSnaType
