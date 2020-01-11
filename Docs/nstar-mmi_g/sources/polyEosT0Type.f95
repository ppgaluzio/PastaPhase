
! =============================================================
!
!                      MODULE polyEosT0Type
!
! =============================================================

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


module polyEosT0Type
  
  use ACC; use CST;

  implicit none

  !
  ! Parameters
  !
  integer, parameter :: polyEos_neparam = 6

  ! verbose mode (=1), noverbose mode (=0)
  logical, parameter :: polyEos_lverb = .true.

  
    
  TYPE polyEos_Coef
     real (kind=pr)                     :: logp1
     real (kind=pr), dimension(3)       :: gam
     real (kind=pr), dimension(2)       :: logrho
     !
     real (kind=pr), dimension(3)       :: k
     real (kind=pr), dimension(3)       :: a
     !
     real (kind=pr), dimension(2)       :: p
     real (kind=pr), dimension(2)       :: trho
     !
     real (kind=pr)                     :: nsat
  END type polyEos_Coef
   
    ! ---------------------------------------
    ! baryons at T=0 or finite T
    ! ---------------------------------------
  TYPE polyEos_eos
     real (kind=pr)                  :: p, rho, enthalpy, cs2, e2a

  END type polyEos_eos
       


  end module polyEosT0Type
