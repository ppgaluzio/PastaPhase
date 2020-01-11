!
! =============================================================
!
!                        MODULE ACC
!
! =============================================================
!
!
!
! --------------------------NOTES------------------------------
!
!  Defines the precision of the real variables
!
!
! -------------------------------------------------------------
! Code by Jerome Margueron
! INT Seattle, USA and IPN Lyon, France
! Date: 22-01-2018
! Latest revision - 22-01-2018
! -------------------------------------------------------------

module ACC

  implicit none

  public

  integer, parameter :: pr = kind( 0.0d0 )
!  integer, parameter :: pr = selected_real_kind( p = 16 )

end module ACC
