! ===================================================
!
!                       MAIN PROGRAM
!
! ===================================================

! --------------------------NOTES------------------------------
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
! ===================================================

Program metaEosMain

  use acc; use cst;
  use metaEosT0; use metaEosT0Type;
  
  implicit none

  !
  !
  !
!  integer, dimension(neparam)  :: eparam
  TYPE (metaEos_coef) :: coef
  !
  logical                      :: withq = .true. 
!  logical                      :: withq = .false. 
  !
  ! CPU time
  real (kind=pr) :: t0, ttot, cpu0, cpu1, cputot
  !
  ! CPU time
  !
  t0=second()
  call cpu_time(cpu0)
  !
  ! Read arguments
  !
  call metaeos_read_inputs(coef)
  !
  ! calculate the meta-EOS for the MCMC
  !
  write(*,*)
  write(*,*)"Calculate metaEOS for the MCMC"
  write(*,*)
  call metaeos_compute_T0_MCMC(coef,withq)
  !
  ! calculate the meta-EOS in SM and NM
  !
  write(*,*)
  write(*,*)"Calculate metaEOS in SM and NM"
  write(*,*)
  call metaeos_compute_T0_SM_NM(coef,withq)
  !
  ! calculate the meta-EOS in beta-eq
  !
  write(*,*)
  write(*,*)"Calculate metaEOS in beta-eq."
  write(*,*)
  call metaeos_compute_T0_beta(coef,withq)
  !
  ttot=second()-t0
  if(ttot.lt.0.0d0) ttot=ttot+86.4d3
  call cpu_time(cpu1)
  cputot=cpu1-cpu0
  write(6,*)
  write(6,'(''  Elapsed time = '',f10.2,'' seconds'')') ttot
  write(6,'(''      CPU time = '',f10.2,'' seconds'')') cputot
  write(6,*)
 
end Program METAEOSMAIN



!  -------------------------------------------------------------
!  -------------------------------------------------------------
!
!    END OF THE MAIN CODE
!
!  -------------------------------------------------------------
!  -------------------------------------------------------------

