! Fortran95 line length stop there --->                                                                                          132

!
! =============================================================
!
!                    crustSnaMAIN.F90
!
! =============================================================
!
!
! --------------------------NOTES------------------------------
!
! Solve the crust EoS assuming single Nucleus approximation (SNA)
!
! Input file: eos.dat
!             Tabulated EoS for the crust and the core of NS
!             density in g/cm3
!             pressure in dynes/cm2
!
! nb is the baryon density
! rho is the mass/energy-density
!
! -------------------------------------------------------------
! Code by Jerome Margueron
! INT Seattle, USA and IPN Lyon, France
! Date: 22-01-2018
! Latest revision - 22-01-2018
! -------------------------------------------------------------

PROGRAM crustSnaMain

  use acc; use cst;
  use crustSnaType; use crustSnaMod;
  use metaEosT0Type; use metaEosT0;
  use polyEosT0Type; use polyEosT0;

  IMPLICIT None
  !
  ! =====================
  ! DEFINE VARIABLES
  ! =====================
  !
  ! ---------------------
  ! define typed variables
  ! ---------------------
  TYPE (crustSna_inputs)      :: crustSna_inp
  TYPE (metaEos_coef)         :: coef_meta
  TYPE (polyEos_coef)         :: coef_poly
  TYPE (crustSna_eos_outputs) :: crust_eos
  !
  ! ---------------------
  ! local variables
  ! ---------------------
  !
  ! CPU time
  real (kind=pr) :: t0, ttot, cpu0, cpu1, cputot
  !
!  logical                      :: withq = .true.
  logical                      :: withq = .false.
  !
  !
  ! =====================
  ! START OF MAIN CODE
  ! =====================
  !
  write(*,'(/40("-"))')
  write(*,*)"Start nsEosMain"
  write(*,'(40("-")/)')
  !
  ! ---------------------
  ! CPU time
  ! ---------------------
  !
  t0=second()
  call cpu_time(cpu0)
  !
  ! ---------------------
  ! create folder to store
  ! the results
  ! ---------------------
  !
  call system("mkdir -p crustSna-res/")
  !
  ! ---------------------
  !  READ INPUT PARAMETERS
  ! ---------------------
  !
  if (crustSna_lverb) write(*,*)"- read parameters for nsEOS:"
  !
  call crustSna_read_inputs(crustSna_inp)
  !
  ! Read arguments for the core EOS
  !
  if (crustSna_inp%model.eq."meos") then
     if (crustSna_lverb) write(*,*)"- read parameters for the metaEOS:"
     call metaEos_read_inputs(coef_meta)
  else if (crustSna_inp%model.eq."poly") then
     if (crustSna_lverb) write(*,*)"- read parameters for the polytropic EOS:"
     call polyEos_read_inputs(coef_poly)
  else
     if (crustSna_lverb) write(*,*)"- crust EOS read from the file",crustSna_inp%model
  endif
  !
  !
  if (crustSna_lverb) write(*,*)"- compute the EOS from the crust to the core:"
  !
  ! ---------------------
  ! DEFINE THE EOS FROM THE CRUST TO THE CORE
  ! ---------------------
  !
  call crustSna_compute_eos(coef_meta,coef_poly,crustSna_inp,withq,crust_eos)
  !
  ! ---------------------
  ! WRITE EOS IN FILE
  ! ---------------------
  !
  if (crustSna_lverb) write(*,*)"- write eos in file:"
  !
  call crustSna_write_eos(crust_eos)
  !
  print*, "FIM"
  !
  ! ---------------------
  ! CPU time
  ! ---------------------
  !
  ttot=second()-t0
  if(ttot.lt.0.0d0) ttot=ttot+86.4d3
  call cpu_time(cpu1)
  cputot=cpu1-cpu0
  write(*,*)
  write(*,'(''  Elapsed time = '',f10.2,'' seconds'')') ttot
  write(*,'(''      CPU time = '',f10.2,'' seconds'')') cputot
  write(*,*)
  !
  ! =====================
  ! END OF MAIN CODE
  ! =====================
  !
  write(*,'(/40("-"))')
  write(*,*)"End crustSnaMain"
  write(*,'(40("-")/)')
  !
end PROGRAM crustSnaMain



!  -------------------------------------------------------------
!  -------------------------------------------------------------
!
!    END OF THE MAIN CODE
!
!  -------------------------------------------------------------
!  -------------------------------------------------------------