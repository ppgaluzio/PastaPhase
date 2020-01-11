! Fortran95 line length stop there --->                                                                                          132

!
! =============================================================
!
!                    NSMAIN.F90
!
! =============================================================
!
!
! --------------------------NOTES------------------------------
!
! Solve the Tolman-Oppenheimer Equations using a tabulated
! equation of state
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

PROGRAM nsEosMain

  use acc; use cst;
  use nsEosType; use nsEosMod;
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
  TYPE (nsEos_inputs)   :: nsEos_inp
  TYPE (metaEos_coef)   :: coef_meta
  TYPE (polyEos_coef)   :: coef_poly
  TYPE (nsEos_eos)      :: nseos
  !
  ! ---------------------
  ! local variables
  ! ---------------------
  !
  ! CPU time
  real (kind=pr) :: t0, ttot, cpu0, cpu1, cputot
  !
  logical                      :: withq = .true. 
!  logical                      :: withq = .false. 
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
  call system("mkdir -p nsEos-res/")
  !
  ! ---------------------
  !  READ INPUT PARAMETERS
  ! ---------------------
  !
  if (nsEos_lverb) write(*,*)"- read parameters for nsEOS:"
  !
  call nsEos_read_inputs(nsEos_inp)
  !
  ! Read arguments for the core EOS
  !
  if (nsEOS_inp%core.eq."meos") then
     if (nsEos_lverb) write(*,*)"- read parameters for the metaEOS:"
     call metaEos_read_inputs(coef_meta)
  else if (nsEOS_inp%core.eq."poly") then
     if (nsEos_lverb) write(*,*)"- read parameters for the polytropic EOS:"
     call polyEos_read_inputs(coef_poly)
  else if (nsEOS_inp%core.eq."same") then
     if (nsEos_lverb) write(*,*)"- core EOS is the same as crust EOS"
  else
     stop "Issue defining core EOS in input file" 
  endif
  !
  !
  if (nsEos_lverb) write(*,*)"- compute the EOS from the crust to the core:"
  !
  ! ---------------------
  ! DEFINE THE EOS FROM THE CRUST TO THE CORE
  ! ---------------------
  !
  call nsEos_compute_eos(coef_meta,coef_poly,nsEos_inp,withq,nseos)
  !
  ! ---------------------
  ! WRITE EOS IN FILE
  ! ---------------------
  !
  if (nsEos_lverb) write(*,*)"- write eos in file:"
  !
  call nsEos_write_eos(nseos)
  !
  if (nsEos_lverb) write(*,*)"- write eos in table:"
  !
  call nsEos_write_table_eos(nseos,coef_meta)
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
  write(*,*)"End nsEosMain"
  write(*,'(40("-")/)')
  !
end PROGRAM nsEosMain



!  -------------------------------------------------------------
!  -------------------------------------------------------------
!
!    END OF THE MAIN CODE
!
!  -------------------------------------------------------------
!  -------------------------------------------------------------



