! Fortran95 line length stop there --->                                                                                          132
!
! =============================================================
!
!                        MODULE  nsEosMod
!
! =============================================================
!
!
! --------------------------NOTES------------------------------
!
!  Build the EOS in the crust and in the core
!  for the TOV solver.
!             nb   baryon number density in fm-3
!             rho  energy-density in g/cm3
!             pre  pressure in dynes/cm2
!
! -------------------------------------------------------------
! Code by Jerome Margueron
! Date: 15-09-2009
! Latest revision - 10-05-2011
! -------------------------------------------------------------

MODULE nsEosMod

  use acc; use cst;
  use nsEosType;

  IMPLICIT none

contains
    

! Instead of reading an input file, the data are hard coded
! This fasten the code
!
subroutine nsEos_read_inputs(nsEos_inp)

  use metaEosT0Type;

  IMPLICIT None
  !
  ! =====================
  ! DEFINE VARIABLES
  ! =====================
  !
  ! ---------------------
  ! OUTPUT VARIABLES
  ! ---------------------
  !
  TYPE (nsEos_inputs), intent(out) :: nsEos_inp
  !
  ! ---------------------
  ! LOCAL VARIABLES
  ! ---------------------
  !
  ! .....................
  ! get argument
  ! .....................
  !
  integer :: num_args
  character(len=20), dimension(:), allocatable :: args
  !
  ! .....................
  ! read input file
  ! .....................
  !
  character (len=20) :: crust, core
  real (kind=pr)     :: temp
  integer            :: dobeta
  real (kind=pr)     :: iden_core_min, iden_core_max
  namelist /nseos/crust, core, temp, dobeta, iden_core_min, iden_core_max
  integer            :: ix
  character (len=20) :: ifile_nsEos
  !
  ! =====================
  ! START ROUTINE
  ! =====================
  !
  ! Read arguments of the code
  !
  num_args = command_argument_count()
  !
  ! Select option according to the arguments
  !   If no argument, then take default values
  !   If argument, then read the parameter file
  !
  if (num_args.eq.0.or.num_args.eq.metaEos_neparam) then
!  if (nsEos_lfast) then
     !
     !  nsEos_inp%crust     = "sly_nd.dat" ! choose between "sly.dat" and "fps.dat" for the crust EoS
     !  nsEos_inp%crust     = "sly_lorene_sd.dat" ! choose between "sly.dat" and "fps.dat" for the crust EoS
!     nsEos_inp%crust     = "fps.dat" ! choose between "sly.dat" and "fps.dat" for the crust EoS
!     nsEos_inp%crust     = "fps.dat" ! choose between "sly.dat" and "fps.dat" for the crust EoS
     nsEos_inp%crust     = "sly.dat" ! choose between "sly.dat" and "fps.dat" for the crust EoS
!     nsEos_inp%core      = "same"    ! choose between "same" (=sly or fps) or "eeos" or "ceos"
     nsEos_inp%core  = "meos"    ! choose between "same" (=sly or fps) or "meos" for meta-EOS
     nsEos_inp%dobeta    =  3        ! for new and the beta eq.: 1-> use Esym, 2-> use chemical potential
                                   ! 3-> include muons
     nsEos_inp%temp      = 0.d0      ! temperature (MeV)
     nsEos_inp%iden_core_min  = 1.0      ! min density for the empirical eos in unit of den0
     nsEos_inp%iden_core_max  = 20.0     ! max density for the empirical eos in unit of den0
!     nsEos_inp%iden_core_step = 0.03     ! step density for the empirical eos in unit of den0
     !
  else if (num_args.eq.1.or.num_args.eq.2) then
     !
     allocate(args(num_args))
     do ix = 1, num_args
        call get_command_argument(ix,args(ix))
     end do
     !
     ! args contains the name of the input files:
     !    - args(1): nsEos.in
     !    - args(2): tov.in
     !
     read(args(1),'(a20)')ifile_nsEos
     !
     write(*,*)" nsEos input file :",ifile_nsEos
     open(unit=10,file=trim(adjustl(ifile_nsEos)),status="unknown")
!     open(unit=10, file="nsEos.in", status="old")
     read( unit = 10, nml = nseos )
     nsEos_inp%crust     = crust
     nsEos_inp%core      = core
     nsEos_inp%dobeta    = dobeta
     nsEos_inp%temp      = temp
     nsEos_inp%iden_core_min = iden_core_min ! in unit of den0
     nsEos_inp%iden_core_max = iden_core_max ! in unit of den0
!     nsEos_inp%iden_core_step = iden_core_step ! in unit of den0
     !
     close(10)
     !
  else
     !
     write(*,*)"setup_input: issue with arguments of the code"
     write(*,*)"Two options:"
     write(*,*)"   - No argument : take default arguments defined in setup_input"
     write(*,*)"   - 1 or 2 argument(s) : name of the input file(s) (less than 20 characters)"
     write(*,*)"Exemples:"
     write(*,*)"   - ./nsEos.e "
     write(*,*)"   - ./nsEos.e nsEos.in"
     write(*,*)"   - ./nsEos.e nsEos.in tov.in"
     write(*,*)"Code stopped"
     stop
     !
  end if

  if (nsEos_lverb) write(*,'(T4,a26,a20)') "nsEOS_inp%crust         = ",nsEOS_inp%crust
  if (nsEos_lverb) write(*,'(T4,a26,a20)') "nsEOS_inp%core          = ",nsEOS_inp%core
  if (nsEos_lverb) write(*,'(T4,a26,i2)')  "nsEOS_inp%dobeta        = ",nsEOS_inp%dobeta
  if (nsEos_lverb) write(*,'(T4,a26,f7.2)')"nsEOS_inp%temp          = ",nsEOS_inp%temp
  if (nsEos_lverb) write(*,'(T4,a26,f7.2)')"nsEOS_inp%iden_core_min  = ",nsEOS_inp%iden_core_min
  if (nsEos_lverb) write(*,'(T4,a26,f7.2)')"nsEOS_inp%iden_core_max  = ",nsEOS_inp%iden_core_max
!  if (nsEos_lverb) write(*,'(T4,a26,f7.2)')"nsEOS_inp%iden_core_step = ",nsEOS_inp%iden_core_step

  !
  ! =====================
  ! END ROUTINE
  ! =====================
  !
end subroutine nsEos_read_inputs









    
! There are two possibilities for the crust EOS which is matched to
! Either the sly EOS or the FPS EOS
! The choice is determined by the value of inp%eos_crust

subroutine nsEos_compute_eos(coef_meta,coef_poly,nsEos_inp,withq,nseos)

  use metaEosT0Type;
  use polyEosT0Type;


  IMPLICIT None
  !
  ! =====================
  ! DEFINE VARIABLES
  ! =====================
  !
  ! ---------------------
  ! INPUT VARIABLES
  ! ---------------------
  !
  TYPE (metaEos_coef), intent(in)   :: coef_meta
  TYPE (polyEos_coef), intent(in)   :: coef_poly
  TYPE (nsEos_inputs), intent(in)   :: nsEos_inp
  logical,             intent(in)   :: withq 
  !
  ! ---------------------
  ! OUTPUT VARIABLES
  ! ---------------------
  !
  TYPE (nsEos_eos),    intent(out)  :: nseos
  !
  ! local variables
  !
  TYPE (nsEos_crust)  :: crust
  TYPE (nsEos_core)   :: core
  !
  ! =====================
  ! START ROUTINE
  ! =====================
  !
  !
  ! Set up the crust EOS
  !
  call nsEos_compute_crust(nsEos_inp,crust)
  !
  ! Set up the core EOS
  !
  if (nsEos_inp%core.eq."meos") then
     call nsEos_compute_core_metaEos(coef_meta,nsEos_inp,withq,core)
  else if (nsEos_inp%core.eq."poly") then
     call nsEos_compute_core_polyEoS(coef_poly,nsEos_inp,core)
  endif
  !
  ! matches crust-core EOS
  !
  if (nsEos_lverb) write(*,*)"Match EOS: begin"
  call nsEos_match_eos(nsEos_inp,crust,core,nseos)
  if (nsEos_lverb) write(*,*)"Match EOS: end"
  !
  ! =====================
  ! END ROUTINE
  ! =====================
  !
end subroutine nsEos_compute_eos




! There are two possibilities for the crust EOS which is matched to
! Either the sly EOS or the FPS EOS
! The choice is determined by the value of inp%eos_crust
! The array crustEos is defined as:
!   -crust%EosLog(1,i): baryon density (fm-3)
!   -crust%EosLog(2,i): total energy density (g cm-3)
!   -crust%EosLog(3,i): total pressure (dyn cm-2)
!   -crust%EosLog(4,i): square of the total sound velocity (in c2)

subroutine nsEos_compute_crust(nsEos_inp,crust)

  use spl;
  
  IMPLICIT None
  !
  ! =====================
  ! DEFINE VARIABLES
  ! =====================
  !
  ! ---------------------
  ! INPUT VARIABLES
  ! ---------------------
  !
  TYPE (nsEos_inputs), intent(in) :: nsEos_inp
  !
  ! ---------------------
  ! OUTPUT VARIABLES
  ! ---------------------
  !
  TYPE (nsEos_crust), intent(out) :: crust 
  !
  ! ---------------------
  ! LOCAL VARIABLES
  ! ---------------------
  !
  integer        :: ios
  real (kind=pr), dimension(:), allocatable :: xcrust_vs2, xcrust_vs2s
  real (kind=pr), dimension(:), allocatable :: res_data
  integer        :: iden, k
  integer        :: idum
  real (kind=pr) :: xden, rho, pre

  ! local variables
  !
  real (kind=pr) :: xy1, xh1, xy2, xh2, xdif
  real (kind=pr) :: xm, ym, xslope, xdeno
  real (kind=pr), dimension(:), allocatable :: help,au
  !
  ! =====================
  ! START ROUTINE
  ! =====================
  !
  if (nsEos_lverb) write(*,*)" Crust:"
  !
  if (nsEos_lfast) then
     !
     ! nd: using numerical derivative for the speed of sound
     !
     if (trim(nsEos_inp%crust).eq."sly_nd.dat") then
!        include "nsEosCrust/inc-sly_nd.dat"
     elseif (trim(nsEos_inp%crust).eq."sly_lorene_nd.dat") then
!        include "nsEosCrust/inc-sly_lorene_nd.dat"
     else if (trim(nsEos_inp%crust).eq."fps_nd.dat") then
!        include "nsEosCrust/inc-fps_nd.dat"
     endif
     !
     ! sd: using spline derivative for the speed of sound
     !
     if (trim(nsEos_inp%crust).eq."sly_sd.dat") then
!        include "nsEosCrust/inc-sly_sd.dat"
     elseif (trim(nsEos_inp%crust).eq."sly_lorene_sd.dat") then
!        include "nsEosCrust/inc-sly_lorene_sd.dat"
     else if (trim(nsEos_inp%crust).eq."fps_sd.dat") then
!        include "nsEosCrust/inc-fps_sd.dat"
     endif
!     do iden = 1, crust%ndata !ndimcu
!        crust%EosLog(1,iden) = xcrust_nb(iden)
!        crust%EosLog(2,iden) = xcrust_rho(iden)
!        crust%EosLog(3,iden) = xcrust_pre(iden)
!        crust%EosLog(4,iden) = xcrust_cs2(iden)
!     enddo
     !
     ! maximum density
     !
     !nbmaxcu=10**xcunb(ndimcu)
     crust%denb_max = 10**crust%EosLog(1,crust%ndata)
     if (nsEos_lverb) write(*,*)'Fast setup of the crust EOS'
     if (nsEos_lverb) write(*,*)'EOS = ',trim(nsEos_inp%crust)
     if (nsEos_lverb) write(*,*)'crust%i_cc = ',crust%i_cc
     if (nsEos_lverb) write(*,'(a22,2d14.4)')'transition density : ',10**crust%EosLog(1,crust%i_cc),&
          &10**crust%EosLog(1,crust%i_cc)/CST_rsat
     if (nsEos_lverb) write(*,'(a22,2d14.4)')'maximum density    : ',10**crust%EosLog(1,crust%ndata),&
          &10**crust%EosLog(1,crust%ndata)/CST_rsat
     !
  else
     !
     !  READ THE CRUST EOS FROM THE INPUT FILE: crusteos/inp%crust
     !  xn  is the baryonic density (in fm-3)
     !  rho is the energy-density including the rest mass energy (in g/cm^3)
     !  pre is the pressure (in dyn/cm2)
     !
!     write(*,*)"  Read crust file : ","../nstar-EosCrust/"//trim(nsEos_inp%crust)//"."
!     OPEN(UNIT=100,FILE="../nstar-EosCrust/"//trim(nsEos_inp%crust))
     if (nsEos_lverb) write(*,*)"  Read crust file : ","nsEos-EosCrust/"//trim(nsEos_inp%crust)//"."
     OPEN(UNIT=100,FILE="nsEos-EosCrust/"//trim(nsEos_inp%crust))
     read(100,*)crust%ndata ! number of points
     read(100,*)crust%rho_cc ! core-crust energy-density in g/cm3
     read(100,*)crust%rho_drip ! drip energy-density in g/cm3 (neutron drip density)
     !
     if (crust%ndata.gt.nsEos_ncrust) stop "tovEOS:setup_crust: increase nsEos_ncrust"
     !
     do iden = 1, crust%ndata
        READ(100,*,iostat=ios)idum,xden,rho,pre
        if (ios.ne.0) stop "tovEOS:setup_crust: exit crust file before the end"
        crust%EosLog(1,iden)  = log10(xden) 
        crust%EosLog(2,iden) = log10(rho) 
        crust%EosLog(3,iden) = log10(pre)
        if (rho.le.crust%rho_cc) crust%i_cc = iden
     enddo
     close(100)
     !
     ! The maximum density is defined as the highest density in the EoS profile, or 12*nsat.
     !  nbmaxcu=min(10**xcunb(ndimcu),12.0*CST_nsat)
     !
     crust%denb_max = 10**crust%EosLog(1,crust%ndata)
     if (nsEos_lverb) write(*,*)'  crust%ndata = ',crust%ndata
     if (nsEos_lverb) write(*,*)'  den_init = ',10**crust%EosLog(1,1),'fm-3',10**crust%EosLog(2,1),'g/cm3'
     if (nsEos_lverb) write(*,'(T4,a32,d10.3,a6)')"Neutron drip density          = ",crust%rho_drip,"g/cm3"
     if (nsEos_lverb) write(*,'(T4,a32,d10.3,a6)')"Core-crust transition density = ",crust%rho_cc,"g/cm3"
     if (nsEos_lverb) write(*,*)'  crust%i_cc = ',crust%i_cc
     if (nsEos_lverb) write(*,'(a31,2(d14.4,a6))')'  transition energy-density : ',&
          &10**crust%EosLog(2,crust%i_cc),'g/cm3',10**crust%EosLog(2,crust%i_cc)/CST_rsat,'nsat'
     if (nsEos_lverb) write(*,'(a31,2(d14.4,a6))')'  maximum energy-density    : ',&
          &10**crust%EosLog(2,crust%ndata),'g/cm3',10**crust%EosLog(2,crust%ndata)/CST_rsat,'nsat'

     ! calculate the sound velocity using spline routine: vs2 = dP / drho
     allocate(help(crust%ndata)); allocate(au(crust%ndata)); allocate(res_data(crust%ndata))
!     call spls3(xcurho,xcupre,ndimcu,xcurho,xcuvs2s,ndimcu,help,au,2,0)
     call spls3(crust%EosLog(2,:),crust%EosLog(3,:),crust%ndata,crust%EosLog(2,:),res_data(:),crust%ndata,help,au,2,0)
     crust%EosLog(4,:) = res_data(:) * 10.**crust%EosLog(3,:) / ( 10.**crust%EosLog(2,:) * CST_clight**2 )
     deallocate(help,au,res_data)

     allocate(xcrust_vs2(crust%ndata))
     allocate(xcrust_vs2s(crust%ndata))

     ! calculate the sound velocity by centered difference interpolation: vs2 = dP / drho
     xcrust_vs2(1)=0._pr
     DO iden = 2, crust%ndata - 1
        xy1 =   10.**crust%EosLog(3,iden+1) - 10.**crust%EosLog(3,iden)
        xh1 = ( 10.**crust%EosLog(2,iden+1) - 10.**crust%EosLog(2,iden) ) * CST_clight**2
        xy2 =   10.**crust%EosLog(3,iden)   - 10.**crust%EosLog(3,iden-1)
        xh2 = ( 10.**crust%EosLog(2,iden)   - 10.**crust%EosLog(2,iden-1) ) * CST_clight**2
        xdif = ( xy1 / xh1 - xy2 / xh2 ) / ( xy1 / xh1 + xy2 / xh2 )
        ! check if there is not discontinuity in the derivatives
        if (dabs(xdif).le.0.5.or.iden.lt.12) then
           xcrust_vs2(iden) = 0.5 * xy1 / xh1 + 0.5 * xy2 / xh2
        else
           ! if there is a discontinuity
           ! use a simple linear extrapolation
           !xcuvs2(i) = xcuvs2(i-1) + ( xcuvs2(i-1) - xcuvs2(i-2) ) * &
           !     & ( 10.**xcurho(i) - 10.**xcurho(i-1) ) / &
           !     & ( 10.**xcurho(i-1) - 10.**xcurho(i-2) ) 
           ! use linear regression to extrapolate
           xm = 0.0
           ym = 0.0
           do k = iden - 4, iden - 1, 1
              ym = ym + xcrust_vs2(k)  
              xm = xm + 10.**crust%EosLog(2,k)
           enddo
           ym = ym / 4.0
           xm = xm / 4.0
           xslope = 0.0
           xdeno = 0.0
           do k = iden - 4, iden - 1, 1
              xslope = xslope + ( 10.**crust%EosLog(2,k) - xm ) * (  xcrust_vs2(k) - ym )  
              xdeno = xdeno + ( 10.**crust%EosLog(2,k) - xm )**2
           enddo
           xslope = xslope / xdeno
           xcrust_vs2(iden) = ym + xslope * ( 10.**crust%EosLog(2,iden) - xm )
           !xcrust_vs2(i) = 0.5 * xy1 / xh1 + 0.5 * xy2 / xh2
        endif
        !write(68,*)10.**crust%EosLog(2,iden),xy1/xh1,xy2/xh2,xcrust_vs2(iden),crust%EosLog(4,iden),xdif
     ENDDO
     xcrust_vs2(crust%ndata) = 2. * xcrust_vs2(crust%ndata-1) - xcrust_vs2(crust%ndata-2)

     ! smooth the spline over 3 points to soften the oscillations
     xcrust_vs2s(1) = 0._pr
     DO iden = 2, crust%ndata - 1
        !xcuvs2m(i) = ( xcuvs2s(i-1) + xcuvs2s(i) + xcuvs2s(i+1) ) / 3.0
        ym = crust%EosLog(4,iden)
        xm = 1.0
        do k = iden - 2, iden + 2
           if (k.eq.iden) cycle
           ym = ym + crust%EosLog(4,k) / dabs( crust%EosLog(2,k) - crust%EosLog(2,iden) )**2
           xm = xm + 1.0 / dabs( crust%EosLog(2,k) - crust%EosLog(2,iden) )**2
        enddo
        xcrust_vs2s(iden) = ym / xm
     ENDDO
     xcrust_vs2s(crust%ndata) = 2.*xcrust_vs2s(crust%ndata-1) - xcrust_vs2s(crust%ndata-2)
     
     OPEN(UNIT=67,FILE="nsEos-res/crust-"//trim(nsEos_inp%crust))
     write(67,*)" den (fm-3), rho (g/cm3), pression (Dyn.cm2), vs/c**2"
     DO iden = 1, crust%ndata
        write(67,'(6d15.5)')10.**crust%EosLog(1,iden),10.**crust%EosLog(2,iden),10.**crust%EosLog(3,iden),&
             &crust%EosLog(4,iden),xcrust_vs2(iden),xcrust_vs2s(iden)
     ENDDO
     CLOSE(67)
     !
     ! write the include file for the fast code
     !
     OPEN(UNIT=100,FILE="nsEos-res/inc-"//trim(nsEos_inp%crust))
     write(100,'("crust%ndata    = ",i5)')crust%ndata
     write(100,'("crust%rho_cc   = ",d10.4)')crust%rho_cc ! core-crust energy-density in g.cm3 neutron drip density
     write(100,'("crust%i_cc     = ",i5)')crust%i_cc
     write(100,'("crust%rho_drip = ",d10.4)')crust%rho_drip ! drip energy-density in g.cm3 neutron drip density
     do iden = 1, crust%ndata
        write(100,'("crust%EosLog(1,",i3,")=",f10.4)')iden,crust%EosLog(1,iden)
        write(100,'("crust%EosLog(2,",i3,")=",f10.4)')iden,crust%EosLog(2,iden)
        write(100,'("crust%EosLog(3,",i3,")=",f10.4)')iden,crust%EosLog(3,iden)
        write(100,'("crust%EosLog(4,",i3,")=",d10.4)')iden,crust%EosLog(4,iden)
     enddo
     close(100)
     
     deallocate(xcrust_vs2,xcrust_vs2s)

  endif
  !
  ! =====================
  ! END ROUTINE
  ! =====================
  !

end subroutine nsEos_compute_crust





!
!  CALCULATE THE DENSE MATTER EOS FROM THE META EOS
!

subroutine nsEos_compute_core_metaEos(coef_meta,nsEos_inp,withq,core)

  use metaEosT0Type; use metaEosT0;
  use nsEosType;

  IMPLICIT None
  !
  ! =====================
  ! DEFINE VARIABLES
  ! =====================
  !
  ! ---------------------
  ! INPUT VARIABLES
  ! ---------------------
  !
  TYPE (metaEos_coef), intent(in)   :: coef_meta
  TYPE (nsEos_inputs), intent(in)   :: nsEos_inp
  logical,             intent(in)   :: withq 
  !
  ! ---------------------
  ! OUTPUT VARIABLES
  ! ---------------------
  !
  TYPE (nsEos_core), intent(out) :: core
  !
  ! ---------------------
  ! LOCAL VARIABLES
  ! ---------------------
  !
  TYPE (metaEos_Densities)     :: eosDen
  TYPE (metaEos_Baryons)       :: eosb0
  TYPE (metaEos_Q1)            :: eosq0
  TYPE (metaEos_Leptons)       :: eosl0
  TYPE (metaEos_Photons)       :: eosp0
  TYPE (metaEos_All)           :: eosa0
  !
  integer, dimension(metaEos_neparam) :: eparam
  !
  real (kind=pr) :: xden, xden_step, xdenb, xd, xt, xmuon
  integer        :: iden, ibeta, i
  !
  real (kind=pr) :: drho, dpre, cs, rho0, pre0
  !
  ! =====================
  ! START ROUTINE
  ! =====================
  !
  if (nsEos_lverb) write(*,*)" Core (metaEos):"
  !
  ! Calculate beta-eq
  !
  xd    = 0.8_pr
  xmuon = 0.0_pr
  xt    = 0.0_pr
  rho0 = 0.0_pr
  pre0 = 0.0_pr
  xden = coef_meta%nsat * nsEos_inp%iden_core_min
  xden_step = 0.1 * coef_meta%nsat 
  iden = 1
  xdenb = xden
  core%lstab = .false.
  core%lcaus = .false.
  !
  do while (xdenb.lt.coef_meta%nsat*nsEos_inp%iden_core_max)
!  do iden = 1, nsEos_ncore
     !
!     xden = ( nsEos_inp%iden_core_min + nsEos_inp%iden_core_step * ( iden - 1 ) ) * coef_meta%nsat
     !write(*,*)iden,xden
     !
     !if (( xden / coef_meta%nsat ).gt.nsEos_inp%iden_core_max) exit
     !
     if (coef_meta%imuon.eq.0) then
     ! calculate beta equilibrium with n, p, e only
        call metaeos_T0_beta_npe(xden,xd,coef_meta,withq,eosDen,ibeta)
     else if (coef_meta%imuon.eq.1) then
     ! calculate beta equilibrium with n, p, e, muon
        call metaeos_T0_beta_npemuon(xden,xd,xmuon,coef_meta,withq,eosDen,ibeta)
     else
        stop "coef%imuon badly defined (different from 0 or 1)"
     endif
     if (ibeta.eq.0) then
!        ndenb_max=iden-1
        exit
     endif
     !
     if (dabs(xdenb-eosDen%den_b).ge.0.1*coef_meta%nsat.and.iden.ne.1) xden_step = xden_step * 0.53
     if (dabs(xdenb-eosDen%den_b).lt.0.1*coef_meta%nsat.and.iden.ne.1) xden_step = xden_step * 2.07
     xdenb = eosDen%den_b
     !
     ! prepare the starting values for the next iteration
     !
     xd = eosDen%xd
     xmuon = eosDen%xmuon
     !
     ! Calculate Baryon properties
     !
     call metaeos_T0_baryons(coef_meta,eosDen,.false.,withq,eosb0,eosq0)
     !
     ! Calculate Lepton properties
     !
     call metaeos_T0_leptons(coef_meta,eosDen,.false.,eosl0)
     !
     ! Calculate Photon properties
     !
     call metaeos_T0_Photons(eosDen,.false.,eosp0)
     !
     ! Calculate global properties
     !
     call metaeos_T0_all(coef_meta,eosDen,.false.,eosb0,eosq0,eosl0,eosp0,eosa0)         
     !
     ! Check stability and causality
     !
     drho = eosa0%rho - rho0
     dpre = eosa0%p   - pre0
     cs   = eosa0%cs
     !
     if (xden.gt.0.2.and.drho.lt.0._pr) core%lstab = .true.
     if (xden.gt.0.2.and.dpre.lt.0._pr) core%lstab = .true.
     if (xden.gt.0.14.and.(cs.lt.0._pr.or.cs.gt.1._pr)) core%lcaus = .true.
     ! symmetry energy
!     if (xden.gt.0.2.and.eosb0%esym_b.lt.0._pr) core%lstab = .true.
     !
     if (core%lstab) exit
     if (core%lcaus) exit
     !
     core%EosLog(1,iden)  = log10(eosDen%den_b) 
!     core%EosLog(1,iden)  = log10(xden) 
     core%EosLog(2,iden)  = log10(eosa0%rho * CONV_MeV_fm3_to_g_cm3) 
     core%EosLog(3,iden)  = log10(eosa0%p * CONV_MeV_fm3_to_dyn_cm2)
     core%EosLog(4,iden)  = eosa0%cs
!     write(*,*)iden,xden,eosDen%den_b,(core%EosLog(i,iden),i=1,4)
     !
     rho0 = eosa0%rho
     pre0 = eosa0%p
     !
     xden = xden + xden_step
     iden = iden + 1
     if (iden.ge.nsEos_ncore) then
        write(*,*)"stop: increase nsEos_ncore"
        stop
     endif
     !
  enddo
  !
  core%ndata = iden - 1
!  if (core%ndata.ge.nsEos_ncore) stop "tovEos:setup_coreEos: increase nsEos_ncore"
  core%denb_max = 10**core%EosLog(1,core%ndata)
!  core%denb_max = xdenb! - nsEos_inp%iden_core_step * coef_meta%nsat
  !
  write(*,*)'  exit condition:'
!  if (.not.core%lstab) then
  if (core%lstab) then
     write(*,*)'  stability violated'
  else
     write(*,*)'  stability fulfilled'
  endif
  if (core%lcaus) then
     write(*,*)'  causality violated'
  else
     write(*,*)'  causality fulfilled'
  endif
  write(*,*)'  core%denb_max : ',core%denb_max
!  if (nsEos_lverb.and.(.not.core%lstab)) write(*,*)'  exit from stability condition'
!  if (nsEos_lverb.and.(.not.core%lcaus)) write(*,*)'  exit from causality condition'
  !
  !
  if (nsEos_lverb) write(*,*)'  den_init     : ',10**core%EosLog(1,1)
  if (nsEos_lverb) write(*,*)'  core%ndata   : ',core%ndata
!  if (nsEos_lverb) write(*,*)'  core%denb_max : ',core%denb_max
  !
  ! =====================
  ! END ROUTINE
  ! =====================
  !

end subroutine nsEos_compute_core_metaEos



!
!  CALCULATE THE DENSE MATTER EOS FROM THE POLYTROPIC EOS
!

subroutine nsEos_compute_core_polyEos(coef_poly,nsEos_inp,core)

  use polyEosT0Type; use polyEosT0;
  use nsEosType;

  IMPLICIT None
  !
  ! =====================
  ! DEFINE VARIABLES
  ! =====================
  !
  ! ---------------------
  ! INPUT VARIABLES
  ! ---------------------
  !
  TYPE (polyEos_coef), intent(in)   :: coef_poly
  TYPE (nsEos_inputs), intent(in)   :: nsEos_inp
  !
  ! ---------------------
  ! OUTPUT VARIABLES
  ! ---------------------
  !
  TYPE (nsEos_core), intent(out) :: core
  !
  ! ---------------------
  ! LOCAL VARIABLES
  ! ---------------------
  !
  TYPE (polyEos_eos)       :: core_eos
  !
  integer, dimension(polyEos_neparam) :: eparam
  !
  real (kind=pr) :: xden, xdenb, rho
  integer :: iden, iden_step
  !
  real (kind=pr) :: drho, dpre, cs, rho0, pre0
  !
  ! =====================
  ! START ROUTINE
  ! =====================
  !
  if (nsEos_lverb) write(*,*)" Core (polyEos):"
  !
  iden_step = 0.1
  !
  do iden = 1, nsEos_ncore
     !
     xden = ( nsEos_inp%iden_core_min + iden_step * ( iden - 1 ) ) * coef_poly%nsat
     rho = xden * CST_mnuc2 * CONV_MeV_fm3_to_g_cm3
     !
     if (( xden / coef_poly%nsat ).gt.nsEos_inp%iden_core_max) exit
     !
     ! calculate poly EOS
     !
     call polyeos_T0(rho,coef_poly,core_eos)
     !
     core%EosLog(1,iden)  = log10(xden) 
     core%EosLog(2,iden)  = log10(core_eos%rho * CONV_MeV_fm3_to_g_cm3) 
     core%EosLog(3,iden)  = log10(core_eos%p * CONV_MeV_fm3_to_dyn_cm2)
     core%EosLog(4,iden)  = core_eos%cs2
     !
  enddo
  !
  core%ndata = iden - 1
  if (core%ndata.ge.nsEos_ncore) stop "tovEos:setup_coreEos: increase nsEos_ncore"
  core%denb_max = xdenb - iden_step * coef_poly%nsat
  !
  if (nsEos_lverb) write(*,*)'  den_init     : ',10**core%EosLog(1,1)
  if (nsEos_lverb) write(*,*)'  core%ndata   : ',core%ndata
  if (nsEos_lverb) write(*,*)'  core%denb_max : ',core%denb_max
  !
  ! =====================
  ! END ROUTINE
  ! =====================
  !

end subroutine nsEos_compute_core_polyEos




!



!
!  MERGE THE CRUST AND THE DENSE MATTER EOS
!

subroutine nsEos_match_eos(nsEos_inp,crust,core,nseos)

  IMPLICIT None
  !
  ! =====================
  ! DEFINE VARIABLES
  ! =====================
  !
  ! ---------------------
  ! INPUT VARIABLES
  ! ---------------------
  !
  TYPE (nsEos_inputs), intent(in)  :: nsEos_inp
  TYPE (nsEos_crust),  intent(in)  :: crust
  TYPE (nsEos_core),   intent(in)  :: core
  !
  ! ---------------------
  ! OUTPUT VARIABLES
  ! ---------------------
  !
  TYPE (nsEos_eos), intent(out)  :: nseos
  !
  ! ---------------------
  ! LOCAL VARIABLES
  ! ---------------------
  !
  integer        :: iden, ii
  real (kind=pr) :: cs, xden, xrho, xpre
  real (kind=pr) :: xrho0, xpre0
  real (kind=pr) :: drho, dpre
  !
  ! =====================
  ! START ROUTINE
  ! =====================
  !
  if (nsEos_lverb) write(*,*)" Match crust-core:"
  !
  nseos%eosLog(:,:) = 0._pr
  if (nsEos_inp%core.eq."same") then
     !
     nseos%eosLog(:,1:crust%ndata) = crust%eosLog(:,1:crust%ndata)
     nseos%ndata = crust%ndata
!     nseos%nmax  = crust%ndata
     nseos%denb_max = crust%denb_max
     !
  else if (nsEos_inp%core.eq."meos") then
     !
     nseos%eosLog(:,1:crust%i_cc) = crust%eosLog(:,1:crust%i_cc)
!     do iden = 1, crust%i_cc
        !
!        nseos%eosLog(:,iden) = crust%eosLog(:,iden)
        !
!     enddo
     !
     ii = crust%i_cc
     do iden = 1, core%ndata
        !
        !write(*,*)"tt1:",ii,10**core%EosLog(1,iden)
        if (10**core%EosLog(1,iden).le.0.12) cycle
        ii = ii + 1
        nseos%eosLog(:,ii) = core%eosLog(:,iden)
!        write(*,*)"tt2:",ii,10**core%EosLog(1,iden)
        !
     enddo
     !
     nseos%ndata = ii
     nseos%denb_max = 10**nsEos%EosLog(1,ii)
     nseos%rho_max = 10**nsEos%EosLog(2,ii)
     !
  endif
  nseos%rho_cc   = crust%rho_cc
  nseos%i_cc     = crust%i_cc
  nseos%rho_drip = crust%rho_drip
  !
  ! Check stability and causality of the entier EOS
  ! -> define nmax
  !
  xrho0=0._pr
  xpre0=0._pr
  do iden = 1, nseos%ndata
!     tov_eos%EosLogNba(:) = nsEos%EosLog(1,:)
!     tov_eos%EosLogRho(:) = nsEos%EosLog(2,:)
!     tov_eos%EosLogPre(:) = nsEos%EosLog(3,:)
!     tov_eos%EosLogCs2(:) = nsEos%EosLog(4,:)
     cs = nsEos%EosLog(4,iden)
     xden = 10**nsEos%EosLog(1,iden)
     xrho = 10**nsEos%EosLog(2,iden)
     xpre = 10**nsEos%EosLog(3,iden)
     drho = xrho - xrho0
     dpre = xpre - xpre0
!     write(*,'(i8,5d14.4)')iden,xden,xrho,xpre,cs
     if (drho.lt.0._pr) exit
     if (dpre.lt.0._pr) exit
     if (cs.lt.0._pr.or.cs.gt.1._pr) exit
     xrho0 = xrho
     xpre0 = xpre
  enddo
  nseos%nmax=iden-2
!  write(*,*)'nmax:',nseos%nmax
!  write(*,*)'data:',nseos%ndata,nseos%denb_max
  if (nseos%ndata.gt.nseos%nmax) then
     nseos%ndata = nseos%nmax
     nseos%denb_max = 10**nsEos%EosLog(1,nseos%nmax)
     nseos%rho_max = 10**nsEos%EosLog(2,nseos%nmax)
  endif
!  write(*,*)'nmax:',nseos%nmax
!  write(*,*)'data:',nseos%ndata,nseos%denb_max
!  stop
  !
  if (nsEos_lverb) write(*,'(T4,a32,d10.3,a6)')"Neutron drip density          = ",nseos%rho_drip,"g/cm3"
  if (nsEos_lverb) write(*,'(T4,a32,d10.3,a6)')"Core-crust transition density = ",nseos%rho_cc,"g/cm3"
  if (nsEos_lverb) write(*,'(T4,a32,f10.7,a5)')"Max value of baryon-density   = ",nseos%denb_max,"fm-3"
  if (nsEos_lverb) write(*,'(T4,a32,d10.3,a6)')"Max value of baryon-density   = ",nseos%rho_max,"g/cm3"
  !
  ! =====================
  ! END ROUTINE
  ! =====================
  !
end subroutine nsEos_match_eos



subroutine nsEos_write_eos(nseos)

  use spl;

  IMPLICIT None
  !
  ! =====================
  ! DEFINE VARIABLES
  ! =====================
  !
  ! ---------------------
  ! INPUT VARIABLES
  ! ---------------------
  !
  TYPE (nsEos_eos), intent(in)  :: nseos
  !
  ! ---------------------
  ! LOCAL VARIABLES
  ! ---------------------
  !
  integer :: iden, ii
  !
  real (kind=pr) :: den
  real (kind=pr), dimension(1) :: denLog, res_den, res_rho, res_pre, res_cs2
  real (kind=pr), dimension(:), allocatable :: help,au
  !
  if (nsEos_lverb) write(*,*)" Write file:"
  !
  ! =====================
  ! START ROUTINE
  ! =====================
  !
  if (nsEos_lverb) write(*,*)"  Write: nsEos-res/eos.out"
  open(unit=20,file="nsEos-res/eos.out",status="unknown")
  write(20,*)"# Neutron drip density          = ",nseos%rho_drip,"g/cm3"
  write(20,*)"# Core-crust transition density = ",nseos%rho_cc,"g/cm3"
  write(20,*)"# Max value of baryon-density   = ",nseos%denb_max,"fm-3"
  write(20,*)"# Max value of baryon-density   = ",nseos%rho_max,"g/cm3"
  write(20,*)"# Number of data                = ",nseos%ndata
  write(20,*)"# i, den (fm-3), rho (g/cm3), pression (Dyn/cm2), vs/c**2"
  !
  do iden = 1, nseos%ndata
     !
     write(20,'(i6,4d15.5)')iden,10**nseos%eosLog(1,iden),10**nseos%eosLog(2,iden),10**nseos%eosLog(3,iden),nseos%eosLog(4,iden)
     !
  enddo
  close(20)
  !
  if (nsEos_lverb) write(*,*)"  Write: nsEos-res/eos-lin.out"
  open(unit=20,file="nsEos-res/eos-lin.out",status="unknown")
  write(20,*)"# Neutron drip density          = ",nseos%rho_drip,"g/cm3"
  write(20,*)"# Core-crust transition density = ",nseos%rho_cc,"g/cm3"
  write(20,*)"# Max value of baryon-density   = ",nseos%denb_max,"fm-3"
  write(20,*)"# Number of data                = ",nseos%ndata
  write(20,*)"# i, den (fm-3), rho (g/cm3), pression (Dyn/cm2), vs/c**2"
  !
  allocate(help(nseos%ndata)); allocate(au(nseos%ndata));
  den = 0.01
  do
!     denLog(1) = log10(939*den*1.78266181e12) ! conversion from MeV.fm-3 to g cm-3
     denLog(1) = log10(den)
     call spls3(nseos%EosLog(1,:),nseos%EosLog(1,:),nseos%ndata,denLog,res_den,1,help,au,1,0)
     call spls3(nseos%EosLog(1,:),nseos%EosLog(2,:),nseos%ndata,denLog,res_rho,1,help,au,1,0)
     call spls3(nseos%EosLog(1,:),nseos%EosLog(3,:),nseos%ndata,denLog,res_pre,1,help,au,1,0)
     call spls3(nseos%EosLog(1,:),nseos%EosLog(4,:),nseos%ndata,denLog,res_cs2,1,help,au,1,0)
     write(20,'(4d15.5)')den,10**res_rho(1),10**res_pre(1),res_cs2(1)
     den = den + 0.01
!     if (nsEos_lverb) write(*,*)den,nseos%denb_max,denLog(1),res_den(1),10**res_rho(1)
     if (den.gt.nseos%denb_max) exit
  enddo
  close(20)
  deallocate(help,au)
  !
  if (nsEos_lverb) write(*,*)"  Write: nsEos-res/eos-log.out"
  open(unit=20,file="nsEos-res/eos-log.out",status="unknown")
  write(20,*)"# Neutron drip density          = ",nseos%rho_drip,"g/cm3"
  write(20,*)"# Core-crust transition density = ",nseos%rho_cc,"g/cm3"
  write(20,*)"# Max value of baryon-density   = ",nseos%denb_max,"fm-3"
  write(20,*)"# Number of data                = ",nseos%ndata
  write(20,*)"# i, den (fm-3), rho (g/cm3), pression (Dyn/cm2), vs/c**2"
  !
  allocate(help(nseos%ndata)); allocate(au(nseos%ndata));
!  denLog(1) = -2.0*939*1.78266181e12 ! conversion from fm-3 to g cm-3
  denLog(1) = -12
  do
     den = 10**denLog(1)
     if (den.gt.nseos%denb_max) exit
     call spls3(nseos%EosLog(1,:),nseos%EosLog(2,:),nseos%ndata,denLog(1),res_rho(1),1,help,au,1,0)
     call spls3(nseos%EosLog(1,:),nseos%EosLog(3,:),nseos%ndata,denLog(1),res_pre(1),1,help,au,1,0)
     call spls3(nseos%EosLog(1,:),nseos%EosLog(4,:),nseos%ndata,denLog(1),res_cs2(1),1,help,au,1,0)


     write(20,'(4d15.5)')den,10**res_rho(1),10**res_pre(1),res_cs2(1)
!     if (nsEos_lverb) write(*,*)den,nseos%denb_max,10**res_den(1),10**res_rho(1)
     denLog(1) = denLog(1) + 0.1
  enddo
  close(20)
  deallocate(help,au)
  !
  ! =====================
  ! END ROUTINE
  ! =====================
  !
end subroutine nsEos_write_eos



subroutine nsEos_write_table_eos(nseos,coef_meta)

  use metaEosT0Type;
  use spl;

  IMPLICIT None
  !
  ! =====================
  ! DEFINE VARIABLES
  ! =====================
  !
  ! ---------------------
  ! INPUT VARIABLES
  ! ---------------------
  !
  TYPE (nsEos_eos), intent(in)      :: nseos
  TYPE (metaEos_coef), intent(in)   :: coef_meta
  !
  ! ---------------------
  ! LOCAL VARIABLES
  ! ---------------------
  !
  integer :: iden, ii
  !
  real (kind=pr) :: den
  real (kind=pr), dimension(1) :: denLog, res_den, res_rho, res_pre, res_cs2
  real (kind=pr), dimension(:), allocatable :: help,au
  !
  ! =====================
  ! START ROUTINE
  ! =====================
  !
!  open(unit=20,file="nsEos-res/eos.out",status="unknown")
!  write(20,*)"# Neutron drip density          = ",nseos%rho_drip,"g/cm3"
!  write(20,*)"# Core-crust transition density = ",nseos%rho_cc,"g/cm3"
!  write(20,*)"# Max value of baryon-density   = ",nseos%denb_max,"fm-3"
!  write(20,*)"# Number of data                = ",nseos%ndata
!  write(20,*)"# i, den (rho-3), rho (g/cm3), pression (Dyn.cm2)"
!  !
!  do iden = 1, nseos%ndata
!     !
!     write(20,*)iden,10**nseos%eosLog(1,iden),10**nseos%eosLog(2,iden),10**nseos%eosLog(3,iden),nseos%eosLog(4,iden)
!     !
!  enddo
!  close(20)
  !
  open(unit=20,file="nsEos-res/table.out",status="unknown")
  write(20,'(a6,13i7,2i3)')"# MM:",coef_meta%eparam(:)
  write(20,*)"# Neutron drip density          = ",nseos%rho_drip,"g/cm3"
  write(20,*)"# Core-crust transition density = ",nseos%rho_cc,"g/cm3"
  write(20,*)"# Max value of baryon-density   = ",nseos%denb_max,"fm-3"
  write(20,*)"# i, den (fm-3), rho (g/cm3), pression (Dyn/cm2), vs/c**2"
  !
  allocate(help(nseos%ndata)); allocate(au(nseos%ndata));
  denLog(1) = -12.0
  iden = 1
  do
     den = 10**denLog(1)
     if (den.gt.nseos%denb_max) exit
     call spls3(nseos%EosLog(1,:),nseos%EosLog(2,:),nseos%ndata,denLog(1),res_rho(1),1,help,au,1,0)
     call spls3(nseos%EosLog(1,:),nseos%EosLog(3,:),nseos%ndata,denLog(1),res_pre(1),1,help,au,1,0)
     call spls3(nseos%EosLog(1,:),nseos%EosLog(4,:),nseos%ndata,denLog(1),res_cs2(1),1,help,au,1,0)
     write(20,'(i6,4d15.5)')iden,den,10**res_rho(1),10**res_pre(1),res_cs2(1)
!     if (nsEos_lverb) write(*,*)den,nseos%denb_max,10**res_den(1),10**res_rho(1)
     denLog(1) = denLog(1) + 0.1
     iden = iden + 1
  enddo
  close(20)
  deallocate(help,au)
  !
  ! =====================
  ! END ROUTINE
  ! =====================
  !
end subroutine nsEos_write_table_eos


end MODULE nsEosMod
