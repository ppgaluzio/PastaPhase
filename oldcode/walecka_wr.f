!EoS no sigma-omega(com autointeracao) + propriedades da materia nuclear
!(energia de ligacao, m_eff, a_sym)

	  program solido5
	  IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	  logical check
	  REAL*8 Y(1), FVEC(1)
	  EXTERNAL F1,F2,F3
	  common/cte/RMA,PI,PI2,Msig
	  COMMON/r/RHO,RKF
	  COMMON/TEST/FVEC1
	  common/gs/SIGMA,gs
	  common/rgs/RGSIGMA

	  open(35,file='tovqhd2.dat',status='unknown')
C----------------------------------------------------
C     Dados de input (if rb and c =0, no self interaction is taken)
c     GM1 set of parameters
	  RMA=939.D0				!MeV
	  hc=197.33d0				!fm*MeV
	  PI=DACOS(-1.D0)
	  PI2=PI*PI
	  GO=10.608D0				!Parametro g_omega without units
	  GS=9.5684D0				!Parametro g_sigma without units
	  Msig = 550.d0				!sigma mass in MeV
	  Momeg = 783.d0			!omega mass in MeV
	  rb=0.						! 0.2948d0/100.d0		        !const. acopl. (glendenning)
	  c=0.						!  -0.1071d0/100.d0  		!const. acopl. (glendenning)


	  RNBINF=0.05D0				!fm^-3
	  RNBSUP=1.D0				!fm^-3
	  NPOINT=200

	  RNBINF=RNBINF*(hc)**3		!MeV^3
	  RNBSUP=RNBSUP*(hc)**3		!MeV^3
	  DNB=(RNBSUP-RNBINF)/(NPOINT-1) !MeV^3

c-----------------------------------------
	  DO i=1,NPOINT
		 RHO=(RNBINF+(i-1)*DNB)	!DENSIDADE EM MeV^3

		 RKF=(1.5d0*PI2*RHO)**(1.D0/3.D0) !fermi momentum em MeV


c-----------------------------------------------------
		 Y(1)=0.1d0

		 CALL BROYDN(Y,1,CHECK)
		 IF(CHECK)THEN
			WRITE(*,*) "Não há raízes"
		 END IF
!     write(6,*)'Debg1:gm,pi2,b,c',gm,pi2,rb,c
c-------------mesons sigma e omega---------------

		 SIGMA=Y(1)				!em MeV

		 RGOMEGA=GO*RHO/Momeg**2 !meson omega em MeV

C-----------------------------------------------------------------
!     WRITE(50,*) RHO/(hc)**3,RKF/(hc),FVEC1


!     WRITE(51,*) RGSIGMA/RMA,RGOMEGA/RMA,(RHO/(hc**3))/0.15D0

c----------------------------------------------------
c     calculo pressao e energia

		 CALL GAUSS(F2,0.D0,RKF,10,RE2,II)


		 ENER=0.5D0*(Msig*SIGMA)**2.D0 + 0.5d0*(Momeg*RGOMEGA)**2.d0
	1		  +(2.D0/PI2)*RE2
!     +
!  &  (1/3.D0)*rb*RMA*(RGSIGMA**3.D0)+			   !MeV^4
!  &	(1/4.D0)*c*(RGSIGMA**4.D0)+(2.D0/PI2)*RE2

		 CALL GAUSS(F3,0.D0,RKF,10,RE3,II)

		 PRES=-0.5D0*(Msig*SIGMA)**2.D0 + 0.5d0*(Momeg*RGOMEGA)**2.d0
	1		  +(2.D0/(3.D0*PI2))*RE3
!-
!   &  (1/3.D0)*rb*RMA*(RGSIGMA**3.D0)-
!   &	(1/4.D0)*c*(RGSIGMA**4.D0)+(2.D0/(3.D0*PI2))*RE3	   !Mev^4

c-------saidas pressão e energia em MeV.fm^-3---------
		 energ=ENER/(hc**3.d0)
		 press=PRES/(hc**3.d0)


!     WRITE(52,*) energ,press

c----------Saida tov---------------------------------

!     write(35,*)RHO/(hc**3.d0),ENER/(hc**4.d0),PRES/(hc**4.d0)

c-------energia de ligação---------------------------

		 rmeff=RMA-RGSIGMA

		 B=(ENER/RHO)-RMA		!MeV

		 WRITE(54,*)RHO/(hc**3.d0),rmeff/rma,energ,press

c-------calulo energia de simetria ------------


!     Esym=(RKF**2.d0)/(6*(RKF*RKF+rmeff*rmeff)**0.5d0)   !MeV

!     write(55,*)RHO/(hc)**3,Esym

c-----------------------------------------------------------------

	  ENDDO
	  close(35)
	  END
!-----------------------------------------------------
	  SUBROUTINE funcv(N,Y,FVEC)
	  IMPLICIT double precision(A-H,O-Z)
	  REAL*8 Y(1), FVEC(1)
	  EXTERNAL F1,F2,F3
	  common/cte/RMA,PI,PI2,Msig
	  common/r/RHO,RKF
	  COMMON/TEST/FVEC1
	  common/gs/SIGMA,gs
	  SIGMA=Y(1)

	  CALL GAUSS(F1,0.D0,RKF,10,RE1,II)

	  gm2=(GS/Msig)**2

	  FVEC(1)=GS*SIGMA-gm2*(2.D0/PI2)*RE1
!     &  -rb*RMA*(GSIGMA**2.D0)
!     &  -c*(GSIGMA**3.D0))

	  FVEC1=FVEC(1)
	  write(*,*)gS,SIGMA,gm2,PI2,RE1,Msig
	  RETURN

	  END
C--------------------------------------------------------
      FUNCTION F1(X)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      common/cte/RMA,PI,PI2,Msig
      common/gs/SIGMA,gs

      E=DSQRT(X*X+(RMA-GS*SIGMA)**2.d0)
      F1=((X**2.d0)*(RMA-GS*SIGMA))/E


      RETURN
      END
C------------------------------------------------------------------
      FUNCTION F2(X)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      common/cte/RMA,PI,PI2,Msig
      common/gs/SIGMA,gs
!     common/rgs/RGSIGMA
      E=DSQRT(X*X+(RMA-GS*SIGMA)**2.D0)
      F2=X*X*E

      RETURN
      END
C--------------------------------------------------------
      FUNCTION F3(X)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      common/cte/RMA,PI,PI2,Msig
      common/gs/SIGMA,gs
!     common/rgs/RGSIGMA
      E=DSQRT(X*X+(RMA-GS*SIGMA)**2.D0)
      F3=(X**4.d0)/E

      RETURN
      END
c----------------------------------------------------
C-----Gauss
      SUBROUTINE GAUSS(F,UG,OG,NN,FXINT,II)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION WG(10),ZG(10)
      DATA WG/0.6667134430869 D-01,
     C     0.1494513491506D00,
     C     0.2190863625160D00,
     C     0.2692667193100D00,
     C     0.2955242247148D00,
     C     0.2955242247148D00,
     C     0.2692667193100D00,
     C     0.2190863625160D00,
     C     0.1494513491506D00,
     C     0.6667134430869D-01/
      DATA ZG/-0.9739065285172D00,
     C     -0.8650633666890D00,
     C     -0.6794095682990D00,
     C     -0.4333953941292D00,
     C     -0.1488743389816D00,
     C     +0.1488743389816D00,
     C     +0.4333953941292D00,
     C     +0.6794095682990D00,
     C     +0.8650633666890D00,
     C     +0.9739065285172D00/

      FXINT=0.D0
      HH=(OG-UG)/DBLE(FLOAT(NN))
      U=UG
      O=U+HH
      KK=1
 24   OU=O+U
      RI=0.D0
      DO 26 I=1,10
         X=0.5D0*(ZG(I)*HH+OU)
         FUNCAO=F(X)
 26      RI=RI+WG(I)*FUNCAO
         FXINT=RI*HH/2.D0+FXINT
         KK=KK+1
         IF(KK-NN)28,28,9999
 28      U=O
         O=O+HH
         GO TO 24
 9999    RETURN
         END
c------------------------------------------

	  include 'brodyn.f'
