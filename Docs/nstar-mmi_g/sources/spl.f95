
module spl

  use acc

  implicit none

contains



!     ----------------------------------------------------------------
!     linear interpolation and extrapolation

!     ----------------------------------------------------------------
!     linear interpolation and extrapolation

  subroutine intp_1D(xvec1,ndim1,fvec,xpoint,isg,rf)
  ! isg = 1 if xvec1 is increasing
  ! isg =-1 if xvec1 is decreasing

    use acc

    implicit none

    integer, intent(in)              :: ndim1, isg
    real (kind=pr), dimension(ndim1) :: xvec1
    real (kind=pr), dimension(ndim1) :: fvec
    real (kind=pr), dimension(1) :: xpoint
    real (kind=pr), dimension(1) :: x0, xh, xp
    integer :: i1
    real (kind=pr) :: rf

    if (isg.eq.1.and.xpoint(1).lt.xvec1(1)) then
      i1=1
      write(91,*)'1D 1 min ',xvec1(1),' not <',xpoint(1)
    else if (isg.eq.1.and.xpoint(1).gt.xvec1(ndim1)) then
      i1=ndim1-1
      write(91,*)'1D 1 max ',xpoint(1),' not <',xvec1(ndim1)
    else if (isg.eq.-1.and.xpoint(1).gt.xvec1(1)) then
      i1=1
      write(91,*)'1D 1 min ',xvec1(1),' not >',xpoint(1)
    else if (isg.eq.-1.and.xpoint(1).lt.xvec1(ndim1)) then
      i1=ndim1-1
      write(91,*)'1D 1 max ',xpoint(1),' not >',xvec1(ndim1)
    else
      i1=1
      do while (isg*(xvec1(i1)-xpoint(1)).lt.0.d0.and.i1.lt.ndim1)
        i1=i1+1
      enddo
      i1=i1-1
    endif
    x0(1)=xvec1(i1)
    xh(1)=xvec1(i1+1)-xvec1(i1)
    xp(1)=1._pr
    if (xh(1).ne.0._pr) xp(1)=(xpoint(1)-x0(1))/xh(1)

    rf = (1. - xp(1)) * fvec(i1) + xp(1) * fvec(i1+1)

  end subroutine intp_1D


  subroutine intp_2D(xvec1,ndim1,xvec2,ndim2,fvec,xpoint,rf)

    use acc

    implicit none

!    integer, parameter :: ndim1=20, ndim2=25
    integer, intent(in)              :: ndim1, ndim2
    real (kind=pr), dimension(ndim1) :: xvec1
    real (kind=pr), dimension(ndim2) :: xvec2
    real (kind=pr), dimension(ndim1,ndim2) :: fvec
    real (kind=pr), dimension(2) :: xpoint
    real (kind=pr), dimension(2) :: x0, xh, xp
    integer :: i1, i2, j1, j2
    real (kind=pr) :: rf

!     test:
!      write(*,*)(xpoint(i),i=1,3)
!      write(*,*)(xvec1(i1),i1=1,ndim1)
!      write(*,*)(xvec2(i2),i2=1,ndim2)
!     end test

    if (xpoint(1).lt.xvec1(1)) then
      i1=1
      write(91,*)'2D 1 min ',xvec1(1),' not <',xpoint(1)
    else if (xpoint(1).gt.xvec1(ndim1)) then
      i1=ndim1-1
      write(91,*)'2D 1 max ',xpoint(1),' not <',xvec1(ndim1)
    else
      i1=1
      do while ((xvec1(i1)-xpoint(1)).lt.0.d0.and.i1.lt.ndim1)
        i1=i1+1
      enddo
      i1=i1-1
    endif
    x0(1)=xvec1(i1)
    xh(1)=xvec1(i1+1)-xvec1(i1)
    xp(1)=1._pr
    if (xh(1).ne.0._pr) xp(1)=(xpoint(1)-x0(1))/xh(1)
    if (xpoint(2).lt.xvec2(1)) then
      i2=1
      write(92,*)'2D 2 min ',xvec2(1),' not <',xpoint(2)
    else if (xpoint(2).gt.xvec2(ndim2)) then
      i2=ndim2-1
      write(92,*)'2D 2 max ',xpoint(2),' not <',xvec2(ndim2)
    else
      i2=1
      do while ((xvec2(i2)-xpoint(2)).lt.0._pr.and.i2.lt.ndim2)
        i2=i2+1
      enddo
      i2=i2-1
    endif
    x0(2)=xvec2(i2)
    xh(2)=xvec2(i2+1)-xvec2(i2)
    xp(2)=1._pr
    if (xh(2).ne.0._pr) xp(2)=(xpoint(2)-x0(2))/xh(2)

!    write(*,*)'x0',x0(:)
!    write(*,*)'xh',xh(:)
!    write(*,*)'xp',xp(:)

    rf=0._pr
    do j1=0,1
      do j2=0,1
        rf=rf+(1-j1-(-1._pr)**j1*xp(1))*(1-j2-(-1._pr)**j2*xp(2))* &
 &         fvec(i1+j1,i2+j2)
!        write(*,*)xvec1(i1+j1),xvec2(i2+j2),fvec(i1+j1,i2+j2)
      enddo
    enddo

  end subroutine intp_2D


!     ----------------------------------------------------------------
  subroutine intp_3D(xvec1,ndim1,xvec2,ndim2,xvec3,ndim3,fvec,xpoint,rf)

    use acc

    implicit none

!    integer, parameter :: ndim1=20, ndim2=25, ndim3=15
    integer, intent(in)              :: ndim1, ndim2, ndim3
    real (kind=pr), dimension(ndim1) :: xvec1
    real (kind=pr), dimension(ndim2) :: xvec2
    real (kind=pr), dimension(ndim3) :: xvec3
    real (kind=pr), dimension(ndim1,ndim2,ndim3) :: fvec
    real (kind=pr), dimension(3) :: xpoint
    real (kind=pr), dimension(3) :: x0, xh, xp
    integer :: i1, i2, i3, j1, j2, j3
    real (kind=pr) :: rf

!     test:
!      write(*,*)(xpoint(i),i=1,3)
!      write(*,*)(xvec1(i1),i1=1,ndim1)
!      write(*,*)(xvec3(i3),i3=1,ndim3)
!     end test

    if (xpoint(1).lt.xvec1(1)) then
      i1=1
      write(91,*)'3D 1 min ',xvec1(1),' not <',xpoint(1)
    else if (xpoint(1).gt.xvec1(ndim1)) then
      i1=ndim1-1
      write(91,*)'3D 1 max ',xpoint(1),' not <',xvec1(ndim1)
    else
      i1=1
      do while ((xvec1(i1)-xpoint(1)).lt.0.d0.and.i1.lt.ndim1)
        i1=i1+1
      enddo
      i1=i1-1
    endif
    x0(1)=xvec1(i1)
    xh(1)=xvec1(i1+1)-xvec1(i1)
    xp(1)=1._pr
    if (xh(1).ne.0._pr) xp(1)=(xpoint(1)-x0(1))/xh(1)
    if (xpoint(2).lt.xvec2(1)) then
      i2=1
      write(92,*)'3D 2 min ',xvec2(1),' not <',xpoint(2)
    else if (xpoint(2).gt.xvec2(ndim2)) then
      i2=ndim2-1
      write(92,*)'3D 2 max ',xpoint(2),' not <',xvec2(ndim2)
    else
      i2=1
      do while ((xvec2(i2)-xpoint(2)).lt.0._pr.and.i2.lt.ndim2)
        i2=i2+1
      enddo
      i2=i2-1
    endif
    x0(2)=xvec2(i2)
    xh(2)=xvec2(i2+1)-xvec2(i2)
    xp(2)=1._pr
    if (xh(2).ne.0.d0) xp(2)=(xpoint(2)-x0(2))/xh(2)
    if (xpoint(3).lt.xvec3(1)) then
      i3=1
      write(93,*)'3D 3 min ',xvec3(1),' not <',xpoint(3)
    else if (xpoint(3).gt.xvec3(ndim3)) then
      i3=ndim3-1
      write(93,*)'3D 3 max ',xpoint(3),' not <',xvec3(ndim3)
    else
      i3=1
      do while ((xvec3(i3)-xpoint(3)).lt.0._pr.and.i3.lt.ndim3)
        i3=i3+1
      enddo
      i3=i3-1
    endif
    x0(3)=xvec3(i3)
    xh(3)=xvec3(i3+1)-xvec3(i3)
    xp(3)=1._pr
    if (xh(3).ne.0.d0) xp(3)=(xpoint(3)-x0(3))/xh(3)

    write(*,*)'x0',x0(:)
    write(*,*)'xh',xh(:)
    write(*,*)'xp',xp(:)

    rf=0._pr
    do j1=0,1
      do j2=0,1
        do j3=0,1
          rf=rf+(1-j1-(-1._pr)**j1*xp(1))* &
          (1-j2-(-1._pr)**j2*xp(2))* &
          (1-j3-(-1._pr)**j3*xp(3))* &
          fvec(i1+j1,i2+j2,i3+j3)
        enddo
      enddo
    enddo

  end subroutine intp_3D


!     ----------------------------------------------------------------
  subroutine intp_4D(xvec1,ndim1,xvec2,ndim2,xvec3,ndim3,xvec4,ndim4,fvec,xpoint,rf)

    use acc

    implicit none

!    integer, parameter :: ndim1=20, ndim2=25, ndim3=15, ndim4=20
    integer, intent(in)              :: ndim1, ndim2, ndim3, ndim4
    real (kind=pr), dimension(ndim1) :: xvec1
    real (kind=pr), dimension(ndim2) :: xvec2
    real (kind=pr), dimension(ndim3) :: xvec3
    real (kind=pr), dimension(ndim4) :: xvec4
    real (kind=pr), dimension(ndim1,ndim2,ndim3,ndim4) :: fvec
    real (kind=pr), dimension(4) :: xpoint
    real (kind=pr), dimension(4) :: x0, xh, xp
    integer :: i1, i2, i3, i4, j1, j2, j3, j4
    real (kind=pr) :: rf

!     test:
!      write(*,*)(xpoint(i),i=1,3)
!      write(*,*)(xvec1(i1),i1=1,ndim1)
!      write(*,*)(xvec3(i3),i3=1,ndim3)
!     end test

    if (xpoint(1).lt.xvec1(1)) then
      i1=1
      write(91,*)'4D 1 min ',xvec1(1),' not <',xpoint(1)
    else if (xpoint(1).gt.xvec1(ndim1)) then
      i1=ndim1-1
      write(91,*)'4D 1 max ',xpoint(1),' not <',xvec1(ndim1)
    else
      i1=1
      do while ((xvec1(i1)-xpoint(1)).lt.0.d0.and.i1.lt.ndim1)
        i1=i1+1
      enddo
      i1=i1-1
    endif
    x0(1)=xvec1(i1)
    xh(1)=xvec1(i1+1)-xvec1(i1)
    xp(1)=1._pr
    if (xh(1).ne.0.d0) xp(1)=(xpoint(1)-x0(1))/xh(1)
    if (xpoint(2).lt.xvec2(1)) then
      i2=1
      write(92,*)'4D 2 min ',xvec2(1),' not <',xpoint(2)
    else if (xpoint(2).gt.xvec2(ndim2)) then
      i2=ndim2-1
      write(92,*)'4D 2 max ',xpoint(2),' not <',xvec2(ndim2)
    else
      i2=1
      do while ((xvec2(i2)-xpoint(2)).lt.0._pr.and.i2.lt.ndim2)
        i2=i2+1
      enddo
      i2=i2-1
    endif
    x0(2)=xvec2(i2)
    xh(2)=xvec2(i2+1)-xvec2(i2)
    xp(2)=1._pr
    if (xh(2).ne.0._pr) xp(2)=(xpoint(2)-x0(2))/xh(2)
    if (xpoint(3).lt.xvec3(1)) then
      i3=1
      write(93,*)'4D 3 min ',xvec3(1),' not <',xpoint(3)
    else if (xpoint(3).gt.xvec3(ndim3)) then
      i3=ndim3-1
      write(93,*)'4D 3 max ',xpoint(3),' not <',xvec3(ndim3)
    else
      i3=1
      do while ((xvec3(i3)-xpoint(3)).lt.0._pr.and.i3.lt.ndim3)
        i3=i3+1
      enddo
      i3=i3-1
    endif
    x0(3)=xvec3(i3)
    xh(3)=xvec3(i3+1)-xvec3(i3)
    xp(3)=1._pr
    if (xh(3).ne.0._pr) xp(3)=(xpoint(3)-x0(3))/xh(3)

    if (xpoint(4).lt.xvec4(1)) then
      i4=1
      write(94,*)'4D 4 min ',xvec4(1),' not <',xpoint(4)
    else if (xpoint(4).gt.xvec4(ndim4)) then
      i4=ndim4-1
      write(94,*)'4D 4 max ',xpoint(4),' not <',xvec4(ndim4)
    else
      i4=1
      do while ((xvec4(i4)-xpoint(4)).lt.0._pr.and.i4.lt.ndim4)
        i4=i4+1
      enddo
      i4=i4-1
    endif
    x0(4)=xvec4(i4)
    xh(4)=xvec4(i4+1)-xvec4(i4)
    xp(4)=1._pr
    if (xh(4).ne.0.d0) xp(4)=(xpoint(4)-x0(4))/xh(4)

!         write(*,*)'x0',(x0(i),i=1,4)
!         write(*,*)'xh',(xh(i),i=1,4)
!         write(*,*)'xp',(xp(i),i=1,4)

    rf=0._pr
    do j1=0,1
      do j2=0,1
        do j3=0,1
          do j4=0,1
            rf=rf+(1-j1-(-1._pr)**j1*xp(1))* &
            (1-j2-(-1._pr)**j2*xp(2))* &
            (1-j3-(-1._pr)**j3*xp(3))* &
            (1-j4-(-1._pr)**j4*xp(4))* &
            fvec(i1+j1,i2+j2,i3+j3,i4+j4)
          enddo
        enddo
      enddo
    enddo

  end subroutine intp_4D




SUBROUTINE SPLS3(X,Y,N,XI,FI,M,Q,AU,IGO,ISPL)
!
!     ******************************************************************
!
!     CUBIC SPLINE, STARTING WITH ZERO SECOND DERIVATIVES AT THE
!     BOUNDARIES OF THE APPROXIMATION INTERVAL.
!
!     IGO = 0      BUILD UP SPLINE ONLY.
!     IGO = 1      BUILD UP SPLINE AND INTERPOLATE FOR M POINTS.
!     IGO = 2      BUILD UP SPLINE AND COMPUTE DERIVATIVES AT M POINTS.
!
!     ISPL = 0     BUILD UP SPLINE
!     ISPL = 1     INTERPOLE ONLY
!
!C     REAL*8 VERSION.        J.GALONSKA, 15.12.1971
!     ******************************************************************
  use acc
  IMPLICIT NONE
  ! inputs
  integer                     :: ISPL,N,M,IGO
  real (kind=pr),dimension(M) :: XI,FI
  real (kind=pr),dimension(N) :: X,Y
  real(kind=pr), dimension(N) :: Q,AU
  !outputs
  integer        :: NN,K,NN2,KK,J,M1,M2,M3
  real (kind=pr) :: ZERO=0.,THREE=3.,SIX=6.,FACT=0.1666666666667
  real (kind=pr) :: HK,HX,HI,HI2
  real (kind=pr) :: YSAVE,AUX,DIVQ,YK,DIJ,DIM1J,DIJ3,DIM1J3

  print*, "N=", N
  print*, "size", size(au), size(Q)

  IF (ISPL.EQ.0) THEN
     AU(1) = ZERO
     AU(N) = ZERO
     Q(1) = ZERO
     HK = X(2) - X(1)
     YSAVE = (Y(2)-Y(1)) / HK
     AUX = ZERO
     NN = N - 1
     DO K = 2,NN
        HX = X(K+1) - X(K-1)
        DIVQ = (HK*Q(K-1)+HX+HX)
        HK = X(K+1) - X(K)
        YK = (Y(K+1)-Y(K)) / HK
        Q(K) = - HK / DIVQ
        AU(K) = (SIX*(YK-YSAVE)-AUX) / DIVQ
        YSAVE = YK
        AUX = AU(K) * HK
     ENDDO
     NN2 = NN + 2
     DO KK = 2,NN
        K = NN2 - KK
        AU(K) = Q(K) * AU(K+1) + AU(K)
        print*, "K = ", K
     ENDDO
     IF (IGO.EQ.0)  RETURN
  ENDIF

!
!     ******************************************************************
!
!     INTERPOLATION OR COMPUTATION OF DERIVATIVES.
!
!     IGO = 1      INTERPOLATE FOR M POINTS.
!     IGO = 2      COMPUTE DERIVATIVES AT M POINTS.
!
!     ******************************************************************
!
  DO J = 1,M
     IF (X(1).GT.XI(J))  THEN
        M1=1
        M2=2
     ELSEIF (XI(J).GT.X(N))  THEN
        M1=N-1
        M2=N
     ELSE
        M1 = 1
        M2 = N
50      M3 = (M2+M1)/2
        IF (XI(J).GE.X(M3)) THEN
           M1 = M3
        ELSE
           M2 = M3
        ENDIF
        IF (M1+1-M2.NE.0)  GO TO 50
     ENDIF
     DIJ = X(M2) - XI(J)
     DIM1J = X(M1) - XI(J)
     HI = X(M2) - X(M1)
     HI2 = HI * HI
     IF (IGO.GE.2) THEN
        FI(J) = FACT * (THREE*(AU(M2)*DIM1J*DIM1J-AU(M1)*DIJ*DIJ) &
 &                      -SIX*(Y(M1)-Y(M2))+HI2*(AU(M1)-AU(M2))) / HI
     ELSE
        DIJ3 = DIJ * DIJ * DIJ
        DIM1J3 = DIM1J * DIM1J * DIM1J
        FI(J) = FACT * (AU(M1)*DIJ3-AU(M2)*DIM1J3+(SIX*Y(M1)-HI2*AU(M1)) &
 &           *DIJ-(SIX*Y(M2)-HI2*AU(M2))*DIM1J) / HI
     ENDIF

     print*, "m1 = ", m1, "m2 = ", m2
  ENDDO

END SUBROUTINE SPLS3


end module spl