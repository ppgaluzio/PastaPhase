!
!
! =============================================================  
!
!                        Module NLE
!
! =============================================================
!
!
!
! --------------------------NOTES------------------------------
!
!   Find the solution of non-linear multidimensional equations.
!
! ---------------------Public subroutines----------------------
!   nles(fcn,x,stp,eps,n):
!     Find the solution of non-linear equations using the simple
!     sequant algorithm. Very pedestrian but robust algorithm.
!         n: number of equations (dimension of the system to solve)
!         x: variables of the minimization
!         stp: step
!         eps: accuracy to reach
!     The system of multidimensional equations is writen in fcn(x,y,n) 
!     which is defined in the module nle_fcn.
!
!   nles2(fcn,x,stp,eps,n):
!     Find the solution of non-linear equations moving in the multidimensional
!     space at each iteration. Still pedestrian, can go faster but less
!     robust than nles.
!         n: number of equations (dimension of the system to solve)
!         x: variables of the minimization
!         stp: step
!         eps: accuracy to reach
!     The system of multidimensional equations is writen in fcn(x,y,n) 
!     which is defined in the module nle_fcn.
!
! --------------------Private subroutines----------------------
!
!   function norm(y,n)
!
!
! -------------------------------------------------------------
! Code by Jerome Margueron
! Date: 10-10-2009
! Latest revision - 19-10-2011
! -------------------------------------------------------------


module NLE

  use acc

  implicit none

  public

!  real (kind=pr), private :: norm
  
contains


  subroutine nles(fcn,x,stp,eps,n,nor,ifail)

    use acc; use cst;

    implicit none

! interface to recognize crust souboutine and function ?
!    interface
!       subroutine crust_sna(x,y,n)
!         use acc;
!         implicit none
!         real (kind=pr), dimension(n), intent(in)  :: x
!         integer, intent(in)                       :: n
!         real (kind=pr), dimension(n), intent(out) :: y
!       end subroutine crust_sna
!    end interface

    interface
       subroutine fcn(x,y,n)
         use acc
         implicit none
         real (kind=pr), dimension(n), intent(in)  :: x
         integer, intent(in)                       :: n
         real (kind=pr), dimension(n), intent(out) :: y
       end subroutine fcn
    end interface

    !input variables
    integer, intent(in)          :: n
    real (kind=pr), intent(in)   :: eps
    ! output variables
    real (kind=pr), intent(out)  :: nor
    integer, intent(out)         :: ifail
    ! input and output variables
    real (kind=pr), dimension(n), intent(inout) :: x,stp
    ! local variables
    integer                      :: i, j, k
    integer                      :: in, in2, icount, ii
    real (kind=pr), dimension(n) :: y0, y1, y2
    real (kind=pr)               :: nor0, nor1
    logical                      :: test

    ifail = 0

    nor1 = 10._pr

    in=0

    do while (nor1 > eps)

       in = in+1

       ! loop over the variables x(1), x(2), ...
       loopi : do i=1,n

          ! the variables stored in x are given in the argument of the
          ! routine for the first iteration, and are incremented in this
          ! loop to reached the zeros of y.
          call fcn(x,y0,n)

          ! calculate the norm of y0
          nor0 = norm(y0,n)

          ii = 1
          test = .true.
          icount = 1
          ! this loop increaments the variable x(i) while y(i) decreases
          ! in absolute value
          do while (test)

             x(i)=x(i)+stp(i)
             call fcn(x,y1,n)
             nor1 = norm(y1,n)

             ! if y(i) increases, change the sign of the step, decreases it,
             ! and go to the next variable
             if ( nor1 > nor0 ) then
                stp(i)=-stp(i)/3._pr
                test=.false.
             endif

             ! acceleration of the step to avoid being stuck into the wild
             if (ii.eq.4) stp(i)=3._pr*stp(i)
             if (ii > 10) stp(i)=(ii-9._pr)*stp(i)

             y0=y1
             nor0=nor1
             ii=ii+1

          enddo

       enddo loopi

       if (in.eq.2000) then
          write(*,*)'Too much iteration, stop 1'
          ifail=1
       endif

    enddo

    nor=nor1

  end subroutine nles




  subroutine nles2(fcn,x,stp,eps,n)

    use acc; use cst;

    implicit none

    interface
       subroutine fcn(x,y,n)
         use acc
         implicit none
         real (kind=pr), dimension(n), intent(in)  :: x
         integer, intent(in)                       :: n
         real (kind=pr), dimension(n), intent(out) :: y
       end subroutine fcn
    end interface

    !input variables
    integer, intent(in)          :: n
    real (kind=pr), intent(in)   :: eps
    ! input and output variables
    real (kind=pr), dimension(n), intent(inout) :: x,stp
    ! local variables
    integer                      :: i, j, k, in, in2, icount
    real (kind=pr), dimension(n) :: y0, y1, y2
    real (kind=pr)               :: nor0, nor1
    logical                      :: test

    in=0

    call fcn(x,y1,n)
    nor1 = norm(y1,n)

    do while (nor1 > eps)

       y0=y1
       nor0=nor1

       in=in+1

       x(:)=x(:)+stp(:)
       call fcn(x,y1,n)
       nor1 = norm(y1,n)

       do i=1,n
          if ( y0(i)*y1(i) < 0._pr ) then
             in2=1+in/(n)
             stp(i)=-stp(i)/real(in2)
          endif
       enddo

       if (in.eq.10000) then
          write(*,*)'Too much iteration, stop'
          stop
       endif

    enddo

  end subroutine nles2




  function norm(y,n)

    use acc;

    implicit none

    real (kind=pr) :: norm 
    !input variables
    integer, intent(in)                      :: n
    real (kind=pr), dimension(n), intent(in) :: y
    ! local variables
    integer :: i

    norm = 0._pr
    do i=1,n
       norm=norm+y(i)**2
    enddo

    norm=sqrt(norm)

  end function norm

end module NLE
