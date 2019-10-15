!====================================================================
    SUBROUTINE vinty(ygr,yv)
!====================================================================
!
!   Computes the vector elements   <B_i, y(r)>
!
!--------------------------------------------------------------------
!
!   on entry
!   --------
!       ygr   array of values of a specific function  y(r) at the
!             gaussian points of each interval, weighted by the
!             gaussian weight
!
!   on exit
!   -------
!       yv    vector of integrals of <B_i, y(r)>, where i=1,..,ns
!
!--------------------------------------------------------------------

    USE spline_param; USE spline_grid

    IMPLICIT NONE

    REAL(8), DIMENSION(nv,ks), INTENT(IN) :: ygr
    REAL(8), DIMENSION(ns), INTENT(INOUT) :: yv

    INTEGER :: ith, i, m

    yv = 0.d0

    do ith = 1,ks
      do i = 1,nv
        do m = 1,ks
          yv(i+ith-1) = yv(i+ith-1) + ygr(i,m)*bsp(i,m,ith)
        end do
      end do
    end do

  END SUBROUTINE vinty
