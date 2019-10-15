!=======================================================================
      Subroutine YVAL (id,iw,mm,yv,ygr)
!=======================================================================
!
!     This routine computes the values of  r^mm f(r), f'(r), f''(r)
!     at the gaussian points of each interval, where f(r) is defined
!     by the spline expansion vector yv.
!
!     on entry
!     --------
!     id      indivates derivatives of f(r)	(=0,1,2)
!     iw      if iw > 0, results are weighted by gaussian coef.s
!     mm      integer defining the power of r
!     yv      the spline expansion vector for the funciton f(r)
!
!     on exit
!     -------
!     ygr     array of values of r^mm f(r) at the gaussian points
!             of each interval
!-----------------------------------------------------------------------

    USE spline_param
    USE spline_grid

    IMPLICIT NONE
    INTEGER, INTENT(in) :: id,iw,mm
    REAL(KIND=8), DIMENSION(ns), INTENT(in) :: yv
    REAL(KIND=8), DIMENSION(nv,ks), INTENT(out) :: ygr

    ! .. local variables

    INTEGER :: m, i, ith
    REAL(KIND=8), DIMENSION(nv,ks) :: gw

    if(mm.eq. 0) then
       gw = 1.d0
    elseif(mm.eq. 1) then
       gw = gr
    elseif(mm.eq.-1) then
       gw = grm
    elseif(mm.gt. 1) then
       gw = gr**mm
    elseif(mm.lt.-1) then
       gw = grm**(-mm)
    end if

    if(iw.ne.0) gw = gw * grw

    ygr = 0.d0

    if(id.eq.0) then

    do m = 1,ks
      do i = 1,nv
        do ith = 1,ks
         ygr(i,m) = ygr(i,m) + yv(i+ith-1)*bsp(i,m,ith)*gw(i,m)
       end do
      end do
     end do

    elseif(id.eq.1.or.id.eq.2) then

     do m = 1,ks
      do i = 1,nv
       do ith = 1,ks
        ygr(i,m) = ygr(i,m) + yv(i+ith-1)*bspd(i,m,ith,id)*gw(i,m)
       end do
      end do
     end do

    end if

    End Subroutine YVAL
