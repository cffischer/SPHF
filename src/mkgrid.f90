!=====================================================================
   SUBROUTINE mkgrid (z)
!=====================================================================
!  Computes the  the knots for spline calculation in the Rho =Zr
!  variable convenient for scaling.
!---------------------------------------------------------------------   
       USE spline_param 
       USE spline_grid
       IMPLICIT NONE
       REAL(kind=8), INTENT(IN) :: z
       !INTEGER, INTRINSIC:: NINT
       INTEGER::  i,nt
       REAL(KIND=8):: hp1
 
       ! .. determine ml, the number of equally spaced steps
       ml = NINT(1.d0/h)
       nt = ns + ks
       nv = ns-ks+1
       me = nv-ml
       hp1 = 1.d0 + h
       
       ! .. establish the grid for z*r     
       
       If (Allocated(t)) Deallocate(t)
       ALLOCATE (t(nt))
 
       ! the multiple knots at the origin
       t(1:ks) = 0.d0
 
       ! the equally spaced points
       DO i = ks+1, ks+ml
         t(i) = t(i-1) + h
       END DO
       
       ! the exponentially spaced points
       DO i = ks+ml+1, ns+1
        t(i) = t(i-1)*hp1
       END DO
 
       t(ns+2:nt) = t(ns+1)
 
       ! final knots for r variable
       t = t/z
       rmax = t(ns+1)
       hmax = rmax - t(ns)
END SUBROUTINE mkgrid

