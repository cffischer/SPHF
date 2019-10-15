!====================================================================
    MODULE spline_param
!====================================================================
!
!   contains basic spline parameters
!
!--------------------------------------------------------------------

    IMPLICIT NONE
    SAVE

    INTEGER(4) :: ks     !   order of B-splines
    INTEGER(4) :: ns     !   number of splines
    INTEGER(4) :: nv     !   number of intervals ( = ns-ks+1 )
    INTEGER(4) :: ml     !   number of intervals from 0 to 1 (=1/h)
    INTEGER(4) :: me     !   number of intervals in the exponential region

    REAL(8) :: h         !   initial step in the knot sequence for z*r
    REAL(8) :: hmax      !   maximum step, t(ns+1) - t(ns) 
    REAL(8) :: rmax      !   border radius, t(ns+1)

    END MODULE spline_param
