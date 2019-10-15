!=======================================================================
  SUBROUTINE reset_int
!=======================================================================
!   Reset the integral parameters
!----------------------------------------------------------------------
!
    USE spline_hl
    USE spline_integrals
    USE spline_moments
   
    lh    = -1   
    itype = 'aaa'
    krk   = -100
    mtype = 'aaa'
    kmk   = -100

  END SUBROUTINE reset_int
