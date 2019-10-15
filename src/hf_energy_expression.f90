!=======================================================================
  MODULE hf_energy_expression
!=======================================================================
!   This module defines the data structure for the energy expression
!----------------------------------------------------------------------
!
    IMPLICIT NONE
    SAVE

    ! energy expression variables
    INTEGER :: kmax
    INTEGER, DIMENSION(5,5) :: ijptr
    REAL(KIND=8), DIMENSION(100) :: coef

 END module hf_energy_expression
