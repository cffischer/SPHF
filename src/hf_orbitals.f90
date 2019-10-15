!=======================================================================
  MODULE hf_orbitals
!=======================================================================
!   This module defines the orbital parameters 
!----------------------------------------------------------------------
!
    IMPLICIT NONE
    SAVE

    ! orbital variables dimensions depend on nwf  and ns
    REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: p

    CONTAINS
    !===================================================================
      SUBROUTINE allocate_orbital_arrays
    !===================================================================
    !   This program allocates arrays associated with the orbitals 
    !   nwf and ns
    !-------------------------------------------------------------------
        USE spline_param, ONLY: ns
        USE hf_atomic_state, ONLY: nwf
        IMPLICIT NONE
        
	if (allocated(p)) Deallocate(p)
  	ALLOCATE( p(ns,nwf) )
       
      END SUBROUTINE allocate_orbital_arrays
 END module hf_orbitals
