!====================================================================
    MODULE spline_hl
!====================================================================
!
!   contains the spline representation of L operator
!
!--------------------------------------------------------------------
   
    Use spline_param

    IMPLICIT NONE
    SAVE

    INTEGER :: lh = -1    ! l-value for current L-operator
    
    REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: hl, vc

!   hl(1:ns,1:ks) - matrix of L-operator in the B-spline basis
!                   (in almost symmetric lower-column mode)
!   vc(1:ns,1:ks) - matrix of mass-velocity correction

   CONTAINS

!====================================================================
      SUBROUTINE allocate_hl
!====================================================================

      USE spline_param

      if (allocated(hl)) deallocate(hl, vc)
      ALLOCATE( hl(ns,ks), vc(ns,ks))
	  
	  hl = 0.d0; vc = 0.d0; lh = -2

      END SUBROUTINE allocate_hl

    END MODULE spline_hl






