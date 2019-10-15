!=======================================================================
  SUBROUTINE coulom(z,l,a, b)  
!=======================================================================
!
!   Builds full a and b matrices for the coulomb problem
!
!   a(i,j) = integration of B_i*(DD+2*z/r-l(l+1)/r**2)*B_j
!   b(i,j) = integration of B_i*B_j
!
!   from the banded symmetric form of these elementary operators
!
!-----------------------------------------------------------------------
!
!   on exit
!   -------
!       a      full matrix for coulomb operator not assumed
!              to be symmetric
!       b      full matrix of the overlap operator
!-------------------------------------------------------------------
!
    USE spline_param
    USE spline_galerkin
    USE spline_hl

    IMPLICIT NONE
    REAL(KIND=8), INTENT(IN) :: z
    INTEGER :: l
    REAL(KIND=8), DIMENSION(ns,ns), INTENT(OUT) :: a, b

    ! ..  LOCAL variables
    INTEGER :: i,j

    ! .. assemble the hl array
    CALL hlm(l)

    ! .. initialize
    a = 0.d0
    b = 0.d0
 
    ! .. lower matrix
    do j = 1,ks  
      do i = ks-j+1,ns  
        a(i,i-ks+j) = hl(i,j)  
        b(i,i-ks+j) = sb(i,j)  
      end do
    end do
 
    ! .. upper matrix by symmetry
    do j = ks+1,2*ks-1  
      do i = j+1-ks,ns  
        a(i-j+ks,i) = hl(i,2*ks-j)  
        b(i-j+ks,i) = sb(i,2*ks-j)
      end do
    end do
 
    ! .. correction for asymmetry in db2
    a(ns-1,ns) = hl(1,1)
    
  END SUBROUTINE coulom
