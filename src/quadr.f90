!======================================================================
      Function quadr(i,j,k)
!======================================================================
!
!       <P_i| r^k |P_j> with inclusion of rel.shift if rel = .true.
!
!----------------------------------------------------------------------

      USE spline_param
      USE hf_atomic_state
      USE hf_orbitals 
      USE spline_galerkin

      IMPLICIT NONE

      INTEGER, INTENT(in) :: i,j,k
      REAL(KIND=8) :: quadr
      REAL(KIND=8), EXTERNAL :: BVMV

      REAL(KIND=8) , DIMENSION(ns,ks) :: array

      ! get the banded array

      if (k == 0) then
         array =sb
      else if (k == 1) then
         array = r1
      else if (k == -1) then
         array = rm1
      else if (k == -2) then
         array = rm2
      else
         Call mrm(k,array)
      end if
     
      ! type of symmetry is symmetric
      quadr = BVMV(ns,ks,array,'s',p(1,i),p(1,j))

    END FUNCTION quadr

