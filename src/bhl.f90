!======================================================================
   Function BHL(i,j)
!======================================================================
!
!       <P_i| L |P_j> with inclusion of rel.shift if rel = .true.
!
!----------------------------------------------------------------------

     USE spline_param
     USE hf_atomic_state
     USE hf_orbitals 
     USE spline_hl

     IMPLICIT NONE

     INTEGER, INTENT(in) :: i,j
     REAL(KIND=8) :: BHL
     REAL(KIND=8), EXTERNAL :: BVMV

     if (iabs(L(i)-L(j)) .NE. 0) Stop ' HL:  LI <> LJ'

     Call HLM (l(i))

     BHL = BVMV(ns,ks,hl,'s',p(1,i),p(1,j))

    END FUNCTION BHL

