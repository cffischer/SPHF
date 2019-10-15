      
!======================================================================
      MODULE factorial
!======================================================================
!
!     Compute and store log factorial i, namely
!     GAM(I) = LOG( GAMMA(I-1) ), WHERE GAMMA(I) = FACTORIAL I-1
!
!----------------------------------------------------------------------
      REAL(KIND=8), DIMENSION(100) :: GAM 


      CONTAINS
      !================================================================
         SUBROUTINE FACTRL(NFACT) 
      !================================================================
         IMPLICIT NONE
         INTEGER , INTENT(IN) :: NFACT 
         INTEGER :: I 
         REAL(KIND=8) :: X, GAMMA

         GAMMA = 1.d0
         GAM(1) = 0.d0
         DO I = 1, NFACT - 1 
            GAMMA = I*GAMMA 
            GAM(I+1) = DLOG(GAMMA) 
         END DO 
         DO I = NFACT + 1, 100 
            X = I - 1 
            GAM(I) = GAM(I-1) + DLOG(X) 
         END DO 
         RETURN  
         END SUBROUTINE FACTRL 
   
      END MODULE factorial 
