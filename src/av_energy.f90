!    Determines average energy.  Calls
!    ..  rme
!======================================================================
   MODULE av_energy
!======================================================================

     REAL(KIND=8), DIMENSION(10) :: CCA 
     REAL(KIND=8), DIMENSION(35) :: CCB 
   
     CONTAINS
   
   !======================================================================
     SUBROUTINE init
   !======================================================================
   !
   !    This routine initializes arrays for the average energy 
        USE factorial
        IMPLICIT NONE
         
   !
         CALL FACTRL (32)
   !
   !  *****  AVERAGE INTERACTIONS FOR EQUIVALENT ELECTRONS
   !
   !  *****  P - P
   !
         CCA(1) = 2.D0/25.D0 
   !
   !  *****  D - D
   !
         CCA(2) = 2.D0/63.D0 
         CCA(3) = 2.D0/63.D0 
   !
   !  *****  F - F
   !
         CCA(4) = 4.D0/195.D0 
         CCA(5) = 2.D0/143.D0 
         CCA(6) = 100.D0/5577.D0 
   !
   !  *****  G - G
   !
         CCA(7) = 20.D0/1309.D0 
         CCA(8) = 162.D0/17017.D0 
         CCA(9) = 20.D0/2431.D0 
         CCA(10) = 4410.D0/371943.D0 
   !
   !
   !  ***** AVERAGE INTERACTIONS FOR NON-EQUIVALENT ELECTRONS
   !
   !  *****  S - ( S, P, D, F, G )
   !
         CCB(1) = 1.D0/2.D0 
         CCB(2) = 1.D0/6.D0 
         CCB(3) = 1.D0/10.D0 
         CCB(4) = 1.D0/14.D0 
         CCB(5) = 1.D0/18.D0 
   !
   !  *****  P - ( P, D, F, G )
   !
         CCB(6) = 1.D0/6.D0 
         CCB(7) = 1.D0/15.D0 
         CCB(8) = 1.D0/15.D0 
         CCB(9) = 3.D0/70.D0 
         CCB(10) = 3.D0/70.D0 
         CCB(11) = 2.D0/63.D0 
         CCB(12) = 2.D0/63.D0 
         CCB(13) = 5.D0/198.D0 
   !
   !  *****  D - ( D, F, G )
   !
         CCB(14) = 1.D0/10.D0 
         CCB(15) = 1.D0/35.D0 
         CCB(16) = 1.D0/35.D0 
         CCB(17) = 3.D0/70.D0 
         CCB(18) = 2.D0/105.D0 
         CCB(19) = 5.D0/231.D0 
         CCB(20) = 1.D0/35.D0 
         CCB(21) = 10.D0/693.D0 
         CCB(22) = 5.D0/286.D0 
   !
   !  *****  F - ( F, G )
   !
         CCB(23) = 1.D0/14.D0 
         CCB(24) = 2.D0/105.D0 
         CCB(25) = 1.D0/77.D0 
         CCB(26) = 50.D0/3003.D0 
         CCB(27) = 2.D0/63.D0 
         CCB(28) = 1.D0/77.D0 
         CCB(29) = 10.D0/1001.D0 
         CCB(30) = 35.D0/2574.D0 
   !
   !  *****  G - ( G )
   !
         CCB(31) = 1.D0/18.D0 
         CCB(32) = 10.D0/693.D0 
         CCB(33) = 9.D0/1001.D0 
         CCB(34) = 10.D0/1287.D0 
         CCB(35) = 245.D0/21879.D0 
         RETURN  
     END SUBROUTINE INIT 
   
   !======================================================================
     REAL(KIND=8) FUNCTION CA (l,k)
   !======================================================================
   !
   !    Computes the average-energy direct contribution from a pair of 
   !    l orbitals    
   !----------------------------------------------------------------------
       IMPLICIT NONE
       INTEGER, EXTERNAL :: rme
       INTEGER, INTENT(IN) :: l,k
   !-----------------------------------------------
   !
       IF (L <= 4) THEN 
          CA = CCA((L*(L-1)+K)/2) 
       ELSE 
          ! corrected according to Prof. P. Bogdanovich 1996.03.18
          CA = RME(L,L,K)**2/((2*L + 1)*(4*L + 1)) 
       ENDIF 
      
     END FUNCTION CA 
   
   !======================================================================
     REAL(KIND=8) FUNCTION CB (l,lp,k)
   !======================================================================
   !
   !    Computes the average-energy exchange contribution from a pair of 
   !    l and lp orbitals    
   !----------------------------------------------------------------------
       IMPLICIT NONE
       INTEGER, EXTERNAL :: rme
       INTEGER, INTENT(IN) :: l,lp,k
       INTEGER , DIMENSION(0:4) :: ICBPTR 
       INTEGER :: L1, L2 
   !-----------------------------------------------
       DATA ICBPTR/ 1, 6, 14, 23, 31/  
   !
       IF (L <= LP) THEN 
         L1 = L 
         L2 = LP 
       ELSE 
         L1 = LP 
         L2 = L 
       ENDIF 

       IF (L2 <= 4) THEN 
        CB = CCB(ICBPTR(L1)+(K+L1-L2)/2+(L1+1)*(L2-L1)) 
       ELSE 
        CB = RME(L,LP,K)**2/(2*(2*L + 1)*(2*LP + 1)) 
       ENDIF 

     END FUNCTION CB 
   
   END MODULE av_energy
!======================================================================
   FUNCTION RME (L, LP, K) 
!======================================================================
!
!--- EVALUATES THE REDUCED MATRIX ELEMENT (L//C(K)//LP)  -  SEE FANO
!    AND RACAH, IRREDUCIBLE TENSORIAL SETS, CHAP. 14, P. 81
!
!
!------------------------------------------------------------------
     USE factorial
     IMPLICIT NONE
     REAL(KIND=8) :: RME
     INTEGER , INTENT(IN) :: L 
     INTEGER , INTENT(IN) :: LP 
     INTEGER , INTENT(IN) :: K 
     INTEGER :: I2G, IG, I1, I2, I3 
     REAL(KIND=8) :: QUSQRT 

     IF (MIN0(L,LP) == 0) THEN 
       RME = 1.D0 
     ELSE IF (K == 0) THEN 
       RME = 2*L + 1 
       RME = DSQRT(RME) 
     ELSE 
       I2G = L + LP + K 
       IG = I2G/2 
       IF (I2G - 2*IG /= 0) THEN 
          RME = 0.D0 
       ELSE 
         I1 = IG - L 
         I2 = IG - LP 
         I3 = IG - K 
         QUSQRT = (2*L + 1)*(2*LP + 1) 
         RME = DSQRT(QUSQRT)*DEXP((GAM(2*I1+1)+GAM(2*I2+1)+GAM(2*I3+1)-GAM(&
            I2G+2))/2.D0+GAM(IG+1)-GAM(I1+1)-GAM(I2+1)-GAM(I3+1)) 
       ENDIF 
     ENDIF 

   END FUNCTION RME 

