! Routine get_energy calls
!  .. array
!  .. enexper  which calls
!     .. looktm, looktmdf, lookup
!     .. dev, devdf which calls
!               devf
!     .. add
!     .. lval

!======================================================================
	SUBROUTINE get_energy(done)
!======================================================================
!
!    Determine the energy data structure
!
!----------------------------------------------------------------------
     Use hf_atomic_state
     Use hf_energy_expression
     IMPLICIT NONE
     LOGICAL :: DONE
 
     ! set kmax.  This could have an upper limit.
      kmax = 2*lmax
      
! .. get contributions from the average energy
      CALL ARRAY 
      ! .. add contributions for the deviation from the averge energy
      CALL ENEXPR(done)

      
      END SUBROUTINE get_energy

!======================================================================
     SUBROUTINE array
!======================================================================
!
!    This routine gets information about the energy for the atomic 
!    state from closed shells.  QUESTION:  Do we want to divide by
!    sum(i)?  
!----------------------------------------------------------------------
     Use hf_atomic_state
     Use hf_energy_expression
     Use av_energy
     IMPLICIT NONE
     INTEGER :: IP, I, ISUMI, J, ISUMJ, K 
     REAL(KIND=8) :: DSUMI, DSUMJ, C 

      !.. define average energy arrays
      CALL init

      !print *, 'array: nclosd,nwf',nclosd, nwf
      !print *, 'qsum(1:nwf)',qsum(1:nwf)
      IP = 0 
      DO I = NCLOSD + 1, NWF 
         ISUMI =qsum(I) 
         DSUMI =qsum(I) - ISUMI 
         DO J = NCLOSD + 1, NWF 
            ISUMJ =qsum(J) 
            DSUMJ =qsum(J) - ISUMJ 
            IF (I /= J) THEN 
               C =qsum(J) 
               IF (DSUMI/=0.d0 .AND. DSUMJ/=0.d0) C = (DSUMI*(ISUMI + 1)*ISUMJ + &
                  DSUMJ*(ISUMJ + 1)*ISUMI)/qSUM(I) 
            ELSE 
               C =qsum(I) - 1.d0
               IF (DSUMI /= 0.d0) C = (ISUMI*(qSUM(I)+DSUMI-1))/qSUM(I) 
            ENDIF 
!
            IJPTR(I-NCLOSD,J-NCLOSD) = IP 
!
!  *****        Direct contribution
!
            DO K = 0, 2*MIN0(L(I),L(J)), 2 
               IP = IP + 1 
               IF (IP > 100) STOP ' COEF array too small: MAX = (100)' 
               COEF(IP) = 0.d0 
               IF (K == 0) THEN 
                  COEF(IP) = C 
               ELSE IF (I == J) THEN 
                  COEF(IP) = -C*CA(L(I),K) 
               ENDIF 
 	    !print *, 'Direct: IP, C, COEF(IP)', ip,c,coef(ip)
            END DO 
!
!  *****        Exchange contribution
!
            IF (I == J) CYCLE  
            DO K = ABS(L(I)-L(J)), L(I) + L(J), 2 
               IP = IP + 1 
               IF (IP > 100) STOP ' COEF array too small: MAX = (100)' 
               COEF(IP) = -C*CB(L(I),L(J),K) 
 	    !print *, 'XCH: IP, C, COEF(IP)', ip,c,coef(ip)
            END DO 
         END DO 
      END DO 
      RETURN  
      END SUBROUTINE ARRAY 

!======================================================================
      SUBROUTINE ENEXPR(DONE) 
!======================================================================
!
!     Determine the deviations to the average energy for the following:
!        i) an open p- , d- or f-shell
!       ii) a single electron or hole, any l
!      iii) an s-electron and a single electron, any l
!       iv) an s-electron and an open p- , d- or f-shell
!        v) an open p-shell and a single electron, any l
!       vi) a single electron in f- and a single electron in d- shells
!
!----------------------------------------------------------------------
      USE hf_atomic_state, ONLY: term, nclosd,qsum, nwf,L
      Use hf_energy_expression
      USE hf_inout
      IMPLICIT NONE
      LOGICAL, INTENT(OUT)  :: DONE 
!-----------------------------------------------
      INTEGER , EXTERNAL :: lval
      INTEGER , DIMENSION(5) :: SUMTAB 
      INTEGER , DIMENSION(11) :: PARTAB, PTRTAB 
      INTEGER , DIMENSION(54) :: LTAB 
      INTEGER , DIMENSION(2) :: NOS 
      INTEGER , DIMENSION(11) :: PLVAL 
      INTEGER :: PACVAL, SP, PS1, PS2 
      INTEGER , DIMENSION(3,54) :: FINT, GINT1, GINT2 
      INTEGER :: IP, IL, IS, J, I, IIS, IIL, NSL, IPM, NSLM, IPP, NSLP, NGGD, &
         NGGF, ISUMP, NL, LP, IPTR1, IPTR2, NOMACH, IND, LV, NP 
      REAL(KIND=8) :: C, CSP, VAL1, VAL2, VAL3 
      CHARACTER :: SL*2, SENOR, PSL*2, SLM*2, SLP*2 
      CHARACTER, DIMENSION(11) :: PARCH 
!-----------------------------------------------
!
!     ... FINT, GINT1, and GINT2 are coefficients of polynomials
!         in l, tabulated by Slater,
!
!
!
!     ... coefficients of F2 integrals for p(n)l(1) configurations
!
      DATA FINT/ 2, -1, 0, -4, -4, 3, 2, 5, 3, 2, -1, 0, -4, -4, 3, 2, 5, 3, -2&
         , 1, 0, 4, 4, -3, -2, -5, -3, -2, 1, 0, 4, 4, -3, -2, -5, -3, 4, -2, 0&
         , -2, -11, 6, -4, -4, 15, -2, 7, 15, 4, 10, 6, 0, 0, 0, 0, 0, 0, 0, 0&
         , 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0&
         , 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0&
         , 0, 0, 0, 2, -1, 0, -4, -4, 3, 2, 5, 3, 2, -1, 0, -4, -4, 3, 2, 5, 3&
         , -4, 2, 0, 2, 11, -6, 4, 4, -15, 2, -7, -15, -4, -10, -6, 0, 0, 0, -2&
         , 1, 0, 4, 4, -3, -2, -5, -3, -2, 1, 0, 4, 4, -3, -2, -5, -3/  
!
!     ... coefficients of G(l-1) integrals
!
      DATA GINT1/ -10, 5, 0, 2, 11, -6, 2, -1, -6, 14, -7, 0, 2, -13, 6, 2, -1&
         , 6, -8, 4, 0, -8, 4, 0, 4, 10, 0, 10, -5, 0, 10, -5, 0, 4, -8, 0, -8&
         , 4, 0, -2, 13, -6, 2, 11, -12, 4, 4, -12, 4, -2, -12, 0, 0, 0, -6, 3&
         , 0, 10, -5, 0, -6, 3, 0, -6, 3, 0, -4, 8, 3, 0, 12, 3, 6, 9, 0, -6, 3&
         , 0, 6, 21, 12, 12, 12, -27, 12, -6, 27, 6, -15, -12, -6, 3, 0, 0, 6, &
         -3, 0, 0, -3, 6, -3, 0, 0, -6, 3, 12, 6, 3, -4, 2, 0, -4, 2, 0, -4, 2&
         , 0, -4, 2, 0, 14, 11, -9, 14, -7, -9, -4, -2, 0, -4, 2, 0, -2, 7, 3, &
         2, 11, 3, 8, 8, 0, 0, 0, 0, -2, 1, 0, -2, 1, 0, -2, 1, 0, -2, 1, 0, -2&
         , 1, 0, 22, 13, 0/  
!
!     ... coefficients of G(l+1) integrals
!
      DATA GINT2/ 2, 5, -3, 2, -7, -15, -10, -25, -15, 2, 5, 9, 2, 17, 21, 14, &
         35, 21, 4, -2, -6, -8, -20, -12, -8, -20, -12, 4, 16, 12, 10, 25, 15, &
         10, 25, 15, 4, 10, 0, 4, 4, -12, 2, -7, -21, -2, -17, -21, -8, -20, &
         -12, 0, 0, 0, -6, -15, -9, 10, 25, 15, 6, 3, -3, 0, -12, -9, -4, -16, &
         -9, -6, -15, -9, -6, -15, -9, 6, 27, 9, 12, 30, -9, 12, 12, -27, 6, -9&
         , -27, -6, -15, -9, 0, 0, -3, 0, -6, -9, -6, -15, -9, 12, 18, 9, 0, 6&
         , 9, 6, 15, 9, -4, -10, -6, -4, -10, -6, -4, -10, -6, 14, 35, 12, 14, &
         17, -6, -4, -10, -6, 8, 8, 0, 2, -7, -6, -2, -11, -6, -4, -10, -6, -4&
         , -10, -6, 0, 0, 0, -2, -5, -3, -2, -5, -3, -2, -5, -3, 22, 31, 9, -2&
         , -5, -3, -2, -5, -3/  
!
!     ... Encoded term value -- S = LTAB/10
!                         Lterm = L + (LTAB mod 10 - 5)
!         Example:  LTAB = 36 with L = 2  is 3F
!
      DATA LTAB/ 36, 35, 34, 16, 15, 14, 46, 45, 44, 26, 25, 24, 27, 26, 25, 24&
         , 23, 25, 55, 35, 37, 36, 35, 34, 33, 17, 16, 15, 14, 13, 36, 35, 34, &
         16, 15, 14, 46, 45, 44, 26, 25, 24, 27, 26, 25, 24, 23, 25, 36, 35, 34&
         , 16, 15, 14/  
      DATA SUMTAB/ 1, 4, 7, 10, 11/  
      DATA PARTAB/ 2, 3, 1, 1, 4, 2, 2, 3, 1, 1, 2/  
      DATA PTRTAB/ 6, 12, 17, 18, 20, 30, 36, 42, 47, 48, 54/  
      DATA PLVAL/ 1, 1, 2, 0, 0, 2, 1, 1, 2, 0, 1/  
      DATA PARCH/ 'P', 'P', 'D', 'S', 'S', 'D', 'P', 'P', 'D', 'S', 'P'/  
!
      !print *, 'Entering enexpr with term =', TERM
      IP = 1 
      DO
        IF (TERM(IP:IP) == ' ') THEN 
           IP = IP + 1 
        ELSE
	   EXIT
        ENDIF 
      END DO
      SL = TERM(IP:IP+1) 
      SENOR = ' ' 
      IF (IP <= 4) SENOR = TERM(IP+2:IP+2) 
      !print *, 'SL=',SL,' Senor =', senor
!
!   ---  convert lowercase L symbol to uppercase
!
      IF (SL(2:2)>'a' .AND. SL(2:2)<'z') SL(2:2) = CHAR(ICHAR(SL(2:2))+ICHAR(&
         'A')-ICHAR('a')) 
!
!  ---  determine if FK or GK data needs to be input
!
      IL = 0 
      IS = 0 
      J = 1 
      IF (SL/='AV' .AND. SL/='aV') THEN 
         DO I = NCLOSD + 1, NWF 
            IF (qsum(I)==4*L(I) + 2 .OR. qsum(I)==0.D0) CYCLE  
            IF (J > 2) THEN 
               DONE = .FALSE. 
               RETURN  
            ENDIF 
            NOS(J) = I 
            J = J + 1 
            IF (L(I)==0 .AND. IS==0) THEN 
               IS = IS + 1 
               IIS = I 
            ELSE 
               IL = IL + 1
               IIL = I 
            ENDIF 
         END DO 
      ELSE 
         DO I = NCLOSD + 1, NWF 
            IF (qsum(I)==4*L(I) + 2 .OR. qsum(I)==0.D0) CYCLE  
            IF (J > 2) THEN 
               DONE = .TRUE. 
               RETURN  
            ENDIF 
            NOS(J) = I 
            J = J + 1 
            IF (L(I)==0 .AND. IS==0) THEN 
               IS = IS + 1 
               IIS = I 
            ELSE 
               IL = IL + 1 
               IIL = I 
            ENDIF 
         END DO 
      ENDIF 
      IF (SL/='AV' .AND. SL/='aV' .AND. IS+IL/=0) THEN 
         DONE = .FALSE. 
         C = 0.D0 
         IF (IS + IL<=2 .AND. IL<=1) THEN 
            IF (IS==0 .AND. IL==1) THEN 
    3          CONTINUE 
      !         print *,'Calling looktm', L(IIL), SL, SENOR, qsum(IIL), IP, NSL
               CALL LOOKTM (L(IIL), SL, SENOR, qsum(IIL), IP, NSL) 
      !         print *,'Returning looktm', L(IIL), SL, SENOR, qsum(IIL), IP, NSL
               IF (NSL > 1) THEN 
                  IF (L(IIL) /= 3) THEN 
                     WRITE (0, *) ' Ambiguous term: enter seniority' 
                  ELSE 
                     WRITE (0, *) ' Ambiguous term: enter Nr for f-subshells' 
                  ENDIF 
                  READ (5, '(A1)') SENOR 
                  GO TO 3 
               ENDIF 
               CALL DEV (IIL, L(IIL), qsum(IIL), IP, DONE) 
            ELSE IF (IS==1 .AND. IL==1) THEN 
               SLM = SL 
               SLP = SL 
               SLM(1:1) = CHAR(ICHAR(SLM(1:1))-1) 
               SLP(1:1) = CHAR(ICHAR(SLP(1:1))+1) 
               CALL LOOKTM (L(IIL), SLM, SENOR, qsum(IIL), IPM, NSLM) 
               CALL LOOKTM (L(IIL), SLP, SENOR, qsum(IIL), IPP, NSLP) 
               IF (NSLM + NSLP == 0) THEN 
                  DONE = .FALSE. 
                  RETURN  
               ELSE IF (NSLM==1 .AND. NSLP==0) THEN 
                  SL = SLM 
                  IP = IPM 
               ELSE IF (NSLM==0 .AND. NSLP==1) THEN 
                  SL = SLP 
                  IP = IPP 
               ELSE IF (NSLM==1 .AND. NSLP==1) THEN 
    4             CONTINUE 
                  WRITE (0, '(A,A3,A,A3)') ' Ambiguous l**n term: enter', SLM, &
                     ' or ', SLP 
                  READ (5, '(A2)') SL 
                  IF (SL == SLM) THEN 
                     IP = IPM 
                  ELSE IF (SL == SLP) THEN 
                     IP = IPP 
                  ELSE 
                     WRITE (0, *) ' Term not allowed: re-enter' 
                     GO TO 4 
                  ENDIF 
               ELSE 
    5             CONTINUE 
                  IF (L(IIL) /= 3) THEN 
                     WRITE (0, '(A,A)') ' Ambiguous l**n parent term:', &
                        'Enter term and seniority' 
                  ELSE 
                     WRITE (0, '(A,A)') ' Ambiguous l**n parent term:', &
                        'Enter term and Nr for f-subshells' 
                  ENDIF 
                  READ (5, '(A2,A1)') SL, SENOR 
                  IF (SENOR == ' ') THEN 
                     IF (L(IIL) /= 3) THEN 
                        WRITE (0, *) 'Seniority is needed' 
                     ELSE 
                        WRITE (0, *) 'Nr for f-subshells is needed' 
                     ENDIF 
                     GO TO 5 
                  ENDIF 
		  ! print *, 'First call to Looktm',L(IIL),SL,SENOR,qSUM(IIL),IP,NSL
                  CALL LOOKTM (L(IIL), SL, SENOR, qsum(IIL), IP, NSL) 
		  ! print *, 'Return from Looktm',L(IIL),SL,SENOR,qSUM(IIL),IP,NSL
                  IF (NSL /= 1) THEN 
                     WRITE (0, '(A,A3,A,A3,A)') ' Allowed terms are ', SLM, &
                        ' or ', SLP, ' plus seniority' 
                     GO TO 5 
                  ENDIF 
               ENDIF 
	       !print *, 'Call DEV', IIL,L(IIL),qSUM(IIL),IP,DONE
               CALL DEV (IIL, L(IIL), qsum(IIL), IP, DONE) 
	       !print *, 'Return DEV Before ADD', IIL,L(IIL),qSUM(IIL),IP,DONE
               IF (DONE) THEN 
                  SP = ICHAR(SL(1:1)) - ICHAR('0') 
                  CSP = (SP - 1)/2. 
		  !! This is where the average energy coefficient is changed
		  !! sp 3P.
		  !print *, 'SP, CSP, SLM', sp, csp,slm
                  IF (SL == SLM) THEN 
                     C = -CSP/(2*L(IIL)+1) 
                  ELSE 
                     C = (CSP + 1)/(2*L(IIL)+1) 
                  ENDIF 
		  !print *, 'ENEXPR: Calling add',C,L(IIL),IIS,IIL
	          !print *,'ENEXPR: Coef', coef(1:14)
                  CALL ADD (C, L(IIL), IIS, IIL, .FALSE.) 
		  !print *, 'ENEXPR: Calling add',C,L(IIL),IIL,IIS
                  CALL ADD (C, L(IIL), IIL, IIS, .FALSE.) 
               ENDIF 
	       !print *,'ENEXPR: Coef', coef(1:16)
            ELSE IF (IS==1 .AND. IL==0) THEN 
               DONE = .TRUE. 
            ENDIF 
!GG
!GG
         ELSE IF (L(NOS(1))==2 .AND. L(NOS(2))==3 .OR. L(NOS(1))==3 .AND. L(NOS&
               (2))==2) THEN 
            IF (qsum(NOS(1))==1 .AND. qsum(NOS(2))==1.D0) THEN 
   11          CONTINUE 
               CALL LOOKTMDF (SL, IP, NSL) 
               IF (NSL == 0) THEN 
                  WRITE (0, *) ' Re-enter the term of configuration' 
                  READ (5, '(A2)') SL 
		  Term = sl
                  GO TO 11 
               ENDIF 
               IF (L(NOS(1)) == 2) THEN 
                  NGGD = NOS(1) 
                  NGGF = NOS(2) 
               ELSE 
                  NGGD = NOS(2) 
                  NGGF = NOS(1) 
               ENDIF 
               CALL DEVDF (NGGD, NGGF, IP, DONE) 
            ENDIF  
!GG
!GG
         ELSE 
            IF (L(NOS(1))==1 .AND. qsum(NOS(2))==1.D0 .OR. L(NOS(2))==1 .AND. &
               qsum(NOS(1))==1.D0) THEN 
               IF (L(NOS(1))==1 .AND. qsum(NOS(2))==1.D0) THEN 
                  ISUMP = qsum(NOS(1)) 
                  NP = NOS(1) 
                  NL = NOS(2) 
               ELSE 
                  ISUMP = qsum(NOS(2)) 
                  NP = NOS(2) 
                  NL = NOS(1) 
               ENDIF 
               SP = ICHAR(SL(1:1)) - ICHAR('0') 
               LP = LVAL(SL(2:2)) 
               PS1 = SP + 1 
               PS2 = SP - 1 
               IF (ISUMP == 1) THEN 
                  IPTR1 = 1 
               ELSE 
                  IPTR1 = SUMTAB(ISUMP-1) + 1 
               ENDIF 
               IPTR2 = SUMTAB(ISUMP) 
               NOMACH = 0 
               CALL LOOKUP (PARTAB, IPTR1, IPTR2, IND, NOMACH, PS1) 
               CALL LOOKUP (PARTAB, IPTR1, IPTR2, IND, NOMACH, PS2) 
               PSL(1:1) = CHAR(PARTAB(IND)+ICHAR('0')) 
               PSL(2:2) = PARCH(IND) 
               IF (NOMACH > 1) THEN 
                  WRITE (0, *) ' AMBIGUOUS PARENT CASE' 
   10             CONTINUE 
                  WRITE (0, *) ' ENTER THE SL TERM FOR p(n) SUBSHELL' 
                  READ (5, '(A)') PSL 
                  IF (PSL(2:2)>'a' .AND. PSL(2:2)<'z') PSL(2:2) = CHAR(ICHAR(&
                     PSL(2:2))+ICHAR('A')-ICHAR('a')) 
                  PS1 = ICHAR(PSL(1:1)) - ICHAR('0') 
                  PS2 = LVAL(PSL(2:2)) 
                  CALL LOOKUP (PLVAL, IPTR1, IPTR2, IND, NOMACH, PS2) 
                  IF (NOMACH/=1 .AND. PARTAB(IND)/=PS1) GO TO 10 
               ENDIF 
               IF (ISUMP == 1) THEN 
                  IPTR1 = 1 
               ELSE 
                  IPTR1 = PTRTAB(IND-1) + 1 
               ENDIF 
               IPTR2 = PTRTAB(IND) 
               LV = L(NL) 
               PACVAL = SP*10 + LP - LV + 5 
               NOMACH = 0 
               CALL LOOKUP (LTAB, IPTR1, IPTR2, IND, NOMACH, PACVAL) 
               IF (NOMACH /= 1) THEN 
                  DONE = .FALSE. 
                  RETURN  
               ENDIF 
               VAL1 = ((FINT(1,IND)*LV+FINT(2,IND))*LV+FINT(3,IND))/(5.D0*(2*LV&
                   - 1)*(2*LV + 3)) 
               VAL2 = ((GINT1(1,IND)*LV+GINT1(2,IND))*LV+GINT1(3,IND))/(2.D0*(2&
                  *LV + 1)*(2*LV - 1)**2) 
               VAL3 = ((GINT2(1,IND)*LV+GINT2(2,IND))*LV+GINT2(3,IND))/(2.D0*(2&
                  *LV + 1)*(2*LV + 3)**2) 
!	       !print *, 'Val1,Val2,Val3', val1, val2, val3
!
!     ...  Add contributions from between p-subshell and l-electron
!
               CALL ADD (VAL1, 2, NP, NL, .TRUE.) 
               CALL ADD (VAL1, 2, NL, NP, .TRUE.) 
               CALL ADD (VAL2, LV - 1, NP, NL, .FALSE.) 
               CALL ADD (VAL2, LV - 1, NL, NP, .FALSE.) 
               CALL ADD (VAL3, LV + 1, NP, NL, .FALSE.) 
               CALL ADD (VAL3, LV + 1, NL, NP, .FALSE.) 
!
!     ... Add deviations for p-subshell
!
               CALL LOOKTM (1, PSL, ' ', qSUM(NP), IP, NSL) 
               CALL DEV (NP, 1, qSUM(NP), IP, DONE) 
            ELSE 
               DONE = .FALSE. 
            ENDIF 
         ENDIF 
      ELSE 
         DONE = .TRUE. 
      ENDIF 
      END SUBROUTINE ENEXPR 

!======================================================================
      SUBROUTINE LOOKTM(L, SL, SEN, Q, IP, NSL) 
!======================================================================
!
!     Add the deviations to the average energy for a partially filled
!       p-, d-, or f-shell
!
!----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER, INTENT(IN)  :: L
      INTEGER, INTENT(INOUT) :: IP, NSL 
      REAL(KIND=8) , INTENT(IN) :: Q 
      CHARACTER(LEN=2)  :: SL 
      CHARACTER  :: SEN 

      INTEGER :: I1, J2 
      INTEGER , DIMENSION(5) :: IPTR 
      INTEGER :: N, IBEGIN, IEND, I 
      CHARACTER(LEN=3), DIMENSION(51) :: TERMS
!-----------------------------------------------
!
      DATA IPTR/ 6, 11, 19, 35, 51/  
      DATA TERMS/ '3P2', '1D2', '1S0', '4S3', '2D3', '2P1', '3F2', '3P2', '1G2'&
         , '1D2', '1S0', '4F3', '4P3', '2H3', '2G3', '2F3', '2D1', '2D3', '2P3'&
         , '5D4', '3H4', '3G4', '3F2', '3F4', '3D4', '3P2', '3P4', '1I4', '1G2'&
         , '1G4', '1F4', '1D2', '1D4', '1S0', '1S4', '6S5', '4G5', '4F3', '4D5'&
         , '4P3', '2I5', '2H3', '2G3', '2G5', '2F3', '2F5', '2D1', '2D3', '2D5'&
         , '2P3', '2S5'/  
 
!
!  --- search for a partially unfilled p- or d-shell
!
      N = Q 
      IF (N > 2*L + 1) N = 4*L + 2 - N 
      IP = 0 
      NSL = 0 
      IF (N>1 .AND. L<=2) THEN 
         IF (L == 1) THEN 
            IBEGIN = 1 
            IEND = 6 
         ELSE 
            IBEGIN = IPTR(N-1) + 1 
            IEND = IPTR(N) 
         ENDIF 
         I = IBEGIN 
         I1 = I 
         J2 = MAX(IEND,I1) 
         DO I = I1, J2 
            IF (SL /= TERMS(I)(1:2)) CYCLE  
            IF (SEN/=' ' .AND. SEN/=TERMS(I)(3:3)) CYCLE  
            NSL = NSL + 1 
            IP = I 
         END DO 
      ELSE IF (N==1 .AND. SL(1:1)=='2') THEN 
         NSL = 1 
!GG
      ELSE IF (L == 3) THEN 
         CALL LOOKF (L, SL, SEN, N, IP, NSL) 
!GG
      ENDIF 
      
      CONTAINS

!======================================================================
      SUBROUTINE LOOKF(L, SL, SEN, N, IP, NSL) 
!======================================================================
!
!     Add the deviations to the average energy for a partially filled
!      f- shell
!
!----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER , INTENT(IN) :: L 
      INTEGER , INTENT(IN) :: N 
      INTEGER , INTENT(OUT) :: IP 
      INTEGER , INTENT(INOUT) :: NSL 
      CHARACTER(LEN=2) , INTENT(IN) :: SL 
      CHARACTER , INTENT(IN) :: SEN 

      INTEGER :: I1, J2 
      INTEGER , DIMENSION(4) :: IPTR 
      INTEGER :: IBEGIN, IEND, I 
      CHARACTER :: TER*3 
      CHARACTER , DIMENSION(144) :: TERMS*3 
      CHARACTER, DIMENSION(119) :: TERMS6*3, TERMS7*3 
!------------------------------------------------------------------------------
!
      DATA IPTR/ 7, 24, 71, 144/  
      DATA TERMS/ '3F1', '3P1', '3H1', '1D1', '1G1', '1S1', '1I1', '4S1', '4D1'&
         , '4F1', '4G1', '4I1', '2P1', '2D1', '2D2', '2F1', '2F2', '2G1', '2G2'&
         , '2H1', '2H2', '2I1', '2K1', '2L1', '5D1', '5F1', '5G1', '5S0', '5I1'&
         , '3F1', '3F2', '3D1', '3D2', '3F3', '3G1', '3G2', '3F4', '3G3', '3P1'&
         , '3P2', '3P3', '3H1', '3H2', '3H3', '3H4', '3I1', '3I2', '3K1', '3K2'&
         , '3L1', '3M1', '1D1', '1D2', '1D3', '1F1', '1G1', '1G2', '1G3', '1D4'&
         , '1G4', '1H1', '1H2', '1S1', '1I1', '1S2', '1I2', '1I3', '1K1', '1L1'&
         , '1L2', '1N1', '6P0', '6F0', '6H0', '4S1', '4P1', '4P2', '4D1', '4D2'&
         , '4D3', '4F1', '4F2', '4F3', '4F4', '4G1', '4G2', '4G3', '4G4', '4H1'&
         , '4H2', '4H3', '4I1', '4I2', '4I3', '4K1', '4K2', '4L1', '4M0', '2P1'&
         , '2P2', '2P3', '2P4', '2D1', '2D2', '2D3', '2D4', '2D5', '2F1', '2F2'&
         , '2F3', '2F4', '2F5', '2F6', '2F7', '2G1', '2G2', '2G3', '2G4', '2G5'&
         , '2G6', '2H1', '2H2', '2H3', '2H4', '2H5', '2H6', '2H7', '2I1', '2I2'&
         , '2I3', '2I4', '2I5', '2K1', '2K2', '2K3', '2K4', '2K5', '2L1', '2L2'&
         , '2L3', '2M1', '2M2', '2N1', '2O0'/  
!            ... f6 terms ...
      DATA TERMS6/ '7F0', '5D1', '5D2', '5D3', '5F1', '5F2', '5G1', '5G2', &
         '5G3', '5P0', '5H1', '5H2', '5S0', '5I1', '5I2', '5K0', '5L0', '3F1', &
         '3F2', '3F6', '3F8', '3D1', '3D2', '3D3', '3D4', '3F3', '3F5', '3G1', &
         '3G2', '3G4', '3G5', '3D5', '3F4', '3F7', '3F9', '3G3', '3G6', '3G7', &
         '3P1', '3P2', '3P3', '3H1', '3H2', '3H3', '3H4', '3P4', '3H5', '3H6', &
         '3P5', '3P6', '3H7', '3H8', '3H9', '3I1', '3I2', '3I3', '3I4', '3I5', &
         '3I6', '3K1', '3K2', '3K3', '3K4', '3K5', '3K6', '3L1', '3L2', '3L3', &
         '3M1', '3M2', '3M3', '3N0', '3O0', '1F2', '1F3', '1F4', '1D1', '1D2', &
         '1D3', '1F1', '1G1', '1G2', '1G3', '1D5', '1G5', '1D6', '1G6', '1G7', &
         '1G8', '1D4', '1G4', '1H1', '1H2', '1P0', '1H3', '1H4', '1S1', '1I1', &
         '1S2', '1I2', '1I3', '1S3', '1I4', '1I5', '1S4', '1I6', '1I7', '1K1', &
         '1K2', '1K3', '1L1', '1L2', '1L3', '1L4', '1M1', '1M2', '1N1', '1N2', &
         '1Q0'/  
!            ... f7 terms ...
      DATA TERMS7/ '8S0', '6P0', '6D0', '6F0', '6G0', '6H0', '6I0', '4S1', &
         '4S2', '4P1', '4P2', '4D1', '4D2', '4D3', '4D4', '4D5', '4D6', '4F1', &
         '4F2', '4F3', '4F4', '4F5', '4G1', '4G2', '4G3', '4G4', '4G5', '4G6', &
         '4G7', '4H1', '4H2', '4H3', '4H4', '4H5', '4I1', '4I2', '4I3', '4I4', &
         '4I5', '4K1', '4K2', '4K3', '4L1', '4L2', '4L3', '4M0', '4N0', '2S1', &
         '2S2', '2P1', '2P2', '2P3', '2P4', '2P5', '2D1', '2D2', '2D3', '2D4', &
         '2D5', '2D6', '2D7', '2F1', '2F2', '2F3', '2F4', '2F5', '2F6', '2F7', &
         '2F8', '2F9', '2FA', '2G1', '2G2', '2G3', '2G4', '2G5', '2G6', '2G7', &
         '2G8', '2G9', '2GA', '2H1', '2H2', '2H3', '2H4', '2H5', '2H6', '2H7', &
         '2H8', '2H9', '2I1', '2I2', '2I3', '2I4', '2I5', '2I6', '2I7', '2I8', &
         '2I9', '2K1', '2K2', '2K3', '2K4', '2K5', '2K6', '2K7', '2L1', '2L2', &
         '2L3', '2L4', '2L5', '2M1', '2M2', '2M3', '2M4', '2N1', '2N2', '2O0', &
         '2Q0'/  
 
!
!  --- search for a partially unfilled f- shell
!
      IF (L == 3) THEN 
         SELECT CASE (N)  
         CASE (2)  
            IBEGIN = 1 
            IEND = 7 
         CASE (6)  
            IBEGIN = 1 
            IEND = 119 
         CASE (7)  
            IBEGIN = 1 
            IEND = 119 
         CASE DEFAULT 
            IBEGIN = IPTR(N-2) + 1 
            IEND = IPTR(N-1) 
         END SELECT 
         I = IBEGIN 
         I1 = I 
         J2 = MAX(IEND,I1) 
         IF (N == 6) THEN 
            DO I = I1, J2 
               TER = TERMS6(I) 
               IF (SL /= TER(1:2)) CYCLE  
               IF (SEN/=' ' .AND. SEN/=TER(3:3)) CYCLE  
               NSL = NSL + 1 
               IP = I + 144 
               CYCLE  
            END DO 
         ELSE 
            IF (N == 7) THEN 
               DO I = I1, J2 
                  TER = TERMS7(I) 
                  IF (SL /= TER(1:2)) CYCLE  
                  IF (SEN/=' ' .AND. SEN/=TER(3:3)) CYCLE  
                  NSL = NSL + 1 
                  IP = I 
                  IP = I + 263 
               END DO 
            ELSE 
               DO I = I1, J2 
                  TER = TERMS(I) 
                  IF (SL /= TER(1:2)) CYCLE  
                  IF (SEN/=' ' .AND. SEN/=TER(3:3)) CYCLE  
                  NSL = NSL + 1 
                  IP = I 
                  CYCLE  
               END DO 
            ENDIF 
         ENDIF 
      ELSE IF (N==1 .AND. SL(1:1)=='2') THEN 
         NSL = 0 
      ENDIF 
      END SUBROUTINE LOOKF 

      END SUBROUTINE LOOKTM 

!======================================================================
      SUBROUTINE LOOKTMDF(SL, IP, NSL) 
!======================================================================
!
!     Terms for the configuration d(1) f(1)
!
!----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER , INTENT(OUT) :: IP 
      INTEGER , INTENT(OUT) :: NSL 
      CHARACTER(LEN=2) , INTENT(IN) :: SL 
      INTEGER :: I1, J2, I 
      CHARACTER(LEN=2), DIMENSION(10) :: TERMS 

!             .. d1 f1 terms
      DATA TERMS/ '1P', '1D', '1F', '1G', '1H', '3P', '3D', '3F', '3G', '3H'/  
!
      IP = 0 
      NSL = 0 
      I = 1 
      I1 = I 
      J2 = MAX(10,I1) 
      DO I = I1, J2 
         IF (SL /= TERMS(I)) CYCLE  
         NSL = NSL + 1 
         IP = I 
      END DO 
      RETURN  
      END SUBROUTINE LOOKTMDF 


!
!     -----------------------------------------------------------------
!                L O O K - U P
!     -----------------------------------------------------------------
!
      SUBROUTINE LOOKUP(TAB, P1, P2, IND, NO, KEY) 
!...Translated by Pacific-Sierra Research 77to90  4.3E  07:08:07  10/22/01  
!...Switches:                     
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: P1 
      INTEGER , INTENT(IN) :: P2 
      INTEGER , INTENT(OUT) :: IND 
      INTEGER , INTENT(INOUT) :: NO 
      INTEGER , INTENT(IN) :: KEY 
      INTEGER , INTENT(IN) :: TAB(*) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: I 
!-----------------------------------------------
      DO I = P1, P2 
         IF (TAB(I) /= KEY) CYCLE  
         NO = NO + 1 
         IND = I 
      END DO 
      RETURN  
      END SUBROUTINE LOOKUP 

!======================================================================
      SUBROUTINE DEV(IEL, L, Q, I, DONE) 
!======================================================================
!
!     Add the deviations to the average energy for a partially filled
!       p-, d-, or f- shell
!
!---------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER, INTENT(IN)  :: IEL 
      INTEGER, INTENT(IN)  :: L
      INTEGER, INTENT(INOUT) :: I
      REAL(KIND=8) , INTENT(IN) :: Q 
      LOGICAL, INTENT(OUT)  :: DONE 

      INTEGER , DIMENSION(6) :: F2PP 
      INTEGER , DIMENSION(45) :: F2DD, F4DD 
      INTEGER :: N 
!
      DATA F2PP/ -3, 3, 12, -9, 0, 6/  
      DATA F2DD/ -58, 77, 50, -13, 140, -93, 42, -12, -57, 123, 105, 69, -12, &
         -105, -69, -24, 66, 12, 39, 21, 57, -51, 30, 48, 84, 219, 111, 210, &
         138, -175, -85, 23, -22, -112, -76, -58, 167, 23, -85, 59, 140, 104, &
         86, 320, 113/  
      DATA F4DD/ 5, -70, 15, 50, 140, -30, -105, 30, 55, -45, 105, -15, 30, &
         -105, 15, -10, 45, -30, -45, 70, -55, 75, 135, 20, 0, 30, -15, 210, &
         -30, -175, -50, -40, -85, 35, 50, 110, -15, -5, 125, -25, 140, 20, -40&
         , -100, -55/  
 
      DONE = .TRUE. 
      N = Q 
      IF (N > 2*L + 1) N = 4*L + 2 - N 
      IF (N > 1) THEN 
         SELECT CASE (L)  
         CASE (1)  
            CALL ADD (2*F2PP(I)/25.D0, 2, IEL, IEL, .TRUE.) 
         CASE (2)  
            I = I - 6 
            CALL ADD (2*F2DD(I)/441.D0, 2, IEL, IEL, .TRUE.) 
            CALL ADD (2*F4DD(I)/441.D0, 4, IEL, IEL, .TRUE.) 
!GG
         CASE (3)  
            CALL DEVF (IEL, L, N, I, DONE)
!GG
         CASE DEFAULT 
            DONE = .FALSE. 
         END SELECT 
      ENDIF 

   END SUBROUTINE DEV 

!======================================================================
      SUBROUTINE DEVDF(IED, IEF, I, DONE) 
!======================================================================
!
!     Add the deviations to the average energy for a pair of
!     d and f electrons
!
!---------------------------------------------------------------------

      INTEGER  :: IED 
      INTEGER  :: IEF 
      INTEGER , INTENT(IN) :: I 
      LOGICAL , INTENT(OUT) :: DONE 
      INTEGER, DIMENSION(10) :: F1DF, F2DF, F3DF, F4DF, F5DF 

      DATA F1DF/ 5, -3, 15, -17, 33, 1, 9, -9, 23, -27/  
      DATA F2DF/ 24, 6, -11, -15, 10, 24, 6, -11, -15, 10/  
      DATA F3DF/ 30, -36, 25, 41, 16, -18, 48, -13, -29, -4/  
      DATA F4DF/ 66, -99, 66, -22, 3, 66, -99, 66, -22, 3/  
      DATA F5DF/ 1815, 990, 440, 220, 170, -1485, -660, -110, 110, 160/  

      DONE = .TRUE. 
      CALL ADD (F2DF(I)/105.D0, 2, IED, IEF, .TRUE.) 
      CALL ADD (F2DF(I)/105.D0, 2, IEF, IED, .TRUE.) 
      CALL ADD (F4DF(I)/693.D0, 4, IED, IEF, .TRUE.) 
      CALL ADD (F4DF(I)/693.D0, 4, IEF, IED, .TRUE.) 
      CALL ADD (F1DF(I)/70.D0, 1, IED, IEF, .FALSE.) 
      CALL ADD (F1DF(I)/70.D0, 1, IEF, IED, .FALSE.) 
      CALL ADD (F3DF(I)/315.D0, 3, IED, IEF, .FALSE.) 
      CALL ADD (F3DF(I)/315.D0, 3, IEF, IED, .FALSE.) 
      CALL ADD (F5DF(I)/7623.D0, 5, IED, IEF, .FALSE.) 
      CALL ADD (F5DF(I)/7623.D0, 5, IEF, IED, .FALSE.) 
      RETURN  
      END SUBROUTINE DEVDF 
   !======================================================================
         SUBROUTINE DEVF(IEL, L, N, I, DONE)
   !======================================================================
   !
   !     Add the deviations to the average energy for a partially filled
   !     f- shell
   !
   !---------------------------------------------------------------------
       IMPLICIT NONE
       INTEGER,INTENT(IN)  :: IEL
       INTEGER, INTENT(IN)  :: L
       INTEGER, INTENT(INOUT) :: I
       LOGICAL, INTENT(OUT)  :: DONE
       INTEGER :: N

       INTEGER , DIMENSION(7) :: F2FF2, F4FF2, F6FF2 
       INTEGER , DIMENSION(17) :: F2FF3, F4FF3, F6FF3 
       INTEGER , DIMENSION(47) :: F2FF4, F4FF4, F6FF4 
       INTEGER , DIMENSION(73) :: F2FF5, F4FF5, F6FF5, F2VV5, F4VV5, F6VV5 
       INTEGER , DIMENSION(119) :: F2FF6, F4FF6, F6FF6, F2VV6, F4VV6, F6VV6, &
          F2FF7, F4FF7, F6FF7, F2VV7, F4VV7, F6VV7 
       INTEGER :: J 
       REAL(KIND=8) :: VV 

       DATA F2FF2/ -2092090, 19277115, -7920055, 9175309, -9862710, 25105080, &
          11506495/  
       DATA F4FF2/ -1426425, 3871725, -2871375, -6724575, 9009325, 17117100, &
          1945125/  
       DATA F6FF2/ -1828750, -13715625, 1413125, 10058125, 2493750, 21945000, &
          1579375/  
 !             ... f3 coefficients
       DATA F2FF3/ -6936930, 16681665, -6936930, 1651650, -21966945, -4789785, &
          -6139419, 12019293, 23123100, 36005970, 7172880, 10731006, 5230225, &
          -13092079, 3798795, -11231220, 5945940/  
       DATA F4FF3/ -4729725, 1126125, -4729725, -2600325, -8456175, 150150, &
          4890600, -2372175, 15765750, 2190825, -1497600, 13953225, 2634450, &
          -2777775, 3521700, 4142775, -1535625/  
       DATA F6FF3/ -6063750, -19201875, -6063750, -10841250, 2296875, 7074375, &
          3320625, 3661875, 20212500, -8452500, 4515000, -10395000, 1500625, &
          12464375, 1194375, 5696250, 4961250/  
 !             ... f4 coefficients
       DATA F2FF4/ 74340, -105840, -40320, -105840, -220500, 123480, -3920, &
          -54828, 216612, 68880, 53280, -29088, -9016, 71064, 159516, -59976, &
          234864, 113652, 11368, 11732, 81984, -23940, 107100, 44856, -145152, &
          -138600, -89460, 362628, 23445, -75420, 156240, 110376, 24264, 307872&
          ,274995, 5544, -17388, 284004, 352800, 146412, 549360, 71883, 50589, &
          -59976, 36120, -136416, 8820/  
       DATA F4FF4/ -55440, -145530, -112770, -145530, -202860, 169785, -5390, &
          73260, -57195, -10395, -71430, 81125, -56350, 98350, 187803, -33418, &
          180250, 164871, 2254, 127925, -74550, 10080, 18270, -120505, 43750, &
          -18585, -80010, 151767, 33363, 234135, 243495, 373443, 71037, -6105, &
          202545, -45885, -84105, 83790, 485100, 181251, 67410, 72324, 239400, &
          93345, 49665, 7350, 7350/  
      DATA F6FF4/ -1732500, -831600, -1159200, -831600, -258300, 970200, -30800&
         , -188100, -492300, -231000, 57600, -297000, 186200, -667800, 790020, &
         249480, -1486800, 1019340, -107240, -605500, -16800, -31500, -459900, &
         277200, 201600, 428400, 522900, 1150380, 564795, 69300, -667800, &
         1035720, 413280, -138600, -864675, 932400, 1152900, -422100, 2772000, &
         1687140, -1159200, 264285, -332325, 579600, 462000, 1150800, 711900/  
 !             ... f5 coefficients
      DATA F2FF5/ 1, -28, -179, 14, 74, -436, 14, 256, 443, 14, -70, -134, 5452&
         , 14, -538, 874, 344, -506, 5881, -436, 14, -4, -101, -1814, 578, -38&
         , -17, 22, 131, 5939, 1789, 1738, 3064, 74, 4742, 10487, 224, 37, 196&
         , 2543, 5413, 1676, 929, 128, 1639, 8, 14251, 1241, 821, 46, 2996, -37&
         , -13426, 3209, 269, 43627, 2*14, 16, -406, 2831, 14, 6566, 422, &
         174706, -3776, 14, -14, -5672, -5, 197, -583, -88/  
      DATA F2VV5/ 195, 117, 585, 195, 1755, 8775, 195, 6825, 20475, 195, 1053, &
         1755, 26325, 195, 4095, 61425, 8775, 5265, 26325, 8775, 195, 65, 585, &
         8775, 8775, 195, 65, 65, 702, 35100, 11700, 6825, 6825, 6825, 45045, &
         32175, 585, 117, 5265, 19305, 26325, 3861, 8775, 455, 20475, 4095, &
         61425, 8775, 2925, 351, 8775, 10530, 289575, 35100, 70200, 154440, 2*&
         195, 585, 14625, 14625, 195, 96525, 8775, 740025, 67275, 195, 143, &
         32175, 39, 2925, 2925, 585/  
      DATA F4FF5/ -4, -70, -848, 7, -71, -331, 7, -269, 942, 7, -175, -94, 2339&
         , 7, 223, -226, -2705, -2237, -2458, 485, 7, -289, -79, 584, -3953, &
         -536, -31, 206, 241, 3989, 5195, 192, 355, 4, 7177, 188, 112, 335, 98&
         , 5557, 3197, 652, 283, 3242, 163, 3226, 4049, 1702, 1954, 1538, 2149&
         , 59, 19957, 1633, 97, 7057, 168, 140, 89, 553, -16, 2149, -22631, 337&
         , -107135, 24194, 7, 2177, 3334, -97, 397, -932, -40/  
      DATA F4VV5/ 39, 429, 4719, 143, 3861, 3861, 143, 3003, 11011, 143, 3861, &
         4719, 42471, 143, 99099, 2457, 42471, 42471, 42471, 14157, 143, 4719, &
         1573, 42471, 42471, 4719, 429, 1287, 3861, 84942, 28314, 1001, 3003, &
         3003, 121121, 51909, 429, 4719, 3861, 51909, 42471, 17303, 4719, 33033&
         , 819, 99099, 297297, 42471, 14157, 14157, 14157, 3861, 467181, 28314&
         , 1716, 207636, 1573, 4719, 3146, 9438, 4719, 14157, 467181, 42471, &
         3581721, 325611, 143, 51909, 51909, 9438, 9438, 14157, 1573/  
      DATA F6FF5/ -1925, -3500, -31675, 350, -5950, -1400, 350, -100, -37525, &
         350, -8750, -42700, -396200, 350, -14450, -200, -42700, -68950, &
         -239225, -74200, 350, -2800, -3325, -35000, -33250, 2800, 175, 350, &
         -5075, -68425, -107975, 3*350, -41300, -91525, 5600, 27475, 4900, &
         -435575, -113225, -679000, -11725, 2800, 175, 7000, -12425, -8575, &
         -27475, 98350, 38150, 16625, 824600, -31675, 6475, -2796325, 12250, &
         350, 6125, 4585, 1015, 9800, 532700, 24500, 687400, 86800, 1400, 42700&
         , 141400, 2450, 1400, 35525, 25900/  
      DATA F6VV5/ 5577, 16731, 184041, 5577, 50193, 50193, 5577, 1859, 184041, &
         5577, 150579, 552123, 1656369, 5577, 184041, 4563, 552123, 1656369, &
         1656369, 552123, 5577, 61347, 184041, 552123, 552123, 61347, 5577, &
         5577, 100386, 2208492, 736164, 3*5577, 2024451, 2024451, 16731, 184041&
         , 150579, 6073353, 1656369, 6073353, 552123, 20449, 1521, 184041, &
         552123, 552123, 184041, 552123, 552123, 301158, 18220059, 2208492, &
         401544, 48586824, 61347, 20449, 184041, 184041, 184041, 61347, 6073353&
         , 552123, 46562373, 4232943, 5577, 674817, 2024451, 20449, 184041, &
         184041, 184041/  
!             ... f6 coefficients
      DATA F2FF6/ -14, -227, -587, -223, -28, -16, -136, -2174, 332, -41, -61, &
         -1031, -28, 7, -331, -16, -184, 154, 476, 3184, 1546, 157, 2713, 983, &
         18869, 488, -14, 128, 3524, 142, -188, 3023, 11558, 3407, 1751, 2458, &
         2042, 614, 209, 8672, 118, 809, 664, -961, 274, 359, -187, 21721, 3797&
         , 1069, 659, 1387, 1013, -1, 59, -41, 53, 1219, -47, 3602, 976, -46048&
         , -1741, -77111, 6686, -58, -1288, 674, -97, -571, -193, -647, -158, &
         98, 1166, 1238, 401, 94541, 961, 8, 274, 13648, 17236, 1543, 974, 287&
         , 692, -28, 358066, 499, 5716, 1511, 209, 887, -23, 341, 112, 679, 4, &
         1037, 233, 98, 259, 8, 1166, 4724, 322903, 646, 484, 854, -58, 2336, &
         3382, 388, -16, -452, -35, -121, -82/  
      DATA F2VV6/ 39, 1755, 61425, 2275, 585, 117, 1755, 12285, 6825, 585, 351&
         , 8775, 585, 1755, 1755, 117, 585, 585, 5265, 19305, 19305, 975, 8775&
         , 20475, 135135, 1755, 1053, 585, 26325, 4095, 184275, 32175, 26325, &
         26325, 8775, 26325, 26325, 2925, 975, 26325, 26325, 2925, 26325, 26325&
         , 975, 5265, 5265, 289575, 26325, 2925, 8775, 8775, 19305, 39, 1755, &
         585, 3510, 9750, 4875, 26325, 26325, 289575, 26325, 740025, 67275, &
         1755, 19305, 10725, 1755, 3510, 1950, 2925, 585, 1755, 8775, 2925, 585&
         , 245700, 4095, 195, 975, 61425, 61425, 8775, 19305, 975, 2925, 15795&
         , 868725, 2340, 8775, 2925, 975, 2925, 2925, 2925, 195, 2925, 13, &
         35100, 11700, 1755, 19305, 195, 8775, 30225, 997425, 8775, 2925, 8775&
         , 1755, 8775, 90675, 3627, 117, 2925, 351, 8775, 325/  
      DATA F4FF6/ -35, -68, -1087, -688, -14, -76, -566, -6133, -1412, -107, &
         -1541, -1216, -14, -280, -1775, -2230, -466, 7, 238, -400, 15601, 64, &
         97, 71, 12916, 379, -35, 422, 179, -1297, 16064, 5755, 446, 772, 370, &
         1514, 95, -85, 359, 778, 1202, 4313, 878, 4867, 158, 43, -619, -11788&
         , 1747, -631, 2*-85, 9257, 236, 370, -89, -19, -379, 467, -293, 166, &
         16942, -86, 261440, -17321, 955, -11270, -4255, 58, -1457, -373, -397&
         , -323, 49, 1501, 5, 137, 6824, 9335, 769, 23027, 67511, 31, 487, &
         10453, 161, 227, 30695, -21961, 652, 2411, 17, 262, 461, 2*213, 56, &
         1351, -34, 5191, 899, 49, 3017, 314, 3139, 3522, -50785, 2021, -787, &
         301, 773, 2242, 2417, 311, 106, 174, 266, -193, -193/  
      DATA F4VV6/ 143, 1287, 9009, 33033, 429, 4719, 14157, 99099, 9009, 1287, &
         14157, 14157, 429, 14157, 14157, 14157, 4719, 39, 3861, 17303, 155727&
         , 429, 14157, 3003, 1090089, 4719, 3861, 14157, 1287, 99099, 297297, &
         51909, 3861, 42471, 14157, 14157, 3861, 14157, 2145, 6435, 14157, &
         23595, 19305, 42471, 4719, 3861, 42471, 467181, 42471, 14157, 2*14157&
         , 155727, 4719, 14157, 4719, 2574, 15730, 7865, 14157, 1287, 467181, &
         42471, 3581721, 325611, 14157, 155727, 51909, 14157, 28314, 9438, &
         14157, 4719, 1287, 14157, 1573, 715, 45045, 33033, 4719, 70785, 495495&
         , 9009, 14157, 155727, 1573, 14157, 127413, 1401543, 4719, 14157, 143&
         , 4719, 14157, 2*1573, 143, 7865, 1573, 70785, 4719, 1287, 155727, &
         4719, 14157, 48763, 1609179, 14157, 14157, 4719, 14157, 14157, 48763, &
         48763, 4719, 1573, 14157, 4719, 4719/  
      DATA F6FF6/ -1750, 175, -8275, -4575, -700, -39550, -14000, -88750, -850&
         , -2975, -66325, -88025, -700, -39025, -55475, -16450, -15050, 350, &
         11900, -358750, -850150, -125, 29875, -1025, -495575, -350, -1750, &
         11800, -50, -4450, -118850, -312025, -12950, -10675, -63175, 62650, &
         -10850, -16450, 1435, -8120, 144550, 40985, 17360, 123725, 2450, -8575&
         , 1225, -491225, -180775, -19775, -24325, -39725, -574175, 2625, 71575&
         , 2975, -175, -5495, -1505, 226100, 5600, -52150, -175, -193025, &
         -43750, 4900, 379750, 18550, 98875, 6475, 5425, 13825, 23450, 2450, &
         -20650, -8050, 3395, 5185, -7375, 8050, 13510, 51910, 9650, 6475, &
         150850, -5425, 1400, -555800, -1272950, 70175, -44450, 175, 10325, &
         -5425, -2975, -10675, 2800, 68635, 4200, 474005, 83825, 2450, 483875, &
         1750, -5950, -29050, 1289575, 76300, 17500, 6650, 139300, 32900, &
         255850, 340900, 23450, 10150, 169225, 78575, 350/  
      DATA F6VV6/ 5577, 2*50193, 20449, 16731, 184041, 552123, 552123, 5577, &
         16731, 2*552123, 16731, 2*552123, 2*184041, 1521, 150579, 6073353, &
         6073353, 5577, 552123, 16731, 6073353, 552123, 150579, 184041, 11583, &
         184041, 1656369, 2024451, 150579, 127413, 552123, 1656369, 150579, &
         184041, 5577, 150579, 1656369, 184041, 150579, 1656369, 61347, 150579&
         , 1656369, 18220059, 1656369, 184041, 2*552123, 6073353, 20449, 552123&
         , 184041, 7722, 122694, 20449, 1656369, 150579, 18220059, 1656369, &
         3581721, 4232943, 42471, 6073353, 674817, 552123, 84942, 122694, 2*&
         184041, 50193, 552123, 184041, 16731, 200772, 184041, 61347, 61347, &
         552123, 50193, 552123, 6073353, 61347, 184041, 4969107, 54660177, &
         736164, 552123, 16731, 61347, 3*184041, 5577, 184041, 20449, 2208492, &
         736164, 50193, 6073353, 61347, 42471, 1901757, 62757981, 552123, &
         184041, 552123, 2*552123, 2*5705271, 2*184041, 2*552123, 1573/  
!             ... f7 coefficients
      DATA F2FF7/ -98, -17, -937, -112, -4, -73, -203, 14, -34, 1, 1793, 67, &
         5323, -3517, -2629, 919, 17, 14, -98, 418, -7328, 14, 158, -82, 6038, &
         688, -109, -3322, -3, -137, -4559, -991, -391, -157, 301, -73, 109, &
         -503, 661, -638, -274, -74, 2, -362, -926, -11, -59, 28, 4144, 821, 2&
         , 241, 557, 959, 22157, 11407, 9067, 13073, 15277, 517, 779, 28, 58, &
         448, 542, 3163, 3220, 4079, 28, 1804, 272, 4304, 4282, 776, 3022, 662&
         , -274, 1556, 374, 7246, 38924, 55, 4391, 526, 73559, 1957, 13807, &
         11801, 623, 779, 29, -73, 179, 8453, 3661, -175, 44, 12488, -9881, 236&
         , 1604, 43, 102031, 7504, -1342, 272, 44, -108, 1136, 3128, -638, -9, &
         -577, -4, 272, 467, -191, -164, -872/  
      DATA F2VV7/ 195, 39, 2925, 585, 39, 585, 585, 117, 585, 117, 2925, 1755, &
         20475, 20475, 122850, 4095, 130, 117, 5265, 1755, 26325, 117, 1755, &
         4095, 20475, 2925, 2457, 20475, 325, 5265, 26325, 8775, 2925, 2925, &
         1755, 585, 585, 3510, 5850, 2925, 2925, 975, 65, 1755, 8775, 65, 195, &
         1755, 8775, 1755, 65, 650, 5850, 8775, 61425, 20475, 20475, 45045, &
         32175, 8775, 2925, 65, 195, 5265, 19305, 26325, 19305, 8775, 1755, &
         8775, 2925, 12285, 61425, 4095, 20475, 2925, 2925, 19305, 8775, 15795&
         , 868725, 351, 8775, 5265, 289575, 17550, 35100, 77220, 2925, 2925, &
         351, 585, 1170, 29250, 14625, 3861, 585, 90675, 332475, 1755, 10725, &
         975, 740025, 67275, 8775, 2925, 585, 715, 10725, 90675, 18135, 130, &
         1950, 39, 2925, 2925, 975, 585, 2925/  
      DATA F4FF7/ -49, -82, -10, -56, -3538, -538, -266, 35, 541, -23, 4219, 79&
         , -41, 1928, -545, -4546, -205, 35, -49, 124, -1891, -434, 1051, 5017&
         , -44, 241, -6046, 514, 284, -617, -5752, 557, 34, 73, 1337, -265, 281&
         , -1993, -1193, 1346, -4751, -1588, -58, -1172, -1237, -45, -1237, 14&
         , -574, 274, 172, 4999, 2689, 476, 2152, 5477, 388, 45739, 3182, 778, &
         320, 42, 137, 224, 1767, 289, -406, 267, 14, 128, 736, 13480, 1957, &
         9892, 4001, 3478, 934, -542, 398, -12290, 224416, 1922, 3001, 2620, &
         64453, 1181, 547, 2347, -84, -4, 1844, 8, 757, 7451, 742, 602, -31, &
         4655, 218216, 2663, -3485, 107, -119963, 44492, 424, 206, 111, 857, &
         8044, -30, 2076, 79, -257, -239, -149, 622, 284, -58, -304/  
      DATA F4VV7/ 143, 429, 143, 429, 14157, 4719, 1573, 429, 4719, 3861, 42471&
         , 1287, 3003, 33033, 18018, 33033, 9438, 429, 3861, 1573, 42471, 4719&
         , 14157, 99099, 27027, 42471, 99099, 9009, 14157, 42471, 42471, 14157&
         , 4719, 4719, 14157, 4719, 4719, 28314, 9438, 42471, 42471, 14157, &
         1573, 14157, 14157, 1573, 14157, 1287, 14157, 1287, 3861, 42471, 14157&
         , 4719, 9009, 33033, 3003, 363363, 51909, 14157, 4719, 143, 1573, 3861&
         , 17303, 3861, 51909, 1573, 1287, 14157, 4719, 99099, 9009, 99099, &
         297297, 42471, 14157, 155727, 4719, 127413, 1401543, 14157, 14157, &
         42471, 467181, 14157, 3146, 103818, 1573, 429, 14157, 4719, 9438, &
         47190, 23595, 155727, 4719, 146289, 1609179, 14157, 467181, 3861, &
         3581721, 325611, 4719, 14157, 1573, 17303, 51909, 48763, 48763, 3146, &
         9438, 4719, 1573, 14157, 14157, 1573, 4719/  
      DATA F6FF7/ -2450, -175, -5075, -2800, -13300, -37625, -37975, 1750, &
         -44450, -175, -20125, 7525, -1825, -5125, -12575, -20725, -7525, 1750&
         , -2450, -101500, 164500, -21700, 66850, -13250, -200, -21700, -49825&
         , -200, -10675, -20125, 239575, -16975, -22225, -28525, 41825, 10675, &
         -27475, -53725, -31325, 700, 5950, -700, -700, -700, -25550, 2975, &
         1225, 700, -51800, 2975, 1750, -27475, -16975, -53725, 3625, 13075, &
         -1825, -110975, -9275, -21875, -22225, 700, 14000, 11200, 330400, 7175&
         , 639800, -6475, 700, -26600, -25900, 92200, 9500, 400, 1900, 350, &
         23450, 30100, -33250, -197050, -7008400, 128275, 27125, 109550, &
         -972475, 47425, -291725, 2017225, 875, -875, 144725, 35875, 11725, &
         -2345, 7805, 247625, 1750, -260050, -1702925, 106750, 224350, 175, &
         6418825, -1400, 350, -2800, 58450, 108850, -4200, 102200, 17150, 19075&
         , 875, 6650, 15050, 875, 1575, 51800, 34300/  
      DATA F6VV7/ 2*5577, 2*16731, 61347, 2*184041, 16731, 184041, 5577, 61347&
         , 50193, 16731, 184041, 100386, 184041, 40898, 16731, 150579, 552123, &
         1656369, 184041, 552123, 184041, 1521, 184041, 552123, 1287, 61347, &
         1656369, 1656369, 552123, 2*184041, 552123, 2*184041, 1104246, 368082&
         , 61347, 61347, 20449, 61347, 42471, 552123, 61347, 61347, 50193, &
         552123, 50193, 16731, 2*368082, 552123, 50193, 184041, 16731, 2024451&
         , 155727, 552123, 184041, 1859, 61347, 150579, 6073353, 150579, &
         6073353, 42471, 50193, 552123, 184041, 552123, 50193, 184041, 20449, &
         61347, 184041, 6073353, 552123, 4969107, 54660177, 2*552123, 1656369, &
         18220059, 1104246, 2208492, 24293412, 184041, 16731, 552123, 184041, &
         368082, 28314, 184041, 6073353, 184041, 5705271, 20919327, 552123, &
         2024451, 1521, 46562373, 325611, 552123, 2*184041, 674817, 224939, &
         5705271, 5705271, 122694, 3146, 61347, 184041, 14157, 20449, 2*184041&
         /  
      DONE = .TRUE. 
      IF (L == 3) THEN 
         SELECT CASE (N)  
         CASE (2)  
            CALL ADD (2*F2FF2(I)/87419475.D0, 2, IEL, IEL, .TRUE.) 
            CALL ADD (2*F4FF2(I)/87419475.D0, 4, IEL, IEL, .TRUE.) 
            CALL ADD (2*F6FF2(I)/87419475.D0, 6, IEL, IEL, .TRUE.) 
         CASE (3)  
            J = I - 7 
            CALL ADD (2*F2FF3(J)/96621525.D0, 2, IEL, IEL, .TRUE.) 
            CALL ADD (2*F4FF3(J)/96621525.D0, 4, IEL, IEL, .TRUE.) 
            CALL ADD (2*F6FF3(J)/96621525.D0, 6, IEL, IEL, .TRUE.) 
         CASE (4)  
            J = I - 24 
            CALL ADD (2*F2FF4(J)/737100.D0, 2, IEL, IEL, .TRUE.) 
            CALL ADD (2*F4FF4(J)/1486485.D0, 4, IEL, IEL, .TRUE.) 
            CALL ADD (2*F6FF4(J)/6625476.D0, 6, IEL, IEL, .TRUE.) 
         CASE (5)  
            J = I - 71 
            VV = F2VV5(J)*1.D0 
            CALL ADD (2*F2FF5(J)/VV, 2, IEL, IEL, .TRUE.) 
            VV = F4VV5(J)*1.D0 
            CALL ADD (2*F4FF5(J)/VV, 4, IEL, IEL, .TRUE.) 
            VV = F6VV5(J)*1.D0 
            CALL ADD (2*F6FF5(J)/VV, 6, IEL, IEL, .TRUE.) 
         CASE (6)  
            J = I - 144 
            VV = F2VV6(J)*1.D0 
            CALL ADD (2*F2FF6(J)/VV, 2, IEL, IEL, .TRUE.) 
            VV = F4VV6(J)*1.D0 
            CALL ADD (2*F4FF6(J)/VV, 4, IEL, IEL, .TRUE.) 
            VV = F6VV6(J)*1.D0 
            CALL ADD (2*F6FF6(J)/VV, 6, IEL, IEL, .TRUE.) 
         CASE (7)  
            J = I - 263 
            VV = F2VV7(J)*1.D0 
            CALL ADD (2*F2FF7(J)/VV, 2, IEL, IEL, .TRUE.) 
            VV = F4VV7(J)*1.D0 
            CALL ADD (2*F4FF7(J)/VV, 4, IEL, IEL, .TRUE.) 
            VV = F6VV7(J)*1.D0 
            CALL ADD (2*F6FF7(J)/VV, 6, IEL, IEL, .TRUE.) 
         CASE DEFAULT 
            DONE = .FALSE. 
         END SELECT 
      ELSE 
         DONE = .FALSE. 
      ENDIF 
      RETURN  
      END SUBROUTINE DEVF 

!======================================================================
      SUBROUTINE ADD(C, K, I, J, FIRST) 
!======================================================================
!
!     Add a Slater integral to the data structure associated with the
!     energy expression
!
!----------------------------------------------------------------------
      USE hf_atomic_state
      USE hf_energy_expression
      IMPLICIT NONE
      INTEGER , INTENT(IN) :: K 
      INTEGER , INTENT(IN) :: I 
      INTEGER , INTENT(IN) :: J 
      REAL(KIND=8) , INTENT(IN) :: C 
      LOGICAL , INTENT(IN) :: FIRST 
      INTEGER :: IP 
!-----------------------------------------------
!
      IP = IJPTR(I-NCLOSD,J-NCLOSD) 
 
      IF (FIRST) THEN 
         COEF(IP+K/2+1) = C/qSUM(I) + COEF(IP+K/2+1) 
      ELSE 
         IP = IP + MIN(L(I),L(J)) + 1 + (K - ABS(L(I)-L(J)))/2 + 1 
	 Print *,'L(i),L(j),IP', L(i),L(j), IP
         COEF(IP) = COEF(IP) + C/qSUM(I) 
      ENDIF 
      RETURN  
      END SUBROUTINE ADD 


!======================================================================
      INTEGER FUNCTION LVAL (SYMBOL) 
!======================================================================
!
!    Look up the l-values associated with the symbol
!----------------------------------------------------------------------
      IMPLICIT NONE
      CHARACTER , INTENT(IN) :: SYMBOL 
      INTEGER :: LOCATE 
      CHARACTER(LEN=40) :: & 
                SET='spdfghiklmnopqrstuvwxSPDFGHIKLMNOPQRSTUVWX' 
 
      LOCATE = INDEX(SET,SYMBOL) 
      IF (LOCATE <= 20) THEN 
         LVAL = LOCATE - 1 
      ELSE 
         LVAL = LOCATE - 20 
      ENDIF 
      RETURN  
      END FUNCTION LVAL 


