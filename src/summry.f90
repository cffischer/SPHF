!     ------------------------------------------------------------------
      SUBROUTINE SUMMRY(orb_diff, scf_diff, tail_c)
!     ------------------------------------------------------------------
!
!       The results of a calculation are summarized.   These include
!   the following for each electron:
!
!          E(NL)   - diagonal energy parameter
!          I(NL)   - -(1/2)<nl|L|nl>
!          KE      - I(NL) + Z <r>
!          !REL     - Relativistic shift (mass-velocity, Darwin term,
!          !          spin-spin contact term)
!          SIGMA   - screening parameter as defined by Eq. (6-  ).
!          AZ(NL)  - starting parameter, P(r)/r**(l+1) as r -> 0.
!          1/R**3  - expected value of <1/r**3>
!          1/R     - expected value of <1/r>
!          R       - expected mean radius
!          R**2    - expected value of <r**2>
!
!   These results are followed by:
!
!          KINETIC ENERGY (EK)
!          POTENTIAL ENERGY (EP) = ETotal - EN
!          RATIO                 = EP/EN
!          TOTAL ENERGY (ETotal)
!
!-----------------------------------------------

      USE spline_param
      USE hf_atomic_state
      USE hf_inout
      USE hf_orbitals 
      USE spline_integrals

      IMPLICIT NONE
      REAL(KIND=8), INTENT(IN) :: orb_diff, scf_diff, tail_c

      INTEGER :: i,lp,ip
      REAL(KIND=8), EXTERNAL :: BVMV, azl, bhl, quadr

      REAL(KIND=8), DIMENSION(nwf) :: r1,hli,rm1
      REAL(KIND=8) :: PI, ekinp, en, rm3, rp2, rz, epotl,ratio, rh, Enl
!
      WRITE(UNIT=log,FMT='(/7X,A)') 'Convergence (latest difference)'
      WRITE(UNIT=log,FMT= '(12X,A,1PD11.3)') '    SCF diff. =', scf_diff
      WRITE(UNIT=log,FMT= '(12X,A,1PD11.3)') 'Orbital diff. =', orb_diff
      WRITE(UNIT=log,FMT= '(12X,A,1PD11.3)') 'Tail cut-off  =', tail_c
      PI = ACOS((-1.d0))
      WRITE (UNIT=log,FMT= 9) ATOM, TERM 
    9 FORMAT(/,/,/,24X,'ATOM ',A6,3X,'TERM ',A6,/,/,2X,'nl',8X,'E(nl)',9X,&
         'I(nl)',9X,'KE(nl)',7X,'S(nl)',6X,'Az(nl)',6X,'MAXR',4X,'MAXP')

      EN = 0.d0
!     end_tol = nwf*1.d-09
!
!  *****  COMPUTE AND PRINT ONE-ELECTRON PARAMETERS
!
      DO ip = 1, NWF 
         i = iord(ip)
         lp = l(i)+1
         az(i) = p(lp+1,i)*azl(z,h,ks,lp)
         R1(I) = QUADR(I,I,1) 
         HLI(I) = -0.5*bhl(i,i)
         RM1(I) = QUADR(I,I,-1) 
         EKINP = HLI(I) + Z*RM1(I) 
         !Print '(10X,A,I5,3F16.9)', 'i,hli,z*rm1(i),ekinp', i,hli(i),z*rm1(i),ekinp
         !EKINP = -0.5*bvmv(ns,ks,db2,'s',p(:,i),p(:,i))
         EN = EN + qSUM(I)*EKINP 
         RH = 3*N(I)*N(I) - L(I)*(L(I)+1) 
         S(I) = Z - 0.5*RH/R1(I) 
         ! Find the end of the function
         !m=ns/2
         !do j = ns, m, -1
         !  if (abs(p(j,i)) > end_tol ) exit
         !end do
         !maxr(i)=min(ns,j+1)
         
         Enl = e(i,i)
         WRITE (UNIT=log,FMT= 15) EL(I),Enl,HLI(I),EKINP,S(I),AZ(i),maxr(i), maxloc(abs(p(:,i)))
   15    FORMAT(1X,A3,F15.8,2F15.8,F9.3,F14.6,2I8) 
      END DO 
!
!  *****  Compute Moments
!
      WRITE (UNIT=log,FMT= 8) 'Delta(r)'
    8 FORMAT(//2X,'nl',7X,A8,8X,'1/R**3',8X,'1/R',10X,'R',10X,'R**2') 
      DO ip = 1, NWF 
         i = iord(ip)
         RM3 = 0 
         IF (L(I) /= 0) RM3 = QUADR(I,I,-3) 
         RP2 = QUADR(I,I,2) 
         RZ = 0.D0 
         IF (L(I) == 0) RZ = AZ(I)**2/(4.*PI)
         WRITE (UNIT=log,FMT= 16) EL(I), RZ, RM3, RM1(I), R1(I), RP2 
   16    FORMAT(1X,A3,F15.6,F15.6,F13.6,F12.8,F12.6) 
      END DO 
      
      EPOTL = ETotal - EN 
      RATIO = EPOTL/EN 
      WRITE (UNIT=err,FMT= 26) ETotal, EN,  EPOTL,  RATIO 
      WRITE (UNIT=log,FMT= 26) ETotal, EN,  EPOTL,  RATIO 
   26 FORMAT(//7X,'TOTAL ENERGY (a.u.)',5X,F25.15 / &
          T22,'Kinetic   ',F25.15,/,10X, &
          T22,'Potential ',F25.15,/,10X, &
          T22,'Ratio     ',F25.15) 
      RETURN  
      END SUBROUTINE SUMMRY 
!======================================================================
