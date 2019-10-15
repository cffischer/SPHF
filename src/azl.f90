!=======================================================================
   Real(8) FUNCTION azl(z,h,ks,lp)
!=======================================================================
!
!   Value of B_(lp+1)/r^lp at r = 0 where lp = l+1
!----------------------------------------------------------------------

     IMPLICIT NONE

     INTEGER, INTENT(in) :: ks,lp
     REAL(8), INTENT(in) :: z,h

     INTEGER(4) :: j
     REAL(8) :: c

     IF (lp < ks ) THEN
       azl = 1.d0
       c = z/h
       do j = 1,lp
         azl = (azl*c*(ks-j))/(j*j)
       end do
     ELSE
       azl = 0.d0
     END IF

   END FUNCTION azl
