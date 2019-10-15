!=======================================================================
   FUNCTION bvalu2 (t, bcoef, ns, ks, x, jderiv)
!=======================================================================
!
!   This routine is a modification of de Boor's bvalue routine, modified
!   to return the value at right endpoint continuous from the left, 
!   instead of 0. It assumes the usual knot multiplicity ks at the right
!   endpoint.
!
!   It calculates the value at  x of the jderiv-th derivative of spline 
!   from b-repr. The spline is taken to be continuous from the right.
!
!   SUBROUTINES contained:
!       interv
!
!-----------------------------------------------------------------------
!
!   on entry
!   --------
!
!       t      knot sequence, of length  ns+ks, assumed nondecreasing.
!       bcoef  coefficient sequence, of length  ns .
!       ns     length of  bcoef, assumed positive.
!       ks     order of the spline .
!
!             . . . W A R N I N G . . .   
!       The restriction  ks <= kmax (=15)  is imposed
!       arbitrarily by the parameter statement defining dimensions
!       for aj, dm, dp  below, but is  NEVER CHECKED.
!
!       x      the point at which to evaluate the spline.
!       jderiv integer giving the order of the derivative to be evaluated
!              ASSUMED to be zero or positive.
!
!   on exit
!   -------
!   bvalu2 - the value of the (jderiv)-th derivative of  f  at  x .
!
!   method
!   ------
!     the nontrivial knot interval  (t(i),t(i+1))  containing  x  is lo-
!   cated with the aid of  INTERV. the  ks  b-coeffs of  f  relevant for
!   this interval are then obtained from  bcoef (or taken to be zero if
!   not explicitly available) and are then differenced  jderiv  times to
!   obtain the b-coeffs of  (d**jderiv)f  relevant for that interval.
!   precisely, with  j = jderiv, we have from x.(12) of the text that
!
!           (d**j)f  =  sum ( bcoef(.,j)*b(.,ks-j,t) )
!
!   where
!                  / bcoef(.),                     ,  j = 0
!                  /
!   bcoef(.,j)  =  / bcoef(.,j-1) - bcoef(.-1,j-1)
!                  / ----------------------------- ,  j > 0
!                  /    (t(.+ks-j) - t(.))/(ks-j)
!
!   then, we use repeatedly the fact that
!
!     sum ( a(.)*b(.,m,t)(x) )  =  sum ( a(.,x)*b(.,m-1,t)(x) )
!   with
!                  (x - t(.))*a(.) + (t(.+m-1) - x)*a(.-1)
!     a(.,x)  =    ---------------------------------------
!                  (x - t(.))      + (t(.+m-1) - x)
!
!   to write  (d**j)f(x)  eventually as a linear combination of b-splines
!   of order  1 , and the coefficient for  b(i,1,t)(x)  must then
!   be the desired number  (d**j)f(x). (see x.(17)-(19) of text).
!
!   current fortran standard makes it impossible to specify the length
!   of t  precisely without the introduction of otherwise superfluous
!   additional arguments.
!
!-----------------------------------------------------------------------

     IMPLICIT NONE
     REAL(KIND=8) :: bvalu2
     INTEGER, INTENT(IN) :: ns, ks, jderiv
     REAL(KIND=8), INTENT(IN) :: x
     REAL(KIND=8), DIMENSION(ns), INTENT(IN) :: bcoef
     REAL(KIND=8), DIMENSION(ns+ks), INTENT(IN) :: t

     ! ..  local variable

     REAL(KIND=8), DIMENSION(ks) :: aj, dm, dp
     REAL(KIND=8) :: fkmj
     INTEGER :: i, mflag, km1, jcmin, imk, nmi, jcmax, jc
     INTEGER :: j, jj, kmj, ilo
  
     bvalu2 = 0.d0
     IF (jderiv >= ks) RETURN

     ! .. find  i  such that 1 <= i < ns+ks  and  t(i) < t(i+1) and
     ! t(i) <= x < t(i+1) . if no such i can be found,  x  lies
     ! outside the support of  the spline  f  and  bvalu2 = 0.
     ! (the asymmetry in this choice of i makes f rightcontinuous)

     IF( x /= t(ns+1) .OR. t(ns+1) /= t(ns+ks) ) THEN
       CALL interv(x, i, mflag)
       IF (mflag /= 0) RETURN
     ELSE
       i = ns
     END IF

     ! .. if ks = 1 (and jderiv = 0), bvalu2 = bcoef(i).

     km1 = ks - 1
     IF ( km1 <= 0 ) THEN
       bvalu2 = bcoef(i)
       RETURN
     END IF
 
     ! .. store the ks b-spline coefficients relevant for the knot interval
     ! (t(i),t(i+1)) in aj(1),...,aj(ks); compute dm(j) = x - t(i+1-j),
     ! dp(j) = t(i+j) - x, j=1,...,ks-1. Set any of the aj not obtainable
     ! from input to zero. set any t.s not obtainable equal to t(1) or
     ! to t(ns+ks) appropriately.
 
     jcmin = 1
     imk = i - ks
     IF (imk < 0) THEN   
       jcmin = 1 - imk
       do j=1,i
         dm(j) = x - t(i+1-j)
       end do
       do j=i,km1
         aj(ks-j) = 0.
         dm(j) = dm(i)
       end do
     ELSE
       do j=1,km1
         dm(j) = x - t(i+1-j)
       end do
     END IF
 
     jcmax = ks
     nmi = ns - i
     IF (nmi < 0) THEN  
       jcmax = ks + nmi
       do j=1,jcmax
         dp(j) = t(i+j) - x
       end do
       do j=jcmax,km1
         aj(j+1) = 0.
         dp(j) = dp(jcmax)
       end do
     ELSE
       do j=1,km1
         dp(j) = t(i+j) - x
       end do
     END IF    
       do jc=jcmin,jcmax
         aj(jc) = bcoef(imk + jc)
       end do
 
     ! .. difference the coefficients  jderiv  times.
 
     IF (jderiv /= 0) THEN
       DO j=1,jderiv
         kmj = ks-j
         fkmj = kmj
         ilo = kmj
         DO jj=1,kmj
           aj(jj) = ((aj(jj+1) - aj(jj))/(dm(ilo) + dp(jj)))*fkmj
           ilo = ilo - 1
         END DO
       END DO
     END IF
 
     ! .. compute value at  x  in (t(i),t(i+1)) of jderiv-th derivative,
     ! given its relevant b-spline coeffs in aj(1),...,aj(ks-jderiv).
 
     IF (jderiv /= km1) THEN
       DO j=jderiv+1,km1
         kmj = ks-j
         ilo = kmj
         DO jj=1,kmj
           aj(jj) = (aj(jj+1)*dm(ilo) + aj(jj)*dp(jj))/(dm(ilo)+dp(jj))
           ilo = ilo - 1
         END DO
       END DO
     END IF
     bvalu2 = aj(1)
 
     CONTAINS
     
     !=======================================================================
       SUBROUTINE interv (x, left, mflag)
     !=======================================================================
     !
     !   Computes  left = max( i ; 1 <= i <= ns+ks  .and.  t(i) <= x )
     !   which is the interval containing x .
     !
     !   A reformatted version of the de Boor routine
     !
     !-----------------------------------------------------------------------
     !
     !   on entry
     !   --------
     !
     !       x   the point whose location with respect to the sequence t is
     !           to be determined.
     !
     !   on exit
     !   -------
     !       left, mflag   integers, whose values are
     !          1     -1   if               x <  t(1)
     !          i      0   if      t(i)  <= x < t(i+1)
     !      ns+ks      1   if   t(ns+ks) <= x
     !
     !   In particular,  mflag = 0 is the 'usual' case.  mflag /= 0
     !   indicates that  x  lies outside the halfopen interval
     !   t(1) <= y < t(ns+ks). the asymmetric treatment of the
     !   interval is due to the decision to make all pp functions cont-
     !   inuous from the right.
     !
     !   ...  m e t h o d  ...
     !   The program is designed to be efficient in the common situation where
     !   it is called repeatedly, with  x  taken from an increasing or decrea-
     !   sing sequence. this will happen, e.g., when a pp function is to be
     !   graphed. the first guess for  left  is therefore taken to be the val-
     !   ue returned at the previous call and stored in the  l o c a l  varia-
     !   ble  ilo . a first check ascertains that  ilo < ns+ks (this is nec-
     !   essary since the present call may have nothing to do with the previ-
     !   ous call). then, if  t(ilo) <= x < t(ilo+1), we set  left = ilo
     !   and are done after just three comparisons.
     !   otherwise, we repeatedly double the difference  istep = ihi - ilo
     !   while also moving  ilo  and  ihi  in the direction of  x , until
     !                      t(ilo) <= x < t(ihi) ,
     !   after which we use bisection to get, in addition, ilo+1 = ihi .
     !   left = ilo  is then returned.
     !
     !------------------------------------------------------------------------
     !
       IMPLICIT NONE
       REAL(KIND=8), INTENT(IN) :: x
       INTEGER, INTENT(OUT) :: left, mflag

       ! .. local variables
 
       INTEGER:: lxt, ihi, istep, middle
       INTEGER:: ilo=1
 
       lxt=ns+ks
       ihi = ilo + 1
       IF (ihi >= lxt) THEN
         IF (x >= t(lxt)) THEN
           mflag = 1
           left = lxt
           RETURN
         END IF      
 
         IF (lxt <= 1) THEN
           mflag = -1
           left = 1
           RETURN       
         END IF
         ilo = lxt - 1
         ihi = lxt
       END IF
 
       IF (x < t(ihi) ) THEN
         IF (x >= t(ilo)) THEN
           mflag = 0
           left = ilo
           RETURN
         END IF        
 
         ! .. now x < t(ilo), decrease  ilo  to capture  x.
 
         istep = 1
         DO 
           ihi = ilo
           ilo = ihi - istep
           IF (ilo <= 1 ) EXIT
           IF ( x >= t(ilo) ) EXIT
           istep = istep*2
         END DO
 
         IF ( ilo <= 1 ) THEN
           ilo = 1
           IF (x < t(1)) THEN
             mflag = -1
             left = 1
             RETURN        
           END IF
         END IF
       ELSE
 
         ! .. now x >= t(ihi), increase ihi to capture x.
 
         istep = 1
         DO
           ilo = ihi
           ihi = ilo + istep
           IF (ihi >= lxt) EXIT
           IF ( x < t(ihi) ) EXIT
           istep = istep*2
         END DO
 
         IF( ihi >= lxt ) THEN
           IF (x >= t(lxt)) THEN
             mflag = 1
             left = lxt
             RETURN      
           END IF
           ihi = lxt
         END IF  
       END IF
 
     ! .. now t(ilo) <= x < t(ihi), narrow the interval.
 
       DO 
         middle = (ilo + ihi)/2
         IF (middle == ilo) THEN 
           mflag = 0
           left = ilo      
           EXIT
         END IF
   
         ! .. note: it is assumed that middle = ilo in case ihi = ilo+1 
   
         IF (x >= t(middle)) THEN
           ilo = middle
         ELSE
           ihi = middle
         END IF
       END DO
     END SUBROUTINE interv
   
   END FUNCTION bvalu2
   
