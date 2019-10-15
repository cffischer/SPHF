! routines that return slater coefficients
!   a(i,j,k) returns coefficient for Fk(i,j,k)
!   b(i,j,k) returns coefficient for Fk(i,j,k)
!=======================================================================
  FUNCTION a(i,j,k)
!=======================================================================
!
!  determine the coefficient of Fk(i,j)
!
!----------------------------------------------------------------------
!
     USE hf_atomic_state
     USE av_energy
     USE hf_energy_expression

     IMPLICIT NONE
     INTEGER, INTENT(IN) :: i,j,k
     REAL(KIND=8) :: a,c
     INTEGER :: istart

     if (  k <= 2*min(l(i),l(j)) .and. mod(k,2)==0 ) then
       if (i > nclosd .and. j > nclosd) then
         istart = ijptr(i-nclosd,j-nclosd) + 1
         a = coef(istart + k/2)
       else if (i == j) then
         c = qsum(i) - 1.d0
         if (k == 0) then
            a = c
         else
            a = -c*ca(l(i),k)
         end if
       else if (k.eq.0) then
         a = qsum(j)
       else
         a = 0.d0
       end if
       ! .. coef is the coefficient in the differential equation
       !    divided by qsum(i). Here we want the coefficient times
       !    2*qsum(i), 
       a = a*qsum(i)
     else
       a = 0
     end if
  END FUNCTION a
!=======================================================================
  Function b(i,j,k) 
!=======================================================================
!
!   determine the coefficient of Gk(i,j)
!
!----------------------------------------------------------------------
!
     USE hf_atomic_state
     USE av_energy
     USE hf_energy_expression
  
     IMPLICIT NONE
     INTEGER, INTENT(IN) :: i,j,k
     REAL(KIND=8) :: b
     INTEGER :: istart,ll,kk
  
     If ( abs(l(i)-l(j)) <= k .and. k <= l(i)+l(j) ) then 
       if (i == j .or. mod(k-abs(l(i)-l(j)),2) /= 0) then
            b = 0.d0

       else if (i > nclosd .and. j > nclosd) then
      
         !  ll is the number of direct terms
         !  istart the beginning of the exchange terms
         ll = min(l(i),l(j)) + 1
         istart = ijptr(i-nclosd,j-nclosd) + 1 + ll
         kk = (k - abs(l(i)-l(j)))/2
         b = coef(istart + kk)

       else
         b = -qsum(j)*cb(l(i),l(j),k)

       end if
       !  This coefficient is the coefficient in a differential equation.
       !  Here we want 2*qsum(i) times that coefficient
       b = qsum(i)*b
     else 
       b = 0
     end if

   END FUNCTION b
