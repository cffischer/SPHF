!=======================================================================
  SUBROUTINE vhwf(n,l,z,nr,r,vh,ierr)  
!=======================================================================
!
!   This program returns the vector of values, vh(i), of a normalized 
!   hydrogenic function with nuclear charge z and quantum numbers (n,l)
!   at values of the radius, r(i),  i=1,nr.
!
!   SUBROUTINE contained:  hnorm
!
!-----------------------------------------------------------------------
!
!   on entry
!   --------
!
!       n     principal quantum number
!       l     angular quantum number
!       z     Nuclear charge
!       nr    number of points in the vector
!       r     values at which radial function is to be evaluated, in
!             increasing order
!
!   on exit
!   -------
!
!       vh    values of the hydrogenic radial functions
!       ierr  error condition; 0 for normal return; 1 otherwise
!
!-----------------------------------------------------------------------
!
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n, l, nr
    INTEGER, INTENT(OUT) :: ierr
    REAL(KIND=8), INTENT(IN) :: z
    REAL(KIND=8), DIMENSION(nr), INTENT(IN) :: r
    REAL(KIND=8), DIMENSION(nr), INTENT(OUT) :: vh

    ! .. Local variables

    INTEGER :: k, i, nm
    REAL(KIND=8) :: a, b, c, factor
    REAL(KIND=8), DIMENSION(nr) :: w, p

    CALL hnorm       ! .. gets the normalization factor
    k = n-l-1  
 
    ! .. store the argument of the exponential factor
 
    w = -2.D0*z*r/n  
 
    ! .. to avoid underlow, determine point at which function will be zero
    ! .. Note that exp(-150) = 0.7175 10**-65
 
    nm = nr  
    do
      if (w(nm) < -150.D0 .AND. nm /= -1) then  
        vh(nm) = 0.D0  
        nm = nm-1  
      else
        EXIT
      end if  
    end do
 
    ! .. Initialize the recurrence relation for each value of r(i)
 
    p(1:nm) = 1.D0  
 
    a = 1.D0  
    b = k  
    c = n+l  
    ierr = 0  
 
    SELECT CASE (k)
 
    CASE(:-1)      ! .. all values below 0
      ierr = 1  
    RETURN  

    CASE (1:)      ! .. all values above 0
 
      ! .. Apply the recurrence relation when k is positive
 
      do I = 1,k  
        p(1:nm) = 1.D0+a/b*p(1:nm)/c*w(1:nm)  
        a = a+1.D0  
        b = b-1.D0  
        c = c-1.D0  
      end do
    END SELECT
 
    ! .. Multiply the factor by the exponential
 
    vh(1:nm) = factor*p(1:nm)*EXP(w(1:nm)/2.D0)*(-w(1:nm))**(l+1)  
 
    CONTAINS  

    !====================================================================
      SUBROUTINE hnorm
    !====================================================================
    !
    !   returns the value of the normalization constant for an 
    !   hydrogenic function with nuclear charge z, and orbital
    !   quantum numbers (n,l)
    !
    !--------------------------------------------------------------------
    !
        IMPLICIT NONE

        ! .. local variables
        INTEGER :: m, i
        REAL(KIND=8) :: a, b, d, t

        m = l+l+1
        a = n+l
        b = m
        t = a
        d = b
        m = m-1

        if ( m > 0) then
          do  i = 1,m
            a = a-1.d0
            b = b-1.d0
            t = t*a
            d = d*b
          end do
        end if
        factor = SQRT(z*t)/(n*d)
      END SUBROUTINE hnorm

  END SUBROUTINE vhwf
