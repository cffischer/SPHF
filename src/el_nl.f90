!======================================================================
      SUBROUTINE el_nl(el,n,l) 
!======================================================================
!
!     Determine 'n' and 'l' from electron string
!----------------------------------------------------------------------
      IMPLICIT NONE
      CHARACTER(LEN=3) , INTENT(INOUT) :: el
      INTEGER, INTENT(OUT) :: n,l
      INTEGER, EXTERNAL :: lval
 
     el = adjustr(el)     ! right justify
     Read (el(1:2), *) n
     l = lval(el(3:3))

     END SUBROUTINE el_nl
      
