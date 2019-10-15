
!======================================================================
      Subroutine density(ns,ks,d,p1,p2,type,ms)
!======================================================================
!
!     d - density of two w.f. p1,p2
!     type = 's' - simmetrical direct case
!            'n' - nonsimmetrical direct
!            'x' - exchange case
!----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER, INTENT(in) :: ns,ks
      REAL(KIND=8), DIMENSION(ns) :: p1,p2
      REAL(KIND=8), DIMENSION(ns,*) :: d
      INTEGER, INTENT(out) :: ms
      CHARACTER(LEN=1) :: type

      INTEGER :: i,j,imin,imax, m1, m2
      
      ! .. determine range of orbital
      !Print '(A/(10F8.4))', 'density: p1:', p1
      do i = ns,ns/2,-1
        if (abs(p1(i)) >1.d-10) exit
      end do
      m1 = i
      !Print '(A/(10F8.4))', 'density: p2:', p2
      do i = ns,ns/2,-1
        if (abs(p2(i)) >1.d-10) exit
      end do
      m2 = i
      ms = min(m1,m2)
      
      if(type.eq.'s') then                    !     o***
      
        d(1:ns,1:ks) = 0.d0
                                              !     o***
        do i =1,ms                            !     o***
          d(i,1) =  p1(i)*p2(i)               !     o***
        end do                                !     o**
                                              !     o*
        do j = 2,ks                           !     o
          do i = 1,ms-j+1
           d(i,j) =  p1(i)*p2(i+j-1) + p1(i+j-1)*p2(i)
          end do
        end do

      elseif(type.eq.'n') then                !     o***
                                              !    *o***
        d(1:ns,1:ks+ks-1) = 0.d0


        do j = 1,ks+ks-1                      !   **o***
          imin=max0( 1, 1+ks-j)               !  ***o***
          imax=min0(ms,ms+ks-j)               !  ***o**
          do i = imin,imax                    !  ***o*
            d(i,j) = p1(i)*p2(i+j-ks)         !  ***o
          end do
        end do

      else

        d(1:ns,1:ns) = 0.d0
        
        do i = 1,ms
          do j = 1,ms
            d(i,j) =  p1(i)*p2(j)
          end do
        end do

      end if

      END SUBROUTINE density
