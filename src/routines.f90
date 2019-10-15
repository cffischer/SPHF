!     ==================================================================
       subroutine add_rkm_d(k,d,hd)
!     ==================================================================
!
!	 Adds RK(.a,.b) to the matrix hd in symmetric banded form given 
!	 the density matrix d in symmetric banded form
!            
!     ------------------------------------------------------------------
         USE spline_param, ONLY: ns, ks
         USE spline_integrals

         IMPLICIT NONE
         INTEGER(4), INTENT (in) :: k
         REAL(KIND=8), DIMENSION(ns,ks), INTENT(in) :: d
         REAL(KIND=8),DIMENSION(ns,ks),INTENT(inout):: hd
         INTEGER :: i,j,jth,ith

         Call mrk_cell(k)
	 ! Form the contribution in band matrix form
         do jth = 1, ks
            do j = 1, ns - jth + 1
               do ith = 1, ks
                  do i = 1, ns - ith + 1
                     hd(i+ith-1,ks-ith+1) = hd(i+ith-1,ks-ith+1) &
                         + d(j,jth)*rkb(i,j,ith,jth)
                  end do
               end do
            end do
         end do

       END Subroutine add_rkm_d

!     ==================================================================
       subroutine add_band_to_full(band,full)
!     ==================================================================
!
!	Add the matrix hd in symmetric band form to 
!       the matrix hfm in the symmetrci full matrix form
!
!     ------------------------------------------------------------------
         USE spline_param, ONLY: ns, ks
         
         IMPLICIT NONE
         REAL(KIND=8), DIMENSION(ns,ks), INTENT(in) :: band
         REAL(KIND=8), DIMENSION(ns,ns), INTENT(inout) :: full
         INTEGER :: j,jth

        ! Add to hm
         Do jth = 1,ks-1
           Do j = ks-jth+1,ns
             full(j,j-ks+jth) = full(j,j-ks+jth) + band(j,jth)
             full(j-ks+jth,j) = full(j-ks+jth,j) + band(j,jth)
           end do
         end do

         do j = 1,ns
            full(j,j) =full(j,j) + band(j,ks)
         end do
        
       END Subroutine add_band_to_full
      

      !==================================================================
      subroutine add_rkm_x(k,dx,hfm)
      !==================================================================
      ! 
      !     Add the full symmetric matrix RK(.,a;b,.) to the  full 
      !     symmetric matrix hfm given the density dx as a symmetric
      !     banded matrix
      !-----------------------------------------------
      USE spline_param
      USE spline_integrals

      IMPLICIT NONE
      INTEGER(4), INTENT (in) :: k
      REAL(KIND=8), DIMENSION(ns,ns), INTENT(in):: dx
      REAL(KIND=8), DIMENSION(ns,ns), INTENT(inout):: hfm
      INTEGER :: i,j,ith,jth, ms
        
      Call mrk_cell(k)

!     do j = ns,ns/2,-1
!       if (abs(dx(j,j)) > 1.d-18) exit
!     end do
!     ms = j
      ms= ns

!  .. j' >= j
      do jth = 1, ks 
!       .. i' >= i
         do ith = 1, ks 
            do j = 1, ms - jth + 1 
               do i = 1, ms - ith + 1 
                hfm(i,j) = hfm(i,j) + &
		           rkb(i,j,ith,jth)*dx(i+ith-1,j+jth-1) 
               end do 
            end do 
         end do 
!       .. i'<i
         do ith = 2, ks 
            do j = 1, ms - jth + 1 
               do i = ith, ms 
                  hfm(i,j) = hfm(i,j) + rkb(i-ith+1,j,ith,jth)* &
                                      dx(i-ith+1,j+jth-1) 
               end do 
            end do 
         end do 
      end do 
!
!  .. j' < j
      do jth = 2, ks 
!       .. i' >= i
         do ith = 1, ks 
            do j = jth, ms 
               do i = 1, ms - ith + 1 
                  hfm(i,j) = hfm(i,j) + rkb(i,j-jth + 1,ith,jth)* &
                                      dx(i+ith-1,j-jth+1) 
               end do 
            end do 
         end do 
!       .. i'<i
         do ith = 2, ks 
            do j = jth, ms 
               do i = ith, ms 
                  hfm(i,j) = hfm(i,j)+rkb(i-ith+1,j-jth+1,ith,jth) &
                                     *dx(i-ith+1,j-jth+1) 
               end do 
            end do 
         end do 
      end do 
!     print *, 'hfm'
!     do 2000 i = 1,ns
!       print '(I5/(6f12.7))', hfm(i:1:ns)
! 2000   continue
      end subroutine add_rkm_x

!==================================================================
      subroutine apply_bc(aa,type)
!==================================================================
! 
!   Apply boundary conditions   
!       
!   type = n : zero, aa(1,1) = aa(ns,ns) =1 
!   type = o : zero, aa(1,1) = aa(ns,ns) =0 
!   type = e : zero, aa(1,1) = 1.d20, aa(ns,ns) = -1.d20
!                    (suitable for an eigenvalue problem)
!------------------------------------------------------------------
      USE spline_param

      IMPLICIT NONE
      REAL(KIND=8), DIMENSION(ns,ns), INTENT(inout) :: aa
      CHARACTER(LEN=1), INTENT(in) :: type
    
      If (type .eq. 'n') then
        aa(1,1:ns) =  0.d0
        aa(1:ns,1) =  0.d0
        aa(ns,1:ns) = 0.d0
        aa(1:ns,ns) = 0.d0
        aa(1,1) =     1.d0
        aa(ns,ns) =   1.d0
        aa(ns-1,1:ns-2) = 0.d0
        aa(1:ns-2,ns-1) = 0.d0
        aa(ns-1,ns-1)=    1.d0
      else if (type .eq. 'o') then
        aa(1,1:ns) =  0.d0
        aa(1:ns,1) =  0.d0
        aa(ns,1:ns) = 0.d0
        aa(1:ns,ns) = 0.d0
        aa(1,1) =     0.d0     
        aa(ns,ns) =   0.d0
      else if (type .eq. 'e') then
        aa(1,1:ns) =  0.d0
        aa(1:ns,1) =  0.d0
        aa(ns,1:ns) = 0.d0
        aa(1:ns,ns) = 0.d0
        aa(1,1) =    -1.0d20     
        aa(ns,ns) =   1.0d20
        aa(ns-1,1:ns-2) = 0.d0
        aa(1:ns-2,ns-1) = 0.d0
        aa(ns-1,ns-1)=    1.d020
      else
        Print *, type, ' type of boundary condition not defined'
        STOP
      end if
      END SUBROUTINE apply_bc
!======================================================================
  SUBROUTINE print_m( m,n, header)
!======================================================================
!   This program prints the matrix m (either band or full)
!----------------------------------------------------------------------
!
    USE spline_param

    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(in) :: header
    INTEGER(4), INTENT(in) :: n
    REAL(KIND=8), DIMENSION(ns,n), INTENT(in) ::m
 
    ! .. local variables
    INTEGER :: j

    PRINT *, header, '  From print_m', size(m,dim=2)

    do j = 1,n
      WRITE(6,'(I3/(8f10.7))') j,m(:,j)
    end do
    PRINT *

  END SUBROUTINE print_m

!======================================================================
  SUBROUTINE normalize(v)
!======================================================================
!   This program normalizes the vector v so the v^t B v = 1
!----------------------------------------------------------------------
!
    USE spline_param
    USE spline_galerkin, ONLY: sb

    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(ns), INTENT(inout) :: v

    REAL(KIND=8), DIMENSION(ns) :: y
   
    call bxv(ks,ns,sb,v,y)
  
    v = v/sqrt(Dot_Product(v,y))
    
    ! Define orbitals to be positive near the nucleus
    if (v(ml) < 0.d0) v = -v

  END SUBROUTINE normalize
!======================================================================
  SUBROUTINE orthonormalize(nw)
!======================================================================
!   This program orthonormalize the set of nw orbitals
!----------------------------------------------------------------------
!
    USE spline_param
    USE spline_galerkin, ONLY: sb
    USE hf_atomic_state
    USE hf_orbitals

    IMPLICIT NONE
    INTEGER(4) :: nw

    REAL(KIND=8), DIMENSION(ns) :: y
    REAL(KIND=8) :: c
    INTEGER(4) :: i,j
   
    do i=1,nw
       call bxv(ks,ns,sb,p(:,i),y)
       do j = 1,i-1
          if (l(i) == l(j)) then
             c = Dot_Product(p(:,j),y)
             p(:,i) = p(:,i) - c*p(:,j)
          end if
       end do
       call normalize(p(:,i))
    end do
    END SUBROUTINE orthonormalize
