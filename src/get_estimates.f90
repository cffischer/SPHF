!=====================================================================
      SUBROUTINE get_estimates
!=====================================================================
!     Get initial estimates i) reading bsw.inp or 2) screened hydrogenic
!  Read bsw.inp for initial estimates of radial function and
!  use screened hydrogenic functions for those not found
!  Atomic number and grid must match
    !------------------------------------------------------------------
      USE spline_param, ONLY: ns,ks
      USE spline_galerkin, ONLY: sb
      USE hf_atomic_state
      USE hf_inout
      USE hf_orbitals
      USE hf_param
      
      IMPLICIT NONE

!
       REAL(KIND=8) :: ss
      INTEGER :: north, i, j, m
      INTEGER, DIMENSION(nwf) :: jorth

      CHARACTER(LEN=6), DIMENSION(nwf) :: ATM, TRM

      in = 0
    
      ! .. Read bwfn.inp for fixed orbitals or initial estimates
      call read_bsw
      ! Assume 
      !  .. fixed orbitals are from 1 to nwf-nit (at the beginning)
      !  .. have the same spline parameters 
      !  .. in(i) = -1 if same z (accepte); in(i) = 1 if scaled.

      ! Determine the overlap matrix
      ! .. generate the b matrix (could be part of a  module)

      Print *
      ss = 0.d0; atm(1:nwf)=atom; trm(1:nwf)=term
      Do i = 1, nwf
         ss = ss + qsum(i)
         IF (in(i) == -1) CYCLE  ! We have an initial estimate
         !  .. Let m be the index of the eigenvalue AFTER orthogonality
         !  .. constraints  have been applied.
         m = n(i) -l(i)
         !.. set orthogonality constraints by forming the list jorth
         north=0
         do j = 1,i-1
           if (l(i) == l(j)) then
              north = north+1
              jorth(north) = j
              m = m-1
           end if
         end do
         if (in(i) == 0) then
            call screened_hydrogenic(i,z-ss+(qsum(i)+1)/2)
            ! this option is for testing with hydrogenic orbitals
            !call screened_hydrogenic(i,z)
            in(i) = 1
         end if
         !.. may wish to change this
         do j = ns, ns/3,-1
            if (abs(p(j,i)) > end_tol ) exit
         end do
         maxr(i)=min(ns,j+1)
         p(j+1:ns,i) = 0.d0
         dpm(i) = 1.d0
         s(i)   = ss-(qsum(i)+1)/2
         e(i,i) = -((z-s(i))/n(i))**2/2
      end do
! We now have an orthonormal set of initial estimates
      ! call energy   ! provides initial estimate of the total energy
      ! write(scr,'(7X,A,F20.12)') &
      !             'Total energy from initial estimates =', etotal
      
     ! Do i = 1,nwf
     !   Write(50,'(10X,A3,F10.5,I4, 2F10.5)') el(i), Z, maxr(i), e(i,i), dpm(i)
     !   Write(50,'(8F13.8)') p(:,i)
     ! End do

    CONTAINS
    !==================================================================
        SUBROUTINE screened_hydrogenic(i,zz)
    !==================================================================
    !   Compute the screened hydrogenic function for orbital i    
    !   that satisfies orthogonality constraints.
    !------------------------------------------------------------------
          INTEGER, INTENT(in) :: i
          REAL(KIND=8), INTENT(in) :: zz

          INTEGER :: j, jp
          REAL(KIND=8), DIMENSION(ns) :: u,y

          !print '(A,3I4,F16.9)', 'Entering screened ..', i, n(i), l(i), zz
          call bhwf(n(i), l(i), zz, u)
          call bxv(ks,ns,sb,u,y)
          do jp = 1, north
             j = jorth(jp)
             u = u - dot_product(y, p(:,j))*p(:,j)
          end do
          call normalize(u)
          p(:,i) = u
        END SUBROUTINE screened_hydrogenic
   END subroutine get_estimates
