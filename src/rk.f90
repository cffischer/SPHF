!  Routines for rk, fk, and gk
!======================================================================
      REAL(8) FUNCTION rk (i1,j1,i2,j2,k)
!======================================================================
!               k
!     Returns  R (i1, j1; i2, j2) base on the assembling the B-spline
!     integrals, which are supposed to be placed in the module
!     spline-integrals. If not, they are calculated by programs
!     mrk_diff or mrk_cell, depending of the parameter 'meth':
!     meth = 'd' - differential equation method
!          = 'c' - cell integration method
!----------------------------------------------------------------------

      USE spline_param
      USE hf_orbitals
      USE spline_integrals
  
      IMPLICIT NONE
      INTEGER(4), INTENT(in) :: i1,j1,i2,j2,k
  
      ! .. local variables
  
      INTEGER(4) :: i,ip, j,jp, msi, msj
      REAL(8), DIMENSION(ns,ks) :: a,b
      REAL(8) :: rkj
  
      ! .. check the B-spline integrals in module spline-integrals
  
      if(k.ne.krk.or.itype.ne.'rk ') then
         Call MRK_cell(k)
      end if
  
      ! .. form cross-products
  
      Call density (ns,ks,a,p(1,i1),p(1,i2),'s',msi)
      Call density (ns,ks,b,p(1,j1),p(1,j2),'s',msj)
  
      ! .. assembling the B-spline integrals

      rk = 0.d0
      do ip = 1,ks
        do i = 1,msi-ip+1
          rkj = 0.d0
          do jp = 1,ks
            do j = 1,msj-jp+1
              rkj = rkj+b(j,jp)*rkb(j,i,jp,ip)
            end do
          end do
          rk = rk + a(i,ip)*rkj
        end do
      end do
  
      END FUNCTION rk
  
!======================================================================
  REAL (KIND=8) FUNCTION fk(i,j,k)
!======================================================================
!                            k
!      Returns the value of F (i,j) based on the assembly of B-spline
!     integrals, which are supposed to be placed in the module
!     spline-integrals. 
!---------------------------------------------------------------------
!

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: i,j,k
    REAL (KIND=8), EXTERNAL :: rk

    fk = rk(i,j,i,j,k)

  END FUNCTION fk
!======================================================================
  REAL (KIND=8) FUNCTION gk(i,j,k)
!======================================================================
!                            k
!      Returns the value of F (i,j) based on the assembly of B-spline
!     integrals, which are supposed to be placed in the module
!     spline-integrals. 
!---------------------------------------------------------------------
!

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: i,j,k
    REAL (KIND=8), EXTERNAL :: rk

    gk = rk(i,j,j,i,k)

  END FUNCTION gk
