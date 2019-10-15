!====================================================================    
    MODULE spline_integrals
!====================================================================
!
!   contains the B-spline representation of two-electron integral
!   rkb(i,j;i',j') in symmetric or non-symmetric column storage mode:
!
!            rkb(1:ns, 1:ns, 1:2*ks-1, 1:2*ks-1)  changed to
!            rkb(1:ns, 1:ns, 1:ks, 1:ks)  changed to
!
!   itype - character (rk, rk1, rk2, ...) which indicates the type
!           of integral and method of calculation for integral storing
!           in the rkb array  at the moment
!
!   krk   - multipole index for the integral
!
!--------------------------------------------------------------------


    IMPLICIT NONE
    SAVE

    INTEGER(4) :: krk = -100, krk_min, krk_max

    CHARACTER(3) :: itype='aaa'
    
    REAL(8), DIMENSION(:,:,:,:), POINTER :: rkb
    INTEGER, ALLOCATABLE :: irka(:)
    REAL(8), ALLOCATABLE, TARGET :: rka(:,:,:,:,:)

    END MODULE spline_integrals


!====================================================================    
      SUBROUTINE allocate_rka(kmin, kmax)
!====================================================================    
!
! ... allocates space for spline integrals
!
!--------------------------------------------------------------------

      USE spline_param, only: ns,ks
      USE spline_integrals

      Implicit none
      INTEGER, INTENT(in) :: kmin,kmax

      if (associated(rkb)) nullify(rkb)

      if (allocated(irka)) Deallocate (irka,rka)
      
      ALLOCATE(rka(ns,ns,ks,ks,kmin:kmax), irka(kmin:kmax))
      irka = -1; krk_min=kmin; krk_max=kmax; krk=-100

      END SUBROUTINE allocate_rka

!====================================================================    
      SUBROUTINE dealloc_integrals
!====================================================================    
!
! ... deallocates arrays in module "spline_integrals"
!
!--------------------------------------------------------------------

      USE spline_integrals

      if (associated(rkb)) nullify(rkb)

      if(allocated(irka)) DEALLOCATE(irka,rka)

      itype='aaa'
      krk = -100

      END SUBROUTINE dealloc_integrals


