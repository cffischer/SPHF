!====================================================================
    MODULE spline_moments
!====================================================================
!                                                          _
!   contains moments defining as <B_i|r^k|B_j> or <B_i|r^k|B_j> 
!   over an interval, where  _
!                            B = B' - B/r
!
!   These moments are used for calculating of two-electron integrals
!   according to the sell algorithm.
!
!--------------------------------------------------------------------
!
!   rkd(1:ks*ks;1:ks*ks,1:nv) - 
!   
!      the two-dimensional array of integrals <B_i B_j|...|B_i' B_j'>
!      over a triangle (or square) diagonal cell
!     
!   rkd[1,2,3,4](1:ks*ks,1:nv) - different moments defining as
!
!     <B_i|r^k|B_j>  and <B_i|r^k|B_j>  over an interval iv
!    
!   rkd, rkd1, ... differ from rkt, rkt1, ... in module spline_moments)
!   only by reduced dimensions, that increases slightly the speed of
!   calculations
!--------------------------------------------------------------------

    IMPLICIT NONE
    SAVE

    INTEGER :: kmk = -100
    CHARACTER(3) :: mtype = 'aaa'

    REAL(8), ALLOCATABLE, DIMENSION(:,:,:) :: rkd
    REAL(8), ALLOCATABLE, DIMENSION(:,:) :: rkd1, rkd2, rkd3, rkd4

    END MODULE spline_moments



!====================================================================
      SUBROUTINE allocate_moments
!====================================================================
!
! ... allocates space for arrays in MODULE spline_moments
!
!--------------------------------------------------------------------

      USE spline_param
      USE spline_moments
   
      INTEGER :: jk

      if(allocated(rkd)) DEALLOCATE(rkd, rkd1,rkd2,rkd3,rkd4)

      jk = ks*ks
      ALLOCATE(rkd(jk,jk,nv), rkd1(jk,nv),rkd2(jk,nv), &
                              rkd3(jk,nv),rkd4(jk,nv))
      rkd = 0.d0; rkd1 = 0.d0; rkd2 = 0.d0; rkd3 = 0.d0; rkd4 = 0.d0
      mtype='bbb'
      kmk=-100

      END SUBROUTINE allocate_moments


!====================================================================
      SUBROUTINE dealloc_moments
!====================================================================
!
!     deallocate space in MODULE spline_moments
!
!--------------------------------------------------------------------

      USE spline_moments
   
      if(allocated(rkd)) DEALLOCATE(rkd, rkd1,rkd2,rkd3,rkd4)
      mtype='aaa'
      kmk=-100

      END SUBROUTINE dealloc_moments


