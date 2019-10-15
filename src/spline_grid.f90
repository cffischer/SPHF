!====================================================================
     MODULE spline_grid
!====================================================================
!
!    defines the values of splines at the gaussian points for each
!    interval of a grid; included in the module is the gaussian data
!    for performing integrations on the grid
!
!--------------------------------------------------------------------

     IMPLICIT NONE
     SAVE

! .. knot sequence, t(1:ns+ks)

     REAL(8), DIMENSION(:), ALLOCATABLE:: t

! .. arrays for spline values in gausian points
!
!    bsp(1:nv+1,1:ks,1:ks), bspd(1:nv+1,1:ks,1:ks,2) 
!
!    bsp(i,m,ith)  - values of the i+ith-1 B-spline in interval i
!                    at gausian point m
!    bspd(i,m,ith,1|2)  - corresponding values of first and second
!                         derivatives 
!    bsp(nv+1,1,.) and bspd(nv+1,1,.,.) - corresponding values at
!                                         last knot point (rmax)

     REAL(8), DIMENSION(:,:,:), ALLOCATABLE::   bsp, bsq
     REAL(8), DIMENSION(:,:,:,:), ALLOCATABLE:: bspd

! .. arrays for gaussian data
!
!    gr(1:nv;1:ks),  grm(1:nv;1:ks),  grw(1:nv;1:ks)  
!
!    gr(i,m)  - gaussian points m in the interval i
!    grm(i,m) - reciprocal value of gr(i,m)
!    grw(i,m) - gaussian weights at corresponding points

     REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE:: gr, grm, grw


     CONTAINS


!====================================================================
      SUBROUTINE allocate_grid
!====================================================================
!
! ... allocates space of the arrays in MODULE spline_grid
!
!--------------------------------------------------------------------

      USE spline_param

      if(Allocated(bsp)) Deallocate(bsp,bsq,bspd,gr,grm,grw)
      ALLOCATE(bsp(nv+1,ks,ks), bsq(nv+1,ks,ks), bspd(nv+1,ks,ks,2), &
               gr(nv,ks), grm(nv,ks), grw(nv,ks))
      bsp = 0.d0; bsq = 0.d0; bspd = 0.d0
	  gr = 0.d0; grm = 0.d0; grw = 0.d0

      END SUBROUTINE allocate_grid

    END MODULE spline_grid

