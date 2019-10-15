!====================================================================
    SUBROUTINE facdyk(ktx,k,ipvtd,dyk)
!====================================================================
!
!   Sets up the matrix
!
!       DD(i,j)-k*(k+1)*rm2
!
!   for solving the following differential equation by the
!   spline-Galerkin method:
!
!       (d^2-k(k+1)/r^2)yk(r) = -(2k+1) B_jf B_js(r)/r,
!
!   with boundary conditions:
!       at r=0, yk(r) = 0;
!       at r=rmax, yk(r) = const, if k=0,
!                  dy/dr + (k/r)y = B_jf(r) B_js(r) if k>0.
!
!   CALL:    dgbfa (LINPACK)  or DGBTRF (LAPACK)
!
!--------------------------------------------------------------------
!
!   on entry
!   --------
!       ktx  the leading dimension of dyk
!       k    the order of y.
!
!   on exit
!   -------
!       dyk  factorized array for the differential with operator
!            (d^2-k(k+1)/r^2). dyk is banded with
!            width 3*ks-2. The first ks-1 rows are zeros, and the
!            following k rows are the elements above the diagnal,
!            the last ks-1 rows are the elements below the diagnal of
!            the original arrays.
!
!--------------------------------------------------------------------

    USE spline_param; USE spline_grid; USE spline_galerkin

    IMPLICIT NONE
    INTEGER(4), INTENT(IN) :: ktx, k
    REAL(8), DIMENSION(ktx,ns) :: dyk
    INTEGER(4), DIMENSION(ns) :: ipvtd

    REAL(8) :: fkk
    INTEGER(4) :: i, j, nn, ier

    fkk = -k*(k+1)
    dyk=0.d0              ! .. clear the dyk array

    ! .. set up dyk

    ! .. lower portion

    do j=1,ks
      do i=ks-j+1,ns
        dyk(3*ks-1-j,i-ks+j)=db2(i,j)+fkk*rm2(i,j)
      end do
    end do

    ! .. upper portion

    do j=2,ks
      do i=1,ns-j+1
        dyk(2*ks-j,i+j-1)=db2(i+j-1,ks-j+1)+fkk*rm2(i+j-1,ks-j+1)
      end do
    end do

    ! .. correct the upper portion of the matrix for the asymmetry
    ! .. of d^2/dr^2

    j = 2*ks-2
    dyk(j,ns) = dyk(j,ns) + db2(1,1) - db2(ns,ks-1)

    ! .. apply the zero boundary condition at the origin and modify the
    ! .. last equation for the boundary condition at rmax

    do i=1,ks
      j=2*ks-i
      dyk(j,i)=0.d0
      j=3*ks-1-i
      dyk(j,ns-ks+i)=0.d0
    end do
    dyk(2*ks-1,1)=1.d0

    ! .. apply the boundary condition at rmax
    ! .. yk(r) = const, if k=0,
    ! .. dy/dr + (k/r)y = B_jf(r) B_js(r) if k>0.

    if (k /= 0) then
      nn = nv+1
      dyk(2*ks,ns-1)=bspd(nn,1,ks-1,1)
      dyk(2*ks-1,ns)=bspd(nn,1,ks,  1) + bsp(nn,1,ks)*k/t(ns+1)
    else
      dyk(2*ks-1,ns) = 1.d0
    end if

    ! .. factorize dyk

!    CALL dgbfa(dyk,ktx,ns,ks-1,ks-1,ipvtd,ier)  
!    IF (ier /= 0) STOP 'FACDYK: dgbtrf from LINPACK failed'

     Call DGBTRF (ns,ns,ks-1,ks-1,dyk,ktx,ipvtd,ier)
     if (ier .ne. 0) Stop 'FACDYK: dgbtrf from LAPACK failed'

  END SUBROUTINE facdyk
