!======================================================================
    SUBROUTINE define_spline
!======================================================================
!
!   initializes the values of the spline and its derivatives
!   and evaluates the spline basic arrays (elementary operators in
!   spline basis).
!
!   SUBROUTINE called:
!       gauss
!       allocate_grid
!       initvb
!       allocate_galerkin
!       initas
!
!   calling sequence:
!                        define_spline
!                   ----------------------
!                  / |          ||      ||
!                 /  |        initvb   initas
!                /   |           |     // \  \\
!           gauss    |        vbsplvd mdb mrm facsb
!                    |           ||
!       allocate_grid,galerkin vbsplvb
!
!----------------------------------------------------------------------

    USE spline_param
    USE spline_grid
    USE spline_galerkin

    IMPLICIT NONE
    REAL(8), DIMENSION(ks) :: gx, gw

    ! .. initializes variables for gaussian integration

    CALL gauss(ks,gx,gw)

    ! .. initializes the values of the spline and its derivatives
    
    CALL allocate_grid
    
    CALL initvb
    
    ! .. initializes the spline array (operators in spline basis)

    Call allocate_galerkin

    CALL initas

    CONTAINS

!=======================================================================
    SUBROUTINE initvb
!=======================================================================
!
!   Sets (or Initializes) the arrays
!       gr      The gaussian points for the nint intervals of [0,Rmax]
!       grm     Reciprocals of the values of gr
!               (to avoid repeated division on the CRAY)
!       grw     Gaussian weights for each of the points in gr
!       bsp     array of B-spline values at the gaussian points
!       bspd    array of values of the first and second derivative of
!               the B-splines at the gaussian points
!
!   Calling sequence:
!        initvb
!          |
!       vbsplvd
!         ||
!       vbsplvb
!
!-----------------------------------------------------------------------
!   on entry
!   --------
!       gx      the gaussian points for the interval [0,1]
!       gw      the gaussian weights for the interval [0,1]
!
!   working arrays
!   --------------
!       dbiatx  working array of dimension (ns,ks,ks) which contains the
!               values and the second derivative values of of B-spline at
!               gaussian points in the interval [0,1]
!
!-----------------------------------------------------------------------

    IMPLICIT NONE
    REAL(8), Allocatable, DIMENSION(:,:,:) :: dbiatx
    INTEGER :: m, i
    Allocate(dbiatx(nv,ks,ks)); dbiatx = 0.d0 

    Do m=1,ks
      gr(1:nv,m)=(t(1+ks:nv+ks)-t(ks:nv+ks-1))*gx(m)+t(ks:nv+ks-1)
      grm(1:nv,m) = 1.d0/gr(1:nv,m)
      grw(1:nv,m)=(t(1+ks:nv+ks)-t(ks:nv+ks-1))*gw(m)
      Call vbsplvd(t,ks,nv,gr(1,m),3,dbiatx)
      bsp(1:nv,m,1:ks)    = dbiatx(1:nv,1:ks,1)
      bspd(1:nv,m,1:ks,1) = dbiatx(1:nv,1:ks,2)
      bspd(1:nv,m,1:ks,2) = dbiatx(1:nv,1:ks,3)
      Do i = 1,ks
       bsq(1:nv,m,i) = bspd(1:nv,m,i,1) - grm(1:nv,m)*bsp(1:nv,m,i)
      End do
    End do

    ! .. store also the values at the last knot

    call vbsplvd(t,ns,1,t(ns+1),3,dbiatx)

    bsp(nv+1,1,1:ks)    = dbiatx(1,1:ks,1)
    bspd(nv+1,1,1:ks,1) = dbiatx(1,1:ks,2)
    bspd(nv+1,1,1:ks,2) = dbiatx(1,1:ks,3)
    bsq(nv+1,1,1:ks)    = bspd(nv+1,1,1:ks,1)-bsp(nv+1,1,1:ks)/t(ns+1)

    Deallocate(dbiatx) 

    END SUBROUTINE initvb


!==================================================================
    SUBROUTINE initas
!==================================================================
!
!   Sets ( or Initializes ) the array in symmetric storage mode:
!
!       db1 --- matrix of integral <B_i,B'_j>
!       db2 --- matrix of integral <B_i,B"_j>
!       sb  --- matrix of integral <B_i,B_j>
!       bb  --- full matrix version of sb
!       r1  --- matrix of integral <B_i,r B_j>
!       rm1 --- matrix of integral <B_i,(1/r)B_j>
!       rm2 --- matrix of integral <B_i,(1/r^2)B_j>
!       rm3 --- matrix of integral <B_i,(1/r^3)B_j>
!               where i=1,..,ns, j=1,...ks
!
!   Calling sequence:
!       initas
!        /  \
!      mdb  mrm
!
!------------------------------------------------------------------


    IMPLICIT NONE


    CALL mdb1      ! .. sets db1 --- matrix of integral <B_i,B'_j>

    CALL mdb2      ! .. sets db2 --- matrix of integral <B_i,B"_j>

    CALL mrm(0, sb)              ! .. sets sb  ---  <B_i,B_j>

    CALL add_band_to_full(sb,bb) !    sets full matrix version of sb

    CALL facsb                   ! .. factorizes sb

    CALL mrm(1, r1)              ! .. sets r1  ---  <B_i,r B_j>

    CALL mrm(-1, rm1)            ! .. sets rm1 ---  <B_i,(1/r)B_j>

    CALL mrm(-2, rm2)            ! .. sets rm2 ---  <B_i,(1/r^2)B_j>

    CALL mrm(-3, rm3)            ! .. sets rm2 ---  <B_i,(1/r^3)B_j>

   END SUBROUTINE initas


!====================================================================
   SUBROUTINE mdb1
!====================================================================
!
!  Computes the matrix elements <B_i,B'_i+j-1> in the B-spline basis
!
!--------------------------------------------------------------------

    IMPLICIT NONE

    ! Local variables

    INTEGER :: i, irow, jcol, ith, jth

    db1 = 0.d0

    do ith = 1,ks
      do jth = 1,ith
        jcol = jth-ith+ks
        do i = 1,nv
          irow = i+ith-1
          db1(irow,jcol) = db1(irow,jcol) &
        + SUM(grw(i,:)*bsp(i,:,ith)*bspd(i,:,jth,1))
        end do
      end do
    end do

  END SUBROUTINE mdb1


!====================================================================
   SUBROUTINE mdb2
!====================================================================
!
!  Computes the matrix elements <B_i,B"_i+j-1> in the B-spline basis
!
!--------------------------------------------------------------------

    IMPLICIT NONE

    ! Local variables

    INTEGER :: i, irow, jcol, ith, jth

    db2 = 0.d0

    do ith = 1,ks
      do jth = 1,ith
        jcol = jth-ith+ks
        do i = 1,nv
          irow = i+ith-1
          db2(irow,jcol) = db2(irow,jcol) &
        + SUM(grw(i,:)*bsp(i,:,ith)*bspd(i,:,jth,2))
        end do
      end do
    end do

    db2(1,1) = db2(1,1) + SUM(grw(nv,:)*bsp(nv,:,ks-1)*bspd(nv,:,ks,2))

  END SUBROUTINE mdb2


  END SUBROUTINE define_spline

!====================================================================
    SUBROUTINE facsb
!====================================================================
!
!   Factorizes bs matrix which is a transpose of overlap matrix sb,
!   <B_i,B_j>,  with the correct boundary condition at r=0 
!
!   SUBROUTINES called:  dpbtrf (from LAPACK)
!
!--------------------------------------------------------------------

    USE spline_param
    USE spline_grid
    USE spline_galerkin

    IMPLICIT NONE
    INTEGER :: m, ierr

    ! .. copy the array, converting to row oriented band storage mode

    bs = TRANSPOSE(sb)

    ! .. apply boundary condition at r=0

    do m = 1,ks-1
      bs(m,ks-m+1)=0.d0
    end do
    bs(ks,1) = 1.d0

    ! .. apply boundary condition at r=rmax
    bs(1:ks-1,ns) = 0.d0
    bs(ks,ns)=1.d0

    Call DPBTRF('U',ns,ks-1,bs,ks,ierr)
    if (ierr.ne.0)  Stop 'facsb: dpbtrf (LAPACK) failed'

  END SUBROUTINE facsb
