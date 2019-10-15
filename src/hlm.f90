!  Routine hlm calls
!  .. mvc  (relativistic correction)
!====================================================================
    SUBROUTINE hlm(L)
!====================================================================
!
!   Sets up the matrix hl for the coulomb operator
!
!       hl(i,j) = DD(i,j) - l*(l+1)*rm2(i,j) + 2 * z*rm1(i,j)
!
!   in symmetric band storage mode, where i=1,ns, j=1,ks and
!   hl(1,1) stores <B_{ns-1},H B_{ns}>.
!
!--------------------------------------------------------------------
!   on entry
!   --------
!       z    nuclear charge
!       l    orbital angular momentum
!--------------------------------------------------------------------

    USE spline_param
    Use hf_atomic_state, ONLY: z, rel, fine
    USE spline_grid
    Use spline_galerkin
    Use spline_hl

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: L

    REAL(KIND=8) :: fl, zz, a
    Real(8), External :: AZL

    if(lh.eq.-1) Call Allocate_hl
    if(lh.eq.l) Return

    fl = L*(L+1.d0)
    zz = 2.d0*z

!.. set up hl matrix

    hl = db2 - fl*rm2 + zz*rm1

!.. store the (ns-1,ns) asymmetric value in hl(1,1)

    hl(1,1) = hl(ns,ks-1) + (db2(1,1)-db2(ns,ks-1))

!.. relativistic shift

    if(rel) then
     Call mvc(l)
     hl = hl + fine*vc
     if(l.eq.0) then
      a = azl(z,h,ks,l+1)
      hl(2,ks) = hl(2,ks) - z*a*a*fine
     end if
    end if

    lh = l
  END SUBROUTINE hlm
!====================================================================
    SUBROUTINE mvc(L)
!====================================================================
!
!   Computes the matrix elements for the mass-velocity correction
!   in the B-spline basis.  The lack of symmetry in the d^2/dr^2
!   operator is ignored.
!
!     VC(i,j) = INT [  (d^2/dr^2 - l(l+1)/r) B_i(r) *
!                      (d^2/dr^2 - l(l+1)/r) B_j(r)  ] dr
!--------------------------------------------------------------------
!
!   on entry
!   --------
!       L    the angular momentum
!
!   on exit
!   -------
!       vc   the mass velocity correction in symmetric storage mode
!
!--------------------------------------------------------------------

    USE spline_param
    USE spline_grid
    !USE hf_atomic_state, ONLY: z, 
    USE spline_hl

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: L

    ! .. local variables

    INTEGER :: m, ith, jth, i, irow, jcol
    REAL(KIND=8) :: fll, y1, y2, SS

    ! .. initialize the vc array

    vc = 0.d0
    fll  =  L*(L+1)

    ! .. compute the matrix elements

    do m = 1,ks
      do i = 1,nv
        SS = fll*grm(i,m)*grm(i,m)

! ... cutoff correction

!        B = gr(i,m)/(gr(i,m)+2*fine*Z)
!        B = B*B*B

        do ith = 1,ks
          irow = i+ith-1
          do jth = 1,ith
          jcol = jth-ith+ks

            y1 = bspd(i,m,ith,2) - SS*bsp(i,m,ith)
            y2 = bspd(i,m,jth,2) - SS*bsp(i,m,jth)
            vc(irow,jcol) = vc(irow,jcol) + grw(i,m)*y1*y2 !*B

          end do
        end do
      end do
    end do

    END SUBROUTINE mvc

