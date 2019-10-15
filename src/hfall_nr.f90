!==================================================================
       SUBROUTINE hfall_nr(nd)
!==================================================================
!    Improve the estimate v by the Newton_Raphson method
!    subject to orthogonality.
!------------------------------------------------------------------
      USE spline_param
      USE spline_galerkin
      USE hf_atomic_state
      USE hf_orbitals
      USE hf_energy_expression
      USE hf_param, ONLY: end_tol
      USE lapack_interface

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: nd

      REAL(KIND=8), DIMENSION(ns,ns,nwf) :: hfm
      REAL(KIND=8), DIMENSION(ns):: v
      REAL(KIND=8), DIMENSION(nwf) :: ek

      REAL(KIND=8), DIMENSION(nd,nd) :: aa
      REAL(KIND=8), DIMENSION(nd)    :: res
      INTEGER, DIMENSION(nd)         :: ipiv
      LOGICAL :: modified

      REAL(KIND=8), EXTERNAL         :: a, b,  bvmv
      REAL(KIND=8), DIMENSION(ns)    :: x, y
      REAL(KIND=8), DIMENSION(ns,ns) :: dx, dsumx, hx
      REAL(KIND=8), DIMENSION(3*nd)  :: work
      REAL(KIND=8)                   :: eij, eji, c
      REAL(KIND=8), DIMENSION(ns,ks) :: dxb ! the band form of dx
      REAL(KIND=8), DIMENSION(ns,ks) :: hxb 
      
      INTEGER :: mm, k, ms, info, i, j, nns, nnd, ioff, joff, iwork

      nns=nwf*ns
      nnd=nns+nwf
      aa=0; res=0; ipiv=0

      !  Create the diagonal block matrices in hfm
      Do i = 1,nwf
        call hf_matrix(i,hfm(:,:,i))
      end do
      
      ! before modifying the diagonal blocks, generate off-diagonal ones
      Call hm_off

      ! now modify the diagonal blocks
      Call hm_diag

      ! Add normalization constraints
      Call hm_norm

      ! Add orthogonality constraints
      Call hm_orthog
 
       ! solve the Newton Raphson hm equations
     Call hm_solve

      ! update and test the solutions
      Call hm_update

   CONTAINS

!==================================================================
      SUBROUTINE hm_off
!==================================================================
!   generate the off-diagonal blocks for Newton Raphson method
!
!   Rules: (for a given pair, both varied)
!     
! Integral  F(i)       H(ii)      Hnr(ii)         Hnr(i,j)
!  F(i,i)   2R(.i;ii)  2R(.i,.i)  2R(.i,.i)
!                                +4R(...ii)
!  
!  F(i,j)   R(.j;ij)   R(.j;.j)    R(.j,.j)(d)    2R(..;ij)(x)
!  
!  G(i,j)   R(.j;ji)   R(.j;j.)    R(.j;j.)(x)     R(..;ji)(x)
!                                                 +R(.j;.i)(d)
!  
!           -eiiBPi     -eiiB      -eiiB
!           -eijBPj                               -eijB
!  
! In this implementation, Hnr(ii) stored in H(ii). Routine needs to 
! be modified for the case where orbital j is not varied.
! Note that R(..;ji) = transpose(R(..;ij) but it is not clear that
! computing the transpose is faster or as accurate.
!------------------------------------------------------------------
!

      mm = nnd
      ! Create off-diagonal blocks directly before modifying hfm
      Do i = 1,nwf
       ioff = (i-1)*ns
       do j= i+1,nwf
         joff = (j-1)*ns
         hx = 0

         ! direct
         modified = .false.
         do k = 0, 2*min(l(i),l(j)), 2
           c = 2*a(i,j,k)
           if (c /= 0.d0) then
              call density(ns,ks,dx,p(:,i),p(:,j),'x',ms)
              dsumx = c*dx
              call add_rkm_x(k,dsumx, hx)
              modified = .true.
            end if
         end do

         ! exchange
         do k = 0,kmax
            if (abs(l(i)-l(j)) <= k .and. k <= l(i)+l(j)) then
                c= b(i,j,k)
                if ( c /= 0.d0) then
                   modified = .true.
                   call density(ns,ns,dx,p(:,j),p(:,i),'x',ms)
                   dsumx = c*dx
                   call add_rkm_x(k,dsumx,hx)
                    call density(ns,ks,dxb,p(:,i),p(:,j),'s',ms)
                    dxb = c*dxb
                    hxb=0
                    call add_rkm_d(k,dxb,hxb)
                    call add_band_to_full(hxb,hx)
                 end if
             end if
          end do
         if (abs(e(i,j)) > 1.D-12 ) then
            ! We have an orthogonality condition
             mm=mm+1
             eij = dot_product(p(:,j),matmul(hfm(:,:,i),p(:,i)))
             eji = dot_product(p(:,i),matmul(hfm(:,:,j),p(:,j)))
             if (dabs(eij-eji) > 1.d-10) then
                eij = (eij+eji)/2  ! use the average
                hx = hx -eij*bb
                ! set the residual
                call bxv(ks,ns,sb,p(:,j),y)
                y(1) = 0.d0
                y(ns-1:ns) = 0.d0
                res(ioff+1:ioff+ns) = +eij*y
                call bxv(ks,ns,sb,p(:,i),y)
                y(1) = 0.d0
                y(ns-1:ns) = 0.d0
                res(joff+1:joff+ns) = +eij*y
                modified=.true.
             else
                eij = 0.d0
             end if
          end if
          if (modified) then
             call apply_bc(hx,'o')
             aa(ioff+1:ioff+ns, joff+1:joff+ns) = hx
             aa(joff+1:joff+ns, ioff+1:ioff+ns) = transpose(hx)
          end if
        end do
      end do
      if (mm /= nd) then
         Print *, 'Something wrong: mm,nd', mm,nd
      end if

      END SUBROUTINE hm_off
       
    !==================================================================
      SUBROUTINE hm_diag
    !==================================================================
    !    Determine the Lagrange multipliers, subtract from diagonal
    !    block. Then determine the residual (-res), and finally 
    !    modify the diagonal block for Newton Raphson
    !    NOTE:  hf_nr uses res and x = x - Deltax
    !------------------------------------------------------------------
    !
       !Now modify the diagonal blocks and move to aa

       Do i = 1,nwf
       ioff = (i-1)*ns
       x = matmul(hfm(:,:,i),p(:,i))
       ek(i)= dot_product(x,p(:,i))

       ! compute the residuals     
       hfm(:,:,i) = hfm(:,:,i) -ek(i)*bb
       call apply_bc(hfm(:,:,i),'n')
       res(ioff+1:ioff+ns)=res(ioff+1:ioff+ns)-matmul(hfm(:,:,i),p(:,i))
       end do
       
       ! add the Hii corrections to the diagonal blocks
       !   (after res has been computed)
       Do i = 1,nwf
       ioff = (i-1)*ns
       Do k = 0, 2*l(i), 2
         c = 2*a(i,i,k)
         if (c /= 0.d0) then
            call density(ns,ks,dx,p(:,i),p(:,i),'x',ms)
            dsumx = c*dx
            call add_rkm_x(k,dsumx,hfm(:,:,i))
         end if
       end do
       call apply_bc(hfm(:,:,i),'n')
   
       aa(ioff+1:ioff+ns,ioff+1:ioff+ns) = hfm(:,:,i)
       End Do

      END SUBROUTINE hm_diag

    !==================================================================
      SUBROUTINE hm_norm
    !==================================================================
    !  Add normalization constraint to assure the change is orthogonal
    !  to current estimate 
    !------------------------------------------------------------------

      ! Add normalization constraint for orbital i
      mm = nns
      do i = 1,nwf
        ioff = (i-1)*ns
        mm = mm+1
        call bxv(ks,ns,sb,p(:,i),y)
        y(1) = 0.d0
        y(ns)= 0.d0
        aa(ioff+1:ioff+ns,mm) = -y
        aa(mm,ioff+1:ioff+ns) =-y
        res(mm) = 0.d0   ! since all current orbitals are normalized
      end do
       
      END SUBROUTINE hm_norm

    !==================================================================
       SUBROUTINE hm_orthog
    !==================================================================
    !   Add orthogonality constraint to assure the change is orthogonal
    !   to a given orbital
    !------------------------------------------------------------------
       
       mm = nnd
       !Add orthogonality constraints
       !  If (l(i) = l(j))
       do i = 1,nwf
        ioff = (i-1)*ns
        do j = i+1,nwf
          joff = (j-1)*ns
          !if (l(i) == l(j)) then
          if (abs(e(i,j))> 1.d-12 ) then
              mm = mm+1
              call bxv(ks,ns,sb,p(:,j),y)
              y(1) = 0.d0
              y(ns-1) = 0.d0
              y(ns)= 0.d0
              aa(ioff+1:ioff+ns,mm) = -y
              aa(mm,ioff+1:ioff+ns) = -y
              call bxv(ks,ns,sb,p(:,i),y)
              y(1) = 0.d0
              y(ns-1) = 0.d0
              y(ns)= 0.d0
              aa(joff+1:joff+ns,mm) = -y
              aa(mm,joff+1:joff+ns) = -y
              res(mm) =  sum(y*p(:,j))
           end if
         end do
       end do

       If ( mm /= nd) then
         Print *, "mm /= nd"
         STOP "Something wrong"
       end if
      END SUBROUTINE hm_orthog

    !==================================================================
       SUBROUTINE hm_solve
    !==================================================================
    !   Solve and test the solution of the Newton Raphson equations
    !------------------------------------------------------------------

       iwork = 3*nd
       CALL DSYTRF('L',nd, aa, nd, ipiv, work, iwork, info)
       If (info /= 0) then
        Print *, ' Error in factorization: DSYTRF in hf_nr'
        STOP
       End if
       CALL DSYTRS('L',nd, 1, aa, nd, ipiv, res, nd, info)
       If (info /= 0) then
        Print *, ' Error in solve routine : DSYTRS in hfall_nr'
        STOP
       End if

       END SUBROUTINE hm_solve    

    !==================================================================
       SUBROUTINE hm_update
    !==================================================================
    !  Solve and test the solution of the Newton Raphson equations
    !------------------------------------------------------------------

       ! save the current set of orbitals
       do i = 1,nwf
        ioff = (i-1)*ns
        v = p(:,i)+res(ioff+1:ioff+ns)
        do  j = ns,ns/2,-1
          if (abs(v(j)) > end_tol) exit
        end do
        maxr(i) = min(ns,j+1)
        v(j+1:ns) = 0.d0
        e(i,i) = (ek(i) +res(nns+i))/qsum(i)
        dpm(i) = maxval( abs(p(:,i)-v(:)))/maxval(abs(p(:,i)))
        p(:,i) = v
        call normalize(p(:,i))
       end do
       END SUBROUTINE hm_update

     END SUBROUTINE hfall_nr 

