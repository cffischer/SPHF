!=====================================================================
   subroutine solve_HF_equations
!=====================================================================
!
!     Solve the HF equations 
!------------------------------------------------------------------
     USE spline_param, ONLY: h,hmax,rmax,ns,ks,ml
     USE spline_galerkin, ONLY: sb,bb
     USE hf_atomic_state
     USE hf_energy_expression, ONLY: kmax
     USE hf_inout
     USE hf_orbitals
     USE hf_param
     
     IMPLICIT NONE

     REAL(KIND=8), DIMENSION(ns,ns) :: hfm
     REAL(KIND=8), DIMENSION(ns) :: v

     REAL(KIND=8), EXTERNAL :: a,b, bhl, bvmv, fk, gk, azl
     INTEGER :: i,j,m, md, it, north
     INTEGER, DIMENSION(nwf) :: jorth

     REAL(KIND=8) ::  et, orb_diff, scf_diff
     INTEGER :: ipr, lp, ip

     Call allocate_rka(0,kmax)  !! needed before call to energy
     call energy 
     write(scr,'(/7X,A,F20.12)') 'Total energy from initial estimates =',&
                                   etotal
     ipr = 0; et = etotal
     Do it = 1,12
       write(6,'(//7X,A,I6/7X,A/)') 'Iteration ',it,'----------------'
       WRITE (6,'(11X,A,11X,A,15X,A,10X,A,8X,A)') &
          'nl', 'E(nl)','AZ', 'DPM', 'MAXR'
       Do ip = nwf-nit+1, nwf
         i = iord(ip)
         m = n(i)-l(i)
         north = 0
         !.. set orthogonality constraints by forming the list jorth
         do j = 1,nwf
           if (i == j .or. l(i) /= l(j)) cycle
           if (clsd(i) .and. clsd(j)) cycle  ! both shells are closed
   	   if (it == 1 .and. (in(i)==1 .or. in(j)==1) .and. j >=i &
                       .and. acc ==1) then
                cycle ! orthogonal only to j < i for screened estimates 
           else
             north = north+1
             jorth(north) = j
             if (j < i) m = m-1   ! index of eiv reduces only for j<i
           end if
           if (( it > 1 .or. ( in(i)==-1 .and. in(j)==-1)) .and. &
                 j > nwf-nit .and.  i < j .and. e(i,j) >1.d-10) then
             !print *, 'calling rotate'
             call rotate(i,j) 
           end if
         end do

         !print *, 'calling hf_matrix for orbital', el(i)
         call hf_matrix(i,hfm)
         if (it == 1 .or. qsum(i) < 2.d0 .or. dpm(i) > 0.5 ) then
           !Print *, 'Calling hf_eiv with north =', north
	    call hf_eiv
         Else 
           md = ns+1+north
           v = p(:,i)                         ! needed for hf_nr
           !Print *, 'Calling hf_nr with north =', north
           call hf_nr(i,north, jorth, hfm, bb, md, v, e(i,i))
         end if
         dpm(i) = maxval( abs(p(:,i)-v(:)))/maxval(abs(p(:,i)))

         do j = ns, ns/3,-1
           if (abs(v(j)) > end_tol ) exit
         end do

         maxr(i)=min(ns,j+1)
         v(j+1:ns) = 0.d0
         p(:,i) =v(:)             
         lp = l(i)+1
         az(i) = p(lp+1,i)*azl(z,h,ks,lp) 
         WRITE(scr,'(10X,A3,F20.12,F18.10,1P,D12.2,I8)') &
                 el(i), e(i,i), az(i), dpm(i), maxr(i)
       end do
       !call orthonormalize(nwf)
       call energy
       ! test convergence
       orb_diff =  maxval(abs(dpm(1:nwf)))
       scf_diff =  abs(et-etotal)/abs(etotal)
        i = maxr(nwf)-1
       print '(/(10X,A,T50,1P,D10.2,D12.2))', &
        'SCF convergence (diff vs tol)', scf_diff, scf_tol, &
        'Orbital convergence (diff vs tol)', orb_diff, orb_tol, &
        'Tail cut-off (orbital nwf)', p(i,nwf), end_tol
       If ( orb_diff < orb_tol .or. scf_diff  < scf_tol) then
         exit
       else
         et = etotal
       End if
       write(scr,'(7X,A,F20.12)') 'Total energy =', etotal
     End do
     call summry(orb_diff, scf_diff, p(i,nwf))
     if (acc == acc_max) call write_bsw

   CONTAINS

   !==================================================================
       SUBROUTINE hf_eiv
   !==================================================================
   !    Find the eigenvector of hfm for the m'th eigenvalue after 
   !    orthogonality has been applied
   !------------------------------------------------------------------

     IMPLICIT NONE

     REAL(KIND=8), DIMENSION(ns,ns) :: aa, ss
     REAL(KIND=8), DIMENSION(3*ns) :: W
     REAL(KIND=8), DIMENSION(ns) :: eigval
     INTEGER :: j, jp, INFO
   
     ! Apply orthogonality conditions for orbitals
     Do jp = 1,north
        j = jorth(jp)
       call apply_orthogonality(hfm,p(:,j))
     End do

     !.. apply boundary conditions for matrices of size mm
     !   This way hfm will not be destroyed by the eigensolver
     aa = hfm(1:ns,1:ns)
     ss = bb
     CALL apply_bc (aa,'e')
     CALL apply_bc (ss,'n')

     ! .. evaluates the eigenvalues and eigenvectors
     call dsygv(1,'V','L',ns,aa,ns,ss,ns,eigval,W,3*ns,INFO)
     ! .. the eigenvectors are now in aa
     if (INFO /= 0) then
       WRITE(UNIT=err, FMT='(A,I6)') 'Error in Eigenvalue routine', info
       STOP
     end if
     v(1:ns) = aa(1:ns,1+m)
     if (v(ml) < 0.d0) v = -v
     e(i,i) = eigval(1+m)/qsum(i)
     END SUBROUTINE  hf_eiv
   !==================================================================
     SUBROUTINE hf_nr(i,north, jorth, hfm, bb, md, v, eval)
   !==================================================================
   !    Improve the estimate v by the Newton_Raphson method
   !    subject to orthogonality.
   !   Rules: (for a given pair, both varied)
   !     
   ! Integral  F(i)       H(ii)      Hnr(ii)       
   !  F(i,i)   2R(.i;ii)  2R(.i,.i)  2R(.i,.i)
   !                                +4R(...ii)
   !------------------------------------------------------------------
     USE lapack_interface
       
     INTEGER, INTENT(in) :: i, north, md
     INTEGER, DIMENSION(:), INTENT(in) :: jorth
     REAL(KIND=8), DIMENSION(ns,ns), INTENT(in) :: hfm, bb
     REAL(KIND=8), DIMENSION(ns), INTENT(inout) :: v
     REAL(KIND=8), INTENT(out) :: eval

     REAL(KIND=8),  DIMENSION(md,md) :: aa, aaa
     REAL(KIND=8),  DIMENSION(md) :: res, rhs
     INTEGER(4),  DIMENSION(md) :: ipiv

     REAL(KIND=8), EXTERNAL  :: a
     REAL(KIND=8), DIMENSION(ns) :: x, y
     REAL(KIND=8), DIMENSION(ns,ns) :: dx, dsumx, hx
     REAL(KIND=8) :: eij
       
     INTEGER :: mm, k,  info, j, jp, ms
     aa = 0.d0
     x = matmul(hfm,v)
     eval= dot_product(x,v)

     ! compute the residuals     
     aa(1:ns,1:ns) = (hfm -eval*bb)
     call apply_bc(aa(1:ns,1:ns),'n')
     res(1:ns) =  matmul(aa(1:ns, 1:ns),v)
     
     mm = ns+1
     ! Add normalization constraint for orbital i
     call bxv(ks,ns,sb,v,y)
     y(1) = 0.d0
     y(ns-1) = 0.d0
     y(ns)= 0.d0
     aa(1:ns,mm) = -y
     aa(mm,1:ns) =-y
     res(mm) = 0.d0
     !Add orthogonality constraints with orbital j
     do jp = 1,north
         mm = mm+1
         j = jorth(jp)
         eij = dot_product(p(:,j),x(:))
         call bxv(ks,ns,sb,p(:,j),y)
         y(1) = 0.d0
         y(ns-1) = 0.d0
         y(ns)= 0.d0
         res(mm) = 0.d0
         aa(1:ns,mm) = -y
         aa(mm,1:ns) = -y
         res(1:ns) = res(1:ns) +eij*y
     end do

     If ( mm /= md) then
        Print *, "mm /= md"
        STOP
     end if
  
     ! Add Haa contributions
       
     hx = 0.d0
     Do k = 0, 2*l(i), 2
        call density(ns,ks,dx,p(:,i),p(:,i),'x',ms)
        dsumx = 2*a(i,i,k)*dx
        call add_rkm_x(k,dsumx,hx)
     end do

     call apply_bc(hx,'n')
     aa(1:ns,1:ns) = aa(1:ns,1:ns) + hx
   
     ! save the equations for testing
     aaa = aa
     rhs = res
     CALL DGETRF(md, md, aa, md, ipiv, info)
     If (info /= 0) then
       WRITE(UNIT=err,FMT='(A)') 'Error in factorization: DSYTRF in hf_nr'
       STOP
     End if
     CALL DGETRS('N', md, 1, aa, md, ipiv, res, md, info)
     If (info /= 0) then
       WRITE(UNIT=err,FMT='(A)') 'Error in solve routine: DSYTRS in hf_nr'
        STOP
     End if
     v = v-res(1:ns)
     call normalize(v)
     eval = (eval -res(ns+1))/qsum(i)
    END SUBROUTINE  hf_nr

   !==================================================================
    SUBROUTINE  write_bsw
   !==================================================================
   !    Ouput the orbitals
   !------------------------------------------------------------------
     Do i = 1,nwf
       Write(ouf) el(i), Z, h, hmax, rmax,ks,ns,maxr(i),e(i,i),dpm(i)
       Write(ouf) p(:,i)
     End do
    End Subroutine write_bsw

   END subroutine solve_HF_equations
