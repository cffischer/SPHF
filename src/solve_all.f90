!=====================================================================
      subroutine solve_HF_equations
!=====================================================================
!
!  Solve the HF equations 
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
      REAL(KIND=8), DIMENSION(ns)    ::  v
      REAL(KIND=8), DIMENSION(ns,ns) :: aa, ss
      REAL(KIND=8), DIMENSION(ns) :: eigval

      ! REAL(KIND=8) :: eti,eval, ans
      REAL(KIND=8), EXTERNAL :: a,b, bhl, bvmv, fk, gk
      ! INTEGER :: ios, eof, i,j,m, mm, md, ll, k, it, nw, north, nd
      INTEGER :: i,j,m,md, it, north, nd
      INTEGER, DIMENSION(nwf) :: jorth

      REAL(KIND=8) ::  et, orb_diff, scf_diff, azl
      INTEGER :: ipr, lp, ip, it_max
      LOGICAL :: converged

      Call Allocate_rka(0,kmax)


      Print '(/7X,A,I2/7X,A)', 'Level ',acc,'-----'
      call energy   ! provides initial estimate of the total energy
      write(scr,'(/7X,A,F20.12)') &
                    'Total energy from initial estimates =', etotal
      ipr = 0; et = etotal
      converged = .false.

      Do it = 1,5
         write(6,'(/7X,A,I2.2/7X,A/)') 'Iteration 1.',it,'----------------'
         WRITE (6,'(11X,A,11X,A,15X,A,10X,A,8X,A)') &
            'nl', 'E(nl)','AZ', 'DPM', 'MAXR'
         Do ip = nwf-nit+1, nwf
            i = iord(ip)
            m = n(i)-l(i)
            north = 0
            !.. set orthogonality constraints by forming the list jorth
            do j = 1,nwf
              if (i == j .or. l(i) /= l(j)) cycle
              if (clsd(i) .and. clsd(j)) cycle  ! if both shells are closed
	      if (it == 1 .and. (in(i)==1 .or. in(j)==1) .and. j >=i &
                          .and. acc ==1) then
                 cycle    ! orthogonal only to j < i
              else
                 north = north+1
                 jorth(north) = j
                 if (j < i) m = m-1         ! index of eiv reduces only for j<i
              end if
              if (( it > 1 .or. ( in(i)==-1 .and. in(j)==-1)) .and. &
                    j > nwf-nit .and.  i < j .and. e(i,j) >1.d-12) then
                call rotate(i,j) 
              end if
            end do
            call hf_matrix(i,hfm)
            if (it == 1 .or. qsum(i) < 2.d0 .or. dpm(i) > 0.5 ) then
	      call hf_eiv
            Else 
              md = ns+1+north
              v = p(:,i)                         ! needed for hf_nr
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
         ! test convergence
         call energy
         orb_diff =  maxval(abs(dpm(1:nwf)))
         scf_diff =  abs(et-etotal)/abs(etotal)
         do i= maxr(nwf),ns/2,-1
           if (abs(p(i,nwf))  > end_tol) exit
         end do
         print '(/(10X,A,T50,1P,D10.2,D12.2))', &
           'SCF convergence (diff vs tol)', scf_diff, scf_tol, &
           'Orbital convergence (diff vs tol)', orb_diff, orb_tol, &
           'Tail cut-off (orbital nwf)', p(i,nwf), end_tol
         write(scr,'(7X,A,F20.12)') 'Total energy =', etotal
         If ( orb_diff < orb_tol .or. scf_diff  < 1000*scf_tol) then
            converged = .true.
            exit
         else
             et = etotal
         End if
      End do
  
      If (converged) then
         if (acc == 1 .and. acc_max == 2 ) then
             it_max=0   !  this is the Level one
         else 
             it_max=2
         end if
         Print '(/10X,A,I2,A)', 'The first phase of Level',acc, ' has converged'
         Print '(10X,A,I2)', 'it_max for the second phase set to', it_max
      else
         it_max = 10
         Print '(/10X,A,I2,A)', 'The first phase of Level',acc, ' has not converged'
      end if
      ! hfall_nr
      ! Compute total number of orthogonality conditions
      north = 0
      if (nit < nwf) it_max = 0 !   <<<  do not vary all
      Do i = nwf-nit+1,nwf
        Do j = i+1,nwf
           if ( l(i) == l(j) .and. abs(e(i,j)) >1.d-12) north = north +1
        End do
      End do
      Do it = 1,it_max
         write(6,'(/7X,A,I2.2/7X,A)') 'Iteration 2.',it,'----------------'
         Do ip = nwf-nit+1, nwf
           i = iord(ip)
           Do j = 1,nwf
              if ( j > nwf-nit .and.  i < j .and. abs(e(i,j)) >1.d-12) then
                call rotate(i,j)
              end if
           End do
         End do
         nd = (ns+1)*nwf + north
         call hfall_nr(nd)
         call energy
         ! test convergence
         orb_diff =  maxval(abs(dpm(1:nwf)))
         scf_diff =  abs(et-etotal)/abs(etotal)
         do i= maxr(nwf),ns/2,-1
           if (abs(p(i,nwf))  > end_tol) exit
         end do
         print '((10X,A,T50,1P,D10.2,D12.2))', &
           'SCF convergence (diff vs tol)', scf_diff, scf_tol, &
           'Orbital convergence (diff vs tol)', orb_diff, orb_tol, &
           'Tail cut-off (orbital nwf)', p(i,nwf), end_tol
         write(scr,'(7X,A,F20.12)') 'Total energy =', etotal
         If ( orb_diff < orb_tol .or. scf_diff  < scf_tol) then
            exit
         else
             et = etotal
         End if
      End do
      call summry(orb_diff, scf_diff, p(i,nwf))
      if (acc == acc_max) then
         call write_bsw
         Call plot_bsw(plt)
      end if
      Call dealloc_integrals

  CONTAINS

    !==================================================================
        SUBROUTINE hf_eiv
    !==================================================================
    !    Find the eigenvector of hfm for the m'th eigenvalue after 
    !    orthogonality has been applied
    !------------------------------------------------------------------

        IMPLICIT NONE

        REAL(KIND=8), DIMENSION(3*ns) :: W
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
      !else
       !  WRITE(99,'(3(I3,F22.16))') (j, eigval(j), j = 2,4)
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
   !------------------------------------------------------------------
        USE lapack_interface
        
        INTEGER, INTENT(in) :: i, north, md
        INTEGER, DIMENSION(:), INTENT(in) :: jorth
        REAL(KIND=8), DIMENSION(ns,ns), INTENT(in) :: hfm, bb
        REAL(KIND=8), DIMENSION(ns), INTENT(inout) :: v
        REAL(KIND=8), INTENT(out) :: eval

        REAL(KIND=8),  DIMENSION(md,md) :: aa
        REAL(KIND=8),  DIMENSION(md) :: res
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
          !res(mm) =  sum(y*v)
          res(mm) = 0.d0
          !y = -eij*y
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
          CALL DGETRF(md, md, aa, md, ipiv, info)
         If (info /= 0) then
           WRITE(UNIT=err,FMT='(A)') ' Error in factorization: DSYTRF in hf_nr'
           STOP
         End if
         CALL DGETRS('N', md, 1, aa, md, ipiv, res, md, info)
         If (info /= 0) then
           WRITE(UNIT=err,FMT='(A)') ' Error in solve routine : DSYTRS in hf_nr'
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
      USE spline_galerkin, ONLY: r1
      IMPLICIT NONE
      INTEGER          :: nel, jj
      REAL(KIND=8)     :: en
      CHARACTER(LEN=3) :: ell = '   '
      
      Do i = 1,nwf
        Write(ouf) el(i), Z, h, hmax, rmax,ks,ns,maxr(i), e(i,i), dpm(i)
        Write(ouf) p(:,i)
      End do

      if (qsum(nwf) == 1) then
         Print '(/7X,A)', 'Rydberg orbitals written to bsw.out'
         Print '(/6X,A)', 'Orbital       Energy         <r>       rmax'
         call hf_matrix(nwf, hfm)
         north =0
         m = n(nwf)-l(nwf)
         do j = 1,nwf-1
           if (l(nwf) /= l(j)) cycle
           north = north+1
           jorth(north) = j
           m = m-1         ! index of eiv reduces only for j<i
         end do
         i=nwf
         call hf_eiv
         nel = n(nwf)
         ell(1:3) = el(nwf)
         Do j = m+1, min(99,ns-3)
            en = eigval(1+j)/qsum(nwf)
            v(1:ns) = aa(1:ns, 1+j)
            do jj = ns, ns/3,-1
              if (abs(v(jj)) > end_tol ) exit
            end do
            write( ell(1:2), FMT='(I2)') j+ n(nwf) -m
            Print '(7X,A4,2F15.6, I8)', ell,en,bvmv(ns,ks,r1,'s',v,v),jj
            Write(ouf) ell, Z, h, hmax, rmax,ks,ns,jj, en, 1.0
            Write(ouf) v
         end do
      end if

      End Subroutine write_bsw

   END subroutine solve_HF_equations
          
