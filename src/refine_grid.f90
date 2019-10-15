!  Routines for refining the grid
!  . refine_grid
!  . reset_int    
!==================================================================
    SUBROUTINE  refine_grid
    !==================================================================
    ! Improve accuracy: h => h/2; ks => ks+2 ; ns =>2*min(maxr(:)+2, ns)
    ! Tolerances decreased by a factor of 10
    ! Change form method 'd' to 'c'
    !------------------------------------------------------------------
    Use spline_param
    USE spline_grid
    USE spline_galerkin
    Use hf_atomic_state
    USE hf_orbitals
    USE hf_inout
    USE hf_param

      IMPLICIT NONE
      REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: t2,v
      REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: p2, yg, ygw, bsf
      REAL(KIND=8), EXTERNAL :: bvalu2
      INTEGER :: ns2,ks2,nt2,i,j,m, ierr
      REAL(KIND=8) :: h2, hmax2, rmax2, ss,rd, sm

      ! save current values
      ns2 = ns; ks2 = ks; nt2=ns2+ks2

      allocate(p2(ns2,nwf),t2(nt2))
      p2 = p; t2=t
      h2 = h; hmax2=hmax; rmax2=rmax
      
      !define new spline parameters
      h = h/2
      ks = ks + 2
      !print '(A,20I6)', 'maxr:', maxr(1:nwf)
      !print *, 'maxval(maxr(:))', maxval(maxr(:))
      ns = 2*min((maxval(maxr(:))+1),ns)   
      !Print *, 'New value of ns:', ns
      call mkgrid(z) 
      scf_tol = scf_tol/1000
      orb_tol = orb_tol/1000
      end_tol = end_tol/1000
      
      ! .. for the current nwf and new ns
      CALL allocate_orbital_arrays
      Call define_spline
      
      ! transform the orbitals to the new basis
      Allocate(yg(nv,ks),ygw(nv,ks), v(ns), bsf(ks,ns))
      
      write(UNIT=scr,FMT='(/7X,A)') 'Refined grid orbitsl'
      do m = 1,nwf
        !Print *, 'Orbital ', el(m)
        Do i=1,nv; do j=1,ks
           if (gr(i,j) > rmax2) then
             yg(i,j) = 0.d0
           else
             yg(i,j) = bvalu2(t2,p2(:,m),ns2,ks2,gr(i,j),0)
           end if
           ygw(i,j)= yg(i,j)*grw(i,j)
        End do; End do

!   ... form the vector of inner products of the radial function 
!   ... and the new spline basis functions:

         Call VINTY (ygw,v)

!   ... solve the system of equations  Sb x = cb:  
   
         !Print '(10X,A,A/(5D14.5))', 'RHS expansion for',el(m),v
         v(1) = 0.d0
         Do i = ns-1,ns/2,-1
            if (abs(v(i)) > 0.d0) exit
         End do
         bsf=transpose(sb)
	! Apply bc at large r
         i = i+1
         bsf(1:ks-1,i:ns) = 0.d0
         bsf(ks,i:ns)     = 1.d0
         ! ... apply bc at the origin
         do j = 1,ks-1
           bsf(j,ks-j+1)=0.d0
         end do
         bsf(ks,1) = 1.d0
         Call dpbtrf('U', ns,ks-1,bsf,ks,ierr)
         if (ierr.ne.0)  then
         write(UNIT=err,FMT='(10X,A/I6)') &
             'refine_grid: dpbtrs (LAPACK) error No: ', ierr
         end if
         Call dpbtrs ('U',ns,ks-1,1,bsf,ks,v,ns,ierr)  
         if(ierr.ne.0) then
         write(UNIT=err,FMT='(10X,A/I6)') &
             'refine_grid: dpbtrs (LAPACK) error No: ', ierr
         end if

!   ... find tail of function and remove too small values
                                                
         Do i=ns,1,-1
          if(abs(v(i)) > 1.d-10) Exit
         End do
         v(i+1:ns) = 0.d0
         p(:,m) = v
         maxr(m) = i

         sm = 0.d0; rd = 0
         Do i=1,nv; Do j=1,ks
           ss = bvalu2(t,p(:,m),ns,ks,gr(i,j),0) 
           if(abs(ss-yg(i,j)).lt.sm) Cycle
           sm = abs(ss-yg(i,j))
           rd = gr(i,j)
         End do; End do

         write(UNIT=scr,FMT='(10X,a4,a,1PD10.2,a,F8.2)') &
           EL(m),'   max. diff.',SM,'  at r =',rd
         dpm(m) = sm
     End do
     Call reset_int

    deallocate(p2,t2,bsf,yg,ygw,v)
    END SUBROUTINE refine_grid

!=======================================================================
  SUBROUTINE reset_int
!=======================================================================
!   Reset the integral parameters
!----------------------------------------------------------------------
!
    USE spline_hl
    USE spline_integrals
    USE spline_moments
   
    lh    = -1   
    itype = 'aaa'
    krk   = -100
    mtype = 'aaa'
    kmk   = -100

  END SUBROUTINE reset_int
