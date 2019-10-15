!
    SUBROUTINE  read_bsw
    !==================================================================
    !  Read bsw.inp for initial estimates
    !  No scaling law is applied at the moment nor expansions
    !      for a different grid
    !------------------------------------------------------------------
    Use spline_param
    Use hf_atomic_state
    USE hf_orbitals
    USE hf_param
    USE hf_inout

         IMPLICIT NONE
         INTEGER      :: ksw, nsw, mw, ios, iel, i
         REAL(KIND=8) :: zw,hw,hmaxw,rmaxw
         CHARACTER(LEN=3) :: elw
         Logical      :: skip, connected, keep

       INQUIRE(FILE='bsw.inp', OPENED=connected)
       If (.not. connected) RETURN
       Do
         READ(UNIT=iuf, IOSTAT=ios) elw,zw,hw,hmaxw,rmaxw,ksw,nsw,mw
         If (ios <0 .or. zw <1.d0 )  EXIT  ! end of file
         if (ios >0 ) STOP 'Read_bsw: Error in reading bsw.inp'
         skip = .false.
         keep = .false.
         iel =0
         elw = adjustr(elw)
         Do i=1,nwf
           if (elw == el(i)) then
              iel=i; exit
           end if
         End do
         if (iel /= 0) then
           if (in(i) == 0) keep = .true.
         end if
         if (iel==0  .or. (.not. keep ) ) then
            READ(IUF, IOSTAT=ios)      ! orbital not in the list
            CYCLE
         end if
         if (zw /= z .or. hw /= h .or. ksw /= ks )  then
            skip = .true.
            write(err, FMT='(/7X,A,A3,1X,A)') 'Read_bsw: Grid for ', elw
            write(err,FMT='(12X,A,F6.1,F8.4,2I6)') 'File:',zw,hw,ksw,nsw
            write(err,FMT='(12X,A,F6.1,F8.4,2I6)') 'Curr:',z,h,ks,ns
         end if
         if (in(i) == -1) cycle    ! We already have estimates
         Call el_nl(elw,n(iel),l(iel))
         if (skip) then
           call different_grid
         else
           If (mw > ns) then
              write(err,FMT='(A,A,A/A)') 'Range for ',el(i),' is too large:', &
                                  ' will truncate the tail'
              mw = ns
           End if
           READ(iuf) p(1:mw,iel)
           If (mw < ns) p(mw+1:ns,iel) = 0.d0
         end if
         maxr(iel) = mw
         in(iel) = -1
      End do
    CONTAINS
    !==================================================================
    SUBROUTINE  different_grid
    !==================================================================
    !  Obtain estimates for the current new grid
    !------------------------------------------------------------------

       USE spline_grid
       USE spline_galerkin, ONLY: sb
       REAL(KIND=8), DIMENSION(nsw+ksw) :: tw
       REAL(KIND=8), DIMENSION(nsw)     :: pw
       REAL(KIND=8), DIMENSION(nv,ks)   :: yg, ygw
       REAL(KIND=8), DIMENSION(ks,ns)   :: bsf
       REAL(KIND=8), DIMENSION(ns)      :: v
       REAL(KIND=8), EXTERNAL           :: bvalu2
       INTEGER :: i, j, ierr, mlw
       REAL(KIND=8) :: rmaxw

       READ(iuf) pw(1:mw)
       pw(mw+1:nsw) = 0.d0

      mlw = NINT(1.d0/hw)

      tw(1:ksw) = 0.d0
      Do i = ksw+1, ksw+mlw
        tw(i) = tw(i-1) +hw
      end do
      Do i = ksw+mlw+1,nsw+1
        tw(i) = tw(i-1)*(1+hw)
      End do
      tw(nsw+2:nsw+ksw) = tw(nsw+1)
      !tw = tw/zw
      tw = tw/z
      rmaxw = tw(nsw+1)

       Do i=1,nv; do j=1,ks
           if (gr(i,j) > rmaxw) then
             yg(i,j) = 0.d0
           else
             yg(i,j) = bvalu2(tw,pw(:),nsw,ksw,gr(i,j),0)
           end if
           ygw(i,j)= yg(i,j)*grw(i,j)
        End do; End do

!   ... form the vector of inner products of the radial function 
!   ... and the new spline basis functions:

         Call VINTY (ygw,v)

!   ... solve the system of equations  Sb x = cb:  

         !Print '(10X,A,A/(5D15.5))', 'RHS expansion for',el(iel),v
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
         if (zw /= z) call normalize(v)
         !   ... find tail of function and remove too small values
         Do i=ns,1,-1
          if(abs(v(i)) > end_tol) Exit
         End do
         v(i+1:ns) = 0.d0
         p(:,iel) = v
        ! maxr(iel) = i
         mw = i
      END SUBROUTINE different_grid

      END SUBROUTINE read_bsw


