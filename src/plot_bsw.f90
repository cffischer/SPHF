!======================================================================
    Subroutine plot_bsw(nu)
!======================================================================
! .. Computes and tabulates the radial functions p(:,i), i=1,nwf at all 
! .. gausian points, plus border  values, for further plots
!
!----------------------------------------------------------------------

    Use spline_param
    Use spline_grid
    Use hf_atomic_state 
    Use hf_orbitals

    IMPLICIT NONE
    Integer, Intent(in) :: nu
    REAL(8) :: y(nv*ks+2,nwf),r(nv*ks+2)
    Integer :: i,j,iv,ith,m,io,mr
    Character(4) :: pel(nwf)
    Character(3) :: el3

    y = 0.d0
    Do io=1,nwf

    do iv = 1,nv                 ! over intervals
     do m = 1,ks                 ! over gausian points
      do ith = 1,ks              ! over B-splines in given interval
        i = iv+ith-1             ! B-spline index
        j = 1 + (iv-1)*ks + m    ! radial point index
        y(j,io) = y(j,io) + p(i,io)*bsp(iv,m,ith)
        end do
      end do
    end do

    m = nv*ks+2 ! last point
    Do i=1,ks
     y(m,io) = y(m,io) + bsp(nv+1,1,i) * p(ns-ks+i,io)
    End do

    m = 1       ! first point
    Do i=1,ks
     y(m,io) = y(m,io) + bsp(nv+1,2,i) * p(i,io)
    End do

    End do ! over orbutals (io)

! .. find last point with neglegible values

    mr=nv*ks+2
	Do i=mr,1,-1; if( SUM(abs(y(i,:))) > 1.d-5) Exit;  End do
    mr = i

! .. prepare columns names (not to begin with integer - important for ORIGIN)

    Do io=1,nwf
	 el3=el(io)
     if(el3(1:1).eq.' ') then
	  pel(io)='p'//el3(2:3)
     else
	  pel(io)='p'//el3
     end if
    End do

! .. prepare radial array

    m=1; R(1)=t(1)
    Do i=1,nv; Do j=1,ks; m=m+1; R(m) = gr(i,j); End do; End do
    m=m+1; R(m) = t(ns+1)

! .. finally, recording ...

    rewind(nu)
    write(nu,'(20(10x,a4))') 'r',pel
    Do i=1,mr
     write(nu,'(20(E14.5))') r(i),y(i,:)
	End do

    End Subroutine plot_bsw


