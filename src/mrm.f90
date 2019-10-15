!=======================================================================
    SUBROUTINE mrm(mm,rm)
!=======================================================================
!
!   Computes the matrix representation of the operator
!
!             r^mm
!
!   in the B-spline basis.
!
!-----------------------------------------------------------------------
!
!   on entry
!   --------
!       mm      integer defining the power of r
!
!   on exit
!   -------
!       rm      <B_i,r^mm B_j> in symmetric storage mode
!
!-----------------------------------------------------------------------

    USE spline_param
    USE spline_grid

    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(ns,ks), INTENT(OUT) :: rm
    INTEGER, INTENT(IN) :: mm

    ! .. local variables

    INTEGER :: i, irow, jcol, ith, jth

    ! .. clear the rm array

    rm = 0.d0

    ! .. assemble the matrix elements

    if (mm == 0) then
      do  ith=1,ks
        do jth=1,ith
          jcol=jth-ith+ks
          do i=1,nv
            irow=i+ith-1
            rm(irow,jcol)=rm(irow,jcol) &
              + SUM(grw(i,:)*bsp(i,:,ith)*bsp(i,:,jth))
          end do
        end do
      end do
    else if (mm == -1) then
      do ith=1,ks
        do jth=1,ith
          jcol=jth-ith+ks
          do i=1,nv
            irow=i+ith-1
            rm(irow,jcol)=rm(irow,jcol) &
              + SUM(grw(i,:)*bsp(i,:,ith)*grm(i,:)*bsp(i,:,jth))
          end do
        end do
      end do
    else if (mm <= -2) then
      do ith=1,ks
        do jth=1,ith
          jcol=jth-ith+ks
          do i=1,nv
            irow=i+ith-1
            rm(irow,jcol)=rm(irow,jcol) &
              + SUM(grw(i,:)*bsp(i,:,ith)*grm(i,:)**(-mm)*bsp(i,:,jth))
          end do
        end do
      end do
    else if (mm == 1) then
      do ith=1,ks
        do jth=1,ith
          jcol=jth-ith+ks
          do i=1,nv
            irow=i+ith-1
            rm(irow,jcol)=rm(irow,jcol) &
             + SUM(grw(i,:)*bsp(i,:,ith)*gr(i,:)*bsp(i,:,jth))
          end do
        end do
      end do
    else if (mm >= 2) then
      do ith=1,ks
        do jth=1,ith
          jcol=jth-ith+ks
          do i=1,nv
            irow=i+ith-1
            rm(irow,jcol)=rm(irow,jcol) &
             + SUM(grw(i,:)*bsp(i,:,ith)*gr(i,:)**mm*bsp(i,:,jth))
          end do
        end do
      end do
    end if

  END SUBROUTINE mrm
