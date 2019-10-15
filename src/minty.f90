!===================================================================
    SUBROUTINE minty(icase,ygr,ym)
!===================================================================
!
!   Computes the array elements   <B_i, y(r) B_j>.
!
!-------------------------------------------------------------------
!
!   on entry
!   --------
!       icase  0 if i >=j (subdiagonals); 1 if i<=j (superdiagonals)
!       ygr    array of values of a specific function  y(r) at the
!              gaussian points of each interval, weighted by the
!              gaussian weights.
!
!   on exit
!   -------
!       ym     <B_i, y(r) B_j> in symmetric lower-band (icase=0) or
!                              upper-band (icase=1) storage mode
!
!--------------------------------------------------------------------

    USE spline_param
    USE spline_grid

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: icase
    REAL(KIND=8), DIMENSION(nv,ks), INTENT(IN) :: ygr
    REAL(KIND=8), DIMENSION(ns,ks), INTENT(OUT) :: ym

    ! .. local variables

    INTEGER :: ith, jth, i, irow, jcol


    ym = 0.d0      ! .. clear the ym array

    if (icase /= 0) then

      do ith = 1,ks
        do jth = ith,ks
          jcol = jth-ith+1
          do i = 1,nv
            irow = i+ith-1
            ym(irow,jcol) = ym(irow,jcol)+   &
                  SUM( ygr(i,:)*bsp(i,:,ith)*bsp(i,:,jth) )
          end do
        end do
      end do

    else

      do ith = 1,ks
        do jth = 1,ith
          jcol = jth-ith+ks
          do i = 1,nv
            irow = i+ith-1
            ym(irow,jcol) = ym(irow,jcol)+   &
                SUM( ygr(i,:)*bsp(i,:,ith)*bsp(i,:,jth) )
          end do
        end do
      end do

    end if

  END SUBROUTINE minty
