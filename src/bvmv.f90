!=====================================================================
    FUNCTION bvmv(n,k,array,type,x,y)
!=====================================================================
!
!   Returns bvmv = <x,array y> for the band array.
!
!   The original size of array is n*n. The bandwidth of array is 2k-1
!   for the non-symmetric case, and is k for symmetric case.
!   The lower-band column storage mode is assumed.
!---------------------------------------------------------------------
!
!   on entry
!   --------
!       n      the leading dimension.
!       k      the band width.
!       array  a banded matrix in the lower-band column storage mode.
!              Its element (i,j) represents the element (i,i+j-k) of
!              the original matrix. For a symmetric matrix, only the
!              lower diagonals are needed.
!       type   'g' for general;
!              's' for symmetric;
!              'a' for "almost"  symmetric, as db2
!       x      vector with length n.
!       y      vector with length n.
!
!----------------------------------------------------------------------

    REAL(KIND=8) :: bvmv
    INTEGER, INTENT(in) :: n, k
    CHARACTER(LEN=1), INTENT(in) :: type
    REAL(KIND=8), DIMENSION(n,*), INTENT(in) :: array
    REAL(KIND=8), DIMENSION(n), INTENT(in) :: x, y

    ! .. Local variables

    INTEGER :: i,j,jp,kp

    bvmv = 0.d0

    if ( type /= 'g') then

      ! .. symmetric

      do jp = 1,k-1
        do i = k+1-jp,n
          j = i+jp-k
          bvmv = bvmv + array(i,jp) * (x(i)*y(j) + x(j)*y(i))
        end do
      end do

      ! .. make correction for 'a' (almost symmetric)

      if (type == 'a')  &
          bvmv = bvmv + x(n-1)*y(n)*(array(1,1)-array(n,k-1))

    else

      do jp = 1,k-1
        do i = k-jp+1,n
          j = i+jp-k
          kp = 2*k-j
          bvmv = bvmv + array(i,jp)*x(i)*y(j)
          bvmv = bvmv + array(j,kp)*x(j)*y(i)
        end do
      end do

    end if

    ! .. add central diagonal

    bvmv = bvmv + SUM(array(:,k)*x(:)*y(:))

  END FUNCTION bvmv

