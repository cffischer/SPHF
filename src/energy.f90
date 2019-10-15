!======================================================================
  SUBROUTINE energy
!======================================================================
!   Determine the total energy of the state and the Virial Theorem                                                                 
!
!----------------------------------------------------------------------
!
    Use spline_param
    Use hf_atomic_state
    USE hf_orbitals
    USE hf_param

    IMPLICIT NONE
    INTEGER :: i,j,k
    REAL(KIND=8), DIMENSION(nwf) :: hli
    REAL (KIND=8), EXTERNAL :: bhl, a, b, fk, gk, rk

    REAL (KIND=8):: c

    DO i = 1,nwf
      hli(i) = -0.5D0*bhl(i,i)
      !print '(I5,2F16.9)', i, qsum(i), hli(i)
    END DO
    
    Etotal = 0.d0
    do i = 1,nwf
        etotal = etotal + qsum(i)*hli(i)
        do j = 1,i
           do  k = 0,2*min0(l(i),l(j)),2
               c = a(i,j,k)
               if (i .eq. j) c = c/2.d0
               if (abs(c).ne.0.d0) then
                  etotal = etotal + c*fk(i,j,k)
                  !print '(3I5,F10.6,F16.9,2X,A)', i,j,k,c,fk(i,j,k), 'Fk'
               end if
            end do
         end do 
         do  j = 1,i-1
            do k = abs(l(i)-l(j)),l(i)+l(j),2
               c = b(i,j,k)
               if (abs(c).ne.0.d0) etotal=etotal+c*gk(i,j,k)
               !print '(3I5,F10.6,F16.9,2X,A)', i,j,k,c,gk(i,j,k), 'Gk'
            end do
         end do
     end do
     end SUBROUTINE energy
    
