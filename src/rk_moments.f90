!======================================================================
    SUBROUTINE rk_moments(k)
!======================================================================
!
!   Defines moments for Rk-integrals in the B-spline cells
!
!   Calling sequence:          rk_moments             
!                              ----------             
!                               /    \\           
!                           moments rk_pdiag      
!                                     ||          
!                                   rk_triang     
!                                    /   \        
!                                 gauss  vbsplvb  
! 
!----------------------------------------------------------------------
!
!   on entry    k  -  multipole index
!   --------
!       
!   on exit     rkd1,rkd2,rkd - off-diagonal and diagonal moments 
!   -------                     (in module spline_moments)
!
!----------------------------------------------------------------------

    USE spline_param
    USE spline_grid
    USE spline_moments
    USE hf_atomic_state

    IMPLICIT NONE
    INTEGER(4), INTENT(in) :: k

    ! .. check the need of calculations

    if(mtype == 'aaa') Call allocate_moments
    if(mtype == 'rk ' .and. kmk == k) Return

    ! .. compute the moments in the spline basis

    CALL moments(  k   , rkd1,'s','b')
    CALL moments(-(k+1), rkd2,'s','b')
    CALL rk_pdiag

    ! .. add the relativistic correction

    if(rel) Call rk_rel      

    mtype='rk '
    kmk=k

    CONTAINS


!======================================================================
    SUBROUTINE rk_pdiag
!======================================================================
!
!   Controls the scaling propeties for diagonal B-spline Rk-interals
!
!   Calls:  rk_triang
!
!----------------------------------------------------------------------

    IMPLICIT NONE

    ! .. local variables

    INTEGER(4) :: iv,ik
    REAL(8) :: hp

    ! .. non-exponential grid

    if(me.eq.0) then
     DO iv=1,nv
      CALL rk_triang(iv)
     END DO
     Return
    end if

    ! .. the first equal step region

    DO iv=1,ml+ks-1
      CALL rk_triang(iv)
    END DO

    ! .. the exponential region - using scaling law

    hp=h+1.d0
    ik = ks*(ks+1)/2
    DO iv=ml+ks,ml+me-ks+2
      rkd(1:ik,1:ik,iv) =  rkd(1:ik,1:ik,iv-1) * hp
    END DO

    ! .. the last equal step region

    DO iv=ml+me-ks+3,nv
      CALL rk_triang(iv)
    END DO

    END SUBROUTINE rk_pdiag


!======================================================================
    SUBROUTINE rk_triang(iv)
!======================================================================
!
!   Returns the two-dimensional array of B-spline integrals 
!               <B_i B_j|r^k/r^(k+1)|B_i' B_j'>
!   over the given triangle diagonal cell 
!
!   Calls:   gauss, vbsplvd
!
!   On entry   iv  -  index of the diagonal cell
!   --------
!
!   On exit    rkd(.,.,iv) - arrays of Rk B-spline integrals for given 
!   --------                 interval iv in the reduced-dimension mode
!                            (in module spline_moments)
!----------------------------------------------------------------------

    IMPLICIT NONE

    INTEGER(4), INTENT(IN) :: iv

    ! .. local variables

    INTEGER(4) :: i,j, ip,jp, ii,jj, m, left,  ik
    REAL(8) :: xbase
    REAL(8), DIMENSION(ks) :: x,w, bi, gx,gw
    REAL(8), DIMENSION(ks,ks) :: bspTmp
    REAL(8), DIMENSION(ks,ks,ks) :: INT
    REAL(8), DIMENSION(nv,ks,ks) :: dbiatx
    REAL(8), DIMENSION(ks*(ks+1)/2,ks*(ks+1)/2) :: a

    left=iv+ks-1
    xbase=t(left)

! .. setup the gaussian points

    CALL gauss(ks,x,w)

    DO m=1,ks

! .. the absolute coordinate at the new gaussian point

      gx(:) = (gr(iv,m)-xbase)*x(:) + xbase

! .. the bspline values at the new gaussian points

      DO i=1,ks
       CALL vbsplvd(t,left,1,gx(i),1,dbiatx)
       bspTmp(i,1:ks)= dbiatx(1,1:ks,1)
      END DO

! .. and the corresponding gaussian weights

      gw(:) = (gr(iv,m)-xbase)*w(:)

      IF(k>1) THEN
        gx(:) = gw(:)*gx(:)**k
      ELSE IF(k==1) THEN
        gx(:) = gw(:)*gx(:)
      ELSE IF(k==0) THEN
        gx(:) = gw(:)
      END IF

!            / r(iv,m)                             k
! .. INT =  |      bsp(iv,:,j)(r) bsp(iv,:,jp)(r) r  dr
!           / r_iv

      DO j=1,ks
       gw(:) = gx(:)*bspTmp(:,j)
       DO jp=j,ks
        INT(j,jp,m)= SUM(gw(:)*bspTmp(:,jp))
       END DO
      END DO
    
    END DO	 !  over m

! .. second integration 

    IF(k/=0) THEN
      gx(:) = grw(iv,:)*grm(iv,:)**(k+1)
    ELSE
      gx(:) = grw(iv,:)*grm(iv,:)
    END IF

    ii = 0
    DO i=1,ks
     DO ip=i,ks
      ii = ii+1
      bi(:) = bsp(iv,:,i)*bsp(iv,:,ip)*gx(:)
      jj = 0
      DO j=1,ks
       DO jp=j,ks
        jj = jj + 1
        a(ii,jj) =  SUM(bi(:)*INT(j,jp,:))
       END DO
      END DO
     END DO
    END DO
    
    ik = ks*(ks+1)/2
    rkd(1:ik,1:ik,iv) = a + TRANSPOSE(a)

    END SUBROUTINE rk_triang


!======================================================================
    SUBROUTINE rk_rel 
!======================================================================
!
!   relativistic corrections to the Rk integrals
!
!   (they have another scaling multiplier, then RK !)
!
!----------------------------------------------------------------------

    Implicit none

    REAl(8) :: C
    REAl(8), DIMENSION(ks) :: a,b
    Integer(4) :: iv, i,j, ip,jp, ii,jj 

     C = fine*(k+k+1)
     Do iv = 1,nv
     a(:) = grm(iv,:)*grm(iv,:)*grw(iv,:) * C
     ii = 0
     DO i=1,ks
      DO ip=i,ks
       ii = ii+1
       b(:) = bsp(iv,:,i)*bsp(iv,:,ip)*a(:)
       jj = 0
       DO j=1,ks
        DO jp=j,ks
         jj = jj + 1
         rkd(ii,jj,iv) = rkd(ii,jj,iv) +  &
                         SUM(b(:)*bsp(iv,:,j)*bsp(iv,:,jp))
        END DO
       END DO
      END DO
     END DO
     END DO

    End SUBROUTINE rk_rel 



    END SUBROUTINE rk_moments

