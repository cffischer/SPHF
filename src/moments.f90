!=====================================================================
    SUBROUTINE  moments (k,rkm,sym,dir)
!=====================================================================
!
!   Computes moments defining as <B_i|r^k|B_j> over an interval
!
!---------------------------------------------------------------------
!
!   on entry
!   --------
!      k      the power of moments
!      sym    's' or 'n' - symmetrical or non-symetrical case
!      div    'b' or 'd' - B_j or B'_j 
!
!   on exit
!   -------
!                                              
!      rkm    array of moments for every interval 'iv'
!
!---------------------------------------------------------------------
  
      USE spline_param
  
      IMPLICIT NONE
      INTEGER(4), INTENT(in) :: k
      REAL(8), DIMENSION(ks*ks,nv), INTENT(out) :: rkm
      CHARACTER(1), INTENT(in) :: sym, dir
  
      ! .. local variables
  
      INTEGER :: iv, jk
      REAL(KIND=8) :: hp
      REAL(KIND=8), DIMENSION(ks*ks) :: rv
  
      jk = ks*(ks+1)/2
      if(sym.eq.'n') jk=ks*ks
  
      if(me.eq.0) then     ! .. there is no exponential grid
  
       Do iv = 1,nv
        Call moment(iv,k,rv,sym,dir)
        rkm(1:jk,iv) = rv(1:jk)
       End do 
  
      else                 ! .. case of exponential grid

       ! .. the first non-exponential region

       DO iv = 1,ml+ks-1
        Call moment(iv,k,rv,sym,dir)
        rkm(1:jk,iv) = rv(1:jk)
       END DO

       ! .. the exponential region --- using scaling law

       hp = (1.d0+h)**k
       if(dir.eq.'b') hp = hp*(1.d0+h)
       DO iv=ml+ks,ml+me-ks+2
        rkm(1:jk,iv) = rkm(1:jk,iv-1)*hp
       END DO

       ! .. the last non-exponential region

       DO iv=ml+me-ks+3,nv
        Call moment(iv,k,rv,sym,dir)
        rkm(1:jk,iv) = rv(1:jk)
       END DO

      end if

      END SUBROUTINE moments



!=====================================================================
    SUBROUTINE  moment (iv,k,rv,sym,dir)
!=====================================================================
!
!   Computes moment defining as <B_i|r^k|B_j> for given interval
!
!---------------------------------------------------------------------
!
!   on entry
!   --------
!      k      the power of moments
!      sym    's' or 'n' - symmetrical or non-symetrical case
!      div    'b' or 'd' - B_j or B'_j 
!      iv     index of interval
!
!   on exit
!   -------
!                                              
!      rv     array of moments for given interval 'iv'
!
!---------------------------------------------------------------------

      USE spline_param
      USE spline_grid
      
      IMPLICIT NONE
      INTEGER(4), INTENT(in) :: k,iv
      REAL(8), DIMENSION(ks*ks), INTENT(out) :: rv
      Character(1), INTENT(in) :: sym, dir
      
      ! .. local variables
      
      INTEGER(4) :: i, j, ii, jj
      REAL(8), DIMENSION(ks) :: gw
      REAL(8), DIMENSION(ks,ks) :: bi,bj

      gw = grw(iv,:)
      if( k > 0 ) gw = grw(iv,:) * gr (iv,:)**(+k)
      if( k < 0 ) gw = grw(iv,:) * grm(iv,:)**(-k)

      bi(:,:) = bsp(iv,:,:)
      Do j = 1,ks
       if(dir.eq.'b') then
        bj(:,j) = bsp(iv,:,j)*gw(:)
       else           
        bj(:,j) = bsq(iv,:,j)*gw(:)
       end if
      End do

      ii = 0
      DO i=1,ks
       jj = 1; if(sym.eq.'s') jj = i
       DO j=jj,ks
        ii = ii + 1
        rv(ii) = SUM(bi(:,i)*bj(:,j))
       END DO
      END DO

      END SUBROUTINE moment
