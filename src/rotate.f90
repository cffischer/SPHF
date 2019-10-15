!======================================================================
  SUBROUTINE rotate(i,j)
!======================================================================
!   Determine the rotation parameter for a  a stationary solution for
!   the (j,i) orbitals, i <= j
!   Rules derived from transformation of densities
!  
!  Let P*i = (Pi - ePj)/(1+e^2)2
!      P*j = (Pj + ePi)/(1+e^2)2
!  
! Integral    coeff of e           coeff of e^2  
!  F(i,i)   -4 R(iiij)             2 F(i,j) +4 G(i,j) -2 F(i,i)
!  
!  F(j,j)    4 R(ij;jj)            2 F(i,j) +4 F(i,j) - 2 F(j,j)
!  
!  F(i,j)    2[R(iiij)-R(ijjj)     F(i,i) + F(jj) -2 F(i,j) -4G(i,j)
!  
!  G(i,j)    2[R(iiij)-R(ijjj)     F(i,i) + F(j,j)-2 F(i,j)
!  
!  F(i,m):   -2 R(im,jm)           F(j,m) - F(i,m)
!  
!  F(j,m):    2 R(im,jm)           F(i,m) - F(j,m)
!  
!  G(i,m)    -2R(ij,mm)            G(j,m) - G(i,m)
!  
!  G(j,m)     2R(ij,mm)            G(i,m) - G(j,m)
!  
!----------------------------------------------------------------------
!
    Use spline_param
    USE spline_galerkin, ONLY: sb
    Use hf_atomic_state
    USE hf_orbitals
    USE hf_param

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: i,j
    INTEGER :: m,k
    LOGICAL :: comp_dg
    REAL (KIND=8), EXTERNAL :: bhl, a, b, fk, gk, rk
    REAL (KIND=8), DIMENSION(ns) :: Pi
    REAL (KIND=8):: eps,dn,g,dg,Hii,Hij,Hjj, Fkii,Fkij,Fkjj,Gkij

    g= 0.d0; dg = 0.d0; comp_dg = .true.
    if (dpm(i) <= 1.d-5 .and. dpm(j) <= 1.d-5 .and. abs(rdg(i,j))>1.d-5) then
       ! Previous is sufficiently accurate
       comp_dg = .false.
       dg = rdg(i,j)
    end if

    ! Contributions from ii, jj, and ij
    Hij = -0.5*bhl(i,j)
    G = 2*(qsum(j) - qsum(i))*Hij
    If (comp_dg) then
      Hii = -0.5*bhl(i,i)
      Hij = -0.5*bhl(i,j)
      Hjj = -0.5*bhl(j,j)
      dg = (qsum(i)-qsum(j))*(Hjj-Hii)
    End if
    !Print '(10X,A,T20, 3F16.9)', 'HL', Hii, Hij, Hjj
    Do k= 0,2*min0(l(i),l(j)),2
      ! Note coefficient is divided by 2 if i=j
      G = G + 2*(a(i,j,k)+b(i,j,k) - a(i,i,k))*rk(i,i,i,j,k)
      G = G - 2*(a(i,j,k)+b(i,j,k) - a(j,j,k))*rk(i,j,j,j,k)
      If (comp_dg) then
        FKii = fk(i,i,k)
        FKij = fk(i,j,k)
        FKjj = fk(j,j,k)
        !Print '(10X,A, T20, 4F16.9)', 'Fk', FKii, FKij, FKjj
        !Print '(10X,A, T20, 3F16.9)', 'A(i,j,k)', a(i,i,k)/2, a(i,j,k), a(j,j,k)/2
        GKij = gk(i,j,k)
        !Print '(10X,A, T20, 4F16.9)', 'Gk', GKii, GKij, GKjj
        !Print '(10X,A, T20, F16.9)', 'B(i,j,k)', b(i,j,k)
        dG = dG + a(i,i,k)*(FKij-FKii +2*GKij) +a(j,j,k)*(FKij-FKjj +2*GKij)
        dG = dG + (a(i,j,k)+b(i,j,k))*(FKii + Fkjj -2*FKij -4*GKij)
        !Print '(10X,A,T20,2F16.9)', 'Rk(iiij)', rk(i,i,i,j,k), rk(i,j,j,j,k)
        ! Print '(10X,(A, T20, 2F16.9))', 'G,dG', g, dg
      End if
    End do
    ! .. Contributions from im, and jm
    do m = 1,nwf
         !Print *, 'Rotate: Loop on m', i,j,m
       if (m ==i .or. m == j) Cycle
       do  k = 0,2*min0(l(i),l(j)),2
          if (a(i,m,k) == a(j,m,k)) Cycle
          g  = g + 2*(a(j,m,k)-a(i,m,k))*rk(i,m,j,m,k)
          if (comp_dg) &
          dg =dg + (a(j,m,k)-a(i,m,k))*(Fk(i,m,k)-Fk(j,m,k))
       end do
       do k = abs(l(i)-l(m)),l(i)+l(m),2
          if (b(i,m,k) == b(j,m,k)) Cycle
          g  = g + 2*(b(j,m,k)-b(i,m,k))*rk(i,m,m,j,k)
          if (comp_dg) &
            dg =dg + (b(j,m,k)-b(i,m,k))*(Gk(i,m,k)-Gk(j,m,k))
       end do
     end do
     !Print *, 'rotate: g,dg', g,dg
     if ( abs(g) <1.d-12 .and. abs(dg) < 1.d-12 ) then
        eps = 0.d0
        e(i,j) = 0.d0
        e(j,i) = 0.d0
     else
        eps = -g/(2*dg)
        if (comp_dg) rdg(i,j) = dg
        dn = 1.d0/sqrt(1+eps*eps)
        Pi = (p(:,i)-eps*p(:,j))*dn
        p(:,j) = (p(:,j)+eps*p(:,i))*dn
        p(:,i) = Pi
      end if
   end SUBROUTINE rotate
       
           
   
 
