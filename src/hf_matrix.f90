!=====================================================================
      subroutine hf_matrix(i,hfm) 
!=====================================================================
!
!     Set up the hf_matrix for orbital i.  Banded matrices will be
!     summed in the array hd, converted to full hfm, and exchange
!     contributions added to hfm.
!
!-----------------------------------------------------------------------
      USE spline_param
      USE hf_atomic_state
      USE hf_orbitals
      USE spline_hl
      USE hf_energy_expression,  ONLY: kmax
      
      IMPLICIT NONE
      REAL(KIND=8), DIMENSION(ns,ns),INTENT(out) :: hfm
      INTEGER(4) :: i

      REAL(KIND=8), DIMENSION(ns,ks) :: hd,d,dsum
      REAL(KIND=8), DIMENSION(ns,ns) :: dx,dxsum
      REAL(KIND=8), DIMENSION(ns) :: y
      REAL(KIND=8) :: sumi,c
      REAL(KIND=8), EXTERNAL :: a,b
      INTEGER :: j, ll, k, nw, ms
      LOGICAL :: found

     !print *, 'Entering hf_matrix for orbital', i
     hfm = 0.d0

     !... copy H0(l(i)) to h1
     ll = l(i) 
     sumi = qsum(i) 
     Call hlm(ll)
     hd = -sumi*hl/2

     !.. add contributions from direct terms to dsum
     !.. first the F0(i,j)
     dsum = 0.d0
     if (in(i) == 0) then
        nw = i-1
     else
        nw = i
     end if
     ! For NOW we want to include all contributions.
     nw =nwf
     k = 0
     do j = 1,nw   !For all orbitals
       c = a(i,j,0)
       !Print *, 'a(i,j,0)', i, j, a(i,j,0)
       if ( c ==0.d0) cycle
       call density(ns,ks,d,p(:,j),p(:,j),'s',ms)
       dsum = dsum+ c*d
     end do
     
     ! Add hf_direct to hd  for k=0
     call add_rkm_d(k,dsum,hd)

     ! Add the direct contribution from Fk(i,i)
     If (qsum(i) > 1 .and. in(i) /= 0) then
       Do k = 2, 2*l(i), 2
          call density(ns,ks,d,p(:,i),p(:,i),'s', ms)
          dsum = a(i,i,k)*d
          call add_rkm_d(k,dsum,hd)
       end do
     End if
     
     Call add_band_to_full(hd,hfm)
     y = matmul(hfm,p(:,i))

     ! Add exchange contribution to x
     If (in(i) /= 0) then
       Do k = 0,kmax
          dxsum = 0.d0
          found = .false.
          Do j = 1,nw
            if (abs(l(i)-l(j)) <= k .and. k <= l(i)+l(j)) then
  	     c= b(i,j,k)
             !Print *, 'b(i,j,0)', i, j,k, b(i,j,k)
  	     if ( c /= 0.d0) then
  	       found = .true.
                 call density(ns,ns,dx,p(:,j),p(:,j),'x',ms)
                 dxsum = dxsum +c*dx
               end if
            end if
          end do
          if (found) call add_rkm_x(k,dxsum,hfm)
       end do
     y = matmul(hfm,p(:,i))
     End if
            
   END SUBROUTINE hf_matrix
