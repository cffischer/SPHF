!====================================================================
   MODULE spline_galerkin
!====================================================================
!
!  contains common arrays used in the application of splines and 
!  the Galerkin method 
!
!--------------------------------------------------------------------

    IMPLICIT NONE
    SAVE

    REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:):: r1, rm1, rm2, rm3
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:):: sb, bs, bb, db1, db2

! ...  sb (1:ns,1:ks)  -  sb (i,j) -->  <B_i|B_i+j-k>
! ...  r1 (1:ns,1:ks)  -  r1 (i,j) -->  <B_i|r|B_i+j-k>
! ...  rm1(1:ns,1:ks)  -  rm1(i,j) -->  <B_i|1/r|B_i+j-k>
! ...  rm2(1:ns,1:ks)  -  rm2(i,j) -->  <B_i|1/r^2|B_i+j-k>
! ...  rm3(1:ns,1:ks)  -  rm2(i,j) -->  <B_i|1/r^3|B_i+j-k>
! ...  db1(1:ns,1:ks)  -  db1(i,j) -->  <B_i|B'_i+j-k>
! ...  db2(1:ns,1:ks)  -  db2(i,j) -->  <B_i|B"_i+j-k>
! ...  bs(1:ks,1:ns)   -  factorization of sb (upper band form)
!      bb(1:ns,1:ns)   -  sb in full-matrix form

! ...  all arrays (except bs) in the symmetric lower-columb storage
! ...  mode, db2 - 'almost' symmetruc, db1 - antisymmetric


    CONTAINS


!====================================================================
    SUBROUTINE allocate_galerkin
!====================================================================
!
!   allocates space for the arrays in the MODULE spline_galerkin
!
!--------------------------------------------------------------------

    USE spline_param

    if(Allocated(r1)) then
      Deallocate(r1,rm1,rm2,rm3,sb,bs,bb,db1,db2)
    end if
    ALLOCATE( r1(ns,ks), rm1(ns,ks), rm2(ns,ks), rm3(ns,ks), &
              sb(ns,ks), bb(ns,ns), bs(ks,ns), db1(ns,ks), db2(ns,ks))
    r1 = 0.d0; rm1 = 0.d0; rm2 = 0.d0; sb = 0.d0; bs = 0.d0; bb=0.d0
	db1 = 0.d0; db2 = 0.d0


    END SUBROUTINE allocate_galerkin

    END MODULE spline_galerkin
