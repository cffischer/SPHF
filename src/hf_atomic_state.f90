!=======================================================================
  MODULE hf_atomic_state 
!=======================================================================
!   This module defines the parameters for the problem to be solved
!----------------------------------------------------------------------
!
    IMPLICIT NONE
    SAVE

    ! atom variables
    REAL(KIND=8) :: z, Etotal
    REAL(KIND=8) :: EC = 0.d0   ! core energy
    REAL(KIND=8) :: fine = 0.25D0/(137.036D0)**2

    LOGICAL :: rel = .FALSE.    ! relativistic corrections

    CHARACTER(LEN=2) :: atom
    CHARACTER(LEN=3) :: term
    CHARACTER(LEN=32) :: configuration
 

    ! orbital variables that depend on nwf 
    INTEGER :: nclosd, nwf, nit
    CHARACTER(LEN=3), DIMENSION(:), ALLOCATABLE :: el
    INTEGER, DIMENSION(:), ALLOCATABLE :: n, l, maxr, iord, in
    LOGICAL, DIMENSION(:), ALLOCATABLE :: clsd, vary
    INTEGER :: lmax
    REAL(KIND=8), DIMENSION(:), ALLOCATABLE ::  dpm, s,  qsum, az
    REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: e, rdg

 END module hf_atomic_state
  
    !===================================================================
      SUBROUTINE allocate_atomic_state
    !===================================================================
    !   This program allocates arrays associated with the orbitals (nwf)
    !-------------------------------------------------------------------
       USE hf_atomic_state
        
        IMPLICIT NONE
	ALLOCATE( el(nwf) )
	ALLOCATE( n(nwf), l(nwf), maxr(nwf), iord(nwf), in(nwf), &
                 clsd(nwf), vary(nwf))
	ALLOCATE( dpm(nwf), s(nwf), qsum(nwf), az(nwf) )
        ALLOCATE( e(nwf,nwf), rdg(nwf,nwf) )
        e=0.d0;; rdg=0.d0
     
      END SUBROUTINE allocate_atomic_state

