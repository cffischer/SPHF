!====================================================================
    MODULE hf_param
!====================================================================
!
!   Parameters for solving the HF equations
!     Tests for convergence:
!     . scf_tol :  relative change in the total energy
!     . orb_tol :  maximum change in an orbital relative to its 
!                  maximum value
!     . end_tol :  values in the tail all less than this value in
!                  magnitude are set to zero.
!     Other parameters:
!     . acc_max :  Levels of accuracy 
!
!--------------------------------------------------------------------
    IMPLICIT NONE
   
    CHARACTER(LEN=4) :: varied ='all'
    CHARACTER(LEN=80):: vstring =' '
    REAL(KIND=8)     :: scf_tol=1.d-12, &
                        orb_tol=1.d-06, &
                        end_tol=1.d-06
    INTEGER          :: acc, acc_max=2

    END MODULE hf_param
