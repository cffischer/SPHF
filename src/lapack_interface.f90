!======================================================================
      MODULE lapack_interface
!======================================================================

      Implicit none

        INTERFACE

         SUBROUTINE DSYTRF( UPLO, N, A, LDA, IPIV, WORK, LWORK, INFO )
           CHARACTER(LEN=1), INTENT(IN) :: UPLO
           INTEGER, INTENT(IN) :: LDA, LWORK, N
           INTEGER, INTENT(OUT) :: INFO, IPIV(*)
           REAL(KIND=8), INTENT(INOUT) :: A( LDA,*)
           REAL(KIND=8), INTENT(OUT) :: WORK( LWORK )
         END SUBROUTINE DSYTRF

         SUBROUTINE DSYTRS( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, INFO )
           CHARACTER(LEN=1), INTENT(IN) :: UPLO
           INTEGER, INTENT(IN) :: LDA, LDB, N, NRHS
           INTEGER, INTENT(OUT) :: INFO
           INTEGER , INTENT(IN) :: IPIV(*)
           REAL(KIND=8), INTENT(IN) :: A( LDA,*)
           REAL(KIND=8), INTENT(INOUT) :: B( LDB,*)
         END SUBROUTINE DSYTRS

         SUBROUTINE DGETRF( M, N, A, LDA, PIV, INFO )
           INTEGER, INTENT(IN) :: LDA, M, N
           INTEGER, INTENT(OUT) :: INFO
           INTEGER, INTENT( OUT ) :: PIV( * )
           REAL(KIND=8), INTENT( INOUT ) :: A( LDA, * )
         END SUBROUTINE DGETRF

         SUBROUTINE DGETRS( TRANS, N, NRHS, A, LDA, PIV, B, LDB, INFO )
           CHARACTER(LEN=1), INTENT(IN) :: TRANS
           INTEGER, INTENT(IN) :: LDA, LDB, NRHS, N
           INTEGER, INTENT(OUT) :: INFO
           INTEGER, INTENT(IN) :: PIV(*)
           REAL(KIND=8), INTENT(IN) :: A(LDA,*)
           REAL(KIND=8), INTENT(INOUT) :: B(LDB,*)
         END SUBROUTINE DGETRS

        END INTERFACE

      End MODULE lapack_interface