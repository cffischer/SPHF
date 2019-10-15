!====================================================================
    MODULE hf_inout
!====================================================================
!
!   contains input/output parameters
!
!--------------------------------------------------------------------

    IMPLICIT NONE
    SAVE

    INTEGER(4) :: err  =  0 !   error (system default for screen)
    INTEGER(4) :: log  =  3 !   hf_log file
    INTEGER(4) :: scr  =  6 !   screen output
    INTEGER(4) :: iuc  = 20 !   input unit for bsw.c (if present)
    INTEGER(4) :: iuf  = 21 !   input unit for bsw.inp (if present)
    INTEGER(4) :: iup  = 22 !   input unit for hf_param (if present)
    INTEGER(4) :: ouf  = 31 !   output unit for bsw.out
    INTEGER(4) :: plt  = 41 !   plot.dat file (for plotting)
    INTEGER(4) :: ryd  = 42 !   file for Rydberg series orbitals

    CONTAINS
!====================================================================
    SUBROUTINE open_files
!====================================================================
!   Open files to be used
!--------------------------------------------------------------------
    IMPLICIT NONE
    Logical :: old

    OPEN(UNIT=log, FILE='sphf.log', STATUS='UNKNOWN', POSITION='asis')
    INQUIRE(FILE='bsw.c', EXIST=OLD) 
    IF (OLD) &
      OPEN(UNIT=IUC, FILE='bsw.c', STATUS='OLD', FORM='FORMATTED') 
    INQUIRE(FILE='bsw.inp', EXIST=OLD) 
    IF (OLD) &
      OPEN(UNIT=IUF, FILE='bsw.inp', STATUS='OLD', FORM='UNFORMATTED') 
    INQUIRE(FILE='hf_param', EXIST=OLD) 
    IF (OLD)  OPEN(UNIT=iup, FILE='hf_param', STATUS='OLD', &
                             FORM='FORMATTED') 
    OPEN(UNIT=OUF, FILE='bsw.out', STATUS='UNKNOWN', FORM='UNFORMATTED') 
    OPEN(UNIT=plt, FILE='plot.dat', STATUS='UNKNOWN')
    END subroutine open_files

!====================================================================
    END MODULE hf_inout
