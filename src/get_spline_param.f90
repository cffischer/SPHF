!======================================================================
     SUBROUTINE get_spline_param(z)
!======================================================================
!
!    This routine gets information about the spline parameters,
!    how the calculation is to be performed, and making the grid in
!    the Z*r variable.
!
!----------------------------------------------------------------------
     Use spline_param
     Use spline_grid
     Use hf_inout

     IMPLICIT NONE
     REAL(KIND=8), INTENT(IN) :: z
     CHARACTER(LEN=80) :: line, line1, line2
     CHARACTER(LEN=7)  :: var
     CHARACTER(LEN=20) :: value
     INTEGER           :: ios
     LOGICAL           :: connected

     ! Set the default values
     h=0.25; ks=4; ns= 30+ 2*sqrt(z)
  
     ! Read optional user defined values
     INQUIRE(FILE='hf_param', OPENED=connected)
     If (connected) then
       Rewind(UNIT=iup)
        Read(iup,'(A)')  line1
        line1 = adjustl(line1)
        If (line1(1:2) /= 'HF') &
            Stop 'First line of hf_param does not start with HF'
        Read(iup,'(A)')  line2
    
        Do 
         call get_next
         if (ios /= 0) EXIT
         select case (trim(var))
          case ('h'); Read(value,*) h
          case ('ks'); Read(value,*) ks
          case ('ns'); Read(value,*) ns
         end select
        END DO
     END IF

     call mkgrid(z)

     
    CONTAINS
 !======================================================================
    SUBROUTINE get_next
 !======================================================================
    IMPLICIT NONE
 
    INTEGER ::  i
 
    Read(iup,'(A)',IOSTAT=ios) line
    IF (ios /= 0) Return
    i = INDEX(line,'=')  
    var = trim(adjustl(line(1:i-1)))
    value = adjustl(line(i+1:))
    
   END subroutine get_next
   END SUBROUTINE get_spline_param


 !======================================================================
    SUBROUTINE write_spline_param
 !======================================================================
    USE spline_param
    USE spline_grid, ONLY: t
    USE hf_param
    USE hf_inout
    IMPLICIT NONE

     WRITE(UNIT=err,FMT='(/7X,A,I2,A)') 'Level',acc,': HF_parameters '
     WRITE(UNIT=err,FMT='(7X,A)')       '----------------------'
     WRITE(UNIT=err,FMT='(10X, A,T42, F10.5)') 'Step-size (h)',h
     WRITE(UNIT=err,FMT='(10X, A,T42, I10)' )  'Spline order (ks)',ks
     WRITE(UNIT=err,FMT='(10X, A,T42, I10)' )  'Size of basis (ns)',ns
     WRITE(UNIT=err,FMT='(10X, A,T42, F10.2)' ) 'Maximum radius', t(ns+1)
     WRITE(UNIT=err,FMT='(10X,A,T40,1PD12.2)') 'SCF convergence tolerance', scf_tol
     WRITE(UNIT=err,FMT='(10X,A,T40,1PD12.2)') 'Orbital convergence tolerance',orb_tol
     WRITE(UNIT=err,FMT='(10X,A,T40,1PD12.2)') 'Orbital tail cut-off', end_tol
     WRITE(UNIT=err,FMT='(10X,A,T48, A4)') 'Orbitals varied', adjustr(varied)

     WRITE(UNIT=log,FMT='(/7X,A,I2,A)') 'Level',acc,': HF_parameters '
     WRITE(UNIT=log,FMT='(7X,A)')       '----------------------'
     WRITE(UNIT=log,FMT='(10X, A,T42, F10.5)') 'Step-size (h) ',h
     WRITE(UNIT=log,FMT='(10X, A,T42, I10)' )  'Spline order (ks)',ks
     WRITE(UNIT=log,FMT='(10X, A,T42, I10)' )  'Size of basis  (ns)',ns
     WRITE(UNIT=log,FMT='(10X, A,T40, F12.2)' ) 'Maximum radius  ', t(ns+1)
     WRITE(UNIT=log,FMT='(10X,A,T40,1PD12.2)') 'SCF convergence tolerance', scf_tol
     WRITE(UNIT=log,FMT='(10X,A,T40,1PD12.2)') 'Orbital convergence tolerance',orb_tol
     WRITE(UNIT=log,FMT='(10X,A,T40,1PD12.2)') 'Orbital tail cut-off', end_tol
     WRITE(UNIT=log,FMT='(10X,A,T48, A4)') 'Orbitals varied', adjustr(varied)

    END SUBROUTINE write_spline_param
