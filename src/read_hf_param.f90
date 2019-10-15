!====================================================================
    SUBROUTINE read_hf_param
!====================================================================
!   Read hf_param if it is connected (exists)
!--------------------------------------------------------------------

    USE hf_inout
    USE hf_param

    IMPLICIT NONE
    INTEGER       :: ios
    CHARACTER(LEN=80) :: line, line1, line2
    CHARACTER(LEN=7)  :: var
    CHARACTER(LEN=20) :: value
    Logical           :: connected

    INQUIRE(FILE='hf_param', OPENED=connected)
    !Print *, 'Is hf_param connected?',connected
    If (.not. connected) return

    Read(iup,'(A)')  line1
    line1 = adjustl(line1)
    If (line1(1:2) /= 'HF') &
       Stop 'First line of hf_param does not start with HF'
    Read(iup,'(A)')  line2
    
    Do 
      call get_next
      if (ios /= 0) EXIT
      select case (trim(var))
       case ('scf_tol'); Read(value,*) scf_tol
       case ('orb_tol'); Read(value,*) orb_tol
       case ('end_tol'); Read(value,*) end_tol 
       case ('acc_max'); Read(value,*) acc_max
       case ('varied'); Read(value,*) varied
       case ('list')  ; READ(value, FMT='(A)') vstring
      end select
    END DO
    !Print *, 'read hf_param:', scf_tol, orb_tol, end_tol, acc_max, varied

   CONTAINS
!======================================================================
   SUBROUTINE get_next
!======================================================================
   IMPLICIT NONE

   INTEGER ::  i

   Read(iup,'(A)',IOSTAT=ios) line
   IF (ios /= 0) Return
   i = INDEX(line,'=')  
   var =trim(adjustl( line(1:i-1)))
   value = adjustl(line(i+1:))
   
  END subroutine get_next

  END subroutine read_hf_param
    
!====================================================================
    SUBROUTINE write_hf_param
!====================================================================
!   Read hf_param if it is connected
!--------------------------------------------------------------------
    USE spline_param
    USE hf_inout
    USE hf_param
    IMPLICIT NONE
    Logical           :: connected


    INQUIRE(FILE='hf_param', OPENED=connected)
    If (.not. connected) then
      OPEN(UNIT=iup, FILE='hf_param', STATUS='UNKNOWN', &
                     FORM='FORMATTED', ACTION='READWRITE')
    Else
        rewind(unit=iup)
    End if
    WRITE(UNIT=iup,FMT='(A)') 'HF_parameters'
    WRITE(UNIT=iup,FMT='(A)') '-------------'
    WRITE(UNIT=iup,FMT='(A,F12.6)') 'h       = ',h
    WRITE(UNIT=iup,FMT='(A,I12)') 'ks      = ',ks
    WRITE(UNIT=iup,FMT='(A,I12)') 'ns      = ',ns
    WRITE(UNIT=iup,FMT='(A,1PD12.2)') 'scf_tol = ',scf_tol
    WRITE(UNIT=iup,FMT='(A,1PD12.2)') 'orb_tol = ',orb_tol
    WRITE(UNIT=iup,FMT='(A,1PD12.2)') 'end_tol = ',end_tol
    WRITE(UNIT=iup,FMT='(A,I12)') 'acc_max = ',  1
    WRITE(UNIT=iup,FMT='(A,8X,A4)') 'varied  = ', adjustr(varied)
    If (varied == 'list') &
      WRITE(UNIT=iup,FMT='(A,A)')   'list    = ', trim(vstring)
    
    
    END subroutine write_hf_param 
!====================================================================
