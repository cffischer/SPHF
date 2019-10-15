!======================================================================
  PROGRAM sphf
!======================================================================
!
!  This program computes the radial functions for simple Hartree-Fock
!  cases.
!
!   SUBROUTINE called:
!       open_files (in hf_inout.f90) 
!       read_hf_param
!       get_case
!       define_spline
!       get_estimates
!       refine_grid
!       write_spline_param  (in get_spline_param.fi90)
!       solve_HF_equations  (in solve_scf.f90 and solve_all.f90
!       write_hf_param ( in read_hf_param.f90)
!
!   Written by:  by Charlotte  Froese Fischer
!   Date:        April, 2010 
!                                                                
!----------------------------------------------------------------------
!
    USE hf_inout
    USE hf_param
    IMPLICIT NONE

    CALL open_files

    ! .. read parameters (optional)
    CALL read_hf_param
    

    ! .. get data about the problem to be solved and how the spline
    ! .. calculation is to be performed, including the grid points
    CALL get_case

    ! .. Initialize the basic B-spline environment and
    ! .. evaluate spline matrix elements for typical operators in
    ! .. in non-relativistic atomic structure theory using the
    ! .. spline-galerkin method
    CALL define_spline

    ! Solve for good and then excellent accuracy
    Do acc = 1,acc_max
      If (acc==1) then
          Call get_estimates
      Else
          Call Refine_grid
      End if
      Call write_spline_param
      Call solve_HF_equations

    End DO

    Call write_hf_param        ! for rerunning the highest level
    
  END PROGRAM sphf
