!======================================================================
      SUBROUTINE get_case
!======================================================================
!
!     This routine obtains information about the problem to be solved
!     and how the spline methods are to be applied

!----------------------------------------------------------------------
      Use hf_atomic_state
      Use hf_orbitals
      Use hf_energy_expression
      
      LOGICAL :: done
      INTEGER :: iselec

      ! .. Header information
      WRITE(6,9)
9     FORMAT(///////22X,'===========================================',&
                   /22X,'S P L I N E  H A R T R E E - F O C K : 2010',&
                   /22X,'==========================================='/)
 
      !WRITE(6,'(//A/A//)') ' START OF CASE',' ============='

      ! .. get the atomic state problem
 get_a: DO
      CALL get_atom

      ! .. determine the energy expression
 get_t: DO
        CALL get_energy(done)
        IF (.NOT.DONE) THEN
           WRITE (0, "(/,' The program could not derive the energy &
               expression'/,&
            ' Select one of the following options and enter:'/,&
            '    1  Re-enter the term and configuration'/,&
            '    2  Use expression from  Eav as input'/,' &
	    3  STOP'/)")
           READ (5, *) ISELEC
           SELECT CASE (iselec)
              CASE (1)
                 ! get new atom and term
                 CYCLE get_a
             CASE (2)
                 ! use average energy
                 IJPTR = 0
                 TERM= 'AV'
                 CYCLE get_t
            CASE (3)
               Stop
           END SELECT
        ELSE
          EXIT get_a
        ENDIF
        END DO get_t
      END DO get_a


      ! .. define the spline step and basis size
      CALL get_spline_param(z)

      ! .. now we know both the nwf and ns
      CALL allocate_orbital_arrays
 
    END SUBROUTINE get_case
