!======================================================================
    SUBROUTINE get_atom
!======================================================================
!
!    This routine gets information about the atomic state -- the
!    closed shells and the configuration from which it determines
!    the number of wave functions, and the electrons
!
!----------------------------------------------------------------------
      Use hf_atomic_state
      Use hf_orbitals
      USE hf_inout
      USE hf_param, ONLY: varied, vstring
      IMPLICIT NONE
      CHARACTER(LEN=80) :: string, form
      CHARACTER(LEN=3),DIMENSION(25) :: elc
      INTEGER :: i, j, start, itemp,ios, ip
      INTEGER, EXTERNAL :: lval
      REAL(KIND=8), DIMENSION(25) :: qsumc
      LOGICAL :: connected
 !-----------------------------------------------
 !
      INQUIRE(Unit=iuc, OPENED=connected)
      If (connected) then
        call read_bswc
      Else
         ! .. get ATOM, TERM, and Z
            WRITE (0, '(A)') ' Enter ATOM,TERM,Z (separated by blanks)'
            READ (5, *) atom, term, z
         ! .. get closed shells
         WRITE (0, *) 
         call get_closed_shells
    
         ! .. get configuration
         WRITE (0, '(/A,A/A)') ' Enter electrons outside CLOSED shells ', &
             '(blank line if none)', ' Example: 2s(1)2p(3)' 
         READ (5, '(A)') configuration
      END IF
 
      nwf = nclosd
      start = 1
      configuration = adjustl(configuration)
      If (configuration(1:1) /= '*') then
      DO
        i = start + index(configuration(start:),'(') - 1
        if (i >= start) then
           nwf = nwf +1
           elc(nwf)=configuration(start:i-1)
           elc(nwf) = adjustr(elc(nwf))
           j = start + index(configuration(start:),')') -1
           READ(configuration(i+1:j-1), *) qsumc(nwf)
           start = j+1
           CYCLE
        end if
        EXIT
      END DO
      End if
      call write_bswc
  
      ! .. now we know nwf, allocate memory
      call allocate_atomic_state
    
      el(1:nwf) = elc(1:nwf)
 
      ! .. store properites of the closed shells
      lmax=0
      DO i = 1,nwf
        Call el_nl(el(i),n(i),l(i))
        lmax = max(l(i),lmax)
        if (i <= nclosd) then
          qsum(I) = 2*(2*L(I)+1) 
          clsd(i) = .true.
        else
          qsum(i) = qsumc(i)
          if (qsum(i) == 2.d0*(2*L(i)+1)) then
             clsd(i)= .true.
          else
             clsd(i) = .false.
          end if
        end if
      END DO
 
      ! Determine Lagrange multiplier structure
      DO i = 1,nwf
          e(i,i) = 0.d0
          DO j = 1, i-1
             e(i,j) = 0.d0
             IF (L(I) .EQ. L(J)) then
               IF (clsd(i) .AND. clsd(j)) then
                  E(I,J) = 1.D-12
                  E(J,I) = 1.d-12
               Else
                  E(I,J) = 1.D-5
                  E(J,I) = 1.D-5
               END IF
             END IF
          End do
      End DO


      form="(16X,'Core =',5(1X,A3,'(',I4,')')/(22X,5(1X,A3,'(',I4,')')))"
      WRITE (UNIT=err,FMT="(//7X,'WAVE FUNCTIONS FOR  ',2A6,' Z =',F5.1,/)") &
                  ATOM, TERM, Z
      WRITE (UNIT=err,FMT=form) (el(i),int(qsum(i)),i=1,nclosd)
      WRITE (UNIT=err,FMT="(7X,A,3(1X,A3,'(',F4.1,')'))") 'Configuration =', &
                 (el(i),qsum(i),i=nclosd+1,nwf)

      ! .. write information to summary file (unit=3)
      WRITE (UNIT=log,FMT="(//7X,'HARTREE-FOCK WAVE FUNCTIONS FOR  ',2A6, &
            ' Z =',F5.1,/)")   ATOM, TERM, Z
      WRITE (UNIT=log,FMT=form) (el(i),int(qsum(i)),i=1,nclosd)
      WRITE (UNIT=log,FMT="(7X,A,3(1X,A3,'(',F4.1,')'))") 'Configuration =',&
                 (el(i),qsum(i),i=nclosd+1,nwf)

      ! .. order the orbitals by n and l
      iord = 0
      DO i = 1, NWF 
         IORD(i) = i
      END DO
      DO i = 1, nwf-1
          Do j = i+1, nwf
          IF (N(iord(i)) < N(iord(j)) .OR. &
             (N(iord(i))==N(iord(j)) .AND. L(iord(i)) < L(iord(j)))) CYCLE
          itemp = iord(i)
          IORD(i) = IORD(j)
          IORD(j) = itemp
          End DO
       END DO
       Do ip = 1, nwf
         i = iord(ip)
      end do

      ! .. determine orbitals to be varied, place at end
      call get_varied
!
!      ... Define an order for the functions to be iterated
!      ... Orbitals varied must be last.
!
     ! .. array of energy parameters

     CONTAINS
!======================================================================
     SUBROUTINE get_closed_shells
!======================================================================
!
!    Determines the closed shells and the number of electrons, nclosd
!
!----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER :: i,length,ib,ie
      CHARACTER(LEN=45) :: He='1s', &
                           Be='1s 2s', &
                           Ne='1s 2s 2p', &
                           Ar='1s 2s 2p 3s 3p', &
                           Kr='1s 2s 2p 3s 3p 3d 4s 4p', &
                           Xe='1s 2s 2p 3s 3p 3d 4s 4p 4d 5s 5p', &
                           Rn='1s 2s 2p 3s 3p 3d 4s 4p 4d 4f 5s 5p 5d 6s 6p'
      
      WRITE (0, '(A)')' List the CLOSED shells (separated by blanks, * for none)'
      READ(5,'(A)') string
      If (string(1:1) == '*') then
        nclosd = 0
        return
      End if
      string= adjustl(string) 
      Length = LEN(trim(string))                  
      ib = index(string(1:Length),'[')+1
      ie = index(string(1:Length),']')-1
      If (ib /= 0) then
         !  Append to defined string
         Select case(string(ib:ie))
           case ('He'); string = trim(He)//' '//string(ie+2:)
           case ('Be'); string = trim(Be)//' '//string(ie+2:)
           case ('Ne'); string = trim(Ne)//' '//string(ie+2:)
           case ('Ar'); string = trim(Ar)//' '//string(ie+2:)
           case ('Kr'); string = trim(Kr)//' '//string(ie+2:)
           case ('Xe'); string = trim(Xe)//' '//string(ie+2:)
           case ('Rn'); string = trim(Rn)//' '//string(ie+2:)
           case DEFAULT; WRITE(err, *) &
                     'Core must be one of He, Be, Ne, Ar, Kr, Xe, Rn '
                WRITE(err, *) 'Example:  [Ne] 3s'
         END select
      End if
      Do i =1,20
         READ(string,*, IOSTAT=ios) elc(1:i)
         If (ios /= 0) Exit
      End do
      nclosd = i-1
      ! right-justify
      do i = 1,nclosd
       elc(i) = adjustr(elc(i))
      end do
   End subroutine get_closed_shells
 
!======================================================================
     SUBROUTINE get_varied
!======================================================================
!
!    Determines the orbitals to be varied, possibly re-ordering
!
!----------------------------------------------------------------------
     IMPLICIT NONE
     LOGICAL :: done
     

     vary(1:nwf) = .false.
     INQUIRE(Unit=iup, OPENED=connected)
     If (.not. connected) then
        WRITE (UNIT=err, FMT='(/A,I3,A/(1X,18(1X,A3)))')  &
        '  There is(are) ', NWF,' orbital(s) as follows:', EL(1:nwf)
        WRITE (UNIT=err, FMT='(/A,A)') ' Orbitals to be varied: ', &
            'ALL/NONE/=i (last i) or "LIST"'
        READ (5, '(A)') STRING
      Else
        string(1:4) = varied
      END IF
      string = adjustl(string)
      IF (STRING(1:3)=='ALL' .OR. STRING(1:3)=='all') THEN 
         NIT = NWF 
         vary(1:nwf) = .true.
      ELSE IF (STRING(1:4)=='NONE' .OR. STRING(1:4)=='none') THEN 
         NIT = 0 
      ELSE IF (INDEX(STRING,'=') /= 0) THEN 
         J = INDEX(STRING,'=') 
         READ (STRING(J+1:), *) NIT 
         vary(nwf-nit+1:nwf) = .true.
      ELSE 
         Do
           If (vstring(1:4) == '    ') then
              WRITE(0, '(A)') &
                   'Enter space delimited list of orbitals to be varied'
              READ(5,'(A)') string
              string= adjustl(string)
              vstring = string
           Else
              string = vstring
           END if
           nit=0; done=.false.
           Do i =1,nwf
             READ(string,*, IOSTAT=ios) elc(1:i)
             IF (ios /= 0) then
               done = .true.
             ELSE
               elc(i) = adjustr(elc(i))
             ! .. orbitals to be varied should be at the end of the list
               Call eptr(el,elc(i),j,nwf)
               IF (j /= 0) THEN 
                  NIT = NIT + 1 
                  If (j < nwf) then
                    itemp = iord(j)
                    iord(j:nwf-1) = iord(j+1:nwf)
                    iord(nwf) = itemp
                  End if
                  vary(i) = .true.
               ELSE 
                 WRITE (0, *) ' Case must match: Re-enter '
               ENDIF
             END IF
             If (done) Exit
           END DO 
           Exit
         END DO
      END IF
     End subroutine get_varied

!======================================================================
     SUBROUTINE read_bswc
!======================================================================
!
!    Store input data into bsw.c, if it does not exist
!
!----------------------------------------------------------------------
     IMPLICIT NONE
     
     INQUIRE(Unit=iuc, OPENED=connected)
     If (connected) then
       READ(IUC, *) atom, term, z, nclosd, nwf 
       READ(IUC, *) elc(1:nclosd)
       READ(IUC, *) configuration
     End if
     END SUBROUTINE read_bswc
!======================================================================
     SUBROUTINE write_bswc
!======================================================================
!
!    Store input data into bsw.c
!
!----------------------------------------------------------------------
     IMPLICIT NONE
     
     INQUIRE(Unit=iuc, OPENED=connected)
     If (.not. connected) then
       OPEN(UNIT=IUC, FILE='bsw.c', STATUS='UNKNOWN', FORM='FORMATTED')
       WRITE(UNIT=IUC, FMT='(1X,A6,A6,F6.1,2I6)') &
             trim(atom), trim(term), Z, nclosd, nwf 
       IF (nclosd > 0) then
        WRITE(IUC, *) elc(1:nclosd)
       Else
         WRITE(iuc, FMT='(A)') '*'
       END if
       WRITE(IUC, FMT='(A)') configuration
     End if
     END SUBROUTINE write_bswc
     END SUBROUTINE get_atom
!=======================================================================
  SUBROUTINE eptr(el, elsymb, iel, nwf) 
!=======================================================================
!
!   Determines the position of the electron in the electron list
!   Zero if not found.
!      
!----------------------------------------------------------------------
!
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: nwf
    CHARACTER(LEN=3), DIMENSION(nwf), INTENT(IN) :: el
    CHARACTER(LEN=3), INTENT(IN) :: elsymb
    INTEGER, INTENT(OUT) :: iel

    INTEGER :: i
    iel = 0
    do  i=1,nwf
      if (el(i) .eq. elsymb ) then
        iel = i
        return
      endif
    end do
    end
