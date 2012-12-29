! MDFF parallel Molecular Dynamics ... For Fun
! Copyright (C) 2011  F. Vasconcelos
!
! This program is free software; you can redistribute it and/or
! modify it under the terms of the GNU General Public License
! as published by the Free Software Foundation; either version 2
! of the License, or (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program; if not, write to the Free Software
! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

! ===== fmV =====

! ======= Hardware =======
! ======= Hardware =======

!*********************** MODULE control ***************************************
!
! main parameters of the code
!
!******************************************************************************
MODULE control 

  implicit none

  integer, SAVE :: myrank                    ! rank of the actual proc   
  integer, SAVE :: numprocs                  ! total number of procs    
  character*8,  SAVE   :: DATE
  character*10, SAVE   :: HOUR 

  logical, SAVE :: ltest
  logical, SAVE :: lstatic                   ! no MD                                                
  logical, SAVE :: lvnlist                   ! verlet list if .true.                            
  logical, SAVE :: lpbc                      ! PBC (periodic ...) if .true.                     
  logical, SAVE :: lminimg                   ! minimum convention image convention if .true.                     
  logical, SAVE :: lrestart                  ! restart or not if true strart from the velocities read in POSFF
  logical, SAVE :: lreduced                  ! print reduced thermo quantites by the number of atoms (natm)
  logical, SAVE :: lshiftpot                 ! shifted potential
  logical, SAVE :: lbmlj                     ! binary mixture lennard-jones potential
  logical, SAVE :: lcoulomb                  ! coulombic potential

  double precision, SAVE :: cutoff           ! small-range cutoff
  double precision, SAVE :: cutlongrange     ! longrange-range cutoff
  double precision, SAVE :: skindiff         ! verlet-list cutoff ( cutoff + skindiff )


  ! type of calculation : md, opt, vib, efg ...              
  character*60, SAVE :: calc                
  character*60, SAVE :: calc_allowed(12)    
  data calc_allowed / 'md' , 'opt' , 'vib' , 'vib+fvib' , 'vib+gmod' , 'vib+band' , &
                      'vib+dos' , 'efg' , 'efg+acf', 'efg+stat' , 'gr' , 'NL_field' /

  ! algorithm for long-range calculation
  character*60, SAVE :: longrange            
  character*60, SAVE :: longrange_allowed(2) 
  data longrange_allowed / 'ewald' , 'direct' /

  ! algorithm for gaussian distribution (no fondamental difference just for fun)  
  character*60, SAVE :: dgauss    
  character*60, SAVE :: dgauss_allowed(3) 
  data dgauss_allowed / 'boxmuller_basic', 'boxmuller_polar' , 'knuth' /

CONTAINS


!*********************** SUBROUTINE control_init ******************************
!
! Set default values, read  and check consistenstency of control parameters
!
!******************************************************************************

SUBROUTINE control_init ( MDFF )

  USE io_file,  ONLY :  ionode , stdin , stdout, kunit_OUTFF

  implicit none

  ! local
  character*132 :: filename
  integer :: ioerr
  CHARACTER (LEN=80) :: MDFF

  namelist /controltag/  lbmlj     , &
                         lcoulomb  , &
                         cutlongrange   , &
                         lvnlist   , &
                         lstatic   , &
                         lpbc      , &
                         lminimg   , &
                         lreduced  , & 
                         lshiftpot , &
                         ltest     , &
                         lrestart  , &
                         calc      , &
                         dgauss    , &
                         longrange , &
                         cutoff    , &
                         skindiff  

  ! ======================
  !  set default values
  ! ======================
  CALL control_default_tag
  ! =======================
  !  read control namelist
  ! =======================
  CALL getarg (1,filename)
  OPEN ( stdin , file = filename)
  READ ( stdin , controltag, iostat=ioerr)
  if ( ioerr .lt. 0 )  then
    if ( ionode ) WRITE ( stdout, '(a)') 'ERROR reading input_file : controltag section is absent'
    STOP
  elseif ( ioerr .gt. 0 )  then
    if ( ionode ) WRITE ( stdout, '(a)') 'ERROR reading input_file : controltag wrong tag',ioerr
    STOP    
  endif
  CLOSE ( stdin )
  ! ======================
  !      check tags
  ! ======================
  CALL control_check_tag
  ! ======================
  !  print info to output
  ! ======================
  CALL control_print_info( stdout , MDFF )
  CALL control_print_info( kunit_OUTFF , MDFF )

  return

END SUBROUTINE control_init



!*********************** SUBROUTINE control_default_tag ***********************
!
! set default values to control tag
!
!******************************************************************************

SUBROUTINE control_default_tag

  implicit none

  ! ================
  !  default values
  ! ================
  lbmlj         = .true.
  lcoulomb      = .false.
  lvnlist       = .true.
  lstatic       = .false.
  lpbc          = .true.
  lreduced      = .true.
  lshiftpot     = .true.
  ltest         = .false.
  lminimg       = .true.
  lrestart      = .false.
  calc          = 'md'
  dgauss        = 'boxmuller_basic'
  longrange     = 'ewald'
  cutoff        = 2.5d0
  cutlongrange  = 400.d0
  skindiff      = 0.1d0

  return 
 
END SUBROUTINE control_default_tag


!*********************** SUBROUTINE control_check_tag ***********************
!
! check control tag values
!
!******************************************************************************

SUBROUTINE control_check_tag

  USE io_file,  ONLY :  stdout , ionode

  implicit none

  ! local
  logical :: allowed
  integer :: i

  if ( .not. lbmlj .and. .not. lcoulomb ) then
   if ( ionode )  WRITE ( stdout , '(a)' ) 'ERROR controltag: lj or coulomb or both. In anyway make a choice !! '
   STOP
  endif

  ! ======
  !  calc
  ! ======
  do i = 1 , size( calc_allowed ) 
   if ( trim(calc) .eq. calc_allowed(i))  allowed = .true.
  enddo
  if ( .not. allowed ) then
      if ( ionode )  WRITE ( stdout , '(a)' ) 'ERROR controltag: calc should be ', calc_allowed
      STOP 
  endif
  allowed = .false.
  ! ===========
  !  longrange
  ! ===========
  do i = 1 , size( longrange_allowed )
   if ( trim(longrange) .eq. longrange_allowed(i))  allowed = .true.
  enddo
  if ( .not. allowed ) then
    if ( ionode ) WRITE ( stdout , '(a)' ) 'ERROR controltag: longrange should be ', longrange_allowed
  endif
  allowed = .false.
  ! =========
  !  dgauss
  ! =========
  do i = 1 , size( dgauss_allowed )
   if ( trim(dgauss) .eq. dgauss_allowed(i))  allowed = .true.
  enddo
  if ( .not. allowed ) then
    if ( ionode ) WRITE ( stdout , '(a)' ) 'ERROR controltag: dgauss should be ', dgauss_allowed
  endif


  return

END SUBROUTINE control_check_tag


!*********************** SUBROUTINE control_print_info ************************
!
! print information to standard output about general parameters
!
! note : should had more information about other control parameters
!
!******************************************************************************

SUBROUTINE control_print_info( kunit , MDFF )

  USE io_file,  ONLY :  ionode 

  implicit none
 
  !local 
  integer :: kunit
  integer :: status
#ifdef GFORTRAN
  character(len=80) :: host
#endif
  CHARACTER (LEN=80) :: MDFF
 
  ! ===============
  !  date time info      
  ! ===============
  CALL DATE_AND_TIME( DATE, HOUR)

  ! ===============
  !  hostname info      
  ! ===============
#ifdef GFORTRAN
  status = hostnm(host)
#endif 
  ! =================
  !  standard output
  ! =================

  if ( ionode ) call dumb_guy(kunit)
  if ( ionode ) then
     WRITE ( kunit ,'(a)')       "          ____    ____  ______   ________  ________  "
     WRITE ( kunit ,'(a)')       "         |_   \  /   _||_   _ `.|_   __  ||_   __  | "
     WRITE ( kunit ,'(a)')       "           |   \/   |    | | `. \ | |_ \_|  | |_ \_| "
     WRITE ( kunit ,'(a)')       "           | |\  /| |    | |  | | |  _|     |  _|    "
     WRITE ( kunit ,'(a)')       "          _| |_\/_| |_  _| |_.' /_| |_     _| |_     "
     WRITE ( kunit ,'(a)')       "         |_____||_____||______.'|_____|   |_____|    "
     WRITE ( kunit ,'(a)')       '  '
     WRITE ( kunit ,'(a)')       '============================================================='
     WRITE ( kunit ,'(a)')       ''
     WRITE ( kunit ,'(a)')       'MOLECULAR DYNAMICS ...for fun                 '
     WRITE ( kunit ,'(a)')       MDFF
     WRITE ( kunit ,'(a)')       'filipe.manuel.vasconcelos@gmail.com  '
     WRITE ( kunit ,'(a,i4,a)')  'Running on',numprocs,' nodes                  '
     if ( status == 0 ) then
#ifdef GFORTRAN
     WRITE ( kunit ,'(a,a)')     'host     : ',trim(host)
#endif
     endif
     WRITE ( kunit ,'(a,a4,a1,a2,a1,a2,a4,a2,a1,a2,a1,a2)') &
                                'date     : ',DATE(1:4),'/',DATE(5:6),'/',DATE(7:8),'   ',HOUR(1:2),&
                                ':',HOUR(3:4),':',HOUR(5:6)
     WRITE ( kunit ,'(a,a)')     'calc     : ', calc 
     WRITE ( kunit ,'(a)')       ''
  endif

  return 
 
END SUBROUTINE control_print_info

END MODULE control 
! ===== fmV =====
