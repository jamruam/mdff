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
#include "symbol.h"
!#define GFORTRAN
! ======= Hardware =======

! *********************** MODULE control ***************************************
!> \brief
!! Module related to main parameters of the code
! ******************************************************************************
MODULE control 

  USE constants , ONLY : dp 
  USE mpimdff
      
  implicit none

  integer,           SAVE :: nprop            !< properties period
  integer                 :: itraj_start      !< write trajectory from step itraj_start
  integer                 :: itraj_period     !< write trajectory each itraj_period steps 
  integer                 :: itraj_format     !< choose the trajectory format ( = 0 BINARY, = 1 FORMATED)
  integer                 :: iscff_format     !< format of optimized configuration 0 = BINARY
  character(len=8),  SAVE :: DATE             !< execution DATE
  character(len=10), SAVE :: HOUR             !< execution HOUR

  logical,           SAVE :: ltraj            !< save trajectory                                    
  logical,           SAVE :: lstatic          !< no MD                                                
  logical,           SAVE :: lvnlist          !< verlet list if .true.                            
  logical,           SAVE :: lpbc             !< PBC (periodic ...) if .true.                     
  logical,           SAVE :: lminimg          !< minimum convention image convention if .true.                     
  logical,           SAVE :: lrestart         !< restart or not if true strart from the velocities read in POSFF
  logical,           SAVE :: lreduced         !< print reduced thermo quantites by the number of atoms (natm)
  logical,           SAVE :: lshiftpot        !< shifted potential
  logical,           SAVE :: lbmlj            !< binary mixture lennard-jones potential
  logical,           SAVE :: lharm            !< harmonic oscilaltor ( to test integration )
  logical,           SAVE :: lmorse           !< morse potential 
  logical,           SAVE :: lcoulomb         !< coulombic potential
  logical,           SAVE :: lsurf            !< add surface contribution in electrostatic quantities  
  logical,           SAVE :: ltest            !< testing flag

  real(kind=dp),     SAVE :: cutlongrange     !< longrange cutoff
  real(kind=dp),     SAVE :: cutshortrange    !< shortrange cutoff
  real(kind=dp),     SAVE :: skindiff         !< verlet-list cutoff ( cutoff + skindiff )
  TYPE(decomposition), SAVE :: kpt_dec


  ! =====================================================
  !   type of calculation : md, opt, vib, efg ...              
  ! =====================================================
  character(len=60), SAVE :: calc             !< type of calculation : md, opt, vib, efg ...  
  character(len=60), SAVE :: calc_allowed(12)    
  data calc_allowed / 'md'       , 'opt'     , 'vib'        , 'vib+fvib'       , 'vib+gmod' , &
                      'vib+band' , 'vib+dos' , 'efg'        , 'efg+acf'        , 'gr'       , &
                      'vois1'    , 'rmc'     /

  ! =====================================================
  ! algorithm for long-range calculation
  ! =====================================================
  character(len=60), SAVE :: longrange        !< algorithm for long-range interaction  
  character(len=60), SAVE :: longrange_allowed(2) 
  data longrange_allowed / 'ewald' , 'direct' /

  ! =====================================================
  !   algorithm for gaussian distribution 
  !   (no fondamental difference just for fun)  
  ! =====================================================
  character(len=60), SAVE :: dgauss           !< algorithm for gaussian distribution
  character(len=60), SAVE :: dgauss_allowed(3) 
  data dgauss_allowed / 'boxmuller_basic', 'boxmuller_polar' , 'knuth' /

  ! =====================================================
  !  format of TRAJFF allowed  
  ! =====================================================
  character(len=3), SAVE :: itraj_save
  character(len=3), SAVE :: iscff_save
  character(len=3), SAVE :: itraj_save_allowed(4)
  data itraj_save_allowed / 'rvf' , 'rnn' , 'rnf' , 'rvn' /
  
CONTAINS

! *********************** SUBROUTINE control_init ******************************
!
!> \brief
!! Initialization of control parameters.
!! Set default values, read  and check consistenstency of control parameters
!
! ******************************************************************************
SUBROUTINE control_init ( MDFF )

  USE io_file,  ONLY :  ionode , stdin , stdout

  implicit none

  ! local
  integer            :: ioerr
  character(len=80)  :: MDFF
  character(len=132) :: filename

  namelist /controltag/  lbmlj        , &
                         lmorse       , &
                         lcoulomb     , &
                         lsurf        , &
                         lharm        , &
                         ltraj        , &
                         cutlongrange , &
                         cutshortrange, &
                         lvnlist      , &
                         lstatic      , &
                         lpbc         , &
                         lminimg      , &
                         lreduced     , & 
                         lshiftpot    , &
                         ltest        , &
                         lrestart     , &
                         calc         , &
                         dgauss       , &
                         longrange    , &
                         itraj_start  , & 
                         itraj_period , & 
                         itraj_format , & 
                         itraj_save   , & 
                         iscff_format , & 
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
    io_node WRITE ( stdout, '(a)') 'ERROR reading input_file : controltag section is absent'
    STOP
  elseif ( ioerr .gt. 0 )  then
    io_node WRITE ( stdout, '(a,i8)') 'ERROR reading input_file : controltag wrong tag',ioerr
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

  return

END SUBROUTINE control_init


! *********************** SUBROUTINE control_default_tag ***********************
!
!> \brief
!! set default values to control tag
!
! ******************************************************************************
SUBROUTINE control_default_tag

  implicit none

  ! ================
  !  default values
  ! ================
  lbmlj         = .false.
  lmorse        = .false.
  lcoulomb      = .false.
  lsurf         = .false.
  lharm         = .false.
  ltraj         = .false.
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
  skindiff      = 0.15_dp
  cutshortrange = 0.0_dp
  cutlongrange  = 0.0_dp
  nprop         = 1
  itraj_start   = 1          
  itraj_period  = 10000
  itraj_format  = 1
  itraj_save    = 'rnn'
  iscff_format  = 1

  return 
 
END SUBROUTINE control_default_tag

! *********************** SUBROUTINE control_check_tag *************************
!
!> \brief
!! check control tag values
!
! ******************************************************************************
SUBROUTINE control_check_tag

  USE io_file,  ONLY :  stdout , ionode

  implicit none

  ! local
  logical :: allowed
  integer :: i

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
    io_node WRITE ( stdout , '(a)' ) 'ERROR controltag: longrange should be ', longrange_allowed
  endif
  allowed = .false.
  ! =========
  !  dgauss
  ! =========
  do i = 1 , size( dgauss_allowed )
   if ( trim(dgauss) .eq. dgauss_allowed(i))  allowed = .true.
  enddo
  if ( .not. allowed ) then
    io_node WRITE ( stdout , '(a)' ) 'ERROR controltag: dgauss should be ', dgauss_allowed
  endif
  ! =========
  ! traj_save  
  ! =========
  do i = 1 , size( itraj_save_allowed )
   if ( trim(itraj_save) .eq. itraj_save_allowed(i))  allowed = .true.
  enddo
  if ( .not. allowed ) then
    io_node WRITE ( stdout , '(a)' ) 'ERROR controltag: itraj_save should be ', itraj_save_allowed
  endif
  do i = 1 , size( itraj_save_allowed )
   if ( trim(iscff_save) .eq. itraj_save_allowed(i))  allowed = .true.
  enddo
  if ( .not. allowed ) then
    io_node WRITE ( stdout , '(a)' ) 'ERROR controltag: iscff_save should be ', itraj_save_allowed
  endif


  if ( calc .ne. 'md' ) return 

  if ( .not. lbmlj .and. .not. lcoulomb .and. .not. lmorse .and. .not. lharm ) then
   if ( ionode )  WRITE ( stdout , '(a)' ) 'ERROR controltag: lj, harm , morse or coulomb or all of them . Anyway make a choice !! '
   STOP
  endif

  if ( ( lbmlj .or. lmorse ) .and. cutshortrange .eq. 0.0_dp ) then
    io_node WRITE ( stdout , '(a)' ) 'controltag: cutshortrange is null', cutshortrange
  endif

  if ( lcoulomb .and. ( cutlongrange .eq. 0.0_dp .or. cutshortrange .eq. 0.0_dp ) ) then
    io_node WRITE ( stdout , '(a)' ) 'controltag: cutlongrange or cutshortrange is null', cutlongrange, cutshortrange
  endif

  return

END SUBROUTINE control_check_tag

! *********************** SUBROUTINE control_print_info ************************
!
!> \brief
!! print information to standard output about general parameters
!
!! \note
!! should have more information about other control parameters
!
! ******************************************************************************
SUBROUTINE control_print_info( kunit , MDFF )

  USE io_file,  ONLY :  ionode 
  USE mpimdff,  ONLY :  numprocs

  implicit none
 
  !local 
  integer :: kunit
  integer :: status
#ifdef GFORTRAN
  character(len=80) :: HOST 
#endif
  character(len=80) :: MDFF
  character(len=30) :: user_name
#ifndef HOST
#define HOST "unknown"
#endif
 
  ! ===============
  !  date time info      
  ! ===============
  CALL DATE_AND_TIME( DATE, HOUR)
  ! ===============
  !  user name 
  ! ===============
  CALL GETENV ( 'USER', user_name )

#ifdef GFORTRAN
  ! ===============
  !  hostname info      
  ! ===============
  status = hostnm(host)
#else   
  status = 0
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
     blankline(kunit)
     separator(kunit)
     blankline(kunit)
     WRITE ( kunit ,'(a)')       'MOLECULAR DYNAMICS ...for fun                 '
     WRITE ( kunit ,'(a)')       MDFF
     WRITE ( kunit ,'(a)')       'filipe.manuel.vasconcelos@gmail.com  '
     WRITE ( kunit ,'(a,i4,a)')  'Running on  : ',numprocs,' nodes                  '
     WRITE ( kunit ,'(a,a)')     'by user     : ',user_name
     if ( status == 0 ) then
#ifdef GFORTRAN
     WRITE ( kunit ,'(a,a)')     'host        : ',trim(HOST)
#else
     WRITE ( kunit ,'(a,a)')     'host        : ',HOST
#endif
     endif
     WRITE ( kunit ,'(a,a4,a1,a2,a1,a2,a4,a2,a1,a2,a1,a2)') &
                                'date        : ',DATE(1:4),'/',DATE(5:6),'/',DATE(7:8),'   ',HOUR(1:2),&
                                ':',HOUR(3:4),':',HOUR(5:6)
     blankline(kunit)
     WRITE ( kunit ,'(a)'  )     'CONTROL MODULE ... WELCOME'
     blankline(kunit)
     WRITE ( kunit ,'(a,a)')     'calc        =  ', calc 
     WRITE ( kunit ,'(a,l2)')    'lbmlj       = ', lbmlj 
     WRITE ( kunit ,'(a,l2)')    'lmorse      = ', lmorse 
     WRITE ( kunit ,'(a,l2)')    'lcoulomb    = ', lcoulomb
     WRITE ( kunit ,'(a,l2)')    'lsurf       = ', lsurf
     WRITE ( kunit ,'(a,l2)')    'lvnlist     = ', lvnlist
     WRITE ( kunit ,'(a,l2)')    'lstatic     = ', lstatic
     WRITE ( kunit ,'(a,l2)')    'lpbc        = ', lpbc
     WRITE ( kunit ,'(a,l2)')    'lminimg     = ', lminimg
     WRITE ( kunit ,'(a,l2)')    'lreduced    = ', lreduced 
     WRITE ( kunit ,'(a,l2)')    'lshiftpot   = ', lshiftpot
     WRITE ( kunit ,'(a,l2)')    'lrestart    = ', lrestart 
     blankline(kunit)
  endif

  return 
 
END SUBROUTINE control_print_info

END MODULE control 
! ===== fmV =====
