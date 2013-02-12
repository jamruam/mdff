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

MODULE md

  implicit none

  logical, SAVE :: ltraj                   ! save trajectory                                    
  logical, SAVE :: lleapequi               ! leap-frog used in the equilibration part together with verlet -> vv + lf 

  integer :: npas                          ! number of time steps
  integer :: nequil                        ! number of equilibration steps
  integer :: nequil_period                 ! equilibration period
  integer :: spas                          ! save configuration each spas step 
  integer :: nprint                        ! print thermo info to standard output
  integer :: fprint                        ! print thermo info to file OSZIFF
  integer :: itraj_start                   ! write trajectory from step itraj_start
  integer :: itraj_period                  ! write trajectory each itraj_period steps 
  integer :: itraj_format                  ! choose the trajectory format ( = 0 BINARY, = 1 FORMATED)
  integer :: updatevnl                     ! number of verlet list update  
  integer :: itime

  double precision :: dt                   ! time step
  double precision :: temp                 ! temperature   
  double precision :: Qnosehoover          ! Q parameter in Nose-Hoover Chain 
  double precision :: tauberendsen         ! characteristic time in berendsen thermostat (simple rescale if tauberendsen = dt )
  double precision :: nuandersen           ! characteristic frequency in andersen thermostat
  double precision :: vxi1, vxi2, xi1, xi2 ! extra variables of Nose-Hoover Chain

  character*60 :: integrator               ! algorithm for dynamic integration             
  character*60 :: integrator_allowed(7) 
  data integrator_allowed / 'nve-vv' , 'nve-lf', 'nve-be' ,  'nvt-and' , 'nvt-nh' , 'nvt-nhc2' , 'nve-lfq'/

  character*60 :: setvel                   ! velocity distribution (Maxwell-Boltzmann or uniform)    
  character*60 :: setvel_allowed(2) 
  data setvel_allowed / 'MaxwBoltz', 'Uniform' /

CONTAINS


!*********************** SUBROUTINE md_init ***********************************
!
! molecular dynamics initialisation of main parameters
!
!******************************************************************************

SUBROUTINE md_init

  USE control,  ONLY :  lpbc , lstatic , calc
  USE io_file,  ONLY :  ionode ,stdin, stdout, kunit_OUTFF

  implicit none

  character*132 :: filename
  integer :: ioerr

  namelist /mdtag/    ltraj         , &
                      integrator    , & 
                      setvel        , & 
                      npas          , & 
                      nequil        , &
                      nequil_period , & 
                      nprint        , & 
                      fprint        , & 
                      itraj_start   , & 
                      itraj_period  , & 
                      spas          , & 
                      dt            , &
                      temp          , & 
                      nuandersen    , & 
                      tauberendsen  , & 
                      Qnosehoover      

  if ( calc .ne. 'md' ) return
  ! ======================
  !  mdtag default values
  ! ======================
  CALL md_default_tag 
  ! =====================
  !  read mdtag namelist
  ! =====================
  CALL getarg (1,filename)
  OPEN ( stdin , file = filename)
  READ ( stdin , mdtag ,iostat =ioerr)
  if ( ioerr .lt. 0 )  then
    if ( ionode ) WRITE ( stdout, '(a)') 'ERROR reading input_file : mdtag section is absent'
    STOP
  elseif ( ioerr .gt. 0 )  then
    if ( ionode ) WRITE ( stdout, '(a)') 'ERROR reading input_file : mdtag wrong tag'
    STOP
  endif
  CLOSE  ( stdin )
  ! ======================
  !  check mdtag namelist
  ! ======================
  CALL md_check_tag
  ! ===================
  !  print mdtag info
  ! ===================
  CALL md_print_info(stdout)
  CALL md_print_info(kunit_OUTFF)

  return

END SUBROUTINE md_init 


!*********************** SUBROUTINE md_default_tag ****************************
!
! set default values to md tag
!
!******************************************************************************

SUBROUTINE md_default_tag

  implicit none
  
  ! =================
  !  default values
  ! =================
  ltraj         = .false.
  lleapequi     = .false.
  integrator    = 'nve-vv'
  setvel        = 'MaxwBoltz'
  npas          = 10
  nequil        = 10
  nequil_period = 1
  nprint        = 1             
  fprint        = 1
  itraj_start   = 1          
  itraj_period  = 10000
  itraj_format  = 1
  spas          = 1000           
  dt            = 0.001d0
  temp          = 1.0d0
  tauberendsen  = 0.0d0

  return

END SUBROUTINE md_default_tag


!*********************** SUBROUTINE md_check_tag ******************************
!
! check md tag values
!
!******************************************************************************
 
SUBROUTINE md_check_tag

  USE io_file,  ONLY :  ionode , stdout

  implicit none

  ! local
  logical :: allowed
  integer :: i

  ! =========================================
  !  scaling velocities berendsen, 
  !  if not defined simple velocity rescale 
  ! =========================================
  if (tauberendsen.eq.0.0D0 ) tauberendsen = dt

  ! ====================
  !  check integrator
  ! ====================
  do i = 1 , size (integrator_allowed)
   if (trim (integrator) .eq. integrator_allowed(i) ) allowed = .true.
  enddo
  if (  .not. allowed ) then
    if ( ionode )  WRITE ( stdout ,'(a)') 'ERROR mdtag: integrator should be ', integrator_allowed
    STOP 
  endif
  allowed = .false.
  ! ===============
  !  check setvel
  ! ===============
  do i = 1 , size (setvel_allowed)
   if (trim (setvel) .eq. setvel_allowed(i) ) allowed = .true.
  enddo
  if (  .not. allowed ) then
    if ( ionode )  WRITE ( stdout ,'(a)') 'ERROR mdtag: setvel should be ', setvel_allowed
    STOP
  endif

  ! ===================
  !  check Qnosehoover
  ! ===================
  if ( integrator .eq. 'nvt-nhc2' .and. Qnosehoover .eq. 0.0d0 ) then
     if ( ionode )  WRITE ( stdout ,'(a,f10.5)') 'ERROR mdtag: with integrator = "nvt-nhc2" Qnoseehoover should be set',Qnosehoover
    STOP
  endif

  if (integrator.eq.'nvt-nh' ) then
    if ( ionode )  WRITE ( stdout ,'(a)') 'ERROR mdtag: integrator = "nvt-nh" not yet implemented try nvt-nhc2 '
    if ( ionode )  WRITE ( stdout ,'(a)') integrator
    STOP 
  endif

  ! ============================
  !  check if leap frog is 
  !  wanted with equilibration 
  !  then we use nve-vv for the 
  !  equilibration period
  ! ============================
  if ( integrator.eq.'nve-lfq' ) then
    lleapequi = .true.
    integrator = 'nve-vv' 
  endif

  return


END SUBROUTINE md_check_tag


!*********************** SUBROUTINE md_print_info *************************
!
! print general information to standard output for md control tag
!
!**************************************************************************

SUBROUTINE md_print_info(kunit)

  USE control,  ONLY :  lpbc , lstatic , lvnlist , lreduced , lminimg 
  USE io_file,  ONLY :  ionode 

  implicit none
  
  !local
  integer :: kunit

  if ( ionode ) then 
                                          WRITE ( kunit ,'(a)')       '=============================================================' 
                                          WRITE ( kunit ,'(a)')       ''
      if (lstatic) then
                                          WRITE ( kunit ,'(a)')       'static  calculation ....boring                 '
      else
        if ( .not. lpbc )                 WRITE ( kunit ,'(a)')       'NO PERIODIC BOUNDARY CONDITIONS (cubic md_cell)'
        if ( lpbc )                       WRITE ( kunit ,'(a)')       'periodic boundary conditions in cubic cell     '  
        if ( lminimg )                    WRITE ( kunit ,'(a)')       'using minimum image convention                 '  
        if ( .not. lminimg )              WRITE ( kunit ,'(a)')       'no minimum image convention                    '  
        if ( lvnlist )                    WRITE ( kunit ,'(a)')       'verlet list used '  
        if ( lreduced )                   WRITE ( kunit ,'(a)')       'units reduced by the number of atom'
        if ( .not.lstatic)                WRITE ( kunit ,'(a)')       'dynamic calculation'
        if ( integrator .eq. 'nve-lf')    WRITE ( kunit ,'(a)')       'NVE ensemble --- leap-frog integrator          '
        if ( integrator .eq. 'nve-be')    WRITE ( kunit ,'(a)')       'NVE ensemble --- beeman integrator             '

        if ( integrator .eq. 'nve-vv' .and. lleapequi )   &
                                          WRITE ( kunit ,'(a)')       'NVE ensemble --- velocity verlet (equil) + leap-frog integrator        '
        if ( integrator .eq. 'nve-vv'.and.  .not. lleapequi )  &
                                          WRITE ( kunit ,'(a)')       'NVE ensemble --- velocity verlet integrator    '
        if ( integrator .eq. 'nvt-and' )  WRITE ( kunit ,'(a)')       'NVT ensemble --- velocity verlet integrator    ' 
        if ( integrator .eq. 'nvt-and' )  WRITE ( kunit ,'(a)')       ' + Andersen thermostat'
        if ( integrator .eq. 'nvt-and' )  WRITE ( kunit ,'(a,f10.5)') 'nuandersen                         = ',nuandersen  
        if ( integrator .eq. 'nvt-nh' )   WRITE ( kunit ,'(a)')       'NVT ensemble --- velocity verlet integrator'
        if ( integrator .eq. 'nvt-nh' )   WRITE ( kunit ,'(a)')       ' + Nose Hoover thermostat'
        if ( integrator .eq. 'nvt-nhc2' ) WRITE ( kunit ,'(a)')       'NVT ensemble --- velocity verlet integrator'
        if ( integrator .eq. 'nvt-nhc2' ) WRITE ( kunit ,'(a)')       ' + Nose Hoover chain 2 thermostat  (see Frenkel and Smit)'
        if ( integrator .eq. 'nvt-nhc2' ) WRITE ( kunit ,'(a,f10.5)') 'Qnosehoover                          = ',Qnosehoover
        if ( ( integrator .ne. 'nvt-and' )  .and. &
             ( integrator .ne. 'nvt-nhc2' ) .and. &
             ( integrator .ne. 'nvt-nh' ) ) then
              if ( nequil .eq. 0 ) then
                                          WRITE ( kunit ,'(a)')       'with no equilibration          '
              else
                                          WRITE ( kunit ,'(a)')       'with equilibration:             '
                                          WRITE ( kunit ,'(a)')       'berendsen scaling ( is not canonical ...and so NVE)'
               if ( integrator .eq. 'nve-lf' ) &
                                          WRITE ( kunit ,'(a)')       'WARNING WITH nve-lf no equilibration possible'
              endif !nequil
        endif !integrator
      endif !static
                                          WRITE ( kunit ,'(a,i10)')   'number of steps                      = ',npas
                                          WRITE ( kunit ,'(a,e12.5)') 'timestep                             = ',dt
                                          WRITE ( kunit ,'(a,e12.5)') 'time range                           = ',dt*npas
                                          WRITE ( kunit ,'(a,f10.5)') 'temperature                          = ',temp
      if ( integrator .eq. 'nve-vv' .and. nequil .ne. 0 ) then           
                                          WRITE ( kunit ,'(a,i10)')   'number of equilibration steps        = ',nequil
                                          WRITE ( kunit ,'(a,i10)')   'equilibration period                 = ',nequil_period
      endif 
      if ( nequil .ne. 0 )                WRITE ( kunit ,'(a,e12.5)') 'Berendsen thermo time scale          = ',tauberendsen
      if ( nequil .ne. 0 .and. tauberendsen .eq. dt )   &
                                          WRITE ( kunit ,'(a)')       'tauberendsen = dt -> simple rescale'
                                          WRITE ( kunit ,'(a,i10)')   'print thermo  periodicity            = ',nprint
      if ( ltraj )                   then    
                                          WRITE ( kunit ,'(a,i10)')   'save trajectory from step            = ',itraj_start
                                          WRITE ( kunit ,'(a,i10)')   'saved trajectory periodicity         = ',itraj_period
      if ( itraj_format .eq. 0 )          WRITE ( kunit ,'(a,I7)')    'trajectory format                    : BINARY'
      if ( itraj_format .ne. 0 )          WRITE ( kunit ,'(a,I7)')    'trajectory format                    : FORMATTED'
      endif       
                                          WRITE ( kunit ,'(a)')       ''
  endif !ionode


  return

END SUBROUTINE md_print_info


!*********************** SUBROUTINE write_traj_xyz ****************************
!
! write trajectory (pos, vel, for) to TRAJFF file
!
!******************************************************************************

SUBROUTINE write_traj_xyz 

  USE io_file,  ONLY :  ionode , kunit_TRAJFF , stdout
  USE config,   ONLY :  system , natm , natmi , ntype , &
                        rx , ry , rz , vx , vy , vz , fx , fy , fz , atype , atypei , simu_cell
  USE cell

  implicit none

  ! local
  integer :: ia , it , i , j

  CALL periodicbc ( natm , rx , ry , rz , simu_cell )

  if ( ionode ) then
    WRITE ( kunit_TRAJFF , '(i6)' ) natm
    WRITE ( kunit_TRAJFF , '(a)' ) system
    do i = 1 , 3
      WRITE ( kunit_TRAJFF , '(3f20.12)' ) (simu_cell%A(i,j),j=1,3) 
    enddo
    WRITE ( kunit_TRAJFF , '(i4)' ) ntype
    WRITE ( kunit_TRAJFF , * ) ( atypei ( it ) , it = 1 , ntype )
    WRITE ( kunit_TRAJFF , * ) ( natmi  ( it ) , it = 1 , ntype )
    
    do ia = 1 , natm
      WRITE ( kunit_TRAJFF , 200 ) atype(ia), rx ( ia ) , ry ( ia ) , rz ( ia ) , &
                                              vx ( ia ) , vy ( ia ) , vz ( ia ) , & 
                                              fx ( ia ) , fy ( ia ) , fz ( ia )
    enddo
  endif

 200 FORMAT(A2,9E20.12)

 return

END SUBROUTINE write_traj_xyz

END MODULE md 
! ===== fmV =====
