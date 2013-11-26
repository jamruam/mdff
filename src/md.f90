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
! ======= Hardware =======

! *********************** MODULE md ********************************************
!> \brief 
!! module related to molecular dynamics calculation ( calc = 'md' )
! ******************************************************************************
MODULE md

  USE constants,                ONLY :  dp

  implicit none

  logical, SAVE :: lleapequi               !< leap-frog used in the equilibration part together with verlet -> vv + lf 

  integer :: npas                          !< number of time steps
  integer :: nequil                        !< number of equilibration steps
  integer :: nequil_period                 !< equilibration period
  integer :: spas                          !< save configuration each spas step 
  integer :: nprint                        !< print thermo info to standard output
  integer :: fprint                        !< print thermo info to file OSZIFF
  integer :: updatevnl                     !< number of verlet list update  
  integer :: itime

  real(kind=dp) :: dt                   !< time step
  real(kind=dp) :: temp                 !< temperature   
  real(kind=dp) :: Qnosehoover          !< Nose-Hoover Chain : Q mass 
  integer       :: nhc_yosh_order       !< Nose-Hoover Chain : order of the yoshida integrator 
  integer       :: nhc_mults            !< Nose-Hoover Chain : number of multiple timesteps 
  integer       :: nhc_n                !< Nose-Hoover Chain : length of the Nose-Hoover chain
  real(kind=dp) :: tauberendsen         !< characteristic time in berendsen thermostat (simple rescale if tauberendsen = dt )
  real(kind=dp) :: nuandersen           !< characteristic frequency in andersen thermostat
  real(kind=dp) , dimension ( : ) , allocatable :: vxi, xi              !< general coordinates of the thermostat (Nose-Hoover Chain)

  ! ================================================
  !     algorithm for dynamic integration
  ! ================================================
  character(len=60) :: integrator               !< integration method   
  character(len=60) :: integrator_allowed(8)    
  data                 integrator_allowed / 'nve-vv' , 'nve-lf', 'nve-be' ,  'nvt-and' , 'nvt-nh' , 'nvt-nhc2' , 'nvt-nhcn', 'nve-lfq'/

  ! ================================================
  !  velocity distribution (Maxwell-Boltzmann or uniform)
  ! ================================================
  character(len=60) :: setvel                   !< velocity distribution 
  character(len=60) :: setvel_allowed(2) 
  data setvel_allowed / 'MaxwBoltz', 'Uniform' /

CONTAINS


! *********************** SUBROUTINE md_init ***********************************
!> \brief
!! molecular dynamics initialisation of main parameters
! ******************************************************************************
SUBROUTINE md_init

  USE control,  ONLY :  lpbc , lstatic , calc
  USE io_file,  ONLY :  ionode ,stdin, stdout

  implicit none

  integer            :: ioerr
  character(len=132) :: filename

  namelist /mdtag/    integrator    , & 
                      setvel        , & 
                      npas          , & 
                      nequil        , &
                      nequil_period , & 
                      nprint        , & 
                      fprint        , & 
                      spas          , & 
                      dt            , &
                      temp          , & 
                      nuandersen    , & 
                      tauberendsen  , &
                      nhc_yosh_order, & 
                      nhc_mults     , & 
                      nhc_n          , & 
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
    io_node &
    WRITE ( stdout, '(a)') 'ERROR reading input_file : mdtag section is absent'
    STOP
  elseif ( ioerr .gt. 0 )  then
    io_node &
    WRITE ( stdout, '(a,i8)') 'ERROR reading input_file : mdtag wrong tag'
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

  return

END SUBROUTINE md_init 


! *********************** SUBROUTINE md_default_tag ****************************
!> \brief
!! set default values to md tag
! ******************************************************************************
SUBROUTINE md_default_tag

  implicit none
  
  ! =================
  !  default values
  ! =================
  lleapequi     = .false.
  integrator    = 'nve-vv'
  setvel        = 'MaxwBoltz'
  npas          = 10
  nequil        = 10
  nequil_period = 1
  nprint        = 1             
  fprint        = 1
  spas          = 1000           
  dt            = 0.0_dp
  temp          = 1.0_dp
  tauberendsen  = 0.0_dp
  nhc_yosh_order= 3 
  nhc_mults     = 2 

  return

END SUBROUTINE md_default_tag


! *********************** SUBROUTINE md_check_tag ******************************
!> \brief
!! check md tag values
! ******************************************************************************
SUBROUTINE md_check_tag

  USE control,  ONLY :  lstatic
  USE io_file,  ONLY :  ionode , stdout

  implicit none

  ! local
  logical :: allowed
  integer :: i

  if (dt.eq.0.0_dp .and. .not. lstatic ) then
    if ( ionode )  WRITE ( stdout ,'(a)') 'ERROR mdtag: timestep dt is zero'
    STOP
  endif

  ! =========================================
  !  scaling velocities berendsen, 
  !  if not defined simple velocity rescale 
  ! =========================================
  if (tauberendsen.eq.0.0_dp ) tauberendsen = dt

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
  if ( integrator .eq. 'nvt-nhc2' .and. Qnosehoover .eq. 0.0_dp ) then
     if ( ionode )  WRITE ( stdout ,'(a,f10.5)') 'ERROR mdtag: with integrator = "nvt-nhc2" Qnosehoover should be set',Qnosehoover
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

  ! allocation of thermostat coordinates
  if ( integrator .eq. 'nvt-nhc2' ) then 
    allocate ( vxi(2) , xi(2) )
  endif 
  if ( integrator .eq. 'nvt-nhcn' ) then 
    allocate ( vxi(nhc_n) , xi(nhc_n) )
  endif 
  ! initial conditions
!   vxi = 1.0_dp
!   xi = 0.0_dp
!  vxi(1) = 1.0_dp      

  return


END SUBROUTINE md_check_tag


! *********************** SUBROUTINE md_print_info *****************************
!> \brief
!! print general information to standard output for md control tag
! ******************************************************************************
SUBROUTINE md_print_info(kunit)

  USE control,          ONLY :  ltraj , lpbc , lstatic , lvnlist , lreduced , lminimg , itraj_start , itraj_period , itraj_format  
  USE io_file,          ONLY :  ionode 

  implicit none
  
  !local
  integer :: kunit

  if ( ionode ) then 
                                          separator(kunit)    
                                          blankline(kunit)    
                                          WRITE ( kunit ,'(a)')       'MD MODULE ... WELCOME'
                                          blankline(kunit)    
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
        if ( integrator .eq. 'nvt-nhc2' .or. &
             integrator .eq. 'nvt-nhcn' ) WRITE ( kunit ,'(a)')       'NVT ensemble --- velocity verlet integrator'
        if ( integrator .eq. 'nvt-nhc2' ) WRITE ( kunit ,'(a)')       ' + Nose Hoover chain 2 thermostat  (see Frenkel and Smit)'
        if ( integrator .eq. 'nvt-nhcn' ) WRITE ( kunit ,'(a)')       ' + Nose Hoover chain N thermostat  (see Martyna et al.)'
        if ( integrator .eq. 'nvt-nhc2' .or. & 
             integrator .eq. 'nvt-nhcn' ) WRITE ( kunit ,'(a,f10.5)') 'Qnosehoover                          = ',Qnosehoover
        if ( ( integrator .ne. 'nvt-and' )  .and. &
             ( integrator .ne. 'nvt-nhc2' ) .and. &
             ( integrator .ne. 'nvt-nhcn' ) .and. &
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
                                          blankline(kunit)    
  endif !ionode


  return

END SUBROUTINE md_print_info


END MODULE md 
! ===== fmV =====
