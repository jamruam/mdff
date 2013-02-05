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

!*********************** MODULE time ******************************************
!
!  this MODULE implements a small timer utility to time
!  SUBROUTINE 
!  (see START_TIMING)
!  from VASP
!
!******************************************************************************
MODULE time

  implicit none
  
  double precision :: timetot        ! engforce SUBROUTINE       
  double precision :: forcetimetot   ! engforce SUBROUTINE       
  double precision :: fcoultimetot1  ! engforce_coul direct space      
  double precision :: fcoultimetot2  ! engforce_coul fourier space      
  double precision :: fcoultimetot3  ! engforce_coul comm only 
  double precision :: efgtimetot1    ! EFG direct space           
  double precision :: efgtimetot2    ! EFG fourier space        
  double precision :: efgtimetot3    ! communication only         
  double precision :: strftimetot    ! facteur de structure
  double precision :: mdsteptimetot  !                            
  double precision :: vacftimetot    !                            
  double precision :: vacftimetot2   !                            
  double precision :: vnlisttimetot  !                            
  double precision :: msdtimetot     !                            
  double precision :: grtimetot      !                            
  double precision :: grtimetot_comm      !                            
  double precision :: opttimetot     !                          
  double precision :: vibtimetot     !                           
  double precision :: fvibtimetot    !                           
  double precision :: timetmp        ! use for opt               

CONTAINS


!*********************** SUBROUTINE time_init *********************************
!
! set time to zero
!
!******************************************************************************

SUBROUTINE time_init

  implicit none

  timetot       = 0.0d0
  efgtimetot1   = 0.0d0
  efgtimetot2   = 0.0d0
  efgtimetot3   = 0.0d0
  strftimetot   = 0.0d0 
  mdsteptimetot = 0.0d0
  vnlisttimetot = 0.0d0
  msdtimetot    = 0.0d0
  grtimetot     = 0.0d0 
  grtimetot_comm= 0.0d0 
  forcetimetot  = 0.0d0
  fcoultimetot1 = 0.0d0
  fcoultimetot2 = 0.0d0
  fcoultimetot3 = 0.0d0
  vibtimetot    = 0.0d0
  opttimetot    = 0.0d0
  vacftimetot   = 0.0d0                             
  vacftimetot2  = 0.0d0                             

  return 
 
END SUBROUTINE time_init

!*********************** SUBROUTINE print_time_info ***************************
!
! print final time info
!
!******************************************************************************

SUBROUTINE print_time_info ( kunit )  

  USE io_file,  ONLY :  ionode 
  USE control,  ONLY :  longrange , lcoulomb , calc

  implicit none

  !global
  integer, intent (in) :: kunit

  ! local
  double precision :: dstime , estime
  double precision :: fcoultime_es , fcoultime_ds

  dstime = efgtimetot1 + efgtimetot3               ! direct
  estime = strftimetot + efgtimetot1 + efgtimetot2 + efgtimetot3 ! ewald
  fcoultime_es = fcoultimetot1 + fcoultimetot2 + fcoultimetot3
  fcoultime_ds = fcoultimetot1 + fcoultimetot3 

  !timing information at the end
  if ( ionode ) then
                                       WRITE ( kunit ,110)             'TOTAL'  , timetot
    if ( mdsteptimetot .ne. 0.0d0 )    WRITE ( kunit ,110)             'MD'     , mdsteptimetot
    if ( opttimetot    .ne. 0.0d0 )    WRITE ( kunit ,110)             'OPT'    , opttimetot
    if ( vibtimetot    .ne. 0.0d0 )    WRITE ( kunit ,110)             'VIB'    , vibtimetot
    if ( fvibtimetot   .ne. 0.0d0 )    WRITE ( kunit ,110)             'FVIB'   , fvibtimetot
    if ( dstime .ne.0.0d0 .and. &
         efgtimetot2   .eq. 0.0d0 )    WRITE ( kunit ,110)             'EFG_DS' , dstime
    if ( efgtimetot2   .ne. 0.0d0 )    WRITE ( kunit ,110)             'EFG_ES' , estime
    if ( fcoultime_ds  .ne. 0.0d0 .and. &
         fcoultimetot2 .eq. 0.0d0 )    WRITE ( kunit ,110)             'FIELD_DS', fcoultime_ds
    if ( fcoultimetot2 .ne. 0.0d0 )    WRITE ( kunit ,110)             'FIELD_ES', fcoultime_es
    if ( grtimetot     .ne. 0.0d0 )    WRITE ( kunit ,110)             'GR'                , grtimetot 
    if ( grtimetot_comm.ne. 0.0d0 )    WRITE ( kunit ,110)             'GR(comm only)'     , grtimetot_comm 
    if ( msdtimetot    .ne. 0.0d0 )    WRITE ( kunit ,110)             'MSD'    , msdtimetot
    if ( vacftimetot   .ne. 0.0d0 )    WRITE ( kunit ,110)             'VACF'   , vacftimetot
    if ( vacftimetot2  .ne. 0.0d0 )    WRITE ( kunit ,110)             'VACF2'  , vacftimetot2
                                       WRITE ( kunit ,'(a)')           ''
                                       WRITE ( kunit ,'(a)')           '======================='
                                       WRITE ( kunit ,'(a)')           'main subroutines:'
                                       WRITE ( kunit ,'(a)')           '======================='
    if ( calc .eq. 'md' ) then
                                       WRITE ( kunit ,'(a)')           'MD:'
                                       WRITE ( kunit ,'(a)')           '======================='
      if ( forcetimetot  .ne. 0.0d0 )  WRITE ( kunit ,110)             'engforce_bmlj      ' , forcetimetot        
      if ( vnlisttimetot .ne. 0.0d0 )  WRITE ( kunit ,110)             'vnlistcheck        ' , vnlisttimetot
    endif
    
    if ( longrange .eq. 'direct' )   then
                                       WRITE ( kunit ,'(a)')           '======================='
      if ( dstime        .ne. 0.0d0 )  WRITE ( kunit ,'(a)')           'EFG:'
                                       WRITE ( kunit ,'(a)')           '======================='
      if ( efgtimetot1   .ne. 0.0d0 )  WRITE ( kunit ,110)             'efg_DS(efg  only)       ', efgtimetot1 
      if ( efgtimetot3   .ne. 0.0d0 )  WRITE ( kunit ,110)             'efg_DS(comm only)       ', efgtimetot3
      if ( fcoultimetot1 .ne. 0.0d0 )  WRITE ( kunit ,110)             'multipole_DS            ', fcoultimetot1
      if ( fcoultimetot3 .ne. 0.0d0 )  WRITE ( kunit ,110)             'multipole_DS(comm only) ', fcoultimetot3
    endif
    if ( longrange .eq. 'ewald')   then
      if ( estime        .ne. 0.0d0 )  WRITE ( kunit ,'(a)')           '======================='
      if ( estime        .ne. 0.0d0 )  WRITE ( kunit ,'(a)')           'EFG:'
      if ( estime        .ne. 0.0d0 )  WRITE ( kunit ,'(a)')           '======================='
      if ( strftimetot   .ne. 0.0d0 )  WRITE ( kunit ,110)             'struct fact             ', strftimetot
      if ( efgtimetot1   .ne. 0.0d0 )  WRITE ( kunit ,110)             'efg_ES(real  part)      ', efgtimetot1
      if ( efgtimetot2   .ne. 0.0d0 )  WRITE ( kunit ,110)             'efg_ES(recip part)      ', efgtimetot2 
      if ( efgtimetot3   .ne. 0.0d0 )  WRITE ( kunit ,110)             'efg_ES(comm  only)      ', efgtimetot3
      if ( fcoultime_es  .ne. 0.0d0 )  WRITE ( kunit ,'(a)')           '======================='
      if ( fcoultime_es  .ne. 0.0d0 )  WRITE ( kunit ,'(a)')           'Coulomb:'
      if ( fcoultime_es  .ne. 0.0d0 )  WRITE ( kunit ,'(a)')           '======================='
      if ( fcoultimetot1 .ne. 0.0d0 )  WRITE ( kunit ,110)             'multipole_ES(real  part)', fcoultimetot1
      if ( fcoultimetot2 .ne. 0.0d0 )  WRITE ( kunit ,110)             'multipole_ES(recip part)', fcoultimetot2 
      if ( fcoultimetot3 .ne. 0.0d0 )  WRITE ( kunit ,110)             'multipole_ES(comm only) ', fcoultimetot3
    endif

    WRITE ( kunit ,'(a)') '============================================================='
  endif

110   FORMAT(2X,A30,' :  cpu time',F9.2)


END SUBROUTINE print_time_info

END MODULE time
! ===== fmV =====
