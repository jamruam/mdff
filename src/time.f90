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

  USE constants, ONLY : dp

  implicit none
  
  real(kind=dp) :: timetot        ! engforce SUBROUTINE       
  real(kind=dp) :: forcetimetot   ! engforce SUBROUTINE       
  real(kind=dp) :: fcoultimetot1  ! engforce_coul direct space      
  real(kind=dp) :: fcoultimetot2  ! engforce_coul fourier space      
  real(kind=dp) :: fcoultimetot3  ! engforce_coul comm only 
  real(kind=dp) :: efgtimetot1    ! EFG direct space           
  real(kind=dp) :: efgtimetot2    ! EFG fourier space        
  real(kind=dp) :: efgtimetot3    ! communication only         
  real(kind=dp) :: strftimetot    ! facteur de structure
  real(kind=dp) :: mdsteptimetot  !                            
  real(kind=dp) :: vacftimetot    !                            
  real(kind=dp) :: vacftimetot2   !                            
  real(kind=dp) :: vnlisttimetot  !                            
  real(kind=dp) :: msdtimetot     !                            
  real(kind=dp) :: grtimetot      !                            
  real(kind=dp) :: grtimetot_comm !                            
  real(kind=dp) :: opttimetot     !                          
  real(kind=dp) :: vibtimetot     !                           
  real(kind=dp) :: hessiantimetot !                           
  real(kind=dp) :: bandtimetot !                           
  real(kind=dp) :: doskpttimetot !                           
  real(kind=dp) :: fvibtimetot    !                           
  real(kind=dp) :: timetmp        ! use for opt               

CONTAINS


!*********************** SUBROUTINE time_init *********************************
!
! set time to zero
!
!******************************************************************************

SUBROUTINE time_init

  implicit none

  timetot       = 0.0_dp
  efgtimetot1   = 0.0_dp
  efgtimetot2   = 0.0_dp
  efgtimetot3   = 0.0_dp
  strftimetot   = 0.0_dp 
  mdsteptimetot = 0.0_dp
  vnlisttimetot = 0.0_dp
  msdtimetot    = 0.0_dp
  grtimetot     = 0.0_dp 
  grtimetot_comm= 0.0_dp 
  forcetimetot  = 0.0_dp
  fcoultimetot1 = 0.0_dp
  fcoultimetot2 = 0.0_dp
  fcoultimetot3 = 0.0_dp
  vibtimetot    = 0.0_dp
  hessiantimetot= 0.0_dp
  bandtimetot   = 0.0_dp
  doskpttimetot = 0.0_dp
  opttimetot    = 0.0_dp
  vacftimetot   = 0.0_dp                             
  vacftimetot2  = 0.0_dp                             

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
  real(kind=dp) :: dstime , estime
  real(kind=dp) :: fcoultime_es , fcoultime_ds

  dstime = efgtimetot1 + efgtimetot3               ! direct
  estime = strftimetot + efgtimetot1 + efgtimetot2 + efgtimetot3 ! ewald
  fcoultime_es = fcoultimetot1 + fcoultimetot2 + fcoultimetot3
  fcoultime_ds = fcoultimetot1 + fcoultimetot3 

  !timing information at the end
  if ( ionode ) then
                                        WRITE ( kunit ,110)             'TOTAL'  , timetot
    if ( mdsteptimetot .ne. 0.0_dp )    WRITE ( kunit ,110)             'MD'     , mdsteptimetot
    if ( opttimetot    .ne. 0.0_dp )    WRITE ( kunit ,110)             'OPT'    , opttimetot
    if ( vibtimetot    .ne. 0.0_dp )    WRITE ( kunit ,110)             'VIB'    , vibtimetot
    if ( fvibtimetot   .ne. 0.0_dp )    WRITE ( kunit ,110)             'FVIB'   , fvibtimetot
    if ( dstime .ne.0.0_dp .and. &
         efgtimetot2   .eq. 0.0_dp )    WRITE ( kunit ,110)             'EFG_DS' , dstime
    if ( efgtimetot2   .ne. 0.0_dp )    WRITE ( kunit ,110)             'EFG_ES' , estime
    if ( fcoultime_ds  .ne. 0.0_dp .and. &
         fcoultimetot2 .eq. 0.0_dp )    WRITE ( kunit ,110)             'FIELD_DS', fcoultime_ds
    if ( fcoultimetot2 .ne. 0.0_dp )    WRITE ( kunit ,110)             'FIELD_ES', fcoultime_es
    if ( grtimetot     .ne. 0.0_dp )    WRITE ( kunit ,110)             'GR'                , grtimetot 
    if ( grtimetot_comm.ne. 0.0_dp )    WRITE ( kunit ,110)             'GR(comm only)'     , grtimetot_comm 
    if ( msdtimetot    .ne. 0.0_dp )    WRITE ( kunit ,110)             'MSD'    , msdtimetot
    if ( vacftimetot   .ne. 0.0_dp )    WRITE ( kunit ,110)             'VACF'   , vacftimetot
    if ( vacftimetot2  .ne. 0.0_dp )    WRITE ( kunit ,110)             'VACF2'  , vacftimetot2
                                        WRITE ( kunit ,'(a)')           ''
                                        WRITE ( kunit ,'(a)')           '======================='
                                        WRITE ( kunit ,'(a)')           'main subroutines:'
                                        WRITE ( kunit ,'(a)')           '======================='
    if ( calc .eq. 'md' ) then
                                        WRITE ( kunit ,'(a)')           'MD:'
                                        WRITE ( kunit ,'(a)')           '======================='
      if ( forcetimetot  .ne. 0.0_dp )  WRITE ( kunit ,110)             'engforce_bmlj      ' , forcetimetot        
      if ( vnlisttimetot .ne. 0.0_dp )  WRITE ( kunit ,110)             'vnlistcheck        ' , vnlisttimetot
                                        WRITE ( kunit ,'(a)')           '======================='
    endif
    
    if ( longrange .eq. 'direct' )   then
      if ( dstime        .ne. 0.0_dp )  WRITE ( kunit ,'(a)')           'EFG:'
                                        WRITE ( kunit ,'(a)')           '======================='
      if ( efgtimetot1   .ne. 0.0_dp )  WRITE ( kunit ,110)             'efg_DS(efg  only)       ', efgtimetot1 
      if ( efgtimetot3   .ne. 0.0_dp )  WRITE ( kunit ,110)             'efg_DS(comm only)       ', efgtimetot3
      if ( fcoultimetot1 .ne. 0.0_dp )  WRITE ( kunit ,110)             'multipole_DS            ', fcoultimetot1
      if ( fcoultimetot3 .ne. 0.0_dp )  WRITE ( kunit ,110)             'multipole_DS(comm only) ', fcoultimetot3
    endif
    if ( longrange .eq. 'ewald')   then
      if ( estime        .ne. 0.0_dp )  WRITE ( kunit ,'(a)')           'EFG:'
      if ( estime        .ne. 0.0_dp )  WRITE ( kunit ,'(a)')           '======================='
      if ( strftimetot   .ne. 0.0_dp )  WRITE ( kunit ,110)             'struct fact             ', strftimetot
      if ( efgtimetot1   .ne. 0.0_dp )  WRITE ( kunit ,110)             'efg_ES(real  part)      ', efgtimetot1
      if ( efgtimetot2   .ne. 0.0_dp )  WRITE ( kunit ,110)             'efg_ES(recip part)      ', efgtimetot2 
      if ( efgtimetot3   .ne. 0.0_dp )  WRITE ( kunit ,110)             'efg_ES(comm  only)      ', efgtimetot3
      if ( fcoultime_es  .ne. 0.0_dp )  WRITE ( kunit ,'(a)')           '======================='
      if ( fcoultime_es  .ne. 0.0_dp )  WRITE ( kunit ,'(a)')           'Coulomb:'
      if ( fcoultime_es  .ne. 0.0_dp )  WRITE ( kunit ,'(a)')           '======================='
      if ( fcoultimetot1 .ne. 0.0_dp )  WRITE ( kunit ,110)             'multipole_ES(real  part)', fcoultimetot1
      if ( fcoultimetot2 .ne. 0.0_dp )  WRITE ( kunit ,110)             'multipole_ES(recip part)', fcoultimetot2 
      if ( fcoultimetot3 .ne. 0.0_dp )  WRITE ( kunit ,110)             'multipole_ES(comm only) ', fcoultimetot3
    endif
    if ( vibtimetot.ne.0.0_dp )         WRITE ( kunit ,'(a)')           '======================='
    if ( vibtimetot.ne.0.0_dp )         WRITE ( kunit ,'(a)')           'VIB:'
    if ( vibtimetot.ne.0.0_dp )         WRITE ( kunit ,'(a)')           '======================='
    if ( hessiantimetot.ne.0.0_dp )     WRITE ( kunit ,110)             'hessian                ', hessiantimetot
    if ( bandtimetot.ne.0.0_dp )        WRITE ( kunit ,110)             'band                   ', bandtimetot
    if ( doskpttimetot.ne.0.0_dp )      WRITE ( kunit ,110)             'doskpttimetot          ', doskpttimetot

    WRITE ( kunit ,'(a)') '============================================================='
  endif

110   FORMAT(2X,A30,' :  cpu time',F9.2)


END SUBROUTINE print_time_info

END MODULE time
! ===== fmV =====
