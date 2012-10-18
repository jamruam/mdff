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

MODULE prop

  implicit none

  logical :: lefg          ! calculate efg on the fly during the md   ( electric field gradient           )
  logical :: lgr           ! calculate gr on the fly during the md    ( radial distribution function      )
  logical :: lstrfac       ! calculate str_av on the fly during the md    ( static structural factor    )
  logical :: lmsd          ! calculate msd on the fly during the md   ( mean square displacment           ) 
  logical :: lvacf         ! calculate ivacf on the fly during the md ( velocity autocorrelation function ) 

  integer :: nprop         ! calculate properties each nprop step for lgr and efg
  integer :: nprop_start   ! calculate properties from step nprop_start
  integer :: nprop_print   ! print properties each nprop_print 
 

CONTAINS


!*********************** SUBROUTINE prop_init ***********************************
!
! init on-the-fly properties
!
!******************************************************************************

SUBROUTINE prop_init

  USE control,  ONLY :  cutoff
  USE io_file,  ONLY :  ionode , stdin, stdout, kunit_OUTFF
  
  implicit none

  ! local
  character * 132 :: filename
  integer :: ioerr       

  namelist /proptag/   lefg        , &
                       lgr         , &
                       lstrfac     , &
                       lmsd        , &
                       lvacf       , &
                       nprop_start , &
                       nprop       , &
                       nprop_print 
 
  ! ==================== 
  !  set default values
  ! ====================  
  CALL prop_default_tag
  ! =================  
  !  reads prop tags
  ! =================  
  CALL getarg (1,filename)
  OPEN ( stdin , file = filename)
  READ ( stdin , proptag, iostat=ioerr)
  if ( ioerr .lt. 0 )  then
    if ( ionode ) WRITE ( stdout, '(a)') 'ERROR reading input_file : proptag section is absent'
    STOP
  elseif ( ioerr .gt. 0 )  then
    if ( ionode ) WRITE ( stdout, '(a)') 'ERROR reading input_file : proptag wrong tag'
    STOP
  endif
  CLOSE ( stdin )
!  CALL prop_check_tag
  ! =================  
  !  print info
  ! =================  
  CALL prop_print_info(stdout)
  CALL prop_print_info(kunit_OUTFF)

  return

END SUBROUTINE prop_init



!*********************** SUBROUTINE prop_default_tag ***************************
!
! set default values to prop tags
!
!******************************************************************************

SUBROUTINE prop_default_tag

  implicit none

! default value
  lefg        = .false.
  lgr         = .false.
  lstrfac     = .false.
  lmsd        = .false.
  lvacf       = .false.
  nprop_start = 0
  nprop       = 1
  nprop_print = nprop

  return

END SUBROUTINE prop_default_tag


!*********************** SUBROUTINE prop_check_tag ***************************
!
! check prop tag values
!
!******************************************************************************

SUBROUTINE prop_check_tag

  implicit none

!something on nprop_print

  return

END SUBROUTINE prop_check_tag

!*********************** SUBROUTINE prop_print_info ***************************
!
! print info 
!
!******************************************************************************

SUBROUTINE prop_print_info(kunit)

  USE control,  ONLY :  lstatic
  USE io_file,  ONLY :  ionode

  implicit none

  ! local
  integer :: kunit
  
  if ( lefg .or. lgr .or. lmsd .or. lvacf ) then
    if ( ionode ) then
      WRITE ( kunit , '(a)'     )        '=============================================================' 
      WRITE ( kunit , '(a)'     )        ''
      WRITE ( kunit , '(a)'     )        'properties on-the-fly'
      WRITE ( kunit , '(a,i6,a)')        'calculate properties each ', nprop, ' steps '
      WRITE ( kunit , '(a,i6)'  )        'calculate properties from step ', nprop_start
      WRITE ( kunit , '(a,i6,a)')        'print properties each ', nprop_print, ' steps '
    endif
  endif

  return

END SUBROUTINE prop_print_info


END MODULE prop
! ===== fmV =====
