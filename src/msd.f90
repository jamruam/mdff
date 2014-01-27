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
!#define debug
! ======= Hardware =======

! *********************** MODULE msd *******************************************
!
!> \brief
!! module related mean square diplacement calculation
!
! ******************************************************************************
MODULE msd

  USE constants,  ONLY : dp
  USE mpimdff 
  implicit none

!  integer ,PARAMETER :: nblock=10, ibmax=20
  integer                                            :: nblock, ibmax

  real(kind=dp) , dimension (:,:,:) , allocatable ::  vxb , vyb , vzb
  real(kind=dp) , dimension (:,:) , allocatable   ::  avb

  integer :: iblm
  integer ,dimension (:,:) , allocatable :: telav
  integer ,dimension (:)   , allocatable :: ibl

  real(kind=dp) :: tdifmax , dtime
  logical :: lmsd

  integer :: nprop

CONTAINS

! *********************** SUBROUTINE msd_init **********************************
!
!> \brief
!! mean square diplacement initialisation
!
! ******************************************************************************
SUBROUTINE msd_init

  USE io,  ONLY :  stdin, stdout, ionode

  implicit none

  ! local
  integer            :: ioerr
  character(len=132) :: filename


  namelist /msdtag/  nblock  , &
                     ibmax   , &
                     tdifmax 

  if ( .not. lmsd ) return

  CALL msd_default_tag
  
  ! ==================
  !  read msdtag tags
  ! ==================
  CALL getarg (1,filename)
  OPEN ( stdin , file = filename)
  READ ( stdin , msdtag, iostat=ioerr)
  if ( ioerr .lt. 0 )  then
   io_node WRITE ( stdout, '(a)') 'ERROR reading input_file : msdtag section is absent'
   STOP
  elseif ( ioerr .gt. 0 )  then
   io_node WRITE ( stdout, '(a,i8)') 'ERROR reading input_file : msdtag wrong tag'
   STOP
  endif

  CLOSE ( stdin )

  CALL msd_check_tag

  CALL msd_print_info(stdout)

  return 
 
END SUBROUTINE msd_init


! *********************** SUBROUTINE msd_default_tag ***************************
!
!> \brief
!! set default values to msd tags
!
! ******************************************************************************
SUBROUTINE msd_default_tag

  implicit none

  tdifmax = 100.0_dp
  nblock  = 10
  ibmax   = 20

  return 
 
END SUBROUTINE msd_default_tag



! *********************** SUBROUTINE msd_check_tag *****************************
!
!> \brief
!! check msd tag values
!
! ******************************************************************************
SUBROUTINE msd_check_tag

  implicit none

  return 
 
END SUBROUTINE msd_check_tag


! *********************** SUBROUTINE msd_print_info ****************************
!
!> \brief
!! print msd information to standard output
!
! ******************************************************************************
SUBROUTINE msd_print_info(kunit)

  USE io,  ONLY :  ionode 

  implicit none
  
  !local
  integer :: kunit

  if ( ionode  ) then
    blankline(kunit)
    WRITE ( kunit ,'(a)')           'mean square displacement:'
    WRITE ( kunit ,'(a,f10.4)')     'tdifmax                              = ',tdifmax
    WRITE ( kunit ,'(a,i10)')       'nblock                               = ',nblock
    WRITE ( kunit ,'(a,i10)')       'ibmax                                = ',ibmax
    WRITE ( kunit ,'(a)')           'output file                          : MSDFF'
  endif

  return 
 
END SUBROUTINE msd_print_info


! *********************** SUBROUTINE msd_alloc *********************************
!
!> \brief
!! allocate and initialize variables
!
! ******************************************************************************
SUBROUTINE msd_alloc

  USE md,         ONLY :  dt
  USE config,     ONLY :  natm

  implicit none

  ! local
  integer :: ib , i , j

  if ( .not. lmsd ) return

  allocate ( vxb ( ibmax , nblock , natm ) , vyb ( ibmax , nblock , natm )  , vzb ( ibmax , nblock , natm )  ) 
  allocate ( avb ( ibmax , nblock) )
  allocate ( telav ( ibmax , nblock ) , ibl ( ibmax ) )

  dtime = nprop * dt
  DO ib = 1, ibmax
    ibl ( ib ) = 0
    DO j = 1, nblock
      telav ( ib , j ) = 0
      avb ( ib , j ) = 0.0_dp
      DO i = 1, natm 
        vxb ( ib , j , i ) = 0.0_dp
        vyb ( ib , j , i ) = 0.0_dp
        vzb ( ib , j , i ) = 0.0_dp
      ENDDO
    ENDDO
  ENDDO

  return 
 
END SUBROUTINE msd_alloc


! *********************** SUBROUTINE msd_dealloc *******************************
!
!> \brief
!! deallocate variables
!
! ******************************************************************************
SUBROUTINE msd_dealloc


  implicit none

  if ( .not. lmsd ) return

    deallocate ( vxb , vyb  , vzb  )
    deallocate ( avb )
    deallocate ( telav , ibl )

  return 
 
END SUBROUTINE msd_dealloc

! *********************** SUBROUTINE msd_main **********************************
!
!> \brief
!! Determine the mean square displacement using Algorithm 9
!
!> \note
!! Adapted from Frenkel and Smit
!
! ******************************************************************************
SUBROUTINE msd_main ( nmsd )
 
  USE config,     ONLY :  natm , vx , vy , vz
  USE md,         ONLY :  dt
  USE time,       ONLY :  msdtimetot
  USE io,    ONLY :  stdout       

  implicit none

  ! global
  integer :: nmsd

  ! local
  integer :: ia , iblock, ib, it,  inp, ii, inmax, ierr
  real(kind=dp) :: delx, dely, delz, xtime
#ifdef debug 
  real(kind=dp) :: r2asum
#endif
  ! timeinfo
  real(kind=dp) :: ttt1 , ttt2 

  ttt1 = MPI_WTIME(ierr)
 
  ! ===================================================
  !  determine current maximum number of blocks: iblm
  ! ===================================================
  iblm = 1
  ii = nmsd / nblock
  do while  (ii.ne.0)
    iblm = iblm + 1
    ii = ii / nblock
  ! ===========================================
  !  test maximu time not longer than tdifmax:
  ! ===========================================
    xtime = dtime * ( nblock ** ( iblm ) )
    if ( xtime .gt. tdifmax ) ii = 0
  enddo
  WRITE ( stdout ,'(2i5,3f8.3)') iblm , nmsd , dtime , xtime , tdifmax
  ! =====================================
  !  limit the maximum number of blocks
  ! =====================================
  iblm = MIN( iblm , ibmax )
  do ib = 1, iblm
    iblock = nblock ** ( ib - 1 )
    ! ==============================
    !  test for blocking operation
    ! ==============================
    if ( MOD ( nmsd , iblock ) .eq. 0 ) then
      ibl ( ib ) = ibl ( ib ) + 1
      ! ==============================
      !  limit to length n (=nblock)
      ! ==============================
      inmax = MIN( ibl ( ib ) , nblock )
#ifdef debug
     r2asum = 0
#endif
      do ia = 1 , natm 
        if (ib.eq.1) then
          ! ================================
          !  zero block: ordinary velocity
          ! ================================
          delx = vx ( ia )
          dely = vy ( ia )
          delz = vz ( ia )
        ELSE
          ! ===============================================
          !  (ib)th block: coarsed velocity previous block
          ! ===============================================
          delx = vxb ( ib-1 , 1 , ia )
          dely = vyb ( ib-1 , 1 , ia )
          delz = vzb ( ib-1 , 1 , ia )
        endif
        do it = 1, inmax
          inp = it
          if ( ibl ( ib ) .gt. nblock ) inp = it + 1
          if ( it .lt. inmax ) then
            vxb ( ib , it , ia ) = vxb ( ib , inp , ia ) + delx
            vyb ( ib , it , ia ) = vyb ( ib , inp , ia ) + dely
            vzb ( ib , it , ia ) = vzb ( ib , inp , ia ) + delz
          else
            vxb ( ib , it , ia ) = delx
            vyb ( ib , it , ia ) = dely
            vzb ( ib , it , ia ) = delz
          endif
        enddo
        do it = 1, inmax
          telav ( ib , it ) = telav ( ib , it ) + 1

          avb ( ib , it )   = avb ( ib , it )    + ( vxb ( ib , inmax - it + 1 , ia ) * dtime ) ** 2 &
                                                 + ( vyb ( ib , inmax - it + 1 , ia ) * dtime ) ** 2 &
                                                 + ( vzb ( ib , inmax - it + 1 , ia ) * dtime ) ** 2
#ifdef debug
          print*,'debug'
          if ( it .eq. 1 ) r2asum = r2asum       + ( vxb ( ib , inmax - it + 1 , ia ) * dtime ) ** 2 &
                                                 + ( vyb ( ib , inmax - it + 1 , ia ) * dtime ) ** 2 &
                                                 + ( vzb ( ib , inmax - it + 1 , ia ) * dtime ) ** 2
#endif
        enddo
      enddo
#ifdef debug
      r2asum = r2asum / natm
      ! ============================================================
      !  print mean square displacement to file for t=1,10,100,etc
      ! ============================================================
      WRITE ( 79 + ib , *) 'debug: ',dtime * ( nblock ** ( ib - 1 ) ), r2asum
#endif
    endif
  enddo
 
  ttt2 = MPI_WTIME(ierr)
  msdtimetot = msdtimetot + ( ttt2 - ttt1 )
 
  return

END SUBROUTINE msd_main


! *********************** SUBROUTINE msd_write_output **************************
!
!> \brief
!! write results to file MSDFF
!
! ******************************************************************************
SUBROUTINE msd_write_output ( quite ) 

  USE io,    ONLY :  ionode , kunit_MSDFF , stdout

  implicit none

  ! global
  integer :: quite

  ! local
  integer :: j, ib, ihbmax
  real(kind=dp) :: thmax

  if ( ionode ) then
    OPEN( kunit_MSDFF , file = 'MSDFF' ) 
    WRITE ( kunit_MSDFF, 99003 ) nblock , ibmax , tdifmax
    thmax = 0
    ihbmax = 0
    do ib = 1 , MIN( ibmax , iblm )
      do j = 2, MIN( ibl ( ib ), nblock)
        if (telav(ib,j).ne.0) then
          WRITE ( kunit_MSDFF , 99002 ) &
          j * dtime * ( nblock ** ( ib - 1 ) ), avb ( ib , j ) / telav ( ib , j ) , telav ( ib , j )
          if ( j * dtime * ( nblock ** ( ib - 1 ) ) .gt. thmax ) then
            ihbmax = telav ( ib , j )
            thmax  = j * dtime * ( nblock ** ( ib - 1 ) )
          endif
        endif
      enddo
    enddo
    CLOSE ( kunit_MSDFF ) 
  endif

  if ( quite .eq. 1 ) then
    if ( ionode ) then
      WRITE ( stdout , '(a)') 'Diffusion calculated with order-n scheme '
      WRITE ( stdout , 99001) 2*dtime, telav(1, 2), thmax, ihbmax
    endif 
  endif


  return 

99001 FORMAT ('Number of samples for tmin = ', f8.3, ' is : ', i10, /, & 
              'Number of samples for tmax = ', f8.3, ' is : ', i10)
99002 FORMAT (2(2x,e14.4), 2x, i10)
99003 FORMAT ('# nblock =', i6, ' ibmax = ', i6 , ' tdifmax = ', f10.4 )

END SUBROUTINE msd_write_output

END MODULE msd
! ===== fmV =====
