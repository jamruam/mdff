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
MODULE vacf 

! ==============================================================================================================
! note the commented statements are related to the convential scheme to compute mean square displacement 
! maybe add a tag in vacftag for that ... but it is really less efficient than the n-order scheme in msd.f90 
! removed 27/06/11
! ===============================================================================================================
! well I do not understand what the code is doing ;)
! tmax , t0max, tdifmax ????????????????? 
! ===============================================================================================================

  implicit none

  integer, PARAMETER  :: tmax=10000, t0max=10000 
  integer             :: it0  
  double precision    :: tdifmax

  integer             :: tvacf , t0 
  double precision    :: dtime 

  integer , dimension (:) , allocatable :: ttv0                           ! t0max
  integer , dimension (:) , allocatable :: nvacf                          ! tmax
  double precision , dimension (:)   , allocatable :: vacff               ! tmax 
  double precision , dimension (:,:) , allocatable :: vxt0 , vyt0 , vzt0  ! natm , t0max

CONTAINS


!*********************** SUBROUTINE vacf_init **********************************
!
! init vacf calc
!
!******************************************************************************

SUBROUTINE vacf_init

  USE prop,     ONLY :  lvacf
  USE io_file,  ONLY :  stdin , stdout, kunit_OUTFF

  implicit none

  ! local
  character * 132 :: filename

  namelist /vacftag/ tdifmax , it0

  if ( .not. lvacf ) return

  CALL vacf_default_tag
  
  ! ============================
  ! reads vacftag tags
  ! ============================
  CALL getarg (1,filename)
  OPEN ( stdin , file = filename)
  READ ( stdin , vacftag)
  CLOSE ( stdin )

  CALL vacf_check_tag

  CALL vacf_print_info(stdout)
  CALL vacf_print_info(kunit_OUTFF)

  return 
 
END SUBROUTINE vacf_init


!*********************** SUBROUTINE vacf_default_tag ***************************
!
! set default values to vacf tag
!
!******************************************************************************

SUBROUTINE vacf_default_tag

  implicit none

  tdifmax = 100.0d0
  it0     = 1

  return 
 
END SUBROUTINE vacf_default_tag

!*********************** SUBROUTINE vacf_check_tag *****************************
!
! check vacf tag values
!
!******************************************************************************

SUBROUTINE vacf_check_tag

  implicit none

  return 
 
END SUBROUTINE vacf_check_tag


!*********************** SUBROUTINE vacf_print_info ****************************
!
! print info to outputs (stdout / OUTFF)
!
!******************************************************************************

SUBROUTINE vacf_print_info(kunit)

  USE io_file,  ONLY :  ionode 

  implicit none

  ! local
  integer :: kunit

  if( ionode  ) then
                               WRITE ( kunit ,'(a)')           ''
                               WRITE ( kunit ,'(a)')           'velocity auto-correlation function   : '
                               WRITE ( kunit ,'(a,f10.4)')     'tdifmax                              = ',tdifmax
                               WRITE ( kunit ,'(a,i5)')        'it0                                  = ',it0
                               WRITE ( kunit ,'(a)')           'output file                          : VACFFF'
  endif

  return 
 
END SUBROUTINE vacf_print_info

!*********************** SUBROUTINE vacf_alloc *********************************
! 
! allocate / deallocate principal arrays for vacf
! 
!******************************************************************************

SUBROUTINE vacf_alloc

  USE md,       ONLY :  dt
  USE prop,     ONLY :  nprop , lvacf
  USE config,   ONLY :  natm

  implicit none

  ! local
  integer :: i

  if ( .not. lvacf ) return

  allocate (    vxt0 ( natm , t0max ) , vyt0 ( natm , t0max ) , vzt0 ( natm , t0max ) )
  allocate (    ttv0 ( t0max ) )                    
  allocate (   vacff ( tmax  ) )
  allocate (   nvacf ( tmax  ) )

  t0 = 0
  tvacf = 0
  dtime = nprop * dt 
  do i = 1, tmax
    nvacf(i) = 0
    vacff(i) = 0
  enddo 

  return 
 
END SUBROUTINE vacf_alloc

SUBROUTINE vacf_dealloc

  USE prop,     ONLY :  lvacf

  implicit none

  if ( .not. lvacf ) return 
 
  deallocate (   vacff ) 
  deallocate (   nvacf )
  deallocate (    ttv0 )                    
  deallocate (    vxt0 , vyt0 , vzt0 )

  return 
 
END SUBROUTINE vacf_dealloc


!*********************** SUBROUTINE vacf_main *********************************
! adapted from Frenkel and Smit :
!
!        velocity autocorrelation function
!
!******************************************************************************

SUBROUTINE vacf_main 

  USE control,          ONLY :  cutoff
  USE config,           ONLY :  natm , rx , ry , rz , vx , vy , vz , box , rho 
  USE md,               ONLY :  dt
  USE io_file,          ONLY :  kunit_VACFFF
  USE constants,        ONLY :  pi
  USE time,             ONLY :  vacftimetot

  implicit none 
  
  INCLUDE 'mpif.h'

  ! local
  integer :: i, dtt , t , ttel , ierr
  ! timeinfo
  double precision :: ttt1 , ttt2 

  ttt1 = MPI_WTIME(ierr)

  tvacf = tvacf + 1
  ! ===============================================
  ! sample velocity auto correlation function and
  ! mean square displacement
  ! ===============================================
  if ( mod ( tvacf , it0 ) .eq. 0 ) then 
  ! ============
  !  new t=0
  ! ============
    t0 = t0 + 1
    ttel = mod ( t0 - 1 , t0max ) + 1
    ttv0 ( ttel ) = tvacf
    do i = 1 , natm 
      vxt0 ( i , ttel ) = vx ( i )
      vyt0 ( i , ttel ) = vy ( i )
      vzt0 ( i , ttel ) = vz ( i )
    enddo 
  endif 
  do t = 1, MIN ( t0 , t0max )
    dtt = tvacf - ttv0( t ) + 1
    if ( dtt .lt. tmax .and. dtt * dtime .le. tdifmax) then 
      nvacf ( dtt ) = nvacf ( dtt ) + 1
      do i = 1 , natm 
        vacff ( dtt ) = vacff ( dtt ) + vx ( i ) * vxt0 ( i , t ) + vy ( i ) * vyt0 ( i , t ) + vz ( i ) * vzt0 ( i , t )
      enddo 
    endif 
  enddo 
 
  ttt2 = MPI_WTIME(ierr)
  vacftimetot = vacftimetot + ( ttt2 - ttt1 )

  return

END SUBROUTINE vacf_main

!*********************** SUBROUTINE vacf_write ********************************
!
! write vacf and fourier transform = > dos
!
!******************************************************************************

SUBROUTINE vacf_write_output

  USE config,   ONLY :  natm ,box  
  USE io_file,  ONLY :  ionode , kunit_VACFFF , kunit_OUTFF 
  USE time,     ONLY :  vacftimetot2

  implicit none
  INCLUDE 'mpif.h'

  ! local
  integer :: ihbmax , i , ierr
  double precision :: vtime , dif
  double precision :: thmax
  double precision :: tauc , tau0 , errvacf 
  double complex   ,dimension (:), allocatable :: in , out 
  double precision ,dimension (:), allocatable :: rout
  double precision :: ttt1 , ttt2

  ttt1 = MPI_WTIME(ierr)
  
  allocate ( in ( tmax ) , out ( tmax ) , rout ( tmax ) ) 



  OPEN ( UNIT = kunit_VACFFF , FILE = 'VACFFF' ) 

  ! =======================================
  ! velocity auto-correlation function
  ! =======================================
  if ( tvacf .ne. 0 ) then 
    ! =======================================
    ! correlation time (for error estimate)
    ! =======================================
    tauc = 0
    do i = 1 , tmax
      if ( nvacf ( i ) * natm .ne. 0 ) then 
        tauc = tauc + ( vacff ( i ) / ( natm * nvacf ( i ) ) ) ** 2 * dtime
      endif
    enddo 
    ! =======================================
    ! normalisation
    ! =======================================
    if ( nvacf ( 1 ) * natm .ne. 0 ) tauc = tauc / ( vacff ( 1 ) / ( natm * nvacf ( 1 ) ) ** 2 )
    dif = 0
    ! =======================================
    ! total averaging time:
    ! =======================================
    tau0 = dtime*it0*t0
    errvacf = DSQRT( 2.0d0 * tauc * ( vacff ( 1 ) / ( natm * nvacf ( 1 ) ) ** 2 / tau0 ) )
    thmax = 0
    ihbmax = 0

    ! =======================
    !  normalisation
    ! =======================
    do i = 1 , tmax
     if ( nvacf ( i ) .ne. 0 ) then
       vacff ( i ) = vacff ( i ) / ( natm * nvacf ( i ) )
     endif
    enddo
    ! ========
    !   FFT
    ! ========
    in  = vacff
    CALL fft_1D_complex ( in , out , tmax )
    rout = dble ( out ) 


    do i = 1 , tmax
      vtime = dtime * ( i - 1 )
      if ( nvacf ( i ) .ne. 0 ) then 
        dif = dif + vacff ( i ) * dtime
        WRITE ( kunit_VACFFF , '(4e16.8)' ) vtime, vacff ( i ) , 1.0d0 / vtime , rout ( i ) 
        if ( vtime .gt. thmax ) then 
          ihbmax = nvacf ( i )
          thmax = vtime
        endif 
      endif 
    enddo 

    if ( ionode ) then
      WRITE ( kunit_OUTFF , 99002) tauc
      WRITE ( kunit_OUTFF , 99004 ) tvacf , t0 , dtime , dtime * it0 , dif/3.0d0
      WRITE ( kunit_OUTFF , '(a)' ) 'Diffusion calculated with conventional scheme '
      WRITE ( kunit_OUTFF , 99001 ) 2 * dtime , nvacf(3) , thmax , ihbmax 
    endif
  endif 

  CLOSE (kunit_VACFFF)

  deallocate ( in , out , rout )

  ttt2 = MPI_WTIME(ierr)
  vacftimetot2 = vacftimetot2 + ( ttt2 - ttt1 )

  RETURN

99001 FORMAT ('Number of samples for tmin = ', f8.3, ' is : ', i10, /, 'Number of samples for tmax = ', f8.3, ' is : ', i10)
99002 FORMAT ('Decorrelation time ', e12.4)
99004 FORMAT ('Velocity auto correlation function and mean square dis.:',&
             /, '   Number of samples       : ', i8, /, &
                '   Number of t=0           : ', i8, /, &
                '   Timestep between samples: ', f8.3, /, &
                '   Timestep between t=0    : ', f8.5, /, &
                '   Diffusion coef.         : ', f8.5)
 
 
END SUBROUTINE vacf_write_output


END MODULE vacf 
! ===== fmV =====