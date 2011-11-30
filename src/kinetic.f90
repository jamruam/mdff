! ===== fmV =====
!*********************** SUBROUTINE init_velocities ***************************
! 
! This routine initialize the velocities. 
!
!******************************************************************************

SUBROUTINE init_velocities

  USE config,   ONLY :  vx , vy , vz , natm 
  USE md,       ONLY :  nequil , setvel , temp
  USE io_file,  ONLY :  ionode , stdout, kunit_OUTFF
  USE control,  ONLY :  restart

  implicit none

  ! local
  double precision :: T, ekin 
  integer :: key , ia

  ! =======================================
  !  set key: 
  !   key = 0 if all velocities are null
  !   key = 1 if at least one is no null
  ! =======================================
  key = 0
  do ia = 1 , natm
    if ( vx (ia) .ne. 0.0d0 ) key = 1
    if ( vy (ia) .ne. 0.0d0 ) key = 1
    if ( vz (ia) .ne. 0.0d0 ) key = 1
  enddo

  
  ! ===============================================
  !  generate velocities from a given distribution
  ! ===============================================
  if ( (key .eq. 0 .or. .not. restart) .and. temp .ne. 0.0d0 .and. (nequil.ne.0)) then

    if( ionode ) then
      WRITE ( kunit_OUTFF ,'(a)') '============================================================='
      WRITE ( kunit_OUTFF ,'(a)') ''
      WRITE ( kunit_OUTFF ,'(a)') 'generate velocities'
      WRITE ( stdout      ,'(a)') '============================================================='
      WRITE ( stdout      ,'(a)') ''
      WRITE ( stdout      ,'(a)') 'generate velocities'
    endif

    ! ================================
    !  Maxwell-Boltzmann distribution
    ! ================================
    if (setvel .eq. 'MaxwBoltz') CALL maxwellboltzmann_velocities

    ! ============================================
    !  Uniform distribution for test purpose only
    ! ============================================
    if (setvel .eq. 'Uniform')   CALL uniform_random_velocities

    ! ====================
    !  rescale velocities 
    ! ====================
    if ( setvel.ne.'Uniform' )   CALL rescale_velocities(1)

  ! ===========
  !  restart
  ! ===========
  elseif ( key .eq. 1 .and. restart ) then

    ! =======================
    !  input temperature  
    ! =======================
    CALL calc_temp(T, ekin)

    if( ionode ) then
      WRITE ( stdout      ,'(a,f10.4)') 'input temperature                    = ',T
      WRITE ( kunit_OUTFF ,'(a,f10.4)') 'input temperature                    = ',T
    endif

  else 

    ! ============================
    !  no initial kinetic energy
    ! ============================
    vx = 0.0d0
    vy = 0.0d0
    vz = 0.0d0

    if( ionode ) then 
      WRITE( stdout      , '(a)' ) 'no initial kinetic energy' 
      WRITE( kunit_OUTFF , '(a)' ) 'no initial kinetic energy' 
    endif

  endif

  return

END SUBROUTINE init_velocities


!*********************** SUBROUTINE rescale_velocities ************************
!
! this subroutine rescale velocities with beredsen thermostat.
! If tauberendsen = dt , this becomes a simple rescale procedure
!
!******************************************************************************

SUBROUTINE rescale_velocities (quite)

  USE config,   ONLY :  natm , vx , vy , vz
  USE md,       ONLY :  dt , temp , tauberendsen
  USE io_file,  ONLY :  ionode , stdout, kunit_OUTFF

  implicit none

  ! global
  integer, intent(in) :: quite

  ! local
  integer :: i 
  double precision :: T, lambda, ekin

  CALL calc_temp(T,ekin)

  lambda = ( 1.0D0 + (dt / tauberendsen) * (  (temp / T) - 1.0D0) ) ** 0.5D0

  do i = 1, natm
     VX(I) = VX(I) * lambda
     VY(I) = VY(I) * lambda
     VZ(I) = VZ(I) * lambda
  enddo

  if( ionode .and. quite .eq. 1) then
    WRITE ( kunit_OUTFF ,'(a,f10.4)') 'Berendsen thermostat'
    WRITE ( kunit_OUTFF ,'(a,f10.4)') 'velocities rescaled by              = ',lambda
    WRITE ( kunit_OUTFF ,'(a,f10.4)') 'wanted temperature         T0       = ',temp
    WRITE ( kunit_OUTFF ,'(a,f10.4)') 'effective temperature      T        = ',T
    WRITE ( kunit_OUTFF ,'(a)')       ''
    WRITE ( stdout ,'(a,f10.4)') 'Berendsen thermostat'
    WRITE ( stdout ,'(a,f10.4)') 'velocities rescaled by              = ',lambda
    WRITE ( stdout ,'(a,f10.4)') 'wanted temperature         T0       = ',temp
    WRITE ( stdout ,'(a,f10.4)') 'effective temperature      T        = ',T
    WRITE ( stdout ,'(a)')       ''
  endif


  return

END SUBROUTINE rescale_velocities

!*********************** SUBROUTINE andersen_velocities ***********************
!
! andersen thermostat 
! based on algorithm 15 by Frenkel and Smit
!
!******************************************************************************

SUBROUTINE andersen_velocities

  USE config,   ONLY :  natm , vx , vy , vz
  USE md,       ONLY :  dt , temp , nuandersen 
  USE control,  ONLY :  myrank , dgauss

  implicit none

  ! local
  integer :: i, iseed
  double precision :: T, ekin, sigma, U, G
  
  
  CALL calc_temp(T,ekin)
  sigma = dsqrt( temp ) 
 
  do i = 1, natm
  CALL RANDOM_SEED(SIZE = ISEED)
  CALL RANDOM_NUMBER(HARVEST = U)
     if ( U .lt. nuandersen * dt) then  
       if ( dgauss .eq. 'boxmuller_basic' ) then
         CALL boxmuller_basic(G,0.0D0,sigma)
         vx(i) = G
         CALL boxmuller_basic(G,0.0D0,sigma)
         vy(i) = G
         CALL boxmuller_basic(G,0.0D0,sigma)
         vz(i) = G
       elseif ( dgauss .eq. 'boxmuller_polar' ) then
         CALL boxmuller_polar(G,0.0D0,sigma)
         vx(i) = G
         CALL boxmuller_polar(G,0.0D0,sigma)
         vy(i) = G
         CALL boxmuller_polar(G,0.0D0,sigma)
         vz(i) = G
       elseif ( dgauss .eq. 'knuth' ) then
         CALL knuth(G,0.0D0,sigma)
         vx(i) = G
         CALL knuth(G,0.0D0,sigma)
         vy(i) = G
         CALL knuth(G,0.0D0,sigma)
         vz(i) = G
       endif
     endif
  enddo


  return

END SUBROUTINE andersen_velocities


!*********************** SUBROUTINE uniform_random_velocities *****************
!
! only used for test
! Frenkel and Smit
!
!******************************************************************************

SUBROUTINE uniform_random_velocities


  USE config,   ONLY :  natm , vx , vy , vz
  USE md,       ONLY :  dt , temp  
  USE control,  ONLY :  dgauss
  USE io_file,  ONLY :  ionode , stdout, kunit_OUTFF

  implicit none
  INCLUDE "mpif.h"

  ! local
  integer :: i
  double precision :: G
  integer :: iseed
  DOUBLE PRECISION v2, vx0, vy0, vz0, Vx0t, Vy0t, Vz0t, f

  if ( ionode ) then
    WRITE ( kunit_OUTFF ,'(a)') 'Velocities from uniform distribution' 
    WRITE ( kunit_OUTFF ,'(a)') 'routine from Frenkel Smit Case Study 4'
    WRITE ( kunit_OUTFF ,'(a)') 'only   USEd to test the code'
  endif

  ! ===========================
  !  give particle a velocity
  ! ===========================
  vx0 = 0.D0
  vy0 = 0.D0
  vz0 = 0.D0
  v2 = 0.D0
  DO i = 1, natm
    CALL RANDOM_SEED(SIZE = iseed)
    CALL RANDOM_NUMBER(HARVEST = G)
    VX(i) = G - 0.5D0
    CALL RANDOM_SEED(SIZE = iseed)
    CALL RANDOM_NUMBER(HARVEST = G)
    VY(i) = G - 0.5D0
    CALL RANDOM_SEED(SIZE = iseed)
    CALL RANDOM_NUMBER(HARVEST = G)
    VZ(i) = G - 0.5D0
    vx0 = vx0 + VX(i)
    vy0 = vy0 + VY(i)
    vz0 = vz0 + VZ(i)
    v2 = v2 + VX(i) ** 2 + VY(i) ** 2 + VZ(i) ** 2
  ENDDO
  ! ======================================
  !   set centre of mass movement to zero
  ! ======================================
  vx0 = vx0/natm
  vy0 = vy0/natm
  vz0 = vz0/natm
  Vx0t = 0.D0
  Vy0t = 0.D0
  Vz0t = 0.D0
  f = SQRT(3 * natm * temp/v2)
  v2 = 0.D0
  DO i = 1, natm
    VX(i) = (VX(i)-vx0) * f
    VY(i) = (VY(i)-vy0) * f
    VZ(i) = (VZ(i)-vz0) * f
    Vx0t = Vx0t + VX(i)
    Vy0t = Vy0t + VY(i)
    Vz0t = Vz0t + VZ(i)
    v2 = v2 + VX(i) ** 2 + VY(i) ** 2 + VZ(i) ** 2
  ENDDO
  v2 = v2/DBLE(3 * natm)
  Vx0t = Vx0t/natm
  Vy0t = Vy0t/natm
  Vz0t = Vz0t/natm
  Temp = v2
  if( ionode ) WRITE ( stdout , 99001) v2
  if( ionode ) WRITE ( stdout , 99002) Vx0t, Vy0t, Vz0t

  return

99001 FORMAT ('Initial temperature     : ', f5.3)
99002 FORMAT ('Velocity centre of mass : ', /, '          x = ', e8.2, /, '          y = ', e8.2, /, '          z = ', e8.2)

END SUBROUTINE uniform_random_velocities

!*********************** SUBROUTINE maxwellboltzmann_velocities ***************
!
!  f (v) = sqrt( m / 2 pi kT) exp [-mv^2/2kT]
!  F.24 Allen-Tildsley
!
!******************************************************************************

SUBROUTINE maxwellboltzmann_velocities

  USE config,   ONLY :  natm , vx , vy , vz
  USE md,       ONLY :  dt , temp  
  USE control,  ONLY :  dgauss
  USE io_file,  ONLY :  ionode , stdout , kunit_OUTFF

  implicit none
  INCLUDE "mpif.h"

  ! local
  integer :: i
  double precision :: RTEMP, SUMX, SUMY, SUMZ
  double precision :: G

#ifdef fun
  if( ionode ) then
    WRITE ( kunit_OUTFF ,'(a)') 'Heat:Hot, as _____:Cold'
    WRITE ( kunit_OUTFF ,'(a)') ''
    WRITE ( kunit_OUTFF ,'(a)') 'a poem by Roald Hoffman'
    WRITE ( kunit_OUTFF ,'(a)') 'from: Chemistry Imagined, Reflections on Science'
    WRITE ( kunit_OUTFF ,'(a)') ''
    WRITE ( kunit_OUTFF ,'(a)') 'Deep in,'
    WRITE ( kunit_OUTFF ,'(a)') "they're there, they're"
    WRITE ( kunit_OUTFF ,'(a)') "at it all the time, it's jai"
    WRITE ( kunit_OUTFF ,'(a)') 'alai on the hot molecular fronton-'
    WRITE ( kunit_OUTFF ,'(a)') 'a bounce off walls onto the packed aleatory'
    WRITE ( kunit_OUTFF ,'(a)') 'dance floor where sideswipes are medium of exchange,'
    WRITE ( kunit_OUTFF ,'(a)') 'momentum trades sealed in swift carom sequences,'
    WRITE ( kunit_OUTFF ,'(a)') 'or just that quick kick in the rear, the haphaz-'
    WRITE ( kunit_OUTFF ,'(a)') 'ard locomotion of the warm, warm world.'
    WRITE ( kunit_OUTFF ,'(a)') 'But spring nights grow cold in Ithaca;'
    WRITE ( kunit_OUTFF ,'(a)') 'the containing walls, glass or metal,'
    WRITE ( kunit_OUTFF ,'(a)') 'are a jagged rough rut of tethered'
    WRITE ( kunit_OUTFF ,'(a)') 'masses, still vibrant, but now'
    WRITE ( kunit_OUTFF ,'(a)') 'retarding, in each collision,'
    WRITE ( kunit_OUTFF ,'(a)') 'the cooling molecules.'
    WRITE ( kunit_OUTFF ,'(a)') "There, they're there,"
    WRITE ( kunit_OUTFF ,'(a)') 'still there,'
    WRITE ( kunit_OUTFF ,'(a)') 'in deep,'
    WRITE ( kunit_OUTFF ,'(a)') 'slow'
    WRITE ( kunit_OUTFF ,'(a)') '.'
  endif
#endif

  if ( ionode ) then
    WRITE ( kunit_OUTFF ,'(a)')     'Velocities from Maxwell-Boltzmann distribution'
    WRITE ( kunit_OUTFF ,'(a,a20)') 'normal distribution method = ', dgauss  
    WRITE ( stdout      ,'(a)')     'Velocities from Maxwell-Boltzmann distribution'
    WRITE ( stdout      ,'(a,a20)') 'normal distribution method = ', dgauss  
  endif      

  RTEMP = DSQRT ( temp )

  do i = 1, natm
     if ( dgauss .eq. 'boxmuller_basic' ) then
       CALL boxmuller_basic(G,0.0D0,1.0D0)
       vx(i) = RTEMP * G
       CALL boxmuller_basic(G,0.0D0,1.0D0)
       vy(i) = RTEMP * G
       CALL boxmuller_basic(G,0.0D0,1.0D0)
       vz(i) = RTEMP * G
     elseif ( dgauss .eq. 'boxmuller_polar' ) then
       CALL boxmuller_polar(G,0.0D0,1.0D0)
       vx(i) = RTEMP * G
       CALL boxmuller_polar(G,0.0D0,1.0D0)
       vy(i) = RTEMP * G
       CALL boxmuller_polar(G,0.0D0,1.0D0)
       vz(i) = RTEMP * G
     elseif ( dgauss .eq. 'knuth' ) then
       CALL knuth(G,0.0D0,1.0D0)
       vx(i) = RTEMP * G 
       CALL knuth(G,0.0D0,1.0D0)
       vy(i) = RTEMP * G
       CALL knuth(G,0.0D0,1.0D0)
       vz(i) = RTEMP * G
     endif
  enddo

  SUMX = 0.0d0
  SUMY = 0.0d0
  SUMZ = 0.0d0

  do i = 1, natm
    SUMX = SUMX + vx(i)
    SUMY = SUMY + vy(i)
    SUMZ = SUMZ + vz(i)
  enddo

  SUMX = SUMX / dble( natm )
  SUMY = SUMY / dble( natm )
  SUMZ = SUMZ / dble( natm )

  do i = 1, natm
     vx(i) = vx(i) - SUMX
     vy(i) = vy(i) - SUMY
     vz(i) = vz(i) - SUMZ
  enddo
  
  return

END SUBROUTINE maxwellboltzmann_velocities

!*********************** SUBROUTINE calc_temp *********************************
!
! calculate temperature and kinetic enegy from velocities
!
!******************************************************************************

SUBROUTINE calc_temp (T, ekin)

  USE config,   ONLY :  natm , vx , vy , vz 

  implicit none

  ! global
  double precision, intent(out) :: ekin , T

  ! local
  integer :: i

  ekin = 0.D0
  do i = 1,natm
    ekin =  ekin + vx(i) ** 2 + vy(i) ** 2 + vz(i) ** 2
  enddo
  ekin = ekin * 0.5d0
  T = (2.0D0/3.0D0) * ekin
  T = T/dble(natm) 

  return

END SUBROUTINE calc_temp
! ===== fmV =====
