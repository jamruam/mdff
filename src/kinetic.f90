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

! *********************** SUBROUTINE init_velocities ***************************
!> \brief 
!!  This routine initialize the velocities. 
! ******************************************************************************
SUBROUTINE init_velocities

  USE constants,        ONLY :  dp 
  USE config,           ONLY :  vx , vy , vz , natm , ntype , ntypemax , atypei , center_of_mass
  USE md,               ONLY :  nequil , setvel , temp
  USE io_file,          ONLY :  ionode , stdout
  USE control,          ONLY :  lrestart

  implicit none

  ! local
  real(kind=dp) :: T, ekin 
  real(kind=dp) :: com ( 0:ntypemax , 3 )
  integer :: key , ia , it 

  ! =======================================
  !  set key: 
  !   key = 0 if all velocities  are null
  !   key = 1 if at least one is not null
  ! =======================================
  key = 0
  do ia = 1 , natm
    if ( vx (ia) .ne. 0.0_dp ) key = 1
    if ( vy (ia) .ne. 0.0_dp ) key = 1
    if ( vz (ia) .ne. 0.0_dp ) key = 1
  enddo

  ! ===============================================
  !  generate velocities from a given distribution
  ! ===============================================
  if ( (key .eq. 0 .or. .not. lrestart) .and. temp .ne. 0.0_dp .and. (nequil.ne.0) ) then

    separator(stdout)    
    io_node blankline(stdout)    
    io_node WRITE ( stdout ,'(a)') 'generate velocities'

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
  !  lrestart
  ! ===========
  elseif ( key .eq. 1 .and. lrestart ) then

    ! =======================
    !  input temperature  
    ! =======================
    CALL calc_temp(T, ekin)

    if ( ionode )  WRITE ( stdout ,'(a,f10.4)') 'input temperature                    = ',T

  else 

    ! ============================
    !  no initial kinetic energy
    ! ============================
    vx = 0.0_dp
    vy = 0.0_dp
    vz = 0.0_dp

    io_node  WRITE ( stdout , '(a)' ) 'no initial kinetic energy' 

  endif

  CALL center_of_mass ( vx , vy , vz , com )
  io_node WRITE ( stdout ,'(a,4e16.6)') 'center of mass vel. ALL ',com( 0 , :)
  do it = 1 , ntype
    io_node WRITE ( stdout ,'(a,a,a,4e16.6)') 'center of mass vel. ', atypei( it ),' ',com( it , :)
  enddo
  io_node blankline(stdout)

  return

END SUBROUTINE init_velocities


! *********************** SUBROUTINE rescale_velocities ************************
!> \brief
!! this subroutine rescale velocities with beredsen thermostat.
!! If tauberendsen = dt , this becomes a simple rescale procedure
!> \param[in] quite make the subroutine quite
! ******************************************************************************
SUBROUTINE rescale_velocities (quite)

  USE constants,                ONLY :  dp 
  USE config,                   ONLY :  natm , vx , vy , vz
  USE md,                       ONLY :  dt , temp , tauberendsen
  USE io_file,                  ONLY :  ionode , stdout

  implicit none

  ! global
  integer, intent(in) :: quite

  ! local
  integer :: ia
  real(kind=dp) :: T, lambda, ekin

  CALL calc_temp(T,ekin)

  lambda = ( 1.0_dp + (dt / tauberendsen) * (  (temp / T) - 1.0_dp) ) ** 0.5_dp

  do ia = 1 , natm
     vx ( ia ) = vx ( ia ) * lambda
     vy ( ia ) = vy ( ia ) * lambda
     vz ( ia ) = vz ( ia ) * lambda
  enddo

  if ( ionode .and. quite .eq. 1) then
    WRITE ( stdout ,'(a,f10.4)') 'Berendsen thermostat'
    WRITE ( stdout ,'(a,f10.4)') 'effective temperature      T        = ',T
    WRITE ( stdout ,'(a,f10.4)') 'wanted temperature         T0       = ',temp
    WRITE ( stdout ,'(a,f10.4)') 'velocities rescaled by              = ',lambda
#ifdef debug    
    CALL calc_temp(T,ekin)
    WRITE ( stdout ,'(a,f10.4)') 'debug : rescaled temperature       = ',T
#endif   
    blankline(stdout) 
    
  endif


  return

END SUBROUTINE rescale_velocities

! *********************** SUBROUTINE andersen_velocities ***********************
!!> \brief
!! andersen thermostat 
!> \note
!! based on algorithm 15 by Frenkel and Smit
! ******************************************************************************
SUBROUTINE andersen_velocities

  USE constants,                ONLY :  dp 
  USE config,                   ONLY :  natm , vx , vy , vz
  USE md,                       ONLY :  dt , temp , nuandersen 
  USE control,                  ONLY :  myrank , dgauss

  implicit none

  ! local
  integer :: ia , iseed
  real(kind=dp) :: T, ekin, sigma, U, G
  
  
  CALL calc_temp(T,ekin)
  sigma = dsqrt( temp ) 
 
  do ia = 1, natm
  CALL RANDOM_SEED(SIZE = ISEED)
  CALL RANDOM_NUMBER(HARVEST = U)
     if ( U .lt. nuandersen * dt) then  
       if ( dgauss .eq. 'boxmuller_basic' ) then
         CALL boxmuller_basic(G,0.0_dp,sigma)
         vx ( ia ) = G
         CALL boxmuller_basic(G,0.0_dp,sigma)
         vy ( ia ) = G
         CALL boxmuller_basic(G,0.0_dp,sigma)
         vz ( ia ) = G
       elseif ( dgauss .eq. 'boxmuller_polar' ) then
         CALL boxmuller_polar(G,0.0_dp,sigma)
         vx ( ia ) = G
         CALL boxmuller_polar(G,0.0_dp,sigma)
         vy ( ia ) = G
         CALL boxmuller_polar(G,0.0_dp,sigma)
         vz ( ia ) = G
       elseif ( dgauss .eq. 'knuth' ) then
         CALL knuth(G,0.0_dp,sigma)
         vx ( ia ) = G
         CALL knuth(G,0.0_dp,sigma)
         vy ( ia ) = G
         CALL knuth(G,0.0_dp,sigma)
         vz ( ia ) = G
       endif
     endif
  enddo

  return

END SUBROUTINE andersen_velocities


! *********************** SUBROUTINE uniform_random_velocities *****************
!> \brief
!! uniform random velocities distribution 
!> \author
!! Frenkel and Smit
!> \note
!! only used for test
! ******************************************************************************
SUBROUTINE uniform_random_velocities

  USE constants,        ONLY :  dp 
  USE config,   ONLY :  natm , vx , vy , vz
  USE md,       ONLY :  dt , temp  
  USE control,  ONLY :  dgauss
  USE io_file,  ONLY :  ionode , stdout

  implicit none
  INCLUDE "mpif.h"

  ! local
  integer :: i
  real(kind=dp) :: G
  integer :: iseed
  DOUBLE PRECISION v2, vx0, vy0, vz0, Vx0t, Vy0t, Vz0t, f

  if ( ionode ) then
    WRITE ( stdout ,'(a)') 'Velocities from uniform distribution' 
    WRITE ( stdout ,'(a)') 'routine from Frenkel Smit Case Study 4'
    WRITE ( stdout ,'(a)') 'only used to test the code'
  endif

  ! ===========================
  !  give particle a velocity
  ! ===========================
  vx0 = 0.0_dp
  vy0 = 0.0_dp
  vz0 = 0.0_dp
  v2 = 0.0_dp
  DO i = 1, natm
    CALL RANDOM_SEED(SIZE = iseed)
    CALL RANDOM_NUMBER(HARVEST = G)
    VX(i) = G - 0.5_dp
    CALL RANDOM_SEED(SIZE = iseed)
    CALL RANDOM_NUMBER(HARVEST = G)
    VY(i) = G - 0.5_dp
    CALL RANDOM_SEED(SIZE = iseed)
    CALL RANDOM_NUMBER(HARVEST = G)
    VZ(i) = G - 0.5_dp
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
  Vx0t = 0.0_dp
  Vy0t = 0.0_dp
  Vz0t = 0.0_dp
  f = SQRT ( 3 * DBLE ( natm ) * temp / v2 )
  v2 = 0.0_dp
  DO i = 1, natm
    VX(i) = (VX(i)-vx0) * f
    VY(i) = (VY(i)-vy0) * f
    VZ(i) = (VZ(i)-vz0) * f
    Vx0t = Vx0t + VX(i)
    Vy0t = Vy0t + VY(i)
    Vz0t = Vz0t + VZ(i)
    v2 = v2 + VX(i) ** 2 + VY(i) ** 2 + VZ(i) ** 2
  ENDDO
  v2 = v2 / DBLE(3 * natm)
  Vx0t = Vx0t/natm
  Vy0t = Vy0t/natm
  Vz0t = Vz0t/natm
  Temp = v2
  io_node WRITE ( stdout , 99001) v2
  io_node WRITE ( stdout , 99002) Vx0t, Vy0t, Vz0t

  return

99001 FORMAT('Initial temperature     : ',f6.3)
99002 FORMAT('Velocity centre of mass : ',/,'          x = ',e9.2,/, '          y = ',e9.2,/, '          z = ',e9.2)

END SUBROUTINE uniform_random_velocities

! *********************** SUBROUTINE maxwellboltzmann_velocities ***************
!
!> \brief
!! maxwell-boltzmann velocities
!! \f$ f(v) = \sqrt{ \frac{m}{2\pi kT}} \exp{-\frac{mv^2}{2kT}}\f$
!
!> \author
!! Allen-Tildsley
!
!> \note
!! adapted from F.24 
!
! ******************************************************************************
SUBROUTINE maxwellboltzmann_velocities

  USE constants,        ONLY :  dp 
  USE config,   ONLY :  natm , vx , vy , vz
  USE md,       ONLY :  dt , temp  
  USE control,  ONLY :  dgauss
  USE io_file,  ONLY :  ionode , stdout 

  implicit none
  INCLUDE "mpif.h"

  ! local
  integer :: ia
  real(kind=dp) :: RTEMP, SUMX, SUMY, SUMZ
  real(kind=dp) :: G

#ifdef fun
  if ( ionode ) then
    WRITE ( stdout ,'(a)') 'Heat:Hot, as _____:Cold'
    blankline(stdout) 
    WRITE ( stdout ,'(a)') 'a poem by Roald Hoffman'
    WRITE ( stdout ,'(a)') 'from: Chemistry Imagined, Reflections on Science'
    blankline(stdout) 
    WRITE ( stdout ,'(a)') 'Deep in,'
    WRITE ( stdout ,'(a)') "they're there, they're"
    WRITE ( stdout ,'(a)') "at it all the time, it's jai"
    WRITE ( stdout ,'(a)') 'alai on the hot molecular fronton-'
    WRITE ( stdout ,'(a)') 'a bounce off walls onto the packed aleatory'
    WRITE ( stdout ,'(a)') 'dance floor where sideswipes are medium of exchange,'
    WRITE ( stdout ,'(a)') 'momentum trades sealed in swift carom sequences,'
    WRITE ( stdout ,'(a)') 'or just that quick kick in the rear, the haphaz-'
    WRITE ( stdout ,'(a)') 'ard locomotion of the warm, warm world.'
    WRITE ( stdout ,'(a)') 'But spring nights grow cold in Ithaca;'
    WRITE ( stdout ,'(a)') 'the containing walls, glass or metal,'
    WRITE ( stdout ,'(a)') 'are a jagged rough rut of tethered'
    WRITE ( stdout ,'(a)') 'masses, still vibrant, but now'
    WRITE ( stdout ,'(a)') 'retarding, in each collision,'
    WRITE ( stdout ,'(a)') 'the cooling molecules.'
    WRITE ( stdout ,'(a)') "There, they're there,"
    WRITE ( stdout ,'(a)') 'still there,'
    WRITE ( stdout ,'(a)') 'in deep,'
    WRITE ( stdout ,'(a)') 'slow'
    WRITE ( stdout ,'(a)') '.'
    blankline(stdout) 
  endif
#endif

  if ( ionode ) then
    WRITE ( stdout      ,'(a)')     'Velocities from Maxwell-Boltzmann distribution'
    WRITE ( stdout      ,'(a,a20)') 'normal distribution method = ', dgauss  
  endif      

  RTEMP = SQRT ( temp )

  do ia = 1 , natm
     if ( dgauss .eq. 'boxmuller_basic' ) then
       CALL boxmuller_basic(G,0.0_dp,1.0_dp)
       vx ( ia ) = RTEMP * G
       CALL boxmuller_basic(G,0.0_dp,1.0_dp)
       vy ( ia ) = RTEMP * G
       CALL boxmuller_basic(G,0.0_dp,1.0_dp)
       vz ( ia ) = RTEMP * G
     elseif ( dgauss .eq. 'boxmuller_polar' ) then
       CALL boxmuller_polar(G,0.0_dp,1.0_dp)
       vx ( ia ) = RTEMP * G
       CALL boxmuller_polar(G,0.0_dp,1.0_dp)
       vy ( ia ) = RTEMP * G
       CALL boxmuller_polar(G,0.0_dp,1.0_dp)
       vz ( ia ) = RTEMP * G
     elseif ( dgauss .eq. 'knuth' ) then
       CALL knuth(G,0.0_dp,1.0_dp)
       vx ( ia ) = RTEMP * G 
       CALL knuth(G,0.0_dp,1.0_dp)
       vy ( ia ) = RTEMP * G
       CALL knuth(G,0.0_dp,1.0_dp)
       vz ( ia ) = RTEMP * G
     endif
  enddo

  SUMX = 0.0_dp
  SUMY = 0.0_dp
  SUMZ = 0.0_dp

  do ia = 1 , natm
    SUMX = SUMX + vx ( ia )
    SUMY = SUMY + vy ( ia )
    SUMZ = SUMZ + vz ( ia )
  enddo

  SUMX = SUMX / DBLE ( natm )
  SUMY = SUMY / DBLE ( natm )
  SUMZ = SUMZ / DBLE ( natm )

  do ia  = 1 , natm
     vx ( ia ) = vx ( ia ) - SUMX
     vy ( ia ) = vy ( ia ) - SUMY
     vz ( ia ) = vz ( ia ) - SUMZ
  enddo
  
  return

END SUBROUTINE maxwellboltzmann_velocities

! *********************** SUBROUTINE calc_temp *********************************
!> \brief
!! calculate temperature and kinetic enegy from velocities
!> \param[out] T temperature
!> \param[out] ekin kinetic energy
! ******************************************************************************
SUBROUTINE calc_temp (T, ekin)

  USE constants,        ONLY :  dp 
  USE config,   ONLY :  natm , vx , vy , vz 

  implicit none

  ! global
  real(kind=dp), intent(out) :: ekin , T

  ! local
  integer :: ia

  ekin = 0.0_dp
  do ia = 1 , natm
    ekin =  ekin + vx ( ia ) ** 2 + vy ( ia ) ** 2 + vz ( ia ) ** 2
  enddo
  ekin = ekin * 0.5_dp
  T = (2.0_dp/3.0_dp) * ekin
  T = T / DBLE (natm) 

  return

END SUBROUTINE calc_temp
! ===== fmV =====
