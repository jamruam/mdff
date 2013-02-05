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


!WARNING
!WARNING
!================================================================================
!
! this subroutines are not clear at all. need more documentation and references.
! example: leap-frog and verlet are not clearly distinguised
! to be rigorous I should test all of them and compare to existing codes
! 
!================================================================================
!WARNING


!*********************** SUBROUTINE prop_leap_frog ****************************
!
! leap-frog algorithm  ( or verlet algorithm ) 
!
!******************************************************************************

SUBROUTINE prop_leap_frog ( iastart , iaend )

  USE control,          ONLY :  lpbc
  USE md,               ONLY :  dt
  USE config,           ONLY :  natm , rx , ry , rz , fx , fy , fz , vx , vy , vz , rxs , rys , rzs
  USE thermodynamic,    ONLY :  temp_r , e_kin
  USE field 

  implicit none

  ! global
  integer, intent(inout) :: iastart , iaend !, list(250 * natm) , point(natm + 1)

  ! local
  integer :: ia 
  double precision, dimension (:), allocatable :: urx , ury , urz
  double precision :: idt , dtsq
  double precision :: tempi , kin
 
  allocate ( urx ( natm ) , ury ( natm ) , urz ( natm ) )
  urx = 0.0D0
  ury = 0.0D0
  urz = 0.0D0

  idt = 0.5D0 / dt 
  dtsq = dt * dt

 
  ! ==========================
  ! force + potential f(t)
  ! ==========================
  CALL engforce_driver ( iastart , iaend )

  ! ================================================= 
  !  r(t+dt) = 2 r(t) - r (t-dt) + f(t) dt*dt
  ! ================================================= 
  do ia = 1 , natm
    urx ( ia ) = 2.0D0 * rx ( ia ) - rxs ( ia ) + fx ( ia ) * dtsq
    ury ( ia ) = 2.0D0 * ry ( ia ) - rys ( ia ) + fy ( ia ) * dtsq
    urz ( ia ) = 2.0D0 * rz ( ia ) - rzs ( ia ) + fz ( ia ) * dtsq 
  enddo

  ! =====================================
  !  v(t) = ( r(t+dt) - r(t) ) / ( 2 dt) 
  ! =====================================
  do ia = 1, natm
    vx ( ia ) = idt * ( urx ( ia ) - rxs ( ia ) )
    vy ( ia ) = idt * ( ury ( ia ) - rys ( ia ) )
    vz ( ia ) = idt * ( urz ( ia ) - rzs ( ia ) )       
  enddo

  ! ==========================================================
  ! calculate kinetic energy of the full time step velocities
  ! ==========================================================
  CALL calc_temp( tempi , kin ) 
  temp_r = tempi
  e_kin  = kin      

  ! =========================================================
  !  updated positions r(t-dt) <= r(t)  and r(t) <= r (t+dt)
  ! =========================================================
  do ia = 1, natm
    rxs ( ia ) = rx  ( ia )
    rys ( ia ) = ry  ( ia )
    rzs ( ia ) = rz  ( ia )
    rx  ( ia ) = urx ( ia )
    ry  ( ia ) = ury ( ia )
    rz  ( ia ) = urz ( ia )
  enddo

  deallocate ( urx , ury , urz )

  return

END SUBROUTINE prop_leap_frog


!*********************** SUBROUTINE prop_velocity_verlet **********************
!
! propagation for velocity-verlet algotithm (found a paper)
!
!******************************************************************************

SUBROUTINE prop_velocity_verlet ( iastart , iaend )

  USE config,           ONLY :  natm, rx, ry, rz, vx, vy, vz, fx, fy, fz
  USE md,               ONLY :  dt
  USE control,          ONLY :  lpbc , lshiftpot
  USE thermodynamic,    ONLY :  temp_r , e_kin
  USE field

  implicit none

  ! global
  integer, intent(inout) :: iastart , iaend !, list(250 * natm) , point(natm + 1)

  ! local
  integer :: ia
  double precision :: dtsq2 , dt2
  double precision, dimension (:), allocatable :: fsx , fsy , fsz
  double precision :: tempi , kin


  allocate (fsx(natm),fsy(natm),fsz(natm))
  fsx = 0.0D0
  fsy = 0.0D0
  fsz = 0.0D0

  dtsq2 = dt * dt * 0.5d0
  dt2 = dt * 0.5d0

  ! =================================================
  !  r(t+dt) = r(t) + v(t)*dt + f(t) * dt*dt/2 
  !  store forces of the previous step  
  ! =================================================
  do ia = 1 , natm
    rx  ( ia ) = rx ( ia ) + vx ( ia ) * dt + fx ( ia ) * dtsq2
    ry  ( ia ) = ry ( ia ) + vy ( ia ) * dt + fy ( ia ) * dtsq2
    rz  ( ia ) = rz ( ia ) + vz ( ia ) * dt + fz ( ia ) * dtsq2
    fsx ( ia ) = fx ( ia )
    fsy ( ia ) = fy ( ia )
    fsz ( ia ) = fz ( ia )
  enddo

  ! ==========================
  ! force + potential f(t+dt)
  ! ==========================
  CALL engforce_driver ( iastart , iaend )

  ! ==============================================
  !  v(t+dt) = v(t) + ( f(t-dt) + f(t) ) * dt / 2
  ! ==============================================
  do ia = 1 , natm
    vx ( ia ) = vx ( ia ) + ( fsx ( ia ) + fx ( ia ) ) * dt2
    vy ( ia ) = vy ( ia ) + ( fsy ( ia ) + fy ( ia ) ) * dt2
    vz ( ia ) = vz ( ia ) + ( fsz ( ia ) + fz ( ia ) ) * dt2
  enddo

  CALL calc_temp(tempi, kin)
  temp_r = tempi      
  e_kin  = kin

  deallocate ( fsx , fsy , fsz )

  return

END SUBROUTINE prop_velocity_verlet

!*********************** SUBROUTINE nose_hoover_chain2 ************************
!
! Nose-Hoover chain (two) see Frenkel-Smit
!
!******************************************************************************

SUBROUTINE nose_hoover_chain2 ( iastart , iaend )

  USE config,           ONLY :  natm , rx , ry , rz , vx , vy , vz , fx , fy , fz 
  USE md,               ONLY :  dt, vxi1, vxi2, xi1, xi2
  USE thermodynamic,    ONLY :  temp_r , e_kin

  implicit none

  ! global
  integer, intent(inout) :: iastart , iaend 
  ! local
  double precision :: kin , tempi
 
  CALL calc_temp ( tempi , kin )

  CALL chain_nh_2 ( kin , vxi1 , vxi2 , xi1 , xi2 )

  CALL prop_pos_vel_verlet ( kin , iastart , iaend )

  CALL chain_nh_2( kin, vxi1, vxi2, xi1, xi2 ) 

  tempi = (2.0D0/3.0D0) * kin
  tempi = tempi/DBLE (natm)

  e_kin = kin
  temp_r = tempi
 
  return

END SUBROUTINE nose_hoover_chain2

!*********************** SUBROUTINE chain_nh_2 ********************************
!
! adapted from Frenkel and Smit
!
!******************************************************************************
SUBROUTINE chain_nh_2 ( kin, vxi1, vxi2, xi1, xi2)

  USE config,   ONLY :  natm , vx , vy , vz 
  USE md,       ONLY :  temp , dt , Qnosehoover

  implicit none

  ! global
  double precision, intent (inout) :: kin, vxi1, vxi2, xi1, xi2

  ! local
  integer :: ia
  double precision :: G1, G2, Q1, Q2  
  double precision :: dt2, dt4, dt8
  double precision :: s, L

  dt2 = dt  * 0.5d0
  dt4 = dt2 * 0.5d0
  dt8 = dt4 * 0.5D0

  Q1 = natm * Qnosehoover
  Q2 = Qnosehoover
  L = DBLE (3 * natm)

  G2   = ( Q1 * vxi1 * vxi1 - temp)
  vxi2 = vxi2 + G2 * dt4
  vxi1 = vxi1 * EXP ( - vxi2 * dt8 )
  G1 = ( 2.0D0 *  kin - L * temp) / Q1
  vxi1 = vxi1 + G1 * dt4
  vxi1 = vxi1 * EXP ( - vxi2 * dt8 )
  xi1 = xi1 + vxi1 * dt2
  xi2 = xi2 + vxi2 * dt2
  s = EXP ( - vxi1 * dt2 )

  do ia = 1, natm
    vx ( ia ) = s * vx ( ia )
    vy ( ia ) = s * vy ( ia )
    vz ( ia ) = s * vz ( ia )
  enddo
  
  kin = kin * s * s
  vxi1 = vxi1 * EXP ( - vxi2 * dt8 )
  G1 = ( 2.0D0 *  kin - L * temp) / Q1
  vxi1 = vxi1 + G1 * dt4
  vxi1 = vxi1 * EXP ( - vxi2 * dt8 )
  G2   = ( Q1 * vxi1 * vxi1 - temp) / Q2
  vxi2 = vxi2 + G2 * dt4

  return
 
END SUBROUTINE chain_nh_2


!*********************** SUBROUTINE prop_pos_vel_verlet ***********************
!
! propagates position and position in the velet algorithm
!
!******************************************************************************

SUBROUTINE prop_pos_vel_verlet ( kin , iastart , iaend )

  USE config,   ONLY :  natm , rx , ry , rz , ry , vx , vy , vz , fx , fy , fz 
  USE md,       ONLY :  dt
  USE control,  ONLY :  lpbc
  USE field 

  implicit none

  ! global
  double precision, intent (out) :: kin
  integer, intent(inout) :: iastart , iaend 

  ! local
  integer :: ia
  double precision :: dt2

  dt2 = dt * 0.5d0

  ! ================================
  !  r(t+dt) = r(t) + v(t) * dt / 2
  ! ================================
  do ia = 1 , natm  
    rx ( ia ) = rx ( ia ) + vx ( ia ) * dt2
    ry ( ia ) = ry ( ia ) + vy ( ia ) * dt2
    rz ( ia ) = rz ( ia ) + vz ( ia ) * dt2 
  enddo

  ! ==========================
  ! force + potential f(t+dt)
  ! ==========================
  CALL engforce_driver ( iastart , iaend )

  kin  = 0.0d0
  do ia = 1 , natm
    vx ( ia ) = vx ( ia ) + fx ( ia ) * dt
    vy ( ia ) = vy ( ia ) + fy ( ia ) * dt
    vz ( ia ) = vz ( ia ) + fz ( ia ) * dt
    rx ( ia ) = rx ( ia ) + vx ( ia ) * dt2
    ry ( ia ) = ry ( ia ) + vy ( ia ) * dt2
    rz ( ia ) = rz ( ia ) + vz ( ia ) * dt2
    kin = kin + vx ( ia ) * vx ( ia ) +  vy ( ia ) * vy ( ia ) + vz ( ia ) * vz ( ia )
  enddo      
  kin = kin * 0.5D0  


  return

END SUBROUTINE prop_pos_vel_verlet

!*********************** SUBROUTINE beeman ************************************
!
! D. Beeman "Some multistep methods for use in molecular dynamics calculations",
! Journal of Computational Physics 20 pp. 130-139 (1976)
!
!******************************************************************************
SUBROUTINE beeman ( iastart , iaend )

  USE config,           ONLY :  natm , rx , ry , rz , ry , vx , vy , vz , fx , fy , fz , fxs , fys , fzs 
  USE thermodynamic,    ONLY :  temp_r , e_kin
  USE md,               ONLY :  dt
  USE control,          ONLY :  lpbc
  USE field

  implicit none

  ! global
  integer, intent(inout) :: iastart , iaend 

  ! local
  integer :: ia
  double precision :: onesix , twothree , onethree , fivesix , dtsq
  double precision , dimension (:) , allocatable :: fxtmp , fytmp, fztmp
  double precision :: kin , tempi


  allocate ( fxtmp (natm) , fytmp (natm) , fztmp (natm) ) 

  ! ================
  !  some constants
  ! ================
  onesix = 1.0d0 / 6.0d0
  onethree = 2.0d0 * onesix
  twothree = 4.0d0 * onesix
  fivesix  = 5.0d0 * onesix
  dtsq = dt * dt
  
  ! ==================================================================
  ! r (t+dt) = r (t) + v (t) dt  + ( 2/3 f (t) - 1/6 f (t-dt) ) dt^2
  ! ==================================================================
  do ia = 1 , natm
    rx ( ia ) = rx ( ia ) + vx ( ia ) * dt + ( twothree * fx ( ia ) - onesix * fxs ( ia ) ) * dtsq 
    ry ( ia ) = ry ( ia ) + vy ( ia ) * dt + ( twothree * fy ( ia ) - onesix * fys ( ia ) ) * dtsq
    rz ( ia ) = rz ( ia ) + vz ( ia ) * dt + ( twothree * fz ( ia ) - onesix * fzs ( ia ) ) * dtsq
  enddo

  ! ========================
  ! save current force f(t)
  ! ========================
  fxtmp = fx
  fytmp = fy
  fztmp = fz

  ! ===========================
  ! force + potential  f(t+dt)
  ! ===========================
  CALL engforce_driver ( iastart , iaend )

  ! ==================================================================
  ! v (t+dt) = v (t) + ( 1/3 f(t+dt) + 5/6 f(t) - 1/6  f(t-dt) )  dt
  ! ==================================================================
  do ia = 1 , natm
    vx ( ia ) = vx ( ia ) + ( onethree * fx ( ia ) + fivesix * fxtmp ( ia ) - onesix * fxs ( ia ) ) * dt     
    vy ( ia ) = vy ( ia ) + ( onethree * fy ( ia ) + fivesix * fytmp ( ia ) - onesix * fys ( ia ) ) * dt
    vz ( ia ) = vz ( ia ) + ( onethree * fz ( ia ) + fivesix * fztmp ( ia ) - onesix * fzs ( ia ) ) * dt
  enddo

  CALL calc_temp(tempi, kin)
  temp_r = tempi      
  e_kin  = kin

  ! ==============
  ! store f(t-dt)
  ! ==============
  fxs = fxtmp
  fys = fytmp
  fzs = fztmp

  deallocate ( fxtmp , fytmp , fztmp ) 

  return

END SUBROUTINE beeman

! ===== fmV =====
