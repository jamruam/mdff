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
!#define debug_nvt_nhcn
! ======= Hardware =======

!================================================================================
!
!> \file
!! all routines related to dynamical integration of phase-space 
!> \brief
!! this subroutines are not clear at all. need more documentation and references.
!! example: leap-frog and verlet are not clearly distinguised
!! to be rigorous I should test all of them and compare to existing codes
! 
!================================================================================


! *********************** SUBROUTINE prop_leap_frog ****************************
!
!> \brief
! leap-frog algorithm  ( or verlet algorithm ) 
!
!> \todo
!! leap-frog and verlet are not clearly distinguished !!
!
! ******************************************************************************
SUBROUTINE prop_leap_frog 

  USE constants,                ONLY :  dp 
  USE control,                  ONLY :  lpbc
  USE md,                       ONLY :  dt
  USE config,                   ONLY :  natm , rx , ry , rz , fx , fy , fz , vx , vy , vz , rxs , rys , rzs
  USE thermodynamic,            ONLY :  temp_r , e_kin
  USE field,                    ONLY :  engforce_driver

  implicit none

  ! local
  integer                                   :: ia 
  real(kind=dp), dimension (:), allocatable :: urx , ury , urz
  real(kind=dp)                             :: idt , dtsq
  real(kind=dp)                             :: tempi , kin
 
  allocate ( urx ( natm ) , ury ( natm ) , urz ( natm ) )
  urx = 0.0_dp
  ury = 0.0_dp
  urz = 0.0_dp

  idt = 0.5_dp / dt 
  dtsq = dt * dt

  ! ==========================
  ! force + potential f(t)
  ! ==========================
  CALL engforce_driver 

  ! ================================================= 
  !  r(t+dt) = 2 r(t) - r (t-dt) + f(t) dt*dt
  ! ================================================= 
  do ia = 1 , natm
    urx ( ia ) = 2.0_dp * rx ( ia ) - rxs ( ia ) + fx ( ia ) * dtsq
    ury ( ia ) = 2.0_dp * ry ( ia ) - rys ( ia ) + fy ( ia ) * dtsq
    urz ( ia ) = 2.0_dp * rz ( ia ) - rzs ( ia ) + fz ( ia ) * dtsq 
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


! *********************** SUBROUTINE prop_velocity_verlet **********************
!
!> \brief
!! propagation for velocity-verlet algotithm (found a paper)
!
!> \todo
!! test it relatively to prop_leap_frog
!
! ******************************************************************************
SUBROUTINE prop_velocity_verlet 

  USE constants,                ONLY :  dp 
  USE config,                   ONLY :  natm, rx, ry, rz, vx, vy, vz, fx, fy, fz
  USE md,                       ONLY :  dt
  USE control,                  ONLY :  lpbc , lshiftpot
  USE thermodynamic,            ONLY :  temp_r , e_kin
  USE io_file,                  ONLY :  ionode
  USE field,                    ONLY :  engforce_driver

  implicit none

  ! local
  integer :: ia
  real(kind=dp) :: dtsq2 , dt2
  real(kind=dp), dimension (:), allocatable :: fsx , fsy , fsz
  real(kind=dp) :: tempi , kin


  ! save previous forces
  allocate (fsx(natm),fsy(natm),fsz(natm))
  fsx = 0.0_dp
  fsy = 0.0_dp
  fsz = 0.0_dp
  fsx = fx
  fsy = fy
  fsz = fz

  dtsq2 = dt * dt * 0.5_dp
  dt2 = dt * 0.5_dp

  ! =================================================
  !  r(t+dt) = r(t) + v(t)*dt + f(t) * dt*dt/2 
  !  store forces of the previous step  
  ! =================================================
  do ia = 1 , natm
    rx  ( ia ) = rx ( ia ) + vx ( ia ) * dt + fx ( ia ) * dtsq2
    ry  ( ia ) = ry ( ia ) + vy ( ia ) * dt + fy ( ia ) * dtsq2
    rz  ( ia ) = rz ( ia ) + vz ( ia ) * dt + fz ( ia ) * dtsq2
  enddo

  ! ==========================
  ! force + potential f(t+dt)
  ! ==========================
  CALL engforce_driver 

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

! *********************** SUBROUTINE nose_hoover_chain2 ************************
!
!> \brief
!! Nose-Hoover two chains
!
!> \note
!! adapted from Frenkel and Smit
!
! ******************************************************************************
SUBROUTINE nose_hoover_chain2 

  USE constants,                ONLY :  dp 
  USE config,                   ONLY :  natm , rx , ry , rz , vx , vy , vz , fx , fy , fz 
  USE md,                       ONLY :  dt, vxi, xi , Qnosehoover
  USE thermodynamic,            ONLY :  temp_r , e_kin , e_nvt

  implicit none

  ! local
  integer                :: inhc
  real(kind=dp)          :: kin , tempi
  real(kind=dp)          :: Q(2)

  ! thermostat mass
  Q(1) = natm * Qnosehoover
  Q(2) = Qnosehoover
 
  CALL calc_temp ( tempi , kin )

  CALL chain_nhc_2 ( kin , vxi , xi , Q )

  CALL prop_pos_vel_verlet ( kin )

  CALL chain_nhc_2( kin, vxi, xi , Q ) 


  tempi = (2.0_dp/3.0_dp) * kin
  tempi = tempi/DBLE (natm)

  e_nvt = 0.0_dp
  do inhc = 1 , 2
    e_nvt = e_nvt + vxi(inhc) * vxi(inhc) * 0.5_dp * Q(inhc)
  enddo
  e_nvt = e_nvt + 3.0_dp * natm * tempi * xi(1)
  e_nvt = e_nvt + tempi * xi(2)

  e_kin = kin
  temp_r = tempi
 
  return

END SUBROUTINE nose_hoover_chain2


! *********************** SUBROUTINE chain_nh_2 ********************************
!
!> \brief
!! intermediate routine used by nose_hoover_chain2
!
!> \note
!! adapted from Frenkel and Smit
!
! ******************************************************************************
SUBROUTINE chain_nhc_2 ( kin, vxi, xi , Q )

  USE constants,                ONLY :  dp 
  USE config,                   ONLY :  natm , vx , vy , vz 
  USE md,                       ONLY :  temp , dt 

  implicit none

  ! global
  real(kind=dp), intent (inout) :: kin
  real(kind=dp), intent (inout) :: vxi(2), xi(2) , Q(2)

  ! local
  integer       :: ia
  real(kind=dp) :: G1, G2  
  real(kind=dp) :: dt2, dt4, dt8
  real(kind=dp) :: s, L

  ! some constants related integrator step dt
  dt2 = dt  * 0.5_dp
  dt4 = dt2 * 0.5_dp
  dt8 = dt4 * 0.5_dp
  ! degree of freedom
  L  = DBLE (3 * natm)

  G2     = ( Q(1) * vxi(1) * vxi(1) - temp) / Q(2)
  vxi(2) = vxi(2) + G2 * dt4
  vxi(1) = vxi(1) * EXP ( - vxi(2) * dt8 )
  G1     = ( 2.0_dp *  kin - L * temp) / Q(1)
  vxi(1) = vxi(1) + G1 * dt4
  vxi(1) = vxi(1) * EXP ( - vxi(2) * dt8 )
  xi(1)  = xi(1) + vxi(1) * dt2
  xi(2)  = xi(2) + vxi(2) * dt2
  s      = EXP ( - vxi(1) * dt2 )
  kin    = kin * s * s

  do ia = 1, natm
    vx ( ia ) = s * vx ( ia )
    vy ( ia ) = s * vy ( ia )
    vz ( ia ) = s * vz ( ia )
  enddo
  
  vxi(1) = vxi(1) * EXP ( - vxi(2) * dt8 )
  G1     = ( 2.0_dp *  kin - L * temp) / Q(1)
  vxi(1) = vxi(1) + G1 * dt4
  vxi(1) = vxi(1) * EXP ( - vxi(2) * dt8 )
  G2     = ( Q(1) * vxi(1) * vxi(1) - temp) / Q(2)
  vxi(2) = vxi(2) + G2 * dt4

  return
 
END SUBROUTINE chain_nhc_2


! *********************** SUBROUTINE prop_pos_vel_verlet ***********************
!
!> \brief
!! propagates position and position in the velet algorithm
!
! ******************************************************************************
SUBROUTINE prop_pos_vel_verlet ( kin )

  USE constants,                ONLY :  dp 
  USE config,                   ONLY :  natm , rx , ry , rz , ry , vx , vy , vz , fx , fy , fz 
  USE md,                       ONLY :  dt
  USE control,                  ONLY :  lpbc
  USE field,                    ONLY :  engforce_driver

  implicit none

  ! global
  real(kind=dp), intent (out) :: kin

  ! local
  integer :: ia
  real(kind=dp) :: dt2

  dt2 = dt * 0.5_dp

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
  CALL engforce_driver 

  kin  = 0.0_dp
  do ia = 1 , natm
    vx ( ia ) = vx ( ia ) + fx ( ia ) * dt
    vy ( ia ) = vy ( ia ) + fy ( ia ) * dt
    vz ( ia ) = vz ( ia ) + fz ( ia ) * dt
    rx ( ia ) = rx ( ia ) + vx ( ia ) * dt2
    ry ( ia ) = ry ( ia ) + vy ( ia ) * dt2
    rz ( ia ) = rz ( ia ) + vz ( ia ) * dt2
    kin = kin + vx ( ia ) * vx ( ia ) +  vy ( ia ) * vy ( ia ) + vz ( ia ) * vz ( ia )
  enddo      
  kin = kin * 0.5_dp  


  return

END SUBROUTINE prop_pos_vel_verlet

! *********************** SUBROUTINE beeman ************************************
!
!> \brief
!! D. Beeman "Some multistep methods for use in molecular dynamics calculations",
!
!> \note
!! Journal of Computational Physics 20 pp. 130-139 (1976)
!
! ******************************************************************************
SUBROUTINE beeman 

  USE constants,                ONLY :  dp 
  USE config,                   ONLY :  natm , rx , ry , rz , ry , vx , vy , vz , fx , fy , fz , fxs , fys , fzs 
  USE thermodynamic,            ONLY :  temp_r , e_kin
  USE md,                       ONLY :  dt
  USE control,                  ONLY :  lpbc
  USE field,                    ONLY :  engforce_driver

  implicit none

  ! local
  integer :: ia
  real(kind=dp) :: onesix , twothree , onethree , fivesix , dtsq
  real(kind=dp) , dimension (:) , allocatable :: fxtmp , fytmp, fztmp
  real(kind=dp) :: kin , tempi


  allocate ( fxtmp (natm) , fytmp (natm) , fztmp (natm) ) 

  ! ================
  !  some constants
  ! ================
  onesix = 1.0_dp / 6.0_dp
  onethree = 2.0_dp * onesix
  twothree = 4.0_dp * onesix
  fivesix  = 5.0_dp * onesix
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
  CALL engforce_driver 

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

SUBROUTINE nose_hoover_chain_n 

  USE constants,                ONLY :  dp
  USE config,                   ONLY :  natm , rx , ry , rz , vx , vy , vz , fx , fy , fz
  USE md,                       ONLY :  dt, vxi, xi , Qnosehoover , nhc_n
  USE thermodynamic,            ONLY :  temp_r , e_kin , e_nvt

  implicit none

  ! local
  integer                :: inhc
  real(kind=dp)          :: kin , tempi
  real(kind=dp), dimension ( : ) , allocatable :: Q

  ! thermostat mass
  allocate ( Q(nhc_n) )
  Q=Qnosehoover
  Q(1) = Q(1) * natm 

  CALL calc_temp ( tempi , kin )
  CALL chain_nhc_n( kin , vxi , xi , Q )

!  CALL prop_pos_vel_verlet ( kin )
  CALL prop_velocity_verlet 

  CALL chain_nhc_n( kin, vxi, xi , Q )
  kin = kin * 0.5_dp

  tempi = (2.0_dp/3.0_dp) * kin 
  tempi = tempi / DBLE (natm)
  temp_r = tempi

  e_nvt = 0.0_dp
  do inhc = 1 , nhc_n 
    e_nvt = e_nvt + vxi(inhc) * vxi(inhc) * 0.5_dp * Q(inhc)
  enddo
  e_nvt = e_nvt + 3.0_dp * natm * tempi * xi(1)
  e_nvt = e_nvt + tempi * xi(2)
#ifdef debug_nvt_nhcn
  print*,'vxi',vxi
  print*,'xi',xi
  print*,'Q',Q
  print*,'tempi',tempi
  print*,'e_nvt',e_nvt / DBLE (natm)
#endif

  e_kin = kin 

  deallocate(Q)

  return

END SUBROUTINE nose_hoover_chain_n



! ref : 
! [1] Molecular Physics, (1996), v87, n5, p1117 Martyna and al.
! [2] Phys Lett. A, (1190), v150 n5,6,7, p262, Yoshida
SUBROUTINE chain_nhc_n ( kin , vxi , xi , Q )

  USE constants,                ONLY :  dp 
  USE config,                   ONLY :  natm , vx , vy , vz 
  USE md,                       ONLY :  temp , dt , nhc_n , nhc_yosh_order

  implicit none

  ! global
  real(kind=dp), intent (inout) :: kin
  real(kind=dp), intent (inout) :: vxi(nhc_n), xi(nhc_n) , Q(nhc_n)

  ! local
  integer :: ia , j , inh
  real(kind=dp), dimension ( : ) , allocatable :: G
  real(kind=dp) :: s , dts , dts2 , dts4 , dts8 , L

  real(kind=dp) , dimension ( : ), allocatable :: yosh_w  ! integrator order as in [2] YOSHIDA 
  real(kind=dp) , dimension ( : ), allocatable :: dt_yosh ! yoshida time 

  allocate ( yosh_w  ( nhc_yosh_order) )       
  allocate ( dt_yosh ( nhc_yosh_order) )       
       
  SELECT CASE ( nhc_yosh_order )
  CASE (1)
     yosh_w(1) = 1.0_dp
  CASE (3)
     yosh_w(1) = 1.0_dp/(2.0_dp-(2.0_dp)**(1.0_dp/3.0_dp))
     yosh_w(2) = 1.0_dp - 2.0_dp*yosh_w(1)
     yosh_w(3) = yosh_w(1)
  CASE (5)
     yosh_w(1) = 1.0_dp/(4.0_dp-(4.0_dp)**(1.0_dp/3.0_dp))
     yosh_w(2) = yosh_w(1)
     yosh_w(3) = yosh_w(1)
     yosh_w(4) = yosh_w(1)
     yosh_w(5) = 1.0_dp - 4.0_dp*yosh_w(1)
  CASE (7)
     yosh_w(1) = .78451361047756_dp
     yosh_w(2) = .235573213359357_dp
     yosh_w(3) = -1.17767998417887_dp
     yosh_w(4) = 1.0_dp - 2.0_dp*(yosh_w(1)+yosh_w(2)+yosh_w(3))
     yosh_w(5) = yosh_w(3)
     yosh_w(6) = yosh_w(2)
     yosh_w(7) = yosh_w(1)
  END SELECT

  kin = kin * 2.0_dp
  allocate ( G ( nhc_n) )
  ! thermostat mass
  L  = DBLE (3 * natm)
  dt_yosh =  yosh_w * dt

  s = 1.0_dp ! scale

  do j=1, nhc_yosh_order 

    dts = dt_yosh( j ) 
    dts2 = dts  * 0.5d0
    dts4 = dts2 * 0.5d0
    dts8 = dts4 * 0.5d0

    vxi ( nhc_n )  = vxi ( nhc_n ) + G ( nhc_n ) * dts4 
    !print*,"vxi(nhc_n)",vxi

    do inh=nhc_n-1,1,-1
      vxi ( inh )= vxi ( inh  ) * EXP ( - vxi ( inh + 1) * dts8 )
      vxi ( inh )    = vxi ( inh ) + G ( inh ) * dts4 
      vxi ( inh )= vxi ( inh ) * EXP ( - vxi ( inh + 1) * dts8 )
    enddo

    ! propagating xi 
    xi = xi + vxi * dts2

    s = s * EXP ( - vxi(1) * dts2 ) ! typo in instructions (35) [1] ?

    G ( 1 )   = ( kin*s - L * temp) / Q ( 1 )

    do inh = 1 , nhc_n - 1
      vxi ( inh )     = vxi ( inh ) * EXP ( - vxi ( inh + 1) * dts8 )  
      vxi ( inh )     = vxi(inh) + G(inh) * dts4
      vxi ( inh )     = vxi ( inh ) * EXP ( - vxi ( inh + 1) * dts8 )  
        G ( inh + 1 ) = ( Q(inh) * vxi(inh) * vxi(inh) - temp) / Q ( inh + 1 )
    enddo

     vxi ( nhc_n ) = vxi ( nhc_n ) + G ( nhc_n ) * dts4

  enddo

    do ia = 1, natm
      vx ( ia ) = s * vx ( ia )
      vy ( ia ) = s * vy ( ia )
      vz ( ia ) = s * vz ( ia )
    enddo
  

  deallocate ( G )
  deallocate ( yosh_w  )       
  deallocate ( dt_yosh )       

  return

END SUBROUTINE chain_nhc_n

! ===== fmV =====
