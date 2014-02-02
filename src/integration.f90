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
!#define debug_nvt_nhc2
!#define debug_nvt_nhcn
!#define debug_nvt_nhcpn
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
  USE config,                   ONLY :  natm , massia , rx , ry , rz , fx , fy , fz , vx , vy , vz , rxs , rys , rzs
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
  !  r(t+dt) = 2 r(t) - r (t-dt) + f(t)/m dt*dt
  ! ================================================= 
  do ia = 1 , natm
    urx ( ia ) = 2.0_dp * rx ( ia ) - rxs ( ia ) + fx ( ia ) * dtsq / massia(ia)
    ury ( ia ) = 2.0_dp * ry ( ia ) - rys ( ia ) + fy ( ia ) * dtsq / massia(ia)
    urz ( ia ) = 2.0_dp * rz ( ia ) - rzs ( ia ) + fz ( ia ) * dtsq / massia(ia)
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
  USE config,                   ONLY :  natm, massia, rx, ry, rz, vx, vy, vz, fx, fy, fz
  USE md,                       ONLY :  dt, integrator
  USE control,                  ONLY :  lpbc , lshiftpot
  USE thermodynamic,            ONLY :  temp_r , e_kin
  USE io,                  ONLY :  ionode
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
  !  r(t+dt) = r(t) + v(t)*dt + f(t)/m * dt*dt/2 
  !  store forces of the previous step  
  ! =================================================
  do ia = 1 , natm
    rx  ( ia ) = rx ( ia ) + vx ( ia ) * dt + (fx ( ia ) * dtsq2 ) / massia(ia)
    ry  ( ia ) = ry ( ia ) + vy ( ia ) * dt + (fy ( ia ) * dtsq2 ) / massia(ia)
    rz  ( ia ) = rz ( ia ) + vz ( ia ) * dt + (fz ( ia ) * dtsq2 ) / massia(ia)
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
SUBROUTINE nhc2 

  USE io,                       ONLY :  ioprint
  USE constants,                ONLY :  dp 
  USE config,                   ONLY :  natm , rx , ry , rz , vx , vy , vz , fx , fy , fz , center_of_mass, ntypemax
  USE md,                       ONLY :  dt, vxi, xi , timesca, nhc_n,temp
  USE thermodynamic,            ONLY :  temp_r , e_kin , e_nvt , e_tot,u_lj_r,h_tot, calc_thermo

  implicit none

  ! local
  integer                :: inhc
  real(kind=dp)          :: kin , tempi, L , com( 0 : ntypemax , 3 )
  real(kind=dp)          :: Q(2)

  ! degrees of freedom
  L = 3.0_dp * ( REAL ( natm, kind=dp) - 1.0_dp )

  ! thermostat mass
  Q(1) = timesca**2.0_dp * temp * L 
  Q(2) = timesca**2.0_dp * temp
 
  CALL calc_temp ( tempi , kin )
  CALL chain_nhc2 ( kin , vxi , xi , Q , L )

  CALL prop_pos_vel_verlet ( kin )

  CALL chain_nhc2( kin, vxi, xi , Q , L ) 

  CALL calc_temp ( tempi , kin )
  e_kin = kin
  temp_r = tempi
#ifdef debug_nvt_nhc2 
  CALL center_of_mass ( vx , vy , vz , com )
  io_print write(*,*) 'vel. com',com(0,:)
#endif

  e_nvt = 0.0_dp
  e_nvt = e_nvt + L * temp * xi(1)
  e_nvt = e_nvt + vxi(1) * vxi(1) * 0.5_dp / Q(1)
  do inhc = 2 , nhc_n 
    e_nvt = e_nvt + vxi(inhc) * vxi(inhc) * 0.5_dp / Q(inhc)
    e_nvt = e_nvt + temp * xi(inhc)
  enddo
#ifdef debug_nvt_nhc2 
  io_print write(*,'(a,4e16.8)') 'fmv xi vxi',xi,vxi
  CALL calc_thermo
  write(*,'(a,<4+nhc_n*2>e16.8)') 'fmv en',e_tot,u_lj_r,kin,h_tot,xi,vxi
#endif

 
  return

END SUBROUTINE nhc2


! *********************** SUBROUTINE chain_nhc2 ********************************
!
!> \brief
!! intermediate routine used by nose_hoover_chain2
!
!> \note
!! adapted from Frenkel and Smit
!
! ******************************************************************************
SUBROUTINE chain_nhc2 ( kin, vxi, xi , Q , L )

  USE io,                       ONLY :  ioprint
  USE constants,                ONLY :  dp 
  USE config,                   ONLY :  natm , vx , vy , vz 
  USE md,                       ONLY :  temp , dt , nhc_n

  implicit none

  ! global
  real(kind=dp), intent (inout) :: kin
  real(kind=dp), intent (inout) :: vxi(2), xi(2) 
  real(kind=dp), intent (in)    :: L, Q(2)

  ! local
  integer       :: ia
  real(kind=dp) :: G1, G2  
  real(kind=dp) :: dt2, dt4, dt8
  real(kind=dp) :: s

 
  ! some constants related integrator step dt
  s = 1._dp
  dt2 = dt  * 0.5_dp
  dt4 = dt2 * 0.5_dp
  dt8 = dt4 * 0.5_dp

  G1     = ( 2.0_dp*kin - L * temp) / Q(1)
  G2     = ( Q(1) * vxi(1) * vxi(1) - temp) / Q(2)
#ifdef debug_nvt_nhc2
  write(*,'(a,e16.8)') "G2",G2
  write(*,'(a,3e16.8)') "G1",G1,Q(1),L
#endif
  vxi(2) = vxi(2) + G2 * dt4
  vxi(1) = vxi(1) * EXP ( - vxi(2) * dt8 )
  vxi(1) = vxi(1) + G1 * dt4
  vxi(1) = vxi(1) * EXP ( - vxi(2) * dt8 )
#ifdef debug_nvt_nhc2
    io_print write(*,'(a,<nhc_n>e16.8)') "vxi(nhc_n)",vxi
    io_print write(*,'(a,<nhc_n>e16.8)') "G",G1,G2
#endif
  s   = s * EXP ( - vxi(1) * dt2 )
  kin = s * s * kin
  io_print print*,'FMV FMV FMV FMV s =',s

  xi  = xi + vxi * dt2
  G1     = ( 2.0_dp*kin - L * temp) / Q(1)
  vxi(1) = vxi(1) * EXP ( - vxi(2) * dt8 )
  vxi(1) = vxi(1) + G1 * dt4
  vxi(1) = vxi(1) * EXP ( - vxi(2) * dt8 )
  G2     = ( Q(1) * vxi(1) * vxi(1) - temp) / Q(2)
  vxi(2) = vxi(2) + G2 * dt4
#ifdef debug_nvt_nhc2
    io_print write(*,'(a,<nhc_n>e16.8)') "xi(nhc_n)",xi
    io_print write(*,'(a,<nhc_n>e16.8)') "G",G1,G2
#endif

  do ia = 1, natm
    vx ( ia ) = s * vx ( ia )
    vy ( ia ) = s * vy ( ia )
    vz ( ia ) = s * vz ( ia )
  enddo

  return
 
END SUBROUTINE chain_nhc2


! *********************** SUBROUTINE prop_pos_vel_verlet ***********************
!
!> \brief
!! propagates position and position in the velet algorithm
!
! ******************************************************************************
SUBROUTINE prop_pos_vel_verlet ( kin )

  USE constants,                ONLY :  dp 
  USE config,                   ONLY :  natm , massia , rx , ry , rz , ry , vx , vy , vz , fx , fy , fz 
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

  ! =========================================
  !  v(t+dt2) = v(t) + f(t) * dt2
  !  r(t+dt)  = r(t) + v(t+dt2) * dt / m 
  ! note : dt2 = dt / 2
  ! =========================================
  do ia = 1 , natm  
    vx ( ia ) = vx ( ia ) + fx ( ia ) * dt2 / massia(ia) 
    vy ( ia ) = vy ( ia ) + fy ( ia ) * dt2 / massia(ia)
    vz ( ia ) = vz ( ia ) + fz ( ia ) * dt2 / massia(ia)
    rx ( ia ) = rx ( ia ) + vx ( ia ) * dt 
    ry ( ia ) = ry ( ia ) + vy ( ia ) * dt  
    rz ( ia ) = rz ( ia ) + vz ( ia ) * dt  
  enddo

  ! ==========================
  ! f(t+dt)
  ! ==========================
  CALL engforce_driver 

  ! =========================================
  !   v(t) = v(t+dt2) + f(t+dt) * dt2
  !   v(t) = v(t) + dt2 * ( f(t) + f(t+ft) ) 
  ! =========================================
  kin  = 0.0_dp
  do ia = 1 , natm
    vx ( ia ) = vx ( ia ) + fx ( ia ) * dt2 / massia(ia)
    vy ( ia ) = vy ( ia ) + fy ( ia ) * dt2 / massia(ia)
    vz ( ia ) = vz ( ia ) + fz ( ia ) * dt2 / massia(ia)
    kin = kin + ( vx ( ia ) * vx ( ia ) +  vy ( ia ) * vy ( ia ) + vz ( ia ) * vz ( ia ) ) * massia (ia)
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
  USE config,                   ONLY :  natm , massia , rx , ry , rz , ry , vx , vy , vz , fx , fy , fz , fxs , fys , fzs 
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
  ! r (t+dt) = r (t) + v (t) dt  + ( 2/3 f (t) - 1/6 f (t-dt) ) / m dt^2
  ! ==================================================================
  do ia = 1 , natm
    rx ( ia ) = rx ( ia ) + vx ( ia ) * dt + ( twothree * fx ( ia ) - onesix * fxs ( ia ) ) * dtsq / massia(ia)
    ry ( ia ) = ry ( ia ) + vy ( ia ) * dt + ( twothree * fy ( ia ) - onesix * fys ( ia ) ) * dtsq / massia(ia)
    rz ( ia ) = rz ( ia ) + vz ( ia ) * dt + ( twothree * fz ( ia ) - onesix * fzs ( ia ) ) * dtsq / massia(ia)
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
    vx ( ia ) = vx ( ia ) + ( onethree * fx ( ia ) + fivesix * fxtmp ( ia ) - onesix * fxs ( ia ) ) * dt / massia(ia)    
    vy ( ia ) = vy ( ia ) + ( onethree * fy ( ia ) + fivesix * fytmp ( ia ) - onesix * fys ( ia ) ) * dt / massia(ia)
    vz ( ia ) = vz ( ia ) + ( onethree * fz ( ia ) + fivesix * fztmp ( ia ) - onesix * fzs ( ia ) ) * dt / massia(ia)
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

! *********************** SUBROUTINE nose_hoover_chain_n ***********************
!
!> \brief
!
! ******************************************************************************
SUBROUTINE nhcn

  USE io,                       ONLY :  ioprint
  USE constants,                ONLY :  dp
  USE config,                   ONLY :  natm , rx , ry , rz , vx , vy , vz , fx , fy , fz , center_of_mass, ntypemax
  USE md,                       ONLY :  dt, vxi, xi , timesca, nhc_n,temp
  USE thermodynamic,            ONLY :  temp_r , e_kin , e_nvt , e_tot,u_lj_r,h_tot, calc_thermo

  implicit none

  ! local
  integer                :: inhc
  real(kind=dp)          :: kin , tempi, L , com( 0 : ntypemax , 3 )
  real(kind=dp)          :: Q(2)

  ! degrees of freedom
  L = 3.0_dp * ( REAL ( natm, kind=dp) - 1.0_dp )

  ! thermostat mass
  Q(1) = timesca**2.0_dp * temp * L
  Q(2) = timesca**2.0_dp * temp

  CALL calc_temp ( tempi , kin )
  CALL chain_nhcn ( kin , vxi , xi , Q , L )

  CALL prop_pos_vel_verlet ( kin )

  CALL chain_nhcn( kin, vxi, xi , Q , L )

  CALL calc_temp ( tempi , kin )
  e_kin = kin
  temp_r = tempi
#ifdef debug_nvt_nhcn 
  CALL center_of_mass ( vx , vy , vz , com )
  io_print write(*,*) 'vel. com',com(0,:)
#endif

  e_nvt = 0.0_dp
  e_nvt = e_nvt + L * temp * xi(1)
  e_nvt = e_nvt + vxi(1) * vxi(1) * 0.5_dp / Q(1)
  do inhc = 2 , nhc_n
    e_nvt = e_nvt + vxi(inhc) * vxi(inhc) * 0.5_dp / Q(inhc)
    e_nvt = e_nvt + temp * xi(inhc)
  enddo
#ifdef debug_nvt_nhcn 
  io_print write(*,'(a,4e16.8)') 'fmv xi vxi',xi,vxi
  CALL calc_thermo
  write(*,'(a,<4+nhc_n*2>e16.8)') 'fmv en',e_tot,u_lj_r,kin,h_tot,xi,vxi
#endif


  return

END SUBROUTINE nhcn

! *********************** SUBROUTINE chain_nhcn ***********************
!
!> \brief
! ref : 
! [1] Molecular Physics, (1996), v87, n5, p1117 Martyna and al.
! [2] Phys Lett. A, (1190), v150 n5,6,7, p262, Yoshida
!
! ******************************************************************************
SUBROUTINE chain_nhcn ( kin , vxi , xi , Q , L )

  USE constants,                ONLY :  dp 
  USE config,                   ONLY :  natm , vx , vy , vz, massia 
  USE md,                       ONLY :  temp , dt , nhc_n , nhc_yosh_order,nhc_mults
  USE io,                       ONLY :  ioprint

  implicit none

  ! global
  real(kind=dp), intent (inout) :: kin
  real(kind=dp), intent (inout) :: vxi(nhc_n), xi(nhc_n) , Q(nhc_n)

  ! local
  integer :: ia , j , k, inh
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

  allocate ( G ( nhc_n) )
  dt_yosh =  yosh_w * dt

  s = 1.0_dp ! scale
!  G(1) = ( 2.0_dp*kin*s - L * temp) / Q ( 1 )
!  do inh=1,nhc_n-1
!    G ( inh + 1 ) = ( Q(inh) * vxi(inh) * vxi(inh) - temp) / Q ( inh + 1 )
!  enddo  

do k=1,nhc_mults
  do j=1, nhc_yosh_order 

    dts = dt_yosh( j ) 
    dts2 = dts  * 0.5d0
    dts4 = dts2 * 0.5d0
    dts8 = dts4 * 0.5d0 

    G   ( nhc_n ) = ( Q(nhc_n-1) * vxi(nhc_n-1) * vxi(nhc_n-1) - temp) / Q ( nhc_n )
    vxi ( nhc_n ) = vxi ( nhc_n ) + G ( nhc_n ) * dts4 
#ifdef debug_nvt_nhcn
    io_print write(*,'(a,<nhc_n>e16.8)') "vxi(nhc_n)",vxi
    write(*,'(a,e16.8)') "G2",G(2)
#endif
    do inh=1,nhc_n-1
      !print*,'get in 1',inh
      vxi ( nhc_n - inh ) = vxi ( nhc_n - inh ) * EXP ( - vxi ( nhc_n + 1 - inh ) * dts8 )
      vxi ( nhc_n - inh ) = vxi ( nhc_n - inh ) + G ( nhc_n - inh ) * dts4 
      vxi ( nhc_n - inh ) = vxi ( nhc_n - inh ) * EXP ( - vxi ( nhc_n + 1 - inh ) * dts8 )
!      G   ( inh - 1 ) = ( Q(inh) * vxi(inh) * vxi(inh) - temp) / Q ( inh - 1 )
    enddo
    G ( 1 )   = ( 2.0_dp*kin - L * temp) / Q ( 1 )
#ifdef debug_nvt_nhcn
    write(*,'(a,3e16.8)') "G1",G(1),Q(1),L
#endif
    vxi ( 1 ) = vxi ( 1 ) * EXP ( - vxi ( 2 ) * dts8 )
    vxi ( 1 ) = vxi ( 1 ) + G ( 1 ) * dts4
    vxi ( 1 ) = vxi ( 1 ) * EXP ( - vxi ( 2 ) * dts8 )

    s = s * EXP ( - vxi(1) * dts2 ) ! typo in instructions (35) [1] ?
    kin = kin * s * s 
#ifdef debug_nvt_nhcn
    io_print print*,'FMV FMV FMV FMV s =',s
#endif
    ! propagating xi 
    xi = xi + vxi * dts2
!====================================================================================================
!====================================================================================================
!  G1     = ( 2.0_dp*kin - L * temp) / Q(1)
!  G2     = ( Q(1) * vxi(1) * vxi(1) - temp) / Q(2)
!  vxi(2) = vxi(2) + G2 * dt4
!  vxi(1) = vxi(1) * EXP ( - vxi(2) * dt8 )
!  vxi(1) = vxi(1) + G1 * dt4
!  vxi(1) = vxi(1) * EXP ( - vxi(2) * dt8 )
!  s   = s * EXP ( - vxi(1) * dt2 )
!  kin = s * s * kin
!  xi  = xi + vxi * dt2
!  G1     = ( 2.0_dp*kin - L * temp) / Q(1)
!  vxi(1) = vxi(1) * EXP ( - vxi(2) * dt8 )
!  vxi(1) = vxi(1) + G1 * dt4
!  vxi(1) = vxi(1) * EXP ( - vxi(2) * dt8 )
!  G2     = ( Q(1) * vxi(1) * vxi(1) - temp) / Q(2)
!  vxi(2) = vxi(2) + G2 * dt4
!====================================================================================================
!====================================================================================================
#ifdef debug_nvt_nhcn
    io_print write(*,'(a,<nhc_n+1>e16.8)') "xi(nhc_n) ",xi,kin
#endif
    G ( 1 )   = ( 2.0_dp*kin - L * temp) / Q ( 1 )
    do inh = 1 , nhc_n - 1
!     print*,'get in 2',inh
      vxi ( inh )     = vxi ( inh ) * EXP ( - vxi ( inh + 1) * dts8 )  
      vxi ( inh )     = vxi(inh) + G(inh) * dts4
      vxi ( inh )     = vxi ( inh ) * EXP ( - vxi ( inh + 1) * dts8 )  
        G ( inh + 1 ) = ( Q(inh) * vxi(inh) * vxi(inh) - temp) / Q ( inh + 1 )
    enddo
    vxi ( nhc_n ) = vxi ( nhc_n ) + G ( nhc_n ) * dts4
  enddo
enddo
  kin = 0.0_dp
  do ia = 1, natm
    vx ( ia ) = s * vx ( ia )
    vy ( ia ) = s * vy ( ia )
    vz ( ia ) = s * vz ( ia )
    kin =  kin + ( vx ( ia ) ** 2 + vy ( ia ) ** 2 + vz ( ia ) ** 2 ) / massia(ia)
  enddo
  kin = kin * 0.5_dp

  deallocate ( G )
  deallocate ( yosh_w  )       
  deallocate ( dt_yosh )       

  return

END SUBROUTINE chain_nhcn

! *********************** SUBROUTINE chain_nhcpn ***********************
!
!> \brief
! ref : 
! [1] Molecular Physics, (1996), v87, n5, p1117 Martyna and al.
! [2] Phys Lett. A, (1190), v150 n5,6,7, p262, Yoshida
!
! ******************************************************************************
SUBROUTINE chain_nhcpn ( kin , vxi , xi , Q , W )

  USE constants,                ONLY :  dp
  USE config,                   ONLY :  simu_cell, natm , vx , vy , vz
  USE md,                       ONLY :  temp , press, dt , nhc_n , nhc_yosh_order
  USE thermodynamic,            ONLY :  pressure_tot, calc_thermo

  implicit none

  ! global
  real(kind=dp), intent (inout) :: kin , W
  real(kind=dp), intent (inout) :: vxi(nhc_n), xi(nhc_n) , Q(nhc_n)

  ! local
  integer :: ia , j , inh
  real(kind=dp), dimension ( : ) , allocatable :: G
  real(kind=dp) :: s , dts , dts2 , dts4 , dts8 , L, Ge , ve , P, coeff

  real(kind=dp) , dimension ( : ), allocatable :: yosh_w  ! integrator order as in [2] YOSHIDA 
  real(kind=dp) , dimension ( : ), allocatable :: dt_yosh ! yoshida time 

  print*,'in'
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

  print*,'in 2'
  CALL calc_thermo
  P = pressure_tot
  kin = kin * 2.0_dp
  allocate ( G ( nhc_n) )
  ! thermostat mass
  L  = DBLE (3 * natm)
  dt_yosh =  yosh_w * dt

  s = 1.0_dp ! scale
  print*,'in 3',dt_yosh( 1 )

  do j=1, nhc_yosh_order

    write(*,*) j
    dts = dt_yosh( j )
    dts2 = dts  * 0.5d0
    dts4 = dts2 * 0.5d0
    dts8 = dts4 * 0.5d0

    vxi ( nhc_n )  = vxi ( nhc_n ) + G ( nhc_n ) * dts4
    write(*,'(a,<nhc_n>e16.8)') "vxi(nhc_n)",vxi
    do inh=nhc_n-1,1,-1
      vxi ( inh )= vxi ( inh  ) * EXP ( - vxi ( inh + 1) * dts8 )
      vxi ( inh )= vxi ( inh ) + G ( inh ) * dts4
      vxi ( inh )= vxi ( inh ) * EXP ( - vxi ( inh + 1) * dts8 )
    enddo

    ve = ve * EXP ( - vxi ( 1) * dts8 )
    ve = ve + Ge * dts4
    ve = ve * EXP ( - vxi ( 1) * dts8 )
      
    ! propagating xi 
    xi = xi + vxi * dts2

    s = s * EXP ( - vxi(1) * dts2 ) ! typo in instructions (35) [1] ?

    G ( 1 )   = ( kin*s - L * temp + W * ve ** 2 ) / Q ( 1 )
    Ge        = ( coeff * kin + 3.0_dp * ( P - press ) * simu_cell%omega ) / W


    ve = ve * EXP ( - vxi ( 1) * dts8 )
    ve = ve + Ge * dts4
    ve = ve * EXP ( - vxi ( 1) * dts8 )

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

  print*,'end'

  deallocate ( G )
  deallocate ( yosh_w  )
  deallocate ( dt_yosh )

  return

END SUBROUTINE chain_nhcpn

! *********************** SUBROUTINE nhcpn ***********************
!
!> \brief
!
! ******************************************************************************
SUBROUTINE nhcpn

  USE constants,                ONLY :  dp
  USE config,                   ONLY :  natm , rx , ry , rz , vx , vy , vz , fx , fy , fz
  USE md,                       ONLY :  dt, vxi, xi , timesca , Wmass, nhc_n,temp
  USE thermodynamic,            ONLY :  temp_r , e_kin , e_npt, pressure_tot, calc_thermo

  implicit none

  ! local
  integer                :: inhc
  real(kind=dp)          :: kin , tempi , pressi, W
  real(kind=dp), dimension ( : ) , allocatable :: Q

  ! thermostat mass
  allocate ( Q(nhc_n) )
  Q=timesca**2.0_dp
  Q(1) = Q(1) * natm
  W = Wmass

  CALL calc_temp ( tempi , kin )

  CALL chain_nhcpn( kin , vxi , xi , Q , W )

  CALL prop_velocity_verlet

  CALL chain_nhcpn( kin, vxi, xi , Q , W )
  kin = kin * 0.5_dp

  tempi = (2.0_dp/3.0_dp) * kin
  tempi = tempi / DBLE (natm)
  temp_r = tempi

  e_npt = 0.0_dp
  do inhc = 1 , nhc_n
    e_npt = e_npt + vxi(inhc) * vxi(inhc) * 0.5_dp * Q(inhc)
  enddo
  e_npt = e_npt + 3.0_dp * natm * temp * xi(1)
  e_npt = e_npt + temp * xi(2)

  CALL calc_thermo
  pressi=pressure_tot

#ifdef debug_nvt_nhcpn
  print*,'vxi',vxi
  print*,'xi',xi
  print*,'Q',Q
  print*,'W',W
  print*,'temp',temp
  print*,'pressi',pressi
  print*,'e_npt',e_npt / DBLE (natm)
#endif

  e_kin = kin

  deallocate(Q)

  return

END SUBROUTINE nhcpn


! ===== fmV =====
