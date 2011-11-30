! ===== fmV =====
!!!!!!!!
!!!!!!!!
!================================================================================
! this subroutines are not clear at all. need more documentation and references.
! example: leap-frog and verlet are not clearly distinguised
! to be rigorous I should test all of them and compare to existing codes
! 
!================================================================================
!!!!!!!!
!!!!!!!!


!*********************** SUBROUTINE prop_leap_frog ****************************
!
! leap-frog algorithm  ( or verlet algorithm ) 
!
!******************************************************************************

SUBROUTINE prop_leap_frog ( iastart , iaend )!, list , point )

  USE control,          ONLY :  lpbc
  USE md,               ONLY :  dt
  USE config,           ONLY :  natm , rx , ry , rz , fx , fy , fz , vx , vy , vz , rxs , rys , rzs
  USE thermodynamic,    ONLY :  temp_r , e_kin
  USE field 

  implicit none

  ! global
  integer, intent(inout) :: iastart , iaend !, list(250 * natm) , point(natm + 1)

  ! local
  integer :: i
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
    CALL engforce ( iastart , iaend )!, list , point )

  ! ================================================= 
  !  r(t+dt) = 2 r(t) - r (t-dt) + f(t) dt*dt
  ! ================================================= 
  do i = 1, natm
    urx(i) = 2.0D0 * rx(i) - rxs(i) + fx(i) * dtsq
    ury(i) = 2.0D0 * ry(i) - rys(i) + fy(i) * dtsq
    urz(i) = 2.0D0 * rz(i) - rzs(i) + fz(i) * dtsq 
  enddo

  ! =====================================
  !  v(t) = ( r(t+dt) - r(t) ) / ( 2 dt) 
  ! =====================================
  do i = 1, natm
    vx(i) = idt * ( urx(i) - rxs(i) )
    vy(i) = idt * ( ury(i) - rys(i) )
    vz(i) = idt * ( urz(i) - rzs(i) )       
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
  do i = 1, natm
    rxs(i) = rx(i)
    rys(i) = ry(i)
    rzs(i) = rz(i)
    rx(i) = urx(i)
    ry(i) = ury(i)
    rz(i) = urz(i)
  enddo

  deallocate ( urx , ury , urz )

  return

END SUBROUTINE prop_leap_frog


!*********************** SUBROUTINE prop_velocity_verlet **********************
!
! propagation for velocity-verlet algotithm (found a paper)
!
!******************************************************************************

SUBROUTINE prop_velocity_verlet ( iastart , iaend )!, list , point )

  USE config,           ONLY :  natm, rx, ry, rz, vx, vy, vz, fx, fy, fz
  USE md,               ONLY :  dt
  USE control,          ONLY :  lpbc , lshiftpot
  USE thermodynamic,    ONLY :  temp_r , e_kin
  USE field

  implicit none

  ! global
  integer, intent(inout) :: iastart , iaend !, list(250 * natm) , point(natm + 1)

  ! local
  integer :: i
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
  do i = 1 , natm
    rx(i) = rx(i) + vx(i) * dt + fx(i) * dtsq2
    ry(i) = ry(i) + vy(i) * dt + fy(i) * dtsq2
    rz(i) = rz(i) + vz(i) * dt + fz(i) * dtsq2
    fsx(i) = fx(i)
    fsy(i) = fy(i)
    fsz(i) = fz(i)
  enddo

  ! ==========================
  ! force + potential f(t+dt)
  ! ==========================
  CALL engforce ( iastart , iaend )!, list , point )

  ! ==============================================
  !  v(t+dt) = v(t) + ( f(t-dt) + f(t) ) * dt / 2
  ! ==============================================
  do i = 1,natm
    vx(i) = vx(i) + ( fsx(i) + fx(i) ) * dt2
    vy(i) = vy(i) + ( fsy(i) + fy(i) ) * dt2
    vz(i) = vz(i) + ( fsz(i) + fz(i) ) * dt2
  enddo

  CALL calc_temp(tempi, kin)
  temp_r = tempi      
  e_kin  = kin

  deallocate ( fsx , fsy , fsz )

  return

END SUBROUTINE prop_velocity_verlet

!*********************** SUBROUTINE prop_velocity_verlet **********************
!
! propagation for velocity-verlet algotithm (found a paper)
! Note for test purpose only
!
!******************************************************************************

SUBROUTINE prop_velocity_verlet_test ( iastart , iaend )!, list , point )

  USE config,           ONLY :  natm, rx, ry, rz, vx, vy, vz, fx, fy, fz
  USE md,               ONLY :  dt
  USE control,          ONLY :  lpbc
  USE thermodynamic,    ONLY :  temp_r , e_kin
  USE field

  implicit none

  ! global
  integer, intent(inout) :: iastart , iaend !, list(250 * natm) , point(natm + 1)

  ! local
  double precision :: dtsq2 , dt2
  double precision, dimension (:), allocatable :: fsx, fsy, fsz
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
    rx(2) = rx(2) + vx(2) * dt + fx(2) * dtsq2
    ry(2) = ry(2) + vy(2) * dt + fy(2) * dtsq2
    rz(2) = rz(2) + vz(2) * dt + fz(2) * dtsq2
    fsx(2) = fx(2)
    fsy(2) = fy(2)
    fsz(2) = fz(2)

  ! ==========================
  ! force + potential f(t+dt)
  ! ==========================
  if (lpbc) then
    CALL engforce_bmlj_pbc_test ( )
  else
    print*,'only pbc allowed for this test'
    STOP
    CALL engforce_bmlj_nopbc ( iastart , iaend )!, list , point )
  endif

  ! ==============================================
  !  v(t+dt) = v(t) + ( f(t-dt) + f(t) ) * dt / 2
  ! ==============================================
  vx(1) = 0.0d0     
  vy(1) = 0.0d0     
  vz(1) = 0.0d0     
  vx(2) = vx(2) + ( fsx(2) + fx(2) ) * dt2
  vy(2) = vy(2) + ( fsy(2) + fy(2) ) * dt2
  vz(2) = vz(2) + ( fsz(2) + fz(2) ) * dt2

  CALL calc_temp(tempi, kin)
  temp_r = tempi      
  e_kin  = kin

  deallocate ( fsx , fsy , fsz )

  return

END SUBROUTINE prop_velocity_verlet_test


!*********************** SUBROUTINE nose_hoover_chain2 ************************
!
! Nose-Hoover chain (two) see Frenkel-Smit
!
!******************************************************************************

SUBROUTINE nose_hoover_chain2 ( iastart , iaend )!, list , point )

  USE config,           ONLY :  natm , rx , ry , rz , vx , vy , vz , fx , fy , fz 
  USE md,               ONLY :  dt, vxi1, vxi2, xi1, xi2
  USE thermodynamic,    ONLY :  temp_r , e_kin

  implicit none

  ! global
  integer, intent(inout) :: iastart , iaend !, list(natm * 250) , point(natm + 1)
  ! local
  double precision :: kin , tempi
 
  CALL calc_temp ( tempi , kin )

  CALL chain_nh_2 ( kin , vxi1 , vxi2 , xi1 , xi2 )

  CALL prop_pos_vel_verlet ( kin , iastart , iaend )!, list , point ) 

  CALL chain_nh_2( kin, vxi1, vxi2, xi1, xi2 ) 

  tempi = (2.0D0/3.0D0) * kin
  tempi = tempi/dble(natm)

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
  integer :: i
  double precision :: G1, G2, Q1, Q2  
  double precision :: dt2, dt4, dt8
  double precision :: s, L

  dt2 = dt  * 0.5d0
  dt4 = dt2 * 0.5d0
  dt8 = dt4 * 0.5D0

  Q1 = natm * Qnosehoover
  Q2 = Qnosehoover
  L = dble(3 * natm)

  G2   = ( Q1 * vxi1 * vxi1 - temp)
  vxi2 = vxi2 + G2 * dt4
  vxi1 = vxi1 * dexp ( - vxi2 * dt8 )
  G1 = ( 2.0D0 *  kin - L * temp) / Q1
  vxi1 = vxi1 + G1 * dt4
  vxi1 = vxi1 * dexp ( - vxi2 * dt8 )
  xi1 = xi1 + vxi1 * dt2
  xi2 = xi2 + vxi2 * dt2
  s = dexp( - vxi1 * dt2 )

  do i = 1, natm
    vx(i) = s * vx(i)
    vy(i) = s * vy(i)
    vz(i) = s * vz(i)
  enddo
  
  kin = kin * s * s
  vxi1 = vxi1 * dexp ( - vxi2 * dt8 )
  G1 = ( 2.0D0 *  kin - L * temp) / Q1
  vxi1 = vxi1 + G1 * dt4
  vxi1 = vxi1 * dexp ( - vxi2 * dt8 )
  G2   = ( Q1 * vxi1 * vxi1 - temp) / Q2
  vxi2 = vxi2 + G2 * dt4

  return
 
END SUBROUTINE chain_nh_2


!*********************** SUBROUTINE prop_pos_vel_verlet ***********************
!
! propagates position and position in the velet algorithm
!
!******************************************************************************

SUBROUTINE prop_pos_vel_verlet ( kin , iastart , iaend )!, list , point )

  USE config,   ONLY :  natm , rx , ry , rz , ry , vx , vy , vz , fx , fy , fz 
  USE md,       ONLY :  dt
  USE control,  ONLY :  lpbc
  USE field 

  implicit none

  ! global
  double precision, intent (out) :: kin
  integer, intent(inout) :: iastart , iaend !, list(250 * natm), point(natm + 1)

  ! local
  integer :: i
  double precision :: dt2

  dt2 = dt * 0.5d0

  ! ================================
  !  r(t+dt) = r(t) + v(t) * dt / 2
  ! ================================
  do i = 1,natm  
    rx(i) = rx(i) + vx(i) * dt2
    ry(i) = ry(i) + vy(i) * dt2
    rz(i) = rz(i) + vz(i) * dt2 
  enddo

  ! ==========================
  ! force + potential f(t+dt)
  ! ==========================
  CALL engforce ( iastart , iaend )!, list , point )

  kin  = 0.0d0
  do i = 1,natm
    vx(i) = vx(i) + fx(i) * dt
    vy(i) = vy(i) + fy(i) * dt
    vz(i) = vz(i) + fz(i) * dt
    rx(i) = rx(i) + vx(i) * dt2
    ry(i) = ry(i) + vy(i) * dt2
    rz(i) = rz(i) + vz(i) * dt2
    kin = kin + vx(i) * vx(i) +  vy(i) * vy(i) + vz(i) * vz(i)
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
SUBROUTINE beeman ( iastart , iaend )!, list , point ) 

  USE config,           ONLY :  natm , rx , ry , rz , ry , vx , vy , vz , fx , fy , fz , fxs , fys , fzs 
  USE thermodynamic,    ONLY :  temp_r , e_kin
  USE md,               ONLY :  dt
  USE control,          ONLY :  lpbc
  USE field

  implicit none

  ! global
  integer, intent(inout) :: iastart , iaend !, list(250 * natm), point(natm + 1)

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
  CALL engforce ( iastart , iaend )!, list , point )

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
