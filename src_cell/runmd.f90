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
!#define debug

! calculate the stress tensor at each time step
!#define stress_t 

! calculates efg at each time step
!#define efg_t

! calculates center of mass at each time step
!#define com_t

!experimental
!calculate block_averaging 
!#define block

!experimental
!multi_tau correlation
!#define multi_tau

! ======= Hardware =======

!*********************** SUBROUTINE md_run **********************************
!
! main SUBROUTINE for bmlj molecular dynamics and static calculation
! properties EFG, GR ... are called here when calculated on the fly.
!
! input : 
!          iastart , iaend : index for atom decomposition (parallel)           
!          list , point    : verlet list variables 
!          offset          : time offset
!
! properties : 
!  
!          EFG             : electric field gradient
!          GR              : radial distribution function
!          MSD             : mean square displacement
!          VACF            : velocity auto-correlation function
!
! Note : previously called bmlj_run
!
!******************************************************************************

SUBROUTINE md_run ( iastart , iaend , offset )


  USE config,   ONLY :  natm , rx , ry , rz , rxs , rys , rzs , vx , vy , vz , &
                        write_CONTFF , center_of_mass , ntypemax , tau_lj , tau_coul 
  USE control,  ONLY :  lpbc , longrange , calc , lstatic , lvnlist , lbmlj , lcoulomb , lmorse
  USE io_file,  ONLY :  ionode , stdout, kunit_OSZIFF, kunit_TRAJFF , kunit_EFGFF , &
                        kunit_EFGALL , kunit_OUTFF , kunit_EQUILFF
  USE prop,     ONLY :  lstrfac , lmsd , lvacf , nprop , nprop_start
  USE md,       ONLY :  npas , ltraj , lleapequi , itraj_period , itraj_start , nequil , nequil_period , nprint, &
                        fprint, spas , dt,  temp , updatevnl , write_traj_xyz , integrator

  USE thermodynamic
  USE time
  USE field 
  USE efg
  USE radial_distrib
  USE msd
  USE vacf 
  USE multi_tau 

  implicit none
  INCLUDE "mpif.h"

  ! global
  integer, intent(inout) :: iastart , iaend 
  integer, intent(in) :: offset

  ! local
  double precision :: tempi , kin 
  double precision :: ttt1, ttt2 
  double precision, dimension(:), allocatable :: xtmp , ytmp , ztmp
  integer :: itime , ierr , ia 
  integer :: nefg , ngr , nmsd , ntau 
  character*60 :: rescale_allowed(2)
  data rescale_allowed / 'nve-vv' , 'nve-be' /
 ! tmp 
#ifdef com_t
  double precision :: com(0:ntypemax,3) !center of mass
  double precision :: ddtt , sumcom1 , sumcom2, sumcom1sq , sumcom2sq 
  double precision :: mcom1 , mcom2 , msqcom1 , msqcom2, scom1 , scom2 , modcom1 , modcom2      
  integer :: comcount

  comcount = 0
  ddtt=0.0d0;sumcom1sq=0.0d0;sumcom2sq=0.0d0;sumcom1=0.0d0;sumcom2=0.0d0
#endif

  ! =========================================
  ! properties counter where is vacf ??? 
  ! =========================================
  nefg = 0
  ngr  = 0
  nmsd = 0
  ntau = 0
  
  CALL init_general_accumulator

#ifdef multi_tau
  ! related to multi_tau
  call alloc
#endif

  if ( ionode ) WRITE ( stdout , '(a,a)' ) &
  '===============================================================================================', &
  '======================================'
  if ( ionode ) WRITE ( kunit_OUTFF , '(a,a)' ) &
  '===============================================================================================', &
  '======================================'
  if ( ionode ) WRITE ( stdout , '(a)' )      'properties at t=0'
  if ( ionode ) WRITE ( kunit_OUTFF , '(a)' ) 'properties at t=0'
 
  allocate( xtmp(natm), ytmp(natm), ztmp(natm) )

  OPEN (unit = kunit_OSZIFF ,file = 'OSZIFF',STATUS = 'UNKNOWN')
  OPEN (unit = kunit_TRAJFF ,file = 'TRAJFF')
#ifdef block
  OPEN (unit = kunit_EQUILFF,file = 'EQUILFF',STATUS = 'UNKNOWN')
#endif

  ! ===================================
  !  calc. kinetic temperature at t=0
  ! ===================================
  CALL calc_temp ( tempi , kin )
  e_kin  = kin
  temp_r = tempi

  
  ! ==============================================
  !  for integrator .ne. than leapfrog algorithm
  ! ==============================================
  if (integrator.ne.'nve-lf') then
   
    ! =========================
    ! force + potential at t=0
    ! =========================
     CALL engforce_driver ( iastart , iaend )

    ! ==================================================
    ! write thermodynamic information of config at t=0
    ! ==================================================
    CALL write_thermo( offset-1 , stdout )
    CALL write_thermo( offset-1 , kunit_OUTFF )
    CALL write_thermo( offset-1 , kunit_OSZIFF )

#ifdef debug
         CALL print_config_sample(0,0)
#endif

  else
  ! ===========================================================================
  !  leapfrog algorithm set the first leap value for r(t-dt)  (approximation)
  ! ===========================================================================
    xtmp = rx
    ytmp = ry
    ztmp = rz
    do ia = 1 , natm
      rxs ( ia )  = rx ( ia ) - vx ( ia ) * dt
      rys ( ia )  = ry ( ia ) - vy ( ia ) * dt
      rzs ( ia )  = rz ( ia ) - vz ( ia ) * dt
    enddo
    rx = rxs 
    ry = rys
    rz = rzs 
    ! ===================
    ! force + potential
    ! ===================
    CALL engforce_driver ( iastart , iaend )
    ! =======================================================
    ! write thermodynamic information of the starting point
    ! =======================================================
    CALL write_thermo( offset-1 , stdout )
    CALL write_thermo( offset-1 , kunit_OUTFF )
    CALL write_thermo( offset-1 , kunit_OSZIFF )
     rx = xtmp
     ry = ytmp
     rz = ztmp
  endif 

  ! =======================
  !  stress tensor at t=0
  ! =======================
  if ( ionode ) WRITE ( stdout , '(a)' ) ' ' 
  if ( ionode ) WRITE ( stdout , '(a)' ) 'stress tensor of initial configuration' 
  if ( ionode ) then
    CALL print_tensor ( tau_lj   , 'TAU_LJ  ' ) 
    CALL print_tensor ( tau_coul , 'TAU_COUL' ) 
  endif

  ! =========================
  !   MAIN LOOP ( TIME )
  ! =========================
  if ( ionode ) WRITE ( stdout , '(a)' ) ''
  if ( ionode ) WRITE ( stdout , '(a)' ) 'starting main loop'
  if ( ionode ) WRITE ( stdout , '(a,a)' ) &
  '===============================================================================================', &
  '======================================'
  if ( ionode ) WRITE ( kunit_OUTFF , '(a)' ) ''
  if ( ionode ) WRITE ( kunit_OUTFF , '(a)' ) 'starting main loop'
  if ( ionode ) WRITE ( kunit_OUTFF , '(a,a)' ) &
  '===============================================================================================', &
  '======================================'
! be careful with the offset !!!!
MAIN:  do itime = offset , npas + (offset-1)        

#ifdef block
          if ( itime .ge. nequil ) then 
            CALL write_thermo(itime, kunit_EQUILFF)
            CALL general_accumulator
          endif
#endif
         ! ==========================================================
         !  leap-frog algorithm set the first leap value for r(t-dt) 
         !  if first steps are equilibrated with a velocity-verlet 
         ! ==========================================================
         if ( lleapequi .and. itime .ge. nequil ) then
           integrator = 'nve-lf'
           lleapequi  = .false.   
           do ia = 1 , natm
             rxs ( ia )  = rx ( ia ) - vx ( ia ) * dt
             rys ( ia )  = ry ( ia ) - vy ( ia ) * dt
             rzs ( ia )  = rz ( ia ) - vz ( ia ) * dt
           enddo
         endif
         
         ! ===========
         !  time info
         ! ===========
         ttt1 = MPI_WTIME ( ierr )

         ! =========================    
         !  integration t -> t + dt 
         ! =========================
         if ( integrator.eq.'nve-lf'  )      CALL prop_leap_frog ( iastart , iaend )
         if ( integrator.eq.'nve-be'  )      CALL beeman ( iastart , iaend )
         if ( integrator.eq.'nve-vv'  )      CALL prop_velocity_verlet ( iastart , iaend )
         if ( integrator.eq.'nvt-and' )      CALL prop_velocity_verlet ( iastart , iaend )
         if ( integrator.eq.'nvt-nhc2')      CALL nose_hoover_chain2 ( iastart , iaend )

#ifdef debug
         CALL print_config_sample(itime,0)
#endif

         ! ===========
         !  time info
         ! ===========
         ttt2 = MPI_WTIME(ierr)
         mdsteptimetot = mdsteptimetot + (ttt2 - ttt1) 

         ! ================================
         !  rescale velocities (NVE equil)
         ! ================================
         if ( ( ANY ( integrator .eq. rescale_allowed ) ) .and.  &
                  ( ( itime .le. nequil.and.MOD ( itime , nequil_period ) .eq. 0 ) .and. &
                    ( itime .ne. npas + ( offset - 1 ) .and. itime .ne. offset ) ) ) then
           CALL rescale_velocities(0)
         endif

         ! ===================================
         !  rescale velocities (NVT Andersen)
         ! ===================================
         if ( (integrator.eq.'nvt-and') ) CALL andersen_velocities

         ! ===================
         !  print trajectory 
         ! ===================
         if ( ltraj .and. (itime .gt. itraj_start ) .and. MOD (itime,itraj_period) .eq. 0 ) then
           xtmp = rx
           ytmp = ry
           ztmp = rz
           CALL write_traj_xyz
           rx = xtmp
           ry = ytmp
           rz = ztmp
         endif
#ifdef stress_t
  if ( ionode ) WRITE ( stdout , '(a)' ) ' ' 
  if ( ionode ) WRITE ( stdout , '(a)' ) 'stress tensor of initial configuration' 
  if ( ionode ) then
    if ( lbmlj )    CALL print_tensor ( tau_lj    , 'TAU_LJ  ' ) 
    if ( lcoulomb ) CALL print_tensor ( tau_coul  , 'TAU_COUL' ) 
    if ( lmorse )   CALL print_tensor ( tau_morse , 'TAU_MORS' ) 
  endif
#endif
       
#ifdef com_t
     if ( MOD ( itime , nprint ) .eq. 0 ) then
       comcount = comcount + 1 
       CALL center_of_mass ( rx , ry , rz , com )
       write(*,'(a3,i10,3e18.6)') 'ALL',itime , com(0,:)
       write(*,'(a3,i10,3e18.6)') 'A  ',itime , com(1,:)
       write(*,'(a3,i10,3e18.6)') 'B  ',itime , com(2,:)
       modcom1=com(1,1)*com(1,1)+com(1,2)*com(1,2)+com(1,3)*com(1,3) 
       modcom2=com(2,1)*com(2,1)+com(2,2)*com(2,2)+com(2,3)*com(2,3) 
       ! sum 
       sumcom1   = sumcom1 + SQRT (modcom1)
       sumcom2   = sumcom2 + SQRT (modcom2)
       ! sum od square
       sumcom1sq = sumcom1sq + modcom1
       sumcom2sq = sumcom2sq + modcom2
       ! counting  
       ddtt      =1.0D0 / DBLE (comcount)
       ! average
       mcom1     = sumcom1*ddtt
       mcom2     = sumcom2*ddtt 
       ! average of square 
       msqcom1   = sumcom1sq*ddtt 
       msqcom2   = sumcom2sq*ddtt 
       !fluctu
       scom1     = msqcom1 - mcom1*mcom1
       scom2     = msqcom2 - mcom2*mcom2
       write(*,'(a3,i10,4e18.6)') 'MOY',itime,mcom1,mcom2,scom1,scom2
    endif
#endif
         ! =======================
         !  properties on-the-fly
         ! =======================
         ! -----------------------------------------------------------------------------------------
         if ( MOD (itime,nprop) .eq. 0 .and. itime.gt.nprop_start) then 

           ! =======
           !  lvacf
           ! =======
           if ( lvacf ) then
             CALL vacf_main 
           endif

           ! =======
           !  lmsd
           ! =======
           if ( lmsd ) then
             nmsd = nmsd + 1
             CALL msd_main ( nmsd )
           endif
  
#ifdef multi_tau
             ntau = ntau + 1
             CALL multi_tau_main ( vx , vy , vz , ntau )
#endif

         endif 
        ! ----------------------------------------------------------------------------------

        ! =============================================
        !  write instanteanous thermodynamic properties
        !  to standard output and OUTFF
        ! =============================================
        if ( MOD ( itime , nprint ) .eq. 0 ) then  
              CALL write_thermo(itime, stdout)
              CALL write_thermo(itime, kunit_OUTFF)
        endif

        ! =============================================
        !  write instanteanous thermodynamic properties
        !  to file OSZIFF 
        ! =============================================
        if (  MOD ( itime , fprint ) .eq. 0 ) then  ! tmp
            CALL write_thermo(itime, kunit_OSZIFF)
        endif
  
        ! =====================
        !  save configuration 
        ! =====================
        if ( MOD ( itime , spas ) .eq. 0 ) then
          CALL write_CONTFF
        endif 
  enddo MAIN
  
  if ( ionode .and. .not. lstatic ) WRITE ( stdout , '(a,a)' ) &
  '===============================================================================================', &
  '======================================'
  if ( ionode ) WRITE ( stdout , '(a)' ) 'end of the main loop'
  if ( ionode .and. .not. lstatic ) WRITE ( kunit_OUTFF , '(a,a)' ) &
  '===============================================================================================', & 
  '======================================'
  if ( ionode ) WRITE ( kunit_OUTFF , '(a)' ) 'end of the main loop'
  if ( ionode ) WRITE ( stdout , '(a)' ) ' '
  if ( ionode ) WRITE ( stdout , '(a)' ) ' ' 
  if ( ionode ) WRITE ( stdout , '(a)' ) 'stress tensor of final configuration' 
  if ( ionode ) then
    CALL print_tensor ( tau_lj   , 'TAU_LJ  ' ) 
    CALL print_tensor ( tau_coul , 'TAU_COUL' ) 
  endif

  CALL  write_average_thermo ( stdout ) 
  CALL  write_average_thermo ( kunit_OUTFF ) 

  ! ===================
  !  properties output
  ! ===================

  ! =======
  ! lstrfac not working yet
  ! =======
!  if ( lstrfac ) then
!    CALL static_struc_fac (ngr)
!  endif 

  ! =======
  !  lmsd
  ! =======
  if ( lmsd ) then
    CALL msd_write_output ( 1 )
  endif

  ! =======
  !  lvacf
  ! =======
  if ( lvacf ) then
   CALL vacf_write_output 
  endif 

#ifdef multi_tau
    CALL multi_tau_write_output
    call dealloc 
#endif

#ifdef block
  CLOSE ( kunit_EQUILFF )
#endif
  CLOSE ( kunit_TRAJFF )
  CLOSE ( kunit_OSZIFF )

  deallocate( xtmp, ytmp, ztmp )

  return 

END SUBROUTINE md_run
! ===== fmV =====
