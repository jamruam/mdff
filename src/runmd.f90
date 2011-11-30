!for debugging purpose
!#define debug

! calculate the stress tensor at each time step
!#define stress_t 

!calculates efg at each time step
!#define efg_t

!===============================
!  experimental dev
!===============================

!calculate block_averaging 
!#define block

!multi_tau correlation
!#define multi_tau
!===============================

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

SUBROUTINE md_run ( iastart , iaend , list , point , offset )


  USE config,   ONLY :  natm , rx , ry , rz , rxs , rys , rzs , vx , vy , vz , write_CONTFF , box
  USE control,  ONLY :  lpbc , longrange , calc , lstatic , lvnlist , lbmlj , lcoulomb
  USE io_file,  ONLY :  ionode , stdout, kunit_OSZIFF, kunit_TRAJFF , kunit_EFGFF , kunit_EFGALL , kunit_OUTFF , kunit_EQUILFF
  USE prop,     ONLY :  lstrfac , lgr , lefg , lmsd , lvacf , nprop , nprop_start
  USE md,       ONLY :  npas , ltraj , lleapequi , itraj_period , itraj_start , nequil , nequil_period , nprint, &
                        fprint, spas , dt,  temp , updatevnl , write_traj_xyz , write_traj_xyz_test , integrator

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
  integer, intent(inout) :: iastart , iaend , list(250 * natm) , point(natm + 1)
  integer, intent(in) :: offset

  ! local
  double precision :: tempi , kin 
  double precision :: ttt1, ttt2 
  double precision, dimension(:), allocatable :: xtmp , ytmp , ztmp
  integer :: itime, i ,ierr
  integer :: nefg , ngr , nmsd , ntau
  character*60 :: rescale_allowed(3)
  data rescale_allowed / 'nve-vv' , 'nve-be' , 'nve-vv_test' /

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

  if( ionode ) WRITE ( stdout , '(a)' ) &
  '====================================================================================================================================='
  if( ionode ) WRITE ( kunit_OUTFF , '(a)' ) &
  '====================================================================================================================================='
  if( ionode ) WRITE ( stdout , '(a)' )      'properties at t=0'
  if( ionode ) WRITE ( kunit_OUTFF , '(a)' ) 'properties at t=0'
 
  allocate( xtmp(natm), ytmp(natm), ztmp(natm) )

  OPEN (unit = kunit_OSZIFF ,file = 'OSZIFF',STATUS = 'UNKNOWN')
  OPEN (unit = kunit_TRAJFF ,file = 'TRAJFF')
#ifdef block
  OPEN (unit = kunit_EQUILFF,file = 'EQUILFF',STATUS = 'UNKNOWN')
#endif

  ! =============
  !  lefg at t=0
  ! =============
  if ( lefg ) then
    if( longrange .eq. 'direct' )  CALL efg_DS ( 0 , 0 , iastart , iaend )
    if( longrange .eq. 'ewald'  )  CALL efg_ES ( 0 , 0 )
  endif
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
     CALL engforce ( iastart , iaend )!, list , point )

    ! ==================================================
    ! write thermodynamic information of config at t=0
    ! ==================================================
    CALL write_thermo( offset-1 , stdout )
    CALL write_thermo( offset-1 , kunit_OUTFF )
    CALL write_thermo( offset-1 , kunit_OSZIFF )

  else
  ! ===========================================================================
  !  leapfrog algorithm set the first leap value for r(t-dt)  (approximation)
  ! ===========================================================================
    xtmp = rx
    ytmp = ry
    ztmp = rz
    do i = 1 , natm
      rxs ( i )  = rx ( i ) - vx ( i ) * dt
      rys ( i )  = ry ( i ) - vy ( i ) * dt
      rzs ( i )  = rz ( i ) - vz ( i ) * dt
    enddo
    rx = rxs 
    ry = rys
    rz = rzs 
    ! ===================
    ! force + potential
    ! ===================
    CALL engforce ( iastart , iaend )!, list , point )
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
  if( ionode ) WRITE ( stdout , '(a)' ) ' ' 
  if( ionode ) WRITE ( stdout , '(a)' ) 'stress tensor of initial configuration' 
  if ( lbmlj )    CALL stress_bmlj ( iastart , iaend , list , point )
  if ( lcoulomb ) CALL stress_coul 

  ! =========================
  !   MAIN LOOP ( TIME )
  ! =========================
  if( ionode ) WRITE ( stdout , '(a)' ) ''
  if( ionode ) WRITE ( stdout , '(a)' ) 'starting main loop'
  if( ionode ) WRITE ( stdout , '(a)' ) &
  '====================================================================================================================================='
  if( ionode ) WRITE ( kunit_OUTFF , '(a)' ) ''
  if( ionode ) WRITE ( kunit_OUTFF , '(a)' ) 'starting main loop'
  if( ionode ) WRITE ( kunit_OUTFF , '(a)' ) &
  '====================================================================================================================================='
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
           do i = 1, natm
             rxs ( i )  = rx ( i ) - vx ( i ) * dt
             rys ( i )  = ry ( i ) - vy ( i ) * dt
             rzs ( i )  = rz ( i ) - vz ( i ) * dt
           enddo
         endif
         
         ! ===========
         !  time info
         ! ===========
         ttt1 = MPI_WTIME(ierr)

         ! =========================    
         !  integration t -> t + dt 
         ! =========================
         if ( integrator.eq.'nve-lf'  )         CALL prop_leap_frog ( iastart , iaend )!, list , point )
         if ( integrator.eq.'nve-be'  )         CALL beeman ( iastart , iaend )!, list , point )
         if ( integrator.eq.'nve-vv'  )         CALL prop_velocity_verlet ( iastart , iaend )!, list , point )
         if ( integrator.eq.'nvt-and' )         CALL prop_velocity_verlet ( iastart , iaend )!, list , point )
         if ( integrator.eq.'nvt-nhc2')         CALL nose_hoover_chain2 ( iastart , iaend )!, list , point )
         if ( integrator.eq.'nve-vv_test' )     CALL prop_velocity_verlet_test ( iastart , iaend )!, list , point )

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
         if( (any(integrator.eq.rescale_allowed)).and.((itime.le.nequil.and.mod(itime,nequil_period).eq.0).and.(itime.ne.npas+(offset-1).and.itime.ne.offset)) ) then
           CALL rescale_velocities(0)
         endif

         ! ===================================
         !  rescale velocities (NVT Andersen)
         ! ===================================
         if( (integrator.eq.'nvt-and') ) CALL andersen_velocities

         ! ===================
         !  print trajectory 
         ! ===================
         if( ltraj .and. (itime .gt. itraj_start ) .and. mod(itime,itraj_period) .eq. 0 ) then
           xtmp = rx
           ytmp = ry
           ztmp = rz
#ifdef debug
           CALL write_traj_xyz_test
#endif 
           CALL write_traj_xyz
           rx = xtmp
           ry = ytmp
           rz = ztmp
         endif
#ifdef stress_t
  if ( lbmlj )    CALL stress_bmlj ( iastart , iaend , list , point )
  if ( lcoulomb ) CALL stress_coul 
#endif
         ! =======================
         !  properties on-the-fly
         ! =======================
         ! -----------------------------------------------------------------------------------------
         if( mod(itime,nprop) .eq. 0 .and. itime.gt.nprop_start) then 

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
  
           ! =======
           !  lefg
           ! =======
           if ( lefg ) then
             nefg = nefg + 1
             ! ==================
             !  direct summation 
             ! ================== 
             if ( longrange .eq. 'direct' )  CALL efg_DS ( itime , nefg , iastart , iaend )
             ! ==================
             !  ewald summation 
             ! ================== 
             if ( longrange .eq. 'ewald' )   CALL efg_ES ( itime , nefg )
#ifdef efg_t 
             ! =============
             !  efg output 
             ! ============= 
             CALL efg_write_output( nefg )
#endif efg_t
           endif 

           ! =======
           !  lgr
           ! =======
           if ( lgr ) then
             ngr = ngr + 1 
             CALL gr_main ( iastart , iaend , ngr ) 
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
        if ( mod( itime , nprint ) .eq. 0 ) then  
              CALL write_thermo(itime, stdout)
              CALL write_thermo(itime, kunit_OUTFF)
        endif

        ! =============================================
        !  write instanteanous thermodynamic properties
        !  to file OSZIFF 
        ! =============================================
        if (  mod( itime , fprint ) .eq. 0 ) then  ! tmp
            CALL write_thermo(itime, kunit_OSZIFF)
        endif
  
        ! =====================
        !  save configuration 
        ! =====================
        if ( mod( itime , spas ) .eq. 0 ) then
          CALL write_CONTFF
        endif 
  enddo MAIN
  
  if( ionode .and. .not. lstatic ) WRITE ( stdout , '(a)' ) &
  '====================================================================================================================================='
  if( ionode ) WRITE ( stdout , '(a)' ) 'end of the main loop'
  if( ionode .and. .not. lstatic ) WRITE ( kunit_OUTFF , '(a)' ) &
  '====================================================================================================================================='
  if( ionode ) WRITE ( kunit_OUTFF , '(a)' ) 'end of the main loop'
  if( ionode ) WRITE ( stdout , '(a)' ) ' '
  if( ionode ) WRITE ( stdout , '(a)' ) 'stress tensor of final configuration'
  if( ionode ) WRITE ( kunit_OUTFF , '(a)' ) ' '
  if( ionode ) WRITE ( kunit_OUTFF , '(a)' ) 'stress tensor of final configuration'

  if ( lbmlj )    CALL stress_bmlj ( iastart , iaend , list , point )
  if ( lcoulomb ) CALL stress_coul 


  CALL  write_average_thermo ( stdout ) 
  CALL  write_average_thermo ( kunit_OUTFF ) 

  ! ===================
  !  properties output
  ! ===================

  ! =======
  ! lstrfac not working yet
  ! =======
  if ( lstrfac ) then
    CALL static_struc_fac (ngr)
  endif 

  ! =======
  !  lefg
  ! =======
  if ( lefg ) then
    if ( ionode ) WRITE ( stdout , * ) 'write efg distribution files of ',nefg,' configurations'
    CALL efg_write_output( nefg ) 
  endif

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
#endif multi_tau

#ifdef block
  CLOSE ( kunit_EQUILFF )
#endif
  CLOSE ( kunit_TRAJFF )
  CLOSE ( kunit_OSZIFF )

  deallocate( xtmp, ytmp, ztmp )

  return 

END SUBROUTINE md_run

! ===== fmV =====
