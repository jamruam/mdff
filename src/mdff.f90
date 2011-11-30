!#define MULTRUN
!#define block
! ===== fmV =====
! ===============================================================
! *      PARALLEL MOLECULAR DYNAMICS ...for fun                 *
! *      Version: 0.2.1 BMLJ  (q-p) LJ potential                *
! *      filipe.manuel.vasconcelos@gmail.com                    *
! *      Based on different code or books by:                   *
! *      F. Affouard University of Lille                        *
! *      S. Sastry   JNCASR                                     *
! *      D. Frenkel and B. Smit   (book)                        *
! *      Allen-Tildesley (book)                                 *
! *      vasp: some features of the best code ever ;)           *
! *      dea-dec2007-janv2011                                   *
! ===============================================================
!
! ===============================================================
! CHECK (optimise default value):
!
! * nuandersen
! * skindiff
! * Qnosehoover
! * tauberendsen  
!
! ===============================================================
! 
! TODO when programming:
!
! * allocate local arrays
! * let some space in the operation at least ( var = var *...)
!  -> better to follow when using grep or any script
!
! ===============================================================
!
! TODO:  
! * add more comments ( in/out of subroutine and some comments)
! * memory usage ?
! * TRAJFF format BINARY 
! * calculate average values (by block) remind use nequil                       
! * implement nvt-nh 
! * add WARNINGS for some tags where default value (???)
! * add default values for secondary tags
! * add some descriptions of gaussian distribution routine 
!  -> (basically some refs to be able to find bugs)
! * make the force field easier to read in control.F
! * parallelization of reciprocal part (ewald sum) in efg at least
!
! not before version 0.3.0:
! * ntype > 2 ( the complicated part will be how to read the force field )
! * masses   
! * write some manual (with theory, ref teaching oriented)
!
! * remove the TODO list ;)                                        
!
! ===============================================================


! ===============================================================
! MAIN SUBROUTINE :                            
! set, read and check control parameters
! timing information at the end
! ===============================================================

PROGRAM main_MDFF

  USE constants
  USE io_file
  USE control 
  USE config
  USE time
  USE md 
  USE field 
  USE vib
  USE opt
  USE prop
  USE efg
  USE radial_distrib 
  USE msd
  USE vacf 
  USE block

  implicit none
  INCLUDE 'mpif.h'

  ! local
  integer :: argc
  integer :: ierr 
  integer :: iastart , iaend   ! atom decomposition
  integer :: offset         
  double precision :: ttt2, ttt1
#ifdef MULTRUN
  integer :: irun          
#endif
  character*60 :: vib_allowed(5)                     
  data vib_allowed / 'vib' , 'vib+fvib' , 'vib+gmod' , 'vib+band' , 'vib+dos' / 
! tmp removed 

  ! ====================================================
  ! MPI initialisation
  ! ====================================================
  CALL MPI_INIT ( ierr )                           
  ! ====================================================
  ! time information init 
  ! ====================================================
  CALL time_init
  ttt1 = MPI_WTIME(ierr)
  ! ====================================================      
  ! random number generator init
  ! ====================================================      
  ! CALL init_random_seed   ! MD is too sensitive to the initial velocities
  ! ====================================================
  ! gives rank index (myrank) 
  ! ====================================================
  CALL MPI_COMM_RANK ( MPI_COMM_WORLD , myrank , ierr )  
  ! ====================================================
  ! gives total proc number  
  ! ====================================================
  CALL MPI_COMM_SIZE ( MPI_COMM_WORLD , numprocs , ierr )  
  ! ====================================================
  ! input/output init 
  ! ====================================================
  CALL io_init                                            
  ! ====================================================
  ! VERSION 
  ! ====================================================
  VERSION = '0.2.0'
  
  ! ========================================
  ! reads command line
  ! ========================================
  argc = iargc()
  if ( argc .ne. 1 ) then
    if (  ionode  ) WRITE ( stdout , * ) '* command: mdff.x <input_file>'
    STOP 
  endif

  ! ========================================
  ! set some default values to control tags 
  ! parameters before to read them 
  ! and check consistency
  ! ========================================
  CALL control_init

  ! =====================================
  ! bmlj_init is the main PROGRAM: 
  ! only binary mixture lennard-jones
  ! particules are avaiable different 
  ! properties EFG, GR ... 
  ! are calculated inside this routine  
  ! =====================================
  if ( .not. ltest )  then

    ! =============================
    ! md initialization
    ! =============================
    CALL md_init
 
    ! =============================
    ! configuration initialization
    ! =============================
    CALL config_init 

    ! =============================
    ! force field initialization 
    ! ============================= 
    if ( calc .eq. 'md' ) CALL field_init

    ! =============================
    ! properties initialization	 
    ! =============================
    CALL prop_init

    ! ===========================================================
    ! parallelization: atom decomposition split atom into
    ! different procs. iastart , iaend are the index of atoms 
    ! for rank = myrank.  only for calc='md'. 
    ! for the other type of calc this is done when natm is known
    ! ===========================================================
    if ( calc .eq. 'md' ) CALL do_split ( natm , myrank , numprocs , iastart , iaend )

    ! =====================================
    ! verlet list generation
    ! mmmmh only md ? 
    ! vnlist should be test with opt too 
    ! =====================================
    if( calc .eq. 'md' ) then 
      if ( ( lvnlist ) .and. lpbc )           CALL vnlist_pbc ( iastart , iaend , list , point )
      if ( ( lvnlist ) .and. ( .not. lpbc) )  CALL vnlist_nopbc   ( iastart , iaend , list , point )       
    endif  
    !IF MD
    MDLOOP: if( calc .eq. 'md') then

    ! ===============================================
    ! some initialization of on-the-fly properties
    ! ===============================================

          ! =====================================
          ! electric field gradient        
          ! =====================================
          CALL efg_init
          CALL efg_alloc

          ! =====================================
          ! radial distribution function
          ! =====================================
          CALL gr_init
          CALL gr_alloc   

          ! =====================================
          ! mean square displacement
          ! =====================================
          CALL msd_init
          CALL msd_alloc 

          ! =====================================
          ! velocity auto-correlation function
          ! =====================================
          CALL vacf_init
          CALL vacf_alloc

    ! ==================================================
    ! OUTSIDE LOOP (not USE yet)
    ! can be useful for some applications #version 0.2
    ! temperature variations or any other parameter
    ! as only one parameters will variate with time
    ! try to find a condensate way to write it
    ! ==================================================
#ifdef MULTRUN
    RUN: DO irun = 1, 1
    tempi = temp
    RUN: DO irun = 0 , 2
        temp = ( 0.1d0 * irun  ) + tempi
        if(  ionode  ) WRITE ( stdout ,'(a)') ''
        if(  ionode  ) WRITE ( stdout , '(a,f8.5)' ) 'T = ' , temp
#endif
  
         ! =======================
         !        static 
         ! =======================
           if( lstatic ) then  
             offset = 1
             npas = 0
             integrator = 'nve-vv' 
             CALL md_run ( iastart , iaend , list , point , offset )
         ! =======================
         ! .... or dynamic   
         ! =======================
           else               
             ! =========================
             ! velocities intialization
             ! =========================
             CALL init_velocities 
             ! ===============================================================
             ! time offset = zero is the static step (first step of dynamic) 
             ! this offset will be USEd if the OUTSIDE LOOP is USEd
             ! ===============================================================
!             offset = ( ( irun - 1 ) * npas ) + 1  
             CALL md_run ( iastart , iaend , list , point , 1 ) 
!                                                           ^ 
!                                                          offset   
           endif
   
#ifdef MULTRUN
        ENDDO RUN
#endif
  
    ENDIF MDLOOP

    ! todo add some comments on opt vib efg here
     
    ! ==============================================
    ! IF OPT :
    ! optimisation of structure read in TRAJFF 
    ! ==============================================
    if ( calc .eq. 'opt' ) then  
      CALL opt_init
      CALL opt_main 
    endif
  
    ! ==============================================
    ! IF VIB :
    ! ==============================================
    if ( any( calc .eq. vib_allowed ) ) then       
      CALL vib_init
      CALL vib_main 
    endif 
  
    ! ==============================================
    ! IF EFG : 
    !  calculates EFG of structures read in TRAJFF   
    ! ==============================================
    if ( calc .eq. 'efg') then
      CALL efg_init
      CALL efgcalc 
    endif
    if ( calc .eq. 'efg+acf') then
      CALL efg_init
      CALL efg_acf 
    endif
   
    ! ==============================================
    ! IF GR : 
    ! calculates radial distribution  of structures 
    ! read in TRAJFF   
    ! ==============================================
    if ( calc .eq. 'gr') then
      CALL grcalc 
    endif

    ! ========================================================
    ! write final config pos and vel (always) only for md run
    ! ========================================================
    if( calc.eq. 'md' ) then 
      CALL write_CONTFF
    endif
 
    ! ==============================================
    ! deallocate prop quantities
    ! ==============================================
    CALL efg_dealloc
    CALL gr_dealloc 
    CALL msd_dealloc
    CALL vacf_dealloc

    ! ==============================================
    ! deallocate coulombic related quantities
    ! ==============================================
    CALL finalize_coulomb

    ! ==============================================
    ! deallocate principal quantites
    ! ==============================================
    CALL config_dealloc

#ifdef block
    CALL block_
#endif

  ! =======================
  !  only for test purpose
  ! =======================
  else

    if ( ionode ) WRITE ( stdout ,'(a)') ' ############################# '
    if ( ionode ) WRITE ( stdout ,'(a)') ' ### TEST PART OF THE CODE ###'
    if ( ionode ) WRITE ( stdout ,'(a)') ' ############################# '


  endif


  if (ionode ) then
    WRITE ( stdout ,'(a)')      '============================================================='
    WRITE ( stdout ,'(a)')      ''
    WRITE ( kunit_OUTFF ,'(a)') '============================================================='
    WRITE ( kunit_OUTFF ,'(a)') ''
  endif
  ttt2 = MPI_WTIME(ierr)
  timetot = ttt2 - ttt1
  CALL print_time_info ( stdout ) 
  CALL print_time_info ( kunit_OUTFF ) 
  CALL io_end

  CALL MPI_FINALIZE(ierr)

  STOP

END PROGRAM main_MDFF
! ===== fmV =====
