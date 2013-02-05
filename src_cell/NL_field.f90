!#define debug
SUBROUTINE Nymand_and_Linse_test

  USE config,   ONLY :  system , natm , ntype , atype , rx , ry , rz , fx , fy ,fz , itype , &
                        atypei , natmi, rho , simu_cell , config_alloc , phi_coul_tot , qia , dipia , ipolar 
                       
  USE control,  ONLY :  longrange , myrank , numprocs , cutlongrange 
  USE io_file      , ONLY : ionode , stdout , kunit_OUTFF , kunit_POSFF, kunit_EFGALL, kunit_NMRFF 
  USE thermodynamic, ONLY : u_coul_tot , vir_coul_tot , u_pol
  USE efg,      ONLY : nmr_convention
  USE field,    ONLY : ef_t , efg_t , qch , dip , lpolar , field_init , rm_coul , km_coul , initialize_coulomb , finalize_coulomb,  &
                       alphaES , multipole_DS , multipole_ES 
  USE cell

  implicit none

  INCLUDE 'mpif.h' 

  ! local
  integer :: ia , it , ifail 
  integer :: iastart , iaend
  integer, parameter :: lwork = 6
  double precision , dimension(:,:)   , allocatable :: ef_tmp
  double precision , dimension(:,:,:) , allocatable :: efg_tmp
  double precision , dimension(:,:)   , allocatable :: mu_tot
  double precision , dimension(:,:)   , allocatable :: mu_ind
  double precision , dimension (:,:)  , allocatable :: nmr
  double precision                                  :: work(3 * lwork)
  double precision                                  :: nmr_conv ( 4 ) , w ( 3 ) , tmp ( 3 , 3 )

  OPEN (UNIT = kunit_POSFF , FILE = 'POSFF')
  OPEN (UNIT = kunit_NMRFF , FILE = 'NMRFF')
  OPEN (UNIT = kunit_EFGALL, FILE = 'EFGALL')

  READ ( kunit_POSFF, * )  natm 
  READ ( kunit_POSFF, * )  system
  READ ( kunit_POSFF ,* ) simu_cell%A ( 1 , 1 ) , simu_cell%A ( 2 , 1 ) , simu_cell%A ( 3 , 1 )
  READ ( kunit_POSFF ,* ) simu_cell%A ( 1 , 2 ) , simu_cell%A ( 2 , 2 ) , simu_cell%A ( 3 , 2 )
  READ ( kunit_POSFF ,* ) simu_cell%A ( 1 , 3 ) , simu_cell%A ( 2 , 3 ) , simu_cell%A ( 3 , 3 )
  READ ( kunit_POSFF, * )  ntype
  READ ( kunit_POSFF ,* ) ( atypei ( it ) , it = 1 , ntype )
  IF ( ionode ) WRITE ( kunit_OUTFF ,'(A,20A3)' ) 'found type information on POSFF : ', atypei ( 1:ntype )
  IF ( ionode ) WRITE ( stdout      ,'(A,20A3)' ) 'found type information on POSFF : ', atypei ( 1:ntype )
  READ( kunit_POSFF ,*)   ( natmi ( it ) , it = 1 , ntype )

  CALL lattice ( simu_cell )
  natm = 0
  do it = 1 , ntype
    natm = natm + natmi(it)
  enddo
  rho = DBLE ( natm ) / simu_cell%omega

  allocate ( ef_tmp    ( natm , 3 ) ) 
  allocate ( efg_tmp    ( natm , 3 , 3 ) ) 
  allocate ( mu_tot ( natm , 3 ) ) 
  allocate ( mu_ind ( natm , 3 ) )  
  allocate ( nmr ( natm , 4 ) ) 

  ef_tmp = 0.0d0
  mu_tot = 0.0d0
  ! ===================================
  !  here we know natm, then alloc 
  !  and decomposition can be applied 
  ! ================================== 
  CALL config_alloc
  CALL do_split ( natm , myrank , numprocs , iastart , iaend )
  ! read charge in fieldtag
  CALL field_init
  CALL initialize_coulomb
  ! =============
  !  print info
  ! =============
  CALL NL_field_print_info(stdout)
  CALL NL_field_print_info(kunit_OUTFF)

  CALL print_general_info( stdout )
  CALL print_general_info( kunit_OUTFF )

  do ia = 1 , natm
    READ ( kunit_POSFF, * ) atype ( ia ) , rx ( ia ) , ry ( ia ) , rz ( ia )
  enddo

  CALL typeinfo_init 

  CALL distance_tab ( stdout )
 
  ! =============================================================================
  !  start calculation of potential energy , field , efg , virial , forces ...
  ! =============================================================================
  u_coul_tot   = 0.0d0
  phi_coul_tot = 0.0d0
  vir_coul_tot = 0.0d0
  ef_t         = 0.0d0
  efg_t        = 0.0d0
  fx           = 0.0d0
  fy           = 0.0d0
  fz           = 0.0d0
 
  ! ===============================================
  !  get moment from polarisabilities
  ! ===============================================
  CALL moment_from_pola ( iastart, iaend , mu_ind )

  mu_tot = dipia + mu_ind
  ! =============================================
  !  coulombic energy , forces (field) and virial
  ! =============================================
  if ( longrange .eq. 'direct' ) CALL multipole_DS ( iastart, iaend , ef_tmp , efg_tmp , mu_tot , u_coul_tot , vir_coul_tot , phi_coul_tot )
  if ( longrange .eq. 'ewald'  ) CALL multipole_ES ( iastart, iaend , ef_tmp , efg_tmp , mu_tot , u_coul_tot , vir_coul_tot , phi_coul_tot )

  ef_t  = ef_t  + ef_tmp
  efg_t = efg_t + efg_tmp


  if ( longrange .eq. 'direct' ) CALL  Nymand_Linse_output ( 'DS' , '  ', u_coul_tot , phi_coul_tot , ef_t , efg_t , fx , vir_coul_tot , natm )
  if ( longrange .eq. 'ewald' )  CALL  Nymand_Linse_output ( 'ES' , ' 1', u_coul_tot , phi_coul_tot , ef_t , efg_t , fx , vir_coul_tot , natm )
  ! =======================================
  ! write efg for each atom in file EFGALL 
  ! =======================================
  if ( ionode  ) then
    WRITE ( kunit_EFGALL , * )  natm
    WRITE ( kunit_EFGALL , * )  system
    WRITE ( kunit_EFGALL , * )  simu_cell%A(1,1) , simu_cell%A(1,2) , simu_cell%A(1,3)
    WRITE ( kunit_EFGALL , * )  simu_cell%A(2,1) , simu_cell%A(2,2) , simu_cell%A(3,3)
    WRITE ( kunit_EFGALL , * )  simu_cell%A(3,1) , simu_cell%A(3,2) , simu_cell%A(3,3)
    WRITE ( kunit_EFGALL , * )  ntype
    WRITE ( kunit_EFGALL , * )  ( atypei ( it ) , it = 1 , ntype )
    WRITE ( kunit_EFGALL , * )  ( natmi  ( it ) , it = 1 , ntype )
    WRITE ( kunit_EFGALL ,'(a)') &
    '      ia type                   vxx                   vyy                vzz                   vxy                   vxz                   vyz'
    do ia = 1 , natm
      WRITE ( kunit_EFGALL ,'(i8,2x,a3,6f26.16)') &
      ia , atype ( ia ) , efg_t ( ia , 1 , 1) , efg_t ( ia , 2 , 2) , &
                          efg_t ( ia , 3 , 3) , efg_t ( ia , 1 , 2) , &
                          efg_t ( ia , 1 , 3) , efg_t ( ia , 2 , 3)

    enddo
  endif

  do ia = 1 , natm 
    ! =================
    !  diagonalisation
    ! =================
    tmp = efg_t(ia , : , : )   
    CALL DSYEV ( 'N' , 'U' , 3 , tmp  , 3 , w , work , 3 * lwork , ifail )
    if ( ifail .ne. 0 ) then
      if ( ionode ) WRITE ( stdout , * ) 'ERROR: DSYEV, STOP in efg MODULE (ewald)'
      STOP
    endif

    ! =============================
    ! NMR conventions test 
    ! =============================
    CALL nmr_convention( w , nmr_conv ( : )  , ia )
    nmr( ia , : )  = nmr_conv ( : )
  enddo
  ! ========================================
  !  output NMR  parameters 
  ! ========================================
  if ( ionode ) then
      WRITE ( kunit_NMRFF , '(a)' ) &
      '#     site           xx          yy          zz               eta'
    do ia = 1 , natm
      WRITE ( kunit_NMRFF , 150 ) ia , atype ( ia ) , &
                                      nmr( ia , 1 ) , nmr( ia , 2 ) , &
                                      nmr( ia , 3 ) , nmr( ia , 4 )
    enddo
  endif

  CLOSE (UNIT = kunit_EFGALL )
  CLOSE (UNIT = kunit_NMRFF )
  CLOSE (UNIT = kunit_POSFF )

  deallocate ( ef_tmp ) 
  deallocate ( efg_tmp ) 
  deallocate ( mu_tot ) 
  deallocate ( mu_ind ) 
  deallocate ( nmr ) 

150 FORMAT(I8,1X,A4,3F20.6,F22.6)

  return

END SUBROUTINE Nymand_and_Linse_test


SUBROUTINE Nymand_Linse_output ( label , labelES , U , pot , Efield , EFG , fx , vir , natm )

  USE io_file,  ONLY : stdout , ionode
  USE constants, ONLY : e_2 , coulombic_factor

  implicit none

  ! global
  integer :: natm
  double precision :: U , vir
  double precision :: pot ( natm ) , fx ( natm )
  double precision :: Efield ( natm , 3 ) , EFG ( natm , 3 , 3 )
  character(len=2) label , labelES

  if ( ionode ) then
    WRITE ( stdout , '(a)' ) ''
    WRITE ( stdout , 490 ) 'Method  ',' surf ',' U ',' phi_1 ',' E1,x ', 'E1,xx ',' E1,yy ',' f1,x ',' PHI '
    WRITE ( stdout , 500 ) label , labelES, U , pot(1) , Efield ( 1 , 1 ) , EFG ( 1 , 1 , 1 ) , EFG ( 1 , 2 , 2 ) , fx ( 1 ) , vir
    WRITE ( stdout , 500 ) ' (eV) ',labelES, U*e_2 , pot(1)*e_2 , Efield ( 1 , 1 )*e_2 , EFG ( 1 , 1 , 1 )*e_2 , EFG ( 1 , 2 , 2 )*e_2 , fx ( 1 ) * coulombic_factor , vir*e_2
    WRITE ( stdout , '(a)' ) ''
  endif

490 FORMAT(A8,A8,7A20)
500 FORMAT(A8,A8,7F20.7)

  return

END SUBROUTINE Nymand_Linse_output

SUBROUTINE NL_field_print_info ( kunit ) 

  USE io_file,  ONLY :  ionode
  USE control,  ONLY :  calc , longrange , cutlongrange
  USE config,   ONLY :  simu_cell 
  USE field,    ONLY :  ncelldirect , kES , alphaES 

  implicit none

  ! global 
  integer :: kunit
  ! local
  double precision :: aaa , i 

  if ( ionode ) then
      WRITE ( kunit ,'(a)')           '============================================================='
      WRITE ( kunit ,'(a)')           ''
      WRITE ( kunit ,'(a)')           'Nymand and Linse test simple systems : '
      WRITE ( kunit ,'(a)')           'atomic charge, dipole and polarisabilities'
      WRITE ( kunit ,'(a)')           'read config from file            : POSFF'
  endif
 
  return

END SUBROUTINE NL_field_print_info

SUBROUTINE moment_from_pola ( iastart , iaend , mu_ind ) 

  USE io_file, ONLY : ionode , stdout 
  USE config,  ONLY : natm , atype , fx , fy , fz , ntype , dipia , qia , ntypemax
  USE control, ONLY : longrange
  USE field,   ONLY : qch , lpolar , multipole_DS , multipole_ES , induced_moment
  USE thermodynamic, ONLY : u_pol

  implicit none

  ! global
  integer :: iastart , iaend 
  double precision :: mu_ind ( natm , 3 ) 

  ! local
  integer :: ia , iscf , it 
  logical :: linduced
  double precision :: Efield_old , diff_efield , conv_tol
  double precision :: u_tmp , vir_tmp
  double precision :: Efield( natm , 3 ) , Efield_stat ( natm , 3 ) , Efield_ind ( natm , 3 ) 
  double precision :: ef_tmp( natm , 3 ) , efg_tmp ( natm , 3 , 3 ) , phi_tmp ( natm ) 
  double precision :: qia_tmp ( natm )  , qch_tmp ( ntypemax ) 

  ! =========================================================
  !  Is there any polarizability ? if yes linduced = .TRUE.
  ! =========================================================
  linduced = .false.
  do it = 1 , ntype
    if ( lpolar ( it ) .eq. 1 ) linduced = .true.
  enddo
  if ( ionode .and. .not. linduced ) then
    write ( stdout , '(a)' ) 'No polarizabilities'
    return
  endif

  ! =============================================
  !  calculate static Efield ( charge + dipoles )
  ! =============================================
  Efield_stat = 0.0d0

  ! =============================================
  !  coulombic energy , forces (field) and virial
  ! =============================================
  if ( longrange .eq. 'direct' ) CALL  multipole_DS ( iastart, iaend , ef_tmp , efg_tmp , dipia , u_tmp , vir_tmp , phi_tmp ) 
  if ( longrange .eq. 'ewald' )  CALL  multipole_ES ( iastart, iaend , ef_tmp , efg_tmp , dipia , u_tmp , vir_tmp , phi_tmp ) 

  Efield_stat = Efield_stat + ef_tmp
  !print*,Efield_stat

  fx      = 0.0d0 ; fy = 0.0d0 ; fz = 0.0d0
  ef_tmp  = 0.0d0
  efg_tmp = 0.0d0
  u_tmp   = 0.0d0
  vir_tmp = 0.0d0
  phi_tmp = 0.0d0

  ! =============================================
  !  init total Efield to static only
  ! =============================================
  Efield = Efield_stat

  if ( ionode ) WRITE ( stdout , '(a)' ) 'calculate induced dipole SCF'

  if ( ionode ) then
    WRITE ( stdout , '(a)' ) 'We start from the static electric field (charge + static dipole)'
    WRITE ( stdout , '(a)' ) ' '
  endif
  diff_efield = 1000.0d0
  Efield_old  = Efield(1,1)
  iscf = 0
  ! =============================
  !           SCF LOOP
  ! =============================
  conv_tol = 0.0007 !hardware
  do while ( diff_efield .gt. conv_tol )

    iscf = iscf + 1

    ! ==========================================================
    !  calculate mu_ind from Efield = Efield_stat + Efield_ind
    ! ==========================================================
    CALL induced_moment ( Efield , mu_ind , u_pol )  ! Efield in ; mu_ind and u_pol out

#ifdef debug
   do ia =1 , natm
     write( stdout , '(a,3f12.5)' ) 'debug : induced moment from pola', mu_ind ( ia , 1 ),  mu_ind ( ia , 2 ) , mu_ind ( ia , 3 )
   enddo
   do ia =1 , natm
     write( stdout , '(a,3f12.5)' ) 'debug : Efield                  ', Efield ( ia , 1 ),  Efield ( ia , 2 ) , Efield ( ia , 3 )
   enddo
#endif

    ! ==========================================================
    !  calculate Efield_ind from mu_ind
    !  Efield_ind out , mu_ind in ==> charges and static dipoles = 0
    ! ==========================================================
    qch_tmp = qch
    qia_tmp = qia
    qch = 0.0d0
    qia = 0.0d0
    if ( longrange .eq. 'direct' ) CALL multipole_DS ( iastart, iaend , Efield_ind , efg_tmp , mu_ind , u_tmp , vir_tmp , phi_tmp ) 
    if ( longrange .eq. 'ewald' )  CALL multipole_ES ( iastart, iaend , Efield_ind , efg_tmp , mu_ind , u_tmp , vir_tmp , phi_tmp ) 


    qch = qch_tmp
    qia = qia_tmp
    fx = 0.0d0 ; fy = 0.0d0 ; fz = 0.0d0
    Efield = Efield_stat + Efield_ind

    ! ===================
    !  stopping criteria
    ! ===================
    diff_efield = ABS ( Efield(1,1) - Efield_old )
    Efield_old = Efield(1,1)

    if ( ionode ) WRITE ( stdout ,'(a,i4,a,3f18.10,2(a,f18.10))') &
    'scf = ',iscf,' Efield atom 1= ',Efield(1,1),Efield(1,2),Efield(1,3),' conv = ',diff_efield,' u_pol = ',u_pol

  enddo ! end of SCF loop



  if ( ionode ) then
    WRITE ( stdout , '(a)' ) ' '
    WRITE ( stdout , '(a,i6,a)') 'scf calculation of the induced electric moment converged in ',iscf, ' iterations '
    WRITE ( stdout , '(a,e10.3)') 'Electric field is converged within ',conv_tol
    WRITE ( stdout , '(a)' ) ' '
    WRITE ( stdout , '(a)' )     'Induced dipoles at atoms : '
    do ia = 1 , natm
      WRITE ( stdout , '(i5,a3,a,3f18.10)' ) &
      ia,atype(ia),' mu_ind = ', mu_ind ( ia , 1 ) , mu_ind ( ia , 2 ) , mu_ind ( ia , 3 )
    enddo
    WRITE ( stdout , '(a)' ) ' '
  endif


  return

END SUBROUTINE moment_from_pola


