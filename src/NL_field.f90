!#define debug
SUBROUTINE Nymand_and_Linse_test

  USE config,   ONLY :  system , natm , ntype , atype , rx , ry , rz , fx , fy ,fz , itype , &
                        atypei , natmi, rho , box , omega , config_alloc , phi_coul_tot , phi_coul_qq , phi_coul_dd , qia , dipia , ipolar 
                       
  USE control,  ONLY :  longrange , myrank , numprocs , cutlongrange 
  USE io_file,  ONLY :  ionode , stdout , kunit_OUTFF , kunit_POSFF
  USE thermodynamic, ONLY : u_coul_tot , vir_coul_tot , u_coul_qq, u_coul_qd, u_coul_dd , vir_coul_qq, vir_coul_qd , vir_coul_dd , u_pol
  USE efg,      ONLY : get_efg_tensor 
  USE field,    ONLY : ef_t , efg_t , qch , dip , lpolar , field_init , rm_coul , km_coul , initialize_coulomb , finalize_coulomb,  &
                       engforce_charge_DS , engforce_charge_ES , engforce_dipole_DS , engforce_dipole_ES , engforce_charge_and_dipole_DS, alphaES , multipole_DS , multipole_ES ! , &
                       !engforce_charge_and_dipole_ES


  implicit none

  ! local
  integer :: ia , it , cc , ccs , imethod
  integer :: iastart , iaend
  double precision , dimension(:,:)   , allocatable :: ef_tmp
  double precision , dimension(:,:,:) , allocatable :: efg_tmp
  double precision , dimension(:,:)   , allocatable :: mu_tot
  double precision , dimension(:,:)   , allocatable :: mu_ind

  OPEN (UNIT = kunit_POSFF , FILE = 'POSFF')

  READ ( kunit_POSFF, * )  natm 
  READ ( kunit_POSFF, * )  system
  READ ( kunit_POSFF, * )  box,ntype
  READ ( kunit_POSFF ,* ) ( atypei ( it ) , it = 1 , ntype )
  IF ( ionode ) WRITE ( kunit_OUTFF ,'(A,20A3)' ) 'found type information on POSFF : ', atypei ( 1:ntype )
  IF ( ionode ) WRITE ( stdout      ,'(A,20A3)' ) 'found type information on POSFF : ', atypei ( 1:ntype )
  READ( kunit_POSFF ,*)   ( natmi ( it ) , it = 1 , ntype )
  omega = box * box * box
  natm = 0
  do it = 1 , ntype
    natm = natm + natmi(it)
  enddo
  rho = DBLE ( natm ) / omega

  allocate ( ef_tmp    ( natm , 3 ) ) 
  allocate ( efg_tmp    ( natm , 3 , 3 ) ) 
  allocate ( mu_tot ( natm , 3 ) ) 
  allocate ( mu_ind ( natm , 3 ) ) 

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

  ! ==========================
  !  set some type parameters
  ! ==========================      
  natmi ( 0 ) = 0
  cc = 0
  do it = 1 , ntype
    ccs = cc
    cc = cc + natmi ( it )
    do ia = ccs + 1 , cc
      atype ( ia ) = atypei ( it )
      itype ( ia ) = it
      qia   ( ia ) = qch (it)
      dipia ( ia , 1 ) = dip ( it , 1 )
      dipia ( ia , 2 ) = dip ( it , 2 )
      dipia ( ia , 3 ) = dip ( it , 3 )
      ipolar ( ia )    = lpolar ( it )
    enddo
  enddo
  natmi  ( 0 ) = natm
  atypei ( 0 ) = 'ALL'
 
#ifdef debug
  write( stdout , '(a,i4,a)' ) 'method ',imethod,longrange
#endif
  ! =============================================================================
  !  start calculation of potential energy , field , efg , virial , forces ...
  ! =============================================================================
  u_coul_tot   = 0.0d0
  phi_coul_tot = 0.0d0
  vir_coul_tot = 0.0d0
  u_coul_qq   = 0.0d0
  phi_coul_qq = 0.0d0
  vir_coul_qq = 0.0d0
  u_coul_dd   = 0.0d0
  phi_coul_dd = 0.0d0
  vir_coul_dd = 0.0d0
  ef_t         = 0.0d0
  efg_t        = 0.0d0
  fx           = 0.0d0
  fy           = 0.0d0
  fz           = 0.0d0
 
  if ( .FALSE. ) then
 
  ! ===============================================
  !  get moment from polarisabilities
  ! ===============================================
  CALL moment_from_pola ( iastart, iaend , mu_ind ) 
 
  mu_tot = dipia + mu_ind

  ! =============================================
  !  coulombic energy , forces (field) and virial
  ! =============================================
  if ( longrange .eq. 'direct' ) CALL engforce_charge_DS ( iastart, iaend , ef_tmp )
  if ( longrange .eq. 'ewald'  ) CALL engforce_charge_ES ( iastart, iaend , ef_tmp )
  ef_t = ef_t + ef_tmp
  ef_tmp = 0.0d0
  ! =============================================
  !  Electric Field from (static) dipoles i.e linduced=.FALSE.
  ! =============================================
  if ( longrange .eq. 'direct' ) CALL engforce_dipole_DS ( iastart , iaend , ef_tmp , mu_tot )
  if ( longrange .eq. 'ewald' )  CALL engforce_dipole_ES ( iastart , iaend , ef_tmp , mu_tot )
  ef_t = ef_t + ef_tmp
  ef_tmp = 0.0d0

  if ( longrange .eq. 'direct' ) CALL engforce_charge_and_dipole_DS ( iastart , iaend , mu_tot )
!  if ( longrange .eq. 'ewald' ) CALL engforce_charge_and_dipole_ES ( iastart , iaend , mu_tot )

  u_coul_tot   = u_coul_qq   + u_coul_dd   + u_coul_qd + u_pol
  vir_coul_tot = vir_coul_qq + vir_coul_dd + vir_coul_qd
  phi_coul_tot = phi_coul_qq + phi_coul_dd

#ifdef debug 
  if ( ionode ) then
    WRITE ( stdout , '(a)' ) ' '
    WRITE ( stdout , '(a)' )     'Electric field at atoms : '
    do ia = 1 , natm
      WRITE ( stdout , '(i5,a3,a,3f18.10)' ) &
      ia,atype(ia),' Efield  = ', ef_t ( ia , 1)  , ef_t ( ia , 2 ) , ef_t ( ia , 3 )
    enddo
    WRITE ( stdout , '(a)' ) ' '

    WRITE ( stdout , '(a)' ) ' '
    WRITE ( stdout , '(a)' ) 'forces at atoms : '
    do ia = 1 , natm
      WRITE ( stdout , '(i5,a3,a,3f18.10)' ) &
      ia,atype(ia),' f       = ', fx ( ia )  , fy ( ia ) , fz ( ia )
    enddo
    WRITE ( stdout , '(a)' ) ' '
    WRITE ( stdout , '(a)' ) ' '
    WRITE ( stdout , '(a)' ) ' potential : '
    do ia = 1 , natm
      WRITE ( stdout , '(i5,a3,3(a,f18.10))' ) &
      ia,atype(ia),' phi_tot = ', phi_coul_tot ( ia )  , ' phi_qq = ', phi_coul_qq (ia) , ' phi_dd = ', phi_coul_dd (ia)
    enddo
    WRITE ( stdout , '(a)' ) ' '
    WRITE ( stdout , '(a)' ) 'Energy and virial : '
    WRITE ( stdout , '(5(a,f18.10))' ) &
    ' u_coul_tot   = ', u_coul_tot  ,' u_coul_qq   = ',u_coul_qq      ,' u_coul_qd   = ', u_coul_qd     ,' u_coul_dd   = ',u_coul_dd   ,' u_pol = ', u_pol
    WRITE ( stdout , '(4(a,f18.10))' ) &
    ' vir_coul_tot = ',vir_coul_tot ,' vir_coul_qq = ',vir_coul_qq    ,' vir_coul_qd = ', vir_coul_qd   ,' vir_coul_dd = ',vir_coul_dd
  endif
#endif

  ! =======
  !   EFG
  ! =======
  CALL get_efg_tensor ( iastart , iaend , efg_t , rm_coul , km_coul , cutlongrange , alphaES , mu_tot )  

  if ( longrange .eq. 'direct' ) CALL  Nymand_Linse_output ( 'DS' , '  ', u_coul_tot , phi_coul_tot , ef_t , efg_t , fx , vir_coul_tot , natm )
  if ( longrange .eq. 'ewald' )  CALL  Nymand_Linse_output ( 'ES' , ' 1', u_coul_tot , phi_coul_tot , ef_t , efg_t , fx , vir_coul_tot , natm ) 


  else ! tmp

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


  endif

  CLOSE (UNIT = kunit_POSFF )

  deallocate ( ef_tmp ) 
  deallocate ( efg_tmp ) 
  deallocate ( mu_tot ) 
  deallocate ( mu_ind ) 

  return

END SUBROUTINE Nymand_and_Linse_test


SUBROUTINE Nymand_Linse_output ( label , labelES , U , pot , Efield , EFG , fx , vir , natm )

  USE io_file,  ONLY : stdout

  implicit none

  ! global
  integer :: natm
  double precision :: U , vir
  double precision :: pot ( natm ) , fx ( natm )
  double precision :: Efield ( natm , 3 ) , EFG ( natm , 3 , 3 )
  character(len=2) label , labelES


  WRITE ( stdout , '(a)' ) ''
  WRITE ( stdout , 490 ) 'Method ',' surf ',' U ',' phi_1 ',' E1,x ', 'E1,xx ',' E1,yy ',' f1,x ',' PHI '
  WRITE ( stdout , 500 ) label , labelES, U , pot(1) , Efield ( 1 , 1 ) , EFG ( 1 , 1 , 1 ) , EFG ( 1 , 2 , 2 ) , fx ( 1 ) , vir
  WRITE ( stdout , '(a)' ) ''

490 FORMAT(A8,A8,7A12)
500 FORMAT(A8,A8,7F12.7)

  return

END SUBROUTINE Nymand_Linse_output

SUBROUTINE NL_field_print_info ( kunit ) 

  USE io_file,  ONLY :  ionode
  USE control,  ONLY :  calc , longrange , cutlongrange
  USE config,   ONLY :  box
  USE field,    ONLY :  ncelldirect , ncellewald , alphaES 

  implicit none

  ! global 
  integer :: kunit
  ! local
  double precision :: aaa

  if ( ionode ) then
      WRITE ( kunit ,'(a)')           '============================================================='
      WRITE ( kunit ,'(a)')           ''
      WRITE ( kunit ,'(a)')           'Nymand and Linse test simple systems : '
      WRITE ( kunit ,'(a)')           'atomic charge, dipole and polarisabilities'
      if ( longrange .eq. 'direct' )  then
        WRITE ( kunit ,'(a)')         'direct summation'
        WRITE ( kunit ,'(a)')         'cubic cutoff in real space'
        WRITE ( kunit ,'(a,f10.5)')   'distance cutoff           cutlongrange = ',cutlongrange
        WRITE ( kunit ,'(a,i10)')     '-ncelldirect ... ncelldirect     = ',ncelldirect
        WRITE ( kunit ,'(a,i10)')     'total number of cells            = ',( 2 * ncelldirect + 1 ) ** 3
      endif
      if ( longrange .eq. 'ewald' )  then
        CALL estimate_alpha( aaa )
        WRITE ( kunit ,'(a)')                   'ewald summation' 
        WRITE ( kunit ,'(a,f10.5,a,f10.5,a)')   'alpha = ',alphaES,' ( ',aaa,' ) '
        WRITE ( kunit ,'(a,f10.5)')             'cutoff in real part = ',cutlongrange
        WRITE ( kunit ,'(a,i10)')               'ncellewald = ',ncellewald
        WRITE ( kunit ,'(a,i10)')               ''
        WRITE ( kunit ,'(a,f10.5)')             'Note:this should hold alpha^2 * box^2 >> 1',alphaES*alphaES*box*box
        WRITE ( kunit ,'(a,f10.5)')             'from estimate_alpha ', aaa
      endif
      WRITE ( kunit ,'(a)')         'read config from file            : POSFF'
  endif
 
     


  return

END SUBROUTINE NL_field_print_info

SUBROUTINE moment_from_pola ( iastart , iaend , mu_ind ) 

  USE io_file, ONLY : ionode , stdout 
  USE config,  ONLY : natm , atype , fx , fy , fz , ntype , dipia , qia , ntypemax
  USE control, ONLY : longrange
  USE field,   ONLY : qch , lpolar , multipole_DS , induced_moment
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
  if ( .not. linduced ) then
    write ( stdout , '(a)' ) 'No polarizabilities'
    return
  endif

  if ( .FALSE. ) then

  ! =============================================
  !  calculate static Efield ( charge + dipoles )
  ! =============================================
  Efield_stat = 0.0d0

  ! =============================================
  !  coulombic energy , forces (field) and virial
  ! =============================================
  if ( longrange .eq. 'direct' ) CALL engforce_charge_DS ( iastart, iaend , ef_tmp )
  if ( longrange .eq. 'ewald' )  CALL engforce_charge_ES ( iastart, iaend , ef_tmp )
  Efield_stat = Efield_stat + ef_tmp
  ef_tmp = 0.0d0
  fx = 0.0d0 ; fy = 0.0d0 ; fz = 0.0d0
  ! =============================================
  !  Electric Field from (static) dipoles
  ! =============================================
  if ( longrange .eq. 'direct' ) CALL engforce_dipole_DS ( iastart , iaend , ef_tmp , dipia )
  if ( longrange .eq. 'ewald' )  CALL engforce_dipole_ES ( iastart , iaend , ef_tmp , dipia  )
  Efield_stat = Efield_stat + ef_tmp
  fx = 0.0d0 ; fy = 0.0d0 ; fz = 0.0d0
  ef_tmp = 0.0d0

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
  conv_tol = 1e-6 !hardware
  do while ( iscf .lt. 3 )
!  do while ( diff_efield .gt. conv_tol .or. iscf .lt. 2 )

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
    !  Efield_ind out , mu_ind in
    ! ==========================================================
    if ( longrange .eq. 'direct' ) CALL engforce_dipole_DS ( iastart , iaend , Efield_ind , mu_ind )
    if ( longrange .eq. 'ewald' )  CALL engforce_dipole_ES ( iastart , iaend , Efield_ind , mu_ind )
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

  else ! tmp

  ! =============================================
  !  calculate static Efield ( charge + dipoles )
  ! =============================================
  Efield_stat = 0.0d0

  ! =============================================
  !  coulombic energy , forces (field) and virial
  ! =============================================
  if ( longrange .eq. 'direct' ) CALL  multipole_DS ( iastart, iaend , ef_tmp , efg_tmp , dipia , u_tmp , vir_tmp , phi_tmp ) 

  Efield_stat = Efield_stat + ef_tmp

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
  conv_tol = 1e-6 !hardware
  do while ( iscf .lt. 3 )

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





  endif



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



