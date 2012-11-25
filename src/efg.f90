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

#define debug

!fix_grid in efg
!#define fix_grid

! ======= Hardware =======

!*********************** MODULE efg  **********************************
!**********************************************************************
MODULE efg 
 
  USE kspace
  USE rspace

  implicit none

  logical :: lefgprintall               ! print ( or not ) all the efg for each atoms and configs to file EFGALL
  logical :: lefg_it_contrib            ! if ones want to get the contributions to EFG separated in types
  logical :: lefg_restart               ! if EFGALL files are ready
  logical :: lcharge_only               ! use only point charge for EFG calculation
  logical :: ldipole_only               ! use only dipoles for EFG calculation
  logical :: lcharge_and_dipole         ! use point charges and dipoles for EFG calculation
  integer :: ncefg                      ! number of configurations READ  for EFG calc (only when calc = 'efg')
  integer :: ntcor                      ! maximum number of steps for the acf calculation (calc = 'efg+acf')

#ifdef fix_grid
  double precision, dimension ( : , : ) , allocatable :: rgrid
#endif fix_grid

  ! ===================
  !  direct summation
  ! ===================
  integer                                           :: ncelldirect
  double precision                                  :: cutefg       ! cut-off distance for efg calc (direct dimension) 
  double precision, dimension(:,:,:)  , allocatable :: efg_t       ! efg_tensor
  double precision, dimension(:,:,:,:), allocatable :: efg_t_it       ! efg_tensor
  double precision, dimension(:,:,:)  , allocatable :: efg_ia       ! efg_tensor
  double precision, dimension(:,:,:,:), allocatable :: efg_ia_it       ! efg_tensor
  TYPE ( rmesh ) :: rm_efg

  ! ==================
  !  ewald summation
  ! ==================
  integer                                           :: ncellewald
  double precision                                  :: alphaES          ! Ewald sum parameter 

  TYPE ( kmesh ) :: km_efg
  TYPE ( kmesh ) :: km_efg_dip

  ! ==============
  !  distribution
  ! ==============
  integer         , dimension(:,:)  , allocatable :: dibvzztot ! vzz distrib. dim. = (ntype , PAN)
  integer         , dimension(:,:)  , allocatable :: dibetatot ! eta distrib. dim. = (ntype , PAN) 
  integer         , dimension(:,:,:), allocatable :: dibUtot   ! U_i distrib. dim. = ( 1:6 , ntype , PAN ) 

  double precision                                :: resvzz    ! resolution in vzz distribution
  double precision                                :: reseta    ! resolution in eta distribution
  double precision                                :: resu      ! resolution in Ui ditribution
  double precision                                :: vzzmin    ! minimum value for vzz (distrib. in [vzzmin, -vzzmin]
  double precision                                :: umin      ! minimum value of umin
  integer                                         :: PANeta    ! nb of bins in eta distrib. (related reseta)
  integer                                         :: PANvzz    ! nb of bins in vzz distrib. (related resvzz)
  integer                                         :: PANU      ! nb of bins in Ui  distrib. (related resu)

CONTAINS

!*********************** SUBROUTINE efg_init **********************************
!
! iniitialisation
!
!******************************************************************************

SUBROUTINE efg_init

  USE io_file 
  USE control,  ONLY :  calc
  USE prop,     ONLY :  lefg
  USE field,    ONLY :  initialize_coulomb
 
  implicit none
  
  ! local
  integer :: ioerr
  character * 132 :: filename

  namelist /efgtag/  lefgprintall     , & 
                     lefg_it_contrib  , &
                     lefg_restart     , &
                     lcharge_only      , &
                     lcharge_and_dipole, &
                     ldipole_only     , &
                     ncefg            , & 
                     ntcor            , &
                     ncellewald       , &
                     ncelldirect      , &
                     cutefg           , &
                     alphaES          , &
                     resvzz           , &
                     reseta           , &
                     resu             , &
                     vzzmin           , &
                     umin           

  if ( .not. lefg .and. calc .ne. 'efg' .and. calc .ne.'efg+acf' ) return 

  ! ===================
  !  set default values
  ! ===================
  CALL efg_default_tag
  ! ===================
  !  reads efgtag tags
  ! ===================
  CALL getarg (1,filename)
  OPEN ( stdin , file = filename)
  READ ( stdin , efgtag,iostat=ioerr)
  if ( ioerr .lt. 0 )  then
   if ( ionode ) WRITE ( stdout, '(a)') 'ERROR reading input_file : efgtag section is absent'
   STOP
  elseif ( ioerr .gt. 0 )  then
   if ( ionode ) WRITE ( stdout, '(a)') 'ERROR reading input_file : efgtag wrong tag'
   STOP
  endif

  CLOSE ( stdin )
  ! ===============
  !  check efgtag
  ! ===============
  CALL efg_check_tag
   if ( calc .eq. 'efg') return
  ! =============
  !  print info
  ! =============
  CALL efg_print_info(stdout)
  CALL efg_print_info(kunit_OUTFF)

  return

END SUBROUTINE efg_init

!*********************** SUBROUTINE efg_default_tag ***************************
!
! set default values for efgtag
!
!******************************************************************************

SUBROUTINE efg_default_tag

  implicit none
  
  ! =================
  !  default values
  ! =================
  lefgprintall       = .false. 
  lefg_it_contrib    = .false.
  lefg_restart       = .false.
  lcharge_only       = .true.
  ldipole_only       = .false.
  lcharge_and_dipole = .false.
  reseta             =   0.1D0
  resvzz             =   0.1D0
  resu               =   0.1D0 
  ncelldirect        =   0
  ncellewald         =   0 
  ncefg              =   0
  cutefg             =  30.0D0
  alphaES            =   1.0D0
  umin               =  -4.0d0
  vzzmin             =  -4.0d0
  ntcor              =  10

  return

END SUBROUTINE efg_default_tag


!*********************** SUBROUTINE efg_check_tag ***************************
!
! check efg parameters ( check longrange and ncell ) 
!
!******************************************************************************

SUBROUTINE efg_check_tag

  USE field,    ONLY :  qch 
  USE control,  ONLY :  calc , longrange
  USE io_file,  ONLY :  ionode , stdout
  USE prop,     ONLY :  lefg
  USE config,   ONLY :  ntype

  implicit none

  if ( calc .eq. 'efg+acf' .and. ntype .gt. 2 ) then
    if ( ionode ) WRITE ( stdout , '(a)' ) 'ERROR the subroutine efg_acf is not implemented for ntype > 2'
    STOP 
  endif

  if ( calc .eq. 'efg+acf' ) return

  if ( qch(1) .eq. 0 .and. qch(2) .ne. 0 ) then 
    if ( ionode ) WRITE ( stdout ,'(a,2f8.3)') &
    'ERROR: who requested EFG calculation without charge definition, something must be wrong' 
    STOP
  endif

  ! ===================
  !  direct summation
  ! ===================
  if (longrange .eq. 'direct' ) then
    if ( ncelldirect .eq. 0 ) then
      if ( ionode ) WRITE ( stdout ,'(a,2f8.3)') 'ERROR eftag: ncelldirect null'
      STOP 
    endif
  endif

  ! ===================
  !  ewald summation
  ! ===================
  if (longrange .eq. 'ewald' ) then
    if ( ncellewald .eq. 0 ) then
      if ( ionode ) WRITE ( stdout ,'(a,2f8.3)') 'ERROR eftag: ncellewald null'
      STOP
   endif
  endif

  ! ==================================
  !  check vzzmin and Umin .ne. 0
  ! ==================================
  if ( vzzmin .eq. 0d0 .or. umin .eq. 0d0 ) then
    if ( ionode ) WRITE ( stdout ,'(a,2f8.3)') 'ERROR efgtag: vzzmin or umin should be set',vzzmin,umin   
    STOP
  endif             
  ! ==========================================
  !  set PAN (nb bins) from resolution values
  ! ==========================================
  PANeta = int(1.0d0/reseta)
  PANvzz = int((2.0d0*dabs(vzzmin))/resvzz)
  PANU = int((2.0D0*dabs(umin))/resu)

  return

END SUBROUTINE efg_check_tag

!*********************** SUBROUTINE efg_print_info ***************************
!
! print information to standard output about efg calculation 
!
!******************************************************************************

SUBROUTINE efg_print_info(kunit)

  USE io_file,  ONLY :  ionode 
  USE control,  ONLY :  calc , longrange 
  USE config,   ONLY :  box      
  USE prop,     ONLY :  lefg

  implicit none

  !local
  integer :: kunit
  integer :: nkcut_full
  integer :: nkcut_half
  double precision :: aaa

  nkcut_full = ( 2 * ncellewald + 1 ) ** 3
  nkcut_half = ( ncellewald + 1 ) ** 3

  if ( ionode ) then
    if ( calc .eq. 'efg' .or. lefg ) then
      WRITE ( kunit ,'(a)')           '=============================================================' 
      WRITE ( kunit ,'(a)')           ''
      WRITE ( kunit ,'(a)')           'electric field gradient:'
      WRITE ( kunit ,'(a)')           'point charges calculation'
      if ( longrange .eq. 'direct' )  then
        WRITE ( kunit ,'(a)')         'direct summation'
        WRITE ( kunit ,'(a)')         'cubic cutoff in real space'
        WRITE ( kunit ,'(a,f10.5)')   'distance cutoff           cutefg = ',cutefg
        WRITE ( kunit ,'(a,i10)')     '-ncelldirect ... ncelldirect     = ',ncelldirect
        WRITE ( kunit ,'(a,i10)')     'total number of cells            = ',( 2 * ncelldirect + 1 ) ** 3
      endif     
      if ( longrange .eq. 'ewald' )  then
        CALL estimate_alpha( aaa )
        WRITE ( kunit ,'(a)')                   'ewald summation'
        WRITE ( kunit ,'(a,f10.5,a,f10.5,a)')   'alpha                            = ',alphaES,' ( ',aaa,' ) '
        WRITE ( kunit ,'(a,f10.5)')             'cutoff in real part              = ',cutefg
        WRITE ( kunit ,'(a,i10)')               'ncellewald                       = ',ncellewald
        WRITE ( kunit ,'(a,i10)')               'nckut (full)                     = ',nkcut_full
        WRITE ( kunit ,'(a,i10)')               'nckut (half)                     = ',nkcut_half
        WRITE ( kunit ,'(a,i10)')               '' 
        WRITE ( kunit ,'(a,f10.5)')             'Note:this should hold alpha^2 * box^2 >> 1',alphaES*alphaES*box*box
        WRITE ( kunit ,'(a,f10.5)')             'from estimate_alpha ', aaa
      endif
      if ( calc .eq. 'efg')  then 
        WRITE ( kunit ,'(a)')         'read config from file            : TRAJFF'
        WRITE ( kunit ,'(a,i10)')     'numbers of config read ncefg     = ',ncefg 
      endif
      WRITE ( kunit ,'(a)')           'write averages values in file    : EFGFF '          

      WRITE ( kunit ,'(a)')           'Distributions:'
      WRITE ( kunit ,'(a)')           'eta distrib                : DTETAFF'
      WRITE ( kunit ,'(a)')           'Vzz distrib                : DTVZZFF'
      WRITE ( kunit ,'(a)')           'Ui components distrib      : DTIBUFF'           
      if ( lefgprintall )  then
        WRITE ( kunit ,'(a)')         'EFG for all atoms          : EFGALL '
      endif
      WRITE ( kunit ,'(a)')           ''
      WRITE ( kunit ,'(a)')           'distributions parameters:'
      WRITE ( kunit ,'(a,f10.5)')     'eta between    0.00000 and    1.00000 with resolution ',reseta
      WRITE ( kunit ,'(a,f10.5,a,f10.5,a,f10.5)') &
                                      'vzz between ',vzzmin,' and ',-vzzmin,' with resolution ',resvzz
      WRITE ( kunit ,'(a,f10.5,a,f10.5,a,f10.5)') &
                                      'Ui  between ',  umin,' and ',  -umin,' with resolution ',resu
    endif
    if ( calc .eq. 'efg+acf' ) then
      WRITE ( kunit ,'(a)')           '============================================================='
      WRITE ( kunit ,'(a)')           ''
      WRITE ( kunit ,'(a)')           'electric field gradient auto-correlation function:'                    
      WRITE ( kunit ,'(a)')           'data from file EFGALL'
      WRITE ( kunit ,'(a,i10)')       'numbers of config read ncefg           = ',ncefg 
      WRITE ( kunit ,'(a,i10)')       'maximum of evaluated correlation step  = ',ntcor
      WRITE ( kunit ,'(a)')           'output file                            : EFGACFFF'
    endif
  endif

  return

END SUBROUTINE efg_print_info

!*********************** SUBROUTINE efg_alloc *********************************
!
! allocate quantities related to efg calculation
! /deallocate 
!
!******************************************************************************

SUBROUTINE efg_alloc

  USE control,  ONLY :  calc , longrange
  USE config,   ONLY :  natm, ntype , qia , itype
  USE prop,     ONLY :  lefg
  USE field,    ONLY :  qch

  implicit none

  ! local
  integer :: nkcut
  integer :: ncmax
  integer :: ia , it

  ! ======
  !  lefg
  ! ======
  if ( .not. lefg .and. ( calc .ne. 'efg' ) ) return 

  if ( ntype .eq.  1 ) then
    do ia = 1 , natm
      qia(ia) = qch(1)
    enddo
  endif

  if ( ntype .ge.  2 ) then
    do ia = 1 , natm
      do it = 1 ,ntype
        if (itype(ia) .eq. it) then
          qia(ia) = qch(it)
        endif
      enddo
    enddo
  endif

  allocate( dibUtot(6,0:ntype,0:PANU) )
  allocate( dibvzztot(0:ntype,0:PANvzz) )
  allocate( dibetatot(0:ntype,0:PANeta) )  
  dibUtot = 0
  dibvzztot = 0
  dibetatot = 0

#ifdef fix_grid
  allocate ( rgrid ( 3 , natm ) )
#endif fix_grid
  allocate( efg_t    ( natm , 3 , 3 ) )
  allocate( efg_t_it ( natm , ntype , 3 , 3 ) )
  allocate( efg_ia    ( natm , 3 , 3 ) )
  allocate( efg_ia_it ( natm , ntype , 3 , 3 ) )
  efg_t = 0.0d0
  efg_t_it = 0.0d0
  efg_ia = 0.0d0
  efg_ia_it = 0.0d0

  ! ============
  !  direct sum
  ! ============
  if ( longrange .eq. 'direct' ) then
    ncmax = ( 2 * ncelldirect + 1 ) ** 3
    rm_efg%ncmax=ncmax
    rm_efg%ncell=ncelldirect
    allocate( rm_efg%boxxyz( 3 , ncmax ) , rm_efg%lcell( ncmax ) )
    CALL direct_sum_init ( rm_efg )
  endif

  ! ============
  !  ewald sum
  ! ============
  if ( longrange .eq. 'ewald')  then

    ! ==================================
    !  Time -Reversal Symmetry is conserved 
    !     k-point mesh -k and +k 
    ! ==================================
    km_efg%meshlabel='km_efg'
    km_efg%ncell     = ncellewald
    nkcut = ( 2 * ncellewald + 1 ) ** 3 
    nkcut = nkcut - 1
    km_efg%nkcut = nkcut
    allocate( km_efg%kptk( nkcut ) , km_efg%kpt(3,nkcut) )
    allocate ( km_efg%strf ( nkcut, ntype) ) 
    CALL kpoint_sum_init ( km_efg )
    ! ==================================
    !  Time -Reversal Symmetry is broken
    !     k-point mesh only +k 
    ! ==================================
    km_efg_dip%ncell = ncellewald
    km_efg_dip%meshlabel='km_efg_dip'
    !nkcut = ( ncellewald + 1 ) ** 3
    nkcut = ( ncellewald + 1 ) ** 3
    nkcut = nkcut - 1
    km_efg_dip%nkcut = nkcut
    allocate( km_efg_dip%kptk( nkcut ) , km_efg_dip%kpt(3,nkcut) )
    allocate ( km_efg_dip%strf ( nkcut, ntype) ) 
    CALL kpoint_sum_init_half ( km_efg_dip )

  endif

  return 

END SUBROUTINE efg_alloc

SUBROUTINE efg_dealloc

  USE control,  ONLY :  calc , longrange 
  USE prop,     ONLY :  lefg

  implicit none

  if ( .not. lefg .and. ( calc .ne. 'efg' ) ) return

  deallocate( dibUtot )
  deallocate( dibvzztot )
  deallocate( dibetatot )
#ifdef fix_grid
  deallocate ( rgrid )
#endif fix_grid
    deallocate( efg_t )  
    deallocate( efg_t_it )  

  if ( longrange .eq. 'direct' ) then 
    deallocate( rm_efg%boxxyz , rm_efg%lcell )
  endif

  if ( longrange .eq. 'ewald' )  then

    deallocate( km_efg%kptk , km_efg%kpt )
    deallocate( km_efg%strf ) 

    deallocate( km_efg_dip%kptk , km_efg_dip%kpt )
    deallocate( km_efg_dip%strf ) 

  endif

  return

END SUBROUTINE efg_dealloc

!*********************** SUBROUTINE efgcalc ***********************************
!
! this SUBROUTINE initialize the calculation of efg when calc = 'efg' 
! from file TRAJFF. This is different from lefg which calculates EFG on-the-fly 
! during MD trajectories. So calc = 'md' + lefg + lgr
! maybe this is a too complex organisation ...
!
!******************************************************************************

SUBROUTINE efgcalc 

  USE io_file
  USE config,   ONLY :  system , natm , ntype , atype , rx , ry , rz , itype , & 
                        atypei , natmi, rho , box , omega , config_alloc , qia, dipia , ipolar
  USE control,  ONLY :  longrange , myrank , numprocs
  USE field,    ONLY :  qch , dip , dip_ind , field_init , Efield_and_scf_induced , initialize_coulomb , lpolar

  implicit none
  INCLUDE 'mpif.h'

  ! local
  integer :: iastart , iaend
  integer :: ia, iconf , na , it , it2 , ik , cc , ccs
  double precision :: aaaa 
  integer :: iiii , ierr
  character * 60 :: cccc 
  logical :: linduced
  integer :: kunit
  integer :: kunit_it(2) ! <--  not general enough should be ntype dependent

#ifdef fix_grid
  double precision, dimension ( : , : ) , allocatable :: rave !average positions
#endif

  ! ==================================
  !  if lefg_restart EFGALL are ready
  ! ==================================
  if ( .not. lefg_restart ) then

    OPEN (unit = kunit_EFGALL  ,file = 'EFGALL', STATUS='REPLACE')

    if ( lefg_it_contrib ) then
      OPEN (unit = kunit_EFGALLIT1  ,file = 'EFGALLIT1', STATUS='REPLACE')
      OPEN (unit = kunit_EFGALLIT2  ,file = 'EFGALLIT2', STATUS='REPLACE')
    endif

    OPEN (UNIT = kunit_TRAJFF , FILE = 'TRAJFF')

    READ ( kunit_TRAJFF, * )  natm 
    READ ( kunit_TRAJFF, * )  system
    READ ( kunit_TRAJFF, * )  box,ntype
    READ ( kunit_TRAJFF ,* ) ( atypei ( it ) , it = 1 , ntype )
    IF ( ionode ) WRITE ( kunit_OUTFF ,'(A,20A3)' ) 'found type information on TRAJFF : ', atypei ( 1:ntype )
    IF ( ionode ) WRITE ( stdout      ,'(A,20A3)' ) 'found type information on TRAJFF : ', atypei ( 1:ntype )
    READ( kunit_TRAJFF ,*)   ( natmi ( it ) , it = 1 , ntype )
    omega = box * box * box
    rho = natm / omega
    ! ===================================
    !  here we know natm, then alloc 
    !  and decomposition can be applied 
    ! ================================== 
    CALL config_alloc 
    CALL do_split ( natm , myrank , numprocs , iastart , iaend )
    ! read charge in fieldtag
    CALL field_init
    CALL efg_alloc
    ! =============
    !  print info
    ! =============
    CALL efg_print_info(stdout)
    CALL efg_print_info(kunit_OUTFF)
  
    CALL print_general_info( stdout )
    CALL print_general_info( kunit_OUTFF )
  
  ! =============================================
  !  not fix_grid: we calculate efg at 
  !  each atom positions which are moving.
  ! =============================================
  
    ! ============================================
    ! LOOP OVER CONFIGURATIONS 
    ! ============================================
    do iconf = 1, ncefg
      na = 0
      if ( iconf .ne. 1 ) READ ( kunit_TRAJFF, * )  iiii 
      if ( iconf .ne. 1 ) READ ( kunit_TRAJFF, * )  cccc
      if ( iconf .ne. 1 ) READ ( kunit_TRAJFF, * )  aaaa,iiii 
      if ( iconf .ne. 1 ) READ ( kunit_TRAJFF , * ) ( cccc , it = 1 , ntype )
      if ( iconf .ne. 1 ) READ ( kunit_TRAJFF , * ) ( iiii , it = 1 , ntype )
      do ia = 1 , natm
        READ ( kunit_TRAJFF, * ) atype ( ia ) , rx ( ia ) , ry ( ia ) , rz ( ia )
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

      efg_t = 0.0d0
      if ( lcharge_only .or. lcharge_and_dipole ) then

        CALL initialize_coulomb
        CALL Efield_and_scf_induced ( iastart , iaend , dip_ind , lfield_only=.TRUE. )
        dip_ind = 0.0d0
         
        ! ===============
        ! use direct sum
        ! ===============
        if ( longrange .eq. 'direct' )  CALL efg_DS ( iastart , iaend )
        ! ===============
        ! use ewald sum
        ! ===============
        if ( longrange .eq. 'ewald' )   CALL efg_ES ( iastart , iaend )

        if ( ionode ) WRITE( stdout , '(a,i6,a)' ) 'efg calculation of config',iconf,' from point charges '

        efg_t  = efg_t + efg_ia 
        efg_t_it = efg_t_it + efg_ia_it 
        CALL print_tensor( efg_ia ( 1 , : , : ) , 'TOTCHG' )

      endif

      if ( ldipole_only .or. lcharge_and_dipole ) then 
          if ( .not. lcharge_and_dipole )  CALL initialize_coulomb
          CALL Efield_and_scf_induced ( iastart , iaend , dip_ind , lfield_only=.FALSE. )

         ! ===========================
         !  mu_tot = mu_stat + mu_ind
         !  we update dipia ( atoms ) 
         !  and       dip   ( types ) 
         ! ===========================
         dipia = dipia + dip_ind
         ik = 0
         do it = 1 , ntype
           ik = ik + natmi(it)
           dip ( it , 1 ) = dip ( it , 1 ) + dip_ind ( ik , 1 )
           dip ( it , 2 ) = dip ( it , 2 ) + dip_ind ( ik , 2 )
           dip ( it , 3 ) = dip ( it , 3 ) + dip_ind ( ik , 3 )
         enddo


        ! ===============
        ! use direct sum
        ! ===============
        if ( longrange .eq. 'direct' )  CALL efg_DS_dip ( iastart , iaend )

        ! ===============
        ! use ewald sum
        ! ===============
        if ( longrange .eq. 'ewald' )   CALL efg_ES_dip ( iastart , iaend )
 
        if ( ionode ) WRITE( stdout , '(a,i6,a)' ) 'efg calculation of config',iconf,' from dipoles '

        efg_t = efg_t + efg_ia 
        efg_t_it = efg_t_it + efg_ia_it 

        CALL print_tensor( efg_ia ( 1 , : , : ) , 'TOTDIP' )

      endif

CALL print_tensor( efg_t( 1 , : , : ) , 'TOTEFG' )

  ! =======================================
  ! write efg for each atom in file EFGALL 
  ! =======================================
  if ( ionode  .and. lefgprintall ) then
    WRITE ( kunit_EFGALL , * )  natm
    WRITE ( kunit_EFGALL , * )  system
    WRITE ( kunit_EFGALL , * )  box,ntype
    WRITE ( kunit_EFGALL , * )  ( atypei ( it ) , it = 1 , ntype )
    WRITE ( kunit_EFGALL , * )  ( natmi  ( it ) , it = 1 , ntype )
    WRITE ( kunit_EFGALL ,'(a)') &
    '      ia type                   vxx                   vyy                   vzz                   vxy                   vxz                   vyz'
    do ia = 1 , natm 
      WRITE ( kunit_EFGALL ,'(i8,2x,a3,6f22.16)') &
      ia , atype ( ia ) , efg_t ( ia , 1 , 1) , efg_t ( ia , 2 , 2) , &
                          efg_t ( ia , 3 , 3) , efg_t ( ia , 1 , 2) , &
                          efg_t ( ia , 1 , 3) , efg_t ( ia , 2 , 3)
  
    enddo
  endif

  if ( lefg_it_contrib ) then
  kunit_it(1) = kunit_EFGALLIT1
  kunit_it(2) = kunit_EFGALLIT2
    do it = 1 , 2 ! <- should be ntype dependent
      kunit = kunit_it(it)
      if ( ionode  .and. lefgprintall ) then
        WRITE ( kunit , * )  natm
        WRITE ( kunit , * )  system
        WRITE ( kunit , * )  box,ntype
        WRITE ( kunit , * )  ( atypei ( it2 ) , it2 = 1 , ntype )
        WRITE ( kunit , * )  ( natmi  ( it2 ) , it2 = 1 , ntype )
        WRITE ( kunit ,'(a)') &
       '      ia type                   vxx                   vyy                   vzz                   vxy                   vxz                   vyz'
        do ia = 1 , natm
          WRITE ( kunit ,'(i8,2x,a3,6f22.16)') &
          ia , atype ( ia ) , efg_t_it ( ia , it , 1 , 1) , efg_t_it ( ia , it , 2 , 2) , &
                              efg_t_it ( ia , it , 3 , 3) , efg_t_it ( ia , it , 1 , 2) , &
                              efg_t_it ( ia , it , 1 , 3) , efg_t_it ( ia , it , 2 , 3)
        enddo
      endif
    enddo
  endif

    enddo ! iconf loop

    CLOSE(kunit_EFGALL)

    if ( lefg_it_contrib ) then
      CLOSE( kunit_EFGALLIT1 )
      CLOSE( kunit_EFGALLIT2 )
    endif

  endif !lefg_restart

  !calculate statistical properties from EFGALL
  CALL io_open  ( kunit_EFGALL , 'EFGALL' , 'unknown' )
  CALL io_open  ( kunit_EFGFF  , 'EFGFF'  , 'unknown' )
  CALL io_open  ( kunit_NMRFF  , 'NMRFF'  , 'unknown' )
  CALL efg_stat( kunit_EFGALL , kunit_EFGFF , kunit_NMRFF )
  CALL io_close ( kunit_NMRFF ) 
  CALL io_close ( kunit_EFGFF ) 
  CALL io_close ( kunit_EFGALL ) 
  ! write average distribution output from EFGALL 
  CALL io_open  ( kunit_DTETAFF , 'DTETAFF' , 'unknown' ) 
  CALL io_open  ( kunit_DTVZZFF , 'DTVZZFF' , 'unknown' ) 
  CALL io_open  ( kunit_DTIBUFF , 'DTIBUFF' , 'unknown' ) 
  CALL efg_write_output( kunit_DTETAFF , kunit_DTVZZFF , kunit_DTIBUFF )
  CALL io_close ( kunit_DTETAFF )
  CALL io_close ( kunit_DTVZZFF )
  CALL io_close ( kunit_DTIBUFF )
  dibUtot = 0
  dibvzztot = 0
  dibetatot = 0

  CALL MPI_BARRIER( MPI_COMM_WORLD , ierr )

  if ( lefg_it_contrib ) then

    ! calculate statistical properties from EFGALLIT1
    CALL io_open  ( kunit_EFGALLIT1 , 'EFGALLIT1' , 'unknown' )
    CALL io_open  ( kunit_EFGFFIT1  , 'EFGFFIT1'  , 'unknown' )
    CALL io_open  ( kunit_NMRFFIT1  , 'NMRFFIT1'  , 'unknown' )
    CALL efg_stat ( kunit_EFGALLIT1 , kunit_EFGFFIT1 , kunit_NMRFFIT1 )
    CALL io_close ( kunit_NMRFFIT1 ) 
    CALL io_close ( kunit_EFGFFIT1 ) 
    CALL io_close ( kunit_EFGALLIT1 ) 
    ! write average distribution output from EFGALLIT1
    CALL io_open  ( kunit_DTETAFFIT , 'DTETAFFIT' , 'unknown' ) 
    CALL io_open  ( kunit_DTVZZFFIT , 'DTVZZFFIT' , 'unknown' ) 
    CALL io_open  ( kunit_DTIBUFFIT , 'DTIBUFFIT' , 'unknown' ) 
    CALL efg_write_output ( kunit_DTETAFFIT , kunit_DTVZZFFIT , kunit_DTIBUFFIT )
    dibUtot = 0
    dibvzztot = 0
    dibetatot = 0

    CALL MPI_BARRIER( MPI_COMM_WORLD , ierr )

    !calculate statistical properties from EFGALLIT2
    CALL io_open  ( kunit_EFGALLIT2 , 'EFGALLIT2' , 'unknown' )
    CALL io_open  ( kunit_EFGFFIT2  , 'EFGFFIT2'  , 'unknown' )
    CALL io_open  ( kunit_NMRFFIT2  , 'NMRFFIT2'  , 'unknown' )
    CALL efg_stat ( kunit_EFGALLIT2 , kunit_EFGFFIT2 , kunit_NMRFFIT2 )
    CALL io_close ( kunit_NMRFFIT2 ) 
    CALL io_close ( kunit_EFGFFIT2 ) 
    CALL io_close ( kunit_EFGALLIT2 ) 
    ! write average distribution output from EFGALLIT1
    CALL efg_write_output ( kunit_DTETAFFIT , kunit_DTVZZFFIT , kunit_DTIBUFFIT )
    CALL io_close ( kunit_DTETAFFIT ) 
    CALL io_close ( kunit_DTVZZFFIT ) 
    CALL io_close ( kunit_DTIBUFFIT ) 

  endif

  CLOSE(kunit_TRAJFF)

  return

END SUBROUTINE efgcalc

!*********************** SUBROUTINE efg_DS ************************************
!
! Direct summation to calculate electric-field-gradient.
! 
! parallelisation atom decomposition
! 
! I can separate the actual caculation of efg and the distribution calculation
!
!******************************************************************************

SUBROUTINE efg_DS ( iastart , iaend )

  USE control,  ONLY :  myrank, numprocs, calc 
  USE config,   ONLY :  system , natm , natmi , atype , atypei , box , itype , rx , ry , rz , ntype , qia 
  USE field,    ONLY :  qch
  USE prop,     ONLY :  nprop_print
  USE time
  USE io_file

  implicit none

  INCLUDE 'mpif.h'

  ! global
  integer, intent(in)           :: iastart , iaend 

  ! local
  integer :: ia, ja, ierr , it 
  integer :: ncell
  double precision :: d , d2 , d5 , dm5
  double precision :: rxi , ryi , rzi , rxj , ryj , rzj , rxij , ryij , rzij
  double precision :: cutefgsq
  ! time
  double precision :: ttt1 , ttt2 

#ifdef debug
  CALL print_config_sample(0,0)
#endif

  ttt1 = MPI_WTIME(ierr)

  cutefgsq = cutefg * cutefg
  efg_ia(:,:,:) = 0.0d0
  efg_ia_it(:,:,:,:) = 0.0d0

! =========================================================
!  MAIN LOOP calculate EFG(i) for each atom i parallelized
! =========================================================

#ifdef debug
     write(*,*) 'iastart iaend',iastart , iaend
     write(*,*) 'rm_efg%ncmax',rm_efg%ncmax
#endif
atom : do ia = iastart , iaend
    rxi = rx(ia)
    ryi = ry(ia)
    rzi = rz(ia)
    ! ==================================================
    ! sum over neighboring cells (see direct_sum_init )
    ! ==================================================
    do ncell = 1 , rm_efg%ncmax
         ! ==============================
         !  ia and ja in different cells
         ! ==============================
         if ( rm_efg%lcell(ncell) .eq. 1) then

          do ja = 1 , natm
            rxj  = rx(ja) + rm_efg%boxxyz(1,ncell)
            ryj  = ry(ja) + rm_efg%boxxyz(2,ncell)
            rzj  = rz(ja) + rm_efg%boxxyz(3,ncell)
            rxij = rxi - rxj
            ryij = ryi - ryj
            rzij = rzi - rzj
            d2   = rxij * rxij + ryij * ryij + rzij * rzij

            if ( d2 .lt. cutefgsq ) then
              d = dsqrt(d2)
              d5 = d2 * d2 * d
              dm5 = 1.0d0/d5
              dm5 = dm5 * qia(ja)

              efg_ia( ia , 1 , 1 )  = efg_ia( ia , 1 , 1 ) - (3.0d0 * rxij * rxij - d2 ) * dm5
              efg_ia( ia , 2 , 2 )  = efg_ia( ia , 2 , 2 ) - (3.0d0 * ryij * ryij - d2 ) * dm5
              efg_ia( ia , 3 , 3 )  = efg_ia( ia , 3 , 3 ) - (3.0d0 * rzij * rzij - d2 ) * dm5
              efg_ia( ia , 1 , 2 )  = efg_ia( ia , 1 , 2 ) -  3.0d0 * rxij * ryij * dm5
              efg_ia( ia , 1 , 3 )  = efg_ia( ia , 1 , 3 ) -  3.0d0 * rxij * rzij * dm5
              efg_ia( ia , 2 , 3 )  = efg_ia( ia , 2 , 3 ) -  3.0d0 * ryij * rzij * dm5

              if ( lefg_it_contrib ) then
                it = itype(ja)  
                efg_ia_it( ia , it , 1 , 1 )  = efg_ia_it( ia , it , 1 , 1 ) - &
                ( 3.0d0 * rxij * rxij - d2 ) * dm5
                efg_ia_it( ia , it , 2 , 2 )  = efg_ia_it( ia , it , 2 , 2 ) - & 
                ( 3.0d0 * ryij * ryij - d2 ) * dm5
                efg_ia_it( ia , it , 3 , 3 )  = efg_ia_it( ia , it , 3 , 3 ) - & 
                ( 3.0d0 * rzij * rzij - d2 ) * dm5
                efg_ia_it( ia , it , 1 , 2 )  = efg_ia_it( ia , it , 1 , 2 ) -  3.0d0 * rxij * ryij * dm5
                efg_ia_it( ia , it , 1 , 3 )  = efg_ia_it( ia , it , 1 , 3 ) -  3.0d0 * rxij * rzij * dm5
                efg_ia_it( ia , it , 2 , 3 )  = efg_ia_it( ia , it , 2 , 3 ) -  3.0d0 * ryij * rzij * dm5
              endif   

            endif ! d2.lt.cutefgsq
 
          enddo ! ja
        endif 
         ! =======================================
         !  ia and ja in the same cell (ia.ne.ja)
         ! =======================================
        if ( rm_efg%lcell(ncell) .eq. 0) then

          do ja = 1,natm

            if (ja.ne.ia) then
              rxj  = rx(ja)
              ryj  = ry(ja)
              rzj  = rz(ja)
              rxij = rxi - rxj
              ryij = ryi - ryj
              rzij = rzi - rzj
              d2   = rxij * rxij + ryij * ryij + rzij * rzij
 
              if (d2.lt.cutefgsq) then
                d = dsqrt(d2)
                d5 = d2 * d2 * d
                dm5 = 1.0d0/d5
                dm5 = dm5 * qia(ja)
                efg_ia( ia , 1 , 1 )  = efg_ia( ia , 1 , 1 ) - (3.0d0 * rxij * rxij - d2 ) * dm5
                efg_ia( ia , 2 , 2 )  = efg_ia( ia , 2 , 2 ) - (3.0d0 * ryij * ryij - d2 ) * dm5
                efg_ia( ia , 3 , 3 )  = efg_ia( ia , 3 , 3 ) - (3.0d0 * rzij * rzij - d2 ) * dm5
                efg_ia( ia , 1 , 2 )  = efg_ia( ia , 1 , 2 ) -  3.0d0 * rxij * ryij * dm5
                efg_ia( ia , 1 , 3 )  = efg_ia( ia , 1 , 3 ) -  3.0d0 * rxij * rzij * dm5
                efg_ia( ia , 2 , 3 )  = efg_ia( ia , 2 , 3 ) -  3.0d0 * ryij * rzij * dm5

              if ( lefg_it_contrib ) then
                it = itype(ja)  
                efg_ia_it( ia , it , 1 , 1 )  = efg_ia_it( ia , it , 1 , 1 ) - &
                ( 3.0d0 * rxij * rxij - d2 ) * dm5
                efg_ia_it( ia , it , 2 , 2 )  = efg_ia_it( ia , it , 2 , 2 ) - & 
                ( 3.0d0 * ryij * ryij - d2 ) * dm5
                efg_ia_it( ia , it , 3 , 3 )  = efg_ia_it( ia , it , 3 , 3 ) - & 
                ( 3.0d0 * rzij * rzij - d2 ) * dm5
                efg_ia_it( ia , it , 1 , 2 )  = efg_ia_it( ia , it , 1 , 2 ) - &
                  3.0d0 * rxij * ryij * dm5
                efg_ia_it( ia , it , 1 , 3 )  = efg_ia_it( ia , it , 1 , 3 ) - &
                  3.0d0 * rxij * rzij * dm5
                efg_ia_it( ia , it , 2 , 3 )  = efg_ia_it( ia , it , 2 , 3 ) - &
                  3.0d0 * ryij * rzij * dm5
              endif   

              endif ! d2.lt.cutefgsq

            endif ! ia.ne.ja
          enddo ! ja
        endif 

     enddo ! ncell

    efg_ia( ia , 2 , 1 ) = efg_ia( ia , 1 , 2 )
    efg_ia( ia , 3 , 1 ) = efg_ia( ia , 1 , 3 )
    efg_ia( ia , 3 , 2 ) = efg_ia( ia , 2 , 3 )
    if ( lefg_it_contrib ) then
      efg_ia_it( ia , : , 2 , 1 ) = efg_ia_it( ia , : , 1 , 2 )
      efg_ia_it( ia , : , 3 , 1 ) = efg_ia_it( ia , : , 1 , 3 )
      efg_ia_it( ia , : , 3 , 2 ) = efg_ia_it( ia , : , 2 , 3 )
    endif

  enddo atom
  !=========================================================
  !      END OF EFG TENSOR CALCULATION
  !=========================================================
#ifdef debug
CALL print_tensor( efg_ia( 1 , : , : ) , 'TOTCHG' )
#endif

  !======================================
  !  MERGE tensor from different proc
  !======================================
  CALL MPI_ALL_REDUCE_DOUBLE ( efg_ia ( : , 1 , 1 ) , natm ) 
  CALL MPI_ALL_REDUCE_DOUBLE ( efg_ia ( : , 2 , 2 ) , natm ) 
  CALL MPI_ALL_REDUCE_DOUBLE ( efg_ia ( : , 3 , 3 ) , natm ) 
  CALL MPI_ALL_REDUCE_DOUBLE ( efg_ia ( : , 1 , 2 ) , natm ) 
  CALL MPI_ALL_REDUCE_DOUBLE ( efg_ia ( : , 1 , 3 ) , natm ) 
  CALL MPI_ALL_REDUCE_DOUBLE ( efg_ia ( : , 2 , 3 ) , natm ) 
  efg_ia( : , 2, 1) = efg_ia( : , 1, 2)
  efg_ia( : , 3, 1) = efg_ia( : , 1, 3)
  efg_ia( : , 3, 2) = efg_ia( : , 2, 3)

  if ( lefg_it_contrib ) then
    do it=1,ntype
      CALL MPI_ALL_REDUCE_DOUBLE ( efg_ia_it ( : , it , 1 , 1 ) , natm ) 
      CALL MPI_ALL_REDUCE_DOUBLE ( efg_ia_it ( : , it , 2 , 2 ) , natm ) 
      CALL MPI_ALL_REDUCE_DOUBLE ( efg_ia_it ( : , it , 3 , 3 ) , natm ) 
      CALL MPI_ALL_REDUCE_DOUBLE ( efg_ia_it ( : , it , 1 , 2 ) , natm ) 
      CALL MPI_ALL_REDUCE_DOUBLE ( efg_ia_it ( : , it , 1 , 3 ) , natm ) 
      CALL MPI_ALL_REDUCE_DOUBLE ( efg_ia_it ( : , it , 2 , 3 ) , natm ) 
    enddo
    efg_ia_it( : , : , 2, 1) = efg_ia_it( : , : , 1, 2)
    efg_ia_it( : , : , 3, 1) = efg_ia_it( : , : , 1, 3)
    efg_ia_it( : , : , 3, 2) = efg_ia_it( : , : , 2, 3)
  endif

  ttt2 = MPI_WTIME(ierr)
  efgtimetot1 = efgtimetot1 + (ttt2-ttt1)

  CALL MPI_BARRIER( MPI_COMM_WORLD , ierr )

  return

END SUBROUTINE efg_DS

!*********************** SUBROUTINE efg_ES ************************************
!
! Ewald sum
! the implementation follows J. Chem. Phys, 112, 14, p 6512 (2000)
! Correction (16-05-11) following the J. Chem. Phys. 129, 074102 (2008)
! and quantum-espresso
!
! think about how to parallelize the real part loop 
! also the reciproc part should be parallized
!
!******************************************************************************
 
SUBROUTINE efg_ES ( iastart , iaend )

  USE control,    ONLY :  myrank , numprocs, calc
  USE config,     ONLY :  system , natm , natmi , atype , atypei , box , &
                          omega , itype , rx , ry , rz , ntype , qia 
  USE constants,  ONLY :  pi , fpi , piroot , imag
  USE field,      ONLY :  qch
  USE kspace,     ONLY :  struc_fact
  USE prop,       ONLY :  nprop_print , nprop
  USE time
  USE io_file 

  implicit none

  INCLUDE 'mpif.h'

  ! global
  integer, intent(in) :: iastart , iaend 

  ! local
  integer :: ia , ja ,  ierr , it 
  integer :: nxij , nyij , nzij
  double precision :: d , d2 , d4 , d3 , d5 , expon 
  double precision :: alpha2 , alpha3
  double precision :: allrealpart
  double precision :: T0 , T1 , T2  ! real part 
  double precision :: rxi , ryi , rzi , rxij , ryij , rzij 
  double precision :: invbox 
  double precision :: ak, kx, ky, kz, kk
  integer :: ik
  double precision :: kri 
  double complex   :: rhon , carg , recarg , recarg_dgg 
  double precision, external :: errfc 
  double precision :: ttt1 , ttt2 , ttt3 , ttt4 , ttt5 
  double precision :: efg_ia_real ( natm , 3 , 3 )
  double precision :: efg_ia_dual_real ( natm , 3 , 3 )
  double complex   :: efg_ia_dual ( natm , 3 , 3 )
  double precision :: efg_ia_real_it ( natm , ntype , 3 , 3 )
  double precision :: efg_ia_dual_real_it ( natm , ntype , 3 , 3 )
  double complex   :: efg_ia_dual_it ( natm , ntype , 3 , 3 )

#ifdef debug
  CALL print_config_sample(0,0)
#endif
  ! ==========================
  !  init some quantities
  ! ==========================
  efg_ia_real         = 0.0d0
  efg_ia_dual         = (0.0d0,0.0d0)
  efg_ia_dual_real    = 0.0d0
  efg_ia              = 0.0d0
  efg_ia_real_it      = 0.0d0
  efg_ia_dual_it      = (0.0d0,0.0d0)
  efg_ia_dual_real_it = 0.0d0
  efg_ia_it           = 0.0d0

  ! =================
  !  some constants 
  ! =================
  !onethird = 1/3.0D0
  invbox = 1.0d0 / box

  ! related ewald parameter
  alpha2 = alphaES * alphaES      
  alpha3 = alpha2  * alphaES


  ttt1 = MPI_WTIME(ierr)

  ! =====================
  ! facteur de structure 
  ! =====================
  CALL struc_fact ( km_efg )

  ttt2 = MPI_WTIME(ierr)
  strftimetot = strftimetot + (ttt2 - ttt1)
! ==============================================
!        direct space part
! ==============================================

atom1: do ia = iastart , iaend 
     rxi = rx(ia)
     ryi = ry(ia)
     rzi = rz(ia)
#ifdef fix_grid
     rxi = rgrid(1,ia)
     ryi = rgrid(2,ia)
     rzi = rgrid(3,ia) 
#endif

     do ja = 1, natm

       if (ja .ne. ia ) then

         rxij = rxi - rx(ja)
         ryij = ryi - ry(ja)
         rzij = rzi - rz(ja)
         nxij = nint( rxij * invbox )
         nyij = nint( ryij * invbox )
         nzij = nint( rzij * invbox )
         rxij = rxij - box * nxij
         ryij = ryij - box * nyij
         rzij = rzij - box * nzij
  
         d2 = rxij * rxij + ryij * ryij + rzij * rzij
         d = dsqrt( d2 )
         d3 = d2 * d
         d4 = d2 * d2
         d5 = d3 * d2
         expon = dexp( - alpha2 * d2 ) / piroot
         T0 = errfc( alphaES * d ) / d5
         T1 = ( 2.0d0 * alphaES ) / d4 
         T2 = ( 4.0d0 * alpha3  ) / d2 / 3.0d0
         allrealpart = qia(ja) * ( T0 + ( T1 + T2 ) *expon )

         efg_ia_real( ia , 1 , 1 ) = efg_ia_real( ia , 1 , 1 ) - &
         ( 3.0D0 * rxij * rxij - d2 ) * allrealpart
         efg_ia_real( ia , 2 , 2 ) = efg_ia_real( ia , 2 , 2 ) - &
         ( 3.0D0 * ryij * ryij - d2 ) * allrealpart
         efg_ia_real( ia , 3 , 3 ) = efg_ia_real( ia , 3 , 3 ) - &
         ( 3.0D0 * rzij * rzij - d2 ) * allrealpart
         efg_ia_real( ia , 1 , 2 ) = efg_ia_real( ia , 1 , 2 ) - &
           3.0d0 * rxij * ryij * allrealpart
         efg_ia_real( ia , 1 , 3 ) = efg_ia_real( ia , 1 , 3 ) - &
           3.0d0 * rxij * rzij * allrealpart
         efg_ia_real( ia , 2 , 3 ) = efg_ia_real( ia , 2 , 3 ) - &
           3.0d0 * ryij * rzij * allrealpart

         if ( lefg_it_contrib ) then
           it = itype(ja)
           efg_ia_real_it( ia , it , 1 , 1 ) = efg_ia_real_it( ia , it , 1 , 1 ) - &
           ( 3.0D0 * rxij * rxij - d2 ) * allrealpart
           efg_ia_real_it( ia , it , 2 , 2 ) = efg_ia_real_it( ia , it , 2 , 2 ) - &
           ( 3.0D0 * ryij * ryij - d2 ) * allrealpart
           efg_ia_real_it( ia , it , 3 , 3 ) = efg_ia_real_it( ia , it , 3 , 3 ) - &
           ( 3.0D0 * rzij * rzij - d2 ) * allrealpart
           efg_ia_real_it( ia , it , 1 , 2 ) = efg_ia_real_it( ia , it , 1 , 2 ) - &
           3.0d0 * rxij * ryij * allrealpart
           efg_ia_real_it( ia , it , 1 , 3 ) = efg_ia_real_it( ia , it , 1 , 3 ) - &
           3.0d0 * rxij * rzij * allrealpart
           efg_ia_real_it( ia , it , 2 , 3 ) = efg_ia_real_it( ia , it , 2 , 3 ) - &
           3.0d0 * ryij * rzij * allrealpart
         endif

       endif

     enddo

  enddo atom1

  ttt3 = MPI_WTIME(ierr)
  efgtimetot1 = efgtimetot1 + (ttt3-ttt2)

! ==============================================
!            reciprocal space part
! ==============================================
  ! sum on atoms 
atom2:  do ia = 1 , natm
    rxi = rx ( ia )
    ryi = ry ( ia )
    rzi = rz ( ia )
#ifdef fix_grid
     rxi = rgrid ( 1 , ia )
     ryi = rgrid ( 2 , ia )
     rzi = rgrid ( 3 , ia ) 
#endif
kpoint: do ik = 1, km_efg%nkcut 
      ! =================
      !   k-space  
      ! =================
      kx = km_efg%kpt ( 1 , ik )
      ky = km_efg%kpt ( 2 , ik )
      kz = km_efg%kpt ( 3 , ik )
      kk = km_efg%kptk( ik )
      Ak = dexp( - kk * 0.25d0 / alpha2 ) 
      if ( km_efg%kptk( ik ) .eq. 0 ) then 
        WRITE ( stdout , * ) 'the sum should be done on k! =  0',ik
        STOP 
      endif

      ! ===============================
      !                              ---
      !  charge density in k-space ( \   q * facteur de structure  )
      !                              /__
      ! ===============================
      rhon = (0.d0, 0.d0)
      do it = 1, ntype
        rhon = rhon + qch(it) * CONJG( km_efg%strf ( ik , it ) )
      enddo

      kri = ( kx * rxi + ky * ryi + kz * rzi ) 
      carg = EXP ( imag * kri )
      recarg_dgg =  rhon * carg * Ak / kk 
      recarg = recarg_dgg * kk 

      efg_ia_dual ( ia , 1 , 1 ) = efg_ia_dual ( ia , 1 , 1 ) +  3.0d0 * kx * kx * recarg_dgg - recarg
      efg_ia_dual ( ia , 2 , 2 ) = efg_ia_dual ( ia , 2 , 2 ) +  3.0d0 * ky * ky * recarg_dgg - recarg
      efg_ia_dual ( ia , 3 , 3 ) = efg_ia_dual ( ia , 3 , 3 ) +  3.0d0 * kz * kz * recarg_dgg - recarg
      efg_ia_dual ( ia , 1 , 2 ) = efg_ia_dual ( ia , 1 , 2 ) +  3.0d0 * kx * ky * recarg_dgg
      efg_ia_dual ( ia , 1 , 3 ) = efg_ia_dual ( ia , 1 , 3 ) +  3.0d0 * kx * kz * recarg_dgg
      efg_ia_dual ( ia , 2 , 3 ) = efg_ia_dual ( ia , 2 , 3 ) +  3.0d0 * ky * kz * recarg_dgg

      if ( lefg_it_contrib ) then
        
        do it =1 , ntype
          ! charge density in k-space it specific
!          rhon = qch(it) * CONJG( km_efg%strf ( ik , it ) )
          rhon = qch(it) * CONJG( km_efg_dip%strf ( ik , it ) )

          recarg_dgg =  rhon*carg*ak / kk
          recarg = recarg_dgg * kk

          efg_ia_dual_it ( ia , it , 1 , 1 ) = efg_ia_dual_it ( ia , it , 1 , 1 ) +  &
          3.0d0 * kx * kx * recarg_dgg - recarg
          efg_ia_dual_it ( ia , it , 2 , 2 ) = efg_ia_dual_it ( ia , it , 2 , 2 ) +  &
          3.0d0 * ky * ky * recarg_dgg - recarg
          efg_ia_dual_it ( ia , it , 3 , 3 ) = efg_ia_dual_it ( ia , it , 3 , 3 ) +  &
          3.0d0 * kz * kz * recarg_dgg - recarg
          efg_ia_dual_it ( ia , it , 1 , 2 ) = efg_ia_dual_it ( ia , it , 1 , 2 ) +  & 
          3.0d0 * kx * ky * recarg_dgg
          efg_ia_dual_it ( ia , it , 1 , 3 ) = efg_ia_dual_it ( ia , it , 1 , 3 ) +  & 
          3.0d0 * kx * kz * recarg_dgg
          efg_ia_dual_it ( ia , it , 2 , 3 ) = efg_ia_dual_it ( ia , it , 2 , 3 ) +  &
          3.0d0 * ky * kz * recarg_dgg
          
        enddo
      endif

    enddo kpoint

  enddo atom2

  ttt4 = MPI_WTIME(ierr)
  efgtimetot2 = efgtimetot2 + (ttt4-ttt3)

!=============================
!  MERGE REAL PART
!=============================
  CALL MPI_ALL_REDUCE_DOUBLE ( efg_ia_real ( : , 1 , 1 ) , natm ) 
  CALL MPI_ALL_REDUCE_DOUBLE ( efg_ia_real ( : , 2 , 2 ) , natm ) 
  CALL MPI_ALL_REDUCE_DOUBLE ( efg_ia_real ( : , 3 , 3 ) , natm ) 
  CALL MPI_ALL_REDUCE_DOUBLE ( efg_ia_real ( : , 1 , 2 ) , natm ) 
  CALL MPI_ALL_REDUCE_DOUBLE ( efg_ia_real ( : , 1 , 3 ) , natm ) 
  CALL MPI_ALL_REDUCE_DOUBLE ( efg_ia_real ( : , 2 , 3 ) , natm ) 

 ! in fact we don't need this as the EFG tensor is completely define from the upper diagonal value
  efg_ia_real( : , 2 , 1) = efg_ia_real( : , 1 , 2) 
  efg_ia_real( : , 3 , 1) = efg_ia_real( : , 1 , 3) 
  efg_ia_real( : , 3 , 2) = efg_ia_real( : , 2 , 3)

  efg_ia_dual( : , 2 , 1) = efg_ia_dual( : , 1 , 2) 
  efg_ia_dual( : , 3 , 1) = efg_ia_dual( : , 1 , 3) 
  efg_ia_dual( : , 3 , 2) = efg_ia_dual( : , 2 , 3)

  if ( lefg_it_contrib ) then

    do it=1,ntype
      CALL MPI_ALL_REDUCE_DOUBLE ( efg_ia_real_it ( : , it , 1 , 1 ) , natm )
      CALL MPI_ALL_REDUCE_DOUBLE ( efg_ia_real_it ( : , it , 2 , 2 ) , natm )
      CALL MPI_ALL_REDUCE_DOUBLE ( efg_ia_real_it ( : , it , 3 , 3 ) , natm )
      CALL MPI_ALL_REDUCE_DOUBLE ( efg_ia_real_it ( : , it , 1 , 2 ) , natm )
      CALL MPI_ALL_REDUCE_DOUBLE ( efg_ia_real_it ( : , it , 1 , 3 ) , natm )
      CALL MPI_ALL_REDUCE_DOUBLE ( efg_ia_real_it ( : , it , 2 , 3 ) , natm )
    enddo

    efg_ia_real_it( : , : , 2, 1) = efg_ia_real_it( : , : , 1, 2) 
    efg_ia_real_it( : , : , 3, 1) = efg_ia_real_it( : , : , 1, 3) 
    efg_ia_real_it( : , : , 3, 2) = efg_ia_real_it( : , : , 2, 3)

    efg_ia_dual_it( : , : , 2, 1) = efg_ia_dual_it( : , : , 1, 2) 
    efg_ia_dual_it( : , : , 3, 1) = efg_ia_dual_it( : , : , 1, 3) 
    efg_ia_dual_it( : , : , 3, 2) = efg_ia_dual_it( : , : , 2, 3)

  endif

  ! =======
  ! 4pi/3V
  ! =======
  efg_ia_dual_real( : , : , :) = dble( efg_ia_dual( : , : , :)  ) 
  efg_ia_dual_real =  efg_ia_dual_real * fpi / omega / 3.0d0

  if ( lefg_it_contrib ) then
    ! =======
    ! 4pi/3V
    ! =======
    efg_ia_dual_real_it( : , : , : , :) = dble( efg_ia_dual_it( : , : , : , :)  ) 
    efg_ia_dual_real_it =  efg_ia_dual_real_it * fpi / omega / 3.0d0
  endif

  ! ==============
  !  total tensor
  ! ==============
  efg_ia = efg_ia_dual_real + efg_ia_real 
  if ( lefg_it_contrib ) then
    efg_ia_it = efg_ia_dual_real_it + efg_ia_real_it 
  endif
!=========================================================
!      END OF EFG TENSOR CALCULATION
!=========================================================

  ttt5 = MPI_WTIME(ierr)
  efgtimetot3 = efgtimetot3 + (ttt5-ttt4)

  CALL MPI_BARRIER( MPI_COMM_WORLD , ierr )

  return

END SUBROUTINE efg_ES

!*********************** SUBROUTINE efg_DS_dip ********************************
!
! calculate EFG from dipoles with direct method
! the implementation follows J. Chem. Phys, 112, 14, p 6512 (2000)
!
!******************************************************************************

SUBROUTINE efg_DS_dip ( iastart , iaend )

  USE control,    ONLY :  calc
  USE config,     ONLY :  ntype , atype , atypei , natm , natmi , box  , system , &
                          rx , ry , rz , qia , dipia , omega , ipolar
  USE constants,  ONLY :  piroot , imag , fpi
  USE kspace,     ONLY :  struc_fact
  USE field,      ONLY :  dip , dip_ind 
  USE time
  USE io_file


  implicit none
  INCLUDE 'mpif.h'

  ! global
  integer, intent(in) :: iastart , iaend 

  ! local
  integer :: ia , ja , ncell 
  integer :: jstart , jend 
  integer :: ierr
  double precision :: d , d2 , d5 , d7 
  double precision :: cutefgsq
  double precision :: rxi , ryi , rzi , rxij , ryij , rzij
  double precision :: rxj , ryj , rzj 
  double precision :: Txxx , Txxy , Txxz 
  double precision :: Tyyx , Tyyy , Tyyz 
  double precision :: Tzzx , Tzzy , Tzzz 
  double precision :: Txyx , Txyy , Txyz 
  double precision :: Txzx , Txzy , Txzz 
  double precision :: Tyzx , Tyzy , Tyzz 
  double precision :: dipx , dipy , dipz 

  double precision :: ttt1 , ttt2 , ttt3 

#ifdef debug
  CALL print_config_sample(0,0)
#endif

  ! ==========================
  !  init some quantities 
  ! ==========================
  cutefgsq = cutefg * cutefg
  efg_ia = 0.0d0
  efg_ia_it = 0.0d0

#ifdef debug
     write(*,*) 'iastart iaend',iastart , iaend
     write(*,*) 'rm_efg%ncmax',rm_efg%ncmax
#endif

ttt1 = MPI_WTIME(ierr)

! ==============================================
!        direct space part
! ==============================================

atom : do ia = iastart , iaend
    rxi = rx(ia)
    ryi = ry(ia)
    rzi = rz(ia)
    ! ==================================================
    ! sum over neighboring cells (see direct_sum_init )
    ! ==================================================
    do ncell = 1 , rm_efg%ncmax
         ! ==============================
         !  ia and ja in different cells
         ! ==============================
         if ( rm_efg%lcell(ncell) .eq. 1) then

          do ja = 1 , natm
            dipx = 3.0d0 * dipia( ja , 1 )
            dipy = 3.0d0 * dipia( ja , 2 )
            dipz = 3.0d0 * dipia( ja , 3 )
            rxj  = rx(ja) + rm_efg%boxxyz(1,ncell)
            ryj  = ry(ja) + rm_efg%boxxyz(2,ncell)
            rzj  = rz(ja) + rm_efg%boxxyz(3,ncell)
            rxij = rxi - rxj
            ryij = ryi - ryj
            rzij = rzi - rzj
            d2   = rxij * rxij + ryij * ryij + rzij * rzij

            if ( d2 .lt. cutefgsq ) then
              d = dsqrt(d2)
              d5 = d2 * d2 * d
              d7 = d5 * d2

         Txxx = 5.0d0 * rxij * rxij * rxij -  3.0d0 * d2 * ( rxij )
         Txxy = 5.0d0 * rxij * rxij * ryij -          d2 * ( ryij )
         Txxz = 5.0d0 * rxij * rxij * rzij -          d2 * ( rzij )
         efg_ia( ia , 1 , 1 ) = efg_ia( ia , 1 , 1 ) -  ( Txxx * dipx + &
                                                          Txxy * dipy + &
                                                          Txxz * dipz ) / d7 

         Tyyx = 5.0d0 * ryij * ryij * rxij -          d2 * ( rxij )
         Tyyy = 5.0d0 * ryij * ryij * ryij -  3.0d0 * d2 * ( ryij )
         Tyyz = 5.0d0 * ryij * ryij * rzij -          d2 * ( rzij )
         efg_ia( ia , 2 , 2 ) = efg_ia( ia , 2 , 2 ) -  ( Tyyx * dipx + &
                                                          Tyyy * dipy + &
                                                          Tyyz * dipz ) / d7 

         Tzzx = 5.0d0 * rzij * rzij * rxij -          d2 * ( rxij )
         Tzzy = 5.0d0 * rzij * rzij * ryij -  3.0d0 * d2 * ( ryij )
         Tzzz = 5.0d0 * rzij * rzij * rzij -          d2 * ( rzij )
         efg_ia( ia , 3 , 3 ) = efg_ia( ia , 3 , 3 ) -  ( Tzzx * dipx + &
                                                          Tzzy * dipy + &
                                                          Tzzz * dipz ) / d7

         Txyx = 5.0d0 * rxij * ryij * rxij -          d2 * ( ryij )
         Txyy = 5.0d0 * rxij * ryij * ryij -          d2 * ( rxij )
         Txyz = 5.0d0 * rxij * ryij * rzij 
         efg_ia( ia , 1 , 2 ) = efg_ia( ia , 1 , 2 ) -  ( Txyx * dipx + &
                                                          Txyy * dipy + &
                                                          Txyz * dipz ) / d7

         Txzx = 5.0d0 * rxij * rzij * rxij -          d2 * ( rzij )
         Txzy = 5.0d0 * rxij * rzij * ryij 
         Txzz = 5.0d0 * rxij * rzij * rzij -          d2 * ( rxij )
         efg_ia( ia , 1 , 3 ) = efg_ia( ia , 1 , 3 ) -  ( Txzx * dipx + &
                                                          Txzy * dipy + &
                                                          Txzz * dipz ) / d7

         Tyzx = 5.0d0 * ryij * rzij * rxij 
         Tyzy = 5.0d0 * ryij * rzij * ryij -          d2 * ( rzij)
         Tyzz = 5.0d0 * ryij * rzij * rzij -          d2 * ( ryij )
         efg_ia( ia , 2 , 3 ) = efg_ia( ia , 2 , 3 ) -  ( Tyzx * dipx + &
                                                          Tyzy * dipy + &
                                                          Tyzz * dipz ) / d7

            endif ! d2.lt.cutefgsq
          enddo ! ja
        endif 
         ! =======================================
         !  ia and ja in the same cell (ia.ne.ja)
         ! =======================================
        if ( rm_efg%lcell(ncell) .eq. 0) then

          do ja = 1 , natm

            dipx = 3.0d0 * dipia( ja , 1 )
            dipy = 3.0d0 * dipia( ja , 2 )
            dipz = 3.0d0 * dipia( ja , 3 )

            if (ja.ne.ia) then
              rxj  = rx(ja)
              ryj  = ry(ja)
              rzj  = rz(ja)
              rxij = rxi - rxj
              ryij = ryi - ryj
              rzij = rzi - rzj
              d2   = rxij * rxij + ryij * ryij + rzij * rzij
 
              if (d2.lt.cutefgsq) then
                d = dsqrt(d2)
                d5 = d2 * d2 * d
              d7 = d5 * d2
         Txxx = 5.0d0 * rxij * rxij * rxij -  3.0d0 * d2 * ( rxij )
         Txxy = 5.0d0 * rxij * rxij * ryij -          d2 * ( ryij )
         Txxz = 5.0d0 * rxij * rxij * rzij -          d2 * ( rzij )
         efg_ia( ia , 1 , 1 ) = efg_ia( ia , 1 , 1 ) -  ( Txxx * dipx + &
                                                          Txxy * dipy + &
                                                          Txxz * dipz ) / d7 

         Tyyx = 5.0d0 * ryij * ryij * rxij -          d2 * ( rxij )
         Tyyy = 5.0d0 * ryij * ryij * ryij -  3.0d0 * d2 * ( ryij )
         Tyyz = 5.0d0 * ryij * ryij * rzij -          d2 * ( rzij )
         efg_ia( ia , 2 , 2 ) = efg_ia( ia , 2 , 2 ) -  ( Tyyx * dipx + &
                                                          Tyyy * dipy + &
                                                          Tyyz * dipz ) / d7 

         Tzzx = 5.0d0 * rzij * rzij * rxij -          d2 * ( rxij )
         Tzzy = 5.0d0 * rzij * rzij * ryij -  3.0d0 * d2 * ( ryij )
         Tzzz = 5.0d0 * rzij * rzij * rzij -          d2 * ( rzij )
         efg_ia( ia , 3 , 3 ) = efg_ia( ia , 3 , 3 ) -  ( Tzzx * dipx + &
                                                          Tzzy * dipy + &
                                                          Tzzz * dipz ) / d7

         Txyx = 5.0d0 * rxij * ryij * rxij -          d2 * ( ryij )
         Txyy = 5.0d0 * rxij * ryij * ryij -          d2 * ( rxij )
         Txyz = 5.0d0 * rxij * ryij * rzij 
         efg_ia( ia , 1 , 2 ) = efg_ia( ia , 1 , 2 ) -  ( Txyx * dipx + &
                                                          Txyy * dipy + &
                                                          Txyz * dipz ) / d7

         Txzx = 5.0d0 * rxij * rzij * rxij -          d2 * ( rzij )
         Txzy = 5.0d0 * rxij * rzij * ryij 
         Txzz = 5.0d0 * rxij * rzij * rzij -          d2 * ( rxij )
         efg_ia( ia , 1 , 3 ) = efg_ia( ia , 1 , 3 ) -  ( Txzx * dipx + &
                                                          Txzy * dipy + &
                                                          Txzz * dipz ) / d7

         Tyzx = 5.0d0 * ryij * rzij * rxij 
         Tyzy = 5.0d0 * ryij * rzij * ryij -          d2 * ( rzij)
         Tyzz = 5.0d0 * ryij * rzij * rzij -          d2 * ( ryij )
         efg_ia( ia , 2 , 3 ) = efg_ia( ia , 2 , 3 ) -  ( Tyzx * dipx + &
                                                          Tyzy * dipy + &
                                                          Tyzz * dipz ) / d7

              endif ! d2.lt.cutefgsq

            endif ! ia.ne.ja
          enddo ! ja
        endif 

     enddo ! ncell

    efg_ia( ia , 2 , 1 ) = efg_ia( ia , 1 , 2 )
    efg_ia( ia , 3 , 1 ) = efg_ia( ia , 1 , 3 )
    efg_ia( ia , 3 , 2 ) = efg_ia( ia , 2 , 3 )

  enddo atom

ttt2 = MPI_WTIME(ierr)
efgtimetot1 = efgtimetot1 + (ttt2-ttt1)

!=============================
!  MERGE REAL PART
!=============================
  CALL MPI_ALL_REDUCE_DOUBLE ( efg_ia ( : , 1 , 1 ) , natm ) 
  CALL MPI_ALL_REDUCE_DOUBLE ( efg_ia ( : , 2 , 2 ) , natm ) 
  CALL MPI_ALL_REDUCE_DOUBLE ( efg_ia ( : , 3 , 3 ) , natm ) 
  CALL MPI_ALL_REDUCE_DOUBLE ( efg_ia ( : , 1 , 2 ) , natm ) 
  CALL MPI_ALL_REDUCE_DOUBLE ( efg_ia ( : , 1 , 3 ) , natm ) 
  CALL MPI_ALL_REDUCE_DOUBLE ( efg_ia ( : , 2 , 3 ) , natm ) 

  ! not needed ... just for consistency
  efg_ia( : , 2 , 1) = efg_ia( : , 1 , 2) 
  efg_ia( : , 3 , 1) = efg_ia( : , 1 , 3) 
  efg_ia( : , 3 , 2) = efg_ia( : , 2 , 3)

!  ! ==============
!  !  total tensor
!  ! ==============
#ifdef debug
  CALL print_tensor( efg_ia( 1 , : , : ) , 'TOTDIP' )
#endif
!!=========================================================
!!      END OF EFG TENSOR CALCULATION
!!=========================================================
!
  ttt3 = MPI_WTIME(ierr)
  efgtimetot3 = efgtimetot3 + (ttt3 - ttt2)

  CALL MPI_BARRIER( MPI_COMM_WORLD , ierr )

  return

END SUBROUTINE efg_DS_dip


!*********************** SUBROUTINE efg_ES_dip  ********************************
!
! calculate EFG from dipoles
! the implementation follows J. Chem. Phys, 112, 14, p 6512 (2000)
!
!
! WARNING : maybe non sense !!
! Observation if we calculate electric field gradient from dipoles, time reversal symmetry 
! is breaking. Also -k and +k symmetry is not conserved. One should only sum on positive k-points 
! in reciprocal space
!
!
!******************************************************************************

SUBROUTINE efg_ES_dip ( iastart , iaend )

  USE control,    ONLY :  calc
  USE config,     ONLY :  ntype , atype , atypei , natm , natmi , box  , system , &
                          rx , ry , rz , qia , dipia , omega
  USE constants,  ONLY :  piroot , imag , fpi , tpi
  USE kspace,     ONLY :  struc_fact , struc_fact_dip 
  USE field,      ONLY :  dip 
  USE time
  USE io_file


  implicit none
  INCLUDE 'mpif.h'

  ! global
  integer, intent(in) :: iastart , iaend 

  ! local
  integer :: ia , ja , it 
  integer :: nxij , nyij , nzij
  integer :: ierr
  double precision :: d , d2 , d4 , d6, d7 , expon
  double precision :: alpha2 , alpha3 , alpha5
  double precision :: allrealpart
  double precision :: invbox 
  double precision :: T0 , T1 , T2 , T3 ! real part
  double precision, external :: errfc
  double precision :: rxi , ryi , rzi , rxij , ryij , rzij
  double precision :: Txxx , Txxy , Txxz 
  double precision :: Tyyx , Tyyy , Tyyz 
  double precision :: Tzzx , Tzzy , Tzzz 
  double precision :: Txyx , Txyy , Txyz 
  double precision :: Txzx , Txzy , Txzz 
  double precision :: Tyzx , Tyzy , Tyzz 
  double precision :: dipx , dipy , dipz 
  double precision :: ak, kx, ky, kz, kk
  integer :: ik
  double precision :: kri 
  double complex   :: rhon , carg , recarg , recarg_dgg
  double precision :: efg_ia_real ( natm , 3 , 3 )
  double precision :: efg_ia_dual_real ( natm , 3 , 3 )
  double complex   :: efg_ia_dual ( natm , 3 , 3 )
  double precision :: efg_ia_real_it ( natm , ntype , 3 , 3 )
  double precision :: efg_ia_dual_real_it ( natm , ntype , 3 , 3 )
  double complex   :: efg_ia_dual_it ( natm , ntype , 3 , 3 )

  double precision :: ttt1 , ttt2 , ttt3 , ttt4 , ttt5 

#ifdef debug
  CALL print_config_sample(0,0)
  write(*,*) 'iastart iaend   ',iastart , iaend
  write(*,*) 'km_efg%nkcut    ',km_efg%nkcut
  write(*,*) 'km_efg_dip%nkcut',km_efg_dip%nkcut
#endif

  ! ==========================
  !  init some quantities 
  ! ==========================
  efg_ia_real         = 0.0d0
  efg_ia_dual         = (0.0d0,0.0d0)
  efg_ia_dual_real    = 0.0d0
  efg_ia              = 0.0d0
  efg_ia_real_it      = 0.0d0
  efg_ia_dual_it      = (0.0d0,0.0d0)
  efg_ia_dual_real_it = 0.0d0
  efg_ia_it           = 0.0d0

  ! =================
  !  some constants 
  ! =================
!  onethird = 1/3.0D0
  invbox = 1.0d0 / box

  ! related ewald parameter
  alpha2 = alphaES * alphaES      
  alpha3 = alpha2  * alphaES
  alpha5 = alpha2  * alpha3


  ttt1 = MPI_WTIME(ierr)

  ! =====================
  ! facteur de structure 
  ! =====================
  CALL struc_fact ( km_efg )

  ttt2 = MPI_WTIME(ierr)
  strftimetot = strftimetot + (ttt2 - ttt1)
! ==============================================
!        direct space part
! ==============================================

atom1: do ia = iastart , iaend 
     rxi = rx(ia)
     ryi = ry(ia)
     rzi = rz(ia)

     do ja = 1, natm

       if (ja .ne. ia ) then

         dipx = 3.0d0 * dipia( ja , 1 )
         dipy = 3.0d0 * dipia( ja , 2 )
         dipz = 3.0d0 * dipia( ja , 3 )

         rxij = rxi - rx(ja)
         ryij = ryi - ry(ja)
         rzij = rzi - rz(ja)
         nxij = nint( rxij * invbox )
         nyij = nint( ryij * invbox )
         nzij = nint( rzij * invbox )
         rxij = rxij - box * nxij
         ryij = ryij - box * nyij
         rzij = rzij - box * nzij
  
         d2 = rxij * rxij + ryij * ryij + rzij * rzij
         d = dsqrt( d2 )
         d4 = d2 * d2
         d6 = d4 * d2
         d7 = d6 * d
         expon = dexp( - alpha2 * d2 ) / piroot

         T0 = errfc( alphaES * d ) / d7
         T1 = ( 2.0d0 * alphaES ) / d6 
         T2 = ( 4.0d0 * alpha3  ) / d4 / 3.0d0
         T3 = ( 8.0d0 * alpha5  ) / d2 / 15.0d0
         allrealpart =  ( T0 +  ( T1 + T2 + T3) * expon )


         Txxx = 5.0d0 * rxij * rxij * rxij -  3.0d0 * d2 * ( rxij )
         Txxy = 5.0d0 * rxij * rxij * ryij -          d2 * ( ryij )
         Txxz = 5.0d0 * rxij * rxij * rzij -          d2 * ( rzij )
         efg_ia_real( ia , 1 , 1 ) = efg_ia_real( ia , 1 , 1 ) -  ( Txxx * dipx + &
                                                                    Txxy * dipy + &
                                                                    Txxz * dipz ) * allrealpart

         Tyyx = 5.0d0 * ryij * ryij * rxij -          d2 * ( rxij )
         Tyyy = 5.0d0 * ryij * ryij * ryij -  3.0d0 * d2 * ( ryij )
         Tyyz = 5.0d0 * ryij * ryij * rzij -          d2 * ( rzij )
         efg_ia_real( ia , 2 , 2 ) = efg_ia_real( ia , 2 , 2 ) -  ( Tyyx * dipx + &
                                                                    Tyyy * dipy + &
                                                                    Tyyz * dipz ) * allrealpart

         Tzzx = 5.0d0 * rzij * rzij * rxij -          d2 * ( rxij )
         Tzzy = 5.0d0 * rzij * rzij * ryij -  3.0d0 * d2 * ( ryij )
         Tzzz = 5.0d0 * rzij * rzij * rzij -          d2 * ( rzij )
         efg_ia_real( ia , 3 , 3 ) = efg_ia_real( ia , 3 , 3 ) -  ( Tzzx * dipx + &
                                                                    Tzzy * dipy + &
                                                                    Tzzz * dipz ) * allrealpart

         Txyx = 5.0d0 * rxij * ryij * rxij -          d2 * ( ryij )
         Txyy = 5.0d0 * rxij * ryij * ryij -          d2 * ( rxij )
         Txyz = 5.0d0 * rxij * ryij * rzij 
         efg_ia_real( ia , 1 , 2 ) = efg_ia_real( ia , 1 , 2 ) -  ( Txyx * dipx + &
                                                                    Txyy * dipy + &
                                                                    Txyz * dipz ) * allrealpart

         Txzx = 5.0d0 * rxij * rzij * rxij -          d2 * ( rzij )
         Txzy = 5.0d0 * rxij * rzij * ryij 
         Txzz = 5.0d0 * rxij * rzij * rzij -          d2 * ( rxij )
         efg_ia_real( ia , 1 , 3 ) = efg_ia_real( ia , 1 , 3 ) -  ( Txzx * dipx + &
                                                                    Txzy * dipy + &
                                                                    Txzz * dipz ) * allrealpart

         Tyzx = 5.0d0 * ryij * rzij * rxij 
         Tyzy = 5.0d0 * ryij * rzij * ryij -          d2 * ( rzij)
         Tyzz = 5.0d0 * ryij * rzij * rzij -          d2 * ( ryij )
         efg_ia_real( ia , 2 , 3 ) = efg_ia_real( ia , 2 , 3 ) -  ( Tyzx * dipx + &
                                                                    Tyzy * dipy + &
                                                                    Tyzz * dipz ) * allrealpart
       endif

     enddo

  enddo atom1

  ttt3 = MPI_WTIME(ierr)
  efgtimetot1 = efgtimetot1 + (ttt3-ttt2)

! ==============================================
!            reciprocal space part
! ==============================================
! sum on atoms 
atom2:  do ia = 1 , natm 
    rxi = rx ( ia )
    ryi = ry ( ia )
    rzi = rz ( ia )
    dipx = 3.0d0 * dipia( ia , 1 )
    dipy = 3.0d0 * dipia( ia , 2 )
    dipz = 3.0d0 * dipia( ia , 3 )
    kpoint : do ik = 1, km_efg%nkcut 
      ! =================
      !   k-space  
      ! =================
      kx = km_efg%kpt ( 1 , ik )
      ky = km_efg%kpt ( 2 , ik )
      kz = km_efg%kpt ( 3 , ik )
      kk = km_efg%kptk( ik )
      Ak = dexp( - kk * 0.25d0 / alpha2 ) 
      if ( km_efg%kptk( ik ) .eq. 0 ) then 
        WRITE ( stdout , * ) 'the sum should be done on k!=  0',ik
        STOP 
      endif

      ! ===============================
      !                              ---
      !  charge density in k-space ( \   q * facteur de structure  )
      !                              /__
      ! ===============================
      rhon = (0.d0, 0.d0)
      do it = 1, ntype
        rhon = rhon + imag * ( dip ( it , 1 )* kx + dip ( it , 2 )* ky + dip ( it , 3 )* kz ) *  CONJG ( km_efg%strf ( ik , it ) )
      enddo

      kri        = ( kx * rxi + ky * ryi + kz * rzi ) 
      carg       = EXP ( imag * kri )
      recarg_dgg = rhon * carg * Ak / kk 
      
#ifdef debug2
      if ( mod ( ik, 20 ) .eq. 0 ) then
        write(*,'(a)') '              recarg_dgg                         rhon                            carg                     Ak            kk                   rhon * carg                         efg'
        write(*,'(a)') '          real           imag            real            imag            real            imag            real           real            real              imag            real           imag' 
      endif
#endif 

      efg_ia_dual( ia , 1 , 1 ) = efg_ia_dual( ia , 1 , 1 ) +  kx * kx * recarg_dgg 
      efg_ia_dual( ia , 2 , 2 ) = efg_ia_dual( ia , 2 , 2 ) +  ky * ky * recarg_dgg 
      efg_ia_dual( ia , 3 , 3 ) = efg_ia_dual( ia , 3 , 3 ) +  kz * kz * recarg_dgg 
      efg_ia_dual( ia , 1 , 2 ) = efg_ia_dual( ia , 1 , 2 ) +  kx * ky * recarg_dgg 
      efg_ia_dual( ia , 1 , 3 ) = efg_ia_dual( ia , 1 , 3 ) +  kx * kz * recarg_dgg 
      efg_ia_dual( ia , 2 , 3 ) = efg_ia_dual( ia , 2 , 3 ) +  ky * kz * recarg_dgg 

#ifdef debug2
      write(*,'(12f16.8)') recarg_dgg , rhon , carg , Ak , kk , rhon * carg , efg_ia_dual( ia , 1 , 1 )
#endif

    enddo kpoint

  enddo atom2

  ttt4 = MPI_WTIME(ierr)
  efgtimetot2 = efgtimetot2 + (ttt4-ttt3)

!=============================
!  MERGE REAL PART
!=============================
  CALL MPI_ALL_REDUCE_DOUBLE ( efg_ia_real ( : , 1 , 1 ) , natm ) 
  CALL MPI_ALL_REDUCE_DOUBLE ( efg_ia_real ( : , 2 , 2 ) , natm ) 
  CALL MPI_ALL_REDUCE_DOUBLE ( efg_ia_real ( : , 3 , 3 ) , natm ) 
  CALL MPI_ALL_REDUCE_DOUBLE ( efg_ia_real ( : , 1 , 2 ) , natm ) 
  CALL MPI_ALL_REDUCE_DOUBLE ( efg_ia_real ( : , 1 , 3 ) , natm ) 
  CALL MPI_ALL_REDUCE_DOUBLE ( efg_ia_real ( : , 2 , 3 ) , natm ) 

  ! not needed ... just for consistency
  efg_ia_real( : , 2 , 1) = efg_ia_real( : , 1 , 2) 
  efg_ia_real( : , 3 , 1) = efg_ia_real( : , 1 , 3) 
  efg_ia_real( : , 3 , 2) = efg_ia_real( : , 2 , 3)

  efg_ia_dual( : , 2 , 1) = efg_ia_dual( : , 1 , 2) 
  efg_ia_dual( : , 3 , 1) = efg_ia_dual( : , 1 , 3) 
  efg_ia_dual( : , 3 , 2) = efg_ia_dual( : , 2 , 3)

  ! =======
  ! 4pi/3V
  ! =======
  efg_ia_dual_real( : , : , :) = dble( efg_ia_dual( : , : , :)  ) 
  efg_ia_dual_real =  efg_ia_dual_real * fpi / omega !/ 3.0d0
  ! ==============
  !  total tensor
  ! ==============
  efg_ia = efg_ia_dual_real + efg_ia_real 
#ifdef debug
CALL print_tensor( efg_ia_dual_real ( 1 , : , : ) , 'DUALRE' )
CALL print_tensor( efg_ia_real      ( 1 , : , : ) , 'REAL  ' )
CALL print_tensor( efg_ia           ( 1 , : , : ) , 'TOTAL ' )
#endif
if ( lefg_it_contrib ) then
  efg_ia_it = efg_ia_dual_real_it + efg_ia_real_it 
endif
!=========================================================
!      END OF EFG TENSOR CALCULATION
!=========================================================

  ttt5 = MPI_WTIME(ierr)
  efgtimetot3 = efgtimetot3 + (ttt5-ttt4)

  CALL MPI_BARRIER( MPI_COMM_WORLD , ierr )

  return

END SUBROUTINE efg_ES_dip


!*********************** SUBROUTINE efg_write_output ********************************
!
! write distributions of EFG parameters and components to files
!
!******************************************************************************

SUBROUTINE efg_write_output ( kunit_eta , kunit_vzz , kunit_u )

  USE io_file,          ONLY :  ionode 
  USE config,           ONLY :  natm, natmi, ntype 
  USE constants,        ONLY :  dzero
  USE time

  implicit none
 
  INCLUDE 'mpif.h'

  ! global
  integer, intent(in)  :: kunit_eta , kunit_vzz , kunit_u 

  ! local
  integer :: i , it 

  if ( ionode ) then
   ! ======================== 
   !  write eta distributions
   ! ======================== 
   WRITE (kunit_eta,'(a,4f15.8)') '#',reseta,dble(natmi(1)),dble(natmi(2)),dble(ncefg)  
   if ( ntype .eq. 1 ) WRITE (kunit_eta,'(4f15.8)') dzero,dzero,dzero
   if ( ntype .eq. 2 ) WRITE (kunit_eta,'(4f15.8)') dzero,dzero,dzero,dzero
   do i = 0 , PANeta-1 
     WRITE (kunit_eta,'(4f15.8)') &
     dble((i+1-0.5D0) * reseta ), &
   ( dble(dibetatot(it,i)) / (reseta * dble(natmi(it)) * dble(ncefg)), it = 0 , ntype )  
   enddo
   ! ========================
   !  write Vzz distribution
   ! ========================
   do i = 0 , PANvzz 
     WRITE (kunit_vzz ,'(4f15.8)') &
     vzzmin + dble(i * resvzz), &
   ( dble(dibvzztot(it,i))/(resvzz * dble(natmi(it)) * dble(ncefg)), it = 0 ,ntype ) 
   enddo

   ! =================================================
   ! write U1 and average Uk (with k>1) distribution 
   ! =================================================
   do i = 0 , PANU 
   WRITE (kunit_u ,'(7f15.8)') umin + dble(i * resu), &
                              ( dble(dibUtot(1,it,i)) / ( resu * dble(natmi(it)) * dble(ncefg)), & 
                                dble(dibUtot(6,it,i)) / ( resu * 4.0d0*dble(natmi(it)) * dble(ncefg)), &
                                it = 0 ,ntype )
   enddo
   WRITE ( kunit_eta , '(a)' ) ''
   WRITE ( kunit_eta , '(a)' ) ''
   WRITE ( kunit_vzz , '(a)' ) ''
   WRITE ( kunit_vzz , '(a)' ) ''
   WRITE ( kunit_u   , '(a)' ) ''
   WRITE ( kunit_u   , '(a)' ) ''

  endif

  return

END SUBROUTINE efg_write_output

!*********************** SUBROUTINE efg_acf ***********************************
!
! calculates the auto-correlation function of the efg principal component 
! based on Allen-Tieldsley
!
!******************************************************************************

SUBROUTINE efg_acf

  USE config,           ONLY  : system , natm , ntype , itype , atype, omega ,box ,rho, config_alloc
  USE io_file,          ONLY  : ionode , stdout , kunit_EFGALL , kunit_EFGACFFF  , kunit_OUTFF

  implicit none

  integer :: ia , it , t , tt0 , t0 , t0max

  integer         , dimension (:,:) , allocatable :: norm
  double precision, dimension (:,:) , allocatable :: acfxx , acfyy , acfzz , acfet ! autocorelation function 
  double precision, dimension (:)   , allocatable :: timeo
  double precision, dimension (:)   , allocatable :: vxx0 , vyy0 , vzz0 , eta0     ! stored value at t=0
  double precision, dimension (:,:) , allocatable :: vxxt , vyyt , vzzt , etat     ! value at any t
  !trash
  integer :: iiii
  double precision :: aaaa      
  character*200 :: XXXX

  OPEN(kunit_EFGALL,FILE='EFGALL')
  READ(kunit_EFGALL,*) natm,box,ntype
  READ(kunit_EFGALL,*) system,iiii
  READ(kunit_EFGALL,*) XXXX
  omega = box * box * box
  rho = natm / omega
  CALL config_alloc 

  CALL print_general_info ( stdout ) 
  CALL print_general_info ( kunit_OUTFF ) 

  allocate ( vxxt(natm,ncefg) , vyyt(natm,ncefg)      , vzzt(natm,ncefg)        , etat(natm,ncefg)                                                            )
  allocate ( vxx0(natm)       , vyy0(natm)            , vzz0(natm)              , eta0(natm)                                                                  )
  allocate ( timeo (ncefg)    , norm(0:ntype,0:ntcor) , &
             acfxx (0:ntype,0:ntcor) , acfyy (0:ntype,0:ntcor) , &
             acfzz (0:ntype,0:ntcor) , acfet (0:ntype,0:ntcor) )
 
  acfxx = 0.0d0
  acfyy = 0.0d0
  acfzz = 0.0d0
  acfet = 0.0d0
  
  do t0 = 1 , ncefg 
    if ( t0 .ne. 1 ) then
      READ(kunit_EFGALL,*) natm,box,ntype
      READ(kunit_EFGALL,*) system,iiii
      READ(kunit_EFGALL,*) XXXX
    endif
    do ia = 1 , natm
      READ(kunit_EFGALL,*) iiii , atype(ia) , vxxt( ia , t0 ) , vyyt( ia , t0 ) &
                                            , vzzt( ia , t0 ) , etat( ia , t0 ) , aaaa , aaaa , aaaa
    enddo
  enddo
  CLOSE(kunit_EFGALL)
 
  if ( ionode ) WRITE ( stdout , '(a)' ) 'EFGALL successfully readed'

  do t0 = 1 , ncefg 
    do ia = 1 ,natm
      vxx0( ia ) = vxxt( ia , t0 )
      vyy0( ia ) = vyyt( ia , t0 )
      vzz0( ia ) = vzzt( ia , t0 )
      eta0( ia ) = etat( ia , t0 )
    enddo
    t0max = min ( ncefg , t0 +ntcor )
    do tt0 = t0 , t0max
      t = tt0 - t0
      do ia = 1 , natm
        if ( atype (ia) .eq. 'A' ) then
          acfxx ( 1, t ) = acfxx ( 1, t ) + vxx0( ia ) * vxxt( ia , tt0 )
          acfyy ( 1, t ) = acfyy ( 1, t ) + vyy0( ia ) * vyyt( ia , tt0 )
          acfzz ( 1, t ) = acfzz ( 1, t ) + vzz0( ia ) * vzzt( ia , tt0 )
          acfet ( 1, t ) = acfet ( 1, t ) + eta0( ia ) * etat( ia , tt0 )
          norm  ( 1, t ) = norm  ( 1, t ) + 1    
        elseif ( atype (ia) .eq. 'B' ) then
          acfxx ( 2, t ) = acfxx ( 2, t ) + vxx0( ia ) * vxxt( ia , tt0 )
          acfyy ( 2, t ) = acfyy ( 2, t ) + vyy0( ia ) * vyyt( ia , tt0 )
          acfzz ( 2, t ) = acfzz ( 2, t ) + vzz0( ia ) * vzzt( ia , tt0 )
          acfet ( 2, t ) = acfet ( 2, t ) + eta0( ia ) * etat( ia , tt0 )
          norm  ( 2, t ) = norm  ( 2, t ) + 1    
        endif
          acfxx ( 0, t ) = acfxx ( 0, t ) + vxx0( ia ) * vxxt( ia , tt0 )
          acfyy ( 0, t ) = acfyy ( 0, t ) + vyy0( ia ) * vyyt( ia , tt0 )
          acfzz ( 0, t ) = acfzz ( 0, t ) + vzz0( ia ) * vzzt( ia , tt0 )
          acfet ( 0, t ) = acfet ( 0, t ) + eta0( ia ) * etat( ia , tt0 )
          norm  ( 0, t ) = norm  ( 0, t ) + 1
      enddo
    enddo      
  enddo

  OPEN(UNIT=kunit_EFGACFFF,FILE='EFGACFFF')
  do t = 0 ,ntcor
    do it = 0 , ntype
      acfxx( it , t ) = acfxx( it , t ) / norm( it , t )
      acfyy( it , t ) = acfyy( it , t ) / norm( it , t )
      acfzz( it , t ) = acfzz( it , t ) / norm( it , t )
      acfet( it , t ) = acfet( it , t ) / norm( it , t )
    enddo
    if ( ionode ) then
    if ( ntype .eq. 2 ) &
    WRITE ( kunit_EFGACFFF , '(i12,12f20.10)' ) t , acfxx(1,t) , acfyy(1,t) , acfzz(1,t) , acfet(1,t) , &
                                                    acfxx(2,t) , acfyy(2,t) , acfzz(2,t) , acfet(2,t) , &
                                                    acfxx(0,t) , acfyy(0,t) , acfzz(0,t) , acfet(0,t)
    if ( ntype .eq. 1 ) &
    WRITE ( kunit_EFGACFFF , '(i12,8f20.10)' )  t , acfxx(1,t) , acfyy(1,t) , acfzz(1,t) , acfet(1,t)
    endif
  enddo

  CLOSE(kunit_EFGACFFF)      

  deallocate(norm,acfxx,acfyy,acfzz,acfet)
  deallocate(vxxt,vyyt,vzzt,etat)
  deallocate(vxx0,vyy0,vzz0,eta0)

  return
        
END SUBROUTINE efg_acf

!*********************** SUBROUTINE efg_stat **********************************
!
!
!******************************************************************************

SUBROUTINE efg_stat ( kunit_input , kunit_output , kunit_nmroutput )

  USE config,           ONLY  : system , natm , natmi , ntype , itype , &
                                atype , atypei , omega , box , rho, config_alloc
  USE io_file,          ONLY  : ionode , stdout , kunit_OUTFF 
  USE prop,              ONLY : nprop_print
  USE constants

  implicit none

  integer            :: i , ic , ia , it , ui
  integer, parameter :: lwork = 6
  integer            :: ifail
  integer            :: kvzz, keta ,ku
  integer            :: kunit_input , kunit_output , kunit_nmroutput
  integer            :: cc , ccs

  double precision                                 :: sq3 , sq32
  double precision, dimension (:,:), allocatable   :: U
  double precision                                 :: w(3) , efgt(3,3) , nmr_conv ( 4 )
  double precision                                 :: work(3 * lwork)
  double precision, dimension (:,:)  , allocatable :: nmr
  double precision, dimension (:)    , allocatable :: vzzmini , vzzmaxi ! ntype
  double precision, dimension (:)    , allocatable :: vzzm, vzzma , etam , pvzz, rho_z, sigmavzz, vzzsq ! ntype
  double precision :: vzzk , etak , uk

  !trash
  integer :: iiii
  character :: XXXX


  ! some constants
  sq3 = dsqrt( 3.0d0 )
  sq3 = 1.0d0 / sq3
  sq32 = sq3 * 0.5d0

  READ ( kunit_input , * )  natm
  READ ( kunit_input , * )  system
  READ ( kunit_input , * )  box,ntype
  READ ( kunit_input ,* ) ( atypei ( it ) , it = 1 , ntype )
  IF ( ionode ) WRITE ( kunit_OUTFF ,'(A,20A3)' ) 'found type information on EFGALL : ', atypei ( 1:ntype )
  IF ( ionode ) WRITE ( stdout      ,'(A,20A3)' ) 'found type information on EFGALL : ', atypei ( 1:ntype )
  READ( kunit_input ,*)   ( natmi ( it ) , it = 1 , ntype )
  READ ( kunit_input , * )  xxxx
  omega = box * box * box
  rho = natm / omega
  ! ===================================
  !  here we know natm, then alloc 
  !  and decomposition can be applied 
  ! ================================== 
  if ( lefg_restart ) CALL config_alloc 
  if ( lefg_restart ) CALL efg_alloc

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
    enddo
  enddo
  natmi  ( 0 ) = natm
  atypei ( 0 ) = 'ALL'
  !  ================
  !   allocation
  !  ================
  allocate ( U        ( natm , 5 ) )                       ! Czjzek U component
  allocate ( nmr      ( natm , 4 ) )                       ! 1 = Vxx ; 2 = Vyy ; 3 = Vzz ; 4 = eta
  allocate ( vzzmini  ( 0:ntype  ), vzzmaxi ( 0:ntype ) )  ! max min value ! not so important
  allocate ( vzzm     ( 0:ntype  ) )                       ! vzz mean value (per type)
  allocate ( vzzma    ( 0:ntype  ) )                       ! mean value of absolute vzz (per type)
  allocate ( vzzsq    ( 0:ntype  ) )                       ! square of vzz (per type)
  allocate ( etam     ( 0:ntype  ) )                       ! eta mean value (per type)
  allocate ( pvzz     ( 0:ntype  ) )                       ! proprtion of positiv vzz (per type)
  allocate ( rho_z    ( 0:ntype  ) )                       ! rho_z (per type)
  allocate ( sigmavzz ( 0:ntype  ) )                       ! variance of vzz distribution (per type)
  ! ==============
  ! set to zero
  ! ==============
  nmr      = 0.0D0 
  pvzz     = 0.0d0
  vzzmini  = 0.0D0
  vzzmaxi  = 0.0D0
  etam     = 0.0d0
  vzzm     = 0.0d0
  vzzma    = 0.0d0
  vzzsq    = 0.0d0
  U        = 0.0D0

  do ic = 1 , ncefg
    if ( ic .ne. 1 ) then
      READ ( kunit_input , * )  natm
      READ ( kunit_input , * )  system
      READ ( kunit_input , * )  box,ntype
      READ ( kunit_input , * )  ( xxxx , it = 1 , ntype )
      READ ( kunit_input , * )  ( iiii  , it = 1 , ntype )
      READ ( kunit_input , * )  xxxx
    endif

    do ia = 1 , natm

      READ( kunit_input , * ) iiii , xxxx , efgt(1,1) , efgt(2,2) , efgt(3,3) , &
                                            efgt(1,2) , efgt(1,3) , efgt(2,3)

      ! ===================================================================== 
      !  Czjzek components (see J. Phys.: Condens. Matter 10 (1998). p10719)
      ! =====================================================================
      U ( ia , 1 ) = efgt ( 3 , 3 ) * 0.5d0
      U ( ia , 2 ) = efgt ( 1 , 3 ) * sq3
      U ( ia , 3 ) = efgt ( 2 , 3 ) * sq3
      U ( ia , 4 ) = efgt ( 1 , 2 ) * sq3
      U ( ia , 5 ) = ( efgt ( 1 , 1 ) - efgt ( 2 , 2 ) ) * sq32
  
      ! =================
      !  diagonalisation
      ! =================
      CALL DSYEV ( 'N' , 'U' , 3 , efgt , 3 , w , work , 3 * lwork , ifail )
      if ( ifail .ne. 0 ) then
        if ( ionode ) WRITE ( stdout , * ) 'ERROR: DSYEV, STOP in efg MODULE (ewald)'
        STOP
      endif

      ! =============================
      ! NMR conventions test 
      ! =============================
      CALL nmr_convention( w , nmr_conv ( : )  , ia )
      nmr( ia , : )  = nmr_conv ( : )
    
      ! ===================================================
      !  statistic for all sites index 0 in array of size ntype
      ! ===================================================
      vzzm  ( 0 )  =  vzzm ( 0 ) +       nmr ( ia , 3 )
      etam  ( 0 )  =  etam ( 0 ) +       nmr ( ia , 4 )
      vzzma ( 0 )  =  vzzma( 0 ) + dabs( nmr ( ia , 3 ) )
      if ( nmr ( ia , 3 ) .le. vzzmini ( 0 ) )  vzzmini ( 0 ) = nmr ( ia , 3 )
      if ( nmr ( ia , 3 ) .ge. vzzmaxi ( 0 ) )  vzzmaxi ( 0 ) = nmr ( ia , 3 )
      if ( nmr ( ia , 3 ) .ge. 0.0d0 ) pvzz ( 0 ) = pvzz ( 0 ) + 1.d0
      vzzsq ( 0 ) = vzzsq ( 0 ) + nmr ( ia , 3 ) * nmr ( ia , 3 )

      ! ===================================================
      !  statistic type specific 
      ! ===================================================
      do it = 1 , ntype
        if (itype(ia).eq.it) then
          vzzm ( it ) = vzzm ( it )  + nmr ( ia , 3 )
          vzzma( it ) = vzzma( it )  + dabs( nmr (ia , 3 ) )
          etam ( it ) = etam ( it )  + nmr ( ia , 4 )
          if (nmr(ia,3).le.vzzmini(it))  vzzmini(it) = nmr(ia,3)
          if (nmr(ia,3).ge.vzzmaxi(it))  vzzmaxi(it) = nmr(ia,3)
          vzzsq(it) = vzzsq(it) + nmr(ia,3) * nmr(ia,3)
          if (nmr(ia,3).ge.0.0d0 ) pvzz(it) = pvzz(it) + 1.d0
        endif ! it 
      enddo
    enddo
    ! =================================
    !  END OF STAT FOR A GIVEN CONFIG
    ! =================================

    if (ntype.eq.1) then
      vzzm  = vzzm  / dble( natm )
      vzzsq = vzzsq / dble( natm )
      vzzma = vzzma / dble( natm )
      etam  = etam  / dble( natm )
      pvzz  = pvzz  / dble( natm )
    else
      vzzm  ( 0 )  = vzzm  ( 0 ) / dble( natm )
      vzzsq ( 0 )  = vzzsq ( 0 ) / dble( natm )
      vzzma ( 0 )  = vzzma ( 0 ) / dble( natm )
      etam  ( 0 )  = etam  ( 0 ) / dble( natm )
      pvzz  ( 0 )  = pvzz  ( 0 ) / dble( natm )
      do it=1,ntype
        vzzm ( it ) = vzzm ( it ) / dble( natmi ( it ) )
        vzzsq( it ) = vzzsq( it ) / dble( natmi ( it ) )
        vzzma( it ) = vzzma( it ) / dble( natmi ( it ) )
        etam ( it ) = etam ( it ) / dble( natmi ( it ) )
        pvzz ( it ) = pvzz ( it ) / dble( natmi ( it ) )
      enddo
    endif

    do it =0 , ntype
      sigmavzz ( it ) = vzzsq ( it ) - ( vzzma ( it ) * vzzma ( it ) )
      sigmavzz ( it ) = dsqrt( sigmavzz ( it ) )
      rho_z( it ) = sigmavzz ( it ) / vzzma ( it )
    enddo

    ! ========================================
    !  output average values ( instantaneous )
    ! ========================================
    do it=1,ntype 
      if ( ( ionode ) .and. (mod(ic,nprop_print) .eq. 0) ) &
      WRITE ( stdout ,100) & 
      ic , atypei(it) , vzzmini(it) , vzzmaxi(it) , vzzm(it) , vzzma(it) , etam(it) , pvzz(it) , rho_z ( it)
      if (   ionode ) WRITE ( kunit_output , 100) &
      ic , atypei(it) , vzzmini(it) , vzzmaxi(it) , vzzm(it) , vzzma(it) , etam(it) , pvzz(it) , rho_z(it)
    enddo
    if ( ionode .and. ntype .ne. 1 .and. (mod(ic,nprop_print) .eq. 0) ) &
      WRITE ( stdout ,100) &
      ic , atypei(0) , vzzmini(0) , vzzmaxi(0) , vzzm(0) , vzzma(0) , etam(0) , pvzz(0) , rho_z(0)
    if ( ionode .and. ntype .ne. 1 ) &
      WRITE ( kunit_output ,100) &
      ic , atypei(0) , vzzmini(0) , vzzmaxi(0) , vzzm(0) , vzzma(0) , etam(0) , pvzz(0) , rho_z(0)

    ! ========================================
    !  output NMR  parameters 
    ! ========================================
    if ( ionode ) then
        WRITE ( kunit_nmroutput , '(a)' ) &
        '#     site           xx          yy          zz               eta'
      do ia = 1 , natm
        WRITE ( kunit_nmroutput , 150 ) ia , atype ( ia ) , & 
                                        nmr( ia , 1 ) , nmr( ia , 2 ) , &
                                        nmr( ia , 3 ) , nmr( ia , 4 )
      enddo
    endif


    ! =======================================================
    !  DISTRIBUTION CALCULATION LOOP 
    ! =======================================================
    do ia = 1 , natm
      ! ========
      !    Ui
      ! ========
      do ui=1,5
        uk = (U(ia,ui)-umin)/resu
        ku = int(uk) + 1
        ! ====================== 
        !  test out of bound
        ! ====================== 
        if (ku.lt.0.or.ku.gt.PANU) then
          if ( ionode ) WRITE ( stdout , * ) 'ERROR: out of bound dibU1'
          if ( ionode ) WRITE ( stdout ,310) i,ku,U(i,ui),umin,dabs(umin)
          STOP
        endif
        ! all types
        dibUtot(ui,0,ku) = dibUtot(ui,0,ku) + 1
        ! type specific
        do it=1,ntype
          if (itype(ia).eq.it) then
            dibUtot(ui,it,ku) = dibUtot(ui,it,ku) + 1
          endif
        enddo
        ! average Uk,k>1
        if ( ui .ne. 1 ) then
          dibUtot(6,0,ku) = dibUtot(6,0,ku) + 1
          do it=1,2
            if (itype(ia).eq.it) then
              dibUtot(6,it,ku) = dibUtot(6,it,ku) + 1
            endif
          enddo
        endif
       enddo
      ! ======================== 
      !  quadrupolar parameters
      ! ======================== 
      vzzk = (nmr(ia,3)-vzzmin) / resvzz
      etak = nmr(ia,4) / reseta
      kvzz = int(vzzk)
      keta = int(etak)
      ! ====================== 
      !  test out of bound
      ! ====================== 
      if ( kvzz .lt. 0 .or. kvzz .gt. PANvzz ) then
        if ( ionode .and. itype(ia) .eq. 1 ) &
        WRITE ( stdout , * ) 'ERROR: out of bound distribvzz A'
        if ( ionode .and. itype(ia) .eq. 2 ) &
        WRITE ( stdout , * ) 'ERROR: out of bound distribvzz B'
        if ( ionode ) WRITE ( stdout ,200) &
        ia , kvzz , nmr ( ia , 3 ) , vzzmini( itype ( ia ) ) , vzzmaxi ( itype ( ia ) ) , &
                    nmr ( ia , 4 ) , nmr ( ia , 1 ) , nmr ( ia , 2 ) , nmr ( ia , 3 )
        STOP
      endif
      if ( keta .lt. 0 .or. keta .gt. PANeta ) then
        if ( ionode ) WRITE ( stdout , * ) 'ERROR: out of bound distribeta'
        if ( ionode ) WRITE ( stdout ,210) ia, keta, nmr(ia,4), nmr(ia,4), &
                                                     nmr(ia,1), nmr(ia,2), &
                                                     nmr(ia,3)
        STOP
      endif
      ! distribution per type 
      dibvzztot(0,kvzz) = dibvzztot(0,kvzz) + 1
      dibetatot(0,keta) = dibetatot(0,keta) + 1
      ! type specific 
      do it=1,ntype
        if ( itype(ia) .eq. it ) then
          dibvzztot(it,kvzz) = dibvzztot(it,kvzz) + 1
          dibetatot(it,keta) = dibetatot(it,keta) + 1
        endif
      enddo
    enddo  !ia 
    !==========================
    ! END of distributions loop
    !==========================

  enddo ! config loop

  lefg_restart = .false. ! only allocate for the first EFGALL file

  !  ================
  !   deallocation
  !  ================
  deallocate ( nmr               )   
  deallocate ( vzzmini , vzzmaxi )
  deallocate ( U                 )
  deallocate ( vzzm              )
  deallocate ( vzzma             )
  deallocate ( etam              )
  deallocate ( pvzz              )
  deallocate ( rho_z             )
  deallocate ( sigmavzz          )
  deallocate ( vzzsq             )

100 FORMAT(I7,1X,' ATOM ',A3,'  minVZZ = ',F7.3,' maxVZZ = ',F7.3,& 
             ' <VZZ> =  ',F10.5,' <|VZZ|> =  ',F10.5,' mETA =  ',F10.5, &
          ' P(Vzz>0) =  ',F10.5, ' RHO_Z = ',F10.5)

150 FORMAT(I7,1X,A3,3F12.4,F18.4)

200 FORMAT('atom = ',I6,' k = ',I12, 'value = ',F14.8,&
           ' min =  ',F7.3,' max =  ',F7.3,' vaa{a = xx,yy,zz}, eta = ',4F14.8)
210 FORMAT('atom = ',I6,' k = ',I12, 'value = ',F14.8,& 
           ' min =  0.0D0    max =  1.0D0',' vaa{a = xx,yy,zz}, eta = ',4F14.8)

310 FORMAT('atom = ',I6,' k = ',I12, 'value = ',F14.8,' min =  ',F7.3,' max =  ',F7.3)

  return

END SUBROUTINE efg_stat

!*********************** SUBROUTINE nmr_convention ****************************
!
!
!******************************************************************************

SUBROUTINE nmr_convention( w , nmr , ia )

  USE io_file,  ONLY :  stdout , ionode

  implicit none

  double precision ::  W(3) , tmpw
  double precision ::  nmr(4)
  integer :: k , ia
  
 ! =======================================
 !  NMR convention: |Vzz| > |Vyy| > |Vxx|
 ! =======================================
  DO K = 1,2
    IF (abs(W(1))>abs(W(2))) THEN
      TMPW = W(2)
      W(2) = W(1)
      W(1) = TMPW
    ENDIF
    IF (abs(W(2))>abs(W(3))) THEN
      TMPW = W(3)
      W(3) = W(2)
      W(2) = TMPW
    ENDIF
  ENDDO

  nmr(3) = W(3) !ZZ
  nmr(2) = W(1) !YY
  nmr(1) = W(2) !XX
  nmr(4) = (nmr(2)-nmr(1))/nmr(3) !ETA

  ! ===================
  !  test rearangement
  ! ===================
  if ( nmr(1) .eq. nmr(2) )  nmr(4) = 0.0D0
  if ( nmr(3) .eq. 0.0 )  nmr(4) = 0.0D0
  if ( nmr(1) .lt. 1.0E-08 .and. nmr(2) .lt. 1.0E-8 .and. nmr(3) .lt. 1.0E-8 ) nmr(4) = 0.0D0
  if (dabs(nmr(3)).lt.dabs(nmr(2))) then
    if ( ionode ) WRITE ( stdout , * ) 'ERROR: |Vzz| < |Vxx|',ia,nmr(1),nmr(2),nmr(3)
    STOP
  endif
  if (dabs(nmr(3)).lt.dabs(nmr(1))) then
    if ( ionode ) WRITE ( stdout , * ) 'ERROR: |Vzz| < |Vyy|',ia,nmr(1),nmr(2),nmr(3)
    STOP
  endif
  if (dabs(nmr(1)).lt.dabs(nmr(2))) then
    if ( ionode ) WRITE ( stdout , * ) 'ERROR: |Vyy| < |Vxx|',ia,nmr(1),nmr(2),nmr(3)
    STOP
  endif
  if ( nmr(4) .gt. 1.0d0 .or. nmr(4) .lt. 0.0d0) then
    if ( ionode ) &
    WRITE ( stdout ,'(a,4f48.24)') 'ERROR: eta > 1.0d0 or eta < 0.0d0',ia,nmr(4),nmr(1),nmr(2),nmr(3)
    STOP
  endif

  return

END SUBROUTINE nmr_convention

END MODULE EFG
! ===== fmV =====
