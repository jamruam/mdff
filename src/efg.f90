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
  double precision, dimension(:,:,:)  , allocatable :: efg_ia       ! efg_tensor
  double precision, dimension(:,:,:,:), allocatable :: efg_ia_it    ! efg_tensor (it contribution)
  TYPE ( rmesh ) :: rm_efg

  ! ==================
  !  ewald summation
  ! ==================
  integer                                           :: ncellewald
  double precision                                  :: alphaES          ! Ewald sum parameter 
  double precision, dimension(:,:,:)  , allocatable :: efg_sum3         ! efg_tensor
  double precision, dimension(:,:,:)  , allocatable :: efg_ia_real         ! efg_tensor
  double precision, dimension(:,:,:,:), allocatable :: efg_ia_real_it      ! efg_tensor
  double complex  , dimension(:,:,:)  , allocatable :: efg_ia_dual         ! efg_tensor
  double complex  , dimension(:,:,:,:), allocatable :: efg_ia_dual_it      ! efg_tensor
  double precision, dimension(:,:,:)  , allocatable :: efg_ia_dual_real    ! efg_tensor
  double precision, dimension(:,:,:,:), allocatable :: efg_ia_dual_real_it ! efg_tensor

  TYPE ( kmesh ) :: km_efg

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

  USE control,  ONLY :  calc
  USE prop,     ONLY :  lefg
  USE io_file,  ONLY :  stdin, stdout, kunit_OUTFF, ionode
  USE field,    ONLY :  initialize_coulomb
 
  implicit none
  
  ! local
  integer :: ioerr
  character * 132 :: filename

  namelist /efgtag/  lefgprintall    , & 
                     lefg_it_contrib , &
                     ncefg           , & 
                     ntcor           , &
                     ncellewald      , &
                     ncelldirect     , &
                     cutefg          , &
                     alphaES         , &
                     resvzz          , &
                     reseta          , &
                     resu            , &
                     vzzmin          , &
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
  lefgprintall  = .false. 
  lefg_it_contrib = .false.
  reseta        =   0.1D0
  resvzz        =   0.1D0
  resu          =   0.1D0 
  ncelldirect   =   0
  ncellewald    =   0 
  ncefg         =   0
  cutefg        =  30.0D0
  alphaES       =   1.0D0
  umin          =  -4.0d0
  vzzmin        =  -4.0d0
  ntcor         =  10

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
        WRITE ( kunit ,'(a)')         'ewald summation'
        WRITE ( kunit ,'(a,f10.5)')   'alpha                            = ',alphaES
        WRITE ( kunit ,'(a,f10.5)')   'cutoff in real part              = ',cutefg
        WRITE ( kunit ,'(a,i10)')     'ncellewald x ncellewald x ncellewald               = ',ncellewald
        WRITE ( kunit ,'(a,i10)')     '' 
        WRITE ( kunit ,'(a,f10.5)')   'Note:this should hold alpha^2 * box^2 >> 1',alphaES*alphaES*box*box
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

  ! ============
  !  direct sum
  ! ============
  if ( longrange .eq. 'direct' ) then
    ncmax = ( 2 * ncelldirect + 1 ) ** 3
    rm_efg%ncmax=ncmax
    rm_efg%ncell=ncelldirect
    allocate( rm_efg%boxxyz( 3 , ncmax ) , rm_efg%lcell( ncmax ) )
    CALL direct_sum_init ( rm_efg )
    allocate( efg_ia   ( natm , 3 , 3 ) )
    allocate( efg_sum3 ( natm , 3 , 3 ) )
    efg_ia    = 0.0d0 
    if ( lefg_it_contrib ) then
      allocate( efg_ia_it( natm , ntype , 3 , 3 ) )
    efg_ia_it = 0.0d0 
    endif
  endif

  ! ============
  !  ewald sum
  ! ============
  if ( longrange .eq. 'ewald')  then

    allocate( efg_ia_real     ( natm , 3 , 3 ) )
    allocate( efg_sum3        ( natm , 3 , 3 ) )
    allocate( efg_ia_dual     ( natm , 3 , 3 ) ) 
    allocate( efg_ia_dual_real( natm , 3 , 3 ) ) 
    allocate( efg_ia    ( natm , 3 , 3 ) )
    efg_ia_real = 0.0d0 ; efg_ia_dual = 0.0d0 ; efg_ia_dual_real =0.0d0 ; efg_ia = 0.0d0  

    if ( lefg_it_contrib ) then
      allocate( efg_ia_real_it     ( natm , ntype , 3 , 3 ) )
      allocate( efg_ia_dual_it     ( natm , ntype , 3 , 3 ) ) 
      allocate( efg_ia_dual_real_it( natm , ntype , 3 , 3 ) ) 
      allocate( efg_ia_it    ( natm , ntype , 3 , 3 ) )
    efg_ia_real_it = 0.0d0 ; efg_ia_dual_it = 0.0d0 ; efg_ia_dual_real_it =0.0d0 ; efg_ia_it = 0.0d0  
    endif

    km_efg%ncell = ncellewald
    nkcut = ( 2 * ncellewald + 1 ) ** 3 
    nkcut = nkcut - 1
    km_efg%nkcut = nkcut
    allocate( km_efg%kptk( nkcut ) , km_efg%kpt(3,nkcut) )
    allocate ( km_efg%strf ( nkcut, ntype) ) 
    CALL kpoint_sum_init ( km_efg )
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

  if ( longrange .eq. 'direct' ) then 
    deallocate( rm_efg%boxxyz , rm_efg%lcell )
    deallocate( efg_ia )  
    deallocate( efg_sum3  )
    if ( lefg_it_contrib ) then
      deallocate( efg_ia_it )  
    endif
  endif

  if ( longrange .eq. 'ewald' )  then

    deallocate( efg_sum3  )
    deallocate( efg_ia_real )
    deallocate( efg_ia_dual ) 
    deallocate( efg_ia_dual_real ) 
    deallocate( efg_ia )

    if ( lefg_it_contrib ) then
      deallocate( efg_ia_real_it )
      deallocate( efg_ia_dual_it ) 
      deallocate( efg_ia_dual_real_it ) 
      deallocate( efg_ia_it )
    endif

    deallocate( km_efg%kptk , km_efg%kpt )
    deallocate( km_efg%strf ) 

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

  USE io_file,  ONLY :  ionode , stdout , kunit_TRAJFF, kunit_EFGFF, kunit_EFGALL, & 
                        kunit_OUTFF , kunit_tmp , kunit_EFGALLIT1 , kunit_EFGALLIT2 , kunit_EFGFFIT1, kunit_EFGFFIT2 ,&
                        kunit_DTETAFF , kunit_DTETAFFIT , kunit_DTVZZFF , kunit_DTVZZFFIT , kunit_DTIBUFF , kunit_DTIBUFFIT

  USE config,   ONLY :  system , natm , ntype , atype , rx , ry , rz , itype , & 
                        atypei , natmi, rho , box , omega , config_alloc , qia 
  USE control,  ONLY :  longrange , myrank , numprocs
  USE field,    ONLY :  qch , field_init

  implicit none

  ! local
  integer :: iastart , iaend
  integer :: ia, iconf , na , it , cc , ccs
  double precision :: aaaa 
  integer :: iiii
  character * 60 :: cccc 

#ifdef fix_grid
  double precision, dimension ( : , : ) , allocatable :: rave !average positions
#endif

  OPEN (unit = kunit_EFGALL  ,file = 'EFGALL')

  if ( lefg_it_contrib ) then
    OPEN (unit = kunit_EFGALLIT1  ,file = 'EFGALLIT1')
    OPEN (unit = kunit_EFGALLIT2  ,file = 'EFGALLIT2')
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

#ifndef fix_grid
    if ( ionode ) WRITE ( stdout , '(a)' ) 'not fix_grid' 
! =============================================
!  not fix_grid: we calculate efg at 
!  each atom positions which are moving.
! =============================================

  ! ============================================
  ! LOOP OVER CONFIGURATIONS 
  ! ============================================
  do iconf = 1, ncefg
    if ( ionode ) WRITE( stdout ,'(a,i6)') 'efg calculation of config',iconf
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
      enddo
    enddo
    natmi  ( 0 ) = natm
    atypei ( 0 ) = 'ALL'

    ! ===============
    ! use direct sum
    ! ===============
    if ( longrange .eq. 'direct' )  CALL efg_DS( iconf , iastart , iaend )
    ! ===============
    ! use ewald sum
    ! ===============
    if ( longrange .eq. 'ewald' )   CALL efg_ES( iconf , iastart , iaend )

  enddo ! iconf loop

  CLOSE(kunit_EFGALL)
  if ( lefg_it_contrib ) then
    CLOSE(kunit_EFGALLIT1)
    CLOSE(kunit_EFGALLIT2)
  endif

  !calculate statistical properties from EFGALL
  OPEN ( unit = kunit_EFGALL  ,file = 'EFGALL')
  OPEN ( unit = kunit_EFGFF  ,file = 'EFGFF')
  CALL efg_stat( kunit_EFGALL , kunit_EFGFF )
  CLOSE( kunit_EFGALL )
  CLOSE( kunit_EFGFF )
  ! write average distribution output from EFGALL 
  OPEN ( unit = kunit_DTETAFF,file = 'DTETAFF')
  OPEN ( unit = kunit_DTVZZFF,file = 'DTVZZFF')
  OPEN ( unit = kunit_DTIBUFF,file = 'DTIBUFF')
  CALL efg_write_output( kunit_DTETAFF , kunit_DTVZZFF , kunit_DTIBUFF )
  CLOSE ( kunit_DTETAFF )
  CLOSE ( kunit_DTVZZFF )
  CLOSE ( kunit_DTIBUFF )
  dibUtot = 0
  dibvzztot = 0
  dibetatot = 0

  if ( lefg_it_contrib ) then
    ! calculate statistical properties from EFGALLIT
    OPEN ( unit = kunit_EFGALLIT1  ,file = 'EFGALLIT1')
    OPEN ( unit = kunit_EFGFFIT1   ,file = 'EFGFFIT1')
    CALL efg_stat ( kunit_EFGALLIT1 , kunit_EFGFFIT1 )
    CLOSE( kunit_EFGALLIT1 )
    CLOSE( kunit_EFGFFIT1 )
    ! write average distribution output from EFGALLIT1
    OPEN ( unit = kunit_DTETAFFIT ,file = 'DTETAFFIT')
    OPEN ( unit = kunit_DTVZZFFIT ,file = 'DTVZZFFIT')
    OPEN ( unit = kunit_DTIBUFFIT ,file = 'DTIBUFFIT')
    CALL efg_write_output( kunit_DTETAFFIT , kunit_DTVZZFFIT , kunit_DTIBUFFIT )
    dibUtot = 0
    dibvzztot = 0
    dibetatot = 0
    !calculate statistical properties from EFGALL
    OPEN ( unit = kunit_EFGALLIT2 ,file = 'EFGALLIT2')
    OPEN ( unit = kunit_EFGFFIT2  ,file = 'EFGFFIT2')
    CALL efg_stat( kunit_EFGALLIT2 , kunit_EFGFFIT2 )
    CLOSE( kunit_EFGALLIT2 )
    CLOSE( kunit_EFGFFIT2 )
    ! write average distribution output from EFGALLIT1
    WRITE ( kunit_DTETAFFIT,'(a)') ''
    WRITE ( kunit_DTETAFFIT,'(a)') ''
    WRITE ( kunit_DTVZZFFIT,'(a)') ''
    WRITE ( kunit_DTVZZFFIT,'(a)') ''
    WRITE ( kunit_DTIBUFFIT,'(a)') ''
    WRITE ( kunit_DTIBUFFIT,'(a)') ''
    CALL efg_write_output( kunit_DTETAFFIT , kunit_DTVZZFFIT , kunit_DTIBUFFIT )
    CLOSE ( kunit_DTETAFFIT )
    CLOSE ( kunit_DTVZZFFIT )
    CLOSE ( kunit_DTIBUFFIT )
  endif

  CLOSE(kunit_TRAJFF)


#else
!its probably wrong since the rewriting work of the efg module
  WRITE ( stdout , '(a)' ) 'fix_grid used' 

  allocate ( rave ( 3 , natm ) )
  rave = 0.0D0
! =================================================
!  fix_grid: we calculate efg at fixed positions
!  which are the average positions of the dynamic
! =================================================

! =================================================
!  so first we calculate the average configuration
! =================================================
  do iconf = 1, ncefg
    na = 0
    if ( iconf .ne. 1 ) READ ( kunit_TRAJFF, * )  iiii
    if ( iconf .ne. 1 ) READ ( kunit_TRAJFF, * )  cccc
    if ( iconf .ne. 1 ) READ ( kunit_TRAJFF, * )  aaaa , iiii
    do ia = 1 , natm
      READ ( kunit_TRAJFF, * ) atype ( ia ) , rx ( ia ) , ry ( ia ) , rz ( ia )
      if ( atype ( ia ) .eq. 'A' ) na = na + 1
      if ( atype ( ia ) .eq. 'A' ) itype ( ia ) = 1
      if ( atype ( ia ) .eq. 'B' ) itype ( ia ) = 2
      if ( iconf .eq. 1) then
        rave ( 1 , ia ) = rave ( 1 , ia ) + rx ( ia ) 
        rave ( 2 , ia ) = rave ( 2 , ia ) + ry ( ia ) 
        rave ( 3 , ia ) = rave ( 3 , ia ) + rz ( ia ) 
      else
        if ( rave ( 1 , ia ) .gt. 0 ) rave ( 1 , ia ) = rave ( 1 , ia ) + abs( rx ( ia ) )
        if ( rave ( 2 , ia ) .gt. 0 ) rave ( 2 , ia ) = rave ( 2 , ia ) + abs( ry ( ia ) )
        if ( rave ( 3 , ia ) .gt. 0 ) rave ( 3 , ia ) = rave ( 3 , ia ) + abs( rz ( ia ) )
        if ( rave ( 1 , ia ) .lt. 0 ) rave ( 1 , ia ) = rave ( 1 , ia ) - abs( rx ( ia ) )
        if ( rave ( 2 , ia ) .lt. 0 ) rave ( 2 , ia ) = rave ( 2 , ia ) - abs( ry ( ia ) )
        if ( rave ( 3 , ia ) .lt. 0 ) rave ( 3 , ia ) = rave ( 3 , ia ) - abs( rz ( ia ) )
      endif
    enddo
  enddo ! ifconf

    rave(1,:) = rave(1,:) / dble ( ncefg )
    rave(2,:) = rave(2,:) / dble ( ncefg )
    rave(3,:) = rave(3,:) / dble ( ncefg )

    rgrid = rave

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
    enddo
  enddo
  natmi  ( 0 ) = natm
  atypei ( 0 ) = 'ALL'

  ! write average configuration
  if ( ionode ) then
  OPEN ( kunit_tmp ,file = 'CAVERFF',STATUS = 'UNKNOWN')
      WRITE ( kunit_tmp,'(a)') system
      WRITE ( kunit_tmp,'(f20.12,i4)') box,ntype
      WRITE ( kunit_tmp,*) ( atypei(it) , it=1,ntype )
      WRITE ( kunit_tmp,*) ( natmi (it) , it=1,ntype )
      WRITE ( kunit_tmp,'(3f20.12)') ( rave ( 1, ia ) , rave ( 2 , ia ) , rave ( 3 , ia ) , ia = 1 , natm )
  CLOSE (kunit_tmp )
  endif
  deallocate ( rave ) 

  CLOSE ( kunit_TRAJFF )

! we read the configurations again 

  OPEN (UNIT = kunit_TRAJFF , FILE = 'TRAJFF')

  READ ( kunit_TRAJFF, * )  natm
  READ ( kunit_TRAJFF, * )  system
  READ ( kunit_TRAJFF, * )  box,ntype

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
      READ ( kunit_TRAJFF, * ) atype(i),rx(i),ry(i),rz(i)
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
    enddo
  enddo
  natmi  ( 0 ) = natm
  atypei ( 0 ) = 'ALL'

  ! ===============
  ! use direct sum
  ! ===============
  if ( longrange .eq. 'direct' )  CALL efg_DS( iconf , iastart , iaend )
  ! ===============
  ! use ewald sum
  ! ===============
  if ( longrange .eq. 'ewald' )   CALL efg_ES( iconf , iastart , iaend)
enddo

  ! write average distribution output 
  CALL efg_write_output( ncefg )

  CLOSE (kunit_EFGALL)
  CLOSE (kunit_EFGFF)
  CLOSE ( kunit_TRAJFF )

#endif fix_grid

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

SUBROUTINE efg_DS ( itime , iastart , iaend  )

  USE control,  ONLY :  myrank, numprocs, calc 
  USE config,   ONLY :  system , natm , natmi , atype , atypei , box , itype , rx , ry , rz , ntype , qia 
  USE io_file,  ONLY :  ionode , stdout , kunit_EFGALL , kunit_EFGALLIT1 , kunit_EFGALLIT2
  USE field,    ONLY :  qch
  USE prop,     ONLY :  nprop_print
  USE time

  implicit none

  INCLUDE 'mpif.h'

  ! global
  integer, intent(in) :: itime , iastart , iaend 
  ! local
  integer :: ia, ja, ierr , sumcount , it
  integer :: ncell
  integer :: kunit
  integer :: kunit_it(2) ! <--  not general enough should be ntype dependent
  double precision :: d , d2 , d5 , dm5
  double precision :: rxi , ryi , rzi , rxj , ryj , rzj , rxij , ryij , rzij
  double precision :: cutefgsq
  ! time
  double precision :: ttt1 , ttt2 , ttt3 

  ttt1 = MPI_WTIME(ierr)

  cutefgsq = cutefg * cutefg
  efg_ia = 0.0d0
  efg_ia_it = 0.0d0

! =========================================================
!  MAIN LOOP calculate EFG(i) for each atom i parallelized
! =========================================================

  sumcount = 0
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
                efg_ia_it( ia , it , 1 , 1 )  = efg_ia_it( ia , it , 1 , 1 ) - (3.0d0 * rxij * rxij - d2 ) * dm5
                efg_ia_it( ia , it , 2 , 2 )  = efg_ia_it( ia , it , 2 , 2 ) - (3.0d0 * ryij * ryij - d2 ) * dm5
                efg_ia_it( ia , it , 3 , 3 )  = efg_ia_it( ia , it , 3 , 3 ) - (3.0d0 * rzij * rzij - d2 ) * dm5
                efg_ia_it( ia , it , 1 , 2 )  = efg_ia_it( ia , it , 1 , 2 ) -  3.0d0 * rxij * ryij * dm5
                efg_ia_it( ia , it , 1 , 3 )  = efg_ia_it( ia , it , 1 , 3 ) -  3.0d0 * rxij * rzij * dm5
                efg_ia_it( ia , it , 2 , 3 )  = efg_ia_it( ia , it , 2 , 3 ) -  3.0d0 * ryij * rzij * dm5
              endif   

              sumcount = sumcount + 1
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
                efg_ia_it( ia , it , 1 , 1 )  = efg_ia_it( ia , it , 1 , 1 ) - (3.0d0 * rxij * rxij - d2 ) * dm5
                efg_ia_it( ia , it , 2 , 2 )  = efg_ia_it( ia , it , 2 , 2 ) - (3.0d0 * ryij * ryij - d2 ) * dm5
                efg_ia_it( ia , it , 3 , 3 )  = efg_ia_it( ia , it , 3 , 3 ) - (3.0d0 * rzij * rzij - d2 ) * dm5
                efg_ia_it( ia , it , 1 , 2 )  = efg_ia_it( ia , it , 1 , 2 ) -  3.0d0 * rxij * ryij * dm5
                efg_ia_it( ia , it , 1 , 3 )  = efg_ia_it( ia , it , 1 , 3 ) -  3.0d0 * rxij * rzij * dm5
                efg_ia_it( ia , it , 2 , 3 )  = efg_ia_it( ia , it , 2 , 3 ) -  3.0d0 * ryij * rzij * dm5
              endif   

                sumcount = sumcount + 1
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

  !======================================
  !  MERGE tensor from different proc
  !======================================
  efg_sum3=0.0d0
  CALL MPI_ALLREDUCE(efg_ia(:,1,1),efg_sum3(:,1,1),natm,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
  efg_ia(:,1,1) = efg_sum3(:,1,1)
  efg_sum3=0.0d0
  CALL MPI_ALLREDUCE(efg_ia(:,2,2),efg_sum3(:,2,2),natm,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
  efg_ia(:,2,2) = efg_sum3(:,2,2)
  efg_sum3=0.0d0
  CALL MPI_ALLREDUCE(efg_ia(:,3,3),efg_sum3(:,3,3),natm,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
  efg_ia(:,3,3) = efg_sum3(:,3,3)
  efg_sum3=0.0d0
  CALL MPI_ALLREDUCE(efg_ia(:,1,2),efg_sum3(:,1,2),natm,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
  efg_ia(:,1,2) = efg_sum3(:,1,2)
  efg_sum3=0.0d0
  CALL MPI_ALLREDUCE(efg_ia(:,1,3),efg_sum3(:,1,3),natm,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
  efg_ia(:,1,3) = efg_sum3(:,1,3)
  efg_sum3=0.0d0
  CALL MPI_ALLREDUCE(efg_ia(:,2,3),efg_sum3(:,2,3),natm,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
  efg_ia(:,2,3) = efg_sum3(:,2,3)

  efg_ia( : , 2, 1) = efg_ia( : , 1, 2)
  efg_ia( : , 3, 1) = efg_ia( : , 1, 3)
  efg_ia( : , 3, 2) = efg_ia( : , 2, 3)

  if ( lefg_it_contrib ) then
    do it=1,ntype
      efg_sum3=0.0d0
      CALL MPI_ALLREDUCE(efg_ia_it(:,it,1,1),efg_sum3(:,1,1),natm,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
      efg_ia_it(:,it,1,1) = efg_sum3(:,1,1)
      efg_sum3=0.0d0
      CALL MPI_ALLREDUCE(efg_ia_it(:,it,2,2),efg_sum3(:,2,2),natm,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
      efg_ia_it(:,it,2,2) = efg_sum3(:,2,2)
      efg_sum3=0.0d0
      CALL MPI_ALLREDUCE(efg_ia_it(:,it,3,3),efg_sum3(:,3,3),natm,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
      efg_ia_it(:,it,3,3) = efg_sum3(:,3,3)
      efg_sum3=0.0d0
      CALL MPI_ALLREDUCE(efg_ia_it(:,it,1,2),efg_sum3(:,1,2),natm,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
      efg_ia_it(:,it,1,2) = efg_sum3(:,1,2)
      efg_sum3=0.0d0
      CALL MPI_ALLREDUCE(efg_ia_it(:,it,1,3),efg_sum3(:,1,3),natm,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
      efg_ia_it(:,it,1,3) = efg_sum3(:,1,3)
      efg_sum3=0.0d0
      CALL MPI_ALLREDUCE(efg_ia_it(:,it,2,3),efg_sum3(:,2,3),natm,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
      efg_ia_it(:,it,2,3) = efg_sum3(:,2,3)
    enddo
    efg_ia_it( : , : , 2, 1) = efg_ia_it( : , : , 1, 2)
    efg_ia_it( : , : , 3, 1) = efg_ia_it( : , : , 1, 3)
    efg_ia_it( : , : , 3, 2) = efg_ia_it( : , : , 2, 3)
  endif

#ifdef debug
ia=1
CALL print_tensor( efg_ia( ia , : , : ) , 'TOTAL' )
#endif
  ! =======================================
  ! write efg for each atom in file EFGALL 
  ! =======================================
  if ( ionode  .and. lefgprintall .and. itime.ne.0 ) then
    WRITE ( kunit_EFGALL ,'(i5,f14.8,i5)') natm,box,ntype 
    WRITE ( kunit_EFGALL ,'(a10,i5)')  system,itime
      WRITE ( kunit_EFGALL ,'(a)') &
      '      ia type                   vxx                   vyy                   vzz                   vxy                   vxz                   vyz'
    do ia = 1 , natm
      WRITE ( kunit_EFGALL ,'(i8,2x,a3,6f22.16)') & 
        ia , atype ( ia ) , efg_ia ( ia , 1 , 1) , efg_ia ( ia , 2 , 2) , efg_ia ( ia , 3 , 3) , &
        efg_ia ( ia , 1 , 2) , efg_ia ( ia , 1 , 3) , efg_ia ( ia , 2 , 3)

    enddo
  endif

  kunit_it(1) = kunit_EFGALLIT1
  kunit_it(2) = kunit_EFGALLIT2
  if ( lefg_it_contrib ) then
    do it = 1 , 2 ! <- should be ntype dependent
      kunit = kunit_it(it)
      WRITE ( kunit ,'(i5,f14.8,i5)') natm,box,ntype
      WRITE ( kunit ,'(a10,i5)')  system,itime
      WRITE ( kunit ,'(a)') &
      '      ia type                   vxx                   vyy                   vzz                   vxy                   vxz                   vyz'
      do ia = 1 , natm
      WRITE ( kunit ,'(i8,2x,a3,6f22.16)') & 
        ia , atype ( ia ) , efg_ia_it ( ia , it , 1 , 1) , efg_ia_it ( ia , it , 2 , 2) , efg_ia_it ( ia , it , 3 , 3) , &
        efg_ia_it ( ia , it , 1 , 2) , efg_ia_it ( ia , it , 1 , 3) , efg_ia_it ( ia , it , 2 , 3)
      enddo
    enddo
  endif
! ====================================
!  no accumulation for the firts step
! ====================================
  if (itime .eq. 0 .and. calc.ne.'efg' ) then
    dibvzztot = 0
    dibetatot = 0
    dibUtot   = 0
  endif

  ttt2 = MPI_WTIME(ierr)
  efgtimetot1 = efgtimetot1 + (ttt2-ttt1)

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
 
SUBROUTINE efg_ES ( itime , iastart , iaend)

  USE control,    ONLY :  myrank , numprocs, calc
  USE config,     ONLY :  system , natm , natmi , atype , atypei , box , &
                          omega , itype , rx , ry , rz , ntype , qia 
  USE io_file,  ONLY :  ionode , stdout , kunit_EFGALL , kunit_EFGALLIT1 , kunit_EFGALLIT2
  USE constants,  ONLY :  pi , fpi , piroot , imag
  USE field,      ONLY :  qch
  USE kspace,     ONLY :  struc_fact
  USE prop,       ONLY :  nprop_print , nprop
  USE time

  implicit none

  INCLUDE 'mpif.h'

  ! global
  integer, intent(in) :: itime , iastart , iaend 

  ! local
  integer :: ia , ja ,  ierr , it 
  integer :: nxij , nyij , nzij
  double precision :: d , d2 , d4 , d3 , d5 , expon 
  double precision :: alpha2 , alpha3
  double precision :: allrealpart
  double precision :: T0 , T1 , T2  ! real part 
  double precision :: rxi , ryi , rzi , rxij , ryij , rzij 
  double precision :: cutoffsq
  double precision :: invbox 
  double precision :: ak, kx, ky, kz, kk
  integer :: ik
  double precision :: kri , onethird
  double precision :: ttt1 , ttt2 , ttt3 , ttt4 , ttt5
  double complex   :: rhon , carg , recarg , recarg_dgg 
  double precision, external :: errfc 
  integer :: kunit
  integer :: kunit_it(2) !<--should be ntype dependent

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
  onethird = 1/3.0D0
  invbox = 1.0d0 / box

  ! cutoff direct space part 
  cutoffsq = box * box * 0.25d0 
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
         T1 = ( 2.0d0 * alphaES * expon ) / d4 
         T2 = ( 4.0d0 * alpha3 * expon ) / d2 / 3.0d0
         allrealpart = qia(ja) * ( T0 + T1 + T2 )

         efg_ia_real( ia , 1 , 1 ) = efg_ia_real( ia , 1 , 1 ) - ( 3.0D0 * rxij * rxij - d2 ) * allrealpart
         efg_ia_real( ia , 2 , 2 ) = efg_ia_real( ia , 2 , 2 ) - ( 3.0D0 * ryij * ryij - d2 ) * allrealpart
         efg_ia_real( ia , 3 , 3 ) = efg_ia_real( ia , 3 , 3 ) - ( 3.0D0 * rzij * rzij - d2 ) * allrealpart
         efg_ia_real( ia , 1 , 2 ) = efg_ia_real( ia , 1 , 2 ) -   3.0d0 * rxij * ryij * allrealpart
         efg_ia_real( ia , 1 , 3 ) = efg_ia_real( ia , 1 , 3 ) -   3.0d0 * rxij * rzij * allrealpart
         efg_ia_real( ia , 2 , 3 ) = efg_ia_real( ia , 2 , 3 ) -   3.0d0 * ryij * rzij * allrealpart

         if ( lefg_it_contrib ) then
           it = itype(ja)
           efg_ia_real_it( ia , it , 1 , 1 ) = efg_ia_real_it( ia , it , 1 , 1 ) - ( 3.0D0 * rxij * rxij - d2 ) * allrealpart
           efg_ia_real_it( ia , it , 2 , 2 ) = efg_ia_real_it( ia , it , 2 , 2 ) - ( 3.0D0 * ryij * ryij - d2 ) * allrealpart
           efg_ia_real_it( ia , it , 3 , 3 ) = efg_ia_real_it( ia , it , 3 , 3 ) - ( 3.0D0 * rzij * rzij - d2 ) * allrealpart
           efg_ia_real_it( ia , it , 1 , 2 ) = efg_ia_real_it( ia , it , 1 , 2 ) -   3.0d0 * rxij * ryij * allrealpart
           efg_ia_real_it( ia , it , 1 , 3 ) = efg_ia_real_it( ia , it , 1 , 3 ) -   3.0d0 * rxij * rzij * allrealpart
           efg_ia_real_it( ia , it , 2 , 3 ) = efg_ia_real_it( ia , it , 2 , 3 ) -   3.0d0 * ryij * rzij * allrealpart
         endif

       endif

     enddo

  enddo atom1

  ttt3 = MPI_WTIME(ierr)
  efgtimetot1 = efgtimetot1 + (ttt3-ttt2)

! ==============================================
!            reciprocal space part
! ==============================================
  ! sum on atoms for one givn ik
atom2:  do ia = 1 , natm
    rxi = rx(ia)
    ryi = ry(ia)
    rzi = rz(ia)
#ifdef fix_grid
     rxi = rgrid(1,ia)
     ryi = rgrid(2,ia)
     rzi = rgrid(3,ia) 
#endif
    kpoint : do ik = 1, km_efg%nkcut 
      ! =================
      !   k-space  
      ! =================
      kx = km_efg%kpt(1,ik)
      ky = km_efg%kpt(2,ik)
      kz = km_efg%kpt(3,ik)
      kk = km_efg%kptk(ik)
      Ak = dexp( - kk * 0.25d0 / alpha2 ) 
      if ( km_efg%kptk(ik) .eq. 0 ) then 
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
      recarg_dgg =  rhon*carg*ak / kk 
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
          rhon = qch(it) * CONJG( km_efg%strf ( ik , it ) )

          recarg_dgg =  rhon*carg*ak / kk
          recarg = recarg_dgg * kk

          efg_ia_dual_it ( ia , it , 1 , 1 ) = efg_ia_dual_it ( ia , it , 1 , 1 ) +  3.0d0 * kx * kx * recarg_dgg - recarg
          efg_ia_dual_it ( ia , it , 2 , 2 ) = efg_ia_dual_it ( ia , it , 2 , 2 ) +  3.0d0 * ky * ky * recarg_dgg - recarg
          efg_ia_dual_it ( ia , it , 3 , 3 ) = efg_ia_dual_it ( ia , it , 3 , 3 ) +  3.0d0 * kz * kz * recarg_dgg - recarg
          efg_ia_dual_it ( ia , it , 1 , 2 ) = efg_ia_dual_it ( ia , it , 1 , 2 ) +  3.0d0 * kx * ky * recarg_dgg
          efg_ia_dual_it ( ia , it , 1 , 3 ) = efg_ia_dual_it ( ia , it , 1 , 3 ) +  3.0d0 * kx * kz * recarg_dgg
          efg_ia_dual_it ( ia , it , 2 , 3 ) = efg_ia_dual_it ( ia , it , 2 , 3 ) +  3.0d0 * ky * kz * recarg_dgg
          
        enddo
      endif

    enddo kpoint

  enddo atom2

  ttt4 = MPI_WTIME(ierr)
  efgtimetot2 = efgtimetot2 + (ttt4-ttt3)

!=============================
!  MERGE REAL PART
!=============================
  efg_sum3=0.0d0
  CALL MPI_ALLREDUCE(efg_ia_real(:,1,1),efg_sum3(:,1,1),natm,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
  efg_ia_real(:,1,1) = efg_sum3(:,1,1)
  efg_sum3=0.0d0
  CALL MPI_ALLREDUCE(efg_ia_real(:,2,2),efg_sum3(:,2,2),natm,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
  efg_ia_real(:,2,2) = efg_sum3(:,2,2)
  efg_sum3=0.0d0
  CALL MPI_ALLREDUCE(efg_ia_real(:,3,3),efg_sum3(:,3,3),natm,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
  efg_ia_real(:,3,3) = efg_sum3(:,3,3)
  efg_sum3=0.0d0
  CALL MPI_ALLREDUCE(efg_ia_real(:,1,2),efg_sum3(:,1,2),natm,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
  efg_ia_real(:,1,2) = efg_sum3(:,1,2)
  efg_sum3=0.0d0
  CALL MPI_ALLREDUCE(efg_ia_real(:,1,3),efg_sum3(:,1,3),natm,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
  efg_ia_real(:,1,3) = efg_sum3(:,1,3)
  efg_sum3=0.0d0
  CALL MPI_ALLREDUCE(efg_ia_real(:,2,3),efg_sum3(:,2,3),natm,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
  efg_ia_real(:,2,3) = efg_sum3(:,2,3)

  efg_ia_real( : , 2, 1) = efg_ia_real( : , 1, 2) 
  efg_ia_real( : , 3, 1) = efg_ia_real( : , 1, 3) 
  efg_ia_real( : , 3, 2) = efg_ia_real( : , 2, 3)

  if ( lefg_it_contrib ) then
    do it=1,ntype
      efg_sum3=0.0d0
      CALL MPI_ALLREDUCE(efg_ia_real_it(:,it,1,1),efg_sum3(:,1,1),natm,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
      efg_ia_real_it(:,it,1,1) = efg_sum3(:,1,1)
      efg_sum3=0.0d0
      CALL MPI_ALLREDUCE(efg_ia_real_it(:,it,2,2),efg_sum3(:,2,2),natm,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
      efg_ia_real_it(:,it,2,2) = efg_sum3(:,2,2)
      efg_sum3=0.0d0
      CALL MPI_ALLREDUCE(efg_ia_real_it(:,it,3,3),efg_sum3(:,3,3),natm,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
      efg_ia_real_it(:,it,3,3) = efg_sum3(:,3,3)
      efg_sum3=0.0d0
      CALL MPI_ALLREDUCE(efg_ia_real_it(:,it,1,2),efg_sum3(:,1,2),natm,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
      efg_ia_real_it(:,it,1,2) = efg_sum3(:,1,2)
      efg_sum3=0.0d0
      CALL MPI_ALLREDUCE(efg_ia_real_it(:,it,1,3),efg_sum3(:,1,3),natm,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
      efg_ia_real_it(:,it,1,3) = efg_sum3(:,1,3)
      efg_sum3=0.0d0
      CALL MPI_ALLREDUCE(efg_ia_real_it(:,it,2,3),efg_sum3(:,2,3),natm,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
      efg_ia_real_it(:,it,2,3) = efg_sum3(:,2,3)
    enddo
    efg_ia_real_it( : , : , 2, 1) = efg_ia_real_it( : , : , 1, 2) 
    efg_ia_real_it( : , : , 3, 1) = efg_ia_real_it( : , : , 1, 3) 
    efg_ia_real_it( : , : , 3, 2) = efg_ia_real_it( : , : , 2, 3)
  endif

  efg_ia_dual ( : , 2 , 1 ) = efg_ia_dual ( : , 1 , 2 )
  efg_ia_dual ( : , 3 , 1 ) = efg_ia_dual ( : , 1 , 3 )
  efg_ia_dual ( : , 3 , 2 ) = efg_ia_dual ( : , 2 , 3 )
  efg_ia_dual_real( : , : , :) = dble( efg_ia_dual( : , : , :)  ) 
  ! =======
  ! 4pi/V
  ! =======
  efg_ia_dual_real =  efg_ia_dual_real * fpi / omega / 3.0d0

  if ( lefg_it_contrib ) then
    efg_ia_dual_it ( : , : , 2 , 1 ) = efg_ia_dual_it ( : , : , 1 , 2 )
    efg_ia_dual_it ( : , : , 3 , 1 ) = efg_ia_dual_it ( : , : , 1 , 3 )
    efg_ia_dual_it ( : , : , 3 , 2 ) = efg_ia_dual_it ( : , : , 2 , 3 )
    efg_ia_dual_real_it( : , : , : , :) = dble( efg_ia_dual_it( : , : , : , :)  ) 
    ! =======
    ! 4pi/V
    ! =======
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
  
  ! =======================================
  ! write efg for each atom in file EFGALL 
  ! =======================================
  if ( ionode  .and. lefgprintall .and. itime.ne.0 ) then
    WRITE ( kunit_EFGALL ,'(i5,f14.8,i5)') natm,box,ntype
    WRITE ( kunit_EFGALL ,'(a10,i5)')  system,itime
      WRITE ( kunit_EFGALL ,'(a)') &
      '      ia type                   vxx                   vyy                   vzz                   vxy                   vxz                   vyz'
    do ia = 1 , natm
      WRITE ( kunit_EFGALL ,'(i8,2x,a3,6f22.16)') &
      ia , atype ( ia ) , efg_ia ( ia , 1 , 1) , efg_ia ( ia , 2 , 2) , efg_ia ( ia , 3 , 3) , &
      efg_ia ( ia , 1 , 2) , efg_ia ( ia , 1 , 3) , efg_ia ( ia , 2 , 3)
  
    enddo
  endif

  kunit_it(1) = kunit_EFGALLIT1
  kunit_it(2) = kunit_EFGALLIT2
  if ( lefg_it_contrib ) then
    do it = 1 , 2 ! <- should be ntype dependent
      kunit = kunit_it(it)
      WRITE ( kunit ,'(i5,f14.8,i5)') natm,box,ntype
      WRITE ( kunit ,'(a10,i5)')  system,itime
      WRITE ( kunit ,'(a)') &
      '      ia type                   vxx                   vyy                   vzz                   vxy                   vxz                   vyz'
      do ia = 1 , natm
        WRITE ( kunit ,'(i8,2x,a3,6f22.16)') &
        ia , atype ( ia ) , efg_ia_it ( ia , it , 1 , 1) , efg_ia_it ( ia , it , 2 , 2) , efg_ia_it ( ia , it , 3 , 3) , &
        efg_ia_it ( ia , it , 1 , 2) , efg_ia_it ( ia , it , 1 , 3) , efg_ia_it ( ia , it , 2 , 3)
      enddo
    enddo
  endif

! ====================================
!  no accumulation for the first step
! ====================================
  if ( itime .eq. 0 .and. calc.ne.'efg') then
    dibvzztot = 0
    dibetatot = 0
    dibUtot   = 0
  endif

  return

END SUBROUTINE efg_ES


!*********************** SUBROUTINE efg_write_output ********************************
!
! write distributions of EFG parameters and components to files
!
!******************************************************************************

SUBROUTINE efg_write_output ( kunit_eta , kunit_vzz , kunit_u )

  USE io_file,          ONLY :  ionode 
  USE config,           ONLY :  natm, natmi, ntype 
  USE constants,        ONLY :  dzero
  USE control,          ONLY :  longrange                
  USE time

  implicit none
 
  INCLUDE 'mpif.h'

  ! global
  integer, intent(in)  :: kunit_eta , kunit_vzz , kunit_u 

  ! local
  integer :: i , ierr , it 

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
             acfxx (0:ntype,0:ntcor) , acfyy (0:ntype,0:ntcor) , acfzz (0:ntype,0:ntcor) , acfet (0:ntype,0:ntcor) )
 
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

SUBROUTINE efg_stat ( kunit_input , kunit_output )

  USE config,           ONLY  : system , natm , natmi , ntype , itype , atype , atypei , omega , box , rho, config_alloc
  USE io_file,          ONLY  : ionode , stdout , kunit_OUTFF 
  USE prop,     ONLY :  nprop_print

  implicit none

  integer            :: i , ic , ia , it , ui
  integer, parameter :: lwork = 6
  integer            :: ifail
  integer            :: kvzz, keta ,ku
  integer            :: kunit_input , kunit_output 

  double precision                                 :: sq3 , sq32
  double precision, dimension (:,:), allocatable   :: U
  double precision                                 ::  w(3) , efgt(3,3)
  double precision                                 :: work(3 * lwork)
  double precision, dimension (:,:)  , allocatable :: nmr
  double precision, dimension (:)    , allocatable :: vzzmini , vzzmaxi ! ntype
  double precision, dimension (:)    , allocatable :: vzzm, vzzma , etam , pvzz, rho_z, sigmavzz, vzzsq ! ntype
  double precision :: vzzk , etak , uk

  !trash
  double precision :: aaaa
  integer :: iiii
  character :: XXXX

  !  ================
  !   allocation
  !  ================
  allocate ( U ( natm , 5 ) )                           ! Czjzek U component
  allocate ( nmr ( natm , 4 ) )                         ! 1 = Vxx ; 2 = Vyy ; 3 = Vzz ; 4 = eta
  allocate ( vzzmini ( 0:ntype ), vzzmaxi ( 0:ntype ) ) ! max min value ! not so important
  allocate ( vzzm ( 0:ntype ) )                         ! vzz mean value (per type)
  allocate ( vzzma ( 0:ntype ) )                        ! mean value of absolute vzz (per type)
  allocate ( vzzsq ( 0:ntype ) )                        ! square of vzz (per type)
  allocate ( etam ( 0:ntype ) )                         ! eta mean value (per type)
  allocate ( pvzz ( 0:ntype ) )                         ! proprtion of positiv vzz (per type)
  allocate ( rho_z ( 0:ntype ) )                        ! rho_z (per type)
  allocate ( sigmavzz ( 0:ntype ) )                     ! variance of vzz distribution (per type)
  ! ==============
  ! set to zero
  ! ==============
  nmr = 0.0D0 
  pvzz = 0.0d0
  vzzmini  = 0.0D0
  vzzmaxi  = 0.0D0
  etam = 0.0d0
  vzzm = 0.0d0
  vzzma = 0.0d0
  vzzsq = 0.0d0
  U = 0.0D0

  ! some constants
  sq3 = dsqrt( 3.0d0 )
  sq3 = 1.0d0 / sq3
  sq32 = sq3 * 0.5d0

  READ( kunit_input  , * ) iiii , aaaa , iiii 
  READ( kunit_input  , * ) XXXX , iiii
  READ( kunit_input  , * ) XXXX

  do ic = 1 , ncefg

    if ( ic .ne. 1 ) then
      READ( kunit_input , *) natm,box,ntype
      READ( kunit_input , *) system,iiii
      READ( kunit_input , *) XXXX
    endif

    do ia = 1 , natm

      READ( kunit_input , * ) iiii , atype(ia) , efgt(1,1) , efgt(2,2) , efgt(3,3) , &
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
      CALL nmr_convention( w , nmr(ia,:) , ia )
    
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
        if ( ionode .and. itype(ia) .eq. 1 ) WRITE ( stdout , * ) 'ERROR: out of bound distribvzz A'
        if ( ionode .and. itype(ia) .eq. 2 ) WRITE ( stdout , * ) 'ERROR: out of bound distribvzz B'
        if ( ionode ) WRITE ( stdout ,200) &
        ia , kvzz , nmr(ia,3) , vzzmini( itype ( ia ) ) , vzzmaxi ( itype ( ia ) ) , &
        nmr ( ia , 4 ), nmr ( ia , 1 ) , nmr ( ia , 2 ) , nmr ( ia , 3 )
        STOP
      endif
      if ( keta .lt. 0 .or. keta .gt. PANeta ) then
        if ( ionode ) WRITE ( stdout , * ) 'ERROR: out of bound distribeta'
        if ( ionode ) WRITE ( stdout ,210) ia, keta, nmr(ia,4), nmr(ia,4), nmr(ia,1), nmr(ia,2), nmr(ia,3)
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

  !  ================
  !   deallocation
  !  ================
  deallocate ( nmr  )   
  deallocate ( vzzmini , vzzmaxi )
  deallocate ( U )
  deallocate ( vzzm )
  deallocate ( vzzma )
  deallocate ( etam )
  deallocate ( pvzz )
  deallocate ( rho_z )
  deallocate ( sigmavzz )
  deallocate ( vzzsq )

100 FORMAT(I7,1X,' ATOM ',A3,'  minVZZ = ',F7.3,' maxVZZ = ',F7.3,& 
             ' <VZZ> =  ',F10.5,' <|VZZ|> =  ',F10.5,' mETA =  ',F10.5, &
          ' P(Vzz>0) =  ',F10.5, ' RHO_Z = ',F10.5)
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

  double precision ::  W(3),tmpw
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
