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
!#define debug
MODULE efg 

  USE kspace
  USE rspace

  implicit none

  logical :: lefgprintall                                               ! whether or not we want to print all the efg for each atoms and configurations in file EFGALL
  logical :: efg_it_contrib                                             ! if ones want to get the contributions to EFG separated in types
  integer :: ncefg                                                      ! number of configurations READ  for EFG calc (only when calc = 'efg')
  integer :: ntcor                                                      ! maximum number of steps for the acf calculation (calc = 'efg+acf')

  ! ===================
  !  direct summation
  ! ===================
  integer                                           :: ncelldirect
  double precision                                  :: cutefg             ! cut-off distance for efg calculation in the direct 
  double precision, dimension(:,:,:), allocatable   :: efg_ia             ! efg_tensor
  double precision, dimension(:,:,:,:), allocatable :: efg_ia_it          ! efg_tensor (it contribution)
  TYPE ( rmesh ) :: rm_efg

  ! ==================
  !  ewald summation
  ! ==================
  integer                                         :: ncellewald
  double precision                                :: alphaES            ! Ewald sum parameter 
  double precision, dimension(:,:,:), allocatable :: efg_ia_real        ! efg_tensor
  double complex  , dimension(:,:,:), allocatable :: efg_ia_dual        ! efg_tensor
  double precision, dimension(:,:,:), allocatable :: efg_ia_dual_real   ! efg_tensor
  double precision, dimension(:,:,:), allocatable :: efg_ia_total       ! efg_tensor

  TYPE ( kmesh ) :: km_efg

  ! ==============
  !  distribution
  ! ==============
  integer         , dimension(:,:)  , allocatable :: dibvzztot          ! vzz distribution dimension = (ntype , PAN)
  integer         , dimension(:,:)  , allocatable :: dibetatot          ! eta distribution dimension = (ntype , PAN) 
  integer         , dimension(:,:,:), allocatable :: dibUtot            ! U_i distribution dimension = ( 1:6 , ntype , PAN ) (Czjzek components)
  double precision                                :: resvzz             ! resolution in vzz distribution
  double precision                                :: reseta             ! resolution in eta distribution
  double precision                                :: resu               ! resolution in Ui ditribution
  double precision                                :: vzzmin             ! minimum value for vzz (distribution between [vzzmin, -vzzmin]
  double precision                                :: umin               ! minimum value of umin
  integer                                         :: PANeta             ! (internal) number of bins in eta distribution (related reseta)
  integer                                         :: PANvzz             ! (internal) number of bins in vzz distribution (related resvzz)
  integer                                         :: PANU               ! (internal) number of bins in Ui distribution (related resu)

!added 19-11-11
  integer         , dimension(:,:,:)  , allocatable :: dibvzztot_it     ! vzz distribution dimension = (ntype , PAN)
  integer         , dimension(:,:,:)  , allocatable :: dibetatot_it     ! eta distribution dimension = (ntype , PAN) 

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

  namelist /efgtag/  lefgprintall  , & 
                     efg_it_contrib, &
                     ncelldirect   , &
                     ncefg         , & 
                     ncellewald    , &
                     cutefg        , &
                     alphaES       , &
                     resvzz        , &
                     reseta        , &
                     resu          , &
                     vzzmin        , &
                     umin          , &
                     ntcor 

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
  if( ioerr .lt. 0 )  then
   if( ionode ) WRITE ( stdout, '(a)') 'ERROR reading input_file : efgtag section is absent'
   STOP
 elseif( ioerr .gt. 0 )  then
   if( ionode ) WRITE ( stdout, '(a)') 'ERROR reading input_file : efgtag wrong tag'
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
  efg_it_contrib = .false.
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

  USE field,    ONLY :  qa ,qb 
  USE control,  ONLY :  calc , longrange
  USE io_file,  ONLY :  ionode , stdout
  USE prop,     ONLY :  lefg
  USE config,   ONLY :  ntype

  implicit none

  if ( calc .eq. 'efg+acf' .and. ntype .gt. 2 ) then
    if ( ionode ) WRITE( stdout , '(a)' ) 'ERROR the subroutine efg_acf is not implemented for ntype > 2'
    STOP 
  endif

  if ( calc .eq. 'efg+acf' ) return

  if ( qa .eq. 0 .and. qb .ne. 0 ) then 
    if ( ionode ) WRITE ( stdout ,'(a,2f8.3)') 'ERROR: who requested EFG calculation without charge definition, something must be wrong' 
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
  if( vzzmin .eq. 0d0 .or. umin .eq. 0d0 ) then
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
  if( ionode ) then
    if( calc .eq. 'efg' .or. lefg ) then
                               WRITE ( kunit ,'(a)')           '=============================================================' 
                               WRITE ( kunit ,'(a)')           ''
                               WRITE ( kunit ,'(a)')           'electric field gradient:'
                               WRITE ( kunit ,'(a)')           'point charges calculation'
      if(longrange .eq. 'direct')  then
                               WRITE ( kunit ,'(a)')           'direct summation'
                               WRITE ( kunit ,'(a)')           'cubic cutoff in real space'
                               WRITE ( kunit ,'(a,f10.5)')     'distance cutoff           cutefg = ',cutefg
                               WRITE ( kunit ,'(a,i10)')       '-ncelldirect ... ncelldirect     = ',ncelldirect
                               WRITE ( kunit ,'(a,i10)')       'total number of cells            = ',( 2 * ncelldirect + 1 ) ** 3
      endif     
      if(longrange .eq. 'ewald')  then
                               WRITE ( kunit ,'(a)')           'ewald summation'
                               WRITE ( kunit ,'(a,f10.5)')     'alpha                            = ',alphaES
                               WRITE ( kunit ,'(a,f10.5)')     'cutoff in real part              = ',cutefg
                               WRITE ( kunit ,'(a,i10)')       'ncellewald x ncellewald x ncellewald               = ',ncellewald
                               WRITE ( kunit ,'(a,i10)')       '' 
                               WRITE ( kunit ,'(a,f10.5)')     'Note:this should hold alpha^2 * box^2 >> 1',alphaES*alphaES*box*box
      endif
      if( calc .eq. 'efg')     WRITE ( kunit ,'(a)')           'read config from file            : TRAJFF'
      if( calc .eq. 'efg')     WRITE ( kunit ,'(a,i10)')       'numbers of config read ncefg     = ',ncefg 
                               WRITE ( kunit ,'(a)')           'write averages values in file    : EFGFF '          

                               WRITE ( kunit ,'(a)')           'Distributions:'
                               WRITE ( kunit ,'(a)')           'eta distrib                : DTETAFF'
                               WRITE ( kunit ,'(a)')           'Vzz distrib                : DTVZZFF'
                               WRITE ( kunit ,'(a)')           'Ui components distrib      : DTIBUFF'           
     if(lefgprintall)          WRITE ( kunit ,'(a)')           'EFG for all atoms          : EFGALL '
                               WRITE ( kunit ,'(a)')           ''
                               WRITE ( kunit ,'(a)')                                 'distributions parameters:'
                               WRITE ( kunit ,'(a,f10.5)')                           'eta between    0.00000 and    1.00000 with resolution ',reseta
                               WRITE ( kunit ,'(a,f10.5,a,f10.5,a,f10.5)')           'vzz between ',vzzmin,' and ',-vzzmin,' with resolution ',resvzz
                               WRITE ( kunit ,'(a,f10.5,a,f10.5,a,f10.5)')           'Ui  between ',  umin,' and ',  -umin,' with resolution ',resu
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

  USE io_file,  ONLY :  kunit_EFGFF , kunit_EFGALL , kunit_EFGFFIT , kunit_EFGALLIT1 , kunit_EFGALLIT2
  USE control,  ONLY :  calc , longrange
  USE config,   ONLY :  natm, ntype , qia , qit , itype
  USE prop,     ONLY :  lefg
  USE field,    ONLY :  qa , qb

  implicit none

  ! local
  integer :: nkcut
  integer :: ncmax
  integer :: ia

  ! ======
  !  lefg
  ! ======
  if( .not. lefg .and. ( calc .ne. 'efg' ) ) return 

  qit(1)=qa
  qit(2)=qb
  
  if ( ntype .eq.  1 ) then
    do ia = 1 , natm
      qia(ia) = qit(1)
    enddo
  endif

  if ( ntype .eq.  2 ) then
    do ia = 1 , natm
      if(itype(ia) .eq. 1) then
        qia(ia) = qit(1)
      endif
      if(itype(ia) .eq. 2) then
        qia(ia) = qit(2)
      endif
    enddo
  endif

 
    OPEN (unit = kunit_EFGFF   ,file = 'EFGFF')
    OPEN (unit = kunit_EFGALL  ,file = 'EFGALL')

    if ( efg_it_contrib ) then
      OPEN (unit = kunit_EFGFFIT ,file = 'EFGFFIT')
      OPEN (unit = kunit_EFGALLIT1  ,file = 'EFGALLIT1')
      OPEN (unit = kunit_EFGALLIT2  ,file = 'EFGALLIT2')
    endif

    allocate( dibUtot(6,0:ntype,0:PANU) )
    allocate( dibvzztot(0:ntype,0:PANvzz) )
    allocate( dibetatot(0:ntype,0:PANeta) )  
    if ( efg_it_contrib ) then
      allocate( dibvzztot_it(0:ntype,ntype,0:PANvzz) )
      allocate( dibetatot_it(0:ntype,ntype,0:PANeta) )  
    endif
    dibvzztot = 0
    dibetatot = 0
    dibUtot = 0
    if ( efg_it_contrib ) then
      dibvzztot_it = 0
      dibetatot_it = 0
    endif

  ! ============
  !  direct sum
  ! ============
  if( longrange .eq. 'direct' ) then
    ncmax = ( 2 * ncelldirect + 1 ) ** 3
    rm_efg%ncmax=ncmax
    rm_efg%ncell=ncelldirect
    allocate( rm_efg%boxxyz( 3 , ncmax ) , rm_efg%lcell( ncmax ) )
    CALL direct_sum_init ( rm_efg )
    allocate( efg_ia( natm , 3 , 3 ) )
    allocate( efg_ia_it( natm , ntype , 3 , 3 ) )
    efg_ia    = 0.0d0 
    efg_ia_it = 0.0d0 
  endif

  ! ============
  !  ewald sum
  ! ============
  if( longrange .eq. 'ewald')  then
    allocate( efg_ia_real( natm , 3 , 3 ) )
    allocate( efg_ia_dual( natm , 3 , 3 ) ) 
    allocate( efg_ia_dual_real( natm , 3 , 3 ) ) 
    allocate( efg_ia_total( natm , 3 , 3 ) )
    efg_ia_real = 0.0d0 ; efg_ia_dual = 0.0d0 ; efg_ia_total = 0.0d0  
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

  if( .not. lefg .and. ( calc .ne. 'efg' ) ) return

  deallocate( dibUtot )
  deallocate( dibvzztot )
  deallocate( dibetatot )
  if (efg_it_contrib) then
    deallocate( dibvzztot_it )
    deallocate( dibetatot_it )
  endif
  if( longrange .eq. 'direct' ) then 
    deallocate( rm_efg%boxxyz , rm_efg%lcell )
    deallocate( efg_ia )  
    deallocate( efg_ia_it )  
  endif
  if( longrange .eq. 'ewald' )  then
    deallocate( efg_ia_real )
    deallocate( efg_ia_dual ) 
    deallocate( efg_ia_total )
    deallocate( km_efg%kptk , km_efg%kpt )
    deallocate ( km_efg%strf ) 
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

  USE io_file,  ONLY :  ionode , stdout , kunit_TRAJFF, kunit_EFGFF, kunit_EFGALL, kunit_OUTFF , kunit_EFGALLIT1 , kunit_EFGALLIT2 , kunit_EFGFFIT
  USE config,   ONLY :  system , natm , ntype , atype , rx , ry , rz , itype , atypei , natmi, rho , box , omega , config_alloc , qia , qit
  USE control,  ONLY :  longrange , myrank , numprocs
  USE field,    ONLY :  qa , qb, field_init

  implicit none

  ! local
  integer :: iastart , iaend
  integer :: i, ia, iconf , na
  double precision :: aaaa 
  integer :: iiii
  character * 60 :: cccc 


  OPEN (UNIT = kunit_TRAJFF,FILE = 'TRAJFF')

  READ ( kunit_TRAJFF, * )  natm 
  READ ( kunit_TRAJFF, * )  system
  READ ( kunit_TRAJFF, * )  box,ntype
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
  if( ionode ) then
    WRITE ( stdout , '(a)'      )    'Remind some parameters of the system:'
    WRITE ( stdout , '(a,i12)'  )    'natm  = ',natm
    WRITE ( stdout , '(a,i12)'  )    'ntype = ',ntype
    WRITE ( stdout , '(a,f12.5)')    'rho   = ',rho
    WRITE ( stdout , '(a,f12.5)')    'box   = ',box
    WRITE ( stdout , '(a,f12.5)')    'vol   = ',omega
    WRITE ( stdout , '(a)'      )    ''
  endif
  qit(1)=qa
  qit(2)=qb
  ! ============================================
  ! LOOP OVER CONFIGURATIONS 
  ! ============================================
  do iconf = 1, ncefg
    na = 0
    if( iconf .ne. 1 ) READ ( kunit_TRAJFF, * )  iiii 
    if( iconf .ne. 1 ) READ ( kunit_TRAJFF, * )  cccc
    if( iconf .ne. 1 ) READ ( kunit_TRAJFF, * )  aaaa,iiii 
    do i = 1,natm
      READ ( kunit_TRAJFF, * ) atype(i),rx(i),ry(i),rz(i)
      if(atype(i) .eq. 'A') na = na + 1  
      if(atype(i) .eq. 'A') itype(i)=1
      if(atype(i) .eq. 'B') itype(i)=2
    enddo
  
  if ( ntype .eq.  1 ) then
    natmi(0)=natm
    natmi(1)=na
    atypei(0)='ALL'
    atypei(1)='A'
    do ia = 1 , natm
      qia(ia) = qit(1)
    enddo
  endif

  if ( ntype .eq.  2 ) then
    natmi(0)=natm
    natmi(1)=na
    natmi(2)=natm-na
    atypei(0)='ALL'
    atypei(1)='A'
    atypei(2)='B'
    do ia = 1 , natm
      if(itype(ia) .eq. 1) then
        qia(ia) = qit(1)
      endif
      if(itype(ia) .eq. 2) then
        qia(ia) = qit(2)
      endif
    enddo
  endif

  ! ===============
  ! use direct sum
  ! ===============
  if ( longrange .eq. 'direct' )  CALL efg_DS( iconf , iconf , iastart , iaend )
  ! ===============
  ! use ewald sum
  ! ===============
  if ( longrange .eq. 'ewald' )   CALL efg_ES( iconf , iconf )

  enddo        
  ! write average distribution output 
  CALL efg_write_output( ncefg )

  CLOSE(kunit_EFGALL)
  CLOSE(kunit_EFGFF)
  CLOSE(kunit_TRAJFF)
  if ( efg_it_contrib ) then
    CLOSE(kunit_EFGALLIT1)
    CLOSE(kunit_EFGALLIT2)
    CLOSE(kunit_EFGFFIT)
  endif

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

SUBROUTINE efg_DS ( itime , nefg , iastart , iaend  )

  USE control,  ONLY :  myrank, numprocs, calc
  USE config,   ONLY :  system , natm , natmi , atype , atypei , box , itype , rx , ry , rz , ntype , qia , qit
  USE io_file,  ONLY :  ionode , stdout , kunit_EFGFF , kunit_EFGALL , kunit_EFGFFIT , kunit_EFGALLIT1 , kunit_EFGALLIT2
  USE field,    ONLY :  qa , qb
  USE prop,     ONLY :  nprop_print
  USE time

  implicit none

  INCLUDE 'mpif.h'

  ! global
  integer, intent(in) :: itime , nefg , iastart , iaend 

  ! local
  integer :: i, ia, ja, nb, ierr , it , ui , na , sumcount , itja
  integer :: ncell
  integer :: kvzz, keta
  double precision :: d , d2 , d5 , dm5
  double precision :: rxi , ryi , rzi , rxj , ryj , rzj , rxij , ryij , rzij
  double precision ::  W(3), EFGT(3,3) 
  double precision, dimension (:,:)  , allocatable :: WIT
  double precision, dimension (:,:,:), allocatable :: EFGTIT
  INTEGER, PARAMETER :: LWORK = 6
  double precision :: WORK(3 * LWORK)
  INTEGER :: IFAIL
  double precision :: cutefgsq
 
  double precision, dimension(:), allocatable :: vxx_sum, vyy_sum, vzz_sum, eta_sum
  double precision, dimension(:), allocatable :: vxx, vyy, vzz, eta
  double precision, dimension (:,:), allocatable :: U
  double precision :: sq3, sq32

  double precision, dimension(:,:)  , allocatable :: nmr
  double precision, dimension(:,:,:), allocatable :: nmrit
  
  double precision, dimension (:),   allocatable :: npvzzmin , npvzzmax ! ntype 
  double precision, dimension (:,:), allocatable :: vzzmini , vzzmaxi ! ntype proc
  double precision, dimension (:),   allocatable :: vzzmini2 , vzzmaxi2 ! proc 
  double precision, dimension (:),   allocatable :: vzzm, vzzma , etam , pvzz, rho_z, sigmavzz, vzzsq ! ntype
  double precision, dimension (:),   allocatable :: vzzm_sum, vzzma_sum , etam_sum , pvzz_sum, vzzsq_sum ! ntype
  double precision, dimension (:,:), allocatable :: vzzmin_sum , vzzmax_sum
  double precision, dimension (:),   allocatable :: vzzmin2_sum , vzzmax2_sum
  !  added 19-11-11
  double precision, dimension (:,:), allocatable :: pvzzit , vzzsqit , etamit , vzzmait , vzzmit , sigmavzzit , rho_zit
  
  double precision :: vzzk , etak
  integer :: ku,uk

  ! time
  double precision :: ttt1 , ttt2 , ttt3 

  ttt1 = MPI_WTIME(ierr)

  na = natmi(1)
  nb = natm - natmi(1)

  ! =================
  !  some constants 
  ! =================
  sq3 = dsqrt(3.0d00)
  sq3 = 1.0d0/sq3
  sq32 = sq3 * 0.5d0
  cutefgsq = cutefg * cutefg

  ! ==========================
  !  allocate some quantities
  ! ==========================
  ! 1 = Vxx ; 2 = Vyy ; 3 = Vzz ; 4 = eta
  allocate ( nmr    (natm , 4 ) ) 
  !new 19-11-11
  allocate ( nmrit  (natm , ntype , 4 ) ) 
  allocate ( EFGTIT ( ntype , 3 , 3 ) )
  allocate ( WIT   ( ntype , 3 ) )
  ! max min value ! not so important
  allocate ( vzzmini ( 0:ntype , 0:numprocs-1 ), vzzmaxi ( 0:ntype , 0:numprocs-1) )
  allocate ( npvzzmin (0:ntype) , npvzzmax (0:ntype) )
  ! Czjzek U component
  allocate ( U ( natm , 5 ) )
  ! vzz mean value (per type)
  allocate ( vzzm ( 0:ntype ) )
  allocate ( vzzmit ( ntype,0:ntype ) )
  ! mean value of absolute vzz (per type) 
  allocate ( vzzma ( 0:ntype ) )
  allocate ( vzzmait ( ntype, 0:ntype ) )
  ! square of vzz (per type)
  allocate ( vzzsq ( 0:ntype ) )
  allocate ( vzzsqit ( ntype, 0:ntype ) )
  ! eta mean value (per type)
  allocate ( etam ( 0:ntype ) ) 
  allocate ( etamit ( ntype, 0:ntype ) ) 
  ! proprtion of positiv vzz (per type)
  allocate ( pvzz ( 0:ntype ) ) 
  allocate ( pvzzit ( ntype , 0:ntype ) ) 
  ! rho_z (per type)
  allocate ( rho_z ( 0:ntype ) )
  allocate ( rho_zit ( ntype,0:ntype ) )
  ! variance of vzz distribution (per type)
  allocate ( sigmavzz ( 0:ntype ) )
  allocate ( sigmavzzit ( ntype,0:ntype ) )

  ! ==============
  ! set to zero
  ! ==============
  nmr = 0.0D0
  nmrit = 0.0D0
  pvzz = 0.0d0
  pvzzit = 0.0d0
  vzzmini  = 0.0D0
  vzzmaxi  = 0.0D0
  npvzzmin = 30000.0D0
  npvzzmax = -30000.0D0
  vzzmini ( : , myrank ) = 30000.0D0
  vzzmaxi ( : , myrank ) = -30000.0D0
  etam = 0.0d0
  etamit = 0.0d0
  vzzm = 0.0d0
  vzzmit = 0.0d0
  vzzma = 0.0d0
  vzzmait= 0.0d0
  vzzsq = 0.0d0
  vzzsqit = 0.0d0
  U = 0.0D0
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
         if( rm_efg%lcell(ncell) .eq. 1) then

          do ja = 1 , natm
            rxj  = rx(ja) + rm_efg%boxxyz(1,ncell)
            ryj  = ry(ja) + rm_efg%boxxyz(2,ncell)
            rzj  = rz(ja) + rm_efg%boxxyz(3,ncell)
            rxij = rxi - rxj
            ryij = ryi - ryj
            rzij = rzi - rzj
            d2   = rxij * rxij + ryij * ryij + rzij * rzij

            if( d2 .lt. cutefgsq ) then
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

              ! added 19-11-11
              if ( efg_it_contrib ) then
                itja = itype(ja)  
                efg_ia_it( ia , itja , 1 , 1 )  = efg_ia_it( ia , itja , 1 , 1 ) - (3.0d0 * rxij * rxij - d2 ) * dm5
                efg_ia_it( ia , itja , 2 , 2 )  = efg_ia_it( ia , itja , 2 , 2 ) - (3.0d0 * ryij * ryij - d2 ) * dm5
                efg_ia_it( ia , itja , 3 , 3 )  = efg_ia_it( ia , itja , 3 , 3 ) - (3.0d0 * rzij * rzij - d2 ) * dm5
                efg_ia_it( ia , itja , 1 , 2 )  = efg_ia_it( ia , itja , 1 , 2 ) -  3.0d0 * rxij * ryij * dm5
                efg_ia_it( ia , itja , 1 , 3 )  = efg_ia_it( ia , itja , 1 , 3 ) -  3.0d0 * rxij * rzij * dm5
                efg_ia_it( ia , itja , 2 , 3 )  = efg_ia_it( ia , itja , 2 , 3 ) -  3.0d0 * ryij * rzij * dm5
              endif   

              sumcount = sumcount + 1
            endif ! d2.lt.cutefgsq
 
          enddo ! ja
        endif 
         ! =======================================
         !  ia and ja in the same cell (ia.ne.ja)
         ! =======================================
        if( rm_efg%lcell(ncell) .eq. 0) then

          do ja = 1,natm

            if(ja.ne.ia) then
              rxj  = rx(ja)
              ryj  = ry(ja)
              rzj  = rz(ja)
              rxij = rxi - rxj
              ryij = ryi - ryj
              rzij = rzi - rzj
              d2   = rxij * rxij + ryij * ryij + rzij * rzij
 
              if(d2.lt.cutefgsq) then
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

              ! added 19-11-11
              if ( efg_it_contrib ) then
                itja = itype(ja)  
                efg_ia_it( ia , itja , 1 , 1 )  = efg_ia_it( ia , itja , 1 , 1 ) - (3.0d0 * rxij * rxij - d2 ) * dm5
                efg_ia_it( ia , itja , 2 , 2 )  = efg_ia_it( ia , itja , 2 , 2 ) - (3.0d0 * ryij * ryij - d2 ) * dm5
                efg_ia_it( ia , itja , 3 , 3 )  = efg_ia_it( ia , itja , 3 , 3 ) - (3.0d0 * rzij * rzij - d2 ) * dm5
                efg_ia_it( ia , itja , 1 , 2 )  = efg_ia_it( ia , itja , 1 , 2 ) -  3.0d0 * rxij * ryij * dm5
                efg_ia_it( ia , itja , 1 , 3 )  = efg_ia_it( ia , itja , 1 , 3 ) -  3.0d0 * rxij * rzij * dm5
                efg_ia_it( ia , itja , 2 , 3 )  = efg_ia_it( ia , itja , 2 , 3 ) -  3.0d0 * ryij * rzij * dm5
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
    ! addes 19-11-11
    if ( efg_it_contrib ) then
      efg_ia_it( ia , : , 2 , 1 ) = efg_ia_it( ia , : , 1 , 2 )
      efg_ia_it( ia , : , 3 , 1 ) = efg_ia_it( ia , : , 1 , 3 )
      efg_ia_it( ia , : , 3 , 2 ) = efg_ia_it( ia , : , 2 , 3 )
    endif

    EFGT=efg_ia(ia,:,:)
    if ( efg_it_contrib ) then 
      EFGTIT(1,:,:)=efg_ia_it(ia,1,:,:)
      EFGTIT(2,:,:)=efg_ia_it(ia,2,:,:)
    endif

#ifdef debug
     if(ia.eq.1) CALL print_tensor(EFGT,'DIREC')
     if ( efg_it_contrib ) then
       if(ia.eq.1) CALL print_tensor(EFGTIT(1,:,:),'DRITA')
       if(ia.eq.1) CALL print_tensor(EFGTIT(2,:,:),'DRITB')
       if(ia.eq.1) CALL print_tensor(EFGTIT(2,:,:)+EFGTIT(1,:,:),'DITAB')
     endif
#endif 

    ! =====================================================================
    !  Czjzek components (see J. Phys.: Condens. Matter 10 (1998). p10719)
    ! =====================================================================
    U( ia, 1) = EFGT(3,3) * 0.5d0
    U( ia, 2) = EFGT(1,3) * sq3
    U( ia, 3) = EFGT(2,3) * sq3
    U( ia, 4) = EFGT(1,2) * sq3
    U( ia, 5) = ( EFGT(1,1) - EFGT(2,2) ) * sq32

    ! =================
    !  diagonalisation
    ! =================
    CALL DSYEV('N','U',3,EFGT,3,w,work,3 * lwork,ifail)
    if (ifail.ne.0) then
      if( ionode ) WRITE ( stdout , * ) 'ERROR: DSYEV, STOP in efg MODULE (EFGT)'
      STOP 
    endif
    ! added 19-11-11
    if ( efg_it_contrib ) then
      work=0.0d0
      CALL DSYEV('N','U',3,EFGTIT(1,:,:),3,wit(1,:),work,3 * lwork,ifail)    
      if (ifail.ne.0) then
        if( ionode ) WRITE ( stdout , * ) 'ERROR: DSYEV, STOP in efg MODULE (EFGTIT1)'
        STOP 
      endif
      work=0.0d0
      CALL DSYEV('N','U',3,EFGTIT(2,:,:),3,wit(2,:),work,3 * lwork,ifail)    
      if (ifail.ne.0) then
        if( ionode ) WRITE ( stdout , * ) 'ERROR: DSYEV, STOP in efg MODULE (EFGTIT2)'
        STOP 
      endif
    endif

   ! =============================
   ! NMR conventions test 
   ! changed and added 19-11-11
   ! =============================
   call nmr_convention( w , nmr(ia,:) , ia )
   if ( efg_it_contrib ) then
     call nmr_convention( wit(1,:) , nmrit(ia,1,:) , ia )
     call nmr_convention( wit(2,:) , nmrit(ia,2,:) , ia )
   endif

#ifdef debug
   print*,nmr(ia,1),nmrit(ia,1,1),nmrit(ia,2,1)   
   print*,nmr(ia,2),nmrit(ia,1,2),nmrit(ia,2,2)
   print*,nmr(ia,3),nmrit(ia,1,3),nmrit(ia,2,3)   
   print*,nmr(ia,4),nmrit(ia,1,4),nmrit(ia,2,4)   
   print*,''
   print*,w(1),w(2),w(3)     
   print*,wit(1,1),wit(1,2),wit(1,3)     
   print*,wit(2,1),wit(2,2),wit(2,3)     
#endif

    ! ========================================================
    !  statistic for all sites index 0 in array of size ntype
    ! ========================================================
    vzzm(0) = vzzm(0) + nmr(ia,3)
    etam(0) = etam(0) + nmr(ia,4)
    vzzma(0) = vzzma(0) + dabs(nmr(ia,3))
    if(nmr(ia,3).le.vzzmini(0,myrank))  vzzmini(0,myrank) = nmr(ia,3)
    if(nmr(ia,3).ge.vzzmaxi(0,myrank))  vzzmaxi(0,myrank) = nmr(ia,3)
    if(nmr(ia,3).ge.0) pvzz(0) = pvzz(0) + 1.d0
    vzzsq(0) = vzzsq(0) + nmr(ia,3) * nmr(ia,3)
        
    do it = 1 , ntype 
      if(itype(ia).eq.it) then
        vzzm(it)  = vzzm(it)  + nmr(ia,3)
        vzzma(it) = vzzma(it) + dabs(nmr(ia,3))
        etam(it)  = etam(it)  + nmr(ia,4)
        if(nmr(ia,3).le.vzzmini(it,myrank))  vzzmini(it,myrank) = nmr(ia,3)
        if(nmr(ia,3).ge.vzzmaxi(it,myrank))  vzzmaxi(it,myrank) = nmr(ia,3)
        vzzsq(it) = vzzsq(it) + nmr(ia,3) * nmr(ia,3)
        if(nmr(ia,3).ge.0.0d0 ) then
          pvzz(it) = pvzz(it) + 1.d0
        endif
      endif ! it 
    enddo

   ! added 19-11-11
    if ( efg_it_contrib ) then
      do itja = 1 , ntype
        vzzmit(itja,0) = vzzmit(itja,0) + nmrit(ia,itja,3)
        etamit(itja,0) = etamit(itja,0) + nmrit(ia,itja,4)
        vzzmait(itja,0) = vzzmait(itja,0) + dabs(nmrit(ia,itja,3))
        if(nmrit(ia,itja,3).ge.0) pvzzit(itja,0) = pvzzit(itja,0) + 1.d0
        vzzsqit(itja,0) = vzzsqit(itja,0) + nmrit(ia,itja,3) * nmrit(ia,itja,3)
            
        do it = 1 , ntype 
          if(itype(ia).eq.it) then
            vzzmit(itja,it)  = vzzmit(itja,it)  + nmrit(ia,itja,3)
            vzzmait(itja,it) = vzzmait(itja,it) + dabs(nmrit(ia,itja,3))
            etamit(itja,it)  = etamit(itja,it)  + nmrit(ia,itja,4)
            vzzsqit(itja,it) = vzzsqit(itja,it) + nmrit(ia,itja,3) * nmrit(ia,itja,3)
            if(nmrit(ia,itja,3).ge.0.0d0 ) then
              pvzzit(itja,it) = pvzzit(itja,it) + 1.d0
            endif
          endif ! it 
        enddo
     enddo
    endif

  enddo atom

  ttt2 = MPI_WTIME(ierr)
  efgtimetot1 = efgtimetot1 + (ttt2-ttt1)

! ==========================================
!  MERGE, CALCULATE DISTRIBUTION AND OUTPUT 
! ==========================================

  allocate( vzzmin_sum(0:ntype,0:numprocs-1) , vzzmax_sum(0:ntype,0:numprocs-1) , vzzmax2_sum(0:numprocs-1) , vzzmin2_sum(0:numprocs-1) )
  allocate( vzzmaxi2(0:numprocs-1) , vzzmini2(0:numprocs-1) )

  do it = 0,ntype 
    vzzmax2_sum=0.0d0
    vzzmin2_sum=0.0d0
    vzzmaxi2 = vzzmaxi(it,:)
    vzzmini2 = vzzmini(it,:)
    CALL MPI_ALLREDUCE(vzzmaxi2,vzzmax2_sum,numprocs,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
    CALL MPI_ALLREDUCE(vzzmini2,vzzmin2_sum,numprocs,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
    vzzmax_sum(it,:)  = vzzmax2_sum
    vzzmin_sum(it,:)  = vzzmin2_sum
  enddo

  do i = 0,numprocs-1
    if(vzzmin_sum(0,i).le.npvzzmin(0)) npvzzmin(0) = vzzmin_sum(0,i)
    if(vzzmax_sum(0,i).ge.npvzzmax(0)) npvzzmax(0) = vzzmax_sum(0,i)
    do it=1,ntype
      if(vzzmin_sum(it,i).le.npvzzmin(it)) npvzzmin(it) = vzzmin_sum(it,i)
      if(vzzmax_sum(it,i).ge.npvzzmax(it)) npvzzmax(it) = vzzmax_sum(it,i)
    enddo
  enddo

  deallocate( vzzmaxi2 , vzzmini2 )
  deallocate (vzzmin_sum , vzzmax_sum , vzzmin2_sum, vzzmax2_sum) 
  allocate ( vzzm_sum ( 0:ntype ) , vzzma_sum ( 0:ntype ) , vzzsq_sum ( 0:ntype ) , etam_sum ( 0:ntype ) , pvzz_sum ( 0:ntype ) )
  pvzz_sum=0.0d0
 
  do it = 0 , ntype
    CALL MPI_ALLREDUCE(vzzm(it),vzzm_sum(it),1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
    CALL MPI_ALLREDUCE(vzzma(it),vzzma_sum(it),1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
    CALL MPI_ALLREDUCE(vzzsq(it),vzzsq_sum(it),1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
    CALL MPI_ALLREDUCE(etam(it),etam_sum(it),1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
    CALL MPI_ALLREDUCE(pvzz(it),pvzz_sum(it),1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
  enddo

  vzzm  = vzzm_sum
  vzzsq = vzzsq_sum
  vzzma = vzzma_sum
  etam  = etam_sum
  pvzz  = pvzz_sum

  !added 19-11-11
  if ( efg_it_contrib ) then
    do itja = 1 , ntype
      vzzm_sum = 0.0d0
      vzzsq_sum = 0.0d0
      vzzma_sum = 0.0d0
      etam_sum = 0.0d0
      pvzz_sum = 0.0d0
      do it = 0 , ntype
        CALL MPI_ALLREDUCE(vzzmit(itja,it),vzzm_sum(it),1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
        CALL MPI_ALLREDUCE(vzzmait(itja,it),vzzma_sum(it),1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
        CALL MPI_ALLREDUCE(vzzsqit(itja,it),vzzsq_sum(it),1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
        CALL MPI_ALLREDUCE(etamit(itja,it),etam_sum(it),1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
        CALL MPI_ALLREDUCE(pvzzit(itja,it),pvzz_sum(it),1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
      enddo
     vzzmit(itja,:)  = vzzm_sum
     vzzsqit(itja,:) = vzzsq_sum
     vzzmait(itja,:) = vzzma_sum
     etamit(itja,:)  = etam_sum
     pvzzit(itja,:)  = pvzz_sum
    enddo
  endif

  deallocate ( vzzm_sum , vzzma_sum , vzzsq_sum  , etam_sum , pvzz_sum )

  if(ntype.eq.1) then
    vzzm  = vzzm  / dble( natm )
    vzzsq = vzzsq / dble( natm )
    vzzma = vzzma / dble( natm )
    etam  = etam  / dble( natm )
    pvzz  = pvzz  / dble( natm )
  else
    vzzm(0)  = vzzm(0)  / dble( natm )
    vzzsq(0) = vzzsq(0) / dble( natm )
    vzzma(0) = vzzma(0) / dble( natm )
    etam(0)  = etam(0)  / dble( natm )
    pvzz(0)  = pvzz(0)  / dble( natm )
    do it=1,ntype
      vzzm(it)  = vzzm(it)  / dble( natmi(it) )
      vzzsq(it) = vzzsq(it) / dble( natmi(it) )
      vzzma(it) = vzzma(it) / dble( natmi(it) )
      etam(it)  = etam(it)  / dble( natmi(it) )
      pvzz(it)  = pvzz(it)  / dble( natmi(it) )
    enddo
  !added 19-11-11
    if (efg_it_contrib) then 
      do itja=1,ntype
        vzzmit(itja,0)  = vzzmit(itja,0)  / dble( natm )
        vzzsqit(itja,0) = vzzsqit(itja,0) / dble( natm )
        vzzmait(itja,0) = vzzmait(itja,0) / dble( natm )
        etamit(itja,0)  = etamit(itja,0)  / dble( natm )
        pvzzit(itja,0)  = pvzzit(itja,0)  / dble( natm )
        do it=1,ntype
          vzzmit(itja,it)  = vzzmit(itja,it)  / dble( natmi(it) )
          vzzsqit(itja,it) = vzzsqit(itja,it) / dble( natmi(it) )
          vzzmait(itja,it) = vzzmait(itja,it) / dble( natmi(it) )
          etamit(itja,it)  = etamit(itja,it)  / dble( natmi(it) )
          pvzzit(itja,it)  = pvzzit(itja,it)  / dble( natmi(it) )
        enddo
      enddo
    endif
  endif

  sigmavzz = vzzsq - ( vzzma * vzzma )
  sigmavzz = dsqrt( sigmavzz )
  rho_z = sigmavzz / vzzma
  if (efg_it_contrib) then
    sigmavzzit = vzzsqit - ( vzzmait * vzzmait )
    sigmavzzit = dsqrt( sigmavzzit )
    rho_zit = sigmavzzit / vzzmait
  endif

! ================================
!  output average values ( time )
! ================================
    if( ionode .and. (mod(nefg,nprop_print).eq.0)) WRITE ( stdout ,'(a)') 'ALL TYPES'
  do it=1,ntype
    if( ionode .and. (mod(nefg,nprop_print).eq.0)) WRITE ( stdout ,100) &
                     itime , atypei(it) , npvzzmin(it) , npvzzmax(it) , vzzm(it) , vzzma(it) , etam(it) , pvzz(it) , rho_z(it)
    if( ionode ) WRITE ( kunit_EFGFF,100) &
                     itime , atypei(it) , npvzzmin(it) , npvzzmax(it) , vzzm(it) , vzzma(it) , etam(it) , pvzz(it) , rho_z(it)
  enddo
  if( ionode  .and. ntype .ne. 1 .and. (mod(nefg,nprop_print).eq.0) ) WRITE ( stdout ,100) &
                                                             itime , atypei(0) , npvzzmin(0) , npvzzmax(0) , vzzm(0) , vzzma(0) , etam(0) , pvzz(0) , rho_z(0)
  if( ionode  .and. ntype .ne. 1 ) WRITE ( kunit_EFGFF ,100) itime , atypei(0) , npvzzmin(0) , npvzzmax(0) , vzzm(0) , vzzma(0) , etam(0) , pvzz(0) , rho_z(0)


! added 19-11-11
  if ( efg_it_contrib ) then
  do itja=1,ntype
    if ( ionode .and. (mod(nefg,nprop_print).eq.0) )  WRITE ( stdout ,'(a,i5)') 'ITJA = ',itja
    do it=1,ntype
      if( ionode .and. (mod(nefg,nprop_print).eq.0)) WRITE ( stdout ,100) &
                        itime , atypei(it) , npvzzmin(it) , npvzzmax(it) , vzzmit(itja,it) , vzzmait(itja,it) , etamit(itja,it) , pvzzit(itja,it) , rho_zit(itja,it)
      if( ionode ) WRITE ( kunit_EFGFFIT,100) &
                        itime , atypei(it) , npvzzmin(it) , npvzzmax(it) , vzzmit(itja,it) , vzzmait(itja,it) , etamit(itja,it) , pvzzit(itja,it) , rho_zit(itja,it)
    enddo
  enddo
  endif


! =======================================================
!  DISTRIBUTION CALCULATION LOOP I PARALLEL (NTYPE = 2)
! =======================================================
do ia=iastart , iaend

  !Ui
  do ui=1,5
    uk = (U(ia,ui)-umin)/resu
    ku = int(uk) + 1
    if(ku.lt.0.or.ku.gt.PANU) then
      if( ionode ) WRITE ( stdout , * ) 'ERROR: out of bound dibU1'
      if( ionode ) WRITE ( stdout ,310) i,ku,U(i,ui),umin,dabs(umin)
      STOP 
    endif
    dibUtot(ui,0,ku) = dibUtot(ui,0,ku) + 1
    do it=1,2
      if (itype(ia).eq.it) then
        dibUtot(ui,it,ku) = dibUtot(ui,it,ku) + 1
      endif
    enddo
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
  
  ! ====================
  !  test out of bound
  ! ====================
  if( keta .eq. PANeta + 1) then
    print*,'str',keta,etak,nmr(ia,4)
    keta = PANeta
  endif
  if( kvzz .lt. 0 .or. kvzz .gt. PANvzz ) then
    if(  ionode .and. itype(ia) .eq. 1 ) WRITE ( stdout , * ) 'ERROR: out of bound distribvzz A'
    if(  ionode .and. itype(ia) .eq. 2 ) WRITE ( stdout , * ) 'ERROR: out of bound distribvzz B'
    if(  ionode ) WRITE ( stdout ,200) ia,kvzz,nmr(ia,3),npvzzmin(itype(ia)),npvzzmax(itype(ia)),nmr(ia,4),nmr(ia,1),nmr(ia,2),nmr(ia,3)
    STOP 
  endif
  if( keta .lt. 0 .or. keta .gt. PANeta + 1) then
    if( ionode ) WRITE ( stdout , * ) 'ERROR: out of bound distribeta'
    if( ionode ) WRITE ( stdout ,210) ia, keta, nmr(ia,4), nmr(ia,4), nmr(ia,1), nmr(ia,2), nmr(ia,3)
    STOP 
  endif

  dibvzztot(0,kvzz) = dibvzztot(0,kvzz) + 1
  dibetatot(0,keta) = dibetatot(0,keta) + 1
  do it=1,ntype
    if( itype(ia) .eq. it ) then
      dibvzztot(it,kvzz) = dibvzztot(it,kvzz) + 1
      dibetatot(it,keta) = dibetatot(it,keta) + 1
    endif
  enddo

!added 19-11-11
  if ( efg_it_contrib ) then

    do itja = 1, ntype
      ! ========================
      !  quadrupolar parameters
      ! ========================
      vzzk = (nmrit(ia,itja,3)-vzzmin) / resvzz
      etak = nmrit(ia,itja,4) / reseta
      kvzz = int(vzzk) 
      keta = int(etak) 
  
      ! ====================
      !  test out of bound
      ! ====================
      if( keta .eq. PANeta + 1) then
        print*,'str',keta,etak,nmr(ia,4)
        keta = PANeta
      endif
      if( kvzz .lt. 0 .or. kvzz .gt. PANvzz ) then
        if(  ionode .and. itype(ia) .eq. 1 ) WRITE ( stdout , * ) 'ERROR: out of bound distribvzz A'
        if(  ionode .and. itype(ia) .eq. 2 ) WRITE ( stdout , * ) 'ERROR: out of bound distribvzz B'
        if(  ionode ) WRITE ( stdout ,200) ia,kvzz,nmr(ia,3),npvzzmin(itype(ia)),npvzzmax(itype(ia)),nmr(ia,4),nmr(ia,1),nmr(ia,2),nmr(ia,3)
        STOP 
      endif
      if( keta .lt. 0 .or. keta .gt. PANeta + 1) then
        if( ionode ) WRITE ( stdout , * ) 'ERROR: out of bound distribeta'
        if( ionode ) WRITE ( stdout ,210) ia, keta, nmr(ia,4), nmr(ia,4), nmr(ia,1), nmr(ia,2), nmr(ia,3)
        STOP 
      endif

      dibvzztot_it(0,itja,kvzz) = dibvzztot_it(0,itja,kvzz) + 1
      dibetatot_it(0,itja,keta) = dibetatot_it(0,itja,keta) + 1
      do it=1,ntype
        if( itype(ia) .eq. it ) then
          dibvzztot_it(it,itja,kvzz) = dibvzztot_it(it,itja,kvzz) + 1
          dibetatot_it(it,itja,keta) = dibetatot_it(it,itja,keta) + 1
        endif
      enddo
    enddo
  endif
enddo  !ia 

! ============================
!  END of distributions loop
! ============================
  allocate( eta_sum ( natm ) , vxx_sum ( natm ), vyy_sum ( natm ), vzz_sum ( natm ) )
  allocate( eta ( natm ) , vxx ( natm ), vyy ( natm ), vzz ( natm ) )

  eta_sum = 0.0D0
  vxx_sum = 0.0D0
  vyy_sum = 0.0D0
  vzz_sum = 0.0D0    

  ! ===================================
  !  merge eta(i),vxx(i),vyy(i),vzz(i) 
  !  where i is the atom index
  ! ===================================
  eta = nmr(:,4)
  vxx = nmr(:,1)
  vyy = nmr(:,2)
  vzz = nmr(:,3)
  CALL MPI_ALLREDUCE(vxx,vxx_sum,natm,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
  CALL MPI_ALLREDUCE(vyy,vyy_sum,natm,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
  CALL MPI_ALLREDUCE(vzz,vzz_sum,natm,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
  CALL MPI_ALLREDUCE(eta,eta_sum,natm,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
  nmr(:,1) = vxx_sum 
  nmr(:,2) = vyy_sum
  nmr(:,3) = vzz_sum
  nmr(:,4) = eta_sum

  if( ionode  .and. lefgprintall .and. itime.ne.0 ) then
    WRITE ( kunit_EFGALL ,'(i5,f14.8,i5)') natm,box,ntype 
    WRITE ( kunit_EFGALL ,'(a10,i5)')  system,itime
    WRITE ( kunit_EFGALL ,'(a)') '      ia  type     vxx           vyy          vzz           eta           rx            ry            rz'
    do i = 1,natm
      WRITE ( kunit_EFGALL ,'(i8,2x,a3,7f14.8)') i,atype(i),nmr(i,1),nmr(i,2),nmr(i,3),nmr(i,4),rx(i),ry(i),rz(i)
    enddo
  endif

  if ( efg_it_contrib ) then
    itja =1
    eta_sum = 0.0D0
    vxx_sum = 0.0D0
    vyy_sum = 0.0D0
    vzz_sum = 0.0D0
    ! ===================================
    !  merge eta(i),vxx(i),vyy(i),vzz(i) 
    !  where i is the atom index
    ! ===================================
    eta = nmrit(:,itja,4)
    vxx = nmrit(:,itja,1)
    vyy = nmrit(:,itja,2)
    vzz = nmrit(:,itja,3)
    CALL MPI_ALLREDUCE(vxx,vxx_sum,natm,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
    CALL MPI_ALLREDUCE(vyy,vyy_sum,natm,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
    CALL MPI_ALLREDUCE(vzz,vzz_sum,natm,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
    CALL MPI_ALLREDUCE(eta,eta_sum,natm,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
    nmrit(:,itja,1) = vxx_sum
    nmrit(:,itja,2) = vyy_sum
    nmrit(:,itja,3) = vzz_sum
    nmrit(:,itja,4) = eta_sum

    WRITE ( kunit_EFGALLIT1 ,'(i5,f14.8,i5)') natm,box,ntype
    WRITE ( kunit_EFGALLIT1 ,'(a10,i5)')  system,itime
    WRITE ( kunit_EFGALLIT1 ,'(a)') '      ia  type     vxx           vyy        vzz           eta           rx            ry            rz'
    do i = 1,natm
      WRITE ( kunit_EFGALLIT1 ,'(i8,2x,a3,7f14.8)') i,atype(i),nmrit(i,itja,1),nmrit(i,itja,2),nmrit(i,itja,3),nmrit(i,itja,4),rx(i),ry(i),rz(i)
    enddo
    itja=2
    eta_sum = 0.0D0
    vxx_sum = 0.0D0
    vyy_sum = 0.0D0
    vzz_sum = 0.0D0
    ! ===================================
    !  merge eta(i),vxx(i),vyy(i),vzz(i) 
    !  where i is the atom index
    ! ===================================
    eta = nmrit(:,itja,4)
    vxx = nmrit(:,itja,1)
    vyy = nmrit(:,itja,2)
    vzz = nmrit(:,itja,3)
    CALL MPI_ALLREDUCE(vxx,vxx_sum,natm,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
    CALL MPI_ALLREDUCE(vyy,vyy_sum,natm,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
    CALL MPI_ALLREDUCE(vzz,vzz_sum,natm,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
    CALL MPI_ALLREDUCE(eta,eta_sum,natm,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
    nmrit(:,itja,1) = vxx_sum
    nmrit(:,itja,2) = vyy_sum
    nmrit(:,itja,3) = vzz_sum
    nmrit(:,itja,4) = eta_sum
    WRITE ( kunit_EFGALLIT2 ,'(i5,f14.8,i5)') natm,box,ntype
    WRITE ( kunit_EFGALLIT2 ,'(a10,i5)')  system,itime
    WRITE ( kunit_EFGALLIT2 ,'(a)') '      ia  type     vxx           vyy        vzz           eta           rx            ry            rz'
    do i = 1,natm
      WRITE ( kunit_EFGALLIT2 ,'(i8,2x,a3,7f14.8)') i,atype(i),nmrit(i,itja,1),nmrit(i,itja,2),nmrit(i,itja,3),nmrit(i,itja,4),rx(i),ry(i),rz(i)
    enddo
  endif



  ! ================
  !  deallocation
  ! ================
  deallocate ( nmrit  ) 
  deallocate ( EFGTIT  )
  deallocate ( WIT    )
  deallocate ( eta , vxx , vyy , vzz )
  deallocate ( eta_sum , vxx_sum , vyy_sum , vzz_sum )
  deallocate ( nmr  ) 
  deallocate ( vzzmini , vzzmaxi )
  deallocate ( U )
  deallocate ( npvzzmin , npvzzmax )
  deallocate ( vzzm ) 
  deallocate ( vzzmit ) 
  deallocate ( vzzma )
  deallocate ( vzzmait )
  deallocate ( etam ) 
  deallocate ( etamit ) 
  deallocate ( pvzz )
  deallocate ( pvzzit )
  deallocate ( rho_z )
  deallocate ( rho_zit )
  deallocate ( sigmavzz )
  deallocate ( sigmavzzit )
  deallocate ( vzzsq )
  deallocate ( vzzsqit )

! ====================================
!  no accumulation for the firts step
! ====================================

  if(itime .eq. 0 .and. calc.ne.'efg' ) then
    dibvzztot = 0
    dibetatot = 0
    dibUtot   = 0
  endif

    ttt3 = MPI_WTIME(ierr)
    efgtimetot3 = efgtimetot3 + (ttt3-ttt2)

  return

100 FORMAT(I7,1X,' ATOM ',A3,'  minVZZ = ',F7.3,' maxVZZ = ',F7.3,' <VZZ> = ',F10.5,' <|VZZ|> =  ',F10.5,' mETA =  ',F10.5,' P(Vzz>0) =  ',F10.5, ' RHO_Z = ',F10.5)
200 FORMAT('atom = ',I6,' k = ',I12, 'value = ',F14.8,' min =  ',F7.3,' max =  ',F7.3,' vaa{a = xx,yy,zz}, eta = ',4F14.8)
210 FORMAT('atom = ',I6,' k = ',I12, 'value = ',F14.8,' min =  0.0D0    max =  1.0D0',' vaa{a = xx,yy,zz}, eta = ',4F14.8)
310 FORMAT('atom = ',I6,' k = ',I12, 'value = ',F14.8,' min =  ',F7.3,' max =  ',F7.3)


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
 
SUBROUTINE efg_ES ( itime , nefg )

  USE control,    ONLY :  myrank , numprocs, calc
  USE config,     ONLY :  system , natm , natmi , atype , atypei , box , omega , itype , rx , ry , rz , ntype , qia , qit
  USE io_file,    ONLY :  ionode , stdout , kunit_EFGFF , kunit_EFGALL
  USE constants,  ONLY :  pi , fpi , piroot , imag
  USE field,      ONLY :  qa , qb
  USE kspace,     ONLY :  struc_fact
  USE prop,       ONLY :  nprop_print , nprop
  USE time

  implicit none

  INCLUDE 'mpif.h'

  ! global
  integer, intent(in) :: itime , nefg 

  ! local
  integer :: i, ia, ja, k, nb, ierr , it , ui , na
  integer :: kvzz, keta
  integer :: nxij , nyij , nzij
  double precision :: d, d2, d4, d3 , d5 , expon, TMPW
  double precision :: alpha2, alpha3
  double precision :: allrealpart
  double precision :: T0 , T1 , T2  ! real part 
  double precision :: rxi , ryi , rzi , rxij , ryij , rzij 
  double precision :: cutoffsq
  double precision ::  W(3), EFGT(3,3)
  INTEGER, PARAMETER :: LWORK = 6
  double precision :: WORK(3 * LWORK)
  INTEGER :: IFAIL
 
  double precision, dimension (:,:), allocatable :: U
  double precision :: sq3 , sq32 , invbox 
  double precision :: ak, kx, ky, kz, kk
  integer :: ik
  double precision, dimension (:,:), allocatable :: nmr
  double precision, dimension (:),   allocatable :: npvzzmin , npvzzmax ! ntype 
  double precision, dimension (:), allocatable :: vzzmini , vzzmaxi ! ntype
  double precision, dimension (:),   allocatable :: vzzm, vzzma , etam , pvzz, rho_z, sigmavzz, vzzsq ! ntype
  double precision :: vzzk , etak 
  integer :: ku,uk
  double precision :: kri , onethird
  double precision :: ttt1 , ttt2 , ttt3 , ttt4 , ttt5
  double complex   :: rhon , carg , recarg , recarg_dgg 
  double precision, external :: errfc 
  integer :: ialpha, ibeta
!  double precision :: eself , tweatpi

  ttt1 = MPI_WTIME(ierr)

#ifdef debug
  CALL print_config_sample(0,0)
#endif
  ! ==========================
  !  allocate some quantities
  ! ==========================
  ! 1 = Vxx ; 2 = Vyy ; 3 = Vzz ; 4 = eta
  allocate ( nmr (natm , 4 ) )
  ! max min value ! not so important
  allocate ( vzzmini ( 0:ntype ), vzzmaxi ( 0:ntype ) )
  allocate ( npvzzmin (0:ntype) , npvzzmax (0:ntype) )
  ! Czjzek U component
  allocate ( U ( natm , 5 ) )
  ! vzz mean value (per type)
  allocate ( vzzm ( 0:ntype ) )
  ! mean value of absolute vzz (per type) 
  allocate ( vzzma ( 0:ntype ) )
  ! square of vzz (per type)
  allocate ( vzzsq ( 0:ntype ) )
  ! eta mean value (per type)
  allocate ( etam ( 0:ntype ) )
  ! proprtion of positiv vzz (per type)
  allocate ( pvzz ( 0:ntype ) )
  ! rho_z (per type)
  allocate ( rho_z ( 0:ntype ) )
  ! variance of vzz distribution (per type)
  allocate ( sigmavzz ( 0:ntype ) )

  ! ==============
  !  set to zero
  ! ==============
  nmr = 0.0D0
  pvzz = 0.0d0
  vzzmini  = 0.0D0
  vzzmaxi  = 0.0D0
  npvzzmin = 30000.0D0
  npvzzmax = -30000.0D0
  vzzmini ( : ) = 30000.0D0
  vzzmaxi ( : ) = -30000.0D0
  etam = 0.0d0
  vzzm = 0.0d0
  vzzma = 0.0d0
  vzzsq = 0.0d0
  U = 0.0D0
  efg_ia_real = 0.0d0
  efg_ia_dual = 0.0d0

  ! =================
  !  some constants 
  ! =================
  na = natmi(1)  
  nb = natm - na
  onethird = 1/3.0D0
  invbox = 1.0d0 / box
  sq3 = dsqrt( 3.0d0 )
  sq3 = 1.0d0 / sq3
  sq32 = sq3 * 0.5d0
  ! cutoff direct space part 
  cutoffsq = box * box * 0.25d0 
  ! related ewald parameter
  alpha2 = alphaES * alphaES      
  alpha3 = alpha2  * alphaES
!  tweatpi=2.0d0*sqrt(alphaES/pi)

  ! =====================
  ! facteur de structure 
  ! =====================
  CALL struc_fact ( km_efg )

  ttt2 = MPI_WTIME(ierr)
  efgtimetot2 = efgtimetot2 + (ttt2-ttt1) 

! ==============================================
!        direct space part
! ==============================================

   do ia = 1 , natm 
     rxi = rx(ia)
     ryi = ry(ia)
     rzi = rz(ia)
     do ja = 1, natm

       if(ja .ne. ia ) then

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
         !if( d2 .lt. cutoffsq ) then
           d = dsqrt( d2 )
           d3 = d2 * d
           d4 = d2 * d2
           d5 = d3 * d2
           expon = dexp( - alpha2 * d2 ) / piroot
           T0 = errfc( alphaES * d ) / d5
           T1 = ( 2.0d0 * alphaES * expon ) / d4 
           T2 = ( 4.0d0 * alpha3 * expon ) / d2 / 3.0d0
           allrealpart = qia(ja) * ( T0 + T1 + T2 )

           efg_ia_real(ia, 1, 1) = efg_ia_real(ia, 1, 1) - ( 3.0D0 * rxij * rxij - d2 ) * allrealpart
           efg_ia_real(ia, 2, 2) = efg_ia_real(ia, 2, 2) - ( 3.0D0 * ryij * ryij - d2 ) * allrealpart
           efg_ia_real(ia, 3, 3) = efg_ia_real(ia, 3, 3) - ( 3.0D0 * rzij * rzij - d2 ) * allrealpart
           efg_ia_real(ia, 1, 2) = efg_ia_real(ia, 1, 2) -   3.0d0 * rxij * ryij * allrealpart
           efg_ia_real(ia, 1, 3) = efg_ia_real(ia, 1, 3) -   3.0d0 * rxij * rzij * allrealpart
           efg_ia_real(ia, 2, 3) = efg_ia_real(ia, 2, 3) -   3.0d0 * ryij * rzij * allrealpart
         !endif

       endif
     enddo
  enddo
  efg_ia_real(:, 2, 1) = efg_ia_real(:, 1, 2) 
  efg_ia_real(:, 3, 1) = efg_ia_real(:, 1, 3) 
  efg_ia_real(:, 3, 2) = efg_ia_real(:, 2, 3)

#ifdef debug  
!  do i=1,natm
  i=1
  call print_tensor(efg_ia_real(i,:,:),' REAL')
!  enddo
#endif   

  ttt3 = MPI_WTIME(ierr)
  efgtimetot1 = efgtimetot1 + (ttt3-ttt2)

! ==============================================
!            reciprocal space part
! ==============================================
  efg_ia_dual = (0.0d0,0.0d0)
  ! sum on atoms for one givn ik
  do ia = 1 , natm
    rxi = rx(ia)
    ryi = ry(ia)
    rzi = rz(ia)
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
        rhon = rhon + qit(it) * CONJG( km_efg%strf ( ik , it ) )
      enddo
      kri = ( kx * rxi + ky * ryi + kz * rzi ) 
      carg = EXP ( imag * kri )
      recarg_dgg =  rhon*carg*ak / kk 
      recarg = recarg_dgg * kk 
      efg_ia_dual ( ia , 1, 1) = efg_ia_dual ( ia , 1 , 1) +  3.0d0 * kx * kx * recarg_dgg - recarg
      efg_ia_dual ( ia , 2, 2) = efg_ia_dual ( ia , 2 , 2) +  3.0d0 * ky * ky * recarg_dgg - recarg
      efg_ia_dual ( ia , 3, 3) = efg_ia_dual ( ia , 3 , 3) +  3.0d0 * kz * kz * recarg_dgg - recarg
      efg_ia_dual ( ia , 1, 2) = efg_ia_dual ( ia , 1 , 2) +  3.0d0 * kx * ky * recarg_dgg
      efg_ia_dual ( ia , 1, 3) = efg_ia_dual ( ia , 1 , 3) +  3.0d0 * kx * kz * recarg_dgg
      efg_ia_dual ( ia , 2, 3) = efg_ia_dual ( ia , 2 , 3) +  3.0d0 * ky * kz * recarg_dgg
    enddo kpoint
  enddo 

  efg_ia_dual ( : , 2, 1) = efg_ia_dual ( : , 1, 2)
  efg_ia_dual ( : , 3, 1) = efg_ia_dual ( : , 1, 3)
  efg_ia_dual ( : , 3, 2) = efg_ia_dual ( : , 2, 3)

  efg_ia_dual_real( : , :, :) = dble( efg_ia_dual( : , :, :)  ) 

  ! =======
  ! 4pi/V
  ! =======
  efg_ia_dual_real =  efg_ia_dual_real * fpi / omega / 3.0d0

#ifdef debug
!  do i=1,natm
i=1
  call print_tensor(efg_ia_dual_real(i,:,:),'RECIP')
!  enddo      
#endif
  
  ! ==============
  !  total tensor
  ! ==============
  efg_ia_total = efg_ia_dual_real + efg_ia_real 

#ifdef debug
!  do i=1,natm
i=1
  call print_tensor(efg_ia_total(i,:,:),'TOTAL')
!  enddo      
#endif

   ttt4 = MPI_WTIME(ierr)
   efgtimetot2 = efgtimetot2 + (ttt4-ttt3)
  
!=========================================================
!      STATISTICS on NMR parameters
!=========================================================
  atom: do ia = 1, natm  

    do ialpha = 1 , 3
      do ibeta = 1 , 3
        EFGT(ialpha,ibeta) = efg_ia_total ( ia , ialpha, ibeta) 
       enddo
    enddo
    EFGT = efg_ia_total ( ia , :, :)
    ! ===================================================================== 
    !  Czjzek components (see J. Phys.: Condens. Matter 10 (1998). p10719)
    ! =====================================================================
    U( ia, 1) = EFGT(3,3) * 0.5d0
    U( ia, 2) = EFGT(1,3) * sq3
    U( ia, 3) = EFGT(2,3) * sq3
    U( ia, 4) = EFGT(1,2) * sq3
    U( ia, 5) = ( EFGT(1,1) - EFGT(2,2) ) * sq32

    ! =================
    !  diagonalisation
    ! =================
    CALL DSYEV('N','U',3,EFGT,3,w,work,3 * lwork,ifail)
    if (ifail.ne.0) then
      if( ionode ) WRITE ( stdout , * ) 'ERROR: DSYEV, STOP in efg MODULE (ewald)'
      STOP 
    endif

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

    nmr(ia,3) = W(3) !ZZ
    nmr(ia,2) = W(1) !YY
    nmr(ia,1) = W(2) !XX
    nmr(ia,4) = (nmr(ia,2)-nmr(ia,1))/nmr(ia,3) !ETA

    ! ===================            
    !  test rearangement
    ! ===============================================================
    ! be aware that this is just a convention 
    ! because if all quantities are very small etaQ can be anything
    ! ===============================================================
    if(dabs(nmr(ia,3)).lt.dabs(nmr(ia,2))) then
      if( ionode ) WRITE ( stdout , * ) 'ERROR: |nmr| < |nmr|',nmr(ia,1),nmr(ia,2),nmr(ia,3)
      STOP 
    endif
    if(dabs(nmr(ia,3)).lt.dabs(nmr(ia,1))) then
      if( ionode ) WRITE ( stdout , * ) 'ERROR: |nmr| < |nmr|',nmr(ia,1),nmr(ia,2),nmr(ia,3)
      STOP 
    endif
    if(dabs(nmr(ia,1)).lt.dabs(nmr(ia,2))) then
      if( ionode ) WRITE ( stdout , * ) 'ERROR: |nmr| < |nmr|',nmr(ia,1),nmr(ia,2),nmr(ia,3)
      STOP 
    endif

    if( nmr(ia,1) .eq. nmr(ia,2) )  nmr(ia,4) = 0.0D0
    if( nmr(ia,3) .eq. 0.0 )  nmr(ia,4) = 0.0D0
    if( nmr(ia,1) .lt. 1.0E-08 .and. nmr(ia,2) .lt. 1.0E-8 .and. nmr(ia,3) .lt. 1.0E-8 ) nmr(ia,4) = 0.0D0
    if( nmr(ia,4) .gt. 1.0d0 .or. nmr(ia,4) .lt. 0.0d0) then
      if( ionode ) WRITE ( stdout ,'(a,4f48.24)') 'ERROR: eta > 1.0d0 or eta < 0.0d0',nmr(ia,4),nmr(ia,1),nmr(ia,2),nmr(ia,3)
      STOP 
    endif

    ! ===================================================
    !  stat for all sites index 0 in array of size ntype
    ! ===================================================

    vzzm (0)  =  vzzm (0) +       nmr(ia,3)
    etam (0)  =  etam (0) +       nmr(ia,4)
    vzzma(0)  =  vzzma(0) + dabs( nmr(ia,3) )
    if( nmr(ia,3) .le. vzzmini(0) )  vzzmini(0) = nmr(ia,3)
    if( nmr(ia,3) .ge. vzzmaxi(0) )  vzzmaxi(0) = nmr(ia,3)
    if( nmr(ia,3) .ge. 0.0d0 ) then
      pvzz(0) = pvzz(0) + 1.d0
    endif 

    vzzsq(0) = vzzsq(0) + nmr(ia,3) * nmr(ia,3)
        
    do it = 1 , ntype 
      if(itype(ia).eq.it) then
        vzzm (it) = vzzm(it)  + nmr(ia,3)
        vzzma(it) = vzzma(it) + dabs(nmr(ia,3))
        etam (it) = etam(it)  + nmr(ia,4)
        if(nmr(ia,3).le.vzzmini(it))  vzzmini(it) = nmr(ia,3)
        if(nmr(ia,3).ge.vzzmaxi(it))  vzzmaxi(it) = nmr(ia,3)
        vzzsq(it) = vzzsq(it) + nmr(ia,3) * nmr(ia,3)
        if(nmr(ia,3).ge.0.0d0 ) then 
          pvzz(it) = pvzz(it) + 1.d0
        endif
      endif ! it 
    enddo

  enddo atom

  if(ntype.eq.1) then
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
    if( ( ionode ) .and. (mod(nefg,nprop_print) .eq. 0) ) & 
    WRITE ( stdout ,100) itime , atypei(it) , vzzmini(it) , vzzmaxi(it) , vzzm(it) , vzzma(it) , etam(it) , pvzz(it) , rho_z ( it)
    if(   ionode ) WRITE ( kunit_EFGFF,100) &
    itime , atypei(it) , vzzmini(it) , vzzmaxi(it) , vzzm(it) , vzzma(it) , etam(it) , pvzz(it) , rho_z(it)
  enddo
  if( ionode .and. ntype .ne. 1 .and. (mod(nefg,nprop_print) .eq. 0) ) & 
  WRITE ( stdout ,100) itime , atypei(0) , vzzmini(0) , vzzmaxi(0) , vzzm(0) , vzzma(0) , etam(0) , pvzz(0) , rho_z(0)
  if( ionode .and. ntype .ne. 1 ) WRITE ( kunit_EFGFF ,100) &
  itime , atypei(0) , vzzmini(0) , vzzmaxi(0) , vzzm(0) , vzzma(0) , etam(0) , pvzz(0) , rho_z(0)

  do ia=1 , natm
    
    ! ====
    !  Ui
    ! ====
    do ui=1,5
      uk = (U(ia,ui)-umin)/resu
      ku = int(uk) + 1
      if(ku.lt.0.or.ku.gt.PANU) then
        if( ionode ) WRITE ( stdout , * ) 'ERROR: out of bound dibU1'
        if( ionode ) WRITE ( stdout ,310) i,ku,U(i,ui),umin,dabs(umin)
        STOP 
      endif
      dibUtot(ui,0,ku) = dibUtot(ui,0,ku) + 1
      do it=1,2
        if (itype(ia).eq.it) then
          dibUtot(ui,it,ku) = dibUtot(ui,it,ku) + 1
        endif
      enddo
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
    if( kvzz .lt. 0 .or. kvzz .gt. PANvzz ) then
      if( ionode .and. itype(ia) .eq. 1 ) WRITE ( stdout , * ) 'ERROR: out of bound distribvzz A'
      if( ionode .and. itype(ia) .eq. 2 ) WRITE ( stdout , * ) 'ERROR: out of bound distribvzz B'
      if( ionode ) WRITE ( stdout ,200) ia,kvzz,nmr(ia,3),npvzzmin(itype(ia)),npvzzmax(itype(ia)),nmr(ia,4),nmr(ia,1),nmr(ia,2),nmr(ia,3)
      STOP 
    endif
    if( keta .lt. 0 .or. keta .gt. PANeta ) then
      if( ionode ) WRITE ( stdout , * ) 'ERROR: out of bound distribeta'
      if( ionode ) WRITE ( stdout ,210) ia, keta, nmr(ia,4), nmr(ia,4), nmr(ia,1), nmr(ia,2), nmr(ia,3)
      STOP 
    endif
    dibvzztot(0,kvzz) = dibvzztot(0,kvzz) + 1
    dibetatot(0,keta) = dibetatot(0,keta) + 1
    do it=1,ntype
      if( itype(ia) .eq. it ) then
        dibvzztot(it,kvzz) = dibvzztot(it,kvzz) + 1
        dibetatot(it,keta) = dibetatot(it,keta) + 1
      endif
    enddo
  enddo  !ia 
!==========================
! END of distributions loop
!==========================
 
  ! =======================================
  ! write efg for each atom in file EFGALL 
  ! =======================================
  if( ionode  .and. lefgprintall .and. itime.ne.0) then
    WRITE ( kunit_EFGALL ,'(i5,f14.8,i5)') natm,box,ntype 
    WRITE ( kunit_EFGALL ,'(a10,i5)')  system,itime
    WRITE ( kunit_EFGALL ,'(a)') '      ia  type     vxx           vyy          vzz           eta           rx            ry            rz'
    do i = 1,natm
      WRITE ( kunit_EFGALL ,'(i8,2x,a3,7f14.8)') i,atype(i),nmr(i,1),nmr(i,2),nmr(i,3),nmr(i,4),rx(i),ry(i),rz(i)
    enddo
  endif

  ! ===============
  !  deallocation
  ! ===============
  deallocate ( nmr  )   
  deallocate ( vzzmini , vzzmaxi )
  deallocate ( U )
  deallocate ( npvzzmin , npvzzmax )
  deallocate ( vzzm ) 
  deallocate ( vzzma )
  deallocate ( etam ) 
  deallocate ( pvzz )
  deallocate ( rho_z )
  deallocate ( sigmavzz )
  deallocate ( vzzsq )

! ====================================
!  no accumulation for the first step
! ====================================
  if( itime .eq. 0 .and. calc.ne.'efg') then
    dibvzztot = 0
    dibetatot = 0
    dibUtot   = 0
  endif

  ! ===========
  ! total time
  ! ===========
  ttt5        = MPI_WTIME(ierr)
  efgtimetot3 = efgtimetot3 + (ttt5-ttt4)

  return

100 FORMAT(I7,1X,' ATOM ',A3,'  minVZZ = ',F7.3,' maxVZZ = ',F7.3,' <VZZ> =  ',F10.5,' <|VZZ|> =  ',F10.5,' mETA =  ',F10.5,' P(Vzz>0) =  ',F10.5, ' RHO_Z = ',F10.5)
200 FORMAT('atom = ',I6,' k = ',I12, 'value = ',F14.8,' min =  ',F7.3,' max =  ',F7.3,' vaa{a = xx,yy,zz}, eta = ',4F14.8)
210 FORMAT('atom = ',I6,' k = ',I12, 'value = ',F14.8,' min =  0.0D0    max =  1.0D0',' vaa{a = xx,yy,zz}, eta = ',4F14.8)
310 FORMAT('atom = ',I6,' k = ',I12, 'value = ',F14.8,' min =  ',F7.3,' max =  ',F7.3)


END SUBROUTINE efg_ES


!*********************** SUBROUTINE efg_write_output ********************************
!
! write distributions of EFG parameters and components to files
!
!******************************************************************************

SUBROUTINE efg_write_output ( iefgcount )

  USE io_file,          ONLY :  ionode , kunit_DTETAFF, kunit_DTVZZFF, kunit_DTIBUFF , kunit_DTETAFFIT, kunit_DTVZZFFIT
  USE config,           ONLY :  natm, natmi, ntype 
  USE constants,        ONLY :  dzero
  USE control,          ONLY :  longrange                
  USE time

  implicit none
 
  INCLUDE 'mpif.h'

  ! global
  integer, intent(in)  :: iefgcount

  ! local
  integer :: i , ierr , it , ui , itja
  integer, dimension(:,:),   allocatable :: dibvzztot_sum, dibetatot_sum
  integer, dimension(:,:,:),   allocatable :: dibUtot_sum
  double precision :: efgtime1 , efgtime2

  allocate( dibUtot_sum(6,0:ntype,0:PANU))
  allocate( dibvzztot_sum(0:ntype,0:PANvzz),dibetatot_sum(0:ntype,0:PANeta))

  OPEN (unit = kunit_DTETAFF,file = 'DTETAFF')
  OPEN (unit = kunit_DTVZZFF,file = 'DTVZZFF')
  OPEN (unit = kunit_DTIBUFF,file = 'DTIBUFF')
  if ( efg_it_contrib ) then
    OPEN (unit = kunit_DTETAFFIT,file = 'DTETAFFIT')
    OPEN (unit = kunit_DTVZZFFIT,file = 'DTVZZFFIT')
  endif

  efgtime1 = MPI_WTIME(ierr)

  dibetatot_sum = 0
  dibvzztot_sum = 0
  dibUtot_sum = 0  

  ! =============================================
  !  only the firect sumamtion is parallelized
  ! =============================================
  if ( longrange .eq. 'direct' ) then
    do it=0,ntype
      CALL MPI_ALLREDUCE(dibvzztot(it,:),dibvzztot_sum(it,:),PANvzz + 1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
      CALL MPI_ALLREDUCE(dibetatot(it,:),dibetatot_sum(it,:),PANeta + 1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
      do ui=1,6
        CALL MPI_ALLREDUCE(dibUtot(ui,it,:),dibUtot_sum(ui,it,:),PANU + 1,  MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
      enddo
    enddo 
    ! added 19-11-11
    if ( efg_it_contrib ) then
      do itja=1,ntype
        dibetatot_sum = 0
        dibvzztot_sum = 0
        do it=0,ntype
          CALL MPI_ALLREDUCE(dibvzztot_it(it,itja,:),dibvzztot_sum(it,:),PANvzz + 1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
          CALL MPI_ALLREDUCE(dibetatot_it(it,itja,:),dibetatot_sum(it,:),PANeta + 1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
        enddo
        dibvzztot_it(:,itja,:) = dibvzztot_sum
        dibetatot_it(:,itja,:) = dibetatot_sum
      enddo
    endif
  else
   do it=0,ntype
     dibvzztot_sum(it,:) = dibvzztot(it,:)
     dibetatot_sum(it,:) = dibetatot(it,:)
     do ui=1,6
       dibUtot_sum(ui,it,:) = dibUtot(ui,it,:) 
     enddo
   enddo
  endif

  if( ionode ) then
   ! ======================== 
   !  write eta distributions
   ! ======================== 
   WRITE (kunit_DTETAFF,'(a,4f15.8)') '#',reseta,dble(natmi(1)),dble(natmi(2)),dble(iefgcount)  
   if( ntype .eq. 1 ) WRITE (kunit_DTETAFF,'(4f15.8)') dzero,dzero,dzero
   if( ntype .eq. 2 ) WRITE (kunit_DTETAFF,'(4f15.8)') dzero,dzero,dzero,dzero
   do i = 0 , PANeta-1 
     WRITE (kunit_DTETAFF,'(4f15.8)') dble((i+1-0.5D0) * reseta ), ( dble(dibetatot_sum(it,i)) / (reseta * dble(natmi(it)) * dble(iefgcount)), it = 0 , ntype )  
   enddo
   ! ========================
   !  write Vzz distribution
   ! ========================
   do i = 0 , PANvzz 
     WRITE (kunit_DTVZZFF,'(4f15.8)') vzzmin + dble(i * resvzz), ( dble(dibvzztot_sum(it,i))/(resvzz * dble(natmi(it)) * dble(iefgcount)), it = 0 ,ntype ) 
   enddo

   ! =================================================
   ! write U1 and average Uk (with k>1) distribution 
   ! =================================================
   do i = 0 , PANU 
   WRITE (kunit_DTIBUFF,'(7f15.8)') umin + dble(i * resu), &
                              ( dble(dibUtot_sum(1,it,i)) / ( resu * dble(natmi(it)) * dble(iefgcount)), & 
                                dble(dibUtot_sum(6,it,i)) / ( resu * 4.0d0*dble(natmi(it)) * dble(iefgcount)), it = 0 ,ntype )
   enddo

 !added 19-11-11
  if ( efg_it_contrib ) then

   do itja = 1 , ntype 

     ! ======================== 
     !  write eta distributions IT
     ! ======================== 
     WRITE (kunit_DTETAFFIT,'(a,4f15.8)') '#',reseta,dble(natmi(1)),dble(natmi(2)),dble(iefgcount)
     WRITE (kunit_DTETAFFIT,'(4f15.8)') dzero,dzero,dzero,dzero
     do i = 0 , PANeta-1
       WRITE (kunit_DTETAFFIT,'(4f15.8)') dble((i+1-0.5D0) * reseta ), ( dble(dibetatot_it(it,itja,i)) / (reseta * dble(natmi(it)) * dble(iefgcount)), it = 0 , ntype )
     enddo
     ! ========================
     !  write Vzz distribution IT
     ! ========================
     do i = 0 , PANvzz
       WRITE (kunit_DTVZZFFIT,'(4f15.8)') vzzmin + dble(i * resvzz), ( dble(dibvzztot_it(it,itja,i))/(resvzz * dble(natmi(it)) * dble(iefgcount)), it = 0 ,ntype )
     enddo

   WRITE (kunit_DTETAFFIT,'(a)') ''
   WRITE (kunit_DTETAFFIT,'(a)') ''
   WRITE (kunit_DTVZZFFIT,'(a)') ''
   WRITE (kunit_DTVZZFFIT,'(a)') ''
   enddo
  endif ! it_contrib

  endif

  if (efg_it_contrib ) then
  CLOSE (kunit_DTVZZFFIT)
  CLOSE (kunit_DTETAFFIT)
  endif
  CLOSE (kunit_DTVZZFF)
  CLOSE (kunit_DTETAFF)
  CLOSE (kunit_DTIBUFF)

  deallocate( dibUtot_sum)
  deallocate( dibvzztot_sum,dibetatot_sum)

  efgtime2 = MPI_WTIME(ierr)
  efgtimetot3 = efgtimetot3 + (efgtime2-efgtime1)

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
  USE io_file,          ONLY  : ionode , stdout , kunit_EFGALL , kunit_EFGACFFF  

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

  if( ionode ) then
    WRITE ( stdout , '(a)'      )    'Remind some parameters of the system:'
    WRITE ( stdout , '(a,a12)'  )    'system    = ',system
    WRITE ( stdout , '(a,i12)'  )    'natm      = ',natm
    WRITE ( stdout , '(a,i12)'  )    'ntype     = ',ntype
    WRITE ( stdout , '(a,f12.5)')    'rho       = ',rho
    WRITE ( stdout , '(a,f12.5)')    'box       = ',box
    WRITE ( stdout , '(a,f12.5)')    'vol       = ',omega
    WRITE ( stdout , '(a)')          ''
  endif

  allocate ( vxxt(natm,ncefg) , vyyt(natm,ncefg)      , vzzt(natm,ncefg)        , etat(natm,ncefg)                                                            )
  allocate ( vxx0(natm)       , vyy0(natm)            , vzz0(natm)              , eta0(natm)                                                                  )
  allocate ( timeo (ncefg)    , norm(0:ntype,0:ntcor) , acfxx (0:ntype,0:ntcor) , acfyy (0:ntype,0:ntcor) , acfzz (0:ntype,0:ntcor) , acfet (0:ntype,0:ntcor) )
 
  acfxx = 0.0d0
  acfyy = 0.0d0
  acfzz = 0.0d0
  acfet = 0.0d0
  
  do t0 = 1 , ncefg 
    if( t0 .ne. 1 ) then
      READ(kunit_EFGALL,*) natm,box,ntype
      READ(kunit_EFGALL,*) system,iiii
      READ(kunit_EFGALL,*) XXXX
    endif
    do ia = 1 , natm
      READ(kunit_EFGALL,*) iiii , atype(ia) , vxxt( ia , t0 ) , vyyt( ia , t0 ) , vzzt( ia , t0 ), etat( ia , t0 ) , aaaa , aaaa , aaaa
    enddo
  enddo
  CLOSE(kunit_EFGALL)
 
  if( ionode ) WRITE ( stdout , '(a)' ) 'EFGALL successfully readed'

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
    if ( ntype .eq. 2 ) WRITE( kunit_EFGACFFF , '(i12,12f20.10)' ) t , acfxx(1,t) , acfyy(1,t) , acfzz(1,t) , acfet(1,t) , &
                                                                       acfxx(2,t) , acfyy(2,t) , acfzz(2,t) , acfet(2,t) , &
                                                                       acfxx(0,t) , acfyy(0,t) , acfzz(0,t) , acfet(0,t)
    if ( ntype .eq. 1 ) WRITE( kunit_EFGACFFF , '(i12,8f20.10)' )  t , acfxx(1,t) , acfyy(1,t) , acfzz(1,t) , acfet(1,t)
    endif
  enddo

  CLOSE(kunit_EFGACFFF)      

  deallocate(norm,acfxx,acfyy,acfzz,acfet)
  deallocate(vxxt,vyyt,vzzt,etat)
  deallocate(vxx0,vyy0,vzz0,eta0)

  return
        
END SUBROUTINE efg_acf

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
    if( nmr(1) .eq. nmr(2) )  nmr(4) = 0.0D0
    if( nmr(3) .eq. 0.0 )  nmr(4) = 0.0D0
    if( nmr(1) .lt. 1.0E-08 .and. nmr(2) .lt. 1.0E-8 .and. nmr(3) .lt. 1.0E-8 ) nmr(4) = 0.0D0
    if(dabs(nmr(3)).lt.dabs(nmr(2))) then
      if( ionode ) WRITE ( stdout , * ) 'ERROR: |Vzz| < |Vxx|',ia,nmr(1),nmr(2),nmr(3)
     STOP
    endif
    if(dabs(nmr(3)).lt.dabs(nmr(1))) then
      if( ionode ) WRITE ( stdout , * ) 'ERROR: |Vzz| < |Vyy|',ia,nmr(1),nmr(2),nmr(3)
      STOP
    endif
    if(dabs(nmr(1)).lt.dabs(nmr(2))) then
      if( ionode ) WRITE ( stdout , * ) 'ERROR: |Vyy| < |Vxx|',ia,nmr(1),nmr(2),nmr(3)
      STOP
    endif
    if( nmr(4) .gt. 1.0d0 .or. nmr(4) .lt. 0.0d0) then
      if( ionode ) WRITE ( stdout ,'(a,4f48.24)') 'ERROR: eta > 1.0d0 or eta < 0.0d0',ia,nmr(4),nmr(1),nmr(2),nmr(3)
      STOP
    endif

  return

END SUBROUTINE nmr_convention


END MODULE EFG
! ===== fmV =====
