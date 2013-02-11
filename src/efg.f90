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
!#define debug2
#define debug_multipole

!fix_grid in efg
!#define fix_grid

! ======= Hardware =======

!*********************** MODULE efg  **********************************
! two implementation of direct and Ewald sum are present for EFG calculation
! the old one ( lefg_old ) does only point charge
! the new one ( .not. lefg_old ) is based on the multipole expansion ( see
! multipole_DS and multipole_ES in field.f90
!**********************************************************************
MODULE efg 
 
  USE kspace
  USE rspace

  implicit none

  logical :: lefgprintall               ! print ( or not ) all the efg for each atoms and configs to file EFGALL
  logical :: lefg_it_contrib            ! if ones want to get the contributions to EFG separated in types
  logical :: lefg_restart               ! if EFGALL files are ready
  logical :: lefg_old                   ! use efg_DS and efg_ES ( old routines ) 
  logical :: lefg_stat                  ! compute statitics distribution on EFG's 
  integer :: ncefg                      ! number of configurations READ  for EFG calc (only when calc = 'efg')
  integer :: ntcor                      ! maximum number of steps for the acf calculation (calc = 'efg+acf')

#ifdef fix_grid
  double precision, dimension ( : , : ) , allocatable :: rgrid
#endif fix_grid

  double precision, dimension ( : , : ) , allocatable :: mu   ! elecctric dipoles

  ! ===================
  !  efg tensors 
  ! ===================
  double precision, dimension(:,:,:)  , allocatable :: efg_t         ! efg_tensor
  double precision, dimension(:,:,:,:), allocatable :: efg_t_it      ! efg_tensor
  double precision, dimension(:,:,:)  , allocatable :: efg_ia        ! efg_tensor
  double precision, dimension(:,:,:,:), allocatable :: efg_ia_it     ! efg_tensor

  ! ===============
  !  distributions
  ! ===============
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
 
  implicit none
  
  ! local
  integer :: ioerr
  character * 132 :: filename

  namelist /efgtag/  lefgprintall     , & 
                     lefg_it_contrib  , &
                     lefg_restart     , &
                     lefg_old         , &
                     lefg_stat        , &
                     ncefg            , & 
                     ntcor            , &
                     resvzz           , &
                     reseta           , &
                     resu             , &
                     vzzmin           , &
                     umin           

  if ( calc .ne. 'efg' .and. calc .ne.'efg+acf' ) return 

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
   if ( ionode ) WRITE ( stderr, '(a)') 'ERROR reading input_file : efgtag section is absent'
   STOP
  elseif ( ioerr .gt. 0 )  then
   if ( ionode ) WRITE ( stderr, '(a)') 'ERROR reading input_file : efgtag wrong tag'
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
  lefgprintall       = .true. 
  lefg_it_contrib    = .false.
  lefg_old           = .false.
  lefg_restart       = .false.
  lefg_stat          = .false.
  reseta             =   0.1D0
  resvzz             =   0.1D0
  resu               =   0.1D0 
  ncefg              =   0
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

  USE control,  ONLY :  calc , longrange
  USE io_file,  ONLY :  ionode , stderr
  USE config,   ONLY :  ntype

  implicit none

  if ( calc .eq. 'efg+acf' .and. ntype .gt. 2 ) then
    if ( ionode ) WRITE ( stderr , '(a)' ) 'ERROR the subroutine efg_acf is not implemented for ntype > 2'
    STOP 
  endif

  if ( calc .eq. 'efg+acf' ) return

  ! ==================================
  !  check vzzmin and Umin .ne. 0
  ! ==================================
  if ( vzzmin .eq. 0d0 .or. umin .eq. 0d0 ) then
    if ( ionode ) WRITE ( stderr ,'(a,2f8.3)') 'ERROR efgtag: vzzmin or umin should be set',vzzmin,umin   
    STOP
  endif             
  ! ==========================================
  !  set PAN (nb bins) from resolution values
  ! ==========================================
  PANeta = int(1.0d0/reseta)
  PANvzz = int((2.0d0*ABS (vzzmin))/resvzz)
  PANU = int((2.0D0*ABS (umin))/resu)

  return

END SUBROUTINE efg_check_tag

!*********************** SUBROUTINE efg_print_info ***************************
!
! print information to standard output about efg calculation 
!
!******************************************************************************

SUBROUTINE efg_print_info(kunit)

  USE io_file,  ONLY :  ionode 
  USE control,  ONLY :  calc , longrange , cutlongrange
  USE config,   ONLY :  ntype , atypei , natmi , simu_cell
  USE constants,ONLY :  tpi
  USE field    ,ONLY :  qch , alphaES , ncelldirect , kES

  implicit none

  ! global 
  integer :: kunit 

  if ( ionode ) then
    if ( calc .eq. 'efg' ) then
      WRITE ( kunit ,'(a)')                     '=============================================================' 
      WRITE ( kunit ,'(a)')                     ''
      WRITE ( kunit ,'(a)')                     'electric field gradient:'
      WRITE ( kunit ,'(a)')                     'point charges calculation'
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

!*********************** SUBROUTINE efgcalc ***********************************
!
! this SUBROUTINE initialize the calculation of efg when calc = 'efg' 
! from file TRAJFF.  
!
!******************************************************************************

SUBROUTINE efgcalc 

  USE io_file
  USE config,   ONLY :  system , natm , ntype , atype , rx , ry , rz , itype , & 
                        atypei , natmi, rho , simu_cell , config_alloc , qia, &
                        dipia , dipia_ind , dipia_wfc , ipolar , fx , fy , fz , phi_coul_tot , config_print_info 
  USE control,  ONLY :  longrange , myrank , numprocs
  USE field,    ONLY :  qch , dip , field_init , lpolar , lwfc , moment_from_pola , moment_from_wfc , rm_coul , km_coul , alphaES , multipole_DS , multipole_ES
  USE cell
  USE constants, ONLY : fpi
  USE thermodynamic, ONLY : u_coul_tot , vir_coul_tot , u_pol

  implicit none
  INCLUDE 'mpif.h'

  ! local
  integer :: iastart , iaend
  integer :: ia, iconf , na , it , it2 
  double precision :: aaaa 
  integer :: iiii , ierr
  character * 60 :: cccc 
  integer :: kunit
  integer :: kunit_it(2) ! <--  not general enough should be ntype dependent
  double precision :: ttt1 , ttt2
  double precision , dimension(:,:)   , allocatable :: ef_tmp
  double precision , dimension(:,:,:) , allocatable :: efg_tmp
 
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

    READ ( kunit_TRAJFF , * ) natm 
    READ ( kunit_TRAJFF , * ) system
    READ ( kunit_TRAJFF , * ) simu_cell%A ( 1 , 1 ) , simu_cell%A ( 2 , 1 ) , simu_cell%A ( 3 , 1 )
    READ ( kunit_TRAJFF , * ) simu_cell%A ( 1 , 2 ) , simu_cell%A ( 2 , 2 ) , simu_cell%A ( 3 , 2 )
    READ ( kunit_TRAJFF , * ) simu_cell%A ( 1 , 3 ) , simu_cell%A ( 2 , 3 ) , simu_cell%A ( 3 , 3 )
    READ ( kunit_TRAJFF , * ) ntype
    READ ( kunit_TRAJFF , * ) ( atypei ( it ) , it = 1 , ntype )
    IF ( ionode ) WRITE ( kunit_OUTFF ,'(A,20A3)' ) 'found type information on TRAJFF : ', atypei ( 1:ntype )
    IF ( ionode ) WRITE ( stdout      ,'(A,20A3)' ) 'found type information on TRAJFF : ', atypei ( 1:ntype )
    READ( kunit_TRAJFF ,*)   ( natmi ( it ) , it = 1 , ntype )

    allocate ( ef_tmp  ( natm , 3     ) )
    allocate ( efg_tmp ( natm , 3 , 3 ) )

    CALL lattice ( simu_cell ) 
    rho = natm / simu_cell%omega
    ! ===================================
    !  here we know natm, then alloc 
    !  and decomposition can be applied 
    ! ================================== 
    CALL config_alloc 
    CALL do_split ( natm , myrank , numprocs , iastart , iaend )
    ! read charge in fieldtag
    CALL field_init
    CALL efg_alloc
    CALL efg_mesh_alloc
    CALL typeinfo_init
    ! =============
    !  print info
    ! =============
    CALL efg_print_info(stdout)
    CALL efg_print_info(kunit_OUTFF)
  
    CALL config_print_info(stdout)
  
    ! =============================================
    !  not fix_grid: we calculate efg at 
    !  each atom positions which are moving.
    ! =============================================
  
    ! ============================================
    ! LOOP OVER CONFIGURATIONS 
    ! ============================================
    do iconf = 1, ncefg

      ttt1 = MPI_WTIME(ierr)

      na = 0
      if ( iconf .ne. 1 ) READ ( kunit_TRAJFF, * )  iiii 
      if ( iconf .ne. 1 ) READ ( kunit_TRAJFF, * )  cccc
      if ( iconf .ne. 1 ) READ ( kunit_TRAJFF, * )  aaaa,iiii 
      if ( iconf .ne. 1 ) READ ( kunit_TRAJFF , * ) ( cccc , it = 1 , ntype )
      if ( iconf .ne. 1 ) READ ( kunit_TRAJFF , * ) ( iiii , it = 1 , ntype )
      do ia = 1 , natm
        READ ( kunit_TRAJFF , * ) atype ( ia ) , rx ( ia ) , ry ( ia ) , rz ( ia ) ,aaaa,aaaa,aaaa,aaaa,aaaa,aaaa
      enddo

#ifdef debug
      ! ============================
      !  print minimum distance 
      !  and distance distribution
      ! ============================
      call distance_tab( kunit )
#endif

      ! =======================
      !  total tensor (efg_t)
      ! =======================
      efg_t = 0.0d0

      ! ======================================
      !     induced moment from polarisation 
      ! ======================================
      CALL moment_from_pola ( iastart , iaend , dipia_ind )

      ! ======================================
      !     induced moment from wannier center  
      ! ======================================
      CALL moment_from_WFc ( dipia_wfc ) 

      mu = dipia + dipia_ind + dipia_wfc

#ifdef debug_multipole
      if ( longrange .eq. 'ewald' )  CALL multipole_ES ( iastart, iaend , ef_tmp , efg_tmp , mu , u_coul_tot , vir_coul_tot , phi_coul_tot ) 
      if ( longrange .eq. 'direct' ) CALL multipole_DS ( iastart, iaend , ef_tmp , efg_tmp , mu , u_coul_tot , vir_coul_tot , phi_coul_tot ) 

      WRITE ( 10000, * ) '#electric field'
      WRITE ( 10001, * ) '#dipoles'
      WRITE (10002, * ) '#forces'
      WRITE (10003, * ) '#efg'
      do ia = 1 , natm
        WRITE (10000, '(3e16.8)' ) ef_tmp(ia,1), ef_tmp(ia,2) , ef_tmp(ia,3)
        WRITE (10001, '(3e16.8)' ) mu    (ia,1), mu    (ia,2) , mu    (ia,3)
        WRITE (10002, '(3e16.8)' ) fx    (ia), fy    (ia) , fz    (ia)
        WRITE (10003, '(6e16.8)' ) efg_tmp(ia,1,1),efg_tmp(ia,1,2),efg_tmp(ia,2,2),efg_tmp(ia,1,3),efg_tmp(ia,2,3),efg_tmp(ia,3,3)
      enddo
#endif

      ! =======================================
      !       efg charge only lefg_old=.true. 
      ! =======================================
      if ( lefg_old ) then 
        ! ===============
        ! use direct sum
        ! ===============
        if ( longrange .eq. 'direct' )  CALL efg_DS ( iastart , iaend , rm_coul )
        ! ===============
        ! use ewald sum
        ! ===============
        if ( longrange .eq. 'ewald' )   CALL efg_ES ( iastart , iaend , km_coul , alphaES )
      else
      ! =======================================
      !       efg charge + dipoles  
      ! =======================================
        ! ===============
        ! use direct sum
        ! ===============
        if ( longrange .eq. 'direct' )  CALL multipole_efg_DS ( iastart , iaend , rm_coul , mu )
        ! ===============
        ! use ewald sum
        ! ===============
        if ( longrange .eq. 'ewald' )   CALL multipole_efg_ES ( iastart , iaend , km_coul , alphaES , mu )
      endif

      efg_t    = efg_t + efg_ia 
      efg_t_it = efg_t_it + efg_ia_it 

#ifdef debug
      CALL print_tensor( efg_t( 1 , : , : ) , 'TOTEFG  ' )
#endif

      ! =======================================
      ! write efg for each atom in file EFGALL 
      ! =======================================
      if ( ionode  .and. lefgprintall ) then
        WRITE ( kunit_EFGALL , * )  natm
        WRITE ( kunit_EFGALL , * )  system
        WRITE ( kunit_EFGALL , * )  simu_cell%A(1,1) , simu_cell%A(2,1) , simu_cell%A(3,1)
        WRITE ( kunit_EFGALL , * )  simu_cell%A(1,2) , simu_cell%A(2,2) , simu_cell%A(3,2)
        WRITE ( kunit_EFGALL , * )  simu_cell%A(1,3) , simu_cell%A(2,3) , simu_cell%A(3,3)
        WRITE ( kunit_EFGALL , * )  ntype
        WRITE ( kunit_EFGALL , * )  ( atypei ( it ) , it = 1 , ntype )
        WRITE ( kunit_EFGALL , * )  ( natmi  ( it ) , it = 1 , ntype )
        WRITE ( kunit_EFGALL ,'(a)') &
        '      ia type                   vxx                   vyy                   vzz                   vxy                   vxz                   vyz'
        do ia = 1 , natm 
          it = itype ( ia ) 
          if ( lwfc( it ) .ge. 0 ) then 
            WRITE ( kunit_EFGALL ,'(i8,2x,a3,6e24.16)') ia , atype ( ia ) , efg_t ( ia , 1 , 1) , efg_t ( ia , 2 , 2) , &
                                                                            efg_t ( ia , 3 , 3) , efg_t ( ia , 1 , 2) , &
                                                                            efg_t ( ia , 1 , 3) , efg_t ( ia , 2 , 3)
          endif
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
            WRITE ( kunit , * )  simu_cell%A(1,1) , simu_cell%A(2,1) , simu_cell%A(3,1)
            WRITE ( kunit , * )  simu_cell%A(1,2) , simu_cell%A(2,2) , simu_cell%A(3,2)
            WRITE ( kunit , * )  simu_cell%A(1,3) , simu_cell%A(2,3) , simu_cell%A(3,3)
            WRITE ( kunit , * )  ntype
            WRITE ( kunit , * )  ( atypei ( it2 ) , it2 = 1 , ntype )
            WRITE ( kunit , * )  ( natmi  ( it2 ) , it2 = 1 , ntype )
            WRITE ( kunit ,'(a)') &
           '      ia type                   vxx                   vyy                   vzz                   vxy                   vxz                   vyz'
            do ia = 1 , natm
              it2 = itype ( ia ) 
              if ( lwfc( it2 ) .ge. 0 ) then 
                WRITE ( kunit ,'(i8,2x,a3,6e24.16)') &
                ia , atype ( ia ) , efg_t_it ( ia , it , 1 , 1) , efg_t_it ( ia , it , 2 , 2) , &
                                    efg_t_it ( ia , it , 3 , 3) , efg_t_it ( ia , it , 1 , 2) , &
                                    efg_t_it ( ia , it , 1 , 3) , efg_t_it ( ia , it , 2 , 3)
              endif
            enddo
          endif
        enddo
      endif

    ttt2 = MPI_WTIME(ierr)
    if ( ionode ) WRITE ( stdout , 110 ) 'config : ',iconf,' EFG  ', ttt2 - ttt1  

    enddo ! iconf loop

    CLOSE(kunit_EFGALL)

    if ( lefg_it_contrib ) then
      CLOSE( kunit_EFGALLIT1 )
      CLOSE( kunit_EFGALLIT2 )
    endif

  endif !lefg_restart

  CALL MPI_BARRIER( MPI_COMM_WORLD , ierr )

  ! ========================================================
  !     calculate statistical properties from EFGALL
  ! ========================================================
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

  CALL efg_dealloc
  CALL efg_mesh_dealloc

  return

110   FORMAT(2X,A8,I4,A20,' :  cpu time',F9.2)

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

SUBROUTINE efg_DS ( iastart , iaend , rm )

  USE control,  ONLY :  myrank , numprocs , calc , cutlongrange
  USE config,   ONLY :  system , natm , natmi , atype , atypei , itype , rx , ry , rz , ntype , qia , simu_cell
  USE prop,     ONLY :  nprop_print
  USE cell
  USE time
  USE io_file

  implicit none

  INCLUDE 'mpif.h'

  ! global
  integer, intent ( in )         :: iastart , iaend 
  TYPE ( rmesh ) , intent ( in ) :: rm

  ! local
  integer :: ia, ja, ierr , it 
  integer :: ncell
  double precision :: d , d2 , d5 , dm5
  double precision :: rxi , ryi , rzi , rxj , ryj , rzj , rxij , ryij , rzij , sxij , syij , szij 
  double precision :: cutefgsq
  double precision :: ttt1 , ttt2 , ttt3
  logical :: lcharge

  ! ===============================================
  !  check if there is any charge otherwise return
  ! ===============================================
  lcharge= .false.
  do ia = 1 , natm
    if ( qia ( ia ) .ne. 0.0d0 ) lcharge = .true.
  enddo
  if ( .not. lcharge ) return


#ifdef debug
  CALL print_config_sample(0,0)
#endif

  ttt1 = MPI_WTIME(ierr)

  cutefgsq = cutlongrange * cutlongrange 
  efg_ia(:,:,:) = 0.0d0
  efg_ia_it(:,:,:,:) = 0.0d0

! =========================================================
!  MAIN LOOP calculate EFG(i) for each atom i parallelized
! =========================================================

#ifdef debug
     WRITE ( stdout ,'(a,2i)')    'debug : iastart iaend',iastart , iaend
     WRITE ( stdout ,'(a,i)')     'debug : rm%ncmax',rm%ncmax
     WRITE ( stdout ,'(a,f16.5)') 'debug : cutefgsq ',cutefgsq
#endif

  ! ======================================
  !         cartesian to direct 
  ! ======================================
  CALL KARDIR(natm,rx,ry,rz,simu_cell%B)

atom : do ia = iastart , iaend
    rxi = rx(ia)
    ryi = ry(ia)
    rzi = rz(ia)
    ! ==================================================
    ! sum over neighboring cells (see direct_sum_init )
    ! ==================================================
    do ncell = 1 , rm%ncmax
         ! ==============================
         !  ia and ja in different cells
         ! ==============================
         if ( rm%lcell(ncell) .eq. 1) then

          do ja = 1 , natm
            rxj  = rx(ja) + rm%boxxyz(1,ncell)
            ryj  = ry(ja) + rm%boxxyz(2,ncell)
            rzj  = rz(ja) + rm%boxxyz(3,ncell)
            sxij = rxi - rxj  
            syij = ryi - ryj  
            szij = rzi - rzj  
            rxij = sxij * simu_cell%A(1,1) + syij * simu_cell%A(1,2) + szij * simu_cell%A(1,3)
            ryij = sxij * simu_cell%A(2,1) + syij * simu_cell%A(2,2) + szij * simu_cell%A(2,3)
            rzij = sxij * simu_cell%A(3,1) + syij * simu_cell%A(3,2) + szij * simu_cell%A(3,3)
            d2   = rxij * rxij + ryij * ryij + rzij * rzij
            if ( d2 .lt. cutefgsq ) then
              d = SQRT (d2)
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
        if ( rm%lcell(ncell) .eq. 0) then

          do ja = 1,natm

            if (ja.ne.ia) then
              rxj  = rx(ja)
              ryj  = ry(ja)
              rzj  = rz(ja)
              sxij = rxi - rxj
              syij = ryi - ryj
              szij = rzi - rzj
              rxij = sxij * simu_cell%A(1,1) + syij * simu_cell%A(1,2) + szij * simu_cell%A(1,3)
              ryij = sxij * simu_cell%A(2,1) + syij * simu_cell%A(2,2) + szij * simu_cell%A(2,3)
              rzij = sxij * simu_cell%A(3,1) + syij * simu_cell%A(3,2) + szij * simu_cell%A(3,3)
              d2   = rxij * rxij + ryij * ryij + rzij * rzij
              if (d2.lt.cutefgsq) then
                d = SQRT (d2)
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
  CALL print_tensor( efg_ia( 1    , : , : ) , 'EFG_1B  ' )
  CALL print_tensor( efg_ia( natm , : , : ) , 'EFG_NB  ' )
#endif

  ttt2 = MPI_WTIME(ierr)
  efgtimetot1 = efgtimetot1 + (ttt2-ttt1)

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
#ifdef debug
  CALL print_tensor( efg_ia( 1    , : , : ) , 'EFG_1A  ' )
  CALL print_tensor( efg_ia( natm , : , : ) , 'EFG_NA  ' )
#endif

  ttt3 = MPI_WTIME(ierr)
  efgtimetot3 = efgtimetot3 + (ttt3-ttt2)

  ! ======================================
  !         direct to cartesian
  ! ======================================
  CALL DIRKAR(natm,rx,ry,rz,simu_cell%A)


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
 
SUBROUTINE efg_ES ( iastart , iaend , km , alphaES )

  USE control,    ONLY :  myrank , numprocs, calc
  USE config,     ONLY :  system , natm , natmi , atype , atypei , &
                          simu_cell , itype , rx , ry , rz , ntype , qia 
  USE constants,  ONLY :  pi , fpi , piroot , imag
  USE field,      ONLY :  qch
  USE kspace,     ONLY :  struc_fact
  USE prop,       ONLY :  nprop_print , nprop
  USE cell 
  USE time
  USE io_file 

  implicit none

  INCLUDE 'mpif.h'

  ! global
  integer, intent(in) :: iastart , iaend 
  TYPE ( kmesh ) :: km
  double precision :: alphaES

  ! local
  integer :: ia , ja ,  ierr , it 
  double precision :: d , d2 , d4 , d3 , d5 , expon 
  double precision :: alpha2 , alpha3
  double precision :: allrealpart
  double precision :: T0 , T1 , T2  ! real part 
  double precision :: rxi , ryi , rzi , rxij , ryij , rzij , sxij , syij , szij
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
  logical :: lcharge 

  ! ===============================================
  !  check if there is any charge otherwise return
  ! ===============================================
  lcharge= .false.
  do ia = 1 , natm
    if ( qia ( ia ) .ne. 0.0d0 ) lcharge = .true.
  enddo
  if ( .not. lcharge ) return


#ifdef debug
  CALL print_config_sample(0,0)
  WRITE ( stdout ,'(a,2i)')    'debug in efg_ES : iastart iaend',iastart , iaend
  WRITE ( stdout ,'(a,i)')     'debug in efg_ES : km%ncmax',km%nkcut
  WRITE ( stdout ,'(a,f16.5)') 'debug in efg_ES : alphaES ',alphaES
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
  alpha2 = alphaES * alphaES      
  alpha3 = alpha2  * alphaES

  ttt1 = MPI_WTIME(ierr)

  ! =====================
  ! facteur de structure 
  ! =====================
  CALL struc_fact ( km )

  ! ======================================
  !         cartesian to direct 
  ! ======================================
  CALL KARDIR(natm,rx,ry,rz,simu_cell%B)

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
         sxij = rxij - nint ( rxij )
         syij = ryij - nint ( ryij )
         szij = rzij - nint ( rzij )
         rxij = sxij * simu_cell%A(1,1) + syij * simu_cell%A(1,2) + szij * simu_cell%A(1,3)
         ryij = sxij * simu_cell%A(2,1) + syij * simu_cell%A(2,2) + szij * simu_cell%A(2,3)
         rzij = sxij * simu_cell%A(3,1) + syij * simu_cell%A(3,2) + szij * simu_cell%A(3,3)
         d2 = rxij * rxij + ryij * ryij + rzij * rzij
         d = SQRT ( d2 )
         d3 = d2 * d
         d4 = d2 * d2
         d5 = d3 * d2
         expon = EXP ( - alpha2 * d2 ) / piroot
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

  ! ======================================
  !         cartesian to direct 
  ! ======================================
  CALL DIRKAR(natm,rx,ry,rz,simu_cell%A)

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
kpoint: do ik = 1, km%nkcut 
      ! =================
      !   k-space  
      ! =================
      kx = km%kpt ( 1 , ik )
      ky = km%kpt ( 2 , ik )
      kz = km%kpt ( 3 , ik )
      kk = km%kptk( ik )
      Ak = EXP ( - kk * 0.25d0 / alpha2 ) 
      if ( km%kptk( ik ) .eq. 0 ) then 
        WRITE ( stdout , * ) 'the sum should be done on k! =  0',ik
        STOP 
      endif

      ! ===============================
      !                              ---
      !  charge density in k-space ( \   q * facteur de structure  )
      !                              /__
      ! ===============================
      rhon = (0.d0, 0.d0)
!      do it = 1, ntype
      do ja = 1, natm
!        rhon = rhon + qch(it) * CONJG( km%strf ( ik , it ) )
        rhon = rhon + qia(ja) * CONJG( km%strf ( ik , ja ) )
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
!          rhon = qch(it) * CONJG( km%strf ( ik , it ) )
          rhon = qch(it) * CONJG( km%strf ( ik , it ) )

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
  efg_ia_dual_real( : , : , :) = DBLE ( efg_ia_dual( : , : , :)  ) 
  efg_ia_dual_real =  efg_ia_dual_real * fpi / simu_cell%omega / 3.0d0

  if ( lefg_it_contrib ) then
    ! =======
    ! 4pi/3V
    ! =======
    efg_ia_dual_real_it( : , : , : , :) = DBLE ( efg_ia_dual_it( : , : , : , :)  ) 
    efg_ia_dual_real_it =  efg_ia_dual_real_it * fpi / simu_cell%omega / 3.0d0
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

#ifdef debug2 
  CALL print_tensor( efg_ia_real  ( 1 , : , : ) , 'EFG_DIRO' )
  CALL print_tensor( efg_ia_dual_real  ( 1 , : , : ) , 'EFG_RECO' )
  CALL print_tensor( efg_ia   ( 1 , : , : ) , 'EFG_TOTO' )
#endif


  return

END SUBROUTINE efg_ES


SUBROUTINE multipole_efg_DS ( iastart, iaend , rm , mu )

  USE config,  ONLY : natm , atype , natmi , ntype , qia , rx , ry , rz , fx , fy , fz , tau_coul , simu_cell 
  USE control, ONLY : cutlongrange , myrank
  USE io_file, ONLY : stdout , ionode 
  USE field,   ONLY :  qch
  USE cell
  USE time

  implicit none

  INCLUDE 'mpif.h'

  ! global 
  integer, intent(in) :: iastart , iaend
  TYPE ( rmesh ) :: rm
  double precision :: mu    ( natm , 3 )

  ! local 
  integer :: ia, ja , ierr , ncell 
  double precision :: cutsq
  double precision :: rxi , ryi , rzi 
  double precision :: rxj , ryj , rzj
  double precision :: rxij , ryij , rzij
  double precision :: sxij , syij , szij
  double precision :: qj 
  double precision :: mujx , mujy , mujz
  double precision :: Txx , Tyy , Tzz , Txy , Txz , Tyz
  double precision :: Txxx,  Tyyy,  Tzzz, Txxy, Txxz, Tyyx, Tyyz, Tzzx, Tzzy, Txyz
  double precision :: d , d2 
  double precision :: dm5 , dm7 

  double precision :: ttt1 , ttt2 , ttt3
  logical          :: lcentralbox

  ttt1 = MPI_WTIME(ierr)

  ! =============================== 
  !         some constants
  ! =============================== 
  cutsq = cutlongrange * cutlongrange

#ifdef debug
  if ( ionode ) then
  WRITE( stdout , '(a)')    'debug: in multipole_efg_DS'
  WRITE( stdout , '(a,i8,a)') 'debug : rm ',rm%ncmax,rm%meshlabel
  do ia = 1 , natm
    WRITE( stdout , '(a,f12.5)')  'debug : charge (atom)  ',qia(ia)
  enddo
  do ia = 1 , natm
    WRITE( stdout , '(a,3f12.5)') 'debug : dipole (atom)  ', mu ( ia , 1 ) , mu ( ia , 2 ) , mu ( ia , 3 )
  enddo
  WRITE( stdout , '(a,2i8)')     'debug : iastart iaend',iastart ,iaend
  WRITE( stdout , '(a,f20.5)')   'debug : cutsq ',cutsq
  endif
  call print_config_sample(0,0)
#endif  

  efg_ia   = 0.0d0

  ! ======================================
  !         cartesian to direct 
  ! ======================================
  CALL KARDIR(natm,rx,ry,rz,simu_cell%B)
  
! ==================================================================================================
!  MAIN LOOP calculate EFG(i)  for each atom i parallelized
! ==================================================================================================
  atom : do ia = iastart , iaend

    rxi  = rx  ( ia )
    ryi  = ry  ( ia )
    rzi  = rz  ( ia )

    ! ==================================================
    ! sum over neighboring cells (see direct_sum_init )
    ! ==================================================
    cells: do ncell = 1 , rm%ncmax

      lcentralbox=rm%lcell(ncell).eq.0   
                   
      do ja = 1 , natm

        if ( ( .not.lcentralbox ) .or. ( lcentralbox .and. ja .ne. ia )  ) then
 
            rxj  = rx ( ja ) + rm%boxxyz( 1 , ncell )
            ryj  = ry ( ja ) + rm%boxxyz( 2 , ncell )
            rzj  = rz ( ja ) + rm%boxxyz( 3 , ncell )
            sxij = rxi - rxj
            syij = ryi - ryj
            szij = rzi - rzj
            rxij = sxij * simu_cell%A(1,1) + syij * simu_cell%A(1,2) + szij * simu_cell%A(1,3)
            ryij = sxij * simu_cell%A(2,1) + syij * simu_cell%A(2,2) + szij * simu_cell%A(2,3)
            rzij = sxij * simu_cell%A(3,1) + syij * simu_cell%A(3,2) + szij * simu_cell%A(3,3)
            d2    = rxij * rxij + ryij * ryij + rzij * rzij

            if ( d2 .lt. cutsq .and. d2 .ne. 0.0d0 ) then 

              qj    = qia ( ja )
              mujx  = mu ( ja , 1 )  
              mujy  = mu ( ja , 2 )  
              mujz  = mu ( ja , 3 )  
              d     = SQRT ( d2 )
              dm5   = 1.0d0 / ( d2 * d2 * d )
              dm7   = dm5 / d2 * 3.0d0
 
              ! multipole interaction tensor rank = 2
              Txx = ( 3.0d0 * rxij * rxij - d2 ) * dm5
              Tyy = ( 3.0d0 * ryij * ryij - d2 ) * dm5
              Tzz = ( 3.0d0 * rzij * rzij - d2 ) * dm5
              Txy = ( 3.0d0 * rxij * ryij      ) * dm5
              Txz = ( 3.0d0 * rxij * rzij      ) * dm5
              Tyz = ( 3.0d0 * ryij * rzij      ) * dm5

              ! multipole interaction tensor rank = 3  
              Txxx = ( 5.0d0 * rxij * rxij * rxij -  3.0d0 * d2 * ( rxij ) ) * dm7 
              Tyyy = ( 5.0d0 * ryij * ryij * ryij -  3.0d0 * d2 * ( ryij ) ) * dm7 
              Tzzz = ( 5.0d0 * rzij * rzij * rzij -  3.0d0 * d2 * ( rzij ) ) * dm7 
              Txxy = ( 5.0d0 * rxij * rxij * ryij -          d2 * ( ryij ) ) * dm7 
              Txxz = ( 5.0d0 * rxij * rxij * rzij -          d2 * ( rzij ) ) * dm7 
              Tyyx = ( 5.0d0 * ryij * ryij * rxij -          d2 * ( rxij ) ) * dm7 
              Tyyz = ( 5.0d0 * ryij * ryij * rzij -          d2 * ( rzij ) ) * dm7 
              Tzzx = ( 5.0d0 * rzij * rzij * rxij -          d2 * ( rxij ) ) * dm7 
              Tzzy = ( 5.0d0 * rzij * rzij * ryij -          d2 * ( ryij ) ) * dm7 
              Txyz = ( 5.0d0 * rxij * ryij * rzij                          ) * dm7 

              ! ===========================================================
              !                  charge-charge interaction
              ! ===========================================================

              ! electric field gradient
              efg_ia ( ia , 1 , 1 )  = efg_ia ( ia , 1 , 1 ) - qj * Txx 
              efg_ia ( ia , 2 , 2 )  = efg_ia ( ia , 2 , 2 ) - qj * Tyy 
              efg_ia ( ia , 3 , 3 )  = efg_ia ( ia , 3 , 3 ) - qj * Tzz 
              efg_ia ( ia , 1 , 2 )  = efg_ia ( ia , 1 , 2 ) - qj * Txy 
              efg_ia ( ia , 1 , 3 )  = efg_ia ( ia , 1 , 3 ) - qj * Txz  
              efg_ia ( ia , 2 , 3 )  = efg_ia ( ia , 2 , 3 ) - qj * Tyz

              ! ===========================================================
              !                  dipole-dipole interaction
              ! ===========================================================

              ! electric field gradient
              efg_ia ( ia , 1 , 1 ) = efg_ia ( ia , 1 , 1 ) - ( Txxx * mujx + Txxy * mujy + Txxz * mujz ) 
              efg_ia ( ia , 2 , 2 ) = efg_ia ( ia , 2 , 2 ) - ( Tyyx * mujx + Tyyy * mujy + Tyyz * mujz ) 
              efg_ia ( ia , 3 , 3 ) = efg_ia ( ia , 3 , 3 ) - ( Tzzx * mujx + Tzzy * mujy + Tzzz * mujz ) 
              efg_ia ( ia , 1 , 2 ) = efg_ia ( ia , 1 , 2 ) - ( Txxy * mujx + Tyyx * mujy + Txyz * mujz ) 
              efg_ia ( ia , 1 , 3 ) = efg_ia ( ia , 1 , 3 ) - ( Txxz * mujx + Txyz * mujy + Tzzx * mujz ) 
              efg_ia ( ia , 2 , 3 ) = efg_ia ( ia , 2 , 3 ) - ( Txyz * mujx + Tyyz * mujy + Tzzy * mujz ) 

            endif

        endif

      enddo

    enddo cells

  enddo atom 

  ttt2 = MPI_WTIME(ierr)
  efgtimetot1 = efgtimetot1 + ( ttt2 - ttt1 )

  CALL MPI_ALL_REDUCE_DOUBLE ( efg_ia ( : , 1 , 1 ) , natm )
  CALL MPI_ALL_REDUCE_DOUBLE ( efg_ia ( : , 2 , 2 ) , natm )
  CALL MPI_ALL_REDUCE_DOUBLE ( efg_ia ( : , 3 , 3 ) , natm )
  CALL MPI_ALL_REDUCE_DOUBLE ( efg_ia ( : , 1 , 2 ) , natm )
  CALL MPI_ALL_REDUCE_DOUBLE ( efg_ia ( : , 1 , 3 ) , natm )
  CALL MPI_ALL_REDUCE_DOUBLE ( efg_ia ( : , 2 , 3 ) , natm )

  ! EFG is symmetric
  ! not needed ... just for consistency
  efg_ia ( : , 2 , 1 ) = efg_ia ( : , 1 , 2 )
  efg_ia ( : , 3 , 1 ) = efg_ia ( : , 1 , 3 )
  efg_ia ( : , 3 , 2 ) = efg_ia ( : , 2 , 3 )

  ttt3 = MPI_WTIME(ierr)
  efgtimetot3 = efgtimetot3 + ( ttt3 - ttt2 )

  ! ======================================
  !         cartesian to direct 
  ! ======================================
  CALL DIRKAR(natm,rx,ry,rz,simu_cell%A)

  return

END SUBROUTINE multipole_efg_DS



SUBROUTINE multipole_efg_ES ( iastart, iaend , km , alphaES , mu )

  USE control,   ONLY : lsurf ,cutoff
  USE config,    ONLY : natm , ntype , natmi , atype , rx , ry , rz , fx , fy , fz , qia , simu_cell 
  USE constants, ONLY : imag , pi , piroot , tpi , fpi
  USE io_file  , ONLY : ionode , stdout 
  USE field,     ONLY : qch
  USE cell
  USE time

  implicit none

  INCLUDE 'mpif.h'

  ! global 
  integer, intent(in) :: iastart , iaend
  TYPE ( kmesh ) :: km
  double precision :: alphaES
  double precision :: mu    ( natm , 3 )

  ! local 
  integer             :: ia , ja , ik , it , ip , ierr
  double precision, dimension(:,:,:), allocatable :: efg_dir , efg_rec , efg_self

  double precision :: mip    ( ntype , 3 )

  double precision :: rxi  , ryi  , rzi
  double precision :: rxj  , ryj  , rzj
  double precision :: kx   , ky   , kz 
  double precision :: rxij , ryij , rzij
  double precision :: sxij , syij , szij
  double precision :: kk   , kri  , Ak 
  double precision :: qj  
  double precision :: mujx , mujy , mujz
  double complex   :: rhon , carg 
  double precision :: recarg 
  double precision :: expon , F1 , F2 , F3 
  double precision :: k_dot_mu
  double precision :: Txx , Tyy , Tzz , Txy , Txz , Tyz
  double precision :: Txxx,  Tyyy,  Tzzz, Txxy, Txxz, Tyyx, Tyyz, Tzzx, Tzzy, Txyz
  double precision :: d , d2 , d3  , d5 
  double precision :: dm1 , dm5 , dm7 
  double precision :: alpha2 , alpha3 , alpha5 
  double precision :: selfa 
  double precision, external :: errfc
  double precision :: fpi_V 
  double precision :: ttt1 , ttt2  , ttt3 , ttt4 
  double precision :: ttt1p , ttt2p 


  ! =============================
  !  init mip ( itype dependent ) 
  ! ==============================
  ip = 0
  do it = 1, ntype
    ip = ip + natmi ( it )
    mip ( it , 1 ) = mu ( ip , 1 )
    mip ( it , 2 ) = mu ( ip , 2 )
    mip ( it , 3 ) = mu ( ip , 3 )
  enddo

#ifdef debug
  if ( ionode ) then
    WRITE( stdout , '(a)')          'debug : in multipole_efg_ES'
    WRITE( stdout , '(a,i8,a,a)')   'debug : km        ',km%nkcut,' ',km%meshlabel
    do ia = 1 , natm
      WRITE( stdout , '(a,f12.5)')  'debug : charge (atom)  ', qia(ia)
    enddo
    do it = 1 , ntype
      WRITE( stdout , '(a,f12.5)')  'debug : charge (type)  ', qch(it)
    enddo
    do ia = 1 , natm
      WRITE( stdout , '(a,3f12.5)') 'debug : dipole (atom)  ', mu ( ia , 1 ) , mu ( ia , 2 ) , mu ( ia , 3 )
    enddo
    do it = 1 , ntype
      WRITE( stdout , '(a,3f12.5)') 'debug : dipole (type)  ', mip ( it , 1 ) , mip ( it , 2 ) , mip ( it , 3 )
    enddo
    WRITE( stdout , '(a,2i8)')      'debug : iastart iaend  ', iastart ,iaend
    WRITE( stdout , '(a,f20.5)')    'debug : alphaES        ', alphaES
  endif
  call print_config_sample(0,0)
#endif 

  allocate( efg_dir ( natm , 3 , 3 ) , efg_rec ( natm , 3 , 3 ) , efg_self ( natm , 3 , 3 ) ) 

  efg_dir  = 0.0d0
  efg_rec  = 0.0d0
  efg_self = 0.0d0

  ttt1p = MPI_WTIME(ierr)
  ! =====================
  ! facteur de structure 
  ! =====================
  CALL struc_fact ( km )

  ttt2p = MPI_WTIME(ierr)
  strftimetot = strftimetot + ( ttt2p - ttt1p ) 

  ! =================
  !  some constants 
  ! =================
  fpi_V  = fpi / simu_cell%omega  ! 4pi / V
  alpha2 = alphaES * alphaES
  alpha3 = alpha2  * alphaES
  alpha5 = alpha3  * alpha2
  selfa  = - 4.0d0 * alpha3 / 3.0d0 / piroot

  ttt1 = MPI_WTIME(ierr)

  ! ======================================
  !         cartesian to direct 
  ! ======================================
  CALL KARDIR(natm,rx,ry,rz,simu_cell%B)

  ! ==============================================
  !        direct space part
  ! ==============================================

  do ia = iastart , iaend
    rxi = rx(ia)
    ryi = ry(ia)
    rzi = rz(ia)

    do ja = 1, natm

      if (ja .ne. ia ) then

        qj   = qia(ja)
        mujx = mu ( ja , 1 )
        mujy = mu ( ja , 2 )
        mujz = mu ( ja , 3 )
        rxj  = rx(ja)
        ryj  = ry(ja)
        rzj  = rz(ja)
        rxij = rxi - rxj
        ryij = ryi - ryj
        rzij = rzi - rzj
        sxij = rxij - nint ( rxij )
        syij = ryij - nint ( ryij )
        szij = rzij - nint ( rzij )
        rxij = sxij * simu_cell%A(1,1) + syij * simu_cell%A(1,2) + szij * simu_cell%A(1,3)
        ryij = sxij * simu_cell%A(2,1) + syij * simu_cell%A(2,2) + szij * simu_cell%A(2,3)
        rzij = sxij * simu_cell%A(3,1) + syij * simu_cell%A(3,2) + szij * simu_cell%A(3,3)
        d2  = rxij * rxij + ryij * ryij + rzij * rzij
        d   = SQRT ( d2 )
        d3  = d2 * d
        d5  = d3 * d2
        dm1 = 1.0d0 / d
        !dm3 = dm1 / d2
        dm5 = 1.0d0 / d2 / d3 
        dm7 = dm5 / d2 * 3.0d0

        expon = EXP ( - alpha2 * d2 )    / piroot
        F1    = errfc( alphaES * d ) + 2.0d0 * alphaES * d  * expon
        F2    = F1 + 4.0d0 * alpha3  * d3 * expon / 3.0d0
        F3    = F2 + 8.0d0 * alpha5  * d5 * expon / 15.0d0

        ! multipole interaction tensor rank = 2
        Txx = ( 3.0d0 * rxij * rxij * F2 - d2 * F1 ) * dm5  
        Tyy = ( 3.0d0 * ryij * ryij * F2 - d2 * F1 ) * dm5 
        Tzz = ( 3.0d0 * rzij * rzij * F2 - d2 * F1 ) * dm5
        Txy = ( 3.0d0 * rxij * ryij * F2           ) * dm5
        Txz = ( 3.0d0 * rxij * rzij * F2           ) * dm5 
        Tyz = ( 3.0d0 * ryij * rzij * F2           ) * dm5 
 
        ! multipole interaction tensor rank = 3  
        Txxx = ( -5.0d0 * rxij * rxij * rxij * F3 +  3.0d0 * d2 * ( rxij ) * F2 ) * dm7 
        Tyyy = ( -5.0d0 * ryij * ryij * ryij * F3 +  3.0d0 * d2 * ( ryij ) * F2 ) * dm7 
        Tzzz = ( -5.0d0 * rzij * rzij * rzij * F3 +  3.0d0 * d2 * ( rzij ) * F2 ) * dm7
        Txxy = ( -5.0d0 * rxij * rxij * ryij * F3 +          d2 * ( ryij ) * F2 ) * dm7 
        Txxz = ( -5.0d0 * rxij * rxij * rzij * F3 +          d2 * ( rzij ) * F2 ) * dm7 
        Tyyx = ( -5.0d0 * ryij * ryij * rxij * F3 +          d2 * ( rxij ) * F2 ) * dm7 
        Tyyz = ( -5.0d0 * ryij * ryij * rzij * F3 +          d2 * ( rzij ) * F2 ) * dm7 
        Tzzx = ( -5.0d0 * rzij * rzij * rxij * F3 +          d2 * ( rxij ) * F2 ) * dm7 
        Tzzy = ( -5.0d0 * rzij * rzij * ryij * F3 +          d2 * ( ryij ) * F2 ) * dm7 
        Txyz = ( -5.0d0 * rxij * ryij * rzij * F3                               ) * dm7 
        ! ===========================================================
        !                  charge-charge interaction
        ! ===========================================================

        ! electric field gradient 
        efg_dir ( ia , 1 , 1 ) = efg_dir( ia , 1 , 1 ) - qj * Txx 
        efg_dir ( ia , 2 , 2 ) = efg_dir( ia , 2 , 2 ) - qj * Tyy
        efg_dir ( ia , 3 , 3 ) = efg_dir( ia , 3 , 3 ) - qj * Tzz
        efg_dir ( ia , 1 , 2 ) = efg_dir( ia , 1 , 2 ) - qj * Txy
        efg_dir ( ia , 1 , 3 ) = efg_dir( ia , 1 , 3 ) - qj * Txz
        efg_dir ( ia , 2 , 3 ) = efg_dir( ia , 2 , 3 ) - qj * Tyz

        ! ===========================================================
        !                  dipole-dipole interaction
        ! ===========================================================
        
        ! electric field gradient 
        efg_dir ( ia , 1 , 1 ) = efg_dir ( ia , 1 , 1 ) + ( Txxx * mujx + Txxy * mujy + Txxz * mujz ) 
        efg_dir ( ia , 2 , 2 ) = efg_dir ( ia , 2 , 2 ) + ( Tyyx * mujx + Tyyy * mujy + Tyyz * mujz )
        efg_dir ( ia , 3 , 3 ) = efg_dir ( ia , 3 , 3 ) + ( Tzzx * mujx + Tzzy * mujy + Tzzz * mujz )
        efg_dir ( ia , 1 , 2 ) = efg_dir ( ia , 1 , 2 ) + ( Txxy * mujx + Tyyx * mujy + Txyz * mujz ) 
        efg_dir ( ia , 1 , 3 ) = efg_dir ( ia , 1 , 3 ) + ( Txxz * mujx + Txyz * mujy + Tzzx * mujz )
        efg_dir ( ia , 2 , 3 ) = efg_dir ( ia , 2 , 3 ) + ( Txyz * mujx + Tyyz * mujy + Tzzy * mujz )
 
      endif

    enddo

  enddo 

  ttt2 = MPI_WTIME(ierr)
  efgtimetot1 = efgtimetot1 + ( ttt2 - ttt1 )

  ! ======================================
  !         direct to cartesian
  ! ======================================
  CALL DIRKAR(natm,rx,ry,rz,simu_cell%A)

  ! ==============================================
  !            reciprocal space part
  ! ==============================================
  kpoint : do ik = 1, km%nkcut
    ! =================
    !   k-space  
    ! =================
    kx   = km%kpt(1,ik)
    ky   = km%kpt(2,ik)
    kz   = km%kpt(3,ik)
    kk   = km%kptk(ik)
    Ak   = EXP ( - kk * 0.25d0 / alpha2 ) / kk

    if (km%kptk(ik) .eq. 0 ) then
      WRITE ( stdout , * ) 'the sum should be done on k! =  0',ik
      STOP
    endif
    ! ===============================
    !                              ---
    !  charge density in k-space ( \   q * facteur de structure  )
    !                              /__
    ! ===============================
    rhon   = (0.d0, 0.d0)
    do ja = 1, natm
      k_dot_mu = ( mu ( ja , 1 ) * kx + mu ( ja , 2 ) * ky + mu ( ja , 3 ) * kz  ) 
      rhon = rhon + ( qia(ja) + imag * k_dot_mu ) * km%strf ( ik , ja ) 
    enddo

    do ia = 1 , natm

      rxi = rx(ia)
      ryi = ry(ia)
      rzi = rz(ia)
      kri = ( kx * rxi + ky * ryi + kz * rzi )
      carg   = EXP  ( imag * kri )
      recarg = DBLE ( CONJG ( rhon ) * carg  * Ak )

      ! elecctric field gradient
      efg_rec ( ia , 1 , 1 ) = efg_rec ( ia , 1 , 1 ) + kx * kx * recarg
      efg_rec ( ia , 2 , 2 ) = efg_rec ( ia , 2 , 2 ) + ky * ky * recarg 
      efg_rec ( ia , 3 , 3 ) = efg_rec ( ia , 3 , 3 ) + kz * kz * recarg 
      efg_rec ( ia , 1 , 2 ) = efg_rec ( ia , 1 , 2 ) + kx * ky * recarg 
      efg_rec ( ia , 1 , 3 ) = efg_rec ( ia , 1 , 3 ) + kx * kz * recarg 
      efg_rec ( ia , 2 , 3 ) = efg_rec ( ia , 2 , 3 ) + ky * kz * recarg 
 
    enddo

  enddo kpoint

  ttt3 = MPI_WTIME(ierr)
  efgtimetot2 = efgtimetot2 + ( ttt3 - ttt2 )

  ! ====================================================== 
  !             merge real part 
  ! ====================================================== 
  CALL MPI_ALL_REDUCE_DOUBLE ( efg_dir ( : , 1 , 1 ) , natm )
  CALL MPI_ALL_REDUCE_DOUBLE ( efg_dir ( : , 2 , 2 ) , natm )
  CALL MPI_ALL_REDUCE_DOUBLE ( efg_dir ( : , 3 , 3 ) , natm )
  CALL MPI_ALL_REDUCE_DOUBLE ( efg_dir ( : , 1 , 2 ) , natm )
  CALL MPI_ALL_REDUCE_DOUBLE ( efg_dir ( : , 1 , 3 ) , natm )
  CALL MPI_ALL_REDUCE_DOUBLE ( efg_dir ( : , 2 , 3 ) , natm )

  ! ======================================================
  ! remark on the unit :
  ! 1/(4*pi*epislon_0) = 1 => epsilon_0 = 1/4pi
  ! ======================================================
  efg_rec =   efg_rec * fpi_V  

  ! ====================================================== 
  !              Self contribution 
  ! ====================================================== 

  ! field gradient 
  do ia = 1 , natm
   efg_self ( ia , 1 , 1 ) = selfa * qia ( ia ) 
   efg_self ( ia , 2 , 2 ) = selfa * qia ( ia ) 
   efg_self ( ia , 3 , 3 ) = selfa * qia ( ia ) 
  enddo

  ! =====================================================
  !                     TOTAL
  ! =====================================================

  efg_ia = ( efg_dir + efg_rec  + efg_self )
  ! only for consistency 
  efg_ia ( : , 2 , 1 ) = efg_ia ( : , 1 , 2 )
  efg_ia ( : , 3 , 1 ) = efg_ia ( : , 1 , 3 )
  efg_ia ( : , 3 , 2 ) = efg_ia ( : , 2 , 3 )

#ifdef debug2
  CALL print_tensor( efg_dir  ( 1 , : , : ) , 'EFG_DIRN' )
  CALL print_tensor( efg_dir  ( 2 , : , : ) , 'EFG_DIRN' )
  CALL print_tensor( efg_rec  ( 1 , : , : ) , 'EFG_RECN' )
  CALL print_tensor( efg_rec  ( 2 , : , : ) , 'EFG_RECN' )
  CALL print_tensor( efg_self ( 1 , : , : ) , 'EFG_SELN' )
  CALL print_tensor( efg_self ( 2 , : , : ) , 'EFG_SELN' )
  CALL print_tensor( efg_ia   ( 1 , : , : ) , 'EFG_TOTN' )
  CALL print_tensor( efg_ia   ( 2 , : , : ) , 'EFG_TOTN' )
#endif

  deallocate( efg_dir , efg_rec , efg_self ) 

  ttt4 = MPI_WTIME (ierr)
  efgtimetot3 = efgtimetot3 + ( ttt4 - ttt3 )

  return

END SUBROUTINE multipole_efg_ES

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
   WRITE (kunit_eta,'(a,<ntype+2>f15.8)') '#',reseta,(DBLE (natmi(it)),it=0,ntype )   
   WRITE (kunit_eta,'(<ntype+2>f15.8)')  dzero,(dzero,it=0,ntype)
   do i = 0 , PANeta-1 
     WRITE (kunit_eta,'(<ntype+2>f15.8)') &
     DBLE ((i+1-0.5D0) * reseta ), &
   ( DBLE (dibetatot(it,i)) / (reseta * DBLE (natmi(it)) * DBLE (ncefg)), it = 0 , ntype )  
   enddo
   ! ========================
   !  write Vzz distribution
   ! ========================
   do i = 0 , PANvzz 
     WRITE (kunit_vzz ,'(<ntype+2>f15.8)') &
     vzzmin + DBLE (i * resvzz), &
   ( DBLE (dibvzztot(it,i))/(resvzz * DBLE (natmi(it)) * DBLE (ncefg)), it = 0 ,ntype ) 
   enddo

   ! =================================================
   ! write U1 and average Uk (with k>1) distribution 
   ! =================================================
   do i = 0 , PANU 
   WRITE (kunit_u ,'(<2*ntype+3>f15.8)') umin + DBLE (i * resu), &
                              ( DBLE (dibUtot(1,it,i)) / ( resu * DBLE (natmi(it)) * DBLE (ncefg)), & 
                                DBLE (dibUtot(6,it,i)) / ( resu * 4.0d0*DBLE (natmi(it)) * DBLE (ncefg)), &
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

  USE config,           ONLY  : system , natm , ntype , itype , atype , simu_cell , rho , config_alloc
  USE io_file,          ONLY  : ionode , stdout , kunit_EFGALL , kunit_EFGACFFF  , kunit_OUTFF
  USE cell

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
  READ(kunit_EFGALL,*) natm
  READ(kunit_EFGALL,*) system
  READ(kunit_EFGALL,*) simu_cell%A(1,1) , simu_cell%A(2,1) , simu_cell%A(3,1)
  READ(kunit_EFGALL,*) simu_cell%A(1,2) , simu_cell%A(2,2) , simu_cell%A(3,2)
  READ(kunit_EFGALL,*) simu_cell%A(1,3) , simu_cell%A(2,3) , simu_cell%A(3,3)
  READ(kunit_EFGALL,*) ntype
  READ(kunit_EFGALL,*) XXXX
  CALL lattice ( simu_cell ) 
  rho = natm / simu_cell%omega
  CALL config_alloc 

  CALL print_general_info ( stdout ) 
  CALL print_general_info ( kunit_OUTFF ) 

  allocate ( vxxt(natm,ncefg) , vyyt(natm,ncefg)      , vzzt(natm,ncefg)        , etat(natm,ncefg)    )
  allocate ( vxx0(natm)       , vyy0(natm)            , vzz0(natm)              , eta0(natm)          )
  allocate ( timeo (ncefg)    , norm(0:ntype,0:ntcor) , &
             acfxx (0:ntype,0:ntcor) , acfyy (0:ntype,0:ntcor) , &
             acfzz (0:ntype,0:ntcor) , acfet (0:ntype,0:ntcor) )
 
  acfxx = 0.0d0
  acfyy = 0.0d0
  acfzz = 0.0d0
  acfet = 0.0d0
  
  do t0 = 1 , ncefg 
    if ( t0 .ne. 1 ) then
      READ(kunit_EFGALL,*) natm
      READ(kunit_EFGALL,*) system
      READ(kunit_EFGALL,*) simu_cell%A(1,1) , simu_cell%A(2,1) , simu_cell%A(3,1)
      READ(kunit_EFGALL,*) simu_cell%A(1,2) , simu_cell%A(2,2) , simu_cell%A(3,2)
      READ(kunit_EFGALL,*) simu_cell%A(1,3) , simu_cell%A(2,3) , simu_cell%A(3,3)
      READ(kunit_EFGALL,*) ntype
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
    t0max = MIN ( ncefg , t0 +ntcor )
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
                                atype , atypei , simu_cell , rho, config_alloc
  USE io_file,          ONLY  : ionode , stdout , stderr , kunit_OUTFF 
  USE prop,             ONLY : nprop_print
  USE field,            ONLY : lwfc         
  USE cell
  USE constants

  implicit none

  integer            :: i , ic , ia , it , ui
  integer, parameter :: lwork = 6
  integer            :: ifail
  integer            :: kvzz, keta ,ku
  integer            :: kunit_input , kunit_output , kunit_nmroutput

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
  sq3 = SQRT ( 3.0d0 )
  sq3 = 1.0d0 / sq3
  sq32 = sq3 * 0.5d0

  READ ( kunit_input , * )  natm
  READ ( kunit_input , * )  system
  READ ( kunit_input , * )  simu_cell%A(1,1) , simu_cell%A(2,1) , simu_cell%A(3,1)
  READ ( kunit_input , * )  simu_cell%A(1,2) , simu_cell%A(2,2) , simu_cell%A(3,2)
  READ ( kunit_input , * )  simu_cell%A(1,3) , simu_cell%A(2,3) , simu_cell%A(3,3)
  READ ( kunit_input , * )  ntype
  READ ( kunit_input ,* ) ( atypei ( it ) , it = 1 , ntype )
  IF ( ionode ) WRITE ( kunit_OUTFF ,'(A,20A3)' ) 'found type information on EFGALL : ', atypei ( 1:ntype )
  IF ( ionode ) WRITE ( stdout      ,'(A,20A3)' ) 'found type information on EFGALL : ', atypei ( 1:ntype )
  READ( kunit_input ,*)   ( natmi ( it ) , it = 1 , ntype )
  READ ( kunit_input , * )  xxxx

  CALL lattice ( simu_cell )
  rho = natm / simu_cell%omega
  ! ===================================
  !  here we know natm, then alloc 
  !  and decomposition can be applied 
  ! ================================== 
  if ( lefg_restart ) CALL config_alloc 
  if ( lefg_restart ) CALL efg_alloc

  CALL typeinfo_init

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
      READ ( kunit_input , * )  simu_cell%A(1,1) , simu_cell%A(2,1) , simu_cell%A(3,1)
      READ ( kunit_input , * )  simu_cell%A(1,2) , simu_cell%A(2,2) , simu_cell%A(3,2)
      READ ( kunit_input , * )  simu_cell%A(1,3) , simu_cell%A(2,3) , simu_cell%A(3,3)
      READ ( kunit_input , * )  ntype
      READ ( kunit_input , * )  ( xxxx , it = 1 , ntype )
      READ ( kunit_input , * )  ( iiii  , it = 1 , ntype )
      READ ( kunit_input , * )  xxxx
    endif

    do ia = 1 , natm
      it = itype ( ia ) 
      if ( lwfc ( it ) .ge. 0 ) then 
      READ( kunit_input , * ) iiii , xxxx , efgt(1,1) , efgt(2,2) , efgt(3,3) , &
                                            efgt(1,2) , efgt(1,3) , efgt(2,3)
      endif

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
        if ( ionode ) WRITE ( stderr , * ) 'ERROR: DSYEV, STOP in efg MODULE (ewald)'
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
      vzzma ( 0 )  =  vzzma( 0 ) + ABS ( nmr ( ia , 3 ) )
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
          vzzma( it ) = vzzma( it )  + ABS ( nmr (ia , 3 ) )
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
      vzzm  = vzzm  / DBLE ( natm )
      vzzsq = vzzsq / DBLE ( natm )
      vzzma = vzzma / DBLE ( natm )
      etam  = etam  / DBLE ( natm )
      pvzz  = pvzz  / DBLE ( natm )
    else
      vzzm  ( 0 )  = vzzm  ( 0 ) / DBLE ( natm )
      vzzsq ( 0 )  = vzzsq ( 0 ) / DBLE ( natm )
      vzzma ( 0 )  = vzzma ( 0 ) / DBLE ( natm )
      etam  ( 0 )  = etam  ( 0 ) / DBLE ( natm )
      pvzz  ( 0 )  = pvzz  ( 0 ) / DBLE ( natm )
      do it=1,ntype
        vzzm ( it ) = vzzm ( it ) / DBLE ( natmi ( it ) )
        vzzsq( it ) = vzzsq( it ) / DBLE ( natmi ( it ) )
        vzzma( it ) = vzzma( it ) / DBLE ( natmi ( it ) )
        etam ( it ) = etam ( it ) / DBLE ( natmi ( it ) )
        pvzz ( it ) = pvzz ( it ) / DBLE ( natmi ( it ) )
      enddo
    endif

    do it =0 , ntype
      sigmavzz ( it ) = vzzsq ( it ) - ( vzzma ( it ) * vzzma ( it ) )
      sigmavzz ( it ) = SQRT ( sigmavzz ( it ) )
      rho_z( it ) = sigmavzz ( it ) / vzzma ( it )
    enddo

    ! ========================================
    !  output average values ( instantaneous )
    ! ========================================
    do it=1,ntype 
      if ( ( ionode ) .and. (MOD (ic,nprop_print) .eq. 0) ) &
      WRITE ( stdout ,100) & 
      ic , atypei(it) , vzzmini(it) , vzzmaxi(it) , vzzm(it) , vzzma(it) , etam(it) , pvzz(it) , rho_z ( it)
      if (   ionode ) WRITE ( kunit_output , 100) &
      ic , atypei(it) , vzzmini(it) , vzzmaxi(it) , vzzm(it) , vzzma(it) , etam(it) , pvzz(it) , rho_z(it)
    enddo
    if ( ionode .and. ntype .ne. 1 .and. (MOD (ic,nprop_print) .eq. 0) ) &
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
        it = itype ( ia )   
        if ( lwfc ( it ) .ge. 0 ) then
        WRITE ( kunit_nmroutput , 150 ) ia , atype ( ia ) , & 
                                        nmr( ia , 1 ) , nmr( ia , 2 ) , &
                                        nmr( ia , 3 ) , nmr( ia , 4 )
        endif
      enddo
    endif

   if ( .not. lefg_stat ) then
     if ( ionode ) WRITE ( stdout , '(a)' ) "No statistic on EFG's"
     return
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
          if ( ionode ) WRITE ( stderr , * ) 'ERROR: out of bound dibU1'
          if ( ionode ) WRITE ( stderr ,310) i,ku,U(i,ui),umin,ABS (umin)
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
        WRITE ( stderr , '(a)' ) 'ERROR: out of bound distribvzz A'
        if ( ionode .and. itype(ia) .eq. 2 ) &
        WRITE ( stderr , '(a)' ) 'ERROR: out of bound distribvzz B'
        if ( ionode ) WRITE ( stderr ,200) &
        ia , kvzz , nmr ( ia , 3 ) , vzzmini( itype ( ia ) ) , vzzmaxi ( itype ( ia ) ) , &
                    nmr ( ia , 4 ) , nmr ( ia , 1 ) , nmr ( ia , 2 ) , nmr ( ia , 3 )
        STOP
      endif
      if ( keta .lt. 0 .or. keta .gt. PANeta ) then
        if ( ionode ) WRITE ( stderr , '(a)' ) 'ERROR: out of bound distribeta'
        if ( ionode ) WRITE ( stderr ,210) ia, keta, nmr(ia,4), nmr(ia,4), &
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

SUBROUTINE nmr_convention( EIG , nmr , ia )

  USE io_file,  ONLY :  stdout , stderr , ionode

  implicit none

  double precision ::  EIG(3) , del11 , del22 , del33 , diso
  double precision ::  nmr(4)
  integer :: ia

  ! =======================================================
  !  NMR convention: |Vzz| > |Vyy| > |Vxx|
  ! =======================================================
  diso=(EIG(1)+EIG(2)+EIG(3))/3.0d0
  if ( ABS ( diso ) .ne. 0.0d0 .and. ABS ( diso ) .gt. 1e-4 ) then
    if ( ionode ) WRITE ( stderr ,'(a,i6,f48.24)') 'ERROR: trace of EFG is not null',ia,diso
     diso = 0.0d0
  !  STOP
  endif

  ! NMR Convention 1 :
  ! ZZ >= YY >= XX  
  ! ZZ = DEL11
  ! YY = DEL22
  ! XX = DEL33

  IF ( ( ABS ( EIG(1) - diso ) .GE. ABS ( EIG(2) - diso ) ) .AND. ( ABS ( EIG(1) - diso ) .GE. ABS ( EIG(3) - diso ) ) ) THEN
    DEL11=EIG(1)
    IF ( ABS ( EIG(2) - diso ) .GE. ABS ( EIG(3) - diso ) ) THEN
      DEL22=EIG(2)
      DEL33=EIG(3)
    ELSE
      DEL22=EIG(3)
      DEL33=EIG(2)
    ENDIF
  ELSEIF ( ( ABS ( EIG(2) - diso ) .GE. ABS ( EIG(1) - diso ) ) .AND. ( ABS ( EIG(2) - diso ) .GE. ABS ( EIG(3) - diso ) ) )  THEN
    DEL11=EIG(2)
    IF ( ABS ( EIG(1) - diso ) .GE. ABS (EIG(3)-diso) ) THEN
      DEL22=EIG(1)
      DEL33=EIG(3)
    ELSE
      DEL22=EIG(3)
      DEL33=EIG(1)
    ENDIF
  ELSEIF (( ABS ( EIG(3) - diso ) .GE. ABS ( EIG(1) - diso ) ) .AND. ( ABS ( EIG(3) - diso ) .GE. ABS ( EIG(2) - diso ) ) ) THEN
    DEL11=EIG(3)
    IF ( ABS ( EIG(1) - diso ) .GE. ABS ( EIG(2) - diso ) ) THEN
      DEL22=EIG(1)
      DEL33=EIG(2)
    ELSE
      DEL22=EIG(2)
      DEL33=EIG(1)
    ENDIF
  ELSE
    WRITE (0,*) 'INTERNAL ERROR DIAGSHIFT, STOP'
    STOP
  ENDIF
  EIG(1) = DEL33 ! XX
  EIG(2) = DEL22 ! YY
  EIG(3) = DEL11 ! ZZ

  nmr(1) = EIG(1)                          ! XX
  nmr(2) = EIG(2)                          ! YY
  nmr(3) = EIG(3)                          ! ZZ

  nmr(4) = ( nmr(1) - nmr(2) ) / ( nmr(3) )   ! ETA
  ! ===================
  !  test rearangement
  ! ===================
  if ( nmr(1) .eq. nmr(2) )  nmr(4) = 0.0D0
  if ( nmr(3) .eq. 0.0 )  nmr(4) = 0.0D0
  if ( ABS ( nmr(1) ) .lt. 1.0E-05 .and. ABS ( nmr(2) ) .lt. 1.0E-5 .and. ABS ( nmr(3) ) .lt. 1.0E-5 ) nmr(4) = 0.0D0
  if (ABS (nmr(3)- diso ).lt.ABS (nmr(1) -diso)) then
    if ( ionode ) WRITE ( stderr , '(a,i4,3e24.16)' ) 'ERROR: |Vzz| < |Vxx|',ia,nmr(1),nmr(2),nmr(3)
    STOP
  endif
  if (ABS (nmr(3) -diso ).lt.ABS (nmr(2) -diso )) then
    if ( ionode ) WRITE ( stderr , '(a,i4,3e24.16)' ) 'ERROR: |Vzz| < |Vyy|',ia,nmr(1),nmr(2),nmr(3)
    STOP
  endif
  if (ABS (nmr(2) -diso ).lt.ABS (nmr(1)-diso)) then
    if ( ionode ) WRITE ( stderr , '(a,i4,3e24.16)' ) 'ERROR: |Vyy| < |Vxx|',ia,nmr(1),nmr(2),nmr(3)
    STOP
  endif
  if ( nmr(4) .gt. 1.0d0 ) then
    if ( ionode ) &
    WRITE ( stderr ,'(a,i6,5f24.16)') 'ERROR1: eta > 1.0d0 ',ia,nmr(4),ABS(nmr(1)-diso),ABS(nmr(2)-diso),ABS(nmr(3)-diso)
    WRITE ( stderr ,'(a,i6,5f24.16)') 'ERROR2: eta > 1.0d0 ',ia,nmr(4),nmr(1),nmr(2),nmr(3)
    !STOP
  !  nmr(4) = 1.0d0
  endif
  if ( nmr(4) .gt. 1.0d0 .or. nmr(4) .lt. 0.0d0) then
    if ( ionode ) &
    WRITE ( stderr ,'(a,i6,4f24.16)') 'ERROR: eta < 0.0d0',ia,nmr(4),nmr(1),nmr(2),nmr(3)
  !  nmr(4) = 0.0d0
    !STOP
  endif

  return

END SUBROUTINE nmr_convention


!*********************** SUBROUTINE efg_alloc *********************************
!
! allocate quantities related to efg calculation
! /deallocate 
!
!******************************************************************************

SUBROUTINE efg_alloc 

  USE control,  ONLY :  calc , longrange
  USE config,   ONLY :  natm, ntype , qia , itype

  implicit none

  if ( calc .ne. 'efg' ) return

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
  allocate( mu ( natm , 3 ) ) 
  mu = 0.0d0
  efg_t = 0.0d0
  efg_t_it = 0.0d0
  efg_ia = 0.0d0
  efg_ia_it = 0.0d0

  return

END SUBROUTINE efg_alloc

SUBROUTINE efg_dealloc

  USE control,  ONLY :  calc , longrange

  implicit none

  if ( calc .ne. 'efg' ) return

  deallocate( dibUtot )
  deallocate( dibvzztot )
  deallocate( dibetatot )
#ifdef fix_grid
  deallocate ( rgrid )
#endif fix_grid
  deallocate( efg_t )
  deallocate( efg_t_it )
  deallocate( efg_ia    )
  deallocate( efg_ia_it )
  deallocate( mu ) 

  return

END SUBROUTINE efg_dealloc


!*********************** SUBROUTINE efg_alloc *********************************
!
! allocate quantities related to efg calculation
! /deallocate 
!
!******************************************************************************

SUBROUTINE efg_mesh_alloc

  USE control,  ONLY :  calc , longrange
  USE config,   ONLY :  natm, ntype , qia , itype
  USE field,    ONLY :  km_coul ,rm_coul , ncelldirect , kES

  implicit none

  ! local
  integer :: nkcut
  integer :: ncmax

  if ( calc .ne. 'efg' ) return

  ! ============
  !  direct sum
  ! ============
  if ( longrange .eq. 'direct' ) then
    rm_coul%meshlabel='rm_efg'
    ncmax = ( 2 * ncelldirect + 1 ) ** 3
    rm_coul%ncmax=ncmax
    rm_coul%ncell=ncelldirect
    allocate( rm_coul%boxxyz( 3 , ncmax ) , rm_coul%lcell( ncmax ) , rm_coul%rr( ncmax ) )
    CALL direct_sum_init ( rm_coul )
  endif

  ! ============
  !  ewald sum
  ! ============
  if ( longrange .eq. 'ewald')  then
    km_coul%meshlabel='km_coul'
    km_coul%kmax(1) = kES(1)
    km_coul%kmax(2) = kES(2)
    km_coul%kmax(3) = kES(3)
    nkcut = ( 2 * km_coul%kmax(1) + 1 ) * ( 2 * km_coul%kmax(2) + 1 ) * ( 2 * km_coul%kmax(3) + 1 )   
    nkcut = nkcut - 1
    km_coul%nkcut = nkcut
    allocate( km_coul%kptk( nkcut ) , km_coul%kpt(3,nkcut) )
!    allocate ( km_coul%strf ( nkcut, ntype) )
    allocate ( km_coul%strf ( nkcut, natm) )
    CALL kpoint_sum_init ( km_coul )
  endif

  return

END SUBROUTINE efg_mesh_alloc

SUBROUTINE efg_mesh_dealloc

  USE control,  ONLY :  calc , longrange
  USE field,    ONLY :  km_coul , rm_coul 

  implicit none

  if ( calc .ne. 'efg' ) return

  if ( longrange .eq. 'direct' ) then
    deallocate( rm_coul%boxxyz , rm_coul%lcell )
  endif

  if ( longrange .eq. 'ewald' )  then

    deallocate( km_coul%kptk , km_coul%kpt )
    deallocate( km_coul%strf )

  endif

  return

END SUBROUTINE efg_mesh_dealloc

END MODULE efg
