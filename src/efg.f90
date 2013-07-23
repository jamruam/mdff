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
#include "symbol.h"
!// general debug flag
!#define debug 
!#define debug_main
!#define debug_es
!#define debug_multipole
!#define debug_efg_stat
!#define fix_grid
!#define debug_non_null_trace_ewald
! ======= Hardware =======

! *********************** MODULE efg  **********************************
!
!> \brief
!! two implementation of direct and Ewald sum are present for EFG calculation
!! the old one ( lefg_old ) does only point charge
!! the new one ( .not. lefg_old ) is based on the multipole expansion 
!! ( see multipole_DS and multipole_ES in field.f90
! 
! **********************************************************************
MODULE efg 
 
  USE constants,                ONLY :  dp 
  USE kspace,                   ONLY :  kmesh  
  USE rspace,                   ONLY :  rmesh
  USE config,                   ONLY :  ntypemax
  USE debug,                    ONLY : intdebug

  implicit none

  logical :: lefgprintall               !< print ( or not ) all the efg for each atoms and configs to file EFGALL
  logical :: lefg_restart               !< if EFGALL files are ready
  logical :: lefg_old                   !< use efg_DS and efg_ES ( old routines ) 
  logical :: lefg_stat                  !< compute statitics distribution on EFG's 
  logical :: lvasp_units               
  logical :: lefg_it_contrib            !< only on kind is contributing too efg (default false)
  integer :: ncefg                      !< number of configurations READ  for EFG calc (only when calc = 'efg')
  integer :: ntcor                      !< maximum number of steps for the acf calculation (calc = 'efg+acf')
  integer :: it_efg                     !< type which is contributing to efg_ia
  real(kind=dp) :: dt                   !< timestep for acf calculation should be the same has the one define in md

#ifdef fix_grid
  real(kind=dp)   , dimension(:,:)  , allocatable :: rgrid
#endif

  real(kind=dp)   , dimension(:,:)  , allocatable :: mu   !< electric dipoles

  ! ===================
  !  efg tensors 
  ! ===================
  real(kind=dp)   , dimension(:,:,:), allocatable :: efg_t         !< efg_tensor
  real(kind=dp)   , dimension(:,:,:), allocatable :: efg_ia        !< efg_tensor

  ! ===============
  !  distributions
  ! ===============
  integer         , dimension(:,:)  , allocatable :: dibvzztot !< vzz distrib. dim. = (ntype , PAN)
  integer         , dimension(:,:)  , allocatable :: dibetatot !< eta distrib. dim. = (ntype , PAN) 
  integer         , dimension(:,:,:), allocatable :: dibUtot   !< U_i distrib. dim. = ( 1:6 , ntype , PAN ) 

  real(kind=dp)                                   :: reseta    !< resolution in eta distribution
  real(kind=dp)                                   :: resu      !< resolution in Ui ditribution
  real(kind=dp)                                   :: resvzz    !< resolution in vzz distribution
  real(kind=dp)                                   :: vzzmin    !< minimum value for vzz (distrib. in [vzzmin, -vzzmin]
  real(kind=dp)                                   :: umin      !< minimum value of umin
  integer                                         :: PANeta    !< nb of bins in eta distrib. (related reseta)
  integer                                         :: PANvzz    !< nb of bins in vzz distrib. (related resvzz)
  integer                                         :: PANU      !< nb of bins in Ui  distrib. (related resu)
  real(kind=dp)                                   :: correct
  logical,SAVE                                    :: lmesh_update = .false. !< true if the kmesh need to be updated

CONTAINS

! *********************** SUBROUTINE efg_init **********************************
!
!> \brief
!! initialisation
!
! ******************************************************************************
SUBROUTINE efg_init

  USE io_file,                  ONLY :  ionode , stdin , stdout , stderr 
  USE control,                  ONLY :  calc
 
  implicit none
  
  ! local
  integer            :: ioerr
  character(len=132) :: filename

  namelist /efgtag/  lefgprintall     , & 
                     lefg_restart     , &
                     lefg_old         , &
                     lefg_it_contrib  , &
                     lefg_stat        , &
                     ncefg            , & 
                     ntcor            , &
                     resvzz           , &
                     reseta           , &
                     resu             , &
                     vzzmin           , &
                     correct          , &
                     lvasp_units      , &
                     dt               , & 
                     it_efg           , &
                     umin           

  if ( calc .ne. 'efg' .and. calc .ne. 'efg+acf' ) return 

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
    io_node WRITE ( stderr, '(a)') 'ERROR reading input_file : efgtag section is absent'
    STOP
  elseif ( ioerr .gt. 0 )  then
    io_node WRITE ( stderr, '(a,i8)') 'ERROR reading input_file : efgtag wrong tag'
    STOP
  endif

  CLOSE ( stdin )
  ! ===============
  !  check efgtag
  ! ===============
  CALL efg_check_tag
   if ( calc .eq. 'efg' ) return
  ! =============
  !  print info
  ! =============
  CALL efg_print_info( stdout )

  return

END SUBROUTINE efg_init

! *********************** SUBROUTINE efg_default_tag ***************************
!
! set default values for efgtag
!
! ******************************************************************************
SUBROUTINE efg_default_tag

  implicit none
  
  ! =================
  !  default values
  ! =================
  lefgprintall       = .true. 
  lefg_old           = .false.
  lefg_restart       = .false.
  lefg_stat          = .false.
  lvasp_units        = .false.
  reseta             =   0.1_dp
  resvzz             =   0.1_dp
  resu               =   0.1_dp 
  ncefg              =   0
  umin               =  -4.0_dp
  vzzmin             =  -4.0_dp
  correct            =  0.0_dp 
  ntcor              =  10

  return

END SUBROUTINE efg_default_tag


! *********************** SUBROUTINE efg_check_tag *****************************
!
! check efg parameters ( check longrange and ncell ) 
!
! ******************************************************************************
SUBROUTINE efg_check_tag

  USE control,                  ONLY :  calc , longrange
  USE io_file,                  ONLY :  ionode , stderr
  USE config,                   ONLY :  ntype

  implicit none

  if ( calc .eq. 'efg+acf' .and. ntype .gt. 2 ) then
    io_node WRITE ( stderr , '(a)' ) 'ERROR the subroutine efg_acf is not implemented for ntype > 2'
    STOP 
  endif
  
  if ( calc .eq. 'efg+acf' .and. dt .eq. 0.0d0 ) then
    io_node WRITE ( stderr , '(a)' ) 'ERROR efgtag : dt should be set for efg_acf calculation'
    STOP
  endif

  if ( calc .eq. 'efg+acf' ) return

  ! ==================================
  !  check vzzmin and Umin .ne. 0
  ! ==================================
  if ( vzzmin .eq. 0_dp .or. umin .eq. 0_dp ) then
    io_node WRITE ( stderr ,'(a,2f8.3)') 'ERROR efgtag: vzzmin or umin should be set',vzzmin,umin
    STOP
  endif             
  ! ==========================================
  !  set PAN (nb bins) from resolution values
  ! ==========================================
  PANeta = int ( 1.0_dp/reseta)
  PANvzz = int ((2.0_dp*ABS (vzzmin))/resvzz)
  PANU   = int ((2.0_dp*ABS (umin))/resu) 
 
  if ( ncefg .eq. 0 ) then
    io_node WRITE ( stderr ,'(a,2f8.3)') 'ERROR efgtag: ncefg is zero but calc=efg requested',ncefg
    STOP
  endif

  return

END SUBROUTINE efg_check_tag

! *********************** SUBROUTINE efg_print_info ****************************
!
! print information to standard output about efg calculation 
!
! ******************************************************************************
SUBROUTINE efg_print_info(kunit)

  USE io_file,                  ONLY :  ionode 
  USE control,                  ONLY :  calc , longrange , cutlongrange
  USE config,                   ONLY :  ntype , atypei , natmi , simu_cell , atype
  USE constants,                ONLY :  tpi ,dzero , done
  USE field,                    ONLY :  qch , alphaES , ncelldirect , kES

  implicit none

  ! global 
  integer :: kunit 

  if ( ionode ) then
    if ( calc .eq. 'efg' ) then
      separator(kunit)
      blankline(kunit) 
      WRITE ( kunit ,'(a)')                     'electric field gradient:'
      WRITE ( kunit ,'(a)')                     'point charges calculation'
      if ( calc .eq. 'efg')  then 
        WRITE ( kunit ,'(a)')                   'read config from file            : TRAJFF'
        WRITE ( kunit ,'(a,i10)')               'numbers of config read ncefg     = ',ncefg 
      endif
      WRITE ( kunit ,'(a)')                     'write averages values in file    : EFGFF '          

      WRITE ( kunit ,'(a)')                     'Distributions:'
      WRITE ( kunit ,'(a)')                     'eta distrib                      : DTETAFF'
      WRITE ( kunit ,'(a)')                     'Vzz distrib                      : DTVZZFF'
      WRITE ( kunit ,'(a)')                     'Ui components distrib            : DTIBUFF'           
      if ( lefgprintall )  then
        WRITE ( kunit ,'(a)')                   'EFG for all atoms                : EFGALL '
      endif
      WRITE ( kunit ,'(a)')                     'distributions parameters:'
      WRITE ( kunit ,500)               'eta between ', dzero,' and ',   done,' with res. ',reseta
      WRITE ( kunit ,500) &
                                                'vzz between ',vzzmin,' and ',-vzzmin,' with res. ',resvzz
      WRITE ( kunit ,500) &
                                                'Ui  between ',  umin,' and ',  -umin,' with res. ',resu
    endif
    if ( calc .eq. 'efg+acf' ) then
      separator(kunit) 
      blankline(kunit) 
      WRITE ( kunit ,'(a)')                     'electric field gradient auto-correlation function:'                    
      WRITE ( kunit ,'(a)')                     'data from file EFGALL'
      WRITE ( kunit ,'(a,i10)')                 'numbers of config read ncefg           = ',ncefg 
      WRITE ( kunit ,'(a,i10)')                 'maximum of evaluated correlation step  = ',ntcor
      WRITE ( kunit ,'(a)')                     'output file                            : EFGACFFF'
    endif
  endif

  return

500 FORMAT(a,e10.3,a,e10.3,a,e10.3) 

END SUBROUTINE efg_print_info

! *********************** SUBROUTINE efgcalc ***********************************
!
! this subroutine initialize the calculation of efg when calc = 'efg' 
! from file TRAJFF.  
!
! ******************************************************************************
SUBROUTINE efgcalc 

  USE io_file,                  ONLY :  ionode , stdout , stderr , kunit_EFGALL , kunit_TRAJFF , kunit_EFGFF , &
                                        kunit_NMRFF , kunit_DTETAFF , kunit_DTVZZFF , kunit_DTIBUFF , io_open , io_close
  USE constants,                ONLY :  fpi , e_2 
  USE config,                   ONLY :  system , natm , ntype , atype , rx , ry , rz , itype , & 
                                        atypei , natmi, rho , simu_cell , config_alloc , qia, &
                                        dipia , dipia_ind , dipia_wfc , ipolar , fx , fy , fz , phi_coul_tot , config_print_info 
  USE control,                  ONLY :  longrange , myrank , numprocs
  USE field,                    ONLY :  qch , dip , field_init , lpolar , lwfc , & 
                                        moment_from_pola , moment_from_wfc , rm_coul , &
                                        km_coul , alphaES , multipole_DS , multipole_ES , field_print_info
  USE thermodynamic,            ONLY :  u_coul_tot , vir_coul_tot , u_pol
  USE cell,                     ONLY :  lattice , dirkar

  implicit none

  INCLUDE 'mpif.h'

  ! local
  integer                                               :: i
  integer                                               :: iastart , iaend 
  integer                                               :: ikstart , ikend
  integer                                               :: ia , iconf , it 
  integer                                               :: iiii , ierr
  character(len=60)                                     :: cccc 
  real(kind=dp)                                         :: aaaa , rho_new
  real(kind=dp)                                         :: ttt1 , ttt2
  real(kind=dp) , dimension ( : , : )     , allocatable :: ef_tmp
  real(kind=dp) , dimension ( : , : , : ) , allocatable :: efg_tmp
#ifdef fix_grid
  real(kind=dp) , dimension ( : , : )     , allocatable :: rave !average positions
#endif
  ! =====================================================
  !   type of positions coordinates 
  ! =====================================================
  logical :: allowed
  character(len=60), SAVE :: cpos
  character(len=60), SAVE :: cpos_allowed(4)
  data cpos_allowed / 'Direct' , 'D' , 'Cartesian' , 'C' /


  ! ==================================
  !  if lefg_restart EFGALL are ready
  ! ==================================
  if ( .not. lefg_restart ) then

    OPEN (unit = kunit_EFGALL  ,file = 'EFGALL', STATUS='REPLACE')

    OPEN (UNIT = kunit_TRAJFF , FILE = 'TRAJFF')

    READ ( kunit_TRAJFF , * ) natm 
    READ ( kunit_TRAJFF , * ) system
    READ ( kunit_TRAJFF , * ) simu_cell%A ( 1 , 1 ) , simu_cell%A ( 2 , 1 ) , simu_cell%A ( 3 , 1 )
    READ ( kunit_TRAJFF , * ) simu_cell%A ( 1 , 2 ) , simu_cell%A ( 2 , 2 ) , simu_cell%A ( 3 , 2 )
    READ ( kunit_TRAJFF , * ) simu_cell%A ( 1 , 3 ) , simu_cell%A ( 2 , 3 ) , simu_cell%A ( 3 , 3 )
    READ ( kunit_TRAJFF , * ) ntype
    READ ( kunit_TRAJFF , * ) ( atypei ( it ) , it = 1 , ntype )
    IF ( ionode ) WRITE ( stdout ,'(A,20A3)' ) 'found type information on TRAJFF : ', atypei ( 1:ntype )
    READ ( kunit_TRAJFF , * )   ( natmi ( it ) , it = 1 , ntype )
    READ ( kunit_TRAJFF , * ) cpos
    ! ======
    !  cpos
    ! ======
    do i = 1 , size( cpos_allowed )
      if ( trim(cpos) .eq. cpos_allowed(i))  allowed = .true.
    enddo
    if ( .not. allowed ) then
      if ( ionode )  WRITE ( stdout , '(a)' ) 'ERROR in POSFF at line 9 should be ', cpos_allowed
      STOP
    endif


    allocate ( ef_tmp  ( natm , 3     ) )
    allocate ( efg_tmp ( natm , 3 , 3 ) )

    CALL lattice ( simu_cell ) 
    rho = natm / simu_cell%omega
    ! ===================================
    !  here we know natm, then alloc 
    !  and decomposition can be applied 
    ! ================================== 
    CALL config_alloc 
    ! read charge in fieldtag
    CALL field_init
    CALL do_split ( natm , myrank , numprocs , iastart , iaend , 'atoms' )
    CALL do_split ( km_coul%nkcut , myrank , numprocs , ikstart , ikend ,'k-pts')
    CALL efg_alloc
    CALL efg_mesh_alloc
    CALL typeinfo_init
    ! =============
    !  print info
    ! =============
    CALL efg_print_info(stdout)
  
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

      if ( iconf .ne. 1 ) READ ( kunit_TRAJFF, * )  natm
      if ( iconf .ne. 1 ) READ ( kunit_TRAJFF, * )  system 
      if ( iconf .ne. 1 ) READ ( kunit_TRAJFF, * )  simu_cell%A ( 1 , 1 ) , simu_cell%A ( 2 , 1 ) , simu_cell%A ( 3 , 1 )
      if ( iconf .ne. 1 ) READ ( kunit_TRAJFF, * )  simu_cell%A ( 1 , 2 ) , simu_cell%A ( 2 , 2 ) , simu_cell%A ( 3 , 2 ) 
      if ( iconf .ne. 1 ) READ ( kunit_TRAJFF, * )  simu_cell%A ( 1 , 3 ) , simu_cell%A ( 2 , 3 ) , simu_cell%A ( 3 , 3 )
      if ( iconf .ne. 1 ) READ ( kunit_TRAJFF, * )  iiii 
      if ( iconf .ne. 1 ) READ ( kunit_TRAJFF , * ) ( cccc , it = 1 , ntype )
      if ( iconf .ne. 1 ) READ ( kunit_TRAJFF , * ) ( iiii , it = 1 , ntype )
      if ( iconf .ne. 1 ) READ ( kunit_TRAJFF , * ) cpos
      ! ======
      !  cpos
      ! ======
      do i = 1 , size( cpos_allowed )
        if ( trim(cpos) .eq. cpos_allowed(i))  allowed = .true.
      enddo
      if ( .not. allowed ) then
        if ( ionode )  WRITE ( stdout , '(a)' ) 'ERROR in POSFF at line 9 should be ', cpos_allowed
        STOP
      endif

      do ia = 1 , natm
        READ ( kunit_TRAJFF , * ) atype ( ia ) , rx ( ia ) , ry ( ia ) , rz ( ia ) ,aaaa,aaaa,aaaa,aaaa,aaaa,aaaa
      enddo
      if ( cpos .eq. 'Direct' ) then
        ! ======================================
        !         direct to cartesian
        ! ======================================
        CALL dirkar ( natm , rx , ry , rz , simu_cell%A )
        if ( ionode .and. iconf .eq. 1 ) WRITE ( stdout      ,'(A,20A3)' ) 'atomic positions in direct coordinates in POSFF'
      else if ( cpos .eq. 'Cartesian' ) then
        if ( ionode .and. iconf .eq. 1 ) WRITE ( stdout      ,'(A,20A3)' ) 'atomic positions in cartesian coordinates in POSFF'
      endif

      CALL lattice ( simu_cell )
      rho_new = natm / simu_cell%omega
      if ( rho_new .ne. rho ) then
        rho = rho_new
        lmesh_update = .true.
      endif
      if ( iconf .ne. 1 .and. lmesh_update ) CALL field_init 
      if ( iconf .ne. 1 .and. lmesh_update ) CALL efg_mesh_alloc

#ifdef debug
      ! ============================
      !  print minimum distance 
      !  and distance distribution
      ! ============================
      call distance_tab 
#endif

      ! =======================
      !  total tensor (efg_t)
      ! =======================
      efg_t    = 0.0_dp

      ! ======================================
      !     induced moment from polarisation 
      ! ======================================
      CALL moment_from_pola ( iastart , iaend , ikstart , ikend , dipia_ind )

      ! ======================================
      !     induced moment from wannier center  
      ! ======================================
      CALL moment_from_WFc ( dipia_wfc ) 

      ! ======================================
      !  total dipole mu :
      !   static +  induced  + "wannier"
      ! ======================================
      mu = dipia + dipia_ind + dipia_wfc

#ifdef debug_multipole
      if ( longrange .eq. 'ewald' )  CALL multipole_ES ( iastart, iaend , ikstart , ikend , ef_tmp , efg_tmp , mu , u_coul_tot , vir_coul_tot , phi_coul_tot ) 
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

      efg_t    = efg_ia 

      if ( lvasp_units ) then
        efg_t    = - efg_t    * e_2 
      endif

#ifdef debug_main
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

    if ( lmesh_update ) CALL efg_mesh_dealloc
    ttt2 = MPI_WTIME(ierr)
    io_node WRITE ( stdout , 110 ) 'config : ',iconf,' EFG  ', ttt2 - ttt1  

    enddo ! iconf loop

    CLOSE(kunit_EFGALL)

  endif !lefg_restart

  CALL MPI_BARRIER( MPI_COMM_WORLD , ierr )
  
  io_node blankline ( stdout ) 
  io_node WRITE ( stdout , '(a)' ) 'calculate statistical properties from EFGALL'

  ! ========================================================
  !     calculate statistical properties from EFGALL
  ! ========================================================
  CALL io_open  ( kunit_EFGALL , 'EFGALL' , 'unknown' )
  CALL io_open  ( kunit_EFGFF  , 'EFGFF'  , 'unknown' )
  CALL io_open  ( kunit_NMRFF  , 'NMRFF'  , 'unknown' )
  CALL efg_stat ( kunit_EFGALL , kunit_EFGFF , kunit_NMRFF )
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

  CLOSE(kunit_TRAJFF)

  CALL efg_dealloc

  return

110   FORMAT(2X,A8,I6,A20,' :  cpu time',F9.2)

END SUBROUTINE efgcalc

! *********************** SUBROUTINE efg_DS ************************************
!
! (only point charges )
! Direct summation to calculate electric-field-gradient.
! parallelisation atom decomposition
! 
! ******************************************************************************
SUBROUTINE efg_DS ( iastart , iaend , rm )

  USE control,                  ONLY :  myrank , numprocs , calc , cutlongrange
  USE config,                   ONLY :  system , natm , natmi , atype , atypei , itype , &
                                        rx , ry , rz , ntype , qia , simu_cell
  USE cell,                     ONLY :  kardir , dirkar
  USE time,                     ONLY :  efgtimetot1 , efgtimetot3
  USE io_file,                  ONLY :  stdout

  implicit none

  INCLUDE 'mpif.h'

  ! global
  integer, intent ( in )         :: iastart , iaend 
  TYPE ( rmesh ) , intent ( in ) :: rm

  ! local
  integer :: ia, ja, ierr 
  integer :: ncell
  real(kind=dp) :: d , d2 , d5 , dm5
  real(kind=dp) :: rxi , ryi , rzi , rxj , ryj , rzj , rxij , ryij , rzij , sxij , syij , szij 
  real(kind=dp) :: cutefgsq
  real(kind=dp) :: ttt1 , ttt2 , ttt3
  logical :: lcharge

  ! ===============================================
  !  check if there is any charge otherwise return
  ! ===============================================
  lcharge= .false.
  do ia = 1 , natm
    if ( qia ( ia ) .ne. 0.0_dp ) lcharge = .true.
  enddo
  if ( .not. lcharge ) return


#ifdef debug
  CALL print_config_sample(0,0)
#endif

  ttt1 = MPI_WTIME(ierr)

  cutefgsq = cutlongrange * cutlongrange 
  efg_ia(:,:,:) = 0.0_dp

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
  CALL kardir ( natm , rx , ry , rz , simu_cell%B )

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
              dm5 = 1.0_dp/d5
              dm5 = dm5 * qia(ja)

              efg_ia( ia , 1 , 1 )  = efg_ia( ia , 1 , 1 ) - (3.0_dp * rxij * rxij - d2 ) * dm5
              efg_ia( ia , 2 , 2 )  = efg_ia( ia , 2 , 2 ) - (3.0_dp * ryij * ryij - d2 ) * dm5
              efg_ia( ia , 3 , 3 )  = efg_ia( ia , 3 , 3 ) - (3.0_dp * rzij * rzij - d2 ) * dm5
              efg_ia( ia , 1 , 2 )  = efg_ia( ia , 1 , 2 ) -  3.0_dp * rxij * ryij * dm5
              efg_ia( ia , 1 , 3 )  = efg_ia( ia , 1 , 3 ) -  3.0_dp * rxij * rzij * dm5
              efg_ia( ia , 2 , 3 )  = efg_ia( ia , 2 , 3 ) -  3.0_dp * ryij * rzij * dm5

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
                dm5 = 1.0_dp/d5
                dm5 = dm5 * qia(ja)

                efg_ia( ia , 1 , 1 )  = efg_ia( ia , 1 , 1 ) - (3.0_dp * rxij * rxij - d2 ) * dm5
                efg_ia( ia , 2 , 2 )  = efg_ia( ia , 2 , 2 ) - (3.0_dp * ryij * ryij - d2 ) * dm5
                efg_ia( ia , 3 , 3 )  = efg_ia( ia , 3 , 3 ) - (3.0_dp * rzij * rzij - d2 ) * dm5
                efg_ia( ia , 1 , 2 )  = efg_ia( ia , 1 , 2 ) -  3.0_dp * rxij * ryij * dm5
                efg_ia( ia , 1 , 3 )  = efg_ia( ia , 1 , 3 ) -  3.0_dp * rxij * rzij * dm5
                efg_ia( ia , 2 , 3 )  = efg_ia( ia , 2 , 3 ) -  3.0_dp * ryij * rzij * dm5

              endif ! d2.lt.cutefgsq

            endif ! ia.ne.ja
          enddo ! ja
        endif 

     enddo ! ncell

    efg_ia( ia , 2 , 1 ) = efg_ia( ia , 1 , 2 )
    efg_ia( ia , 3 , 1 ) = efg_ia( ia , 1 , 3 )
    efg_ia( ia , 3 , 2 ) = efg_ia( ia , 2 , 3 )

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

#ifdef debug
  CALL print_tensor( efg_ia( 1    , : , : ) , 'EFG_1A  ' )
  CALL print_tensor( efg_ia( natm , : , : ) , 'EFG_NA  ' )
#endif

  ttt3 = MPI_WTIME(ierr)
  efgtimetot3 = efgtimetot3 + (ttt3-ttt2)

  ! ======================================
  !         direct to cartesian
  ! ======================================
  CALL dirkar ( natm , rx , ry , rz , simu_cell%A )


  return

END SUBROUTINE efg_DS

! *********************** SUBROUTINE efg_ES ************************************
!
! Ewald sum (only point charge)
! the implementation follows J. Chem. Phys, 112, 14, p 6512 (2000)
! Correction (16-05-11) following the J. Chem. Phys. 129, 074102 (2008)
! and quantum-espresso
!
! the reciproc part could be ( should be ) parallized
!
! ******************************************************************************
SUBROUTINE efg_ES ( iastart , iaend , km , alphaES )

  USE control,                  ONLY :  myrank , numprocs , calc
  USE config,                   ONLY :  system , natm , natmi , atype , atypei , &
                                        simu_cell , itype , rx , ry , rz , ntype , qia 
  USE constants,                ONLY :  pi , fpi , piroot , imag
  USE field,                    ONLY :  qch
  USE kspace,                   ONLY :  struc_fact
  USE cell,                     ONLY :  kardir , dirkar 
  USE time,                     ONLY :  strftimetot , efgtimetot1 , efgtimetot2 , efgtimetot3
  USE io_file,                  ONLY :  stdout

  implicit none

  INCLUDE 'mpif.h'

  ! global
  integer, intent(in) :: iastart , iaend 
  TYPE ( kmesh ) :: km
  real(kind=dp) :: alphaES

  ! local
  !integer            :: it 
  integer            :: ia , ja ,  ierr 
  real(kind=dp)      :: d , d2 , d4 , d3 , d5 , expon 
  real(kind=dp)      :: alpha2 , alpha3
  real(kind=dp)      :: allrealpart
  real(kind=dp)      :: T0 , T1 , T2  ! real part 
  real(kind=dp)      :: rxi , ryi , rzi , rxij , ryij , rzij , sxij , syij , szij
  real(kind=dp)      :: ak, kx , ky , kz , kk
  integer            :: ik
  real(kind=dp)      :: kri 
  complex(kind=dp)   :: rhon , carg , recarg , recarg_dgg 
  real(kind=dp), external :: errfc 
  real(kind=dp)      :: ttt1 , ttt2 , ttt3 , ttt4 , ttt5 
  real(kind=dp)      :: efg_ia_real ( natm , 3 , 3 )
  real(kind=dp)      :: efg_ia_dual_real ( natm , 3 , 3 )
  complex(kind=dp)   :: efg_ia_dual ( natm , 3 , 3 )
  real(kind=dp)      :: efg_ia_real_it ( natm , ntype , 3 , 3 )
  real(kind=dp)      :: efg_ia_dual_real_it ( natm , ntype , 3 , 3 )
  complex(kind=dp)   :: efg_ia_dual_it ( natm , ntype , 3 , 3 )
  logical            :: lcharge 

  ! ===============================================
  !  check if there is any charge otherwise return
  ! ===============================================
  lcharge= .false.
  do ia = 1 , natm
    if ( qia ( ia ) .ne. 0.0_dp ) lcharge = .true.
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
  efg_ia_real         = 0.0_dp
  efg_ia_dual         = (0.0_dp,0.0_dp)
  efg_ia_dual_real    = 0.0_dp
  efg_ia              = 0.0_dp
  efg_ia_real_it      = 0.0_dp
  efg_ia_dual_it      = (0.0_dp,0.0_dp)
  efg_ia_dual_real_it = 0.0_dp

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
  ttt2 = MPI_WTIME(ierr)
  strftimetot = strftimetot + (ttt2 - ttt1)

  ! ======================================
  !         cartesian to direct 
  ! ======================================
  CALL kardir ( natm , rx , ry , rz , simu_cell%B )

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
         T1 = ( 2.0_dp * alphaES ) / d4 
         T2 = ( 4.0_dp * alpha3  ) / d2 / 3.0_dp
         allrealpart = qia(ja) * ( T0 + ( T1 + T2 ) *expon )

         efg_ia_real( ia , 1 , 1 ) = efg_ia_real( ia , 1 , 1 ) - ( 3.0_dp * rxij * rxij - d2 ) * allrealpart
         efg_ia_real( ia , 2 , 2 ) = efg_ia_real( ia , 2 , 2 ) - ( 3.0_dp * ryij * ryij - d2 ) * allrealpart
         efg_ia_real( ia , 3 , 3 ) = efg_ia_real( ia , 3 , 3 ) - ( 3.0_dp * rzij * rzij - d2 ) * allrealpart
         efg_ia_real( ia , 1 , 2 ) = efg_ia_real( ia , 1 , 2 ) -   3.0_dp * rxij * ryij        * allrealpart
         efg_ia_real( ia , 1 , 3 ) = efg_ia_real( ia , 1 , 3 ) -   3.0_dp * rxij * rzij        * allrealpart
         efg_ia_real( ia , 2 , 3 ) = efg_ia_real( ia , 2 , 3 ) -   3.0_dp * ryij * rzij        * allrealpart

       endif

     enddo

  enddo atom1

  ttt3 = MPI_WTIME(ierr)
  efgtimetot1 = efgtimetot1 + (ttt3-ttt2)

  ! ======================================
  !         cartesian to direct 
  ! ======================================
  CALL dirkar ( natm , rx , ry , rz , simu_cell%A )

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
    Ak   = EXP ( - kk * 0.25_dp / alpha2 ) 

    if (km%kptk(ik) .eq. 0 ) then
      WRITE ( stdout , * ) 'the sum should be done on k! =  0',ik
      STOP
    endif
    ! ===============================
    !                              ---
    !  charge density in k-space ( \   q * facteur de structure  )
    !                              /__
    ! ===============================
    rhon   = (0.0_dp, 0.0_dp)
    do ja = 1, natm
      rhon = rhon +  qia(ja) * CONJG( km%strf ( ik , ja ) )
    enddo

    do ia = 1 , natm
#ifdef fix_grid
     rxi = rgrid ( 1 , ia )
     ryi = rgrid ( 2 , ia )
     rzi = rgrid ( 3 , ia ) 
#endif
      rxi = rx(ia)
      ryi = ry(ia)
      rzi = rz(ia)
      kri = ( kx * rxi + ky * ryi + kz * rzi )
      carg = EXP ( imag * kri )
      recarg_dgg =  rhon * carg * Ak / kk
      recarg = recarg_dgg * kk
      ! =========================
      ! electric field gradient
      ! =========================
      efg_ia_dual ( ia , 1 , 1 ) = efg_ia_dual ( ia , 1 , 1 ) +  3.0_dp * kx * kx * recarg_dgg - recarg
      efg_ia_dual ( ia , 2 , 2 ) = efg_ia_dual ( ia , 2 , 2 ) +  3.0_dp * ky * ky * recarg_dgg - recarg
      efg_ia_dual ( ia , 3 , 3 ) = efg_ia_dual ( ia , 3 , 3 ) +  3.0_dp * kz * kz * recarg_dgg - recarg
      efg_ia_dual ( ia , 1 , 2 ) = efg_ia_dual ( ia , 1 , 2 ) +  3.0_dp * kx * ky * recarg_dgg
      efg_ia_dual ( ia , 1 , 3 ) = efg_ia_dual ( ia , 1 , 3 ) +  3.0_dp * kx * kz * recarg_dgg
      efg_ia_dual ( ia , 2 , 3 ) = efg_ia_dual ( ia , 2 , 3 ) +  3.0_dp * ky * kz * recarg_dgg

    enddo

  enddo kpoint

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

  ! =======
  ! 4pi/3V
  ! =======
  efg_ia_dual_real( : , : , :) = REAL ( efg_ia_dual( : , : , :) , kind = dp ) 
  efg_ia_dual_real =  efg_ia_dual_real * fpi / simu_cell%omega / 3.0_dp

  ! ==============
  !  total tensor
  ! ==============
  efg_ia = efg_ia_dual_real + efg_ia_real 

!=========================================================
!      END OF EFG TENSOR CALCULATION
!=========================================================

  ttt5 = MPI_WTIME(ierr)
  efgtimetot3 = efgtimetot3 + (ttt5-ttt4)

#ifdef debug 
  CALL print_tensor( efg_ia_real  ( 1 , : , : ) , 'EFG_DIRO' )
  CALL print_tensor( efg_ia_dual_real  ( 1 , : , : ) , 'EFG_RECO' )
  CALL print_tensor( efg_ia   ( 1 , : , : ) , 'EFG_TOTO' )
#endif


  return

END SUBROUTINE efg_ES

! *********************** SUBROUTINE multipole_efg_DS **************************
!
! New version of the Direct summation to calculate electric-field-gradient.
! from point charges and dipoles
! parallelisation atom decomposition
! Based on the multipole_DS subroutine in field.f90
! 
! ******************************************************************************
SUBROUTINE multipole_efg_DS ( iastart, iaend , rm , mu )

  USE config,                   ONLY :  natm , atype , natmi , ntype , qia , &
                                        rx , ry , rz , fx , fy , fz , tau_coul , simu_cell ,itype 
  USE control,                  ONLY :  cutlongrange , myrank
  USE io_file,                  ONLY :  stdout , ionode 
  USE field,                    ONLY :  qch
  USE cell,                     ONLY :  kardir , dirkar
  USE time,                     ONLY :  efgtimetot1 , efgtimetot3

  implicit none

  INCLUDE 'mpif.h'

  ! global 
  integer, intent(in) :: iastart , iaend
  TYPE ( rmesh ) :: rm
  real(kind=dp) :: mu    ( natm , 3 )

  ! local 
  integer :: ia, ja , ierr , ncell 
  real(kind=dp) :: cutsq
  real(kind=dp) :: rxi , ryi , rzi 
  real(kind=dp) :: rxj , ryj , rzj
  real(kind=dp) :: rxij , ryij , rzij
  real(kind=dp) :: sxij , syij , szij
  real(kind=dp) :: qj 
  real(kind=dp) :: mujx , mujy , mujz
  real(kind=dp) :: Txx , Tyy , Tzz , Txy , Txz , Tyz
  real(kind=dp) :: Txxx,  Tyyy,  Tzzz, Txxy, Txxz, Tyyx, Tyyz, Tzzx, Tzzy, Txyz
  real(kind=dp) :: d , d2 
  real(kind=dp) :: dm5 , dm7 

  real(kind=dp) :: ttt1 , ttt2 , ttt3
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

  efg_ia   = 0.0_dp

  ! ======================================
  !         cartesian to direct 
  ! ======================================
  CALL kardir ( natm , rx , ry , rz , simu_cell%B )
  
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

        if ( lefg_it_contrib .and. itype(ja).ne.it_efg ) cycle  
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

            if ( d2 .lt. cutsq .and. d2 .ne. 0.0_dp ) then 

              qj    = qia ( ja )
              mujx  = mu ( ja , 1 )  
              mujy  = mu ( ja , 2 )  
              mujz  = mu ( ja , 3 )  
              d     = SQRT ( d2 )
              dm5   = 1.0_dp / ( d2 * d2 * d )
              dm7   = dm5 / d2 * 3.0_dp
 
              ! multipole interaction tensor rank = 2
              Txx = ( 3.0_dp * rxij * rxij - d2 ) * dm5
              Tyy = ( 3.0_dp * ryij * ryij - d2 ) * dm5
              Tzz = ( 3.0_dp * rzij * rzij - d2 ) * dm5
              Txy = ( 3.0_dp * rxij * ryij      ) * dm5
              Txz = ( 3.0_dp * rxij * rzij      ) * dm5
              Tyz = ( 3.0_dp * ryij * rzij      ) * dm5

              ! multipole interaction tensor rank = 3  
              Txxx = ( 5.0_dp * rxij * rxij * rxij -  3.0_dp * d2 * ( rxij ) ) * dm7 
              Tyyy = ( 5.0_dp * ryij * ryij * ryij -  3.0_dp * d2 * ( ryij ) ) * dm7 
              Tzzz = ( 5.0_dp * rzij * rzij * rzij -  3.0_dp * d2 * ( rzij ) ) * dm7 
              Txxy = ( 5.0_dp * rxij * rxij * ryij -           d2 * ( ryij ) ) * dm7 
              Txxz = ( 5.0_dp * rxij * rxij * rzij -           d2 * ( rzij ) ) * dm7 
              Tyyx = ( 5.0_dp * ryij * ryij * rxij -           d2 * ( rxij ) ) * dm7 
              Tyyz = ( 5.0_dp * ryij * ryij * rzij -           d2 * ( rzij ) ) * dm7 
              Tzzx = ( 5.0_dp * rzij * rzij * rxij -           d2 * ( rxij ) ) * dm7 
              Tzzy = ( 5.0_dp * rzij * rzij * ryij -           d2 * ( ryij ) ) * dm7 
              Txyz = ( 5.0_dp * rxij * ryij * rzij                           ) * dm7 

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
  CALL dirkar ( natm , rx , ry , rz , simu_cell%A )

  return

END SUBROUTINE multipole_efg_DS

! *********************** SUBROUTINE multipole_efg_ES **************************
!
! New version of the Ewald summation to calculate electric-field-gradient.
! from point charges and dipoles
! parallelisation atom decomposition
! Based on the multipole_ES subroutine in field.f90
! 
! ******************************************************************************
SUBROUTINE multipole_efg_ES ( iastart , iaend , km , alphaES , mu )

  USE control,                  ONLY :  lsurf , cutshortrange 
  USE config,                   ONLY :  natm , ntype , natmi , atype , &
                                        rx , ry , rz , fx , fy , fz ,  &
                                        qia , simu_cell , itype
  USE constants,                ONLY :  imag , pi , piroot , tpi , fpi
  USE io_file,                  ONLY :  ionode , stdout 
  USE field,                    ONLY :  qch
  USE kspace,                   ONLY :  struc_fact        
  USE cell,                     ONLY :  kardir , dirkar
  USE time,                     ONLY :  strftimetot , efgtimetot1 , efgtimetot2 , efgtimetot3           

  implicit none

  INCLUDE 'mpif.h'

  ! global 
  integer, intent(in) :: iastart , iaend
  TYPE ( kmesh ) :: km
  real(kind=dp) :: alphaES
  real(kind=dp) :: mu    ( natm , 3 )

  ! local 
  integer             :: ia , ja , ik , it , ip , ierr
  real(kind=dp), dimension(:,:,:), allocatable :: efg_dir , efg_rec , efg_self

  real(kind=dp) :: mip    ( ntype , 3 )

  real(kind=dp) :: rxi  , ryi  , rzi
  real(kind=dp) :: rxj  , ryj  , rzj
  real(kind=dp) :: kx   , ky   , kz 
  real(kind=dp) :: rxij , ryij , rzij
  real(kind=dp) :: sxij , syij , szij
  real(kind=dp) :: kk   , kri  , Ak 
  real(kind=dp) :: qj  
  real(kind=dp) :: mujx , mujy , mujz
  complex(kind=dp)   :: rhon , carg 
  real(kind=dp) :: recarg 
  real(kind=dp) :: expon , F1 , F2 , F3 
  real(kind=dp) :: k_dot_mu
  real(kind=dp) :: Txx , Tyy , Tzz , Txy , Txz , Tyz
  real(kind=dp) :: Txxx,  Tyyy,  Tzzz, Txxy, Txxz, Tyyx, Tyyz, Tzzx, Tzzy, Txyz
  real(kind=dp) :: d , d2 , d3  , d5 
  real(kind=dp) :: dm1 , dm5 , dm7 
  real(kind=dp) :: alpha2 , alpha3 , alpha5 
  real(kind=dp) :: selfa 
  real(kind=dp), external :: errfc
  real(kind=dp) :: fpi_V 
  real(kind=dp) :: ttt1 , ttt2  , ttt3 , ttt4 
  real(kind=dp) :: ttt1p , ttt2p 


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

  efg_ia   = 0.0_dp
  efg_dir  = 0.0_dp
  efg_rec  = 0.0_dp
  efg_self = 0.0_dp

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
  selfa  = - 4.0_dp * alpha3 / 3.0_dp / piroot

  ttt1 = MPI_WTIME(ierr)

  ! ======================================
  !         cartesian to direct 
  ! ======================================
  CALL kardir ( natm , rx , ry , rz , simu_cell%B )

  ! ==============================================
  !        direct space part
  ! ==============================================

  do ia = iastart , iaend
    rxi = rx(ia)
    ryi = ry(ia)
    rzi = rz(ia)

    do ja = 1, natm

      if ( lefg_it_contrib .and. itype(ja).ne.it_efg ) cycle  
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
        dm1 = 1.0_dp / d
        !dm3 = dm1 / d2
        dm5 = 1.0_dp / d2 / d3 
        dm7 = dm5 / d2 * 3.0_dp

        expon = EXP ( - alpha2 * d2 )    / piroot
        F1    = errfc( alphaES * d ) + 2.0_dp * alphaES * d  * expon
        F2    = F1 + 4.0_dp * alpha3  * d3 * expon / 3.0_dp
        F3    = F2 + 8.0_dp * alpha5  * d5 * expon / 15.0_dp

        ! multipole interaction tensor rank = 2
        Txx = ( 3.0_dp * rxij * rxij * F2 - d2 * F1 ) * dm5  
        Tyy = ( 3.0_dp * ryij * ryij * F2 - d2 * F1 ) * dm5 
        Tzz = ( 3.0_dp * rzij * rzij * F2 - d2 * F1 ) * dm5
        Txy = ( 3.0_dp * rxij * ryij * F2           ) * dm5
        Txz = ( 3.0_dp * rxij * rzij * F2           ) * dm5 
        Tyz = ( 3.0_dp * ryij * rzij * F2           ) * dm5 
 
        ! multipole interaction tensor rank = 3  
        Txxx = ( -5.0_dp * rxij * rxij * rxij * F3 +  3.0_dp * d2 * ( rxij ) * F2 ) * dm7 
        Tyyy = ( -5.0_dp * ryij * ryij * ryij * F3 +  3.0_dp * d2 * ( ryij ) * F2 ) * dm7 
        Tzzz = ( -5.0_dp * rzij * rzij * rzij * F3 +  3.0_dp * d2 * ( rzij ) * F2 ) * dm7
        Txxy = ( -5.0_dp * rxij * rxij * ryij * F3 +           d2 * ( ryij ) * F2 ) * dm7 
        Txxz = ( -5.0_dp * rxij * rxij * rzij * F3 +           d2 * ( rzij ) * F2 ) * dm7 
        Tyyx = ( -5.0_dp * ryij * ryij * rxij * F3 +           d2 * ( rxij ) * F2 ) * dm7 
        Tyyz = ( -5.0_dp * ryij * ryij * rzij * F3 +           d2 * ( rzij ) * F2 ) * dm7 
        Tzzx = ( -5.0_dp * rzij * rzij * rxij * F3 +           d2 * ( rxij ) * F2 ) * dm7 
        Tzzy = ( -5.0_dp * rzij * rzij * ryij * F3 +           d2 * ( ryij ) * F2 ) * dm7 
        Txyz = ( -5.0_dp * rxij * ryij * rzij * F3                                ) * dm7 
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
  CALL dirkar ( natm , rx , ry , rz , simu_cell%A )

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
    Ak   = EXP ( - kk * 0.25_dp / alpha2 ) / kk

    if (km%kptk(ik) .eq. 0 ) then
      WRITE ( stdout , * ) 'the sum should be done on k! =  0',ik
      STOP
    endif
    ! ===============================
    !                              ---
    !  charge density in k-space ( \   q * facteur de structure  )
    !                              /__
    ! ===============================
    rhon   = (0.0_dp, 0.0_dp)
    do ja = 1, natm
      if ( lefg_it_contrib .and. (itype(ja) .ne. it_efg) ) cycle
      k_dot_mu = ( mu ( ja , 1 ) * kx + mu ( ja , 2 ) * ky + mu ( ja , 3 ) * kz  ) 
      rhon = rhon + ( qia(ja) + imag * k_dot_mu ) * km%strf ( ik , ja ) 
    enddo

    do ia = 1 , natm

      rxi = rx(ia)
      ryi = ry(ia)
      rzi = rz(ia)
      kri = ( kx * rxi + ky * ryi + kz * rzi )
      carg   = EXP  ( imag * kri )
      recarg = REAL ( CONJG ( rhon ) * carg  * Ak , kind = dp )

      ! electric field gradient
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
   if ( lefg_it_contrib .and. (itype(ia) .ne. it_efg) ) cycle
   efg_self ( ia , 1 , 1 ) = selfa * qia ( ia ) 
   efg_self ( ia , 2 , 2 ) = selfa * qia ( ia ) 
   efg_self ( ia , 3 , 3 ) = selfa * qia ( ia ) 
  enddo

  ! =====================================================
  !                     TOTAL
  ! =====================================================

  efg_ia = ( efg_dir + efg_rec  + efg_self )
  ! =====================
  ! only for consistency 
  ! =====================
  efg_ia ( : , 2 , 1 ) = efg_ia ( : , 1 , 2 )
  efg_ia ( : , 3 , 1 ) = efg_ia ( : , 1 , 3 )
  efg_ia ( : , 3 , 2 ) = efg_ia ( : , 2 , 3 )

  efg_ia ( : , 1 , 1 ) = efg_ia ( : , 1 , 1 ) - correct
  efg_ia ( : , 2 , 2 ) = efg_ia ( : , 2 , 2 ) - correct
  efg_ia ( : , 3 , 3 ) = efg_ia ( : , 3 , 3 ) - correct
  

#ifdef debug_es
  CALL print_tensor( efg_dir  ( 1 , : , : ) , 'EFG_DIRN' )
  CALL print_tensor( efg_rec  ( 1 , : , : ) , 'EFG_RECN' )
  CALL print_tensor( efg_self ( 1 , : , : ) , 'EFG_SELN' )
  CALL print_tensor( efg_ia   ( 1 , : , : ) , 'EFG_TOTN' )
#endif

  deallocate( efg_dir , efg_rec , efg_self ) 

  ttt4 = MPI_WTIME (ierr)
  efgtimetot3 = efgtimetot3 + ( ttt4 - ttt3 )

  return

END SUBROUTINE multipole_efg_ES

! *********************** SUBROUTINE efg_write_output **************************
!
! write distributions of EFG parameters and components to files
!
! ******************************************************************************
SUBROUTINE efg_write_output ( kunit_eta , kunit_vzz , kunit_u )

  USE io_file,                  ONLY :  ionode 
  USE config,                   ONLY :  natm, natmi, ntype 
  USE constants,                ONLY :  dzero
  USE field,                    ONLY :  lwfc

  implicit none
 
  INCLUDE 'mpif.h'

  ! global
  integer, intent(in)  :: kunit_eta , kunit_vzz , kunit_u 

  ! local
  integer :: i , it , saveit0 , totions
  ! make some quantities real to print them 
  ! it also to avoid REAL( ..., kind = dp ) everywhere
  
  real(kind=dp), dimension ( : , :  )     , allocatable :: r_dibetatot 
  real(kind=dp), dimension ( : , :  )     , allocatable :: r_dibvzztot 
  real(kind=dp), dimension ( : , : , :  ) , allocatable :: r_dibUtot 
  real(kind=dp), dimension ( :  )         , allocatable :: r_natmi
  real(kind=dp) :: r_ncefg    

  allocate( r_dibUtot(6,0:ntype,0:PANU) )
  allocate( r_dibvzztot(0:ntype,0:PANvzz) )
  allocate( r_dibetatot(0:ntype,0:PANeta) )
  allocate( r_natmi(0:ntype) )
  r_dibUtot   = 0.0_dp
  r_dibvzztot = 0.0_dp
  r_dibetatot = 0.0_dp
  r_natmi     = 0.0_dp



  r_dibetatot = REAL ( dibetatot , kind = dp ) 
  r_dibvzztot = REAL ( dibvzztot , kind = dp ) 
  r_dibUtot   = REAL ( dibUtot   , kind = dp ) 
  r_natmi     = REAL ( natmi     , kind = dp )
  r_ncefg     = REAL ( ncefg     , kind = dp )

  if ( ionode ) then
    ! ======================== 
    !  write eta distributions
    ! ======================== 
    saveit0=natmi(0)
    totions = 0
    do it = 1 , ntype
      if ( lwfc ( it ) .eq. -1 ) cycle  
      totions = totions + natmi (it)   
    enddo    
    natmi(0)=totions
    WRITE (kunit_eta,'(a,<ntype+2>f15.8)') '#', reseta , ( r_natmi ( it ), it = 0 , ntype )   
    WRITE (kunit_eta,'(<ntype+2>f15.8)')  dzero, ( dzero , it = 0 , ntype )
    do i = 0 , PANeta-1
      WRITE (kunit_eta,'(<ntype+2>f15.8)') &
       REAL ( ( (i+1) - 0.5_dp )  * reseta , kind = dp ), (  r_dibetatot( it , i )  / ( reseta * r_natmi(it) * r_ncefg ) , it = 0 , ntype )   
    enddo
    ! ========================
    !  write Vzz distribution
    ! ========================
    do i = 0 , PANvzz 
      WRITE (kunit_vzz ,'(<ntype+2>f15.8)') &
      vzzmin  + REAL ( i * resvzz , kind = dp ), ( r_dibvzztot(it,i) / ( resvzz * r_natmi(it) * r_ncefg ) , it = 0 , ntype )
    enddo
 
    ! =================================================
    ! write U1 and average Uk (with k>1) distribution 
    ! =================================================
    do i = 0 , PANU 
      WRITE (kunit_u ,'(<2*ntype+3>f15.8)') umin  + REAL ( i * resu , kind = dp ) , &
                                            ( r_dibUtot(1,it,i) / ( resu * r_natmi(it) * r_ncefg ) , &
                                              r_dibUtot(6,it,i) / ( resu * 4.0_dp * r_natmi(it) * r_ncefg) , it = 0 , ntype )
    enddo
    blankline(kunit_eta) 
    blankline(kunit_eta) 
    blankline(kunit_vzz) 
    blankline(kunit_vzz) 
    blankline(kunit_u) 
    blankline(kunit_u) 
 
  natmi(0) = saveit0
  endif

  deallocate( r_dibUtot )
  deallocate( r_dibvzztot )
  deallocate( r_dibetatot )
  deallocate( r_natmi )

  return

END SUBROUTINE efg_write_output

! *********************** SUBROUTINE efg_acf ***********************************
!
! calculate the auto-correlation function of the efg principal component 
! based on Allen-Tieldsley
!
! ******************************************************************************
SUBROUTINE efg_acf

  USE config,                   ONLY :  system , natm , ntype , itype , atype , atypei, natmi , simu_cell , rho , config_alloc
  USE io_file,                  ONLY :  ionode , stdout , stderr , kunit_EFGALL , kunit_EFGACFFF , kunit_NMRACFFF
  USE cell,                     ONLY :  lattice

  implicit none

  integer, parameter                             :: lwork = 6
  integer                                        :: ifail
  integer                                        :: ia , it , t , tt0 , t0 , t0max
!  integer                                        :: ui commented in the code
  integer,       dimension (:,:)   , allocatable :: norm
  real(kind=dp), dimension (:,:)   , allocatable :: acfxx , acfyy , acfzz          ! autocorelation function ( EFG tensor )
  real(kind=dp), dimension (:,:)   , allocatable :: acfxy , acfxz , acfyz          ! autocorelation function ( EFG tensor )
  real(kind=dp), dimension (:,:)   , allocatable :: acf11 , acf22 , acf33 , acfeta ! autocorelation function ( principal components )
!  real(kind=dp), dimension (:,:,:) , allocatable :: acfU                           ! autocorelation function ( principal components )

  real(kind=dp), dimension (:)     , allocatable :: vxx0 , vyy0 , vzz0             ! stored value at t=0 ( EFG tensor )
  real(kind=dp), dimension (:)     , allocatable :: vxy0 , vxz0 , vyz0             ! stored value at t=0 ( EFG tensor )
  real(kind=dp), dimension (:)     , allocatable :: v110 , v220 , v330 , eta0      ! stored value at t=0 ( principal components )
!  real(kind=dp), dimension (:,:)   , allocatable :: U0                             ! stored value at t=0 ( U vector )

  real(kind=dp), dimension (:,:)   , allocatable :: vxxt , vyyt , vzzt             ! value at any t      ( EFG tensor )   
  real(kind=dp), dimension (:,:)   , allocatable :: vxyt , vxzt , vyzt             ! value at any t      ( EFG tensor ) 
  real(kind=dp), dimension (:,:)   , allocatable :: v11t , v22t , v33t , etat      ! value at any t      ( principal components )
!  real(kind=dp), dimension (:,:,:) , allocatable :: Ut                             ! value at any t      ( U vector ) 
  real(kind=dp)                                  :: r_norm , nmr ( 4 ) , sq3 , sq32 , efgt(3,3)
  real(kind=dp)                                  :: w(3) 
  real(kind=dp)                                  :: work(3 * lwork)

  !trash
  integer            :: iiii
  character(len=200) :: XXXX

  ! ===================
  !   some constants
  ! ===================
  sq3 = SQRT ( 3.0_dp )
  sq3 = 1.0_dp / sq3
  sq32 = sq3 * 0.5_dp

  OPEN ( kunit_EFGALL,FILE='EFGALL')
  READ ( kunit_EFGALL , * ) natm
  READ ( kunit_EFGALL , * ) system
  READ ( kunit_EFGALL , * ) simu_cell%A(1,1) , simu_cell%A(2,1) , simu_cell%A(3,1)
  READ ( kunit_EFGALL , * ) simu_cell%A(1,2) , simu_cell%A(2,2) , simu_cell%A(3,2)
  READ ( kunit_EFGALL , * ) simu_cell%A(1,3) , simu_cell%A(2,3) , simu_cell%A(3,3)
  READ ( kunit_EFGALL , * ) ntype
  READ ( kunit_EFGALL , * )  ( atypei ( it ) , it = 1 , ntype )
  READ ( kunit_EFGALL , * )  ( natmi  ( it ) , it = 1 , ntype )
  READ ( kunit_EFGALL , * ) XXXX
  CALL lattice ( simu_cell ) 
  rho = natm / simu_cell%omega
  CALL config_alloc 
  CALL typeinfo_init

  CALL print_general_info ( stdout ) 
  WRITE ( stdout , '(a,e12.5)' ) 'time step ', dt

  allocate ( norm   ( 0:ntype , 0:ntcor )                                                             )
  allocate ( vxxt   ( natm    , ncefg   ) , vyyt  ( natm , ncefg )      , vzzt  ( natm , ncefg )      )
  allocate ( vxyt   ( natm    , ncefg   ) , vxzt  ( natm , ncefg )      , vyzt  ( natm , ncefg )      )
  allocate ( v11t   ( natm    , ncefg   ) , v22t  ( natm , ncefg )      , v33t  ( natm , ncefg )      )
  allocate ( etat   ( natm    , ncefg   )                                                             )
!  allocate ( Ut    ( natm , 5 , ncefg )                                                               )
  allocate ( vxx0   ( natm              ) , vyy0  ( natm )              , vzz0  ( natm )              )
  allocate ( vxy0   ( natm              ) , vxz0  ( natm )              , vyz0  ( natm )              )
  allocate ( v110   ( natm              ) , v220  ( natm )              , v330  ( natm )              )
  allocate ( eta0   ( natm              )                                                             )
!  allocate ( U0    ( natm , 5 )                                                                       )
  allocate ( acfxx  ( 0:ntype , 0:ntcor ) , acfyy ( 0:ntype , 0:ntcor ) , acfzz ( 0:ntype , 0:ntcor ) )
  allocate ( acfxy  ( 0:ntype , 0:ntcor ) , acfxz ( 0:ntype , 0:ntcor ) , acfyz ( 0:ntype , 0:ntcor ) )
  allocate ( acf11  ( 0:ntype , 0:ntcor ) , acf22 ( 0:ntype , 0:ntcor ) , acf33 ( 0:ntype , 0:ntcor ) )
  allocate ( acfeta ( 0:ntype , 0:ntcor )                                                             )
!  allocate ( acfU  ( 0:ntype , 5 , 0:ntcor )                           )
 
  acfxx  = 0.0_dp
  acfyy  = 0.0_dp
  acfzz  = 0.0_dp
  acfxy  = 0.0_dp
  acfxz  = 0.0_dp
  acfyz  = 0.0_dp
  acf11  = 0.0_dp
  acf22  = 0.0_dp
  acf33  = 0.0_dp
  acfeta = 0.0_dp
!  acfU   = 0.0_dp
  
  do t0 = 1 , ncefg 
    if ( t0 .ne. 1 ) then
      READ ( kunit_EFGALL , * ) natm
      READ ( kunit_EFGALL , * ) system
      READ ( kunit_EFGALL , * ) simu_cell%A(1,1) , simu_cell%A(2,1) , simu_cell%A(3,1)
      READ ( kunit_EFGALL , * ) simu_cell%A(1,2) , simu_cell%A(2,2) , simu_cell%A(3,2)
      READ ( kunit_EFGALL , * ) simu_cell%A(1,3) , simu_cell%A(2,3) , simu_cell%A(3,3)
      READ ( kunit_EFGALL , * ) ntype
      READ ( kunit_EFGALL , * )  ( atypei ( it ) , it = 1 , ntype )
      READ ( kunit_EFGALL , * )  ( natmi  ( it ) , it = 1 , ntype )
      READ ( kunit_EFGALL , * ) XXXX
    endif
    do ia = 1 , natm
      READ(kunit_EFGALL,*) iiii , atype(ia) , vxxt( ia , t0 ) , vyyt( ia , t0 ) , vzzt( ia , t0 ) , vxyt( ia , t0 ) , vxzt( ia , t0 ) , vyzt( ia , t0 )
    ! ===================================================================== 
    !       NMR parameters : principal components and assymetry
    ! ===================================================================== 
    efgt = 0.0_dp
    efgt ( 1 , 1 ) = vxxt( ia , t0 )
    efgt ( 2 , 2 ) = vyyt( ia , t0 )
    efgt ( 3 , 3 ) = vzzt( ia , t0 )
    efgt ( 1 , 2 ) = vxyt( ia , t0 )
    efgt ( 1 , 3 ) = vxzt( ia , t0 )
    efgt ( 2 , 3 ) = vyzt( ia , t0 )
#ifdef debug
    CALL print_tensor ( efgt , 'EFGTACFF' )
#endif
    ! ===================================================================== 
    !  Czjzek components (see J. Phys.: Condens. Matter 10 (1998). p10719)
    ! =====================================================================
    !Ut ( ia , 1 , t0 ) =   efgt ( 3 , 3 ) * 0.5_dp
    !Ut ( ia , 2 , t0 ) =   efgt ( 1 , 3 ) * sq3
    !Ut ( ia , 3 , t0 ) =   efgt ( 2 , 3 ) * sq3
    !Ut ( ia , 4 , t0 ) =   efgt ( 1 , 2 ) * sq3
    !Ut ( ia , 5 , t0 ) = ( efgt ( 1 , 1 ) - efgt ( 2 , 2 ) ) * sq32
    ! =================
    !  diagonalisation
    ! =================
    CALL DSYEV ( 'N' , 'U' , 3 , efgt , 3 , w , work , 3 * lwork , ifail )
    if ( ifail .ne. 0 ) then
      io_node WRITE ( stderr , * ) 'ERROR: DSYEV, STOP in efg MODULE (ewald)'
      STOP
    endif
#ifdef debug
    print*,w
#endif
    ! =============================
    ! NMR conventions test 
    ! =============================
    nmr = 0.0_dp
    CALL nmr_convention( w , nmr ( : )  , ia )
    v11t( ia , t0 ) = nmr ( 1 )
    v22t( ia , t0 ) = nmr ( 2 )
    v33t( ia , t0 ) = nmr ( 3 )
    etat( ia , t0 ) = nmr ( 4 )
#ifdef debug
    write ( stdout , '(4e12.5)') ,v11t( ia , t0 ),v22t( ia , t0 ),v33t( ia , t0 ),etat( ia , t0 )
#endif
    enddo

  enddo
  CLOSE(kunit_EFGALL)
 
  io_node WRITE ( stdout , '(a)' ) 'EFGALL successfully readed'

  do t0 = 1 , ncefg 
    if ( MOD ( t0 , ncefg / 20 ) .eq. 0 ) then
      WRITE ( stdout , '(a,i9)') 'time = ',t0
    endif
    do ia = 1 ,natm
      vxx0( ia )   = vxxt( ia , t0 )
      vyy0( ia )   = vyyt( ia , t0 )
      vzz0( ia )   = vzzt( ia , t0 )
      vxy0( ia )   = vxyt( ia , t0 )
      vxz0( ia )   = vxzt( ia , t0 )
      vyz0( ia )   = vyzt( ia , t0 )
      v110( ia )   = v11t( ia , t0 )
      v220( ia )   = v22t( ia , t0 )
      v330( ia )   = v33t( ia , t0 )
      eta0( ia )   = etat( ia , t0 )
  !    U0( ia , : ) = Ut  ( ia , : , t0 )
    enddo
    t0max = MIN ( ncefg , t0 + ntcor )
    do tt0 = t0 , t0max
      t = tt0 - t0
      do ia = 1 , natm
          it = itype ( ia ) 
          acfxx ( it , t ) = acfxx ( it , t ) + vxx0 ( ia ) * vxxt ( ia , tt0 )
          acfyy ( it , t ) = acfyy ( it , t ) + vyy0 ( ia ) * vyyt ( ia , tt0 )
          acfzz ( it , t ) = acfzz ( it , t ) + vzz0 ( ia ) * vzzt ( ia , tt0 )
          acfxy ( it , t ) = acfxy ( it , t ) + vxy0 ( ia ) * vxyt ( ia , tt0 )
          acfxz ( it , t ) = acfxz ( it , t ) + vxz0 ( ia ) * vxzt ( ia , tt0 )
          acfyz ( it , t ) = acfyz ( it , t ) + vyz0 ( ia ) * vyzt ( ia , tt0 )
          acf11 ( it , t ) = acf11 ( it , t ) + v110 ( ia ) * v11t ( ia , tt0 )
          acf22 ( it , t ) = acf22 ( it , t ) + v220 ( ia ) * v22t ( ia , tt0 )
          acf33 ( it , t ) = acf33 ( it , t ) + v330 ( ia ) * v33t ( ia , tt0 )
          acfeta( it , t ) = acfeta( it , t ) + eta0 ( ia ) * etat ( ia , tt0 )
  !        acfU  ( it ,  : , t ) = acfU  ( it , : , t ) + U0   ( ia , : ) * Ut   ( ia , : , tt0 )
          norm  ( it , t ) = norm  ( it , t ) + 1    
    !      acfxx ( 0  , t ) = acfxx ( 0  , t ) + vxx0 ( ia ) * vxxt ( ia , tt0 )
    !      acfyy ( 0  , t ) = acfyy ( 0  , t ) + vyy0 ( ia ) * vyyt ( ia , tt0 )
    !      acfzz ( 0  , t ) = acfzz ( 0  , t ) + vzz0 ( ia ) * vzzt ( ia , tt0 )
    !      acfxy ( 0  , t ) = acfxy ( 0  , t ) + vxy0 ( ia ) * vxyt ( ia , tt0 )
    !      acfxz ( 0  , t ) = acfxz ( 0  , t ) + vxz0 ( ia ) * vxzt ( ia , tt0 )
    !      acfyz ( 0  , t ) = acfyz ( 0  , t ) + vyz0 ( ia ) * vyzt ( ia , tt0 )
    !      acf11 ( 0  , t ) = acf11 ( 0  , t ) + v110 ( ia ) * v11t ( ia , tt0 )
    !      acf22 ( 0  , t ) = acf22 ( 0  , t ) + v220 ( ia ) * v22t ( ia , tt0 )
    !      acf33 ( 0  , t ) = acf33 ( 0  , t ) + v330 ( ia ) * v33t ( ia , tt0 )
    !      acfeta( 0  , t ) = acfeta( 0  , t ) + eta0 ( ia ) * etat ( ia , tt0 )
    !      acfU  ( 0  , : , t ) = acfU  ( 0  , : , t ) + U0   ( ia , : ) * Ut   ( ia , : , tt0 )
    !      norm  ( 0  , t ) = norm  ( 0  , t ) + 1
      enddo
    enddo      
  enddo

  OPEN(UNIT=kunit_EFGACFFF,FILE='EFGACFFF')
  OPEN(UNIT=kunit_NMRACFFF,FILE='NMRACFFF')

  io_node WRITE ( kunit_EFGACFFF , '(a,<ntype+1>a120)' )    '#     ',( atypei ( it ) , it = 1 , ntype ) 
  io_node WRITE ( kunit_EFGACFFF , '(a,<6*(ntype+1)>a20)' ) '#time ',(' vxx ',' vyy ',' vzz ',' vxy ',' vxz ',' vyz ', it = 1 , ntype ) 
  io_node WRITE ( kunit_NMRACFFF , '(a,<ntype+1>a120)' )    '#     ',( atypei ( it ) , it = 1 , ntype ) 
!  io_node WRITE ( kunit_NMRACFFF , '(a,<6*(ntype+1)>a20)' ) '#time ',(' VXX ',' VYY ',' VZZ ',' ETA ',' U1 ',' U2 ', ' U3 ',' U4 ',' U5 ', it = 1 , ntype ) 
  io_node WRITE ( kunit_NMRACFFF , '(a,<6*(ntype+1)>a20)' ) '#time ',(' VXX ',' VYY ',' VZZ ',' ETA ', it = 1 , ntype ) 
  do t = 0 , ntcor
    do it = 1 , ntype
      r_norm = 1.0_dp / REAL ( norm ( it , t ) , kind = dp )
      acfxx ( it , t ) = acfxx ( it , t ) * r_norm 
      acfyy ( it , t ) = acfyy ( it , t ) * r_norm 
      acfzz ( it , t ) = acfzz ( it , t ) * r_norm 
      acfxy ( it , t ) = acfxy ( it , t ) * r_norm 
      acfxz ( it , t ) = acfxz ( it , t ) * r_norm 
      acfyz ( it , t ) = acfyz ( it , t ) * r_norm 
      acf11 ( it , t ) = acf11 ( it , t ) * r_norm 
      acf22 ( it , t ) = acf22 ( it , t ) * r_norm 
      acf33 ( it , t ) = acf33 ( it , t ) * r_norm 
      acfeta( it , t ) = acfeta( it , t ) * r_norm 
   !   acfU  ( it , : , t ) = acfU  ( it , : , t ) * r_norm 
    enddo
    io_node WRITE ( kunit_EFGACFFF , '(<6*(ntype+1)+1>f20.10)' ) REAL ( t , kind = dp ) * dt , ( acfxx(it,t) , acfyy(it,t) , acfzz(it,t) , acfxy (it,t) ,  acfxz(it,t)  , acfyz(it,t)  , it = 1 , ntype ) 
!    io_node WRITE ( kunit_NMRACFFF , '(<9*(ntype+1)+1>f20.10)' ) REAL ( t , kind = dp ) * dt , ( acf11(it,t) , acf22(it,t) , acf33(it,t) , acfeta(it,t) ,( acfU(it,ui,t), ui = 1 , 5 ) , it = 1 , ntype ) 
    io_node WRITE ( kunit_NMRACFFF , '(<9*(ntype+1)+1>f20.10)' ) REAL ( t , kind = dp ) * dt , ( acf11(it,t) , acf22(it,t) , acf33(it,t) , acfeta(it,t) , it = 1 , ntype ) 
  enddo

  CLOSE ( kunit_NMRACFFF )      
  CLOSE ( kunit_EFGACFFF )      

  deallocate ( norm )
  deallocate ( acfxx , acfyy , acfzz )
  deallocate ( acfxy , acfxz , acfyz )
  deallocate ( vxxt  , vyyt  , vzzt  )
  deallocate ( vxyt  , vxzt  , vyzt  )
  deallocate ( vxx0  , vyy0  , vzz0  )
  deallocate ( vxy0  , vxz0  , vyz0  )
  deallocate ( v11t  , v22t  , v33t  )
  deallocate ( v110  , v220  , v330  )
  deallocate ( acf11 , acf22 , acf33 )
  deallocate ( eta0                  )
  deallocate ( etat                  )
  deallocate ( acfeta                )
!  deallocate ( U0                    )
!  deallocate ( Ut                    )
!  deallocate ( acfU                  )

  return
        
END SUBROUTINE efg_acf

! *********************** SUBROUTINE efg_stat **********************************
!
! Statitics on EFG tensors ( principal components, U vector and other
! quantities )
! Calculate distribution arrays of main EFG quantities 
!
! It reads EFGALL files and print EFGFF and NMRFF files
! 
! EFGFF is not really useful
!
! ******************************************************************************
SUBROUTINE efg_stat ( kunit_input , kunit_output , kunit_nmroutput )

  USE config,                   ONLY :  system , natm , natmi , ntype , itype , &
                                        atype , atypei , simu_cell , rho, config_alloc , quadia
  USE io_file,                  ONLY :  ionode , stdout , stderr 
  USE field,                    ONLY :  lwfc , field_init        
  USE constants,                ONLY :  eh 
  USE cell,                     ONLY :  lattice

  implicit none

  integer, parameter :: lwork = 6
  integer            :: ic , ia , it , ui
  integer            :: ifail
  integer            :: kvzz, keta ,ku
  integer            :: kunit_input , kunit_output , kunit_nmroutput

  real(kind=dp)                                 :: sq3 , sq32
  real(kind=dp)                                 :: w(3) , efgt(3,3) , nmr_conv ( 4 )
  real(kind=dp)                                 :: work(3 * lwork)
  real(kind=dp)                                 :: vzzk , etak , uk
  real(kind=dp)                                 :: CQ_UNIT
  real(kind=dp), dimension (:,:)  , allocatable :: U
  real(kind=dp), dimension (:,:)  , allocatable :: nmr
  real(kind=dp), dimension (:)    , allocatable :: vzzmini , vzzmaxi ! ntype
  real(kind=dp), dimension (:)    , allocatable :: vzzm, vzzma , etam , pvzz, rho_z, sigmavzz, vzzsq ! ntype

  !trash
  integer   :: iiii
  character :: XXXX

  ! some constants
  sq3 = SQRT ( 3.0_dp )
  sq3 = 1.0_dp / sq3
  sq32 = sq3 * 0.5_dp

#ifdef debug_efg_stat
  WRITE ( stdout , '(a)' ) 'debug : in efg_stat'
#endif

  READ ( kunit_input , * )  natm
  READ ( kunit_input , * )  system
  READ ( kunit_input , * )  simu_cell%A(1,1) , simu_cell%A(2,1) , simu_cell%A(3,1)
  READ ( kunit_input , * )  simu_cell%A(1,2) , simu_cell%A(2,2) , simu_cell%A(3,2)
  READ ( kunit_input , * )  simu_cell%A(1,3) , simu_cell%A(2,3) , simu_cell%A(3,3)
  READ ( kunit_input , * )  ntype
  READ ( kunit_input ,* ) ( atypei ( it ) , it = 1 , ntype )
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
  if ( lefg_restart ) CALL field_init 
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
  nmr      = 0.0_dp 
  pvzz     = 0.0_dp
  vzzmini  = 0.0_dp
  vzzmaxi  = 0.0_dp
  etam     = 0.0_dp
  vzzm     = 0.0_dp
  vzzma    = 0.0_dp
  vzzsq    = 0.0_dp
  U        = 0.0_dp

  do ic = 1 , ncefg
    if ( ic .ne. 1 ) then
      READ ( kunit_input , * )  natm
      READ ( kunit_input , * )  system
      READ ( kunit_input , * )  simu_cell%A(1,1) , simu_cell%A(2,1) , simu_cell%A(3,1)
      READ ( kunit_input , * )  simu_cell%A(1,2) , simu_cell%A(2,2) , simu_cell%A(3,2)
      READ ( kunit_input , * )  simu_cell%A(1,3) , simu_cell%A(2,3) , simu_cell%A(3,3)
      READ ( kunit_input , * )  ntype
      READ ( kunit_input , * )  ( xxxx , it = 1 , ntype )
      READ ( kunit_input , * )  ( iiii , it = 1 , ntype )
      READ ( kunit_input , * )  xxxx
    endif
    CALL lattice ( simu_cell )
    rho = natm / simu_cell%omega


    do ia = 1 , natm
      it = itype ( ia )
      if ( lwfc ( it ) .ge. 0 ) then
        READ( kunit_input , * ) iiii , xxxx , efgt(1,1) , efgt(2,2) , efgt(3,3) , &
                                              efgt(1,2) , efgt(1,3) , efgt(2,3)
      endif

#ifdef debug_efg_stat
  WRITE ( stdout , '(a,2i8,a,4f12.6)' ) 'debug : reading tensor',ia, iiii , xxxx , efgt(1,1) , efgt(2,2) , efgt(3,3),efgt(1,1) + efgt(2,2) + efgt(3,3) 
#endif
      ! ===================================================================== 
      !  Czjzek components (see J. Phys.: Condens. Matter 10 (1998). p10719)
      ! =====================================================================
      U ( ia , 1 ) = efgt ( 3 , 3 ) * 0.5_dp
      U ( ia , 2 ) = efgt ( 1 , 3 ) * sq3
      U ( ia , 3 ) = efgt ( 2 , 3 ) * sq3
      U ( ia , 4 ) = efgt ( 1 , 2 ) * sq3
      U ( ia , 5 ) = ( efgt ( 1 , 1 ) - efgt ( 2 , 2 ) ) * sq32
  
!#ifdef debug_efg_stat
! STOP
!#endif
      ! =================
      !  diagonalisation
      ! =================
      CALL DSYEV ( 'N' , 'U' , 3 , efgt , 3 , w , work , 3 * lwork , ifail )
      if ( ifail .ne. 0 ) then
        io_node WRITE ( stderr , * ) 'ERROR: DSYEV, STOP in efg MODULE (ewald)'
        STOP
      endif

      ! =============================
      ! NMR conventions test 
      ! =============================
      CALL nmr_convention( w , nmr_conv ( : )  , ia )
      nmr( ia , : )  = nmr_conv ( : )
    
!#ifdef debug_efg_stat
! STOP
!#endif
      ! ===================================================
      !  statistic for all sites index 0 in array of size ntype
      ! ===================================================
      vzzm  ( 0 )  =  vzzm ( 0 ) +       nmr ( ia , 3 )
      etam  ( 0 )  =  etam ( 0 ) +       nmr ( ia , 4 )
      vzzma ( 0 )  =  vzzma( 0 ) + ABS ( nmr ( ia , 3 ) )
      if ( nmr ( ia , 3 ) .le. vzzmini ( 0 ) )  vzzmini ( 0 ) = nmr ( ia , 3 )
      if ( nmr ( ia , 3 ) .ge. vzzmaxi ( 0 ) )  vzzmaxi ( 0 ) = nmr ( ia , 3 )
      if ( nmr ( ia , 3 ) .ge. 0.0_dp ) pvzz ( 0 ) = pvzz ( 0 ) + 1.0_dp
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
          if (nmr(ia,3).ge.0.0_dp ) pvzz(it) = pvzz(it) + 1.0_dp
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
      if ( ( ionode ) ) &
      WRITE ( stdout ,100) & 
      ic , atypei(it) , vzzmini(it) , vzzmaxi(it) , vzzm(it) , vzzma(it) , etam(it) , pvzz(it) , rho_z ( it)
      if (   ionode ) WRITE ( kunit_output , 100) &
      ic , atypei(it) , vzzmini(it) , vzzmaxi(it) , vzzm(it) , vzzma(it) , etam(it) , pvzz(it) , rho_z(it)
    enddo
    if ( ionode .and. ntype .ne. 1 ) &
      WRITE ( stdout ,100) &
      ic , atypei(0) , vzzmini(0) , vzzmaxi(0) , vzzm(0) , vzzma(0) , etam(0) , pvzz(0) , rho_z(0)
    if ( ionode .and. ntype .ne. 1 ) &
      WRITE ( kunit_output ,100) &
      ic , atypei(0) , vzzmini(0) , vzzmaxi(0) , vzzm(0) , vzzma(0) , etam(0) , pvzz(0) , rho_z(0)

    ! ========================================
    !  output NMR  parameters 
    ! ========================================
    CQ_UNIT = eh / 1000.0_dp  
    if ( ionode ) then
        WRITE ( kunit_nmroutput , '(a)' ) &
        '#     site           xx          yy          zz           Cq(MHz)               eta'
      do ia = 1 , natm
        it = itype ( ia )   
        if ( lwfc ( it ) .ge. 0 ) then
        WRITE ( kunit_nmroutput , 150 ) ia , atype ( ia ) , & 
                                        nmr( ia , 1 ) , nmr( ia , 2 ) , &
                                        nmr( ia , 3 ) , nmr( ia , 3 ) * CQ_UNIT * quadia ( ia ) , nmr( ia , 4 )
        endif
      enddo
    endif

   if ( .not. lefg_stat ) then
     io_node WRITE ( stdout , '(a)' ) "No statistic on EFG's"
     return
   endif

    ! =======================================================
    !  DISTRIBUTION CALCULATION LOOP 
    ! =======================================================
    do ia = 1 , natm
      it = itype ( ia )
      if ( lwfc ( it ) .eq. -1 ) cycle
      ! ========
      !    Ui
      ! ========
      do ui=1,5
        uk = (U(ia,ui)-umin )/resu
        ku = int(uk) + 1
        ! ====================== 
        !  test out of bound
        ! ====================== 
        if (ku.lt.0.or.ku.gt.PANU) then
          io_node WRITE ( stderr , * ) 'ERROR: out of bound dibU1'
          io_node WRITE ( stderr ,310) ia,ku,U(ia,ui),umin,ABS (umin)
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
      kvzz = int(vzzk) + 1
      keta = int(etak)
      ! ====================== 
      !  test out of bound
      ! ====================== 
      if ( kvzz .lt. 0 .or. kvzz .gt. PANvzz ) then
        if ( ionode .and. itype(ia) .eq. 1 ) &
        WRITE ( stderr , '(a,a)' ) 'ERROR: out of bound distribvzz ',atype(ia)
        if ( ionode .and. itype(ia) .eq. 2 ) &
        WRITE ( stderr , '(a,a)' ) 'ERROR: out of bound distribvzz ',atype(ia)
        io_node WRITE ( stderr ,200) &
        ia , kvzz , nmr ( ia , 3 ) , vzzmini( itype ( ia ) ) , vzzmaxi ( itype ( ia ) ) , &
                    nmr ( ia , 4 ) , nmr ( ia , 1 ) , nmr ( ia , 2 ) , nmr ( ia , 3 )
        STOP
      endif
      if ( keta .lt. 0 .or. keta .gt. PANeta ) then
        io_node WRITE ( stderr , '(a)' ) 'ERROR: out of bound distribeta'
        io_node WRITE ( stderr ,210) ia, keta, nmr(ia,4), nmr(ia,4), &
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

100 FORMAT(I7,1X,' ATOM ',A3,'  minVZZ = ',F10.5,' maxVZZ = ',F10.5,& 
             ' <VZZ> =  ',F10.5,' <|VZZ|> =  ',F10.5,' mETA =  ',F10.5, &
          ' P(Vzz>0) =  ',F10.5, ' RHO_Z = ',F10.5)

150 FORMAT(I7,1X,A3,3F12.4,2F18.4)

200 FORMAT('atom = ',I6,' k = ',I12, 'value = ',F14.8,&
           ' min =  ',F10.5,' max =  ',F10.5,' vaa{a = xx,yy,zz}, eta = ',4F14.8)
210 FORMAT('atom = ',I6,' k = ',I12, 'value = ',F14.8,& 
           ' min =  0.0_dp    max =  1.0_dp',' vaa{a = xx,yy,zz}, eta = ',4F14.8)

310 FORMAT('atom = ',I6,' k = ',I12, 'value = ',F14.8,' min =  ',F10.5,' max =  ',F10.5)

#ifdef debug_efg_stat 
      lseparator(sdterr) 
      WRITE(stderr, '(a)' ) 'check distribution counting (sum arrays)'
      lseparator(sdterr) 
      do it = 0 , ntype 
        WRITE(stderr, '(a,2i9)' ) 'number of atoms for type (it=0 => ALL) : ',it,natmi(it)
      enddo
      do it = 0 , ntype 
        WRITE(stderr, '(a,i6)' ) 'it = ',it
        WRITE(stderr, '(a,a)' )    '     ','  actual sum      expected sum'
        WRITE(stderr, '(a,2i12)' ) 'Vzz  ' , SUM(dibvzztot(it,:)),ncefg*natmi(it)
        WRITE(stderr, '(a,2i12)' ) 'eta  ' , SUM(dibetatot(it,:)),ncefg*natmi(it)
        WRITE(stderr, '(a,2i12)' ) 'U1   ' , SUM(dibUtot(1,it,:)),ncefg*natmi(it)
        WRITE(stderr, '(a,2i12)' ) 'Uk>1 ' , SUM(dibUtot(6,it,:)),ncefg*natmi(it)*4
        lseparator(sdterr) 
      enddo
#endif


  return

END SUBROUTINE efg_stat

! *********************** SUBROUTINE nmr_convention ****************************
!
! This routine determine the order of principal components from the  
! NMR convention: |Vzz| > |Vyy| > |Vxx|
! the asymmetry is then define from 
!
!          Vxx - Vyy
!  eta = ------------
!             Vzz
!
!  In input EIG are the eigenvectors of the EFG tensor
! 
!  In output the array nmr(4) gives :
!    nmr ( 1 ) = Vxx
!    nmr ( 2 ) = Vyy
!    nmr ( 3 ) = Vzz
!    nmr ( 4 ) = eta
! 
! ******************************************************************************
SUBROUTINE nmr_convention( EIG , nmr , ia )

  USE io_file,                  ONLY :  stdout , stderr , ionode

  implicit none

  ! global
  real(kind=dp) ::  EIG(3)
  real(kind=dp) ::  nmr(4)
  integer :: ia
  ! local 
  real(kind=dp) ::  del11 , del22 , del33 , diso

  ! =======================================================
  ! =======================================================
  diso=(EIG(1)+EIG(2)+EIG(3))/3.0_dp
  if ( ABS ( diso ) .ne. 0.0_dp .and. ABS ( diso ) .gt. 1e-3 ) then
!    io_node WRITE ( stderr ,'(a,i6,f48.24)') 'ERROR: trace of EFG is not null',ia,diso
     diso = 0.0_dp
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
  if ( nmr(1) .eq. nmr(2) )  nmr(4) = 0.0_dp
  if ( nmr(3) .eq. 0.0 )  nmr(4) = 0.0_dp
  if ( ABS ( nmr(1) ) .lt. 1.0E-05 .and. ABS ( nmr(2) ) .lt. 1.0E-5 .and. ABS ( nmr(3) ) .lt. 1.0E-5 ) nmr(4) = 0.0_dp
  if (ABS (nmr(3)- diso ).lt.ABS (nmr(1) -diso)) then
    io_node WRITE ( stderr , '(a,i4,3e24.16)' ) 'ERROR: |Vzz| < |Vxx|',ia,nmr(1),nmr(2),nmr(3)
    STOP
  endif
  if (ABS (nmr(3) -diso ).lt.ABS (nmr(2) -diso )) then
    io_node WRITE ( stderr , '(a,i4,3e24.16)' ) 'ERROR: |Vzz| < |Vyy|',ia,nmr(1),nmr(2),nmr(3)
    STOP
  endif
  if (ABS (nmr(2) -diso ).lt.ABS (nmr(1)-diso)) then
    io_node WRITE ( stderr , '(a,i4,3e24.16)' ) 'ERROR: |Vyy| < |Vxx|',ia,nmr(1),nmr(2),nmr(3)
    STOP
  endif
  if ( nmr(4) .gt. 1.0_dp ) then
    !if ( ionode ) &
    !WRITE ( stderr ,'(a,i6,5f24.16)') 'ERROR1: eta > 1.0_dp ',ia,nmr(4),ABS(nmr(1)-diso),ABS(nmr(2)-diso),ABS(nmr(3)-diso)
    !WRITE ( stderr ,'(a,i6,5f24.16)') 'ERROR2: eta > 1.0_dp ',ia,nmr(4),nmr(1),nmr(2),nmr(3)
    !STOP
  !  nmr(4) = 1.0_dp
  endif
  if ( nmr(4) .gt. 1.0_dp .or. nmr(4) .lt. 0.0_dp) then
    !if ( ionode ) &
    !WRITE ( stderr ,'(a,i6,4f24.16)') 'ERROR: eta < 0.0_dp',ia,nmr(4),nmr(1),nmr(2),nmr(3)
  !  nmr(4) = 0.0_dp
    !STOP
  endif

  return

END SUBROUTINE nmr_convention


! *********************** SUBROUTINE efg_alloc *********************************
!
! allocate quantities related to efg calculation
! /deallocate 
!
! ******************************************************************************
SUBROUTINE efg_alloc 

  USE control,                  ONLY :  calc , longrange
  USE config,                   ONLY :  natm , ntype , qia , itype

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
#endif
  allocate( efg_t    ( natm , 3 , 3 ) )
  allocate( efg_ia    ( natm , 3 , 3 ) )
  allocate( mu ( natm , 3 ) ) 
  mu = 0.0_dp
  efg_t = 0.0_dp
  efg_ia = 0.0_dp

  return

END SUBROUTINE efg_alloc

! *********************** SUBROUTINE efg_dealloc *******************************
!
! deallocate quantities related to efg calculation
!
! ******************************************************************************
SUBROUTINE efg_dealloc

  USE control,                  ONLY :  calc , longrange

  implicit none

  if ( calc .ne. 'efg' ) return

  deallocate( dibUtot )
  deallocate( dibvzztot )
  deallocate( dibetatot )
#ifdef fix_grid
  deallocate ( rgrid )
#endif
  deallocate( efg_t )
  deallocate( efg_ia    )
  deallocate( mu ) 

  return

END SUBROUTINE efg_dealloc


! *********************** SUBROUTINE efg_mesh_alloc ****************************
!
! allocate quantities related to the real space or reciprocal space mesh used to
! calculate efg's
! /deallocate 
!
! ******************************************************************************
SUBROUTINE efg_mesh_alloc

  USE control,                  ONLY :  calc , longrange
  USE config,                   ONLY :  natm , ntype , qia , itype
  USE field,                    ONLY :  km_coul , rm_coul , ncelldirect , kES
  USE kspace,                   ONLY :  kpoint_sum_init
  USE rspace,                   ONLY :  direct_sum_init

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
    allocate ( km_coul%strf ( nkcut, natm) )
    CALL kpoint_sum_init ( km_coul )
  endif

  return

END SUBROUTINE efg_mesh_alloc

! *********************** SUBROUTINE efg_mesh_dealloc **************************
!
! allocate quantities related to the real space or reciprocal space mesh used to
! calculate efg's
! /deallocate 
!
! ******************************************************************************
SUBROUTINE efg_mesh_dealloc

  USE control,                  ONLY :  calc , longrange
  USE field,                    ONLY :  km_coul , rm_coul 

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
