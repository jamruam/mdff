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
! along with this program; if not, WRITE to the Free Software
! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
! ===== fmV =====

! ======= Hardware =======
#include "symbol.h"
!#define debug
!#define debug_units
!#define debug_multipole_ES
!#define debug2
!#define debug_wfc
!#define debug_morse
!#define debug_bmlj
!#define debug_bmlj_pbc
!#define debug_quadratic
!#define debug_para
! ======= Hardware =======

! *********************** MODULE field *****************************************
!> \brief
!! Module related to force-field calculation
! ******************************************************************************
MODULE field 

  USE constants,                        ONLY :  dp 
  USE config,                           ONLY :  ntypemax 
  USE kspace,                           ONLY :  kmesh 
  USE rspace,                           ONLY :  rmesh
  USE io,                          ONLY :  ionode
  USE mpimdff

  implicit none

  integer :: cccc=0
  real(kind=dp)     :: utail               !< long-range correction of short-range interaction 
  character(len=60) :: units               !< energy units currently used 
  character(len=60) :: units_allowed(5)
  data                 units_allowed / 'eV' , 'K' , 'kcal' , 'kJ' , 'internal' /  ! see define_units 
  character(len=60) :: ctrunc                                                     !< truncation of bmlj
  character(len=60) :: ctrunc_allowed(3)                                          !< truncation of bmlj 
  data                 ctrunc_allowed / 'notrunc', 'linear' , 'quadratic' /       !< see initialize_param_bmlj
  integer           :: trunc                                                      !< integer definition of truncation 

  logical, SAVE     :: lKA               !< use Kob-Andersen model for BMLJ                        
  logical, SAVE     :: lautoES           !< auto-determination of Ewald parameter from epsw ( accuracy)
  logical, SAVE     :: lwrite_dip_wfc    !< write dipoles from wannier centers to file
  logical, SAVE     :: ldip_wfc          !< calculate electrostatic contribution from dipolar momemt coming from wfc
  logical, SAVE     :: lquiet            !< unknown


  ! ============================================================  
  !                         Lennard - Jones
  ! ============================================================  
  !
  !              eps    /    / sigma*\ q         / sigma*\ p  \
  !     V  =   ------- |  p | ------- |   -  q  | ------- |    |      sigma* = 2^(1/6)*sigma
  !             q - p   \    \   r   /           \   r   /    /
  !
  ! main parameters
  real(kind=dp) :: qlj     ( ntypemax , ntypemax )
  real(kind=dp) :: plj     ( ntypemax , ntypemax )
  real(kind=dp) :: epslj   ( ntypemax , ntypemax )
  real(kind=dp) :: sigmalj ( ntypemax , ntypemax )
  real(kind=dp) :: rcutsq  ( ntypemax , ntypemax )  
  real(kind=dp) :: sigsq   ( ntypemax , ntypemax ) 
  real(kind=dp) :: epsp    ( ntypemax , ntypemax )
  real(kind=dp) :: fc      ( ntypemax , ntypemax )
  real(kind=dp) :: uc      ( ntypemax , ntypemax )
  real(kind=dp) :: uc1     ( ntypemax , ntypemax )
  real(kind=dp) :: uc2     ( ntypemax , ntypemax )
  real(kind=dp) :: testtab ( ntypemax , ntypemax )

  ! ============================================================  
  !                            Morse
  ! ============================================================  
  !
  !     V  =  
  !
  !
  ! main parameters
  real(kind=dp) :: rhomor  ( ntypemax , ntypemax )
  real(kind=dp) :: epsmor  ( ntypemax , ntypemax )
  real(kind=dp) :: sigmamor( ntypemax , ntypemax )
  real(kind=dp) :: rs      ( ntypemax , ntypemax )
  real(kind=dp) :: fm      ( ntypemax , ntypemax )


  ! ============================================================  
  !                    FORCE FIELD 
  ! ============================================================  

  ! type dependent properties
  real(kind=dp)    :: mass     ( ntypemax )            !< masses ( not yet tested one everywhere )
  real(kind=dp)    :: qch      ( ntypemax )            !< charges 
  real(kind=dp)    :: quad_efg ( ntypemax )            !< quadrupolar moment
  real(kind=dp)    :: dip      ( ntypemax , 3 )        !< dipoles 
  real(kind=dp)    :: pol      ( ntypemax , 3 , 3 )    !< polarizability if lpolar( it ) .ne. 0  
  integer          :: lpolar   ( ntypemax )            !< induced moment from pola 
  integer          :: lwfc     ( ntypemax )            !< moment from wannier centers 
  real(kind=dp)    :: conv_tol_ind                     !< convergence tolerance of the scf induced dipole calculation
  real(kind=dp)    :: rcut_wfc                         !< radius cut-off for WFs searching
  
  ! ewald sum related 
  real(kind=dp)    :: epsw                             !< accuracy of the ewald sum 
  real(kind=dp)    :: alphaES                          !< Ewald sum parameter 
  real(kind=dp)    :: cutshortrange                    !< Ewald sum parameter cutoff shortrange 
  integer          :: kES(3)                           !< kmax of ewald sum in reciprocal space
  TYPE ( kmesh )   :: km_coul                          !< kpoint mesh ( see kspace.f90 )
  ! direct sum
  integer          :: ncelldirect                      !< number of cells  in the direct summation
  TYPE ( rmesh )   :: rm_coul                          !< real space mesh ( see rspace.f90 )

  
  real(kind=dp), dimension ( : , : )     , allocatable :: ef_t  !< electric field vector
  real(kind=dp), dimension ( : , : , : ) , allocatable :: efg_t !< electric field gradient tensor


CONTAINS

! *********************** SUBROUTINE field_default_tag *************************
!> \brief
!! set default values to field tags
! ******************************************************************************
SUBROUTINE field_default_tag

  implicit none

  ! =================
  !  default values
  ! =================

  units         = 'internal'

  ! BMLJ
  lKA           = .false.
  ctrunc        = 'linear' 
  epslj         = 1.0_dp
  sigmalj       = 1.0_dp
  epsmor        = 1.0_dp 
  sigmamor      = 1.0_dp 
  rhomor        = 3.0_dp 
  qlj           = 12.0_dp
  plj           = 6.0_dp
  mass          = 1.0_dp
  ! Coulomb
  ncelldirect   =  2
  kES(1)        = 10
  kES(2)        = 10
  kES(3)        = 10
  alphaES       = 1.0_dp
  qch           = 0.0_dp
  quad_efg      = 0.0_dp
  dip           = 0.0_dp
  pol           = 0.0_dp
  lpolar        = 0
  lwfc          = 0
  lwrite_dip_wfc= .false.            
  ldip_wfc      = .true.            
  rcut_wfc      = 0.5_dp

  conv_tol_ind  = 1e-5
  epsw = 1e-7
  lautoES = .false.

  return

END SUBROUTINE field_default_tag


! *********************** SUBROUTINE field_check_tag ***************************
!> \brief
!! check field tag values
! ******************************************************************************
SUBROUTINE field_check_tag

  USE config,                   ONLY :  ntype
  USE io,                  ONLY :  ionode , stdout

  implicit none

  integer :: i
  logical :: allowed

  allowed = .false.
  ! ========
  !  units 
  ! ========
  do i = 1 , size ( units_allowed )
    if ( trim ( units ) .eq. units_allowed ( i ) )  allowed = .true.
  enddo

  if ( .not. allowed ) then
    if ( ionode )  WRITE ( stdout , '(a)' ) 'ERROR fieldtag: units should be ', units_allowed
    STOP
  endif
  allowed = .false.

  ! ========
  !  ctrunc
  ! ========
  do i = 1 , size ( ctrunc_allowed )
    if ( trim ( ctrunc ) .eq. ctrunc_allowed ( i ) )  allowed = .true.
  enddo

  if ( .not. allowed ) then
    if ( ionode )  WRITE ( stdout , '(a)' ) 'ERROR fieldtag: ctrunc should be ', ctrunc_allowed
    STOP
  endif

  if ( ctrunc .eq. 'notrunc'   ) trunc = 0
  if ( ctrunc .eq. 'linear'    ) trunc = 1
  if ( ctrunc .eq. 'quadratic' ) trunc = 2

  ! =================================================
  !  KOB-ANDERSEN MODEL --- PhysRevE 51-4626 (1995) 
  ! =================================================
  if ( lKA .and. ntype .ne. 2 ) then
    io_node WRITE ( stdout , '(a)' ) 'ERROR fieldtag lKA should be used with 2 differents types'
    STOP 
  endif
  if ( lKA ) then
    sigmalj ( 1 , 1 ) = 1.0_dp
    sigmalj ( 2 , 2 ) = 0.88_dp
    sigmalj ( 1 , 2 ) = 0.8_dp
    sigmalj ( 2 , 1 ) = 0.8_dp
    epslj   ( 1 , 1 ) = 1.0_dp
    epslj   ( 2 , 2 ) = 0.5_dp
    epslj   ( 1 , 2 ) = 1.5_dp
    epslj   ( 2 , 1 ) = 1.5_dp
  endif

  return 

END SUBROUTINE field_check_tag

! *********************** SUBROUTINE field_init ********************************
!> \brief
!! force field initialisation
! ******************************************************************************
SUBROUTINE field_init

  USE control,                  ONLY :  calc , lbmlj , lcoulomb , lmorse , longrange
  USE io,                  ONLY :  ionode, stdin, stdout 
  USE config,                   ONLY :  dipia

  implicit none

  ! local
  character(len=132)    :: filename
  integer               :: ioerr

  namelist /fieldtag/    lKA           , &       
                         ctrunc        , &
                         units         , &
                         ncelldirect   , &
                         kES           , &
                         alphaES       , &
                         qlj           , & 
                         plj           , & 
                         sigmalj       , &
                         epslj         , &
                         sigmamor      , &
                         epsmor        , &
                         rhomor        , &
                         mass          , &
                         qch           , &
                         quad_efg      , &
                         dip           , &
                         dipia         , &
                         pol           , &  
                         conv_tol_ind  , &  
                         epsw          , &  
                         lautoES       , &  
                         lwfc          , &            
                         lwrite_dip_wfc, &            
                         ldip_wfc      , &            
                         rcut_wfc      , &            
                         lpolar           
                 

  ! ================================
  ! defaults values for field tags 
  ! ================================
  CALL field_default_tag

  ! ================================
  ! read field tags
  ! ================================
  CALL getarg (1,filename)
  OPEN ( stdin , file = filename)
    READ ( stdin , fieldtag, iostat=ioerr)
    if ( ioerr .lt. 0 )  then
      io_node WRITE ( stdout, '(a)') 'ERROR reading input_file : fieldtag section is absent'
      STOP
    elseif ( ioerr .gt. 0 )  then
      io_node WRITE ( stdout, '(a,i8)') 'ERROR reading input_file : fieldtag wrong tag'
      STOP
    endif
  CLOSE ( stdin )

  ! ================================
  ! check field tags values
  ! ================================
  CALL field_check_tag

  CALL define_units

  ! ===============================================
  !  this routines generates the ewald parameters
  ! ===============================================
  if ( longrange .eq. 'ewald' ) CALL ewald_param
 
  ! =====================================  
  !  if efg print field info and return 
  ! =====================================  
  if ( calc .eq. 'efg' .or. calc .eq. 'rmc' ) then 
    CALL field_print_info(stdout,quiet=.true.)
    return
  endif

  ! ================================
  ! initialize constant parameters
  ! ================================
  if ( lbmlj .or. lmorse )    then
    CALL initialize_param_bmlj_morse
  endif

  if ( lcoulomb ) then
    CALL initialize_coulomb
  endif

  ! ================================
  !  pint field info
  ! ================================
  CALL field_print_info(stdout,quiet=.false.)

  return

END SUBROUTINE field_init

! *********************** SUBROUTINE field_print_info **************************
!> \brief
!! print force field information to standard output
! ******************************************************************************
SUBROUTINE field_print_info ( kunit , quiet )

  USE config,           ONLY :  natm , ntype , atypei , natmi , simu_cell 
  USE control,          ONLY :  calc , cutshortrange , lbmlj , lmorse , lcoulomb , longrange , lreduced , cutlongrange
  USE io,          ONLY :  ionode 
  USE constants,        ONLY :  pi , pisq

  implicit none

  !local
  logical , optional :: quiet
  integer :: kunit, it , it1 , it2 , i , j
  real(kind=dp) :: rcut2 , kmax2 , alpha2 , ereal , ereci(3) , ereci2(3) , qtot , qtot2
  logical :: linduced

  if ( ( present ( quiet ) .and. quiet ) .and. .not. lquiet ) then 
    lquiet = .true.
  else if ( ( present ( quiet ) .and. quiet ) .and. lquiet ) then
    return
  endif

  qtot   = 0.0_dp
  qtot2  = 0.0_dp
  do it = 1 , ntype
      qtot  = qtot  +   qch(it) * natmi ( it )
      qtot2 = qtot2 + ( qch(it) * natmi ( it ) ) * ( qch(it) * natmi ( it ) )
  enddo
  linduced = .false.
  do it = 1 , ntype
    if ( lpolar(it) .eq. 1 )  linduced = .true.
  enddo

  if ( ionode ) then
    separator(kunit)    
    blankline(kunit)
    WRITE ( kunit ,'(a)')               'FIELD MODULE ... WELCOME'
    blankline(kunit)
    WRITE ( kunit ,'(a)')               'energy units in input/output ',units 
    blankline(kunit)
    lseparator(kunit) 
    WRITE ( kunit ,'(a)')               'point charges: '
    lseparator(kunit) 
    do it = 1 , ntype 
      WRITE ( kunit ,'(a,a,a,e10.3)')   'q'   ,atypei(it),'                    = ',qch(it)
      WRITE ( kunit ,'(a,a,a,e10.3,a)') 'quad',atypei(it),'                 = ',quad_efg(it),' mb'
    enddo
    WRITE ( kunit ,'(a,e10.3)')         'total charge            = ',  qtot
    WRITE ( kunit ,'(a,e10.3)')         'second moment of charge = ',  qtot2
    blankline(kunit)
    lseparator(kunit) 
    WRITE ( kunit ,'(a)')               'static dipoles: '
    lseparator(kunit) 
    do it = 1 , ntype
      WRITE ( kunit ,'(a,a,a,3f10.5)')  'mu',atypei(it),'      = ',dip(it,1),dip(it,2),dip(it,3)
    enddo
    blankline(kunit)
    if ( linduced ) then 
      lseparator(kunit) 
      WRITE ( kunit ,'(a)')             'polarizabilities on atoms'
      lseparator(kunit) 
      do it = 1 , ntype
        if ( lpolar( it ) .eq. 1 ) then
          WRITE ( kunit ,'(a,a2,a,f12.4)')'polarizability on type ', atypei(it),' : ' 
          WRITE ( kunit ,'(3f12.4)')      ( pol ( it , 1 , j ) , j = 1 , 3 ) 
          WRITE ( kunit ,'(3f12.4)')      ( pol ( it , 2 , j ) , j = 1 , 3 ) 
          WRITE ( kunit ,'(3f12.4)')      ( pol ( it , 3 , j ) , j = 1 , 3 ) 
          blankline(kunit)
        else
          WRITE ( kunit ,'(a,a2)')      'no polarizability on type ', atypei(it)
        endif
        lseparator(kunit) 
        blankline(kunit)
      enddo
    endif
    ! =================================
    !      LONG RANGE INTERACTIONS 
    ! =================================
    if ( lcoulomb )    then 
      lseparator(kunit) 
      WRITE ( kunit ,'(a)')             'coulombic interaction : '
      lseparator(kunit) 
      blankline(kunit)
      WRITE ( kunit ,'(a)')             '        qi qj   '
      WRITE ( kunit ,'(a)')             ' Vij = -------  '
      WRITE ( kunit ,'(a)')             '         rij    '          
      blankline(kunit)
      ! =================================
      !         Direct Summation
      ! =================================
      if ( longrange .eq. 'direct' )  then
        WRITE ( kunit ,'(a)')           'direct summation'
        WRITE ( kunit ,'(a)')           'cubic cutoff in real space'
        WRITE ( kunit ,'(a,i10)')       '-ncelldirect ... ncelldirect     = ',ncelldirect
        WRITE ( kunit ,'(a,i10)')       'total number of cells            = ',( 2 * ncelldirect + 1 ) ** 3
        WRITE ( kunit ,'(a,f10.5)')     'radial cutoff                    = ',cutlongrange
      endif     
      ! =================================
      !         Ewald Summation
      ! =================================
      if ( longrange .eq. 'ewald' )  then
        if ( lautoES ) then
          WRITE ( kunit ,'(a)')         'ewald summation parameters ( automatic see field.f90 tmp construction )'
          WRITE ( kunit ,'(a,f10.5)')   'alpha                            = ',alphaES
          WRITE ( kunit ,'(a,f10.5)')   'cut-off (short range)            = ',cutshortrange
          WRITE ( kunit ,'(a,3i10)')    'kmax                             = ',(kES(i),i=1,3)
          blankline(kunit)
          WRITE ( kunit ,'(a,e12.5)')   'relative error (user defined)    : ',epsw
        else
          alpha2 = alphaES * alphaES
          do i=1,3
            kmax2 = pi * kES(i) / simu_cell%ANORM(i) / alphaES
            kmax2 = kmax2 * kmax2    
            ereci(i)  = EXP ( - kmax2 ) / kmax2 * ( SQRT ( DBLE(kES(i)) ) / alphaES / simu_cell%ANORM(i) / simu_cell%ANORM(i) )
            ereci2(i) = qtot2 * alphaES / pisq / ( SQRT ( DBLE(kES(i)) * DBLE(kES(i)) * DBLE(kES(i)) ) ) &
            * EXP ( - ( pi * DBLE(kES(i)) / alphaES / simu_cell%ANORM(i) ) ** 2 )
          enddo
          rcut2  = cutshortrange * cutshortrange
          ereal  = EXP ( - alpha2 * rcut2 ) / alpha2  / rcut2 * SQRT ( cutshortrange / 2.0_dp / simu_cell%omega )
          WRITE ( kunit ,'(a)')         'ewald summation parameters ( from input file )'
          WRITE ( kunit ,'(a,f10.5)')   'alpha                            = ',alphaES
          WRITE ( kunit ,'(a,f10.5)')   'cut-off (short range)            = ',cutshortrange
          WRITE ( kunit ,'(a,3i10)')    'kmax                             = ',(kES(i),i=1,3)
          blankline(kunit)
          WRITE ( kunit ,'(a,e12.5)')   'relative error in real space       with alphaES from input : ',ereal
          WRITE ( kunit ,'(a,3e12.5)')  'relative error in reciprocal space with alphaES from input : ',(ereci(i),i=1,3)
          WRITE ( kunit ,'(a,3e12.5)')  'relative error in reciprocal space with alphaES from input : ',(ereci2(i),i=1,3)
          blankline(kunit)
          blankline(kunit)
        endif
      endif
    endif
    ! =================================
    !     SHORT RANGE INTERACTIONS 
    ! =================================
    !       LENNARD-JONES
    ! =================================
    if ( lbmlj )       then     
      lseparator(kunit) 
      WRITE ( kunit ,'(a)')             'lennard-jones           '
      lseparator(kunit) 
      blankline(kunit)
      WRITE ( kunit ,'(a)')             '       eps    /    / sigma* \ q       / sigma*  \ p \'
      WRITE ( kunit ,'(a)')             ' V = ------- |  p | ------- |    - q | -------- |   |'   
      WRITE ( kunit ,'(a)')             '      q - p   \    \   r    /         \    r    /   /'
      blankline(kunit)
      WRITE ( kunit ,'(a,f10.5)')       'cutoff      = ',cutshortrange
      WRITE ( kunit ,'(a,a)')           'truncation  = ',ctrunc
      if ( .not.lreduced ) &
      WRITE ( kunit ,'(a,2f20.9)')      'long range correction : ',utail
      blankline(kunit)
      blankline(kunit)
      do it1 = 1 , ntype
        do it2 = it1 , ntype 
          lseparator(kunit) 
          WRITE ( kunit ,'(a,a,a,a)')   atypei(it1),'-',atypei(it2),' interactions:'    
          lseparator(kunit) 
          if ( it1 .eq. it2 ) then
            WRITE ( kunit ,100)         'sigma                                = ',sigmalj ( it1 , it2 )
            WRITE ( kunit ,100)         'eps                                  = ',epslj   ( it1 , it2 )
            WRITE ( kunit ,100)         'q                                    = ',qlj     ( it1 , it2 )
            WRITE ( kunit ,100)         'p                                    = ',plj     ( it1 , it2 )
          else
            WRITE ( kunit ,110)         'sigma                                = ',sigmalj ( it1 , it2 ) , '( ',sigmalj ( it2 , it1 ), ' )' 
            WRITE ( kunit ,110)         'eps                                  = ',epslj   ( it1 , it2 ) , '( ',epslj   ( it2 , it1 ), ' )'
            WRITE ( kunit ,110)         'q                                    = ',qlj     ( it1 , it2 ) , '( ',qlj     ( it2 , it1 ), ' )'
            WRITE ( kunit ,110)         'p                                    = ',plj     ( it1 , it2 ) , '( ',plj     ( it2 , it1 ), ' )'
          endif
          if ( trunc .eq. 1 ) then
            WRITE ( kunit ,'(a,f10.5)')'shift correction (linear)            : ',uc(it1,it2)
          else if ( trunc .eq. 2 ) then
            WRITE ( kunit ,'(a,2f10.5)')'shift correction (quadratic)         : ',uc(it1,it2),uc2(it1,it2)
          endif
        enddo
      enddo
    endif
    ! =================================
    !      BUCKINGHAM - MORSE 
    ! =================================
    if ( lmorse )       then
      blankline(kunit)
      lseparator(kunit) 
      WRITE ( kunit ,'(a)')             'morse  potential'
      lseparator(kunit) 
      blankline(kunit)
      do it1 = 1 , ntype
        do it2 = it1 , ntype
          lseparator(kunit) 
          WRITE ( kunit ,'(a,a,a,a)')   atypei(it1),'-',atypei(it2),' interactions:'
          lseparator(kunit) 
          if ( it1 .eq. it2 ) then
            WRITE ( kunit ,100)         'sigma                                = ',sigmamor ( it1 , it2 )
            WRITE ( kunit ,100)         'eps                                  = ',epsmor   ( it1 , it2 )
            WRITE ( kunit ,100)         'rho                                  = ',rhomor   ( it1 , it2 )
          else
            WRITE ( kunit ,110)         'sigma                                = ',sigmamor ( it1 , it2 ) , '( ',sigmamor ( it2 , it1 ), ' )'
            WRITE ( kunit ,110)         'eps                                  = ',epsmor   ( it1 , it2 ) , '( ',epsmor   ( it2 , it1 ), ' )'
            WRITE ( kunit ,110)         'rho                                  = ',rhomor   ( it1 , it2 ) , '( ',rhomor   ( it2 , it1 ), ' )'
          endif
        enddo
      enddo

    endif
    blankline(kunit)
    blankline(kunit)
    separator(kunit)    
    blankline(kunit)
  endif

  return
100 FORMAT(a,f10.5) 
110 FORMAT(a,f10.5,a,f10.5,a) 

END SUBROUTINE field_print_info


! *********************** SUBROUTINE gen_ewald_param ***************************
!> \brief
!! automatic determination of ewald parameter from epsw 
!> \note
!! there is several methods ( we follow dl_poly )
! ******************************************************************************
SUBROUTINE ewald_param
   
  USE constants,                ONLY :  pi , pisq
  USE config,                   ONLY :  simu_cell 
  USE control,                  ONLY :  cutshortrange
  USE io,                  ONLY :  stdout

  implicit none

  !local
  real(kind=dp) :: aaa , aaa2 , rcut 
  real(kind=dp) :: alpha , eps , tol , tol1
  integer       :: nc ( 3 )
  integer       :: kmax_1 , kmax_2 , kmax_3

  rcut =  0.5_dp * MIN(simu_cell%WA,simu_cell%WB,simu_cell%WC) - 0.1_dp
  CALL estimate_alpha( aaa , epsw ,rcut )
  CALL accur_ES_frenkel_smit( epsw , aaa2 , rcut , nc )
  ! dl_poly like
  if ( min (cutshortrange,rcut) .ne. cutshortrange ) &
  WRITE ( stdout , '(a)') 'WARNING : cutshortrange will be changed according to simu_cell%W*'
  cutshortrange = min (cutshortrange,rcut)
  eps=min(abs(epsw),0.5_dp)
  tol=sqrt(abs(log(eps*cutshortrange)))
  alpha=sqrt(abs(log(eps*cutshortrange*tol)))/rcut
  tol1=sqrt(-log(eps*cutshortrange*(2.0_dp*tol*alpha)**2))
  kmax_1=nint(0.25_dp+simu_cell%ANORM(1)*alpha*tol1/pi)
  kmax_2=nint(0.25_dp+simu_cell%ANORM(2)*alpha*tol1/pi)
  kmax_3=nint(0.25_dp+simu_cell%ANORM(3)*alpha*tol1/pi)
  if ( lautoES ) then
    ! dl_poly like
    io_node write ( stdout , '(a)' ) 'auto ES dl_poly like'  
    alphaES = alpha
    kES(1)=kmax_1
    kES(2)=kmax_2
    kES(3)=kmax_3
    ! frenkel_smit like
    !write ( * , * ) 'auto ES frenkel_smit like'  
    !alphaES = aaa2
    !kES=nc
    ! qe like
    !alphaES = aaa
    !kES=nc
  endif

  return

END SUBROUTINE ewald_param

! *********************** SUBROUTINE define_units ******************************
!> \brief
!! define conversion factor from wanted units
!> \warning
!! the units are not completely consistent in MDFF
! ******************************************************************************
SUBROUTINE define_units

  USE constants,                ONLY :  dp , boltz
  USE thermodynamic,            ONLY :  engunit

  implicit none

    engunit = 1.0_dp
!   eV => internal 
   if ( units.eq.'eV')   engunit = 9648.530821_dp
!  kcal ==> internal
   if ( units.eq.'kcal') engunit = 418.4_dp
!  kJ ==> internal
   if ( units.eq.'kJ')   engunit = 100.0_dp
!  K ==> internal
   if ( units.eq.'kJ')   engunit = boltz

   epslj  = epslj  * engunit
   epsmor = epsmor * engunit

#ifdef debug_units
   io_node WRITE ( stderr , * ) 'engunit',engunit
#endif

  return

END SUBROUTINE

! *********************** SUBROUTINE initialize_param_bmlj_morse ***************
!> \brief
!! initialisation of principal parameters for lennard-jones and morse potentials
! ******************************************************************************
SUBROUTINE initialize_param_bmlj_morse 

  USE constants,                ONLY :  tpi
  USE config,                   ONLY :  ntypemax , natm , natmi , rho , atype , itype  , ntype , simu_cell
  USE control,                  ONLY :  skindiff , cutshortrange , lreduced, calc
  USE io,                  ONLY :  ionode, stdout

  implicit none

  ! local
  integer :: it,jt
  real(kind=dp) :: one13, one16, two16, rskinmax
  real(kind=dp) :: rcut3 ( ntypemax , ntypemax )
  real(kind=dp) :: rskin ( ntypemax , ntypemax ) 
  real(kind=dp) :: rskinsq ( ntypemax , ntypemax ) 
  real(kind=dp) :: ut ( ntypemax , ntypemax ) 

  real(kind=dp) :: rcut ( ntype , ntype ) 
  real(kind=dp) :: ppqq ( ntype , ntype )
  real(kind=dp) :: pp ( ntype , ntype )
  real(kind=dp) :: qq ( ntype , ntype )
  real(kind=dp) :: pp3 ( ntype , ntype ) 
  real(kind=dp) :: qq3 ( ntype , ntype ) 
  real(kind=dp) :: sr2 ( ntype , ntype ) 
  real(kind=dp) :: sr ( ntype , ntype ) 
  real(kind=dp) :: srp ( ntype , ntype ) 
  real(kind=dp) :: srq ( ntype , ntype ) 

  rskinmax = 0.0_dp
  utail    = 0.0_dp
  rcut3    = 0.0_dp
  rcut     = 0.0_dp
  rskin    = 0.0_dp
  rskinsq  = 0.0_dp
  ut       = 0.0_dp
  ppqq     = 0.0_dp
  pp       = 0.0_dp
  qq       = 0.0_dp
  sr2      = 0.0_dp
  sr       = 0.0_dp
  srp      = 0.0_dp
  srq      = 0.0_dp
  
  do it = 1 , ntype
    do jt = 1 , ntype
      pp ( it , jt ) = plj ( it , jt ) 
      qq ( it , jt ) = qlj ( it , jt ) 
    enddo
  enddo
  pp3 = pp - 3.0_dp
  qq3 = qq - 3.0_dp

  one13 = (1.0_dp / 3.0_dp)
  one16 = (1.0_dp / 6.0_dp)
  two16 = 2.0_dp **  one16

  ! ==================================================================================================
  ! TAIL ( checked september 2011) :
  ! be careful two minus are vanishing
  !
  !     2 pi rc^3 espilon NA NB    /     p         / sigma* \ q           q         / sigma*  \ p  \
  !    -------------------------- |  ----------   | ------- |    --   ----------   | -------- |    |
  !         ( q - p )   V          \ ( q - 3 )     \   rc   /          ( p - 3 )    \  rc     /    /
  !
  !
  ! rskinmax ??????
  ! ==================================================================================================
  
  
  ! ==================================================================================================
  !  simple truncation  
  !
  !  ctrunc = 'linear'
  !  trunc = 1 
  !
  !
  !           eps    /    / sigma* \ q         / sigma* \ p  \
  !  V  =   ------- |  p | ------- |    -  q  | --------|    |   -  c1   with  sigma* = 2^(1/6)*sigma
  !          q - p   \    \   r    /           \    r   /    /
  !
  ! 
  !          eps      /     /  sigma* \ q        /  sigma* \ p  \
  !  c1 = ---------  |   p |  -------- |    - q |  -------- |   |      with rc = cutshortrange 
  !         q - p     \     \   rc    /          \   rc    /    / 
  ! 
  ! ==================================================================================================
  

  ! ==================================================================================================
  !  truncation presented in J. Chem. Phys. 135 (2011) , Sengupta, Vasconcelos, Affouard, Sastry
  ! 
  !  ctrunc = 'quadratic'
  !  trunc  = 2    
  ! 
  !
  !           eps    /    / sigma* \ q         / sigma* \ p  \
  !  V  =   ------- |  p | ------- |    -  q  | --------|    |   +  c1 r^2  -  c2    with  sigma* = 2^(1/6)*sigma
  !          q - p   \    \   r   /            \    r   /    /
  !
  !
  ! Sage (January 2013) 
  !     
  !              eps p  q           /  / sigma* \ q       /  sigma* \ p \
  !    c1 =  --------------------- |  | --------|    -   | --------- |  |
  !          2  ( q - p ) * rc^2    \  \   rc   /         \   rc    /   /
  ! 
  !    and
  !     
  !           epsilon     /           / sigma* \ q               / sigma* \ p  \
  !    c2 =  ----------- |  (2p+qp)  | --------|     - (2q+pq)  | --------|   |
  !          2 ( q - p )  \           \   rc   /                 \   rc   /   /       
  !
  ! ==================================================================================================
  do it = 1 , ntype 
    do jt = 1 , ntype

      ! intermediate terms 
          rcut   ( it , jt ) = cutshortrange                                      ! rc
          rcutsq ( it , jt ) = rcut ( it , jt ) * rcut   ( it , jt )       ! rc^2 
          rcut3  ( it , jt ) = rcut ( it , jt ) * rcutsq ( it , jt )       ! rc^3
          ppqq   ( it , jt ) = pp   ( it , jt ) * qq     ( it , jt )       ! p x q
          rskin  ( it , jt ) = rcut ( it , jt ) + skindiff

         if ( rskin ( it , jt ) .gt. rskinmax ) rskinmax = rskin ( it , jt )
           rskinsq ( it , jt ) = rskin ( it , jt ) * rskin ( it , jt )
           ! eps / q - p   
           epsp ( it , jt )= epslj( it , jt )/( qq ( it , jt )-pp ( it , jt ) )
           ! sigma*^2 
           sigsq ( it , jt )= two16 * two16 * sigmalj ( it , jt ) * sigmalj ( it , jt )
           ! sigma*^2 / rc ^2
           sr2 ( it , jt ) = sigsq ( it , jt )/ rcutsq ( it , jt )
           ! sigma* / rc 
           sr( it ,jt )    = SQRT ( sr2( it , jt )  )
           ! (sigma* / rc ) ^ p
           srp( it , jt )  = sr ( it , jt ) ** pp ( it , jt )
           ! (sigma* / rc ) ^ q 
           srq( it , jt )  = sr ( it , jt ) ** qq ( it , jt )
           ! trunc = 1
           uc ( it , jt )  = epsp ( it , jt ) * ( pp ( it , jt ) * srq( it , jt ) - qq ( it , jt ) * srp( it , jt ) )
           ! trunc = 2
           ! c1
           uc1 ( it , jt ) = epsp ( it , jt ) *  ppqq ( it , jt ) / ( 2.0_dp * rcutsq ( it , jt ) ) * ( srq( it , jt ) - srp( it , jt ) ) 
           ! c2
           uc2 ( it , jt ) = 0.5_dp * epsp( it , jt )  * (  &
                         ( 2.0_dp * pp ( it , jt ) + ppqq ( it , jt )  ) * srq( it , jt ) - &
                         ( 2.0_dp * qq ( it , jt ) + ppqq ( it , jt )  ) * srp( it , jt ) )
           ! for the virial
           fc ( it , jt ) =  ppqq ( it , jt ) * epsp ( it , jt ) /  sigsq ( it , jt ) 
           ! morse
           fm ( it , jt ) = - 2.0_dp * epsmor ( it , jt ) * rhomor ( it , jt ) * EXP ( rhomor ( it , jt ) * sigmamor ( it , jt ) )  
           rs ( it , jt ) = EXP ( rhomor ( it , jt ) * sigmamor( it , jt ) ) 
           ! tail energy
           ut ( it , jt ) = epsp ( it , jt ) * ( pp ( it , jt ) * srq ( it , jt ) / qq3 ( it , jt ) - &
                                                 qq ( it , jt ) * srp ( it , jt ) / pp3 ( it , jt )  )       
           ut ( it , jt ) = ut ( it , jt ) * rcut3 ( it , jt ) * tpi 
           if ( ( natmi ( it ) .ne. 0 ) .and. ( natmi ( jt ) .ne. 0 ) ) &
           utail = utail + ut ( it , jt ) * natmi ( it ) * natmi ( jt ) / simu_cell%omega

           io_node write(stdout,*) 'long range correction (init)',utail

!#ifdef debug
!  WRITE ( stdout , '(2i6,7e16.6)' ) it , jt , uc ( it , jt )  , epsp ( it , jt ) , pp ( it , jt ) , qq ( it , jt ) , srq( it , jt ) , srp( it , jt ) ,rcutsq ( it , jt ) 
!#endif
    enddo
  enddo
 
#ifdef debug_quadratic
  do it = 1 , ntype 
    WRITE ( stdout , '(a,<ntype>e16.6)' ) 'uc1 ',(uc1(it,jt),jt=1,ntype)
    WRITE ( stdout , '(a,<ntype>e16.6)' ) 'uc2 ',(uc2(it,jt),jt=1,ntype)
  enddo
#endif

  return

END SUBROUTINE initialize_param_bmlj_morse


! *********************** SUBROUTINE engforce_driver ***************************
!
!> \brief
!! this subroutine is main driver to select the different potential and forces
!! subroutines
!
!> \todo
!! make it more clean
!
! ******************************************************************************
SUBROUTINE engforce_driver 

  USE config,                   ONLY :  natm , dipia , fx,rx,vx
  USE control,                  ONLY :  lpbc , lminimg , lbmlj , lcoulomb , lmorse , lharm , lshiftpot , longrange
  USE io,                  ONLY :  ionode , stdout 

  implicit none

  ! local 
  logical :: lj_and_coul
  real(kind=dp) :: eftmp( natm , 3 ) , efgtmp ( natm , 3 , 3 ) , u_coultmp , vir_coultmp , phi_coultmp ( natm ) 

  ! harmonic oscillator ( test purpose )
  if ( lharm ) then
    CALL engforce_harm
  endif


  ! ====================================
  !  check if LJ and COULOMBIC 
  !  interactions are needed together
  ! ====================================
  lj_and_coul = .false.
  if ( lbmlj .and. lcoulomb ) then
    lj_and_coul = .true.
  endif

  if ( lmorse ) then
    lbmlj=.false.
    lcoulomb =.false.
    CALL engforce_morse_pbc         
  endif

  !  io_node WRITE ( stdout ,'(a)') 'debug ' 
  ! ==============================
  !  LJ + COULOMBIC INTERACTIONS 
  ! ==============================
  if ( lj_and_coul ) then
#ifdef debug  
    WRITE ( stdout ,'(a)') 'debug : lj + coulomb' 
#endif  
    ! =======
    !   PBC
    ! =======
    if ( lpbc ) then
      if ( lshiftpot .and. lminimg ) CALL engforce_bmlj_pbc         
      if ( .not.lshiftpot )          CALL engforce_bmlj_pbc_noshift 
#ifdef debug  
      print*,'before coul'
      CALL print_config_sample(0,0)
#endif
      if ( longrange .eq. 'ewald'  ) CALL multipole_ES ( eftmp , efgtmp , dipia , u_coultmp , vir_coultmp , phi_coultmp )
      if ( longrange .eq. 'direct' ) CALL multipole_DS ( eftmp , efgtmp , dipia , u_coultmp , vir_coultmp , phi_coultmp )
#ifdef debug  
      print*,'before after coul'
      CALL print_config_sample(0,0)
#endif
      ! =======
      !  NO PBC
      ! =======
    else
      WRITE ( stdout ,'(a)') 'ERROR sick job: coulomb + nopbc'
      if ( lshiftpot )               CALL engforce_bmlj_nopbc       
!      if ( .not.lshiftpot )          CALL engforce_bmlj_pbc_noshift 
    endif
    ! ================================
    !      ONLY LENNARD-JONES
    ! ================================
    elseif ( lbmlj ) then
#ifdef debug  
      WRITE ( stdout , '(a)' ) 'debug : lj only' 
#endif  
     ! =======
     !   PBC
     ! =======
     if ( lpbc ) then
#ifdef debug  
       WRITE ( stdout , '(a)' ) 'debug : pbc ' 
#endif  
!   if ( ionode ) then
!    print*,'before lj'
!    write(*,'(a,f60.40)') 'ri before prop',rx(1)
!    write(*,'(a,f60.40)') 'rj before prop',rx(167)
!    write(*,'(a,f60.40)') 'vi before prop',vx(1)
!    write(*,'(a,f60.40)') 'vj before prop',vx(167)
!    write(*,'(a,f60.40)') 'fi before prop',fx(1)
!    write(*,'(a,f60.40)') 'fj before prop',fx(167)
!  endif
       if ( lshiftpot       )     CALL engforce_bmlj_pbc         
       if ( .not. lshiftpot )     CALL engforce_bmlj_pbc_noshift 
!   if ( ionode ) then
!    print*,'after lj'
!    write(*,'(a,f60.40)') 'ri before prop',rx(1)
!    write(*,'(a,f60.40)') 'rj before prop',rx(167)
!    write(*,'(a,f60.40)') 'vi before prop',vx(1)
!    write(*,'(a,f60.40)') 'vj before prop',vx(167)
!    write(*,'(a,f60.40)') 'fi before prop',fx(1)
!    write(*,'(a,f60.40)') 'fj before prop',fx(167)
!  endif

       ! =======
       !  NO PBC
       ! =======
     else
#ifdef debug  
       WRITE ( stdout , '(a)' ) 'debug : no pbc ' 
#endif  
       if ( lshiftpot )          CALL engforce_bmlj_nopbc          
!      if ( .not. lshiftpot )   CALL engforce_bmlj_nopbc_noshift 
       if ( .not. lshiftpot ) then
         WRITE ( stdout , '(a)' ) 'not yet : nopbc + no shift pot'
         STOP
       endif
     endif

     ! ================================
     !      ONLY COULOMB 
     ! ================================
  elseif ( lcoulomb ) then
#ifdef debug  
    WRITE ( stdout ,'(a)') 'debug : coulomb only ' 
#endif  
    ! =======
    !   PBC
    ! =======
    if ( lpbc ) then
      if ( longrange .eq. 'ewald'  ) CALL multipole_ES    ( eftmp , efgtmp , dipia , u_coultmp , vir_coultmp , phi_coultmp )
      if ( longrange .eq. 'direct' ) CALL multipole_DS    ( eftmp , efgtmp , dipia , u_coultmp , vir_coultmp , phi_coultmp )
      ! =======
      !  NO PBC
      ! =======
    else
      WRITE ( stdout ,'(a)') 'ERROR sick job: coulomb + nopbc'
      STOP 
    endif

  endif

  return

END SUBROUTINE engforce_driver

! *********************** SUBROUTINE engforce_bmlj_pbc *************************
!
!> \brief
!! total potential energy forces for each atoms for a bmlj potential with 
!! periodic boundaries conditions, with or without vnlist.
!
!> \author
!! F.Affouard / FMV
!
!> \note
!! adapted from F. Affouard code. Parallelized in december 2008
!
! ******************************************************************************
SUBROUTINE engforce_bmlj_pbc 

  USE config,                   ONLY :  natm , rx , ry , rz , fx , fy , fz, vx, vy , vz , tau_nonb ,  &
                                        atype , itype , list , point , ntype , simu_cell , atom_dec
  USE control,                  ONLY :  lvnlist 
  USE thermodynamic,            ONLY :  u_lj , vir_lj , write_thermo
  USE io,                  ONLY :  stdout
  USE time,                     ONLY :  forcetimetot 
  USE cell,                     ONLY :  kardir , dirkar

  implicit none

  ! local
  integer          :: ia , ja , it , jt , j1 , je , jb , ierr 
  integer          :: p1 , p2
  real(kind=dp) :: rxi , ryi , rzi
  real(kind=dp) :: rxij , ryij , rzij 
  real(kind=dp) :: sxij , syij , szij 
  real(kind=dp) :: sr2 , rijsq , srp , srq
  real(kind=dp) :: wij , fxij , fyij , fzij
  real(kind=dp) :: ptwo ( ntype , ntype )
  real(kind=dp) :: qtwo ( ntype , ntype )
  real(kind=dp) :: ttt1 , ttt2 
  real(kind=dp) :: u , vir 

#ifdef debug_bmlj
  WRITE ( stdout , '(a,2i6)' ) 'debug : atom decomposition ',atom_dec%istart,atom_dec%iend
  do ia=1,natm
    WRITE ( stdout , '(a,i6,a,a)' )  'debug : atype ',ia,'',atype(ia)
    WRITE ( stdout , '(a,i6,a,i4)' ) 'debug : itype ',ia,'',itype(ia)
  enddo
#endif

  ttt1 = MPI_WTIME(ierr) ! timing info
  cccc = cccc + 1

  u   = 0.0_dp
  vir = 0.0_dp
  fx  = 0.0_dp
  fy  = 0.0_dp
  fz  = 0.0_dp
  tau_nonb = 0.0d0
 
  do jt = 1 , ntype
    do it = 1 , ntype
      ptwo ( it , jt ) = plj ( it , jt ) * 0.5_dp
      qtwo ( it , jt ) = qlj ( it , jt ) * 0.5_dp
    enddo
  enddo

  if ( lvnlist ) CALL vnlistcheck 

  ! ======================================
  !         cartesian to direct 
  ! ======================================
  CALL kardir ( natm , rx , ry , rz , simu_cell%B )

  do ia = atom_dec%istart , atom_dec%iend
    rxi = rx ( ia )
    ryi = ry ( ia )
    rzi = rz ( ia )
    ! =====================================
    !  verlet list : ja index in point arrays
    ! =====================================
    if ( lvnlist ) then
      jb = point( ia )
      je = point( ia + 1 ) - 1
    else
    ! ====================================
    !         else all ja   
    ! ====================================
      jb = 1 
      je = natm
    endif
    do j1 = jb, je
      if ( lvnlist ) then
        ja = list ( j1 )
      else 
        ja = j1
      endif
      if ( ( lvnlist .and. ja .ne. ia ) .or. ( .not. lvnlist .and. ja .gt. ia )  ) then
!      if ( ( lvnlist .and. ja .ne. ia ) .or. ( .not. lvnlist .and.  ( ( ia .gt. ja .and. ( MOD ( ia + ja , 2 ) .eq. 0 ) ) .or. &
!                                                                      ( ia .lt. ja .and. ( MOD ( ia + ja , 2 ) .ne. 0 ) ) ) ) ) then
        rxij = rxi - rx ( ja )
        ryij = ryi - ry ( ja )
        rzij = rzi - rz ( ja )
        sxij = rxij - nint ( rxij )
        syij = ryij - nint ( ryij )
        szij = rzij - nint ( rzij )
        rxij = sxij * simu_cell%A(1,1) + syij * simu_cell%A(1,2) + szij * simu_cell%A(1,3)
        ryij = sxij * simu_cell%A(2,1) + syij * simu_cell%A(2,2) + szij * simu_cell%A(2,3)
        rzij = sxij * simu_cell%A(3,1) + syij * simu_cell%A(3,2) + szij * simu_cell%A(3,3)
        rijsq = rxij * rxij + ryij * ryij + rzij * rzij
        p1 = itype ( ia )
        p2 = itype ( ja )
        if ( rijsq .lt. rcutsq(p1,p2) ) then
          sr2 = sigsq(p1,p2) / rijsq
          srp = sr2 ** (ptwo(p1,p2))
          srq = sr2 ** (qtwo(p1,p2))
          ! potential energy ( truncated )  
          if ( trunc .eq. 0 ) then
            u = u  + epsp(p1,p2) * ( plj(p1,p2) * srq -qlj(p1,p2) * srp )
          endif
          if ( trunc .eq. 1 ) then
            u =  u + epsp(p1,p2) * ( plj(p1,p2) * srq -qlj(p1,p2) * srp ) - uc(p1,p2)
          endif
          if ( trunc .eq. 2 ) then
            u =  u + epsp(p1,p2) * ( plj(p1,p2) * srq -qlj(p1,p2) * srp ) + uc1(p1,p2) * rijsq - uc2(p1,p2)
          endif
          wij = fc(p1,p2) * (srq-srp) * sr2
          fxij = wij * rxij
          fyij = wij * ryij
          fzij = wij * rzij
          !virial 
          vir = vir + wij * rijsq
          ! forces 
          fx ( ia ) = fx ( ia ) + fxij
          fy ( ia ) = fy ( ia ) + fyij
          fz ( ia ) = fz ( ia ) + fzij
          fx ( ja ) = fx ( ja ) - fxij
          fy ( ja ) = fy ( ja ) - fyij
          fz ( ja ) = fz ( ja ) - fzij
          ! stress tensor
          tau_nonb(1,1) = tau_nonb(1,1) + rxij * fxij
          tau_nonb(1,2) = tau_nonb(1,2) + rxij * fyij
          tau_nonb(1,3) = tau_nonb(1,3) + rxij * fzij
          tau_nonb(2,1) = tau_nonb(2,1) + ryij * fxij
          tau_nonb(2,2) = tau_nonb(2,2) + ryij * fyij
          tau_nonb(2,3) = tau_nonb(2,3) + ryij * fzij
          tau_nonb(3,1) = tau_nonb(3,1) + rzij * fxij
          tau_nonb(3,2) = tau_nonb(3,2) + rzij * fyij
          tau_nonb(3,3) = tau_nonb(3,3) + rzij * fzij
        endif
      endif
    enddo
  enddo
  tau_nonb = tau_nonb / simu_cell%omega 
  vir = vir/3.0_dp

  ttt2 = MPI_WTIME(ierr) ! timing info
  forcetimetot = forcetimetot + ( ttt2 - ttt1 )

  CALL MPI_ALL_REDUCE_DOUBLE_SCALAR ( u   ) 
  CALL MPI_ALL_REDUCE_DOUBLE_SCALAR ( vir ) 
  
  CALL MPI_ALL_REDUCE_DOUBLE ( fx , natm ) 
  CALL MPI_ALL_REDUCE_DOUBLE ( fy , natm ) 
  CALL MPI_ALL_REDUCE_DOUBLE ( fz , natm ) 

  CALL MPI_ALL_REDUCE_DOUBLE ( tau_nonb( 1, : ) , 3  )
  CALL MPI_ALL_REDUCE_DOUBLE ( tau_nonb( 2, : ) , 3  )
  CALL MPI_ALL_REDUCE_DOUBLE ( tau_nonb( 3, : ) , 3  )

  u_lj = u
  vir_lj = vir

  ! ======================================
  !         direct to cartesian
  ! ======================================
  CALL dirkar ( natm , rx , ry , rz , simu_cell%A )

  return

END SUBROUTINE engforce_bmlj_pbc

! *********************** SUBROUTINE engforce_bmlj_nopbc ***********************
!
!> \brief
!! total potential energy forces for each atoms for a bmlj potential with 
!! *NO* periodic boundaries conditions, with or without vnlist (lvnlist=.TRUE.OR.FALSE.)
!
! ******************************************************************************
SUBROUTINE engforce_bmlj_nopbc 

  USE config,                   ONLY :  natm , rx , ry , rz , fx , fy , fz , itype , list , point , ntype , tau_nonb , simu_cell , atom_dec
  USE control,                  ONLY :  lvnlist
  USE thermodynamic,            ONLY :  u_lj , vir_lj

  implicit none

  ! local
  integer :: ia , ja , it, jt, j1, je, jb !, ierr
  integer :: p1, p2
  real(kind=dp) :: rxi , ryi , rzi
  real(kind=dp) :: rxij , ryij , rzij , sr2 , rijsq , srp , srq
  real(kind=dp) :: wij , fxij , fyij , fzij
  real(kind=dp) :: ptwo ( ntype , ntype )
  real(kind=dp) :: qtwo ( ntype , ntype )
  real(kind=dp) :: u, vir

  u = 0.0_dp
  vir = 0.0_dp
  fx = 0.0_dp
  fy = 0.0_dp
  fz = 0.0_dp
  tau_nonb = 0.0_dp

  do jt = 1 , ntype
    do it = 1 , ntype
      ptwo ( it , jt ) = plj ( it , jt ) * 0.5_dp
      qtwo ( it , jt ) = qlj ( it , jt ) * 0.5_dp
    enddo
  enddo

  if ( lvnlist ) CALL vnlistcheck 
  do ia = atom_dec%istart, atom_dec%iend
    rxi = rx ( ia )
    ryi = ry ( ia )
    rzi = rz ( ia )
    if ( lvnlist ) then
      jb = point ( ia )
      je = point ( ia + 1 ) - 1
    else
      jb = ia 
      je = atom_dec%iend
    endif
    do j1 = jb, je
      if ( lvnlist ) then
        ja = list(j1)
      else
        ja = j1
      endif
      if ( ( lvnlist .and. ja .ne. ia ) .or. ( .not. lvnlist .and. ja .gt. ia )  ) then
        rxij = rxi - rx ( ja )
        ryij = ryi - ry ( ja )
        rzij = rzi - rz ( ja )
        rijsq = rxij * rxij + ryij * ryij + rzij * rzij
        p1 = itype ( ia )
        p2 = itype ( ja )
        if ( rijsq .lt. rcutsq(p1,p2) ) then
          sr2 = sigsq(p1,p2) / rijsq
          srp = sr2 ** (ptwo(p1,p2))
          srq = sr2 ** (qtwo(p1,p2))
          if ( trunc .eq. 1 ) then
            u =  u + epsp(p1,p2) * ( plj(p1,p2) * srq - qlj(p1,p2) * srp ) - uc(p1,p2)
          endif
          if ( trunc .eq. 2 ) then
            u =  u + epsp(p1,p2) * ( plj(p1,p2) * srq -qlj(p1,p2) * srp ) + uc1(p1,p2) * rijsq - uc2(p1,p2)
          endif
          wij = fc(p1,p2) * (srq-srp) * sr2
          vir = vir + wij * rijsq
          fxij = wij * rxij
          fyij = wij * ryij
          fzij = wij * rzij
          fx ( ia ) = fx ( ia ) + fxij
          fy ( ia ) = fy ( ia ) + fyij
          fz ( ia ) = fz ( ia ) + fzij
          fx ( ja ) = fx ( ja ) - fxij
          fy ( ja ) = fy ( ja ) - fyij
          fz ( ja ) = fz ( ja ) - fzij
          ! stress tensor
          tau_nonb(1,1) = tau_nonb(1,1) + rxij * fxij
          tau_nonb(1,2) = tau_nonb(1,2) + rxij * fyij
          tau_nonb(1,3) = tau_nonb(1,3) + rxij * fzij
          tau_nonb(2,1) = tau_nonb(2,1) + ryij * fxij
          tau_nonb(2,2) = tau_nonb(2,2) + ryij * fyij
          tau_nonb(2,3) = tau_nonb(2,3) + ryij * fzij
          tau_nonb(3,1) = tau_nonb(3,1) + rzij * fxij
          tau_nonb(3,2) = tau_nonb(3,2) + rzij * fyij
          tau_nonb(3,3) = tau_nonb(3,3) + rzij * fzij
        endif
      endif
    enddo
  enddo
  tau_nonb = tau_nonb / simu_cell%omega 
  vir = vir / 3.0_dp

  CALL MPI_ALL_REDUCE_DOUBLE_SCALAR ( u ) 
  CALL MPI_ALL_REDUCE_DOUBLE_SCALAR ( vir ) 
  
  CALL MPI_ALL_REDUCE_DOUBLE ( fx , natm )
  CALL MPI_ALL_REDUCE_DOUBLE ( fy , natm )
  CALL MPI_ALL_REDUCE_DOUBLE ( fz , natm )
  CALL MPI_ALL_REDUCE_DOUBLE ( tau_nonb ( 1 , : )  , 3  )
  CALL MPI_ALL_REDUCE_DOUBLE ( tau_nonb ( 2 , : )  , 3  )
  CALL MPI_ALL_REDUCE_DOUBLE ( tau_nonb ( 3 , : )  , 3  )

  u_lj = u
  vir_lj = vir

  return

END SUBROUTINE engforce_bmlj_nopbc

! *********************** SUBROUTINE engforce_bmlj_pbc_noshift *****************
!
!> \brief
!! same as engforce_bmlj_pbc but with no shift in the potential 
!
!> \note
!! if ok it should be merged !
!
!> \todo
!! test it !
!
! ******************************************************************************
SUBROUTINE engforce_bmlj_pbc_noshift 

  USE config,                   ONLY :  natm , rx , ry , rz , fx , fy , fz , itype , list , point , ntype , simu_cell , atom_dec
  USE control,                  ONLY :  lvnlist 
  USE thermodynamic,            ONLY :  u_lj , vir_lj
  USE time,                     ONLY :  forcetimetot

  implicit none

  ! local
  integer       :: ia , ja , it , jt , j1 , je , jb , ierr
  integer       :: p1 , p2
  real(kind=dp) :: rxi , ryi, rzi
  real(kind=dp) :: rxij , ryij , rzij
  real(kind=dp) :: sxij , syij , szij
  real(kind=dp) :: sr2 , rijsq , srp , srq
  real(kind=dp) :: wij , fxij , fyij , fzij
  real(kind=dp) :: ptwo ( ntype , ntype )
  real(kind=dp) :: qtwo ( ntype , ntype )
  real(kind=dp) :: ttt1 , ttt2
  real(kind=dp) :: u , vir
  
  ttt1 = MPI_WTIME(ierr) ! timing info
  
  u   = 0.0_dp
  vir = 0.0_dp
  fx  = 0.0_dp
  fy  = 0.0_dp
  fz  = 0.0_dp
  

  do jt = 1 , ntype
    do it = 1 , ntype
      ptwo ( it , jt ) = plj ( it , jt ) * 0.5_dp
      qtwo ( it , jt ) = qlj ( it , jt ) * 0.5_dp
    enddo
  enddo

  if ( lvnlist ) CALL vnlistcheck 
  do ia = atom_dec%istart , atom_dec%iend
    rxi = rx ( ia )
    ryi = ry ( ia )
    rzi = rz ( ia )
    if ( lvnlist ) then
      jb = point ( ia )
      je = point ( ia + 1 ) - 1
    else
      jb = ia 
      je = atom_dec%iend 
    endif
    do j1 = jb, je
      if ( lvnlist ) then
        ja = list ( j1 )
      else
        ja = j1
      endif 
      if ( ( lvnlist .and. ja .ne. ia ) .or. ( .not. lvnlist .and. ja .gt. ia )  ) then
        rxij  = rxi - rx ( ja )
        ryij  = ryi - ry ( ja )
        rzij  = rzi - rz ( ja )
        sxij = rxij - nint ( rxij )
        syij = ryij - nint ( ryij )
        szij = rzij - nint ( rzij )
        rxij = sxij * simu_cell%A(1,1) + syij * simu_cell%A(1,2) + szij * simu_cell%A(1,3)
        ryij = sxij * simu_cell%A(2,1) + syij * simu_cell%A(2,2) + szij * simu_cell%A(2,3)
        rzij = sxij * simu_cell%A(3,1) + syij * simu_cell%A(3,2) + szij * simu_cell%A(3,3)
        rijsq = rxij * rxij + ryij * ryij + rzij * rzij
        p1    = itype ( ia )
        p2    = itype ( ja )
        if ( rijsq .lt. rcutsq(p1,p2) ) then
          sr2 = sigsq(p1,p2) / rijsq
          srp = sr2 ** (ptwo(p1,p2))
          srq = sr2 ** (qtwo(p1,p2))
          u =  u + epsp(p1,p2) * ( plj(p1,p2) * srq - qlj(p1,p2) * srp ) 
          wij = fc(p1,p2) * (srq-srp) * sr2
          vir = vir + wij * rijsq
          fxij = wij * rxij
          fyij = wij * ryij
          fzij = wij * rzij
          fx ( ia ) = fx ( ia ) + fxij
          fy ( ia ) = fy ( ia ) + fyij
          fz ( ia ) = fz ( ia ) + fzij
          fx ( ja ) = fx ( ja ) - fxij
          fy ( ja ) = fy ( ja ) - fyij
          fz ( ja ) = fz ( ja ) - fzij
        endif
      endif
    enddo
  enddo
  vir = vir/3.0_dp

  ttt2 = MPI_WTIME(ierr) ! timing info
  forcetimetot = forcetimetot + ( ttt2 - ttt1 )

  CALL MPI_ALL_REDUCE_DOUBLE_SCALAR ( u )
  CALL MPI_ALL_REDUCE_DOUBLE_SCALAR ( vir )

  CALL MPI_ALL_REDUCE_DOUBLE ( fx , natm )
  CALL MPI_ALL_REDUCE_DOUBLE ( fy , natm )
  CALL MPI_ALL_REDUCE_DOUBLE ( fz , natm )

  u_lj = u
  vir_lj = vir

  return

END SUBROUTINE engforce_bmlj_pbc_noshift

! *********************** SUBROUTINE initialize_coulomb ************************
!> \brief
!! this subroutine initialize the common quantities for charged particules.
!! real space and reciprocal space summation
!> \note
!! It is used by EFG and Coulombic subroutines 
!> \warning
!! EFG routines and + Coulombic forces routines has be used together 
! ******************************************************************************
SUBROUTINE initialize_coulomb

  USE config,   ONLY  : natm , natmi , ntype , qia , itype
  USE control,  ONLY  : longrange
  USE rspace,   ONLY  : direct_sum_init
  USE kspace,   ONLY  : kpoint_sum_init 

  implicit none

  ! local
  integer :: nk , ncmax 

  allocate ( ef_t ( natm , 3 ) )
  allocate ( efg_t ( natm , 3 , 3 ) )

  ! ============
  !  direct sum
  ! ============
  if ( longrange .eq. 'direct' ) then
   rm_coul%meshlabel='rm_coul'
   ncmax = ( 2 * ncelldirect + 1 ) ** 3
   rm_coul%ncmax=ncmax
   rm_coul%ncell=ncelldirect
   allocate( rm_coul%boxxyz( 3 , ncmax ) , rm_coul%lcell( ncmax ) , rm_coul%rr ( ncmax ) )
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
    nk = ( 2 * km_coul%kmax(1) + 1 ) * ( 2 * km_coul%kmax(2) + 1 ) * ( 2 * km_coul%kmax(3) + 1 )   
    km_coul%nk = nk
    allocate ( km_coul%kptk( nk ) , km_coul%kpt(3,nk) )
    allocate ( km_coul%strf ( nk, natm ) ) 
    CALL kpoint_sum_init ( km_coul , alphaES )
  endif

  return

END SUBROUTINE initialize_coulomb

! *********************** SUBROUTINE finalize_coulomb **************************
!> \brief
!! deallocate main quanties used during coulombic calculation
! ******************************************************************************
SUBROUTINE finalize_coulomb

  USE control,  ONLY :  longrange , lcoulomb , calc

  implicit none

  if ( .not. lcoulomb .or. calc.ne.'md') return

  deallocate ( ef_t )
  deallocate ( efg_t )
  ! ============
  !  direct sum
  ! ============
  if ( longrange .eq. 'direct' ) then
    deallocate( rm_coul%boxxyz , rm_coul%lcell , rm_coul%rr )
  endif

  ! ============
  !  ewald sum
  ! ============
  if ( longrange .eq. 'ewald')  then
    deallocate( km_coul%kptk , km_coul%kpt )
    deallocate ( km_coul%strf )
  endif


END SUBROUTINE finalize_coulomb

! *********************** SUBROUTINE induced_moment ****************************
!> \brief
!! this subroutine calculates the induced moment from the total Electric field and the
!! polarizability tensor
!! Basically used in the SCF loop to get induced_moment
!! \f$ \mu_{i,\alpha} =  p_{i,\alpha,\beta} * E_{i,\beta} \f$
!> \param[in]  Efield electric field vector define at ion position 
!> \param[out] mu_ind induced electric dipole define at ion position 
!> \param[in]  u_pol potential energy of polarizability
!> \note
!! polia is the polarizability tensor
! ******************************************************************************
SUBROUTINE induced_moment ( Efield , mu_ind , u_pol )

  USE config, ONLY : natm , itype , atypei, ntype , polia
  USE io, ONLY : stdout

  implicit none

  ! global
  real(kind=dp) , intent ( in  ) :: Efield ( natm , 3 ) 
  real(kind=dp) , intent ( out ) :: mu_ind ( natm , 3 ) 
  real(kind=dp) , intent ( out ) :: u_pol   

  ! local 
  integer :: alpha , beta
  integer :: ia , it 
  real(kind=dp) :: invpol ( 3 , 3 )   
  integer, parameter :: LWORK=1000  
  real(kind=dp) :: WORK ( LWORK ) 
  integer :: ipiv ( 3 ) 
  integer :: ierr

  mu_ind = 0.0_dp
  do ia = 1 , natm 
    do alpha = 1 , 3 
      do beta = 1 , 3  
        mu_ind ( ia , alpha ) = mu_ind ( ia , alpha ) + polia ( ia , alpha , beta ) * Efield ( ia , beta )  
      enddo
    enddo
  enddo

  u_pol = 0.0_dp 
  do ia = 1 , natm
    it = itype(ia) 
    if ( lpolar ( it ) .eq. 1 ) then 
      invpol ( : , : )  = polia ( ia , : , : )
      CALL DGETRF( 3, 3, invpol, 3, ipiv, ierr )
      if ( ierr.lt.0 ) then
          WRITE( stdout , '(a,i6)' ) 'ERROR call to DGETRF failed in induced_moment',ierr
          STOP
      endif
      CALL DGETRI( 3 , invpol , 3 ,  ipiv , WORK, LWORK, ierr )
      if ( ierr.lt.0 ) then
          WRITE( stdout , '(a,i6)' ) 'ERROR call to DGETRI failed in induced_moment',ierr
          STOP
      endif

      do alpha = 1 , 3
        do beta = 1 , 3
           u_pol = u_pol + mu_ind ( ia , alpha ) * invpol ( alpha , beta ) * mu_ind ( ia , beta ) 
        enddo
      enddo

    endif
  enddo

  u_pol = u_pol * 0.5_dp

  return

END SUBROUTINE induced_moment

! *********************** SUBROUTINE multipole_DS ******************************
!
!> \brief 
!! This subroutine calculates electric field, electric field gradient, 
!! potential energy, virial, electric potential and forces at ions in
!! multipole expansion by Direct summation
!
!> \param[in]  mu electric dipole at ions
!> \param[out] ef electric field
!> \param[out] efg electric field gradient
!> \param[out] u_coul coulombic potential energy 
!> \param[out] virial_coul coulombic virial 
!> \param[out] phi_coul electric potential at ions 
!
!> \note
!! it is not efficient for large systems
!! some quantities are not even converges at all
!
! ******************************************************************************
SUBROUTINE multipole_DS ( ef , efg , mu , u_coul , vir_coul , phi_coul )

  USE config,                   ONLY :  natm , atype , natmi , ntype , qia , rx , ry , rz , fx , fy , fz , tau_coul , simu_cell , atom_dec 
  USE control,                  ONLY :  cutlongrange
  USE io,                  ONLY :  stdout , ionode 
  USE thermodynamic,            ONLY :  u_pol, u_coul_tot , vir_coul_tot 
  USE cell,                     ONLY :  kardir , dirkar
  USE time,                     ONLY :  fcoultimetot1 , fcoultimetot3

  implicit none

  ! global 
  real(kind=dp) , intent(in)    :: mu     ( natm , 3 )
  real(kind=dp) , intent(out)   :: ef     ( natm , 3 )
  real(kind=dp) , intent(out)   :: efg    ( natm , 3 , 3 )
  real(kind=dp) , intent(out)   :: u_coul 
  real(kind=dp) , intent(out)   :: vir_coul
  real(kind=dp) , intent(out)   :: phi_coul ( natm ) 

  ! local 
  integer :: ia, ja , ierr , ncell 
  real(kind=dp), dimension(:), allocatable :: fx_coul , fy_coul , fz_coul
  real(kind=dp), dimension(:), allocatable :: phi_coul_qq , phi_coul_dd 
  real(kind=dp) :: u_coul_qq , u_coul_dd , u_coul_qd 
  real(kind=dp) :: vir_coul_qq , vir_coul_dd , vir_coul_qd
  real(kind=dp) :: cutsq
  real(kind=dp) :: rxi , ryi , rzi 
  real(kind=dp) :: rxj , ryj , rzj
  real(kind=dp) :: rxij , ryij , rzij
  real(kind=dp) :: sxij , syij , szij
  real(kind=dp) :: fxij , fyij , fzij
  real(kind=dp) :: qi , qj , qij
  real(kind=dp) :: muix , muiy , muiz
  real(kind=dp) :: mujx , mujy , mujz
  real(kind=dp) :: T
  real(kind=dp) :: Tx , Ty , Tz
  real(kind=dp) :: Txx , Tyy , Tzz , Txy , Txz , Tyz
  real(kind=dp) :: Txxx,  Tyyy,  Tzzz, Txxy, Txxz, Tyyx, Tyyz, Tzzx, Tzzy, Txyz
  real(kind=dp) :: d , d2 
  real(kind=dp) :: dm1 , dm3 , dm5 , dm7 

  real(kind=dp) :: ttt1 , ttt2  , ttt3
  logical          :: lcentralbox
  !debug 
  integer          :: it , ip
  real(kind=dp) :: mip    ( ntype , 3 )

  ttt1 = MPI_WTIME(ierr)

  ! =============================== 
  !         some constants
  ! =============================== 
  cutsq = cutlongrange * cutlongrange

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
  WRITE ( stdout , '(a)')          'debug: in multipole_DS '
  WRITE ( stdout , '(a,i8,a)')     'debug : rm_coul        ',rm_coul%ncmax,rm_coul%meshlabel
  do ia = 1 , natm
    WRITE ( stdout , '(a,f12.5)')  'debug : charge (atom)  ',qia(ia)
  enddo
  do it = 1 , ntype
    WRITE ( stdout , '(a,f12.5)')  'debug : charge (type ) ',qch(it)
  enddo
  do ia = 1 , natm
    WRITE ( stdout , '(a,3f12.5)') 'debug : dipole (atom)  ', mu ( ia , 1 ) , mu ( ia , 2 ) , mu ( ia , 3 )
  enddo
  do it = 1 , ntype
    WRITE ( stdout , '(a,3f12.5)') 'debug : dipole (type ) ',  mip ( it , 1 ) , mip ( it , 2 ) , mip ( it , 3 )
  enddo
  WRITE ( stdout , '(a,2i8)')      'debug : atom_dec ',atom_dec%istart , atom_dec%iend
  WRITE ( stdout , '(a,f20.5)')    'debug : cutsq          ',cutsq
  endif
  call print_config_sample(0,0)
#endif  


  allocate( fx_coul ( natm ), fy_coul ( natm ), fz_coul ( natm ) )
  allocate( phi_coul_qq ( natm ), phi_coul_dd ( natm ) )
  ef          = 0.0_dp
  efg         = 0.0_dp
  vir_coul    = 0.0_dp
  vir_coul_qq = 0.0_dp
  vir_coul_qd = 0.0_dp
  vir_coul_dd = 0.0_dp
  u_coul      = 0.0_dp
  u_coul_qq   = 0.0_dp
  u_coul_qd   = 0.0_dp
  u_coul_dd   = 0.0_dp
  fx_coul     = 0.0_dp
  fy_coul     = 0.0_dp
  fz_coul     = 0.0_dp
  phi_coul    = 0.0_dp
  phi_coul_qq = 0.0_dp
  phi_coul_dd = 0.0_dp
  
  ! ======================================
  !         cartesian to direct 
  ! ======================================
  CALL kardir ( natm , rx , ry , rz , simu_cell%B )

! ==================================================================================================
!  MAIN LOOP calculate EF(i) , EFG(i) , virial , potential, forces ... for each atom i parallelized
! ==================================================================================================
  atom : do ia = atom_dec%istart , atom_dec%iend

    rxi  = rx  ( ia )
    ryi  = ry  ( ia )
    rzi  = rz  ( ia )
    qi   = qia ( ia )
    muix = mu ( ia , 1 )  
    muiy = mu ( ia , 2 )  
    muiz = mu ( ia , 3 )  

    ! ==================================================
    ! sum over neighboring cells (see direct_sum_init )
    ! ==================================================
    cells: do ncell = 1 , rm_coul%ncmax

      lcentralbox=rm_coul%lcell(ncell).eq.0   
                   
      do ja = 1 , natm

        if ( ( .not.lcentralbox ) .or. ( lcentralbox .and. ja .ne. ia )  ) then
 
            rxj   = rx ( ja ) + rm_coul%boxxyz( 1 , ncell )
            ryj   = ry ( ja ) + rm_coul%boxxyz( 2 , ncell )
            rzj   = rz ( ja ) + rm_coul%boxxyz( 3 , ncell )
            sxij = rxi - rxj
            syij = ryi - ryj
            szij = rzi - rzj
            rxij = sxij * simu_cell%A(1,1) + syij * simu_cell%A(1,2) + szij * simu_cell%A(1,3)
            ryij = sxij * simu_cell%A(2,1) + syij * simu_cell%A(2,2) + szij * simu_cell%A(2,3)
            rzij = sxij * simu_cell%A(3,1) + syij * simu_cell%A(3,2) + szij * simu_cell%A(3,3)
            d2   = rxij * rxij + ryij * ryij + rzij * rzij

            if ( d2 .lt. cutsq .and. d2 .ne. 0.0_dp ) then 

              qj    = qia ( ja )
              qij   = qi * qj
              mujx  = mu ( ja , 1 )  
              mujy  = mu ( ja , 2 )  
              mujz  = mu ( ja , 3 )  
              d     = SQRT ( d2 )
              dm1   = 1.0_dp / d 
              dm3   = dm1 / d2
              dm5   = dm3 / d2
              dm7   = dm5 / d2
 
              ! multipole interaction tensor rank = 0 
              T  = dm1

              ! multipole interaction tensor rank = 1
              Tx = - rxij * dm3 
              Ty = - ryij * dm3 
              Tz = - rzij * dm3 

              ! multipole interaction tensor rank = 2
              Txx = ( 3.0_dp * rxij * rxij - d2 ) * dm5
              Tyy = ( 3.0_dp * ryij * ryij - d2 ) * dm5
              Tzz = ( 3.0_dp * rzij * rzij - d2 ) * dm5
              Txy = ( 3.0_dp * rxij * ryij      ) * dm5
              Txz = ( 3.0_dp * rxij * rzij      ) * dm5
              Tyz = ( 3.0_dp * ryij * rzij      ) * dm5

              ! multipole interaction tensor rank = 3  
              Txxx = ( - 5.0_dp * rxij * rxij * rxij +  3.0_dp * d2 * ( rxij ) ) * dm7 * 3.0_dp
              Tyyy = ( - 5.0_dp * ryij * ryij * ryij +  3.0_dp * d2 * ( ryij ) ) * dm7 * 3.0_dp
              Tzzz = ( - 5.0_dp * rzij * rzij * rzij +  3.0_dp * d2 * ( rzij ) ) * dm7 * 3.0_dp
              Txxy = ( - 5.0_dp * rxij * rxij * ryij +          d2 * ( ryij ) ) * dm7 * 3.0_dp
              Txxz = ( - 5.0_dp * rxij * rxij * rzij +          d2 * ( rzij ) ) * dm7 * 3.0_dp
              Tyyx = ( - 5.0_dp * ryij * ryij * rxij +          d2 * ( rxij ) ) * dm7 * 3.0_dp
              Tyyz = ( - 5.0_dp * ryij * ryij * rzij +          d2 * ( rzij ) ) * dm7 * 3.0_dp
              Tzzx = ( - 5.0_dp * rzij * rzij * rxij +          d2 * ( rxij ) ) * dm7 * 3.0_dp
              Tzzy = ( - 5.0_dp * rzij * rzij * ryij +          d2 * ( ryij ) ) * dm7 * 3.0_dp
              Txyz = ( - 5.0_dp * rxij * ryij * rzij                          ) * dm7 * 3.0_dp

              ! ===========================================================
              !                  charge-charge interaction
              ! ===========================================================

              ! potential 
              phi_coul_qq ( ia ) = phi_coul_qq ( ia ) + qj * T

              ! electrostatic energy
              u_coul_qq          = u_coul_qq          + qij * T

              ! electric field
              ef ( ia , 1 )      = ef ( ia , 1 ) - qj * Tx
              ef ( ia , 2 )      = ef ( ia , 2 ) - qj * Ty
              ef ( ia , 3 )      = ef ( ia , 3 ) - qj * Tz 
 
              ! forces 
              fxij = qij * Tx
              fyij = qij * Ty
              fzij = qij * Tz

              fx_coul ( ia ) = fx_coul ( ia ) - fxij
              fy_coul ( ia ) = fy_coul ( ia ) - fyij
              fz_coul ( ia ) = fz_coul ( ia ) - fzij

              ! virial
              vir_coul_qq =  vir_coul_qq - ( fxij * rxij + fyij * ryij + fzij * rzij )

              ! electric field gradient
              efg ( ia , 1 , 1 )  = efg ( ia , 1 , 1 ) - qj * Txx 
              efg ( ia , 2 , 2 )  = efg ( ia , 2 , 2 ) - qj * Tyy 
              efg ( ia , 3 , 3 )  = efg ( ia , 3 , 3 ) - qj * Tzz 
              efg ( ia , 1 , 2 )  = efg ( ia , 1 , 2 ) - qj * Txy 
              efg ( ia , 1 , 3 )  = efg ( ia , 1 , 3 ) - qj * Txz  
              efg ( ia , 2 , 3 )  = efg ( ia , 2 , 3 ) - qj * Tyz

              ! stress tensor
              tau_coul(1,1) = tau_coul(1,1) - rxij * fxij
              tau_coul(1,2) = tau_coul(1,2) - rxij * fyij
              tau_coul(1,3) = tau_coul(1,3) - rxij * fzij
              tau_coul(2,1) = tau_coul(2,1) - ryij * fxij
              tau_coul(2,2) = tau_coul(2,2) - ryij * fyij
              tau_coul(2,3) = tau_coul(2,3) - ryij * fzij
              tau_coul(3,1) = tau_coul(3,1) - rzij * fxij
              tau_coul(3,2) = tau_coul(3,2) - rzij * fyij
              tau_coul(3,3) = tau_coul(3,3) - rzij * fzij

              ! ===========================================================
              !                  dipole-dipole interaction
              ! ===========================================================

              ! potential 
              phi_coul_dd ( ia ) = phi_coul_dd ( ia ) - ( Tx * mujx +  Ty * mujy + Tz * mujz )

              ! electrostatic energy
              u_coul_dd = u_coul_dd - ( muix * Txx * mujx + &
                                        muix * Txy * mujy + &
                                        muix * Txz * mujz + &
                                        muiy * Txy * mujx + &
                                        muiy * Tyy * mujy + &
                                        muiy * Tyz * mujz + &
                                        muiz * Txz * mujx + &
                                        muiz * Tyz * mujy + &
                                        muiz * Tzz * mujz )

              ! electric field
              ef ( ia , 1 ) = ef ( ia , 1 ) + ( Txx * mujx + Txy * mujy + Txz * mujz ) 
              ef ( ia , 2 ) = ef ( ia , 2 ) + ( Txy * mujx + Tyy * mujy + Tyz * mujz ) 
              ef ( ia , 3 ) = ef ( ia , 3 ) + ( Txz * mujx + Tyz * mujy + Tzz * mujz ) 

              ! forces 
              fxij = ( muix * Txxx * mujx         + &
                       muix * Txxy * mujy * 2.0_dp + &
                       muix * Txxz * mujz * 2.0_dp + &
                       muiy * Txyz * mujz * 2.0_dp + &   
                       muiy * Tyyx * mujy         + &    
                       muiz * Tzzx * mujz )

              fyij = ( muiy * Tyyy * mujy         + &
                       muix * Tyyx * mujy * 2.0_dp + &
                       muiy * Tyyz * mujz * 2.0_dp + &
                       muiy * Txyz * mujz * 2.0_dp + &     
                       muix * Txxy * mujx         + &
                       muiz * Tzzy * mujz )

              fzij = ( muiz * Tzzz * mujz         + &
                       muiz * Tzzx * mujx * 2.0_dp + &
                       muiy * Tzzy * mujz * 2.0_dp + &
                       muiy * Txyz * mujx * 2.0_dp + &
                       muix * Txxz * mujx         + &
                       muiy * Tyyz * mujy )

              fx_coul ( ia ) = fx_coul ( ia ) + fxij
              fy_coul ( ia ) = fy_coul ( ia ) + fyij
              fz_coul ( ia ) = fz_coul ( ia ) + fzij

              ! virial
              vir_coul_dd =  vir_coul_dd + ( fxij * rxij + fyij * ryij + fzij * rzij )
 
              ! electric field gradient
              efg ( ia , 1 , 1 ) = efg ( ia , 1 , 1 ) + ( Txxx * mujx + Txxy * mujy + Txxz * mujz ) 
              efg ( ia , 2 , 2 ) = efg ( ia , 2 , 2 ) + ( Tyyx * mujx + Tyyy * mujy + Tyyz * mujz ) 
              efg ( ia , 3 , 3 ) = efg ( ia , 3 , 3 ) + ( Tzzx * mujx + Tzzy * mujy + Tzzz * mujz ) 
              efg ( ia , 1 , 2 ) = efg ( ia , 1 , 2 ) + ( Txxy * mujx + Tyyx * mujy + Txyz * mujz ) 
              efg ( ia , 1 , 3 ) = efg ( ia , 1 , 3 ) + ( Txxz * mujx + Txyz * mujy + Tzzx * mujz ) 
              efg ( ia , 2 , 3 ) = efg ( ia , 2 , 3 ) + ( Txyz * mujx + Tyyz * mujy + Tzzy * mujz ) 

              ! stress tensor
              tau_coul(1,1) = tau_coul(1,1) + rxij * fxij
              tau_coul(1,2) = tau_coul(1,2) + rxij * fyij
              tau_coul(1,3) = tau_coul(1,3) + rxij * fzij
              tau_coul(2,1) = tau_coul(2,1) + ryij * fxij
              tau_coul(2,2) = tau_coul(2,2) + ryij * fyij
              tau_coul(2,3) = tau_coul(2,3) + ryij * fzij
              tau_coul(3,1) = tau_coul(3,1) + rzij * fxij
              tau_coul(3,2) = tau_coul(3,2) + rzij * fyij
              tau_coul(3,3) = tau_coul(3,3) + rzij * fzij


              ! ===========================================================
              !                  charge-dipole interaction
              ! ===========================================================

              ! electrostatic energy
              u_coul_qd = u_coul_qd - ( qi * ( Tx * mujx + Ty * mujy + Tz * mujz ) ) + ( qj * ( Tx * muix + Ty * muiy + Tz * muiz ) ) 

              ! forces 
              fxij = qi * ( Txx * mujx + Txy * mujy + Txz * mujz ) - &
                     qj * ( Txx * muix + Txy * muiy + Txz * muiz ) 

              fyij = qi * ( Txy * mujx + Tyy * mujy + Tyz * mujz ) - &
                     qj * ( Txy * muix + Tyy * muiy + Tyz * muiz ) 

              fzij = qi * ( Txz * mujx + Tyz * mujy + Tzz * mujz ) - &
                     qj * ( Txz * muix + Tyz * muiy + Tzz * muiz ) 


              fx_coul ( ia ) = fx_coul ( ia ) + fxij
              fy_coul ( ia ) = fy_coul ( ia ) + fyij
              fz_coul ( ia ) = fz_coul ( ia ) + fzij

              ! virial
              vir_coul_qd =  vir_coul_qd + ( fxij * rxij + fyij * ryij + fzij * rzij )

              ! stress tensor
              tau_coul(1,1) = tau_coul(1,1) + rxij * fxij
              tau_coul(1,2) = tau_coul(1,2) + rxij * fyij
              tau_coul(1,3) = tau_coul(1,3) + rxij * fzij
              tau_coul(2,1) = tau_coul(2,1) + ryij * fxij
              tau_coul(2,2) = tau_coul(2,2) + ryij * fyij
              tau_coul(2,3) = tau_coul(2,3) + ryij * fzij
              tau_coul(3,1) = tau_coul(3,1) + rzij * fxij
              tau_coul(3,2) = tau_coul(3,2) + rzij * fyij
              tau_coul(3,3) = tau_coul(3,3) + rzij * fzij

            endif

        endif

      enddo

    enddo cells

  enddo atom 

  ttt2 = MPI_WTIME(ierr)
  fcoultimetot1 = fcoultimetot1 + ( ttt2 - ttt1 )

  CALL MPI_ALL_REDUCE_DOUBLE_SCALAR ( u_coul_qq )
  CALL MPI_ALL_REDUCE_DOUBLE_SCALAR ( u_coul_dd )
  CALL MPI_ALL_REDUCE_DOUBLE_SCALAR ( u_coul_qd )
  CALL MPI_ALL_REDUCE_DOUBLE_SCALAR ( vir_coul_qq )
  CALL MPI_ALL_REDUCE_DOUBLE_SCALAR ( vir_coul_qd )
  CALL MPI_ALL_REDUCE_DOUBLE_SCALAR ( vir_coul_dd )

  CALL MPI_ALL_REDUCE_DOUBLE ( fx_coul , natm )
  CALL MPI_ALL_REDUCE_DOUBLE ( fy_coul , natm )
  CALL MPI_ALL_REDUCE_DOUBLE ( fz_coul , natm )

  CALL MPI_ALL_REDUCE_DOUBLE ( phi_coul_qq , natm )
  CALL MPI_ALL_REDUCE_DOUBLE ( phi_coul_dd , natm )

  CALL MPI_ALL_REDUCE_DOUBLE ( ef ( : , 1 )  , natm )
  CALL MPI_ALL_REDUCE_DOUBLE ( ef ( : , 2 )  , natm )
  CALL MPI_ALL_REDUCE_DOUBLE ( ef ( : , 3 )  , natm )

  CALL MPI_ALL_REDUCE_DOUBLE ( efg ( : , 1 , 1 ) , natm )
  CALL MPI_ALL_REDUCE_DOUBLE ( efg ( : , 2 , 2 ) , natm )
  CALL MPI_ALL_REDUCE_DOUBLE ( efg ( : , 3 , 3 ) , natm )
  CALL MPI_ALL_REDUCE_DOUBLE ( efg ( : , 1 , 2 ) , natm )
  CALL MPI_ALL_REDUCE_DOUBLE ( efg ( : , 1 , 3 ) , natm )
  CALL MPI_ALL_REDUCE_DOUBLE ( efg ( : , 2 , 3 ) , natm )

  CALL MPI_ALL_REDUCE_DOUBLE ( tau_coul ( 1 , : ) , 3 )
  CALL MPI_ALL_REDUCE_DOUBLE ( tau_coul ( 2 , : ) , 3 )
  CALL MPI_ALL_REDUCE_DOUBLE ( tau_coul ( 3 , : ) , 3 )

  ! EFG is symmetric
  ! not needed ... just for consistency
  efg ( : , 2 , 1 ) = efg ( : , 1 , 2 )
  efg ( : , 3 , 1 ) = efg ( : , 1 , 3 )
  efg ( : , 3 , 2 ) = efg ( : , 2 , 3 )

  CALL MPI_BARRIER( MPI_COMM_WORLD , ierr )

  u_coul_qq   = u_coul_qq * 0.5_dp
  u_coul_dd   = u_coul_dd * 0.5_dp
  u_coul_qd   = u_coul_qd * 0.5_dp
  vir_coul_qq = - vir_coul_qq  * 0.5_dp
  vir_coul_qd = - vir_coul_qd  * 0.5_dp
  vir_coul_dd = - vir_coul_dd  * 0.5_dp

  fx = fx + fx_coul
  fy = fy + fy_coul
  fz = fz + fz_coul

  u_coul   = u_coul_qq   + u_coul_dd   + u_coul_qd + u_pol
  vir_coul = vir_coul_qq + vir_coul_dd + vir_coul_qd
  phi_coul = phi_coul_qq + phi_coul_dd
  tau_coul = - tau_coul / simu_cell%omega * 0.5_dp

#ifdef debug2 
  if ( ionode ) then
    blankline(stdout)
    WRITE ( stdout , '(a)' )     'Electric field at atoms : '
    blankline(stdout)
    do ia = 1 , natm
      WRITE ( stdout , '(i5,a3,a,3f18.10)' ) &
      ia,atype(ia),' Efield  = ', ef ( ia , 1)  , ef ( ia , 2 ) , ef ( ia , 3 )
    enddo
    blankline(stdout)

    blankline(stdout)
    WRITE ( stdout , '(a)' ) 'forces at atoms : '
    do ia = 1 , natm
      WRITE ( stdout , '(i5,a3,a,3f18.10)' ) &
      ia,atype(ia),' f       = ', fx ( ia )  , fy ( ia ) , fz ( ia )
    enddo
    blankline(stdout)
    blankline(stdout)
    WRITE ( stdout , '(a)' ) ' potential : '
    do ia = 1 , natm
      WRITE ( stdout , '(i5,a3,3(a,f18.10))' ) &
      ia,atype(ia),' phi_tot = ', phi_coul ( ia )  , ' phi_qq = ', phi_coul_qq (ia) , ' phi_dd = ', phi_coul_dd (ia)
    enddo
    blankline(stdout)
    WRITE ( stdout , '(a)' ) 'Energy and virial : '
    WRITE ( stdout , '(5(a,f18.10))' ) &
    ' u_coul_tot   = ', u_coul  ,' u_coul_qq   = ',u_coul_qq      ,' u_coul_qd   = ', u_coul_qd     ,' u_coul_dd   = ',u_coul_dd   ,' u_pol = ', u_pol
    WRITE ( stdout , '(4(a,f18.10))' ) &
    ' vir_coul_tot = ',vir_coul ,' vir_coul_qq = ',vir_coul_qq    ,' vir_coul_qd = ', vir_coul_qd   ,' vir_coul_dd = ',vir_coul_dd
  endif
  CALL print_tensor( tau_coul ( : , : )     , 'TAU_COUL' )
  CALL print_tensor( efg ( 1 , : , : ) , 'TOTAL   ' )
#endif

  deallocate( fx_coul , fy_coul , fz_coul )
  deallocate( phi_coul_qq , phi_coul_dd  )

  u_coul_tot = u_coul
  vir_coul_tot = vir_coul

  ttt3 = MPI_WTIME(ierr)
  fcoultimetot3 = fcoultimetot3 + ( ttt3 - ttt2 )

  ! ======================================
  !         cartesian to direct 
  ! ======================================
  CALL dirkar ( natm , rx , ry , rz , simu_cell%A )

  return

END SUBROUTINE multipole_DS

! *********************** SUBROUTINE multipole_ES ******************************
!> \brief
!! This subroutine calculates electric field, electric field gradient, 
!! potential energy, virial, electric potential and forces at ions in
!! a multipole expansion by Ewald summation
!
!> \param[in]  mu electric dipole at ions
!> \param[out] ef electric field
!> \param[out] efg electric field gradient
!> \param[out] u_coul coulombic potential energy 
!> \param[out] virial_coul coulombic virial 
!> \param[out] phi_coul electric potential at ions 
!
!> \todo
!! make it more condensed 
! ******************************************************************************
SUBROUTINE multipole_ES ( ef , efg , mu , u_coul , vir_coul , phi_coul )

  USE control,                  ONLY :  lsurf , cutshortrange 
  USE config,                   ONLY :  natm , ntype , natmi , atype , &
                                        rx , ry , rz , fx , fy , fz , tau_coul , qia , simu_cell , atom_dec
  USE constants,                ONLY :  imag , pi , piroot , tpi , fpi
  USE io,                  ONLY :  ionode , stdout 
  USE thermodynamic,            ONLY :  u_pol, u_coul_tot , vir_coul_tot 
  USE cell,                     ONLY :  kardir , dirkar 
  USE time,                     ONLY :  strftimetot , fcoultimetot1 , fcoultimetot2 , fcoultimetot3 , fcoultimetot2_2
  USE kspace,                   ONLY :  struc_fact

  implicit none

  ! global 
  real(kind=dp)    :: ef     ( natm , 3 )
  real(kind=dp)    :: efg    ( natm , 3 , 3 )
  real(kind=dp)    :: mu     ( natm , 3 )
  real(kind=dp)    :: tau_dir( 3 , 3 )
  real(kind=dp)    :: tau_rec( 3 , 3 )
  real(kind=dp)    :: u_coul 
  real(kind=dp)    :: vir_coul 
  real(kind=dp)    :: phi_coul ( natm )  

  ! local 
  integer             :: ia , ja , ik , it , ip , ierr
  real(kind=dp), dimension(:)    , allocatable :: fx_coul , fy_coul , fz_coul
  real(kind=dp), dimension(:)    , allocatable :: fx_dir , fy_dir , fz_dir
  real(kind=dp), dimension(:)    , allocatable :: fx_rec , fy_rec , fz_rec
  real(kind=dp), dimension(:)    , allocatable :: fx_surf , fy_surf , fz_surf
  real(kind=dp), dimension(:)    , allocatable :: phi_dir , phi_rec , phi_surf , phi_self
  real(kind=dp), dimension(:,:)  , allocatable :: ef_dir , ef_rec , ef_surf , ef_self
  real(kind=dp), dimension(:,:,:), allocatable :: efg_dir , efg_rec , efg_self

  real(kind=dp) :: u_dir , u_rec , u_surf , u_self
  real(kind=dp) :: u_surf_qq , u_surf_qd , u_surf_dd
  real(kind=dp) :: vir_dir , vir_rec , vir_surf , vir_self

  real(kind=dp) :: mip    ( ntype , 3 )

  real(kind=dp) :: cutsq
  real(kind=dp) :: rxi  , ryi  , rzi
  real(kind=dp) :: rxj  , ryj  , rzj
  real(kind=dp) :: kx   , ky   , kz 
  real(kind=dp) :: rxij , ryij , rzij
  real(kind=dp) :: sxij , syij , szij
  real(kind=dp) :: fxij , fyij , fzij
  real(kind=dp) :: kk   , kri  , Ak , kcoe
  real(kind=dp) :: qi   , qj   , qij
  real(kind=dp) :: muix , muiy , muiz
  real(kind=dp) :: mujx , mujy , mujz
  complex(kind=dp)   :: rhon , rhonmu , carg 
  real(kind=dp) :: str  , recarg , recargi1 
  real(kind=dp) :: expon , F0 , F1 , F2 , F3 
  real(kind=dp) :: k_dot_mu
  real(kind=dp) :: T
  real(kind=dp) :: Tx , Ty , Tz
  real(kind=dp) :: Txx , Tyy , Tzz , Txy , Txz , Tyz
  real(kind=dp) :: Txxx,  Tyyy,  Tzzz, Txxy, Txxz, Tyyx, Tyyz, Tzzx, Tzzy, Txyz
  real(kind=dp) :: d , d2 , d3  , d5 
  real(kind=dp) :: dm1 , dm3 , dm5 , dm7 
  real(kind=dp) :: alpha2 , alpha3 , alpha5 
  real(kind=dp) :: selfa , selfa2 
  real(kind=dp), external :: errfc
  real(kind=dp) :: qtot ( 3 ) , qsq , mutot ( 3 ) , musq , qmu_sum ( 3 ) 
  real(kind=dp) :: tpi_V , tpi_3V , fpi_V , fpi_3V 
  real(kind=dp) :: ttt1 , ttt2  , ttt3 , ttt4 , ttt5 , ttt6 , ttt7 , ttt_2

  ttt1 = MPI_WTIME(ierr)

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

  ! =============================
  !  total charge / moment / square
  ! =============================
  mutot = 0.0_dp
  musq  = 0.0_dp
  qtot  = 0.0_dp
  qsq   = 0.0_dp
  do ia = 1 , natm
    mutot ( 1 ) = mutot ( 1 ) + mu ( ia , 1 )
    mutot ( 2 ) = mutot ( 2 ) + mu ( ia , 2 )
    mutot ( 3 ) = mutot ( 3 ) + mu ( ia , 3 )
    musq = musq + ( mu ( ia , 1 ) * mu ( ia , 1 ) +  mu ( ia , 2 ) * mu ( ia , 2 ) +  mu ( ia , 3 ) * mu ( ia , 3 ) )
    qtot ( 1 ) = qtot ( 1 ) + qia ( ia ) * rx ( ia ) 
    qtot ( 2 ) = qtot ( 2 ) + qia ( ia ) * ry ( ia ) 
    qtot ( 3 ) = qtot ( 3 ) + qia ( ia ) * rz ( ia ) 
    qsq = qsq + qia ( ia ) * qia ( ia )
  enddo
  qmu_sum ( 1 ) = qtot ( 1 ) + mutot ( 1 ) 
  qmu_sum ( 2 ) = qtot ( 2 ) + mutot ( 2 ) 
  qmu_sum ( 3 ) = qtot ( 3 ) + mutot ( 3 ) 


#ifdef debug_multipole_ES
  if ( ionode ) then
  WRITE( stdout , '(a)')          'debug : in engforce_charge_ES'
  WRITE( stdout , '(a,i8,a,a)')   'debug : km_coul       ',km_coul%nk,' ',km_coul%meshlabel
  do ia = 1 , natm
    WRITE( stdout , '(a,f12.5)')  'debug : charge (atom) ',qia(ia)
  enddo
  do it = 1 , ntype
    WRITE( stdout , '(a,f12.5)')  'debug : charge (type) ',qch(it)
  enddo
  do ia = 1 , natm
    WRITE( stdout , '(a,3f12.5)') 'debug : dipole (atom) ', mu ( ia , 1 ) , mu ( ia , 2 ) , mu ( ia , 3 )
  enddo
  do it = 1 , ntype
    WRITE( stdout , '(a,3f12.5)') 'debug : dipole (type) ', mip ( it , 1 ) , mip ( it , 2 ) , mip ( it , 3 )
  enddo
  WRITE( stdout , '(a,2i8)')      'debug : atom decompositon istart,iend ', atom_dec%istart , atom_dec%iend
  WRITE( stdout , '(a,f20.5)')    'debug : alphaES       ', alphaES
  endif
  call print_config_sample(0,0)
#endif 


  allocate( ef_dir ( natm , 3 ) , ef_rec ( natm , 3 ) , ef_surf ( natm , 3 ) , ef_self ( natm , 3 ) ) 
  allocate( efg_dir ( natm , 3 , 3 ) , efg_rec ( natm , 3 , 3 ) , efg_self ( natm , 3 , 3 ) ) 
  allocate( fx_coul (natm) , fy_coul (natm) , fz_coul (natm) )
  allocate( fx_dir (natm) , fy_dir (natm) , fz_dir (natm) )
  allocate( fx_rec (natm) , fy_rec (natm) , fz_rec (natm) )
  allocate( fx_surf (natm) , fy_surf (natm) , fz_surf (natm) )
  allocate( phi_dir ( natm ) , phi_rec ( natm ) , phi_surf ( natm ) , phi_self ( natm ) )

#ifdef debug_multipole_ES
  WRITE( stdout , '(a)')   'debug : in multipole_ES after allocating '
#endif

  ef_dir   = 0.0_dp
  ef_rec   = 0.0_dp
  ef_surf  = 0.0_dp
  efg_dir  = 0.0_dp
  efg_rec  = 0.0_dp
  efg_self = 0.0_dp
  !fx   = 0.0_dp
  !fy   = 0.0_dp
  !fz   = 0.0_dp
  fx_dir   = 0.0_dp
  fy_dir   = 0.0_dp
  fz_dir   = 0.0_dp
  fx_rec   = 0.0_dp
  fy_rec   = 0.0_dp
  fz_rec   = 0.0_dp
  fx_surf  = 0.0_dp
  fy_surf  = 0.0_dp
  fz_surf  = 0.0_dp
  phi_dir  = 0.0_dp
  phi_rec  = 0.0_dp
  phi_surf = 0.0_dp
  phi_self = 0.0_dp
  u_dir    = 0.0_dp 
  u_rec    = 0.0_dp 
  u_self   = 0.0_dp 
  u_surf   = 0.0_dp 
  vir_dir  = 0.0_dp 
  vir_rec  = 0.0_dp 
  vir_surf = 0.0_dp 
  vir_self = 0.0_dp 
  tau_dir  = 0.0_dp
  tau_rec  = 0.0_dp

  ! =====================
  ! facteur de structure 
  ! =====================
  CALL struc_fact ( km_coul )
  ttt2 = MPI_WTIME(ierr)
  strftimetot = strftimetot + (ttt2 - ttt1)

#ifdef debug_multipole_ES
  WRITE( stdout , '(a)')   'debug : in multipole_ES after struct fact '
#endif
  ! ======================================
  !         cartesian to direct 
  ! ======================================
  CALL kardir ( natm , rx , ry , rz , simu_cell%B )

  ! =================
  !  some constants 
  ! =================
  tpi_V  = tpi    / simu_cell%omega  ! 2pi / V
  tpi_3V = tpi_V  / 3.0_dp  ! 2pi / 3V 
  fpi_V  = tpi_V  * 2.0_dp  ! 4pi / V
  fpi_3V = tpi_3V * 2.0_dp  ! 4pi / 3V
  alpha2 = alphaES * alphaES
  alpha3 = alpha2  * alphaES
  alpha5 = alpha3  * alpha2
  selfa  = alphaES / piroot
  selfa2 = 2.0_dp * selfa * alpha2 / 3.0_dp
  cutsq  = cutshortrange * cutshortrange

  ! ==============================================
  !        direct space part
  ! ==============================================

  do ia = atom_dec%istart , atom_dec%iend
    rxi = rx(ia)
    ryi = ry(ia)
    rzi = rz(ia)
    qi  = qia(ia)
    muix = mu ( ia , 1 )
    muiy = mu ( ia , 2 )
    muiz = mu ( ia , 3 )

    do ja = 1, natm

      if (ja .ne. ia ) then

        qj   = qia(ja)
        mujx = mu ( ja , 1 )
        mujy = mu ( ja , 2 )
        mujz = mu ( ja , 3 )
        qij  = qi * qj 
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
      !if ( d2 .lt. cutsq ) then
        d   = SQRT ( d2 )
        d3  = d2 * d 
        d5  = d3 * d2 
        dm1 = 1.0_dp / d
        dm3 = dm1 / d2 
        dm5 = dm3 / d2 
        dm7 = dm5 / d2 
 
        expon = EXP ( - alpha2 * d2 )    / piroot
        F0    = errfc( alphaES * d )
        F1    = F0 + 2.0_dp * alphaES * d  * expon 
        F2    = F1 + 4.0_dp * alpha3  * d3 * expon / 3.0_dp
        F3    = F2 + 8.0_dp * alpha5  * d5 * expon / 15.0_dp

        ! multipole interaction tensor rank = 0 
        T  = dm1 * F0

        ! multipole interaction tensor rank = 1
        Tx = - rxij * F1 * dm3
        Ty = - ryij * F1 * dm3
        Tz = - rzij * F1 * dm3

        ! multipole interaction tensor rank = 2
        Txx = ( 3.0_dp * rxij * rxij * F2 - d2 * F1 ) * dm5 
        Tyy = ( 3.0_dp * ryij * ryij * F2 - d2 * F1 ) * dm5 
        Tzz = ( 3.0_dp * rzij * rzij * F2 - d2 * F1 ) * dm5
        Txy = ( 3.0_dp * rxij * ryij * F2           ) * dm5
        Txz = ( 3.0_dp * rxij * rzij * F2           ) * dm5 
        Tyz = ( 3.0_dp * ryij * rzij * F2           ) * dm5 
 
        ! multipole interaction tensor rank = 3  
        Txxx = ( - 5.0_dp * rxij * rxij * rxij * F3 +  3.0_dp * d2 * ( rxij ) * F2 ) * dm7 * 3.0_dp
        Tyyy = ( - 5.0_dp * ryij * ryij * ryij * F3 +  3.0_dp * d2 * ( ryij ) * F2 ) * dm7 * 3.0_dp
        Tzzz = ( - 5.0_dp * rzij * rzij * rzij * F3 +  3.0_dp * d2 * ( rzij ) * F2 ) * dm7 * 3.0_dp
        Txxy = ( - 5.0_dp * rxij * rxij * ryij * F3 +           d2 * ( ryij ) * F2 ) * dm7 * 3.0_dp
        Txxz = ( - 5.0_dp * rxij * rxij * rzij * F3 +           d2 * ( rzij ) * F2 ) * dm7 * 3.0_dp
        Tyyx = ( - 5.0_dp * ryij * ryij * rxij * F3 +           d2 * ( rxij ) * F2 ) * dm7 * 3.0_dp
        Tyyz = ( - 5.0_dp * ryij * ryij * rzij * F3 +           d2 * ( rzij ) * F2 ) * dm7 * 3.0_dp
        Tzzx = ( - 5.0_dp * rzij * rzij * rxij * F3 +           d2 * ( rxij ) * F2 ) * dm7 * 3.0_dp
        Tzzy = ( - 5.0_dp * rzij * rzij * ryij * F3 +           d2 * ( ryij ) * F2 ) * dm7 * 3.0_dp
        Txyz = ( - 5.0_dp * rxij * ryij * rzij * F3                                ) * dm7 * 3.0_dp

        ! ===========================================================
        !                  charge-charge interaction
        ! ===========================================================

        ! potential 
        phi_dir ( ia ) = phi_dir ( ia ) + qj * T 
 
        ! energy
        u_dir = u_dir + qij * T 
 
        ! electric field
        ef_dir ( ia , 1 ) = ef_dir ( ia , 1 ) - qj * Tx 
        ef_dir ( ia , 2 ) = ef_dir ( ia , 2 ) - qj * Ty 
        ef_dir ( ia , 3 ) = ef_dir ( ia , 3 ) - qj * Tz 
 
        ! forces
        fxij = qij * Tx
        fyij = qij * Ty
        fzij = qij * Tz

        fx_dir ( ia ) = fx_dir ( ia ) - fxij
        fy_dir ( ia ) = fy_dir ( ia ) - fyij
        fz_dir ( ia ) = fz_dir ( ia ) - fzij

        ! electric field gradient 
        efg_dir ( ia , 1 , 1 ) = efg_dir( ia , 1 , 1 ) - qj * Txx 
        efg_dir ( ia , 2 , 2 ) = efg_dir( ia , 2 , 2 ) - qj * Tyy
        efg_dir ( ia , 3 , 3 ) = efg_dir( ia , 3 , 3 ) - qj * Tzz
        efg_dir ( ia , 1 , 2 ) = efg_dir( ia , 1 , 2 ) - qj * Txy
        efg_dir ( ia , 1 , 3 ) = efg_dir( ia , 1 , 3 ) - qj * Txz
        efg_dir ( ia , 2 , 3 ) = efg_dir( ia , 2 , 3 ) - qj * Tyz

        ! stress tensor
        tau_dir(1,1) = tau_dir(1,1) - rxij * fxij
        tau_dir(1,2) = tau_dir(1,2) - rxij * fyij
        tau_dir(1,3) = tau_dir(1,3) - rxij * fzij
        tau_dir(2,1) = tau_dir(2,1) - ryij * fxij
        tau_dir(2,2) = tau_dir(2,2) - ryij * fyij
        tau_dir(2,3) = tau_dir(2,3) - ryij * fzij
        tau_dir(3,1) = tau_dir(3,1) - rzij * fxij
        tau_dir(3,2) = tau_dir(3,2) - rzij * fyij
        tau_dir(3,3) = tau_dir(3,3) - rzij * fzij

        ! virial
        vir_dir = vir_dir - ( fxij * rxij + fyij * ryij + fzij * rzij )

        ! ===========================================================
        !                  dipole-dipole interaction
        ! ===========================================================
        
        ! potential 
        phi_dir ( ia ) = phi_dir ( ia ) - ( mujx * Tx + mujy * Ty + mujz * Ty  ) 
 
        ! energy
        u_dir = u_dir - ( muix * Txx * mujx + &   
                          muix * Txy * mujy + &
                          muix * Txz * mujz + &
                          muiy * Txy * mujx + &
                          muiy * Tyy * mujy + &
                          muiy * Tyz * mujz + &
                          muiz * Txz * mujx + &
                          muiz * Tyz * mujy + &
                          muiz * Tzz * mujz )

        ! electric field
        ef_dir ( ia , 1 ) = ef_dir ( ia , 1 ) + ( Txx * mujx + Txy * mujy + Txz * mujz ) 
        ef_dir ( ia , 2 ) = ef_dir ( ia , 2 ) + ( Txy * mujx + Tyy * mujy + Tyz * mujz ) 
        ef_dir ( ia , 3 ) = ef_dir ( ia , 3 ) + ( Txz * mujx + Tyz * mujy + Tzz * mujz ) 

        ! forces
        fxij = ( muix * Txxx * mujx + &
                 muix * Txxy * mujy + &
                 muix * Txxz * mujz + &
                 muiy * Txxy * mujx + &
                 muiy * Tyyx * mujy + &
                 muiy * Txyz * mujz + &
                 muiz * Txxz * mujx + &
                 muiz * Txyz * mujy + &
                 muiz * Tzzx * mujz ) 

        fyij = ( muix * Txxy * mujx + &
                 muix * Tyyx * mujy + &
                 muix * Txyz * mujz + &
                 muiy * Tyyx * mujx + &
                 muiy * Tyyy * mujy + &
                 muiy * Tyyz * mujz + &
                 muiz * Txyz * mujx + &
                 muiz * Tyyz * mujy + &
                 muiz * Tzzy * mujz ) 

        fzij = ( muix * Txxz * mujx + &
                 muix * Txyz * mujy + &
                 muix * Tzzx * mujz + &
                 muiy * Txyz * mujx + &
                 muiy * Tyyz * mujy + &
                 muiy * Tzzy * mujz + &
                 muiz * Tzzx * mujx + &
                 muiz * Tzzy * mujy + &
                 muiz * Tzzz * mujz ) 

        fx_dir ( ia ) = fx_dir ( ia ) + fxij
        fy_dir ( ia ) = fy_dir ( ia ) + fyij
        fz_dir ( ia ) = fz_dir ( ia ) + fzij
 
        ! electric field gradient 
        efg_dir ( ia , 1 , 1 ) = efg_dir ( ia , 1 , 1 ) + ( Txxx * mujx + Txxy * mujy + Txxz * mujz ) 
        efg_dir ( ia , 2 , 2 ) = efg_dir ( ia , 2 , 2 ) + ( Tyyx * mujx + Tyyy * mujy + Tyyz * mujz )
        efg_dir ( ia , 3 , 3 ) = efg_dir ( ia , 3 , 3 ) + ( Tzzx * mujx + Tzzy * mujy + Tzzz * mujz )
        efg_dir ( ia , 1 , 2 ) = efg_dir ( ia , 1 , 2 ) + ( Txxy * mujx + Tyyx * mujy + Txyz * mujz ) 
        efg_dir ( ia , 1 , 3 ) = efg_dir ( ia , 1 , 3 ) + ( Txxz * mujx + Txyz * mujy + Tzzx * mujz )
        efg_dir ( ia , 2 , 3 ) = efg_dir ( ia , 2 , 3 ) + ( Txyz * mujx + Tyyz * mujy + Tzzy * mujz )
 
        ! stress tensor
        tau_dir(1,1) = tau_dir(1,1) - rxij * fxij
        tau_dir(1,2) = tau_dir(1,2) - rxij * fyij
        tau_dir(1,3) = tau_dir(1,3) - rxij * fzij
        tau_dir(2,1) = tau_dir(2,1) - ryij * fxij
        tau_dir(2,2) = tau_dir(2,2) - ryij * fyij
        tau_dir(2,3) = tau_dir(2,3) - ryij * fzij
        tau_dir(3,1) = tau_dir(3,1) - rzij * fxij
        tau_dir(3,2) = tau_dir(3,2) - rzij * fyij
        tau_dir(3,3) = tau_dir(3,3) - rzij * fzij

        ! virial
        vir_dir = vir_dir + ( fxij * rxij + fyij * ryij + fzij * rzij )

        ! ===========================================================
        !                  charge-dipole interaction
        ! ===========================================================
 
        ! electrostatic energy
        u_dir = u_dir - ( qi * ( Tx * mujx + Ty * mujy + Tz * mujz ) ) + ( qj * ( Tx * muix + Ty * muiy + Tz * muiz ) )
 
        ! forces 
        fxij = qi * ( Txx * mujx + Txy * mujy + Txz * mujz ) - &
               qj * ( Txx * muix + Txy * muiy + Txz * muiz )

        fyij = qi * ( Txy * mujx + Tyy * mujy + Tyz * mujz ) - &
               qj * ( Txy * muix + Tyy * muiy + Tyz * muiz )

        fzij = qi * ( Txz * mujx + Tyz * mujy + Tzz * mujz ) - &
               qj * ( Txz * muix + Tyz * muiy + Tzz * muiz )

        fx_dir ( ia ) = fx_dir ( ia ) + fxij
        fy_dir ( ia ) = fy_dir ( ia ) + fyij
        fz_dir ( ia ) = fz_dir ( ia ) + fzij

        ! stress tensor
        tau_dir(1,1) = tau_dir(1,1) + rxij * fxij
        tau_dir(1,2) = tau_dir(1,2) + rxij * fyij
        tau_dir(1,3) = tau_dir(1,3) + rxij * fzij
        tau_dir(2,1) = tau_dir(2,1) + ryij * fxij
        tau_dir(2,2) = tau_dir(2,2) + ryij * fyij
        tau_dir(2,3) = tau_dir(2,3) + ryij * fzij
        tau_dir(3,1) = tau_dir(3,1) + rzij * fxij
        tau_dir(3,2) = tau_dir(3,2) + rzij * fyij
        tau_dir(3,3) = tau_dir(3,3) + rzij * fzij

        ! virial
        vir_dir = vir_dir + ( fxij * rxij + fyij * ryij + fzij * rzij )
!       endif

      endif

    enddo

  enddo 

  ttt3 = MPI_WTIME(ierr)
  fcoultimetot1 = fcoultimetot1 + ( ttt3 - ttt2 )
  ttt_2 = 0.0_dp
  ! ======================================
  !         direct to cartesian
  ! ======================================
  CALL dirkar ( natm , rx , ry , rz , simu_cell%A )

  ! ==============================================
  !            reciprocal space part
  ! ==============================================
!  kpoint : do ik = 1 , km_coul%nk
  kpoint : do ik = km_coul%kpt_dec%istart, km_coul%kpt_dec%iend 
  
    ttt6 = MPI_WTIME(ierr)
    ! =================
    !   k-space  
    ! =================
    kx   = km_coul%kpt(1,ik)
    ky   = km_coul%kpt(2,ik)
    kz   = km_coul%kpt(3,ik)
    kk   = km_coul%kptk(ik)
    Ak   = EXP ( - kk * 0.25_dp / alpha2 ) / kk
    kcoe = 2.0_dp * ( 1.0_dp / alpha2  + 1.0_dp / kk )

    if (km_coul%kptk(ik) .eq. 0 ) cycle 
    !  WRITE ( stdout , * ) 'the sum should be done on k! =  0',ik
    ! ===============================
    !                              ---
    !  charge density in k-space ( \   q * facteur de structure  )
    !                              /__
    ! ===============================
    ! check multipole_efg_ES 
    rhon   = (0.0_dp, 0.0_dp)
    rhonmu = (0.0_dp, 0.0_dp)
    do ja = 1, natm
      k_dot_mu = ( mu ( ja , 1 ) * kx + mu ( ja , 2 ) * ky + mu ( ja , 3 ) * kz  ) 
      rhon = rhon + ( qia(ja) + imag * k_dot_mu ) * km_coul%strf ( ik , ja ) 
      rhonmu = rhonmu + imag * k_dot_mu * km_coul%strf ( ik , ja ) 
    enddo

    str = rhon * CONJG(rhon)

    u_rec   = u_rec   + Ak * str 
    vir_rec = vir_rec + Ak * str * ( 3.0_dp - kcoe * kk ) / 3.0_dp

    do ia = 1 , natm
      rxi = rx(ia)
      ryi = ry(ia)
      rzi = rz(ia)
      qi  = qia ( ia )
      muix = mu ( ia , 1 ) 
      muiy = mu ( ia , 2 ) 
      muiz = mu ( ia , 3 ) 
      kri = ( kx * rxi + ky * ryi + kz * rzi )
      k_dot_mu  =( muix * kx + muiy * ky + muiz * kz  )
      carg        = EXP  ( imag * kri )
      recarg      = DBLE ( CONJG ( rhon ) * carg )
      recargi1    = DBLE ( CONJG ( rhon ) * carg * imag )
      phi_rec ( ia ) = phi_rec ( ia ) + Ak * recarg 

      fxij = Ak * kx * recargi1 
      fyij = Ak * ky * recargi1 
      fzij = Ak * kz * recargi1  

      ef_rec ( ia , 1 ) = ef_rec ( ia , 1 ) - fxij 
      ef_rec ( ia , 2 ) = ef_rec ( ia , 2 ) - fyij
      ef_rec ( ia , 3 ) = ef_rec ( ia , 3 ) - fzij

      fx_rec ( ia ) = fx_rec ( ia ) - qi * fxij + Ak * k_dot_mu * recarg * kx 
      fy_rec ( ia ) = fy_rec ( ia ) - qi * fyij + Ak * k_dot_mu * recarg * ky 
      fz_rec ( ia ) = fz_rec ( ia ) - qi * fzij + Ak * k_dot_mu * recarg * kz 

      ! electric field gradient
      efg_rec ( ia , 1 , 1 ) = efg_rec ( ia , 1 , 1 ) +  Ak * kx * kx * recarg
      efg_rec ( ia , 2 , 2 ) = efg_rec ( ia , 2 , 2 ) +  Ak * ky * ky * recarg 
      efg_rec ( ia , 3 , 3 ) = efg_rec ( ia , 3 , 3 ) +  Ak * kz * kz * recarg 
      efg_rec ( ia , 1 , 2 ) = efg_rec ( ia , 1 , 2 ) +  Ak * kx * ky * recarg 
      efg_rec ( ia , 1 , 3 ) = efg_rec ( ia , 1 , 3 ) +  Ak * kx * kz * recarg 
      efg_rec ( ia , 2 , 3 ) = efg_rec ( ia , 2 , 3 ) +  Ak * ky * kz * recarg 

     ! stress tensor
      tau_rec(1,1) = tau_rec(1,1) + ( 1.0_dp - kcoe * kx * kx ) * str * ak
      tau_rec(1,2) = tau_rec(1,2) -            kcoe * kx * ky   * str * ak
      tau_rec(1,3) = tau_rec(1,3) -            kcoe * kx * kz   * str * ak
      tau_rec(2,1) = tau_rec(2,1) -            kcoe * ky * kx   * str * ak 
      tau_rec(2,2) = tau_rec(2,2) + ( 1.0_dp - kcoe * ky * ky ) * str * ak
      tau_rec(2,3) = tau_rec(2,3) -            kcoe * ky * kz   * str * ak
      tau_rec(3,1) = tau_rec(3,1) -            kcoe * kz * kx   * str * ak
      tau_rec(3,2) = tau_rec(3,2) -            kcoe * kz * ky   * str * ak
      tau_rec(3,3) = tau_rec(3,3) + ( 1.0_dp - kcoe * kz * kz ) * str * ak

    enddo

    ttt7 = MPI_WTIME(ierr)
    ttt_2 = ttt_2 + ( ttt7 - ttt6 )
  enddo kpoint
!  print*,ttt_2
  ttt_2 = ttt_2 / km_coul%nk
  fcoultimetot2_2 = fcoultimetot2_2 + ttt_2

  ttt4 = MPI_WTIME(ierr)
  fcoultimetot2 = fcoultimetot2 + ( ttt4 - ttt3 )

  CALL MPI_ALL_REDUCE_DOUBLE_SCALAR ( u_dir )
  CALL MPI_ALL_REDUCE_DOUBLE_SCALAR ( vir_dir )
  CALL MPI_ALL_REDUCE_DOUBLE_SCALAR ( u_rec )
  CALL MPI_ALL_REDUCE_DOUBLE_SCALAR ( vir_rec )

  CALL MPI_ALL_REDUCE_DOUBLE ( fx_dir , natm )
  CALL MPI_ALL_REDUCE_DOUBLE ( fy_dir , natm )
  CALL MPI_ALL_REDUCE_DOUBLE ( fz_dir , natm )
  CALL MPI_ALL_REDUCE_DOUBLE ( fx_rec , natm )
  CALL MPI_ALL_REDUCE_DOUBLE ( fy_rec , natm )
  CALL MPI_ALL_REDUCE_DOUBLE ( fz_rec , natm )

  CALL MPI_ALL_REDUCE_DOUBLE ( phi_dir , natm )
  CALL MPI_ALL_REDUCE_DOUBLE ( phi_rec , natm )

  CALL MPI_ALL_REDUCE_DOUBLE ( ef_dir ( : , 1 )  , natm )
  CALL MPI_ALL_REDUCE_DOUBLE ( ef_dir ( : , 2 )  , natm )
  CALL MPI_ALL_REDUCE_DOUBLE ( ef_dir ( : , 3 )  , natm )
  CALL MPI_ALL_REDUCE_DOUBLE ( ef_rec ( : , 1 )  , natm )
  CALL MPI_ALL_REDUCE_DOUBLE ( ef_rec ( : , 2 )  , natm )
  CALL MPI_ALL_REDUCE_DOUBLE ( ef_rec ( : , 3 )  , natm )

  CALL MPI_ALL_REDUCE_DOUBLE ( efg_dir ( : , 1 , 1 ) , natm )
  CALL MPI_ALL_REDUCE_DOUBLE ( efg_dir ( : , 2 , 2 ) , natm )
  CALL MPI_ALL_REDUCE_DOUBLE ( efg_dir ( : , 3 , 3 ) , natm )
  CALL MPI_ALL_REDUCE_DOUBLE ( efg_dir ( : , 1 , 2 ) , natm )
  CALL MPI_ALL_REDUCE_DOUBLE ( efg_dir ( : , 1 , 3 ) , natm )
  CALL MPI_ALL_REDUCE_DOUBLE ( efg_dir ( : , 2 , 3 ) , natm )
  CALL MPI_ALL_REDUCE_DOUBLE ( efg_rec ( : , 1 , 1 ) , natm )
  CALL MPI_ALL_REDUCE_DOUBLE ( efg_rec ( : , 2 , 2 ) , natm )
  CALL MPI_ALL_REDUCE_DOUBLE ( efg_rec ( : , 3 , 3 ) , natm )
  CALL MPI_ALL_REDUCE_DOUBLE ( efg_rec ( : , 1 , 2 ) , natm )
  CALL MPI_ALL_REDUCE_DOUBLE ( efg_rec ( : , 1 , 3 ) , natm )
  CALL MPI_ALL_REDUCE_DOUBLE ( efg_rec ( : , 2 , 3 ) , natm )

  CALL MPI_ALL_REDUCE_DOUBLE ( tau_dir ( 1 , : ) , 3 )
  CALL MPI_ALL_REDUCE_DOUBLE ( tau_dir ( 2 , : ) , 3 )
  CALL MPI_ALL_REDUCE_DOUBLE ( tau_dir ( 3 , : ) , 3 )
  CALL MPI_ALL_REDUCE_DOUBLE ( tau_rec ( 1 , : ) , 3 )
  CALL MPI_ALL_REDUCE_DOUBLE ( tau_rec ( 2 , : ) , 3 )
  CALL MPI_ALL_REDUCE_DOUBLE ( tau_rec ( 3 , : ) , 3 )

  ! ======================================================
  ! remark on the unit :
  ! 1/(4*pi*epislon_0) = 1 => epsilon_0 = 1/4pi
  ! ======================================================

  tau_dir =   tau_dir * 0.5_dp 
  u_dir   =   u_dir   * 0.5_dp
  vir_dir = - vir_dir * 0.5_dp

  tau_rec =   tau_rec * tpi_V / simu_cell%omega
  u_rec   =   u_rec   * tpi_V  
  vir_rec =   vir_rec * tpi_V 
  phi_rec =   phi_rec * fpi_V 
  ef_rec  =   ef_rec  * fpi_V 
  fx_rec  =   fx_rec  * fpi_V 
  fy_rec  =   fy_rec  * fpi_V 
  fz_rec  =   fz_rec  * fpi_V 
  efg_rec =   efg_rec * fpi_V  

  ! ====================================================== 
  !              Surface contribution 
  ! ====================================================== 
  
  ! electrostatic energy and virial
  ! qq
  u_surf_qq = qtot ( 1 ) * qtot ( 1 ) + qtot ( 2 ) * qtot ( 2 ) + qtot ( 3 ) * qtot ( 3 )
  u_surf_dd = mutot ( 1 ) * mutot ( 1 ) + mutot ( 2 ) * mutot ( 2 ) + mutot ( 3 ) * mutot ( 3 )
  u_surf_qd = 2.0_dp * ( qtot ( 1 ) * mutot ( 1 ) + qtot ( 2 ) * mutot ( 2 ) + qtot ( 3 ) * mutot ( 3 ) )
  u_surf = u_surf_qq + u_surf_qd + u_surf_dd

  u_surf = u_surf * tpi_3V 

  vir_surf = - ( u_surf_qq + 4.0_dp * u_surf_qd + 3.0_dp * u_surf_dd ) 
  vir_surf = vir_surf * tpi_3V 

  ! potential, field , forces ( no contrib to efg ) 
  do ia = 1 , natm
    phi_surf ( ia ) = rx ( ia ) * qmu_sum ( 1 ) + ry ( ia ) * qmu_sum ( 2 ) + rz ( ia ) * qmu_sum ( 3 ) 
    ef_surf ( ia , 1 ) = qmu_sum ( 1 )
    ef_surf ( ia , 2 ) = qmu_sum ( 2 )
    ef_surf ( ia , 3 ) = qmu_sum ( 3 )
    fx_surf ( ia ) = qia ( ia ) * qmu_sum ( 1 )
    fy_surf ( ia ) = qia ( ia ) * qmu_sum ( 2 )
    fz_surf ( ia ) = qia ( ia ) * qmu_sum ( 3 )
  enddo
  ef_surf  = - ef_surf  * fpi_3V
  phi_surf =   phi_surf * fpi_3V
  fx_surf  = - fx_surf  * fpi_3V
  fy_surf  = - fy_surf  * fpi_3V
  fz_surf  = - fz_surf  * fpi_3V

  ! ====================================================== 
  !              Self contribution 
  ! ====================================================== 

  ! electrostic energy 
  u_self   = - selfa * qsq - selfa2 * musq
!  vir_self = - 9.0_dp * pi / 2.0_dp * qsq / alpha2 / simu_cell%omega 
 
  ! potential , field , field gradient ( no contrib to forces ) 
  do ia = 1 , natm
   phi_self ( ia ) = - 2.0_dp * selfa * qia ( ia )
   ef_self( ia , 1 ) = 2.0_dp * selfa2 * mu ( ia , 1 ) 
   ef_self( ia , 2 ) = 2.0_dp * selfa2 * mu ( ia , 2 ) 
   ef_self( ia , 3 ) = 2.0_dp * selfa2 * mu ( ia , 3 ) 
   efg_self ( ia , 1 , 1 ) = - 2.0_dp * selfa2 * qia ( ia )
   efg_self ( ia , 2 , 2 ) = - 2.0_dp * selfa2 * qia ( ia ) 
   efg_self ( ia , 3 , 3 ) = - 2.0_dp * selfa2 * qia ( ia ) 
  enddo

  ! =====================================================
  !                     TOTAL
  ! =====================================================

  if ( lsurf ) then
    u_coul   = ( u_dir   + u_rec   + u_surf   + u_self  + u_pol  )
    vir_coul = ( vir_dir + vir_rec + vir_surf + vir_self         )
    phi_coul = ( phi_dir + phi_rec + phi_surf + phi_self         )
    ef       = ( ef_dir  + ef_rec  + ef_surf  + ef_self          )
    efg      = ( efg_dir + efg_rec            + efg_self         )
    tau_coul = ( tau_dir + tau_rec )
    fx       = fx + ( fx_rec + fx_dir + fx_surf )
    fy       = fy + ( fy_rec + fy_dir + fy_surf )
    fz       = fz + ( fz_rec + fz_dir + fz_surf )
  else
#ifdef debug 
    io_node WRITE ( stdout , '(a)' ) 'Surface contribution is not added to the different quantities ( see multipole_ES inf field.f90 ) '
#endif
    u_coul   = ( u_dir   + u_rec   + u_self  + u_pol  )
    vir_coul = ( vir_dir + vir_rec + vir_self         )
    phi_coul = ( phi_dir + phi_rec + phi_self         )
    ef       = ( ef_dir  + ef_rec  + ef_self          )
    efg      = ( efg_dir + efg_rec + efg_self         )
    tau_coul = ( tau_dir + tau_rec )
    fx       = fx + ( fx_rec + fx_dir                 )
    fy       = fy + ( fy_rec + fy_dir                 )
    fz       = fz + ( fz_rec + fz_dir                 )
  endif

  ! only for consistency 
  efg ( : , 2 , 1 ) = efg ( : , 1 , 2 )
  efg ( : , 3 , 1 ) = efg ( : , 1 , 3 )
  efg ( : , 3 , 2 ) = efg ( : , 2 , 3 )

#ifdef debug2 
  if ( ionode ) then
    blankline(stdout)
    WRITE ( stdout , '(a)' )     'Electric field at atoms : '
    do ia = 1 , natm
      WRITE ( stdout , '(i5,a3,a,3f18.10)' ) &
      ia,atype(ia),' Efield   = ', ef ( ia , 1)  , ef ( ia , 2 ) , ef ( ia , 3 )
    enddo
    do ia = 1 , natm
      WRITE ( stdout , '(i5,a3,a,3f18.10)' ) &
      ia,atype(ia),' ef_dir   = ', ef_dir ( ia , 1)  , ef_dir ( ia , 2 ) , ef_dir ( ia , 3 )
    enddo
    do ia = 1 , natm
      WRITE ( stdout , '(i5,a3,a,3f18.10)' ) &
      ia,atype(ia),' ef_rec   = ', ef_rec ( ia , 1)  , ef_rec ( ia , 2 ) , ef_rec ( ia , 3 )
    enddo
    do ia = 1 , natm
      WRITE ( stdout , '(i5,a3,a,3f18.10)' ) &
      ia,atype(ia),' ef_surf  = ', ef_surf ( ia , 1)  , ef_surf ( ia , 2 ) , ef_surf ( ia , 3 )
    enddo
    do ia = 1 , natm
      WRITE ( stdout , '(i5,a3,a,3f18.10)' ) &
      ia,atype(ia),' ef_self  = ', ef_self ( ia , 1)  , ef_self ( ia , 2 ) , ef_self ( ia , 3 )
    enddo

    blankline(stdout)
    blankline(stdout)

    WRITE ( stdout , '(a)' ) 'forces at atoms : '
    do ia = 1 , natm
      WRITE ( stdout , '(i5,a3,a,3f18.10)' ) &
      ia,atype(ia),' f       = ', fx ( ia )  , fy ( ia ) , fz ( ia )
    enddo
    blankline(stdout)
    WRITE ( stdout , '(a)' ) 'potential at atoms : '
    do ia = 1 , natm
      WRITE ( stdout , '(i5,a3,a,f18.10)' )  ia,atype(ia),' phi         = ', phi_coul ( ia )  
    enddo
    do ia = 1 , natm
      WRITE ( stdout , '(i5,a3,a,f18.10)' )  ia,atype(ia),' phi_dir     = ', phi_dir ( ia )  
    enddo
    do ia = 1 , natm
      WRITE ( stdout , '(i5,a3,a,f18.10)' )  ia,atype(ia),' phi_rec     = ', phi_rec ( ia )  
    enddo
    do ia = 1 , natm
      WRITE ( stdout , '(i5,a3,a,f18.10)' )  ia,atype(ia),' phi_surf    = ', phi_surf ( ia )  
    enddo
    do ia = 1 , natm
      WRITE ( stdout , '(i5,a3,a,f18.10)' )  ia,atype(ia),' phi_self    = ', phi_self ( ia )  
    enddo

    blankline(stdout)
    blankline(stdout)
    WRITE ( stdout , '(a)' ) 'Energy and virial : '
    WRITE ( stdout , '(5(a,f16.8))' ) ,'vir_dir    = ', vir_dir    ,' vir_rec    = ', vir_rec    ,' vir_surf    = ',vir_surf   ,' vir_self    = '   ,vir_self   ,'                           vir_coul    = ',vir_coul
    WRITE ( stdout , '(5(a,f16.8))' ) ,'phi_dir(1) = ', phi_dir(1) ,' phi_rec(1) = ', phi_rec(1) ,' phi_surf(1) = ',phi_surf(1),' phi_self(1) = '   ,phi_self(1),'                           phi_coul(1) = ',phi_coul (1)
    WRITE ( stdout , '(6(a,f16.8))' ) ,'u_dir      = ', u_dir      ,' u_rec      = ', u_rec      ,' u_surf      = ',u_surf     ,' u_self      = '   ,u_self     ,' u_pol  = ',u_pol,' u_coul      = ',u_coul

  endif

  CALL print_tensor( efg_dir  ( 1 , : , : ) , 'EFG_DIR ' )
  CALL print_tensor( efg_rec  ( 1 , : , : ) , 'EFG_REC ' )
  CALL print_tensor( efg_self ( 1 , : , : ) , 'EFG_SELF' )
  CALL print_tensor( efg      ( 1 , : , : ) , 'EFG_TOT ' )
  CALL print_tensor( tau_dir  ( : , : )     , 'TAU_DIR ' )
  CALL print_tensor( tau_rec  ( : , : )     , 'TAU_REC ' )
  CALL print_tensor( tau_coul ( : , : )     , 'TAU_COUL' )
#endif

  deallocate( ef_dir  , ef_rec , ef_surf , ef_self ) 
  deallocate( efg_dir , efg_rec , efg_self ) 
  deallocate( fx_coul , fy_coul , fz_coul )
  deallocate( fx_dir  , fy_dir  , fz_dir )
  deallocate( fx_rec  , fy_rec  , fz_rec  )
  deallocate( fx_surf , fy_surf , fz_surf )
  deallocate( phi_dir , phi_rec , phi_surf , phi_self )

  u_coul_tot = u_coul
  vir_coul_tot = vir_coul

  ttt5 = MPI_WTIME(ierr)
  fcoultimetot3 = fcoultimetot3 + ( ttt5 - ttt4 )

  return

END SUBROUTINE multipole_ES

! *********************** SUBROUTINE engforce_morse_pbc ************************
!
!> \brief
!! total potential energy forces for each atoms for a morse potential with 
!! periodic boundaries conditions, with or without vnlist
!! (lvnlist=.TRUE.OR.FALSE.)
!
!> \warning
!! not fully tested
!
! ******************************************************************************
SUBROUTINE engforce_morse_pbc 

  USE config,                   ONLY :  natm , rx , ry , rz , fx , fy , fz, atype , itype , list , point , ntype , simu_cell , atom_dec
  USE control,                  ONLY :  lvnlist  
  USE thermodynamic,            ONLY :  u_morse , vir_morse
  USE time,                     ONLY :  forcetimetot
  USE io,                  ONLY :  stdout , ionode
  USE cell,                     ONLY :  kardir , dirkar

  implicit none

  ! local
  integer          :: ia , ja , j1 , je , jb , ierr
  integer          :: p1 , p2
  real(kind=dp) :: rxi , ryi , rzi
  real(kind=dp) :: rxij , ryij , rzij
  real(kind=dp) :: sxij , syij , szij
  real(kind=dp) :: rijsq , rij , erh , erh2 
  real(kind=dp) :: wij , fxij , fyij , fzij
  real(kind=dp) :: forcetime1 , forcetime2 
  real(kind=dp) :: u , vir , u2
  real(kind=dp) :: expon 

#ifdef debug_morse
  if ( ionode ) then
  do ia=1, natm
    WRITE ( stdout , '(a,a)' )  'debug : atype ',atype(ia)
    WRITE ( stdout , '(a,i4)' ) 'debug : itype ',itype(ia)
  enddo
  endif
#endif

  forcetime1 = MPI_WTIME(ierr) ! timing info

  u   = 0.0_dp
  vir = 0.0_dp
  fx  = 0.0_dp
  fy  = 0.0_dp
  fz  = 0.0_dp

  if ( lvnlist ) CALL vnlistcheck 
  ! ======================================
  !         cartesian to direct 
  ! ======================================
  CALL kardir ( natm , rx , ry , rz , simu_cell%B )

  do ia = atom_dec%istart , atom_dec%iend
    rxi = rx ( ia )
    ryi = ry ( ia )
    rzi = rz ( ia )
    if ( lvnlist ) then
      jb = point( ia )
      je = point( ia + 1 ) - 1
    else
      jb = 1 
      je = natm !atom_dec%iend
    endif
    do j1 = jb, je
      if ( lvnlist ) then
        ja = list ( j1 )
      else
        ja = j1
      endif
      if ( ja .ne. ia ) then
        rxij = rxi - rx ( ja )
        ryij = ryi - ry ( ja )
        rzij = rzi - rz ( ja )
        sxij = rxij - nint ( rxij )
        syij = ryij - nint ( ryij )
        szij = rzij - nint ( rzij )
        rxij = sxij * simu_cell%A(1,1) + syij * simu_cell%A(1,2) + szij * simu_cell%A(1,3)
        ryij = sxij * simu_cell%A(2,1) + syij * simu_cell%A(2,2) + szij * simu_cell%A(2,3)
        rzij = sxij * simu_cell%A(3,1) + syij * simu_cell%A(3,2) + szij * simu_cell%A(3,3)
        rijsq= rxij * rxij + ryij * ryij + rzij * rzij
        p1   = itype ( ja )
        p2   = itype ( ia )
        if ( rijsq .lt. rcutsq(p1,p2) ) then
           rij   = SQRT ( rijsq ) 
           erh   = EXP ( - rij * rhomor(p1,p2) ) 
           expon = erh * rs(p1,p2) - 1.0_dp
           expon = expon * expon - 1.0_dp
           u     = u  + epsmor(p1,p2) * expon 
           expon = EXP( rhomor(p1,p2) * ( 1.0_dp - rij / sigmamor(p1,p2) )  ) 
           u2    = u2 + ( expon * ( expon - 2.0_dp ) ) * epsmor(p1,p2)
        !  if ( ia .eq. 1 )  print*,rij,erh,erh2,expon,u
!          if ( trunc .eq. 1 ) then
!            u =  u + epsp(p1,p2) * ( plj(p1,p2) * srq -qlj(p1,p2) * srp ) - uc(p1,p2)
!          endif
!          if ( trunc .eq. 2 ) then
!            u =  u + epsp(p1,p2) * ( plj(p1,p2) * srq -qlj(p1,p2) * srp ) + uc1(p1,p2) * rijsq - uc2(p1,p2)
!          endif
          wij  = fm(p1,p2) * ( erh - rs(p1,p2) * erh2 ) 
          vir  = vir + wij * rijsq
          fxij = wij * rxij
          fyij = wij * ryij
          fzij = wij * rzij
          fx ( ia ) = fx ( ia ) + fxij
          fy ( ia ) = fy ( ia ) + fyij
          fz ( ia ) = fz ( ia ) + fzij
          fx ( ja ) = fx ( ja ) - fxij
          fy ( ja ) = fy ( ja ) - fyij
          fz ( ja ) = fz ( ja ) - fzij
!TODO
!!! stresss tensor tau_nonb 
        endif
      endif
    enddo
  enddo
!  vir = vir/3.0_dp

  forcetime2 = MPI_WTIME(ierr) ! timing info
  forcetimetot = forcetimetot+(forcetime2-forcetime1)

  CALL MPI_ALL_REDUCE_DOUBLE_SCALAR ( u   )
  CALL MPI_ALL_REDUCE_DOUBLE_SCALAR ( u2  )
  CALL MPI_ALL_REDUCE_DOUBLE_SCALAR ( vir )

  CALL MPI_ALL_REDUCE_DOUBLE ( fx , natm )
  CALL MPI_ALL_REDUCE_DOUBLE ( fy , natm )
  CALL MPI_ALL_REDUCE_DOUBLE ( fz , natm )

  u_morse = u
  vir_morse = vir

  WRITE ( stdout ,'(2f24.12)') u_morse , u2 

  ! ======================================
  !         direct to cartesian
  ! ======================================
  CALL dirkar ( natm , rx , ry , rz , simu_cell%A )

  return

END SUBROUTINE engforce_morse_pbc

! *********************** SUBROUTINE moment_from_pola **************************
!
!> \brief
!! This routines evaluates the dipole moment induced by polarizabilities on atoms. 
!! The evaluation is done self-consistently starting from the the field due the
!! point charges only.
!
!> \param[out] mu_ind induced electric dipole from polarizabilities
!
!> \note
!! The stopping criteria is governed by conv_tol_ind
!
! ******************************************************************************
SUBROUTINE moment_from_pola ( mu_ind ) 

  USE io, ONLY : ionode , stdout 
  USE config,  ONLY : natm , atype , fx , fy , fz , ntype , dipia , qia , ntypemax
  USE control, ONLY : longrange
  USE thermodynamic, ONLY : u_pol

  implicit none

  ! global
  real(kind=dp) , intent (out) :: mu_ind ( natm , 3 ) 

  ! local
  integer :: ia , iscf , it 
  logical :: linduced
  real(kind=dp) :: Efield_old , diff_efield
  real(kind=dp) :: u_tmp , vir_tmp
  real(kind=dp) :: Efield( natm , 3 ) , Efield_stat ( natm , 3 ) , Efield_ind ( natm , 3 ) 
  real(kind=dp) :: ef_tmp( natm , 3 ) , efg_tmp ( natm , 3 , 3 ) , phi_tmp ( natm ) 
  real(kind=dp) :: qia_tmp ( natm )  , qch_tmp ( ntypemax ) 

  ! =========================================================
  !  Is there any polarizability ? if yes linduced = .TRUE.
  ! =========================================================
  linduced = .false.
  do it = 1 , ntype
    if ( lpolar ( it ) .eq. 1 ) linduced = .true.
  enddo
  if ( .not. linduced ) then
    return
  endif

  ! =============================================
  !  calculate static Efield ( charge + dipoles )
  ! =============================================
  Efield_stat = 0.0_dp

  ! =============================================
  !  coulombic energy , forces (field) and virial
  ! =============================================
  if ( longrange .eq. 'direct' ) CALL  multipole_DS ( ef_tmp , efg_tmp , dipia , u_tmp , vir_tmp , phi_tmp ) 
  if ( longrange .eq. 'ewald' )  CALL  multipole_ES ( ef_tmp , efg_tmp , dipia , u_tmp , vir_tmp , phi_tmp ) 

  Efield_stat = Efield_stat + ef_tmp

  fx      = 0.0_dp ; fy = 0.0_dp ; fz = 0.0_dp
  ef_tmp  = 0.0_dp
  efg_tmp = 0.0_dp
  u_tmp   = 0.0_dp
  vir_tmp = 0.0_dp
  phi_tmp = 0.0_dp

  ! =============================================
  !  init total Efield to static only
  ! =============================================
  Efield = Efield_stat

  io_node WRITE ( stdout , '(a)' ) 'calculate induced dipole SCF'

  if ( ionode ) then
    WRITE ( stdout , '(a)' ) 'We start from the static electric field (charge + static dipole)'
    blankline(stdout)
  endif
  diff_efield = 1000.0_dp
  Efield_old  = Efield(1,1)
  iscf = 0
  ! =========================
  !  charges are set to zero 
  ! =========================
  qch_tmp = qch
  qia_tmp = qia
  qch = 0.0_dp
  qia = 0.0_dp
  ! =============================
  !           SCF LOOP
  ! =============================
  do while ( diff_efield .gt. conv_tol_ind )

    iscf = iscf + 1

    ! ==========================================================
    !  calculate mu_ind from Efield = Efield_stat + Efield_ind
    ! ==========================================================
    CALL induced_moment ( Efield , mu_ind , u_pol )  ! Efield in ; mu_ind and u_pol out

#ifdef debug
   do ia =1 , natm
     WRITE ( stdout , '(a,3f12.5)' ) 'debug : induced moment from pola', mu_ind ( ia , 1 ),  mu_ind ( ia , 2 ) , mu_ind ( ia , 3 )
   enddo
   do ia =1 , natm
     WRITE ( stdout , '(a,3f12.5)' ) 'debug : Efield                  ', Efield ( ia , 1 ),  Efield ( ia , 2 ) , Efield ( ia , 3 )
   enddo
#endif

    ! ==========================================================
    !  calculate Efield_ind from mu_ind
    !  Efield_ind out , mu_ind in ==> charges and static dipoles = 0
    ! ==========================================================
    if ( longrange .eq. 'direct' ) CALL multipole_DS ( Efield_ind , efg_tmp , mu_ind , u_tmp , vir_tmp , phi_tmp ) 
    if ( longrange .eq. 'ewald' )  CALL multipole_ES ( Efield_ind , efg_tmp , mu_ind , u_tmp , vir_tmp , phi_tmp ) 
    fx = 0.0_dp ; fy = 0.0_dp ; fz = 0.0_dp

    Efield = Efield_stat + Efield_ind

    ! ===================
    !  stopping criteria
    ! ===================
    diff_efield = ABS ( Efield(1,1) - Efield_old )
    Efield_old = Efield(1,1)

    io_node WRITE ( stdout ,'(a,i4,a,3f18.10,2(a,f18.10))') &
    'scf = ',iscf,' Efield atom 1= ',Efield(1,1),Efield(1,2),Efield(1,3),' conv = ',diff_efield,' u_pol = ',u_pol

  enddo ! end of SCF loop

  ! ===========================
  !  charge info is recovered
  ! ===========================
  qch = qch_tmp
  qia = qia_tmp


  if ( ionode ) then
    blankline(stdout)
    WRITE ( stdout , '(a,i6,a)') 'scf calculation of the induced electric moment converged in ',iscf, ' iterations '
    WRITE ( stdout , '(a,e10.3)') 'Electric field is converged within ',conv_tol_ind
    blankline(stdout)
    WRITE ( stdout , '(a)' )     'Induced dipoles at atoms : '
    do ia = 1 , natm
      WRITE ( stdout , '(i5,a3,a,3f18.10)' ) &
      ia,atype(ia),' mu_ind = ', mu_ind ( ia , 1 ) , mu_ind ( ia , 2 ) , mu_ind ( ia , 3 )
    enddo
    blankline(stdout)
  endif


  return

END SUBROUTINE moment_from_pola

! *********************** SUBROUTINE moment_from_WFc ***************************
!
!> \brief
!!  This routines evaluates the dipole moment induced by Wannier centers (Wfc).
!!  the listing of Wfc is done as for the verlet list
!
!> \description  
!!           wfc_point         : array of size natm+1 
!!                               gives the starting and finishing index of array
!!                               list for a given atom i
!!                               jbegin = point(i)  jend = point(i+1) - 1
!!           wfc_list          : index list of neighboring atoms
!!
!!  how to use it :
!!                          do ia = 1, natm
!!                            jbegin = wfc_point(i)
!!                            jend = wfc_point(i+1) - 1
!!                            do jvnl = jbegin , jend
!!                              ja = wfc_list ( jvnl ) 
!!                              then ia en ja are neighboors   
!!                            enddo
!!                          enddo
!!
!> \param[out] mu electric dipole 
!
!> \author
!! FMV
!
!> \date 
!! January 2013
!
! ******************************************************************************
SUBROUTINE moment_from_WFc ( mu )

  USE constants,                ONLY :  Debye_unit
  USE config,                   ONLY :  natm , ntype , itype , atype , simu_cell , rx , ry , rz 
  USE cell,                     ONLY :  kardir , dirkar 
  USE io,                  ONLY :  ionode , kunit_DIPWFC , stdout 

  implicit none

  ! global
  real(kind=dp), intent ( out )  :: mu ( natm ,  3 ) 

  ! local
  integer :: it, jt , ia , ja , icount , k , jb , je , jwfc , ia_last
  integer, dimension ( : ) , allocatable :: cwfc 
  integer, dimension ( : ) , allocatable :: wfc_list, wfc_point         ! wfc neighbour list info
  real(kind=dp) :: rijsq, cutsq
  real(kind=dp) :: rxi , ryi , rzi
  real(kind=dp) :: rxij , ryij , rzij
  real(kind=dp) :: sxij , syij , szij
  real(kind=dp) :: dmu
  logical          :: lwannier 

  ! ==========================================================================
  !  quick return 
  !  Is there any Wannier center ? if yes lwannier = .TRUE.
  !  redondant : cette condition est egalement evalué avant ( fast anyway )
  ! ==========================================================================
  lwannier = .false.
  do it = 1 , ntype
    if ( lwfc ( it ) .gt. 0 ) lwannier = .true.
  enddo
  if ( .not. lwannier ) then
    return
  endif

#ifdef debug_wfc
  WRITE ( stdout , '(a)' ) 'debug : in moment_from_WFs'
#endif

  cutsq = rcut_wfc * rcut_wfc
  mu = 0.0_dp

  allocate( wfc_list ( natm * 250 ) , wfc_point(  natm + 1 ) )
  allocate ( cwfc ( natm ) )
  cwfc = 0
  wfc_list      = 0
  wfc_point     = 0

  ! ======================================
  !         cartesian to direct 
  ! ======================================
  CALL kardir ( natm , rx , ry , rz , simu_cell%B )

  icount = 1
  do ia = 1 , natm 
    it = itype ( ia ) 
    if ( lwfc ( it ) .gt. 0 )  then
      rxi = rx ( ia )
      ryi = ry ( ia )
      rzi = rz ( ia )
      k = 0    
      do ja = 1 , natm
        jt = itype (ja) 
        if ( lwfc ( jt ) .lt. 0 )  then
          rxij  = rxi - rx ( ja )
          ryij  = ryi - ry ( ja )
          rzij  = rzi - rz ( ja )
          sxij  = rxij - nint ( rxij )
          syij  = ryij - nint ( ryij )
          szij  = rzij - nint ( rzij )
          rxij  = sxij * simu_cell%A(1,1) + syij * simu_cell%A(1,2) + szij * simu_cell%A(1,3)
          ryij  = sxij * simu_cell%A(2,1) + syij * simu_cell%A(2,2) + szij * simu_cell%A(2,3)
          rzij  = sxij * simu_cell%A(3,1) + syij * simu_cell%A(3,2) + szij * simu_cell%A(3,3)
          rijsq = rxij * rxij + ryij * ryij + rzij * rzij
          if ( rijsq .lt. rcut_wfc ) then
            cwfc ( ia ) = cwfc ( ia ) + 1
            icount = icount+1
            k = k + 1
            wfc_list( icount - 1 ) = ja
            ia_last = ia
#ifdef debug_wfc
          WRITE ( stdout ,'(a4,i4,a4,i4,2e17.8)') atype(ia),ia,atype(ja),ja,SQRT(rijsq),rcut_wfc
          WRITE ( stdout ,'(a,5e17.8)')          'distance ',rxij,ryij,rzij,SQRT(rijsq),rcut_wfc
#endif
          endif

        endif
 
      enddo
      wfc_point( ia ) = icount-k
    endif
  enddo
  wfc_point ( ia_last + 1 ) = icount


  ! ===========================================
  !  check if some wannier centers is missing
  ! ===========================================
  do ia = 1 , natm
    it = itype ( ia ) 
    if ( ( lwfc(it) .gt. 0 ) .and. ( cwfc ( ia ) .ne. lwfc ( it ) ) ) then
      WRITE ( * ,* ) 'ERROR in moment_from_WFc : wannier center is missing for atom ',ia , cwfc ( ia ) , lwfc ( it ) 
      STOP 
    endif
  enddo

  ! =============================================
  !  compute dipole moments from Wannier centers
  ! =============================================
  do ia = 1 , natm
    it = itype ( ia ) 
    jb = wfc_point ( ia )
    je = wfc_point ( ia + 1 ) - 1
    if ( lwfc ( it ) .gt. 0 ) then
      rxi = rx ( ia )
      ryi = ry ( ia )
      rzi = rz ( ia )
      jb = wfc_point ( ia )
      je = wfc_point ( ia + 1 ) - 1
      do jwfc = jb , je
        ja = wfc_list ( jwfc )
        rxij  = rx ( ja ) - rxi
        ryij  = ry ( ja ) - ryi
        rzij  = rz ( ja ) - rzi
        sxij  = rxij - nint ( rxij )
        syij  = ryij - nint ( ryij )
        szij  = rzij - nint ( rzij )
        rxij  = sxij * simu_cell%A(1,1) + syij * simu_cell%A(1,2) + szij * simu_cell%A(1,3)
        ryij  = sxij * simu_cell%A(2,1) + syij * simu_cell%A(2,2) + szij * simu_cell%A(2,3)
        rzij  = sxij * simu_cell%A(3,1) + syij * simu_cell%A(3,2) + szij * simu_cell%A(3,3)
        mu ( ia , 1 ) = mu ( ia , 1 ) + rxij 
        mu ( ia , 2 ) = mu ( ia , 2 ) + ryij
        mu ( ia , 3 ) = mu ( ia , 3 ) + rzij
      enddo
    endif
  enddo

  ! charges of WFC
  mu = -2.0_dp * mu 

  if ( ionode .and. lwrite_dip_wfc ) then 
    OPEN ( UNIT = kunit_DIPWFC , FILE='DIPWFC' )     
    do ia= 1 , natm
      WRITE ( kunit_DIPWFC , '(a,3e16.8)' ) atype( ia ) , mu ( ia , 1 ) , mu ( ia , 2 ) , mu ( ia , 3 )  
    enddo
    CLOSE ( kunit_DIPWFC ) 
  endif

  ! see constants.f90
  ! Debye_unit = 2.54174618782479816355_dp / bohr

  if ( ionode ) then
    WRITE ( stdout , '(a)' ) 'Dipoles from WFCs'
    WRITE ( stdout , '(a)' ) 'atom         |mu| [D]      |mu| [eA]'
    do ia = 1 , natm 
      it = itype ( ia )  
      if ( lwfc ( it ) .eq. -1 ) cycle
      dmu = mu ( ia , 1 ) * mu ( ia , 1 ) + mu ( ia , 2 ) * mu ( ia , 2 ) + mu ( ia , 3 ) *  mu ( ia , 3 )
      dmu = SQRT ( dmu ) 
      WRITE ( stdout , '(a,4x,2e16.8)' ) atype(ia), dmu * Debye_unit, dmu 
    enddo 
  endif

  deallocate ( wfc_list , wfc_point )
  deallocate ( cwfc ) 

  ! ======================================
  !         direct to cartesian
  ! ======================================
  CALL dirkar ( natm , rx , ry , rz , simu_cell%A )

  return

END SUBROUTINE moment_from_WFc

END MODULE field 
! ===== fmV =====
