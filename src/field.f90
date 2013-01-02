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
#define debug
#define debug3
! ======= Hardware =======

MODULE field 

  USE config , ONLY : ntypemax , npolmax
  USE kspace 
  USE rspace

  implicit none

  character*60 :: ctrunc
  character*60 :: ctrunc_allowed(2)
  data ctrunc_allowed / 'linear' , 'quadratic' /    ! see initialize_param_bmlj
  integer   :: trunc

  logical, SAVE :: lKA                              ! Kob-Andersen model for BMLJ                        


  !              eps    /    / sigma*\ q         / sigma*\ p  \
  !     V  =   ------- |  p | ------- |   -  q  | ------- |    |      sigma* = 2^(1/6)*sigma
  !             q - p   \    \   r   /           \   r   /    /

  double precision :: qlj     ( ntypemax , ntypemax )
  double precision :: plj     ( ntypemax , ntypemax )
  double precision :: epslj   ( ntypemax , ntypemax )
  double precision :: sigmalj ( ntypemax , ntypemax )

  double precision :: mass  ( ntypemax )    ! masses ( not yet )
  double precision :: qch   ( ntypemax )    ! charges only for efg (no electrostatic interaction and forces)
  double precision :: dip   ( ntypemax , 3 )! dipole for efg (no electrostatic interaction and forces)
  double precision :: pol   ( npolmax  , 3 , 3 )! polarizability ( is defined at particular position )

  double precision :: utail

  !todo : make the following quantities dependent on ntypemax ( hardware parameter) 
  ! or find a way to read them directly from control.F 
  double precision :: rcutsq ( ntypemax , ntypemax )  
  double precision :: sigsq  ( ntypemax , ntypemax ) 
  double precision :: epsp   ( ntypemax , ntypemax )
  double precision :: fc     ( ntypemax , ntypemax )
  double precision :: uc     ( ntypemax , ntypemax )
  double precision :: uc1    ( ntypemax , ntypemax )
  double precision :: uc2    ( ntypemax , ntypemax )
  double precision :: testtab ( ntypemax , ntypemax )

  ! about the charge:
  ! the array  q (natm) should be defined in the config module even if it is only
  ! called outside ( coulomb + EFG ) -----> one solution would be to read q(natm)
  ! from the POSFF file 
  ! alphaES should be read in control.F
  !
  ! one should be careful where the force are set to zero. 
  ! fx , fy , fz are the total force (coulumb + LJ) 

  ! problem defined twice ( is it a problem ?)
  double precision                                :: alphaES            ! Ewald sum parameter 
  integer                                         :: ncellewald
  integer                                         :: ncelldirect
  integer                                         :: lpolar(ntypemax)

  TYPE ( kmesh ) :: km_coul
  TYPE ( kmesh ) :: km_coul_dip
  TYPE ( rmesh ) :: rm_coul

  double precision, dimension ( : , : ) , allocatable :: ef_t
  double precision, dimension ( : , : , : ) , allocatable :: efg_t

CONTAINS

!*********************** SUBROUTINE field_default_tag *************************
!
! set default values to field tags
!
!******************************************************************************

SUBROUTINE field_default_tag

  implicit none

  ! =================
  !  default values
  ! =================

  ! BMLJ
  lKA           = .false.
  ctrunc        = 'linear' 
  epslj         = 1.0d0
  sigmalj       = 1.0d0
  qlj           = 12.0d0
  plj           = 6.0d0
  mass          =  1.0d0
  ! Coulomb
  ncelldirect   =  2
  ncellewald    = 10
  alphaES       =  1.0d0
  qch           =  0.0d0
  dip           =  0.0d0
  pol           =  0.0d0
  lpolar        =  0

  return

END SUBROUTINE field_default_tag


!*********************** SUBROUTINE field_check_tag ***************************
!
! check field tag values
!
!******************************************************************************

SUBROUTINE field_check_tag

  USE config,   ONLY :  xn , ntype
  USE io_file,  ONLY :  ionode , stdout

  implicit none

  integer :: i
  logical :: allowed

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

  if ( ctrunc .eq. 'linear' ) then
    trunc = 1
  endif

  if ( ctrunc .eq. 'quadratic' ) then
    trunc = 2
  endif

  ! =====
  !  lKA
  ! =====
  if ( lKA .and. ntype .ne. 2 ) then
    if ( ionode ) WRITE ( stdout , '(a)' ) 'ERROR fieldtag lKA should be used with 2 differents types'
    STOP 
  endif
  if ( lKA ) then
    sigmalj ( 1 , 1 ) = 1.0d0
    sigmalj ( 2 , 2 ) = 0.88d0
    sigmalj ( 1 , 2 ) = 0.8d0
    sigmalj ( 2 , 1 ) = 0.8d0
    epslj   ( 1 , 1 ) = 1.0d0
    epslj   ( 2 , 2 ) = 0.5d0
    epslj   ( 1 , 2 ) = 1.5d0
    epslj   ( 2 , 1 ) = 1.5d0
    xn    ( 1 )       = 0.8d0
    xn    ( 2 )       = 0.2d0
  endif

  return 

END SUBROUTINE field_check_tag

!*********************** SUBROUTINE field_init ********************************
!
! force field initialisation
!
!******************************************************************************

SUBROUTINE field_init

  USE control,  ONLY :  calc , lbmlj , lcoulomb
  USE io_file,  ONLY :  ionode, stdin, stdout , kunit_OUTFF

  implicit none

  ! local
  character * 132 :: filename
  integer :: ioerr

  namelist /fieldtag/    lKA           , &       
                         ctrunc        , &
                         ncelldirect   , &
                         ncellewald    , &
                         alphaES       , &
                         qlj           , & 
                         plj           , & 
                         sigmalj       , &
                         epslj         , &
                         mass          , &
                         qch           , &
                         dip           , &
                         pol           , &  
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
      if ( ionode ) WRITE ( stdout, '(a)') 'ERROR reading input_file : fieldtag section is absent'
      STOP
    elseif ( ioerr .gt. 0 )  then
      if ( ionode ) WRITE ( stdout, '(a)') 'ERROR reading input_file : fieldtag wrong tag'
      STOP
    endif
  CLOSE ( stdin )

  ! ================================
  ! check field tags values
  ! ================================
  CALL field_check_tag

  ! ================================
  !  pint field info
  ! ================================
  CALL field_print_info(stdout)
  CALL field_print_info(kunit_OUTFF)
  
  if ( calc .eq. 'efg' .or. calc .eq. 'NL_field' ) return

  ! ================================
  ! initialize constant parameters
  ! ================================
  if ( lbmlj )    then
    CALL initialize_param_bmlj
    if ( ionode ) WRITE ( stdout ,'(a)' ) 'bmlj quantities initialized' 
  endif

  if ( lcoulomb ) then
    CALL initialize_coulomb
    if ( ionode ) WRITE ( stdout ,'(a)') 'coulombic quantities initialized' 
  endif

  return

END SUBROUTINE field_init

!*********************** SUBROUTINE field_print_info **************************
!
! print force field informationto standard output
!
!******************************************************************************

SUBROUTINE field_print_info(kunit)

  USE config,   ONLY :  ntype , box, atypei
  USE control,  ONLY :  calc , cutoff , lbmlj , lcoulomb , longrange
  USE io_file,  ONLY :  ionode 

  implicit none

  !local
  integer :: kunit, it , it1 , it2
  double precision :: aaa
  logical :: linduced

  if ( ionode ) then
    WRITE ( kunit ,'(a)')       '=============================================================' 
    WRITE ( kunit ,'(a)')       ''
    if ( .not. lcoulomb .and. calc .eq. 'efg' )  &
    WRITE ( kunit ,'(a)')       'point charges: (no forces yet !)'
    do it = 1 , ntype 
      WRITE ( kunit ,'(a,a,a,f10.5)') 'q',atypei(it),'      = ',qch(it) 
    enddo
    WRITE ( kunit ,'(a)')       ''
    WRITE ( kunit ,'(a)')       'static dipoles: (no forces yet !)'
    do it = 1 , ntype
      WRITE ( kunit ,'(a,a,a,3f10.5)') 'mu',atypei(it),'      = ',dip(it,1),dip(it,2),dip(it,3)
    enddo
    WRITE ( kunit ,'(a)')       ''
    if ( lcoulomb )    then 
      WRITE ( kunit ,'(a)')       'force field information :                      '
      WRITE ( kunit ,'(a)')       'discret charges:'
      WRITE ( kunit ,'(a)')       ''
      WRITE ( kunit ,'(a)')       '        qi qj   '
      WRITE ( kunit ,'(a)')       ' Vij = -------  '
      WRITE ( kunit ,'(a)')       '         rij    '          
      WRITE ( kunit ,'(a)')       ''
      if ( longrange .eq. 'direct' )  then
        WRITE ( kunit ,'(a)')       'direct summation'
        WRITE ( kunit ,'(a)')       'cubic cutoff in real space'
        WRITE ( kunit ,'(a,i10)')   '-ncelldirect ... ncelldirect     = ',ncelldirect
        WRITE ( kunit ,'(a,i10)')   'total number of cells            = ',( 2 * ncelldirect + 1 ) ** 3
      endif     
      if ( longrange .eq. 'ewald' )  then
        CALL estimate_alpha( aaa )
        WRITE ( kunit ,'(a)')       'ewald summation'
        WRITE ( kunit ,'(a,f10.5)') 'alpha                            = ',alphaES
        WRITE ( kunit ,'(a,i10)')   'ncellewald                       = ',ncellewald
        WRITE ( kunit ,'(a,i10)')   '' 
        WRITE ( kunit ,'(a,f10.5)') 'Note:this should hold alpha^2 * box^2 >> 1',alphaES*alphaES*box*box
        WRITE ( kunit ,'(a,f10.5)') 'from estimate_alpha alpha should be ', aaa
      endif
    endif

    ! RETURN ?
    if ( calc .ne. 'efg' .and. calc .ne. 'NL_field' ) then

      if ( lbmlj )       then     
        WRITE ( kunit ,'(a)')       'force field information :                      '
        WRITE ( kunit ,'(a)')       'no masses are implemented                      '
        WRITE ( kunit ,'(a)')       'LENNARD-JONES                  '
        if ( ntype .eq. 2 )  &
        WRITE ( kunit ,'(a)')       'BINARY MIXTURE LENNARD-JONES           '
        WRITE ( kunit ,'(a)')       '' 
        WRITE ( kunit ,'(a)')       '       eps    /    / sigma* \ q       / sigma*  \ p \'
        WRITE ( kunit ,'(a)')       ' V = ------- |  p | ------- |    - q | -------- |   |'   
        WRITE ( kunit ,'(a)')       '      q - p   \    \   r    /         \    r    /   /'
        WRITE ( kunit ,'(a)')       ''
        if ( ntype .eq. 2 .and. lKA )  &
        WRITE ( kunit ,'(a)')       'KOB-ANDERSEN MODEL --- PhysRevE 51-4626 (1995) '
        if ( .not. lKA ) &
        WRITE ( kunit ,'(a)')       'USER DEFINED MODEL                             '
        WRITE ( kunit ,'(a,f10.5)') 'cutoff      = ',cutoff
        WRITE ( kunit ,'(a,a)')     'truncation  = ',ctrunc
        WRITE ( kunit ,'(a)')       ''

        do it1 = 1 , ntype
          do it2 = it1 , ntype          
            WRITE ( kunit ,'(a)')       '--------------------------------------------------------' 
            WRITE ( kunit ,'(a1,a1,a1,a)')       atypei(it1),'-',atypei(it2),' interactions:'    
            WRITE ( kunit ,'(a)')       '--------------------------------------------------------' 
            WRITE ( kunit ,'(a,f10.5)') 'sigma                                = ',sigmalj ( it1 , it2 )
            WRITE ( kunit ,'(a,f10.5)') 'eps                                  = ',epslj   ( it1 , it2 )
            WRITE ( kunit ,'(a,f10.5)') 'q                                    = ',qlj     ( it1 , it2 )
            WRITE ( kunit ,'(a,f10.5)') 'p                                    = ',plj     ( it1 , it2 )
          enddo
        enddo
      endif
    endif ! efg

    WRITE ( kunit ,'(a)') '' 

    linduced = .false.
    do it = 1 , ntype
      if ( lpolar(it) .eq. 1 )  linduced = .true.
    enddo
    if ( linduced ) then 
      WRITE ( kunit ,'(a)')          '--------------------------------------------------------'
      WRITE ( kunit ,'(a)')   ' '
      WRITE ( kunit ,'(a)')          'polarizabilities on atoms'
      WRITE ( kunit ,'(a)')   ' '
      WRITE ( kunit ,'(a)')          '--------------------------------------------------------'
      do it = 1 , ntype
        if ( lpolar(it) .eq. 1 ) then
          WRITE ( kunit ,'(a,a1,a)') 'polarizabiliy tensor of type ', atypei(it),' : ' 
          WRITE ( kunit ,'(a)')   ' '
          WRITE ( kunit ,'(3f12.5)')  pol( it , 1 , 1 ) , pol( it , 1 , 2 ) , pol( it , 1 , 3 ) 
          WRITE ( kunit ,'(3f12.5)')  pol( it , 2 , 1 ) , pol( it , 2 , 2 ) , pol( it , 2 , 3 ) 
          WRITE ( kunit ,'(3f12.5)')  pol( it , 3 , 1 ) , pol( it , 3 , 2 ) , pol( it , 3 , 3 ) 
          WRITE ( kunit ,'(a)')   ' '
        else
          WRITE ( kunit ,'(a,a1)') 'no polarizabiliy on type ', atypei(it)
        endif
        WRITE ( kunit ,'(a)')   ' '
      enddo
    endif

  endif

  return 

END SUBROUTINE field_print_info

!*********************** SUBROUTINE initialize_param_bmlj *********************
!
! iniitialisation of BMLJ parameters
!
!******************************************************************************

SUBROUTINE initialize_param_bmlj 

  USE constants,        ONLY :  pi
  USE config,           ONLY :  ntypemax , natm , box , omega , natmi , rho , atype , itype  , ntype
  USE control,          ONLY :  skindiff , cutoff , lreduced, calc
  USE io_file,          ONLY :  ionode, stdout, kunit_OUTFF

  implicit none

  ! local
  integer :: it,jt
  double precision :: one13, one16, two16, rskinmax, utail
  double precision :: sr2, sr, srp, srq
  double precision :: sig ( ntypemax , ntypemax )  
  double precision :: rcut3 ( ntypemax , ntypemax )
  double precision :: epsl ( ntypemax , ntypemax )
  double precision :: rskin ( ntypemax , ntypemax ) 
  double precision :: rskinsq ( ntypemax , ntypemax ) 
  double precision :: ut ( ntypemax , ntypemax ) 

  double precision :: rcut ( ntype , ntype ) 
  double precision :: ppqq ( ntype , ntype )
  double precision :: pp ( ntype , ntype )
  double precision :: qq ( ntype , ntype )
  double precision :: pp3 ( ntype , ntype ) 
  double precision :: qq3 ( ntype , ntype ) 

  do it = 1 , ntype
    do jt = 1 , ntype
      pp ( it , jt ) = plj ( it , jt ) 
      qq ( it , jt ) = qlj ( it , jt ) 
    enddo
  enddo
  pp3 = pp - 3.0d0
  qq3 = qq - 3.0d0

  one13 = (1.0D0 / 3.0D0)
  one16 = (1.0D0 / 6.0D0)
  two16 = 2.D0 **  one16

  ! =========================================
  !  calculate potential, force parameters   
  ! =========================================
  sig = sigmalj
  epsl = epslj 
  
  rskinmax = 0.0D0
  utail = 0.0d0
  
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
  !  c1 = ---------  |   p |  -------- |    - q |  -------- |   |      with rc = cutoff
  !         q - p     \     \   rc    /          \   rc    /    / 
  ! 
  ! ==================================================================================================
  

  ! ==================================================================================================
  !  truncation presented in J. Chem. Phys. 135 (2011) , Sengupta, Vasconcelos, Affouard, Sastry
  ! 
  !  ctrunc = 'quadratic'
  !  trunc  = 2    
  ! 
  !  Not sure !!!!!! 
  !
  !           eps    /    / sigma* \ q         / sigma* \ p  \
  !  V  =   ------- |  p | ------- |    -  q  | --------|    |   +  c1 r^2  -  c2    with  sigma* = 2^(1/6)*sigma
  !          q - p   \    \   r   /            \    r   /    /
  !
  !
  ! Maple (november 2011)
  !     
  !               eps p  q          /  / sigma* \ q       /  sigma* \ p \
  !    c1 =  --------------------- |  | --------|    -   | --------- |  |
  !          2  ( q - p ) * rc^2    \  \   rc   /         \   rc    /   /
  ! 
  !    and
  !     
  !          epsilon     /           / sigma* \ q               / sigma* \ p  \
  !    c2 = ----------- |  (2p+qp)  | --------|     - (2q+pq)  | --------|   |
  !         2 ( q - p )  \           \   rc   /                 \   rc   /   /       
  !
  ! ==================================================================================================

  do jt = 1 , ntype 
    do it = 1 , ntype

      ! intermediate terms 
          rcut   ( it , jt ) = cutoff                                      ! rc
          rcutsq ( it , jt ) = rcut ( it , jt ) * rcut   ( it , jt )       ! rc^2 
          rcut3  ( it , jt ) = rcut ( it , jt ) * rcutsq ( it , jt )       ! rc^3
          ppqq   ( it , jt ) = pp   ( it , jt ) * qq     ( it , jt )       ! p x q
          rskin  ( it , jt ) = rcut ( it , jt ) + skindiff

         if ( rskin ( it , jt ) .gt. rskinmax ) rskinmax = rskin ( it , jt )
           rskinsq ( it , jt ) = rskin ( it , jt ) * rskin ( it , jt )
           ! eps / q - p   
           epsp ( it , jt ) = epsl( it , jt )/( qq ( it , jt )-pp ( it , jt ) )
           !                                  sigma*^2 
           sigsq( it , jt ) = two16 * two16 * sig ( it , jt ) * sig ( it , jt )
           ! sigma*^2 / rc ^2
           sr2 = sigsq ( it , jt )/rcutsq ( it , jt )
           ! sigma* / rc 
           sr = SQRT (sr2)
           ! (sigma* / rc ) ^ p
           srp = sr ** pp ( it , jt )
           ! (sigma* / rc ) ^ q 
           srq = sr ** qq ( it , jt )
           ! trunc = 1
           uc ( it , jt )  = epsp ( it , jt ) * ( pp ( it , jt ) * srq - qq ( it , jt ) * srp )
           ! trunc = 2
           ! c1
           uc1 ( it , jt ) = epsp ( it , jt ) *  ppqq ( it , jt ) / ( 2.0d0 * rcutsq ( it , jt ) ) * ( srq - srp ) 
           ! c2
           uc2 ( it , jt ) = 0.5d0 * epsp( it , jt )  * (  &
                         ( 2.0d0 * qq ( it , jt ) + ppqq ( it , jt )  ) * srq - &
                         ( 2.0d0 * pp ( it , jt ) + ppqq ( it , jt )  ) * srp )
           ! for the virial
           fc ( it , jt ) =  ppqq ( it , jt ) * epsp ( it , jt ) /  sigsq ( it , jt ) 
           ! tail energy
           ut ( it , jt ) = epsp ( it , jt ) * ( pp ( it , jt ) * srq / qq3 ( it , jt ) - &
                                                qq ( it , jt ) * srp / pp3 ( it , jt )  )       
           ut ( it , jt ) = ut ( it , jt ) * rcut3 ( it , jt ) * 2.0d0 * pi 
           if ( ( natmi ( it ) .ne. 0 ) .and. ( natmi ( jt ) .ne. 0 ) ) &
           utail = utail + ut ( it , jt ) * natmi ( it ) * natmi ( jt ) / omega
    enddo
  enddo
 
  if ( calc .ne. 'md' ) return
  if (ionode.and..not.lreduced) WRITE( stdout     ,'(a,2f20.9)') 'long range correction : ',utail
  if (ionode.and.lreduced)      WRITE( stdout     ,'(a,2f20.9)') 'long range correction : ',utail/natm
  if (ionode.and..not.lreduced) WRITE( kunit_OUTFF,'(a,2f20.9)') 'long range correction : ',utail
  if (ionode.and.lreduced)      WRITE( kunit_OUTFF,'(a,2f20.9)') 'long range correction : ',utail/natm
  if (ionode) WRITE( stdout     ,'(a)') ' '
  if (ionode) WRITE( kunit_OUTFF,'(a)') ' '

  return

END SUBROUTINE initialize_param_bmlj


!*********************** SUBROUTINE engforce_driver ***************************
!
! this subroutine is used as an interfaced to the different potential + forces
! subroutines
!
!******************************************************************************

SUBROUTINE engforce_driver ( iastart , iaend )!, list , point )

  USE config,   ONLY :  natm
  USE control,  ONLY :  lpbc , lminimg , lbmlj , lcoulomb , lshiftpot , longrange
  USE io_file,  ONLY :  stdout 

  implicit none

  ! global
  integer, intent(inout)  :: iastart , iaend 
  ! local 
  logical :: lj_and_coul

  ! ====================================
  !  check if LJ and COULOMBIC 
  !  interactions are needed together
  ! ====================================
  lj_and_coul = .false.
  if ( lbmlj .and. lcoulomb ) then
    lj_and_coul = .true.
  endif

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
      if ( lshiftpot .and. lminimg ) CALL engforce_bmlj_pbc         ( iastart , iaend )
      if ( .not.lshiftpot )          CALL engforce_bmlj_pbc_noshift ( iastart , iaend )
  !    if ( longrange .eq. 'ewald'  ) CALL engforce_charge_ES       ( iastart , iaend , tmp ) 
  !    if ( longrange .eq. 'direct' ) CALL engforce_charge_DS       ( iastart , iaend ) 
      ! =======
      !  NO PBC
      ! =======
    else
      WRITE ( stdout ,'(a)') 'ERROR sick job: coulomb + nopbc'
      STOP 
      if ( lshiftpot )               CALL engforce_bmlj_nopbc       ( iastart , iaend )
!      if ( .not.lshiftpot )          CALL engforce_bmlj_pbc_noshift ( iastart , iaend ) 
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
       if ( lshiftpot       )     CALL engforce_bmlj_pbc         ( iastart , iaend )
       if ( .not. lshiftpot )     CALL engforce_bmlj_pbc_noshift ( iastart , iaend )
       ! =======
       !  NO PBC
       ! =======
     else
#ifdef debug  
       WRITE ( stdout , '(a)' ) 'debug : no pbc ' 
#endif  
       if ( lshiftpot )          CALL engforce_bmlj_nopbc          ( iastart , iaend )
!      if ( .not. lshiftpot )   CALL engforce_bmlj_nopbc_noshift ( iastart , iaend ) 
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
 !     if ( longrange .eq. 'ewald'  )  CALL engforce_charge_ES ( )
 !     if ( longrange .eq. 'direct' )  CALL engforce_charge_DS ( iastart , iaend ) 
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



!*********************** SUBROUTINE engforce_bmlj_pbc *************************
!
! total potential energy forces for each atoms for a bmlj potential with 
! periodic boundaries conditions, with or without vnlist (lvnlist=.TRUE.OR.FALSE.)
!
!******************************************************************************

SUBROUTINE engforce_bmlj_pbc ( iastart , iaend )

  USE config,           ONLY : natm , box , rx , ry , rz , fx , fy , fz, atype , itype , list , point , ntype
  USE control,          ONLY : lvnlist , myrank
  USE thermodynamic,    ONLY : u_lj , vir_lj
  USE time
  USE io_file,          ONLY : stdout

  implicit none
  INCLUDE 'mpif.h'

  ! global
  integer, intent(in)  :: iastart , iaend 

  ! local
  integer          :: ia , ja , it , jt , j1 , je , jb , ierr
  integer          :: p1 , p2
  integer          :: nxij , nyij , nzij
  double precision :: rxi , ryi , rzi
  double precision :: rxij , ryij , rzij , sr2 , rijsq , srp , srq
  double precision :: wij , fxij , fyij , fzij
  double precision :: ptwo ( ntype , ntype )
  double precision :: qtwo ( ntype , ntype )
  double precision :: one13 , two13 , forcetime1 , forcetime2 , invbox
  double precision :: u , vir 

#ifdef debug
  write( stdout , '(a,a)' )  'debug : atype',atype
  write( stdout , '(a,i4)' ) 'debug : itype',itype
#endif

  forcetime1 = MPI_WTIME(ierr) ! timing info

  u   = 0.0D0
  vir = 0.0D0
  fx  = 0.0D0
  fy  = 0.0D0
  fz  = 0.0D0
 
  invbox = 1.0d0 / box

  one13 = ( 1.0D0 / 3.0D0 )
  two13 = 2.D0 ** one13

  do jt = 1 , ntype
    do it = 1 , ntype
      ptwo ( it , jt ) = plj ( it , jt ) * 0.5d0
      qtwo ( it , jt ) = qlj ( it , jt ) * 0.5d0
    enddo
  enddo

  if ( lvnlist ) CALL vnlistcheck ( iastart , iaend )
  do ia = iastart , iaend
    rxi = rx ( ia )
    ryi = ry ( ia )
    rzi = rz ( ia )
    if ( lvnlist ) then
      jb = point( ia )
      je = point( ia + 1 ) - 1
    else
      jb = ia 
      je = iaend
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
        nxij = NINT ( rxij * invbox )
        nyij = NINT ( ryij * invbox )
        nzij = NINT ( rzij * invbox )
        rxij = rxij - box * nxij
        ryij = ryij - box * nyij
        rzij = rzij - box * nzij
        rijsq = rxij * rxij + ryij * ryij + rzij * rzij
        p1 = itype ( ja )
        p2 = itype ( ia )
        if ( rijsq .lt. rcutsq(p1,p2) ) then
          sr2 = sigsq(p1,p2) / rijsq
          srp = sr2 ** (ptwo(p1,p2))
          srq = sr2 ** (qtwo(p1,p2))
          if ( trunc .eq. 1 ) then
            u =  u + epsp(p1,p2) * ( plj(p1,p2) * srq -qlj(p1,p2) * srp ) - uc(p1,p2)
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
        endif
      endif
    enddo
  enddo
  vir = vir/3.0D0

  forcetime2 = MPI_WTIME(ierr) ! timing info
  forcetimetot = forcetimetot+(forcetime2-forcetime1)

  CALL MPI_ALL_REDUCE_DOUBLE_SCALAR ( u   ) 
  CALL MPI_ALL_REDUCE_DOUBLE_SCALAR ( vir ) 

  CALL MPI_ALL_REDUCE_DOUBLE ( fx , natm ) 
  CALL MPI_ALL_REDUCE_DOUBLE ( fy , natm ) 
  CALL MPI_ALL_REDUCE_DOUBLE ( fz , natm ) 

  u_lj = u
  vir_lj = vir

  return

END SUBROUTINE engforce_bmlj_pbc

!*********************** SUBROUTINE engforce_bmlj_nopbc ***********************
!
! total potential energy forces for each atoms for a bmlj potential with 
! *NO* periodic boundaries conditions, with or without vnlist (lvnlist=.TRUE.OR.FALSE.)
!
!******************************************************************************

SUBROUTINE engforce_bmlj_nopbc ( iastart , iaend )!, list , point )

  USE config,           ONLY : natm , box , rx , ry , rz , fx , fy , fz , itype , list , point , ntype 
  USE control,          ONLY : lvnlist , myrank
  USE thermodynamic,    ONLY : u_lj , vir_lj
  USE time

  implicit none
  INCLUDE 'mpif.h'

  ! global
  integer, intent(inout) :: iastart , iaend !, list(250 * natm),point(natm+1)

  ! local
  integer :: ia , ja , it, jt, j1, je, jb !, ierr
  integer :: p1, p2
  double precision :: rxi , ryi , rzi
  double precision :: rxij , ryij , rzij , sr2 , rijsq , srp , srq
  double precision :: wij , fxij , fyij , fzij
  double precision :: ptwo ( ntype , ntype )
  double precision :: qtwo ( ntype , ntype )
  double precision :: one13, two13
  double precision :: u, vir

  u = 0.0D0
  vir = 0.0D0
  fx = 0.0D0
  fy = 0.0D0
  fz = 0.0D0

  one13 = (1.0D0/3.0D0)
  two13 = 2.D0 ** one13

  do jt = 1 , ntype
    do it = 1 , ntype
      ptwo ( it , jt ) = plj ( it , jt ) * 0.5d0
      qtwo ( it , jt ) = qlj ( it , jt ) * 0.5d0
    enddo
  enddo

  if ( lvnlist ) CALL vnlistcheck ( iastart , iaend ) !, list , point )
  do ia = iastart, iaend
    rxi = rx ( ia )
    ryi = ry ( ia )
    rzi = rz ( ia )
    if ( lvnlist ) then
      jb = point ( ia )
      je = point ( ia + 1 ) - 1
    else
      jb = ia 
      je = iaend
    endif
    do j1 = jb, je
      if ( lvnlist ) then
        ja = list(j1)
      else
        ja = j1
      endif
      if ( ja .ne. ia ) then
        rxij = rxi - rx ( ja )
        ryij = ryi - ry ( ja )
        rzij = rzi - rz ( ja )
        rijsq = rxij * rxij + ryij * ryij + rzij * rzij
        p1 = itype ( ja )
        p2 = itype ( ia )
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
        endif
      endif
    enddo
  enddo
  vir = vir/3.0D0

  CALL MPI_ALL_REDUCE_DOUBLE_SCALAR ( u ) 
  CALL MPI_ALL_REDUCE_DOUBLE_SCALAR ( vir ) 
  
  CALL MPI_ALL_REDUCE_DOUBLE ( fx , natm )
  CALL MPI_ALL_REDUCE_DOUBLE ( fy , natm )
  CALL MPI_ALL_REDUCE_DOUBLE ( fz , natm )

  u_lj = u
  vir_lj = vir

  return

END SUBROUTINE engforce_bmlj_nopbc

!*********************** SUBROUTINE engforce_bmlj_pbc_noshift *****************
!
!FOR TEST PURPOSE
!
! same as engforce_bmlj_pbc but with no shift in the potential 
! if ok this should be merged !!
!
!FOR TEST PURPOSE
!
!******************************************************************************

SUBROUTINE engforce_bmlj_pbc_noshift ( iastart , iaend )!, list , point )

  USE config,           ONLY : natm , box , rx , ry , rz , fx , fy , fz , itype , list , point , ntype
  USE control,          ONLY : lvnlist , myrank
  USE thermodynamic,    ONLY : u_lj , vir_lj
  USE time

  implicit none
  INCLUDE 'mpif.h'

  ! global
  integer, intent(in)  :: iastart , iaend 
  !  integer, intent(out) :: list( 250 * natm ) , point( natm+1 )

  ! local
  integer :: ia , ja , it , jt , j1 , je , jb , ierr
  integer :: p1, p2
  integer :: nxij , nyij , nzij
  double precision :: rxi, ryi, rzi
  double precision :: rxij,ryij,rzij,sr2,rijsq,srp,srq
  double precision :: wij,fxij,fyij,fzij
  double precision :: ptwo ( ntype , ntype )
  double precision :: qtwo ( ntype , ntype )
  double precision :: one13, two13, forcetime1, forcetime2, invbox
  double precision :: u, vir
  
  forcetime1 = MPI_WTIME(ierr) ! timing info
  
  u   = 0.0D0
  vir = 0.0D0
  fx  = 0.0D0
  fy  = 0.0D0
  fz  = 0.0D0
  
  invbox = 1.0d0/box

  one13 = (1.0D0/3.0D0)
  two13 = 2.D0 ** one13

  do jt = 1 , ntype
    do it = 1 , ntype
      ptwo ( it , jt ) = plj ( it , jt ) * 0.5d0
      qtwo ( it , jt ) = qlj ( it , jt ) * 0.5d0
    enddo
  enddo

  if ( lvnlist ) CALL vnlistcheck ( iastart , iaend )
  do ia = iastart , iaend
    rxi = rx ( ia )
    ryi = ry ( ia )
    rzi = rz ( ia )
    if ( lvnlist ) then
      jb = point ( ia )
      je = point ( ia + 1 ) - 1
    else
      jb = ia 
      je = iaend 
    endif
    do j1 = jb, je
      if ( lvnlist ) then
        ja = list ( j1 )
      else
        ja = j1
      endif 
      if ( ja .ne. ia ) then
        rxij  = rxi - rx ( ja )
        ryij  = ryi - ry ( ja )
        rzij  = rzi - rz ( ja )
        nxij  = NINT ( rxij * invbox )
        nyij  = NINT ( ryij * invbox )
        nzij  = NINT ( rzij * invbox )
        rxij  = rxij - box * nxij
        ryij  = ryij - box * nyij
        rzij  = rzij - box * nzij
        rijsq = rxij * rxij + ryij * ryij + rzij * rzij
        p1    = itype ( ja )
        p2 = itype ( ia )
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
  vir = vir/3.0D0

  forcetime2 = MPI_WTIME(ierr) ! timing info
  forcetimetot = forcetimetot+(forcetime2-forcetime1)

  CALL MPI_ALL_REDUCE_DOUBLE_SCALAR ( u )
  CALL MPI_ALL_REDUCE_DOUBLE_SCALAR ( vir )

  CALL MPI_ALL_REDUCE_DOUBLE ( fx , natm )
  CALL MPI_ALL_REDUCE_DOUBLE ( fy , natm )
  CALL MPI_ALL_REDUCE_DOUBLE ( fz , natm )

  u_lj = u
  vir_lj = vir

  return

END SUBROUTINE engforce_bmlj_pbc_noshift

!*********************** SUBROUTINE engforce_bmlj_pbc_test *************************
!
! total potential energy forces for each atoms for a bmlj potential with 
! periodic boundaries conditions, with or without vnlist (lvnlist=.TRUE.OR.FALSE.)
!
! Note : for test only 
! here we only consider the potential/force due to a single particule at the
! origin 
!******************************************************************************

SUBROUTINE engforce_bmlj_pbc_test 

  USE config,           ONLY : natm , box , rx , ry , rz , fx , fy , fz , itype , ntype
  USE control,          ONLY : lvnlist , myrank
  USE thermodynamic,    ONLY : u_lj , vir_lj
  USE time

  implicit none
  INCLUDE 'mpif.h'

  ! local
  integer :: ja , it , jt , ierr
  integer :: p1, p2
  integer :: nxij , nyij , nzij
  double precision :: rxij,ryij,rzij,sr2,rijsq,srp,srq
  double precision :: wij,fxij,fyij,fzij
  double precision :: ptwo ( ntype , ntype )
  double precision :: qtwo ( ntype , ntype )
  double precision :: one13, two13, forcetime1, forcetime2, invbox
  double precision :: u , vir
  
  forcetime1 = MPI_WTIME(ierr) ! timing info
  
  u = 0.0D0
  vir = 0.0D0
  fx = 0.0D0
  fy = 0.0D0
  fz = 0.0D0

  invbox = 1.0d0/box

  one13 = (1.0D0/3.0D0)
  two13 = 2.D0 ** one13

  do jt = 1 , ntype
    do it = 1 , ntype
      ptwo ( it , jt ) = plj ( it , jt ) * 0.5d0
      qtwo ( it , jt ) = qlj ( it , jt ) * 0.5d0
    enddo
  enddo

  if ( lvnlist ) then
    print*,'no vnlist with this test purpose'
    STOP 
  else ! lvnlist .false.
    ja = 2
    rxij = rx ( ja )
    ryij = ry ( ja )
    rzij = rz ( ja )
    nxij = NINT ( rxij * invbox )
    nyij = NINT ( ryij * invbox )
    nzij = NINT ( rzij * invbox )
    rxij = rxij - box * nxij
    ryij = ryij - box * nyij
    rzij = rzij - box * nzij
    rijsq = rxij  * rxij + ryij * ryij + rzij * rzij
    p1 = 1
    p2 = 1 
    if ( rijsq .lt. rcutsq(p1,p2) ) then
      sr2 = sigsq(p1,p2)/rijsq
      srp = sr2 ** (ptwo(p1,p2))
      srq = sr2 ** (qtwo(p1,p2))
      u =  u + epsp(p1,p2) * ( plj(p1,p2) * srq -qlj(p1,p2)  * srp ) !- uc(p1,p2)
      wij = fc(p1,p2) * (srq-srp) * sr2
      vir = vir + wij * rijsq
      fxij = wij * rxij
      fyij = wij * ryij
      fzij = wij * rzij
      fx ( ja ) = fx ( ja ) + fxij
      fy ( ja ) = fy ( ja ) + fyij
      fz ( ja ) = fz ( ja ) + fzij
    endif
    vir = vir/3.0D0
  endif    

  forcetime2 = MPI_WTIME(ierr) ! timing info
  forcetimetot = forcetimetot+(forcetime2-forcetime1)

  CALL MPI_ALL_REDUCE_DOUBLE_SCALAR ( u )
  CALL MPI_ALL_REDUCE_DOUBLE_SCALAR ( vir )

  CALL MPI_ALL_REDUCE_DOUBLE ( fx , natm ) 
  CALL MPI_ALL_REDUCE_DOUBLE ( fy , natm )
  CALL MPI_ALL_REDUCE_DOUBLE ( fz , natm )

  u_lj = u
  vir_lj = vir

  return

END SUBROUTINE engforce_bmlj_pbc_test


!*********************** SUBROUTINE initialize_charge *************************
!
! this subroutine initialize the common quantities for charged particules.
! together with the k-space for the Ewald summation
! It is used by EFG and Coulombic subroutines 
!
! Problem : cannot be used together (i.e EFG + Coulomb ) or then converged the 
! Coulomb summation and used the same parameters for both quantities 
!
!******************************************************************************
SUBROUTINE initialize_coulomb

  USE config,   ONLY  : natm , natmi , ntype , qia , itype
  USE control,  ONLY  : longrange
  USE rspace,   ONLY  : direct_sum_init
  USE kspace,   ONLY  : kpoint_sum_init 

  implicit none

  ! local
  integer :: ia , it , nkcut , ncmax , ccs , cc

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
    km_coul%ncell = ncellewald
    nkcut = ( 2 * ncellewald + 1 ) ** 3 
    nkcut = nkcut - 1
    km_coul%nkcut = nkcut
    allocate( km_coul%kptk( nkcut ) , km_coul%kpt(3,nkcut) )
    allocate ( km_coul%strf ( nkcut, ntype ) ) 
    allocate ( km_coul%strf2 ( nkcut, ntype ) ) 
    CALL kpoint_sum_init ( km_coul )
    km_coul_dip%ncell = ncellewald
    km_coul_dip%meshlabel='km_coul_dip'
    !      nkcut = ( ncellewald + 1 ) ** 3
    nkcut = ( ncellewald + 1 ) ** 3
    nkcut = nkcut - 1
    km_coul_dip%nkcut = nkcut
    allocate( km_coul_dip%kptk( nkcut ) , km_coul_dip%kpt(3,nkcut) )
    allocate ( km_coul_dip%strf ( nkcut, ntype ) )
    allocate ( km_coul_dip%strf2 ( nkcut, ntype ) )
    CALL kpoint_sum_init_half ( km_coul_dip )

  endif

  return

END SUBROUTINE initialize_coulomb


SUBROUTINE finalize_coulomb

  USE control,  ONLY :  longrange , lcoulomb , calc

  implicit none

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
    deallocate ( km_coul%strf2 )
    deallocate ( km_coul_dip%strf )
    deallocate ( km_coul_dip%strf2 )
  endif


END SUBROUTINE finalize_coulomb

!*********************** SUBROUTINE engforce_coulomb_ES ***********************
!
! this subroutine calculates the potential energy and the force acting on each
! charge due to the coulombic interaction. The Ewald method is used
! 
! Note :  same consruction as efg_ES 
!
!******************************************************************************
SUBROUTINE engforce_charge_ES  ( iastart, iaend , Efield )

  USE control,          ONLY :  myrank , numprocs
  USE config,           ONLY :  system , natm , natmi , atype , atypei , box , omega , &
                        itype , rx , ry , rz , fx , fy , fz , ntype , qia , phi_coul_qq
  USE io_file,          ONLY :  ionode , stdout 
  USE constants,        ONLY :  pi , fpi , piroot , imag , tpi
  USE kspace,           ONLY :  struc_fact
  USE thermodynamic,    ONLY :  u_coul_qq , vir_coul_qq
  USE time

  implicit none

  INCLUDE 'mpif.h'

  ! global
  integer, intent(in)  :: iastart , iaend
  double precision, optional  :: Efield ( natm , 3 )

  ! local
  integer :: ia, ja , it , ierr , ip
  integer :: ik
  integer :: nxij , nyij , nzij
  double precision :: rij , rijsq , expon , invbox
  double precision :: alpha2
  double precision :: qi , qj , qij, qijf 
  double precision :: rxi , ryi , rzi , rxj , ryj , rzj , rxij , ryij , rzij 
  double precision :: wij0 , wij , fxij , fyij , fzij 
  double precision :: ak, kx, ky, kz, kk
  double precision :: kri 
  double precision :: pot , pot3 
  double precision :: u_dir , u_rec , u_surf , u_self , selfa
  double precision :: vir_dir , vir_rec , vir_surf
  double complex   :: rhon , carg 
  double precision   :: str
  double precision, external :: errfc 
  double precision :: ttt1 , ttt2 , ttt3 
  logical :: lcharge
  double precision  :: Efield_dir ( natm , 3 )
  double precision  :: Efield_rec ( natm , 3 )
  double precision  :: Efield_surf ( natm , 3 )
  double precision :: qsum ( 3 ) , qsq
  double precision, dimension(:), allocatable :: fx_dir, fy_dir, fz_dir
  double precision, dimension(:), allocatable :: fx_rec, fy_rec, fz_rec
  double precision, dimension(:), allocatable :: fx_surf, fy_surf, fz_surf
  double precision, dimension(:), allocatable :: phi_dir , phi_rec , phi_surf , phi_self 

  lcharge = .false.
  do it = 1 , ntype
    if ( qch(it) .ne. 0.0d0 ) lcharge = .true.
  enddo

  if ( .not. lcharge ) WRITE ( stdout , '(a)') 'No point charges'
  if ( .not. lcharge ) return

#ifdef debug
  write( stdout , '(a)')    'debug: in engforce_charge_ES'
  write( stdout , '(a,i8,a,a)') 'debug : km_coul ',km_coul%nkcut,' ',km_coul%meshlabel
  do ia = 1 , natm
    write( stdout , '(a,f12.5)') 'debug : charge (atom) ',qia(ia)
  enddo
  do it = 1 , ntype
    write( stdout , '(a,f12.5)') 'debug : charge (type ) ',qch(it)
  enddo
  write( stdout , '(a,2i8)')     'debug : iastart iaend',iastart ,iaend
  write( stdout , '(a,f20.5)')   'debug : alphaES ', alphaES
  call print_config_sample(0,0)
#endif 
 
  ttt1 = MPI_WTIME(ierr)

  u_rec   = 0.0d0
  u_dir   = 0.0d0
  u_surf  = 0.0d0
  u_self  = 0.0d0
  vir_rec = 0.0d0
  vir_dir = 0.0d0
  phi_coul_qq = 0.0d0
  Efield  = 0.0d0
  Efield_dir  = 0.0d0
  Efield_rec  = 0.0d0
  Efield_surf  = 0.0d0

  allocate( fx_dir(natm), fy_dir(natm), fz_dir(natm) )
  allocate( fx_rec(natm), fy_rec(natm), fz_rec(natm) )
  allocate( fx_surf(natm), fy_surf(natm), fz_surf(natm) )
  allocate( phi_rec(natm) , phi_dir(natm) , phi_surf(natm) , phi_self ( natm )  )
  fx_dir  = 0.0D0
  fy_dir  = 0.0D0
  fz_dir  = 0.0D0
  fx_rec  = 0.0D0
  fy_rec  = 0.0D0
  fz_rec  = 0.0D0
  fx_surf  = 0.0D0
  fy_surf  = 0.0D0
  fz_surf = 0.0D0
  phi_rec = 0.0d0
  phi_dir = 0.0d0
  phi_surf = 0.0d0
  phi_self = 0.0d0

  ! ==================
  !  total charge sum
  ! ==================
  do ia = 1 , natm
    qsum ( 1 ) = qsum ( 1 ) + qia ( ia ) * rx ( ia ) 
    qsum ( 2 ) = qsum ( 2 ) + qia ( ia ) * ry ( ia ) 
    qsum ( 3 ) = qsum ( 3 ) + qia ( ia ) * rz ( ia ) 
    qsq = qsq + qia ( ia ) * qia ( ia ) 
  enddo


  ! =====================
  ! facteur de structure 
  ! =====================
  CALL struc_fact ( km_coul )

  ! =================
  !  some constants 
  ! =================
  invbox = 1.0d0 / box
  ! related ewald parameter
  alpha2 = alphaES * alphaES     
  selfa = alphaES / piroot 

  ! ==============================================
  !        direct space part
  ! ==============================================

  do ia = 1 , natm 
    rxi = rx(ia)
    ryi = ry(ia)
    rzi = rz(ia)
    qi = qia(ia)
    do ja = 1, natm

      rxj = rx(ja)
      ryj = ry(ja) 
      rzj = rz(ja) 
      qj = qia(ja)
      rxij = rxi - rxj
      ryij = ryi - ryj
      rzij = rzi - rzj
      nxij = NINT ( rxij * invbox )
      nyij = NINT ( ryij * invbox )
      nzij = NINT ( rzij * invbox )
      rxij = rxij - box * nxij
      ryij = ryij - box * nyij
      rzij = rzij - box * nzij

      rijsq = rxij * rxij + ryij * ryij + rzij * rzij
      if (ja .ne. ia ) then

        rij = SQRT ( rijsq )
 
        pot  = qj  / rij
        qij  = qi  * pot
        pot3 = pot / rijsq 
        qijf = qi  * pot3         ! qi * qj / r^3

        wij0  = errfc( alphaES * rij ) 
        expon = EXP ( - alpha2 * rijsq ) / piroot
        wij   = ( wij0 + 2.0d0 * rij * alphaES * expon )

        u_dir          = u_dir + qij * wij0
        phi_dir ( ia ) = phi_dir ( ia ) + pot * wij0

        Efield_dir ( ia , 1 ) = Efield_dir ( ia , 1 ) + rxij * wij * pot3
        Efield_dir ( ia , 2 ) = Efield_dir ( ia , 2 ) + ryij * wij * pot3
        Efield_dir ( ia , 3 ) = Efield_dir ( ia , 3 ) + rzij * wij * pot3
 
        fxij = wij * rxij * qijf 
        fyij = wij * ryij * qijf
        fzij = wij * rzij * qijf 

        fx_dir(ia) = fx_dir(ia) + fxij
        fy_dir(ia) = fy_dir(ia) + fyij
        fz_dir(ia) = fz_dir(ia) + fzij

        vir_dir =  vir_dir + ( fxij * rxij + fyij * ryij + fzij * rzij )

      endif

       Efield_surf ( ia , 1 ) = Efield_surf ( ia , 1 ) + qj * rxj
       Efield_surf ( ia , 2 ) = Efield_surf ( ia , 2 ) + qj * ryj
       Efield_surf ( ia , 3 ) = Efield_surf ( ia , 3 ) + qj * rzj

       fx_surf( ia ) = fx_surf( ia ) + qi * qj * rxj
       fy_surf( ia ) = fy_surf( ia ) + qi * qj * ryj
       fz_surf( ia ) = fz_surf( ia ) + qi * qj * rzj

       phi_surf ( ia ) = phi_surf ( ia ) + qj * rijsq 

    enddo
  enddo

  ttt2 = MPI_WTIME(ierr)
  fcoultimetot1 = fcoultimetot1 + ( ttt2 - ttt1 ) 

  ! ==============================================
  !            reciprocal space part
  ! ==============================================
  kpoint : do ik = 1, km_coul%nkcut
    ! =================
    !   k-space  
    ! =================
    kx = km_coul%kpt(1,ik)
    ky = km_coul%kpt(2,ik)
    kz = km_coul%kpt(3,ik)
    kk = km_coul%kptk(ik)
    ak = EXP ( - kk * 0.25d0 / alpha2 ) / kk


    if (km_coul%kptk(ik) .eq. 0 ) then 
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
      rhon = rhon + qch(it) * CONJG( km_coul%strf ( ik , it ) )
    enddo

    str = rhon * CONJG(rhon)

    u_rec   = u_rec   + DBLE ( ak * str )
    vir_rec = vir_rec + DBLE ( ak * str * 0.5d0 * ( 1.0d0 - ( kk * 0.5d0 / alpha2 ) )  )


    do ia = 1 , natm 
      rxi = rx(ia)
      ryi = ry(ia)
      rzi = rz(ia)
      qi  = qia ( ia ) 

      kri = ( kx * rxi + ky * ryi + kz * rzi )

      carg = EXP ( imag * kri )
      wij = DBLE ( rhon * carg * ak * imag )

      phi_rec ( ia ) = phi_rec ( ia ) + DBLE ( ak * carg * rhon ) 

      Efield_rec ( ia , 1 ) = Efield_rec ( ia , 1 ) - kx * wij
      Efield_rec ( ia , 2 ) = Efield_rec ( ia , 2 ) - ky * wij
      Efield_rec ( ia , 3 ) = Efield_rec ( ia , 3 ) - kz * wij

      fxij = wij * kx
      fyij = wij * ky
      fzij = wij * kz

      fx_rec( ia ) = fx_rec( ia ) - qi * fxij
      fy_rec( ia ) = fy_rec( ia ) - qi * fyij
      fz_rec( ia ) = fz_rec( ia ) - qi * fzij
     
    enddo


  enddo kpoint

  ttt3 = MPI_WTIME(ierr)
  fcoultimetot2 = fcoultimetot2 + ( ttt3 - ttt2 ) 

  ! ======================================================
  ! remark on the unit :
  ! 1/(4*pi*epislon_0) = 1 => epsilon_0 = 1/4pi
  ! ======================================================
  u_dir   =   u_dir   * 0.5d0
  vir_dir = - vir_dir * 0.5d0

  u_rec   =   u_rec   * tpi / omega 
  vir_rec = - vir_rec * fpi / omega
  phi_rec =   phi_rec * fpi / omega

  u_surf = qsum ( 1 ) * qsum ( 1 ) + qsum ( 2 ) * qsum ( 2 ) + qsum ( 3 ) * qsum ( 3 )
  u_surf   =   u_surf * tpi / 3.0d0 / omega
  vir_surf = - u_surf
  phi_surf = - phi_surf * tpi / omega / 3.0d0

  u_self  = - selfa * qsq 
  ! self
  do ia = 1 , natm
   phi_self ( ia ) = - 2.0d0 * selfa * qia ( ia ) 
  enddo

  u_coul_qq   = ( u_dir + u_rec + u_surf + u_self )
  vir_coul_qq = ( vir_dir + vir_rec + vir_surf )
  phi_coul_qq = ( phi_dir + phi_rec + phi_surf + phi_self )

#ifdef debug
  write( stdout , '(4(a,f16.8))' ) ,'vir_dir    = ', vir_dir    ,' vir_rec    = ', vir_rec    ,' vir_surf    = ',vir_surf   ,' vir_coul_qq = '   ,vir_coul_qq
  write( stdout , '(5(a,f16.8))' ) ,'u_dir      = ', u_dir      ,' u_rec      = ', u_rec      ,' u_surf      = ',u_surf     ,' u_self      = '   ,u_self     ,' u_coul_qq      = ',u_coul_qq
  write( stdout , '(5(a,f16.8))' ) ,'phi_dir(1) = ', phi_dir(1) ,' phi_rec(1) = ', phi_rec(1) ,' phi_surf(1) = ',phi_surf(1),' phi_self(1) = '   ,phi_self(1),' phi_coul_qq(1) = ',phi_coul_qq(1)
#endif

  Efield_rec  = Efield_rec  * fpi / omega
  Efield_surf = Efield_surf * fpi / omega / 3.0d0

  fx_rec = fx_rec * fpi / omega
  fy_rec = fy_rec * fpi / omega
  fz_rec = fz_rec * fpi / omega
  fx_surf = fx_surf * fpi / omega / 3.0d0
  fy_surf = fy_surf * fpi / omega / 3.0d0
  fz_surf = fz_surf * fpi / omega / 3.0d0

  fx = fx_rec + fx_dir - fx_surf
  fy = fy_rec + fy_dir - fy_surf 
  fz = fz_rec + fz_dir - fz_surf



#ifdef debug
  do ip=1,natm
    WRITE ( stdout ,'(a,i4,a,3f12.8)') 'atom = ',ip,' Efield_dir         = ', &
    Efield_dir      ( ip , 1 ) , Efield_dir      ( ip , 2 ) , Efield_dir ( ip , 3 )
  enddo
  do ip=1,natm
    WRITE ( stdout ,'(a,i4,a,3f12.8)') 'atom = ',ip,' Efield_rec         = ', &
    Efield_rec ( ip , 1 ) , Efield_rec ( ip , 2 ) , Efield_rec ( ip , 3 )
  enddo
  do ip=1,natm
    WRITE ( stdout ,'(a,i4,a,3f12.8)') 'atom = ',ip,' Efield_surf         = ', &
    Efield_surf      ( ip , 1 ) , Efield_surf      ( ip , 2 ) , Efield_surf ( ip , 3 )
  enddo
#endif
 
   Efield = Efield_dir + Efield_rec - Efield_surf

  deallocate( fx_dir , fy_dir , fz_dir  )
  deallocate( fx_rec , fy_rec , fz_rec  )
  deallocate( fx_surf, fy_surf, fz_surf )
  deallocate( phi_rec , phi_dir , phi_surf , phi_self )

  return


END SUBROUTINE engforce_charge_ES


!*********************** SUBROUTINE engforce_coulomb_DS ***********************
!
! this subroutine calculates the potential energy and the force acting on each
! charge due to the coulombic interaction. The Direct Summation is used
! 
! Note :  same consruction as efg_DS
!
!******************************************************************************

SUBROUTINE engforce_charge_DS (  iastart , iaend , Efield ) 

  USE control,          ONLY :  myrank , numprocs, calc , longrange , cutlongrange
  USE config,           ONLY :  system , natm , natmi , atype , atypei , box , omega , &
                        itype , rx , ry , rz , fx , fy , fz , ntype , qia , phi_coul_qq
  USE io_file,          ONLY :  ionode , stdout 
  USE constants,        ONLY :  pi , fpi , piroot , imag , tpi
  USE thermodynamic,    ONLY :  u_coul_qq , vir_coul_qq
  USE time

  implicit none

  INCLUDE 'mpif.h'

  ! global
  integer, intent(in)  :: iastart , iaend 
  double precision, optional  :: Efield ( natm , 3 )

  ! local
  integer :: ia, ja , ierr , ncell , it
  double precision :: d , d2 
  double precision :: qi , qj , qij, qijf 
  double precision :: rxi , ryi , rzi , rxij , ryij , rzij 
  double precision :: rxj , ryj , rzj 
  double precision :: fxij , fyij , fzij 
  double precision :: u_dir , vir_dir 
  double precision :: ttt1 , ttt2  , ttt3
  double precision :: cutsq
  double precision, dimension(:), allocatable :: fx_dir, fy_dir, fz_dir
  double precision ::  phi ( natm ) , pot , pot3
  logical          :: lcharge

  ttt1 = MPI_WTIME(ierr)

  lcharge = .false.
  do it = 1 , ntype
    if ( qch(it) .ne. 0.0d0 ) lcharge = .true.
  enddo

  if ( .not. lcharge ) WRITE ( stdout , '(a)') 'No point charges'
  if ( .not. lcharge ) return


  ! =============================== 
  !         some constants
  ! =============================== 
  cutsq = cutlongrange * cutlongrange


#ifdef debug2
  write( stdout , '(a)')    'debug: in engforce_charge_DS'
  write( stdout , '(a,i8,a)') 'debug : rm_coul ',rm_coul%ncmax,rm_coul%meshlabel
  do ia = 1 , natm
    write( stdout , '(a,f12.5)') 'debug : charge',qia(ia)
  enddo
  write( stdout , '(a,2i8)')     'debug : iastart iaend',iastart ,iaend
  write( stdout , '(a,f20.5)')   'debug : cutsq ',cutsq
  call print_config_sample(0,0)
#endif  


  Efield (:,:) = 0.0d0
  vir_dir    = 0.0d0
  u_dir      = 0.0d0 
  phi (:)    = 0.0d0

  allocate( fx_dir(natm), fy_dir(natm), fz_dir(natm) )
  fx_dir  = 0.0D0
  fy_dir  = 0.0D0
  fz_dir  = 0.0D0

  ! =========================================================
  !  MAIN LOOP calculate EFG(i) for each atom i parallelized
  ! =========================================================
  atom : do ia = iastart , iaend

    rxi = rx  ( ia )
    ryi = ry  ( ia )
    rzi = rz  ( ia )
    qi  = qia ( ia )

  ! ==================================================
  ! sum over neighboring cells (see direct_sum_init )
  ! ==================================================
      do ncell = 1 , rm_coul%ncmax
      ! ==============================
      !  ia and ja in different cells
      ! ==============================
        if ( rm_coul%lcell ( ncell ) .eq. 1 ) then

          do ja = 1 , natm

            rxj   = rx ( ja ) + rm_coul%boxxyz( 1 , ncell )
            ryj   = ry ( ja ) + rm_coul%boxxyz( 2 , ncell )
            rzj   = rz ( ja ) + rm_coul%boxxyz( 3 , ncell )

            rxij  = rxi - rxj
            ryij  = ryi - ryj
            rzij  = rzi - rzj

            d2    = rxij * rxij + ryij * ryij + rzij * rzij

            if ( d2 .lt. cutsq ) then

              qj    = qia ( ja )
              d     = SQRT ( d2 )
              pot   = qj / d
              pot3  = pot / d2
              qij   = qi * pot 
              qijf  = qi * pot3 

              ! electrostatic energy
              u_dir      = u_dir      + qij 
              phi ( ia ) = phi ( ia ) + pot 

              Efield ( ia , 1 ) = Efield ( ia , 1 ) + rxij * pot3
              Efield ( ia , 2 ) = Efield ( ia , 2 ) + ryij * pot3
              Efield ( ia , 3 ) = Efield ( ia , 3 ) + rzij * pot3  

              fxij = qijf * rxij
              fyij = qijf * ryij
              fzij = qijf * rzij

              fx_dir ( ia ) = fx_dir ( ia ) + fxij
              fy_dir ( ia ) = fy_dir ( ia ) + fyij
              fz_dir ( ia ) = fz_dir ( ia ) + fzij

              ! virial
              vir_dir =  vir_dir + ( fxij * rxij + fyij * ryij + fzij * rzij )


            endif ! d2.lt.cutsq

          enddo ! ja

        endif 
        ! =======================================
        !  ia and ja in the same cell (ia.ne.ja)
        ! =======================================
        if ( rm_coul%lcell(ncell) .eq. 0) then
          do ja = 1,natm
 
            if ( ja .ne. ia ) then

              rxj  = rx ( ja )
              ryj  = ry ( ja )
              rzj  = rz ( ja )

              rxij = rxi - rxj
              ryij = ryi - ryj
              rzij = rzi - rzj

              d2  = rxij * rxij + ryij * ryij + rzij * rzij

              if ( d2 .lt. cutsq ) then

                qj    = qia ( ja )
                d     = SQRT ( d2 )
                pot   = qj / d 
                pot3  = pot / d2 
                qij   = qi * pot
                qijf  = qi * pot3

                ! electrostatic energy
                u_dir = u_dir + qij
                phi ( ia ) = phi ( ia ) + pot

                Efield ( ia , 1 ) = Efield ( ia , 1 ) + rxij * pot3
                Efield ( ia , 2 ) = Efield ( ia , 2 ) + ryij * pot3
                Efield ( ia , 3 ) = Efield ( ia , 3 ) + rzij * pot3

                fxij = qijf * rxij
                fyij = qijf * ryij
                fzij = qijf * rzij

                fx_dir ( ia ) = fx_dir ( ia ) + fxij
                fy_dir ( ia ) = fy_dir ( ia ) + fyij
                fz_dir ( ia ) = fz_dir ( ia ) + fzij

                ! virial
                vir_dir =  vir_dir + ( fxij * rxij + fyij * ryij + fzij * rzij )

              endif

            endif ! ia.ne.ja

          enddo ! ja

        endif 

      enddo ! ncell

  enddo atom

  ttt2 = MPI_WTIME(ierr)

  CALL MPI_ALL_REDUCE_DOUBLE_SCALAR ( u_dir )
  CALL MPI_ALL_REDUCE_DOUBLE_SCALAR ( vir_dir ) 

  CALL MPI_ALL_REDUCE_DOUBLE ( fx_dir , natm )
  CALL MPI_ALL_REDUCE_DOUBLE ( fy_dir , natm )
  CALL MPI_ALL_REDUCE_DOUBLE ( fz_dir , natm )

  u_coul_qq   =   u_dir    * 0.5d0
  vir_coul_qq = - vir_dir  * 0.5d0
  phi_coul_qq = phi
  
  fx = fx + fx_dir
  fy = fy + fy_dir
  fz = fz + fz_dir

#ifdef debug
  write( stdout , '(5(A,f12.8))') 'U_qq  = ', u_coul_qq,' phi_qq,1   = ', phi(1),' E1,x = ',Efield(1,1),' f1,x = ', fx(1) , ' vir_qq = ',vir_coul_qq 
  write( stdout , '(5(A,f12.8))') 'U_dir = ', u_dir,    ' vir_dir    = ',vir_dir
#endif

  deallocate( fx_dir, fy_dir, fz_dir )

  ttt3 = MPI_WTIME(ierr)

  return

END SUBROUTINE engforce_charge_DS

!*********************** SUBROUTINE field_dipole_DS **********************
!
! this subroutine calculates the electric field at atoms position from all the 
! static dipoles of the system. The Direct Summation is used
! 
! Note : same consruction as efg_DS
!                          
!  Ei,alpha  =  sum_j!=i  T^alpha,beta_ij * mu_j,beta  
!
!  T^alpha,alpha_ij  = ( 3 * r_ij,alpha * r_ij,beta - r^2 delta_alpha,beta ) / r_ij^5 
!
!  
! when linduced=.TRUE., the mu in input is used 
! however when linduced=.FALSE., the mu has some value in output (static dipoles)
!
!******************************************************************************

SUBROUTINE engforce_dipole_DS ( iastart , iaend , Efield , mu )

  USE config,   ONLY : natm , ntype , rx , ry , rz , fx , fy, fz , dipia , phi_coul_dd
  USE control,  ONLY : cutlongrange
  USE thermodynamic,    ONLY :  u_coul_dd , vir_coul_dd
  USE io_file

  implicit none

  ! global 
  integer, intent(in) :: iastart , iaend 
  double precision :: Efield ( natm , 3) 
  double precision :: mu    ( natm , 3) 

  ! local
  integer :: ia , ja  , it , ncell
  double precision :: d , d2 , dm3 , dm5 , dm7
  double precision :: Txx , Tyy , Tzz , Txy , Txz , Tyz
  double precision :: rxi , ryi , rzi , rxij , ryij , rzij
  double precision :: rxj , ryj , rzj 
  double precision :: fxij , fyij , fzij 
  double precision ::  Txxx,  Tyyy,  Tzzz, Txxy, Txxz, Tyyx, Tyyz, Tzzx, Tzzy, Txyz
  double precision :: cutsq , u_dir , vir_dir
  double precision, dimension(:), allocatable :: fx_dir, fy_dir, fz_dir
  double precision ::  phi ( natm ) 
  logical :: ldipole

  ldipole = .false.
  do it = 1 , ntype 
    if ( mu(it,1) .ne. 0.0d0 .or. & 
    mu(it,2) .ne. 0.0d0 .or. &
    mu(it,3) .ne. 0.0d0 ) ldipole = .true.
  enddo
  if ( .not. ldipole .and. ionode ) WRITE ( stdout , '(a)') 'No static dipoles' 
  if ( .not. ldipole ) return

  cutsq = cutlongrange * cutlongrange

#ifdef debug2
  write ( stdout , '(a)') 'debug : in engforce_dipole_DS'
  write ( stdout , '(a,i)') 'debug : rm_coul',rm_coul%ncmax
  do ia = 1 , natm
  write( * , '(a,3f12.5)' ) 'debug : dipole ', mu (ia,1), mu (ia,2), mu (ia,3)
  enddo
#endif  

  Efield ( : , : ) = 0.0d0
  vir_dir = 0.0d0
  u_dir   = 0.0d0
  phi ( : ) = 0.0d0

  allocate( fx_dir(natm), fy_dir(natm), fz_dir(natm) )
  fx_dir  = 0.0D0
  fy_dir  = 0.0D0
  fz_dir  = 0.0D0


  atom : do ia = iastart , iaend
           rxi = rx(ia)
           ryi = ry(ia)
           rzi = rz(ia)
          ! ==================================================
          ! sum over neighboring cells (see direct_sum_init )
          ! ==================================================
          do ncell = 1 , rm_coul%ncmax
          ! ==============================
          !  ia and ja in different cells
          ! ==============================
          if ( rm_coul%lcell(ncell) .eq. 1) then
 
            do ja = 1 , natm
              rxj  = rx(ja) + rm_coul%boxxyz(1,ncell)
              ryj  = ry(ja) + rm_coul%boxxyz(2,ncell)
              rzj  = rz(ja) + rm_coul%boxxyz(3,ncell)
              rxij = rxi - rxj
              ryij = ryi - ryj
              rzij = rzi - rzj
              d2   = rxij * rxij + ryij * ryij + rzij * rzij

              if ( d2 .lt. cutsq ) then
                d  = SQRT (d2)
                dm3 = 1.0d0 / ( d2 * d )  
                dm5 = dm3 / d2
                dm7 = dm5 / d2
          
                Txx = 3.0d0 * rxij * rxij - d2 
                Tyy = 3.0d0 * ryij * ryij - d2 
                Tzz = 3.0d0 * rzij * rzij - d2 
                Txy = 3.0d0 * rxij * ryij  
                Txz = 3.0d0 * rxij * rzij  
                Tyz = 3.0d0 * ryij * rzij  
          
                u_dir = u_dir - ( mu ( ia , 1 ) * Txx * mu ( ja , 1 ) + &
                                  mu ( ia , 1 ) * Txy * mu ( ja , 2 ) + &
                                  mu ( ia , 1 ) * Txz * mu ( ja , 3 ) + &
                                  mu ( ia , 2 ) * Txy * mu ( ja , 1 ) + &
                                  mu ( ia , 2 ) * Tyy * mu ( ja , 2 ) + &
                                  mu ( ia , 2 ) * Tyz * mu ( ja , 3 ) + &
                                  mu ( ia , 3 ) * Txz * mu ( ja , 1 ) + &
                                  mu ( ia , 3 ) * Tyz * mu ( ja , 2 ) + &
                                  mu ( ia , 3 ) * Tzz * mu ( ja , 3 ) ) * dm5
              
                phi ( ia ) = phi ( ia ) + ( rxij * mu ( ja , 1 ) +  ryij * mu ( ja , 2 ) + rzij * mu ( ja , 3 ) ) * dm3

                Efield ( ia , 1 ) = Efield ( ia , 1 ) + ( Txx * mu ( ja , 1 ) + Txy * mu ( ja , 2 ) + Txz * mu ( ja , 3 ) ) * dm5
                Efield ( ia , 2 ) = Efield ( ia , 2 ) + ( Txy * mu ( ja , 1 ) + Tyy * mu ( ja , 2 ) + Tyz * mu ( ja , 3 ) ) * dm5
                Efield ( ia , 3 ) = Efield ( ia , 3 ) + ( Txz * mu ( ja , 1 ) + Tyz * mu ( ja , 2 ) + Tzz * mu ( ja , 3 ) ) * dm5

                Txxx = 5.0d0 * rxij * rxij * rxij -  3.0d0 * d2 * ( rxij )
                Tyyy = 5.0d0 * ryij * ryij * ryij -  3.0d0 * d2 * ( ryij )
                Tzzz = 5.0d0 * rzij * rzij * rzij -  3.0d0 * d2 * ( rzij )

                Txxy = 5.0d0 * rxij * rxij * ryij -          d2 * ( ryij )
                Txxz = 5.0d0 * rxij * rxij * rzij -          d2 * ( rzij )
                Tyyx = 5.0d0 * ryij * ryij * rxij -          d2 * ( rxij )
                Tyyz = 5.0d0 * ryij * ryij * rzij -          d2 * ( rzij )
                Tzzx = 5.0d0 * rzij * rzij * rxij -          d2 * ( rxij )
                Tzzy = 5.0d0 * rzij * rzij * ryij -          d2 * ( ryij )

                Txyz = 5.0d0 * rxij * ryij * rzij
          
                fxij = ( mu ( ia , 1 ) * Txxx * mu ( ja , 1 )         + &
                         mu ( ia , 1 ) * Txxy * mu ( ja , 2 ) * 2.0d0 + &
                         mu ( ia , 1 ) * Txxz * mu ( ja , 3 ) * 2.0d0 + &
                         mu ( ia , 2 ) * Txyz * mu ( ja , 3 ) * 2.0d0 + &       ! yxz = xyz
                         mu ( ia , 2 ) * Tyyx * mu ( ja , 2 )         + &       ! yxy = xyy
                         mu ( ia , 3 ) * Tzzx * mu ( ja , 3 ) ) * 3.0d0 * dm7   ! zxz = xzz
          
          
                fyij = ( mu ( ia , 2 ) * Tyyy * mu ( ja , 2 )         + &
                         mu ( ia , 1 ) * Tyyx * mu ( ja , 2 ) * 2.0d0 + &
                         mu ( ia , 2 ) * Tyyz * mu ( ja , 3 ) * 2.0d0 + &  
                         mu ( ia , 2 ) * Txyz * mu ( ja , 3 ) * 2.0d0 + &       ! yxz = xyz
                         mu ( ia , 1 ) * Txxy * mu ( ja , 1 )         + &
                         mu ( ia , 3 ) * Tzzy * mu ( ja , 3 ) ) * 3.0d0 * dm7 ! zyz = yzz
          
                fzij = ( mu ( ia , 3 ) * Tzzz * mu ( ja , 3 )         + & 
                         mu ( ia , 3 ) * Tzzx * mu ( ja , 1 ) * 2.0d0 + & 
                         mu ( ia , 2 ) * Tzzy * mu ( ja , 3 ) * 2.0d0 + &  
                         mu ( ia , 2 ) * Txyz * mu ( ja , 1 ) * 2.0d0 + &  
                         mu ( ia , 1 ) * Txxz * mu ( ja , 1 )         + &
                         mu ( ia , 2 ) * Tyyz * mu ( ja , 2 ) ) * 3.0d0 * dm7 
          
                fx_dir ( ia ) = fx_dir ( ia ) - fxij
                fy_dir ( ia ) = fy_dir ( ia ) - fyij
                fz_dir ( ia ) = fz_dir ( ia ) - fzij
          
                vir_dir =  vir_dir + ( fxij * rxij + fyij * ryij + fzij * rzij )
          
          
              endif ! d2.lt.cutsq
          
            enddo ! ja

         endif
         ! =======================================
         !  ia and ja in the same cell (ia.ne.ja)
         ! =======================================
         if ( rm_coul%lcell(ncell) .eq. 0) then

            do ja = 1,natm

              if (ja.ne.ia) then
                rxj  = rx(ja)
                ryj  = ry(ja)
                rzj  = rz(ja)
                rxij = rxi - rxj
                ryij = ryi - ryj
                rzij = rzi - rzj
                d2   = rxij * rxij + ryij * ryij + rzij * rzij
          
                if (d2.lt.cutsq) then
                  d = SQRT (d2)
                  dm3 = 1.0d0 / ( d2 * d ) 
                  dm5 = dm3 / d2
                  dm7 = dm5 / d2 
          
                  Txx = 3.0d0 * rxij * rxij - d2 
                  Tyy = 3.0d0 * ryij * ryij - d2 
                  Tzz = 3.0d0 * rzij * rzij - d2 
                  Txy = 3.0d0 * rxij * ryij  
                  Txz = 3.0d0 * rxij * rzij  
                  Tyz = 3.0d0 * ryij * rzij  
          
                  u_dir = u_dir - ( mu ( ia , 1 ) * Txx * mu ( ja , 1 ) + &
                                    mu ( ia , 1 ) * Txy * mu ( ja , 2 ) + &
                                    mu ( ia , 1 ) * Txz * mu ( ja , 3 ) + &
                                    mu ( ia , 2 ) * Txy * mu ( ja , 1 ) + &
                                    mu ( ia , 2 ) * Tyy * mu ( ja , 2 ) + &
                                    mu ( ia , 2 ) * Tyz * mu ( ja , 3 ) + & 
                                    mu ( ia , 3 ) * Txz * mu ( ja , 1 ) + &
                                    mu ( ia , 3 ) * Tyz * mu ( ja , 2 ) + &
                                    mu ( ia , 3 ) * Tzz * mu ( ja , 3 ) ) * dm5
          
                  phi ( ia ) = phi ( ia ) + ( rxij * mu ( ja , 1 ) +  ryij * mu ( ja , 2 ) + rzij * mu ( ja , 3 ) ) * dm3
          
                  Efield ( ia , 1 ) = Efield ( ia , 1 ) + ( Txx * mu ( ja , 1 ) + Txy * mu ( ja , 2 ) + Txz * mu ( ja , 3 ) ) * dm5
                  Efield ( ia , 2 ) = Efield ( ia , 2 ) + ( Txy * mu ( ja , 1 ) + Tyy * mu ( ja , 2 ) + Tyz * mu ( ja , 3 ) ) * dm5
                  Efield ( ia , 3 ) = Efield ( ia , 3 ) + ( Txz * mu ( ja , 1 ) + Tyz * mu ( ja , 2 ) + Tzz * mu ( ja , 3 ) ) * dm5
          
                  Txxx = 5.0d0 * rxij * rxij * rxij -  3.0d0 * d2 * ( rxij )
                  Tyyy = 5.0d0 * ryij * ryij * ryij -  3.0d0 * d2 * ( ryij )
                  Tzzz = 5.0d0 * rzij * rzij * rzij -  3.0d0 * d2 * ( rzij )

                  Txxy = 5.0d0 * rxij * rxij * ryij -          d2 * ( ryij )
                  Txxz = 5.0d0 * rxij * rxij * rzij -          d2 * ( rzij )
                  Tyyx = 5.0d0 * ryij * ryij * rxij -          d2 * ( rxij )
                  Tyyz = 5.0d0 * ryij * ryij * rzij -          d2 * ( rzij )
                  Tzzx = 5.0d0 * rzij * rzij * rxij -          d2 * ( rxij )
                  Tzzy = 5.0d0 * rzij * rzij * ryij -          d2 * ( ryij )

                  Txyz = 5.0d0 * rxij * ryij * rzij
         

                  fxij = ( mu ( ia , 1 ) * Txxx * mu ( ja , 1 )         + &
                           mu ( ia , 1 ) * Txxy * mu ( ja , 2 ) * 2.0d0 + &
                           mu ( ia , 1 ) * Txxz * mu ( ja , 3 ) * 2.0d0 + &
                           mu ( ia , 2 ) * Txyz * mu ( ja , 3 ) * 2.0d0 + & 
                           mu ( ia , 2 ) * Tyyx * mu ( ja , 2 )         + & 
                           mu ( ia , 3 ) * Tzzx * mu ( ja , 3 ) ) * 3.0d0 * dm7 
  
  
                  fyij = ( mu ( ia , 2 ) * Tyyy * mu ( ja , 2 )         + &
                           mu ( ia , 1 ) * Tyyx * mu ( ja , 2 ) * 2.0d0 + &
                           mu ( ia , 2 ) * Tyyz * mu ( ja , 3 ) * 2.0d0 + &        
                           mu ( ia , 2 ) * Txyz * mu ( ja , 3 ) * 2.0d0 + &
                           mu ( ia , 1 ) * Txxy * mu ( ja , 1 )         + &
                           mu ( ia , 3 ) * Tzzy * mu ( ja , 3 ) ) * 3.0d0 * dm7
  
                  fzij = ( mu ( ia , 3 ) * Tzzz * mu ( ja , 3 )         + &
                           mu ( ia , 3 ) * Tzzx * mu ( ja , 1 ) * 2.0d0 + &
                           mu ( ia , 2 ) * Tzzy * mu ( ja , 3 ) * 2.0d0 + & 
                           mu ( ia , 2 ) * Txyz * mu ( ja , 1 ) * 2.0d0 + &
                           mu ( ia , 1 ) * Txxz * mu ( ja , 1 )         + &
                           mu ( ia , 2 ) * Tyyz * mu ( ja , 2 ) ) * 3.0d0 * dm7
 
          
                  fx_dir ( ia ) = fx_dir ( ia ) - fxij
                  fy_dir ( ia ) = fy_dir ( ia ) - fyij
                  fz_dir ( ia ) = fz_dir ( ia ) - fzij
          
                  vir_dir =  vir_dir + ( fxij * rxij + fyij * ryij + fzij * rzij )
          
                endif ! d2.lt.cutsq
          
              endif ! ia.ne.ja

            enddo ! ja

          endif

       enddo ! ncell

  enddo atom


  CALL MPI_ALL_REDUCE_DOUBLE_SCALAR ( u_dir )
  CALL MPI_ALL_REDUCE_DOUBLE_SCALAR ( vir_dir )

  CALL MPI_ALL_REDUCE_DOUBLE ( fx_dir , natm )
  CALL MPI_ALL_REDUCE_DOUBLE ( fy_dir , natm )
  CALL MPI_ALL_REDUCE_DOUBLE ( fz_dir , natm )

  u_coul_dd   =   u_dir    * 0.5d0
  vir_coul_dd = - vir_dir  * 0.5d0
  phi_coul_dd = phi

  fx = fx + fx_dir
  fy = fy + fy_dir
  fz = fz + fz_dir

#ifdef debug
  write( stdout , '(5(A,f12.8))') 'U_dd  = ', u_coul_dd,' phi_dd,1   = ', phi(1),' E1,x = ',Efield(1,1),' f1,x = ', fx(1) , ' vir_dd = ',vir_coul_dd
  write( stdout , '(5(A,f12.8))') 'U_dir = ', u_dir,    ' vir_dir    = ',vir_dir
#endif

  return

END SUBROUTINE engforce_dipole_DS

!*********************** SUBROUTINE field_chargeanddipole_DS **********************
!
! this subroutine calculates the electric field at atoms position from all the 
! static dipoles of the system. The Direct Summation is used
! 
! Note : same consruction as efg_DS
!                          
!  Ei,alpha  =  sum_j!=i  T^alpha,beta_ij * mu_j,beta  
!
!  T^alpha,alpha_ij  = ( 3 * r_ij,alpha * r_ij,beta - r^2 delta_alpha,beta ) / r_ij^5 
!
!  
!******************************************************************************

SUBROUTINE engforce_charge_and_dipole_DS ( iastart , iaend , mu )

  USE config,   ONLY : natm , ntype , rx , ry , rz , fx , fy, fz , dipia , phi_coul_dd , qia
  USE control,  ONLY : cutlongrange
  USE thermodynamic,    ONLY :  u_coul_qd , vir_coul_qd
  USE io_file

  implicit none

  ! global 
  integer, intent(in) :: iastart , iaend 
  double precision :: mu    ( natm , 3) 

  ! local
  integer :: ia , ja  , it , ncell
  double precision :: d , d2 , dm3 , dm5 , qi , qj 
  double precision :: Txx , Tyy , Tzz , Txy , Txz , Tyz
  double precision :: rxi , ryi , rzi , rxij , ryij , rzij
  double precision :: rxj , ryj , rzj 
  double precision :: fxij , fyij , fzij 
!  double precision :: Txxx,  Tyyy,  Tzzz, Txxy, Txxz, Tyyx, Tyyz, Tzzx, Tzzy, Txyz
  double precision :: cutsq , u_dir , vir_dir
  double precision, dimension(:), allocatable :: fx_dir, fy_dir, fz_dir
  !double precision ::  phi ( natm ) 
  logical :: ldipole

  ldipole = .false.
  do it = 1 , ntype 
    if ( mu(it,1) .ne. 0.0d0 .or. & 
    mu(it,2) .ne. 0.0d0 .or. &
    mu(it,3) .ne. 0.0d0 ) ldipole = .true.
  enddo
  if ( .not. ldipole .and. ionode ) WRITE ( stdout , '(a)') 'No static dipoles' 
  if ( .not. ldipole ) return

  cutsq = cutlongrange * cutlongrange
#ifdef debug
  write ( stdout , '(a,i6)') 'debug : rm_coul',rm_coul%ncmax
  do ia = 1 , natm
  write( * , '(a,3f12.5)' )  'debug : dipole ', mu (ia,1), mu (ia,2), mu (ia,3)
  enddo
#endif  

  vir_dir = 0.0d0
  u_dir   = 0.0d0

  allocate( fx_dir(natm), fy_dir(natm), fz_dir(natm) )
  fx_dir  = 0.0D0
  fy_dir  = 0.0D0
  fz_dir  = 0.0D0


  atom : do ia = iastart , iaend
           rxi = rx(ia)
           ryi = ry(ia)
           rzi = rz(ia)
           qi = qia(ia)
          ! ==================================================
          ! sum over neighboring cells (see direct_sum_init )
          ! ==================================================
          do ncell = 1 , rm_coul%ncmax
          ! ==============================
          !  ia and ja in different cells
          ! ==============================
          if ( rm_coul%lcell(ncell) .eq. 1) then
 
            do ja = 1 , natm
              rxj  = rx(ja) + rm_coul%boxxyz(1,ncell)
              ryj  = ry(ja) + rm_coul%boxxyz(2,ncell)
              rzj  = rz(ja) + rm_coul%boxxyz(3,ncell)
              rxij = rxi - rxj
              ryij = ryi - ryj
              rzij = rzi - rzj
              d2   = rxij * rxij + ryij * ryij + rzij * rzij

              if ( d2 .lt. cutsq ) then
                qj = qia ( ja ) 
                d  = SQRT (d2)
                dm3 = 1.0d0/ ( d2 * d )
                dm5 = dm3 / d2 


                u_dir = u_dir + ( qi * ( rxij * mu ( ja , 1 ) + ryij * mu ( ja , 2 ) + rzij * mu ( ja , 3 ) ) - &
                                  qj * ( rxij * mu ( ia , 1 ) + ryij * mu ( ia , 2 ) + rzij * mu ( ia , 3 ) ) ) * dm3
                        

                Txx = 3.0d0 * rxij * rxij - d2 
                Tyy = 3.0d0 * ryij * ryij - d2 
                Tzz = 3.0d0 * rzij * rzij - d2 
                Txy = 3.0d0 * rxij * ryij      
                Txz = 3.0d0 * rxij * rzij     
                Tyz = 3.0d0 * ryij * rzij      

                fxij = qi * ( Txx * mu ( ja , 1 ) + Txy * mu ( ja , 2 ) + Txz * mu ( ja , 3 ) ) * dm5 - &
                       qj * ( Txx * mu ( ia , 1 ) + Txy * mu ( ia , 2 ) + Txz * mu ( ia , 3 ) ) * dm5

                fyij = qi * ( Txy * mu ( ja , 1 ) + Tyy * mu ( ja , 2 ) + Tyz * mu ( ja , 3 ) ) * dm5 - &
                       qj * ( Txy * mu ( ia , 1 ) + Tyy * mu ( ia , 2 ) + Tyz * mu ( ia , 3 ) ) * dm5

                fzij = qi * ( Txz * mu ( ja , 1 ) + Tyz * mu ( ja , 2 ) + Tzz * mu ( ja , 3 ) ) * dm5 - &
                       qj * ( Txz * mu ( ia , 1 ) + Tyz * mu ( ia , 2 ) + Tzz * mu ( ia , 3 ) ) * dm5
 
          
                fx_dir ( ia ) = fx_dir ( ia ) + fxij 
                fy_dir ( ia ) = fy_dir ( ia ) + fyij 
                fz_dir ( ia ) = fz_dir ( ia ) + fzij 
          
                vir_dir =  vir_dir + ( fxij * rxij + fyij * ryij + fzij * rzij )
          
          
              endif ! d2.lt.cutsq
          
            enddo ! ja

         endif
         ! =======================================
         !  ia and ja in the same cell (ia.ne.ja)
         ! =======================================
         if ( rm_coul%lcell(ncell) .eq. 0) then

            do ja = 1,natm

              if (ja.ne.ia) then
                rxj  = rx(ja)
                ryj  = ry(ja)
                rzj  = rz(ja)
                rxij = rxi - rxj
                ryij = ryi - ryj
                rzij = rzi - rzj
                d2   = rxij * rxij + ryij * ryij + rzij * rzij
          
                if (d2.lt.cutsq) then
                  qj = qia ( ja ) 
                  d = SQRT (d2)
                  dm3 = 1.0d0 / ( d2 * d )
                  dm5 = dm3 / d2

                  u_dir = u_dir + ( qi * ( rxij * mu ( ja , 1 ) + ryij * mu ( ja , 2 ) + rzij * mu ( ja , 3 ) ) - &
                                    qj * ( rxij * mu ( ia , 1 ) + ryij * mu ( ia , 2 ) + rzij * mu ( ia , 3 ) ) ) * dm3

                  Txx = 3.0d0 * rxij * rxij - d2 
                  Tyy = 3.0d0 * ryij * ryij - d2 
                  Tzz = 3.0d0 * rzij * rzij - d2 
                  Txy = 3.0d0 * rxij * ryij      
                  Txz = 3.0d0 * rxij * rzij      
                  Tyz = 3.0d0 * ryij * rzij      

                  fxij = qi * ( Txx * mu ( ja , 1 ) + Txy * mu ( ja , 2 ) + Txz * mu ( ja , 3 ) ) * dm5 - &
                         qj * ( Txx * mu ( ia , 1 ) + Txy * mu ( ia , 2 ) + Txz * mu ( ia , 3 ) ) * dm5

                  fyij = qi * ( Txy * mu ( ja , 1 ) + Tyy * mu ( ja , 2 ) + Tyz * mu ( ja , 3 ) ) * dm5 - &
                         qj * ( Txy * mu ( ia , 1 ) + Tyy * mu ( ia , 2 ) + Tyz * mu ( ia , 3 ) ) * dm5

                  fzij = qi * ( Txz * mu ( ja , 1 ) + Tyz * mu ( ja , 2 ) + Tzz * mu ( ja , 3 ) ) * dm5 - &
                         qj * ( Txz * mu ( ia , 1 ) + Tyz * mu ( ia , 2 ) + Tzz * mu ( ia , 3 ) ) * dm5


                  fx_dir ( ia ) = fx_dir ( ia ) + fxij
                  fy_dir ( ia ) = fy_dir ( ia ) + fyij 
                  fz_dir ( ia ) = fz_dir ( ia ) + fzij 

                  vir_dir =  vir_dir + ( fxij * rxij + fyij * ryij + fzij * rzij )

                endif ! d2.lt.cutsq
          
              endif ! ia.ne.ja

            enddo ! ja

          endif

       enddo ! ncell

  enddo atom


  CALL MPI_ALL_REDUCE_DOUBLE_SCALAR ( u_dir )
  CALL MPI_ALL_REDUCE_DOUBLE_SCALAR ( vir_dir )

  CALL MPI_ALL_REDUCE_DOUBLE ( fx_dir , natm )
  CALL MPI_ALL_REDUCE_DOUBLE ( fy_dir , natm )
  CALL MPI_ALL_REDUCE_DOUBLE ( fz_dir , natm )

  u_coul_qd   =   u_dir    * 0.5d0
  vir_coul_qd = - vir_dir  * 0.5d0

  fx = fx + fx_dir
  fy = fy + fy_dir
  fz = fz + fz_dir

  return

END SUBROUTINE engforce_charge_and_dipole_DS

!*********************** SUBROUTINE field_dipole_ES ************************************
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

SUBROUTINE engforce_dipole_ES ( iastart , iaend , Efield , mu )

  USE control,    ONLY :  myrank , numprocs, calc
  USE config,     ONLY :  system , natm , natmi , atype , atypei , box , &
                  omega , itype , rx , ry , rz , ntype , qia , dipia , fx , fy ,fz , phi_coul_dd
  USE constants,  ONLY :  pi , tpi , fpi , piroot , imag , mimag 
  USE kspace,     ONLY :  struc_fact , struc_fact_dip
  USE prop,       ONLY :  nprop_print , nprop
  USE thermodynamic,    ONLY :  u_coul_dd , vir_coul_dd
  USE time
  USE io_file 

  implicit none

  INCLUDE 'mpif.h'

  ! global
  integer, intent(in) :: iastart , iaend 
  double precision :: Efield ( natm , 3 ) 
  double precision :: mu    ( natm , 3 ) 

  ! local
  integer :: ia , ja ,  ierr , it , ip , ik 
  integer :: nxij , nyij , nzij
  double precision :: mip    ( ntype , 3 ) 
  double precision :: mutot  ( 3 ) , musq
  double precision :: d , d2 , d3 , d4 , d5 , d7 , expon , qi
  double precision :: alpha2 , alpha3 , alpha4 , alpha5
  double precision :: T0 , T1 , T2 , T3 ! real part 
  double precision :: rxi , ryi , rzi , rxij , ryij , rzij 
  double precision :: invbox 
  double precision :: ak, kx, ky, kz, kk
  double precision :: Txx , Tyy , Tzz , Txy , Txz , Tyz
  double precision :: Dxx , Dyy , Dzz , Dxy , Dxz , Dyz
  double precision :: Txxx,  Tyyy,  Tzzz, Txxy, Txxz, Tyyx, Tyyz, Tzzx, Tzzy, Txyz
  double complex   :: Tx , Ty , Tz 
  double precision :: Txf , Tyf , Tzf 
  double precision :: kri , kmui , k_dot_mu
  double precision :: muix , muiy , muiz 
  double precision :: mujx , mujy , mujz 
  double precision :: muir , mujr , muij
  double complex   :: rhon 
  double precision :: u_dir , u_rec , u_surf , u_self 
  double precision :: str , recarg_ef
  double complex   :: carg , carg_ef , carg_f
  double precision, external :: errfc 
  double precision :: ttt1 , ttt2 , ttt3 , ttt4 , ttt5 
  double precision :: Efield_real  ( natm , 3 )
  double precision :: Efield_surf  ( natm , 3 )
  double precision :: Efield_dual  ( natm , 3 )
  double precision :: fxij , fyij , fzij
  double precision :: C_r , D_r 
  double precision :: vir_dir , vir_rec , vir_surf, vir_self 
  double precision, dimension(:), allocatable :: fx_dir, fy_dir, fz_dir
  double precision, dimension(:), allocatable :: fx_rec, fy_rec, fz_rec
  double precision, dimension(:), allocatable :: fx_surf, fy_surf, fz_surf
  double precision, dimension(:), allocatable :: fx_self, fy_self, fz_self
  double precision, dimension(:), allocatable :: phi_dir , phi_rec , phi_surf , phi_self
  logical :: ldipole

  ldipole = .false.
  do it = 1 , ntype
    if ( mu(it,1) .ne. 0.0d0 .or. &
         mu(it,2) .ne. 0.0d0 .or. &
         mu(it,3) .ne. 0.0d0 ) ldipole = .true.
  enddo
  if ( .not. ldipole .and. ionode ) WRITE ( stdout , '(a)') 'No static dipoles'
  if ( .not. ldipole ) return
  ! =============================
  !  init mip ( itype dependent ) 
  ! ==============================
  ip = 0
  do it = 1, ntype
    ip = ip + natmi ( it )
    print*,it,natmi(it),ip
    mip ( it , 1 ) = mu ( ip , 1 )
    mip ( it , 2 ) = mu ( ip , 2 )
    mip ( it , 3 ) = mu ( ip , 3 )
  enddo

  ! =============================
  !  total moment / square
  ! =============================
  mutot = 0.0d0
  musq  = 0.0d0
  do ia = 1 , natm
    mutot ( 1 ) = mutot ( 1 ) + mu ( ia , 1 )
    mutot ( 2 ) = mutot ( 2 ) + mu ( ia , 2 )
    mutot ( 3 ) = mutot ( 3 ) + mu ( ia , 3 )
    musq = musq + ( mu ( ia , 1 ) * mu ( ia , 1 ) +  mu ( ia , 2 ) * mu ( ia , 2 ) +  mu ( ia , 3 ) * mu ( ia , 3 ) )
  enddo


#ifdef debug
  CALL print_config_sample(0,0)
  write(stdout , '(a,i6)' ) 'debug in engforce_dipole_ES : km_coul%nkcut',km_coul%nkcut
  write(stdout , '(a,i6)' ) 'debug in engforce_dipole_ES : km_coul_dip%nkcut',km_coul_dip%nkcut
  do ia= 1,natm
    write ( stdout ,'(a,i4,3f16.8)')  'debug in engforce_dipole_ES : ia = ', ia , mu  ( ia , 1 ) , mu  ( ia , 2 ) ,  mu  ( ia , 3 )
  enddo
  do it=1,ntype
    write ( stdout ,'(a,i4,3f16.8)')  'debug in engforce_dipole_ES : it = ', it , mip ( it , 1 ) , mip ( it , 2 ) ,  mip ( it , 3 )
  enddo
#endif

  ! ======================
  ! some local allocation
  ! ======================
  allocate( fx_dir(natm), fy_dir(natm), fz_dir(natm) )
  allocate( fx_rec(natm), fy_rec(natm), fz_rec(natm) )
  allocate( fx_surf(natm), fy_surf(natm), fz_surf(natm) )
  allocate( fx_self(natm), fy_self(natm), fz_self(natm) )
  allocate( phi_rec(natm) , phi_dir(natm) , phi_surf(natm) , phi_self ( natm ) )
  ! ==========================
  !  init some quantities
  ! ==========================
  Efield           = 0.0d0
  Efield_surf      = 0.0d0
  Efield_real      = 0.0d0
  Efield_dual      = 0.0d0 
  u_dir   = 0.0d0
  u_rec   = 0.0d0
  u_surf  = 0.0d0
  u_self  = 0.0d0
  vir_dir = 0.0d0
  vir_rec = 0.0d0
  vir_surf= 0.0d0
  phi_dir = 0.0d0
  phi_rec = 0.0d0
  phi_surf= 0.0d0
  phi_self = 0.0d0
  fx_dir  = 0.0D0
  fy_dir  = 0.0D0
  fz_dir  = 0.0D0
  fx_rec  = 0.0D0
  fy_rec  = 0.0D0
  fz_rec  = 0.0D0
  fx_surf  = 0.0D0
  fy_surf  = 0.0D0
  fz_surf = 0.0D0
  fx_self  = 0.0D0
  fy_self  = 0.0D0
  fz_self = 0.0D0

  ! =================
  !  some constants 
  ! =================
  invbox = 1.0d0 / box
  alpha2 = alphaES * alphaES
  alpha3 = alpha2  * alphaES
  alpha4 = alpha2  * alpha2
  alpha5 = alpha2  * alpha3

  ttt1 = MPI_WTIME(ierr)
  ! =====================
  ! facteur de structure 
  ! =====================
  CALL struc_fact ( km_coul )


  ttt2 = MPI_WTIME(ierr)
  strftimetot = strftimetot + (ttt2 - ttt1)
  ! ==============================================
  !        direct space part
  ! ==============================================

  atom1: do ia = iastart , iaend 

    rxi = rx(ia)
    ryi = ry(ia)
    rzi = rz(ia)
    qi = qia ( ia )
    muix = mu ( ia , 1 ) 
    muiy = mu ( ia , 2 ) 
    muiz = mu ( ia , 3 ) 

    do ja = 1, natm

      rxij = rxi - rx(ja)
      ryij = ryi - ry(ja)
      rzij = rzi - rz(ja)
      if (ja .ne. ia ) then
   
        mujx = mu ( ja , 1 ) 
        mujy = mu ( ja , 2 ) 
        mujz = mu ( ja , 3 ) 
         
        muij = muix * mujx + muiy * mujy + muiz * mujz 

        nxij = NINT ( rxij * invbox )
        nyij = NINT ( ryij * invbox )
        nzij = NINT ( rzij * invbox )
        rxij = rxij - box * nxij
        ryij = ryij - box * nyij
        rzij = rzij - box * nzij

        mujr = mujx * rxij + mujy * ryij + mujz * rzij

        muir = muix * rxij + muiy * ryij + muiz * rzij
 
        d2 = rxij * rxij + ryij * ryij + rzij * rzij
        d  = SQRT ( d2 )
        d3 = d2 * d 
        d4 = d2 * d2
        d5 = d4 * d
        d7 = d5 * d2

        expon = EXP ( - alpha2 * d2   ) / piroot
        T0 = errfc   ( alphaES * d     ) 
        T1 =         ( 2.0d0 * alphaES ) * d
        T2 =         ( 4.0d0 * alpha3  ) / d2 
        T3 =         ( 8.0d0 * alpha5  ) / 15.0d0 * d5 

        C_r = ( 3.0d0  * T0 + T1 * expon * ( 3.0d0 + 2.0d0 * alpha2 * d2 )  )  / d5
        D_r = ( 15.0d0 * T0 + T1 * expon * ( 15.0d0 + 10.0d0 * alpha2 * d2 + 4.0d0 * alpha4 * d4) ) / d7

        Dxx = rxij * rxij
        Dyy = ryij * ryij
        Dzz = rzij * rzij
        Dxy = rxij * ryij
        Dxz = rxij * rzij
        Dyz = ryij * rzij

        Txx = 3.0d0 * rxij * rxij - d2
        Tyy = 3.0d0 * ryij * ryij - d2
        Tzz = 3.0d0 * rzij * rzij - d2
        Txy = 3.0d0 * rxij * ryij
        Txz = 3.0d0 * rxij * rzij
        Tyz = 3.0d0 * ryij * rzij

        u_dir = u_dir - ( mu ( ia , 1 ) * Txx * mu ( ja , 1 ) + &
                          mu ( ia , 1 ) * Txy * mu ( ja , 2 ) + &
                          mu ( ia , 1 ) * Txz * mu ( ja , 3 ) + &
                          mu ( ia , 2 ) * Txy * mu ( ja , 1 ) + &
                          mu ( ia , 2 ) * Tyy * mu ( ja , 2 ) + &
                          mu ( ia , 2 ) * Tyz * mu ( ja , 3 ) + &
                          mu ( ia , 3 ) * Txz * mu ( ja , 1 ) + &
                          mu ( ia , 3 ) * Tyz * mu ( ja , 2 ) + &
                          mu ( ia , 3 ) * Tzz * mu ( ja , 3 ) ) * ( T0 + T1 * expon ) / d5

        u_dir = u_dir -  ( mu ( ia , 1 ) * Dxx * mu ( ja , 1 ) + &
                           mu ( ia , 1 ) * Dxy * mu ( ja , 2 ) + &
                           mu ( ia , 1 ) * Dxz * mu ( ja , 3 ) + &
                           mu ( ia , 2 ) * Dxy * mu ( ja , 1 ) + &
                           mu ( ia , 2 ) * Dyy * mu ( ja , 2 ) + &
                           mu ( ia , 2 ) * Dyz * mu ( ja , 3 ) + &
                           mu ( ia , 3 ) * Dxz * mu ( ja , 1 ) + &
                           mu ( ia , 3 ) * Dyz * mu ( ja , 2 ) + &
                           mu ( ia , 3 ) * Dzz * mu ( ja , 3 ) ) * T2 * expon
 
        phi_dir ( ia ) = phi_dir ( ia ) +  ( mu ( ja , 1 ) * rxij + mu ( ja , 2 ) * ryij + mu ( ja , 3 ) * rzij ) * ( T0  + T1 * expon ) / d3  

        Efield_real ( ia , 1 ) = Efield_real ( ia , 1 ) + &
        ( Txx * mu ( ja , 1 ) + Txy * mu ( ja , 2 ) + Txz * mu ( ja , 3 ) ) * ( T0 + T1 * expon ) / d5 + &
        ( Dxx * mu ( ja , 1 ) + Dxy * mu ( ja , 2 ) + Dxz * mu ( ja , 3 ) ) * ( T2 * expon ) 
        Efield_real ( ia , 2 ) = Efield_real ( ia , 2 ) + &
        ( Txy * mu ( ja , 1 ) + Tyy * mu ( ja , 2 ) + Tyz * mu ( ja , 3 ) ) * ( T0 + T1 * expon ) / d5 + &
        ( Dxy * mu ( ja , 1 ) + Dyy * mu ( ja , 2 ) + Dyz * mu ( ja , 3 ) ) * ( T2 * expon )
        Efield_real ( ia , 3 ) = Efield_real ( ia , 3 ) + &
        ( Txz * mu ( ja , 1 ) + Tyz * mu ( ja , 2 ) + Tzz * mu ( ja , 3 ) ) * ( T0 + T1 * expon ) / d5 + &
        ( Dxz * mu ( ja , 1 ) + Dyz * mu ( ja , 2 ) + Dzz * mu ( ja , 3 ) ) * ( T2 * expon )

        Txxx = 15.0d0 * rxij * rxij * rxij -  9.0d0 * d2 * ( rxij )
        Tyyy = 15.0d0 * ryij * ryij * ryij -  9.0d0 * d2 * ( ryij )
        Tzzz = 15.0d0 * rzij * rzij * rzij -  9.0d0 * d2 * ( rzij )

        Txxy = 15.0d0 * rxij * rxij * ryij -  3.0d0 * d2 * ( ryij )
        Txxz = 15.0d0 * rxij * rxij * rzij -  3.0d0 * d2 * ( rzij )
        Tyyx = 15.0d0 * ryij * ryij * rxij -  3.0d0 * d2 * ( rxij )
        Tyyz = 15.0d0 * ryij * ryij * rzij -  3.0d0 * d2 * ( rzij )
        Tzzx = 15.0d0 * rzij * rzij * rxij -  3.0d0 * d2 * ( rxij )
        Tzzy = 15.0d0 * rzij * rzij * ryij -  3.0d0 * d2 * ( ryij )

        Txyz = 15.0d0 * rxij * ryij * rzij

        fxij = ( muij * rxij + muix * mujr + muir * mujx ) * C_r - &
                 muir * mujr * D_r * rxij    

        fyij = ( muij * ryij + muiy * mujr + muir * mujy ) * C_r - &
                 muir * mujr * D_r * ryij    

        fzij = ( muij * rzij + muiz * mujr + muir * mujz ) * C_r - &
                 muir * mujr * D_r * rzij    

        fx_dir ( ia ) = fx_dir ( ia ) + fxij
        fy_dir ( ia ) = fy_dir ( ia ) + fyij
        fz_dir ( ia ) = fz_dir ( ia ) + fzij

        vir_dir =  vir_dir + ( fxij * rxij + fyij * ryij + fzij * rzij )

      endif

      phi_surf ( ia ) = phi_surf ( ia ) + 2.0d0 * ( mu ( ia , 1 ) * rxij + mu ( ia , 2 ) * ryij + mu ( ia , 3 ) * rzij  ) 

      Efield_surf ( ia , 1 ) = Efield_surf ( ia , 1 ) - mu ( ja , 1 ) / 3.0D0
      Efield_surf ( ia , 2 ) = Efield_surf ( ia , 2 ) - mu ( ja , 2 ) / 3.0D0
      Efield_surf ( ia , 3 ) = Efield_surf ( ia , 3 ) - mu ( ja , 3 ) / 3.0D0

    enddo

  enddo atom1


  ttt3 = MPI_WTIME(ierr)
  !  efgtimetot1 = efgtimetot1 + (ttt3-ttt2)

  ! ==============================================
  !            reciprocal space part
  ! ==============================================
  kpoint : do ik = 1, km_coul%nkcut

    kx = km_coul%kpt ( 1 , ik )
    ky = km_coul%kpt ( 2 , ik )
    kz = km_coul%kpt ( 3 , ik )
    kk = km_coul%kptk( ik )
    Ak = EXP ( - kk * 0.25d0 / alpha2 ) / kk
    if ( km_coul%kptk( ik ) .eq. 0 ) then 
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
      k_dot_mu = ( mip ( it , 1 ) * kx + mip ( it , 2 ) * ky + mip ( it , 3 ) * kz  )
      rhon = rhon + imag * CONJG ( k_dot_mu * km_coul%strf ( ik , it ) )
    enddo

    str = rhon * CONJG ( rhon ) 
    u_rec = u_rec + DBLE ( str * ak )

    do ia = 1 , natm
      rxi  = rx ( ia )
      ryi  = ry ( ia )
      rzi  = rz ( ia )
      muix = mu ( ia , 1 ) 
      muiy = mu ( ia , 2 ) 
      muiz = mu ( ia , 3 ) 
      kri  = ( kx * rxi + ky * ryi + kz * rzi ) 
      kmui = ( kx * muix + ky * muiy + kz * muiz )   
      carg = EXP ( imag * kri )

      carg_ef = imag * carg * CONJG ( rhon )
      carg_f  =        carg * CONJG ( rhon )

      phi_rec ( ia ) = phi_rec ( ia ) + DBLE ( ak * carg * CONJG ( rhon ) )

      Tx = DBLE ( kx * ak * ( carg_ef - kmui )     ) 
      Ty = DBLE ( ky * ak * ( carg_ef - kmui )    )
      Tz = DBLE ( kz * ak * ( carg_ef - kmui )      )

      Txf = DBLE ( kx * ak * kmui * carg_f )
      Tyf = DBLE ( ky * ak * kmui * carg_f )
      Tzf = DBLE ( kz * ak * kmui * carg_f )

      Efield_dual ( ia , 1 ) = Efield_dual ( ia , 1 ) - Tx 
      Efield_dual ( ia , 2 ) = Efield_dual ( ia , 2 ) - Ty 
      Efield_dual ( ia , 3 ) = Efield_dual ( ia , 3 ) - Tz

      fx_rec ( ia ) = fx_rec ( ia ) + Txf 
      fy_rec ( ia ) = fy_rec ( ia ) + Tyf 
      fz_rec ( ia ) = fz_rec ( ia ) + Tzf 

      vir_rec = vir_rec + DBLE ( ak * str * 0.5d0 * ( 1.0d0 - ( kk * 0.5d0 / alpha2 ) )  + str * ak )

    enddo

  enddo kpoint

  ttt4 = MPI_WTIME(ierr)
  !  efgtimetot2 = efgtimetot2 + (ttt4-ttt3)

  !=============================
  !  MERGE REAL PART
  !=============================
  CALL MPI_ALL_REDUCE_DOUBLE ( Efield_real ( : , 1 ) , natm ) 
  CALL MPI_ALL_REDUCE_DOUBLE ( Efield_real ( : , 2 ) , natm ) 
  CALL MPI_ALL_REDUCE_DOUBLE ( Efield_real ( : , 3 ) , natm ) 

  ! =======
  ! 4pi/V
  ! =======
  ! self 
  do ia = 1 , natm 
    Efield_dual ( ia , 1 ) = Efield_dual ( ia , 1 ) + mu ( ia , 1 ) / 3.0d0
    Efield_dual ( ia , 2 ) = Efield_dual ( ia , 2 ) + mu ( ia , 2 ) / 3.0D0
    Efield_dual ( ia , 3 ) = Efield_dual ( ia , 3 ) + mu ( ia , 3 ) / 3.0D0
  enddo

  Efield_dual = Efield_dual  * fpi / omega 
  Efield_surf = Efield_surf  * fpi / omega 

#ifdef debug
  do ip=1,natm
    WRITE ( stdout ,'(a,i4,a,3f12.8)') 'atom = ',ip,' Efield_real         = ', &
    Efield_real      ( ip , 1 ) , Efield_real      ( ip , 2 ) , Efield_real      ( ip , 3 ) 
  enddo
  do ip=1,natm
    WRITE ( stdout ,'(a,i4,a,3f12.8)') 'atom = ',ip,' Efield_dual         = ', &
    Efield_dual ( ip , 1 ) , Efield ( ip , 2 ) , Efield ( ip , 3 ) 
  enddo
  do ip=1,natm
    WRITE ( stdout ,'(a,i4,a,3f12.8)') 'atom = ',ip,' Efield_surf         = ', &
    Efield_surf      ( ip , 1 ) , Efield_surf      ( ip , 2 ) , Efield_surf      ( ip , 3 ) 
  enddo
#endif
  ! =======================
  !  total electric field 
  ! =======================
  Efield  = Efield_dual + Efield_real 
#ifdef debug
  do ip=1,natm
    WRITE ( stdout ,'(a,i4,a,3f12.8)') 'atom = ',ip,' Efield_without_surf = ',   Efield   ( ip , 1 ) , Efield  ( ip , 2 ) , Efield  ( ip , 3 )   
  enddo
#endif
  Efield  = Efield_dual + Efield_real + Efield_surf
#ifdef debug
  do ip=1,natm
    WRITE ( stdout ,'(a,i4,a,3f12.8)') 'atom = ',ip,' Efield_with_surf    = ',   Efield   ( ip , 1 ) , Efield  ( ip , 2 ) , Efield  ( ip , 3 )  
  enddo
#endif
  fx_rec = fx_rec * fpi / omega
  fy_rec = fy_rec * fpi / omega
  fz_rec = fz_rec * fpi / omega
 
  ! f_surf = 0.0d0 without charge
  fx_surf = 0.0d0
  fy_surf = 0.0d0
  fz_surf = 0.0d0

  fx = fx_rec + fx_dir - fx_surf
  fy = fy_rec + fy_dir - fy_surf
  fz = fz_rec + fz_dir - fz_surf
 
  u_dir  = u_dir * 0.5d0
  u_rec  = u_rec * tpi / omega
  u_surf = mutot ( 1 ) * mutot ( 1 ) + mutot ( 2 ) * mutot ( 2 ) + mutot ( 3 ) * mutot ( 3 ) 
  u_surf = u_surf * tpi / omega / 3.0d0 
  u_self = - musq * 2.0d0 * alpha3 / 3.0d0 / piroot

  phi_rec = phi_rec * fpi / omega
  phi_self = 0.0d0
  phi_surf = phi_surf * tpi / omega / 3.0d0

  vir_dir = vir_dir * 0.5d0
  vir_rec = - vir_rec * fpi / omega
  vir_surf = mutot ( 1 ) * mutot ( 1 ) + mutot ( 2 ) * mutot ( 2 ) + mutot ( 3 ) * mutot ( 3 ) 
  vir_surf = - vir_surf * tpi / omega 

  u_coul_dd   = ( u_dir + u_rec + u_surf + u_self )
  phi_coul_dd = ( phi_dir + phi_rec + phi_surf + phi_self )
  vir_coul_dd = ( vir_dir + vir_rec + vir_surf )

#ifdef debug
  write( stdout , '(4(a,f16.8))' ) ,'vir_dir    = ', vir_dir    ,' vir_rec    = ', vir_rec    ,' vir_surf    = ',vir_surf   ,' vir_self    = '   ,vir_self   ,' vir_coul_dd    = ',vir_coul_dd   
  write( stdout , '(5(a,f16.8))' ) ,'u_dir      = ', u_dir      ,' u_rec      = ', u_rec      ,' u_surf      = ',u_surf     ,' u_self      = '   ,u_self     ,' u_coul_dd      = ',u_coul_dd
  write( stdout , '(5(a,f16.8))' ) ,'phi_dir(1) = ', phi_dir(1) ,' phi_rec(1) = ', phi_rec(1) ,' phi_surf(1) = ',phi_surf(1),' phi_self(1) = '   ,phi_self(1),' phi_coul_dd(1) = ',phi_coul_dd(1)
  write( stdout , '(5(a,f16.8))' ) ,'fx_dir(1)  = ', fx_dir(1)  ,' fx_rec(1)  = ', fx_rec(1)  ,' fx_surf(1)  = ',fx_surf(1) ,' fx_self(1)  = '   ,fx_self(1) ,' fx(1)          = ',fx(1)
#endif


  ttt5 = MPI_WTIME(ierr)
!  efgtimetot3 = efgtimetot3 + (ttt5-ttt4)

  CALL MPI_BARRIER( MPI_COMM_WORLD , ierr )

  deallocate( fx_dir , fy_dir  , fz_dir  )
  deallocate( fx_rec , fy_rec  , fz_rec  )
  deallocate( fx_surf, fy_surf , fz_surf  )
  deallocate( fx_self, fy_self , fz_self  )
  deallocate( phi_rec, phi_dir , phi_surf , phi_self )

  return

END SUBROUTINE engforce_dipole_ES


!*********************** SUBROUTINE induced_moment *****************************
!
! this subroutine calculates the induced moment from the total Electric field and the
! polarizability tensor
! Basically used in the SCF loop to get induced_moment
! 
! Formula : 
!                          
!  mu_i,alpha =  POL_i,alpha,beta * E_i,beta   
!
! POL is the polarizability tensor
!
!******************************************************************************

SUBROUTINE induced_moment ( Efield , mu_ind , u_pol )

  USE config, ONLY : natm , itype , atypei, ntype
  USE io_file, ONLY : stdout

  implicit none

  ! global
  double precision , intent ( in  ) :: Efield ( natm , 3 ) 
  double precision , intent ( out ) :: mu_ind ( natm , 3 ) 
  double precision , intent ( out ) :: u_pol   

  ! local 
  integer :: alpha , beta
  integer :: ia , it 
  double precision :: invpol ( 3 , 3 )   
  integer, parameter :: LWORK=1000  
  double precision :: WORK ( LWORK ) 
  integer :: ipiv ( 3 ) 
  integer :: ierr

  mu_ind = 0.0d0
  do ia = 1 , natm 
    do alpha = 1 , 3 
      do beta = 1 , 3  
        mu_ind ( ia , alpha ) = mu_ind ( ia , alpha ) + pol ( itype(ia) , alpha , beta ) * Efield ( ia , beta )  
      enddo
    enddo
  enddo

  u_pol = 0.0d0 
  do ia = 1 , natm
    it = itype(ia) 
    if ( lpolar ( it ) .eq. 1 ) then 
      invpol ( : , : )  = pol ( it , : , : )
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
  u_pol = u_pol * 0.5d0

  return

END SUBROUTINE induced_moment


SUBROUTINE multipole_DS ( iastart, iaend , ef , efg , mu , u_coul , vir_coul , phi_coul )

  USE config,  ONLY : natm , atype , qia , rx , ry , rz , fx , fy , fz 
  USE control, ONLY : cutlongrange
  USE io_file, ONLY : stdout , ionode 
  USE thermodynamic , ONLY : u_pol

  implicit none

  INCLUDE 'mpif.h'

  ! global 
  integer, intent(in) :: iastart , iaend
  double precision    :: ef     ( natm , 3 )
  double precision    :: efg    ( natm , 3 , 3 )
  double precision    :: mu     ( natm , 3 )
  double precision    :: u_coul 
  double precision    :: vir_coul
  double precision    :: phi_coul ( natm ) 

  ! local 
  integer :: ia, ja , ierr , ncell , it
  double precision, dimension(:), allocatable :: fx_coul , fy_coul , fz_coul
  double precision, dimension(:), allocatable :: phi_coul_qq , phi_coul_dd 
  double precision :: u_coul_qq , u_coul_dd , u_coul_qd 
  double precision :: vir_coul_qq , vir_coul_dd , vir_coul_qd
  double precision :: cutsq
  double precision :: rxi , ryi , rzi 
  double precision :: rxj , ryj , rzj
  double precision :: rxij , ryij , rzij
  double precision :: fxij , fyij , fzij
  double precision :: qi , qj , qij
  double precision :: muix , muiy , muiz
  double precision :: mujx , mujy , mujz
  double precision :: T
  double precision :: Tx , Ty , Tz
  double precision :: Txx , Tyy , Tzz , Txy , Txz , Tyz
  double precision :: Txxx,  Tyyy,  Tzzz, Txxy, Txxz, Tyyx, Tyyz, Tzzx, Tzzy, Txyz
  double precision :: d , d2 
  double precision :: dm1 , dm3 , dm5 , dm7 

  double precision :: ttt1 , ttt2  , ttt3
  logical          :: lcentralbox

  ttt1 = MPI_WTIME(ierr)

  ! =============================== 
  !         some constants
  ! =============================== 
  cutsq = cutlongrange * cutlongrange


#ifdef debug3
  write( stdout , '(a)')    'debug: in multipole_DS'
  write( stdout , '(a,i8,a)') 'debug : rm_coul ',rm_coul%ncmax,rm_coul%meshlabel
  do ia = 1 , natm
    write( stdout , '(a,f12.5)') 'debug : charge',qia(ia)
  enddo
  do ia = 1 , natm
    write( stdout , '(a,3f12.5)') 'debug : dipole', mu ( ia , 1 ),  mu ( ia , 2 ) ,  mu ( ia , 3 )
  enddo
  write( stdout , '(a,2i8)')     'debug : iastart iaend',iastart ,iaend
  write( stdout , '(a,f20.5)')   'debug : cutsq ',cutsq
  call print_config_sample(0,0)
#endif  


  allocate( fx_coul ( natm ), fy_coul ( natm ), fz_coul ( natm ) )
  allocate( phi_coul_qq ( natm ), phi_coul_dd ( natm ) )
  ef          = 0.0d0
  efg         = 0.0d0
  vir_coul    = 0.0d0
  vir_coul_qq = 0.0d0
  vir_coul_qd = 0.0d0
  vir_coul_dd = 0.0d0
  u_coul      = 0.0d0
  u_coul_qq   = 0.0d0
  u_coul_qd   = 0.0d0
  u_coul_dd   = 0.0d0
  fx_coul     = 0.0D0
  fy_coul     = 0.0D0
  fz_coul     = 0.0D0
  phi_coul    = 0.0d0
  phi_coul_qq = 0.0d0
  phi_coul_dd = 0.0d0
  
! ==================================================================================================
!  MAIN LOOP calculate EF(i) , EFG(i) , virial , potential, forces ... for each atom i parallelized
! ==================================================================================================
  atom : do ia = iastart , iaend

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

            rxij  = rxi - rxj
            ryij  = ryi - ryj
            rzij  = rzi - rzj

            d2    = rxij * rxij + ryij * ryij + rzij * rzij

            if ( d2 .lt. cutsq .and. d2 .ne. 0.0d0 ) then 

              qj    = qia ( ja )
              qij   = qi * qj
              mujx  = mu ( ja , 1 )  
              mujy  = mu ( ja , 2 )  
              mujz  = mu ( ja , 3 )  
              d     = SQRT ( d2 )
              dm1   = 1.0d0 / d 
              dm3   = dm1 / d2
              dm5   = dm3 / d2
              dm7   = dm5 / d2
 
              ! multipole interaction tensor rank = 0 
              T  = dm1

              ! multipole interaction tensor rank = 1
              Tx = rxij * dm3 
              Ty = ryij * dm3 
              Tz = rzij * dm3 

              ! multipole interaction tensor rank = 2
              Txx = ( 3.0d0 * rxij * rxij - d2 ) * dm5
              Tyy = ( 3.0d0 * ryij * ryij - d2 ) * dm5
              Tzz = ( 3.0d0 * rzij * rzij - d2 ) * dm5
              Txy = ( 3.0d0 * rxij * ryij      ) * dm5
              Txz = ( 3.0d0 * rxij * rzij      ) * dm5
              Tyz = ( 3.0d0 * ryij * rzij      ) * dm5

              ! multipole interaction tensor rank = 3  
              Txxx = ( 5.0d0 * rxij * rxij * rxij -  3.0d0 * d2 * ( rxij ) ) * dm7 * 3.0d0
              Tyyy = ( 5.0d0 * ryij * ryij * ryij -  3.0d0 * d2 * ( ryij ) ) * dm7 * 3.0d0
              Tzzz = ( 5.0d0 * rzij * rzij * rzij -  3.0d0 * d2 * ( rzij ) ) * dm7 * 3.0d0
              Txxy = ( 5.0d0 * rxij * rxij * ryij -          d2 * ( ryij ) ) * dm7 * 3.0d0
              Txxz = ( 5.0d0 * rxij * rxij * rzij -          d2 * ( rzij ) ) * dm7 * 3.0d0
              Tyyx = ( 5.0d0 * ryij * ryij * rxij -          d2 * ( rxij ) ) * dm7 * 3.0d0
              Tyyz = ( 5.0d0 * ryij * ryij * rzij -          d2 * ( rzij ) ) * dm7 * 3.0d0
              Tzzx = ( 5.0d0 * rzij * rzij * rxij -          d2 * ( rxij ) ) * dm7 * 3.0d0
              Tzzy = ( 5.0d0 * rzij * rzij * ryij -          d2 * ( ryij ) ) * dm7 * 3.0d0
              Txyz = ( 5.0d0 * rxij * ryij * rzij                          ) * dm7 * 3.0d0

              ! ===========================================================
              !                  charge-charge interaction
              ! ===========================================================

              ! potential 
              phi_coul_qq ( ia ) = phi_coul_qq ( ia ) + qj * T

              ! electrostatic energy
              u_coul_qq       = u_coul_qq       + qij * T

              ! electric field
              ef ( ia , 1 ) = ef ( ia , 1 ) + qj * Tx
              ef ( ia , 2 ) = ef ( ia , 2 ) + qj * Ty
              ef ( ia , 3 ) = ef ( ia , 3 ) + qj * Ty 
 
              ! forces 
              fxij = qij * Tx
              fyij = qij * Ty
              fzij = qij * Tz

              fx_coul ( ia ) = fx_coul ( ia ) + fxij
              fy_coul ( ia ) = fy_coul ( ia ) + fyij
              fz_coul ( ia ) = fz_coul ( ia ) + fzij

              ! virial
              vir_coul_qq =  vir_coul_qq + ( fxij * rxij + fyij * ryij + fzij * rzij )

              ! electric field gradient
              efg ( ia , 1 , 1 )  = efg ( ia , 1 , 1 ) - qj * Txx 
              efg ( ia , 2 , 2 )  = efg ( ia , 2 , 2 ) - qj * Tyy 
              efg ( ia , 3 , 3 )  = efg ( ia , 3 , 3 ) - qj * Tzz 
              efg ( ia , 1 , 2 )  = efg ( ia , 1 , 2 ) - qj * Txy 
              efg ( ia , 1 , 3 )  = efg ( ia , 1 , 3 ) - qj * Txz  
              efg ( ia , 2 , 3 )  = efg ( ia , 2 , 3 ) - qj * Tyz


              ! ===========================================================
              !                  dipole-dipole interaction
              ! ===========================================================

              ! potential 
              phi_coul_dd ( ia ) = phi_coul_dd ( ia ) + Tx * mujx +  Ty * mujy + Tz * mujz 

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
                       muix * Txxy * mujy * 2.0d0 + &
                       muix * Txxz * mujz * 2.0d0 + &
                       muiy * Txyz * mujz * 2.0d0 + &   
                       muiy * Tyyx * mujy         + &    
                       muiz * Tzzx * mujz )

              fyij = ( muiy * Tyyy * mujy         + &
                       muix * Tyyx * mujy * 2.0d0 + &
                       muiy * Tyyz * mujz * 2.0d0 + &
                       muiy * Txyz * mujz * 2.0d0 + &     
                       muix * Txxy * mujx         + &
                       muiz * Tzzy * mujz )

              fzij = ( muiz * Tzzz * mujz         + &
                       muiz * Tzzx * mujx * 2.0d0 + &
                       muiy * Tzzy * mujz * 2.0d0 + &
                       muiy * Txyz * mujx * 2.0d0 + &
                       muix * Txxz * mujx         + &
                       muiy * Tyyz * mujy )

              fx_coul ( ia ) = fx_coul ( ia ) - fxij
              fy_coul ( ia ) = fy_coul ( ia ) - fyij
              fz_coul ( ia ) = fz_coul ( ia ) - fzij

              ! virial
              vir_coul_dd =  vir_coul_dd + ( fxij * rxij + fyij * ryij + fzij * rzij )
 
              ! electric field gradient
              efg ( ia , 1 , 1 ) = efg ( ia , 1 , 1 ) - ( Txxx * mujx + Txxy * mujy + Txxz * mujz ) 
              efg ( ia , 2 , 2 ) = efg ( ia , 2 , 2 ) - ( Tyyx * mujx + Tyyy * mujy + Tyyz * mujz ) 
              efg ( ia , 3 , 3 ) = efg ( ia , 3 , 3 ) - ( Tzzx * mujx + Tzzy * mujy + Tzzz * mujz ) 
              efg ( ia , 1 , 2 ) = efg ( ia , 1 , 2 ) - ( Txxy * mujx + Tyyx * mujy + Txyz * mujz ) 
              efg ( ia , 1 , 3 ) = efg ( ia , 1 , 3 ) - ( Txxz * mujx + Txyz * mujy + Tzzx * mujz ) 
              efg ( ia , 2 , 3 ) = efg ( ia , 2 , 3 ) - ( Txyz * mujx + Tyyz * mujy + Tzzy * mujz ) 


              ! ===========================================================
              !                  charge-dipole interaction
              ! ===========================================================

              ! electrostatic energy
              u_coul_qd = u_coul_qd + ( qi * ( Tx * mujx + Ty * mujy + Tz * mujz ) - &
                                        qj * ( Tx * muix + Ty * muiy + Tz * muiz ) ) 

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


            endif



        endif

      enddo
   

    enddo cells


  enddo atom 


  ttt2 = MPI_WTIME(ierr)


  CALL MPI_ALL_REDUCE_DOUBLE_SCALAR ( u_coul_qq )
  CALL MPI_ALL_REDUCE_DOUBLE_SCALAR ( u_coul_dd )
  CALL MPI_ALL_REDUCE_DOUBLE_SCALAR ( u_coul_qd )
  CALL MPI_ALL_REDUCE_DOUBLE_SCALAR ( vir_coul_qq )
  CALL MPI_ALL_REDUCE_DOUBLE_SCALAR ( vir_coul_qd )
  CALL MPI_ALL_REDUCE_DOUBLE_SCALAR ( vir_coul_dd )

  CALL MPI_ALL_REDUCE_DOUBLE ( fx_coul , natm )
  CALL MPI_ALL_REDUCE_DOUBLE ( fy_coul , natm )
  CALL MPI_ALL_REDUCE_DOUBLE ( fz_coul , natm )

  CALL MPI_ALL_REDUCE_DOUBLE ( phi_coul , natm )

  CALL MPI_ALL_REDUCE_DOUBLE ( ef ( : , 1 )  , natm )
  CALL MPI_ALL_REDUCE_DOUBLE ( ef ( : , 2 )  , natm )
  CALL MPI_ALL_REDUCE_DOUBLE ( ef ( : , 3 )  , natm )

  CALL MPI_ALL_REDUCE_DOUBLE ( efg ( : , 1 , 1 ) , natm )
  CALL MPI_ALL_REDUCE_DOUBLE ( efg ( : , 2 , 2 ) , natm )
  CALL MPI_ALL_REDUCE_DOUBLE ( efg ( : , 3 , 3 ) , natm )
  CALL MPI_ALL_REDUCE_DOUBLE ( efg ( : , 1 , 2 ) , natm )
  CALL MPI_ALL_REDUCE_DOUBLE ( efg ( : , 1 , 3 ) , natm )
  CALL MPI_ALL_REDUCE_DOUBLE ( efg ( : , 2 , 3 ) , natm )

  ! EFG is symmetric
  ! not needed ... just for consistency
  efg ( : , 2 , 1 ) = efg ( : , 1 , 2 )
  efg ( : , 3 , 1 ) = efg ( : , 1 , 3 )
  efg ( : , 3 , 2 ) = efg ( : , 2 , 3 )


  u_coul_qq   = u_coul_qq * 0.5d0
  u_coul_dd   = u_coul_dd * 0.5d0
  u_coul_qd   = u_coul_qd * 0.5d0
  vir_coul_qq = - vir_coul_qq  * 0.5d0
  vir_coul_qd = - vir_coul_qd  * 0.5d0
  vir_coul_dd = - vir_coul_dd  * 0.5d0

  fx = fx + fx_coul
  fy = fy + fy_coul
  fz = fz + fz_coul

  u_coul   = u_coul_qq   + u_coul_dd   + u_coul_qd + u_pol
  vir_coul = vir_coul_qq + vir_coul_dd + vir_coul_qd
  phi_coul = phi_coul_qq + phi_coul_dd


#ifdef debug3 
  if ( ionode ) then
    WRITE ( stdout , '(a)' ) ' '
    WRITE ( stdout , '(a)' )     'Electric field at atoms : '
    do ia = 1 , natm
      WRITE ( stdout , '(i5,a3,a,3f18.10)' ) &
      ia,atype(ia),' Efield  = ', ef ( ia , 1)  , ef ( ia , 2 ) , ef ( ia , 3 )
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
      ia,atype(ia),' phi_tot = ', phi_coul ( ia )  , ' phi_qq = ', phi_coul_qq (ia) , ' phi_dd = ', phi_coul_dd (ia)
    enddo
    WRITE ( stdout , '(a)' ) ' '
    WRITE ( stdout , '(a)' ) 'Energy and virial : '
    WRITE ( stdout , '(5(a,f18.10))' ) &
    ' u_coul_tot   = ', u_coul  ,' u_coul_qq   = ',u_coul_qq      ,' u_coul_qd   = ', u_coul_qd     ,' u_coul_dd   = ',u_coul_dd   ,' u_pol = ', u_pol
    WRITE ( stdout , '(4(a,f18.10))' ) &
    ' vir_coul_tot = ',vir_coul ,' vir_coul_qq = ',vir_coul_qq    ,' vir_coul_qd = ', vir_coul_qd   ,' vir_coul_dd = ',vir_coul_dd
  endif
#endif


  deallocate( fx_coul , fy_coul , fz_coul )
  deallocate( phi_coul_qq , phi_coul_dd  )

  ttt3 = MPI_WTIME(ierr)




  return

END SUBROUTINE



SUBROUTINE multipole_ES ( iastart, iaend , ef , efg , mu , u_coul , vir_coul , phi_coul )

  USE config,    ONLY : natm , ntype , natmi , rx , ry , rz , fx , fy , fz , qia , omega
  USE constants, ONLY : imag , piroot , tpi , fpi
  USE io_file  , ONLY : stdout 
  USE time

  implicit none

  INCLUDE 'mpif.h'

  ! global 
  integer, intent(in) :: iastart , iaend
  double precision    :: ef     ( natm , 3 )
  double precision    :: efg    ( natm , 3 , 3 )
  double precision    :: mu     ( natm , 3 )
  double precision    :: u_coul 
  double precision    :: vir_coul 
  double precision    :: phi_coul ( natm )  

  ! local 
  integer          :: ia , ja , ik , it , ip , ierr
  double precision, dimension(:)    , allocatable :: fx_coul , fy_coul , fz_coul
  double precision, dimension(:)    , allocatable :: fx_dir , fy_dir , fz_dir
  double precision, dimension(:)    , allocatable :: fx_rec , fy_rec , fz_rec
  double precision, dimension(:)    , allocatable :: fx_surf , fy_surf , fz_surf
  double precision, dimension(:)    , allocatable :: phi_dir , phi_rec , phi_surf , phi_self
  double precision, dimension(:,:)  , allocatable :: ef_dir , ef_rec , ef_surf
  double precision, dimension(:,:,:), allocatable :: efg_dir , efg_rec , efg_self

  double precision :: u_dir , u_rec , u_surf , u_self
  double precision :: vir_dir , vir_rec , vir_surf

  double precision :: mip    ( ntype , 3 )

  double precision :: rxi , ryi , rzi
  double precision :: rxj , ryj , rzj
  double precision :: kx , ky , kz 
  double precision :: kk , kri , Ak
  double precision :: rxij , ryij , rzij
  double precision :: fxij , fyij , fzij
  double precision :: qi , qj , qij
  double precision :: muix , muiy , muiz
  double precision :: mujx , mujy , mujz
  double complex   :: rhon , carg
  double precision :: str , recarg , recarg_dgg , recargi
  double precision :: expon , F0 , F1 , F2 , Fdamp
  double precision :: T
  double precision :: Tx , Ty , Tz
  double precision :: Txx , Tyy , Tzz , Txy , Txz , Tyz
  double precision :: Dxx , Dyy , Dzz , Dxy , Dxz , Dyz
  double precision :: Txxx,  Tyyy,  Tzzz, Txxy, Txxz, Tyyx, Tyyz, Tzzx, Tzzy, Txyz
  double precision :: d , d2 , d3  
  double precision :: dm1 , dm3 , dm5 , dm7 
  double precision :: alpha2 , alpha3 , selfa , selfa3 
  double precision, external :: errfc
  double precision :: qtot ( 3 ) , qsq , mutot ( 3 ) , musq , qmu_sum ( 3 ) 

  double precision :: wij , old !tmp


  double precision :: ttt1 , ttt2  , ttt3

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
  mutot = 0.0d0
  musq  = 0.0d0
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
  qmu_sum = qtot + mutot



#ifdef debug
  write( stdout , '(a)')    'debug: in engforce_charge_ES'
  write( stdout , '(a,i8,a,a)') 'debug : km_coul ',km_coul%nkcut,' ',km_coul%meshlabel
  do ia = 1 , natm
    write( stdout , '(a,f12.5)') 'debug : charge (atom) ',qia(ia)
  enddo
  do it = 1 , ntype
    write( stdout , '(a,f12.5)') 'debug : charge (type ) ',qch(it)
  enddo
  do ia = 1 , natm
    write( stdout , '(a,3f12.5)') 'debug : dipole (atom) ', mu ( ia , 1 ) , mu ( ia , 2 ) , mu ( ia , 3 )
  enddo
  do it = 1 , ntype
    write( stdout , '(a,3f12.5)') 'debug : dipole (type ) ',  mip ( it , 1 ) , mip ( it , 2 ) , mip ( it , 3 )
  enddo
  write( stdout , '(a,2i8)')     'debug : iastart iaend',iastart ,iaend
  write( stdout , '(a,f20.5)')   'debug : alphaES ', alphaES
  call print_config_sample(0,0)
#endif 

  ttt1 = MPI_WTIME(ierr)
  allocate( ef_dir ( natm , 3 ) , ef_rec ( natm , 3 ) , ef_surf ( natm , 3 ) ) 
  allocate( efg_dir ( natm , 3 , 3 ) , efg_rec ( natm , 3 , 3 ) , efg_self ( natm , 3 , 3 ) ) 
  allocate( fx_coul (natm) , fy_coul (natm) , fz_coul (natm) )
  allocate( fx_dir (natm) , fy_dir (natm) , fz_dir (natm) )
  allocate( fx_rec (natm) , fy_rec (natm) , fz_rec (natm) )
  allocate( fx_surf (natm) , fy_surf (natm) , fz_surf (natm) )
  allocate( phi_dir ( natm ) , phi_rec ( natm ) , phi_surf ( natm ) , phi_self ( natm ) )

  ef_dir   = 0.0d0
  ef_rec   = 0.0d0
  ef_surf  = 0.0d0
  efg_dir  = 0.0d0
  efg_rec  = 0.0d0
  efg_self = 0.0d0
  fx_dir   = 0.0D0
  fy_dir   = 0.0D0
  fz_dir   = 0.0D0
  fx_rec   = 0.0D0
  fy_rec   = 0.0D0
  fz_rec   = 0.0D0
  fx_surf  = 0.0D0
  fy_surf  = 0.0D0
  fz_surf  = 0.0D0
  phi_dir  = 0.0D0
  phi_rec  = 0.0D0
  phi_surf = 0.0D0
  phi_self = 0.0D0

  ! =====================
  ! facteur de structure 
  ! =====================
  CALL struc_fact ( km_coul )

  ! =================
  !  some constants 
  ! =================
  !invbox = 1.0d0 / box
  alpha2 = alphaES * alphaES
  alpha3 = alpha2  * alphaES
  selfa  = alphaES / piroot
  !selfa3 = - 4.0d0 * alpha3 / piroot / 3.0d0 
  selfa3 = - fpi / 3.0d0 / omega  

  ! ==============================================
  !        direct space part
  ! ==============================================

  do ia = 1 , natm
    rxi = rx(ia)
    ryi = ry(ia)
    rzi = rz(ia)
    qi  = qia(ia)
    do ja = 1, natm

      if (ja .ne. ia ) then

        qj   = qia(ja)
        qij  = qi * qj 
        rxj  = rx(ja)
        ryj  = ry(ja)
        rzj  = rz(ja)
        rxij = rxi - rxj
        ryij = ryi - ryj
        rzij = rzi - rzj

!      nxij = NINT ( rxij * invbox )
!      nyij = NINT ( ryij * invbox )
!      nzij = NINT ( rzij * invbox )
!      rxij = rxij - box * nxij
!      ryij = ryij - box * nyij
!      rzij = rzij - box * nzij

       d2  = rxij * rxij + ryij * ryij + rzij * rzij
       d   = SQRT ( d2 )
       d3  = d2 * d 
       dm1 = 1.0d0 / d
       dm3 = dm1 / d2 
       dm5 = dm3 / d2 
       dm7 = dm5 / d2 

       !pot  = qj  / rij
       !qij  = qi  * pot
       !pot3 = pot / rijsq
       !qijf = qi  * pot3         ! qi * qj / r^3

       expon = EXP ( - alpha2 * d2 )    / piroot
       F0    = errfc( alphaES * d )
       F1    = F0 + 2.0d0 * d * alphaES * expon 
       F2    = 4.0d0          * alpha3  * expon * d3 / 3.0d0
       Fdamp = F1 + F2 
       ! multipole interaction tensor rank = 0 
       T  = dm1 * F0

       ! multipole interaction tensor rank = 1
       Tx = rxij * dm3 * F1
       Ty = ryij * dm3 * F1
       Tz = rzij * dm3 * F1

       ! multipole interaction tensor rank = 2
       Txx = ( 3.0d0 * rxij * rxij - d2 ) * dm5 
       Tyy = ( 3.0d0 * ryij * ryij - d2 ) * dm5 
       Tzz = ( 3.0d0 * rzij * rzij - d2 ) * dm5
       Txy = ( 3.0d0 * rxij * ryij      ) * dm5
       Txz = ( 3.0d0 * rxij * rzij      ) * dm5 
       Tyz = ( 3.0d0 * ryij * rzij      ) * dm5 

       ! multipole interaction tensor rank = 3  
       Txxx = ( 5.0d0 * rxij * rxij * rxij -  3.0d0 * d2 * ( rxij ) ) * dm7 * 3.0d0
       Tyyy = ( 5.0d0 * ryij * ryij * ryij -  3.0d0 * d2 * ( ryij ) ) * dm7 * 3.0d0
       Tzzz = ( 5.0d0 * rzij * rzij * rzij -  3.0d0 * d2 * ( rzij ) ) * dm7 * 3.0d0
       Txxy = ( 5.0d0 * rxij * rxij * ryij -          d2 * ( ryij ) ) * dm7 * 3.0d0
       Txxz = ( 5.0d0 * rxij * rxij * rzij -          d2 * ( rzij ) ) * dm7 * 3.0d0
       Tyyx = ( 5.0d0 * ryij * ryij * rxij -          d2 * ( rxij ) ) * dm7 * 3.0d0
       Tyyz = ( 5.0d0 * ryij * ryij * rzij -          d2 * ( rzij ) ) * dm7 * 3.0d0
       Tzzx = ( 5.0d0 * rzij * rzij * rxij -          d2 * ( rxij ) ) * dm7 * 3.0d0
       Tzzy = ( 5.0d0 * rzij * rzij * ryij -          d2 * ( ryij ) ) * dm7 * 3.0d0
       Txyz = ( 5.0d0 * rxij * ryij * rzij                          ) * dm7 * 3.0d0


       ! ===========================================================
       !                  charge-charge interaction
       ! ===========================================================

       phi_dir ( ia ) = phi_dir ( ia ) + qj * T 

       u_dir = u_dir + qij * T

       ef_dir ( ia , 1 ) = ef_dir ( ia , 1 ) + qj * Tx 
       ef_dir ( ia , 2 ) = ef_dir ( ia , 2 ) + qj * Ty 
       ef_dir ( ia , 3 ) = ef_dir ( ia , 3 ) + qj * Tz 

       fxij = qij * Tx
       fyij = qij * Ty
       fzij = qij * Tz 

       fx_dir ( ia ) = fx_dir ( ia ) + fxij
       fy_dir ( ia ) = fy_dir ( ia ) + fyij
       fz_dir ( ia ) = fz_dir ( ia ) + fzij

       vir_dir =  vir_dir + ( fxij * rxij + fyij * ryij + fzij * rzij )

       efg_dir ( ia , 1 , 1 ) = efg_dir ( ia , 1 , 1 ) - qj * Txx * Fdamp
       efg_dir ( ia , 2 , 2 ) = efg_dir ( ia , 2 , 2 ) - qj * Tyy * Fdamp 
       efg_dir ( ia , 3 , 3 ) = efg_dir ( ia , 3 , 3 ) - qj * Tzz * Fdamp
       efg_dir ( ia , 1 , 2 ) = efg_dir ( ia , 1 , 2 ) - qj * Txy * Fdamp 
       efg_dir ( ia , 1 , 3 ) = efg_dir ( ia , 1 , 3 ) - qj * Txz * Fdamp 
       efg_dir ( ia , 2 , 3 ) = efg_dir ( ia , 2 , 3 ) - qj * Tyz * Fdamp


      endif

    enddo

  enddo 

  ttt2 = MPI_WTIME(ierr)
  fcoultimetot1 = fcoultimetot1 + ( ttt2 - ttt1 )

  ! ==============================================
  !            reciprocal space part
  ! ==============================================
  kpoint : do ik = 1, km_coul%nkcut
    ! =================
    !   k-space  
    ! =================
    kx = km_coul%kpt(1,ik)
    ky = km_coul%kpt(2,ik)
    kz = km_coul%kpt(3,ik)
    kk = km_coul%kptk(ik)
    Ak = EXP ( - kk * 0.25d0 / alpha2 ) / kk


    if (km_coul%kptk(ik) .eq. 0 ) then
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
      rhon = rhon + qch(it) * CONJG( km_coul%strf ( ik , it ) )
    enddo

    str = rhon * CONJG(rhon)

    u_rec   = u_rec   + DBLE ( Ak * str )
    vir_rec = vir_rec + DBLE ( Ak * str * 0.5d0 * ( 1.0d0 - ( kk * 0.5d0 / alpha2 ) )  )

    do ia = 1 , natm

      rxi = rx(ia)
      ryi = ry(ia)
      rzi = rz(ia)
      qi  = qia ( ia )
      kri = ( kx * rxi + ky * ryi + kz * rzi )

      carg       = EXP ( imag * kri )
      recarg     = DBLE ( rhon * carg * Ak )
      recarg_dgg = DBLE ( rhon * carg * Ak * kk / 3.0d0 )  
      recargi    = DBLE ( rhon * carg * Ak * imag )
      
      phi_rec ( ia ) = phi_rec ( ia ) + recarg 

      fxij = kx * recargi 
      fyij = ky * recargi 
      fzij = kz * recargi 

      ef_rec ( ia , 1 ) = ef_rec ( ia , 1 ) - fxij
      ef_rec ( ia , 2 ) = ef_rec ( ia , 2 ) - fyij
      ef_rec ( ia , 3 ) = ef_rec ( ia , 3 ) - fzij


      fx_rec ( ia ) = fx_rec ( ia ) - qi * fxij
      fy_rec ( ia ) = fy_rec ( ia ) - qi * fyij
      fz_rec ( ia ) = fz_rec ( ia ) - qi * fzij

      efg_rec ( ia , 1 , 1 ) = efg_rec ( ia , 1 , 1 ) +  kx * kx * recarg - recarg_dgg 
      efg_rec ( ia , 2 , 2 ) = efg_rec ( ia , 2 , 2 ) +  ky * ky * recarg - recarg_dgg
      efg_rec ( ia , 3 , 3 ) = efg_rec ( ia , 3 , 3 ) +  kz * kz * recarg - recarg_dgg
      efg_rec ( ia , 1 , 2 ) = efg_rec ( ia , 1 , 2 ) +  kx * ky * recarg 
      efg_rec ( ia , 1 , 3 ) = efg_rec ( ia , 1 , 3 ) +  kx * kz * recarg 
      efg_rec ( ia , 2 , 3 ) = efg_rec ( ia , 2 , 3 ) +  ky * kz * recarg 

    enddo


  enddo kpoint

  ttt3 = MPI_WTIME(ierr)
  fcoultimetot2 = fcoultimetot2 + ( ttt3 - ttt2 )

  ! ======================================================
  ! remark on the unit :
  ! 1/(4*pi*epislon_0) = 1 => epsilon_0 = 1/4pi
  ! ======================================================
  u_dir   =   u_dir   * 0.5d0
  vir_dir = - vir_dir * 0.5d0

  u_rec   =   u_rec   * tpi / omega
  vir_rec = - vir_rec * fpi / omega
  phi_rec =   phi_rec * fpi / omega

  
  ! surface_contribution 
  ! qq
  u_surf = qtot ( 1 ) * qtot ( 1 ) + qtot ( 2 ) * qtot ( 2 ) + qtot ( 3 ) * qtot ( 3 )
  ! dd
  u_surf = u_surf +  mutot ( 1 ) * mutot ( 1 ) + mutot ( 2 ) * mutot ( 2 ) + mutot ( 3 ) * mutot ( 3 )
  ! qd
  u_surf = u_surf + 2.0d0 * ( qtot ( 1 ) * mutot ( 1 ) + qtot ( 2 ) * mutot ( 2 ) + qtot ( 3 ) * mutot ( 3 ) )
  u_surf = u_surf * tpi / 3.0d0 / omega
  vir_surf = - u_surf
  do ia = 1 , natm
    phi_surf ( ia ) = rx ( ia ) * qmu_sum ( 1 ) + ry ( ia ) * qmu_sum ( 1 ) + rz ( ia ) * qmu_sum ( 1 )
    ef_surf ( ia , 1 ) = qmu_sum ( 1 )
    ef_surf ( ia , 2 ) = qmu_sum ( 2 )
    ef_surf ( ia , 3 ) = qmu_sum ( 3 )
    fx_surf ( ia ) = qia ( ia ) * qmu_sum ( 1 )
    fy_surf ( ia ) = qia ( ia ) * qmu_sum ( 2 )
    fz_surf ( ia ) = qia ( ia ) * qmu_sum ( 3 )
  enddo

  phi_surf = phi_surf * fpi / omega / 3.0d0

  ! self
  u_self  = - selfa * qsq
  do ia = 1 , natm
   phi_self ( ia ) = - 2.0d0 * selfa * qia ( ia )
  ! efg_self ( ia , 1 , 1 ) = selfa3 * qia ( ia ) 
  ! efg_self ( ia , 2 , 2 ) = selfa3 * qia ( ia ) 
  ! efg_self ( ia , 3 , 3 ) = selfa3 * qia ( ia ) 
  enddo

  u_coul   = ( u_dir + u_rec + u_surf + u_self )
  vir_coul = ( vir_dir + vir_rec + vir_surf )
  phi_coul = ( phi_dir + phi_rec + phi_surf + phi_self )

#ifdef debug
  write( stdout , '(4(a,f16.8))' ) ,'vir_dir    = ', vir_dir    ,' vir_rec    = ', vir_rec    ,' vir_surf    = ',vir_surf   ,' vir_coul    = '   ,vir_coul
  write( stdout , '(5(a,f16.8))' ) ,'u_dir      = ', u_dir      ,' u_rec      = ', u_rec      ,' u_surf      = ',u_surf     ,' u_self      = '   ,u_self     ,' u_coul      = ',u_coul
  write( stdout , '(5(a,f16.8))' ) ,'phi_dir(1) = ', phi_dir(1) ,' phi_rec(1) = ', phi_rec(1) ,' phi_surf(1) = ',phi_surf(1),' phi_self(1) = '   ,phi_self(1),' phi_coul(1) = ',phi_coul (1)
#endif

  ef_rec  = ef_rec  * fpi / omega
  ef_surf = ef_surf * fpi / omega / 3.0d0
  ef  = ef_dir + ef_rec - ef_surf

  fx_rec = fx_rec   * fpi / omega
  fy_rec = fy_rec   * fpi / omega
  fz_rec = fz_rec   * fpi / omega

  fx_surf = fx_surf * fpi / omega / 3.0d0
  fy_surf = fy_surf * fpi / omega / 3.0d0
  fz_surf = fz_surf * fpi / omega / 3.0d0

  fx = fx + fx_rec + fx_dir - fx_surf
  fy = fy + fy_rec + fy_dir - fy_surf
  fz = fz + fz_rec + fz_dir - fz_surf


  efg_rec = efg_rec * fpi / omega 

  efg = efg_dir + efg_rec + efg_self


#ifdef debug
  do ia=1,natm
    WRITE ( stdout ,'(a,i4,a,3f12.8)') 'atom = ',ia,' Efield_dir         = ', ef_dir ( ia , 1 ) , ef_dir ( ia , 2 ) , ef_dir ( ia , 3 )
  enddo
  do ia=1,natm
    WRITE ( stdout ,'(a,i4,a,3f12.8)') 'atom = ',ia,' Efield_rec         = ', ef_rec ( ia , 1 ) , ef_rec ( ia , 2 ) , ef_rec ( ia , 3 )
  enddo
  do ia=1,natm
    WRITE ( stdout ,'(a,i4,a,3f12.8)') 'atom = ',ia,' Efield_surf        = ', ef_surf ( ia , 1 ) , ef_surf ( ia , 2 ) , ef_surf ( ia , 3 )
  enddo
  do ia=1,natm
    WRITE ( stdout ,'(a,i4,a,3f12.8)') 'atom = ',ia,' Efield_tot         = ', ef ( ia , 1 ) , ef ( ia , 2 ) , ef ( ia , 3 )
  enddo
    WRITE ( stdout ,'(a)') '' 
  do ia=1,natm
    WRITE ( stdout ,'(a,i4,a,3f12.8)') 'atom = ',ia,' f_dir         = ', fx_dir ( ia ) , fy_dir ( ia ) , fz_dir ( ia )
  enddo
  do ia=1,natm
    WRITE ( stdout ,'(a,i4,a,3f12.8)') 'atom = ',ia,' f_rec         = ', fx_rec ( ia ) , fy_rec ( ia ) , fz_rec ( ia )
  enddo
  do ia=1,natm
    WRITE ( stdout ,'(a,i4,a,3f12.8)') 'atom = ',ia,' f_surf        = ', fx_surf ( ia ) , fy_surf ( ia ) , fz_surf ( ia )
  enddo
  do ia=1,natm
    WRITE ( stdout ,'(a,i4,a,3f12.8)') 'atom = ',ia,' f_tot         = ', fx ( ia ) , fy ( ia ) , fz ( ia )
  enddo

CALL print_tensor( efg_dir ( 1 , : , : ) , 'DIRECT' )
CALL print_tensor( efg_rec ( 1 , : , : ) , 'RECIPR' )
CALL print_tensor( efg     ( 1 , : , : ) , 'TOTAL ' )

#endif


  deallocate( ef_dir  , ef_rec , ef_surf ) 
  deallocate( efg_dir , efg_rec , efg_self ) 
  deallocate( fx_coul , fy_coul , fz_coul )
  deallocate( fx_dir  , fy_dir  , fz_dir )
  deallocate( fx_rec  , fy_rec  , fz_rec  )
  deallocate( fx_surf , fy_surf , fz_surf )
  deallocate( phi_dir , phi_rec , phi_surf , phi_self )

  return

END SUBROUTINE







END MODULE field 
! ===== fmV =====
