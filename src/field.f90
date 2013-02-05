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
!#define debug
!#define debug2
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

  ! Lennard - Jones
  double precision :: qlj     ( ntypemax , ntypemax )
  double precision :: plj     ( ntypemax , ntypemax )
  double precision :: epslj   ( ntypemax , ntypemax )
  double precision :: sigmalj ( ntypemax , ntypemax )
  
  ! Morse
  double precision :: rhomor  ( ntypemax , ntypemax )
  double precision :: epsmor  ( ntypemax , ntypemax )
  double precision :: sigmamor( ntypemax , ntypemax )

  double precision :: mass  ( ntypemax )    ! masses ( not yet )
  double precision :: qch   ( ntypemax )    ! charges only for efg (no electrostatic interaction and forces)
  double precision :: dip   ( ntypemax , 3 )! dipole for efg (no electrostatic interaction and forces)
  double precision :: pol   ( npolmax  , 3 , 3 )! polarizability ( is defined at particular position )
  
  double precision ::  conv_tol_ind

  double precision :: utail

  !todo : make the following quantities dependent on ntypemax ( hardware parameter) 
  ! or find a way to read them directly from control.F 

  ! lennard-jones
  double precision :: rcutsq ( ntypemax , ntypemax )  
  double precision :: sigsq  ( ntypemax , ntypemax ) 
  double precision :: epsp   ( ntypemax , ntypemax )
  double precision :: fc     ( ntypemax , ntypemax )
  double precision :: uc     ( ntypemax , ntypemax )
  double precision :: uc1    ( ntypemax , ntypemax )
  double precision :: uc2    ( ntypemax , ntypemax )
  double precision :: testtab ( ntypemax , ntypemax )

  !morse
  double precision :: rs     ( ntypemax , ntypemax )
  double precision :: fm     ( ntypemax , ntypemax )

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
  epsmor        = 1.0d0 
  sigmamor      = 1.0d0 
  rhomor        = 3.0d0 
  qlj           = 12.0d0
  plj           = 6.0d0
  mass          = 1.0d0
  ! Coulomb
  ncelldirect   =  2
  ncellewald    = 10
  alphaES       =  1.0d0
  qch           =  0.0d0
  dip           =  0.0d0
  pol           =  0.0d0
  lpolar        =  0

  conv_tol_ind  = 1e-5

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

  USE control,  ONLY :  calc , lbmlj , lcoulomb , lmorse
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
                         sigmamor      , &
                         epsmor        , &
                         rhomor        , &
                         mass          , &
                         qch           , &
                         dip           , &
                         pol           , &  
                         conv_tol_ind  , &  
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
  if ( lbmlj .or. lmorse )    then
    CALL initialize_param_bmlj_morse
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

  USE config,   ONLY :  ntype , box, atypei , natmi
  USE control,  ONLY :  calc , cutoff , lbmlj , lcoulomb , longrange
  USE io_file,  ONLY :  ionode 
  USE constants,ONLY :  tpi

  implicit none

  !local
  integer :: kunit, it , it1 , it2
  double precision :: aaa , rcut2 , kmax2 , alpha2 , ereal , ereci , qtot
  logical :: linduced

  qtot = 0.0d0
  kmax2 = tpi / box * ncellewald
  kmax2 = kmax2 * kmax2 
  rcut2 = box * box * 0.25d0
  alpha2 = alphaES * alphaES
  ereal = 2.0D0 * EXP ( - alpha2 * rcut2 ) / box
  ereci = EXP ( - kmax2 / alpha2 * 0.25d0 ) / kmax2

  if ( ionode ) then
    WRITE ( kunit ,'(a)')       '=============================================================' 
    WRITE ( kunit ,'(a)')       ''
    if ( .not. lcoulomb .and. calc .eq. 'efg' )  &
    WRITE ( kunit ,'(a)')       'point charges: (no forces yet !)'
    do it = 1 , ntype 
      WRITE ( kunit ,'(a,a,a,f10.5)') 'q',atypei(it),'      =',qch(it)
      qtot = qtot + qch(it) * natmi ( it ) 
    enddo
    WRITE ( kunit ,'(a,f10.5)')       'total charge of the system = ',  qtot
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
        WRITE ( kunit ,'(a,e12.5)') 'relative error in real space with alphaES from input : ',ereal
        WRITE ( kunit ,'(a,e12.5)') 'relative error in reciprocal space with alphaES from input : ',ereci
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
            if ( it1 .eq. it2 ) then
              WRITE ( kunit ,'(a,f10.5)') 'sigma                                = ',sigmalj ( it1 , it2 )
              WRITE ( kunit ,'(a,f10.5)') 'eps                                  = ',epslj   ( it1 , it2 )
              WRITE ( kunit ,'(a,f10.5)') 'q                                    = ',qlj     ( it1 , it2 )
              WRITE ( kunit ,'(a,f10.5)') 'p                                    = ',plj     ( it1 , it2 )
            else
              WRITE ( kunit ,'(a,f10.5,a,f10.5,a)') 'sigma                                = ',sigmalj ( it1 , it2 ) , '( ',sigmalj ( it2 , it1 ), ' )'
              WRITE ( kunit ,'(a,f10.5,a,f10.5,a)') 'eps                                  = ',epslj   ( it1 , it2 ) , '( ',epslj   ( it2 , it1 ), ' )'
              WRITE ( kunit ,'(a,f10.5,a,f10.5,a)') 'q                                    = ',qlj     ( it1 , it2 ) , '( ',qlj     ( it2 , it1 ), ' )'
              WRITE ( kunit ,'(a,f10.5,a,f10.5,a)') 'p                                    = ',plj     ( it1 , it2 ) , '( ',plj     ( it2 , it1 ), ' )'
            endif
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
          WRITE ( kunit ,'(a,a2,a)') 'polarizability tensor of type ', atypei(it),' : ' 
          WRITE ( kunit ,'(a)')   ' '
          WRITE ( kunit ,'(3f12.5)')  pol( it , 1 , 1 ) , pol( it , 1 , 2 ) , pol( it , 1 , 3 ) 
          WRITE ( kunit ,'(3f12.5)')  pol( it , 2 , 1 ) , pol( it , 2 , 2 ) , pol( it , 2 , 3 ) 
          WRITE ( kunit ,'(3f12.5)')  pol( it , 3 , 1 ) , pol( it , 3 , 2 ) , pol( it , 3 , 3 ) 
          WRITE ( kunit ,'(a)')   ' '
        else
          WRITE ( kunit ,'(a,a2)') 'no polarizability on type ', atypei(it)
        endif
        WRITE ( kunit ,'(a)')   ' '
      enddo
    endif

  endif

  return 

END SUBROUTINE field_print_info

!*********************** SUBROUTINE initialize_param_bmlj_morse ***************
!
! iniitialisation of BMLJ parameters
!
!******************************************************************************

SUBROUTINE initialize_param_bmlj_morse 

  USE constants,        ONLY :  pi
  USE config,           ONLY :  ntypemax , natm , box , omega , natmi , rho , atype , itype  , ntype
  USE control,          ONLY :  skindiff , cutoff , lreduced, calc
  USE io_file,          ONLY :  ionode, stdout, kunit_OUTFF

  implicit none

  ! local
  integer :: it,jt
  double precision :: one13, one16, two16, rskinmax, utail
  double precision :: rcut3 ( ntypemax , ntypemax )
  double precision :: rskin ( ntypemax , ntypemax ) 
  double precision :: rskinsq ( ntypemax , ntypemax ) 
  double precision :: ut ( ntypemax , ntypemax ) 

  double precision :: rcut ( ntype , ntype ) 
  double precision :: ppqq ( ntype , ntype )
  double precision :: pp ( ntype , ntype )
  double precision :: qq ( ntype , ntype )
  double precision :: pp3 ( ntype , ntype ) 
  double precision :: qq3 ( ntype , ntype ) 
  double precision :: sr2 ( ntype , ntype ) 
  double precision :: sr ( ntype , ntype ) 
  double precision :: srp ( ntype , ntype ) 
  double precision :: srq ( ntype , ntype ) 

  rskinmax = 0.0D0
  utail    = 0.0d0
  rcut3    = 0.0d0
  rcut     = 0.0d0
  rskin    = 0.0d0
  rskinsq  = 0.0d0
  ut       = 0.0d0
  ppqq     = 0.0d0
  pp       = 0.0d0
  qq       = 0.0d0
  sr2      = 0.0d0
  sr       = 0.0d0
  srp      = 0.0d0
  srq      = 0.0d0
  
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
  do it = 1 , ntype 
    do jt = 1 , ntype

      ! intermediate terms 
          rcut   ( it , jt ) = cutoff                                      ! rc
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
           uc1 ( it , jt ) = epsp ( it , jt ) *  ppqq ( it , jt ) / ( 2.0d0 * rcutsq ( it , jt ) ) * ( srq( it , jt ) - srp( it , jt ) ) 
           ! c2
           uc2 ( it , jt ) = 0.5d0 * epsp( it , jt )  * (  &
                         ( 2.0d0 * qq ( it , jt ) + ppqq ( it , jt )  ) * srq( it , jt ) - &
                         ( 2.0d0 * pp ( it , jt ) + ppqq ( it , jt )  ) * srp( it , jt ) )
           ! for the virial
           fc ( it , jt ) =  ppqq ( it , jt ) * epsp ( it , jt ) /  sigsq ( it , jt ) 
           ! morse
           fm ( it , jt ) = - 2.0d0 * epsmor ( it , jt ) * rhomor ( it , jt ) * EXP ( rhomor ( it , jt ) * sigmamor ( it , jt ) )  
           rs ( it , jt ) = EXP ( rhomor ( it , jt ) * sigmamor( it , jt ) ) 
           ! tail energy
           ut ( it , jt ) = epsp ( it , jt ) * ( pp ( it , jt ) * srq ( it , jt ) / qq3 ( it , jt ) - &
                                                 qq ( it , jt ) * srp ( it , jt ) / pp3 ( it , jt )  )       
           ut ( it , jt ) = ut ( it , jt ) * rcut3 ( it , jt ) * 2.0d0 * pi 
           if ( ( natmi ( it ) .ne. 0 ) .and. ( natmi ( jt ) .ne. 0 ) ) &
           utail = utail + ut ( it , jt ) * natmi ( it ) * natmi ( jt ) / omega
#ifdef debug 
  write ( * , '(2i6,7e16.6)' ) it , jt , uc ( it , jt )  , epsp ( it , jt ) , pp ( it , jt ) , qq ( it , jt ) , srq( it , jt ) , srp( it , jt ) ,rcutsq ( it , jt ) 
#endif
    enddo
  enddo
 
  if (ionode.and..not.lreduced) WRITE( stdout     ,'(a,2f20.9)') 'long range correction : ',utail
  if (ionode.and.lreduced)      WRITE( stdout     ,'(a,2f20.9)') 'long range correction : ',utail/natm
  if (ionode.and..not.lreduced) WRITE( kunit_OUTFF,'(a,2f20.9)') 'long range correction : ',utail
  if (ionode.and.lreduced)      WRITE( kunit_OUTFF,'(a,2f20.9)') 'long range correction : ',utail/natm
  if (ionode) WRITE( stdout     ,'(a)') ' '
  if (ionode) WRITE( kunit_OUTFF,'(a)') ' '

  return

END SUBROUTINE initialize_param_bmlj_morse


!*********************** SUBROUTINE engforce_driver ***************************
!
! this subroutine is used as an interfaced to the different potential + forces
! subroutines
!
!******************************************************************************

SUBROUTINE engforce_driver ( iastart , iaend )!, list , point )

  USE config,   ONLY :  natm , dipia
  USE control,  ONLY :  lpbc , lminimg , lbmlj , lcoulomb , lmorse , lshiftpot , longrange
  USE io_file,  ONLY :  stdout 

  implicit none

  ! global
  integer, intent(inout)  :: iastart , iaend 
  ! local 
  logical :: lj_and_coul
  double precision :: eftmp( natm , 3 ) , efgtmp ( natm , 3 , 3 ) , mutmp ( natm , 3 ) , u_coultmp , vir_coultmp , phi_coultmp ( natm ) 

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
    CALL engforce_morse_pbc         ( iastart , iaend )
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
      if ( longrange .eq. 'ewald'  ) CALL multipole_ES              ( iastart , iaend , eftmp , efgtmp , dipia , u_coultmp , vir_coultmp , phi_coultmp )
      if ( longrange .eq. 'direct' ) CALL multipole_DS              ( iastart , iaend , eftmp , efgtmp , dipia , u_coultmp , vir_coultmp , phi_coultmp )
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
      if ( longrange .eq. 'ewald'  ) CALL multipole_ES              ( iastart , iaend , eftmp , efgtmp , dipia , u_coultmp , vir_coultmp , phi_coultmp )
      if ( longrange .eq. 'direct' ) CALL multipole_DS              ( iastart , iaend , eftmp , efgtmp , dipia , u_coultmp , vir_coultmp , phi_coultmp )
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

  USE config,           ONLY : natm , box , rx , ry , rz , fx , fy , fz, tau_lj , atype , itype , list , point , ntype , omega
  USE control,          ONLY : lvnlist , myrank
  USE thermodynamic,    ONLY : u_lj , vir_lj , write_thermo
  USE time
  USE io_file,          ONLY : stdout

  implicit none
  INCLUDE 'mpif.h'

  ! global
  integer, intent(in)  :: iastart , iaend 

  ! local
  integer          :: ia , ja , it , jt , j1 , je , jb , ierr , i
  integer          :: p1 , p2
  integer          :: nxij , nyij , nzij
  double precision :: rxi , ryi , rzi
  double precision :: rxij , ryij , rzij , sr2 , rijsq , srp , srq
  double precision :: wij , fxij , fyij , fzij
  double precision :: ptwo ( ntype , ntype )
  double precision :: qtwo ( ntype , ntype )
  double precision :: forcetime1 , forcetime2 , invbox
  double precision :: u , vir 

#ifdef debug
  !do ia = 1 , natm
  ia = natm 
    write( stdout , '(a,i6,a,a)' )  'debug : atype ',ia,'',atype(ia)
    write( stdout , '(a,i6,a,i4)' ) 'debug : itype ',ia,'',itype(ia)
  !enddo
#endif

  forcetime1 = MPI_WTIME(ierr) ! timing info

  u   = 0.0D0
  vir = 0.0D0
  fx  = 0.0D0
  fy  = 0.0D0
  fz  = 0.0D0
 
  invbox = 1.0d0 / box

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
    ! =====================================
    !  if verlet list : ja in point arrays
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
        p1 = itype ( ia )
        p2 = itype ( ja )
        if ( rijsq .lt. rcutsq(p1,p2) ) then
          sr2 = sigsq(p1,p2) / rijsq
          srp = sr2 ** (ptwo(p1,p2))
          srq = sr2 ** (qtwo(p1,p2))
          ! potential energy ( truncated )  
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
          tau_lj(1,1) = tau_lj(1,1) + rxij * fxij
          tau_lj(1,2) = tau_lj(1,2) + rxij * fyij
          tau_lj(1,3) = tau_lj(1,3) + rxij * fzij
          tau_lj(2,1) = tau_lj(2,1) + ryij * fxij
          tau_lj(2,2) = tau_lj(2,2) + ryij * fyij
          tau_lj(2,3) = tau_lj(2,3) + ryij * fzij
          tau_lj(3,1) = tau_lj(3,1) + rzij * fxij
          tau_lj(3,2) = tau_lj(3,2) + rzij * fyij
          tau_lj(3,3) = tau_lj(3,3) + rzij * fzij
        endif
      endif
    enddo
  enddo
  tau_lj = tau_lj / omega
  vir = vir/3.0D0

  forcetime2 = MPI_WTIME(ierr) ! timing info
  forcetimetot = forcetimetot+(forcetime2-forcetime1)

  CALL MPI_ALL_REDUCE_DOUBLE_SCALAR ( u   ) 
  CALL MPI_ALL_REDUCE_DOUBLE_SCALAR ( vir ) 

  CALL MPI_ALL_REDUCE_DOUBLE ( fx , natm ) 
  CALL MPI_ALL_REDUCE_DOUBLE ( fy , natm ) 
  CALL MPI_ALL_REDUCE_DOUBLE ( fz , natm ) 

  CALL MPI_ALL_REDUCE_DOUBLE ( tau_lj( 1, : ) , 3  )
  CALL MPI_ALL_REDUCE_DOUBLE ( tau_lj( 2, : ) , 3  )
  CALL MPI_ALL_REDUCE_DOUBLE ( tau_lj( 3, : ) , 3  )

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

SUBROUTINE engforce_bmlj_nopbc ( iastart , iaend )

  USE config,           ONLY : natm , box , rx , ry , rz , fx , fy , fz , itype , list , point , ntype , tau_lj , omega
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
  double precision :: u, vir

  u = 0.0D0
  vir = 0.0D0
  fx = 0.0D0
  fy = 0.0D0
  fz = 0.0D0

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
          tau_lj(1,1) = tau_lj(1,1) + rxij * fxij
          tau_lj(1,2) = tau_lj(1,2) + rxij * fyij
          tau_lj(1,3) = tau_lj(1,3) + rxij * fzij
          tau_lj(2,1) = tau_lj(2,1) + ryij * fxij
          tau_lj(2,2) = tau_lj(2,2) + ryij * fyij
          tau_lj(2,3) = tau_lj(2,3) + ryij * fzij
          tau_lj(3,1) = tau_lj(3,1) + rzij * fxij
          tau_lj(3,2) = tau_lj(3,2) + rzij * fyij
          tau_lj(3,3) = tau_lj(3,3) + rzij * fzij
        endif
      endif
    enddo
  enddo
  tau_lj = tau_lj / omega 
  vir = vir / 3.0D0

  CALL MPI_ALL_REDUCE_DOUBLE_SCALAR ( u ) 
  CALL MPI_ALL_REDUCE_DOUBLE_SCALAR ( vir ) 
  
  CALL MPI_ALL_REDUCE_DOUBLE ( fx , natm )
  CALL MPI_ALL_REDUCE_DOUBLE ( fy , natm )
  CALL MPI_ALL_REDUCE_DOUBLE ( fz , natm )
  CALL MPI_ALL_REDUCE_DOUBLE ( tau_lj ( 1 , : )  , 3  )
  CALL MPI_ALL_REDUCE_DOUBLE ( tau_lj ( 2 , : )  , 3  )
  CALL MPI_ALL_REDUCE_DOUBLE ( tau_lj ( 3 , : )  , 3  )

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
  double precision :: forcetime1, forcetime2, invbox
  double precision :: u, vir
  
  forcetime1 = MPI_WTIME(ierr) ! timing info
  
  u   = 0.0D0
  vir = 0.0D0
  fx  = 0.0D0
  fy  = 0.0D0
  fz  = 0.0D0
  
  invbox = 1.0d0/box

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
      if ( ( lvnlist .and. ja .ne. ia ) .or. ( .not. lvnlist .and. ja .gt. ia )  ) then
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

!*********************** SUBROUTINE initialize_coulomb ************************
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
  integer :: nkcut , ncmax 

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
  endif

  return

END SUBROUTINE initialize_coulomb


SUBROUTINE finalize_coulomb

  USE control,  ONLY :  longrange , lcoulomb , calc

  implicit none

  if ( .not. lcoulomb ) return

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
  endif


END SUBROUTINE finalize_coulomb

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

  USE config,  ONLY : natm , atype , natmi , ntype , qia , rx , ry , rz , fx , fy , fz , tau_coul , omega
  USE control, ONLY : cutlongrange,myrank
  USE io_file, ONLY : stdout , ionode 
  USE thermodynamic , ONLY : u_pol, u_coul_tot , vir_coul_tot 
  USE time

  implicit none

  INCLUDE 'mpif.h'

  ! global 
  integer, intent(in) :: iastart , iaend
  double precision , intent(in)    :: mu     ( natm , 3 )
  double precision , intent(out)   :: ef     ( natm , 3 )
  double precision , intent(out)   :: efg    ( natm , 3 , 3 )
  double precision , intent(out)   :: u_coul 
  double precision , intent(out)   :: vir_coul
  double precision , intent(out)   :: phi_coul ( natm ) 

  ! local 
  integer :: ia, ja , ierr , ncell 
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
  !debug 
  integer          :: it , ip
  double precision :: mip    ( ntype , 3 ) , tmp ( 3 , 3 ) ! tmp

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
  write( stdout , '(a)')          'debug: in multipole_DS '
  write( stdout , '(a,i8,a)')     'debug : rm_coul        ',rm_coul%ncmax,rm_coul%meshlabel
  do ia = 1 , natm
    write( stdout , '(a,f12.5)')  'debug : charge (atom)  ',qia(ia)
  enddo
  do it = 1 , ntype
    write( stdout , '(a,f12.5)')  'debug : charge (type ) ',qch(it)
  enddo
  do ia = 1 , natm
    write( stdout , '(a,3f12.5)') 'debug : dipole (atom)  ', mu ( ia , 1 ) , mu ( ia , 2 ) , mu ( ia , 3 )
  enddo
  do it = 1 , ntype
    write( stdout , '(a,3f12.5)') 'debug : dipole (type ) ',  mip ( it , 1 ) , mip ( it , 2 ) , mip ( it , 3 )
  enddo
  write( stdout , '(a,2i8)')      'debug : iastart iaend  ',iastart ,iaend
  write( stdout , '(a,f20.5)')    'debug : cutsq          ',cutsq
  endif
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
              Tx = - rxij * dm3 
              Ty = - ryij * dm3 
              Tz = - rzij * dm3 

              ! multipole interaction tensor rank = 2
              Txx = ( 3.0d0 * rxij * rxij - d2 ) * dm5
              Tyy = ( 3.0d0 * ryij * ryij - d2 ) * dm5
              Tzz = ( 3.0d0 * rzij * rzij - d2 ) * dm5
              Txy = ( 3.0d0 * rxij * ryij      ) * dm5
              Txz = ( 3.0d0 * rxij * rzij      ) * dm5
              Tyz = ( 3.0d0 * ryij * rzij      ) * dm5

              ! multipole interaction tensor rank = 3  
              Txxx = ( - 5.0d0 * rxij * rxij * rxij +  3.0d0 * d2 * ( rxij ) ) * dm7 * 3.0d0
              Tyyy = ( - 5.0d0 * ryij * ryij * ryij +  3.0d0 * d2 * ( ryij ) ) * dm7 * 3.0d0
              Tzzz = ( - 5.0d0 * rzij * rzij * rzij +  3.0d0 * d2 * ( rzij ) ) * dm7 * 3.0d0
              Txxy = ( - 5.0d0 * rxij * rxij * ryij +          d2 * ( ryij ) ) * dm7 * 3.0d0
              Txxz = ( - 5.0d0 * rxij * rxij * rzij +          d2 * ( rzij ) ) * dm7 * 3.0d0
              Tyyx = ( - 5.0d0 * ryij * ryij * rxij +          d2 * ( rxij ) ) * dm7 * 3.0d0
              Tyyz = ( - 5.0d0 * ryij * ryij * rzij +          d2 * ( rzij ) ) * dm7 * 3.0d0
              Tzzx = ( - 5.0d0 * rzij * rzij * rxij +          d2 * ( rxij ) ) * dm7 * 3.0d0
              Tzzy = ( - 5.0d0 * rzij * rzij * ryij +          d2 * ( ryij ) ) * dm7 * 3.0d0
              Txyz = ( - 5.0d0 * rxij * ryij * rzij                          ) * dm7 * 3.0d0

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
  tau_coul = - tau_coul / omega * 0.5d0

#ifdef debug2 
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
  CALL print_tensor( tau_coul ( : , : )     , 'TAU_COUL' )
  CALL print_tensor( efg ( 1 , : , : ) , 'TOTAL   ' )
#endif

  deallocate( fx_coul , fy_coul , fz_coul )
  deallocate( phi_coul_qq , phi_coul_dd  )

  u_coul_tot = u_coul
  vir_coul_tot = vir_coul

  ttt3 = MPI_WTIME(ierr)
  fcoultimetot2 = fcoultimetot2 + ( ttt3 - ttt2 )

  return

END SUBROUTINE multipole_DS



SUBROUTINE multipole_ES ( iastart, iaend , ef , efg , mu , u_coul , vir_coul , phi_coul )

  USE control,   ONLY : lsurf
  USE config,    ONLY : natm , ntype , natmi , atype , rx , ry , rz , fx , fy , fz , tau_coul , qia , omega , box
  USE constants, ONLY : imag , pi , piroot , tpi , fpi
  USE io_file  , ONLY : ionode , stdout 
  USE time
  USE thermodynamic , ONLY : u_pol, u_coul_tot , vir_coul_tot 

  implicit none

  INCLUDE 'mpif.h'

  ! global 
  integer, intent(in) :: iastart , iaend
  double precision    :: ef     ( natm , 3 )
  double precision    :: efg    ( natm , 3 , 3 )
  double precision    :: mu     ( natm , 3 )
  double precision    :: tau_dir( 3 , 3 )
  double precision    :: tau_rec( 3 , 3 )
  double precision    :: u_coul 
  double precision    :: vir_coul 
  double precision    :: phi_coul ( natm )  

  ! local 
  integer             :: ia , ja , ik , it , ip , ierr
  integer             :: nxij , nyij , nzij
  double precision, dimension(:)    , allocatable :: fx_coul , fy_coul , fz_coul
  double precision, dimension(:)    , allocatable :: fx_dir , fy_dir , fz_dir
  double precision, dimension(:)    , allocatable :: fx_rec , fy_rec , fz_rec
  double precision, dimension(:)    , allocatable :: fx_surf , fy_surf , fz_surf
  double precision, dimension(:)    , allocatable :: phi_dir , phi_rec , phi_surf , phi_self
  double precision, dimension(:,:)  , allocatable :: ef_dir , ef_rec , ef_surf , ef_self
  double precision, dimension(:,:,:), allocatable :: efg_dir , efg_rec , efg_self

  double precision :: u_dir , u_rec , u_surf , u_self
  double precision :: u_surf_qq , u_surf_qd , u_surf_dd
  double precision :: vir_dir , vir_rec , vir_surf , vir_self

  double precision :: mip    ( ntype , 3 )

  double precision :: rxi  , ryi  , rzi
  double precision :: rxj  , ryj  , rzj
  double precision :: kx   , ky   , kz 
  double precision :: rxij , ryij , rzij
  double precision :: fxij , fyij , fzij
  double precision :: kk   , kri  , Ak , kcoe
  double precision :: qi   , qj   , qij
  double precision :: muix , muiy , muiz
  double precision :: mujx , mujy , mujz
  double complex   :: rhon , rhonmu , carg 
  double precision :: str  , recarg , recargi1 
  double precision :: expon , F0 , F1 , F2 , F3 
  double precision :: k_dot_mu
  double precision :: T
  double precision :: Tx , Ty , Tz
  double precision :: Txx , Tyy , Tzz , Txy , Txz , Tyz
  double precision :: Txxx,  Tyyy,  Tzzz, Txxy, Txxz, Tyyx, Tyyz, Tzzx, Tzzy, Txyz
  double precision :: d , d2 , d3  , d5 
  double precision :: dm1 , dm3 , dm5 , dm7 
  double precision :: alpha2 , alpha3 , alpha5 
  double precision :: selfa , selfa2 
  double precision, external :: errfc
  double precision :: qtot ( 3 ) , qsq , mutot ( 3 ) , musq , qmu_sum ( 3 ) 
  double precision :: invbox
  double precision :: tpi_V , tpi_3V , fpi_V , fpi_3V 
  double precision :: ttt1 , ttt2  , ttt3 , ttt4

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
  qtot  = 0.0d0
  qsq   = 0.0d0
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


#ifdef debug
  if ( ionode ) then
  write( stdout , '(a)')          'debug : in engforce_charge_ES'
  write( stdout , '(a,i8,a,a)')   'debug : km_coul       ',km_coul%nkcut,' ',km_coul%meshlabel
  do ia = 1 , natm
    write( stdout , '(a,f12.5)')  'debug : charge (atom) ',qia(ia)
  enddo
  do it = 1 , ntype
    write( stdout , '(a,f12.5)')  'debug : charge (type) ',qch(it)
  enddo
  do ia = 1 , natm
    write( stdout , '(a,3f12.5)') 'debug : dipole (atom) ', mu ( ia , 1 ) , mu ( ia , 2 ) , mu ( ia , 3 )
  enddo
  do it = 1 , ntype
    write( stdout , '(a,3f12.5)') 'debug : dipole (type) ',  mip ( it , 1 ) , mip ( it , 2 ) , mip ( it , 3 )
  enddo
  write( stdout , '(a,2i8)')      'debug : iastart iaend ',iastart ,iaend
  write( stdout , '(a,f20.5)')    'debug : alphaES       ', alphaES
  endif
  call print_config_sample(0,0)
#endif 

  ttt1 = MPI_WTIME(ierr)

  allocate( ef_dir ( natm , 3 ) , ef_rec ( natm , 3 ) , ef_surf ( natm , 3 ) , ef_self ( natm , 3 ) ) 
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
  fx   = 0.0D0
  fy   = 0.0D0
  fz   = 0.0D0
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
  u_dir    = 0.0d0 
  u_rec    = 0.0d0 
  u_self   = 0.0d0 
  u_surf   = 0.0d0 
  vir_dir  = 0.0d0 
  vir_rec  = 0.0d0 
  vir_surf = 0.0d0 
  vir_self = 0.0d0 
  tau_dir  = 0.0d0
  tau_rec  = 0.0d0

  ! =====================
  ! facteur de structure 
  ! =====================
  CALL struc_fact ( km_coul )

  ! =================
  !  some constants 
  ! =================
  tpi_V  = tpi    / omega  ! 2pi / V
  tpi_3V = tpi_V  / 3.0d0  ! 2pi / 3V 
  fpi_V  = tpi_V  * 2.0d0  ! 4pi / V
  fpi_3V = tpi_3V * 2.0d0  ! 4pi / 3V
  invbox = 1.0d0 / box
  alpha2 = alphaES * alphaES
  alpha3 = alpha2  * alphaES
  alpha5 = alpha3  * alpha2
  selfa  = alphaES / piroot
  selfa2 = 2.0d0 * selfa * alpha2 / 3.0d0

  ! ==============================================
  !        direct space part
  ! ==============================================

  do ia = iastart , iaend
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

        nxij = NINT ( rxij * invbox )
        nyij = NINT ( ryij * invbox )
        nzij = NINT ( rzij * invbox )
        rxij = rxij - box * nxij
        ryij = ryij - box * nyij
        rzij = rzij - box * nzij

        d2  = rxij * rxij + ryij * ryij + rzij * rzij
        d   = SQRT ( d2 )
        d3  = d2 * d 
        d5  = d3 * d2 
        dm1 = 1.0d0 / d
        dm3 = dm1 / d2 
        dm5 = dm3 / d2 
        dm7 = dm5 / d2 
 
        expon = EXP ( - alpha2 * d2 )    / piroot
        F0    = errfc( alphaES * d )
        F1    = F0 + 2.0d0 * alphaES * d  * expon 
        F2    = F1 + 4.0d0 * alpha3  * d3 * expon / 3.0d0
        F3    = F2 + 8.0d0 * alpha5  * d5 * expon / 15.0d0

        ! multipole interaction tensor rank = 0 
        T  = dm1 * F0

        ! multipole interaction tensor rank = 1
        Tx = - rxij * F1 * dm3
        Ty = - ryij * F1 * dm3
        Tz = - rzij * F1 * dm3

        ! multipole interaction tensor rank = 2
        Txx = ( 3.0d0 * rxij * rxij * F2 - d2 * F1 ) * dm5 
        Tyy = ( 3.0d0 * ryij * ryij * F2 - d2 * F1 ) * dm5 
        Tzz = ( 3.0d0 * rzij * rzij * F2 - d2 * F1 ) * dm5
        Txy = ( 3.0d0 * rxij * ryij * F2           ) * dm5
        Txz = ( 3.0d0 * rxij * rzij * F2           ) * dm5 
        Tyz = ( 3.0d0 * ryij * rzij * F2           ) * dm5 
 
        ! multipole interaction tensor rank = 3  
        Txxx = ( - 5.0d0 * rxij * rxij * rxij * F3 +  3.0d0 * d2 * ( rxij ) * F2 ) * dm7 * 3.0d0
        Tyyy = ( - 5.0d0 * ryij * ryij * ryij * F3 +  3.0d0 * d2 * ( ryij ) * F2 ) * dm7 * 3.0d0
        Tzzz = ( - 5.0d0 * rzij * rzij * rzij * F3 +  3.0d0 * d2 * ( rzij ) * F2 ) * dm7 * 3.0d0
        Txxy = ( - 5.0d0 * rxij * rxij * ryij * F3 +          d2 * ( ryij ) * F2 ) * dm7 * 3.0d0
        Txxz = ( - 5.0d0 * rxij * rxij * rzij * F3 +          d2 * ( rzij ) * F2 ) * dm7 * 3.0d0
        Tyyx = ( - 5.0d0 * ryij * ryij * rxij * F3 +          d2 * ( rxij ) * F2 ) * dm7 * 3.0d0
        Tyyz = ( - 5.0d0 * ryij * ryij * rzij * F3 +          d2 * ( rzij ) * F2 ) * dm7 * 3.0d0
        Tzzx = ( - 5.0d0 * rzij * rzij * rxij * F3 +          d2 * ( rxij ) * F2 ) * dm7 * 3.0d0
        Tzzy = ( - 5.0d0 * rzij * rzij * ryij * F3 +          d2 * ( ryij ) * F2 ) * dm7 * 3.0d0
        Txyz = ( - 5.0d0 * rxij * ryij * rzij * F3                               ) * dm7 * 3.0d0

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
    kx   = km_coul%kpt(1,ik)
    ky   = km_coul%kpt(2,ik)
    kz   = km_coul%kpt(3,ik)
    kk   = km_coul%kptk(ik)
    Ak   = EXP ( - kk * 0.25d0 / alpha2 ) / kk
    kcoe = 2.0d0 * ( 1.0d0 / alpha2  + 1.0d0 / kk )

    if (km_coul%kptk(ik) .eq. 0 ) then
      WRITE ( stdout , * ) 'the sum should be done on k! =  0',ik
      STOP
    endif
    ! ===============================
    !                              ---
    !  charge density in k-space ( \   q * facteur de structure  )
    !                              /__
    ! ===============================
    rhon   = (0.d0, 0.d0)
    rhonmu = (0.d0, 0.d0)
    do it = 1, ntype
      k_dot_mu = ( mip ( it , 1 ) * kx + mip ( it , 2 ) * ky + mip ( it , 3 ) * kz  ) 
      rhon = rhon + ( qch(it) + imag * k_dot_mu ) * km_coul%strf ( ik , it ) 
      rhonmu = rhonmu + imag * k_dot_mu * km_coul%strf ( ik , it ) 
    enddo

    str = rhon * CONJG(rhon)

    u_rec   = u_rec   + Ak * str 
    vir_rec = vir_rec + Ak * str * ( 3.0d0 - kcoe * kk ) / 3.0d0

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

      ! elecctric field gradient
      efg_rec ( ia , 1 , 1 ) = efg_rec ( ia , 1 , 1 ) +  Ak * kx * kx * recarg
      efg_rec ( ia , 2 , 2 ) = efg_rec ( ia , 2 , 2 ) +  Ak * ky * ky * recarg 
      efg_rec ( ia , 3 , 3 ) = efg_rec ( ia , 3 , 3 ) +  Ak * kz * kz * recarg 
      efg_rec ( ia , 1 , 2 ) = efg_rec ( ia , 1 , 2 ) +  Ak * kx * ky * recarg 
      efg_rec ( ia , 1 , 3 ) = efg_rec ( ia , 1 , 3 ) +  Ak * kx * kz * recarg 
      efg_rec ( ia , 2 , 3 ) = efg_rec ( ia , 2 , 3 ) +  Ak * ky * kz * recarg 
 
      ! stress tensor
      tau_rec(1,1) = tau_rec(1,1) + ( 1.0d0 - kcoe * kx * kx ) * str * ak
      tau_rec(1,2) = tau_rec(1,2) -           kcoe * kx * ky   * str * ak
      tau_rec(1,3) = tau_rec(1,3) -           kcoe * kx * kz   * str * ak
      tau_rec(2,1) = tau_rec(2,1) -           kcoe * ky * kx   * str * ak 
      tau_rec(2,2) = tau_rec(2,2) + ( 1.0d0 - kcoe * ky * ky ) * str * ak
      tau_rec(2,3) = tau_rec(2,3) -           kcoe * ky * kz   * str * ak
      tau_rec(3,1) = tau_rec(3,1) -           kcoe * kz * kx   * str * ak
      tau_rec(3,2) = tau_rec(3,2) -           kcoe * kz * ky   * str * ak
      tau_rec(3,3) = tau_rec(3,3) + ( 1.0d0 - kcoe * kz * kz ) * str * ak

    enddo

  enddo kpoint

  ttt3 = MPI_WTIME(ierr)
  fcoultimetot2 = fcoultimetot2 + ( ttt3 - ttt2 )

  CALL MPI_ALL_REDUCE_DOUBLE_SCALAR ( u_dir )
  CALL MPI_ALL_REDUCE_DOUBLE_SCALAR ( vir_dir )

  CALL MPI_ALL_REDUCE_DOUBLE ( fx_dir , natm )
  CALL MPI_ALL_REDUCE_DOUBLE ( fy_dir , natm )
  CALL MPI_ALL_REDUCE_DOUBLE ( fz_dir , natm )

  CALL MPI_ALL_REDUCE_DOUBLE ( phi_dir , natm )

  CALL MPI_ALL_REDUCE_DOUBLE ( ef_dir ( : , 1 )  , natm )
  CALL MPI_ALL_REDUCE_DOUBLE ( ef_dir ( : , 2 )  , natm )
  CALL MPI_ALL_REDUCE_DOUBLE ( ef_dir ( : , 3 )  , natm )

  CALL MPI_ALL_REDUCE_DOUBLE ( efg_dir ( : , 1 , 1 ) , natm )
  CALL MPI_ALL_REDUCE_DOUBLE ( efg_dir ( : , 2 , 2 ) , natm )
  CALL MPI_ALL_REDUCE_DOUBLE ( efg_dir ( : , 3 , 3 ) , natm )
  CALL MPI_ALL_REDUCE_DOUBLE ( efg_dir ( : , 1 , 2 ) , natm )
  CALL MPI_ALL_REDUCE_DOUBLE ( efg_dir ( : , 1 , 3 ) , natm )
  CALL MPI_ALL_REDUCE_DOUBLE ( efg_dir ( : , 2 , 3 ) , natm )

  CALL MPI_ALL_REDUCE_DOUBLE ( tau_dir ( 1 , : ) , 3 )
  CALL MPI_ALL_REDUCE_DOUBLE ( tau_dir ( 2 , : ) , 3 )
  CALL MPI_ALL_REDUCE_DOUBLE ( tau_dir ( 3 , : ) , 3 )

  ! ======================================================
  ! remark on the unit :
  ! 1/(4*pi*epislon_0) = 1 => epsilon_0 = 1/4pi
  ! ======================================================

  tau_dir =   tau_dir * 0.5d0 !/ omega
  u_dir   =   u_dir   * 0.5d0
  vir_dir = - vir_dir * 0.5d0

  tau_rec =   tau_rec * tpi_V / omega
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
  u_surf_qd = 2.0d0 * ( qtot ( 1 ) * mutot ( 1 ) + qtot ( 2 ) * mutot ( 2 ) + qtot ( 3 ) * mutot ( 3 ) )
  u_surf = u_surf_qq + u_surf_qd + u_surf_dd

  u_surf = u_surf * tpi_3V 

  vir_surf = - ( u_surf_qq + 4.0d0 * u_surf_qd + 3.0d0 * u_surf_dd ) 
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
!  vir_self = - 9.0d0 * pi / 2.0d0 * qsq / alpha2 / omega 
 
  ! potential , field , field gradient ( no contrib to forces ) 
  do ia = 1 , natm
   phi_self ( ia ) = - 2.0d0 * selfa * qia ( ia )
   ef_self( ia , 1 ) = 2.0d0 * selfa2 * mu ( ia , 1 ) 
   ef_self( ia , 2 ) = 2.0d0 * selfa2 * mu ( ia , 2 ) 
   ef_self( ia , 3 ) = 2.0d0 * selfa2 * mu ( ia , 3 ) 
   efg_self ( ia , 1 , 1 ) = - 2.0d0 * selfa2 * qia ( ia ) 
   efg_self ( ia , 2 , 2 ) = - 2.0d0 * selfa2 * qia ( ia ) 
   efg_self ( ia , 3 , 3 ) = - 2.0d0 * selfa2 * qia ( ia ) 
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
    if ( ionode ) WRITE ( stdout , '(a)' ) 'Surface contribution is not added to the different quantities ( see multipole_ES inf field.f90 ) '
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

#ifdef debug 
  if ( ionode ) then
    WRITE ( stdout , '(a)' ) ' '
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

    WRITE ( stdout , '(a)' ) ' '

    WRITE ( stdout , '(a)' ) ' '
    WRITE ( stdout , '(a)' ) 'forces at atoms : '
    do ia = 1 , natm
      WRITE ( stdout , '(i5,a3,a,3f18.10)' ) &
      ia,atype(ia),' f       = ', fx ( ia )  , fy ( ia ) , fz ( ia )
    enddo
    WRITE ( stdout , '(a)' ) ' '
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

    WRITE ( stdout , '(a)' ) ' '
    WRITE ( stdout , '(a)' ) ' '
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
    WRITE ( stdout , '(a)' ) 'Energy and virial : '
    WRITE ( stdout , '(5(a,f16.8))' ) ,'vir_dir    = ', vir_dir    ,' vir_rec    = ', vir_rec    ,' vir_surf    = ',vir_surf   ,' vir_self    = '   ,vir_self   ,'                           vir_coul    = ',vir_coul
    WRITE ( stdout , '(5(a,f16.8))' ) ,'phi_dir(1) = ', phi_dir(1) ,' phi_rec(1) = ', phi_rec(1) ,' phi_surf(1) = ',phi_surf(1),' phi_self(1) = '   ,phi_self(1),'                           phi_coul(1) = ',phi_coul (1)
    WRITE ( stdout , '(6(a,f16.8))' ) ,'u_dir      = ', u_dir      ,' u_rec      = ', u_rec      ,' u_surf      = ',u_surf     ,' u_self      = '   ,u_self     ,' u_pol  = ',u_pol,' u_coul      = ',u_coul
  CALL print_tensor( tau_dir  ( : , : )     , 'TAU_DIR ' )
  CALL print_tensor( tau_rec  ( : , : )     , 'TAU_REC ' )
  CALL print_tensor( tau_coul ( : , : )     , 'TAU_COUL' )

  deallocate( ef_dir  , ef_rec , ef_surf , ef_self ) 
  deallocate( efg_dir , efg_rec , efg_self ) 
  deallocate( fx_coul , fy_coul , fz_coul )
  deallocate( fx_dir  , fy_dir  , fz_dir )
  deallocate( fx_rec  , fy_rec  , fz_rec  )
  deallocate( fx_surf , fy_surf , fz_surf )
  deallocate( phi_dir , phi_rec , phi_surf , phi_self )

  u_coul_tot = u_coul
  vir_coul_tot = vir_coul

  ttt4 = MPI_WTIME(ierr)
  fcoultimetot3 = fcoultimetot3 + ( ttt4 - ttt3 )

  return

END SUBROUTINE multipole_ES

!*********************** SUBROUTINE engforce_morse_pbc *************************
!
! total potential energy forces for each atoms for a morse potential with 
! periodic boundaries conditions, with or without vnlist
! (lvnlist=.TRUE.OR.FALSE.)
!
!******************************************************************************

SUBROUTINE engforce_morse_pbc ( iastart , iaend )

  USE config,           ONLY : natm , box , rx , ry , rz , fx , fy , fz, atype , itype , list , point , ntype
  USE control,          ONLY : lvnlist , myrank 
  USE thermodynamic,    ONLY : u_morse , vir_morse
  USE time
  USE io_file,          ONLY : stdout , ionode

  implicit none
  INCLUDE 'mpif.h'

  ! global
  integer, intent(in)  :: iastart , iaend

  ! local
  integer          :: ia , ja , j1 , je , jb , ierr
  integer          :: p1 , p2
  integer          :: nxij , nyij , nzij
  double precision :: rxi , ryi , rzi
  double precision :: rxij , ryij , rzij , rijsq , rij , erh , erh2 
  double precision :: wij , fxij , fyij , fzij
  double precision :: forcetime1 , forcetime2 , invbox
  double precision :: u , vir
  double precision :: expon 

#ifdef debug
  if ( ionode ) then
  do ia=1, natm
    write( stdout , '(a,a)' )  'debug : atype',atype(ia)
    write( stdout , '(a,i4)' ) 'debug : itype',itype(ia)
  enddo
  endif
#endif

  forcetime1 = MPI_WTIME(ierr) ! timing info

  u   = 0.0D0
  vir = 0.0D0
  fx  = 0.0D0
  fy  = 0.0D0
  fz  = 0.0D0

  invbox = 1.0d0 / box

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
           rij   = SQRT ( rijsq ) 
           erh   = EXP ( - rij * rhomor (p1,p2) ) 
           erh2  = erh * erh 
           expon = erh * rs (p1,p2) - 1.0d0
           expon = expon * expon - 1.0d0 
           u     = u + epsmor(p1,p2) * expon 
!          if ( trunc .eq. 1 ) theni
!            u =  u + epsp(p1,p2) * ( plj(p1,p2) * srq -qlj(p1,p2) * srp ) - uc(p1,p2)
!          endif
!          if ( trunc .eq. 2 ) then
!            u =  u + epsp(p1,p2) * ( plj(p1,p2) * srq -qlj(p1,p2) * srp ) + uc1(p1,p2) * rijsq - uc2(p1,p2)
!          endif
          wij = fm(p1,p2) * ( erh - rs(p1,p2) * erh2 ) 
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
!  vir = vir/3.0D0

  forcetime2 = MPI_WTIME(ierr) ! timing info
  forcetimetot = forcetimetot+(forcetime2-forcetime1)

  CALL MPI_ALL_REDUCE_DOUBLE_SCALAR ( u   )
  CALL MPI_ALL_REDUCE_DOUBLE_SCALAR ( vir )

  CALL MPI_ALL_REDUCE_DOUBLE ( fx , natm )
  CALL MPI_ALL_REDUCE_DOUBLE ( fy , natm )
  CALL MPI_ALL_REDUCE_DOUBLE ( fz , natm )

  u_morse = u
  vir_morse = vir

  write(*,'(f24.12)') u_morse

  return

END SUBROUTINE engforce_morse_pbc

!*********************** SUBROUTINE engforce_morse_pbc *************************
!
! This routines evaluates the dipole moment induced by polarizabilities on atoms. 
! The evaluation is done self-consistently starting from the the field due the
! point charges only.
! The stopping criteria is governed by conv_tol_ind
!
!******************************************************************************

SUBROUTINE moment_from_pola ( iastart , iaend , mu_ind ) 

  USE io_file, ONLY : ionode , stdout 
  USE config,  ONLY : natm , atype , fx , fy , fz , ntype , dipia , qia , ntypemax
  USE control, ONLY : longrange
  USE thermodynamic, ONLY : u_pol

  implicit none

  ! global
  integer :: iastart , iaend 
  double precision :: mu_ind ( natm , 3 ) 

  ! local
  integer :: ia , iscf , it 
  logical :: linduced
  double precision :: Efield_old , diff_efield
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
  ! =========================
  !  charges are set to zero 
  ! =========================
  qch_tmp = qch
  qia_tmp = qia
  qch = 0.0d0
  qia = 0.0d0
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
    if ( longrange .eq. 'direct' ) CALL multipole_DS ( iastart, iaend , Efield_ind , efg_tmp , mu_ind , u_tmp , vir_tmp , phi_tmp ) 
    if ( longrange .eq. 'ewald' )  CALL multipole_ES ( iastart, iaend , Efield_ind , efg_tmp , mu_ind , u_tmp , vir_tmp , phi_tmp ) 
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

  ! ===========================
  !  charge info is recovered
  ! ===========================
  qch = qch_tmp
  qia = qia_tmp


  if ( ionode ) then
    WRITE ( stdout , '(a)' ) ' '
    WRITE ( stdout , '(a,i6,a)') 'scf calculation of the induced electric moment converged in ',iscf, ' iterations '
    WRITE ( stdout , '(a,e10.3)') 'Electric field is converged within ',conv_tol_ind
    WRITE ( stdout , '(a)' ) ' '
#ifdef debug
    WRITE ( stdout , '(a)' )     'Induced dipoles at atoms : '
    do ia = 1 , natm
      WRITE ( stdout , '(i5,a3,a,3f18.10)' ) &
      ia,atype(ia),' mu_ind = ', mu_ind ( ia , 1 ) , mu_ind ( ia , 2 ) , mu_ind ( ia , 3 )
    enddo
    WRITE ( stdout , '(a)' ) ' '
#endif 
  endif


  return

END SUBROUTINE moment_from_pola



END MODULE field 
! ===== fmV =====
