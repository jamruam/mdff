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

  double precision :: qlj     ( ntypemax , ntypemax )
  double precision :: plj     ( ntypemax , ntypemax )
  double precision :: epslj   ( ntypemax , ntypemax )
  double precision :: sigmalj ( ntypemax , ntypemax )

  double precision :: mass  ( ntypemax )    ! masses ( not yet )
  double precision :: qch   ( ntypemax )    ! charges only for efg (no electrostatic interaction and forces)
  double precision :: dip   ( ntypemax , 3 )! dipole for efg (no electrostatic interaction and forces)
  double precision :: dip_ind( ntypemax , 3 )! dipole for efg (no electrostatic interaction and forces)
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

  if ( calc .eq. 'efg' ) return

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
    if ( calc .ne. 'efg' ) then
    
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
       ! sigma*^2 
       sigsq( it , jt ) = two16 * two16 * sig ( it , jt ) * sig ( it , jt )
       ! sigma*^2 / rc ^2
       sr2 = sigsq ( it , jt )/rcutsq ( it , jt )
       ! sigma* / rc 
       sr = dsqrt(sr2)
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


!*********************** SUBROUTINE engforce **********************************
!
! this subroutine is used as an interfaced to the different potential + forces
! subroutines
!
!******************************************************************************

SUBROUTINE engforce ( iastart , iaend )!, list , point )

  USE config,   ONLY :  natm
  USE control,  ONLY :  lpbc , lminimg , lbmlj , lcoulomb , lshiftpot , longrange
  USE io_file,  ONLY :  stdout 

  implicit none

  ! global
  integer, intent(inout)  :: iastart , iaend 
!  integer, intent(out) :: list( 250 * natm ) , point( natm+1 )
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
    WRITE ( stdout ,'(a)') ' lj + coulomb' 
#endif  
    ! =======
    !   PBC
    ! =======
    if ( lpbc ) then
      if ( lshiftpot .and. lminimg ) CALL engforce_bmlj_pbc         ( iastart , iaend )
      if ( .not.lshiftpot )          CALL engforce_bmlj_pbc_noshift ( iastart , iaend )
      if ( longrange .eq. 'ewald'  ) CALL engforce_coulomb_ES       (  )
      if ( longrange .eq. 'direct' ) CALL engforce_coulomb_DS       ( iastart , iaend ) 
    ! =======
    !  NO PBC
    ! =======
    else
      WRITE ( stdout ,'(a)') 'not yet : coulomb + nopbc'
      STOP 
      if ( lshiftpot )               CALL engforce_bmlj_nopbc       ( iastart , iaend )
!      if ( .not.lshiftpot )          CALL engforce_bmlj_pbc_noshift ( iastart , iaend ) 
!      if ( longrange .eq. 'ewald' )  CALL engforce_coulomb_ES_nopbc (  )
!      if ( longrange .eq. 'direct')  CALL engforce_coulomb_DS_nopbc (  )
    endif
  ! ================================
  !      ONLY LENNARD-JONES
  ! ================================
  elseif ( lbmlj ) then
#ifdef debug  
    WRITE ( stdout , '(a)' ) ' lj only' 
#endif  
    ! =======
    !   PBC
    ! =======
    if ( lpbc ) then
#ifdef debug  
    WRITE ( stdout , '(a)' ) ' pbc ' 
#endif  
      if ( lshiftpot       )     CALL engforce_bmlj_pbc         ( iastart , iaend )
      if ( .not. lshiftpot )     CALL engforce_bmlj_pbc_noshift ( iastart , iaend )
    ! =======
    !  NO PBC
    ! =======
    else
#ifdef debug  
    WRITE ( stdout , '(a)' ) ' no pbc ' 
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
    WRITE ( stdout ,'(a)') ' coulomb only ' 
#endif  
    ! =======
    !   PBC
    ! =======
    if ( lpbc ) then
      if ( longrange .eq. 'ewald'  )  CALL engforce_coulomb_ES ( )
      if ( longrange .eq. 'direct' )  CALL engforce_coulomb_DS ( iastart , iaend ) 
    ! =======
    !  NO PBC
    ! =======
    else
      WRITE ( stdout ,'(a)') 'not yet : coulomb + nopbc'
      STOP 
!      if ( longrange .eq. 'ewald'  )  CALL engforce_coulomb_ES_nopbc ( )
!      if ( longrange .eq. 'direct' )  CALL engforce_coulomb_DS_nopbc ( )
    endif

  endif

  return

END SUBROUTINE engforce



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
  print*,'atype',atype
  print*,'itype',itype
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
          nxij = nint( rxij * invbox )
          nyij = nint( ryij * invbox )
          nzij = nint( rzij * invbox )
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
        rxij = rxi - rx ( ja )
        ryij = ryi - ry ( ja )
        rzij = rzi - rz ( ja )
        nxij = nint( rxij * invbox )
        nyij = nint( ryij * invbox )
        nzij = nint( rzij * invbox )
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
    nxij = nint( rxij * invbox )
    nyij = nint( ryij * invbox )
    nzij = nint( rzij * invbox )
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

  ! ==========================
  !  set some type parameters
  ! ==========================      
  natmi ( 0 ) = 0
  cc = 0
  do it = 1 , ntype
      ccs = cc
      cc = cc + natmi ( it )
    do ia = ccs + 1 , cc
      qia ( ia ) = qch ( it )
    enddo
  enddo
  natmi ( 0 )  = natm

  ! ============
  !  direct sum
  ! ============
  if ( longrange .eq. 'direct' ) then
    ncmax = ( 2 * ncelldirect + 1 ) ** 3
    rm_coul%ncmax=ncmax
    rm_coul%ncell=ncelldirect
    allocate( rm_coul%boxxyz( 3 , ncmax ) , rm_coul%lcell( ncmax ) )
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
      allocate ( km_coul%strf ( nkcut, ntype) ) 
      CALL kpoint_sum_init ( km_coul )
      km_coul_dip%ncell = ncellewald
      km_coul_dip%meshlabel='km_coul_dip'
!      nkcut = ( ncellewald + 1 ) ** 3
      nkcut = ( ncellewald + 1 ) ** 3
      nkcut = nkcut - 1
      km_coul_dip%nkcut = nkcut
      allocate( km_coul_dip%kptk( nkcut ) , km_coul_dip%kpt(3,nkcut) )
      allocate ( km_coul_dip%strf ( nkcut, ntype) )
      CALL kpoint_sum_init_half ( km_coul_dip )

  endif

  return

END SUBROUTINE initialize_coulomb


SUBROUTINE finalize_coulomb

  USE control,  ONLY :  longrange , lcoulomb , calc
  USE prop,     ONLY :  lefg

  implicit none

  if ( .not. lcoulomb .or. calc .ne. 'efg' .or. .not. lefg ) return

  ! ============
  !  direct sum
  ! ============
  if ( longrange .eq. 'direct' ) then
    deallocate( rm_coul%boxxyz , rm_coul%lcell )
  endif

  ! ============
  !  ewald sum
  ! ============
  if ( longrange .eq. 'ewald')  then
      deallocate( km_coul%kptk , km_coul%kpt )
      deallocate ( km_coul%strf )
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
SUBROUTINE engforce_coulomb_ES  

  USE control,          ONLY :  myrank , numprocs, calc , longrange
  USE config,           ONLY :  system , natm , natmi , atype , atypei , box , omega , &
                                itype , rx , ry , rz , fx , fy , fz , ntype , qia 
  USE io_file,          ONLY :  ionode , stdout 
  USE constants,        ONLY :  pi , fpi , piroot , imag , tpi
  USE kspace,           ONLY :  struc_fact
  USE thermodynamic,    ONLY :  u_coul , vir_coul
  USE time

  implicit none

  INCLUDE 'mpif.h'

  ! local
  integer :: ia, ja , it , ierr
  integer :: ik
  integer :: nxij , nyij , nzij
  double precision :: rij , rijsq , expon , invbox
  double precision :: alpha2
  double precision :: qi , qj , qij, qijf 
  double precision :: rxi , ryi , rzi , rxij , ryij , rzij 
  double precision :: wij0 , wij , fxij , fyij , fzij 
  double precision :: ak, kx, ky, kz, kk
  double precision :: kri , str
  double precision :: u_dir , u_rec 
  double precision :: vir_dir , vir_rec 
  double complex   :: rhon , carg 
  double precision, external :: errfc 
  double precision :: ttt1 , ttt2 , ttt3 


  if (longrange .eq. 'direct' ) return

  ttt1 = MPI_WTIME(ierr)

  u_rec   = 0.0d0
  u_dir   = 0.0d0
  vir_rec = 0.0d0
  vir_dir = 0.0d0
  
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

! ==============================================
!        direct space part
! ==============================================

   do ia = 1 , natm 
     rxi = rx(ia)
     ryi = ry(ia)
     rzi = rz(ia)
     qi = qia(ia)
     do ja = 1, natm

       if (ja .ne. ia ) then
         qj = qia(ja)
         rxij = rxi - rx(ja)
         ryij = ryi - ry(ja)
         rzij = rzi - rz(ja)

         nxij = nint( rxij * invbox )
         nyij = nint( ryij * invbox )
         nzij = nint( rzij * invbox )

         rxij = rxij - box * nxij
         ryij = ryij - box * nyij
         rzij = rzij - box * nzij

         rijsq = rxij * rxij + ryij * ryij + rzij * rzij
         rij = dsqrt( rijsq )

         qij  = qi * qj / rij
         qijf = qij / rijsq         ! qi * qj / r^3

         wij0 = errfc( alphaES * rij ) 
         expon = dexp( - alpha2 * rijsq ) / piroot
         wij  = qijf * ( wij0 + 2.0d0 * rij * alphaES * expon )

         u_dir = u_dir + qij * wij0

         fxij = wij * rxij
         fyij = wij * ryij
         fzij = wij * rzij

         fx(ia) = fx(ia) + fxij
         fy(ia) = fy(ia) + fyij
         fz(ia) = fz(ia) + fzij

       endif
     enddo
  enddo

  ttt2 = MPI_WTIME(ierr)
  fcoultimetot1 = fcoultimetot1 + ( ttt2 - ttt1 ) 
  
! ==============================================
!            reciprocal space part
! ==============================================
  ! sum on atoms for one givn ik
  do ia = 1 , natm
    rxi = rx(ia)
    ryi = ry(ia)
    rzi = rz(ia)
    qi = fpi * qia(ia) / omega
    kpoint : do ik = 1, km_coul%nkcut 
      ! =================
      !   k-space  
      ! =================
      kx = km_coul%kpt(1,ik)
      ky = km_coul%kpt(2,ik)
      kz = km_coul%kpt(3,ik)
      kk = km_coul%kptk(ik)
      ak = dexp( - kk * 0.25d0 / alpha2 ) / kk

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
        rhon = rhon + qch(it) * CONJG(km_coul%strf ( ik , it ) )
      enddo
      kri = ( kx * rxi + ky * ryi + kz * rzi ) 
      carg = EXP ( imag * kri )
      wij = rhon * carg * ak * imag 

      str = rhon * CONJG(rhon) 
!     vir_rec = vir_rec + wij * rijsq
      u_rec = u_rec + ak * str

      fxij = wij * kx 
      fyij = wij * ky
      fzij = wij * kz

      fx(ia) = fx(ia) - qi * fxij
      fy(ia) = fy(ia) - qi * fyij
      fz(ia) = fz(ia) - qi * fzij

    enddo kpoint

  enddo 

  ! ======================================================
  ! remark on the unit :
  ! 1/(4*pi*epislon_0) = 1 => epsilon_0 = 1/4pi
  ! ======================================================
  u_rec = u_rec * tpi / omega 

  vir_coul = ( vir_dir + vir_rec )
  u_coul   = ( u_dir   + u_rec   )
  print*,'u_coul_dir =',u_dir,'u_coul_rec =',u_rec,'u_coul_tot =',u_coul

  ttt3 = MPI_WTIME(ierr)
  fcoultimetot2 = fcoultimetot2 + ( ttt3 - ttt2 ) 
  
  return


END SUBROUTINE engforce_coulomb_ES


!*********************** SUBROUTINE engforce_coulomb_DS ***********************
!
! this subroutine calculates the potential energy and the force acting on each
! charge due to the coulombic interaction. The Direct Summation is used
! 
! Note :  same consruction as efg_DS
!
!******************************************************************************
SUBROUTINE engforce_coulomb_DS (  iastart , iaend ) 

  USE control,          ONLY :  myrank , numprocs, calc , longrange
  USE config,           ONLY :  system , natm , natmi , atype , atypei , box , omega , &
                                itype , rx , ry , rz , fx , fy , fz , ntype , qia 
  USE io_file,          ONLY :  ionode , stdout 
  USE constants,        ONLY :  pi , fpi , piroot , imag , tpi
  USE thermodynamic,    ONLY :  u_coul , vir_coul
  USE time

  implicit none

  INCLUDE 'mpif.h'

  ! global
  integer, intent(in)  :: iastart , iaend 

  ! local
  integer :: ia, ja , ierr , ncell , it
  double precision :: rij , rijsq 
  double precision :: qi , qj , qij, qijf 
  double precision :: rxi , ryi , rzi , rxij , ryij , rzij 
  double precision :: rxj , ryj , rzj 
  double precision :: fxij , fyij , fzij 
  double precision :: u_dir , vir_dir 
  double precision :: ttt1 , ttt2  , ttt3
  double precision :: magd_it, ene_it 
  double precision, dimension(:), allocatable :: fx_dir, fy_dir, fz_dir
  double precision, dimension(:), allocatable :: u_dir_magd


  ttt1 = MPI_WTIME(ierr)

  allocate( fx_dir(natm), fy_dir(natm), fz_dir(natm) )
  allocate( u_dir_magd (ntype) )

#ifdef debug
  print*,qia
#endif  

  vir_dir = 0.0d0
  u_dir   = 0.0d0 
 
  fx_dir  = 0.0D0
  fy_dir  = 0.0D0
  fz_dir  = 0.0D0

  u_dir_magd = 0.0d0
! =========================================================
!  MAIN LOOP calculate EFG(i) for each atom i parallelized
! =========================================================
atom : do ia = iastart , iaend

    rxi = rx(ia)
    ryi = ry(ia)
    rzi = rz(ia)
    qi  = qia(ia)

    ! ==================================================
    ! sum over neighboring cells (see direct_sum_init )
    ! ==================================================
    do ncell = 1 , rm_coul%ncmax
         ! ==============================
         !  ia and ja in different cells
         ! ==============================
         if ( rm_coul%lcell ( ncell ) .eq. 1) then
          do ja = 1 , natm

            qj    = qia(ja)

            rxj   = rx(ja) + rm_coul%boxxyz(1,ncell)
            ryj   = ry(ja) + rm_coul%boxxyz(2,ncell)
            rzj   = rz(ja) + rm_coul%boxxyz(3,ncell)

            rxij  = rxi - rxj
            ryij  = ryi - ryj
            rzij  = rzi - rzj

            rijsq = rxij * rxij + ryij * ryij + rzij * rzij

            rij   = dsqrt( rijsq )
            qij   = qi * qj / rij
            qijf  = qij / rijsq

            u_dir = u_dir + qij 
!	    print*,'out',atype(ia),atype(ja),qij,u_dir,ncell
            vir_dir = vir_dir - qij
            u_dir_magd(itype(ia)) = u_dir_magd(itype(ia)) + qij  ! for madlung

            fxij = qijf * rxij
            fyij = qijf * ryij
            fzij = qijf * rzij

            fx_dir(ia) = fx_dir(ia) + fxij
            fy_dir(ia) = fy_dir(ia) + fyij
            fz_dir(ia) = fz_dir(ia) + fzij

          enddo ! ja
 
        endif 
         ! =======================================
         !  ia and ja in the same cell (ia.ne.ja)
         ! =======================================
        if ( rm_coul%lcell(ncell) .eq. 0) then

          do ja = 1,natm

            if ( ja .ne. ia ) then

              qj = qia(ja)

              rxj  = rx(ja)
              ryj  = ry(ja)
              rzj  = rz(ja)

              rxij = rxi - rxj
              ryij = ryi - ryj
              rzij = rzi - rzj

              rijsq   = rxij * rxij + ryij * ryij + rzij * rzij
 
              rij = dsqrt( rijsq )
              qij = qi * qj / rij
              qijf = qij / rijsq

              u_dir = u_dir + qij 
!	      print*,'in',atype(ia),atype(ja),qij,u_dir,0
              vir_dir = vir_dir - qij
              u_dir_magd(itype(ia)) = u_dir_magd(itype(ia)) + qij  ! for madlung

              fxij = qijf * rxij
              fyij = qijf * ryij
              fzij = qijf * rzij

              fx_dir(ia) = fx_dir(ia) + fxij
              fy_dir(ia) = fy_dir(ia) + fyij
              fz_dir(ia) = fz_dir(ia) + fzij

            endif ! ia.ne.ja

          enddo ! ja

        endif 

     enddo ! ncell
  enddo atom

  do it=1,ntype
    magd_it=u_dir_magd(it)/(qch(it))
    ene_it=magd_it*qch(it)
    print*,'MAGD',it,magd_it,ene_it
  enddo

  ttt2 = MPI_WTIME(ierr)

  CALL MPI_ALL_REDUCE_DOUBLE_SCALAR ( u_dir )
  CALL MPI_ALL_REDUCE_DOUBLE_SCALAR ( vir_dir ) 
    
  CALL MPI_ALL_REDUCE_DOUBLE ( fx_dir , natm )
  CALL MPI_ALL_REDUCE_DOUBLE ( fy_dir , natm )
  CALL MPI_ALL_REDUCE_DOUBLE ( fz_dir , natm )
  
  u_coul   = u_dir 
  vir_coul = vir_dir

  fx = fx + fx_dir
  fy = fy + fy_dir
  fz = fz + fz_dir


  ! ================================================
  !  WARNING  Empirical coding !!!!!!! WARNING
  !  Where is the sign error !!! (tested on NaCl and compared to ewald and GULP)
  !  TODO:  test on another system
  ! ================================================
  u_coul   = - u_coul
  vir_coul = - vir_coul


  deallocate( fx_dir, fy_dir, fz_dir )
  deallocate( u_dir_magd )

  ttt3 = MPI_WTIME(ierr)

  return

END SUBROUTINE engforce_coulomb_DS

!*********************** SUBROUTINE field_charge_DS ***********************
!
! this subroutine calculates the electric field at atoms position from all the 
! points charges of the system. The Direct Summation is used
! 
! Note : same consruction as efg_DS
!                          
!  Ei,alpha  =  sum_j!=i  T^alpha_ij * q_j  
!
!  T^alpha_ij  = r_ij,alpha / r_ij^3 
!
!******************************************************************************

SUBROUTINE field_charge_DS ( iastart , iaend , Efield )

  USE config,	ONLY : natm , ntype , rx , ry, rz , qia
  USE control,  ONLY : cutlongrange
  USE io_file

  implicit none

  ! global 
  integer, intent(in) :: iastart , iaend 
  double precision :: Efield ( natm , 3) 

  ! local
  integer :: ia , ja  , it , ncell
  double precision :: d , d2 , d3 , dm3
  double precision :: rxi , ryi , rzi , rxij , ryij , rzij
  double precision :: rxj , ryj , rzj 
  double precision :: cutsq 
  logical :: lcharge

  lcharge = .false.
  do it = 1 , ntype
    if ( qch(it) .ne. 0.0d0 ) lcharge = .true.
  enddo

  if ( .not. lcharge ) WRITE ( stdout , '(a)') 'No point charges' 
  if ( .not. lcharge ) return

#ifdef debug
    print*,'rm_coul',rm_coul%ncmax
#endif

  Efield = 0.0d0
  cutsq = cutlongrange * cutlongrange

atom : do ia = iastart , iaend

         rxi = rx(ia)
         ryi = ry(ia)
         rzi = rz(ia)
#ifdef debug
    WRITE(*,'(a,i4,3f12.6)') 'pol pos ',ia,rxi,ryi,rzi
#endif
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
              d = dsqrt(d2)
              d3 = d2 * d
              dm3 = 1.0d0/ d3
              dm3 = dm3 * qia(ja)

              Efield ( ia , 1 ) = Efield ( ia , 1 ) + rxij * dm3 
              Efield ( ia , 2 ) = Efield ( ia , 2 ) + ryij * dm3 
              Efield ( ia , 3 ) = Efield ( ia , 3 ) + rzij * dm3 

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
                d = dsqrt(d2)
                d3 = d2 * d
                dm3 = 1.0d0/d3
                dm3 = dm3 * qia(ja)

                Efield ( ia , 1 ) = Efield ( ia , 1 ) + rxij * dm3 
                Efield ( ia , 2 ) = Efield ( ia , 2 ) + ryij * dm3 
                Efield ( ia , 3 ) = Efield ( ia , 3 ) + rzij * dm3 

              endif ! d2.lt.cutsq

            endif ! ia.ne.ja
          enddo ! ja
        endif

     enddo ! ncell

  enddo atom

  return

END SUBROUTINE field_charge_DS


!*********************** SUBROUTINE field_charge_ES ************************************
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
 
SUBROUTINE field_charge_ES ( iastart , iaend , Efield )

  USE control,    ONLY :  myrank , numprocs, calc
  USE config,     ONLY :  system , natm , natmi , atype , atypei , box , &
                          omega , itype , rx , ry , rz , ntype , qia 
  USE constants,  ONLY :  pi , fpi , piroot , imag
  USE kspace,     ONLY :  struc_fact
  USE prop,       ONLY :  nprop_print , nprop
  USE time
  USE io_file 

  implicit none

  INCLUDE 'mpif.h'

  ! global
  integer, intent(in) :: iastart , iaend 
  double precision :: Efield ( natm , 3 ) 

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
  double complex   :: rhon , carg , recarg_dgg 
  double precision :: recarg 
  double precision, external :: errfc 
  double precision :: ttt1 , ttt2 , ttt3 , ttt4 , ttt5 
  double precision :: Efield_real      ( natm , 3 )
  double precision :: Efield_surf      ( natm , 3 )
  double precision :: Efield_dual_real ( natm , 3 )
  double complex   :: Efield_dual      ( natm , 3 )
  logical :: lcharge

  lcharge = .false.
  do it = 1 , ntype
    if ( qch(it) .ne. 0.0d0 ) lcharge = .true.
  enddo
  
  if ( .not. lcharge ) WRITE ( stdout , '(a)') 'No point charges'
  if ( .not. lcharge ) return



#ifdef debug
  CALL print_config_sample(0,0)
  print*,'km_coul%nkcut',km_coul%nkcut
#endif
  ! ==========================
  !  init some quantities
  ! ==========================
  Efield      = 0.0d0
  Efield_real = 0.0d0
  Efield_dual = ( 0.0d0 , 0.0d0 ) 

  ! =================
  !  some constants 
  ! =================
  invbox = 1.0d0 / box
  alpha2 = alphaES * alphaES

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

     do ja = 1, natm

        rxij = rxi - rx(ja)
        ryij = ryi - ry(ja)
        rzij = rzi - rz(ja)
       if (ja .ne. ia ) then

         nxij = nint( rxij * invbox )
         nyij = nint( ryij * invbox )
         nzij = nint( rzij * invbox )
         rxij = rxij - box * nxij
         ryij = ryij - box * nyij
         rzij = rzij - box * nzij
  
         d2 = rxij * rxij + ryij * ryij + rzij * rzij
         d = dsqrt( d2 )
         d3 = d2 * d
         expon = dexp( - alpha2 * d2 ) / piroot

         T0 = errfc( alphaES * d ) / d3
         T1 = ( 2.0d0 * alphaES ) / d2 
         allrealpart = qia(ja) * ( T0 + T1*expon )

         Efield_real ( ia , 1 ) = Efield_real ( ia , 1 ) + rxij * allrealpart
         Efield_real ( ia , 2 ) = Efield_real ( ia , 2 ) + ryij * allrealpart
         Efield_real ( ia , 3 ) = Efield_real ( ia , 3 ) + rzij * allrealpart

       endif
         Efield_surf ( ia , 1 ) = Efield_surf ( ia , 1 ) - qia(ja) * rx(ja)
         Efield_surf ( ia , 2 ) = Efield_surf ( ia , 2 ) - qia(ja) * ry(ja)
         Efield_surf ( ia , 3 ) = Efield_surf ( ia , 3 ) - qia(ja) * rz(ja)

     enddo

  enddo atom1

  ttt3 = MPI_WTIME(ierr)
!  efgtimetot1 = efgtimetot1 + (ttt3-ttt2)

! ==============================================
!            reciprocal space part
! ==============================================
  ! sum on atoms 
atom2:  do ia = 1 , natm
    rxi = rx ( ia )
    ryi = ry ( ia )
    rzi = rz ( ia )
    kpoint : do ik = 1, km_coul%nkcut 
      ! =================
      !   k-space  
      ! =================
      kx = km_coul%kpt ( 1 , ik )
      ky = km_coul%kpt ( 2 , ik )
      kz = km_coul%kpt ( 3 , ik )
      kk = km_coul%kptk( ik )
      Ak = dexp( - kk * 0.25d0 / alpha2 ) 
      !write(*, '(5f12.8)')  kx , ky , kz , Ak , kk
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
        rhon = rhon + imag * qch(it) * CONJG( km_coul%strf ( ik , it ) )
      enddo

      kri = ( kx * rxi + ky * ryi + kz * rzi ) 
      carg = EXP ( imag * kri )
      recarg_dgg =  rhon * carg * Ak / kk 
     ! write(*, '(8f12.8)') recarg_dgg,rhon,carg,Ak,kk 

      Efield_dual ( ia , 1 )  =  Efield_dual ( ia , 1 ) - kx * recarg_dgg
      Efield_dual ( ia , 2 )  =  Efield_dual ( ia , 2 ) - ky * recarg_dgg
      Efield_dual ( ia , 3 )  =  Efield_dual ( ia , 3 ) - kz * recarg_dgg

    enddo kpoint

  enddo atom2

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
  Efield_dual_real( : , : ) = dble( Efield_dual( : , :)  ) 
  Efield_dual_real    =  Efield_dual_real * fpi / omega 
  Efield_surf = Efield_surf * fpi / omega / 3.0D0
  write(*,'(a,3f12.8)') ' Efield_real         = ',Efield_real      ( 1 , 1 ) , Efield_real      ( 1 , 2 ) , Efield_real ( 1 , 3 ) 
  write(*,'(a,3f12.8)') ' Efield_dual         = ',Efield_dual_real ( 1 , 1 ) , Efield_dual_real ( 1 , 2 ) , Efield_dual_real ( 1 , 3 ) 
  write(*,'(a,3f12.8)') ' Efield_surf         = ',Efield_surf      ( 1 , 1 ) , Efield_surf      ( 1 , 2 ) , Efield_surf      ( 1 , 3 ) 

  ! =======================
  !  total electric field 
  ! =======================
  Efield  = Efield_dual_real + Efield_real 
  write(*,'(a,3f12.8)') ' Efield_without_surf = ',   Efield   ( 1 , 1 ) , Efield  ( 1 , 2 ) , Efield  ( 1 , 3 )  
  Efield  = Efield_dual_real + Efield_real + Efield_surf
  write(*,'(a,3f12.8)') ' Efield_with_surf    = ',   Efield   ( 1 , 1 ) , Efield  ( 1 , 2 ) , Efield  ( 1 , 3 )  

  ttt5 = MPI_WTIME(ierr)
!  efgtimetot3 = efgtimetot3 + (ttt5-ttt4)

  CALL MPI_BARRIER( MPI_COMM_WORLD , ierr )

  return

END SUBROUTINE field_charge_ES

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

SUBROUTINE field_dipole_DS ( iastart , iaend , Efield , mu , linduced )

  USE config,	ONLY : natm , ntype , rx , ry , rz , dipia
  USE control,  ONLY : cutlongrange
  USE io_file

  implicit none

  ! global 
  integer, intent(in) :: iastart , iaend 
  double precision :: Efield ( natm , 3) 
  double precision :: mu    ( natm , 3) 
  logical :: linduced

  ! local
  integer :: ia , ja  , it , ncell
  double precision :: d , d2 , d4 , dm5 
  double precision :: Txx , Tyy , Tzz , Txy , Txz , Tyz
  double precision :: rxi , ryi , rzi , rxij , ryij , rzij
  double precision :: rxj , ryj , rzj 
  double precision :: cutsq
  logical :: ldipole

  if ( .not. linduced ) then
    ldipole = .false.
    do it = 1 , ntype 
      if ( dipia(it,1) .ne. 0.0d0 .or. & 
           dipia(it,2) .ne. 0.0d0 .or. &
           dipia(it,3) .ne. 0.0d0 ) ldipole = .true.
    enddo
    if ( .not. ldipole .and. ionode ) WRITE ( stdout , '(a)') 'No static dipoles' 
    if ( .not. ldipole ) return
    mu = dipia
  else
#ifdef debug2
   WRITE ( stdout , '(a)') 'induced moment in input will be used'
#endif 
  endif

#ifdef debug
    print*,'rm_coul',rm_coul%ncmax
#endif
 
  cutsq = cutlongrange * cutlongrange
  Efield = 0.0d0
 

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
              d = dsqrt(d2)
              d4 = d2 * d2
              dm5 = 1.0d0/ ( d4 * d )

              Txx = 3.0d0 * rxij * rxij - d2 
              Tyy = 3.0d0 * ryij * ryij - d2 
              Tzz = 3.0d0 * rzij * rzij - d2 
              Txy = 3.0d0 * rxij * ryij  
              Txz = 3.0d0 * rxij * rzij  
              Tyz = 3.0d0 * ryij * rzij  

              Efield ( ia , 1 ) = Efield ( ia , 1 ) + ( Txx * mu ( ja , 1 ) + Txy * mu ( ja , 2 ) + Txz * mu ( ja , 3 ) ) * dm5
              Efield ( ia , 2 ) = Efield ( ia , 2 ) + ( Txy * mu ( ja , 1 ) + Tyy * mu ( ja , 2 ) + Tyz * mu ( ja , 3 ) ) * dm5
              Efield ( ia , 3 ) = Efield ( ia , 3 ) + ( Txz * mu ( ja , 1 ) + Tyz * mu ( ja , 2 ) + Tzz * mu ( ja , 3 ) ) * dm5


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
                d = dsqrt(d2)
                d4 = d2 * d2
                dm5 = 1.0d0/ ( d4 * d )
  
                Txx = 3.0d0 * rxij * rxij - d2 
                Tyy = 3.0d0 * ryij * ryij - d2 
                Tzz = 3.0d0 * rzij * rzij - d2 
                Txy = 3.0d0 * rxij * ryij  
                Txz = 3.0d0 * rxij * rzij  
                Tyz = 3.0d0 * ryij * rzij  

                Efield ( ia , 1 ) = Efield ( ia , 1 ) + ( Txx * mu ( ja , 1 ) + Txy * mu ( ja , 2 ) + Txz * mu ( ja , 3 ) ) * dm5
                Efield ( ia , 2 ) = Efield ( ia , 2 ) + ( Txy * mu ( ja , 1 ) + Tyy * mu ( ja , 2 ) + Tyz * mu ( ja , 3 ) ) * dm5
                Efield ( ia , 3 ) = Efield ( ia , 3 ) + ( Txz * mu ( ja , 1 ) + Tyz * mu ( ja , 2 ) + Tzz * mu ( ja , 3 ) ) * dm5

              endif ! d2.lt.cutsq

            endif ! ia.ne.ja
          enddo ! ja
        endif

     enddo ! ncell

  enddo atom

  return

END SUBROUTINE field_dipole_DS

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
 
SUBROUTINE field_dipole_ES ( iastart , iaend , Efield , mu , linduced )

  USE control,    ONLY :  myrank , numprocs, calc
  USE config,     ONLY :  system , natm , natmi , atype , atypei , box , &
                          omega , itype , rx , ry , rz , ntype , qia , dipia 
  USE constants,  ONLY :  pi , fpi , piroot , imag
  USE kspace,     ONLY :  struc_fact , struc_fact_dip
  USE prop,       ONLY :  nprop_print , nprop
  USE time
  USE io_file 

  implicit none

  INCLUDE 'mpif.h'

  ! global
  integer, intent(in) :: iastart , iaend 
  double precision :: Efield ( natm , 3 ) 
  double precision :: mu    ( natm , 3 ) 
  logical :: linduced

  ! local
  integer :: ia , ja ,  ierr , it , ip 
  integer :: nxij , nyij , nzij
  double precision :: mip    ( ntype , 3 ) 
  double precision :: mutot  ( 3 ) 
  double precision :: d , d2 , d4 , d3 , d5 , expon 
  double precision :: alpha2 , alpha3
  double precision :: allrealpart
  double precision :: T0 , T1 , T2  ! real part 
  double precision :: rxi , ryi , rzi , rxij , ryij , rzij 
  double precision :: invbox 
  double precision :: ak, kx, ky, kz, kk
  double precision :: Txx , Tyy , Tzz , Txy , Txz , Tyz
  double complex   :: Tx , Ty , Tz , sTx 
  !double precision   :: Tx , Ty , Tz , sTx
  integer :: ik
  double precision :: kri , kmu
  double complex   :: rhon , carg , recarg_dgg 
  double precision :: recarg 
  double precision :: cdual  
  double precision, external :: errfc 
  double precision :: ttt1 , ttt2 , ttt3 , ttt4 , ttt5 
  double precision :: Efield_real      ( natm , 3 )
  double precision :: Efield_surf      ( natm , 3 )
  double precision :: Efield_dual_real ( natm , 3 )
  double complex   :: Efield_dual      ( natm , 3 )
  logical :: ldipole

  if ( .not. linduced ) then
    ldipole = .false.
    do it = 1 , ntype
      if ( dipia(it,1) .ne. 0.0d0 .or. &
           dipia(it,2) .ne. 0.0d0 .or. &
           dipia(it,3) .ne. 0.0d0 ) ldipole = .true.
    enddo
    if ( .not. ldipole .and. ionode ) WRITE ( stdout , '(a)') 'No static dipoles'
    if ( .not. ldipole ) return
    mu = dipia
  else
#ifdef debug2
   WRITE ( stdout , '(a)') 'induced moment in input will be used'
#endif 
  endif
sTx = 0.0d0
  ! ==============================
  !  init mip ( itype dependent ) 
  ! ==============================
  ip = 0
  do it = 1, ntype
    ip = ip + natmi ( it )
    mip ( it , 1 ) = mu ( ip , 1 )
    mip ( it , 2 ) = mu ( ip , 2 )
    mip ( it , 3 ) = mu ( ip , 3 )
  enddo

  do ia = 1 , natm
    mutot ( 1 ) = mutot ( 1 ) + mu ( ia , 1 )
    mutot ( 2 ) = mutot ( 2 ) + mu ( ia , 2 )
    mutot ( 3 ) = mutot ( 3 ) + mu ( ia , 3 )
  enddo

#ifdef debug
  CALL print_config_sample(0,0)
  print*,'km_coul%nkcut',km_coul%nkcut
#endif
  ! ==========================
  !  init some quantities
  ! ==========================
  Efield           = 0.0d0
  Efield_surf      = 0.0d0
  Efield_real      = 0.0d0
  Efield_dual      = ( 0.0d0 , 0.0d0 ) 
  Efield_dual_real = 0.0d0

  ! =================
  !  some constants 
  ! =================
  invbox = 1.0d0 / box
  alpha2 = alphaES * alphaES
  alpha3 = alpha2  * alphaES
  cdual  = omega * alpha3 / 3.0d0 / piroot / pi 

  ttt1 = MPI_WTIME(ierr)
  ! =====================
  ! facteur de structure 
  ! =====================
  CALL struc_fact ( km_coul )
!  CALL struc_fact ( km_coul_dip )

  do ia= 1,natm
    write ( stdout ,'(3f16.8)')  mu  ( ia , 1 ) , mu  ( ia , 2 ) ,  mu  ( ia , 3 ) 
  enddo
  do it=1,ntype
    write ( stdout ,'(3f16.8)')  mip ( it , 1 ) , mip ( it , 2 ) ,  mip ( it , 3 )
  enddo

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

        rxij = rxi - rx(ja)
        ryij = ryi - ry(ja)
        rzij = rzi - rz(ja)
       if (ja .ne. ia ) then

         nxij = nint( rxij * invbox )
         nyij = nint( ryij * invbox )
         nzij = nint( rzij * invbox )
         rxij = rxij - box * nxij
         ryij = ryij - box * nyij
         rzij = rzij - box * nzij
  
         d2 = rxij * rxij + ryij * ryij + rzij * rzij
         d = dsqrt( d2 )
         d4 = d2 * d2
         d5 = d4 * d

         expon = dexp ( - alpha2 * d2   ) / piroot
         T0 = errfc   ( alphaES * d     ) / d5
         T1 =         ( 2.0d0 * alphaES ) / d4
         T2 =         ( 4.0d0 * alpha3  ) / d2 / 3.0d0
         allrealpart = ( T0 + ( T1 + T2 ) * expon )
      
         Txx = 3.0d0 * rxij * rxij - d2
         Tyy = 3.0d0 * ryij * ryij - d2
         Tzz = 3.0d0 * rzij * rzij - d2
         Txy = 3.0d0 * rxij * ryij
         Txz = 3.0d0 * rxij * rzij
         Tyz = 3.0d0 * ryij * rzij

         Efield_real ( ia , 1 ) = Efield_real ( ia , 1 ) + &
         ( Txx * mu ( ja , 1 ) + Txy * mu ( ja , 2 ) + Txz * mu ( ja , 3 ) ) * allrealpart 
         Efield_real ( ia , 2 ) = Efield_real ( ia , 2 ) + &
         ( Txy * mu ( ja , 1 ) + Tyy * mu ( ja , 2 ) + Tyz * mu ( ja , 3 ) ) * allrealpart
         Efield_real ( ia , 3 ) = Efield_real ( ia , 3 ) + &
         ( Txz * mu ( ja , 1 ) + Tyz * mu ( ja , 2 ) + Tzz * mu ( ja , 3 ) ) * allrealpart


       endif

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
  ! sum on atoms 
atom2:  do ia = 1 , natm
    rxi = rx ( ia )
    ryi = ry ( ia )
    rzi = rz ( ia )
    sTx = (0.0d0,0.0d0)
    kpoint : do ik = 1, km_coul%nkcut 
!    kpoint : do ik = 1, km_coul_dip%nkcut 
      ! =================
      !   k-space  
      ! =================
      kx = km_coul%kpt ( 1 , ik )
      ky = km_coul%kpt ( 2 , ik )
      kz = km_coul%kpt ( 3 , ik )
      kk = km_coul%kptk( ik )
!      kx = km_coul_dip%kpt ( 1 , ik )
!      ky = km_coul_dip%kpt ( 2 , ik )
!      kz = km_coul_dip%kpt ( 3 , ik )
!      kk = km_coul_dip%kptk( ik )
      Ak = dexp( - kk * 0.25d0 / alpha2 ) 
      Ak = Ak / kk  
      !write(*, '(5f12.8)')  kx , ky , kz , Ak , kk
      if ( km_coul%kptk( ik ) .eq. 0 ) then 
        WRITE ( stdout , * ) 'the sum should be done on k! =  0',ik
        STOP 
      endif
!      if ( km_coul_dip%kptk( ik ) .eq. 0 ) then 
!        WRITE ( stdout , * ) 'the sum should be done on k! =  0',ik
!        STOP 
!      endif

      ! ===============================
      !                              ---
      !  charge density in k-space ( \   q * facteur de structure  )
      !                              /__
      ! ===============================
      rhon = (0.d0, 0.d0)
      do it = 1, ntype
        rhon = rhon + CONJG ( imag * ( mip ( it , 1 ) * kx + mip ( it , 2 ) * ky + mip ( it , 3 ) * kz ) &
             * km_coul%strf ( ik , it ) ) 
!        rhon = rhon + CONJG ( imag * ( mip ( it , 1 ) * kx + mip ( it , 2 ) * ky + mip ( it , 3 ) * kz ) &
!             * km_coul_dip%strf ( ik , it ) ) 
      enddo

      kri        = ( kx * rxi + ky * ryi + kz * rzi ) 
      carg       = EXP ( imag * kri ) 
      recarg_dgg =  imag * rhon * carg 

      kmu = kx * mu ( ia , 1 ) + ky * mu ( ia , 2 ) + kz * mu ( ia , 3 )
!      kmu = kmu / 3.0D0
!      Tx  = Ak * kx * ( dble(recarg_dgg) - kmu ) 
!      Ty  = Ak * ky * ( dble(recarg_dgg) - kmu )
!      Tz  = Ak * kz * ( dble(recarg_dgg) - kmu )

      Tx  = Ak * kx * ( recarg_dgg - kmu ) 
      Ty  = Ak * ky * ( recarg_dgg - kmu )
      Tz  = Ak * kz * ( recarg_dgg - kmu )

!      Tx  = Ak * kx * ( dble ( recarg_dgg ))
!      Ty  = Ak * ky * ( dble ( recarg_dgg ))
!      Tz  = Ak * kz * ( dble ( recarg_dgg ))

!      Tx  = Ak * kx * recarg_dgg 
!      Ty  = Ak * ky * recarg_dgg 
!      Tz  = Ak * kz * recarg_dgg 

!      Tx  = Ak * kx * recarg_dgg - kmu
!     Ty  = Ak * ky * recarg_dgg - kmu 
!     Tz  = Ak * kz * recarg_dgg - kmu

      sTx = STx + Tx
  !    print*,sTx,kmu

      Efield_dual ( ia , 1 ) = Efield_dual ( ia , 1 ) - Tx 
      Efield_dual ( ia , 2 ) = Efield_dual ( ia , 2 ) - Ty 
      Efield_dual ( ia , 3 ) = Efield_dual ( ia , 3 ) - Tz 
!      write(*, '(2i5,6f14.8)') ia,ik,rhon,carg,Ak,kk 
!      write(*, '(2i5,6f14.8)') ia,ik, Tx , Ty , Tz
!      write(*, '(2i5,6f14.8)') ia,ik,Efield_dual ( ia , 1 ),Efield_dual ( ia , 2 ),Efield_dual ( ia , 3 )

    enddo kpoint

!    Efield_dual ( ia , 1 ) = Efield_dual ( ia , 1 ) - ( mutot ( 1 ) - mu( ia , 1 ) ) / 3.0d0 !+ cdual * mu ( ia , 1 )
!    Efield_dual ( ia , 2 ) = Efield_dual ( ia , 2 ) - ( mutot ( 2 ) - mu( ia , 2 ) ) / 3.0d0 !+ cdual * mu ( ia , 2 )
!    Efield_dual ( ia , 3 ) = Efield_dual ( ia , 3 ) - ( mutot ( 3 ) - mu( ia , 3 ) ) / 3.0d0 !+ cdual * mu ( ia , 3 )

!    Efield_dual ( ia , 1 ) = Efield_dual ( ia , 1 ) - ( mutot ( 1 ) - mu( ia , 1 ) ) / 3.0d0 + cdual * mu ( ia , 1 )
!    Efield_dual ( ia , 2 ) = Efield_dual ( ia , 2 ) - ( mutot ( 2 ) - mu( ia , 2 ) ) / 3.0d0 + cdual * mu ( ia , 2 )
!    Efield_dual ( ia , 3 ) = Efield_dual ( ia , 3 ) - ( mutot ( 3 ) - mu( ia , 3 ) ) / 3.0d0 + cdual * mu ( ia , 3 )

    Efield_dual ( ia , 1 ) = Efield_dual ( ia , 1 ) + mu( ia , 1 ) / 6.0d0 
    Efield_dual ( ia , 2 ) = Efield_dual ( ia , 2 ) + mu( ia , 2 ) / 6.0d0 
    Efield_dual ( ia , 3 ) = Efield_dual ( ia , 3 ) + mu( ia , 3 ) / 6.0d0 

!    Efield_dual ( ia , 1 ) = Efield_dual ( ia , 1 ) + mu( ia , 1 ) / 3.0d0 + cdual * mu ( ia , 1 )
!    Efield_dual ( ia , 2 ) = Efield_dual ( ia , 2 ) + mu( ia , 2 ) / 3.0d0 + cdual * mu ( ia , 2 )
!    Efield_dual ( ia , 3 ) = Efield_dual ( ia , 3 ) + mu( ia , 3 ) / 3.0d0 + cdual * mu ( ia , 3 )

  enddo atom2

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
  Efield_dual_real ( : , : ) = dble( Efield_dual ( : , : )  ) 
  Efield_dual_real           = Efield_dual_real  * 2.0D0  * fpi / omega 
!  Efield_dual_real           = Efield_dual_real  * fpi / omega 
  Efield_surf                = Efield_surf       * fpi / omega 
  do ip=1,natm
    WRITE ( stdout ,'(a,i4,a,3f12.8)') 'atom = ',ip,' Efield_real         = ', &
    Efield_real      ( ip , 1 ) , Efield_real      ( ip , 2 ) , Efield_real      ( ip , 3 ) 
  enddo
  do ip=1,natm
    WRITE ( stdout ,'(a,i4,a,3f12.8)') 'atom = ',ip,' Efield_dual         = ', &
    Efield_dual_real ( ip , 1 ) , Efield_dual_real ( ip , 2 ) , Efield_dual_real ( ip , 3 ) 
  enddo
  do ip=1,natm
    WRITE ( stdout ,'(a,i4,a,3f12.8)') 'atom = ',ip,' Efield_surf         = ', &
    Efield_surf      ( ip , 1 ) , Efield_surf      ( ip , 2 ) , Efield_surf      ( ip , 3 ) 
  enddo
  ! =======================
  !  total electric field 
  ! =======================
  Efield  = Efield_dual_real + Efield_real 
  do ip=1,natm
    WRITE ( stdout ,'(a,i4,a,3f12.8)') 'atom = ',ip,' Efield_without_surf = ',   Efield   ( ip , 1 ) , Efield  ( ip , 2 ) , Efield  ( ip , 3 )   
  enddo
  Efield  = Efield_dual_real + Efield_real + Efield_surf
  do ip=1,natm
    WRITE ( stdout ,'(a,i4,a,3f12.8)') 'atom = ',ip,' Efield_with_surf    = ',   Efield   ( ip , 1 ) , Efield  ( ip , 2 ) , Efield  ( ip , 3 )  
  enddo

  ttt5 = MPI_WTIME(ierr)
!  efgtimetot3 = efgtimetot3 + (ttt5-ttt4)

  CALL MPI_BARRIER( MPI_COMM_WORLD , ierr )

  return

END SUBROUTINE field_dipole_ES

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

SUBROUTINE induced_moment ( Efield , mu_ind )

  USE config, 	ONLY : natm

  implicit none

  ! global
  double precision :: Efield ( natm , 3 ) 
  double precision :: mu_ind ( natm , 3 ) 
 
  ! local 
  integer :: alpha , beta
  integer ia

  mu_ind = 0.0d0
  do ia = 1 , natm 
    do alpha = 1 , 3 
      do beta = 1 , 3  
        mu_ind ( ia , alpha ) = mu_ind ( ia , alpha ) + pol ( ia , alpha , beta ) * Efield ( ia , beta )  
      enddo
    enddo
  enddo

  return

END SUBROUTINE induced_moment


!*********************** SUBROUTINE Efield_and_scf_induced *********************
!
! this subroutine is an interface to the calculation of the electric field from 
! point charges and dipoles. it also calculate the induced moment in SCF procedure.
! 
!  if lfield_only = .TRUE. : we calculate only electric field 	
! 
!******************************************************************************

SUBROUTINE Efield_and_scf_induced ( iastart, iaend , mu_ind , lfield_only ) 

  USE io_file , ONLY : stdout , ionode
  USE config, 	ONLY : natm , natmi , ntype , atype
  USE control,  ONLY : longrange

  implicit none

  ! global 
  ! global
  integer, intent(in)           :: iastart , iaend 
  double precision :: mu_ind ( natm , 3 ) 
  logical :: lfield_only

  ! local 
  integer :: ia , iscf , it , ik , npol
  double precision :: efield_old , diff_efield
  double precision :: Efield ( natm , 3 ) , Efield_stat ( natm , 3 ) , Efield_ind ( natm , 3 ) , tmp( natm , 3 ) 
  logical :: linduced

  ! =============================================
  !  calculate static Efield ( charge + dipoles)
  ! =============================================
  Efield_stat = 0.0d0

  ! =============================================
  !  Electric Field from point charges
  ! =============================================
  if ( longrange .eq. 'direct' ) CALL field_charge_DS ( iastart , iaend , tmp )
  if ( longrange .eq. 'ewald' )  CALL field_charge_ES ( iastart , iaend , tmp )

  Efield_stat = Efield_stat + tmp
  tmp = 0.0d0
  mu_ind = 0.0d0
  ! =============================================
  !  Electric Field from (static) dipoles i.e linduced=.FALSE.
  ! =============================================
  if ( longrange .eq. 'direct' ) CALL field_dipole_DS ( iastart , iaend , tmp , mu_ind , linduced=.FALSE. ) 
  if ( longrange .eq. 'ewald' )  CALL field_dipole_ES ( iastart , iaend , tmp , mu_ind , linduced=.FALSE. )
  Efield_stat = Efield_stat + tmp
  mu_ind = 0.0d0

  ! =============================================
  !  init total Efield to static only 
  ! =============================================
  Efield = Efield_stat 

  linduced = .false.
  do it = 1 , ntype
    if ( lpolar ( it ) .eq. 1 ) linduced = .true.
  enddo
  
  if ( linduced .and. .not. lfield_only ) then
  if ( ionode ) WRITE ( stdout , '(a)' ) 'calculate induced dipole SCF'

  
    if ( ionode ) then
      WRITE ( stdout , '(a)' ) 'We start from the static electric field (charge + static dipole)'
      WRITE ( stdout , '(a)' ) ' '
    endif
    diff_efield = 1000 
    Efield_old  = Efield(1,1)
    iscf = 0
    ! =============================
    !           SCF LOOP
    ! =============================
    do while ( diff_efield .gt. 1e-8 .or. iscf .lt. 2 ) 
    
      iscf = iscf + 1 

      ! ==========================================================
      !  calculate mu_ind from Efield = Efield_stat + Efield_ind
      ! ==========================================================
      CALL induced_moment ( Efield , mu_ind )  ! Efield in ; mu_ind out 
  
      ! ==========================================================
      !  calculate Efield_ind from mu_ind 
      !  Efield_ind out , mu_ind in
      ! ==========================================================
      if ( longrange .eq. 'direct' ) CALL field_dipole_DS ( iastart , iaend , Efield_ind , mu_ind , linduced = .TRUE. )
      if ( longrange .eq. 'ewald' )  CALL field_dipole_ES ( iastart , iaend , Efield_ind , mu_ind , linduced = .TRUE. )

      Efield = Efield_stat + Efield_ind

      ! ===================
      !  stopping criteria
      ! ===================
      diff_efield = dabs ( Efield(1,1) - Efield_old ) 
      Efield_old = Efield(1,1)
  
      if ( ionode ) WRITE ( stdout ,'(a,i4,a,3f18.10,a,f18.10)') &
      'scf = ',iscf,' Efield atom 1= ',Efield(1,1),Efield(1,2),Efield(1,3),' conv = ',diff_efield

    enddo ! end of SCF loop

    dip_ind = mu_ind

    if ( ionode ) then
      WRITE ( stdout , '(a)' ) ' '
      WRITE ( stdout , '(a,i6,a)') 'scf calculation of the induced electric moment converged in ',iscf, ' iterations '
      WRITE ( stdout , '(a,i6,a)') 'Electric field is converged within 1e-6'
      WRITE ( stdout , '(a)' ) ' '
      WRITE ( stdout , '(a)' )     'Induced dipoles at atoms : '
      do ia = 1 , natm
        WRITE ( stdout , '(i5,a3,a,3f18.10)' ) &
        ia,atype(ia),' mu_ind = ', dip_ind ( ia , 1 ) , dip_ind ( ia , 2 ) , dip_ind ( ia , 3 )
      enddo
      WRITE ( stdout , '(a)' ) ' '
    endif

  endif

  if ( ionode ) then
    WRITE ( stdout , '(a)' ) ' '
    WRITE ( stdout , '(a)' )     'Electric field at atoms : '
    do ia = 1 , natm
      WRITE ( stdout , '(i5,a3,a,3f18.10)' ) &
      ia,atype(ia),' Efield = ', Efield ( ia , 1)  , Efield ( ia , 2 ) , Efield ( ia , 3 )
    enddo
    WRITE ( stdout , '(a)' ) ' '
  endif

  return

END SUBROUTINE Efield_and_scf_induced 


END MODULE field 
! ===== fmV =====
