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
! ======= Hardware =======

MODULE field 

  USE config , ONLY : ntypemax
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

  TYPE ( kmesh ) :: km_coul
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
  qch(1)        = -1.0d0            
  qch(2)        =  1.0d0 

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
    xn    ( 1 )     = 0.8d0
    xn    ( 2 )     = 0.2d0
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
                         qch            
                          
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

  if ( ionode ) then
               WRITE ( kunit ,'(a)')       '=============================================================' 
               WRITE ( kunit ,'(a)')       ''
    if ( .not. lcoulomb )  &
               WRITE ( kunit ,'(a)')       'discret charges: (no forces!)'
             do it = 1 , ntype 
               WRITE ( kunit ,'(a,a,a,f10.5)') 'q',atypei(it),'      = ',qch(it) 
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
               WRITE ( kunit ,'(a)')       'ewald summation'
               WRITE ( kunit ,'(a,f10.5)') 'alpha                            = ',alphaES
               WRITE ( kunit ,'(a,i10)')   'ncellewald x ncellewald x ncellewald               = ',ncellewald
               WRITE ( kunit ,'(a,i10)')   '' 
               WRITE ( kunit ,'(a,f10.5)') 'Note:this should hold alpha^2 * box^2 >> 1',alphaES*alphaES*box*box
      endif
    endif

    ! RETURN ?
    if ( calc .eq. 'efg' ) return 
    
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
              WRITE ( kunit ,'(a)') '' 
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
  if (ionode.and..not.lreduced) write( stdout     ,'(a,2f20.9)') 'long range correction : ',utail
  if (ionode.and.lreduced)      write( stdout     ,'(a,2f20.9)') 'long range correction : ',utail/natm
  if (ionode.and..not.lreduced) write( kunit_OUTFF,'(a,2f20.9)') 'long range correction : ',utail
  if (ionode.and.lreduced)      write( kunit_OUTFF,'(a,2f20.9)') 'long range correction : ',utail/natm
  if (ionode) write( stdout     ,'(a)') ' '
  if (ionode) write( kunit_OUTFF,'(a)') ' '
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
      if ( lshiftpot .and. .not. lminimg ) CALL engforce_new_bmlj         ( iastart , iaend )
      if ( lshiftpot .and.       lminimg ) CALL engforce_bmlj_pbc         ( iastart , iaend )!, list , point )
      if ( .not.lshiftpot )                CALL engforce_bmlj_pbc_noshift ( iastart , iaend )!, list , point )
      if ( longrange .eq. 'ewald' )        CALL engforce_coulomb_ES       (  )
      if ( longrange .eq. 'direct')        CALL engforce_coulomb_DS       ( iastart , iaend ) 
    ! =======
    !  NO PBC
    ! =======
    else
      WRITE ( stdout ,'(a)') 'not yet : coulomb + nopbc'
      STOP 
      if ( lshiftpot )               CALL engforce_bmlj_nopbc       ( iastart , iaend )!, list , point )
!      if ( .not.lshiftpot )          CALL engforce_bmlj_pbc_noshift ( iastart , iaend , list , point ) 
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
      if ( lshiftpot       )     CALL engforce_bmlj_pbc         ( iastart , iaend )!, list , point )
      if ( .not. lshiftpot )     CALL engforce_bmlj_pbc_noshift ( iastart , iaend )!, list , point )
    ! =======
    !  NO PBC
    ! =======
    else
#ifdef debug  
    WRITE ( stdout , '(a)' ) ' no pbc ' 
#endif  
      if ( lshiftpot )          CALL engforce_bmlj_nopbc          ( iastart , iaend )! , list , point )
!      if ( .not. lshiftpot )   CALL engforce_bmlj_nopbc_noshift ( iastart , iaend , list , point ) 
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
!  integer, intent(out) :: list( 250 * natm ) , point( natm+1 )

  ! local
  integer :: ia , ja , it , jt , j1 , je , jb , ierr
  integer :: p1 , p2
  integer :: nxij , nyij , nzij
  double precision :: rxi , ryi , rzi
  double precision :: rxij , ryij , rzij , sr2 , rijsq , srp , srq
  double precision :: wij , fxij , fyij , fzij
  double precision :: ptwo ( ntype , ntype )
  double precision :: qtwo ( ntype , ntype )
  double precision :: one13 , two13 , forcetime1 , forcetime2 , invbox
  double precision :: pot_sum , vir_sum
  double precision :: u , vir 
  double precision, dimension(:), allocatable :: fx_sum , fy_sum , fz_sum

#ifdef debug
  print*,'atype',atype
  print*,'itype',itype
#endif

  forcetime1 = MPI_WTIME(ierr) ! timing info
  
  allocate( fx_sum (natm), fy_sum (natm), fz_sum (natm) )

  fx_sum = 0.0D0
  fy_sum = 0.0D0
  fz_sum = 0.0D0

  u = 0.0D0
  vir = 0.0D0
  fx = 0.0D0
  fy = 0.0D0
  fz = 0.0D0

  invbox = 1.0d0 / box

  one13 = ( 1.0D0 / 3.0D0 )
  two13 = 2.D0 ** one13

  do jt = 1 , ntype
    do it = 1 , ntype
      ptwo ( it , jt ) = plj ( it , jt ) * 0.5d0
      qtwo ( it , jt ) = qlj ( it , jt ) * 0.5d0
    enddo
  enddo

  if ( lvnlist ) then
    CALL vnlistcheck ( iastart , iaend )!, list , point )
    do ia = iastart , iaend
      rxi = rx ( ia )
      ryi = ry ( ia )
      rzi = rz ( ia )
      jb = point( ia )
      je = point( ia + 1 ) - 1
      do j1 = jb, je
        ja = list ( j1 )
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
    else ! lvnlist .false.
      do ia = iastart , iaend
        rxi = rx ( ia )
        ryi = ry ( ia )
        rzi = rz ( ia )
        do ja =  1 , natm
          if ( ja .gt. ia ) then    
              rxij = rxi - rx ( ja )
              ryij = ryi - ry ( ja )
              rzij = rzi - rz ( ja )
              nxij = nint( rxij * invbox )
              nyij = nint( ryij * invbox )
              nzij = nint( rzij * invbox )
              rxij = rxij - box * nxij
              ryij = ryij - box * nyij
              rzij = rzij - box * nzij
              rijsq = rxij  * rxij + ryij * ryij + rzij * rzij
              p1 = itype ( ja )
              p2 = itype ( ia )
              if ( rijsq .lt. rcutsq(p1,p2) ) then
                sr2 = sigsq(p1,p2)/rijsq
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
    endif    

  forcetime2 = MPI_WTIME(ierr) ! timing info
  forcetimetot = forcetimetot+(forcetime2-forcetime1)

  pot_sum = 0.0d0
  vir_sum = 0.0d0 
  CALL MPI_ALLREDUCE(u,pot_sum,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
  CALL MPI_ALLREDUCE(vir,vir_sum,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)

  CALL MPI_ALLREDUCE(fx,fx_sum,natm,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
  CALL MPI_ALLREDUCE(fy,fy_sum,natm,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
  CALL MPI_ALLREDUCE(fz,fz_sum,natm,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)

  fx = fx_sum
  fy = fy_sum
  fz = fz_sum
  fx_sum = 0.0d0
  fy_sum = 0.0d0
  fz_sum = 0.0d0

  u_lj = pot_sum
  vir_lj = vir_sum

  deallocate( fx_sum, fy_sum, fz_sum )

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
  integer :: ia , ja , it, jt, j1, je, jb , ierr
  integer :: p1, p2
  double precision :: rxi , ryi , rzi
  double precision :: rxij , ryij , rzij , sr2 , rijsq , srp , srq
  double precision :: wij , fxij , fyij , fzij
  double precision :: ptwo ( ntype , ntype )
  double precision :: qtwo ( ntype , ntype )
  double precision :: one13, two13
  double precision :: pot_sum, vir_sum
  double precision, dimension(:), allocatable :: fx_sum, fy_sum, fz_sum
  double precision :: u, vir


  allocate( fx_sum(natm), fy_sum(natm), fz_sum(natm) )

  fx_sum = 0.0D0
  fy_sum = 0.0D0
  fz_sum = 0.0D0

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

  if ( lvnlist ) then
    CALL vnlistcheck ( iastart , iaend ) !, list , point )
    do ia = iastart, iaend
      rxi = rx ( ia )
      ryi = ry ( ia )
      rzi = rz ( ia )
      jb = point ( ia )
      je = point ( ia + 1 ) - 1
      do j1 = jb, je
        ja = list(j1)
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

    else ! lvnlist .false.
    do ia = iastart , iaend
        rxi = rx ( ia )
        ryi = ry ( ia )
        rzi = rz ( ia )
        do ja = 1 , natm
          if ( ja .gt. ia ) then
            rxij = rxi - rx ( ja )
            ryij = ryi - ry ( ja )
            rzij = rzi - rz ( ja )
            rijsq = rxij  * rxij + ryij * ryij + rzij * rzij
            p1 = itype ( ja )
            p2 = itype ( ia )
            if ( rijsq .lt. rcutsq(p1,p2) ) then
              sr2 = sigsq(p1,p2) / rijsq
              srp = sr2 **  (ptwo(p1,p2))
              srq = sr2 **  (qtwo(p1,p2))
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
    endif


  CALL MPI_ALLREDUCE(u,pot_sum,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
  CALL MPI_ALLREDUCE(vir,vir_sum,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)

  CALL MPI_ALLREDUCE(fx,fx_sum,natm,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
  CALL MPI_ALLREDUCE(fy,fy_sum,natm,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
  CALL MPI_ALLREDUCE(fz,fz_sum,natm,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)

  fx = fx_sum
  fy = fy_sum
  fz = fz_sum
  fx_sum = 0.0d0
  fy_sum = 0.0d0
  fz_sum = 0.0d0

  u_lj = pot_sum
  vir_lj = vir_sum
  pot_sum = 0.0d0
  vir_sum = 0.0d0 

  deallocate( fx_sum, fy_sum, fz_sum )

  return

END SUBROUTINE engforce_bmlj_nopbc

!*********************** SUBROUTINE engforce_new_bmlj *************************
!
! total potential energy forces for each atoms for a bmlj potential with 
! periodic boundaries conditions, with lvnlist=.TRUE.
! However, this version avoid the use of minimum image convention.
! The neighbour list is created in the infinite repeated unit cell system  
!
!******************************************************************************

SUBROUTINE engforce_new_bmlj ( iastart , iaend )

  USE config,           ONLY : natm , box , rx , ry , rz , fx , fy , fz, atype , itype , list , point , ntype
  USE control,          ONLY : lvnlist , myrank
  USE thermodynamic,    ONLY : u_lj , vir_lj
  USE time

  implicit none
  INCLUDE 'mpif.h'

  ! global
  integer, intent(in)  :: iastart , iaend

  ! local
  integer :: ia , ja , it , jt , je , jb , ierr
  integer :: p1 , p2
  double precision :: rxi , ryi , rzi
  double precision :: rxij , ryij , rzij , sr2 , rijsq , srp , srq
  double precision :: wij , fxij , fyij , fzij
  double precision :: ptwo ( ntype , ntype )
  double precision :: qtwo ( ntype , ntype )
  double precision :: one13 , two13 , forcetime1 , forcetime2 , invbox
  double precision :: pot_sum , vir_sum
  double precision :: u , vir
  double precision, dimension(:), allocatable :: fx_sum , fy_sum , fz_sum

#ifdef debug
  print*,'atype',atype
  print*,'itype',itype
#endif

  forcetime1 = MPI_WTIME(ierr) ! timing info

  allocate( fx_sum(natm), fy_sum(natm), fz_sum(natm) )

  fx_sum = 0.0D0
  fy_sum = 0.0D0
  fz_sum = 0.0D0

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

  !if ( lvnlist ) then
    CALL vnlistcheck ( iastart , iaend )
    do ia = iastart , iaend
      rxi = rx ( ia )
      ryi = ry ( ia )
      rzi = rz ( ia )
      jb = point ( ia )
      je = point ( ia + 1 ) - 1
!      do j1 = jb, je
      do ja = 1, natm
!        ja = list(j1)
        if ( ja .ne. ia ) then
          rxij = rxi - rx ( ja )
          ryij = ryi - ry ( ja )
          rzij = rzi - rz ( ja )
!          nxij = nint( rxij * invbox )
!          nyij = nint( ryij * invbox )
!          nzij = nint( rzij * invbox )
!          rxij = rxij - box * nxij
!          ryij = ryij - box * nyij
!          rzij = rzij - box * nzij
          rijsq = rxij * rxij + ryij * ryij + rzij * rzij
          p1 = itype ( ja )
          p2 = itype ( ia )
          if ( rijsq .lt. rcutsq(p1,p2) ) then
            sr2 = sigsq(p1,p2) / rijsq
            srp = sr2 ** (ptwo(p1,p2))
            srq = sr2 ** (qtwo(p1,p2))
            if ( trunc .eq. 1 ) then
              u =  u + epsp ( p1 , p2 ) * ( plj ( p1 , p2 ) * srq - &
                                            qlj ( p1 , p2 ) * srp ) & 
                   - uc ( p1 , p2 )
            endif
            if ( trunc .eq. 2 ) then
              u =  u + epsp ( p1 , p2 ) * ( plj ( p1 , p2 ) * srq   - &
                                            qlj ( p1 , p2 ) * srp ) + & 
                                            uc1 (p1 , p2 ) * rijsq - uc2 (p1 , p2 )
            endif
            wij  = fc ( p1 , p2 ) * ( srq - srp ) * sr2
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
          endif
        endif
      enddo
    enddo
    vir = vir/3.0D0
    !endif

  forcetime2 = MPI_WTIME(ierr) ! timing info
  forcetimetot = forcetimetot+(forcetime2-forcetime1)

  pot_sum = 0.0d0
  vir_sum = 0.0d0
  CALL MPI_ALLREDUCE(u,pot_sum,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
  CALL MPI_ALLREDUCE(vir,vir_sum,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)

  CALL MPI_ALLREDUCE(fx,fx_sum,natm,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
  CALL MPI_ALLREDUCE(fy,fy_sum,natm,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
  CALL MPI_ALLREDUCE(fz,fz_sum,natm,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)

  fx = fx_sum
  fy = fy_sum
  fz = fz_sum
  fx_sum = 0.0d0
  fy_sum = 0.0d0
  fz_sum = 0.0d0

  u_lj = pot_sum
  vir_lj = vir_sum

  deallocate( fx_sum, fy_sum, fz_sum )

  return

END SUBROUTINE engforce_new_bmlj


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
  double precision :: pot_sum, vir_sum
  double precision, dimension(:), allocatable :: fx_sum, fy_sum, fz_sum
  double precision :: u, vir

  forcetime1 = MPI_WTIME(ierr) ! timing info
  
  allocate( fx_sum(natm), fy_sum(natm), fz_sum(natm) )

  fx_sum = 0.0D0
  fy_sum = 0.0D0
  fz_sum = 0.0D0

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
    CALL vnlistcheck ( iastart , iaend )!, list , point )
    do ia = iastart , iaend
      rxi = rx ( ia )
      ryi = ry ( ia )
      rzi = rz ( ia )
      jb = point ( ia )
      je = point ( ia + 1 ) - 1
      do j1 = jb, je
        ja = list(j1)
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
   
    else ! lvnlist .false.
      do ia = iastart , iaend
        rxi = rx ( ia )
        ryi = ry ( ia )
        rzi = rz ( ia )
        do ja =  1, natm
          if ( ja .gt. ia ) then    
              rxij = rxi - rx ( ja )
              ryij = ryi - ry ( ja )
              rzij = rzi - rz ( ja )
              nxij = nint( rxij * invbox )
              nyij = nint( ryij * invbox )
              nzij = nint( rzij * invbox )
              rxij = rxij - box * nxij
              ryij = ryij - box * nyij
              rzij = rzij - box * nzij
              rijsq = rxij  * rxij + ryij * ryij + rzij * rzij
              p1 = itype ( ja )
              p2 = itype ( ia )
              if ( rijsq .lt. rcutsq(p1,p2) ) then
                sr2 = sigsq(p1,p2)/rijsq
                srp = sr2 ** (ptwo(p1,p2))
                srq = sr2 ** (qtwo(p1,p2))
                u =  u + epsp(p1,p2) * ( plj(p1,p2) * srq -qlj(p1,p2)  * srp ) 
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
    endif    

  forcetime2 = MPI_WTIME(ierr) ! timing info
  forcetimetot = forcetimetot+(forcetime2-forcetime1)

  pot_sum = 0.0d0
  vir_sum = 0.0d0 
  CALL MPI_ALLREDUCE(u,pot_sum,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
  CALL MPI_ALLREDUCE(vir,vir_sum,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)

  CALL MPI_ALLREDUCE(fx,fx_sum,natm,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
  CALL MPI_ALLREDUCE(fy,fy_sum,natm,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
  CALL MPI_ALLREDUCE(fz,fz_sum,natm,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)

  fx = fx_sum
  fy = fy_sum
  fz = fz_sum
  fx_sum = 0.0d0
  fy_sum = 0.0d0
  fz_sum = 0.0d0

  u_lj   = pot_sum
  vir_lj = vir_sum

  deallocate( fx_sum, fy_sum, fz_sum )

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
  double precision :: pot_sum, vir_sum
  double precision, dimension(:), allocatable :: fx_sum, fy_sum, fz_sum
  double precision :: u , vir

  forcetime1 = MPI_WTIME(ierr) ! timing info
  
  allocate( fx_sum(natm), fy_sum(natm), fz_sum(natm) )

  fx_sum = 0.0D0
  fy_sum = 0.0D0
  fz_sum = 0.0D0

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

  pot_sum = 0.0d0
  vir_sum = 0.0d0 
  CALL MPI_ALLREDUCE(u,pot_sum,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
  CALL MPI_ALLREDUCE(vir,vir_sum,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)

  CALL MPI_ALLREDUCE(fx,fx_sum,natm,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
  CALL MPI_ALLREDUCE(fy,fy_sum,natm,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
  CALL MPI_ALLREDUCE(fz,fz_sum,natm,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)

  fx = fx_sum
  fy = fy_sum
  fz = fz_sum
  fx_sum = 0.0d0
  fy_sum = 0.0d0
  fz_sum = 0.0d0

  u_lj = pot_sum
  vir_lj = vir_sum

  deallocate( fx_sum, fy_sum, fz_sum )

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
      km_coul%ncell = ncellewald
      nkcut = ( 2 * ncellewald + 1 ) ** 3 
      nkcut = nkcut - 1
      km_coul%nkcut = nkcut
      allocate( km_coul%kptk( nkcut ) , km_coul%kpt(3,nkcut) )
      allocate ( km_coul%strf ( nkcut, ntype) ) 
      CALL kpoint_sum_init ( km_coul )
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
  double precision :: pot_sum, vir_sum
  double precision :: magd_it, ene_it 
  double precision, dimension(:), allocatable :: fx_dir, fy_dir, fz_dir
  double precision, dimension(:), allocatable :: fx_sum, fy_sum, fz_sum
  double precision, dimension(:), allocatable :: u_dir_magd


  ttt1 = MPI_WTIME(ierr)

  allocate( fx_sum(natm), fy_sum(natm), fz_sum(natm) )
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
 
  fx_sum  = 0.0D0
  fy_sum  = 0.0D0
  fz_sum  = 0.0D0

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

  pot_sum = 0.0d0
  vir_sum = 0.0d0 
  
  CALL MPI_ALLREDUCE(u_dir,pot_sum,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
  CALL MPI_ALLREDUCE(vir_dir,vir_sum,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)

  CALL MPI_ALLREDUCE(fx_dir,fx_sum,natm,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
  CALL MPI_ALLREDUCE(fy_dir,fy_sum,natm,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
  CALL MPI_ALLREDUCE(fz_dir,fz_sum,natm,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)

  fx = fx + fx_sum
  fy = fy + fy_sum
  fz = fz + fz_sum

  fx_sum = 0.0d0
  fy_sum = 0.0d0
  fz_sum = 0.0d0

  u_coul = pot_sum
  vir_coul = vir_sum

  ! ================================================
  !  WARNING  Empirical coding !!!!!!! WARNING
  !  Where is the sign error !!! (tested on NaCl and compared to ewald and GULP)
  !  TODO:  test on another system
  ! ================================================
  u_coul   = - u_coul
  vir_coul = - vir_coul


  deallocate( fx_sum, fy_sum, fz_sum )
  deallocate( fx_dir, fy_dir, fz_dir )
  deallocate( u_dir_magd )

  ttt3 = MPI_WTIME(ierr)

  return

END SUBROUTINE engforce_coulomb_DS


END MODULE field 
! ===== fmV =====
