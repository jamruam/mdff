! MDFF parallel Molecular Dynamics ... For Fun
! Copyright (C) 2014  F. Vasconcelos
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
! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301,
! USA.
! ===== fmV =====

! ======= Hardware =======
#include "symbol.h"
! ======= Hardware =======


! *********************** MODULE stochio ***************************************
!> \brief 
!! module related to stochiometry calculation ( calc = 'stochio' )
! ******************************************************************************
MODULE stochio

  USE constants,        ONLY :  dp, g_to_am
  USE oxyde
  USE io

  implicit none

  ! =====================================================
  !   type of data : pct (i.e %) or num (i.e 0 < pox < 1 )              
  ! =====================================================
  character(len=3) , SAVE :: def          
  character(len=60), SAVE :: def_allowed(2)
  data def_allowed / 'pct' , 'num' /

  integer             :: i , j , ie
  integer             :: ioerr
  integer             :: minel 
  
  integer             :: target_nions , ntot , save_target
  integer             :: nel    ( nelem )              ! number of ions per element has to be calculated
!  real(kind=dp)       :: pox ( noxyde )                ! proportion en oxyde (input)
  real(kind=dp)       :: rat ( noxyde )                ! proportion en oxyde (output)
  real(kind=dp)       :: sumox , charg , density , a_o_b , a_o_c 
  real(kind=dp)       :: volume , totmass , acell , bcell , ccell , volume2
  character(len=135)  :: filename
  logical             :: allowed
  logical             :: lexact , lcubic
  integer             :: sum_sto, numbands


CONTAINS
  
! *********************** MODULE stochio_calc **********************************
!> \brief 
! ******************************************************************************
SUBROUTINE stochio_calc

  implicit none

  integer :: iox

  ! ========================
  !  define 0 < pox < 1
  ! ========================
  if ( def .eq. 'num' ) then 
    oxydes%relcon = oxydes%relcon * 100._dp
  endif
  ! ============================
  ! check that the sum is == 1
  ! if not STOP 
  ! ============================
  sumox=0._dp
  do iox=1,noxyde
    sumox=sumox+oxydes(iox)%relcon
  enddo
  if ( REAL(sumox,kind=dp) .ne. 100._dp ) then
    WRITE ( stdout , '(a,f24.16)' ) 'ERROR sum of oxydes not 100 % : ', REAL(sumox,kind=dp)
    STOP
  endif
  oxydes%relcon = oxydes%relcon / 100._dp

  ! ================================ 
  !  check if the number of ions 
  !  requested is large enough
  ! ================================ 
  minel=0
  do iox=1,noxyde
    if ( oxydes(iox)%relcon .eq. 0.0_dp ) then
      do j=1,2
        minel = minel + oxydes(iox)%nel_ox(j)
      enddo
    endif
  enddo
  if ( minel .gt. target_nions ) then 
    WRITE ( stdout , '(a,2f24.16)' ) 'ERROR target_nions too small : ',target_nions, minel
    STOP
  endif

  ! ==========================================================
  ! get target_nions
  ! starting from the wanted number of ions (save_target)
  ! we calculate the total number of ions :
  !         ___
  !ntot =   \    target_nions x oxyde(%) x nions_per_oxyde
  !         /__
  !       oxyde,eleme
  ! 
  ! ==========================================================
  save_target = target_nions
  target_nions = target_nions + 1
  ntot = 0
  do while ( ntot .ne. save_target .and. target_nions .gt. 0 ) 
    target_nions = target_nions - 1
    ntot = 0
    do iox = 1 , noxyde
      if ( oxydes(iox)%relcon .ne.0._dp ) then
        do j = 1 , 2 
          ntot = ntot + int(target_nions*oxydes(iox)%relcon)*oxydes(iox)%nel_ox(j)
        enddo
      endif
    enddo
  enddo
  if ( target_nions .le. 0 ) then
    WRITE ( stdout ,'(a,2i)' ) 'could not found configuration for the requested number ions :',target_nions,save_target
    STOP
  endif
  do iox = 1 , noxyde
    if ( oxydes(iox)%relcon .ne.0._dp ) then
      if ( int(target_nions*oxydes(iox)%relcon) .eq. 0 ) STOP
    endif
  enddo
  WRITE ( stdout ,'(a,i)' ) 'target found :',target_nions

  ! ===========================
  ! count the number of ions for  
  ! elements from target_nions
  ! ===========================
  nel = 0
  do iox = 1, noxyde
    if ( oxydes(iox)%relcon.ne.0._dp ) then
      do j=1,2
        do ie = 1 , nelem
          if ( oxydes(iox)%ele_ox( j ) .eq. tabper(ie)%elename ) then
            nel(ie) = nel(ie) + int(target_nions*oxydes(iox)%relcon)*oxydes(iox)%nel_ox( j )
          endif
        enddo
      enddo
    endif
  enddo

  ! ==============================
  ! check that the charg is null
  ! ==============================
  charg = 0._dp
  do ie=1,nelem
    if ( nel(ie) .ne. 0 ) then
      charg = charg + nel(ie) * tabper(ie)%numoxyd
    endif
  enddo
  if ( charg .ne. 0._dp ) then
    WRITE ( stdout , '(a,2f24.16)' ) 'ERROR total charge is not zero : ',charg
    STOP
  endif

  ! ========================================
  !  WRITE info to standard output
  ! ========================================
  WRITE( stdout , '(a)' ) '#elem     nsites' 
  ntot = 0
  do ie=1,nelem
    if ( nel(ie) .ne. 0 ) then
      ntot = ntot + nel(ie)
      WRITE( stdout, '(A5,I10)' ) tabper(ie)%elename,nel(ie) 
    endif
  enddo
  WRITE( stdout, '(A5,I10)' )   'total',ntot


  ! ===============================================================
  !  calculate the maybe approximated oxyde proportion : rat(iox)
  ! ===============================================================
  ! first sum total number of oxyde molecules (per mole)
  sum_sto = 0
  do iox = 1,noxyde
    do ie = 1 , nelem
      if ( tabper(ie)%elename .eq. oxydes(iox)%ele_ox( 1 ) ) then
        sum_sto =sum_sto + (nel(ie)/oxydes(iox)%nel_ox( 1 ) )
      endif
    enddo
  enddo
  ! then calculates the ration ratio "rat( iox )"
  ! and check if it is the exact wanted ratio comparing to all pox( iox )
  lexact = .true.
  do iox = 1 , noxyde
    if ( oxydes(iox)%relcon .ne. 0._dp ) then
      do ie = 1 , nelem
        if ( tabper(ie)%elename .eq. oxydes(iox)%ele_ox( 1 ) ) then
          rat(iox)=REAL(nel(ie)/oxydes(iox)%nel_ox( 1 ),kind=dp)/REAL(sum_sto,kind=dp)
          if ( rat( iox ) .ne. oxydes(iox)%relcon ) then
            lexact =.false.
          endif
        endif
      enddo
    endif
  enddo
  ! finally prints out
  WRITE(stdout , * ) ''
  if ( lexact ) then
    WRITE(stdout , '(a)' ) '#stoichiometry is EXACT'
  else
    WRITE(stdout , '(a)' ) '#stoichiometry is APPROXIMATED'
  endif
  WRITE(stdout , '(a)' ) '#oxyde    current(%)      wanted(%)           nb of oxyde (mole)'
  do iox = 1 , noxyde
    if ( oxydes(iox)%relcon .ne. 0.0_dp ) then
      do ie = 1 , nelem
        if ( tabper(ie)%elename .eq. oxydes(iox)%ele_ox(1) ) then
          WRITE(stdout , '(A4,2F14.3,I20)' ) oxydes(iox)%nameox,rat( iox )*100._dp, oxydes(iox)%relcon *100._dp, nel(ie)/oxydes(iox)%nel_ox(1)
        endif
      enddo
    endif
  enddo
  WRITE(stdout , * ) ''

  ! ==========================================================
  !  calculate the cell parameters cubic or not orthorhombic
  !  default : cubic 
  ! if orthorhombic the ratio a/b and a/c should be given
  ! ==========================================================
  if ( density .eq. 0._dp ) then
    WRITE(stdout , '(a)' ) 'No density in input file'
  else
    totmass= 0.0_dp
    numbands = 0
    WRITE( stdout , '(a)' ) '           Z     valence    mass    element       nb ele'
    do ie = 1 , nelem
      if ( nel(ie) .ne. 0 ) then
        totmass = totmass + tabper(ie)%massele * nel(ie) 
        numbands = numbands + tabper(ie)%valence * nel(ie)
        WRITE( stdout , '(i,i,f8.3,10x,a,i)' ) ie,tabper(ie)%valence,tabper(ie)%massele,tabper(ie)%elename,nel(ie)
      endif
    enddo
    blankline(stdout)
    volume = totmass /  density * g_to_am   
    acell =( volume /  REAL(a_o_b,kind=dp) / REAL(a_o_c,kind=dp) )**(1._dp/3._dp) 
    WRITE(stdout , '(a,i12)' )       'NBANDS        = ',numbands/2
    WRITE(stdout , '(a,f12.4,a)' )   'density       = ',density,' g/cm^3 '
    WRITE(stdout , '(a,f12.4,a)' )   'total mass    = ',totmass,' [a.m]  '
    WRITE(stdout , '(a,f12.4,a)' )   'volume        = ',volume ,' A '
    if ( lcubic ) then 
      WRITE(stdout , '(a,f12.4,a)' ) 'cell param. a = ', acell, ' A ' 
    else
      bcell = acell*a_o_b
      ccell = acell*a_o_c
      WRITE(stdout , '(a,f16.8,a)' ) 'cell param. a = ', acell, ' A ' 
      WRITE(stdout , '(a,f16.8,a)' ) 'cell param. b = ', bcell, ' A ' 
      WRITE(stdout , '(a,f16.8,a)' ) 'cell param. c = ', ccell, ' A ' 
      volume2 = acell*bcell*ccell
      if ( abs(volume2 - volume ) >= 1e-8_dp ) then
        WRITE(stdout , '(a,2f42.24)' ) 'ERROR in volume bug found',volume2 ,volume
        STOP
      endif
    endif
  endif

  CALL tabper_oxydes_dalloc 

  return

END SUBROUTINE stochio_calc

! *********************** MODULE stochio_init **********************************
!> \brief 
! ******************************************************************************
SUBROUTINE stochio_init

  USE control,          ONLY :  calc
  USE io,               ONLY :  ionode ,stdin, stdout

  implicit none

  integer :: iox

  namelist /stochiotag/    def           , &
                           density       , &
                           lcubic        , &
                           a_o_b         , &
                           a_o_c         , &
                           target_nions  , &
                           sio2          , &
                           geo2          , &
                           na2o          , &
                           b2o3          , &
                           la2o3         , &
                           cao           , &
                           al2o3         , &
                           p2o5

  if ( calc .ne. 'stochio' ) return



  CALL stochio_default_tag

  CALL getarg (1,filename)
  OPEN ( stdin , file = filename)
  READ ( stdin , stochiotag, iostat=ioerr)
  if ( ioerr .lt. 0 )  then
    WRITE ( stdout, '(a)') 'ERROR reading input_file : input section is absent'
    STOP
  elseif ( ioerr .gt. 0 )  then
    WRITE ( stdout, '(a,i6)') 'ERROR reading input_file : input wrong tag',ioerr
    WRITE ( stdout, '(a,i6)' ) 'Number of available oxydes :' ,noxyde
    WRITE ( stdout, '(<noxyde>a6)' ) ( oxydes(iox)%nameox , iox = 1 , noxyde )
    WRITE ( stdout, '(a)' ) 'add some oxyde in oxyde.f90'
    STOP
  endif
  CLOSE ( stdin )

  CALL gen_tab_period
  CALL gen_oxydes

  ! =======================
  !  check stochiotag info
  ! =======================
  CALL stochio_check_tag 

  ! =======================
  !  print stochiotag info
  ! =======================
  CALL stochio_print_info(stdout)

  return

END SUBROUTINE stochio_init


! *********************** MODULE stochio_default_tag ***************************
!> \brief 
! ******************************************************************************
SUBROUTINE stochio_default_tag

  implicit none

  ! ======================
  !  set default values
  ! ======================
  lcubic  = .true.
  a_o_b   = 1.0_dp
  a_o_c   = 1.0_dp
  sio2    = 0._dp
  na2o    = 0._dp
  b2o3    = 0._dp
  la2o3   = 0._dp
  cao     = 0._dp
  p2o5    = 0._dp
  al2o3   = 0._dp
  geo2    = 0._dp
  def     = 'pct'
  density = 0.0d0

  return

END SUBROUTINE stochio_default_tag

! *********************** MODULE stochio_check_tag *****************************
!> \brief 
! ******************************************************************************
SUBROUTINE stochio_check_tag
  

  implicit none

 integer :: ntype, iox


  if ( lcubic ) then
    if ( ( a_o_b .eq. 0._dp) .or. ( a_o_c .eq. 0._dp) ) then
      WRITE ( stdout, '(a,2f8.3)') 'with lcubic cell parameters ration should be equal to 1.0 ',a_o_b,a_o_c
      a_o_b = 1.0_dp
      a_o_c = 1.0_dp
    endif
  endif

  ! =================
  !  check def param
  ! =================
  allowed = .false.
  do i = 1 , size( def_allowed )
   if ( trim(def) .eq. def_allowed(i))  allowed = .true.
  enddo
  if ( .not. allowed ) then
      WRITE ( stdout , '(a)' ) 'ERROR controltag: calc should be ', def_allowed
      STOP
  endif


  do iox = 1 , noxyde
    if ( oxydes(iox)%relcon .ne. 0.0d0) ntype=ntype+1
  enddo
  WRITE ( stdout , '(a,i)' ) 'ntypes ', ntype

  return

END SUBROUTINE stochio_check_tag

SUBROUTINE stochio_print_info ( kunit )

  implicit none

  !local
  integer :: kunit, iox

  if ( ionode ) then
    separator(kunit)
    blankline(kunit)
    WRITE( kunit , '(a)' ) 'STOCHIO MODULE ... WELCOME'
    WRITE ( stdout , '(a)'            ) 'Composition :'
    WRITE ( stdout , '(<noxyde>a8)'   ) (oxydes(iox)%nameox,iox=1,noxyde)
    WRITE ( stdout , '(<noxyde>f8.2)' ) (oxydes(iox)%relcon,iox=1,noxyde)
  endif
  

  return

END SUBROUTINE stochio_print_info

END MODULE stochio
