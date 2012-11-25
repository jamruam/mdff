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

MODULE radial_distrib 

  integer :: PANGR             ! (internal) number of bins in g(r) distribution
  integer :: nskip
  integer :: nc
 
  double precision :: resg     ! resolution in g(r) distribution 

  ! gr(4,0:PANGR)
  ! gr(1,:) : AA
  ! gr(2,:) : BB
  ! gr(3,:) : AB
  ! gr(4,:) : all
  integer, dimension(:,:), allocatable :: gr 

CONTAINS

!*********************** SUBROUTINE grcalc_init *******************************
!
!
!******************************************************************************

SUBROUTINE gr_init

  USE control,  ONLY :  calc
  USE prop,     ONLY :  lgr
  USE io_file,  ONLY :  stdin, stdout, kunit_OUTFF, ionode

  implicit none

  ! local
  integer :: ioerr
  character * 132 :: filename

  namelist /grtag/   nc    , &
                     nskip , &
                     resg  

  if ( .not. lgr .and. calc .ne. 'gr' ) return

  CALL gr_default_tag
  
  ! =================
  !  read grtag tags
  ! =================
  CALL getarg (1,filename)
  OPEN ( stdin , file = filename)
  READ ( stdin , grtag , iostat=ioerr )
  if ( ioerr .lt. 0 )  then
   if ( ionode ) WRITE ( stdout, '(a)') 'ERROR reading input_file : vacftag section is absent'
   STOP
  elseif ( ioerr .gt. 0 )  then
   if ( ionode ) WRITE ( stdout, '(a)') 'ERROR reading input_file : vacftag wrong tag'
   STOP
  endif

  CLOSE ( stdin )

  CALL gr_check_tag

  CALL gr_print_info(stdout)
  CALL gr_print_info(kunit_OUTFF)

  return 
 
END SUBROUTINE gr_init

!*********************** SUBROUTINE gr_alloc **********************************
!
!
!******************************************************************************

SUBROUTINE gr_alloc

  USE control,  ONLY :  calc
  USE prop,     ONLY :  lgr

  implicit none

  if ( .not. lgr .and. calc .ne. 'gr' ) return

  allocate(gr(4,0:PANGR))
  gr = 0      

  return 
 
END SUBROUTINE gr_alloc


!*********************** SUBROUTINE gr_dealloc **********************************
!
!
!******************************************************************************

SUBROUTINE gr_dealloc

  USE control,  ONLY :  calc
  USE prop,     ONLY :  lgr

  implicit none

  if ( .not. lgr .and. calc .ne. 'gr' ) return
  
  deallocate( gr )

  return 
 
END SUBROUTINE gr_dealloc


!*********************** SUBROUTINE gr_default_tag ***************************
!
! set default values to gr tag
!
!******************************************************************************

SUBROUTINE gr_default_tag

  implicit none

  ! ===============
  !  default value
  ! ===============
  resg = 0.1d0
  nskip = 0
  nc = 0

  return

END SUBROUTINE gr_default_tag


!*********************** SUBROUTINE gr_check_tag ***************************
!
! check gr tag values
!
!******************************************************************************

SUBROUTINE gr_check_tag

  USE config,   ONLY :  box

  implicit none

  ! local
  double precision :: cutgr
  
  cutgr=box*0.5D0
  PANGR=int(cutgr/resg)+1

  return

END SUBROUTINE gr_check_tag

!*********************** SUBROUTINE gr_print_info ***************************
!
!
!******************************************************************************

SUBROUTINE gr_print_info(kunit)

  USe control,  ONLY :  calc
  USE io_file,  ONLY :  ionode 

  implicit none
 
  ! local
  integer :: kunit

   if ( ionode ) then
                  WRITE ( kunit ,'(a,f10.5)')           'resolution of g(r) function resg     = ',resg
                  WRITE ( kunit ,'(a)')                 'save radial_distribution in file     :    GRTFF' 
      if ( calc .eq. 'gr' )     then 
                  WRITE ( kunit ,'(a)')                 'read configuration from file          :   TRAJFF'
                  WRITE ( kunit ,'(a)')                 ''
                  WRITE ( kunit ,'(a,i5)')              'number of config. in TRAJFF         =     ',nc        
                  WRITE ( kunit ,'(a,i5)')              'number of config. to be skipped     =     ',nskip
                  WRITE ( kunit ,'(a)')                 ''
      endif
   endif 
  return

END SUBROUTINE gr_print_info

!*********************** SUBROUTINE grcalc ************************************
!
!
!******************************************************************************
SUBROUTINE grcalc

  USE config,           ONLY :  system , natm , ntype , rx , ry , rz , atype , &
                                rho , config_alloc , box , omega , atypei , itype, natmi
  USE control,          ONLY :  myrank , numprocs
  USE io_file,          ONLY :  ionode , stdout , kunit_TRAJFF , kunit_GRTFF , kunit_OUTFF
  USE time

  implicit none
  INCLUDE 'mpif.h'

  ! local 
  integer :: ia , ic , it , ngr , na
  integer :: iastart , iaend 
  ! trash 
  double precision :: aaaa
  integer :: iiii
  character * 60 :: cccc

  OPEN (UNIT = kunit_TRAJFF ,FILE = 'TRAJFF') 
    
  READ ( kunit_TRAJFF , * ) natm
  READ ( kunit_TRAJFF , * ) system
  READ ( kunit_TRAJFF , * ) box , ntype
  READ ( kunit_TRAJFF ,* ) ( atypei ( it ) , it = 1 , ntype )
  IF ( ionode ) WRITE ( kunit_OUTFF ,'(A,20A3)' ) 'found type information on TRAJFF : ', atypei ( 1:ntype )
  IF ( ionode ) WRITE ( stdout      ,'(A,20A3)' ) 'found type information on TRAJFF : ', atypei ( 1:ntype )
  READ( kunit_TRAJFF ,*)   ( natmi ( it ) , it = 1 , ntype )
  omega = box * box * box
  rho = dble ( natm )  / omega 
  CALL gr_init

  CALL print_general_info( stdout )
  CALL print_general_info( kunit_OUTFF )

  ! ===================================
  !  here we know natm, then alloc 
  !  and decomposition can be applied 
  ! ================================== 
  CALL config_alloc 
  CALL do_split ( natm , myrank , numprocs , iastart , iaend )
  CALL gr_alloc
#ifdef debug
  print*,'debug',iastart , iaend
#endif
  ! ==========================================   
  !  skip the first nskip configurations 
  ! ==========================================
  if (nskip.gt.0) then
    do ic = 1,nskip
      if ( ic .ne. 1 ) READ ( kunit_TRAJFF , * ) iiii   
      if ( ic .ne. 1 ) READ ( kunit_TRAJFF , * ) cccc
      if ( ic .ne. 1 ) READ ( kunit_TRAJFF , * ) aaaa   ,iiii
      if ( ic .ne. 1 ) READ ( kunit_TRAJFF , * ) ( cccc , it = 1 , ntype )
      if ( ic .ne. 1 ) READ ( kunit_TRAJFF , * ) ( iiii , it = 1 , ntype )
      do ia = 1 , natm 
        READ ( kunit_TRAJFF , * ) atype ( ia ) , rx ( ia ) , ry ( ia ) , rz ( ia ) , aaaa,aaaa,aaaa,aaaa,aaaa,aaaa
      enddo      
    enddo
  endif

  ngr = 0
  do ic = nskip + 1, nc
    if ( ionode ) print*,ic
    na = 0
    ! ===================================
    !  read config from trajectory file
    ! ===================================
    if ( ic .ne. (nskip + 1) .or. nskip.ne.0) READ ( kunit_TRAJFF , * ) iiii
    if ( ic .ne. (nskip + 1) .or. nskip.ne.0) READ ( kunit_TRAJFF , * ) cccc
    if ( ic .ne. (nskip + 1) .or. nskip.ne.0) READ ( kunit_TRAJFF , * ) aaaa , iiii
    if ( ic .ne. (nskip + 1) .or. nskip.ne.0 ) READ ( kunit_TRAJFF , * ) ( cccc , it = 1 , ntype )
    if ( ic .ne. (nskip + 1) .or. nskip.ne.0 ) READ ( kunit_TRAJFF , * ) ( iiii , it = 1 , ntype )
    do ia = 1 , natm
      READ ( kunit_TRAJFF , * ) atype ( ia ) , rx ( ia) , ry ( ia ) , rz ( ia ) , aaaa,aaaa,aaaa,aaaa,aaaa,aaaa
      if ( atype ( ia ) .eq. 'A') na = na + 1
      if ( atype ( ia ) .eq. 'A') itype(ia)=1
      if ( atype ( ia ) .eq. 'B') itype(ia)=2
    enddo
    if ( ntype .eq.  1 ) then
      natmi(0)=natm
      natmi(1)=na
      atypei(0)='ALL'
      atypei(1)='A'
    endif
    if ( ntype .eq.  2 ) then
      natmi(0)=natm
      natmi(1)=na
      natmi(2)=natm-na
      atypei(0)='ALL'
      atypei(1)='A'
      atypei(2)='B'
    endif

      ngr=ngr+1 
      ! ==========================
      !  calc radial_distribution 
      ! ==========================  
      call gr_main ( iastart , iaend , ngr )

  enddo !nc 

  CLOSE( kunit_TRAJFF )

  return

END SUBROUTINE grcalc

!*********************** SUBROUTINE gr_main ***********************************
!
! based on Frenkel and Smit
! the ntype differenciation could be better
!
!******************************************************************************

SUBROUTINE gr_main ( iastart , iaend , ngr )

  USE config,           ONLY :  natm , natmi , rx , ry , rz , atype , box , omega , ntype
  USE prop,             ONLY :  nprop, nprop_start
  USE constants,        ONLY :  pi 
  USE io_file,          ONLY :  ionode , stdout , kunit_GRTFF , kunit_NRTFF
  USE time

  implicit none
  INCLUDE "mpif.h"

  ! global
  integer, intent(in) :: iastart , iaend
  integer, intent(in) :: ngr  

  ! local
  integer :: i , ia , ja , k , ierr , ii 
  integer :: igr , na , nb 
  integer :: nxij , nyij , nzij
  double precision :: cut2 , rr , rijsq , vol
  double precision :: rxi , ryi , rzi
  double precision :: rxij , ryij , rzij
  double precision :: grr1 , grr2 , grr3 , grrtot
  double precision :: grtime1 , grtime2      

  grtime1 = MPI_WTIME(ierr)

  OPEN ( kunit_GRTFF , file = 'GRTFF' )
  OPEN ( kunit_NRTFF , file = 'NRTFF' )

  na = natmi(1)
  nb = natm - na

  do ia = 1 , natm
  if ( ia .gt. iaend ) then
    gr=0
    !write(stdout,*) 'clean part of the array'
  endif
  enddo

  cut2 = box * box * 0.25D0    
  do ia = iastart , iaend
    rxi = rx ( ia )
    ryi = ry ( ia )
    rzi = rz ( ia )
    do ja = 1, natm
      if ( ja .ne. ia ) then  
      rxij = rxi - rx ( ja )
      ryij = ryi - ry ( ja )
      rzij = rzi - rz ( ja )
      nxij = nint( rxij / box )
      nyij = nint( ryij / box )
      nzij = nint( rzij / box )
      rxij = rxij - box *nxij
      ryij = ryij - box *nyij
      rzij = rzij - box *nzij
      rijsq = rxij * rxij + ryij * ryij + rzij * rzij
      if ( rijsq.lt.cut2 ) then
        rr = dsqrt(rijsq)
        igr = int(rr/resg)
        gr(4,igr)  = gr(4,igr)+1  ! all average
        if ( atype ( ia ) .eq. 'A' .and. atype ( ja ) .eq. 'A' ) then
          gr ( 1 , igr )  = gr ( 1 , igr ) + 1   ! AA average
        endif
        if ( atype ( ia ) .eq. 'B' .and. atype ( ja ) .eq. 'B' ) then
          gr ( 2 , igr )  = gr ( 2 , igr ) + 1   ! BB average
        endif
        if ( atype ( ia ) .eq.'A' .and. atype ( ja ) .eq. 'B' ) then!.or.(atype(i).eq.'B'.and.atype(j).eq.'A' )) then 
          gr ( 3 , igr )  = gr ( 3 , igr ) + 1  ! AB or BA average
        endif
      endif
     endif
    enddo
  enddo

  do ii=1,4
    CALL MPI_ALL_REDUCE_INTEGER ( gr(ii,:), PANGR+1 )
  enddo

  do i = 1,PANGR-1
    rr = resg*(dble(i)+0.5D0)
    k = i+1
    k = k*k*k
    k = k-(i*i*i)
    vol = dble(k)*resg*resg*resg ! r^3
    vol = (4.0D0/3.0D0)*pi*vol/omega !4/3pir^3
    grr1   = dble(gr(1,i))/(ngr*vol*na*na)
    if (ntype .ne. 1 ) then
      grr2   = dble(gr(2,i))/(ngr*vol*nb*nb)
      grr3   = dble(gr(3,i))/(ngr*vol*na*nb)
      grrtot = dble(gr(4,i))/(ngr*vol*natm*natm)
      if ( ionode ) then  
        WRITE (kunit_GRTFF,'(5f15.10)') rr , grr1 , grr2 , grr3 , grrtot 
        WRITE (kunit_NRTFF,'(5f20.10)') &
        rr , dble(gr(1,i))*4.0d0*pi*rr*rr/dble(ngr*na) , &
             dble(gr(2,i))*4.0d0*pi*rr*rr/dble(ngr*nb) , &
             dble(gr(3,i))*4.0d0*pi*rr*rr/dble(ngr*nb) , &
             dble(gr(4,i))*4.0d0*pi*rr*rr/dble(ngr*natm) 
      endif
    else
      if ( ionode ) then  
        WRITE (kunit_GRTFF,'(2f15.10,i10)') rr , grr1
        WRITE (kunit_NRTFF,'(f15.10,4i8)')  rr , gr(1,i) 
      endif
     
    endif
  enddo   

#ifdef debug
  do i=1, PANGR-1
    WRITE (kunit_GRTFF,'(a,5i6)') 'debug ',i,gr(:,i)
  enddo
#endif 
 
  CLOSE ( kunit_NRTFF )
  CLOSE ( kunit_GRTFF )

  grtime2 = MPI_WTIME(ierr)
  grtimetot = grtimetot + ( grtime2 - grtime1 )

  return
 
END SUBROUTINE gr_main

!*********************** SUBROUTINE static_struc_fac **************************
!
!
!******************************************************************************

SUBROUTINE static_struc_fac ( ngr )

  USE io_file,  ONLY :  ionode , kunit_STRFACFF
  USE config,   ONLY :  natm , rho , natmi , omega
  USE constants,ONLY :  pi        

  implicit none
  ! global
  integer :: ngr 

  ! local
  integer :: i , k 
  double precision :: q , grr1 , vol , rr 
  double precision , dimension (:) , allocatable  :: stat_str
  double precision   ,dimension (:), allocatable :: in 
  double complex     ,dimension (:), allocatable :: out

  do i = 1,PANGR-1
    rr = resg*(dble(i)+0.5D0)
    k = i+1
    k = k**3
    k = k-(i*i*i)
    vol = dble(k)*resg*resg*resg
    vol = (4.0D0/3.0D0)*pi*vol/omega
    grr1   = dble(gr(1,i))/(ngr*vol*natmi(1)*natm)
  !  gr(1,i) = grr1
  enddo

  allocate ( stat_str (PANGR) )
  allocate ( in(PANGR) , out(PANGR/2 +1) )
  ! ========
  !   FFT
  ! ========
  in  = gr(1,:) 
  CALL fft_1D_real ( in , out , pangr )
  stat_str = dble ( out ) 
  !stat_str = 1.0d0 + rho * stat_str

  do i = 1,PANGR-1
    if ( ionode ) then
      WRITE ( kunit_STRFACFF ,'(2f15.10,i10)') q , stat_str(i) 
    endif
  enddo

  deallocate ( stat_str )

  
  return

END SUBROUTINE static_struc_fac


END MODULE radial_distrib 
! ===== fmV =====
