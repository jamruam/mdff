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
#define debug2
! ======= Hardware =======

! *********************** SUBROUTINE grcalc_init *******************************
!
!> \brief
!! Module related to radial function distribution calculation and/or static
!! factor structure
!
! ******************************************************************************
MODULE radial_distrib 

  USE constants,                ONLY :  dp

  implicit none

  integer :: PANGR            !< (internal) number of bins in g(r) distribution
  integer :: nskip            !< number of configurations skipped 
  integer :: nconf            !< number of configurations used in g(r) calculation
  real(kind=dp) :: cutgr      !< radial cut-off 
  real(kind=dp) :: resg       !< resolution in g(r) distribution 
  !> g(r) function ( bin x ntype x ntype )
  integer, dimension(:,:,:), allocatable :: gr  

CONTAINS

! *********************** SUBROUTINE grcalc_init *******************************
!
!> \brief
!! initialize radial distribution calculation parameters
!
! ******************************************************************************
SUBROUTINE gr_init

  USE config,                   ONLY :  simu_cell
  USE control,                  ONLY :  calc
  USE io_file,                  ONLY :  stdin , stdout , ionode

  implicit none

  ! local
  integer            :: ioerr 
!28/05/13  integer            :: npangr, i
  character(len=132) :: filename

  namelist /grtag/   nconf , &
                     nskip , &
                     cutgr , &
                     resg  

  if ( calc .ne. 'gr' ) return

  CALL gr_default_tag
  
  ! =================
  !  read grtag tags
  ! =================
  CALL getarg (1,filename)
  OPEN ( stdin , file = filename)
  READ ( stdin , grtag , iostat=ioerr )
  if ( ioerr .lt. 0 )  then
   io_node WRITE ( stdout, '(a)') 'ERROR reading input_file : grtag section is absent'
   STOP
  elseif ( ioerr .gt. 0 )  then
   io_node WRITE ( stdout, '(a,i8)') 'ERROR reading input_file : grtag wrong tag'
   STOP
  endif
  CLOSE ( stdin )

  ! ==========================================
  ! define a new resolution to be 2^N points
  ! ==========================================
 PANGR=int(cutgr/resg)+1
!28/05/13 ! i = 1
!28/05/13 ! do while ( 2**i .lt. PANGR )
!28/05/13 !    i = i + 1
!28/05/13 ! enddo
!28/05/13 ! npangr = i
!28/05/13 ! PANGR = 2** npangr
!28/05/13 ! resg = cutgr / DBLE ( PANGR  )

  CALL gr_print_info(stdout)

  return 
 
END SUBROUTINE gr_init

! *********************** SUBROUTINE gr_alloc **********************************
!
!> \brief
!! allocate g(r) function 
!
! ******************************************************************************
SUBROUTINE gr_alloc

  USE control,                  ONLY :  calc
  USE config,                   ONLY :  ntype

  implicit none

  if ( calc .ne. 'gr' ) return

  allocate(gr(0:PANGR,0:ntype,0:ntype))
  gr = 0      

  return 
 
END SUBROUTINE gr_alloc


! *********************** SUBROUTINE gr_dealloc ********************************
!
!> \brief
!! deallocate g(r) function 
!
! ******************************************************************************
SUBROUTINE gr_dealloc

  USE control,                  ONLY :  calc

  implicit none

  if ( calc .ne. 'gr' ) return
  
  deallocate( gr )

  return 
 
END SUBROUTINE gr_dealloc


! *********************** SUBROUTINE gr_default_tag ****************************
!
!> \brief
!! set default values to gr tag
!
! ******************************************************************************
SUBROUTINE gr_default_tag

  USE config,           ONLY : simu_cell

  implicit none

  ! ===============
  !  default value
  ! ===============
  resg = 0.1_dp
  nskip = 0
  nconf = 0
  cutgr=0.5_dp * MIN(simu_cell%WA,simu_cell%WB,simu_cell%WC)-0.01_dp

  return

END SUBROUTINE gr_default_tag


! *********************** SUBROUTINE gr_print_info *****************************
!
!> \brief
!! print infog on g(r) calculation
!
! ******************************************************************************
SUBROUTINE gr_print_info(kunit)

  USe control,                  ONLY :  calc
  USE io_file,                  ONLY :  ionode 

  implicit none
 
  ! local
  integer :: kunit

   if ( ionode ) then
                  WRITE ( kunit ,'(a,f10.5,a)')         'resolution of g(r) function resg     = ',resg,' new value to have 2^N points in g(r)'
                  WRITE ( kunit ,'(a,i5)')              'number of points in g(r)             = ',PANGR
                  WRITE ( kunit ,'(a)')                 'save radial_distribution in file     :   GRTFF' 
      if ( calc .eq. 'gr' )     then 
                  WRITE ( kunit ,'(a)')                 'read configuration from file         :   TRAJFF'
                  blankline(kunit)
                  WRITE ( kunit ,'(a,i5)')              'number of config. in TRAJFF          = ',nconf        
                  WRITE ( kunit ,'(a,i5)')              'number of config. to be skipped      = ',nskip
                  blankline(kunit)
      endif
   endif 
  return

END SUBROUTINE gr_print_info

! *********************** SUBROUTINE grcalc ************************************
!
!> \brief
!! main driver of radial distribution function calculation
!! this subroutine read the trajectory, allocate, call the  
!
! ******************************************************************************
SUBROUTINE grcalc

  USE config,                   ONLY :  system , natm , ntype , rx , ry , rz , atype , &
                                        rho , config_alloc , simu_cell , atypei , itype, natmi
  USE control,                  ONLY :  myrank , numprocs
  USE io_file,                  ONLY :  ionode , stdout , stderr , kunit_TRAJFF , kunit_GRTFF , kunit_NRTFF
  USE constants,                ONLY :  pi 
  USE cell,                     ONLY :  lattice , dirkar
  USE time,                     ONLY :  grtimetot_comm

  implicit none
  INCLUDE 'mpif.h'

  ! local 
  integer                                              :: ia , ic , it , ngr , i  
  integer                                              :: pairs , it1 , it2 , mp , ierr 
  integer                                              :: iastart , iaend 
  real(kind=dp),     dimension ( : , : ) , allocatable :: grr 
  integer,           dimension ( : )     , allocatable :: nr 
  character(len=15), dimension ( : )     , allocatable :: cint
  real(kind=dp)                                        :: rr , vol
  real(kind=dp)                                        :: ttt1 , ttt2      
  ! =====================================================
  !   type of positions coordinates 
  ! =====================================================
  logical :: allowed
  character(len=60), SAVE :: cpos
  character(len=60), SAVE :: cpos_allowed(4)
  data cpos_allowed / 'Direct' , 'D' , 'Cartesian' , 'C' /


  ! trash 
  integer            :: iiii
  real(kind=dp)      :: aaaa
  character(len=60)  :: cccc


  OPEN (UNIT = kunit_TRAJFF ,FILE = 'TRAJFF') 
  OPEN ( kunit_GRTFF , file = 'GRTFF' )
  OPEN ( kunit_NRTFF , file = 'NRTFF' )

  READ ( kunit_TRAJFF , * ) natm
  READ ( kunit_TRAJFF , * ) system
  READ ( kunit_TRAJFF , * ) simu_cell%A ( 1 , 1 ) , simu_cell%A ( 2 , 1 ) , simu_cell%A ( 3 , 1 )
  READ ( kunit_TRAJFF , * ) simu_cell%A ( 1 , 2 ) , simu_cell%A ( 2 , 2 ) , simu_cell%A ( 3 , 2 )
  READ ( kunit_TRAJFF , * ) simu_cell%A ( 1 , 3 ) , simu_cell%A ( 2 , 3 ) , simu_cell%A ( 3 , 3 )
  READ ( kunit_TRAJFF , * ) ntype
  READ ( kunit_TRAJFF ,* ) ( atypei ( it ) , it = 1 , ntype )
  IF ( ionode ) WRITE ( stdout      ,'(A,20A3)' ) 'found type information on TRAJFF : ', atypei ( 1:ntype )
  READ( kunit_TRAJFF ,*)   ( natmi ( it ) , it = 1 , ntype )
  READ( kunit_TRAJFF ,*) cpos
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

  CALL lattice ( simu_cell ) 
  rho = DBLE ( natm )  / simu_cell%omega 

  CALL gr_init

  CALL print_general_info( stdout )

  ! ===================================
  !  here we know natm, then alloc 
  !  and decomposition can be applied 
  ! ================================== 
  CALL config_alloc 
  CALL do_split ( natm , myrank , numprocs , iastart , iaend , 'atoms' )
  CALL gr_alloc

  pairs =  ntype * ( ntype + 1 ) / 2
  allocate ( grr ( 0 : PANGR , 0 : pairs ) , nr ( 0 : pairs ) , cint ( 0 : pairs  ))
  grr  = 0.0_dp
  nr   = 0
  cint = ''

#ifdef debug2
  if ( ionode ) then 
    WRITE ( stdout , '(a,2i6)' ) 'debug : iastart, iaend ',iastart , iaend
    WRITE ( stdout , '(a,2i6)' ) 'debug : number of type pairs ', pairs
  endif
#endif

  CALL typeinfo_init

  ! ==========================================   
  !  skip the first nskip configurations 
  ! ==========================================
  if (nskip.gt.0) then
    do ic = 1,nskip
      if ( ic .ne. 1 ) READ ( kunit_TRAJFF , * ) iiii   
      if ( ic .ne. 1 ) READ ( kunit_TRAJFF , * ) cccc
      if ( ic .ne. 1 ) READ ( kunit_TRAJFF , * ) aaaa   ,  aaaa ,  aaaa
      if ( ic .ne. 1 ) READ ( kunit_TRAJFF , * ) aaaa   ,  aaaa ,  aaaa
      if ( ic .ne. 1 ) READ ( kunit_TRAJFF , * ) aaaa   ,  aaaa ,  aaaa
      if ( ic .ne. 1 ) READ ( kunit_TRAJFF , * ) iiii
      if ( ic .ne. 1 ) READ ( kunit_TRAJFF , * ) ( cccc , it = 1 , ntype )
      if ( ic .ne. 1 ) READ ( kunit_TRAJFF , * ) ( iiii , it = 1 , ntype )
      if ( ic .ne. 1 ) READ ( kunit_TRAJFF , * ) cpos
      do ia = 1 , natm 
        READ ( kunit_TRAJFF , * ) atype ( ia ) , rx ( ia ) , ry ( ia ) , rz ( ia ) , aaaa,aaaa,aaaa,aaaa,aaaa,aaaa
      enddo      
    enddo
  endif

  ngr = 0
  do ic = nskip + 1, nconf
    io_node WRITE ( stdout , '(a,i6,a,i6,a)' ) 'config : [ ',ic,' / ',nconf,' ] '
    ! ===================================
    !  read config from trajectory file
    ! ===================================
    if ( ic .ne. (nskip + 1) .or. nskip .ne. 0 ) then
      READ ( kunit_TRAJFF , * ) iiii
      READ ( kunit_TRAJFF , * ) cccc
      READ ( kunit_TRAJFF , * ) aaaa   ,  aaaa ,  aaaa
      READ ( kunit_TRAJFF , * ) aaaa   ,  aaaa ,  aaaa
      READ ( kunit_TRAJFF , * ) aaaa   ,  aaaa ,  aaaa
      READ ( kunit_TRAJFF , * ) iiii
      READ ( kunit_TRAJFF , * ) ( cccc , it = 1 , ntype )
      READ ( kunit_TRAJFF , * ) ( iiii , it = 1 , ntype )
      READ ( kunit_TRAJFF , * ) cpos
    endif
    do ia = 1 , natm
      READ ( kunit_TRAJFF , * ) atype ( ia ) , rx ( ia) , ry ( ia ) , rz ( ia ) , aaaa,aaaa,aaaa,aaaa,aaaa,aaaa
    enddo
    if ( cpos .eq. 'Direct' ) then
      ! ======================================
      !         direct to cartesian
      ! ======================================
      CALL dirkar ( natm , rx , ry , rz , simu_cell%A )
      if ( ionode .and. ic .eq. 1 ) WRITE ( stdout      ,'(A,20A3)' ) 'atomic positions in direct coordinates in POSFF'
    else if ( cpos .eq. 'Cartesian' ) then
      if ( ionode .and. ic .eq. 1 ) WRITE ( stdout      ,'(A,20A3)' ) 'atomic positions in cartesian coordinates in POSFF'
    endif

    ngr=ngr+1 
    ! ==========================
    !  calc radial_distribution 
    ! ==========================  
    call gr_main ( iastart , iaend )

  enddo !nconf 

  ! ===========================================
  !        merge results  
  ! ===========================================
  ttt1 = MPI_WTIME(ierr)

  CALL MPI_ALL_REDUCE_INTEGER ( gr(:,0,0), PANGR )
  do it1 = 1 , ntype
    do it2 = it1 , ntype
      CALL MPI_ALL_REDUCE_INTEGER ( gr(:,it1,it2), PANGR )
    enddo  
  enddo
  ttt2 = MPI_WTIME(ierr)
  grtimetot_comm = grtimetot_comm + ( ttt2 - ttt1 )


#ifdef debug2
  do i=0, PANGR
    io_node WRITE (stdout , '(a,5i6)') 'debug ( total ) : ',i,gr(i,1,1)
  enddo
#endif

  ! ======================================================= 
  !  write output files GRTFF , NRTFF
  ! ======================================================= 
  cint ( 0) = atypei ( 0 )//' - '//atypei ( 0 )
  mp = 1
  do it1 = 1 , ntype
    do it2 = it1 , ntype
      cint(mp) = atypei(it1)//' - '//atypei(it2)
      mp = mp + 1 
    enddo
  enddo

  WRITE ( kunit_GRTFF , '(<pairs+2>a)' ) '#       rr           ', ( cint ( mp ) , mp = 0 , pairs ) 
  WRITE ( kunit_NRTFF , '(<pairs+2>a)' ) '#       rr           ', ( cint ( mp ) , mp = 0 , pairs ) 

  nr ( 0 ) = 0 

  do i = 1 , PANGR
    rr  = ( dble(i)-0.5d0)*resg
    vol = 4.d0*pi*(resg*rr*rr+(resg**3)/12.d0)
    vol = vol/simu_cell%omega 
!    rr = resg * ( DBLE (i) + 0.5_dp )
!    k  = i+1
!    k  = k*k*k
!    k  = k-(i*i*i)
!    vol = DBLE (k)*resg*resg*resg ! r^3
!    vol = (4.0_dp/3.0_dp)*pi*vol/simu_cell%omega !4/3pir^3
    ! all - all pairs 
    grr ( i , 0 ) = DBLE (gr(i,0,0))/(ngr*vol*natm*natm)
    ! type pairs
    mp = 1
    do it1 = 1 , ntype
      do it2 = it1 , ntype
#ifdef debug2
        WRITE ( stdout , '(a,3i5)' ) 'debug ( pair ) : ', mp , it1 , it2
#endif        
        if ( mp .lt. 0 .and. mp .gt. pairs ) then
          WRITE ( stderr , '(a)' ) 'ERROR out of bound of gr in gr_main'
          STOP
        endif 
        nr ( mp ) = it1          
        grr ( i , mp ) = DBLE ( gr ( i , it1 , it2 ) ) / DBLE ( ngr * vol * natmi ( it1 ) * natmi ( it2 ) )  
        mp = mp + 1
      enddo
    enddo
    !if ( pairs .ne. 1 ) then
      if ( ionode ) then
        WRITE ( kunit_GRTFF ,'(<pairs+2>f15.10)') rr , ( grr ( i , mp ) , mp = 0 , pairs ) 
        WRITE ( kunit_NRTFF ,'(<pairs+2>f20.10)') rr , ( DBLE ( grr ( i , mp) ) * 4.0_dp * pi * rr * rr * DBLE ( natmi(nr(mp)) * vol ) , mp = 0 , pairs )
      endif
    !else
    !  if ( ionode ) then
    !    WRITE (kunit_GRTFF,'(2f15.10)') rr , grr( i , 0 ) 
    !    WRITE (kunit_NRTFF,'(f15.10,4i8)')  rr , gr(i,0,0)
    !  endif
    !endif
  enddo


  CLOSE ( kunit_NRTFF )
  CLOSE ( kunit_GRTFF )

  CLOSE( kunit_TRAJFF )

  !CALL static_struc_fac ( grr , PANGR , pairs ) 

  deallocate ( grr , nr , cint )
  CALL gr_dealloc

  return

END SUBROUTINE grcalc

! *********************** SUBROUTINE gr_main ***********************************
!
!> \brief
!! based on Frenkel and Smit
!
! ******************************************************************************
SUBROUTINE gr_main ( iastart , iaend )

  USE control,                  ONLY :  myrank
  USE config,                   ONLY :  natm , natmi , rx , ry , rz , atype , simu_cell , ntype , itype
  USE io_file,                  ONLY :  ionode , stdout  , stderr
  USE time,                     ONLY :  grtimetot
  USE cell,                     ONLY :  kardir , dirkar

  implicit none
  INCLUDE 'mpif.h'

  ! global
  integer, intent(in) :: iastart , iaend

  ! local
  integer :: ia , ja , ierr , ita , jta 
  integer :: igr 
  real(kind=dp) :: cut2 , rijsq , rr 
  real(kind=dp) :: rxi , ryi , rzi
  real(kind=dp) :: rxij , ryij , rzij
  real(kind=dp) :: sxij , syij , szij
  real(kind=dp) :: ttt1 , ttt2      

  ttt1 = MPI_WTIME(ierr)

  ! =======================
  !  cut-off half box
  ! =======================
  cut2 = cutgr * cutgr 

  ! ======================================
  !         cartesian to direct 
  ! ======================================
  CALL kardir ( natm , rx , ry , rz , simu_cell%B )

!  print*,iastart,iaend
  do ia = iastart , iaend
    rxi = rx ( ia )
    ryi = ry ( ia )
    rzi = rz ( ia )
    ita = itype( ia )
    do ja = 1, natm
      if ( ja .ne. ia ) then  
        jta = itype( ja )
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
        if ( rijsq.lt.cut2 ) then
          rr = SQRT ( rijsq )
          igr = INT ( rr / resg ) + 1
!          print*,igr
          if ( igr .lt. 0 .and. igr .gt. PANGR-1 ) then
            WRITE ( stderr , '(a)' ) 'ERROR out of bound of gr in gr_main'
            STOP
          endif 
          ! all pairs
          gr ( igr , 0 , 0 )     = gr ( igr , 0 , 0 ) + 1 
          gr ( igr , ita , jta ) = gr ( igr , ita , jta ) + 1
        endif
      endif
    enddo
  enddo
  
#ifdef debug2
  do igr=0, PANGR-1
    WRITE (stdout , '(a,5i6)') 'debug: ',myrank,igr,gr(igr,1,1)
  enddo
#endif 

  ttt2 = MPI_WTIME(ierr)
  grtimetot = grtimetot + ( ttt2 - ttt1 )

  ! ======================================
  !         direct to cartesian
  ! ======================================
  CALL dirkar ( natm , rx , ry , rz , simu_cell%A )

  return
 
END SUBROUTINE gr_main

! *********************** SUBROUTINE static_struc_fac **************************
!
!
! ******************************************************************************
!SUBROUTINE static_struc_fac ( gr , PANGR , pairs )

!  USE io_file,                  ONLY :  ionode , kunit_STRFACFF , stdout 
!  USE config,                   ONLY :  rho
!  USE constants,                ONLY :  pi , tpi , imag

!  implicit none
!  ! global
!  integer :: PANGR , pairs
!  real(kind=dp) :: gr ( 0 : PANGR-1 , 0 : pairs ) , rr
!
!  ! local
!  integer :: i , j , is , NN
!  real(kind=dp) :: q , qj , ri , rip
!  complex(kind=dp) :: Sk
!  real(kind=dp) , dimension (:) , allocatable  :: stat_str
!  real(kind=dp) , dimension (:) , allocatable :: in 
!  real(kind=dp) , dimension (:,:) , allocatable :: Uji 
!  complex(kind=dp)   , dimension (:) , allocatable :: out
!  real(kind=dp) :: res , shift
!  real(kind=dp) :: x , k
!
!  io_node WRITE ( stdout , '(a)' ) 'in static_struc_fac'
!
!  allocate ( in ( PANGR ) , out ( PANGR /2 + 1 ) )
!
!  in  = ( 0.0,0.0)
!  out = ( 0.0,0.0)
!  ! ========
!  !   FFT
!  ! ========
!  do i=0,PANGR-1
!    in ( i + 1 ) = gr ( i , 0 ) 
!  enddo
!
!!  CALL fft_1D_complex ( in , out , PANGR )
!  CALL fft_1D_real(in,out,PANGR)
!
!  do i= 1 , PANGR/2+1
!    q = ( dble ( i )  + 0.5_dp ) / DBLE ( PANGR ) / resg
!    Sk = 1.0_dp + rho * out( i + 1 )  
!    io_node WRITE ( 20000 , '(3e16.8)' )  q , Sk  
!  enddo
!
!  deallocate ( in , out )
!
!! other version
!  ! Uji (eq 12) J Phys Cndens Matter V17 2005 )
!!  allocate ( Uji ( PANGR , PANGR ) ) 
!
!  do i = 1 , PANGR
!!    ri  = ( dble ( i )  + 0.5_dp )  * res
!    rip = ( dble ( i + 1 )  + 0.5_dp )  * res
!!    do j = 1 , PANGR
!      qj = tpi * DBLE ( j ) + 0.5_dp / DBLE ( PANGR / 2 + 1 ) / resg
!      Uji ( j , i ) = SR ( ri , qj ) - SR ( rip , qj ) 
! !     Uji ( j , i ) = Uji ( j , i ) / qj  
!    enddo
! ! enddo
!  Uji = 2.0_dp * tpi * Uji
!!
!  do i= 1 , PANGR/2+1
!    q = tpi * ( dble ( i )  + 0.5_dp ) / DBLE ( PANGR / 2 + 1) / resg
!!    do j = 1 , PANGR
!      Sk =  Sk + Uji ( j , i ) * ( gr ( j , 0 ) -1.0_dp )  
!    enddo
!    io_node WRITE ( 30000 , '(3e16.8)' )  q , Sk
!  enddo
!
!
!  deallocate ( Uji )

! test purpose because I'm dumb
! I was enable to get k in "real life" 
! I did a simple discret case  ( N = 4 ) 
! to find the relation between k and q 

!  NN = 4
!  allocate ( in ( NN ) , out  ( NN ) )
!
!  in = ( 0.0,0.0)
!  is = 3 
!  in(is) = ( 1.0,0.0)
!  print*,'in in '
!  CALL fft_1D_complex ( in , out , NN )
!
!  res = 2.0_dp
!  shift= (is-1) * res
!  write( stdout , '(8a)' ) '            x       Re in(i)        Im in(i)            k          Re out(i)       Im out(i)        Re phase        Im phase'
!  do i=1,NN
!    x = DBLE(i-1)*res
!    k = DBLE(i-1) / DBLE ( NN ) 
!    q = DBLE(i-1) / DBLE ( NN ) / res
!    write( stdout , '(8e16.8)' ) x , in(i) , k , out(i) , exp ( -imag * tpi * k * (is-1) )
!    write( stdout , '(8e16.8)' ) x , in(i) , q , out(i) , exp ( -imag * tpi * q * shift )
!  enddo
!
!  deallocate ( in ,out )

!  return
!CONTAINS

!real(kind=dp) function SR(r,qj)
!  implicit none
!  real(kind=dp) :: r , qj  
!  SR = SIN ( qj * r ) / qj / qj 
!  SR = SR - r * COS ( qj * r ) / qj
!end function 

!END SUBROUTINE static_struc_fac


END MODULE radial_distrib 
! ===== fmV =====
