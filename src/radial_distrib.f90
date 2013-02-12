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
!#define debug2
! ======= Hardware =======

MODULE radial_distrib 

  integer :: PANGR             ! (internal) number of bins in g(r) distribution
  integer :: nskip
  integer :: nconf
  double precision :: cutgr
 
  double precision :: resg     ! resolution in g(r) distribution 

  integer, dimension(:,:,:), allocatable :: gr 

CONTAINS

!*********************** SUBROUTINE grcalc_init *******************************
!
!
!******************************************************************************

SUBROUTINE gr_init

  USE config,   ONLY :  simu_cell
  USE control,  ONLY :  calc
  USE io_file,  ONLY :  stdin, stdout, kunit_OUTFF, ionode

  implicit none

  ! local
  integer :: ioerr , npangr , i
  character * 132 :: filename

  namelist /grtag/   nconf , &
                     nskip , &
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
   if ( ionode ) WRITE ( stdout, '(a)') 'ERROR reading input_file : vacftag section is absent'
   STOP
  elseif ( ioerr .gt. 0 )  then
   if ( ionode ) WRITE ( stdout, '(a)') 'ERROR reading input_file : vacftag wrong tag'
   STOP
  endif

  CLOSE ( stdin )

! define a new resolution 2^N points
  cutgr=0.5d0 * MIN(simu_cell%WA,simu_cell%WB,simu_cell%WC)-0.1d0
  PANGR=int(cutgr/resg)
  i = 1
  do while ( 2**i .lt. PANGR )
     i = i + 1
  enddo
  npangr = i
  PANGR = 2** npangr
  resg = cutgr / DBLE ( PANGR  )




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
  USE config ,  ONLY : ntype

  implicit none

  if ( calc .ne. 'gr' ) return

  allocate(gr(0:PANGR-1,0:ntype,0:ntype))
  gr = 0      

  return 
 
END SUBROUTINE gr_alloc


!*********************** SUBROUTINE gr_dealloc **********************************
!
!
!******************************************************************************

SUBROUTINE gr_dealloc

  USE control,  ONLY :  calc

  implicit none

  if ( calc .ne. 'gr' ) return
  
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
  nconf = 0

  return

END SUBROUTINE gr_default_tag


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
                  WRITE ( kunit ,'(a,f10.5,a)')         'resolution of g(r) function resg     = ',resg,' new value to get 2^N points in g(r)'
                  WRITE ( kunit ,'(a,i5)')              'number of points in g(r)             = ',PANGR
                  WRITE ( kunit ,'(a)')                 'save radial_distribution in file     :   GRTFF' 
      if ( calc .eq. 'gr' )     then 
                  WRITE ( kunit ,'(a)')                 'read configuration from file         :   TRAJFF'
                  WRITE ( kunit ,'(a)')                 ''
                  WRITE ( kunit ,'(a,i5)')              'number of config. in TRAJFF          = ',nconf        
                  WRITE ( kunit ,'(a,i5)')              'number of config. to be skipped      = ',nskip
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
                                rho , config_alloc , simu_cell , atypei , itype, natmi
  USE control,          ONLY :  myrank , numprocs
  USE io_file,          ONLY :  ionode , stdout , kunit_TRAJFF , kunit_OUTFF , kunit_GRTFF , kunit_NRTFF
  USE constants,        ONLY :  pi 
  USE cell
  USE time

  implicit none
  INCLUDE 'mpif.h'

  ! local 
  integer :: ia , ic , it , ngr , i , k , pairs , it1 , it2 , mp , ierr 
  integer :: iastart , iaend 
  double precision, dimension ( : , : ) , allocatable :: grr 
  integer, dimension ( : ) , allocatable :: nr 
  character(len=15), dimension ( : ) , allocatable :: cint
  double precision :: rr , vol
  ! trash 
  double precision :: aaaa
  integer :: iiii
  character * 60 :: cccc
  double precision :: grtime1 , grtime2      


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
  IF ( ionode ) WRITE ( kunit_OUTFF ,'(A,20A3)' ) 'found type information on TRAJFF : ', atypei ( 1:ntype )
  IF ( ionode ) WRITE ( stdout      ,'(A,20A3)' ) 'found type information on TRAJFF : ', atypei ( 1:ntype )
  READ( kunit_TRAJFF ,*)   ( natmi ( it ) , it = 1 , ntype )

  CALL lattice ( simu_cell ) 
  rho = DBLE ( natm )  / simu_cell%omega 

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
  pairs =  ntype * ( ntype + 1 ) / 2
  allocate ( grr ( 0 : PANGR-1 , 0 : pairs ) , nr ( 0 : pairs ) , cint ( 0 : pairs  ))
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
      do ia = 1 , natm 
        READ ( kunit_TRAJFF , * ) atype ( ia ) , rx ( ia ) , ry ( ia ) , rz ( ia ) , aaaa,aaaa,aaaa,aaaa,aaaa,aaaa
      enddo      
    enddo
  endif

  ngr = 0
  do ic = nskip + 1, nconf
    if ( ionode ) WRITE ( stdout , '(a,i6,a,i6,a)' ) 'config : [ ',ic,' / ',nconf,' ] '
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
    endif
    do ia = 1 , natm
      READ ( kunit_TRAJFF , * ) atype ( ia ) , rx ( ia) , ry ( ia ) , rz ( ia ) , aaaa,aaaa,aaaa,aaaa,aaaa,aaaa
    enddo
    ngr=ngr+1 
    ! ==========================
    !  calc radial_distribution 
    ! ==========================  
    call gr_main ( iastart , iaend )

  enddo !nconf 


  ! ===========================================
  !        merge results  
  ! ===========================================
  grtime1 = MPI_WTIME(ierr)

  CALL MPI_ALL_REDUCE_INTEGER ( gr(:,0,0), PANGR )
  do it1 = 1 , ntype
    do it2 = it1 , ntype
      CALL MPI_ALL_REDUCE_INTEGER ( gr(:,it1,it2), PANGR )
    enddo  
  enddo
  grtime2 = MPI_WTIME(ierr)
  grtimetot_comm = grtimetot_comm + ( grtime2 - grtime1 )


#ifdef debug2
  do i=0, PANGR
    if ( ionode ) WRITE (stdout , '(a,5i6)') 'debug ( total ) : ',i,gr(i,1,1)
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

  do i = 0 , PANGR-1
    rr = resg * ( DBLE (i) + 0.5D0 )
    k  = i+1
    k  = k*k*k
    k  = k-(i*i*i)
    vol = DBLE (k)*resg*resg*resg ! r^3
    vol = (4.0D0/3.0D0)*pi*vol/simu_cell%omega !4/3pir^3
    ! all - all pairs 
    grr ( i , 0 ) = DBLE (gr(i,0,0))/(ngr*vol*natm*natm)
    ! type pairs
    mp = 1
    do it1 = 1 , ntype
      do it2 = it1 , ntype
#ifdef debug2
       WRITE ( stdout , '(a,3i5)' ) 'debug ( pair ) : ', mp , it1 , it2
#endif        
        nr ( mp ) = it1          
        grr ( i , mp ) = DBLE ( gr ( i , it1 , it2 ) ) / DBLE ( ngr * vol * natmi ( it1 ) * natmi ( it2 ) )  
        mp = mp + 1
      enddo
    enddo
    if ( pairs .ne. 1 ) then
      if ( ionode ) then
        WRITE ( kunit_GRTFF ,'(<pairs+2>f15.10)') rr , ( grr ( i , mp ) , mp = 0 , pairs ) 
        WRITE ( kunit_NRTFF ,'(<pairs+2>f20.10)') rr , ( DBLE ( grr ( i , mp) ) * 4.0d0 * pi * rr * rr * DBLE ( natmi(nr(mp)) * vol ) , mp = 0 , pairs )
      endif
    else
      if ( ionode ) then
        WRITE (kunit_GRTFF,'(2f15.10)') rr , grr( i , 0 ) 
        WRITE (kunit_NRTFF,'(f15.10,4i8)')  rr , gr(i,0,0)
      endif
    endif
  enddo


  CLOSE ( kunit_NRTFF )
  CLOSE ( kunit_GRTFF )

  CLOSE( kunit_TRAJFF )

  CALL static_struc_fac ( grr , PANGR , pairs ) 
  

  deallocate ( grr , nr , cint )
  CALL gr_dealloc

  return

END SUBROUTINE grcalc

!*********************** SUBROUTINE gr_main ***********************************
!
! based on Frenkel and Smit
!
!******************************************************************************

SUBROUTINE gr_main ( iastart , iaend )

  USE control, ONLY : myrank
  USE config,           ONLY :  natm , natmi , rx , ry , rz , atype , simu_cell , ntype , itype
  USE prop,             ONLY :  nprop, nprop_start
  USE io_file,          ONLY :  ionode , stdout 
  USE time

  implicit none
  INCLUDE 'mpif.h'

  ! global
  integer, intent(in) :: iastart , iaend

  ! local
  integer :: ia , ja , ierr , ita , jta  
  integer :: igr 
  integer :: nxij , nyij , nzij
  double precision :: cut2 , rijsq , rr 
  double precision :: rxi , ryi , rzi
  double precision :: rxij , ryij , rzij
  double precision :: grtime1 , grtime2      

  grtime1 = MPI_WTIME(ierr)

  ! =======================
  !  cut-off half box
  ! =======================
  cut2 = cutgr * cutgr 

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
      nxij = NINT ( rxij * simu_cell%BNORM(1) )
      nyij = NINT ( ryij * simu_cell%BNORM(2) )
      nzij = NINT ( rzij * simu_cell%BNORM(3) )
      rxij = rxij - simu_cell%ANORM(1) * nxij
      ryij = ryij - simu_cell%ANORM(2) * nyij
      rzij = rzij - simu_cell%ANORM(3) * nzij
      rijsq = rxij * rxij + ryij * ryij + rzij * rzij
      if ( rijsq.lt.cut2 ) then
        rr = SQRT ( rijsq )
        igr = INT ( rr / resg )
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

  grtime2 = MPI_WTIME(ierr)
  grtimetot = grtimetot + ( grtime2 - grtime1 )

  return
 
END SUBROUTINE gr_main

!*********************** SUBROUTINE static_struc_fac **************************
!
!
!******************************************************************************

SUBROUTINE static_struc_fac ( gr , PANGR , pairs )

  USE io_file,  ONLY :  ionode , kunit_STRFACFF , stdout 
  USE config ,  ONLY :  rho
  USE constants,ONLY :  tpi        

  implicit none
  ! global
  integer :: PANGR , pairs
  double precision :: gr ( 0 : PANGR-1 , 0 : pairs ) , rr

  ! local
  integer :: i , k 
  double precision :: q 
  double complex :: Sk
!  double precision , dimension (:) , allocatable  :: stat_str
  double complex   , dimension (:) , allocatable :: out , in 

  if ( ionode ) WRITE ( stdout , '(a)' ) 'in static_struc_fac'

  allocate ( in ( PANGR ) , out ( PANGR ) )

  ! ========
  !   FFT
  ! ========
  do i=0,PANGR-1
    in ( i + 1 ) = gr ( i , 0 )
  enddo

  CALL fft_1D_complex ( in , out , PANGR )

  do i= 0 , PANGR-1
    q = ( dble ( i )  + 0.5d0 ) * tpi / resg
    Sk = 1.0d0 + rho * out( i + 1)  
    rr = resg * ( DBLE (i) + 0.5D0 )
    if ( ionode ) WRITE ( 20000 , '(5e16.8)' )  rr , q , Sk , gr ( i , 0 )
  enddo

  deallocate ( in , out )

  return

END SUBROUTINE static_struc_fac


END MODULE radial_distrib 
! ===== fmV =====
