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
!*********************** MODULE CONF ******************************************
!
! This module should deal with everything related to the actual configuration
! positions, velocities, forces, density , box .... 
!
!******************************************************************************

MODULE config

  implicit none

  integer, PARAMETER :: ntypemax=10

  character*60, SAVE :: system                                       ! system name                                              

  logical :: lfcc   
  logical :: lsc   ! not tested
  logical :: lbcc  ! not yet

  integer :: natm                                                    ! nb of atoms
  integer :: ntype                                                   ! number of types
  integer :: ncell                                                   ! number of cell (with lfcc)
  integer, dimension(:), allocatable :: itype                        ! type of atome i array 
  integer, dimension(:), allocatable :: list, point                  ! vnlist info
  integer, dimension(0:ntypemax)     :: natmi                        ! number of atoms (per type)

  double precision :: box                                            ! cell dimension
  double precision :: omega                                          ! volume
  double precision :: rho                                            ! density  
  double precision :: xna, xnb                                       ! relative composition 

  double precision , dimension(:), allocatable :: rx  , ry  , rz     ! positions
  double precision , dimension(:), allocatable :: vx  , vy  , vz     ! velocities
  double precision , dimension(:), allocatable :: fx  , fy  , fz     ! forces
  double precision , dimension(:), allocatable :: fxs , fys , fzs    ! forces (previous step t-dt) beeman
  double precision , dimension(:), allocatable :: rxs , rys , rzs    ! previous positions for leap-frog integrator
  double precision , dimension(:), allocatable :: xs  , ys  , zs     ! last positions in verlet list
  double precision , dimension(:), allocatable :: rix , riy , riz    ! positions in the center of mass reference 

  double precision , dimension(:), allocatable :: qia                ! charge of ion 
  double precision , dimension(0:ntypemax)     :: qit                ! charge of ion type  

  character*10, dimension(:), allocatable      :: atype              ! atom type A or B 
  character*10, dimension(0:ntypemax)          :: atypei             ! type of atoms (per type)

  character*60 :: struct                                             ! crystal structure for fcc
  character*60 :: struct_allowed(2)                     
  data struct_allowed / 'NaCl' , 'random' /

  logical :: lgenconf                                                ! if .true. generate the configuration
                                                                     ! if not all config parameters are readed in
                                                                     !   md       opt/efg      vib
                                                                     !  POSFF  /  TRAJFF  /   ISCFF 

CONTAINS


!*********************** SUBROUTINE config_init *******************************
!
! set default values, read and check consistenstency of conifig parameters
!
!******************************************************************************

SUBROUTINE config_init 

  USE control,  ONLY :  calc
  USE io_file,  ONLY :  ionode , stdin, stdout , kunit_OUTFF

  implicit none

  ! local
  character * 132 :: filename
  integer :: ioerr

  namelist /configtag/       system   , & 
                             rho      , &
                             box      , &
                             ntype    , & 
                             natm     , &
                             struct   , & 
                             lfcc     , &
                             lsc      , &
                             lbcc     , & 
                             ncell    , &  
                             xna      , &
                             xnb      , &
                             lgenconf 

  ! ===============================
  ! defaults values for config tags
  ! ===============================
  CALL config_default_tag 
  ! ===============================
  ! reads config tags
  ! ===============================
  CALL getarg (1,filename)
  OPEN  ( stdin , file = filename)
  READ  ( stdin , configtag , iostat= ioerr)
  if( ioerr .lt. 0 )  then
    if( ionode ) WRITE ( stdout, '(a)') 'ERROR reading input_file : configtag section is absent'
    STOP
  elseif( ioerr .gt. 0 )  then
    if( ionode ) WRITE ( stdout, '(a)') 'ERROR reading input_file : configtag wrong tag'
    STOP
  endif
  CLOSE ( stdin )

  ! ===============================
  !  check values for config tags
  ! ===============================
  CALL config_check_tag

  if ( lgenconf ) then
    ! ===============================
    !  allocate princiapl arrays
    ! ===============================
    CALL config_alloc
    ! ===========================================================
    !  generate initial configuration only for calc = md
    !  for opt, vib and efg configuations are read in other files
    ! ===========================================================
    CALL gen_pos 
  else
    ! ===========================================================
    ! read initial configuration only for calc = md
    ! for opt, vib and efg configuations are read in other files
    ! ===========================================================
    if( calc .eq. 'md' ) then       
      CALL read_pos
    endif
  endif

  if( calc .ne. 'md' ) return
  ! ===============================
  ! print info for config tags
  ! ===============================
  CALL config_print_info(stdout)
  CALL config_print_info(kunit_OUTFF)
  ! ===============================
  ! write the starting configuration 
  ! ===============================
  CALL write_CONTFF

  return
 
END SUBROUTINE config_init 


!*********************** SUBROUTINE config_default_tag ************************
!
! set default values to config tag
!
!******************************************************************************

SUBROUTINE config_default_tag

  implicit none

  ! =================
  !  default values
  ! =================
  system   = 'UNKNOWN'
  lgenconf = .true.     
  xna      = 1.0D0
  xnb      = 0.0D0
  ntype    = 1
  ncell    = 4
  lfcc     = .true. 
  lsc      = .false.
  lbcc     = .false.
  rho      = 1.0d0
  box      = 0.0d0 
  struct   = 'random'
  natm     = 0
      
  return

END SUBROUTINE config_default_tag


!*********************** SUBROUTINE config_check_tag **************************
!
! check config tag values
!
!******************************************************************************

SUBROUTINE config_check_tag

  USE control,  ONLY :  cutoff
  USE io_file,  ONLY :  ionode , stdout


  implicit none
  
  ! local
  logical :: allowed
  integer :: i
 
  if ( .not. lgenconf ) then
    if( ionode ) WRITE ( stdout ,'(a)') '============================================================='
    if( ionode ) WRITE ( stdout ,'(a)') 'WARNING: with lgenconf=.false.'
    if( ionode ) WRITE ( stdout ,'(a)') 'all other config flags are not used '
  endif
  ! ===================
  !   generate config
  ! ===================
  if ( .not. lgenconf ) return 

  ! ============
  !  struct check
  ! ============
  do i = 1 , size (struct_allowed)
   if (trim (struct) .eq. struct_allowed(i) ) allowed = .true.
  enddo
  if(  .not. allowed ) then
    if ( ionode )  WRITE ( stdout ,'(a)') 'ERROR configtag: struct should be ', struct_allowed
    STOP 
  endif
  ! ============
  !  lfcc check
  ! ============
  if(lfcc .and. ncell.eq.0) then
    if( ionode ) WRITE ( stdout ,'(a)') 'ERROR configtag: with lfcc, ncell is needed'
    STOP 
  endif 
  if(struct.eq.'NaCl') then
     lfcc = .true.
     xna = 0.5D0
     xnb = 0.5D0
  endif
  if(struct.eq.'random'.and.ntype.eq.2) then
     if(xna.eq.0.0d0 .or. xnb .eq. 0.0d0) then
        if( ionode ) WRITE ( stdout ,'(a)') 'ERROR configtag: with struct = "random" and 2 types. xna and xnb are needed'
        STOP 
     endif 
  endif
  if(struct.eq.'NaCl'.and.ntype.ne.2) then
    if( ionode ) WRITE ( stdout ,'(a)') 'ERROR configtag: with struct = "NaCl" 2 types are needed'
    STOP 
  endif
  if ( lfcc ) then
     natm = 4*(ncell*ncell*ncell)
     if ( ntype .eq. 2 .and. struct.eq.'NaCl') then 
       natm = 8*(ncell*ncell*ncell)
     endif
  endif
  if ( lsc ) then
    natm = ncell * ncell * ncell    
  endif      
  ! ===================
  !  box and rho check
  ! ===================
  if(box.ne.0.0D0.and.rho.ne.0.0D0) then
    if( ionode ) WRITE ( stdout ,'(a)') 'ERROR configtag: box and rho cannot be set together (no mass) --- choose one !!'
    STOP 
  endif
  if(rho.ne.0.0D0.and.box.eq.0.0d0) then
    box = dble(natm) / rho
    box = box**(1.0D0 / 3.0D0 )
  else if(box.ne.0.0d0.and.rho.eq.0.0d0)  then
    rho = dble(natm)/(box**3.0D0)
  endif
  ! ===================
  !  check cutoff
  ! ===================
  if(cutoff.gt.box*0.5D0) then
    if( ionode )   WRITE ( stdout ,'(a)') ' '
    if( ionode )   WRITE ( stdout ,'(a)') 'WARNING configtag: cutoff is greater than half of cell size'
    if( ionode )   WRITE ( stdout ,'(a,f10.6,a2)') 'cutoff should be changed to a safer value, ex. ',box*0.499D0,' ?'
    if( ionode )   WRITE ( stdout ,'(a,f10.6,a,f10.6)') 'box =  ',box,'  cutoff =  ',cutoff
    cutoff = box*0.499D0
    if( ionode )   WRITE ( stdout ,'(a,f10.6)') 'Cutoff have been changed to continue',box*0.499D0
    if( ionode )   WRITE ( stdout ,'(a)') ' '
  endif
  ! ============
  !  natm check
  ! ============
  if ( natm .eq. 0 ) then
    if( ionode )   WRITE ( stdout ,'(a)') 'ERROR configtag: natm is 0'      
  endif

  return

END SUBROUTINE config_check_tag


!*********************** SUBROUTINE config_print_info *************************
!
! print information to standard output about the starting configuration
!
!******************************************************************************

SUBROUTINE config_print_info(kunit)

  USE io_file,  ONLY :  ionode 

  implicit none

  ! local
  integer :: kunit

  omega = box * box * box

  if( ionode ) then
                       WRITE ( kunit ,'(a,a)')          'system                               : ',system
                       WRITE ( kunit ,'(a,i10)')        'natm                                 = ',natm
       if(ntype.eq.2)  WRITE ( kunit ,'(a,i10,f8.2,a1)')'nA                                   = ',natmi(1),xna * 100,'%'
       if(ntype.eq.2)  WRITE ( kunit ,'(a,i10,f8.2,a1)')'nB                                   = ',natm-natmi(1),xnb * 100,'%'
                       WRITE ( kunit ,'(a,f10.4)')      'cell parameter                       = ',box
                       WRITE ( kunit ,'(a,f10.4)')      'volume                               = ',omega
                       WRITE ( kunit ,'(a,f10.4)')      'density                              = ',rho
  endif 
  ! ============================
  !  print minimum distance 
  !  and distance distribution
  ! ============================
  call distance_tab( kunit )

  return
  
END SUBROUTINE config_print_info


!*********************** SUBROUTINE gen_pos ***********************************
!
!  here we read configuration (pos,vel) from file or generate from a given
!  structure lfcc, lsc or lbcc ... (face centered cubic, simple cubic or ...)
!  confusing organisation !!!!!
!  copie a revoir ;)
!
!******************************************************************************

SUBROUTINE gen_pos

  USE io_file,  ONLY :  ionode , stdout , kunit_POSFF, kunit_OUTFF

  implicit none

  ! local
  integer :: i , natmtmp , na 
  double precision ,dimension (:) ,allocatable :: rxtmp,rytmp,rztmp

  IF ( ionode ) then
    WRITE ( kunit_OUTFF ,'(a)')       '=============================================================' 
    WRITE ( kunit_OUTFF ,'(a)')       ''
    WRITE ( kunit_OUTFF ,'(a)')       'structure generated from the code'
    WRITE ( stdout ,'(a)')            '=============================================================' 
    WRITE ( stdout ,'(a)')            ''
    WRITE ( stdout ,'(a)')            'structure generated from the code'
  endif

  allocate(rxtmp(natm),rytmp(natm),rztmp(natm))

  rxtmp = 0.0D0
  rytmp = 0.0D0
  rztmp = 0.0D0

  ! ================================      
  ! lfcc face centered cubic
  ! ================================      
  IF( lfcc ) THEN
    ! ======================================      
    ! if 2 types how to generate A/B 
    ! configuration according to xna and xnb
    ! ======================================      
    if(ntype.eq.2) then 


      ! ========
      !  random
      ! ========
      if (struct.eq.'random') then 
        CALL init_fcc('random', natm , ntype , ncell , rx , ry , rz , box , struct )
        CALL random_config ( natm , atype , xna )


      ! ========
      !  NaCl 
      ! ========
      elseif(struct.eq.'NaCl') then
         
        ! store the number of atom 
        natmtmp = natm
        ! only half of the atom for site A
        natm = natm * 0.5d0 
        ! generate sites A 
        CALL init_fcc('AAAAAA', natm , ntype , ncell , rx , ry , rz , box , struct ) 
        ! store sites A 
        do i = 1, natm
          rxtmp(i) = rx(i)
          rytmp(i) = ry(i)
          rztmp(i) = rz(i)  
        enddo
        ! generate sites B (half of the total number of atoms) 
        CALL init_fcc('BBBBBB', natm , ntype , ncell , rx , ry , rz , box , struct )
        ! store sites B 
        do i = 1, natm
          rxtmp(i + natm) = rx(i) + (box/ncell) * 0.5D0
          rytmp(i + natm) = ry(i) + (box/ncell) * 0.5D0
          rztmp(i + natm) = rz(i) + (box/ncell) * 0.5D0
        enddo
        do i = 1, natm 
          atype(i) = 'A'
          itype(i) = 1
        enddo
        do i = natm + 1, natmtmp
          atype(i) = 'B'
          itype(i) = 2
        enddo
        natm = natmtmp
        do i = 1, natm
          rx(i) = rxtmp(i)
          ry(i) = rytmp(i)
          rz(i) = rztmp(i)   
        enddo   
      endif
    ! ======================================      
    ! else all atoms are A ;)
    ! ======================================      
    else 
      CALL init_fcc('AAAAAA' , natm ,  ntype , ncell , rx , ry , rz , box , struct )      
      do i = 1, natm
        atype(i) = 'A'
        itype(i) = 1
      enddo
    endif
  ! ================================      
  ! lsc: simple cubic
  ! ================================      
  ELSEIF (lsc) THEN 

    ! ======================================      
    ! if 2 types how to generate A/B 
    ! configuration according to xna and xnb
    ! ======================================      
    if(ntype.eq.2) then 
      ! ========
      ! random 
      ! ========
      if (struct.eq.'random') then 
        CALL init_sc
        CALL random_config ( natm , atype , xna )
        print*,'check init_sc compare to fcc'
        STOP 
      endif
    ! ======================================      
    ! else all atoms are A ;)
    ! ======================================      
    else 
      CALL init_sc
      do i = 1, natm
        atype(i) = 'A'
        itype(i) = 1 
      enddo
    endif
  ! ================================      
  ! lbcc: body-centered cubic
  ! ================================      
  ELSEIF (lbcc) THEN ! lbcc

    if( ionode ) then
      WRITE ( stdout , * ) 'not yet'  
      STOP 
    endif
  ENDIF ! config     


  ! ==============================
  !  count A atoms and initialize 
  !  some type arrays
  ! ==============================      
  natmi = 0
  na = 0
  do i = 1 , natm
    if( atype(i) .eq. 'A' ) then 
      na = na + 1
      itype ( i ) = 1
    else
      itype ( i ) = 2 
    endif
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
   
  ! ======================      
  ! recalculate xna,xnb  
  ! ======================
  xna = dble(na)/dble(natm)
  xnb = 1.0D0-xna

  deallocate(rxtmp,rytmp,rztmp)

  return

END SUBROUTINE gen_pos

!*********************** SUBROUTINE read_pos **********************************
!
!  here we read configuration (pos,vel) from file or generate from a given
!  structure lfcc, lsc or lbcc ... (face centered cubic, simple cubic or ...)
!  organisation to confusing !!!!!
!
!******************************************************************************

SUBROUTINE read_pos 

  USE io_file,  ONLY :  ionode , stdout , kunit_POSFF, kunit_OUTFF

  implicit none

  ! local
  integer :: i , it , ia,na

  IF ( ionode ) then
    WRITE ( kunit_OUTFF ,'(a)')       '=============================================================' 
    WRITE ( kunit_OUTFF ,'(a)')       ''
    WRITE ( kunit_OUTFF ,'(a)')       'config from file POSFF'
    WRITE ( stdout ,'(a)')            '=============================================================' 
    WRITE ( stdout ,'(a)')            ''
    WRITE ( stdout ,'(a)')            'config from file POSFF'
  endif

  OPEN ( kunit_POSFF , file = 'POSFF' ) 
  READ ( kunit_POSFF ,*) system
  READ ( kunit_POSFF ,*) box , ntype 
  READ( kunit_POSFF ,*) (atypei(it),it=1,ntype)
  IF ( ionode ) WRITE( kunit_OUTFF ,'(A,20A3)') 'found type information on POSFF : ',atypei(1:ntype)
  IF ( ionode ) WRITE( stdout      ,'(A,20A3)') 'found type information on POSFF : ',atypei(1:ntype)
  READ( kunit_POSFF ,*) (natmi(it),it=1,ntype)
      
  omega = box * box * box
  natm=0
  do it=1,ntype
    natm = natm + natmi(it)
  enddo
  rho = dble(natm)/ omega

  call config_alloc 

  ! =========================================      
  ! read positions and velocities from disk 
  ! =========================================      
  READ  ( kunit_POSFF , * ) ( rx ( ia ) , ry ( ia ) , rz ( ia ) , vx( ia ) , vy ( ia ) , vz ( ia ) , ia = 1 , natm )
  CLOSE ( kunit_POSFF )

  ! ==========================
  !  set some type parameters
  ! ==========================      
  do i = 1 , natm
    if ( i .le. natmi(1) )  then
      atype ( i ) = 'A' 
      itype ( i ) =  1
    else
      atype ( i ) = 'B' 
      itype ( i ) =  2 
     endif
  enddo

  ! ======================      
  ! recalculate xna,xnb  
  ! ======================
  xna = dble(natmi(1))/dble(natm)
  xnb = 1.0D0-xna

  atypei(0)='ALL'
  natmi(0)=natm


  return

END SUBROUTINE read_pos

!*********************** SUBROUTINE write_CONTFF ******************************
!
! write configuration (pos,vel) to CONTFF file
!
!******************************************************************************

SUBROUTINE write_CONTFF

  USE io_file,  ONLY :  kunit_CONTFF, ionode

  implicit none

  ! local
  integer :: ia , it
  double precision, dimension (:) , allocatable :: xxx , yyy , zzz

  allocate ( xxx ( natm ) , yyy ( natm ) , zzz ( natm ) )

  xxx = rx
  yyy = ry
  zzz = rz

  CALL  periodicbc ( natm , xxx , yyy , zzz , box )
  
  if ( ionode ) then
  OPEN ( kunit_CONTFF ,file = 'CONTFF',STATUS = 'UNKNOWN')
      WRITE ( kunit_CONTFF,'(a)') system
      WRITE ( kunit_CONTFF,'(f20.12,i4)') box,ntype 
      WRITE ( kunit_CONTFF,*) ( atypei(it) , it=1,ntype ) 
      WRITE ( kunit_CONTFF,*) ( natmi (it) , it=1,ntype ) 
      WRITE ( kunit_CONTFF,'(6f20.12)') ( xxx ( ia ) , yyy ( ia ) , zzz ( ia ) , vx ( ia ) , vy ( ia ) , vz ( ia ) , ia = 1 , natm )
  CLOSE (kunit_CONTFF)
  endif

  deallocate ( xxx , yyy , zzz ) 

  return

END SUBROUTINE write_CONTFF


!*********************** SUBROUTINE config_alloc ******************************
!
! allocation of principal arrays of the calculation:
! 
! positions of atom ia                      : rx(ia)  , ry(ia)   , rz(ia) 
! velocities of atom ia                     : vx(ia)  , vy(ia)   , vz(ia) 
! forces on atom ia                         : fx(ia)  , fy(ia)   , fz(ia) 
! previous pos. leap-frog algo.             : rxs(ia) , rys(ia)  , rzs(ia) 

! atom type information (still not nice) 
!
! type of atom ia (character)               : atype(ia) = 'A' or 'B'
! type of atom ia (integer)                 : itype(ia) = 1 or 2
! number of atoms of type it                : natmi(it) .le. natm
! name of type of type it                   : atypei(it) = 'A' or 'B' or ALL
!
!******************************************************************************

SUBROUTINE config_alloc

  USE control,  ONLY :  calc , lvnlist

  implicit none

  allocate( rx  ( natm ) , ry ( natm )  , rz ( natm ) )
  allocate( vx  ( natm ) , vy ( natm )  , vz ( natm ) )
  allocate( fx  ( natm ) , fy ( natm )  , fz ( natm ) )
  allocate( fxs ( natm ) , fys ( natm ) , fzs ( natm ) )
  allocate( rxs ( natm ) , rys ( natm ) , rzs ( natm ) )
  allocate( atype ( natm ) , itype ( natm ) )
  allocate( list ( natm * 250 ) , point(  natm + 1 ) )
  allocate( xs ( natm ) , ys ( natm ) , zs ( natm ) ) 
  allocate( qia ( natm ) )

  rx = 0.0D0
  ry = 0.0D0
  rz = 0.0D0
  vx = 0.0D0
  vy = 0.0D0
  vz = 0.0D0
  fx = 0.0D0
  fy = 0.0D0
  fz = 0.0D0
  fxs = 0.0d0
  fys = 0.0d0
  fzs = 0.0d0
  xs = 0.0D0
  ys = 0.0D0
  zs = 0.0D0
  rxs = 0.0D0
  rys = 0.0D0
  rzs = 0.0D0
  atype = ''
  list = 0
  point = 0
  qia = 0.0d0

  return 
 
END SUBROUTINE config_alloc


!*********************** SUBROUTINE config_dealloc ****************************
!
! deallocate config quantiteies (see config_alloc)
!
!******************************************************************************

SUBROUTINE config_dealloc

  USE control, ONLY : lvnlist 

  implicit none 

  deallocate( rx  , ry  , rz )
  deallocate( vx  , vy  , vz )
  deallocate( fx  , fy  , fz )
  deallocate( fxs , fys , fzs )
  deallocate( rxs , rys , rzs )
  deallocate( atype )
  deallocate( itype )
  deallocate( list , point )
  deallocate( xs , ys , zs )
  deallocate( qia ) 

  return 

END SUBROUTINE config_dealloc

!*********************** SUBROUTINE center_of_mass ****************************
!  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  ! should depend on mass
!  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! for the moment still here but should be moved somewhere else
! input ax
!******************************************************************************

SUBROUTINE center_of_mass ( ax , ay , az , com )

  implicit none

  ! global
  double precision :: ax ( natm ) , ay ( natm ) , az ( natm )
  double precision :: com ( 3 )

  ! local
  integer :: i


  do i = 1, natm
    com ( 1 ) = com ( 1 )  + ax ( i ) ! * m 
    com ( 2 ) = com ( 2 )  + ay ( i ) ! * m
    com ( 3 ) = com ( 3 )  + az ( i ) ! * m 
  enddo

  com ( 1 )  = com ( 1 ) / dble(natm)
  com ( 2 )  = com ( 2 ) / dble(natm)
  com ( 3 )  = com ( 3 ) / dble(natm)

  return

END SUBROUTINE center_of_mass


!*********************** SUBROUTINE linear_momentum ***************************
!
! Calculate the linear momentum (should be conserved along nve traj)
! NOTE : not used so far
!
!******************************************************************************

SUBROUTINE linear_momentum

  implicit none

  ! local
  integer :: i
  double precision :: Px, Py, Pz, normP

  
  do i = 1,natm
  Px = Px + vx(i)
  Py = Py + vy(i) 
  Pz = Pz + vz(i)
  enddo

  normP = dsqrt(Px*Px + Py*Py + Pz*Pz)

  return  

END SUBROUTINE linear_momentum

!*********************** SUBROUTINE angular_momentum **************************
!
! Calculate the angular momentum (not conserved with pbc)
! NOTE : not used so far
!
!******************************************************************************

SUBROUTINE angular_momentum ( Lx , Ly , Lz , normL )

  implicit none

  ! local
  integer :: i
  double precision :: Lx, Ly, Lz, normL

  Lx = 0.0D0
  Ly = 0.0D0
  Lz = 0.0D0      
  do i = 1,natm
   Lx = Lx + ry(i)*vz(i)- rz(i)*vy(i) 
   Ly = Ly + rz(i)*vx(i)- rx(i)*vz(i) 
   Lz = Lz + rx(i)*vy(i)- ry(i)*vx(i)
  enddo

  normL = dsqrt(Lx*Lx + Ly*Ly + Lz*Lz)

  return

END SUBROUTINE angular_momentum


!TMP
!*********************** SUBROUTINE ions_reference_positions ******************
!
! Calculate the real position of atoms relative to the center of mass (cdm)
! and store them in taui
! cdmi: initial position of the center of mass (cdm) in cartesian coor.  
! NOTE : not used so far
!
!******************************************************************************

SUBROUTINE ions_reference_positions

  implicit none
  double precision :: com ( 3 )
  integer  :: ia

  CALL center_of_mass ( rx , ry , rz , com )

  do ia = 1, natm
    rix ( ia ) = rx ( ia ) - com ( 1 )
    riy ( ia ) = ry ( ia ) - com ( 1 )
    riz ( ia ) = rz ( ia ) - com ( 1 )
  enddo

  return 
 
END SUBROUTINE ions_reference_positions


!*********************** SUBROUTINE ions_displacement *************************
!
! Calculate the sum of the quadratic displacements of the atoms in the ref.
! of cdm respect to the initial positions.
! taui: initial positions in real units in the ref. of cdm
! -------------------------------------------------------------------------
!  att!     tau_ref: starting position in center-of-mass ref. in real units
! -------------------------------------------------------------------------
! NOTE : not used so far
!******************************************************************************

SUBROUTINE ions_displacement( dis, ax , ay , az )

  implicit none
 
  ! global  
  double precision, INTENT(OUT) :: dis(:)
  double precision, intent(in) :: ax(:) , ay(:) , az(:)

  ! local
  double precision :: rdist(3), r2, com(3)
  INTEGER  :: it, ia, isa

 
  ! =========================================================
  !  Compute the current value of cdm "Centro Di Massa"
  ! =========================================================

  CALL center_of_mass ( ax , ay , az , com )
 
  isa = 0
  do it = 1, ntype
    dis(it) = 0.0d0
    r2      = 0.0d0
    do ia = 1 , natmi(it)
      isa = isa + 1
      rdist ( 1 ) = rx (isa) - com ( 1 ) 
      rdist ( 2 ) = ry (isa) - com ( 2 )
      rdist ( 3 ) = rz (isa) - com ( 3 )
      r2 = r2 + ( rdist( 1 ) - rix(isa) )**2  + ( rdist( 2 ) - riy(isa) )**2 + ( rdist( 3 ) - riz(isa) )**2 
    enddo 
    dis(it) = dis(it) + r2 / DBLE(natmi(it))
  enddo
  
  return
 
END SUBROUTINE ions_displacement

END MODULE config
! ===== fmV =====
