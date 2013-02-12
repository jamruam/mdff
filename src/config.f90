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

!*********************** MODULE CONF ******************************************
!
! This module should deal with everything related to the current configuration
! positions, velocities, forces, density , box .... 
!
!******************************************************************************

MODULE config

  USE cell 

  implicit none

  integer, PARAMETER :: natmmax  = 16   ! tmp not controling all the arrays
  integer, PARAMETER :: ntypemax = 16
  integer, PARAMETER :: npolmax  = 16

  character*60, SAVE :: system                                        ! system name                                              

  logical :: lfcc   
  logical :: lsc   ! not tested
  logical :: lbcc  ! not yet

  integer :: natm                                                     ! nb of atoms
  integer :: ntype                                                    ! number of types
  integer :: ncell                                                    ! number of cell (with lfcc)
  integer, dimension(:), allocatable :: itype                         ! type of atome i array 
  integer, dimension(:), allocatable :: ipolar                        ! .eq. 1 if polar 
  integer, dimension(:), allocatable :: list, point                   ! vnlist info
  integer, dimension(0:ntypemax)     :: natmi                         ! number of atoms (per type)

  TYPE ( celltype ) :: simu_cell
  double precision :: tau_nonb ( 3 , 3 )                              ! stress tensor ( lennard-jones , morse ... )
  double precision :: tau_coul ( 3 , 3 )                              ! stress tensor coulombic
  double precision :: rho                                             ! density  

  double precision , dimension(:)    , allocatable :: rx  , ry  , rz  ! positions
  double precision , dimension(:)    , allocatable :: vx  , vy  , vz  ! velocities
  double precision , dimension(:)    , allocatable :: fx  , fy  , fz  ! forces
  double precision , dimension(:)    , allocatable :: fxs , fys , fzs ! forces (previous step t-dt) beeman
  double precision , dimension(:)    , allocatable :: rxs , rys , rzs ! previous positions for leap-frog integrator
  double precision , dimension(:)    , allocatable :: xs  , ys  , zs  ! last positions in verlet list
  double precision , dimension(:)    , allocatable :: rix , riy , riz ! positions in the center of mass reference 

  double precision , dimension(:)    , allocatable :: qia             ! charge on ion 
  double precision , dimension(:,:)  , allocatable :: dipia           ! dipole on ion 
  double precision , dimension(:,:)  , allocatable :: dipia_ind       ! induced dipole on ion 
  double precision , dimension(:,:)  , allocatable :: dipia_wfc       ! induced dipole on ion from Wannier centers
  double precision , dimension(:,:,:), allocatable :: polia           ! polarisation on ion

  double precision , dimension(:), allocatable     :: phi_coul_tot    ! coulombic potential 

  character*3, dimension(:), allocatable           :: atype           ! atom type A or B 
  character*3, dimension(0:ntypemax)               :: atypei          ! type of atoms (per type)

  character*60 :: struct                                              ! crystal structure for fcc
  character*60 :: struct_allowed(2)                     
  data struct_allowed / 'NaCl' , 'random' /


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

  ! ===========================================================
  ! read initial configuration only for calc = md
  ! for opt, vib and efg configuations are read in other routines 
  ! ===========================================================
  if ( calc .eq. 'md' ) then       
    CALL read_pos
  endif

#ifdef debug
  CALL print_config_sample(0,0)
#endif

  return
 
END SUBROUTINE config_init 


!*********************** SUBROUTINE config_print_info *************************
!
! print information to standard output about the starting configuration
!
!******************************************************************************

SUBROUTINE config_print_info(kunit)

  USE io_file,  ONLY :  ionode 

  implicit none

  ! global 
  integer , intent ( in ) :: kunit 

  ! local
  integer :: it , i

  if ( ionode ) then
    WRITE ( kunit ,'(a)')            ''
    WRITE ( kunit ,'(a,a)')          'system                             : ',system
    WRITE ( kunit ,'(a,i16)')        'natm                               = ',natm
    WRITE ( kunit ,'(a,i16)')        'ntype                              = ',ntype
    do it = 1 , ntype     
      WRITE ( kunit ,'(a,a,a,i16,f8.2,a1)') &
                          'n',atypei(it),'                               = ',natmi(it),DBLE(natmi(it))/DBLE(natm) * 100.0d0,'%'
    enddo
    WRITE ( kunit ,'(a)')            ''
    WRITE ( kunit ,'(a)')            '---------------------------------------------------------------------------------------------'
    WRITE ( kunit ,'(a)')            ''
    WRITE ( kunit ,'(a)')            'direct     basis : '
    WRITE ( kunit ,'(a,3f16.4)')     'a_vector                           = ',simu_cell%A(1,1),simu_cell%A(2,1),simu_cell%A(3,1) 
    WRITE ( kunit ,'(a,3f16.4)')     'b_vector                           = ',simu_cell%A(1,2),simu_cell%A(2,2),simu_cell%A(3,2) 
    WRITE ( kunit ,'(a,3f16.4)')     'c_vector                           = ',simu_cell%A(1,3),simu_cell%A(2,3),simu_cell%A(3,3) 
    WRITE ( kunit ,'(a)')            ''
    WRITE ( kunit ,'(a,3f16.4)')     'cell parameters     (direct)       = ',(simu_cell%ANORM(i),i=1,3)
    WRITE ( kunit ,'(a,3f16.4)')     'perpendicular width (direct)       = ',simu_cell%WA,simu_cell%WB,simu_cell%WC
    WRITE ( kunit ,'(a,3f16.4)')     'angles              (direct)       = ',simu_cell%ALPH,simu_cell%BET,simu_cell%GAMM
    WRITE ( kunit ,'(a,f16.4)')      'volume              (direct)       = ',simu_cell%omega
    WRITE ( kunit ,'(a)')            ''
    WRITE ( kunit ,'(a)')            '---------------------------------------------------------------------------------------------'
    WRITE ( kunit ,'(a)')            ''
    WRITE ( kunit ,'(a)')            'reciprocal basis : '
    WRITE ( kunit ,'(a,3f16.4)')     'a*_vector                          = ',simu_cell%B(1,1),simu_cell%B(2,1),simu_cell%B(3,1) 
    WRITE ( kunit ,'(a,3f16.4)')     'b*_vector                          = ',simu_cell%B(1,2),simu_cell%B(2,2),simu_cell%B(3,2) 
    WRITE ( kunit ,'(a,3f16.4)')     'c*_vector                          = ',simu_cell%B(1,3),simu_cell%B(2,3),simu_cell%B(3,3) 
    WRITE ( kunit ,'(a)')            ''
    WRITE ( kunit ,'(a,3f16.4)')     'cell parameters     (reciprocal)   = ',(simu_cell%BNORM(i),i=1,3)
    WRITE ( kunit ,'(a,3f16.4)')     'perpendicular width (reciprocal)   = ',simu_cell%RWA,simu_cell%RWB,simu_cell%RWC
    WRITE ( kunit ,'(a,3f16.4)')     'angles              (reciprocal)   = ',simu_cell%RALPH,simu_cell%RBET,simu_cell%RGAMM
    WRITE ( kunit ,'(a,f16.4)')      'volume              (reciprocal)   = ',simu_cell%romega
    WRITE ( kunit ,'(a)')            ''
    WRITE ( kunit ,'(a)')            '---------------------------------------------------------------------------------------------'
    WRITE ( kunit ,'(a)')            ''
    WRITE ( kunit ,'(a,f16.4)')      'density                              = ',rho
    WRITE ( kunit ,'(a)')            ''
    
  endif 

  return
  
END SUBROUTINE config_print_info

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

  !CALL  periodicbc ( natm , xxx , yyy , zzz )
  
  if ( ionode ) then
  OPEN ( kunit_CONTFF ,file = 'CONTFF',STATUS = 'UNKNOWN')
      WRITE ( kunit_CONTFF,'(i)') natm 
      WRITE ( kunit_CONTFF,'(a)') system
      WRITE ( kunit_CONTFF,'(3f20.12)') simu_cell%A ( 1 , 1 ) , simu_cell%A ( 2 , 1 ) , simu_cell%A ( 3 , 1 )
      WRITE ( kunit_CONTFF,'(3f20.12)') simu_cell%A ( 1 , 2 ) , simu_cell%A ( 2 , 2 ) , simu_cell%A ( 3 , 2 )
      WRITE ( kunit_CONTFF,'(3f20.12)') simu_cell%A ( 1 , 3 ) , simu_cell%A ( 2 , 3 ) , simu_cell%A ( 3 , 3 )
      WRITE ( kunit_CONTFF,'(i4)') ntype 
      WRITE ( kunit_CONTFF,*) ( atypei(it) , it=1,ntype ) 
      WRITE ( kunit_CONTFF,*) ( natmi (it) , it=1,ntype ) 
      WRITE ( kunit_CONTFF,'(a,9e20.12)') ( atype ( ia ) , xxx ( ia ) , yyy ( ia ) , zzz ( ia ) , & 
                                                           vx  ( ia ) , vy  ( ia ) , vz  ( ia ) , &
                                                           fx  ( ia ) , fy  ( ia ) , fz ( ia )  , ia = 1 , natm )
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
! forces on atom ia previous step (beeman)  : fxs(ia) , fys(ia)  , fzs(ia) 
! previous pos. leap-frog algo.             : rxs(ia) , rys(ia)  , rzs(ia) 
! last config. list verlet update           : xs(ia)  , ys(ia)   , zs(ia) 
! charge on atom                            : qia 
! static dipole on atom                     : dipia
! induced dipole on atom                    : dipia_ind
! induced dipole on atom from wannier c.    : dipia_wfc
! polarisation tensor on atom               : polia

! atom type information (still not nice) 
!
! type of atom ia (character)               : atype(ia)  
! type of atom ia (integer)                 : itype(ia) = 1 .. ntype
! number of atoms of type it                : natmi(it) .le. natm
! name of type of type it                   : atypei(it) => atypei(0) = 'ALL'
! verlet list                               : list , point 
! is the atom polarized ?                   : ipolar = 1 if yes 
! coulombic potential                       : phi_coul_tot
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
  allocate( list ( natm * 1000 ) , point (  natm + 1 ) )
  allocate( xs ( natm ) , ys ( natm ) , zs ( natm ) ) 
  allocate( qia ( natm ) )
  allocate( dipia ( natm , 3 ) )
  allocate( dipia_ind ( natm , 3 ) )
  allocate( dipia_wfc ( natm , 3 ) )
  allocate( polia ( natm , 3 , 3 ) )
  allocate( ipolar ( natm ) )
  allocate( phi_coul_tot ( natm ) ) ! only if we calculated coulombic interactions

  rx    = 0.0D0
  ry    = 0.0D0
  rz    = 0.0D0
  vx    = 0.0D0
  vy    = 0.0D0
  vz    = 0.0D0
  fx    = 0.0D0
  fy    = 0.0D0
  fz    = 0.0D0
  fxs   = 0.0d0
  fys   = 0.0d0
  fzs   = 0.0d0
  xs    = 0.0D0
  ys    = 0.0D0
  zs    = 0.0D0
  rxs   = 0.0D0
  rys   = 0.0D0
  rzs   = 0.0D0
  atype = ''
  list      = 0
  point     = 0
  qia       = 0.0d0
  dipia     = 0.0d0
  dipia_ind = 0.0d0
  dipia_wfc = 0.0d0
  polia     = 0.0d0
  ipolar    = 0
  phi_coul_tot = 0.0d0

  return 
 
END SUBROUTINE config_alloc


!*********************** SUBROUTINE config_dealloc ****************************
!
! deallocate config quantities (see config_alloc)
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
  deallocate( dipia ) 
  deallocate( dipia_ind ) 
  deallocate( dipia_wfc ) 
  deallocate( polia ) 
  deallocate( ipolar ) 
  deallocate( phi_coul_tot ) ! well only if we calculated coulombic interactions

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
  double precision , intent ( in  ) :: ax ( natm ) , ay ( natm ) , az ( natm )
  double precision , intent ( out ) :: com ( 0 : ntypemax , 3 )

  ! local
  integer :: ia , it 

  com = 0.0D0

  do ia = 1 , natm
    it = itype ( ia )    
    com ( it , 1 ) = com ( it , 1 )  + ax ( ia ) ! * m 
    com ( it , 2 ) = com ( it , 2 )  + ay ( ia ) ! * m
    com ( it , 3 ) = com ( it , 3 )  + az ( ia ) ! * m 
    com ( 0  , 1 ) = com ( 0  , 1 )  + ax ( ia ) ! * m 
    com ( 0  , 2 ) = com ( 0  , 2 )  + ay ( ia ) ! * m
    com ( 0  , 3 ) = com ( 0  , 3 )  + az ( ia ) ! * m
  enddo

  do it = 0 , ntype
    com ( it , 1 )  = com ( it , 1 ) / DBLE ( natmi ( it ) )
    com ( it , 2 )  = com ( it , 2 ) / DBLE ( natmi ( it ) )
    com ( it , 3 )  = com ( it , 3 ) / DBLE ( natmi ( it ) )
  enddo

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
  integer :: ia
  double precision :: Px, Py, Pz, normP

  do ia = 1 , natm
    Px = Px + vx ( ia )
    Py = Py + vy ( ia ) 
    Pz = Pz + vz ( ia )
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
  integer :: ia 
  double precision :: Lx, Ly, Lz, normL

  Lx = 0.0D0
  Ly = 0.0D0
  Lz = 0.0D0      
  do ia = 1 , natm
   Lx = Lx + ry ( ia ) * vz ( ia ) - rz ( ia ) * vy ( ia ) 
   Ly = Ly + rz ( ia ) * vx ( ia ) - rx ( ia ) * vz ( ia ) 
   Lz = Lz + rx ( ia ) * vy ( ia ) - ry ( ia ) * vx ( ia )
  enddo

  normL = dsqrt( Lx * Lx + Ly * Ly + Lz * Lz)

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
  double precision :: com ( 0:ntypemax, 3 )
  integer  :: ia

  CALL center_of_mass ( rx , ry , rz , com )

  do ia = 1 , natm
    rix ( ia ) = rx ( ia ) - com ( 0 , 1 )
    riy ( ia ) = ry ( ia ) - com ( 0 , 1 )
    riz ( ia ) = rz ( ia ) - com ( 0 , 1 )
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
  double precision, intent ( out ) :: dis(:)
  double precision, intent ( in  ) :: ax (:) , ay(:) , az(:)

  ! local
  double precision :: rdist(3), r2, com(0:ntypemax,3)
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
      rdist ( 1 ) = rx (isa) - com ( 0 , 1 ) 
      rdist ( 2 ) = ry (isa) - com ( 0 , 2 )
      rdist ( 3 ) = rz (isa) - com ( 0 , 3 )
      r2 = r2 + ( rdist( 1 ) - rix(isa) )**2 + &
                ( rdist( 2 ) - riy(isa) )**2 + &
                ( rdist( 3 ) - riz(isa) )**2 
    enddo 
    dis(it) = dis(it) + r2 / DBLE(natmi(it))
  enddo
  
  return
 
END SUBROUTINE ions_displacement

END MODULE config
! ===== fmV =====
