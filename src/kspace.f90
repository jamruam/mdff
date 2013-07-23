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
!#define debug
! ======= Hardware =======

! *********************** MODULE rspace ****************************************
!> @brief
!> module related to reciprocal space summation 
!> @author
!> FMV
! ******************************************************************************
MODULE kspace 

  USE constants , ONLY : dp 

  implicit none

  TYPE :: kmesh
    integer                                        :: nk        !< total number of kpoints
    integer                                        :: nkcut     !< (internal) number of k-points in each direction 
    integer                                        :: kmax(3)   !< nb of kpts in each recip. direction in ewald sum.
    real(kind=dp)   , dimension(:,:) , allocatable :: kpt       !< kpoints mesh : dimension ( 3 , nk )
    real(kind=dp)   , dimension(:)   , allocatable :: kptk      !< k module
    complex(kind=dp), dimension(:,:) , allocatable :: strf      !< facteur de structure
    character(len=15)                              :: meshlabel !< giving a name to kmesh
  END TYPE

  TYPE :: kpath
    integer                                        :: nk        !< number of kpoints in the path
    real(kind=dp), dimension(:,:) , allocatable    :: kpt       !< kpoint path : dimension ( 3 , nk )
    real(kind=dp)                                  :: kdir(3)   !< path vector  
    character(len=3)                               :: pathlabel !< giving a name to kpath
  END TYPE


CONTAINS

! *********************** SUBROUTINE get_path **********************************
!> \brief
!! this subroutine initialized kpoint path for band structure calculation
!> \author
!! FMV
! ******************************************************************************
SUBROUTINE get_kpath ( ks , ke , nk , path , kp )

  USE config,           ONLY :  simu_cell
  USE io_file,          ONLY :  ionode , stdout
  USE cell,             ONLY :  dirkar

  implicit none

  ! global 
  real(kind=dp)    :: ks(3) , ke(3)
  integer          :: nk 
  TYPE ( kpath )   :: kp
  character(len=3) :: path
  ! local  
  integer        :: ik , i
  real(kind=dp)  :: step , r_ik
  
  ! initialize
  kp%kpt = 0.0d0
  ! first kpoint
  kp%kpt ( 1 , 1  ) = ks ( 1 )
  kp%kpt ( 2 , 1  ) = ks ( 2 )
  kp%kpt ( 3 , 1  ) = ks ( 3 )
  ! last kpoint
  kp%kpt ( 1 , nk ) = ke ( 1 )
  kp%kpt ( 2 , nk ) = ke ( 2 )
  kp%kpt ( 3 , nk ) = ke ( 3 )
  ! path direction
  kp%kdir ( 1 ) = kp%kpt ( 1 , nk ) - kp%kpt ( 1 , 1 )
  kp%kdir ( 2 ) = kp%kpt ( 2 , nk ) - kp%kpt ( 2 , 1 )
  kp%kdir ( 3 ) = kp%kpt ( 3 , nk ) - kp%kpt ( 3 , 1 )
  kp%nk = nk
  kp%pathlabel = path

  if ( ionode ) then
    WRITE ( stdout , '(a,a     )' ) 'path name : ', kp%pathlabel
    WRITE ( stdout , '(a,3f12.5)' ) 'start     : ', ( kp%kpt ( i , 1  ) , i = 1 , 3 )
    WRITE ( stdout , '(a,3f12.5)' ) 'end       : ', ( kp%kpt ( i , nk ) , i = 1 , 3 )
    WRITE ( stdout , '(a,3f12.5)' ) 'direction : ', ( kp%kdir( i )      , i = 1 , 3 )
 endif

  step = 1.0_dp / REAL ( nk - 1 )

  do ik = 2 , nk
    r_ik = REAL ( ik - 1 , kind = dp )
    kp%kpt(1,ik) = kp%kpt(1,1) + r_ik * kp%kdir ( 1 ) * step
    kp%kpt(2,ik) = kp%kpt(2,1) + r_ik * kp%kdir ( 2 ) * step
    kp%kpt(3,ik) = kp%kpt(3,1) + r_ik * kp%kdir ( 3 ) * step
  enddo

  io_node WRITE ( stdout , '(a,a,a)' ) 'path ',kp%pathlabel,' generated'

  return

END SUBROUTINE get_kpath

! *********************** SUBROUTINE kpoint_sum_init ***************************
!> \brief
!! this subroutine initialized the k vectors for the ewald summation
!> \param[in,out] km kpoint mesh being initialized
!> \note
!! kpt(3,:) are the three components in 2pi/box units
!! kptk is the square module
!> \author
!! FMV
! ******************************************************************************

SUBROUTINE kpoint_sum_init( km ) 

  USE config,           ONLY :  simu_cell 
  USE io_file,          ONLY :  ionode , stdout
  USE constants,        ONLY :  tpi

  implicit none
  
  ! global
  TYPE ( kmesh ), intent(inout) :: km

  ! local
  integer :: nx , ny , nz , nk 
  real(kind=dp) :: kx, ky, kz, kk

  io_node WRITE ( stdout      ,'(a,a,a)') 'generate k-points arrays (full) ',km%meshlabel,' mesh'

  nk = 0
  do nx =  - km%kmax(1) , km%kmax(1)
    do ny = - km%kmax(2)  , km%kmax(2)
      do nz = - km%kmax(3) , km%kmax(3)
        if ( ( nx .ne. 0 ) .or. ( ny .ne. 0) .or. ( nz .ne. 0) ) then
          nk = nk + 1
          kx = tpi * ( DBLE (nx) * simu_cell%B(1,1) +  DBLE (ny) * simu_cell%B(1,2) + DBLE (nz) * simu_cell%B(1,3) )
          ky = tpi * ( DBLE (nx) * simu_cell%B(2,1) +  DBLE (ny) * simu_cell%B(2,2) + DBLE (nz) * simu_cell%B(2,3) )
          kz = tpi * ( DBLE (nx) * simu_cell%B(3,1) +  DBLE (ny) * simu_cell%B(3,2) + DBLE (nz) * simu_cell%B(3,3) )
          kk = kx * kx + ky * ky + kz * kz
          km%kpt ( 1 , nk ) = kx
          km%kpt ( 2 , nk ) = ky
          km%kpt ( 3 , nk ) = kz
          km%kptk ( nk )    = kk
        endif
      enddo
    enddo
  enddo
 
  if ( nk .ne. km%nkcut ) then
    io_node WRITE ( stdout ,'(a,2i7,x,a)') 'number of k-points do not match in kpoint_sum_init', nk , km%nkcut , km%meshlabel
    STOP
  endif

  km%nk = nk 

  ! ======================
  !  organized kpt arrays 
  ! ======================
  call reorder_kpt ( km ) 
  io_node WRITE ( stdout      ,'(a)') '(full) kpt arrays sorted'

  return

END SUBROUTINE kpoint_sum_init

! *********************** SUBROUTINE kpoint_sum_init_half **********************
!> \brief
!! this subroutine initialized the k vectors for the ewald summation
! kpt(3,:) are the three components in 2pi/box units
! kptk is the square module
!> \param[in,out] km kpoint mesh being initialized
!> \note
!! kpt(3,:) are the three components in 2pi/box units
!! kptk is the square module
!> \author
!! FMV
!> \note
!! this subroutine generate only half of the brillouin zone 
! ******************************************************************************
SUBROUTINE kpoint_sum_init_half ( km )

  USE config,           ONLY :  simu_cell 
  USE io_file,          ONLY :  ionode , stdout
  USE constants,        ONLY :  tpi

  implicit none
  
  ! global
  TYPE ( kmesh ), intent(inout) :: km
  
  ! local
  integer :: nx , ny , nz , nk 
  real(kind=dp) :: kx, ky, kz, kk

  io_node WRITE ( stdout      ,'(a,a,a)') 'generate k-points arrays (half) ',km%meshlabel,' mesh'
  
  nk = 0
  do nx =  0 , km%kmax(1) 
    do ny = 0  , km%kmax(2)
      do nz = 0 , km%kmax(3)
        if ( ( nx .ne. 0 ) .or. ( ny .ne. 0) .or. ( nz .ne. 0) ) then
          nk = nk + 1
          kx = tpi * DBLE (nx) * simu_cell%BNORM(1) 
          ky = tpi * DBLE (ny) * simu_cell%BNORM(2)
          kz = tpi * DBLE (nz) * simu_cell%BNORM(3)
          kk = kx * kx + ky * ky + kz * kz
          km%kpt(1,nk) = kx
          km%kpt(2,nk) = ky
          km%kpt(3,nk) = kz
          km%kptk(nk) = kk
        endif
      enddo
    enddo
  enddo

  if ( nk .ne. km%nkcut ) then
    io_node WRITE ( stdout ,'(a,2i7,x,a)') 'number of k-points do not match in kpoint_sum_init_half ', nk , km%nkcut , km%meshlabel
    STOP
  endif

  ! ======================
  !  organized kpt arrays 
  ! ======================
  call reorder_kpt ( km )
  io_node WRITE ( stdout      ,'(a)')     '(half) kpt arrays sorted'

  return

END SUBROUTINE kpoint_sum_init_half

! *********************** SUBROUTINE reorder_kpt *******************************
!> \brief
!! this subroutine reorder kpt arrays (increasing k^2)
!> \param[in,out] km kpoint mesh being initialized
!> \author
!! FMV
! ******************************************************************************
SUBROUTINE reorder_kpt ( km )

  USE io_file,  ONLY :  ionode , stdout 

  implicit none

  !global
  TYPE ( kmesh ), intent(inout) :: km

  !local
  integer :: ik, lk
  real(kind=dp), dimension (:), allocatable :: tkpt
  real(kind=dp), dimension (:), allocatable :: tmpkx , tmpky , tmpkz
  integer, dimension (:), allocatable :: labelkpt, labelt

  allocate ( tkpt ((km%nkcut+1)/2) , labelkpt(km%nkcut), labelt((km%nkcut+1)/2) )
  allocate ( tmpkx(km%nkcut) , tmpky(km%nkcut) , tmpkz(km%nkcut)  )

  ! ==============================
  !  set the initial array labels
  ! ==============================
  do ik=1,km%nkcut
    labelkpt(ik)=ik
  enddo

  ! ===========================================
  !  arrays are sorted for increasing k^2 
  !  (see tools.f90 for more details )
  !  the old labels are stored in labelkpt
  !  that will be used to reorganized kx,ky,kz
  ! ===========================================
  call merge_sort ( km%kptk , km%nkcut , tkpt , labelkpt , labelt )

  ! ==============================
  !  save previous order
  ! ==============================
  tmpkx=km%kpt(1,:)
  tmpky=km%kpt(2,:)
  tmpkz=km%kpt(3,:) 
  ! ===============================================
  !  change kptx , kpty , kptz following kptk sort
  ! ===============================================
  do ik=1,km%nkcut
    lk=labelkpt(ik)
    km%kpt(1,ik) = tmpkx(lk) 
    km%kpt(2,ik) = tmpky(lk) 
    km%kpt(3,ik) = tmpkz(lk) 
  enddo

  deallocate ( tkpt , labelkpt , labelt )
  deallocate ( tmpkx , tmpky , tmpkz    )

  return

END SUBROUTINE reorder_kpt

! *********************** SUBROUTINE struct_fact *******************************
!!> \brief
!! calculate the structure factors for each type of atoms in the unit cell
!! Basically used to calculate k-space charge density. 
!> \param[in,out] km kpoint mesh being initialized
!> \author
!! FMV
! ******************************************************************************
SUBROUTINE struc_fact ( km ) 
  
  USE config,           ONLY :  natm , itype , ntype , rx , ry , rz 
  USE io_file,          ONLY :  ionode , kunit_STRFACFF
  USE constants,        ONLY :  imag , mimag

  implicit none

  ! global
  TYPE ( kmesh ), intent(inout) :: km

  ! local
  !integer :: it
  integer :: ia, ik
  real(kind=dp) :: arg , rxi , ryi , rzi 

  !  exp ( i k . r ) 
  km%strf(:,:) = (0.0_dp,0.0_dp)
!  do it = 1, ntype
     do ia = 1, natm
        rxi = rx ( ia ) 
        ryi = ry ( ia ) 
        rzi = rz ( ia ) 
 
!        if ( itype (ia) .eq. it ) then
           do ik = 1, km%nkcut 
              arg = ( km%kpt ( 1 , ik ) * rxi + km%kpt ( 2 , ik ) * ryi + km%kpt ( 3 , ik ) * rzi ) 
              km%strf  ( ik , ia ) = km%strf  ( ik , ia ) + EXP( imag * arg ) 
           enddo
!        endif
     enddo
!  enddo

  return

END SUBROUTINE struc_fact

END MODULE kspace 
! ===== fmV =====


