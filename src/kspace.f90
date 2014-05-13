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

  USE constants ,       ONLY :  dp 
  USE mpimdff,          ONLY :  decomposition        

  implicit none

  TYPE :: kmesh
    integer                                        :: nk        !< total number of kpoints
    integer                                        :: kmax(3)   !< nb of kpts in each recip. direction in ewald sum.
    real(kind=dp)   , dimension(:)   , allocatable :: kptx      !< kpoints mesh : dimension ( nk )
    real(kind=dp)   , dimension(:)   , allocatable :: kpty      !< kpoints mesh : dimension ( nk )
    real(kind=dp)   , dimension(:)   , allocatable :: kptz      !< kpoints mesh : dimension ( nk )
    real(kind=dp)   , dimension(:)   , allocatable :: kptk      !< k module
    real(kind=dp)   , dimension(:)   , allocatable :: Ak        !< Ak in ewald
    real(kind=dp)   , dimension(:)   , allocatable :: kcoe      !< kcoe in ewald
    complex(kind=dp), dimension(:)   , allocatable :: rhon      !< facteur de structure
    complex(kind=dp), dimension(:,:) , allocatable :: rhon_dk   !< facteur de structure
    complex(kind=dp), dimension(:,:) , allocatable :: expikr    !< facteur de structure
    real(kind=dp)   , dimension(:,:) , allocatable :: ckr, skr  !< facteur de structure
    complex(kind=dp), dimension(:,:) , allocatable :: expikm    !< facteur de structure
    character(len=15)                              :: meshlabel !< giving a name to kmesh
    TYPE ( decomposition )                         :: kpt_dec   
  END TYPE

  TYPE :: kpath
    integer                                        :: nk        !< number of kpoints in the path
    real(kind=dp), dimension(:) , allocatable      :: kptz      !< kpoint path : dimension ( nk )
    real(kind=dp), dimension(:) , allocatable      :: kpty      !< kpoint path : dimension ( nk )
    real(kind=dp), dimension(:) , allocatable      :: kptx      !< kpoint path : dimension ( nk )
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
  USE io,               ONLY :  ionode , stdout
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
  kp%kptx = 0.0d0
  kp%kpty = 0.0d0
  kp%kptz = 0.0d0
  ! first kpoint
  kp%kptx ( 1 ) = ks ( 1 )
  kp%kpty ( 1 ) = ks ( 2 )
  kp%kptz ( 1 ) = ks ( 3 )
  ! last kpoint
  kp%kptx ( nk ) = ke ( 1 )
  kp%kpty ( nk ) = ke ( 2 )
  kp%kptz ( nk ) = ke ( 3 )
  ! path direction
  kp%kdir ( 1 ) = kp%kptx ( nk ) - kp%kptx ( 1 )
  kp%kdir ( 2 ) = kp%kpty ( nk ) - kp%kpty ( 1 )
  kp%kdir ( 3 ) = kp%kptz ( nk ) - kp%kptz ( 1 )
  kp%nk = nk
  kp%pathlabel = path

  if ( ionode ) then
    WRITE ( stdout , '(a,a     )' ) 'path name : ', kp%pathlabel
    WRITE ( stdout , '(a,3f12.5)' ) 'start     : ',   kp%kptx ( 1  ) , kp%kpty ( 1  ),kp%kptz ( 1  )
    WRITE ( stdout , '(a,3f12.5)' ) 'end       : ',   kp%kptx ( nk ) , kp%kpty ( nk ),kp%kptz ( nk )
    WRITE ( stdout , '(a,3f12.5)' ) 'direction : ', ( kp%kdir( i )      , i = 1 , 3 )
 endif

  step = 1.0_dp / REAL ( nk - 1 )

  do ik = 2 , nk
    r_ik = REAL ( ik - 1 , kind = dp )
    kp%kptx(ik) = kp%kptx(1) + r_ik * kp%kdir ( 1 ) * step
    kp%kpty(ik) = kp%kpty(1) + r_ik * kp%kdir ( 2 ) * step
    kp%kptz(ik) = kp%kptz(1) + r_ik * kp%kdir ( 3 ) * step
  enddo

  io_node WRITE ( stdout , '(a,a,a)' ) 'path ',kp%pathlabel,' generated'

  return

END SUBROUTINE get_kpath

! *********************** SUBROUTINE kpoint_sum_init ***************************
!> \brief
!! this subroutine initialized the k vectors for the ewald summation
!> \param[in,out] km kpoint mesh being initialized
!> \note
!! kptx,kpty,kptz(:) are the three components in 2pi/box units
!! kptk is the square module
!> \author
!! FMV
! ******************************************************************************

SUBROUTINE kpoint_sum_init( km , alpha ) 

  USE config,           ONLY :  simu_cell 
  USE io,               ONLY :  ionode , stdout
  USE constants,        ONLY :  tpi

  implicit none
  
  ! global
  TYPE ( kmesh ), intent(inout) :: km

  ! local
  integer :: nx , ny , nz , nk , nymin , nzmin 
  real(kind=dp) :: kx, ky, kz, kk , alpha , alpha2

  io_node WRITE ( stdout      ,'(a,a,a)') 'generate k-points arrays (full) ',km%meshlabel,' mesh'

  alpha2 = alpha * alpha

  nk = 0
  do nx =  0 , km%kmax(1)
!  do nx =  -km%kmax(1) , km%kmax(1)
    if ( nx .eq. 0 ) then
      nymin = 0
    else
      nymin = - km%kmax(2)
    endif
    do ny = nymin  , km%kmax(2)
!    do ny = -km%kmax(2)  , km%kmax(2)
      if ( nx .eq. 0 .and. ny .eq. 0 ) then
        nzmin = 1
      else
        nzmin = - km%kmax(3)
      endif
      do nz = nzmin , km%kmax(3)
!        do nz = -km%kmax(3)  , km%kmax(3)
          nk = nk + 1
          kx = tpi * ( REAL (nx) * simu_cell%B(1,1) +  REAL (ny) * simu_cell%B(1,2) + REAL (nz) * simu_cell%B(1,3) )
          ky = tpi * ( REAL (nx) * simu_cell%B(2,1) +  REAL (ny) * simu_cell%B(2,2) + REAL (nz) * simu_cell%B(2,3) )
          kz = tpi * ( REAL (nx) * simu_cell%B(3,1) +  REAL (ny) * simu_cell%B(3,2) + REAL (nz) * simu_cell%B(3,3) )
          kk = kx * kx + ky * ky + kz * kz
          km%kptx ( nk ) = kx
          km%kpty ( nk ) = ky
          km%kptz ( nk ) = kz
          km%kptk ( nk ) = kk
          km%Ak   ( nk ) = EXP ( - kk * 0.25_dp / alpha2 ) / kk 
          km%kcoe ( nk ) = 2.0_dp * ( 1.0_dp / kk + 1.0_dp / alpha2 / 4.0_dp )
      enddo
    enddo
  enddo
 
  if ( nk .ne. km%nk ) then
    io_node WRITE ( stdout ,'(a,2i7,x,a)') 'number of k-points do not match in kpoint_sum_init', nk , km%nk , km%meshlabel
  endif

  km%nk = nk 

  ! ======================
  !  organized kpt arrays 
  ! ======================
  call reorder_kpt ( km ) 
  io_node WRITE ( stdout      ,'(a,i)') '(full) kpt arrays sorted',nk

  return

END SUBROUTINE kpoint_sum_init

! *********************** SUBROUTINE kpoint_sum_init_BZ **********************
!> \brief
!! this subroutine initialized the k vectors of the Brillouin zone 
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
SUBROUTINE kpoint_sum_init_BZ ( km , alpha )

  USE config,           ONLY :  simu_cell 
  USE io,               ONLY :  ionode , stdout
  USE constants,        ONLY :  tpi

  implicit none
  
  ! global
  TYPE ( kmesh ), intent(inout) :: km
  real(kind=dp) :: alpha
  
  ! local
  integer       :: nx , ny , nz , nk 
  real(kind=dp) :: kx, ky, kz, kk
  real(kind=dp) :: kxi , kyi , kzi , alpha2

  alpha2 = alpha * alpha

  io_node WRITE ( stdout      ,'(a,a,a)') 'generate k-points arrays (half) ',km%meshlabel,' mesh'
  
  nk = 0
  do nx =  0 , km%kmax(1) 
    kxi = REAL( nx , kind = dp ) / REAL( km%kmax(1) , kind = dp )
    do ny = 0  , km%kmax(2)
      kyi = REAL( ny , kind = dp ) / REAL( km%kmax(2) , kind = dp )
      do nz = 0 , km%kmax(3)
        kzi = REAL( nz , kind = dp ) / REAL( km%kmax(3) , kind = dp )
        nk = nk + 1
        kx = tpi * ( kxi * simu_cell%B(1,1) +  kyi * simu_cell%B(1,2) + kzi * simu_cell%B(1,3) )
        ky = tpi * ( kxi * simu_cell%B(2,1) +  kyi * simu_cell%B(2,2) + kzi * simu_cell%B(2,3) )
        kz = tpi * ( kxi * simu_cell%B(3,1) +  kyi * simu_cell%B(3,2) + kzi * simu_cell%B(3,3) )
        kk = kx * kx + ky * ky + kz * kz
        km%kptx(nk) = kx 
        km%kpty(nk) = ky 
        km%kptz(nk) = kz 
        km%kptk(nk) = kk 
      enddo
    enddo
  enddo

  if ( nk .ne. km%nk ) then
    io_node WRITE ( stdout ,'(a,2i7,x,a)') 'number of k-points do not match in kpoint_sum_init_BZ ', nk , km%nk , km%meshlabel
    STOP
  endif
  io_node WRITE ( stdout ,'(a,a,i8)') 'number of k-points in ', km%meshlabel , km%nk

  ! ======================
  !  organized kpt arrays 
  ! ======================
  call reorder_kpt ( km )
  io_node WRITE ( stdout      ,'(a)')     '(half) kpt arrays sorted'

  return

END SUBROUTINE kpoint_sum_init_BZ

! *********************** SUBROUTINE reorder_kpt *******************************
!> \brief
!! this subroutine reorder kpt arrays (increasing k^2)
!> \param[in,out] km kpoint mesh being initialized
!> \author
!! FMV
! ******************************************************************************
SUBROUTINE reorder_kpt ( km )

  USE io,  ONLY :  ionode , stdout 

  implicit none

  !global
  TYPE ( kmesh ), intent(inout) :: km

  !local
  integer :: ik, lk
  real(kind=dp), dimension (:), allocatable :: tkpt
  real(kind=dp), dimension (:), allocatable :: tmpkx , tmpky , tmpkz , tmpak , tmpkcoe
  integer, dimension (:), allocatable :: labelkpt, labelt

  allocate ( tkpt ((km%nk+1)/2) , labelkpt(km%nk), labelt((km%nk+1)/2) )
  allocate ( tmpkx(km%nk) , tmpky(km%nk) , tmpkz(km%nk)  , tmpak(km%nk) , tmpkcoe (km%nk) )

  ! ==============================
  !  set the initial array labels
  ! ==============================
  do ik=1,km%nk
    labelkpt(ik)=ik
  enddo

  ! ===========================================
  !  arrays are sorted for increasing k^2 
  !  (see tools.f90 for more details )
  !  the old labels are stored in labelkpt
  !  that will be used to reorganized kx,ky,kz,Ak,kcoe
  ! ===========================================
  call merge_sort ( km%kptk , km%nk , tkpt , labelkpt , labelt )

  ! ==============================
  !  save previous order
  ! ==============================
  tmpkx=km%kptx(:)
  tmpky=km%kpty(:)
  tmpkz=km%kptz(:) 
  if ( ALLOCATED(km%Ak)   ) tmpak=km%Ak
  if ( ALLOCATED(km%kcoe) ) tmpkcoe=km%kcoe
  ! ===============================================
  !  change kptx , kpty , kptz following kptk sort
  ! ===============================================
  do ik=1,km%nk
    lk=labelkpt(ik)
    km%kptx(ik) = tmpkx(lk) 
    km%kpty(ik) = tmpky(lk) 
    km%kptz(ik) = tmpkz(lk) 
    if ( ALLOCATED(km%Ak)   ) km%Ak(ik)   = tmpak(lk) 
    if ( ALLOCATED(km%kcoe) ) km%kcoe(ik) = tmpkcoe(lk) 
  enddo

  deallocate ( tkpt , labelkpt , labelt )
  deallocate ( tmpkx , tmpky , tmpkz, tmpak , tmpkcoe    )

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
! deprecated not working anymore
!SUBROUTINE struc_fact ( km ) 
!  
!  USE config,           ONLY :  natm , itype , ntype , rx , ry , rz 
!  USE io,               ONLY :  ionode , kunit_STRFACFF
!  USE constants,        ONLY :  imag , mimag

!  implicit none

!  ! global
!  TYPE ( kmesh ), intent(inout) :: km

!  ! local
!  !integer :: it
!  integer :: ia, ik
!  real(kind=dp) :: arg , rxi , ryi , rzi 

!  !  exp ( i k . r ) 
!  km%strf(:,:) = (0.0_dp,0.0_dp)
!  do it = 1, ntype
!     do ia = 1, natm
!        rxi = rx ( ia ) 
!        ryi = ry ( ia ) 
!        rzi = rz ( ia ) 
 
!        if ( itype (ia) .eq. it ) then
!           do ik = 1, km%nk 
!              arg = ( km%kpt ( ik , 1) * rxi + km%kpt ( ik , 2 ) * ryi + km%kpt ( ik , 3) * rzi ) 
!              km%strf  ( ik , ia ) = km%strf  ( ik , ia ) + EXP( imag * arg ) !!! check the order of index
!           enddo
!        endif
!     enddo
!  enddo

!  return

!END SUBROUTINE struc_fact

! similar as in kspace.f90 but with type conditions for efg calculation
SUBROUTINE charge_density_k ( km , mu , ldip , update_mu )

  USE constants,        ONLY :  imag
  USE config,           ONLY :  natm , rx , ry , rz , qia , itype

  implicit none

  ! global
  TYPE ( kmesh ), intent(inout) :: km
  real(kind=dp) , intent(in)    :: mu    ( natm , 3 )
  logical                       :: update_mu
  logical                       :: ldip

  ! local
  integer :: ia ,ik
  real(kind=dp) :: rxi , ryi , rzi , k_dot_r , k_dot_mu , mux , muy , muz, qi
  real(kind=dp) :: kx , ky , kz
  complex(kind=dp) :: expikr , expikm , sumia
        
  km%rhon  = (0.0_dp,0.0_dp)
  do ik = km%kpt_dec%istart, km%kpt_dec%iend
    kx = km%kptx ( ik )
    ky = km%kpty ( ik )
    kz = km%kptz ( ik )

    sumia = (0.0_dp,0.0_dp)
    do ia = 1 , natm
      qi  = qia (ia )
      rxi = rx ( ia )
      ryi = ry ( ia )
      rzi = rz ( ia )
      k_dot_r   = ( kx * rxi + ky * ryi + kz * rzi )
      expikr    = EXP ( imag * k_dot_r )
      km%ckr(ia,ik) = COS( k_dot_r )
      km%skr(ia,ik) = SIN( k_dot_r )
      sumia = sumia + qi * expikr 
      if ( .not. ldip ) cycle
      mux = mu ( ia , 1 )
      muy = mu ( ia , 2 )
      muz = mu ( ia , 3 )
      k_dot_mu  = ( mux * kx + muy * ky + muz * kz )
      expikm    = k_dot_mu * expikr
      !km%expikm(ia,ik) = expikm
      sumia = sumia + imag * expikm 
    enddo
    km%rhon(ik) = sumia
  enddo

  return

END SUBROUTINE charge_density_k


END MODULE kspace 
! ===== fmV =====


