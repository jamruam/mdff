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

! module related to real space summation 
MODULE rspace

  ! ===================
  !  direct summation
  ! ===================

  TYPE rmesh
    integer                                         :: ncell              ! numberof cell in the direct calculation in each direction ncell x ncell x ncell 
    integer                                         :: ncmax              ! (internal) number of cell in direct summation
    double precision, dimension(:,:)  , allocatable :: boxxyz
    integer         , dimension(:)    , allocatable :: lcell
  END TYPE rmesh



CONTAINS

!*********************** SUBROUTINE direct_sum_init ***************************
!
! generate vectors in real space of the neigboring cells in a cubic cutoff 
! for the direct summation
! cubic cutoff:  -ncelldirect to ncelldirect in each direction
! vectors are stored in boxxyz(alpha,nc) alpha = 1,2,3 (x,y,z)
! and nc is the index of the given neighboring cell 
! lcell (nc) is set to 1 if nc is not the central box (ncell !=0)
!
!******************************************************************************

SUBROUTINE direct_sum_init ( rm )

  USE config,   ONLY : box
  USE io_file,  ONLY : ionode , stdout , kunit_OUTFF 

  implicit none
 
  ! global
  TYPE ( rmesh ) :: rm

  ! local
  integer :: nc , ncellx , ncelly , ncellz , ncelldirect

  ncelldirect = rm%ncell

  nc = 0
  rm%lcell = 0
  do ncellx = -ncelldirect,ncelldirect
    do ncelly = -ncelldirect,ncelldirect
      do ncellz = -ncelldirect,ncelldirect
        nc = nc + 1
        rm%boxxyz(1,nc) = box * dble(ncellx)
        rm%boxxyz(2,nc) = box * dble(ncelly)
        rm%boxxyz(3,nc) = box * dble(ncellz)
        if(ncellx .ne. 0 .or. ncelly .ne. 0 .or. ncellz .ne. 0) then
          rm%lcell(nc) = 1
        endif
      enddo
    enddo
  enddo


  if ( nc .ne. rm%ncmax ) then
    if( ionode ) WRITE ( stdout ,'(a,3i7)') 'number of ncells do not match in direct_sum_init', nc , rm%ncmax
  endif

  return

END SUBROUTINE direct_sum_init




END MODULE rspace
