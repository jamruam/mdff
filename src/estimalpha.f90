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

PROGRAM estimate_alpha
  
  implicit none

  integer :: k
  double precision :: alpha , alpha2 , rcut , rcut2 , rcut3 , e

  read(*,*) rcut
  alpha = 4.0d0
  alpha2 = alpha*alpha
  rcut2=rcut*rcut
  rcut3=rcut2*rcut

  e = dexp( -alpha2*rcut2 )
  e = e / alpha2 / rcut3
  e = e * 0.56d0  
  print*,alpha,e  
  do while ( e .le. 1e-7 )
    alpha = alpha - 0.05d0
    alpha2 = alpha*alpha
    e = dexp( -alpha2*rcut2 )
    e = e / alpha2 / rcut3
    e = e * 0.56d0
    print*,alpha,e  
  enddo

END PROGRAM estimate_alpha

