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
!*********************** MODULE CONSTANTS *************************************
!
! Main constants of the code
!
!******************************************************************************
MODULE constants 
  
  double precision, PARAMETER :: dzero  = 0.0d0
  double precision, PARAMETER :: pi     = 3.14159265358979323846264338327950288419716939937510d0
  double precision, PARAMETER :: tpi    = 2.d0*pi
  double complex  , PARAMETER :: citpi  = (0.d0,1.d0)*tpi
  double complex  , PARAMETER :: imag   = (0.d0,1.d0)
  double precision, PARAMETER :: fpi    = 2.d0*tpi
  double precision, PARAMETER :: piroot = dsqrt(pi)

END MODULE constants
! ===== fmV =====
