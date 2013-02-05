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
! ======= Hardware =======

! complementary error function
! W. Press et al. : Numerical Recipes, p. 214
! built-in erfc() is buggy e.g. on some Linux distributions
  double precision FUNCTION errfc(x)
  implicit none
  ! local
  double precision :: x , z , t    
  z=ABS(x)
  t=1d0/(1d0+0.5d0*z)
  errfc=t*EXP(-z*z-1.26551223d0+t*(1.00002368d0+t*(.37409196d0+ &
  t*(.09678418d0+t*(-.18628806d0+t*(.27886807d0+t*(-1.13520398d0+ &
  t*(1.48851587d0+t*(-.82215223d0+t*.17087277d0)))))))))
  if(x.lt.0d0) errfc=2d0-errfc
  return
  END FUNCTION

! error function
  double precision FUNCTION errf(x)
  implicit none
  ! local
  double precision :: x , errfc   
  errf=1.0d0-errfc(x)
  END FUNCTION
! ===== fmV =====
