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

!*********************** MODULE CONSTANTS *************************************
!
! Main constants of the code
!
!******************************************************************************
MODULE constants 
  
  double precision, PARAMETER :: dzero  = 0.0d0
  double precision, PARAMETER :: pi     = 3.14159265358979323846264338327950288419716939937510d0
  double precision, PARAMETER :: pisq   = pi * pi 
  double precision, PARAMETER :: tpi    = 2.d0*pi
  double precision, PARAMETER :: fpi    = 2.d0*tpi    ! 4pi
  double precision, PARAMETER :: piroot = SQRT ( pi ) ! sqrt(pi)
  double precision, PARAMETER :: radian = 180.0d0 / pi 
  double complex  , PARAMETER :: imag   = (0.d0,1.d0) ! imaginary number 
  double complex  , PARAMETER :: mimag  = (0.d0,-1.d0)! negative imaginary number 
  double complex  , PARAMETER :: citpi  = imag*tpi    ! 2*i*pi

  double precision, PARAMETER :: hart   = 27.2113838565563D0  ! ( in eV ) 
  double precision, PARAMETER :: bohr   = 0.52917720859000D0  ! ( in angstrom ) 
  double precision, PARAMETER :: coulombic_factor = 138935.4835d0 ! ( for Dl_POLY comparison )
  ! if length  are in angstom and energy in eV 
  ! the square of the unit charge is :
  double precision, PARAMETER :: e_2    = hart * bohr 

END MODULE constants
! ===== fmV =====
