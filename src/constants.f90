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

  implicit none

  integer, PARAMETER ,  PUBLIC    :: dp               = selected_real_kind(15,300)  
  real(kind=dp),        PARAMETER :: dzero            = 0.0_dp
  real(kind=dp),        PARAMETER :: pi               = 3.14159265358979323846264338327950288419716939937510_dp
  real(kind=dp),        PARAMETER :: pisq             = pi * pi          ! pi^2
  real(kind=dp),        PARAMETER :: tpi              = 2.0_dp*pi        ! 2pi
  real(kind=dp),        PARAMETER :: fpi              = 2.0_dp*tpi       ! 4pi
  real(kind=dp),        PARAMETER :: piroot           = SQRT ( pi )      ! sqrt(pi)
  real(kind=dp),        PARAMETER :: radian           = 180.0_dp / pi    ! radian 
  complex(kind=dp),     PARAMETER :: imag             = (0.0_dp, 1.0_dp) ! imaginary number 
  complex(kind=dp),     PARAMETER :: mimag            = (0.0_dp,-1.0_dp) ! negative imaginary number 
  complex(kind=dp),     PARAMETER :: citpi            = imag*tpi         ! 2*i*pi

  real(kind=dp),        PARAMETER :: rytoev           = 13.605826_dp     ! ( a.u to eV )
  real(kind=dp),        PARAMETER :: hart             = rytoev*2.0_dp    ! ( Hartree ) 
  real(kind=dp),        PARAMETER :: bohr             = 0.529177249_dp   ! ( a.u in angstrom ) 
  real(kind=dp),        PARAMETER :: coulombic_factor = 138935.4835_dp   ! ( for Dl_POLY comparison )
  real(kind=dp),        PARAMETER :: boltz            = 0.831451115_dp   ! boltzmann constant
  ! if length  are in angstom and energy in eV 
  ! the square of the unit charge is :
  real(kind=dp),        PARAMETER :: e_2              = hart * bohr               ! ( lenght of angstom and energy in eV )
  real(kind=dp),        PARAMETER :: evtoj            = 1.602176487e-19_dp        ! ( eV to J )
  real(kind=dp),        PARAMETER :: hplanck          = 6.62617636e-34_dp         ! ( Planck constant )
  real(kind=dp),        PARAMETER :: eh               = evtoj / 1e14_dp / hplanck !  e / h 
  ! conversion from a.u to V/A^2 ( if length  are in angstom and energy in eV )   !  Cq unit (  a.u to V/A^2 )
  real(kind=dp),        PARAMETER :: CQ_UNIT          = eh / 1000.0_dp 

END MODULE constants
! ===== fmV =====
