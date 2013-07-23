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

! *********************** SUBROUTINE init_random_seed **************************
!> \brief
!> initialisation of seed for the fortran random number generator RANDOM_NUMBER
! ******************************************************************************
SUBROUTINE init_random_seed

  implicit none

  ! local
  integer :: i , n , clock
  integer, dimension(:), allocatable :: seed
  
  i=1       
  CALL RANDOM_SEED(size = n)
  allocate(seed(n))
          
  CALL SYSTEM_CLOCK(COUNT = clock)
          
  seed = clock + 48 * (/ (i - 1, i = 1, n) /)
  CALL RANDOM_SEED(PUT = seed)
          
  deallocate(seed)

  return

END SUBROUTINE

! *********************** SUBROUTINE knuth *************************************
!> \brief
!! Random variate from the standard normal distribution.
!! The distribution is gaussian with zero mean and unit variance         
!> \note
!! Knuth D, The Art of Computer Programming , (2nd edition Addison-Wesley), 1978                                     
!> \note
!! Adapted from F.24 (Allen-Tildesley)
!> \author
!! Allen-Tildesley
!> \param[in] mean mean value of the gaussian distribution
!> \param[in] sigma variance value of the gaussian distribution
!> \param[out] G random number from the gaussian distribution define by mean and sigma
! ******************************************************************************
SUBROUTINE knuth ( G , mean , sigma )

  USE constants, ONLY : dp
  implicit none

  ! global
  real(kind=dp), intent (out) :: G
  real(kind=dp), intent (in)  :: mean, sigma

  ! local
  real(kind=dp) :: A1, A3, A5, A7, A9
  PARAMETER ( A1 = 3.949846138_dp, A3 = 0.252408784_dp )
  PARAMETER ( A5 = 0.076542912_dp, A7 = 0.008355968_dp )
  PARAMETER ( A9 = 0.029899776_dp                   )
  real(kind=dp) :: x
  real(kind=dp) :: summ, r, r2
  integer :: i
  integer :: iseed

  summ = 0.0_dp
  do i = 1, 12
     CALL RANDOM_SEED(SIZE = iseed)
     CALL RANDOM_NUMBER(HARVEST = x)
     summ = summ + x
  enddo

  r  = ( summ - 6.0_dp ) / 4.0_dp
  r2 = r * r

  G = ((((A9*R2+A7)*R2+A5)*R2+A3)*R2+A1)*R

  G =  mean + G*sigma

  return

END SUBROUTINE knuth 

! *********************** SUBROUTINE boxmuller_polar ***************************
!> \brief
!! Box-Muller polar method
!> \note
!! http://en.wikipedia.org/wiki/Box-Muller_transform
!> \param[in] mean mean value of the gaussian distribution
!> \param[in] sigma variance value of the gaussian distribution
!> \param[out] G random number from the gaussian distribution define by mean and sigma
! ******************************************************************************

SUBROUTINE boxmuller_polar (G, mean, sigma)

  USE constants,        ONLY : dp 
  implicit none

  ! global
  real(kind=dp), intent (out) :: G
  real(kind=dp), intent (in)  :: mean, sigma

  ! local
  real(kind=dp) :: G1,U,V,S
  integer :: iseed
  
  CALL RANDOM_SEED(SIZE = ISEED)
  CALL RANDOM_NUMBER(HARVEST = U)
  CALL RANDOM_SEED(SIZE = ISEED)
  CALL RANDOM_NUMBER(HARVEST = V)
  U = (2.0_dp*U)-1.0_dp
  V = (2.0_dp*V)-1.0_dp
  S = U*U+V*V
  do while ( S .eq. 0 .or. S .ge. 1) 
    CALL RANDOM_SEED(SIZE = iseed)
    CALL RANDOM_NUMBER(HARVEST = U)
    CALL RANDOM_SEED(SIZE = iseed)
    CALL RANDOM_NUMBER(HARVEST = V)
    U = (2.0_dp*U)-1.0_dp
    V = (2.0_dp*V)-1.0_dp
    S = U * U + V * V 
  enddo
  G1 = -2.0_dp * LOG ( S )
  G = U * SQRT ( G1 / S ) 
  G = mean + G * sigma

  return

END SUBROUTINE boxmuller_polar

! *********************** SUBROUTINE boxmuller_basic ***************************
!> \brief
!! Box-Muller cartesian method
!> \note
!! http://en.wikipedia.org/wiki/Box-Muller_transform
!> \param[in] mean mean value of the gaussian distribution
!> \param[in] sigma variance value of the gaussian distribution
!> \param[out] G random number from the gaussian distribution define by mean and sigma
! ******************************************************************************
SUBROUTINE boxmuller_basic (G, mean, sigma)
  
  USE constants,        ONLY : dp , tpi
  implicit none

  ! global
  real(kind=dp), intent (out) :: G
  real(kind=dp), intent (in)  :: mean, sigma

  ! local
  real(kind=dp) :: C,U,V,R
  integer :: iseed

  CALL RANDOM_SEED(SIZE = iseed)
  CALL RANDOM_NUMBER(HARVEST = U)
  CALL RANDOM_SEED(SIZE = iseed)
  CALL RANDOM_NUMBER(HARVEST = V)
  
  R = SQRT ( -2.0_dp * LOG ( U ) )
  C = COS ( tpi * V )
  G = R * C
  G = mean + G * sigma

  return

END SUBROUTINE boxmuller_basic

! ===== fmV =====
