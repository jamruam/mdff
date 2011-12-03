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
!*********************** SUBROUTINE init_random_seed **************************
!
! initialisation of seed 
!
!******************************************************************************

SUBROUTINE init_random_seed

  implicit none

  ! local
  integer :: i, n, clock
  integer, dimension(:), allocatable :: seed
         
  CALL RANDOM_SEED(size = n)
  allocate(seed(n))
          
  CALL SYSTEM_CLOCK(COUNT = clock)
          
  seed = clock + 48 * (/ (i - 1, i = 1, n) /)
  CALL RANDOM_SEED(PUT = seed)
          
  deallocate(seed)

  return

END SUBROUTINE

!*********************** SUBROUTINE knuth *************************************
!
! adapted from F.24 (Allen-Tildesley)
!
! RANDOM VARIATE FROM THE STANDARD NORMAL DISTRIBUTION.         
! THE DISTRIBUTION IS GAUSSIAN WITH ZERO MEAN AND UNIT VARIANCE.
! REFERENCE:                                                    
!     KNUTH D, THE ART OF COMPUTER PROGRAMMING, (2ND EDITION    
!     ADDISON-WESLEY), 1978                                     
!
!******************************************************************************

SUBROUTINE knuth ( G , mean , sigma )

  implicit none

  ! global
  double precision, intent (out) :: G
  double precision, intent (in)  :: mean, sigma

  ! local
  double precision :: A1, A3, A5, A7, A9
  PARAMETER ( A1 = 3.949846138d0, A3 = 0.252408784d0 )
  PARAMETER ( A5 = 0.076542912d0, A7 = 0.008355968d0 )
  PARAMETER ( A9 = 0.029899776d0                   )
  double precision :: x
  double precision :: summ, r, r2
  integer :: i
  integer :: iseed

  summ = 0.0d0
  do i = 1, 12
     CALL RANDOM_SEED(SIZE = iseed)
     CALL RANDOM_NUMBER(HARVEST = x)
     summ = summ + x
  enddo

  r  = ( summ - 6.0d0 ) / 4.0d0
  r2 = r * r

  G = ((((A9*R2+A7)*R2+A5)*R2+A3)*R2+A1)*R

  G =  mean + G*sigma

  return

END SUBROUTINE knuth 

!*********************** SUBROUTINE boxmuller_polar ***************************
!
!
!******************************************************************************

SUBROUTINE boxmuller_polar (G, mean, sigma)

  implicit none

  ! global
  double precision, intent (out) :: G
  double precision, intent (in)  :: mean, sigma

  ! local
  double precision :: G1,U,V,S
  integer :: iseed
  
  CALL RANDOM_SEED(SIZE = ISEED)
  CALL RANDOM_NUMBER(HARVEST = U)
  CALL RANDOM_SEED(SIZE = ISEED)
  CALL RANDOM_NUMBER(HARVEST = V)
  U = (2.0d0*U)-1.0d0
  V = (2.0d0*V)-1.0d0
  S = U*U+V*V
  do while ( S .eq. 0 .or. S .ge. 1) 
  CALL RANDOM_SEED(SIZE = iseed)
  CALL RANDOM_NUMBER(HARVEST = U)
  CALL RANDOM_SEED(SIZE = iseed)
  CALL RANDOM_NUMBER(HARVEST = V)
  U = (2.0d0*U)-1.0d0
  V = (2.0d0*V)-1.0d0
  S = U*U+V*V
  enddo
  G1 = -2.0d0*dlog(S)
  G = U*DSQRT(G1/S) 
  G = mean + G*sigma

  return

END SUBROUTINE boxmuller_polar

!*********************** SUBROUTINE boxmuller_basic ***************************
!
!
!******************************************************************************

SUBROUTINE boxmuller_basic (G, mean, sigma)
  
  USE constants,        ONLY : pi

  implicit none

  ! global
  double precision, intent (out) :: G
  double precision, intent (in)  :: mean, sigma

  ! local
  double precision :: C,U,V,R
  integer :: iseed

  CALL RANDOM_SEED(SIZE = iseed)
  CALL RANDOM_NUMBER(HARVEST = U)
  CALL RANDOM_SEED(SIZE = iseed)
  CALL RANDOM_NUMBER(HARVEST = V)
  
  R = dsqrt(-2.0d0*log(U))
  C = dcos(2.0d0*pi*V)
  G = R*C
  G = mean + G*sigma

  return

END SUBROUTINE boxmuller_basic

! ===== fmV =====
