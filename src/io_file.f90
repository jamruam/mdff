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
!*********************** MODULE IO_FILE ***************************************
!
! definition of the units for output files
!
!******************************************************************************
MODULE io_file

  ! rank for output (if true myrank.eq.0)
  logical :: ionode        

  ! standard output
  integer, PARAMETER :: stdout          = 6  

  ! standard input file ( code argument control.F )
  integer, PARAMETER :: stdin           = 1001

  ! tmp working file
  integer, PARAMETER :: kunit_tmp       =  8

  ! OUTFF file 
  integer, PARAMETER :: kunit_OUTFF     = 9

  ! thermodynamic info
  integer, PARAMETER :: kunit_OSZIFF    = 10

  ! trajectory (positions, velocities , forces )
  integer, PARAMETER :: kunit_TRAJFF    = 20

  ! input configuration
  integer, PARAMETER :: kunit_POSFF     = 30

  ! end configuration 
  integer, PARAMETER :: kunit_CONTFF    = 40

  ! electric field gradient trajectory
  integer, PARAMETER :: kunit_EFGFF     = 50  
  ! electric field gradient trajectory
  integer, PARAMETER :: kunit_EFGFFIT   = 51  

  ! all efg for each atoms (if lefgprintall)
  integer, PARAMETER :: kunit_EFGALL    = 60 
  ! all efg for each atoms (if lefgprintall)
  integer, PARAMETER :: kunit_EFGALLIT1 = 61
  ! all efg for each atoms (if lefgprintall)
  integer, PARAMETER :: kunit_EFGALLIT2 = 62 

  ! EFG eta distribution (average) 
  integer, PARAMETER :: kunit_DTETAFF   = 70
  ! EFG eta distribution (average) 
  integer, PARAMETER :: kunit_DTETAFFIT = 71

  ! EFG vzz distribution (average)
  integer, PARAMETER :: kunit_DTVZZFF   = 80
  ! EFG vzz distribution (average)
  integer, PARAMETER :: kunit_DTVZZFFIT = 81

  ! EFG tensor component distribution Ui  (average)
  integer, PARAMETER :: kunit_DTIBUFF   = 90

  ! GR radial distribution (average)
  integer, PARAMETER :: kunit_GRTFF     = 100

  ! OPT output configuration after optimisation
  integer, PARAMETER :: kunit_ISCFF     = 110

  ! OPT thermodynamics properties of optimized structure
  integer, PARAMETER :: kunit_ISTHFF    = 120 

  ! VIB eigenvalues (frequencies) of the hessian matrix
  integer, PARAMETER :: kunit_EIGFF     = 130

  ! VIB eigenvector (normal modes) of the hessian matrix 
  integer, PARAMETER :: kunit_VECTFF    = 140

  ! VIB density of states 
  integer, PARAMETER :: kunit_DOSFF     = 150

  ! VIB generated configuration of a given mode
  integer, PARAMETER :: kunit_MODFF     = 160

  ! VIB kpoint mesh for the complete dos
  integer, PARAMETER :: kunit_IBZKPTFF  = 170

  ! VIB "complete" density of states
  integer, PARAMETER :: kunit_DOSKFF    = 180
  integer, PARAMETER :: kunit_DKFF      = 181

  ! VIB fvibcalc output
  integer, PARAMETER :: kunit_VIBFF     = 190

  ! MSD output file
  integer, PARAMETER :: kunit_MSDFF     = 200

  ! stress tensor output file
  integer, PARAMETER :: kunit_STRESSFF  = 210

  ! velocity auto-correlation output file
  integer, PARAMETER :: kunit_VACFFF    = 220
 
  ! EFG auto-correlation output file
  integer, PARAMETER :: kunit_EFGACFFF  = 230

  ! static structure factor
  integer, PARAMETER :: kunit_STRFACFF  = 240

  ! static structure factor
  integer, PARAMETER :: kunit_EQUILFF   = 250

  ! mean number of atoms in a shell of width at distance r
  integer, PARAMETER :: kunit_NRTFF     = 260

CONTAINS

!*********************** SUBROUTINE io_init ***********************************
!
! initialize the ionode logical variable and open OUTFF file
!
! if ( ionode ) WRITE(unit,*)  
!
!******************************************************************************
SUBROUTINE io_init

  implicit none
  INCLUDE 'mpif.h'

  ! local 
  integer :: myrank , ierr
  
  CALL MPI_COMM_RANK ( MPI_COMM_WORLD , myrank , ierr )    ! numero de processus (output myrank .eq. 0 )

   if ( myrank .eq. 0 ) then
     ionode = .true.
   else
     ionode = .false.
   endif
  ! ===================
  !  open main file
  ! ===================
  OPEN (unit = kunit_OUTFF  ,file = 'OUTFF',STATUS = 'UNKNOWN')

  return

END SUBROUTINE io_init

!*********************** SUBROUTINE io_end ***********************************
!
! close OUTFF file
!
!******************************************************************************

SUBROUTINE io_end

  implicit none

  CLOSE(kunit_OUTFF)

  return
 
END SUBROUTINE io_end

END MODULE io_file
! ===== fmV =====
