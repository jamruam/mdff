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
        PROGRAM gene_kpt

        implicit none

        integer :: ikx , iky , ikz
        integer :: nk , wi , nktot
        double precision :: ak , pi

        wi = 1
        pi = 3.14159265358979323846d0

        READ (*,*) nk
        nktot = ( ( nk + 1 ) *( nk + 1 ) * ( nk + 1 ) )

        WRITE (*,*) 'Automatically generated mesh FF'
        WRITE (*,*) nktot
        WRITE (*,*) 'Reciprocal lattice'


        ak=(1.0d0)/(nk) 


        do ikx = 0 , nk
          do iky = 0 , nk        
            do ikz = 0 , nk
              WRITE (*,'(3f16.12,i6)') ikx * ak , iky * ak , ikz * ak, wi
            enddo
          enddo
        enddo

        END PROGRAM gene_kpt

! ===== fmV =====
