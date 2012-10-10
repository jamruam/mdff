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
        PROGRAM disttab

        implicit none

        integer N, PAN
        double precision, dimension (:),allocatable :: a , b , c
        double precision, dimension (:),allocatable :: distrib_a, distrib_b, distrib_c, distrib_all
        double precision amin , amax , da , ia , aaaa
        integer i , ak , ka , iiii

        READ ( * , * ) N , PAN , amin , amax

        WRITE (6, * ) N ,'points'
        WRITE (6, * ) PAN ,'intervals between' , amin , 'and' , amax

        allocate( a(N) , b(N) , c(N) , distrib_a(PAN) , distrib_b(PAN), distrib_c(PAN) , distrib_all(PAN) )

        do i=1,N
        READ ( * , * ) iiii,aaaa,aaaa,aaaa, a(i) , b(i) , c(i)
        enddo


        da = amax - amin
        ia = dble( da / PAN )

        distrib_a = 0

        do i = 1, N
          ak = ( a(i) - amin ) * ( 1.0d0 / ia )
          ka = int( ak ) + 1
          if ( ka .gt. PAN+1 ) then
             WRITE (6, * ) 'ERROR out of bound in distrib_a'
             WRITE (6, * ) i , ka , PAN+1     
             STOP 
          endif
          distrib_a(ka) = distrib_a(ka) + 1
          distrib_all(ka) = distrib_all(ka) + 1
          ak = ( b(i) - amin ) * ( 1.0d0 / ia )
          ka = int( ak ) + 1
          if ( ka .gt. PAN+1 ) then
             WRITE (6, * ) 'ERROR out of bound in distrib_a'
             WRITE (6, * ) i , ka , PAN+1     
             STOP 
          endif
          distrib_b(ka) = distrib_b(ka) + 1
          distrib_all(ka) = distrib_all(ka) + 1
          ak = ( c(i) - amin ) * ( 1.0d0 / ia )
          ka = int( ak ) + 1
          if ( ka .gt. PAN+1 ) then
             WRITE (6, * ) 'ERROR out of bound in distrib_a'
             WRITE (6, * ) i , ka , PAN+1     
             STOP 
          endif
          distrib_c(ka) = distrib_c(ka) + 1
          distrib_all(ka) = distrib_all(ka) + 1
        enddo

        WRITE (6, * ) 'END TAB'

        OPEN ( unit=10000 , file='distrib.out' )

        do i = 0, PAN
          WRITE (10000,'(5f15.8)') amin + dble( i * ia ) , distrib_all(i) / ( ia * dble( 3 * N ) ) , &
                                                           distrib_a(i) / ( ia * dble( 3 * N ) )       , & 
                                                           distrib_b(i) / ( ia * dble( 3 * N ) )       , & 
                                                           distrib_c(i) / ( ia * dble( 3 * N ) ) 
        enddo

        CLOSE (10000)

        WRITE (6, * ) 'file distrib.out'

        deallocate(a , b , c , distrib_all, distrib_a ,distrib_b , distrib_c )

        STOP 

        END PROGRAM disttab
! ===== fmV =====
