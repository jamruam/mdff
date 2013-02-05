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

! ======= Hardware =======
! ======= Hardware =======

MODULE block

  TYPE accu
    double precision :: accval
    double precision :: accvalsq
    double precision :: average
    integer :: counter 
  END TYPE


CONTAINS

!*********************** SUBROUTINE block *************************************
!
! calculates block averages using the method of Flyberg and Petersen
! JCP 91 (1989) pg 461
!
! it reads the file EQUILFF wich store all the thermodynamic quantities after
! equilibration step , or all if ( nequil = 0 )
!
!******************************************************************************

SUBROUTINE block_

  USE md,       ONLY :  npas , nequil
  USE io_file,  ONLY :  ionode , stdout , kunit_EQUILFF

  implicit none

  integer :: i, ii, j, nblock, idum, opbl, nb20 , nstep 

  integer, PARAMETER :: nblockmax = 150000
  integer, PARAMETER :: nquan = 8                    !  number of thermo quantities
  integer, PARAMETER :: div = 20 

  double precision , dimension ( :,: ), allocatable :: bdata
  double precision , dimension ( : )  , allocatable :: av  
  double precision , dimension ( : )  , allocatable :: sav 
  double precision , dimension ( : )  , allocatable :: sum 
  double precision , dimension ( : )  , allocatable :: ssum
  double precision , dimension ( : )  , allocatable :: svv 
  double precision , dimension ( : )  , allocatable :: avv 
  integer          , dimension ( : )  , allocatable :: istep 
  character*6      , dimension ( : )  , allocatable :: dname

  !trash
  integer :: iiii
  double precision :: aaaa 
  character*1 :: cccc 

  nstep = npas - nequil + 1
  if ( nstep .le. 1 ) return

  allocate ( bdata ( nblockmax , nquan )     )
  allocate ( av ( nquan )  , sav   ( nquan ) )
  allocate ( sum ( nquan ) , ssum  ( nquan ) )
  allocate ( svv ( nquan ) , avv   ( nquan ) )
  allocate ( istep (nstep) , dname ( nquan ) )

  if ( ionode ) WRITE ( stdout , '(a)' ) ' '
  if ( ionode ) WRITE ( stdout , '(a)' ) ' '
  if ( ionode ) WRITE ( stdout , '(a)' ) ' ***** Calculate block averages ************'

  ! =======================================================
  !  read thermodynamic quantities at each time in EQUILFF
  ! =======================================================
  OPEN ( UNIT = kunit_EQUILFF , FILE='EQUILFF') 
  do i = 1 , nstep 
    nblock = nblock + 1
    READ( kunit_EQUILFF , * )  iiii     , aaaa , ( dname ( j ) , cccc , bdata( i , j ) , j = 1, 5 ) 
    READ( kunit_EQUILFF , * )  istep(i) , aaaa , aaaa , cccc , aaaa , & 
                               ( dname ( j ) , cccc , bdata( i , j ) , j = 6, 8 )
    do j = 1, nquan
      av(j)  = av(j)  + bdata(i, j)
      sav(j) = sav(j) + bdata(i, j)*bdata(i, j)
    enddo
  enddo
  CLOSE ( kunit_EQUILFF )
  
  if ( ionode ) WRITE ( stdout , *) '  number of blocks ', nblock

  idum = 0 
  nb20 = nblock/div

  do j = 1, nquan 
    avv(j) = 0.D0
    svv(j) = 0.D0
  enddo

  do ii = 1, div

    do j = 1, nquan
      sum(j)  = 0.D0
      ssum(j) = 0.D0
    enddo

    do i = 1, nb20
      idum = idum + 1
      if ( idum .gt. nstep ) then
        if ( ionode ) WRITE ( stdout , '(a)' ) 'ERROR in block_ : out of bound bdata',idum,nstep
        STOP 
      endif  
      do j = 1, nquan
        sum(j)  = sum(j)  + bdata(idum, j)
        ssum(j) = ssum(j) + bdata(idum, j)**2
      enddo

    enddo

    do j = 1, nquan
      avv(j) = avv(j) + sum(j)  / DBLE ( nb20 )
      svv(j) = svv(j) + ssum(j) / DBLE ( nb20 )
    enddo

    if ( ionode ) WRITE ( stdout , 99004) ii, ( sum( j ) / DBLE ( nb20 ) , j=1, nquan )
    if ( ionode ) WRITE ( stdout , 99004) ii, ( SQRT ( ssum ( j ) / DBLE ( nb20 ) - &
                                          sum( j ) * sum( j ) / DBLE ( nb20 ) / DBLE ( nb20 ) )  , j=1, nquan )

  enddo

  if ( ionode ) WRITE ( stdout , '(a)' ) '' 

  do j = 1, nquan
    avv(j) = avv(j)/ DBLE ( div )
    svv(j) = svv(j)/ DBLE ( div ) - avv(j)*avv(j)
    if ( ionode ) WRITE ( stdout , 99003) dname(j), avv(j), SQRT (svv(j)/DBLE (div))
  enddo
  if ( ionode ) WRITE ( stdout , '(a)' ) '' 
 

  do j = 1, nquan
    av(j)  =  av(j) / DBLE ( nblock )
    sav(j) = (sav(j)/ DBLE ( nblock ) ) - av(j)*av(j)
! ==========================================
!           estimate: <c0/(n-1)>
! ==========================================
    av(j) = SQRT (sav(j)/DBLE ( (nblock-1) ) )
    sav(j) = av(j)      / SQRT (2.*DBLE ( (nblock-1.)) )
  END DO

  opbl = 0
  if ( ionode ) WRITE ( stdout , 99002) opbl, (av(j), sav(j), j=1, nquan )

! ================================
!  perform block transformations
! ================================
  do while (nblock.ge.4)

    nblock = nblock/2
    i = 1
    do j = 1, nquan
      av(j) = 0
      sav(j) = 0
    enddo

    do ii = 1, nblock

      do j = 1, nquan
        bdata(ii, j) = (bdata(i,j)+bdata(i+1,j))/2.D0
        av(j)        = av(j)  + bdata(ii, j)
        sav(j)       = sav(j) + bdata(ii, j)*bdata(ii, j)
      enddo
      i = i + 2

    enddo

    do j = 1, nquan
      av(j)  =   av (j) / DBLE ( nblock )                        
      sav(j) = ( sav(j) / DBLE ( nblock ) ) - av(j) * av(j)          ! c0 eq.8  
    enddo

!  if ( ionode )  WRITE ( stdout , 99002) opbl, (av(j), sav(j), j=1, nquan)

    do j = 1, nquan
      av(j) = SQRT ( sav(j)/DBLE ( (nblock-1)) )
      sav(j) = av(j) / SQRT (2.*( DBLE (nblock) -1.))
    enddo

    opbl = opbl + 1

  if ( ionode )  WRITE ( stdout , 99002) opbl, (av(j), sav(j), j=1, nquan)

  enddo

  deallocate ( bdata         )
  deallocate ( av    , sav   )
  deallocate ( sum   , ssum  )
  deallocate ( svv   , avv   )
  deallocate ( istep , dname )

  return
 
!99001 FORMAT (i10,30(e12.4))
99002 FORMAT (i4, 30(f10.6))
99003 FORMAT (' ######', 1x, a6, 2x, f16.8, '   ', f16.8)
99004 FORMAT (' block data ', 1x, i4, 30(1x,f10.6))

END SUBROUTINE block_

END MODULE block
! ===== fmV =====
