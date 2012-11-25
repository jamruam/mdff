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

!*********************** SUBROUTINE estimate_alpha ***********************
!
! This routine estiate the alpha parameter in the Ewald summation.
! This estimation is based on the cell size
! The parameter should be anyway optimised for each problem
!
!******************************************************************************

SUBROUTINE estimate_alpha(alpha)
  
  USE config,   ONLY : box
  implicit none

  double precision :: alpha , alpha2 , rcut , rcut2 , rcut3 , e

  rcut = box * 0.5d0
  alpha = 4.0d0
  alpha2 = alpha*alpha
  rcut2=rcut*rcut
  rcut3=rcut2*rcut

  e = dexp( -alpha2*rcut2 )
  e = e / alpha2 / rcut3
  e = e * 0.56d0  
  do while ( e .le. 1e-7 )
    alpha = alpha - 0.05d0
    alpha2 = alpha*alpha
    e = dexp( -alpha2*rcut2 )
    e = e / alpha2 / rcut3
    e = e * 0.56d0
  enddo

  return

END SUBROUTINE estimate_alpha



!*********************** SUBROUTINE do_split ***********************************
!
! this SUBROUTINE split the number of atoms in for each np procs
! iastart and iaend are the atom index for proc myrank
! WARNING : c'est moi qui l'ai fait ;)
!
! input : 
!          *  n           = number of atoms
!          *  mrank       = local index proc
!          *  np          = number of procs  
! output :
!          *  iastart , iaend = starting and ending atom index for proc mrank
!******************************************************************************

SUBROUTINE do_split ( n , mrank , np , iastart , iaend )

  USE io_file,  ONLY :  ionode , stdout, kunit_OUTFF

  implicit none

  ! global
  integer :: np , mrank
  integer :: iastart , iaend , n

  ! local 
  integer :: imin,imax
  integer :: istartV(0:np-1),iendV(0:np-1)
  integer :: splitnumberV(0:np-1),isteps,x,y,me

  imin = 1
  imax = n      

  isteps = (imax-imin)+1
  x = isteps/np
  y = mod(isteps,np)

  do me = 0,np-1
     if ((me.eq.0).or.(me.gt.y)) then
        splitnumberV(me) = x
     else if ((me.gt.0).or.(me.lt.(y+1))) then
        splitnumberV(me) = x+1
     endif
  enddo
  do me = 0,np-1
     if (me.eq.0) then
        istartV(0) = imin
        iendV(0)   = imin + ( x -1 )
     else if (me.gt.0) then
        istartV(me) = iendV(me-1) + 1
        iendV(me)   = istartV(me) + splitnumberV(me) - 1
     endif
  enddo

  iastart = istartV(mrank)
  iaend = iendV(mrank)

  if ( ionode ) then
    WRITE ( stdout ,'(a)') 'paralelisation - atom decomposition'
    WRITE ( kunit_OUTFF ,'(a)') 'paralelisation - atom decomposition'
    do me = 0,np-1
      WRITE ( stdout ,'(a5,i4,a5,i5,a3,i5)') 'rank = ',me,'atom',istartV(me),'to',iendV(me)
      WRITE ( kunit_OUTFF ,'(a5,i4,a5,i5,a3,i5)') 'rank = ',me,'atom',istartV(me),'to',iendV(me)
    enddo     
  endif

  return

END SUBROUTINE do_split


!*********************** SUBROUTINE distance_tab ******************************
!
! this subroutine calculates distance between atoms and check if the smallest distance 
! is not too small ( i.e < sigmaAA*0.001d0) 
!
!******************************************************************************

SUBROUTINE distance_tab ( kunit )

  USE config,   ONLY :  natm , box , rx , ry , rz
  USE field,    ONLY :  sigmalj
  USE io_file,  ONLY :  ionode , stdout

  implicit none
  
  ! global
  integer :: kunit

  ! local
  double precision :: rxi, ryi, rzi
  double precision :: rxij, ryij, rzij, rij, rijsq, norm
  integer :: nxij, nyij, nzij 
  integer :: ia , ja , PANdis, kdis
  integer, dimension (:) ,allocatable :: dist
  double precision :: resdis,mindis 

  resdis = 0.5D0 ! should still local no need to be controled

  PANdis = box/resdis

  allocate ( dist ( 0:PANdis ) )
  dist = 0
  mindis = 100000.0D0

  do ia = 1 , natm - 1
    rxi = rx ( ia )
    ryi = ry ( ia )
    rzi = rz ( ia )
    do ja = ia + 1 , natm
      if ( ia .ne. ja ) then
        rxij = rxi - rx ( ja )
        ryij = ryi - ry ( ja )
        rzij = rzi - rz ( ja )
        nxij = nint( rxij / box )
        nyij = nint( ryij / box )
        nzij = nint( rzij / box )
        rxij = rxij - box * nxij
        ryij = ryij - box * nyij
        rzij = rzij - box * nzij
        rijsq = rxij *rxij + ryij * ryij + rzij * rzij
        rij = dsqrt ( rijsq )  
        mindis = min( mindis , rij )
        if ( rij .lt. sigmalj(1,1) * 0.001d0 ) then
          if ( ionode ) &
          WRITE ( stdout ,'(a,i5,a,i5,a,f12.6)') &
          'ERROR: DISTANCE between atoms', ia ,' and ', ja ,' is very small',rij   
          STOP 
        endif
        rij = rij / resdis
        kdis = int( rij )
        if ( kdis .lt. 0 .or. kdis .gt. PANdis ) then
          if ( ionode ) WRITE ( stdout ,*) 'ERROR: out of bound dist in SUBROUTINE distance_tab'
        endif
        dist ( kdis ) = dist ( kdis ) + 1
      endif
    enddo 
  enddo

  norm = natm * ( natm - 1 ) / 2
 
  if ( ionode ) then
    WRITE ( kunit ,'(a)')       'distance check subroutine'
    WRITE ( kunit ,'(a,f6.2)')  'smallest distance = ',mindis
!    WRITE ( kunit ,'(a)')       'distance ditribution:'
!    WRITE ( kunit ,'(a)')       'dist  pct(%)'      
!    do i = 0, PANdis
!      WRITE ( kunit ,'(2f6.2)') i * resdis , dist( i ) / norm * 100.0D0
!    enddo
    WRITE ( kunit ,'(a)')       ''
  endif
  
  deallocate(dist)

  return

END SUBROUTINE distance_tab


!*********************** SUBROUTINE vnlist_pbc ********************************
!
! verlet list subroutine :
! Periodic boundaries condition version ( minimum image convention ) 
! Parallelized and tested
! 
! input :
!          iastart , iaend : index of atom decomposition ( parallel )
!
! output:
!  
!          point           : array of size natm+1 
!                            gives the starting and finishing index of array list for a given atom i
!                            jbegin = point(i)  jend = point(i+1) - 1
!          list            : index list of neighboring atoms
!
! how to use it :
!                          do ia = 1, natm
!                            jbegin = point(i)
!                            jend = point(i+1) - 1
!                            do jvnl = jbegin , jend
!                              ja = list ( jvnl ) 
!                              then ia en ja are neighboors   
!                            enddo
!                          enddo
!
!******************************************************************************

SUBROUTINE vnlist_pbc ( iastart , iaend )!, list , point )

  USE config,   ONLY :  natm , natmi , rx , ry , rz , box , itype, list , point, ntype 
  USE control,  ONLY :  skindiff , cutoff

  implicit none

  ! global
  integer , intent (in)  :: iastart , iaend
!  integer , intent (out) :: list(natm*250), point(natm+1)
 
  ! local
  integer :: icount , ia , ja , it , jt , k
  integer :: p1 , p2
  integer :: nxij , nyij , nzij
  double precision  :: rskinsq(ntype,ntype) , rcut(ntype,ntype) , rskin(ntype,ntype)
  double precision :: rxi , ryi , rzi , rxij , ryij , rzij , rijsq

  do jt = 1, ntype 
    do it = 1, ntype
       rcut    ( it , jt ) = cutoff
       rskin   ( it , jt ) = rcut  ( it , jt ) + skindiff
       rskinsq ( it , jt ) = rskin ( it , jt ) * rskin ( it , jt )
    enddo
  enddo

  icount = 1
  do ia = iastart , iaend
    rxi = rx ( ia )
    ryi = ry ( ia )
    rzi = rz ( ia )
    k = 0
    do ja = 1 , natm
      if ( ( ia .gt. ja .and. ( mod ( ia + ja , 2 ) .eq. 0 ) ) .or. &
           ( ia .lt. ja .and. ( mod (ia + ja , 2 ) .ne.0 ) ) ) then
        rxij = rxi - rx ( ja )
        ryij = ryi - ry ( ja )
        rzij = rzi - rz ( ja )
        nxij = nint( rxij / box )
        nyij = nint( ryij / box )
        nzij = nint( rzij / box )
        rxij = rxij - box * nxij
        ryij = ryij - box * nyij
        rzij = rzij - box * nzij
        rijsq = rxij * rxij + ryij * ryij + rzij * rzij
        p1 = itype ( ia )
        p2 = itype ( ja )
        if ( rijsq .le. rskinsq(p1,p2)) then
          icount = icount+1
          k = k+1
          list(icount-1) = ja
        endif
      endif
    enddo
    point(ia) = icount-k
  enddo
  point (iaend + 1 ) = icount

  return

END SUBROUTINE vnlist_pbc

!*********************** SUBROUTINE vnlist_nopbc ******************************
!
! no periodic boundaries condition version
! same as vnlist_pbc but no periodic boundaries
!
!******************************************************************************

SUBROUTINE vnlist_nopbc ( iastart , iaend )!, list , point )

  USE config,   ONLY :  natm , natmi , rx , ry , rz , box , itype , list , point , ntype
  USE control,  ONLY :  skindiff , cutoff

  implicit none

  ! global
  integer :: iastart , iaend
!  integer , intent (out) :: list(natm*250), point(natm+1)

  ! local
  integer :: icount , ia , ja , it , jt , k
  integer :: p1,p2
  double precision :: rskinsq ( ntype , ntype ) , rcut ( ntype , ntype ) , rskin ( ntype , ntype )
  double precision :: rxi,ryi,rzi,rxij,ryij,rzij,rijsq

  do jt = 1, ntype
    do it = 1, ntype
       rcut    ( it , jt ) = cutoff
       rskin   ( it , jt ) = rcut  ( it , jt ) + skindiff
       rskinsq ( it , jt ) = rskin ( it , jt ) * rskin ( it , jt )
    enddo
  enddo

  icount = 1
  do ia = iastart , iaend
    rxi = rx ( ia )
    ryi = ry ( ia )
    rzi = rz ( ia )
    k = 0
    do ja = 1 , natm
      if ( ( ia .gt. ja .and. ( mod ( ia + ja , 2 ) .eq. 0 ) ) .or. &
           ( ia .lt. ja .and. ( mod ( ia + ja , 2 ) .ne. 0 ) ) ) then
        rxij = rxi - rx ( ja )
        ryij = ryi - ry ( ja )
        rzij = rzi - rz ( ja )
        rijsq = rxij * rxij + ryij * ryij + rzij * rzij
        p1 = itype ( ia )
        p2 = itype ( ja )
        if (rijsq .le. rskinsq(p1,p2)) then
          icount = icount + 1
          k = k + 1
          list ( icount - 1 ) = ja
        endif
      endif
    enddo
    point (ia ) = icount-k
  enddo
  point( iaend + 1 ) = icount

  return

END SUBROUTINE vnlist_nopbc

!*********************** SUBROUTINE vnlist_noimg ********************************
!
! verlet list subroutine :
! Periodic boundaries condition version ( no minimum image convention ) 
! !!!!test version 
! 
! input :
!          iastart , iaend : index of atom decomposition ( parallel )
!
! output:
!  
!          point           : array of size natm+1 
!                            gives the starting and finishing index of array list for a given atom i
!                            jbegin = point(i)  jend = point(i+1) - 1
!          list            : index list of neighboring atoms
!
! how to use it :
!                          do ia = 1, natm
!                            jbegin = point(i)
!                            jend = point(i+1) - 1
!                            do jvnl = jbegin , jend
!                              ja = list ( jvnl ) 
!                              then ia en ja are neighboors   
!                            enddo
!                          enddo
!
!******************************************************************************

SUBROUTINE vnlist_noimg ( iastart , iaend )

  USE config,   ONLY :  natm , natmi , rx , ry , rz , box , itype, list , point, ntype
  USE control,  ONLY :  skindiff , cutoff

  implicit none

  ! global
  integer , intent (in)  :: iastart , iaend

  ! local
  integer :: icount, ia , ja , it , jt ,k
  integer :: p1,p2
  integer :: nxij,nyij,nzij
  double precision  :: rskinsq ( ntype , ntype ) , rcut ( ntype , ntype ) , rskin ( ntype , ntype )
  double precision :: rxi,ryi,rzi,rxij,ryij,rzij,rijsq

  do jt = 1, ntype 
    do it = 1, ntype
       rcut    ( it , jt ) = cutoff
       rskin   ( it , jt ) = rcut  ( it , jt ) + skindiff
       rskinsq ( it , jt ) = rskin ( it , jt ) * rskin ( it , jt )
    enddo
  enddo

  icount = 1
  do ia = iastart , iaend
    rxi = rx ( ia )
    ryi = ry ( ia )
    rzi = rz ( ia )
    k = 0
    do ja = 1 , natm
      if ( ( ia .gt. ja .and. ( mod ( ia + ja , 2 ) .eq. 0 ) ) .or. &
           ( ia .lt. ja .and. ( mod ( ia + ja , 2 ) .ne. 0 ) ) ) then
        rxij = rxi - rx ( ja )
        ryij = ryi - ry ( ja )
        rzij = rzi - rz ( ja )
        nxij = nint ( rxij / box )
        nyij = nint ( ryij / box )
        nzij = nint ( rzij / box )
        rxij = rxij - box * nxij
        ryij = ryij - box * nyij
        rzij = rzij - box * nzij
        rijsq = rxij * rxij + ryij * ryij + rzij * rzij
        p1 = itype ( ia )
        p2 = itype ( ja )
        if ( rijsq .le. rskinsq(p1,p2)) then
          icount = icount+1
          k = k+1
          list(icount-1) = ja
        endif
      endif
    enddo
    point ( ia ) = icount-k
  enddo
  point (iaend + 1 ) = icount

  return

END SUBROUTINE vnlist_noimg


!*********************** SUBROUTINE vnlistcheck *******************************
!
! check wether verlet list should be updated
!
!******************************************************************************

SUBROUTINE vnlistcheck ( iastart , iaend ) !, list , point )

  USE config,   ONLY :  natm , rx , ry , rz , xs , ys , zs , list , point
  USE control,  ONLY :  lpbc , lminimg , skindiff 
  USE md,       ONLY :  updatevnl
  USE time,     ONLY :  vnlisttimetot

  implicit none
  INCLUDE 'mpif.h'

  ! global
  integer, intent (in) :: iastart , iaend 
!  integer, intent (inout) :: list(natm*250), point(natm+1)

  ! local
  integer :: ia , ierr
  double precision :: displ
  double precision :: rxsi,rysi,rzsi
  double precision :: ttt1 , ttt2


  ttt1 = MPI_WTIME(ierr)

  displ = 0.0D0
  do ia = 1, natm
    rxsi = dabs( rx ( ia ) - xs ( ia ) )
    rysi = dabs( ry ( ia ) - ys ( ia ) )
    rzsi = dabs( rz ( ia ) - zs ( ia ) )
    if ( rxsi .gt. DISPL ) DISPL = rxsi
    if ( rysi .gt. DISPL ) DISPL = rysi
    if ( rzsi .gt. DISPL ) DISPL = rzsi
  enddo 
        
  ! ========================================
  !  DISPL = 2.0 * SQRT  ( 3.0 * DISPL ** 2 ) 
  !  I don't know where this come from, 
  !  It was in the very old version
  ! ========================================

  if ( displ .ge. skindiff * 0.5d0 ) then
    updatevnl = updatevnl + 1
    if ( lpbc ) then
!      if ( lminimg ) then 
        CALL vnlist_pbc ( iastart, iaend )!, list , point )
!      else
!        CALL vnlist_noimg ( iastart, iaend )
!      endif
    else
    CALL vnlist_nopbc ( iastart, iaend )!, list , point )
    endif
     
    do ia = 1, natm 
      xs ( ia ) = rx ( ia )
      ys ( ia ) = ry ( ia )
      zs ( ia ) = rz ( ia )
    END do
  endif

  ttt2 = MPI_WTIME(ierr)
  vnlisttimetot = vnlisttimetot + ( ttt2 - ttt1 ) 

  return

END SUBROUTINE vnlistcheck

!*********************** SUBROUTINE print_tensor ******************************
!
! subroutine which print an (3,3) array in a tensor format 
! the trace is also given in output 
!
!******************************************************************************

SUBROUTINE print_tensor( tens , key )

  USE io_file,  ONLY :  ionode , stdout

  implicit none

  ! global
  double precision :: tens(3,3)
  character*6 :: key

  ! local
  integer :: i

  if ( ionode ) then
    WRITE ( stdout ,'(a)') ''
    WRITE ( stdout ,'(a)') key
    do i = 1 , 3
      WRITE ( stdout ,'(3f15.8)') tens(i,1) , tens(i,2) , tens(i,3)
    enddo
    WRITE ( stdout ,'(a,f15.8)') 'trace=',(tens(1,1) + tens(2,2) + tens(3,3))/3.0d0
    WRITE ( stdout , '(a)' ) ''
  endif

  return

END SUBROUTINE print_tensor

!*********************** SUBROUTINE merge_sort ********************************
!
!  adapted from :
!  http://rosettacode.org/wiki/Sorting_algorithms/Merge_sort#Fortran
! 
!  I changed the routine to keep the initial labels during the sort process
!
!******************************************************************************

SUBROUTINE merge_1(A,NA,B,NB,C,NC,labela,labelb,labelc)

  implicit none
    
  ! global
  integer, intent(in)              :: NA,NB,NC                  ! Normal usage: NA+NB = NC
  double precision, intent(in out) :: A(NA)        ! B overlays C(NA+1:NC)
  double precision, intent(in)     :: B(NB)
  double precision, intent(in out) :: C(NC)
  integer, intent(in out)          :: labelA(NA)       
  integer, intent(in)              :: labelB(NB)
  integer, intent(in out)          :: labelC(NC)

  ! local
  integer :: i,j,k

  i = 1; j = 1; k = 1;
  do while(i .le. NA .and. j .le. NB)
    if (A ( i ) .le. B ( j )) then
      C ( k ) = A ( i )
      labelc ( k ) = labela ( i ) 
      i = i + 1
    else
      C ( k ) = B ( j )
      labelc ( k ) = labelb ( j )
      J = J+1
    endif
    k = k + 1
  enddo
  do while (i .le. NA)
    C ( k ) = A ( i )
    labelc ( k ) = labela ( i )
    i = i + 1
    k = k + 1
  enddo

  return
 
END SUBROUTINE merge_1
 
RECURSIVE SUBROUTINE merge_sort(A,N,T,labela,labelt)

  implicit none 

  ! global
  integer, intent(in)                                :: N
  double precision, dimension(N), intent(in out)     :: A
  double precision, dimension((N+1)/2), intent (out) :: T
  integer, dimension(N), intent(in out)              :: labelA
  integer, dimension((N+1)/2), intent(out)           :: labelt

  ! local
  double precision :: V
  integer :: NA,NB,labelv
 
  if (N .lt. 2) return
  if (N .eq. 2) then
    if (A(1) .gt. A(2)) then
      V = A(1)
      labelv=labelA(1)
      A(1) = A(2)
      labela(1)=labelA(2)
      A(2) = V
      labela(2)=labelv
    endif
    return
  endif      
  NA=(N+1)/2
  NB=N-NA

  call merge_sort(A,NA,T,labela,labelt)
  call merge_sort(A(NA+1),NB,T,labela(NA+1),labelt)
 
  if (A(NA) .gt. A(NA+1)) then
    T(1:NA)=A(1:NA)
    labelt(1:NA)=labelA(1:NA)
    call merge_1(T,NA,A(NA+1),NB,A,N,labelt,labela(NA+1),labela)
  endif

  return
 
end subroutine merge_sort
 
!*********************** SUBROUTINE print_config_sample ***********************
!
! print a sample of the configurations (pos , vel , force ) 
! essentially for debug purpose
!
! input : 
!          time ,rank :  just to have some information in output
!
!******************************************************************************

SUBROUTINE print_config_sample ( time , rank )

  USE config,   ONLY :  natm , atype , itype , rx , vx , fx , qia , dipia , ipolar
  USE control,  ONLY :  myrank
  USE io_file,  ONLY :  stdout

  implicit none

  ! global
  integer, intent(in) :: time

  ! local 
  integer :: ia , rank

  if ( myrank.eq.rank ) then
       WRITE ( stdout ,'(a)') '==================================================================================================================================='
       WRITE ( stdout ,'(a)') 'debug SAMPLE OF THE CONFIGIURATION debug     '
       WRITE ( stdout ,'(a5,i10)') 'time = ',time
       WRITE ( stdout ,'(a5,i10)') 'rank = ',rank
       WRITE ( stdout ,'(a)') '     i    atype       itype      ipolar      q      mu_x    mu_y    mu_z             rx                 vx                  fx'
    if ( natm .ge. 8)   &
       WRITE ( stdout ,'(i6,a10,2i10,4x,4f8.3,3f20.10)') &
       ( ia , atype ( ia ) , itype ( ia ) , ipolar ( ia ) , qia ( ia ) , dipia ( ia , 1 ), dipia ( ia , 2) ,dipia ( ia , 3 ), &
        rx ( ia ) , vx ( ia ) , fx ( ia ) , ia = 1 , 4 )
    if ( natm .ge. 8)   &
       WRITE ( stdout ,'(i6,a10,2i10,4x,4f8.3,3f20.10)') &
       ( ia , atype ( ia ) , itype ( ia ) , ipolar ( ia ) , qia ( ia ) , dipia ( ia , 1 ), dipia ( ia , 2) ,dipia ( ia , 3 ), &
        rx ( ia ) , vx ( ia ) , fx ( ia ) , ia = natm - 4  , natm )
    if ( natm .lt. 8)   &
       WRITE ( stdout ,'(i6,a10,2i10,4x,4f8.3,3f20.10)') &
       ( ia , atype ( ia ) , itype ( ia ) , ipolar ( ia ) , qia ( ia ) , dipia ( ia , 1 ), dipia ( ia , 2) ,dipia ( ia , 3 ), &
        rx ( ia ) , vx ( ia ) , fx ( ia ) , ia = 1 , natm )
       WRITE ( stdout ,'(a)') ' '
       WRITE ( stdout ,'(a)') '==================================================================================================================================='
  endif

  return

END SUBROUTINE print_config_sample

 
!*********************** SUBROUTINE  ********************************
!
!******************************************************************************

SUBROUTINE print_general_info(kunit)

  USE io_file,  ONLY :  ionode 
  USE config,   ONLY : natm , ncell , ntype , rho , box , omega

  implicit none

  integer :: kunit

  if ( ionode ) then
    WRITE ( kunit ,'(a)')          ''
    WRITE ( kunit ,'(a)')          'Remind some parameters of the system:'
    WRITE ( kunit ,'(a,i5)')       'natm  = ', natm
    if ( ncell .ne. 0 ) WRITE ( kunit ,'(a,i5,a)')       'ncell = ', ncell
    WRITE ( kunit ,'(a,i5)')       'ntype = ', ntype
    WRITE ( kunit ,'(a,f10.3)')    'rho   = ', rho
    WRITE ( kunit ,'(a,f10.3)')    'box   = ', box
    WRITE ( kunit ,'(a,f10.3)')    'vol   = ', omega
    WRITE ( kunit ,'(a)')          ''
  endif

  return

END SUBROUTINE print_general_info 


!*********************** SUBROUTINE dumb_guy ********************************
!
!  This subroutine permits to print out the dumb guy!
!  Here is the guy ...
!
!******************************************************************************

SUBROUTINE dumb_guy(kunit)

  USE io_file,  ONLY :  ionode 

  implicit none

  integer :: kunit

  if ( ionode ) then
     WRITE ( kunit ,'(a)') '                            \\|//                    '
     WRITE ( kunit ,'(a)') '                           -(o o)-                           '
     WRITE ( kunit ,'(a)') '========================oOO==(_)==OOo========================'
  endif

  return

END SUBROUTINE dumb_guy

! **********************************************************************
!  to remove leading and trailing spaces
! **********************************************************************
character(30) FUNCTION sweep_blanks(in_str)

  implicit none

  character(*), intent(in) :: in_str
  character(30)            :: out_str
  character                :: ch
  integer                  :: j

  out_str = " "

  do j=1, len_trim(in_str)

  ! get j-th char
    ch = in_str(j:j)
    if (ch .ne. " ") then
      out_str = trim(out_str) // ch
    endif
    sweep_blanks = out_str
  enddo

END FUNCTION sweep_blanks

SUBROUTINE MPI_ALL_REDUCE_DOUBLE( vec_result , ndim )

  implicit none
  INCLUDE 'mpif.h'

  ! global
  integer :: ndim
  double precision, dimension ( ndim ) :: vec_result 
  ! local 
  integer :: ierr
  double precision, dimension ( : ) , allocatable :: vec_sum

  allocate ( vec_sum ( ndim ) )
  vec_sum=0.0d0

  CALL MPI_ALLREDUCE( vec_result , vec_sum , ndim , MPI_DOUBLE_PRECISION , MPI_SUM , MPI_COMM_WORLD , ierr )
  vec_result = vec_sum 

  deallocate ( vec_sum ) 

  return

END SUBROUTINE MPI_ALL_REDUCE_DOUBLE

SUBROUTINE MPI_ALL_REDUCE_INTEGER( vec_result , ndim )

  implicit none
  INCLUDE 'mpif.h'

  ! global
  integer :: ndim
  integer , dimension ( ndim ) :: vec_result
  ! local 
  integer :: ierr
  integer , dimension ( : ) , allocatable :: vec_sum

  allocate ( vec_sum ( ndim ) )
  vec_sum=0.0d0

  CALL MPI_ALLREDUCE( vec_result , vec_sum , ndim , MPI_INTEGER , MPI_SUM , MPI_COMM_WORLD , ierr )
  vec_result = vec_sum

  deallocate ( vec_sum )

  return

END SUBROUTINE MPI_ALL_REDUCE_INTEGER

SUBROUTINE MPI_ALL_REDUCE_DOUBLE_SCALAR ( sresult )

  implicit none
  INCLUDE 'mpif.h'

  double precision, intent (inout) :: sresult
  ! local 
  integer          :: ierr
  double precision :: ssum

  ssum=0.0d0
  CALL MPI_ALLREDUCE( sresult , ssum , 1 , MPI_DOUBLE_PRECISION , MPI_SUM , MPI_COMM_WORLD , ierr )
  sresult = ssum

  return

END SUBROUTINE
! ===== fmV =====
