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

!*********************** SUBROUTINE read_pos **********************************
!
!  here we read configuration (pos,vel,force) from file.
!
!******************************************************************************

SUBROUTINE read_pos 
  
  USE control,                  ONLY :  calc
  USE config,                   ONLY :  rx , ry , rz , vx , vy , vz , fx , fy , fz , atype , atypei , itype , &
                                        natmi , natm , dipia , qia , ipolar , rho , system , ntype , config_alloc , &
                                        simu_cell , config_print_info
  USE field,                    ONLY :  qch , dip , lpolar , field_init
  USE io_file,                  ONLY :  ionode , stdout , kunit_POSFF
  USE cell,                     ONLY :  lattice, periodicbc , dirkar

  implicit none

  ! local
  integer :: it , ia , i
  ! =====================================================
  !   type of positions coordinates 
  ! =====================================================
  logical :: allowed
  character(len=60), SAVE :: cpos
  character(len=60), SAVE :: cpos_allowed(4)
  data cpos_allowed / 'Direct' , 'D' , 'Cartesian' , 'C' /


  IF ( ionode ) then
    WRITE ( stdout ,'(a)')            '=============================================================' 
    WRITE ( stdout ,'(a)')            ''
    WRITE ( stdout ,'(a)')            'config from file POSFF'
  endif

  OPEN ( kunit_POSFF , file = 'POSFF' ) 
  READ ( kunit_POSFF ,* ) natm 
  READ ( kunit_POSFF ,* ) system
  READ ( kunit_POSFF ,* ) simu_cell%A ( 1 , 1 ) , simu_cell%A ( 2 , 1 ) , simu_cell%A ( 3 , 1 )
  READ ( kunit_POSFF ,* ) simu_cell%A ( 1 , 2 ) , simu_cell%A ( 2 , 2 ) , simu_cell%A ( 3 , 2 )
  READ ( kunit_POSFF ,* ) simu_cell%A ( 1 , 3 ) , simu_cell%A ( 2 , 3 ) , simu_cell%A ( 3 , 3 )
  READ ( kunit_POSFF ,* ) ntype
  READ ( kunit_POSFF ,* ) ( atypei ( it ) , it = 1 , ntype )
  IF ( ionode ) WRITE ( stdout      ,'(A,20A3)' ) 'found type information on POSFF : ', atypei ( 1:ntype )
  READ( kunit_POSFF ,*) ( natmi ( it ) , it = 1 , ntype )
  READ( kunit_POSFF ,*) cpos
  ! ======
  !  cpos
  ! ======
  do i = 1 , size( cpos_allowed )
   if ( trim(cpos) .eq. cpos_allowed(i))  allowed = .true.
  enddo
  if ( .not. allowed ) then
    if ( ionode )  WRITE ( stdout , '(a)' ) 'ERROR in POSFF at line 9 should be ', cpos_allowed
    STOP
  endif
 
  CALL lattice ( simu_cell )
  rho = DBLE ( natm ) / simu_cell%omega

  CALL config_alloc 

  ! ===============================
  ! config info
  ! ===============================
  CALL config_print_info(stdout)
  ! ===============================
  ! init force_field 
  ! ===============================
  CALL field_init

  ! =========================================      
  ! read positions and velocities from disk 
  ! =========================================      
  READ  ( kunit_POSFF , * ) ( atype ( ia ) , rx ( ia ) , ry ( ia ) , rz ( ia ) , &
                                             vx ( ia ) , vy ( ia ) , vz ( ia ) , &
                                             fx ( ia ) , fy ( ia ) , fz ( ia ) , ia = 1 , natm )

  if ( cpos .eq. 'Direct' ) then
    ! ======================================
    !         direct to cartesian
    ! ======================================
    CALL dirkar ( natm , rx , ry , rz , simu_cell%A )
    if ( ionode ) WRITE ( stdout      ,'(A,20A3)' ) 'atomic positions in direct coordinates in POSFF'
  else if ( cpos .eq. 'Cartesian' ) then
    if ( ionode ) WRITE ( stdout      ,'(A,20A3)' ) 'atomic positions in cartesian coordinates in POSFF'
  endif 


  CLOSE ( kunit_POSFF )

  CALL typeinfo_init

  CALL distance_tab ( stdout )

  CALL periodicbc ( natm , rx , ry , rz , simu_cell )

  return

END SUBROUTINE read_pos

!*********************** SUBROUTINE typeinfo_init *****************************
!
!
!******************************************************************************

SUBROUTINE typeinfo_init

  USE config ,  ONLY : atype , atypei , itype , natmi , natm , ntype , dipia , qia , quadia , ipolar , polia  
  USE field ,   ONLY : qch , quad_efg , dip , lpolar , pol
  
  implicit none


  ! local  
  integer :: ccs , cc
  integer :: ia , it, i , j 

  ! ==========================
  !  set some type parameters
  ! ==========================      
  natmi ( 0 ) = 0
  cc = 0
  do it = 1 , ntype
    ccs = cc
    cc = cc + natmi ( it )
    do ia = ccs + 1 , cc
      atype ( ia )     = atypei ( it )
      itype ( ia )     = it
      qia   ( ia )     = qch ( it )
      quadia( ia )     = quad_efg ( it )   
      dipia ( ia , 1 ) = dip ( it , 1 )
      dipia ( ia , 2 ) = dip ( it , 2 )
      dipia ( ia , 3 ) = dip ( it , 3 )
      ipolar ( ia )    = lpolar ( it )
      do i = 1 , 3
        do j = 1 , 3
          polia ( ia , i , j ) = pol ( it , i , j )
        enddo
      enddo 
    enddo
  enddo
  natmi  ( 0 ) = natm
  atypei ( 0 ) = 'ALL'

  return

END SUBROUTINE typeinfo_init



