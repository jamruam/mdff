! ===== fmV =====
!*********************** SUBROUTINE periodicpbc *******************************
!
! replace atoms inside the box
! note:
! pbc are used in calculation of most properties but the positions are not
! effectively reajust to prevent big jump when using verlet list. 
! 
! input :
!          natm      : number of atoms
!          ralpha : positions component alpha
!          box    : size of the cubic box
!
!******************************************************************************

SUBROUTINE periodicbc ( natm , rx , ry , rz , box )

   implicit none

  ! global
  integer :: natm
  double precision :: box
  double precision :: rx ( natm ) , ry ( natm ) , rz ( natm )

  ! local
  integer :: ia

  do ia = 1 , natm
     rx ( ia ) = rx ( ia ) - nint( rx ( ia ) / box ) * box
     ry ( ia ) = ry ( ia ) - nint( ry ( ia ) / box ) * box
     rz ( ia ) = rz ( ia ) - nint( rz ( ia ) / box ) * box
  enddo

  return

END SUBROUTINE periodicbc

! The folloowing subroutines: init-fcc, init-sc and init-bcc should be merged in a more efficient way !
! TODO test LJ crystal and compare to litterature data

        
!*********************** SUBROUTINE init_fcc **********************************
!
! readapted from Allen-Tildsley  
! FICHE F.23.  ROUTINE TO SET UP ALPHA FCC LATTICE OF LINEAR MOLECULES
! This FORTRAN code is intended to illustrate points made in the text.
! To our knowledge it works correctly.  However it is the responsibility of
! the   USEr to test it, if it is to be   USEd in a research application.
!
! SETS UP THE ALPHA FCC LATTICE FOR N LINEAR MOLECULES.         
!                                                                   
! THE SIMULATION BOX IS A UNIT CUBE CENTRED AT THE ORIGIN.      
! N SHOULD BE AN INTEGER OF THE FORM ( 4 * ( NC ** 3 ) ),       
! WHERE NC IS THE NUMBER OF FCC UNIT CELLS IN EACH DIRECTION.   
! SEE FIGURE 5.10 FOR A DIAGRAM OF THE LATTICE AND A            
! DEFINITION OF THE FOUR ORIENTATIONAL SUBLATTICES.             
!                                                                   
! PRINCIPAL VARIABLES:                                          
!                                                                   
! INTEGER N                    NUMBER OF MOLECULES              
! REAL    RX(N),RY(N),RZ(N)    MOLECULAR POSITIONS              
! REAL    EX(N),EY(N),EZ(N)    UNIT VECTORS GIVING ORIENTATIONS 
! REAL    RROOT3               1.0 / SQRT ( 3.0 )              
! 
!******************************************************************************

SUBROUTINE init_fcc ( key , natm , ntype , ncell , rx , ry , rz , box , struct )

  USE io_file,  ONLY :  ionode , stdout, kunit_OUTFF


  ! global
  integer :: natm , ncell , ntype
  double precision , dimension ( natm ) :: rx , ry , rz       ! positions
  double precision :: box
  character*6      :: key
  character*60     :: struct     ! structure type for fcc and 2 types ('NaCl' or 'random')

  ! local
  double precision :: rroot3
  double precision :: cell, cell2, cellp
  integer :: i, ix, iy, iz, iref, m  

        if( ionode .and. ( key .eq. 'AAAAAA' ) ) then
           WRITE ( stdout ,'(a,i3,a1,i3,a1,i3,a)')                                  'Face-centered cubic structure :',ncell,'x',ncell,'x',ncell,'  '
           if( struct .eq. 'random' .and. ntype .eq. 2 ) WRITE ( stdout ,'(a)')     'random disposition of atoms                  '
           if( struct .eq. 'NaCl'   ) WRITE ( stdout ,'(a)')                        'NaCl: natm = 8*(ncell*ncell*ncell)           ' 
           WRITE ( kunit_OUTFF ,'(a,i3,a1,i3,a1,i3,a)')                             'Face-centered cubic structure :',ncell,'x',ncell,'x',ncell,'  '
           if( struct .eq. 'random' .and. ntype .eq. 2 ) WRITE ( kunit_OUTFF,'(a)') 'random disposition of atoms                  '
           if( struct .eq. 'NaCl'   ) WRITE ( kunit_OUTFF,'(a)')                    'NaCl: natm = 8*(ncell*ncell*ncell)           ' 
         endif

  ! ======================================
  !  CALCULATE THE SIDE OF THE UNIT CELL 
  ! ======================================

  cellp = box 

  rroot3 = 1.0d0 / SQRT ( 3.0d0 )
  cell  = 1.0d0 / dble( ncell )
  cell2 = 0.5d0 * cell 

  ! ======================
  !  BUILD THE UNIT CELL 
  ! ======================
  !  ** SUBLATTICE A **
  ! ======================

  rx(1) =  0.0d0 
  ry(1) =  0.0d0 
  rz(1) =  0.0d0 

  ! ======================
  !   ** SUBLATTICE B **
  ! ======================

  rx(2) =  cell2
  ry(2) =  cell2
  rz(2) =  0.0d0

  ! ======================
  !  ** SUBLATTICE C **
  ! ======================

  rx(3) =  0.0d0
  ry(3) =  cell2
  rz(3) =  cell2

  ! ======================
  !  ** SUBLATTiCE D **
  ! ======================

  rx(4) =  cell2
  ry(4) =  0.0d0
  rz(4) =  cell2

  ! ==================================================
  !   ** CONSTRUCT THE LATTiCE FROm THE UNiT cell **
  ! ==================================================

  m = 0

  do iz = 1, ncell
    do iy = 1, ncell
      do ix = 1, ncell
        do iref = 1, 4
          rx ( iref + m ) = rx ( iref ) + cell * DBLE ( ix - 1 ) 
          ry ( iref + m ) = ry ( iref ) + cell * DBLE ( iy - 1 ) 
          rz ( iref + m ) = rz ( iref ) + cell * DBLE ( iz - 1 ) 
        enddo
        m = m + 4
      enddo
    enddo
  enddo

  ! ============================================
  !   ** SHiFT CENTRE OF BOX TO THE ORiGiN **
  ! ============================================

  do i = 1, natm
    rx ( i ) = rx ( i ) - 0.5d0
    ry ( i ) = ry ( i ) - 0.5d0
    rz ( i ) = rz ( i ) - 0.5d0
    rx ( i ) = cellp * rx ( i )
    ry ( i ) = cellp * ry ( i )
    rz ( i ) = cellp * rz ( i ) 
  enddo

  return

END SUBROUTINE init_fcc


!*********************** SUBROUTINE init_sc ***********************************
!
! from frenkel and Smit
!
!******************************************************************************

SUBROUTINE init_sc

  USE config,   ONLY :  natm , ncell , rx , ry , rz , box 
  USE io_file,  ONLY :  ionode , stdout, kunit_OUTFF

  implicit none

  ! local
  integer ::  i , j , k , itel
  double precision :: del , dx , dy , dz

  if( ionode ) then
    WRITE ( stdout ,'(a)') 'simple cubic structure'
    WRITE ( kunit_OUTFF,'(a)') 'simple cubic structure'
  endif

  del = box / DBLE( ncell )
  itel = 0 
  dx = -del
  do i = 1 , ncell
    dx = dx + del
    dy = -del
    do j = 1 , ncell
      dy = dy + del
      dz = -del
      do k = 1 , ncell
        dz = dz + del
        if ( itel .lt. natm ) then
          itel = itel + 1
          rx ( itel ) = dx 
          ry ( itel ) = dy 
          rz ( itel ) = dz 
        endif
      enddo
    enddo
  enddo

  return

END SUBROUTINE init_sc

!*********************** SUBROUTINE init_bcc **********************************
!!TEST
! base centered cubic 
!!TEST
!******************************************************************************

SUBROUTINE init_bcc

  USE config,   ONLY :  natm , ncell , rx , ry , rz , box
  USE io_file,  ONLY :  ionode , stdout, kunit_OUTFF

  implicit none

  if( ionode ) then
    WRITE ( stdout ,'(a)')     'body centered cubic structure'
    WRITE ( kunit_OUTFF,'(a)') 'body centered cubic structure'
  endif

  


  return

END SUBROUTINE init_bcc

!*********************** SUBROUTINE random_config *****************************
!
! this subroutine attributes randomly the atom type to a given configuration
! when the starting configuration is completely generated by the code 
!
!******************************************************************************

SUBROUTINE random_config ( natm , atype , xna )

  USE config,   ONLY :  rx, ry ,rz 

  implicit none

  ! global
  integer :: natm
  double precision :: xna 
  character*10 , dimension ( natm )  :: atype , atmp

  ! local
  integer :: iseed
  integer :: i, k, nna, nnb
  double precision :: X
  double precision, dimension (:) , allocatable :: xtmp , ytmp , ztmp
 
  allocate( xtmp(natm) , ytmp(natm) , ztmp(natm) ) 
 
  nna = int( xna * natm )
  nnb = natm - nna
  ! ==================
  !  initialize atype
  ! ==================
  do i = 1, natm
   atype(i) = 'C'
  enddo

  do i = 1 , nna
    CALL RANDOM_SEED(SIZE = iseed)
    CALL RANDOM_NUMBER(HARVEST = x)
    k = int( X * natm ) + 1
    do while ( atype(k) .eq. 'A' )
      CALL RANDOM_SEED(SIZE = iseed)
      CALL RANDOM_NUMBER(HARVEST = x)
      k = int( X * natm ) + 1
    enddo
    atype(k) = 'A'
  enddo

  do i = 1 , natm
    if( atype(i) .ne. 'A' ) then
      atype(i) = 'B'
    endif
  enddo

  ! ====================
  !  reorganized atoms   
  ! ====================
  k = 0
  do i = 1 , natm
   if ( atype(i) .eq. 'A' ) then
     k = k + 1 
     if ( k .gt. nna ) stop 
     atmp(k) = atype(i)
     xtmp(k) = rx(i)
     ytmp(k) = ry(i)
     ztmp(k) = rz(i)
   endif 
  enddo

  do i = 1 , natm
   if ( atype(i) .eq. 'B' ) then
     k = k + 1
     if ( k .gt. natm ) stop
     atmp(k) = atype(i)
     xtmp(k) = rx(i)
     ytmp(k) = ry(i)
     ztmp(k) = rz(i)
   endif
  enddo

  do i=1, natm
   atype(i) = atmp(i)
   rx(i) = xtmp(i)
   ry(i) = ytmp(i)
   rz(i) = ztmp(i)
  enddo
  
  deallocate( xtmp , ytmp , ztmp ) 


  return


END SUBROUTINE random_config



! ===== fmV =====
