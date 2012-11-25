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

!*********************** SUBROUTINE stress **********************************
!
! this subroutine calculates the stress tensor (Cauchy tensor probably)
! not the virial stress tensor Zhou et al.(Proc. R. soc. lond. A (2003) 459,
! 2347-2392)
!
! atom decomposition version
! 
!******************************************************************************

SUBROUTINE stress_bmlj ( iastart , iaend )

  USE control,  ONLY :  lvnlist
  USE config,   ONLY :  natm , fx , fy , fz , rx , ry , rz , box , omega , atype , itype , list , point , ntype
  USE field,    ONLY :  rcutsq ,  sigsq , epsp , fc , plj , qlj
  USE io_file,  ONLY :  ionode , stdout , kunit_OUTFF

  implicit none
  INCLUDE 'mpif.h'

  ! global
  integer, intent(in)  :: iastart , iaend

  !local 
  integer :: i , j , ia , ja , j1 ,jb , je , ierr 
  integer :: p1, p2
  integer :: nxij , nyij , nzij
  double precision :: rijsq , invbox
  double precision :: rxi , ryi , rzi
  double precision :: rxij, ryij , rzij , sr2 , srp , srq
  double precision :: wij , fxij , fyij , fzij
  double precision, dimension ( 3 , 3 ) :: tau
  double precision, dimension ( 3 )   :: taux , tauxx
  double precision :: ptwo ( ntype , ntype ) , qtwo ( ntype , ntype )

  do j = 1 , ntype 
    do i = 1 , ntype
      ptwo(i,j) = plj(i,j) * 0.5d0
      qtwo(i,j) = qlj(i,j) * 0.5d0
    enddo
  enddo
 
  invbox = 1.0d0/box 

  ! =======================
  !  with verlet list
  ! =======================
  if ( lvnlist ) then
    tau = 0.0d0
    do ia = iastart , iaend
      rxi = rx ( ia )
      ryi = ry ( ia )
      rzi = rz ( ia )
      jb = point( ia )
      je = point( ia + 1 ) - 1
      do j1 = jb , je
        ja = list ( j1 )
        if ( ja .ne. ia ) then
          rxij = rxi - rx ( ja )
          ryij = ryi - ry ( ja )
          rzij = rzi - rz ( ja )
          nxij = nint( rxij * invbox )
          nyij = nint( ryij * invbox )
          nzij = nint( rzij * invbox )
          rxij = rxij - box * nxij
          ryij = ryij - box * nyij
          rzij = rzij - box * nzij
          rijsq = rxij * rxij + ryij * ryij + rzij * rzij
          p1 = itype ( ja )
          p2 = itype ( ia )
          if ( rijsq .lt. rcutsq(p1,p2) ) then
            sr2 = sigsq(p1,p2) / rijsq
            srp = sr2 ** ( ptwo(p1,p2) )
            srq = sr2 ** ( qtwo(p1,p2) )
            wij = fc(p1,p2) * (srq-srp) * sr2
            fxij = wij * rxij
            fyij = wij * ryij
            fzij = wij * rzij
            tau(1,1) = tau(1,1) + rxij * fxij 
            tau(1,2) = tau(1,2) + rxij * fyij
            tau(1,3) = tau(1,3) + rxij * fzij
            tau(2,1) = tau(2,1) + ryij * fxij
            tau(2,2) = tau(2,2) + ryij * fyij
            tau(2,3) = tau(2,3) + ryij * fzij 
            tau(3,1) = tau(3,1) + rzij * fxij
            tau(3,2) = tau(3,2) + rzij * fyij 
            tau(3,3) = tau(3,3) + rzij * fzij
          endif
        endif
      enddo
    enddo
  tau = tau / omega
  ! =================
  !  no verlet list
  ! =================
  else
    tau = 0.0d0    
    do ia = iastart , iaend
      rxi = rx ( ia )
      ryi = ry ( ia )
      rzi = rz ( ia )
      do ja = 1, natm
        if ( ja .gt. ia ) then
          rxij = rxi - rx ( ja )
          ryij = ryi - ry ( ja )
          rzij = rzi - rz ( ja )
          nxij = nint( rxij * invbox )
          nyij = nint( ryij * invbox )
          nzij = nint( rzij * invbox )
          rxij = rxij - box * nxij
          ryij = ryij - box * nyij
          rzij = rzij - box * nzij
          rijsq = rxij  * rxij + ryij * ryij + rzij * rzij
          p1 = itype (j)
          p2 = itype (i)
          if ( rijsq .lt. rcutsq(p1,p2) ) then
            sr2 = sigsq(p1,p2) / rijsq
            srp = sr2 ** ( ptwo(p1,p2) )
            srq = sr2 ** ( qtwo(p1,p2) )
            wij = fc(p1,p2) * (srq-srp) * sr2
            fxij = wij * rxij
            fyij = wij * ryij
            fzij = wij * rzij
            tau(1,1) = tau(1,1) + rxij * fxij 
            tau(1,2) = tau(1,2) + rxij * fyij
            tau(1,3) = tau(1,3) + rxij * fzij
            tau(2,1) = tau(2,1) + ryij * fxij
            tau(2,2) = tau(2,2) + ryij * fyij
            tau(2,3) = tau(2,3) + ryij * fzij 
            tau(3,1) = tau(3,1) + rzij * fxij
            tau(3,2) = tau(3,2) + rzij * fyij 
            tau(3,3) = tau(3,3) + rzij * fzij
          endif
        endif
      enddo
    enddo
  tau = tau / omega
  endif

  do i = 1 , 3
    CALL MPI_ALL_REDUCE_DOUBLE ( tau ( i , : ) , 3  )
  enddo

  if ( ionode ) then
    WRITE ( stdout ,'(a)' )               'LJ : '
    WRITE ( stdout ,'(a,2a15)' )          '                 x' , 'y' , 'z'
    WRITE ( stdout ,'(a,3f15.8)')         'x            ' , tau(1,1) , tau(1,2) , tau(1,3)
    WRITE ( stdout ,'(a,3f15.8)')         'y            ' , tau(2,1) , tau(2,2) , tau(2,3)
    WRITE ( stdout ,'(a,3f15.8)')         'z            ' , tau(3,1) , tau(3,2) , tau(3,3)
    WRITE ( stdout ,'(a,f15.8,a6,f15.8,a1)') &
                                          'Press. iso = ', (tau(1,1) + tau(2,2) + tau(3,3)) / 3.0d0 , &
                                                       '(',(tau(1,1) + tau(2,2) + tau(3,3)) / 3.0d0 / dble(natm),')'
    WRITE ( stdout ,'(a)')                ''
    WRITE ( kunit_OUTFF ,'(a)' )          'LJ : '
    WRITE ( kunit_OUTFF ,'(a,2a15)' )     '                 x' , 'y' , 'z'
    WRITE ( kunit_OUTFF ,'(a,3f15.8)')    'x            ' , tau(1,1) , tau(1,2) , tau(1,3)
    WRITE ( kunit_OUTFF ,'(a,3f15.8)')    'y            ' , tau(2,1) , tau(2,2) , tau(2,3)
    WRITE ( kunit_OUTFF ,'(a,3f15.8)')    'z            ' , tau(3,1) , tau(3,2) , tau(3,3)
    WRITE ( kunit_OUTFF ,'(a,f15.8,a6,f15.8,a1)')  &
                                          'Press. iso = ', (tau(1,1) + tau(2,2) + tau(3,3)) / 3.0d0 , &
                                                       '(',(tau(1,1) + tau(2,2) + tau(3,3)) / 3.0d0 / dble(natm),')'
    WRITE ( kunit_OUTFF ,'(a)')                ''
  endif

  return

END SUBROUTINE stress_bmlj


!*********************** SUBROUTINE stress_coul_ewald *************************
!
!
!******************************************************************************
SUBROUTINE stress_coul_ewald

  USE control,          ONLY :  longrange
  USE config,           ONLY :  qia , ntype , natm , omega , rx , ry , rz , box 
  USE io_file,          ONLY :  ionode , stdout , kunit_OUTFF
  USE constants,        ONLY :  imag , tpi , piroot 
  USE kspace,           ONLY :  struc_fact     
  USE field,            ONLY :  alphaES , km_coul , qch
  USE time

  implicit none

  INCLUDE 'mpif.h'

  ! local
  integer :: ia, ja , it , ierr
  integer :: ik
  integer :: nxij , nyij , nzij
  double precision :: rij , rijsq , expon , invbox
  double precision :: alpha2
  double precision :: qi , qj , qij , qijf
  double precision :: rxi , ryi , rzi , rxij , ryij , rzij 
  double precision :: wij0 , wij , fxij , fyij , fzij 
  double precision :: ak, kx, ky, kz, kk , kcoe
  double complex   :: str
  double complex   :: rhon 
  double precision, dimension(3,3) :: tau , tau_dir , tau_rec
  double precision, external :: errfc 
  double precision :: ttt1 , ttt2 , ttt3 


  if (longrange .eq. 'direct' ) return

  ttt1 = MPI_WTIME(ierr)

  ! =====================
  ! facteur de structure 
  ! =====================
  CALL struc_fact ( km_coul )  

  ! =================
  !  some constants 
  ! =================
  invbox = 1.0d0 / box
  ! related ewald parameter
  alpha2 = alphaES * alphaES      
  tau_dir = 0.0d0
  tau_rec = 0.0d0


! ==============================================
!        direct space part
! ==============================================

   do ia = 1 , natm 

     rxi = rx(ia)
     ryi = ry(ia)
     rzi = rz(ia)
     qi = qia(ia)

     do ja = 1, natm

       if ( ja .ne. ia ) then

         qj = qia(ja)

         rxij = rxi - rx(ja)
         ryij = ryi - ry(ja)
         rzij = rzi - rz(ja)
         nxij = nint( rxij * invbox )
         nyij = nint( ryij * invbox )
         nzij = nint( rzij * invbox )
         rxij = rxij - box * nxij
         ryij = ryij - box * nyij
         rzij = rzij - box * nzij
         rijsq = rxij * rxij + ryij * ryij + rzij * rzij

         rij = dsqrt( rijsq )

         qij  = qi * qj / rij
         qijf = qij / rijsq      ! qi * qj / rij^3

         wij0 = errfc( alphaES * rij ) 
         expon = dexp( - alpha2 * rijsq ) / piroot
         wij  = qijf * ( wij0 + 2.0d0 * rij * alphaES * expon )

         fxij = wij * rxij
         fyij = wij * ryij
         fzij = wij * rzij

         tau_dir(1,1) = tau_dir(1,1) + rxij * fxij 
         tau_dir(1,2) = tau_dir(1,2) + rxij * fyij
         tau_dir(1,3) = tau_dir(1,3) + rxij * fzij
         tau_dir(2,1) = tau_dir(2,1) + ryij * fxij
         tau_dir(2,2) = tau_dir(2,2) + ryij * fyij
         tau_dir(2,3) = tau_dir(2,3) + ryij * fzij 
         tau_dir(3,1) = tau_dir(3,1) + rzij * fxij
         tau_dir(3,2) = tau_dir(3,2) + rzij * fyij 
         tau_dir(3,3) = tau_dir(3,3) + rzij * fzij

       endif
     enddo
  enddo

  ttt2 = MPI_WTIME(ierr)

! ==============================================
!            reciprocal space part
! ==============================================

  do ia = 1 , natm
    rxi = rx(ia)
    ryi = ry(ia)
    rzi = rz(ia)

    kpoint : do ik = 1, km_coul%nkcut 
      ! =================
      !   k-space  
      ! =================
      kx = km_coul%kpt(1,ik)
      ky = km_coul%kpt(2,ik)
      kz = km_coul%kpt(3,ik)
      kk = km_coul%kptk(ik)
      ak = dexp( - kk * 0.25d0 / alpha2 ) / kk
      kcoe     = 2.0d0 * ( 0.25d0 / alpha2  + 1.0d0 / kk )

      if ( km_coul%kptk(ik) .eq. 0.0d0 ) then 
        WRITE ( stdout , * ) 'the sum should be done on k! =  0',ik
        STOP 
      endif

      ! ===============================
      !                              ---
      !  charge density in k-space ( \   q * facteur de structure  )
      !                              /__
      ! ===============================
      rhon = (0.d0, 0.d0)
      do it = 1, ntype
        rhon = rhon + qch(it) * CONJG( km_coul%strf ( ik , it ) )
      enddo
      str = rhon * CONJG(rhon) * ak 
      !print*,ak,kk,- kk * 0.25d0 / alpha2,dexp( - kk * 0.25d0 / alpha2 )

      tau_rec(1,1) = tau_rec(1,1) + ( 1.0d0 - kcoe * kx * kx ) * str
      tau_rec(1,2) = tau_rec(1,2) -           kcoe * kx * ky   * str
      tau_rec(1,3) = tau_rec(1,3) -           kcoe * kx * kz   * str
      tau_rec(2,1) = tau_rec(2,1) -           kcoe * ky * kx   * str
      tau_rec(2,2) = tau_rec(2,2) + ( 1.0d0 - kcoe * ky * ky ) * str
      tau_rec(2,3) = tau_rec(2,3) -           kcoe * ky * kz   * str
      tau_rec(3,1) = tau_rec(3,1) -           kcoe * kz * kx   * str
      tau_rec(3,2) = tau_rec(3,2) -           kcoe * kz * ky   * str
      tau_rec(3,3) = tau_rec(3,3) + ( 1.0d0 - kcoe * kz * kz ) * str

    enddo kpoint

  enddo

  ! ======================================================
  ! remark on the unit :
  ! 1/(4*pi*epislon_0) = 1 => epsilon_0 = 1/4pi
  ! ======================================================
  print*,tau_rec
  tau_rec = tau_rec * tpi / omega 
  print*,tau_rec

  tau = tau_dir + tau_rec

   if ( ionode ) then
    WRITE ( stdout      ,'(a)'       )    'Coulomb : '
    WRITE ( stdout      ,'(a,f15.8)'       )    'here !!!!! : ',tpi / omega
    WRITE ( stdout      ,'(a,2a15)'  )    '                 x' , 'y' , 'z' 
    WRITE ( stdout      ,'(a,3f15.8)')    'x            ' , tau(1,1) , tau(1,2) , tau(1,3)
    WRITE ( stdout      ,'(a,3f15.8)')    'y            ' , tau(2,1) , tau(2,2) , tau(2,3)
    WRITE ( stdout      ,'(a,3f15.8)')    'z            ' , tau(3,1) , tau(3,2) , tau(3,3)
    WRITE ( stdout ,'(a,f15.8,a6,f15.8,a1)') &
                                          'Press. iso = ', (tau(1,1) + tau(2,2) + tau(3,3)) / 3.0d0 , &
                                                       '(',(tau(1,1) + tau(2,2) + tau(3,3)) / 3.0d0 / dble(natm),')'
    WRITE ( stdout      ,'(a)'       )    ''
    WRITE ( kunit_OUTFF ,'(a)'       )    'Coulomb : '
    WRITE ( kunit_OUTFF ,'(a,2a15)'  )    '                 x' , 'y' , 'z'
    WRITE ( kunit_OUTFF ,'(a,3f15.8)')    'x            ' , tau(1,1) , tau(1,2) , tau(1,3)
    WRITE ( kunit_OUTFF ,'(a,3f15.8)')    'y            ' , tau(2,1) , tau(2,2) , tau(2,3)
    WRITE ( kunit_OUTFF ,'(a,3f15.8)')    'z            ' , tau(3,1) , tau(3,2) , tau(3,3)
    WRITE ( kunit_OUTFF ,'(a,f15.8,a6,f15.8,a1)') &
                                          'Press. iso = ', (tau(1,1) + tau(2,2) + tau(3,3)) / 3.0d0 , &
                                                       '(',(tau(1,1) + tau(2,2) + tau(3,3)) / 3.0d0 / dble(natm),')'
    WRITE ( kunit_OUTFF ,'(a)'       )    ''
  endif

  ttt3 = MPI_WTIME(ierr)
  
  return


END SUBROUTINE stress_coul_ewald

!*********************** SUBROUTINE stress_coul_direct *************************
!
!
!******************************************************************************
SUBROUTINE stress_coul_direct ( iastart , iaend )

  USE control,          ONLY :  longrange
  USE config,           ONLY :  qia , ntype , natm , omega , rx , ry , rz , box 
  USE io_file,          ONLY :  ionode , stdout , kunit_OUTFF
  USE constants,        ONLY :  imag , tpi , piroot 
  USE field,            ONLY :  rm_coul , qch
  USE time

  implicit none

  INCLUDE 'mpif.h'

  ! global
  integer, intent(in)  :: iastart , iaend

  ! local
  integer :: ia, ja , ierr
  integer :: ncell
  double precision :: rij , rijsq 
  double precision :: qi , qj , qij , qijf
  double precision :: rxi , ryi , rzi , rxj, ryj, rzj, rxij , ryij , rzij 
  double precision :: fxij, fyij , fzij 
  double precision, dimension(3,3) :: tau , tau_dir
  double precision, external :: errfc 
  double precision :: ttt1 , ttt2 , ttt3 


  if ( longrange .eq. 'ewald' ) return

  ttt1 = MPI_WTIME(ierr)
  tau_dir = 0.0d0

! =========================================================
!  MAIN LOOP calculate EFG(i) for each atom i parallelized
! =========================================================
atom : do ia = iastart , iaend

    rxi = rx(ia)
    ryi = ry(ia)
    rzi = rz(ia)
    qi  = qia(ia)

    ! ==================================================
    ! sum over neighboring cells (see direct_sum_init )
    ! ==================================================
    do ncell = 1 , rm_coul%ncmax
         ! ==============================
         !  ia and ja in different cells
         ! ==============================
         if ( rm_coul%lcell ( ncell ) .eq. 1) then
          do ja = 1 , natm

            qj    = qia(ja)

            rxj   = rx(ja) + rm_coul%boxxyz(1,ncell)
            ryj   = ry(ja) + rm_coul%boxxyz(2,ncell)
            rzj   = rz(ja) + rm_coul%boxxyz(3,ncell)

            rxij  = rxi - rxj
            ryij  = ryi - ryj
            rzij  = rzi - rzj

            rijsq = rxij * rxij + ryij * ryij + rzij * rzij

            rij   = dsqrt( rijsq )
            qij   = qi * qj / rij
            qijf  = qij / rijsq

            fxij = qijf * rxij
            fyij = qijf * ryij
            fzij = qijf * rzij

            tau_dir(1,1) = tau_dir(1,1) + rxij * fxij 
            tau_dir(1,2) = tau_dir(1,2) + rxij * fyij
            tau_dir(1,3) = tau_dir(1,3) + rxij * fzij
            tau_dir(2,1) = tau_dir(2,1) + ryij * fxij
            tau_dir(2,2) = tau_dir(2,2) + ryij * fyij
            tau_dir(2,3) = tau_dir(2,3) + ryij * fzij 
            tau_dir(3,1) = tau_dir(3,1) + rzij * fxij
            tau_dir(3,2) = tau_dir(3,2) + rzij * fyij 
            tau_dir(3,3) = tau_dir(3,3) + rzij * fzij

          enddo ! ja
 
        endif 
         ! =======================================
         !  ia and ja in the same cell (ia.ne.ja)
         ! =======================================
        if ( rm_coul%lcell(ncell) .eq. 0) then

          do ja = 1,natm

            if ( ja .ne. ia ) then

              qj = qia(ja)

              rxj  = rx(ja)
              ryj  = ry(ja)
              rzj  = rz(ja)

              rxij = rxi - rxj
              ryij = ryi - ryj
              rzij = rzi - rzj

              rijsq   = rxij * rxij + ryij * ryij + rzij * rzij
 
              rij = dsqrt( rijsq )
              qij = qi * qj / rij
              qijf = qij / rijsq

              fxij = qijf * rxij
              fyij = qijf * ryij
              fzij = qijf * rzij

              tau_dir(1,1) = tau_dir(1,1) + rxij * fxij 
              tau_dir(1,2) = tau_dir(1,2) + rxij * fyij
              tau_dir(1,3) = tau_dir(1,3) + rxij * fzij
              tau_dir(2,1) = tau_dir(2,1) + ryij * fxij
              tau_dir(2,2) = tau_dir(2,2) + ryij * fyij
              tau_dir(2,3) = tau_dir(2,3) + ryij * fzij 
              tau_dir(3,1) = tau_dir(3,1) + rzij * fxij
              tau_dir(3,2) = tau_dir(3,2) + rzij * fyij 
              tau_dir(3,3) = tau_dir(3,3) + rzij * fzij

            endif ! ia.ne.ja

          enddo ! ja

        endif 

     enddo ! ncell
  enddo atom

  ttt2 = MPI_WTIME(ierr)

  tau = tau_dir 

   if ( ionode ) then
    WRITE ( stdout      ,'(a)'       )    ' WARNING not complete test version'
    WRITE ( stdout      ,'(a)'       )    'Coulomb : '
    WRITE ( stdout      ,'(a,2a15)'  )    '                 x' , 'y' , 'z' 
    WRITE ( stdout      ,'(a,3f15.8)')    'x            ' , tau(1,1) , tau(1,2) , tau(1,3)
    WRITE ( stdout      ,'(a,3f15.8)')    'y            ' , tau(2,1) , tau(2,2) , tau(2,3)
    WRITE ( stdout      ,'(a,3f15.8)')    'z            ' , tau(3,1) , tau(3,2) , tau(3,3)
    WRITE ( stdout ,'(a,f15.8,a6,f15.8,a1)') &
                                          'Press. iso = ', (tau(1,1) + tau(2,2) + tau(3,3)) / 3.0d0 , &
                                                       '(',(tau(1,1) + tau(2,2) + tau(3,3)) / 3.0d0 / dble(natm),')'
    WRITE ( stdout      ,'(a)'       )    ''
    WRITE ( kunit_OUTFF ,'(a)'       )    ' WARNING not complete test version'
    WRITE ( kunit_OUTFF ,'(a)'       )    'Coulomb : '
    WRITE ( kunit_OUTFF ,'(a,2a15)'  )    '                 x' , 'y' , 'z'
    WRITE ( kunit_OUTFF ,'(a,3f15.8)')    'x            ' , tau(1,1) , tau(1,2) , tau(1,3)
    WRITE ( kunit_OUTFF ,'(a,3f15.8)')    'y            ' , tau(2,1) , tau(2,2) , tau(2,3)
    WRITE ( kunit_OUTFF ,'(a,3f15.8)')    'z            ' , tau(3,1) , tau(3,2) , tau(3,3)
    WRITE ( kunit_OUTFF ,'(a,f15.8,a6,f15.8,a1)') &
                                          'Press. iso = ', (tau(1,1) + tau(2,2) + tau(3,3)) / 3.0d0 , &
                                                       '(',(tau(1,1) + tau(2,2) + tau(3,3)) / 3.0d0 / dble(natm),')'
    WRITE ( kunit_OUTFF ,'(a)'       )    ''
  endif

  ttt3 = MPI_WTIME(ierr)
  
  return


END SUBROUTINE stress_coul_direct
! ===== fmV =====
