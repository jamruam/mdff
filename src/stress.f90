! ===== fmV =====
!*********************** SUBROUTINE stress **********************************
!
! this subroutine calculates the stress tensor (Cauchy tensor probably)
! not the virial stress tensor Zhou et al.(Proc. R. soc. lond. A (2003) 459,
! 2347-2392)
!
! atom decomposition version
! 
!******************************************************************************

SUBROUTINE stress_bmlj ( iastart , iaend )!, list , point )

  USE control,  ONLY :  lvnlist
  USE config,   ONLY :  natm , fx , fy , fz , rx , ry , rz , box , omega , atype , itype , list , point
  USE field,    ONLY :  rcutsq ,  sigsq , epsp , fc , pp , qq
  USE io_file,  ONLY :  ionode , stdout , kunit_OUTFF

  implicit none
  INCLUDE 'mpif.h'

  ! global
  integer, intent(in)  :: iastart , iaend
!  integer, intent(in) :: list( 250 * natm ) , point( natm + 1 )

  !local 
  integer :: i , j , ia , ja , j1 ,jb , je , ierr 
  integer :: p1, p2
  integer :: nxij , nyij , nzij
  double precision :: rijsq , invbox
  double precision :: rxi , ryi , rzi
  double precision :: rxij, ryij , rzij , sr2 , srp , srq
  double precision :: wij , fxij , fyij , fzij
  double precision, dimension(3,3) :: tau
  double precision, dimension(3)   :: taux , tauxx
  double precision :: ptwo(2,2),qtwo(2,2)

  do j = 1,2
    do i = 1,2
      ptwo(i,j) = pp(i,j) * 0.5d0
      qtwo(i,j) = qq(i,j) * 0.5d0
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
        if( ja .ne. ia ) then
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
          if( rijsq .lt. rcutsq(p1,p2) ) then
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
        if( ja .gt. ia ) then
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
          if( rijsq .lt. rcutsq(p1,p2) ) then
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

  taux = 0.0d0
  tauxx = tau(1,:)
  CALL MPI_ALLREDUCE(tauxx,taux,3,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
  tau(1,:) = taux

  taux = 0.0d0
  tauxx = tau(2,:)
  CALL MPI_ALLREDUCE(tauxx,taux,3,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
  tau(2,:) = taux

  taux = 0.0d0
  tauxx = tau(3,:)
  CALL MPI_ALLREDUCE(tauxx,taux,3,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
  tau(3,:) = taux

  if ( ionode ) then
    WRITE ( stdout ,'(a)' )               'LJ : '
    WRITE ( stdout ,'(a,2a15)' )          '                 x' , 'y' , 'z'
    WRITE ( stdout ,'(a,3f15.8)')         'x            ' , tau(1,1) , tau(1,2) , tau(1,3)
    WRITE ( stdout ,'(a,3f15.8)')         'y            ' , tau(2,1) , tau(2,2) , tau(2,3)
    WRITE ( stdout ,'(a,3f15.8)')         'z            ' , tau(3,1) , tau(3,2) , tau(3,3)
    WRITE ( stdout ,'(a,f15.8,a6,f15.8,a1)')    'Press. iso = ', (tau(1,1) + tau(2,2) + tau(3,3)) / 3.0d0 , '(',(tau(1,1) + tau(2,2) + tau(3,3)) / 3.0d0 / dble(natm),')'
    WRITE ( stdout ,'(a)')                ''
    WRITE ( kunit_OUTFF ,'(a)' )          'LJ : '
    WRITE ( kunit_OUTFF ,'(a,2a15)' )     '                 x' , 'y' , 'z'
    WRITE ( kunit_OUTFF ,'(a,3f15.8)')    'x            ' , tau(1,1) , tau(1,2) , tau(1,3)
    WRITE ( kunit_OUTFF ,'(a,3f15.8)')    'y            ' , tau(2,1) , tau(2,2) , tau(2,3)
    WRITE ( kunit_OUTFF ,'(a,3f15.8)')    'z            ' , tau(3,1) , tau(3,2) , tau(3,3)
    WRITE ( kunit_OUTFF ,'(a,f15.8,a6,f15.8,a1)')    'Press. iso = ', (tau(1,1) + tau(2,2) + tau(3,3)) / 3.0d0 , '(',(tau(1,1) + tau(2,2) + tau(3,3)) / 3.0d0 / dble(natm),')'
    WRITE ( kunit_OUTFF ,'(a)')                ''
  endif

  return

END SUBROUTINE stress_bmlj


!*********************** SUBROUTINE stress_coulomb ****************************
!
!
!******************************************************************************
SUBROUTINE stress_coul

  USE control,          ONLY :  longrange
  USE config,           ONLY :  qit , qia , ntype , natm , omega , rx , ry , rz , box 
  USE io_file,          ONLY :  ionode , stdout , kunit_OUTFF
  USE constants,        ONLY :  imag , tpi , piroot 
  USE kspace,           ONLY :  struc_fact     
  USE field,            ONLY :  alphaES , km_coul
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
  double precision :: kri , str
  double complex   :: rhon , carg 
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


! ==============================================
!        direct space part
! ==============================================

   do ia = 1 , natm 
     rxi = rx(ia)
     ryi = ry(ia)
     rzi = rz(ia)
     qi = qia(ia)
     do ja = 1, natm

       if(ja .gt. ia ) then
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
         qijf = qij / rijsq      

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
  ! sum on atoms for one givn ik
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
      if ( km_coul%kptk(ik) .eq. 0 ) then 
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
        rhon = rhon + qit(it) * CONJG( km_coul%strf ( ik , it ) )
      enddo
      kri = ( kx * rxi + ky * ryi + kz * rzi ) 
      carg = EXP ( imag * kri )
      wij = rhon * carg * ak * imag 

      str = rhon * CONJG(rhon) * ak 

      tau_rec(1,1) = tau_rec(1,1) + ( 1.0d0 - kcoe * kx * kx ) * str
      tau_rec(1,2) = tau_rec(1,2) +           kcoe * kx * ky   * str
      tau_rec(1,3) = tau_rec(1,3) +           kcoe * kx * kz   * str
      tau_rec(2,1) = tau_rec(2,1) +           kcoe * ky * kx   * str
      tau_rec(2,2) = tau_rec(2,2) + ( 1.0d0 - kcoe * ky * ky ) * str
      tau_rec(2,3) = tau_rec(2,3) +           kcoe * ky * kz   * str
      tau_rec(3,1) = tau_rec(3,1) +           kcoe * kz * kx   * str
      tau_rec(3,2) = tau_rec(3,2) +           kcoe * kz * ky   * str
      tau_rec(3,3) = tau_rec(3,3) + ( 1.0d0 - kcoe * kz * kz ) * str

    enddo kpoint

  enddo 
  tau_rec = tau_rec * tpi / omega

  tau = tau_dir + tau_rec

   if ( ionode ) then
    WRITE ( stdout      ,'(a)'       )    'Coulomb : '
    WRITE ( stdout      ,'(a,2a15)'  )    '                 x' , 'y' , 'z' 
    WRITE ( stdout      ,'(a,3f15.8)')    'x            ' , tau(1,1) , tau(1,2) , tau(1,3)
    WRITE ( stdout      ,'(a,3f15.8)')    'y            ' , tau(2,1) , tau(2,2) , tau(2,3)
    WRITE ( stdout      ,'(a,3f15.8)')    'z            ' , tau(3,1) , tau(3,2) , tau(3,3)
    WRITE ( stdout ,'(a,f15.8,a6,f15.8,a1)')    'Press. iso = ', (tau(1,1) + tau(2,2) + tau(3,3)) / 3.0d0 , '(',(tau(1,1) + tau(2,2) + tau(3,3)) / 3.0d0 / dble(natm),')'
    WRITE ( stdout      ,'(a)'       )    ''
    WRITE ( kunit_OUTFF ,'(a)'       )    'Coulomb : '
    WRITE ( kunit_OUTFF ,'(a,2a15)'  )    '                 x' , 'y' , 'z'
    WRITE ( kunit_OUTFF ,'(a,3f15.8)')    'x            ' , tau(1,1) , tau(1,2) , tau(1,3)
    WRITE ( kunit_OUTFF ,'(a,3f15.8)')    'y            ' , tau(2,1) , tau(2,2) , tau(2,3)
    WRITE ( kunit_OUTFF ,'(a,3f15.8)')    'z            ' , tau(3,1) , tau(3,2) , tau(3,3)
    WRITE ( kunit_OUTFF ,'(a,f15.8,a6,f15.8,a1)')    'Press. iso = ', (tau(1,1) + tau(2,2) + tau(3,3)) / 3.0d0 , '(',(tau(1,1) + tau(2,2) + tau(3,3)) / 3.0d0 / dble(natm),')'
    WRITE ( kunit_OUTFF ,'(a)'       )    ''
  endif

  ttt3 = MPI_WTIME(ierr)
  
  return


END SUBROUTINE stress_coul
! ===== fmV =====
