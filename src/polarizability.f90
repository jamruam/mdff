MODULE polarizability

  USE constants,                ONLY :  dp
  implicit none

  real(kind=dp) :: efield (3,3)
  real(kind=dp) :: epsw



CONTAINS

! *********************** SUBROUTINE pola_init ***********************************
!> \brief
!! initialisation of main parameters
! ******************************************************************************
SUBROUTINE pola_init

  USE control,  ONLY :  lstatic , calc
  USE io,       ONLY :  ionode ,stdin, stdout

  implicit none

  integer            :: ioerr
  character(len=132) :: filename

  namelist /polatag/    efield,&
                        epsw  

  if ( calc .ne. 'polarizability' ) return
  ! ======================
  !  polatag default values
  ! ======================
  CALL pola_default_tag
  ! =====================
  !  read polatag namelist
  ! =====================
  CALL getarg (1,filename)
  OPEN ( stdin , file = filename)
  READ ( stdin , polatag ,iostat =ioerr)
  if ( ioerr .lt. 0 )  then
    io_node &
    WRITE ( stdout, '(a)') 'ERROR reading input_file : polatag section is absent'
    STOP
  elseif ( ioerr .gt. 0 )  then
    io_node &
    WRITE ( stdout, '(a,i8)') 'ERROR reading input_file : polatag wrong tag',ioerr
    STOP
  endif
  CLOSE  ( stdin )
  ! ======================
  !  check polatag namelist
  ! ======================
  CALL pola_check_tag

  ! ===================
  !  print polatag info
  ! ===================
  CALL pola_print_info(stdout)

  return

END SUBROUTINE pola_init


SUBROUTINE pola_main

  USE cell
  USE kspace 

  implicit none

  TYPE ( celltype )                            :: simu_cell          !< simulation cell
  real(kind=dp) :: dmu(16,3,3), T2(3,3)
  real(kind=dp) :: f(16,3,3), invf(16,3,3), pola(16,3,3)
  real(kind=dp) :: f_rec(16,3,3) , f_dir(16,3,3) , f_self(16,3,3)
  real(kind=dp) :: e(3,3)
  real(kind=dp) :: rx(16), ry(16), rz(16) 
  real(kind=dp) :: rxi, ryi, rzi 
  real(kind=dp) :: rxj, ryj, rzj 
  real(kind=dp) :: exij , eyij , ezij
  real(kind=dp) :: mutot(3,3) 
  real(kind=dp) :: rij(3) , sij(3)
  character(len=2) :: atype(16),atypei(16)
  character(len=60) :: system
  integer :: natmi(16) , nions , ntype , kmax
  integer :: i, j, k, ia , ja , it, ik, nk
  real(kind=dp) :: d, d2 ,d3 , d5, dm1, dm3, dm5, F0, F1 , F2 , expon
  real(kind=dp) :: alphaES , alpha2, alpha3, rcut , cutsq,recarg,k_dot_mu,muix,muiy,muiz,k_dot_r
  real(kind=dp) :: rhonk_R(3),rhonk_I(3) ,Ak, selfa, selfa2
  character(len=60) :: cpos

  integer, parameter :: LWORK=1000
  real(kind=dp) :: WORK ( LWORK )
  integer :: ipiv ( 3 )
  integer :: ierr
  integer(kind=4):: lda,ldb,lde
  TYPE ( kmesh )   :: km_coul                          !< kpoint mesh ( see kspace.f90 )
  real(kind=dp) :: kx,ky,kz,kk,cutlongrange
  real(kind=dp) :: ckr(16) , skr(16)
  real(kind=dp) :: eps,alpha,tol,tol1, epsw,fij(3),fpi_V,tpi_V
  integer :: kmax_1, kmax_2, kmax_3,kES(3), ii


  CALL read_pos

  CALL ewald_param

  km_coul%meshlabel='km_coul'
  km_coul%kmax(1) = kES(1)
  km_coul%kmax(2) = kES(2)
  km_coul%kmax(3) = kES(3)

!   with symmetry
  nk = km_coul%kmax(3) + km_coul%kmax(2) * ( 2 * km_coul%kmax(3) + 1 ) + km_coul%kmax(1) * ( 2 * km_coul%kmax(2) + 1 ) * ( 2 * km_coul%kmax(3) + 1 )
  km_coul%nk = nk
  allocate ( km_coul%kptk( nk ) , km_coul%kptx(nk), km_coul%kpty(nk), km_coul%kptz(nk) )
  allocate ( km_coul%Ak      ( nk ) )
  allocate ( km_coul%kcoe    ( nk ) )
  CALL kpoint_sum_init ( km_coul , alphaES )




  rcut =  0.5_dp * MIN(simu_cell%WA,simu_cell%WB,simu_cell%WC) + 1.0d0
  ! dl_poly like
  if ( min (cutlongrange,rcut) .ne. cutlongrange ) &
  WRITE ( 6 , '(a)') 'WARNING : cutlongrange will be changed according to simu_cell%W*'
  cutlongrange = min (cutlongrange,rcut)
  eps=min(abs(epsw),0.5_dp)
  tol=sqrt(abs(log(eps*cutlongrange)))
  alpha=sqrt(abs(log(eps*cutlongrange*tol)))/rcut
  tol1=sqrt(-log(eps*cutlongrange*(2.0_dp*tol*alpha)**2))
  kmax_1=nint(0.25_dp+simu_cell%ANORM(1)*alpha*tol1/pi)
  kmax_2=nint(0.25_dp+simu_cell%ANORM(2)*alpha*tol1/pi)
  kmax_3=nint(0.25_dp+simu_cell%ANORM(3)*alpha*tol1/pi)
 ! dl_poly like
  write ( 6 , '(a)' ) 'automatic Ewald Sum '
  alphaES = alpha
  kES(1)=kmax_1
  kES(2)=kmax_2
  kES(3)=kmax_3
  alpha2 = alphaES * alphaES
  alpha3 = alpha2  * alphaES
  cutsq = cutlongrange * cutlongrange 
  write(6,'(a,f12.6,a,f12.6,a,3i)') 'alpha = ',alphaES,' rcut = ',cutlongrange,' k = ',kES

  km_coul%meshlabel='km_coul'
  km_coul%kmax(1) = kES(1)
  km_coul%kmax(2) = kES(2)
  km_coul%kmax(3) = kES(3)
    ! full kpt
    nk = ( 2 * km_coul%kmax(1) + 1 ) * ( 2 * km_coul%kmax(2) + 1 ) * ( 2 * km_coul%kmax(3) + 1 )
!   with symmetry
!  nk = km_coul%kmax(3) + km_coul%kmax(2) * ( 2 * km_coul%kmax(3) + 1 ) + km_coul%kmax(1) * ( 2 * km_coul%kmax(2) + 1 ) * ( 2 * km_coul%kmax(3) + 1 )
  km_coul%nk = nk
  allocate ( km_coul%kptk( nk ) , km_coul%kptx(nk), km_coul%kpty(nk), km_coul%kptz(nk) )
  allocate ( km_coul%Ak      ( nk ) )
  CALL kpoint_sum_init ( km_coul , alphaES ,simu_cell )

  ! =========================================      
  ! read positions, velocities, forces from disk 
  ! =========================================      
  do ia=1,nions
  READ  ( 5000 , * ) atype ( ia ) , rx ( ia ) , ry ( ia ) , rz ( ia )
  enddo

  CLOSE(5000)

  CALL distance_tab
  ! ======================================
  !         cartesian to direct 
  ! ======================================
  CALL kardir ( nions , rx , ry , rz , simu_cell%B )


  OPEN(UNIT=1000,FILE='dx')
    do ia=1,nions
      read(1000,*) dmu(ia,1,1),dmu(ia,2,1),dmu(ia,3,1)
    enddo 
  CLOSE(1000)
  OPEN(UNIT=1000,FILE='dy')
    do ia=1,nions
      read(1000,*) dmu(ia,1,2),dmu(ia,2,2),dmu(ia,3,2)
    enddo 
  CLOSE(1000)
  OPEN(UNIT=1000,FILE='dz')
    do ia=1,nions
      read(1000,*) dmu(ia,1,3),dmu(ia,2,3),dmu(ia,3,3)
    enddo 
  CLOSE(1000)


  f = 0.0_dp
  e(1,1) = 0.01_dp
  e(2,1) =  0.0_dp
  e(3,1) =  0.0_dp

  e(1,2) =  0.0_dp
  e(2,2) = 0.01_dp
  e(3,2) =  0.0_dp

  e(1,3) =  0.0_dp
  e(2,3) =  0.0_dp
  e(3,3) = 0.01_dp

  do ia=1,nions
    rxi = rx(ia)
    ryi = ry(ia)
    rzi = rz(ia)

      do ja=1 , nions

        if ( ja .le. ia ) cycle

        rxj  = rx(ja)
        ryj  = ry(ja)
        rzj  = rz(ja)
        rij(1) = rxi - rxj
        rij(2) = ryi - ryj
        rij(3) = rzi - rzj
        sij = rij - nint ( rij )
        rij=0.0_dp
        do j=1, 3
          do k=1, 3
            rij(j) = rij(j) + sij(k) * simu_cell%A(j,k)
          enddo
        enddo
        d2  = rij(1) * rij(1) + rij(2) * rij(2) + rij(3) * rij(3)
        if ( d2 .gt. cutsq ) cycle
  !      print*,ia,ja,d2,cutlongrange**2.0_dp

        d   = SQRT ( d2 )
        d3  = d2 * d
        d5  = d3 * d2
        dm1 = 1.0_dp / d
        dm3 = dm1 / d2
        dm5 = dm3 / d2

        expon = EXP ( - alpha2 * d2 )    / piroot
        F0    = errfc( alphaES * d )
        F1    = F0 + 2.0_dp * alphaES * d  * expon
        F2    = F1 + 4.0_dp * alpha3  * d3 * expon / 3.0_dp



        T2 = 0.0_dp
        do j = 1 , 3
          do k = 1 , 3
            if ( j .gt. k ) cycle
              T2(j,k) = ( 3.0_dp * rij(j) * rij(k) * F2 ) * dm5
            if ( j .eq. k ) T2(j,j) = T2 (j,j) - F1 * dm3
         enddo
       enddo
       T2(2,1) = T2(1,2)
       T2(3,1) = T2(1,3)
       T2(3,2) = T2(2,3)

       ! field on x,y,z
       do ii=1,3 
         do k = 1 , 3
           f_dir(ia,:,ii) = f_dir(ia,:,ii) +  T2(:,k) * dmu(ja,k,ii)
           f_dir(ja,:,ii) = f_dir(ja,:,ii) +  T2(:,k) * dmu(ia,k,ii)
         enddo
       enddo

    enddo    

  enddo

  do ia=1, nions
    do ii=1, 3
    f_dir(ia,:,ii) = f_dir(ia,:,ii) + e(:,ii) 
    enddo
  enddo

  
  f_rec=0.0_dp
  ! ==============================================
  !            reciprocal space part
  ! ==============================================
  kpoint : do ik = 1,  km_coul%nk

    if (km_coul%kptk(ik) .eq. 0.0_dp ) cycle

    kx     = km_coul%kptx(ik)
    ky     = km_coul%kpty(ik)
    kz     = km_coul%kptz(ik)
    kk     = km_coul%kptk(ik)
    Ak     = km_coul%Ak( ik )
    !write(6,'(a,i,4e16.8)') 'Ak',ik,Ak,kx,ky,kz
    
    rhonk_R = 0.0_dp
    rhonk_I = 0.0_dp
    do ia = 1, nions 
      rxi = rx(ia)
      ryi = ry(ia)
      rzi = rz(ia)
      k_dot_r  = ( kx * rxi + ky * ryi + kz * rzi )
      ckr(ia)  = COS(k_dot_r)
      skr(ia)  = SIN(k_dot_r)
      do k=1, 3
        muix = dmu ( ia , 1 , k)
        muiy = dmu ( ia , 2 , k)
        muiz = dmu ( ia , 3 , k)
        k_dot_mu = ( muix * kx + muiy * ky + muiz * kz )
        rhonk_R(k) = rhonk_R(k) - k_dot_mu * skr(ia)
        rhonk_I(k) = rhonk_I(k) + k_dot_mu * ckr(ia) ! rhon_R + i rhon_I
      enddo
    enddo
  !  write(6,'(i,a,3e16.8)') ik,' rhonk_R ',rhonk_R
  !  write(6,'(i,a,3e16.8)') ik,' rhonk_I ',rhonk_I

    do ia = 1 , nions 
      do k=1,3
        recarg  = Ak * ( rhonk_I(k) * ckr(ia) - rhonk_R(k) * skr(ia) )
        exij = kx * recarg
        eyij = ky * recarg
        ezij = kz * recarg
        f_rec ( ia , 1 , k) = f_rec ( ia , 1 ,k) - exij
        f_rec ( ia , 2 , k) = f_rec ( ia , 2 ,k) - eyij
        f_rec ( ia , 3 , k) = f_rec ( ia , 3 ,k) - ezij
   !     write(6,'(i,4e16.8)') ia,recarg,exij,eyij,ezij
   !     write(6,'(i,3e16.8)') ia,f_rec ( ia , : , k)
      enddo
    enddo

  enddo kpoint

  selfa  = alphaES / piroot
  selfa2 = 2.0_dp * selfa * alpha2 / 3.0_dp
  tpi_V  = tpi    / simu_cell%omega  ! 2pi / V
  fpi_V  = tpi_V  * 2.0_dp           ! 4pi / V
  !f_rec  = f_rec  * 2.0_dp
  f_rec  = f_rec  * fpi_V

  ! ====================================================== 
  !              Self contribution 
  ! ====================================================== 
  ! electrostic energy 
  do k=1,3
    do ia = 1 , nions 
      f_self( ia , 1 ,k) = 2.0_dp * selfa2 * dmu ( ia , 1 ,k)
      f_self( ia , 2 ,k) = 2.0_dp * selfa2 * dmu ( ia , 2 ,k)
      f_self( ia , 3 ,k) = 2.0_dp * selfa2 * dmu ( ia , 3 ,k)
    enddo
  enddo

  f = 0.0d0
  !f = f_dir + f_rec + f_self
  f = f_dir + f_self
  !f = f_dir 
  
  do ia=1,nions
    write(6,'(a)') '------------------------------'
    write(6,'(a,i)') 'ion : ', ia
    write(6,'(a)') '------------------------------'
    write(6,'(a)') ''
    write(6,'(a)') 'f =  '
    write(6,'(a,3e16.8)') 'x',( f(ia,1,j) , j=1,3)
    write(6,'(a,3e16.8)') 'y',( f(ia,2,j) , j=1,3)
    write(6,'(a,3e16.8)') 'z',( f(ia,3,j) , j=1,3)
    write(6,'(a)') ''
    write(6,'(a)') 'f (rec)=  '
    write(6,'(a,3e16.8)') 'x',( f_rec(ia,1,j) , j=1,3)
    write(6,'(a,3e16.8)') 'y',( f_rec(ia,2,j) , j=1,3)
    write(6,'(a,3e16.8)') 'z',( f_rec(ia,3,j) , j=1,3)
    write(6,'(a)') ''
    write(6,'(a)') 'f (dir)=  '
    write(6,'(a,3e16.8)') 'x',( f_dir(ia,1,j) , j=1,3)
    write(6,'(a,3e16.8)') 'y',( f_dir(ia,2,j) , j=1,3)
    write(6,'(a,3e16.8)') 'z',( f_dir(ia,3,j) , j=1,3)
    write(6,'(a)') ''
    write(6,'(a)') 'f (self)=  '
    write(6,'(a,3e16.8)') 'x',( f_self(ia,1,j) , j=1,3)
    write(6,'(a,3e16.8)') 'y',( f_self(ia,2,j) , j=1,3)
    write(6,'(a,3e16.8)') 'z',( f_self(ia,3,j) , j=1,3)


    invf ( ia , : , : )  = f ( ia , : , : )
    CALL DGETRF( 3, 3, invf(ia,:,:), 3, ipiv, ierr )
    if ( ierr.lt.0 ) then
      WRITE( 6 , '(a,i6)' ) 'ERROR call to DGETRF failed in induced_moment',ierr
      STOP
    endif
    CALL DGETRI( 3 , invf(ia,:,:) , 3 ,  ipiv , WORK, LWORK, ierr )
    if ( ierr.lt.0 ) then
      WRITE( 6, '(a,i6)' ) 'ERROR call to DGETRI failed in induced_moment',ierr
      STOP
    endif
    write(6,'(a)') 'f^-1 =  '
    write(6,'(a,3e16.8)') 'x',( invf(ia,1,j) , j=1,3)
    write(6,'(a,3e16.8)') 'y',( invf(ia,2,j) , j=1,3)
    write(6,'(a,3e16.8)') 'z',( invf(ia,3,j) , j=1,3)
    write(6,'(a)') ''
    write(6,'(a)') 'dmu =  '
    write(6,'(a,3e16.8)') 'x',( dmu(ia,1,j) , j=1,3)
    write(6,'(a,3e16.8)') 'y',( dmu(ia,2,j) , j=1,3)
    write(6,'(a,3e16.8)') 'z',( dmu(ia,3,j) , j=1,3)
    write(6,'(a)') ''

  pola(ia,:,:)=0.0_dp
  call dgemm("N","N",3,3,3,1.0_dp,invf(ia,:,:),3,dmu(ia,:,:),3,0.0_dp,pola(ia,:,:),3)

  write(6,'(a)') 'pola =  '
  write(6,'(a,3e16.8)') 'x',( pola(ia,1,j) , j=1,3) 
  write(6,'(a,3e16.8)') 'y',( pola(ia,2,j) , j=1,3) 
  write(6,'(a,3e16.8)') 'z',( pola(ia,3,j) , j=1,3) 
  write(6,'(a)') ''
  write(6,'(a,e16.8)') 'alpha = ',pola(ia,1,1)+pola(ia,2,2)+pola(ia,3,3) / 3.0_dp

  enddo


END MODULE polarizability
