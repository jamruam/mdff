!e ===== fmV =====
MODULE vib

  implicit none

  integer :: ncvib        ! number of configurations in ISCFF to be analysed
  integer :: ngconf       ! number of configurations that would be generated 
  integer :: imod           
  integer :: nkphon       ! number of kpoint between ks and kf used when calc = 'vib+band'

  double precision :: ksx, ksy, ksz ! start kpoint
  double precision :: kfx, kfy, kfz ! final kpoint    
  double precision :: resdos        ! resolution in density of states     
  double precision :: omegamax      ! maximum value in dos

  logical :: lwrite_vectff ! write the vector field

CONTAINS

!*********************** SUBROUTINE vib_init **********************************
!
! initialize vib calculation
!
!******************************************************************************

SUBROUTINE vib_init

  USE io_file,  ONLY :  stdin, stdout ,kunit_OUTFF

  implicit none

  character*132 :: filename

  namelist /vibtag/ ncvib         , & 
                    ngconf        , & 
                    lwrite_vectff , &
                    imod          , & 
                    nkphon        , &
                    resdos        , &
                    omegamax      , &
                    ksx , ksy , ksz   , &
                    kfx , kfy , kfz   
            
  ! ======================
  !  vibtag default values
  ! ======================
  CALL vib_default_tag
  ! ======================
  !  read vibtag namelist
  ! ======================
  CALL getarg (1,filename)
  OPEN ( stdin , file = filename)
  READ ( stdin , vibtag)
  CLOSE  ( stdin )
  ! CALL vib_check_tag
  ! ======================
  !  vibtag print info
  ! ======================
  CALL vib_print_info(stdout)
  CALL vib_print_info(kunit_OUTFF)


  return
  
END SUBROUTINE vib_init

!*********************** SUBROUTINE vib_default_tag ***************************
!
! set default values to vib tag
!
!******************************************************************************

SUBROUTINE vib_default_tag

  implicit none

  ! =================
  !  default values
  ! =================
  ncvib        = 0
  ngconf       = 0
  imod         = 4
  resdos       = 1.0d0             
  lwrite_vectff = .false.
  omegamax = 100.0d0

  return

END SUBROUTINE vib_default_tag


!*********************** SUBROUTINE vib_check_tag ***************************
!
! check vib tag values
!
!******************************************************************************

!SUBROUTINE vib_check_tag
!
!  USE io_file,  ONLY :  stdout
!
!  implicit none
!
!  return
!
!END SUBROUTINE vib_check_tag

!*********************** SUBROUTINE vib_print_info ***************************
!
! print information about vib calculation
!
!******************************************************************************

SUBROUTINE vib_print_info(kunit)

  USE control,  ONLY :  lpbc , calc 
  USE config,   ONLY :  natm , box, ntype, rho
  USE io_file,  ONLY :  ionode

  implicit none

  ! local
  integer :: kunit

  if( ionode ) then
                               WRITE ( kunit ,'(a)')       ''
                               WRITE ( kunit ,'(a)')       'NORMAL MODE ANALYSIS             '
                               WRITE ( kunit ,'(a)')       'Hessian calc. and generation of config.' 
                               WRITE ( kunit ,'(a)')       'for a given mode'
                               WRITE ( kunit ,'(a)')       ''
     if(calc.eq.'vib + fvib')  then  
                               WRITE ( kunit ,'(a)')       'Fvibcalc:'
                               WRITE ( kunit ,'(a)')       'This program reads in the vibrational frequencies '
                               WRITE ( kunit ,'(a)')       'and calculates the frequency dependent part of    ' 
                               WRITE ( kunit ,'(a)')       'the vibrational free energy fvib.                 ' 
                               WRITE ( kunit ,'(a)')       'It also histograms the fvib according to the energy'
                               WRITE ( kunit ,'(a)')       'of the minima. ener are energies that are READ  in,'
                               WRITE ( kunit ,'(a)')       'eigval are eigen values, fvib is calculated in the' 
                               WRITE ( kunit ,'(a)')       'program, and so is fvibhist.'
                               WRITE ( kunit ,'(a)')       'Original author S. Sastry JNCASR and probably F. Affouard'
     endif
     if(.not.lpbc)             WRITE ( kunit ,'(a)')       'NO PERIODIC BOUNDARY CONDITIONS (cubic box)'
     if(lpbc)                  WRITE ( kunit ,'(a)')       'periodic boundary conditions in cubic cell'
                               WRITE ( kunit ,'(a)')       ''
                               WRITE ( kunit ,'(a)')       'configuration (at equilibrium) file : ISCFF'
                               WRITE ( kunit ,'(a,i5)')    'number of configurations in file    = ',ncvib
                               WRITE ( kunit ,'(a)')       ''
                               WRITE ( kunit ,'(a)')       'START HESSIAN CALCULATION            '
                               WRITE ( kunit ,'(a)')       'save eingenvalues                     : EIGFF'
    if ( lwrite_vectff )       WRITE ( kunit ,'(a)')       'save eingenvectors                    : VECTFF'
                               WRITE ( kunit ,'(a)')       ''   
  endif



  return

END SUBROUTINE vib_print_info

!*********************** SUBROUTINE VIB_MAIN  *********************************
!
! main program of the vib calculation:
!
! this subroutine reads configuration from ISCFF (usually optimazed structures)
!
!******************************************************************************

SUBROUTINE vib_main 

  USE config,           ONLY :  system , natm , rx , ry , rz , fx , fy , fz , atype , box , omega , rho , ntype, config_alloc , list , point
  USE md,               ONLY :  temp
  USE control,          ONLY :  calc , myrank , numprocs
  USE io_file,          ONLY :  ionode , stdout , kunit_ISCFF , kunit_EIGFF , kunit_VECTFF , kunit_DOSFF , kunit_MODFF, kunit_OUTFF
  USE thermodynamic,    ONLY :  u_tot , pressure_tot , calc_thermo
  USE field

  implicit none
  INCLUDE "mpif.h"

  ! local
  integer :: i , ia , ja , imod , ic , PANdos , ka 
  double precision :: pressure0, pot0, ak
  double precision, dimension(:), allocatable :: fx_sum, fy_sum, fz_sum
  integer, dimension(:), allocatable :: dostab
  integer, dimension(:), allocatable :: dostabtot
  double precision, dimension(:,:), allocatable :: hess
  double precision, dimension(:), allocatable :: work, deig, ipiv !,deig_reorder(3 * n)
  CHARACTER * 1 jobz, uplo
  INTEGER * 4   info, lwork
  integer :: iastart , iaend

  ! trash 
  double precision :: aaaa
  integer :: iiii
  character * 60 :: cccc

  PANdos = dble(omegamax/resdos) + 1

  allocate(dostabtot(0:PANdos + 1))
  allocate(dostab(0:PANdos + 1))
  dostabtot = 0
 
  OPEN ( UNIT = kunit_ISCFF , FILE = 'ISCFF'  )
  OPEN ( UNIT = kunit_EIGFF , FILE = 'EIGFF'  )
  OPEN ( UNIT = kunit_VECTFF, FILE = 'VECTFF' )
  OPEN ( UNIT = kunit_DOSFF , FILE = 'DOSFF'  )

  if(calc.ne.'vib+band'.and.calc.ne.'vib+dos') then 
  if ( ionode ) WRITE ( kunit_VECTFF ,'(a)') '       #rx               ry             dx              dy              rx              rz              dx             dz              ry              rz              dy             dz           eigenvalue'   
  endif

  ! ==============================
  !  read main config parameters
  ! ==============================
  READ ( kunit_ISCFF , * ) natm  
  READ ( kunit_ISCFF , * ) system 
  READ ( kunit_ISCFF , * ) box,ntype
  omega = box * box * box      
  rho = dble ( natm ) / omega 
  ! ===================================
  !  here we know natm, then alloc 
  !  and decomposition can be applied 
  ! ================================== 
  CALL config_alloc 
  CALL do_split ( natm , myrank , numprocs , iastart , iaend )
  allocate(fx_sum(natm),fy_sum(natm),fz_sum(natm))
  allocate(hess(3 * natm,3 * natm))
  allocate(work(9 * natm))
  allocate(deig(3 * natm),ipiv(3 * natm))
  WRITE ( stdout ,'(a)')         'Remind some parameters of the system:'
  WRITE ( stdout ,'(a,i5)')      'natm  = ',natm
  WRITE ( stdout ,'(a,i5)')      'ntype = ',ntype
  WRITE ( stdout ,'(a,f6.3)')    'rho   = ',rho
  WRITE ( stdout ,'(a,f6.3)')    'box   = ',box
  WRITE ( stdout ,'(a,f6.3)')    'vol   = ',omega
  WRITE ( stdout ,'(a)')         ''
  WRITE ( stdout ,'(a,i5,a,i5)') 'hessian dimension',3 * natm,'x',3 * natm

  ! ===========================================
  !  LOOP OVER CONFIGURATIONS (AT EQUILIBRIUM)
  ! ===========================================
  do ic = 1, ncvib
    write( stdout ,'(a,i5)') 'read config',ic
    write( kunit_OUTFF ,'(a,i5)') 'read config',ic
    if ( ic .ne. 1 ) READ ( kunit_ISCFF , * ) iiii
    if ( ic .ne. 1 ) READ ( kunit_ISCFF , * ) cccc
    if ( ic .ne. 1 ) READ ( kunit_ISCFF , * ) aaaa 
    do ia = 1,natm
      READ ( kunit_ISCFF , * ) atype(ia),rx(ia),ry(ia),rz(ia)
    enddo

    ! ==============================
    ! force ???? should not be here 
    ! there is no need ... to check!
    ! ==============================
    CALL engforce ( iastart , iaend )!, list , point )

    CALL calc_thermo
    pot0      = u_tot
    pressure0 = pressure_tot 

    ! ========================
    !  get the hessian matrix
    ! ========================
    CALL hessian(hess)

    ! ===========================
    !  real space dos: 3N modes
    ! ===========================
    if(calc.ne.'vib+band'.and.calc.ne.'vib+dos') then

      ! ===========================
      !  diagonalisation of hess
      ! ===========================
      write( stdout ,'(a)') 'start diagonalization'
      write( kunit_OUTFF ,'(a)') 'start diagonalization'
      lwork = 9 * natm
      jobz = 'V'
      uplo = 'U'
      CALL DSYEV( jobz , uplo , 3 * natm , hess ,3 * natm , deig , work , lwork , info )
      if(info.ne.0)then
        if ( ionode ) WRITE ( stdout ,'(a,i5)') 'ERROR: Improper termination. of DSYEV info = ',info
        STOP 
      else
        if ( ionode ) WRITE ( stdout ,'(a)') 'hessian diagonalisation ok'
        if ( ionode ) WRITE ( kunit_OUTFF ,'(a)') 'hessian diagonalisation ok'
      endif

      ! ===========================================
      !  write frequencies^2 (hessian eigenvalues) 
      ! ===========================================
      do imod = 1,3 * natm
        if ( ionode ) WRITE ( kunit_EIGFF ,'(f24.10)') deig(imod)
      enddo
      ! =========================================
      !  DOS: density of states  of the 3N modes
      ! =========================================
      dostab = 0
      do i = 4,3 * natm
        ak = (dsqrt(deig(i)))/resdos
        ka = int(ak) + 1
        if(ka.gt.PANdos + 1.or.ka.lt.0) then
           WRITE ( stdout , * ) 'ERROR out of bound in dostab'
           WRITE ( stdout , * ) i,ka,PANdos + 1,ak,deig(i)
        endif
        dostab(ka) = dostab(ka) + 1
        dostabtot(ka) = dostabtot(ka) + 1
      enddo
      dostab(0) = 3
      dostabtot(0) = dostabtot(0) + 3

      do i = 0,PANdos + 1
        if ( ionode ) WRITE ( kunit_DOSFF ,'(3f16.8)') dble(i) * resdos,dostab(i) / ( 3.0 * natm * resdos ) , dostabtot(i) / ( 3.0 * natm * resdos * ic)
      enddo
      if ( ionode ) WRITE ( kunit_DOSFF ,'(a)') ''  
      if ( ionode ) WRITE ( kunit_DOSFF ,'(a)') ''  
      if ( ionode ) WRITE ( stdout ,'(a)') 'dos ok' 
      if ( ionode ) WRITE ( kunit_OUTFF ,'(a)') 'dos ok' 


      ! ====================
      !  write vector field
      ! ====================
      if ( lwrite_vectff ) then
        do imod = 1,3 * natm
          do ja = 1, natm
            if ( ionode ) WRITE  ( kunit_VECTFF ,'(13f16.8)') rx(ja),ry(ja),hess(ja,imod),hess(ja + natm,imod), &
                                                              ry(ja),rz(ja),hess(ja + natm,imod),hess(ja + 2 * natm,imod), &
                                                              rx(ja),rz(ja),hess(ja,imod),hess(ja + 2 * natm,imod), &
                                                              deig(imod)
          enddo
          if ( ionode ) WRITE  ( kunit_VECTFF , * ) ' '
          if ( ionode ) WRITE  ( kunit_VECTFF , * ) ' '
        enddo
      endif

    endif 
    ! ==============================================
    !  beware hess is changing when diagonalized!!! 
    ! ==============================================
    !  using reciprocal space we get the complete : 
    ! ==============================================
    !  band structure    
    ! ==============================================
    if( calc .eq.'vib+band' ) CALL band(hess)
    ! ==============================================
    ! total dos for more than 3N k - vector
    ! ==============================================
    if( calc .eq.'vib+dos'  ) CALL doskpt(hess)


  enddo ! loop over configurations

  CLOSE ( kunit_DOSFF ) 
  CLOSE ( kunit_VECTFF )
  CLOSE ( kunit_EIGFF )
  CLOSE ( kunit_ISCFF )

  ! ===================================
  !  generate config from a given mode
  ! ===================================
  if(calc.eq.'vib+gmod') then
    if( ionode ) then
      WRITE ( kunit_OUTFF ,'(a)')            ''
      WRITE ( kunit_OUTFF ,'(a,i4,a,f8.3)')  'generate ', ngconf ,' configurations at temp  = ', temp
      WRITE ( kunit_OUTFF ,'(a)')            'save configurations in file          : MODFF'
      WRITE ( kunit_OUTFF ,'(a)')            ''
      WRITE ( kunit_OUTFF ,'(a)')            'Method:'
      WRITE ( kunit_OUTFF ,'(a)')            'gaussian distributions of normal coordinates'
      WRITE ( kunit_OUTFF ,'(a)')            'equation 11 of J.Phys.Chem.B 2005, 109, 7245'
      WRITE ( kunit_OUTFF ,'(a)')            ''
    endif

    !generation of modes imod
    OPEN (unit = kunit_MODFF ,FILE = 'MODFF')
    CALL generate_modes(deig,hess,2040)
    CLOSE ( kunit_MODFF )
  endif 

  ! ====================================
  !  start vibrational free energy calc
  ! ====================================
  if(calc.eq.'vib+fvib') then
    if( ionode ) then
      WRITE ( kunit_OUTFF ,'(a)')       ''    
      WRITE ( kunit_OUTFF ,'(a)')       'read thermo (Fvibcalc)              : ISTHFF'
      WRITE ( kunit_OUTFF ,'(a)')       'WARNING !!!! DEV. STATE !!!! WARNING !!! NEED TO BE TESTED!!!'
    endif
    CALL fvibcalc
    STOP
  endif

  deallocate(dostab)
  deallocate(dostabtot)
  deallocate(deig)
  deallocate(work)
  deallocate(hess) 
  deallocate(fx_sum,fy_sum,fz_sum)

  return

END SUBROUTINE vib_main



!*********************** SUBROUTINE hessian ***********************************
!
!
!******************************************************************************

SUBROUTINE hessian ( hess )

  USE config,   ONLY :  natm, rx, ry, rz, box, itype 
  USE io_file,  ONLY :  ionode , stdout
  USE field,    ONLY :  rcutsq , sigsq , epsp , fc , uc , pp , qq

  implicit none
  ! global
  double precision :: hess(3 * natm,3 * natm)

  ! local
  integer i, j, p1, p2, ia, ja
  integer :: nxij,nyij,nzij
  double precision :: rxi, ryi, rzi, rxij, ryij, rzij, rijsq
  double precision :: sr2,sr,srp2,srp4,srq2,srq4,c1,c2        
  double precision :: txx,txy,tyy,tyz,tzz,txz
  integer np2(2,2),nq2(2,2),np4(2,2),nq4(2,2)


  do i = 1,2
    do j = 1,2
      np2(i,j) = pp(i,j) 
      nq2(i,j) = qq(i,j) 
      np4(i,j) = pp(i,j) 
      nq4(i,j) = qq(i,j) 
      np2(i,j) = np2(i,j) + 2
      nq2(i,j) = nq2(i,j) + 2
      np4(i,j) = np4(i,j) + 4
      nq4(i,j) = nq4(i,j) + 4
    enddo
  enddo

  hess = 0.0D0


  do ia = 1,natm - 1
    rxi = rx(ia)
    ryi = ry(ia)
    rzi = rz(ia)
    do ja = ia + 1,natm
      rxij = rxi - rx(ja)
      ryij = ryi - ry(ja)
      rzij = rzi - rz(ja)
      nxij = nint( rxij / box )
      nyij = nint( ryij / box )
      nzij = nint( rzij / box )
      rxij = rxij - box * nxij
      ryij = ryij - box * nyij
      rzij = rzij - box * nzij
      rijsq = rxij * rxij + ryij * ryij + rzij * rzij

      p1 = itype(ia)
      p2 = itype(ja)

      if (rijsq .lt. rcutsq(p1,p2)) then
        sr2 = sigsq(p1,p2)/rijsq
        sr = dsqrt(sr2)
        srp2 = sr ** np2(p1,p2)
        srp4 = sr ** np4(p1,p2)
        srq2 = sr ** nq2(p1,p2)
        srq4 = sr ** nq4(p1,p2)

        ! ===========
        !  Hessian 
        ! ===========
        c1 = (fc(p1,p2)/sigsq(p1,p2)) * ((pp(p1,p2) + 2.0D0) * srp4 - (qq(p1,p2) + 2.0D0) * srq4)
        c2 = fc(p1,p2) * (srq2 - srp2)


        hess( ia , ja ) = rxij * rxij * c1 + c2                        ! x,x
        hess( ia , natm + ja ) = rxij * ryij * c1                      ! x,y
        hess( ia , 2 * natm + ja ) = rxij * rzij * c1                  ! x,z  

        hess( natm + ia,ja) = ryij * rxij * c1                         ! x,y
        hess( natm + ia,natm + ja) = ryij * ryij * c1 + c2             ! y,y  
        hess( natm + ia,2 * natm + ja) = ryij * rzij * c1              ! y,z

        hess( 2 * natm + ia , ja ) = rzij * rxij * c1                  ! x,z
        hess( 2 * natm + ia , natm + ja ) = rzij * ryij * c1           ! y,z
        hess( 2 * natm + ia , 2 * natm + ja ) = rzij * rzij * c1 + c2  ! z,z 

        ! ====================== 
        !  hessian is symmetric 
        ! ====================== 

        hess(ja,ia) = hess(ia,ja) 
        hess(natm + ja,ia) = hess(ia,natm + ja) 
        hess(2 * natm + ja,ia) = hess(ia,2 * natm + ja) 

        hess(ja,natm + ia) = hess(natm + ia,ja) 
        hess(natm + ja,natm + ia) = hess(natm + ia,natm + ja) 
        hess(2 * natm + ja,natm + ia) = hess(natm + ia,2 * natm + ja) 

        hess(ja,2 * natm + ia) = hess(2 * natm + ia,ja) 
        hess(natm + ja,2 * natm + ia) = hess(2 * natm + ia,natm + ja) 
        hess(2 * natm + ja,2 * natm + ia) = hess(2 * natm + ia,2 * natm + ja) 

      endif
    enddo
  enddo

!=================================================================================================================
! elements diagonaux:
! "appelle abusivement self - force. il traduit simplement le fait que
! lorsqu'un atome est deplace, il ressent une force exercee par l'ensemble des autres atom du cristal.   
! La forme de ce terme se clarifie quand on s'aperercoit que, deplacer un atome dans une direction
! donnee, est equivqlent a deplacer l'ensemble du cristal a l'exception de cet atome dans la dition opposee
!
! Notes de Cours de Physique des materiaux
! P. GHOSEZ, J - Y. RATY
! Universite de Liege
! http://www.phythema.ulg.ac.be/Teaching/Cours/Cours - Phys_Mat/Notes - Phys_Mat.pdf
!=================================================================================================================

  do ia = 1, natm
    txx = 0.0D0
    txy = 0.0D0
    txz = 0.0D0
    tyy = 0.0D0
    tyz = 0.0D0
    tzz = 0.0D0
    do ja = 1, natm
      if((ia.eq.ja).and.(hess(ia,ja).ne.0.0D0))then
        if ( ionode ) WRITE ( stdout , * )'Error! ',ia,ja,hess(ia,ja)
      endif
      txx = txx + hess(ia,ja)
      txy = txy + hess(ia,natm + ja)
      txz = txz + hess(ia,2 * natm + ja)
      tyy = tyy + hess(natm + ia,natm + ja)
      tyz = tyz + hess(natm + ia,2 * natm + ja)
      tzz = tzz + hess(2 * natm + ia,2 * natm + ja)
       
    enddo
    ! diagonal = - sum i,j
    hess(ia,ia) = - txx
    hess(ia,natm + ia) = - txy
    hess(ia,2 * natm + ia) = - txz
    hess(natm + ia,natm + ia) = - tyy
    hess(natm + ia,2 * natm + ia) = - tyz
    hess(2 * natm + ia,2 * natm + ia) = - tzz
   
    !symetry
    hess(natm + ia,ia) = hess(ia,natm + ia)
    hess(2 * natm + ia,ia) = hess(i,2 * natm + ia)
    hess(2 * natm + ia,natm + ia) = hess(natm + ia,2 * natm + ia)
        

    if(hess(ia,ia).eq.0.0D0.and. ionode )                       WRITE ( stdout , * )'in forcehes zero at ',ia
    if(hess(ia + natm,ia + natm).eq.0.0D0.and. ionode )         WRITE ( stdout , * )'in forcehes zero at ',natm + ia
    if(hess(ia + 2 * natm,ia + 2 * natm).eq.0.0D0.and. ionode ) WRITE ( stdout , * )'in forcehes zero at ',2 * natm + ia
  enddo


  return

END SUBROUTINE hessian

!*********************** SUBROUTINE fvibcalc **********************************
!
! This program reads in the vibrational frequencies and calculates 
! the frequency dependent part of the vibrational free energy fvib. 
! It also histograms the fvib according to the energy of the minima.
! ener are energies that are READ  in, eigval are eigen values, fvib 
! is calculated in the program, and so is fvibhist.
!
!******************************************************************************

SUBROUTINE fvibcalc

  USE config, 	ONLY : 	natm
  USE io_file,	ONLY :	ionode , stdout , kunit_ISTHFF , kunit_EIGFF , kunit_VIBFF

  implicit none

!lGocal
  integer :: i, j
  integer,parameter :: nbin = 100, ncmax = 100000
  double precision :: ener(ncmax), eigval(3 * natm), fvib(ncmax), fvibhist(nbin,3)
  character * 100  header
  integer :: ncs , nskip ,ifv, ieg, ind
  double precision :: emin, emax
  double precision :: enerind, pres, de, dum, fvibav, pofet,atmp 

  ! trash  
  integer iiii
  double precision :: xxxx 

!# config                 eIS                grad                Pres  Iter  Neng
!u_initial     Press_initial
!       1     - 3.947162344773      0.000000000037      6.708240000000    20    70 - 2.603606441227   13.556241500514
!       2     - 3.947162344773      0.000000000049      6.708240000000    21    74 - 2.328527224939   15.013373656527

  ncs = 1
  emin = 999999.0d0
  emax = - 999999.0d0
  OPEN (unit = kunit_ISTHFF ,file = 'ISTHFF',status = 'old')
  READ ( kunit_ISTHFF , * ) header
  READ ( kunit_ISTHFF , * ) header
  DO i = 1,ncvib
    READ ( kunit_ISTHFF , * ) iiii, enerind, xxxx, pres, iiii, iiii, xxxx
    if(enerind .gt. emax) emax = enerind
    if(enerind .lt. emin) emin = enerind
  ENDDO
  CLOSE ( kunit_ISTHFF )

  de = abs(emax - emin)/real(nbin)

!  PRINT  * ,' *  - ISthermo file = '
!  PRINT  * ,' *  - vib      file = ',F_IN2
!  PRINT  * ,' - fvibn     file = ',F_IN3
  print  * ,'Nconf = ',ncvib
  print  * ,'emin  = ',emin
  print  * ,'emax  = ',emax
  print  * ,'de    = ',de

  OPEN (unit = kunit_ISTHFF ,file = 'ISTHFF',status = 'old')
  READ ( kunit_ISTHFF , * ) header
  READ ( kunit_ISTHFF , * ) header
  OPEN (unit = kunit_EIGFF ,file = 'EIGFF',status = 'old')
  OPEN (unit = kunit_VIBFF ,file = 'VIBFF',status = 'unknown')

! =========================================
!  specify minimum energy and the bin size 
! =========================================

  do i = 1, ncvib
    fvib(i) = 0.0
  enddo

  do i = 1, nbin
    fvibhist(i,1) = emin + (i - .5) * de  
    fvibhist(i,2) = 0.0
    fvibhist(i,3) = 0.0
  enddo

  dum = 10.0d0
  ifv = 0

  do i = 1, ncvib
    READ ( kunit_ISTHFF , * ) iiii, enerind, xxxx, pres, iiii, iiii, xxxx
    if(mod(i,ncs).eq.0)then
      ifv = ifv + 1
      ener(ifv) = enerind 

      do ieg = 1,3 * natm 
      READ ( kunit_EIGFF , * ) eigval(ieg)
      enddo

      nskip = 3
      do j = 4, 3 * natm 
        if(eigval(j).lt.0.0D0)then
          nskip = nskip + 1 
        else 
          fvib(ifv) = fvib(ifv) + log(eigval(j))
        endif
      enddo 

      if(nskip.gt.3.and. ionode ) WRITE ( stdout , * ) ifv, nskip 

      fvib(ifv) = fvib(ifv) * (3.0d0 * natm)/((2.0d0 * natm) * (3.0d0 * natm - nskip * 1.0d0))
    endif 
  enddo

  CLOSE ( kunit_EIGFF )
  CLOSE ( kunit_ISTHFF )
  

  do i = 1, ncvib
    atmp = (ener(i) - emin)/de
    ind = int(atmp)  
    ind = ind + 1
    fvibhist(ind,2) = fvibhist(ind,2) + 1.0d0
    fvibhist(ind,3) = fvibhist(ind,3) + fvib(i)
  enddo 

  do i = 1, nbin 
    if(fvibhist(i,2).ne.0.0) then
      if( ionode ) WRITE ( kunit_VIBFF ,'(3(2x,e14.6))') fvibhist(i,1),fvibhist(i,3)/fvibhist(i,2),fvibhist(i,2)/(de * ncvib)
    endif
  enddo


  fvibav = 0.0d0
  pofet = 0.0d0

  do i = 1, nbin 

    if(fvibhist(i,2).ne.0.0) then
      fvibav = fvibav + fvibhist(i,3)
      pofet = pofet + fvibhist(i,2)
    endif 
  enddo

  if( ionode ) WRITE ( stdout , * ) fvibav/pofet, pofet

  return

END SUBROUTINE fvibcalc

!*********************** SUBROUTINE generate_modes ****************************
!
! generate configuration from a given mode
!
!******************************************************************************

SUBROUTINE generate_modes ( deig , hess , kunit )

  USE config,   ONLY :  system , natm , rx , ry , rz , atype
  USE md,       ONLY :  temp 
  USE io_file,  ONLY :  ionode , stdout
  USE field 

  implicit none

  ! global
  double precision :: hess(3 * natm,3 * natm)
  double precision :: deig(3 * natm)
  integer :: kunit
  ! local
  integer :: ia, ja, igconf
  double precision, dimension(:), allocatable :: rrx , rry , rrz
  double precision, dimension(:), allocatable :: qx , qxd , xsq
  double precision :: dexpo, omeg, G1

  allocate( rrx ( natm ) , rry ( natm ) , rrz ( natm ) )
  allocate( qx ( 3 * natm ) , qxd ( 3 * natm ) , xsq ( 3 * natm ) )  ! normal coordinates

  ! ===============================
  !  inverse of eigenvtors matrix
  ! ===============================
  !CALL DGETRF( 3 * natm, 3 * natm, hess, 3 * natm, ipiv, info )
  !CALL DGETRI( 3 * natm, hess, 3 * natm, ipiv, work, lwork, info )

  ! ===========================================
  !  transformations to get normal coordinates
  ! ===========================================
!  qx = 0.0D0
!  do i = 1,3 * natm
    do ja = 1,natm
      qx(imod) = qx(imod) + hess(ja,imod) * rx(ja) + hess(ja + natm,imod) * ry(ja) + hess(ja + 2 * natm,imod) * rz(ja)
    enddo
!  enddo

  ! ===========================================
  !  following MgO mauri paper 
  !  xsq(i): width of gaussian distribution 
  !  on q(i) coordinates (normal coordinates) 
  !  of mode i
  ! ===========================================

!  do i = 1,3 * natm ! loop over modes
    omeg = dsqrt(deig(imod))
    dexpo = dexp( - omeg/temp)
    xsq(imod) = 1.0d0 + dexpo
    xsq(imod) = xsq(imod) / (2.0d0 * omeg * (1.0d0 - dexpo) )
!  enddo

  ! =============================================
  !  generate ngconf configurations of mode imod
  ! =============================================
  do igconf = 1,ngconf

    if( ionode ) WRITE ( stdout ,'(a,i4,a,i4)') 'conf = ',igconf,' of mode',imod

        ! ========================================================================
        ! gaussian distribution: 
        !   mean value: equilibrium position rx,ry,rz - > qx, qy, qz 
        !   width x^2: equation (11) of   J. Phys. Chem. B 2005, 109, 7245 - 7250
        ! ========================================================================
        CALL boxmuller_basic(G1, 0.0d0, xsq(imod)) 
        qxd(imod) = G1

      ! ==========================
      ! inverse transformations 
      ! to get usual coordinates
      ! ==========================
      rrx = 0.0d0
      rry = 0.0d0
      rrz = 0.0d0
      do ia = 1,natm
          rrx(ia) = rrx(ia) + hess(ia,imod) * qxd(imod)
          rry(ia) = rry(ia) + hess(ia + natm,imod) * qxd(imod)
          rrz(ia) = rrz(ia) + hess(ia + 2 * natm,imod) * qxd(imod)
      enddo

      if( ionode ) then
        WRITE (kunit, * ) natm
        WRITE (kunit, * ) system  
        do ia  = 1,natm
          WRITE (kunit,'(a,3f16.10)') atype(ia),rx(ia) + rrx(ia),ry(ia) + rry(ia),rz(ia) + rrz(ia)
        enddo
      endif

  enddo 

  deallocate(rrx,rry,rrz)
  deallocate(qx,qxd,xsq)

  return 

END SUBROUTINE generate_modes



!*********************** SUBROUTINE band **************************************
!
! calculates dispersion curve (band) in a given direction
!
!******************************************************************************

SUBROUTINE band ( hess )

  USE config,       ONLY :  natm , rx , ry , rz , box , ncell , itype
  USE control,      ONLY :  calc
  USE constants,    ONLY :  tpi , pi
  USE io_file,      ONLY :  stdout , kunit_DOSKFF
  USE field,        ONLY :  rcutsq , sigsq , epsp , fc , uc , pp , qq

  implicit none

  ! global  
  double precision :: hess(3 * natm,3 * natm)
  ! local
  integer :: nxij,nyij,nzij
  double precision :: rxi, ryi, rzi, rxij, ryij, rzij, rijsq
  double precision :: kxt, kyt, kzt
  integer :: ik, ia, ja, p1, p2
  double precision :: ak, akn
  double precision :: hessij(3,3), sphase
  double precision :: kxi, kyi, kzi, krx, kry, krz
  CHARACTER * 1 jobz, uplo
  INTEGER * 4   info, lwork
  double precision :: tmpxx,tmpxy,tmpyy,tmpyz,tmpzz,tmpxz
  double precision :: ww(3), work(9)


  OPEN (UNIT = kunit_DOSKFF ,FILE = 'DOSKFF')

  lwork = 9
  jobz = 'V'
  uplo = 'U'

  ak = (tpi * ncell)/(box)
  akn = (tpi * ncell)/(box * nkphon)

  kxt = kfx - ksx
  kyt = kfy - ksy
  kzt = kfz - ksz

  WRITE ( stdout ,'(a)') ''
  WRITE ( stdout ,'(a)') 'dispersion curves:'
  WRITE ( stdout ,'(a)') '' 
  WRITE ( stdout ,'(a,f5.2,a,f5.2,a,f5.2,a)') 'from point : [ ',ksx,' , ',ksy,' , ',ksz,' ]'
  WRITE ( stdout ,'(a,f5.2,a,f5.2,a,f5.2,a)') 'to         : [ ',kfx,' , ',kfy,' , ',kfz,' ]' 
  WRITE ( stdout ,'(a,f5.2,a,f5.2,a,f5.2,a)') 'direction  : [ ',kxt,' , ',kyt,' , ',kzt,' ]'
  WRITE ( stdout ,'(a)') ''      

  do ik = 0,nkphon - 1
     kxi = (kxt * akn * ik)  + ksx * ak
     kyi = (kyt * akn * ik)  + ksy * ak 
     kzi = (kzt * akn * ik)  + ksz * ak

      tmpxx = 0.0d0
      tmpyy = 0.0d0
      tmpzz = 0.0d0
      tmpxy = 0.0d0
      tmpxz = 0.0d0
      tmpyz = 0.0d0
      hessij = 0.D0
      work = 0.0d0
      ww = 0.0d0

    do ia = 1,natm
      rxi = rx(ia)
      ryi = ry(ia)
      rzi = rz(ia)
      do ja = 1,natm
        rxij = rxi - rx(ja)
        ryij = ryi - ry(ja)
        rzij = rzi - rz(ja)

        nxij = nint( rxij / box )
        nyij = nint( ryij / box )
        nzij = nint( rzij / box )

        rxij = rxij - box * nxij
        ryij = ryij - box * nyij
        rzij = rzij - box * nzij

        rijsq = rxij * rxij + ryij * ryij + rzij * rzij

        p1 = itype(ia)
        p2 = itype(ja)

        if (rijsq .lt. rcutsq(p1,p2)) then
          krx = (kxi * rxij)
          kry = (kyi * ryij)
          krz = (kzi * rzij)
          sphase = dsin( 0.5d0 * (krx + kry + krz) )
          sphase = sphase * sphase

          tmpxx = tmpxx + hess(ja,ia)               * sphase 
          tmpyy = tmpyy + hess(ja + natm,ia + natm)     * sphase
          tmpzz = tmpzz + hess(ja + 2 * natm,ia + 2 * natm) * sphase
          tmpxy = tmpxy + hess(ja,ia + natm)          * sphase
          tmpxz = tmpxz + hess(ja,ia + 2 * natm)        * sphase
          tmpyz = tmpyz + hess(ja + natm,ia + 2 * natm)   * sphase
        endif
      enddo  ! j atom loop
    enddo ! i atom loop
       hessij(1,1) = - tmpxx * 2.0d0 / natm
       hessij(2,2) = - tmpyy * 2.0d0 / natm
       hessij(3,3) = - tmpzz * 2.0d0 / natm

       hessij(1,2) = - tmpxy * 2.0d0 / natm
       hessij(1,3) = - tmpxz * 2.0d0 / natm
       hessij(2,3) = - tmpyz * 2.0d0 / natm
 
       hessij(2,1) = hessij(1,2)
       hessij(3,1) = hessij(1,3)
       hessij(3,2) = hessij(2,3)

       CALL DSYEV(jobz,uplo,3,hessij,3,ww,work,lwork,info)

       if(mod(ik + 1,10).eq.0) WRITE ( * ,'(i7,3f12.4,3f18.6)')  ik, kxi, kyi, kzi, dsqrt(ww(1)),dsqrt(ww(2)), dsqrt(ww(3))
       WRITE ( kunit_DOSKFF ,'(i7,3f12.4,3f18.6)')  ik, kxi, kyi, kzi, dsqrt(ww(1)),dsqrt(ww(2)), dsqrt(ww(3))

  enddo

  CLOSE ( kunit_DOSKFF )

  return  

END SUBROUTINE band 


!*********************** SUBROUTINE doskpt ************************************
!
! calculates the DOS for a given set of k-points
!
!******************************************************************************

SUBROUTINE doskpt ( hess )

  USE config,           ONLY :  natm , rx , ry , rz , box , ncell , itype
  USE control,          ONLY :  calc
  USE constants,        ONLY :  tpi , pi
  USE io_file,          ONLY :  stdout , kunit_IBZKPTFF , kunit_DOSKFF
  USE field,            ONLY :  rcutsq , sigsq , epsp , fc , uc , pp , qq

  implicit none

  ! global  
  double precision :: hess(3 * natm,3 * natm)
  ! local
  integer :: nxij,nyij,nzij
  double precision :: rxi, ryi, rzi, rxij, ryij, rzij, rijsq
  integer :: ck ,kl, wi
  integer :: ik, i, ia, ja, p1, p2
  double precision :: ak, akn
  double precision :: hessij(3,3), sphase
  double precision :: kxi, kyi, kzi, krx, kry, krz
  CHARACTER * 1 jobz, uplo
  INTEGER * 4   info, lwork
  double precision :: tmpxx,tmpxy,tmpyy,tmpyz,tmpzz,tmpxz
  double precision :: ww(3), work(9)
  ! trash
  character * 60 :: cccc

  OPEN (UNIT = kunit_DOSKFF ,FILE = 'DOSKFF')

  lwork = 9
  jobz = 'V'
  uplo = 'U'

  wi = 1
  ak = (tpi * ncell)/(box)
  akn = (tpi * ncell)/(box * nkphon)


  OPEN (UNIT = kunit_IBZKPTFF ,FILE = 'IBZKPTFF')
  READ ( kunit_IBZKPTFF , * ) cccc
  READ ( kunit_IBZKPTFF , * ) nkphon
  READ ( kunit_IBZKPTFF , * ) cccc

  WRITE ( kunit_DOSKFF , * ) nkphon,'20  0.0 30.0'

  wi = 1
  WRITE ( stdout ,'(a)') 'read kpoints in IBZKPTFF'

  kl = 1
  ck = 1
  do ik = 0,nkphon - 1

    READ ( kunit_IBZKPTFF , * ) kxi , kyi , kzi , wi

    kxi = (kxi * ak)
    kyi = (kyi * ak)
    kzi = (kzi * ak)


    tmpxx = 0.0d0
    tmpyy = 0.0d0
    tmpzz = 0.0d0
    tmpxy = 0.0d0
    tmpxz = 0.0d0
    tmpyz = 0.0d0
    hessij = 0.D0
    work = 0.0d0
    ww = 0.0d0

    do ia = 1,natm

      rxi = rx(ia)
      ryi = ry(ia)
      rzi = rz(ia)
      do ja = 1,natm

        rxij = rxi - rx(ja)
        ryij = ryi - ry(ja)
        rzij = rzi - rz(ja)
        nxij = nint( rxij / box )
        nyij = nint( ryij / box )
        nzij = nint( rzij / box )
        rxij = rxij - box * nxij
        ryij = ryij - box * nyij
        rzij = rzij - box * nzij
        rijsq = rxij * rxij + ryij * ryij + rzij * rzij

        p1 = itype(ia)
        p2 = itype(ja)

        if (rijsq .lt. rcutsq(p1,p2)) then
          krx = (kxi * rxij)
          kry = (kyi * ryij)
          krz = (kzi * rzij)
          sphase = dsin( 0.5d0 * (krx + kry + krz) )
          sphase = sphase * sphase

          tmpxx = tmpxx + hess(ja,ia)                       * sphase
          tmpyy = tmpyy + hess(ja + natm,ia + natm)         * sphase
          tmpzz = tmpzz + hess(ja + 2 * natm,ia + 2 * natm) * sphase
          tmpxy = tmpxy + hess(ja,ia + natm)                * sphase
          tmpxz = tmpxz + hess(ja,ia + 2 * natm)            * sphase
          tmpyz = tmpyz + hess(ja + natm,ia + 2 * natm)     * sphase
        endif
      enddo  ! j atom loop
    enddo ! i atom loop

    hessij(1,1) = - tmpxx * 2.0d0 / natm
    hessij(2,2) = - tmpyy * 2.0d0 / natm
    hessij(3,3) = - tmpzz * 2.0d0 / natm

    hessij(1,2) = - tmpxy * 2.0d0 / natm
    hessij(1,3) = - tmpxz * 2.0d0 / natm
    hessij(2,3) = - tmpyz * 2.0d0 / natm

    hessij(2,1) = hessij(1,2)
    hessij(3,1) = hessij(1,3)
    hessij(3,2) = hessij(2,3)

    CALL DSYEV(jobz,uplo,3,hessij,3,ww,work,lwork,info)

    do i = 1,wi
      if(mod(ik + 1,kl).eq.0) then
        WRITE ( stdout ,'(i7,3f12.4,3f18.6)')  ck, kxi, kyi, kzi, dsqrt(ww(1)),dsqrt(ww(2)), dsqrt(ww(3))
        kl = kl * 2
      endif
      WRITE ( kunit_DOSKFF ,'(i7,3f12.4,3f18.6)')  ck, kxi, kyi, kzi, dsqrt(ww(1)),dsqrt(ww(2)), dsqrt(ww(3))
      ck = ck + 1
    enddo

  enddo
  CLOSE ( kunit_IBZKPTFF )
  CLOSE ( kunit_DOSKFF )

  return  

END SUBROUTINE doskpt 


END MODULE vib

! ===== fmV =====
