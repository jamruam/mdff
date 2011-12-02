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
MODULE opt

  implicit none

  integer :: ncopt             ! number of configurations in TRAJFF  
  integer :: nskipopt          ! number of configurations skipped in the beginning 
  integer :: nmaxopt           ! number of configurations optimized  
  integer :: nperiodopt

  character*60 :: optalgo      ! choose optimization algorithm
  character*60 :: optalgo_allowed(3)
  data optalgo_allowed  / 'sastry' , 'lbfgs' , 'm1qn3' /
  

CONTAINS

!*********************** SUBROUTINE opt_init **********************************
!
! optimisation initialisation
!
!******************************************************************************

SUBROUTINE opt_init

  USE io_file,  ONLY :  stdin, stdout, kunit_OUTFF, ionode

  implicit none
 
  integer :: ioerr 
  character*132 :: filename

  namelist /opttag/ ncopt         , & 
                    nskipopt      , & 
                    nmaxopt       , & 
                    optalgo           

  ! ===============================
  !  set default values to opttag
  ! ===============================
  CALL opt_default_tag
  ! =======================
  !  read opttag namelist
  ! =======================
  CALL getarg ( 1 , filename )
  OPEN ( stdin , file = filename )
  READ ( stdin , opttag ,iostat=ioerr)
  if( ioerr .lt. 0 )  then
   if( ionode ) WRITE ( stdout, '(a)') 'ERROR reading input_file : opttag section is absent'
   STOP
  elseif( ioerr .gt. 0 )  then
   if( ionode ) WRITE ( stdout, '(a)') 'ERROR reading input_file : opttag wrong tag'
   STOP
  endif

  CLOSE  ( stdin )
  ! ======================
  !  check readed opttag
  ! ======================
  CALL opt_check_tag
  ! ======================
  !  print info on opttag
  ! ======================
  CALL opt_print_info(stdout)
  CALL opt_print_info(kunit_OUTFF)

  return

END SUBROUTINE opt_init



!*********************** SUBROUTINE opt_default_tag ***************************
!
! set default values to vacf tag
!
!******************************************************************************

SUBROUTINE opt_default_tag

  implicit none

  ! ================
  !  default values
  ! ================
  nmaxopt      = 1 
  nskipopt     = 0 
  optalgo      = 'step' 
  ncopt        = 0

  return 
 
END SUBROUTINE opt_default_tag


!*********************** SUBROUTINE opt_check_tag ***************************
!
! check opt tag values
!
!******************************************************************************

SUBROUTINE opt_check_tag

  USE io_file,  ONLY :  ionode , stdout

  implicit none

  ! local
  logical :: allowed
  integer :: i

  do i = 1 , size ( optalgo_allowed )
    if( TRIM(optalgo) == optalgo_allowed(i) ) allowed = .TRUE.
  enddo
  if ( .not. allowed ) then
    if ( ionode ) WRITE ( stdout ,'(a,a)') 'ERROR optag: optalgo should be ',optalgo_allowed
    STOP 
  endif

  nperiodopt = (ncopt - nskipopt)/nmaxopt
  if( nperiodopt .lt. 1) nperiodopt = 1

  return

END SUBROUTINE opt_check_tag

!*********************** SUBROUTINE opt_print_info ***************************
!
! print optimisation information to standard output
!
!******************************************************************************

SUBROUTINE opt_print_info(kunit) 

  USE control,  ONLY :  lpbc 
  USE io_file,  ONLY :  ionode 

  implicit none
 
  ! local
  integer :: kunit       

  if( ionode ) then
                               WRITE ( kunit ,'(a)')       '=============================================================' 
                               WRITE ( kunit ,'(a)')       ''
                               WRITE ( kunit ,'(a)')       'minimisation (Inherent Structures)                           ' 
                               WRITE ( kunit ,'(a)')       ''
       if( optalgo .eq. 'sastry' )  then
                               WRITE ( kunit ,'(a)')       'Line mini. via cubic extrapolation of potential '
                               WRITE ( kunit ,'(a)')       'Determines IS via unconstrained optimization    '
                               WRITE ( kunit ,'(a)')       'Original author S. Sastry JNCASR '
                               WRITE ( kunit ,'(a)')       'readapted by FMV'
                               WRITE ( kunit ,'(a)')       ''
                               WRITE ( kunit ,'(a)')       'Steppest Method '
       endif
       if( optalgo .eq. 'lbfgs' )  then
                               WRITE ( kunit ,'(a)')       'Limited memory BFGS method for large scale optimisation'
                               WRITE ( kunit ,'(a)')       'Author: J. Nocedal   *** July 1990 ***'
                               WRITE ( kunit ,'(a)')       'driver by FMV'
       endif
       if( optalgo .eq. 'm1qn3' )  then
                               WRITE ( kunit ,'(a)')       'M1QN3, Version 3.3, October 2009'
                               WRITE ( kunit ,'(a)')       'Authors: Jean Charles Gilbert, Claude Lemarechal, INRIA.'
                               WRITE ( kunit ,'(a)')       'Copyright 2008, 2009, INRIA.'
                               WRITE ( kunit ,'(a)')       'M1QN3 is distributed under the terms of the GNU General Public  License.'
                               WRITE ( kunit ,'(a)')       ''
                               WRITE ( kunit ,'(a)')       'reverse communication and DIS ( see m1qn3 driver and documentation) ' 
       endif

                               WRITE ( kunit ,'(a)')       ''  
     if(.not.lpbc)             WRITE ( kunit ,'(a)')       'no periodic boundary conditions (cubic box wall)'
     if(lpbc)                  WRITE ( kunit ,'(a)')       'periodic boundary conditions in cubic cell'
                               WRITE ( kunit ,'(a)')       ''
                               WRITE ( kunit ,'(a)')       'read configuration from file          :   TRAJFF'
                               WRITE ( kunit ,'(a)')       'save IS thermo. properties into file  :   ISTHFF' 
                               WRITE ( kunit ,'(a)')       'save IS configurations into file      :   ISCFF '  
                               WRITE ( kunit ,'(a)')       ''
                               WRITE ( kunit ,'(a,i5)')    'number of config. in TRAJFF         =     ',ncopt           
                               WRITE ( kunit ,'(a,i5)')    'number of config. to be skipped     =     ',nskipopt
                               WRITE ( kunit ,'(a,i5)')    'maximum of config. to be minimized  =     ',nmaxopt           
                               WRITE ( kunit ,'(a,i5)')    'minimized every config.             =     ',nperiodopt
                               WRITE ( kunit ,'(a)')       ''
                               WRITE ( kunit ,'(a)')       '=============================================================' 
  endif 


  return

END SUBROUTINE opt_print_info


!*********************** SUBROUTINE opt_main  *********************************
!
! main subroutine for optimisation
!
!******************************************************************************

SUBROUTINE opt_main 

  USE config,           ONLY :  system , natm , ntype , rx , ry , rz , fx , fy , fz , atype , rho , config_alloc , list , point , box , omega , atypei , itype, natmi
  USE control,          ONLY :  myrank , numprocs
  USE io_file,          ONLY :  ionode , stdout , kunit_TRAJFF , kunit_ISTHFF , kunit_ISCFF, kunit_OUTFF
  USE thermodynamic,    ONLY :  u_tot , pressure_tot , calc_thermo
  USE constants,        ONLY :  dzero
  USE field 
  USE time

  implicit none
  INCLUDE 'mpif.h'

  ! local 
  integer :: i, ic, neng, iter, na, nopt
  integer :: iastart , iaend , ierr
  double precision :: phigrad, pressure0, pot0, Eis
  double precision :: ttt1,ttt2
  ! trash 
  double precision :: aaaa
  integer :: iiii
  character * 60 :: cccc

  ttt1 = MPI_WTIME(ierr)

  nopt=0
 
  OPEN (UNIT = kunit_TRAJFF ,FILE = 'TRAJFF') 
  OPEN (UNIT = kunit_ISTHFF ,FILE = 'ISTHFF')
  OPEN (UNIT = kunit_ISCFF  ,FILE = 'ISCFF') 

    
  if( ionode ) WRITE ( kunit_ISTHFF , '(a)' )                '#neng: evaluation of force'
  if( ionode ) WRITE ( kunit_ISTHFF , '(a8,3a20,2a6,2a18)') "# config","eIS","grad","Pres","Iter","Neng","u_initial","Press_initial"
       
  READ ( kunit_TRAJFF , * ) natm
  READ ( kunit_TRAJFF , * ) system
  READ ( kunit_TRAJFF , * ) box , ntype
  omega = box * box * box
  rho = dble ( natm )  / omega 
 if( ionode ) then
    WRITE ( stdout , '(a)'      )    'Remind some parameters of the system:'
    WRITE ( stdout , '(a,i12)'  )    'natm  = ',natm
    WRITE ( stdout , '(a,i12)'  )    'ntype = ',ntype
    WRITE ( stdout , '(a,f12.5)')    'rho   = ',rho
    WRITE ( stdout , '(a,f12.5)')    'box   = ',box
    WRITE ( stdout , '(a,f12.5)')    'vol   = ',omega
    WRITE ( stdout , '(a)'      )    ''
  endif

  ! ===================================
  !  here we know natm, then alloc 
  !  and decomposition can be applied 
  ! ================================== 
  CALL config_alloc 
  CALL do_split ( natm , myrank , numprocs , iastart , iaend )
  CALL field_init
#ifdef debug
  print*,iastart , iaend
#endif
  ! ==========================================   
  !  skip the first nskipopt configurations 
  ! ==========================================
  if(nskipopt.gt.0) then
    do ic = 1,nskipopt
      if(ic.ne.1 ) READ ( kunit_TRAJFF , * ) iiii   
      if(ic.ne.1 ) READ ( kunit_TRAJFF , * ) cccc
      if(ic.ne.1 ) READ ( kunit_TRAJFF , * ) aaaa   ,iiii
      do i = 1,natm 
        READ ( kunit_TRAJFF , * ) atype(i),rx(i),ry(i),rz(i),aaaa,aaaa,aaaa,aaaa,aaaa,aaaa
      enddo      
    enddo
  endif

  do ic = nskipopt + 1, ncopt

    ! ===================================
    !  read config from trajectory file
    ! ===================================
    if( ic .ne. (nskipopt + 1) .or. nskipopt.ne.0) READ ( kunit_TRAJFF , * ) iiii
    if( ic .ne. (nskipopt + 1) .or. nskipopt.ne.0) READ ( kunit_TRAJFF , * ) cccc
    if( ic .ne. (nskipopt + 1) .or. nskipopt.ne.0) READ ( kunit_TRAJFF , * ) aaaa , iiii
    do i = 1,natm
      READ ( kunit_TRAJFF , * ) atype(i),rx(i),ry(i),rz(i),aaaa,aaaa,aaaa,aaaa,aaaa,aaaa
      if(atype(i) .eq. 'A') na = na + 1
      if(atype(i) .eq. 'A') itype(i)=1
      if(atype(i) .eq. 'B') itype(i)=2
    enddo
    if ( ntype .eq.  1 ) then
      natmi(0)=natm
      natmi(1)=na
      atypei(0)='ALL'
      atypei(1)='A'
    endif
    if ( ntype .eq.  2 ) then
      natmi(0)=natm
      natmi(1)=na
      natmi(2)=natm-na
      atypei(0)='ALL'
      atypei(1)='A'
      atypei(2)='B'
    endif

    if( ((mod(ic,nperiodopt) .eq. 0) .or. (ic .eq. nskipopt + 1) ) .and. nopt.lt.nmaxopt ) then 
      nopt=nopt+1 
      ! =======================
      !  calc initial thermo  
      ! =======================  
      neng = 0
      CALL engforce( iastart , iaend )
      CALL calc_thermo
      pot0 = u_tot 
      pressure0 = pressure_tot

      ! ====================================
      !  main routine for the optimisation
      ! ====================================
      if(optalgo.eq.'sastry') then 
        if( ionode ) then
          WRITE ( stdout ,'(a)')             ''
          WRITE ( stdout ,'(a,2f16.8)')      'initial energy&pressure  = ',pot0,pressure0
          WRITE ( stdout ,'(a)')             '    its  nstep          grad              ener'
          WRITE (kunit_OUTFF,'(a,2f16.8)')   'initial energy&pressure  = ',pot0,pressure0
        endif
        CALL sastry ( iter , Eis , phigrad , neng , iastart , iaend )
      endif


      if(optalgo.eq.'lbfgs') then
        if( ionode ) then
          WRITE ( stdout ,'(a)')             ''
          WRITE ( stdout ,'(a,2f16.8)')      'initial energy&pressure  = ',pot0,pressure0
          WRITE ( stdout ,'(a)')             '    its       grad              ener'
          WRITE (kunit_OUTFF,'(a,2f16.8)')   'initial energy&pressure  = ',pot0,pressure0
        endif
        CALL lbfgs_driver ( iter, Eis , phigrad , iastart , iaend )
        neng = iter ! the number of function call is iter
      endif


       if(optalgo.eq.'m1qn3') then
        if( ionode ) then
          WRITE ( stdout ,'(a)')             ''
          WRITE ( stdout ,'(a,2f16.8)')      'initial energy&pressure  =',pot0,pressure0
          WRITE ( stdout ,'(a)')             '    its       grad    ener'
          WRITE (kunit_OUTFF,'(a,2f16.8)')   'initial energy&pressure  =',pot0,pressure0
        endif
        CALL m1qn3_driver ( iter, Eis , phigrad , iastart , iaend )
        neng = iter ! the number of function call is iter
      endif

      CALL calc_thermo  
      ! ================================
      !  write final thermodynamic info 
      ! ================================
      if( ionode ) then
        WRITE ( stdout , '(a,2f16.8)' )      'final energy&Pressure = ',Eis,pressure_tot
        WRITE ( kunit_OUTFF , '(a,2f16.8)' ) 'final energy&Pressure = ',Eis,pressure_tot
        WRITE ( kunit_ISTHFF , '(i8,3f20.12,2i8,2f20.12)' ) ic , Eis , phigrad , pressure_tot , iter , neng , pot0 , pressure0
        WRITE ( stdout,'(a)') ''
        WRITE ( stdout,'(a)') ''  
      endif
      
      ! ===========================================
      !  calculated forces (they should be small) 
      ! ===========================================
      CALL engforce( iastart , iaend )
      ! =============================================
      !  write IS structures
      !  new configurations are stored in rx ,ry ,rz, fx , fy ,fz
      ! =============================================         
      if( ionode ) then   
        WRITE ( kunit_ISCFF , * ) natm 
        WRITE ( kunit_ISCFF , * ) system 
        WRITE ( kunit_ISCFF , * ) box,ntype 
        do i = 1,natm 
          WRITE ( kunit_ISCFF ,'(A2,9F20.12)') atype(i),rx(i),ry(i),rz(i),dzero,dzero,dzero,fx(i),fy(i),fz(i)
        enddo       
      endif

    endif ! nperiod         
  enddo !ncopt 

  CLOSE( kunit_ISCFF )
  CLOSE( kunit_ISTHFF )
  CLOSE( kunit_TRAJFF )

  ttt2 = MPI_WTIME(ierr)
  opttimetot = opttimetot + (ttt2-ttt1) 


  return

END SUBROUTINE opt_main


!*********************** SUBROUTINE satry *************************************
!
! Line minimizations via cubic extrapolation of potential. 
! determines inherent structure via unconstrained optimization
!
!******************************************************************************

SUBROUTINE sastry ( iter , Eis , phigrad , neng , iastart , iaend )

  USE config,           ONLY :  natm , rx , ry , rz , fx , fy , fz , list , point
  USE io_file,          ONLY :  ionode , stdout
  USE thermodynamic,    ONLY :  u_tot      
  USE field 

  implicit none

  ! global
  integer :: iastart , iaend 
  integer :: iter, neng
  double precision :: Eis, phigrad, vir
  
  ! local 
  integer :: i , j , kl
  integer :: itmax , nskp , ik , its , nstep
  double precision ftol , epsilon
  double precision umag ,ds , dsk , dutol , ukeep , x0 , u0 , f1d0
  double precision x1 , u1 , f1d1 , uprev , dx , dx2 , dx3 , xsol , f1ds
  double precision uppcub , u3pcub , x1sol , x2sol , curv1 , curv2 , usol , u2
  double precision u2pdis , adx , dd1 , dd2 , gg , dgg , gam

  double precision, dimension(:), allocatable :: gx  , gy  , gz
  double precision, dimension(:), allocatable :: hx  , hy  , hz
  double precision, dimension(:), allocatable :: xix , xiy , xiz
  double precision, dimension(:), allocatable :: fx_sum, fy_sum, fz_sum
 
  allocate( gx(natm) , gy(natm)  ,gz(natm) )
  allocate( hx(natm) , hy(natm) , hz(natm) )
  allocate( xix(natm) , xiy(natm) , xiz(natm) ) 
  allocate( fx_sum(natm), fy_sum(natm), fz_sum(natm) )

  
  ftol = 1.0E-14
  epsilon = 1.0E-14
  itmax = 10000
  nskp = 1
  nstep = 0
  CALL engforce ( iastart , iaend )!, list , point )
  neng = neng + 1

  umag = 0.0D0
  do j = 1,natm
    xix(j) = - fx(j)
    xiy(j) = - fy(j)
    xiz(j) = - fz(j)
    gx(j) = - xix(j)          
    gy(j) = - xiy(j)
    gz(j) = - xiz(j)
    hx(j) = gx(j)
    hy(j) = gy(j)
    hz(j) = gz(j)
    xix(j) = hx(j)
    xiy(j) = hy(j)
    xiz(j) = hz(j)
    umag = umag + xix(j) * xix(j) + xiy(j) * xiy(j) + xiz(j) * xiz(j)
  enddo

  umag = dsqrt(umag)
  do j = 1, natm
    xix(j) = xix(j)/umag
    xiy(j) = xiy(j)/umag
    xiz(j) = xiz(j)/umag
  END do

  ik = 0
  ds = 0.2D0
  dsk = ds
  dutol = 1.0e-12

  kl = 1
  do its = 1,itmax

    phigrad = 0.0D0
    do i = 1, natm
      phigrad = phigrad + fx(i) * fx(i) + fy(i) * fy(i) + fz(i) * fz(i)
    enddo!

    if( ionode .and.mod(its,kl).eq.0) then
      WRITE (stdout,'(2i6,2E20.8,i6)') its , nstep , phigrad , u_tot
      kl = 2 * kl 
    endif
#ifdef debug
  call print_config_sample(its,0)
#endif

    iter = its
    ukeep = u_tot
    x0 = 0.0D0
    u0 = u_tot
    f1d0 = 0.0D0
    do i = 1, natm
      f1d0 = f1d0 + fx(i) * xix(i) + fy(i) * xiy(i) + fz(i) * xiz(i)
    enddo
    f1d0 = - f1d0


!C  reset search direction to steepest descent direction 
    if(f1d0.ge.0.0) then
      if( ionode ) WRITE (stdout, * ) 'RESETTING SEARCH DIRECTION'
      umag = 0.0D0
      do j = 1,natm
        xix(j) = - fx(j)
        xiy(j) = - fy(j)
        xiz(j) = - fz(j)
        gx(j) = - xix(j)          
        gy(j) = - xiy(j)
        gz(j) = - xiz(j)
        hx(j) = gx(j)
        hy(j) = gy(j)
        hz(j) = gz(j)
        xix(j) = hx(j)
        xiy(j) = hy(j)
        xiz(j) = hz(j)
        umag = umag + xix(j) * xix(j) + xiy(j) * xiy(j) + xiz(j) * xiz(j)
      enddo
      umag = dsqrt(umag)
      do j = 1, natm
        xix(j) = xix(j)/umag
        xiy(j) = xiy(j)/umag
        xiz(j) = xiz(j)/umag
      enddo
      ukeep = u_tot
      x0 = 0.0D0
      u0 = u_tot
      f1d0 = 0.0D0
      do i = 1, natm
        f1d0 = f1d0 + fx(i) * xix(i) + fy(i) * xiy(i) + fz(i) * xiz(i)
      enddo
      f1d0 = - f1d0
    endif

    x1 = ds 

    CALL eforce1d ( x1 , u1 , vir , iastart , iaend , f1d1 , xix , xiy , xiz , neng )
    nstep = 0
    uprev = min( u1 , u0 )
777 nstep = nstep + 1

    dx = x1 - x0
    dx2 = dx * dx
    dx3 = dx * dx2

    if(nstep.gt.100) then !nstep

      if( min ( u0 , u1 ) .gt. ukeep ) then
        Eis = u0
        deallocate( fx_sum, fy_sum, fz_sum )
        return
      else 
        if( u0 .lt. u1 ) then
          CALL eforce1d ( x0 , u0 , vir , iastart , iaend , f1d0 , xix , xiy , xiz , neng )
          xsol = x0 
          u_tot = u0 
          f1ds = f1d0
        else 
          xsol = x1
          u_tot = u1
          f1ds = f1d1 
        endif

        do i = 1, natm
          rx(i) = rx(i) + xsol * xix(i)
          ry(i) = ry(i) + xsol * xiy(i)
          rz(i) = rz(i) + xsol * xiz(i)
        enddo

      endif
      goto 676
    endif


    if(f1d0 * f1d1.lt.0.0D0) then !f1d0 * f1d1.lt.0.0D0
      uppcub = (2.0D0/dx2) * (3.0D0 * (u1 - u0) - (f1d1 + 2.0D0 * f1d0) * dx)
      u3pcub = ( - 12.0D0/(dx3)) * ((u1 - u0) - (f1d1 + f1d0) * (dx/2.0D0))

      x1sol = - uppcub + dsqrt(uppcub * uppcub - 2.0D0 * u3pcub * f1d0)
      x1sol = x1sol/u3pcub
      x2sol = - uppcub -  dsqrt(uppcub * uppcub - 2.0D0 * u3pcub * f1d0)
      x2sol = x2sol/u3pcub
              
      curv1 = uppcub + u3pcub * x1sol 
      curv2 = uppcub + u3pcub * x2sol 

      if(curv1.gt.0.0D0) xsol = x1sol + x0
      if(curv2.gt.0.0D0) xsol = x2sol + x0

      dx = xsol - x0 
      dx2 = dx * dx
      dx3 = dx * dx2

      usol = u0 + f1d0 * dx + (uppcub * dx2/2.0D0) + (u3pcub * dx3/6.0D0)

!C    this is to make sure u2 is not > ukeep
524   CALL eforce1d ( xsol , u2 , vir , iastart , iaend , f1ds , xix , xiy , xiz , neng ) 


      if(u2.gt.ukeep) then
        xsol = x0 + (xsol - x0)/2.0D0
        goto 524
      endif

      if((abs((uprev - u2)/u2)).lt.dutol) then

        do i = 1, natm
          rx(i) = rx(i) + xsol * xix(i)
          ry(i) = ry(i) + xsol * xiy(i)
          rz(i) = rz(i) + xsol * xiz(i)
        enddo

        u_tot = u2 
        ds = xsol
        if(xsol.le.1.0e-8) ds = 1.0e-8
        goto 676
      else 

        if(f1d0 * f1ds.le.0.0D0) then
          x1 = xsol 
          u1 = u2 
          f1d1 = f1ds
        else
          x0 = xsol 
          u0 = u2
          f1d0 = f1ds 
        endif
        uprev = min(u0,u1)
        goto 777
      endif
              
    else !f1d0 * f1d1.lt.0.0D0

      u2pdis = (f1d1 - f1d0)/(x1 - x0)
      adx = abs(x1 - x0)

      if(u2pdis.gt.0.0D0) then
        xsol = x0 -  f1d0/u2pdis
        if(abs(xsol - x0).gt.3.0D0 * adx)then
          xsol = x0 + 3.0D0 * adx
        endif
      else
!c      xsol = x0 -  (f1d0/abs(f1d0)) * 1.50D0 * adx
        xsol = x0 + 1.50D0 * adx
      endif

!C    upto here define xsol 

!C    now make sure usol < ukeep
525   CALL eforce1d ( xsol , u2 , vir , iastart , iaend , f1ds , xix , xiy , xiz , neng )


      if(u2.gt.ukeep)then
        xsol = x0 + (xsol - x0)/2.0D0
        goto 525
      endif


      if((abs((uprev - u2)/u2).lt.dutol)) then

        do i = 1, natm
          rx(i) = rx(i) + xsol * xix(i)
          ry(i) = ry(i) + xsol * xiy(i)
          rz(i) = rz(i) + xsol * xiz(i)
        END do

        u_tot = u2
        ds = xsol

        if(xsol.le.1.0e-8) ds = 1.0e-8

        goto 676
      END if

!C    next pick x0 to be lowest x with negative f 

      if(u1.lt.u0)then
        x0 = x1
        u0 = u1
        f1d0 = f1d1
      END if

      x1 = xsol
      u1 = u2
      f1d1 = f1ds 

      uprev = min(u0,u1)

      goto 777 
      
    endif
              
676 continue 
           

    DD1 = 2.0D0 * abs(ukeep - u_tot) 
    DD2 = ftol * (abs(ukeep) + abs(u_tot) + epsilon)

    if ( DD1 .LE. DD2) then 
      CALL engforce ( iastart , iaend )!, list , point )
      neng = neng + 1

      Eis = u_tot
      phigrad = 0.0D0
      do i = 1, natm
        phigrad = phigrad + fx(i) * fx(i) + fy(i) * fy(i) + fz(i) * fz(i)
      enddo

      if( ionode ) then
        WRITE ( stdout,'(2i6,2E20.8,i6)') its,nstep,phigrad,u_tot
        WRITE ( stdout , '(a,i5,a)' ) 'minimum reached in ',its,' iterations'
      endif
      deallocate( fx_sum, fy_sum, fz_sum )
      return
    else
      if(abs(1.0D0 - u_tot/ukeep).lt.dutol)then
        dutol = abs(1.0D0 - u_tot/ukeep)
      endif
    endif
           
    gg = 0.0D0
    dgg = 0.0D0
    do j = 1,natm
      xix(j) = - fx(j)
      xiy(j) = - fy(j)
      xiz(j) = - fz(j)
      gg = gg + gx(j) * gx(j) + gy(j) * gy(j) + gz(j) * gz(j)
      dgg = dgg + (xix(j) + gx(j)) * xix(j) +  &
              (xiy(j) + gy(j)) * xiy(j) +  &
              (xiz(j) + gz(j)) * xiz(j)
    enddo
    if (gg .eq. 0.0D0) then 
      deallocate( fx_sum, fy_sum, fz_sum )
      return
    endif
    gam = dgg/gg
    umag = 0.0D0
    do j = 1,natm
      gx(j) = - xix(j)
      gy(j) = - xiy(j)
      gz(j) = - xiz(j)
      hx(j) = gx(j) + gam * hx(j)
      hy(j) = gy(j) + gam * hy(j)
      hz(j) = gz(j) + gam * hz(j)
      xix(j) = hx(j)
      xiy(j) = hy(j)
      xiz(j) = hz(j)
      umag = umag + xix(j) * xix(j) + xiy(j) * xiy(j) + xiz(j) * xiz(j)
    enddo
           
    umag = dsqrt(umag)
    do j = 1, natm
      xix(j) = xix(j)/umag
      xiy(j) = xiy(j)/umag
      xiz(j) = xiz(j)/umag
    enddo

  enddo 


  deallocate( fx_sum, fy_sum, fz_sum )
  deallocate( gx  , gy  , gz )
  deallocate( hx  , hy  , hz )
  deallocate( xix , xiy , xiz )


  return

END SUBROUTINE sastry     


!*********************** SUBROUTINE eforce1d **********************************
!
! one dimensional minimal search for the Sastry otpmisation routine (sastry)
!
!******************************************************************************

SUBROUTINE eforce1d( x , pot , vir , iastart , iaend , f1d , xix , xiy , xiz , neng ) 

  USE config,           ONLY :  natm, box, rx, ry, rz, fx, fy, fz
  USE thermodynamic,    ONLY : vir_tot , u_tot , calc_thermo
  USE field 

  implicit none

  ! global
  double precision :: pot ,vir, x, f1d
  integer          :: iastart , iaend
!  integer          :: list(natm * 250), point(natm + 1)
  double precision :: xix(natm),xiy(natm),xiz(natm)
  
  ! local
  integer :: i
  integer :: neng
  double precision, dimension(:), allocatable :: fx_sum, fy_sum, fz_sum
  double precision, dimension(:), allocatable :: rxt, ryt, rzt
  double precision, dimension(:), allocatable :: tmpx, tmpy ,tmpz

  allocate( fx_sum(natm), fy_sum(natm), fz_sum(natm) )
  allocate( rxt(natm), ryt(natm), rzt(natm) )
  allocate( tmpx(natm), tmpy(natm) ,tmpz(natm) )

  fx_sum = 0.0D0
  fy_sum = 0.0D0
  fz_sum = 0.0D0
  rxt = 0.0D0
  ryt = 0.0D0
  rzt = 0.0D0
  tmpx = 0.0D0
  tmpy = 0.0D0
  tmpz = 0.0D0

  ! ===============
  !  search in 1D
  ! ===============
  do i = 1, natm
    rxt(i) = rx(i) + x * xix(i)
    ryt(i) = ry(i) + x * xiy(i)
    rzt(i) = rz(i) + x * xiz(i)
  END do
   
  tmpx = rx
  tmpy = ry
  tmpz = rz
  rx = rxt
  ry = ryt
  rz = rzt

  ! ================
  !  warning ! only pbc
  ! ================
  CALL engforce ( iastart , iaend )!, list , point )
  CALL calc_thermo
  vir = vir_tot
  pot = u_tot  
  neng = neng + 1

  rx = tmpx
  ry = tmpy
  rz = tmpz

  f1d = 0.0D0
  do i = 1, natm
    f1d = f1d + fx(i) * xix(i) + fy(i) * xiy(i) + fz(i) * xiz(i)
  enddo
  f1d = - f1d


  deallocate( fx_sum, fy_sum, fz_sum )
  deallocate( rxt, ryt, rzt )
  deallocate( tmpx, tmpy, tmpz )


  return

END SUBROUTINE eforce1d 

!*********************** SUBROUTINE lbfgs_driver *******************************
!
! this subroutine is a driver for tje lbfgs subroutine (see file lbfgs.f90)
!
! Description:
!
!        LIMITED MEMORY BFGS METHOD FOR LARGE SCALE OPTIMIZATION
!                          JORGE NOCEDAL
!                        *** July 1990 ***
!
! 
!     This subroutine solves the unconstrained minimization problem
! 
!                      min F(x),    x= (x1,x2,...,xN),
!
!      using the limited memory BFGS method. The routine is especially
!      effective on problems involving a large number of variables. In
!      a typical iteration of this method an approximation Hk to the
!      inverse of the Hessian is obtained by applying M BFGS updates to
!      a diagonal matrix Hk0, using information from the previous M steps.
!      The user specifies the number M, which determines the amount of
!      storage required by the routine. The user may also provide the
!      diagonal matrices Hk0 if not satisfied with the default choice.
!      The algorithm is described in "On the limited memory BFGS method
!      for large scale optimization", by D. Liu and J. Nocedal,
!      Mathematical Programming B 45 (1989) 503-528.
! 
!      the routine CSRCH written by More' and Thuente.
!
!******************************************************************************

SUBROUTINE lbfgs_driver ( icall, Eis , phigrad , iastart , iaend  )

  USE config,           ONLY :  natm, rx, ry, rz , fx , fy , fz , list, point
  USE thermodynamic,    ONLY :  u_tot , u_lj_r , calc_thermo
  USE field
  USE io_file,          ONLY :  ionode , stdout

  implicit none

  ! global
  integer :: iastart , iaend , icall
  double precision :: Eis, phigrad

  ! local 
  integer :: ia , kl , itmax , its
  integer :: NWORK
  double precision :: F,EPS,XTOL,GTOL,STPMIN,STPMAX
  integer :: IPRINT(2),IFLAG,N,M
  logical :: DIAGCO
  double precision, dimension (:), allocatable :: X , G , DIAG , W

  !==============================================================
  ! The driver for LBFGS must always declare LB2 as EXTERNAL
  !==============================================================
  external LB2
  common  /LB3/ GTOL,STPMIN,STPMAX

  !=====================================
  !  dimension of optimization problem
  !=====================================
  N=3*natm
  ALLOCATE ( X (N) ) 
  ALLOCATE ( G (N) )
  ALLOCATE ( DIAG (N) )
  !=================================================================================
  !  M is an INTEGER variable that must be set by the user to
  !  the number of corrections used in the BFGS update. It
  !  is not altered by the routine. Values of M less than 3 are
  !  not recommended; large values of M will result in excessive
  !  computing time. 3<= M <=7 is recommended. Restriction: M>0
  !=================================================================================
  M=4
  NWORK = N * ( 2 * M + 1 ) + 2*M
  ALLOCATE ( W (NWORK) )
  !=================================================================================
  !  IPRINT  is an INTEGER array of length two which must be set by the
  !             user.
  ! 
  !             IPRINT(1) specifies the frequency of the output:
  !                IPRINT(1) < 0 : no output is generated,
  !                IPRINT(1) = 0 : output only at first and last iteration,
  !                IPRINT(1) > 0 : output every IPRINT(1) iterations.
  ! 
  !             IPRINT(2) specifies the type of output generated:
  !                IPRINT(2) = 0 : iteration count, number of function 
  !                                evaluations, function value, norm of the
  !                                gradient, and steplength,
  !                IPRINT(2) = 1 : same as IPRINT(2)=0, plus vector of
  !                                variables and  gradient vector at the
  !                                initial point,
  !                IPRINT(2) = 2 : same as IPRINT(2)=1, plus vector of
  !                                variables,
  !                IPRINT(2) = 3 : same as IPRINT(2)=2, plus gradient vector.
  !=================================================================================
  IPRINT(1)= -1
  IPRINT(2)= 0

  !=================================================================================
  ! We do not wish to provide the diagonal matrices Hk0, and 
  ! therefore set DIAGCO to FALSE.
  !=================================================================================
  DIAGCO= .FALSE.
  !=================================================================================
  ! EPS     is a positive DOUBLE PRECISION variable that must be set by
  !         the user, and determines the accuracy with which the solution
  !         is to be found. The subroutine terminates when
  !
  !             ||G|| < EPS max(1,||X||),
  !
  !             where ||.|| denotes the Euclidean norm.
  !=================================================================================
  EPS = 1.0D-5
  !=================================================================================
  ! XTOL    is a  positive DOUBLE PRECISION variable that must be set by
  !         the user to an estimate of the machine precision (e.g.
  !         10**(-16) on a SUN station 3/60). The line search routine will
  !         terminate if the relative width of the interval of uncertainty
  !         is less than XTOL.
  !=================================================================================
  XTOL = 1.0D-14
  
  ICALL = 0
  !=================================================================================
  !  IFLAG   is an INTEGER variable that must be set to 0 on initial entry
  !             to the subroutine. A return with IFLAG<0 indicates an error,
  !             and IFLAG=0 indicates that the routine has terminated without
  !             detecting errors. On a return with IFLAG=1, the user must
  !             evaluate the function F and gradient G. On a return with
  !             IFLAG=2, the user must provide the diagonal matrix Hk0.
  ! 
  !             The following negative values of IFLAG, detecting an error,
  !             are possible:
  ! 
  !              IFLAG=-1  The line search routine MCSRCH failed. The
  !                        parameter INFO provides more detailed information
  !                        (see also the documentation of MCSRCH):
  !
  !                       INFO = 0  IMPROPER INPUT PARAMETERS.
  !
  !                       INFO = 2  RELATIVE WIDTH OF THE INTERVAL OF
  !                                 UNCERTAINTY IS AT MOST XTOL.
  !
  !                       INFO = 3  MORE THAN 20 FUNCTION EVALUATIONS WERE
  !                                 REQUIRED AT THE PRESENT ITERATION.
  !
  !                       INFO = 4  THE STEP IS TOO SMALL.
  !
  !                       INFO = 5  THE STEP IS TOO LARGE.
  !
  !                       INFO = 6  ROUNDING ERRORS PREVENT FURTHER PROGRESS. 
  !                                 THERE MAY NOT BE A STEP WHICH SATISFIES
  !                                 THE SUFFICIENT DECREASE AND CURVATURE
  !                                 CONDITIONS. TOLERANCES MAY BE TOO SMALL.
  !
  ! 
  !              IFLAG=-2  The i-th diagonal element of the diagonal inverse
  !                        Hessian approximation, given in DIAG, is not
  !                        positive.
  !           
  !              IFLAG=-3  Improper input parameters for LBFGS (N or M are
  !C                        not positive).
  !=================================================================================
  IFLAG=0

  ! ========================================
  !  set X ( 3 *natm ) to rx , ry , rz
  ! ========================================
  do ia = 1, natm 
    X( ia )            = rx( ia ) 
    X( natm + ia )     = ry( ia ) 
    X( 2 * natm + ia ) = rz( ia ) 
  enddo  

  kl = 1
  itmax = 10000

  !==============================================================
  ! We allow at most 10000 evaluations of F and G
  !==============================================================
  do its=1,itmax

    CALL engforce ( iastart , iaend )!, list , point )
    CALL calc_thermo

#ifdef debug
     call print_config_sample(icall,0)
#endif

    ! ===============================
    !  set the gradient to fx,fy,fz
    ! ===============================
    DO ia=1,natm
      G( ia           ) =   - fx( ia ) 
      G( natm + ia    ) =   - fy( ia ) 
      G( 2* natm + ia ) =   - fz( ia ) 
    ENDDO

    F=u_tot

    ! ====================
    !  main call to LBFGS
    ! ====================
    ICALL=ICALL + 1
    CALL LBFGS(N,M,X,F,G,DIAGCO,DIAG,IPRINT,EPS,XTOL,W,IFLAG)
    IF(IFLAG.LE.0) THEN 
      Eis=u_tot 
      IF(IFLAG.EQ.0) THEN
        if( ionode ) then
          WRITE ( stdout,'(i6,2E20.8,i6)') icall , phigrad , u_tot
          WRITE ( stdout , '(a,i5,a)' ) 'minimum reached in ',icall,' iterations'
        endif
      ENDIF
      RETURN 
    ENDIF
    ! ===============================================
    !  reset positions for energy/forces calculation
    ! ==============================================
    phigrad=0.d0
    do ia = 1, natm  
      rx( ia )   = X ( ia )
      ry( ia )   = X ( natm + ia )
      rz( ia )   = X ( 2* natm + ia )
      phigrad=phigrad+ G ( ia ) * G( ia ) + G( natm + ia    ) * G( natm + ia    )  + G( 2* natm + ia ) * G( 2* natm + ia )
    enddo  

    if( ionode .and.mod(icall,kl).eq.0) then
      WRITE (stdout,'(i6,2E20.8,i6)') icall , phigrad , u_tot
      kl = 2 * kl
    endif

  enddo 


  DEALLOCATE ( X  )
  DEALLOCATE ( G  )
  DEALLOCATE ( DIAG )
  DEALLOCATE ( W )



  RETURN 

END SUBROUTINE lbfgs_driver


!*********************** SUBROUTINE  m1qn3_driver *****************************
!
! driver routine for the optimsation algorithm m1qn3 (see m1qn3.f90)
!
!******************************************************************************
SUBROUTINE m1qn3_driver ( icall, Eis , phigrad, iastart , iaend )

  USE control,          ONLY :  myrank , numprocs
  USE config,           ONLY :  natm , rx , ry , rz, fx ,fy ,fz, list ,point
  USE io_file,          ONLY :  stdout, ionode
  USE thermodynamic,    ONLY :  u_tot , u_lj_r , calc_thermo
  USE field

  implicit none

  ! global
  integer :: iastart , iaend , icall
  double precision :: Eis , phigrad

  !local
  integer :: ia , n , kl 
  character*3 normtype
  integer :: imp,io,imode(3),omode,niter,nsim,iz(5),ndz,izs(1),indic,reverse
  real rzs(1)
  double precision :: f,dx,df1,epsrel,dzs(1)
  double precision, dimension (:), allocatable :: x , g , dz

  ! external 
  external euclid,ctonbe,ctcabe,simul_rc

  ! =====================
  !   initialization
  ! =====================
  n=3*natm
  ndz=40000
  icall = 0
  kl = 1
  ALLOCATE ( x (N) )
  ALLOCATE ( g (N) )
  ALLOCATE ( dz (ndz) )

  ! ========================================
  !  set X ( 3 *natm ) to rx , ry , rz
  ! ========================================
  do ia = 1, natm
    X( ia )            = rx( ia )
    X( natm + ia )     = ry( ia )
    X( 2 * natm + ia ) = rz( ia )
  enddo

  CALL engforce ( iastart , iaend )
  CALL calc_thermo

  f=u_tot

  ! ===============================
  !  set the gradient to fx,fy,fz
  ! ===============================
  DO ia=1,natm
    G( ia           ) =  - fx( ia )
    G( natm + ia    ) =  - fy( ia )
    G( 2* natm + ia ) =  - fz( ia )
  ENDDO

  ! =========================================================
  !  call the optimization code
  !  normtype can be 'sup' for sup-norm,
  !                  'two' for 2-norm, 
  !                  'dfn' for the norm defined by prosca
  ! =========================================================
  dx=1.e-14            ! dxmin  
  df1=1000.0d0         ! df1
  epsrel=1.e-6         ! epsg 
  niter=6000           ! niter
  nsim=6000            ! nsim 
  normtype = 'dfn'     ! normtype   
  io=stdout            ! io  
  imp= 0               ! impress  no print
  imode(1)=0           ! imode(1) = 0 DIS  
                       ! imode(1) = 1 SIS
  imode(2)=0           ! starting mode : = 0 cold start  direction -g
                       !                 = 1 warm start  direction -Wg 
  imode(3)=0           ! imode 
  reverse=1            ! reverse

  100 continue

    icall = icall + 1

    call m1qn3 (simul_rc,euclid,ctonbe,ctcabe,n,x,f,g,dx,df1, &
                epsrel,normtype,imp,io,imode,omode,niter,nsim,iz, & 
                dz,ndz,reverse,indic,izs,rzs,dzs)
    if (reverse.lt.0) goto 101
    ! ===============================================
    !  reset positions for energy/forces calculation
    ! ==============================================
    do ia = 1, natm
      rx( ia )   = X ( ia )
      ry( ia )   = X ( natm + ia )
      rz( ia )   = X ( 2* natm + ia )
    enddo

    CALL engforce ( iastart , iaend )!, list , point )
    CALL calc_thermo

    f=u_tot
 
    ! ===============================
    !  set the gradient to fx,fy,fz
    ! ===============================
    phigrad=0.d0
    DO ia=1,natm
      G( ia           ) =  - fx( ia )
      G( natm + ia    ) =  - fy( ia )
      G( 2* natm + ia ) =  - fz( ia )
      phigrad=phigrad+ G ( ia ) * G( ia ) + G( natm + ia    ) * G( natm + ia    )  + G( 2* natm + ia ) * G( 2* natm + ia )
    ENDDO

    ! write step information
    if( ionode .and.mod(icall,kl).eq.0) then
      WRITE (stdout,'(i6,2E20.8,i6)') icall , phigrad , u_tot
      kl = 2 * kl
    endif


    goto 100
  101 continue

  if( ionode ) then
    WRITE ( stdout,'(i6,2E20.8,i6)') icall , phigrad , u_tot
    WRITE ( stdout , '(a,i5,a)' ) 'minimum reached in ',icall,' iterations'
  endif
  if( omode .ne. 0 .and. omode .ne. 1 ) then
    if( ionode) WRITE ( stdout , '(a,i5)' ) 'm1qn3 did not properly terminate ',omode
    stop
  endif  

  ! final energy
  Eis=u_tot

  deallocate ( x )
  deallocate ( g )
  deallocate ( dz )

  return

END SUBROUTINE m1qn3_driver


END MODULE opt
! ===== fmV =====
