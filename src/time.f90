! ===== fmV =====
!*********************** MODULE time ******************************************
!
!  this MODULE implements a small timer utility to time
!  SUBROUTINE 
!  (see START_TIMING)
!  from VASP
!
!******************************************************************************
MODULE time

  implicit none
  
  double precision :: timetot        ! engforce SUBROUTINE       
  double precision :: forcetimetot   ! engforce SUBROUTINE       
  double precision :: fcoultimetot1  ! engforce_coul direct space      
  double precision :: fcoultimetot2  ! engforce_coul fourier space      
  double precision :: efgtimetot1    ! EFG direct space           
  double precision :: efgtimetot2    ! EFG fourier space        
  double precision :: efgtimetot3    ! ditribution  only         
  double precision :: mdsteptimetot  !                            
  double precision :: vacftimetot    !                            
  double precision :: vacftimetot2   !                            
  double precision :: vnlisttimetot  !                            
  double precision :: msdtimetot     !                            
  double precision :: grtimetot      !                            
  double precision :: opttimetot     !                          
  double precision :: vibtimetot     !                           
  double precision :: fvibtimetot    !                           
  double precision :: timetmp        ! use for opt               

CONTAINS


!*********************** SUBROUTINE time_init *********************************
!
! set time to zero
!
!******************************************************************************

SUBROUTINE time_init

  implicit none

  timetot       = 0.0d0
  efgtimetot1   = 0.0d0
  efgtimetot2   = 0.0d0
  efgtimetot3   = 0.0d0
  mdsteptimetot = 0.0d0
  vnlisttimetot = 0.0d0
  msdtimetot    = 0.0d0
  grtimetot     = 0.0d0 
  forcetimetot  = 0.0d0
  fcoultimetot1 = 0.0d0
  fcoultimetot2 = 0.0d0
  vibtimetot    = 0.0d0
  opttimetot    = 0.0d0
  vacftimetot   = 0.0d0                             
  vacftimetot2  = 0.0d0                             

  return 
 
END SUBROUTINE time_init

!*********************** SUBROUTINE print_time_info ***************************
!
! print final time info
!
!******************************************************************************

SUBROUTINE print_time_info ( kunit )  

  USE io_file,  ONLY :  ionode 
  USE control,  ONLY :  longrange

  implicit none

  ! local
  double precision :: dstime , estime
  double precision :: fcoultime_es , fcoultime_ds
  integer :: kunit

  dstime = efgtimetot1 + efgtimetot3               ! direct
  estime = efgtimetot1 + efgtimetot2 + efgtimetot3 ! ewald
  fcoultime_es = fcoultimetot1 + fcoultimetot2
  fcoultime_ds = fcoultimetot1 

  !timing information at the end
  if( ionode ) then
                                    WRITE ( kunit ,110)             'TOTAL'  , timetot
    if( mdsteptimetot .ne. 0.0d0 )  WRITE ( kunit ,110)             'MD'     , mdsteptimetot
    if( opttimetot    .ne. 0.0d0 )  WRITE ( kunit ,110)             'OPT'    , opttimetot
    if( vibtimetot    .ne. 0.0d0 )  WRITE ( kunit ,110)             'VIB'    , vibtimetot
    if( fvibtimetot   .ne. 0.0d0 )  WRITE ( kunit ,110)             'FVIB'   , fvibtimetot
    if( dstime .ne.0.0d0 .and. &
        efgtimetot2   .eq. 0.0d0 )  WRITE ( kunit ,110)             'EFG_DS' , dstime
    if( efgtimetot2   .ne. 0.0d0 )  WRITE ( kunit ,110)             'EFG_ES' , estime
    if( grtimetot     .ne. 0.0d0 )  WRITE ( kunit ,110)             'GR'     , grtimetot 
    if( msdtimetot    .ne. 0.0d0 )  WRITE ( kunit ,110)             'MSD'    , msdtimetot
    if( vacftimetot   .ne. 0.0d0 )  WRITE ( kunit ,110)             'VACF'   , vacftimetot
    if( vacftimetot2   .ne. 0.0d0 )  WRITE ( kunit ,110)            'VACF2'  , vacftimetot2
                                    WRITE ( kunit ,'(a)')           ''
                                    WRITE ( kunit ,'(a)')           'main subroutines:'
                                    WRITE ( kunit ,'(a)')           'MD:'
    if( forcetimetot  .ne. 0.0d0 )  WRITE ( kunit ,110)             'engforce_bmlj     ' , forcetimetot         
    if( fcoultimetot2 .ne. 0.0d0 )  WRITE ( kunit ,110)             'engforce_coul_ES  ' , fcoultime_es         
    if( fcoultimetot2 .eq. 0.0d0 )  WRITE ( kunit ,110)             'engforce_coul_DS  ' , fcoultime_ds         
    if( vnlisttimetot .ne. 0.0d0 )  WRITE ( kunit ,110)             'vnlistcheck       ' , vnlisttimetot
    
    if( longrange .eq. 'direct' )   then
    if( dstime        .ne. 0.0d0 )  WRITE ( kunit ,'(a)')           'EFG:'
    if( efgtimetot1   .ne. 0.0d0 )  WRITE ( kunit ,110)             'efg_DS(efg  only) ', efgtimetot1 
    if( efgtimetot3   .ne. 0.0d0 )  WRITE ( kunit ,110)             'efg_DS(dist only) ', efgtimetot3
    endif
    if( longrange .eq. 'ewald')   then
    if( estime        .ne. 0.0d0 )  WRITE ( kunit ,'(a)')           'EFG:'
    if( efgtimetot1   .ne. 0.0d0 )  WRITE ( kunit ,110)             'efg_ES(real  part)', efgtimetot1
    if( efgtimetot2   .ne. 0.0d0 )  WRITE ( kunit ,110)             'efg_ES(recip part)', efgtimetot2 
    if( efgtimetot3   .ne. 0.0d0 )  WRITE ( kunit ,110)             'efg_ES(dist  only)', efgtimetot3
    if( fcoultime_es  .ne. 0.0d0 )  WRITE ( kunit ,'(a)')           'Coulomb:'
    if( fcoultimetot1 .ne. 0.0d0 )  WRITE ( kunit ,110)             'coul_ES(real  part)', fcoultimetot1
    if( fcoultimetot2 .ne. 0.0d0 )  WRITE ( kunit ,110)             'coul_ES(recip part)', fcoultimetot2 
    endif

    WRITE ( kunit ,'(a)') '============================================================='
  endif

110   FORMAT(2X,A20,' :  cpu time',F9.2)


END SUBROUTINE print_time_info

END MODULE time
! ===== fmV =====
