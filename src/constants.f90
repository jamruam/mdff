! ===== fmV =====
!*********************** MODULE CONSTANTS *************************************
!
! Main constants of the code
!
!******************************************************************************
MODULE constants 
  
  double precision, PARAMETER :: dzero  = 0.0d0
  double precision, PARAMETER :: pi     = 3.14159265358979323846264338327950288419716939937510d0
  double precision, PARAMETER :: tpi    = 2.d0*pi
  double complex  , PARAMETER :: citpi  = (0.d0,1.d0)*tpi
  double complex  , PARAMETER :: imag   = (0.d0,1.d0)
  double precision, PARAMETER :: fpi    = 2.d0*tpi
  double precision, PARAMETER :: piroot = dsqrt(pi)

END MODULE constants
! ===== fmV =====
