! ===== fmV =====
!======================================
! driver for fftw routines
!======================================
SUBROUTINE fft_1D_complex(in,out,N)

  implicit none
  INCLUDE "fftw3.f"
  INCLUDE "mpif.h"

  integer :: N
  double complex, dimension(N) :: in, out
  integer plan

  CALL dfftw_plan_dft_1d ( plan , N  , in , out , FFTW_FORWARD , FFTW_ESTIMATE )
  CALL dfftw_execute_dft ( plan , in , out )
  CALL dfftw_destroy_plan( plan )

  return 

END SUBROUTINE fft_1D_complex

SUBROUTINE fft_1D_real(in,out,N)
 
  implicit none
  INCLUDE "fftw3.f"

  integer :: N
  double precision, dimension(N)      :: in
  double complex  ,dimension(N/2 + 1) :: out
  integer*8 plan
     
  call dfftw_plan_dft_r2c_1d(plan,N,in,out,FFTW_ESTIMATE)
  call dfftw_execute_dft_r2c(plan, in, out)
  call dfftw_destroy_plan(plan)

END SUBROUTINE fft_1D_real
! ===== fmV =====
