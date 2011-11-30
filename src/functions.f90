! ===== fmV =====
! complementary error function
! W. Press et al. : Numerical Recipes, p. 214
! built-in erfc() is buggy e.g. on some Linux distributions
  double precision FUNCTION errfc(x)
  implicit none
  ! local
  double precision :: x , z , t    
  z=dabs(x)
  t=1d0/(1d0+0.5d0*z)
  errfc=t*exp(-z*z-1.26551223d0+t*(1.00002368d0+t*(.37409196d0+ &
  t*(.09678418d0+t*(-.18628806d0+t*(.27886807d0+t*(-1.13520398d0+ &
  t*(1.48851587d0+t*(-.82215223d0+t*.17087277d0)))))))))
  if(x.lt.0d0) errfc=2d0-errfc
  return
  END FUNCTION

  double precision FUNCTION errf(x)
  implicit none
  ! local
  double precision :: x , errfc   
  errf=1.0d0-errfc(x)
  END FUNCTION
! ===== fmV =====
