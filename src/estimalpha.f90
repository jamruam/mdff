PROGRAM estimate_alpha
  
  implicit none

  integer :: k
  double precision :: alpha , alpha2 , rcut , rcut2 , rcut3 , e

  read(*,*) rcut
  alpha = 4.0d0
  alpha2 = alpha*alpha
  rcut2=rcut*rcut
  rcut3=rcut2*rcut

  e = dexp( -alpha2*rcut2 )
  e = e / alpha2 / rcut3
  e = e * 0.56d0  
  print*,alpha,e  
  do while ( e .le. 1e-7 )
    alpha = alpha - 0.05d0
    alpha2 = alpha*alpha
    e = dexp( -alpha2*rcut2 )
    e = e / alpha2 / rcut3
    e = e * 0.56d0
    print*,alpha,e  
  enddo

END PROGRAM estimate_alpha

