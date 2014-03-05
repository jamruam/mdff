#define debug_TT
MODULE tt_damp

  USE constants,        ONLY :  dp
  USE io,               ONLY :  stdout

  integer,            PARAMETER :: maximum_of_TT_expansion=8
  real(kind=dp), dimension (maximum_of_TT_expansion) :: E_TT 

CONTAINS

SUBROUTINE get_TT_damp 

  implicit none

  integer :: k 

  E_TT(0) = 1.0_dp
  do k = 1 , maximum_of_TT_expansion
    E_TT(k) = E_TT(k-1) / REAL ( k ,kind=dp )
  enddo
#ifdef debug_TT  
  write( stdout , '(<maximum_of_TT_expansion>e16.8)') E_TT
#endif

  return

END SUBROUTINE get_TT_damp


END MODULE tt_damp
