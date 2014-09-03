MODULE tensors_rk

  USE constants,        ONLY :  dp

  implicit none

  TYPE :: tensor_rank0
    real( kind=dp ) :: sca 
    real( kind=dp ) :: sca_damp      
    real( kind=dp ) :: sca_damp2      
  END TYPE
  TYPE :: tensor_rank1
    real( kind=dp ) :: a(3)      
    real( kind=dp ) :: a_damp(3)      
    real( kind=dp ) :: a_damp2(3)      
  END TYPE
  TYPE :: tensor_rank2
    real( kind=dp ) :: ab(3,3)      
    real( kind=dp ) :: ab_damp(3,3)      
    real( kind=dp ) :: ab_damp2(3,3)      
  END TYPE
  TYPE :: tensor_rank3
    real( kind=dp ) :: abc(3,3,3)      
    real( kind=dp ) :: abc_damp(3,3,3)      
    real( kind=dp ) :: abc_damp2(3,3,3)      
  END TYPE
  TYPE :: tensor_rank4
    real( kind=dp ) :: abcd(3,3,3,3)      
    real( kind=dp ) :: abcd_damp(3,3,3,3)      
    real( kind=dp ) :: abcd_damp2(3,3,3,3)      
  END TYPE
  TYPE :: tensor_rank5
    real( kind=dp ) :: abcde(3,3,3,3,3)      
    real( kind=dp ) :: abcde_damp(3,3,3,3,3)      
    real( kind=dp ) :: abcde_damp2(3,3,3,3,3)      
  END TYPE
 

END MODULE tensors_rk
