program trajff_to_xyz

  implicit none

  integer*8         :: ia , ic , it 
  integer*8         :: nconf , natm , ntype , iiii
  character(len=35) :: xxxx , cpos
  character(len=35) :: filename 
  character(len=2)  :: atype
  real*8            :: box , aaaa
  real*8            :: x , y , z  , A ( 3 , 3 ), V1 ,V2 ,V3 

  READ ( * , * ) nconf
  READ ( * , * ) filename
  WRITE( 0 , '(a,a)' ) 'reading file',filename  
  OPEN ( UNIT = 1000,FILE=filename)
  do ic = 1 , nconf
    WRITE( 0 , '(a,i6,a,i6,a)' ) 'reading [',ic,' / ',nconf,' ] '
    READ ( 1000 , * ) natm
    READ ( 1000 , * ) xxxx
    READ ( kunit_POSFF ,* ) A ( 1 , 1 ) , A ( 2 , 1 ) , A ( 3 , 1 )
    READ ( kunit_POSFF ,* ) A ( 1 , 2 ) , A ( 2 , 2 ) , A ( 3 , 2 )
    READ ( kunit_POSFF ,* ) A ( 1 , 3 ) , A ( 2 , 3 ) , A ( 3 , 3 )

    READ ( 1000 , * ) box ,box , box
    READ ( 1000 , * ) box ,box , box
    READ ( 1000 , * ) box ,box , box
    READ ( 1000 , * ) ntype 
    READ ( 1000 , * ) ( atype , it = 1 , ntype )
    READ ( 1000 , * ) ( iiii , it = 1 , ntype )
    READ ( 1000 , * ) cpos 

    do ia = 1 , natm      
      READ ( 1000 , * ) atype , x , y , z , aaaa , aaaa , aaaa , aaaa , aaaa , aaaa
       if ( cpos .eq. 'Direct' ) then
          V1=X*A(1,1)+Y*A(1,2)+Z*A(1,3)
          V2=X*A(2,1)+Y*A(2,2)+Z*A(2,3)
          V3=X*A(3,1)+Y*A(3,2)+Z*A(3,3)
          X=V1
          Y=V2
          Z=V3
         if ( ionode ) WRITE ( stdout      ,'(A,20A3)' ) 'atomic positions in direct coordinates in POSFF'
       endif

      WRITE( * , * ) atype , x , y , z 
    enddo
  enddo 

  WRITE ( 0 , * ) 'output complete'
 
  CLOSE ( 1000 )

end program trajff_to_xyz
