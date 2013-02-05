program trajff_to_xyz

  implicit none

  integer*8 :: ia , ic , it 
  integer*8 :: nconf , natm , ntype , iiii
  character(len=35) :: xxxx
  character(len=35) :: filename 
  character(len=2) :: atype
  real*8   :: box , aaaa
  real*8   :: x , y , z  


  read ( * , * ) nconf
  print*,nconf
  read ( * , * ) filename
  print*,filename
  OPEN ( unit = 1000,FILE=filename)
  do ic = 1 , nconf
    read ( 1000 , * ) natm
    write( * , * ) natm
    read ( 1000 , * ) xxxx
    write ( * , * ) xxxx
    read ( 1000 , * ) box , ntype
    READ ( 1000 , * ) ( atype , it = 1 , ntype )
    READ ( 1000 , * ) ( iiii , it = 1 , ntype )
    do ia = 1 , natm      
      read ( 1000 , * ) atype , x , y , z , aaaa , aaaa , aaaa , aaaa , aaaa , aaaa
      write( * , * ) atype , x , y , z 
    enddo
  enddo 
 
  CLOSE ( 1000 )

end program trajff_to_xyz
