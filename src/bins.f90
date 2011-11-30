! ===== fmV =====
        PROGRAM disttab

        implicit none

        integer N , PAN
        double precision, dimension (:),allocatable :: a
        double precision, dimension (:),allocatable :: distrib_a
        double precision amin , amax , da , ia
        integer i , ak , ka

        READ ( * , * ) N , PAN , amin , amax

        da = amax - amin
        ia = dble( da / PAN )

        WRITE (6,'(a,i5,a)') 'N = ', N ,' points'
        WRITE (6,'(a,i5,a,f10.6,a,f10.6)') 'bin = ', PAN ,'intervals between', amin ,' and ', amax
        WRITE (6,'(a,f10.6)') 'resolution =', ia

        allocate( a(0:N) , distrib_a(0:PAN+1) )

        distrib_a = 0
        a = 0

        do i = 1, N
          READ ( * , * ) a(i)
        enddo


        do i = 1, N
          ak = ( a(i) - amin ) * ( 1.0d0 / ia )
          ka = int( ak ) + 1
          if( ka .gt. PAN+1 ) then
             WRITE (6, * ) 'ERROR out of bound in distrib_a'
             WRITE (6, * ) i , ka , PAN+1     
             STOP 
          endif
          distrib_a(ka) = distrib_a(ka) + 1
        enddo


        OPEN ( unit=10 , file='distrib.out' )
        do i = 0, PAN
          WRITE (6,'(2f15.8)')  amin + dble( i * ia ) , distrib_a(i) / ( ia * dble(N) )
          WRITE (10,'(2f15.8)') amin + dble( i * ia ) , distrib_a(i) / ( ia * dble(N) )
        enddo
        CLOSE (10)

        WRITE (6, * ) 'file distrib.out'

        deallocate( a )
        deallocate( distrib_a )
        
        STOP 
         
        END PROGRAM disttab
! ===== fmV =====
