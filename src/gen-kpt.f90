! ===== fmV =====
        PROGRAM gene_kpt

        implicit none

        integer :: ikx , iky , ikz
        integer :: nk , wi , nktot
        double precision :: ak , pi

        wi = 1
        pi = 3.14159265358979323846d0

        READ (*,*) nk
        nktot = ( ( nk + 1 ) *( nk + 1 ) * ( nk + 1 ) )

        WRITE (*,*) 'Automatically generated mesh FF'
        WRITE (*,*) nktot
        WRITE (*,*) 'Reciprocal lattice'


        ak=(1.0d0)/(nk) 


        do ikx = 0 , nk
          do iky = 0 , nk        
            do ikz = 0 , nk
              WRITE (*,'(3f16.12,i6)') ikx * ak , iky * ak , ikz * ak, wi
            enddo
          enddo
        enddo

        END PROGRAM gene_kpt

! ===== fmV =====
