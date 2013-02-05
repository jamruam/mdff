SUBROUTINE structure_factor_evaluate ( delta, npts, ex, ey, ez )

    INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND ( 14, 200 )
    REAL(KIND=dp), DIMENSION(:), INTENT(in)  :: delta
    INTEGER, DIMENSION(:), INTENT(IN)        :: npts
    COMPLEX(KIND=dp), DIMENSION(:), &
      INTENT(out)                            :: ex
    COMPLEX(KIND=dp), DIMENSION(:), &
      INTENT(out)                            :: ey
    COMPLEX(KIND=dp), DIMENSION(:), &
      INTENT(out)                            :: ez

    COMPLEX(KIND=dp)                         :: fm, fp
    INTEGER                                  :: j, l0, l1, m0, m1, n0, n1
    REAL(KIND=dp)                            :: vec( 3 )

    l0 = LBOUND ( ex, 1 )
    l1 = UBOUND ( ex, 1 )
    m0 = LBOUND ( ey, 1 )
    m1 = UBOUND ( ey, 1 )
    n0 = LBOUND ( ez, 1 )
    n1 = UBOUND ( ez, 1 )

    ! delta is in scaled coordinates
    vec ( : ) = twopi * ( delta ( : ) + 0.5_dp  )

    ex ( l0 ) = 1.0_dp
    ey ( m0 ) = 1.0_dp
    ez ( n0 ) = 1.0_dp
    ex ( l1 ) = 1.0_dp
    ey ( m1 ) = 1.0_dp
    ez ( n1 ) = 1.0_dp

    fp = CMPLX ( COS ( vec ( 1 ) ), -SIN ( vec ( 1 ) ),KIND=dp)
    fm = CONJG ( fp )
    DO j = 1, -l0
       ex (  j + l0 ) = ex (  j + l0 - 1 ) * fp
       ex ( -j + l1 ) = ex ( -j + l1 + 1 ) * fm
    END DO

    fp = CMPLX ( COS ( vec ( 2 ) ), -SIN ( vec ( 2 ) ),KIND=dp)
    fm = CONJG ( fp )
    DO j = 1, -m0
       ey (  j + m0 ) = ey (  j + m0 - 1 ) * fp
       ey ( -j + m1 ) = ey ( -j + m1 + 1 ) * fm
    END DO

    fp = CMPLX ( COS ( vec ( 3 ) ), -SIN ( vec ( 3 ) ),KIND=dp)
    fm = CONJG ( fp )
    DO j = 1, -n0
       ez (  j + n0 ) = ez (  j + n0 - 1 ) * fp
       ez ( -j + n1 ) = ez ( -j + n1 + 1 ) * fm
    END DO

  END SUBROUTINE structure_factor_evaluate

