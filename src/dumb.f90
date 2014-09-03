module dumb

  USE constants,        ONLY : dp
  real(kind=dp) ::  tdt21, tdt32,tdt54,tdt65,tdt76,tdt83,tdt98

CONTAINS

SUBROUTINE print_dumb

  implicit none

  write(*,'(a,e16.8)') "tdt21 = ",tdt21
  write(*,'(a,e16.8)') "tdt32 = ",tdt32
  write(*,'(a,e16.8)') "tdt54 = ",tdt54
  write(*,'(a,e16.8)') "tdt65 = ",tdt65
  write(*,'(a,e16.8)') "tdt76 = ",tdt76
  write(*,'(a,e16.8)') "tdt83 = ",tdt83
  write(*,'(a,e16.8)') "tdt98 = ",tdt98
  write(*,'(a,2e16.8)') "sum loop = ",tdt54+tdt65+tdt76,tdt83
  write(*,'(a,e16.8)') "all",tdt21+tdt32+tdt83+tdt98 
  return
end subroutine

end module dumb
