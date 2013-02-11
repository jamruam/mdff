! MDFF parallel Molecular Dynamics ... For Fun
! Copyright (C) 2011  F. Vasconcelos
!
! This program is free software; you can redistribute it and/or
! modify it under the terms of the GNU General Public License
! as published by the Free Software Foundation; either version 2
! of the License, or (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program; if not, write to the Free Software
! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
! ===== fmV =====

! ======= Hardware =======
! ======= Hardware =======

MODULE thermodynamic

  USE block

  implicit none

  double precision :: e_tot          ! total energy  ( potential + kinetic )
  double precision :: u_tot          ! total potential energy 
  
  double precision :: e_kin          ! kinetic energy
  double precision :: u_lj           ! potential energy from lennard_jones interaction
  double precision :: u_morse        ! potential energy from morse interaction
  double precision :: u_coul_tot
!  double precision :: u_coul_qq      ! potential energy from coulombic interaction 
!  double precision :: u_coul_dd      ! potential energy from coulombic interaction 
!  double precision :: u_coul_qd      ! potential energy from coulombic interaction
  double precision :: u_pol 

  double precision :: e_kin_r        ! kinetic energy
  double precision :: temp_r         ! temperature ( from kinetic energy )  
  double precision :: u_lj_r         ! potential energy from lennard_jones interaction
  double precision :: u_morse_r      ! potential energy from morse interaction
  double precision :: u_coul_r       ! potential energy from coulombic interaction 


  double precision :: vir_tot        ! total virial
  double precision :: vir_lj         ! virial of lj interaction
  double precision :: vir_morse      ! virial of morse interaction
  double precision :: vir_coul_tot   ! virial of coulombic interaction
!  double precision :: vir_coul_qq    ! virial of coulombic interaction
!  double precision :: vir_coul_dd    ! virial of coulombic interaction
!  double precision :: vir_coul_qd    ! virial of coulombic interaction

  double precision :: pvirial_tot    ! virial correction to the pressure
  double precision :: pvirial_lj     ! virial correction to the pressure
  double precision :: pvirial_coul   ! virial correction to the pressure

  double precision :: pressure_tot   ! total pressure
  double precision :: pressure_lj    ! pressure from lj interactions
  double precision :: pressure_coul  ! pressure from coulombic interactions

  TYPE (accu) acc_e_tot          ! total energy  ( potential + kinetic ) 
  TYPE (accu) acc_u_tot          ! total potential energy 
  
  TYPE (accu) acc_e_kin_r        ! kinetic energy
  TYPE (accu) acc_u_lj_r         ! potential energy from lennard_jones interaction
  TYPE (accu) acc_u_coul_r       ! potential energy from coulombic interaction 
  TYPE (accu) acc_temp_r             ! temperature ( from kinetic energy )  

  TYPE (accu) acc_vir_tot        ! total virial
  TYPE (accu) acc_vir_lj         ! virial of lj interaction
  TYPE (accu) acc_vir_coul       ! virial of coulombic interaction

  TYPE (accu) acc_pressure_tot   ! total pressure
  TYPE (accu) acc_pressure_lj    ! pressure from lj interactions
  TYPE (accu) acc_pressure_coul  ! pressure from coulombic interactions

CONTAINS

!*********************** SUBROUTINE calc_thermo *******************************
!
! calculation of the main thermodynamic quantities 
!
!******************************************************************************
SUBROUTINE calc_thermo

  USE control,  ONLY :  lreduced 
  USE config,   ONLY :  natm , rho , simu_cell 

  implicit none
  double precision :: omega

  omega = simu_cell%omega

  if (lreduced) then
    u_lj_r   = u_lj                      / DBLE ( natm )
    u_morse_r= u_morse                   / DBLE ( natm )
    u_coul_r = u_coul_tot                / DBLE ( natm )
    e_kin_r  = e_kin                     / DBLE ( natm )
  else
    u_lj_r   = u_lj
    u_morse_r= u_morse
    u_coul_r = u_coul_tot
    e_kin_r  = e_kin
  endif
  
  u_tot    = u_lj_r + u_coul_r + u_morse_r
  e_tot    = u_tot  + e_kin_r
  vir_tot  = vir_lj + vir_coul_tot + vir_morse

  if (lreduced) then
    pvirial_lj    = vir_lj   / omega / DBLE ( natm ) 
    pvirial_coul  = ( vir_coul_tot ) / omega / DBLE ( natm ) 
    pvirial_tot   = pvirial_lj + pvirial_coul  

    pressure_tot  = pvirial_tot + temp_r / omega
    pressure_lj   = pvirial_lj   
    pressure_coul = pvirial_coul 
  else
    pvirial_lj    = vir_lj   / omega 
    pvirial_coul  = ( vir_coul_tot )  / omega 
    pvirial_tot   = pvirial_lj + pvirial_coul  

    pressure_tot  = pvirial_tot + rho * temp_r 
    pressure_lj   = pvirial_lj   
    pressure_coul = pvirial_coul 
  endif

  return

END SUBROUTINE calc_thermo


!*********************** SUBROUTINE init_general_accumulator ******************
!
!
!******************************************************************************
SUBROUTINE init_general_accumulator

  implicit none

  acc_e_tot%accval           = 0.0d0 
  acc_e_tot%accvalsq         = 0.0d0
  acc_e_tot%counter          = 0

  acc_u_tot%accval           = 0.0d0
  acc_u_tot%accvalsq         = 0.0d0
  acc_u_tot%counter          = 0

  acc_e_kin_r%accval         = 0.0d0
  acc_e_kin_r%accvalsq       = 0.0d0
  acc_e_kin_r%counter        = 0

  acc_u_lj_r%accval          = 0.0d0
  acc_u_lj_r%accvalsq        = 0.0d0
  acc_u_lj_r%counter         = 0
 
  acc_u_coul_r%accval        = 0.0d0
  acc_u_coul_r%accvalsq      = 0.0d0
  acc_u_coul_r%counter       = 0

  acc_temp_r%accval          = 0.0d0
  acc_temp_r%accvalsq        = 0.0d0
  acc_temp_r%counter         = 0

  acc_vir_tot%accval         = 0.0d0
  acc_vir_tot%accvalsq       = 0.0d0
  acc_vir_tot%counter        = 0

  acc_vir_lj%accval          = 0.0d0
  acc_vir_lj%accvalsq        = 0.0d0
  acc_vir_lj%counter         = 0

  acc_vir_coul%accval        = 0.0d0
  acc_vir_coul%accvalsq      = 0.0d0
  acc_vir_coul%counter       = 0

  acc_pressure_tot%accval    = 0.0d0
  acc_pressure_tot%accvalsq  = 0.0d0
  acc_pressure_tot%counter   = 0

  acc_pressure_lj%accval     = 0.0d0
  acc_pressure_lj%accvalsq   = 0.0d0
  acc_pressure_lj%counter    = 0

  acc_pressure_coul%accval   = 0.0d0
  acc_pressure_coul%accvalsq = 0.0d0
  acc_pressure_coul%counter  = 0

  return

END SUBROUTINE init_general_accumulator

!*********************** SUBROUTINE general_accumulator ***********************
!
!
!******************************************************************************
SUBROUTINE general_accumulator

  implicit none

  acc_e_tot%accval           = acc_e_tot%accval           + e_tot
  acc_e_tot%accvalsq         = acc_e_tot%accvalsq         + e_tot * e_tot
  acc_e_tot%counter          = acc_e_tot%counter          + 1

  acc_u_tot%accval           = acc_u_tot%accval           + u_tot
  acc_u_tot%accvalsq         = acc_u_tot%accvalsq         + u_tot * u_tot
  acc_u_tot%counter          = acc_u_tot%counter          + 1

  acc_e_kin_r%accval         = acc_e_kin_r%accval         + e_kin_r
  acc_e_kin_r%accvalsq       = acc_e_kin_r%accvalsq       + e_kin_r * e_kin_r
  acc_e_kin_r%counter        = acc_e_kin_r%counter        + 1

  acc_u_lj_r%accval          = acc_u_lj_r%accval          + u_lj_r
  acc_u_lj_r%accvalsq        = acc_u_lj_r%accvalsq        + u_lj_r * u_lj_r
  acc_u_lj_r%counter         = acc_u_lj_r%counter         + 1
 
  acc_u_coul_r%accval        = acc_u_coul_r%accval        + u_coul_r
  acc_u_coul_r%accvalsq      = acc_u_coul_r%accvalsq      + u_coul_r * u_coul_r
  acc_u_coul_r%counter       = acc_u_coul_r%counter       + 1

  acc_temp_r%accval          = acc_temp_r%accval          + temp_r
  acc_temp_r%accvalsq        = acc_temp_r%accvalsq        + temp_r * temp_r
  acc_temp_r%counter         = acc_temp_r%counter         + 1

  acc_vir_tot%accval         = acc_vir_tot%accval         + vir_tot
  acc_vir_tot%accvalsq       = acc_vir_tot%accvalsq       + vir_tot * vir_tot 
  acc_vir_tot%counter        = acc_vir_tot%counter        + 1

  acc_vir_lj%accval          = acc_vir_lj%accval          + vir_lj
  acc_vir_lj%accvalsq        = acc_vir_lj%accvalsq        + vir_lj * vir_lj
  acc_vir_lj%counter         = acc_vir_lj%counter         + 1

  acc_vir_coul%accval        = acc_vir_coul%accval        + vir_coul_tot
  acc_vir_coul%accvalsq      = acc_vir_coul%accvalsq      + vir_coul_tot * vir_coul_tot
  acc_vir_coul%counter       = acc_vir_coul%counter       + 1

  acc_pressure_tot%accval    = acc_pressure_tot%accval    + pressure_tot
  acc_pressure_tot%accvalsq  = acc_pressure_tot%accvalsq  + pressure_tot * pressure_tot
  acc_pressure_tot%counter   = acc_pressure_tot%counter   + 1

  acc_pressure_lj%accval     = acc_pressure_lj%accval     + pressure_lj
  acc_pressure_lj%accvalsq   = acc_pressure_lj%accvalsq   + pressure_lj * pressure_lj
  acc_pressure_lj%counter    = acc_pressure_lj%counter    + 1

  acc_pressure_coul%accval   = acc_pressure_coul%accval   + pressure_coul
  acc_pressure_coul%accvalsq = acc_pressure_coul%accvalsq + pressure_coul * pressure_coul
  acc_pressure_coul%counter  = acc_pressure_coul%counter  + 1

  return

END SUBROUTINE general_accumulator

!*********************** SUBROUTINE write_thermo **********************
!
! SUBROUTINE write_thermo
! write thermodynamic quantities to file OSZIFF or print to standard output 
!
!******************************************************************************

SUBROUTINE write_thermo ( step , kunit , dummy )

  USE config,   ONLY :  simu_cell 
  USE md,       ONLY :  dt
  USE io_file,  ONLY :  ionode

  implicit none

  ! global
  integer, intent(in) :: kunit , step
  double precision, intent(in), optional :: dummy

  ! local 
  double precision :: omega
  
  omega = simu_cell%omega

  call calc_thermo

  if ( ionode ) then
    if ( PRESENT ( dummy ) ) then
      WRITE ( kunit , 200 ) &
      step , DBLE (step * dt) , e_tot   , e_kin_r      , u_tot        , u_lj_r      , u_coul_r  , u_morse_r    , dummy
      WRITE ( kunit , 201 ) &
      step , DBLE (step * dt) , temp_r  , pressure_tot , pressure_lj , pressure_coul , omega , dummy  
    else
      WRITE ( kunit , 100 ) &
      step , DBLE (step * dt) , e_tot   , e_kin_r      , u_tot        , u_lj_r      , u_coul_r , u_morse_r       
      WRITE ( kunit , 101 ) &
      step , DBLE (step * dt) , temp_r  , pressure_tot , pressure_lj , pressure_coul, omega 
    endif
  endif

 100 FORMAT(I9,2X,E15.8,'  Etot = ',E15.8,'  Ekin  = ',E15.8,'  Utot  = ',&
                E15.8,'  U_lj   = ',E15.8,'  U_coul   = ',E15.8,'  U_morse   = ',E15.8)
 101 FORMAT(I9,2X,E15.8,'  Temp = ',E15.8,'  Press = ',E15.8,'  P_lj  = ',&
                E15.8,'  P_coul = ',E15.8,'  Volume   = ',E15.8)
 200 FORMAT(I9,2X,E15.8,'  Etot = ',E15.8,'  Ekin  = ',E15.8,'  Utot  = ',&
                E15.8,'  U_lj   = ',E15.8,'  U_coul   = ',E15.8,'  U_morse   = ',E15.8,'  dummy = ',E15.8)
 201 FORMAT(I9,2X,E15.8,'  Temp = ',E15.8,'  Press = ',E15.8,'  P_lj  = ',&
                E15.8,'  P_coul = ',E15.8,'  Volume   = ',E15.8,'  dummy = ',E15.8)

  return

END SUBROUTINE write_thermo

!*********************** SUBROUTINE write_average_thermo **********************
!
! SUBROUTINE write_average_thermo
! write time average thermodynamic quantities to file OSZIFF or print to standard output 
!
!******************************************************************************

SUBROUTINE write_average_thermo ( kunit )

  USE md,       ONLY :  dt
  USE io_file,  ONLY :  ionode

  implicit none

  ! global
  integer, intent(in) :: kunit

  !local
  double precision :: e_tot_av , e_kin_r_av , u_tot_av , u_lj_r_av
  double precision :: u_coul_r_av , temp_r_av  , pressure_tot_av , pressure_lj_av  , pressure_coul_av               

  double precision :: e_tot_avsq , e_kin_r_avsq , u_tot_avsq , u_lj_r_avsq
  double precision :: u_coul_r_avsq , temp_r_avsq  , pressure_tot_avsq , pressure_lj_avsq  , pressure_coul_avsq

  double precision :: e_tot_sig , e_kin_r_sig , u_tot_sig , u_lj_r_sig
  double precision :: u_coul_r_sig , temp_r_sig  , pressure_tot_sig , pressure_lj_sig  , pressure_coul_sig
 
  if ( acc_e_tot%counter         .eq. 0 ) return
  if ( acc_e_kin_r%counter       .eq. 0 ) return
  if ( acc_u_tot%counter         .eq. 0 ) return 
  if ( acc_u_lj_r%counter        .eq. 0 ) return
  if ( acc_u_coul_r%counter      .eq. 0 ) return
  if ( acc_temp_r%counter        .eq. 0 ) return
  if ( acc_pressure_tot%counter  .eq. 0 ) return
  if ( acc_pressure_lj%counter   .eq. 0 ) return
  if ( acc_pressure_coul%counter .eq. 0 ) return
 
  ! ========
  !  < A >
  ! ========
  e_tot_av           = acc_e_tot%accval           / DBLE ( acc_e_tot%counter )
  e_kin_r_av         = acc_e_kin_r%accval         / DBLE ( acc_e_kin_r%counter ) 
  u_tot_av           = acc_u_tot%accval           / DBLE ( acc_u_tot%counter )
  u_lj_r_av          = acc_u_lj_r%accval          / DBLE ( acc_u_lj_r%counter )
  u_coul_r_av        = acc_u_coul_r%accval        / DBLE ( acc_u_coul_r%counter )
  temp_r_av          = acc_temp_r%accval          / DBLE ( acc_temp_r%counter )
  pressure_tot_av    = acc_pressure_tot%accval    / DBLE ( acc_pressure_tot%counter )
  pressure_lj_av     = acc_pressure_lj%accval     / DBLE ( acc_pressure_lj%counter )
  pressure_coul_av   = acc_pressure_coul%accval   / DBLE ( acc_pressure_coul%counter )

  ! ========
  !  < A² >
  ! ========
  e_tot_avsq         = acc_e_tot%accvalsq         / DBLE ( acc_e_tot%counter )
  e_kin_r_avsq       = acc_e_kin_r%accvalsq       / DBLE ( acc_e_kin_r%counter ) 
  u_tot_avsq         = acc_u_tot%accvalsq         / DBLE ( acc_u_tot%counter )
  u_lj_r_avsq        = acc_u_lj_r%accvalsq        / DBLE ( acc_u_lj_r%counter )
  u_coul_r_avsq      = acc_u_coul_r%accvalsq      / DBLE ( acc_u_coul_r%counter )
  temp_r_avsq        = acc_temp_r%accvalsq        / DBLE ( acc_temp_r%counter )
  pressure_tot_avsq  = acc_pressure_tot%accvalsq  / DBLE ( acc_pressure_tot%counter )
  pressure_lj_avsq   = acc_pressure_lj%accvalsq   / DBLE ( acc_pressure_lj%counter )
  pressure_coul_avsq = acc_pressure_coul%accvalsq / DBLE ( acc_pressure_coul%counter )

  ! =============================
  !  sqrt ( < A² > - < A > ²)
  ! =============================
  e_tot_sig         = SQRT ( e_tot_avsq         - e_tot_av         * e_tot_av         )
  e_kin_r_sig       = SQRT ( e_kin_r_avsq       - e_kin_r_av       * e_kin_r_av       )
  u_tot_sig         = SQRT ( u_tot_avsq         - u_tot_av         * u_tot_av         )
  u_lj_r_sig        = SQRT ( u_lj_r_avsq        - u_lj_r_av        * u_lj_r_av        )
  u_coul_r_sig      = SQRT ( u_coul_r_avsq      - u_coul_r_av      * u_coul_r_av      )
  temp_r_sig        = SQRT ( temp_r_avsq        - temp_r_av        * temp_r_av        )
  pressure_tot_sig  = SQRT ( pressure_tot_avsq  - pressure_tot_av  * pressure_tot_av  )
  pressure_lj_sig   = SQRT ( pressure_lj_avsq   - pressure_lj_av   * pressure_lj_av   )
  pressure_coul_sig = SQRT ( pressure_coul_avsq - pressure_coul_av * pressure_coul_av )


  if ( ionode ) then
  !   WRITE ( kunit , '(30i6)' )  acc_e_tot%counter , acc_e_kin_r%counter, acc_u_tot%counter , acc_u_lj_r%counter , acc_u_coul_r%counter , acc_temp_r%counter , acc_pressure_tot%counter , acc_pressure_lj%counter , acc_pressure_coul%counter 

      WRITE ( kunit , 100 ) e_tot_av    , e_kin_r_av       , u_tot_av         , u_lj_r_av         , u_coul_r_av        
      WRITE ( kunit , 102 ) e_tot_sig   , e_kin_r_sig      , u_tot_sig        , u_lj_r_sig        , u_coul_r_sig
      WRITE ( kunit , 101 ) temp_r_av   , pressure_tot_av  , pressure_lj_av   , pressure_coul_av  
      WRITE ( kunit , 103 ) temp_r_sig  , pressure_tot_sig , pressure_lj_sig  , pressure_coul_sig  
  endif

 100 FORMAT(2X,'Aver. values:  <Etot>     = ',E15.8,'  <Ekin>      = ',&
                     E15.8,'  <Utot>      = ',E15.8,'  <U_lj>      = ',E15.8,'  <U_coul>       = ',E15.8)
 101 FORMAT(2X,'Aver. values:  <Temp>     = ',E15.8,'  <Press>     = ',&
                     E15.8,'  <P_lj>      = ',E15.8,'  <P_coul>    = ',E15.8)
 102 FORMAT(2X,'Var.  estim :  sig2(Etot) = ',E15.8,'  sig2(Ekin)  = ',&
                     E15.8,'  sig2(Utot)  = ',E15.8,'  sig2(U_lj)  = ',E15.8,'  sig2(U_coul)   = ',E15.8)
 103 FORMAT(2X,'Var.  estim :  sig2(Temp) = ',E15.8,'  sig2(Press) = ',&
                     E15.8,'  sig2(P_lj)  = ',E15.8,'  sig2(P_coul)= ',E15.8)

  return

END SUBROUTINE write_average_thermo

END MODULE thermodynamic
! ===== fmV =====
