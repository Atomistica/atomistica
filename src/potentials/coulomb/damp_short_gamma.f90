!! ==========================================================
!!
!! damping function for short-range gamma for DFTB3 Hamiltonian
!! 
!! Gaus et al. J. Chem. Theory Comput 7, 931 (2001).
!!
!! ==========================================================
module damp_short_gamma
  use supplib

  implicit none

  private
  
  public :: hij                                   !> hij         = h(rij,Ui,Uj)  
  public :: part_deriv_hij_wrt_Ui                 !> dhij/dUi    = dh(rij,Ui,Uj)/dUi   
  public :: part_deriv_hij_wrt_r                  !> dhij/dr     = dh(rij,Ui,Uj)/dr 
  public :: second_part_deriv_hij_wrt_Ui_and_r    !> d2hij/dUidr = d2h(rij,Ui,Uj)/dUidr

contains
 
  function hij(abs_rij, U_i, U_j, zeta) result(res)
  real(DP), intent(in) :: abs_rij
  real(DP), intent(in) :: U_i, U_j
  real(DP), intent(in) :: zeta
  real(DP)             :: res
  real(DP)             :: efact
  
  efact = - (0.50_DP*(U_i + U_j))**zeta
  res = exp(efact*abs_rij**2)

  endfunction hij

  function part_deriv_hij_wrt_Ui(abs_rij, U_i, U_j, zeta) result(res)
  real(DP), intent(in) :: abs_rij
  real(DP), intent(in) :: U_i, U_j
  real(DP), intent(in) :: zeta
  real(DP)             :: fact
  real(DP)             :: res

  fact = - (0.50_DP*(U_i + U_j))**(zeta - 1.0_DP)
  res = - 0.50_DP*zeta*abs_rij**2*fact*hij(abs_rij, U_i, U_j, zeta)

  endfunction part_deriv_hij_wrt_Ui

  function part_deriv_hij_wrt_r(abs_rij, U_i, U_j, zeta) result(res)
  real(DP), intent(in) :: abs_rij
  real(DP), intent(in) :: U_i, U_j
  real(DP), intent(in) :: zeta
  real(DP)             :: fact
  real(DP)             :: res

  fact = - (0.50_DP*(U_i + U_j))**zeta
  res = - 2.0_DP*abs_rij*fact*hij(abs_rij, U_i, U_j, zeta)

  endfunction part_deriv_hij_wrt_r

  function second_part_deriv_hij_wrt_Ui_and_r(abs_rij, U_i, U_j, zeta) result(res)
  real(DP), intent(in) :: abs_rij
  real(DP), intent(in) :: U_i, U_j
  real(DP), intent(in) :: zeta
  real(DP)             :: fact1, fact2
  real(DP)             :: res

  fact1 = - (0.50_DP*(U_i + U_j))**(zeta - 1.0_DP)
  fact2 = - (0.50_DP*(U_i + U_j))**zeta

  res = zeta*abs_rij*fact1*(abs_rij**2*fact2 - 1.0_DP)*hij(abs_rij, U_i, U_j, zeta)

  endfunction second_part_deriv_hij_wrt_Ui_and_r

endmodule damp_short_gamma
