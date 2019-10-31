!! ======================================================================
!! Atomistica - Interatomic potential library and molecular dynamics code
!! https://github.com/Atomistica/atomistica
!!
!! Copyright (2005-2020) Lars Pastewka <lars.pastewka@imtek.uni-freiburg.de>
!! and others. See the AUTHORS file in the top-level Atomistica directory.
!!
!! This program is free software: you can redistribute it and/or modify
!! it under the terms of the GNU General Public License as published by
!! the Free Software Foundation, either version 2 of the License, or
!! (at your option) any later version.
!!
!! This program is distributed in the hope that it will be useful,
!! but WITHOUT ANY WARRANTY; without even the implied warranty of
!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!! GNU General Public License for more details.
!!
!! You should have received a copy of the GNU General Public License
!! along with this program.  If not, see <http://www.gnu.org/licenses/>.
!! ======================================================================
!! ============================================================
!!
!! short-range gamma function for DFTB3 Hamiltonian and forces
!! 
!! Gaus et al. J. Chem. Theory Comput 7, 931 (2001).
!!
module coulomb_short_gamma
  use supplib

  use damp_short_gamma

  implicit none

  private
  
  public :: capital_short_gamma
  public :: derivative_capital_short_gamma
  public :: Sfij, Sgij 

contains
 
  function capital_short_gamma(abs_rij, dU_i, U_i, U_j, zeta) result(res)
  real(DP), intent(in)           :: abs_rij
  real(DP), intent(in)           :: dU_i
  real(DP), intent(in)           :: U_i, U_j
  real(DP), intent(in), optional :: zeta
  real(DP)                       :: res
  
  if (present(zeta)) then

    res = part_deriv_sgamma_wrt_Ui(abs_rij, U_i, U_j, zeta)*dU_i

  else

    res = part_deriv_sgamma_wrt_Ui(abs_rij, U_i, U_j)*dU_i

  endif

  endfunction capital_short_gamma  

  function derivative_capital_short_gamma(abs_rij, dU_i, U_i, U_j, zeta) result(res)
  real(DP), intent(in)           :: abs_rij
  real(DP), intent(in)           :: dU_i
  real(DP), intent(in)           :: U_i, U_j
  real(DP), intent(in), optional :: zeta
  real(DP)                       :: res
  
  if (present(zeta)) then

     res = second_part_deriv_csgamma_wrt_Ui_and_r(abs_rij, U_i, U_j, zeta)*dU_i

  else

     res = second_part_deriv_csgamma_wrt_Ui_and_r(abs_rij, U_i, U_j)*dU_i

  endif

  endfunction derivative_capital_short_gamma  

  function short_gamma(abs_rij, U_i, U_j, zeta) result(res)
  real(DP), intent(in)           :: abs_rij
  real(DP), intent(in)           :: U_i, U_j
  real(DP), intent(in), optional :: zeta
  real(DP)                       :: res

  if (present(zeta)) then
    
     if (abs(U_i - U_j) < 1.0d-06) then
    
        res = - Sgij(abs_rij, U_i)*hij(abs_rij, U_i, U_j, zeta) 
   
     else
     
        res = - Sfij(abs_rij, U_i, U_j)*hij(abs_rij, U_i, U_j, zeta) 
     
     endif

  else 

     if (abs(U_i - U_j) < 1.0d-06) then
    
        res = - Sgij(abs_rij, U_i)
     
     else
     
        res = - Sfij(abs_rij, U_i, U_j)
     
     endif
     
  endif

  endfunction short_gamma   

  function second_part_deriv_csgamma_wrt_Ui_and_r(abs_rij, U_i, U_j, zeta) result(res)
  real(DP), intent(in)           :: abs_rij
  real(DP), intent(in)           :: U_i, U_j
  real(DP), intent(in), optional :: zeta
  real(DP)                       :: res
 
  if (present(zeta)) then
    
     if (abs(U_i - U_j) < 1.0d-06) then
     
        res = - (second_part_deriv_Sgij_wrt_Ui_and_r(abs_rij, U_i)*hij(abs_rij, U_i, U_j, zeta) &
               + part_deriv_Sgij_wrt_Ui(abs_rij, U_i)*part_deriv_hij_wrt_r(abs_rij, U_i, U_j, zeta) &
               + part_deriv_Sgij_wrt_r(abs_rij, U_i)*part_deriv_hij_wrt_Ui(abs_rij, U_i, U_j, zeta) &
               + Sgij(abs_rij, U_i)*second_part_deriv_hij_wrt_Ui_and_r(abs_rij, U_i, U_j, zeta))
     
     else
     
        res = - (second_part_deriv_Sfij_wrt_Ui_and_r(abs_rij, U_i, U_j)*hij(abs_rij, U_i, U_j, zeta) &
               + part_deriv_Sfij_wrt_Ui(abs_rij, U_i, U_j)*part_deriv_hij_wrt_r(abs_rij, U_i, U_j, zeta) &
               + part_deriv_Sfij_wrt_r(abs_rij, U_i, U_j)*part_deriv_hij_wrt_Ui(abs_rij, U_i, U_j, zeta) &
               + Sfij(abs_rij, U_i, U_j)*second_part_deriv_hij_wrt_Ui_and_r(abs_rij, U_i, U_j, zeta))
     
     endif
  
  else 

     if (abs(U_i - U_j) < 1.0d-06) then
    
        res =  - second_part_deriv_Sgij_wrt_Ui_and_r(abs_rij, U_i)
   
     else
    
        res =  - second_part_deriv_Sfij_wrt_Ui_and_r(abs_rij, U_i, U_j)
    
     endif

  endif

  endfunction second_part_deriv_csgamma_wrt_Ui_and_r

  function part_deriv_sgamma_wrt_Ui(abs_rij, U_i, U_j, zeta) result(res)
  real(DP), intent(in)           :: abs_rij
  real(DP), intent(in)           :: U_i, U_j
  real(DP), intent(in), optional :: zeta
  real(DP)                       :: res

  if (present(zeta)) then
    
     if (abs(U_i - U_j) < 1.0d-06) then
    
        res = - 3.20_DP*part_deriv_Sgij_wrt_Ui(abs_rij, U_i)*hij(abs_rij, U_i, U_j, zeta) &
              - Sgij(abs_rij, U_i)*part_deriv_hij_wrt_Ui(abs_rij, U_i, U_j, zeta)
     
     else
    
        res = - 3.20_DP*part_deriv_Sfij_wrt_Ui(abs_rij, U_i, U_j)*hij(abs_rij, U_i, U_j, zeta) &
              - Sfij(abs_rij, U_i, U_j)*part_deriv_hij_wrt_Ui(abs_rij, U_i, U_j, zeta)
    
     endif
  
  else 

     if (abs(U_i - U_j) < 1.0d-06) then
    
        res =  - 3.20_DP*part_deriv_Sgij_wrt_Ui(abs_rij, U_i)
   
     else
 
        res =  - 3.20_DP*part_deriv_Sfij_wrt_Ui(abs_rij, U_i, U_j)

     endif

  endif

  endfunction part_deriv_sgamma_wrt_Ui

  function part_deriv_sgamma_wrt_r(abs_rij, U_i, U_j, zeta) result(res)
  real(DP), intent(in)           :: abs_rij
  real(DP), intent(in)           :: U_i, U_j
  real(DP), intent(in), optional :: zeta
  real(DP)                       :: res

  if (present(zeta)) then
    
     if (abs(U_i - U_j) < 1.0d-06) then
   
        res = - part_deriv_Sgij_wrt_r(abs_rij, U_i)*hij(abs_rij, U_i, U_j, zeta) &
              - Sgij(abs_rij, U_i)*part_deriv_hij_wrt_r(abs_rij, U_i, U_j, zeta)
    
     else 
    
        res = - part_deriv_Sfij_wrt_r(abs_rij, U_i, U_j)*hij(abs_rij, U_i, U_j, zeta) &
              - Sfij(abs_rij, U_i, U_j)*part_deriv_hij_wrt_r(abs_rij, U_i, U_j, zeta)
    
     endif

  else

     if (abs(U_i - U_j) < 1.0d-06) then
     
        res = - part_deriv_Sgij_wrt_r(abs_rij, U_i)
   
     else
    
        res = - part_deriv_Sfij_wrt_r(abs_rij, U_i, U_j)
    
     endif

  endif

  endfunction part_deriv_sgamma_wrt_r

  function second_part_deriv_Sfij_wrt_Ui_and_r(abs_rij, U_i, U_j) result(res)
  real(DP), intent(in) :: abs_rij
  real(DP), intent(in) :: U_i, U_j
  real(DP)             :: expi, expj
  real(DP)             :: res

  expi =exp(-U_i*abs_rij)
  expj =exp(-U_j*abs_rij)

  res = 3.20_DP*(expi*((U_i*abs_rij - 1.0_DP)*fij(abs_rij, U_i, U_j) &
                    - U_i*part_deriv_fij_wrt_Ui(abs_rij, U_i, U_j) &
                    + second_part_deriv_fij_wrt_Ui_and_r(abs_rij, U_i, U_j) &
                    - abs_rij*part_deriv_fij_wrt_r(abs_rij, U_i, U_j)) &
              + expj*(second_part_deriv_fij_wrt_Uj_and_r(abs_rij, U_j, U_i) &
                    - U_j*part_deriv_fij_wrt_Uj(abs_rij, U_j, U_i)))
 
  endfunction second_part_deriv_Sfij_wrt_Ui_and_r

  function second_part_deriv_Sgij_wrt_Ui_and_r(abs_rij, U_i) result(res)
  real(DP), intent(in) :: abs_rij
  real(DP), intent(in) :: U_i
  real(DP)             :: expi
  real(DP)             :: res

  expi =exp(-U_i*abs_rij)

  res = 3.2_DP*expi*((U_i*abs_rij - 1.0_DP)*gij(abs_rij, U_i) &
                    - U_i*part_deriv_gij_wrt_Ui(abs_rij, U_i) &
                    + second_part_deriv_gij_wrt_Ui_and_r(abs_rij, U_i) &
                    - abs_rij*part_deriv_gij_wrt_r(abs_rij, U_i))

  endfunction second_part_deriv_Sgij_wrt_Ui_and_r

  function part_deriv_Sfij_wrt_Ui(abs_rij, U_i, U_j) result(res)
  real(DP), intent(in) :: abs_rij
  real(DP), intent(in) :: U_i, U_j
  real(DP)             :: expi, expj
  real(DP)             :: res
  
  expi = exp(-U_i*abs_rij)
  expj = exp(-U_j*abs_rij)
  
  res = expi*part_deriv_fij_wrt_Ui(abs_rij, U_i, U_j) &
      - abs_rij*expi*fij(abs_rij, U_i, U_j) &
      + expj*part_deriv_fij_wrt_Uj(abs_rij, U_j, U_i)

  endfunction part_deriv_Sfij_wrt_Ui

  function part_deriv_Sgij_wrt_Ui(abs_rij, U_i) result(res)
  real(DP), intent(in) :: abs_rij
  real(DP), intent(in) :: U_i
  real(DP)             :: expi
  real(DP)             :: res

  expi = exp(-U_i*abs_rij)
   
  res = expi*part_deriv_gij_wrt_Ui(abs_rij, U_i) - abs_rij*expi*gij(abs_rij, U_i)

  endfunction part_deriv_Sgij_wrt_Ui

  function part_deriv_Sfij_wrt_r(abs_rij, U_i, U_j) result(res)
  real(DP), intent(in) :: abs_rij
  real(DP), intent(in) :: U_i, U_j
  real(DP)             :: expi, expj
  real(DP)             :: res
  
  expi = exp(-U_i*abs_rij)
  expj = exp(-U_j*abs_rij)
  
  res = expi*(part_deriv_fij_wrt_r(abs_rij, U_i, U_j) - U_i*fij(abs_rij, U_i, U_j)) &
      + expj*(part_deriv_fij_wrt_r(abs_rij, U_j, U_i) - U_j*fij(abs_rij, U_j, U_i))

  endfunction part_deriv_Sfij_wrt_r

  function part_deriv_Sgij_wrt_r(abs_rij, U_i) result(res)
  real(DP), intent(in) :: abs_rij
  real(DP), intent(in) :: U_i
  real(DP)             :: expi
  real(DP)             :: res
  
  expi = exp(-U_i*abs_rij)
  res = expi*(part_deriv_gij_wrt_r(abs_rij, U_i) - U_i*gij(abs_rij, U_i))

  endfunction part_deriv_Sgij_wrt_r

  function part_deriv_fij_wrt_r(abs_rij, U_i, U_j) result(res)
  real(DP), intent(in) :: abs_rij
  real(DP), intent(in) :: U_i, U_j
  real(DP)             :: res
  real(DP)             :: term1, term2
  
  term1 = U_j**6 - 3.0_DP*U_i**2*U_j**4
  term2 = (U_i**2 - U_j**2)**3*abs_rij**2

  res = term1/term2

  endfunction  part_deriv_fij_wrt_r

  function part_deriv_gij_wrt_r(abs_rij, U_i) result(res)
  real(DP), intent(in) :: abs_rij
  real(DP), intent(in) :: U_i
  real(DP)             :: res

  res = - 1.0_DP/abs_rij**2 + 3.0_DP/16.0_DP*U_i**2 + 1.0_DP/24.0_DP*U_i**3*abs_rij

  endfunction part_deriv_gij_wrt_r

  function Sfij(abs_rij, U_i, U_j) result(res)
  real(DP), intent(in) :: abs_rij
  real(DP), intent(in) :: U_i, U_j
  real(DP)             :: expi, expj
  real(DP)             :: res
  
  expi = exp(-U_i*abs_rij)
  expj = exp(-U_j*abs_rij)
  
  res = expi*fij(abs_rij, U_i, U_j) + expj*fij(abs_rij, U_j, U_i)

  return

  endfunction Sfij

  function Sgij(abs_rij, U_i) result(res)
  real(DP), intent(in) :: abs_rij
  real(DP), intent(in) :: U_i
  real(DP)             :: expi
  real(DP)             :: res
  
  expi = exp(-U_i*abs_rij)
  res = expi*gij(abs_rij,U_i)  

  endfunction Sgij

  function fij(abs_rij, U_i, U_j) result(res)
  real(DP), intent(in) :: abs_rij
  real(DP), intent(in) :: U_i, U_j
  real(DP)             :: res

  real(DP)             :: term1, term2, term3, term4
  
  term1 = U_i*U_j**4
  term2 = 2.0_DP*(U_i**2 - U_j**2)**2
  
  term3 = U_j**6 - 3.0_DP*U_i**2*U_j**4
  term4 = (U_i**2 - U_j**2)**3*abs_rij

  res = term1/term2 - term3/term4 

  endfunction fij

  function gij(abs_rij, U_i) result(res)
  real(DP), intent(in) :: abs_rij
  real(DP), intent(in) :: U_i
  real(DP)             :: res
  real(DP)             :: uir
   
  uir = U_i*abs_rij
  res = 48.0_DP + (33.0_DP + (9.0_DP + uir)*uir)*uir
  res = res/(48.0_DP*abs_rij)

  endfunction gij

  function part_deriv_fij_wrt_Ui(abs_rij, U_i, U_j) result(res)
  real(DP), intent(in) :: abs_rij
  real(DP), intent(in) :: U_i, U_j
  real(DP)             :: res
  real(DP)             :: term1, term2, term3, term4  

  term1 = U_j**6 + 3.0_DP*U_i**2*U_j**4
  term2 = 2.0_DP*(U_i**2 - U_j**2)**3

  term3 = 12.0_DP*U_i**3*U_j**4
  term4 = (U_i**2 - U_j**2)**4*abs_rij

  res = - term1/term2 - term3/term4

  return

  endfunction part_deriv_fij_wrt_Ui

  function part_deriv_fij_wrt_Uj(abs_rij, U_i, U_j) result(res)
  real(DP), intent(in) :: abs_rij
  real(DP), intent(in) :: U_i, U_j
  real(DP)             :: res
  real(DP)             :: term1, term2, term3, term4

  term1 = 2.0_DP*U_i**3*U_j**3
  term2 = (U_i**2 - U_j**2)**3

  term3 = 12.0_DP*U_i**4*U_j**3
  term4 = (U_i**2 - U_j**2)**4*abs_rij

  res = term1/term2 + term3/term4

  return

  endfunction part_deriv_fij_wrt_Uj

  function part_deriv_gij_wrt_Ui(abs_rij, U_i) result(res)
  real(DP), intent(in) :: abs_rij
  real(DP), intent(in) :: U_i 
  real(DP)             :: res
  real(DP)             :: uir

  uir = U_i*abs_rij
  res = 33.0_DP + (18.0_DP + 3.0_DP*uir)*uir
  res = res/48.0_DP

  endfunction part_deriv_gij_wrt_Ui

  function second_part_deriv_gij_wrt_Ui_and_r(abs_rij, U_i) result(res)
  real(DP), intent(in) :: abs_rij
  real(DP), intent(in) :: U_i 
  real(DP)             :: res

  res = (3.0_DP + U_i*abs_rij)*U_i/8.0_DP

  endfunction second_part_deriv_gij_wrt_Ui_and_r

  function second_part_deriv_fij_wrt_Ui_and_r(abs_rij, U_i, U_j) result(res)
  real(DP), intent(in) :: abs_rij
  real(DP), intent(in) :: U_i, U_j
  real(DP)             :: term1, term2
  real(DP)             :: res
 
  term1 = 12.0_DP*U_i**3*U_j**4 
  term2 = (U_i**2 - U_j**2)**4*abs_rij**2

  res = term1/term2

  endfunction second_part_deriv_fij_wrt_Ui_and_r

  function second_part_deriv_fij_wrt_Uj_and_r(abs_rij, U_i, U_j) result(res)
  real(DP), intent(in) :: abs_rij
  real(DP), intent(in) :: U_i, U_j
  real(DP)             :: term1, term2
  real(DP)             :: res
 
  term1 = - 12.0_DP*U_i**4*U_j**3 
  term2 = (U_i**2 - U_j**2)**4*abs_rij**2

  res = term1/term2

  endfunction second_part_deriv_fij_wrt_Uj_and_r

endmodule coulomb_short_gamma

