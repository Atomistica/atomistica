module test_linearalgebra
  use fruit
  use system_module
  use linearalgebra

  implicit none

contains

  subroutine test_det
    real(DP), parameter :: tol = 1e-6
    real(DP) :: a(2, 2)
    a = 0.0_DP
    a(1, 1) = 1.0_DP
    a(2, 2) = 2.0_DP
    call assert_equals(det(a), 2.0_DP, tol, "diagonal matrix")
    a(1, 2) = 3.0_DP
    call assert_equals(det(a), 2.0_DP, tol, "triangular matrix")
    a(2, 1) = 4.0_DP
    call assert_equals(det(a), -10.0_DP, tol, "full matrix")
  endsubroutine test_det

endmodule test_linearalgebra
