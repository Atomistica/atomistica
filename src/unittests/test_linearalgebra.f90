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

  subroutine test_sqrtm
    real(DP), parameter :: tol = 1e-6
    real(DP) :: a(3, 3), b(3, 3)
    a = 0.0_DP
    a(1, 1) = 1.0_DP
    a(2, 2) = 2.0_DP
    a(3, 3) = 3.0_DP
    a(1, 2) = 1.2_DP
    a(2, 1) = 1.2_DP
    a(1, 3) = 0.9_DP
    a(3, 1) = 0.9_DP
    b = sqrtm(a)
    call assert_equals(a, matmul(b, b), 3, 3, tol, "sqrtm")
  endsubroutine

endmodule test_linearalgebra
