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

  subroutine test_gauss
    real(DP), parameter :: tol = 1e-6
    real(DP) :: A(3, 4), B(3, 4), x(3), y(3), C(3, 3), z(3, 1)
    A = 0.0_DP
    A(1, 1) = 1.0_DP
    A(2, 2) = 1.0_DP
    A(1, 2) = 2.0_DP
    A(2, 1) = 0.3_DP
    A(3, 2) = 1.5_DP
    A(3, 3) = 4.0_DP
    A(2, 3) = 3.5_DP
    A(1, 4) = 1.0_DP
    A(2, 4) = 1.0_DP
    A(3, 4) = 2.0_DP
    B = A
    call gauss(3, A, x)
    y = matmul(B(:, :3), x)-B(:, 4)
    call assert_equals(y(1), 0.0_DP, tol, "gauss")
    call assert_equals(y(2), 0.0_DP, tol, "gauss")
    call assert_equals(y(3), 0.0_DP, tol, "gauss")
    C = B(1:3, 1:3)
    z(1:3, 1) = B(1:3, 4)
    call gaussn(3, C, 1, z)
    y = matmul(B(:, :3), Z(1:3, 1))-B(:, 4)    
    call assert_equals(y(1), 0.0_DP, tol, "gauss")
    call assert_equals(y(2), 0.0_DP, tol, "gauss")
    call assert_equals(y(3), 0.0_DP, tol, "gauss")
  endsubroutine test_gauss

  subroutine test_gauss_inverse
    real(DP), parameter :: tol = 1e-6
    real(DP) :: A(3, 3), B(3, 3), C(3, 3)
    integer :: error
    A = 0.0_DP
    A(1, 1) = 3
    A(2, 1) = 2
    A(1, 3) = 2
    A(2, 3) = -2
    A(3, 2) = 1
    A(3, 3) = 1
    call identity(3, B)
    call gaussn(3, A, 3, B)
    C = matmul(A, B)
    call assert_equals(C(1, 1), 1.0_DP, tol, "gauss_inverse")
    call assert_equals(C(2, 2), 1.0_DP, tol, "gauss_inverse")
    call assert_equals(C(3, 3), 1.0_DP, tol, "gauss_inverse")
    call assert_equals(C(1, 2), 0.0_DP, tol, "gauss_inverse")
    call assert_equals(C(1, 3), 0.0_DP, tol, "gauss_inverse")
    call assert_equals(C(2, 3), 0.0_DP, tol, "gauss_inverse")
    call assert_equals(C(2, 1), 0.0_DP, tol, "gauss_inverse")
    call assert_equals(C(3, 1), 0.0_DP, tol, "gauss_inverse")
    call assert_equals(C(3, 2), 0.0_DP, tol, "gauss_inverse")
    ! The following matrix does not have an inverse. Gauss should fail.
    A = 0.0_DP
    A(1, 1) = 1
    A(1, 2) = 6
    A(1, 3) = 4
    A(2, 1) = 2
    A(2, 2) = 4
    A(2, 3) = -1
    A(3, 1) = -1
    A(3, 2) = 2
    A(3, 3) = 5
    call identity(3, B)
    call gaussn(3, A, 3, B, error=error)
    call assert_true(error /= 0)
  endsubroutine test_gauss_inverse

endmodule test_linearalgebra
