module test_cutoff
  use fruit
  use system_module
  use cutoff

  implicit none

contains

#define MAKE_CUTOFF_TEST(name, cutoff_t)                                 \
  subroutine name ;                                                      \
    type(cutoff_t) :: cutoff ;                                           \
    real(DP), parameter :: r1 = 1.5_DP ;                                 \
    real(DP), parameter :: r2 = 2.75_DP ;                                \
    real(DP), parameter :: tol = 1d-6 ;                                  \
    real(DP), parameter :: dr = 1e-6 ;                                   \
    real(DP) :: val, dval, val2, dval2 ;                                 \
    integer :: i ;                                                       \
    call init(cutoff, r1, r2) ;                                          \
    call fc(cutoff, r1, val, dval) ;                                     \
    call assert_equals(1.0_DP, val, tol) ;                               \
    call assert_equals(0.0_DP, dval, tol) ;                              \
    call fc(cutoff, r2, val, dval) ;                                     \
    call assert_equals(0.0_DP, val, tol) ;                               \
    call assert_equals(0.0_DP, dval, tol) ;                              \
    do i = 0, 100 ;                                                      \
       call fc(cutoff, r1+i*(r2-r1)/100.0_DP, val, dval) ;               \
       call assert_true(val >= 0.0_DP .and. val <= 1.0_DP) ;             \
       call assert_true(dval <= 0.0_DP) ;                                \
    enddo ;                                                              \
    do i = 0, 99 ;                                                       \
       call fc(cutoff, r1+i*(r2-r1)/1000.0_DP, val, dval) ;              \
       call fc(cutoff, r1+i*(r2-r1)/1000.0_DP+dr, val2, dval2) ;         \
       call assert_true(abs((val2-val)/dr-0.5_DP*(dval+dval2)) < tol) ;  \
    enddo ;                                                              \
  endsubroutine name

MAKE_CUTOFF_TEST(test_exp_cutoff, exp_cutoff_t)
MAKE_CUTOFF_TEST(test_trig_cutoff, trig_cutoff_t)

endmodule test_cutoff
