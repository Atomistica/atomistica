module test_table2d
  use fruit
  use system_module
  use table2d

  implicit none

contains

  subroutine test_table2d_f_and_df
    reaL(DP), parameter :: tol = 1e-6_DP
    type(table2d_t)     :: tab
    real(DP)            :: vals(0:2, 0:2) = reshape([1.0_DP, 2.0_DP, 3.0_DP, &
                                                     4.0_DP, 5.0_DP, 6.0_DP, &
                                                     7.0_DP, 8.0_DP, 9.0_DP], &
                                                     [3,3])
    integer             :: i, j
    real(DP)            :: val, dvaldi, dvaldj

    call init(tab, 2, 2, vals)
    do i = 0, 2
       do j = 0, 2
          call eval(tab, 1.0_DP*i, 1.0_DP*j, val, dvaldi, dvaldj)
          call assert_equals(vals(i, j), val, tol, "val")
          call assert_equals(0.0_DP, dvaldi, tol, "dvaldi")
          call assert_equals(0.0_DP, dvaldj, tol, "dvaldj")
       enddo
    enddo
    call del(tab)

    call init(tab, 2, 2, vals, vals)
    do i = 0, 2
       do j = 0, 2
          call eval(tab, 1.0_DP*i, 1.0_DP*j, val, dvaldi, dvaldj)
          call assert_equals(vals(i, j), val, tol, "val")
          call assert_equals(vals(i, j), dvaldi, tol, "dvaldi")
          call assert_equals(0.0_DP, dvaldj, tol, "dvaldj")
       enddo
    enddo
    call del(tab)
  endsubroutine test_table2d_f_and_df

endmodule test_table2d
