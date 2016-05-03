module test_table3d
  use fruit
  use system_module
  use table3d

  implicit none

contains

  subroutine test_table3d_f_and_df
    reaL(DP), parameter :: tol = 1e-6_DP
    type(table3d_t)     :: tab
    real(DP)            :: vals(0:2, 0:2, 0:2) = &
        reshape([1.0_DP, 2.0_DP, 3.0_DP, 4.0_DP, 5.0_DP, 6.0_DP, &
                 7.0_DP, 8.0_DP, 9.0_DP, 10.0_DP, 11.0_DP, 12.0_DP, &
                 13.0_DP, 14.0_DP, 15.0_DP, 16.0_DP, 17.0_DP, 18.0_DP, &
                 19.0_DP, 20.0_DP, 21.0_DP, 22.0_DP, 23.0_DP, 24.0_DP, &
                 25.0_DP, 26.0_DP, 27.0_DP, 28.0_DP, 29.0_DP], &
                [3,3,3])
    integer             :: i, j, k
    real(DP)            :: val, dvaldi, dvaldj, dvaldk

    call init(tab, 2, 2, 2, vals)
    do i = 0, 2
       do j = 0, 2
          do k = 0, 2
             call eval(tab, 1.0_DP*i, 1.0_DP*j, 1.0_DP*k, val, dvaldi, dvaldj, &
                       dvaldk)
             call assert_equals(vals(i, j, k), val, tol, "val")
             call assert_equals(0.0_DP, dvaldi, tol, "dvaldi")
             call assert_equals(0.0_DP, dvaldj, tol, "dvaldj")
             call assert_equals(0.0_DP, dvaldk, tol, "dvaldk")
          enddo
       enddo
    enddo
    call del(tab)

    call init(tab, 2, 2, 2, vals, vals)
    do i = 0, 2
       do j = 0, 2
          do k = 0, 2
             call eval(tab, 1.0_DP*i, 1.0_DP*j, 1.0_DP*k, val, dvaldi, dvaldj, &
                       dvaldk)
             call assert_equals(vals(i, j, k), val, tol, "val")
             call assert_equals(vals(i, j, k), dvaldi, tol, "dvaldi")
             call assert_equals(0.0_DP, dvaldj, tol, "dvaldj")
             call assert_equals(0.0_DP, dvaldk, tol, "dvaldk")
          enddo
       enddo
    enddo
    call del(tab)
  endsubroutine test_table3d_f_and_df

endmodule test_table3d
