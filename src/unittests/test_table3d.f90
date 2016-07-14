module test_table3d
  use fruit
  use system_module
  use table3d

  implicit none

contains

  subroutine test_table3d_f_and_df
    reaL(DP), parameter :: tol = 1e-6_DP
    type(table3d_t)     :: tab
    real(DP)            :: vals(0:1, 0:2, 0:3) = &
        reshape([1.0_DP, 2.0_DP, 3.0_DP, 4.0_DP, 5.0_DP, 6.0_DP, &
                 7.0_DP, 8.0_DP, 9.0_DP, 10.0_DP, 11.0_DP, 12.0_DP, &
                 13.0_DP, 14.0_DP, 15.0_DP, 16.0_DP, 17.0_DP, 18.0_DP, &
                 19.0_DP, 20.0_DP, 21.0_DP, 22.0_DP, 23.0_DP, 24.0_DP], &
                [2,3,4])
    integer             :: i, j, k
    real(DP)            :: val, val1, val2, dvaldi, dvaldj, dvaldk, ti, tj, tk

    call init(tab, 1, 2, 3, vals)
    do i = 0, 1
       do j = 0, 2
          do k = 0, 3
             call eval(tab, 1.0_DP*i, 1.0_DP*j, 1.0_DP*k, val, dvaldi, dvaldj, &
                       dvaldk)
             call assert_equals(vals(i, j, k), val, tol, "table3d|val")
             call assert_equals(0.0_DP, dvaldi, tol, "table3d|dvaldi")
             call assert_equals(0.0_DP, dvaldj, tol, "table3d|dvaldj")
             call assert_equals(0.0_DP, dvaldk, tol, "table3d|dvaldk")
          enddo
       enddo
    enddo
    call del(tab)

    call init(tab, 1, 2, 3, vals, 2*vals, 3*vals, 4*vals)
    do i = 0, 1
       do j = 0, 2
          do k = 0, 3
             call eval(tab, 1.0_DP*i, 1.0_DP*j, 1.0_DP*k, val, dvaldi, dvaldj, &
                       dvaldk)
             call assert_equals(vals(i, j, k), val, tol, "table3d|val")
             call assert_equals(2*vals(i, j, k), dvaldi, tol, "table3d|dvaldi")
             call assert_equals(3*vals(i, j, k), dvaldj, tol, "table3d|dvaldj")
             call assert_equals(4*vals(i, j, k), dvaldk, tol, "table3d|dvaldk")
          enddo
       enddo
    enddo

    do i = 0, 11
       do j = 0, 22
          do k = 0, 33
             call eval(tab, 0.1_DP*i, 0.1_DP*j, 0.1_DP*k, val, dvaldi, dvaldj, &
                       dvaldk)
             call eval(tab, 0.1_DP*i-tol, 0.1_DP*j, 0.1_DP*k, val1, ti, tj, tk)
             call eval(tab, 0.1_DP*i+tol, 0.1_DP*j, 0.1_DP*k, val2, ti, tj, tk)
             call assert_equals((val2-val1)/(2*tol), dvaldi, tol*1000, "table3d|dvaldi")
             call eval(tab, 0.1_DP*i, 0.1_DP*j-tol, 0.1_DP*k, val1, ti, tj, tk)
             call eval(tab, 0.1_DP*i, 0.1_DP*j+tol, 0.1_DP*k, val2, ti, tj, tk)
             call assert_equals((val2-val1)/(2*tol), dvaldj, tol*1000, "table3d|dvaldj")
             call eval(tab, 0.1_DP*i, 0.1_DP*j, 0.1_DP*k-tol, val1, ti, tj, tk)
             call eval(tab, 0.1_DP*i, 0.1_DP*j, 0.1_DP*k+tol, val2, ti, tj, tk)
             call assert_equals((val2-val1)/(2*tol), dvaldk, tol*1000, "table3d|dvaldk")
          enddo
       enddo
    enddo
    call del(tab)
  endsubroutine test_table3d_f_and_df

endmodule test_table3d
