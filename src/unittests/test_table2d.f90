module test_table2d
  use fruit
  use system_module
  use table2d

  implicit none

contains

  subroutine test_table2d_f_and_df
    reaL(DP), parameter :: tol = 1e-6_DP
    type(table2d_t)     :: tab
    real(DP)            :: vals(0:2, 0:3) = reshape([1.0_DP, 2.0_DP, 3.0_DP, &
                                                     4.0_DP, 5.0_DP, 6.0_DP, &
                                                     7.0_DP, 8.0_DP, 9.0_DP, &
                                                     10.0_DP, 11.0_DP, &
                                                     12.0_DP], &
                                                     [3,4])
    integer             :: i, j
    real(DP)            :: val, val1, val2, dvaldi, dvaldj, ti, tj

    call init(tab, 2, 3, vals)
    do i = 0, 2
       do j = 0, 3
          call eval(tab, 1.0_DP*i, 1.0_DP*j, val, dvaldi, dvaldj)
          call assert_equals(vals(i, j), val, tol, "val")
          call assert_equals(0.0_DP, dvaldi, tol, "dvaldi")
          call assert_equals(0.0_DP, dvaldj, tol, "dvaldj")
       enddo
    enddo
    call del(tab)

    call init(tab, 2, 3, vals, 2*vals, 3*vals)
    do i = 0, 2
       do j = 0, 3
          call eval(tab, 1.0_DP*i, 1.0_DP*j, val, dvaldi, dvaldj)
          call assert_equals(vals(i, j), val, tol, "val")
          call assert_equals(2*vals(i, j), dvaldi, tol, "dvaldi")
          call assert_equals(3*vals(i, j), dvaldj, tol, "dvaldj")
       enddo
    enddo

    do i = 0, 22
       do j = 0, 33
          call eval(tab, 0.1_DP*i, 0.1_DP*j, val, dvaldi, dvaldj)
          call eval(tab, 0.1_DP*i-tol, 0.1_DP*j, val1, ti, tj)
          call eval(tab, 0.1_DP*i+tol, 0.1_DP*j, val2, ti, tj)
          call assert_equals((val2-val1)/(2*tol), dvaldi, tol*100, "dvaldi")
          call eval(tab, 0.1_DP*i, 0.1_DP*j-tol, val1, ti, tj)
          call eval(tab, 0.1_DP*i, 0.1_DP*j+tol, val2, ti, tj)
          call assert_equals((val2-val1)/(2*tol), dvaldj, tol*100, "dvaldj")
       enddo
    enddo
    call del(tab)
  endsubroutine test_table2d_f_and_df

endmodule test_table2d
