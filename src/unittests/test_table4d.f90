module test_table4d
  use fruit
  use system_module
  use table4d

  implicit none

contains

  subroutine test_table4d_f_and_df
    reaL(DP), parameter :: tol = 1e-7_DP
    type(table4d_t)     :: tab
    real(DP)            :: vals(0:1, 0:2, 0:3, 0:1) = &
        reshape([1.0_DP, 2.0_DP, 3.0_DP, 4.0_DP, 5.0_DP, 6.0_DP, &
                 7.0_DP, 8.0_DP, 9.0_DP, 10.0_DP, 11.0_DP, 12.0_DP, &
                 13.0_DP, 14.0_DP, 15.0_DP, 16.0_DP, 17.0_DP, 18.0_DP, &
                 19.0_DP, 20.0_DP, 21.0_DP, 22.0_DP, 23.0_DP, 24.0_DP, &
                 41.0_DP, 42.0_DP, 43.0_DP, 44.0_DP, 45.0_DP, 46.0_DP, &
                 47.0_DP, 48.0_DP, 49.0_DP, 410.0_DP, 411.0_DP, 412.0_DP, &
                 413.0_DP, 414.0_DP, 415.0_DP, 416.0_DP, 417.0_DP, 418.0_DP, &
                 419.0_DP, 420.0_DP, 421.0_DP, 422.0_DP, 423.0_DP, 424.0_DP], &
                [2,3,4,2])
    integer  :: i, j, k, l
    real(DP) :: val, val1, val2, dvaldi, dvaldj, dvaldk, dvaldl, ti, tj, tk, tl

    call init(tab, 1, 2, 3, 1, vals)
    do i = 0, 1
       do j = 0, 2
          do k = 0, 3
             do l = 0, 1
                call eval(tab, 1.0_DP*i, 1.0_DP*j, 1.0_DP*k, 1.0_DP*l, &
                          val, dvaldi, dvaldj, dvaldk, dvaldl)
                call assert_equals(vals(i, j, k, l), val, tol, "table4d|val")
                call assert_equals(0.0_DP, dvaldi, tol, "table4d|dvaldi")
                call assert_equals(0.0_DP, dvaldj, tol, "table4d|dvaldj")
                call assert_equals(0.0_DP, dvaldk, tol, "table4d|dvaldk")
                call assert_equals(0.0_DP, dvaldl, tol, "table4d|dvaldl")
             enddo
          enddo
       enddo
    enddo
    call del(tab)

    call init(tab, 1, 2, 3, 1, vals, 2*vals, 3*vals, 0.25_DP*vals, 0.5_DP*vals)
    do i = 0, 1
       do j = 0, 2
          do k = 0, 3
             do l = 0, 1
                call eval(tab, 1.0_DP*i, 1.0_DP*j, 1.0_DP*k, 1.0_DP*l, &
                          val, dvaldi, dvaldj, dvaldk, dvaldl)
                call assert_equals(vals(i, j, k, l), val, tol, "table4d|val")
                call assert_equals(2*vals(i, j, k, l), dvaldi, tol, "table4d|dvaldi")
                call assert_equals(3*vals(i, j, k, l), dvaldj, tol, "table4d|dvaldj")
                call assert_equals(0.25_DP*vals(i, j, k, l), dvaldk, tol, "table4d|dvaldk")
                call assert_equals(0.5_DP*vals(i, j, k, l), dvaldl, tol, "table4d|dvaldl")
             enddo
          enddo
       enddo
    enddo

    do i = 0, 11
       do j = 0, 22
          do k = 0, 33
             do l = 0, 44
                 call eval(tab, 0.1_DP*i, 0.1_DP*j, 0.1_DP*k, 0.1_DP*l, &
                           val, dvaldi, dvaldj, dvaldk, dvaldl)
                 call eval(tab, 0.1_DP*i-tol, 0.1_DP*j, 0.1_DP*k, 0.1_DP*l, val1, ti, tj, tk, tl)
                 call eval(tab, 0.1_DP*i+tol, 0.1_DP*j, 0.1_DP*k, 0.1_DP*l, val2, ti, tj, tk, tl)
                 call assert_equals((val2-val1)/(2*tol), dvaldi, 0.1_DP, "table4d|dvaldi")
                 call eval(tab, 0.1_DP*i, 0.1_DP*j-tol, 0.1_DP*k, 0.1_DP*l, val1, ti, tj, tk, tl)
                 call eval(tab, 0.1_DP*i, 0.1_DP*j+tol, 0.1_DP*k, 0.1_DP*l, val2, ti, tj, tk, tl)
                 call assert_equals((val2-val1)/(2*tol), dvaldj, 0.1_DP, "table4d|dvaldj")
                 call eval(tab, 0.1_DP*i, 0.1_DP*j, 0.1_DP*k-tol, 0.1_DP*l, val1, ti, tj, tk, tl)
                 call eval(tab, 0.1_DP*i, 0.1_DP*j, 0.1_DP*k+tol, 0.1_DP*l, val2, ti, tj, tk, tl)
                 call assert_equals((val2-val1)/(2*tol), dvaldk, 0.1_DP, "table4d|dvaldk")
                 call eval(tab, 0.1_DP*i, 0.1_DP*j, 0.1_DP*k, 0.1_DP*l-tol, val1, ti, tj, tk, tl)
                 call eval(tab, 0.1_DP*i, 0.1_DP*j, 0.1_DP*k, 0.1_DP*l+tol, val2, ti, tj, tk, tl)
                 call assert_equals((val2-val1)/(2*tol), dvaldl, 0.1_DP, "table4d|dvaldl")
              enddo
          enddo
       enddo
    enddo
    call del(tab)
  endsubroutine test_table4d_f_and_df

endmodule test_table4d
