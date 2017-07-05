module test_table4d
  use fruit
  use system_module
  use table4d

  implicit none

contains

  subroutine test_table4d_f_and_df_single_dimension
    real(DP), parameter :: tol = 1e-7_DP
    real(DP), parameter :: eps = 1e-7_DP
    type(table4d_t)     :: tab3111, tab1311, tab1131, tab1113
    real(DP)            :: vals(0:3) = [1.0_DP, 2.0_DP, 3.0_DP, 4.0_DP]
    integer  :: i, j, k, l
    real(DP) :: val, val1, val2, dvaldi, dvaldj, dvaldk, dvaldl, ti, tj, tk, tl

    real(DP) :: vals3111(0:3, 0:1, 0:1, 0:1)
    real(DP) :: vals1311(0:1, 0:3, 0:1, 0:1)
    real(DP) :: vals1131(0:1, 0:1, 0:3, 0:1)
    real(DP) :: vals1113(0:1, 0:1, 0:1, 0:3)

    do i = 0, 1
       do j = 0, 1
          do k = 0, 1
             vals3111(0:3, i, j, k) = vals
             vals1311(i, 0:3, j, k) = vals
             vals1131(i, j, 0:3, k) = vals
             vals1113(i, j, k, 0:3) = vals
          enddo
       enddo
    enddo

    call init(tab3111, 3, 1, 1, 1, vals3111)
    call init(tab1311, 1, 3, 1, 1, vals1311)
    call init(tab1131, 1, 1, 3, 1, vals1131)
    call init(tab1113, 1, 1, 1, 3, vals1113)
!    do i = 0, 1
!       do j = 0, 1
!          do k = 0, 1
!             do l = 0, 3
!                call eval(tab1113, 1.0_DP*i, 1.0_DP*j, 1.0_DP*k, 1.0_DP*l, &
!                          val, dvaldi, dvaldj, dvaldk, dvaldl)
!                call assert_equals(vals(i, j, k, l), val, tol, "table4d|val")
!                call assert_equals(0.0_DP, dvaldi, tol, "table4d|dvaldi")
!                call assert_equals(0.0_DP, dvaldj, tol, "table4d|dvaldj")
!                call assert_equals(0.0_DP, dvaldk, tol, "table4d|dvaldk")
!                call assert_equals(0.0_DP, dvaldl, tol, "table4d|dvaldl")
!             enddo
!          enddo
!       enddo
!    enddo
!    call del(tab)

    call init(tab3111, 3, 1, 1, 1, vals3111, 2*vals3111, 3*vals3111, 0.25_DP*vals3111, 0.5_DP*vals3111)
    call init(tab1311, 1, 3, 1, 1, vals1311, 2*vals1311, 3*vals1311, 0.25_DP*vals1311, 0.5_DP*vals1311)
    call init(tab1131, 1, 1, 3, 1, vals1131, 2*vals1131, 3*vals1131, 0.25_DP*vals1131, 0.5_DP*vals1131)
    call init(tab1113, 1, 1, 1, 3, vals1113, 2*vals1113, 3*vals1113, 0.25_DP*vals1113, 0.5_DP*vals1113)
    do i = 0, 1
       do j = 0, 1
          do k = 0, 1
             do l = 0, 3
                call eval(tab1113, 1.0_DP*i, 1.0_DP*j, 1.0_DP*k, 1.0_DP*l, &
                          val, dvaldi, dvaldj, dvaldk, dvaldl)
                call assert_equals(vals1113(i, j, k, l), val, tol, "table4d|val")
                call assert_equals(2*vals1113(i, j, k, l), dvaldi, tol, "table4d|dvaldi")
                call assert_equals(3*vals1113(i, j, k, l), dvaldj, tol, "table4d|dvaldj")
                call assert_equals(0.25_DP*vals1113(i, j, k, l), dvaldk, tol, "table4d|dvaldk")
                call assert_equals(0.5_DP*vals1113(i, j, k, l), dvaldl, tol, "table4d|dvaldl")
             enddo
          enddo
       enddo
    enddo

    do i = 0, 11
       do j = 0, 11
          do k = 0, 11
             do l = 0, 33
                 call eval(tab1113, 0.1_DP*i, 0.1_DP*j, 0.1_DP*k, 0.1_DP*l, &
                           val, dvaldi, dvaldj, dvaldk, dvaldl)
                 call eval(tab1113, 0.1_DP*i-eps, 0.1_DP*j, 0.1_DP*k, 0.1_DP*l, val1, ti, tj, tk, tl)
                 call eval(tab1113, 0.1_DP*i+eps, 0.1_DP*j, 0.1_DP*k, 0.1_DP*l, val2, ti, tj, tk, tl)
                 call assert_equals((val2-val1)/(2*tol), dvaldi, tol, "table4d|dvaldi")
                 call eval(tab1113, 0.1_DP*i, 0.1_DP*j-eps, 0.1_DP*k, 0.1_DP*l, val1, ti, tj, tk, tl)
                 call eval(tab1113, 0.1_DP*i, 0.1_DP*j+eps, 0.1_DP*k, 0.1_DP*l, val2, ti, tj, tk, tl)
                 call assert_equals((val2-val1)/(2*tol), dvaldj, tol, "table4d|dvaldj")
                 call eval(tab1113, 0.1_DP*i, 0.1_DP*j, 0.1_DP*k-eps, 0.1_DP*l, val1, ti, tj, tk, tl)
                 call eval(tab1113, 0.1_DP*i, 0.1_DP*j, 0.1_DP*k+eps, 0.1_DP*l, val2, ti, tj, tk, tl)
                 call assert_equals((val2-val1)/(2*tol), dvaldk, tol, "table4d|dvaldk")
                 call eval(tab1113, 0.1_DP*i, 0.1_DP*j, 0.1_DP*k, 0.1_DP*l-eps, val1, ti, tj, tk, tl)
                 call eval(tab1113, 0.1_DP*i, 0.1_DP*j, 0.1_DP*k, 0.1_DP*l+eps, val2, ti, tj, tk, tl)
                 call assert_equals((val2-val1)/(2*tol), dvaldl, tol, "table4d|dvaldl")
              enddo
          enddo
       enddo
    enddo

    do i = 0, 11
       do j = 0, 11
          do k = 0, 33
             do l = 0, 11
                 call eval(tab1131, 0.1_DP*i, 0.1_DP*j, 0.1_DP*k, 0.1_DP*l, &
                           val, dvaldi, dvaldj, dvaldk, dvaldl)
                 call eval(tab1131, 0.1_DP*i-eps, 0.1_DP*j, 0.1_DP*k, 0.1_DP*l, val1, ti, tj, tk, tl)
                 call eval(tab1131, 0.1_DP*i+eps, 0.1_DP*j, 0.1_DP*k, 0.1_DP*l, val2, ti, tj, tk, tl)
                 call assert_equals((val2-val1)/(2*tol), dvaldi, tol, "table4d|dvaldi")
                 call eval(tab1131, 0.1_DP*i, 0.1_DP*j-eps, 0.1_DP*k, 0.1_DP*l, val1, ti, tj, tk, tl)
                 call eval(tab1131, 0.1_DP*i, 0.1_DP*j+eps, 0.1_DP*k, 0.1_DP*l, val2, ti, tj, tk, tl)
                 call assert_equals((val2-val1)/(2*tol), dvaldj, tol, "table4d|dvaldj")
                 call eval(tab1131, 0.1_DP*i, 0.1_DP*j, 0.1_DP*k-eps, 0.1_DP*l, val1, ti, tj, tk, tl)
                 call eval(tab1131, 0.1_DP*i, 0.1_DP*j, 0.1_DP*k+eps, 0.1_DP*l, val2, ti, tj, tk, tl)
                 call assert_equals((val2-val1)/(2*tol), dvaldk, tol, "table4d|dvaldk")
                 call eval(tab1131, 0.1_DP*i, 0.1_DP*j, 0.1_DP*k, 0.1_DP*l-eps, val1, ti, tj, tk, tl)
                 call eval(tab1131, 0.1_DP*i, 0.1_DP*j, 0.1_DP*k, 0.1_DP*l+eps, val2, ti, tj, tk, tl)
                 call assert_equals((val2-val1)/(2*tol), dvaldl, tol, "table4d|dvaldl")
              enddo
          enddo
       enddo
    enddo

    do i = 0, 11
       do j = 0, 33
          do k = 0, 11
             do l = 0, 11
                 call eval(tab1311, 0.1_DP*i, 0.1_DP*j, 0.1_DP*k, 0.1_DP*l, &
                           val, dvaldi, dvaldj, dvaldk, dvaldl)
                 call eval(tab1311, 0.1_DP*i-eps, 0.1_DP*j, 0.1_DP*k, 0.1_DP*l, val1, ti, tj, tk, tl)
                 call eval(tab1311, 0.1_DP*i+eps, 0.1_DP*j, 0.1_DP*k, 0.1_DP*l, val2, ti, tj, tk, tl)
                 call assert_equals((val2-val1)/(2*tol), dvaldi, tol, "table4d|dvaldi")
                 call eval(tab1311, 0.1_DP*i, 0.1_DP*j-eps, 0.1_DP*k, 0.1_DP*l, val1, ti, tj, tk, tl)
                 call eval(tab1311, 0.1_DP*i, 0.1_DP*j+eps, 0.1_DP*k, 0.1_DP*l, val2, ti, tj, tk, tl)
                 call assert_equals((val2-val1)/(2*tol), dvaldj, tol, "table4d|dvaldj")
                 call eval(tab1311, 0.1_DP*i, 0.1_DP*j, 0.1_DP*k-eps, 0.1_DP*l, val1, ti, tj, tk, tl)
                 call eval(tab1311, 0.1_DP*i, 0.1_DP*j, 0.1_DP*k+eps, 0.1_DP*l, val2, ti, tj, tk, tl)
                 call assert_equals((val2-val1)/(2*tol), dvaldk, tol, "table4d|dvaldk")
                 call eval(tab1311, 0.1_DP*i, 0.1_DP*j, 0.1_DP*k, 0.1_DP*l-eps, val1, ti, tj, tk, tl)
                 call eval(tab1311, 0.1_DP*i, 0.1_DP*j, 0.1_DP*k, 0.1_DP*l+eps, val2, ti, tj, tk, tl)
                 call assert_equals((val2-val1)/(2*tol), dvaldl, tol, "table4d|dvaldl")
              enddo
          enddo
       enddo
    enddo

    do i = 0, 33
       do j = 0, 11
          do k = 0, 11
             do l = 0, 11
                 call eval(tab3111, 0.1_DP*i, 0.1_DP*j, 0.1_DP*k, 0.1_DP*l, &
                           val, dvaldi, dvaldj, dvaldk, dvaldl)
                 call eval(tab3111, 0.1_DP*i-eps, 0.1_DP*j, 0.1_DP*k, 0.1_DP*l, val1, ti, tj, tk, tl)
                 call eval(tab3111, 0.1_DP*i+eps, 0.1_DP*j, 0.1_DP*k, 0.1_DP*l, val2, ti, tj, tk, tl)
                 call assert_equals((val2-val1)/(2*eps), dvaldi, tol, "table4d|dvaldi")
                 call eval(tab3111, 0.1_DP*i, 0.1_DP*j-eps, 0.1_DP*k, 0.1_DP*l, val1, ti, tj, tk, tl)
                 call eval(tab3111, 0.1_DP*i, 0.1_DP*j+eps, 0.1_DP*k, 0.1_DP*l, val2, ti, tj, tk, tl)
                 call assert_equals((val2-val1)/(2*eps), dvaldj, tol, "table4d|dvaldj")
                 call eval(tab3111, 0.1_DP*i, 0.1_DP*j, 0.1_DP*k-eps, 0.1_DP*l, val1, ti, tj, tk, tl)
                 call eval(tab3111, 0.1_DP*i, 0.1_DP*j, 0.1_DP*k+eps, 0.1_DP*l, val2, ti, tj, tk, tl)
                 call assert_equals((val2-val1)/(2*eps), dvaldk, tol, "table4d|dvaldk")
                 call eval(tab3111, 0.1_DP*i, 0.1_DP*j, 0.1_DP*k, 0.1_DP*l-eps, val1, ti, tj, tk, tl)
                 call eval(tab3111, 0.1_DP*i, 0.1_DP*j, 0.1_DP*k, 0.1_DP*l+eps, val2, ti, tj, tk, tl)
                 call assert_equals((val2-val1)/(2*eps), dvaldl, tol, "table4d|dvaldl")
              enddo
          enddo
       enddo
    enddo

    call del(tab3111)
    call del(tab1311)
    call del(tab1131)
    call del(tab1113)
  endsubroutine test_table4d_f_and_df_single_dimension

  subroutine test_table4d_f_and_df
    reaL(DP), parameter :: tol = 1e-7_DP
    reaL(DP), parameter :: tol2 = 1e-3_DP
    reaL(DP), parameter :: eps = 5e-7_DP
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
                 if (dvaldi > 1e-9_DP) then
                    call eval(tab, 0.1_DP*i-eps, 0.1_DP*j, 0.1_DP*k, 0.1_DP*l, val1, ti, tj, tk, tl)
                    call eval(tab, 0.1_DP*i+eps, 0.1_DP*j, 0.1_DP*k, 0.1_DP*l, val2, ti, tj, tk, tl)
                    call assert_equals(1.0_DP, (val2-val1)/(2*eps) / dvaldi, tol2, "table4d|dvaldi")
                 endif
                 if (dvaldj > 1e-9_DP) then
                    call eval(tab, 0.1_DP*i, 0.1_DP*j-eps, 0.1_DP*k, 0.1_DP*l, val1, ti, tj, tk, tl)
                    call eval(tab, 0.1_DP*i, 0.1_DP*j+eps, 0.1_DP*k, 0.1_DP*l, val2, ti, tj, tk, tl)
                    call assert_equals(1.0_DP, (val2-val1)/(2*eps) / dvaldj, tol2, "table4d|dvaldj")
                 endif
                 if (dvaldk > 1e-9_DP) then
                    call eval(tab, 0.1_DP*i, 0.1_DP*j, 0.1_DP*k-eps, 0.1_DP*l, val1, ti, tj, tk, tl)
                    call eval(tab, 0.1_DP*i, 0.1_DP*j, 0.1_DP*k+eps, 0.1_DP*l, val2, ti, tj, tk, tl)
                    call assert_equals(1.0_DP, (val2-val1)/(2*eps) / dvaldk, tol2, "table4d|dvaldk")
                 endif
                 if (dvaldl > 1e-9_DP) then
                    call eval(tab, 0.1_DP*i, 0.1_DP*j, 0.1_DP*k, 0.1_DP*l-eps, val1, ti, tj, tk, tl)
                    call eval(tab, 0.1_DP*i, 0.1_DP*j, 0.1_DP*k, 0.1_DP*l+eps, val2, ti, tj, tk, tl)
                    call assert_equals(1.0_DP, (val2-val1)/(2*eps) / dvaldl, tol2, "table4d|dvaldl")
                 endif
              enddo
          enddo
       enddo
    enddo
    call del(tab)
  endsubroutine test_table4d_f_and_df

endmodule test_table4d
