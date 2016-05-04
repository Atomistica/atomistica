program run_tests
  use fruit
  use test_cutoff
  use test_linearalgebra
  use test_table2d
  use test_table3d
  use test_table4d

  implicit none

  call init_fruit

  ! test_cutoff
  call test_exp_cutoff
  call test_trig_on
  call test_trig_off

  ! test_linearalgebra
  call test_det
  call test_sqrtm

  ! test_table2d
  call test_table2d_f_and_df

  ! test_table3d
  call test_table3d_f_and_df

  ! test_table4d
  call test_table4d_f_and_df

  call fruit_summary
  call fruit_finalize
endprogram
