program run_tests
  use fruit
  use test_cutoff
  use test_table2d

  implicit none

  call init_fruit
  call test_exp_cutoff
  call test_trig_cutoff
  call test_table2d_f_and_df
  call fruit_summary
  call fruit_finalize
endprogram
