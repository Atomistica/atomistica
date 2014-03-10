program run_tests
  use fruit
  use test_cutoff

  implicit none

  call init_fruit
  call test_exp_cutoff
  call test_trig_cutoff
  call fruit_summary
  call fruit_finalize
endprogram
