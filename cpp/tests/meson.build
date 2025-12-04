# Tests for Atomistica C++

catch2_dep = dependency('catch2-with-main',
  fallback: ['catch2', 'catch2_with_main_dep'],
  required: true,
)

test_sources = files(
  'test_atomic_system.cpp',
  'test_neighbor_list.cpp',
  'test_lj.cpp',
  'test_spline.cpp',
  'test_cutoff_functions.cpp',
  'test_tersoff.cpp',
  'test_eam.cpp',
  'test_coulomb.cpp',
  'test_tightbinding.cpp',
  'test_integrators.cpp',
)

test_exe = executable('atomistica_cpp_tests',
  test_sources,
  include_directories: inc,
  dependencies: [atomistica_cpp_dep, catch2_dep],
)

test('atomistica_cpp_tests', test_exe)
