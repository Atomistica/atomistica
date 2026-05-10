project('atomistica-cpp', 'cpp',
  version: '0.1.0',
  license: 'GPL-2.0-or-later',
  meson_version: '>= 1.1.0',
  default_options: [
    'cpp_std=c++17',
    'warning_level=2',
    'buildtype=debugoptimized',
    'b_ndebug=if-release',
  ],
)

# Dependencies from wrapdb
eigen_dep = dependency('eigen3',
  fallback: ['eigen', 'eigen_dep'],
  required: true,
)

# Optional OpenMP
openmp_dep = dependency('openmp', required: get_option('enable_openmp'))

# Optional LAPACK (for tight-binding solver)
lapack_dep = dependency('lapack', required: false)

# Optional s-dftd3 (for DFT-D3 dispersion)
sdftd3_dep = dependency('s-dftd3', required: false)

# Optional NetCDF (for trajectory output in standalone mdcore)
netcdf_dep = dependency('netcdf', required: false)

# Collect all dependencies
atomistica_deps = [eigen_dep]
if openmp_dep.found()
  atomistica_deps += openmp_dep
endif
if lapack_dep.found()
  atomistica_deps += lapack_dep
endif
if sdftd3_dep.found()
  atomistica_deps += sdftd3_dep
  add_project_arguments('-DHAVE_SDFTD3', language: 'cpp')
  message('s-dftd3 found: DFT-D3 dispersion will be available')
else
  message('s-dftd3 not found: DFT-D3 dispersion will NOT be available (install s-dftd3 to enable)')
endif

if netcdf_dep.found()
  message('NetCDF found: NetCDF trajectory output will be available in standalone mdcore')
else
  message('NetCDF not found: NetCDF trajectory output will NOT be available (install netcdf-c to enable)')
endif

# Include directories
inc = include_directories('include')

# Core library sources
core_sources = files(
  'src/core/atomic_system.cpp',
  'src/core/neighbor_list.cpp',
  'src/math/spline.cpp',
  'src/math/cutoff_functions.cpp',
  'src/potentials/pair/lj.cpp',
  'src/tightbinding/anderson_mixer.cpp',
  'src/tightbinding/bond_analysis.cpp',
  'src/tightbinding/dftb.cpp',
  'src/tightbinding/hamiltonian.cpp',
  'src/tightbinding/materials.cpp',
  'src/tightbinding/solver.cpp',
)

# Build the core library
atomistica_lib = library('atomistica',
  core_sources,
  include_directories: inc,
  dependencies: atomistica_deps,
  install: true,
)

# Declare dependency for consumers
atomistica_dep = declare_dependency(
  link_with: atomistica_lib,
  include_directories: inc,
  dependencies: atomistica_deps,
)

# Python bindings (optional)
if get_option('enable_python')
  subdir('python')
endif

# Tests
if get_option('enable_tests')
  subdir('tests')
endif

# Standalone mdcore executable (optional)
if get_option('enable_standalone')
  subdir('standalone')
endif

# Install headers
install_subdir('include/atomistica',
  install_dir: get_option('includedir'),
)
