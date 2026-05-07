project('atomistica',
  'cpp',
  # Version determined by setuptools-scm via pyproject.toml
  version : run_command('python3', 'discover_version.py', check : true).stdout().strip(),
  license: 'GPL-2.0-or-later',
  meson_version: '>= 1.1.0',
  default_options: [
    'buildtype=debugoptimized',
  ],
)

py = import('python').find_installation(pure: false)
py_dep = py.dependency()

# Optional LAPACK (used by tight-binding LAPACK eigensolver)
lapack_dep = dependency('lapack', required: false)
if not lapack_dep.found()
  lapack_dep = dependency('openblas', required: false)
endif

# ===========================================================================
# atomistica: Modern C++17 implementation of interatomic potentials
# ===========================================================================

eigen_dep = dependency('eigen3',
  fallback: 'eigen',
  required: false,
)

pybind11_dep = dependency('pybind11',
  fallback: ['pybind11', 'pybind11_dep'],
  required: false,
)

# Optional s-dftd3 for DFT-D3 dispersion
sdftd3_dep = dependency('s-dftd3', required: false)

if eigen_dep.found() and pybind11_dep.found()

  atomistica_deps = [eigen_dep]
  if lapack_dep.found()
    atomistica_deps += lapack_dep
  endif
  if sdftd3_dep.found()
    atomistica_deps += sdftd3_dep
    add_project_arguments('-DHAVE_SDFTD3', language: 'cpp')
    message('s-dftd3 found: DFT-D3 dispersion will be available')
  else
    message('s-dftd3 not found: DFT-D3 dispersion will NOT be available')
  endif

  cpp_inc = include_directories('lib/include')

  cpp_lib_sources = files(
    'lib/src/core/atomic_system.cpp',
    'lib/src/core/neighbor_list.cpp',
    'lib/src/math/spline.cpp',
    'lib/src/math/cutoff_functions.cpp',
    'lib/src/potentials/pair/lj.cpp',
    'lib/src/tightbinding/anderson_mixer.cpp',
    'lib/src/tightbinding/bond_analysis.cpp',
    'lib/src/tightbinding/dftb.cpp',
    'lib/src/tightbinding/hamiltonian.cpp',
    'lib/src/tightbinding/materials.cpp',
    'lib/src/tightbinding/solver.cpp',
  )

  atomistica_lib = static_library('atomistica',
    cpp_lib_sources,
    include_directories: cpp_inc,
    dependencies: atomistica_deps,
    override_options: ['cpp_std=c++17'],
    install: false,
  )

  atomistica_dep = declare_dependency(
    link_with: atomistica_lib,
    include_directories: cpp_inc,
    dependencies: atomistica_deps,
  )

  # Python extension: _atomistica (installed as atomistica package)
  py.extension_module('_atomistica',
    'lib/python/bindings.cpp',
    include_directories: cpp_inc,
    dependencies: [py_dep, pybind11_dep, atomistica_dep],
    override_options: ['cpp_std=c++17'],
    install: true,
    subdir: 'atomistica',
  )

  py.install_sources(
    [
      'lib/python/__init__.py',
      'lib/python/ase_calculator.py',
      'lib/python/parameters.py',
    ],
    subdir: 'atomistica',
  )

  message('atomistica Python extension will be built')

else
  message('Eigen3 or pybind11 not found; skipping atomistica')
endif
