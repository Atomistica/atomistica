pybind11_dep = dependency('pybind11',
  fallback: ['pybind11', 'pybind11_dep'],
  required: true,
)

py = import('python').find_installation(pure: false)

py.extension_module('_atomistica_cpp',
  'bindings.cpp',
  include_directories: inc,
  dependencies: [atomistica_cpp_dep, pybind11_dep],
  install: true,
  subdir: 'atomistica_cpp',
)

# Install Python package files
py.install_sources(
  '__init__.py',
  'ase_calculator.py',
  'parameters.py',
  subdir: 'atomistica_cpp',
)
