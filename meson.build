project('atomistica',
  'c', 'cpp', 'fortran',
  # Version will be determined by setuptools-scm via pyproject.toml
  version: '0.0.0',
  license: 'GPL-2.0-or-later',
  meson_version: '>= 1.1.0',
  default_options: [
    'buildtype=debugoptimized',
    'c_std=c99',
    'fortran_std=legacy',
  ],
)

py = import('python').find_installation(pure: false)
py_dep = py.dependency()

cc = meson.get_compiler('c')
cpp = meson.get_compiler('cpp')
fc = meson.get_compiler('fortran')

# NumPy dependency - use Python to find it rather than pkg-config
numpy_inc = run_command(py,
  ['-c', 'import numpy; print(numpy.get_include())'],
  check: true
).stdout().strip()

# NumPy include must be passed as compiler argument since it's an absolute path
# Create a dependency object that adds the include path as a compiler argument
numpy_dep = declare_dependency(
  compile_args: ['-I' + numpy_inc]
)

# LAPACK dependency
# Try pkg-config first, fall back to direct library search
lapack_dep = dependency('lapack', required: false)
if not lapack_dep.found()
  # Direct library search with explicit directories for Windows
  lapack_dirs = []
  if host_machine.system() == 'windows'
    # Add MSYS2 MinGW library directory if it exists
    msys2_lib = 'C:/msys64/mingw64/lib'
    if fc.has_argument('-L' + msys2_lib)
      lapack_dirs += [msys2_lib]
    endif
  endif

  # On Windows/MinGW, LAPACK is often provided as liblapack
  lapack_dep = fc.find_library('lapack', dirs: lapack_dirs, required: false)
  if not lapack_dep.found()
    # Try alternative names (libopenblas includes LAPACK on some systems)
    lapack_dep = fc.find_library('openblas', dirs: lapack_dirs, required: false)
  endif
  if not lapack_dep.found()
    error('LAPACK library not found. Please install LAPACK development files.')
  endif
endif

# Add -fPIC for all compiled code
add_project_arguments('-fPIC', language: ['c', 'cpp', 'fortran'])
add_project_arguments('-DHAVE_LAPACK', language: ['c', 'cpp'])
add_project_arguments('-DNO_BIND_C_OPTIONAL', language: ['c', 'cpp', 'fortran'])
add_project_arguments('-DPYTHON', language: ['c', 'cpp', 'fortran'])

# Enable Fortran preprocessing for files with preprocessor directives
add_project_arguments('-cpp', language: 'fortran')

# Allow free-form line length beyond 132 characters (GNU Fortran default)
if fc.get_id() == 'gcc'
  add_project_arguments('-ffree-line-length-none', language: 'fortran')
  # Allow argument mismatches for optional arguments (needed for GCC 10+)
  add_project_arguments('-fallow-argument-mismatch', language: 'fortran')
elif fc.get_id() == 'intel' or fc.get_id() == 'intel-cl'
  add_project_arguments('-extend-source', language: 'fortran')
endif

# Include directories
inc_dirs = include_directories(
  'src',
  'src/support',
  'src/potentials',
  'src/notb',
  'src/notb/dense',
)

# Run code generation scripts
message('Running code generation scripts...')

# Generate versioninfo
versioninfo_gen = custom_target('versioninfo',
  output: 'versioninfo.f90',
  command: ['sh', '@SOURCE_ROOT@/src/gen_versioninfo.sh',
            '@SOURCE_ROOT@/src', '@OUTDIR@', 'Python'],
  capture: false,
)

# Scan metadata and generate factories
# Note: We need to generate potential_sources.txt as well for dynamic source discovery
gen_factories = custom_target('gen_factories',
  output: ['coulomb_factory_f90.f90', 'coulomb_factory_c.c', 'coulomb_factory_c.h',
           'potentials_factory_f90.f90', 'potentials_factory_c.c',
           'potentials_factory_c.h', 'have.inc', 'potential_sources.txt'],
  input: files('tools/meta.py', 'tools/listclasses.py'),
  command: [py, '@SOURCE_ROOT@/build_helpers/generate_factories.py',
            '@SOURCE_ROOT@', '@OUTDIR@'],
  depends: [],
)

# Get the list of potential source files from the generated file
# This will be populated after the custom_target runs
# For now, we need to manually specify potential sources that we know exist
# Meson doesn't support reading files generated at build time for source lists

# First, compile coulomb modules separately to avoid parallel build race condition
coulomb_lib_sources = [
  'src/potentials/coulomb/coulomb_short_gamma.f90',
  'src/potentials/coulomb/cutoff_coulomb.f90',
  'src/potentials/coulomb/damp_short_gamma.f90',
  'src/potentials/coulomb/direct_coulomb.f90',
  'src/potentials/coulomb/gaussian_charges.f90',
  'src/potentials/coulomb/pme.f90',
  'src/potentials/coulomb/pme_kernel.f90',
  'src/potentials/coulomb/slater_charges.f90',
  'src/potentials/coulomb/fft_wrap.f',
  'src/potentials/coulomb/fft3-public.f',
]

# Coulomb modules depend on supplib and particles/neighbors modules, so need basic support compiled first
support_lib_sources = [
  versioninfo_gen,
  'src/support/c_f.f90',
  'src/support/error.f90',
  'src/support/System.f90',
  'src/support/MPI_context.f90',
  'src/support/Units.f90',
  'src/support/PeriodicTable.f90',
  'src/support/c_linearalgebra.cpp',
  'src/support/f_linearalgebra.f90',
  'src/support/f_ptrdict.f90',
  'src/support/c_ptrdict.c',
  'src/support/io.f90',
  'src/support/f_logging.f90',
  'src/support/c_logging.c',
  'src/support/timer.f90',
  'src/support/tls.f90',
  'src/support/misc.f90',
  'src/support/data.f90',
  'src/support/simple_spline.f90',
  'src/support/nonuniform_spline.f90',
  'src/support/cutoff.f90',
  'src/support/histogram1d.f90',
  'src/support/supplib.f90',
  'src/support/atomistica.f90',
  # Core particle and neighbor modules needed by potentials
  'src/python/f90/python_particles.f90',
  'src/python/f90/python_neighbors.f90',
  # Filter module needed by coulomb and potentials
  'src/core/filter.f90',
]

# Build support library first
support_lib = static_library('support',
  support_lib_sources,
  include_directories: [inc_dirs],
  dependencies: [lapack_dep, numpy_dep],
  fortran_args: ['-DHAVE_LAPACK'],
  install: false,
)

# Build coulomb library with dependency on support
coulomb_lib = static_library('coulomb',
  coulomb_lib_sources,
  include_directories: [inc_dirs],
  dependencies: [lapack_dep, numpy_dep],
  link_with: support_lib,
  fortran_args: ['-DHAVE_LAPACK'],
  install: false,
)

# Build the main atomisticalib library with static sources
# Note: This now depends on support_lib and coulomb_lib being built first
lib_sources = [
  'src/python/f90/particles_wrap.f90',
  'src/python/f90/neighbors_wrap.f90',
  'src/python/f90/python_helper.f90',
  'src/special/table2d.f90',
  'src/special/table3d.f90',
  'src/special/table4d.f90',
  'src/special/anderson_mixer.f90',
  'src/special/extrapolation.f90',
  gen_factories[0],  # coulomb_factory_f90.f90
  gen_factories[1],  # coulomb_factory_c.c
  gen_factories[3],  # potentials_factory_f90.f90
  gen_factories[4],  # potentials_factory_c.c
]

# Add potential source files
# These need to be added manually as Meson requires sources at configure time
lib_sources += files(
  # Pair potentials
  'src/potentials/pair_potentials/born_mayer.f90',
  'src/potentials/pair_potentials/double_harmonic.f90',
  'src/potentials/pair_potentials/harmonic.f90',
  'src/potentials/pair_potentials/lj_cut.f90',
  'src/potentials/pair_potentials/r6.f90',
  # EAM
  'src/potentials/eam/tabulated_alloy_eam.f90',
  'src/potentials/eam/tabulated_eam.f90',
  # BOP - Main module files only (other files are included via #include directives)
  'src/potentials/bop/tersoff/tersoff.f90',
  'src/potentials/bop/tersoff/tersoff_scr.f90',
  'src/potentials/bop/brenner/brenner.f90',
  'src/potentials/bop/brenner/brenner_scr.f90',
  'src/potentials/bop/kumagai/kumagai.f90',
  'src/potentials/bop/kumagai/kumagai_scr.f90',
  'src/potentials/bop/juslin/juslin.f90',
  'src/potentials/bop/juslin/juslin_scr.f90',
  # REBO2 needs its default tables module compiled first
  'src/potentials/bop/rebo2/rebo2_default_tables.f90',
  'src/potentials/bop/rebo2/rebo2.f90',
  'src/potentials/bop/rebo2/rebo2_scr.f90',
  # Dispersion
  'src/potentials/dispersion/dispdftd3.f90',
  # NOTB (Non-orthogonal Tight Binding)
  'src/notb/materials.f90',
  'src/notb/dense/dense_hamiltonian_type.f90',
  'src/notb/dense/c_dense_hamiltonian.cpp',
  'src/notb/dense/dense_hamiltonian.f90',
  'src/notb/dense/dense_hs.f90',
  'src/notb/dense/dense_forces.f90',
  'src/notb/dense/dense_repulsion.f90',
  'src/notb/dense/dense_scc.f90',
  'src/notb/dense/solver/dense_occupation.f90',
  'src/notb/dense/solver/dense_solver_lapack.f90',
  'src/notb/dense/solver/dense_solver_cp.f90',
  'src/notb/dense/solver/dense_solver_dispatch.f90',
  'src/notb/dense/analysis/dense_bonds.f90',
  'src/notb/dense/dense_notb.f90',
  # coulomb_dispatch depends on coulomb_lib being built - put it last
  'src/python/f90/coulomb_dispatch.f90',
)

# Create a build directory include for generated files
build_inc = include_directories('.')

# Create dependencies to ensure proper module file access
support_dep = declare_dependency(link_with: support_lib)
coulomb_dep = declare_dependency(link_with: coulomb_lib, dependencies: support_dep)

# Build the Python extension module directly with all sources
ext_sources = [
  'src/python/c/py_f.c',
  'src/python/c/particles.c',
  'src/python/c/neighbors.c',
  'src/python/c/coulomb.c',
  'src/python/c/coulomb_callback.c',
  'src/python/c/potential.c',
  'src/python/c/analysis.c',
  'src/python/c/atomisticamodule.c',
]

# Add all lib_sources to ext_sources for direct compilation
ext_sources += lib_sources

# Add include path for generated headers via compiler args
gen_headers_inc = declare_dependency(
  compile_args: ['-I' + meson.current_build_dir()],
  sources: [gen_factories[2], gen_factories[5], gen_factories[6]]  # coulomb_factory_c.h, potentials_factory_c.h, have.inc
)

py.extension_module('_atomistica',
  ext_sources,
  include_directories: [inc_dirs, build_inc],
  link_whole: [support_lib, coulomb_lib],
  dependencies: [py_dep, numpy_dep, lapack_dep, gen_headers_inc, coulomb_dep],
  fortran_args: ['-DHAVE_LAPACK'],
  install: true,
  subdir: 'atomistica',
)

# Install Python package
py.install_sources(
  [
    'src/python/atomistica/__init__.py',
    'src/python/atomistica/analysis.py',
    'src/python/atomistica/aseinterface.py',
    'src/python/atomistica/atomic_strain.py',
    'src/python/atomistica/deformation.py',
    'src/python/atomistica/hardware.py',
    'src/python/atomistica/io.py',
    'src/python/atomistica/join_calculators.py',
    'src/python/atomistica/logger.py',
    'src/python/atomistica/mdcore_io.py',
    'src/python/atomistica/native.py',
    'src/python/atomistica/parameters.py',
    'src/python/atomistica/snippets.py',
    'src/python/atomistica/tests.py',
  ],
  subdir: 'atomistica',
)

# Install scripts from tools directory
py.install_sources(
  [
    'src/python/tools/a_angle_distribution.py',
    'src/python/tools/a_convert.py',
    'src/python/tools/a_fire.py',
    'src/python/tools/a_g2.py',
    'src/python/tools/a_run.py',
    'src/python/tools/a_voro.py',
  ],
  subdir: 'atomistica' / 'tools',
)
