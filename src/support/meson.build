supplib_sources = [
    'c_f.f90',
    'error.f90',
    'System.f90',
    'MPI_context.f90',
    'Units.f90',
    'c_linearalgebra.cpp',
    'f_linearalgebra.f90',
    'PeriodicTable.f90',
    'f_ptrdict.f90',
    'c_ptrdict.c',
    'io.f90',
    'f_logging.f90',
    'c_logging.c',
    'misc.f90',
    'data.f90',
    'timer.f90',
    'tls.f90',
    'simple_spline.f90',
    'nonuniform_spline.f90',
    'cutoff.f90',
    'histogram1d.f90',
    'supplib.f90'
]

incdirs = include_directories('..')

supplib = static_library(
    'supplib',
    supplib_sources,
    include_directories: incdirs,
    dependencies: atomistica_dependencies,
    install: true  # needed, otherwise it is not copied to wheel
)
