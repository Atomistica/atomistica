project(
    'Atomistica', # Project name
    'c', 'cpp', 'fortran', # Project type. We need a C, C++ and Fortran compiler.
    default_options: ['cpp_std=c++11'],
    version: run_command('python3', 'discover_version.py', check: true).stdout().strip()
)

add_global_arguments('-ffree-form', language: 'fortran')
add_global_arguments('-ffree-line-length-none', language: 'fortran')
add_global_arguments('-x', 'f95-cpp-input', language: 'fortran')

add_global_arguments('-DNO_BIND_C_OPTIONAL', language: 'fortran')  # Only gfortran

pymod = import('python')
python = pymod.find_installation('python3',
    required: true,
)

version = run_command(python, 'discover_version.py', '--full', check: true).stdout().strip().split()
version_dirty = version[0]
version_str = version[1]
version_hash = version[2]

fftw3 = dependency('fftw3', required: false)
if fftw3.found()
    atomistica_with_fftw3 = true
    message('Atomistica FFTW: *** YES ***')
    add_global_arguments('-DWITH_FFTW', language : ['fortran'])
else
    atomistica_with_fftw3 = false
    message('Atomistica FFTW: no')
endif

netcdf = dependency('netcdf', required: false)
if netcdf.found()
    atomistica_with_netcdf = true
    message('Atomistica netCDF: *** YES ***')
    add_global_arguments('-DWITH_NETCDF', language : ['fortran'])
else
    atomistica_with_netcdf = false
    message('Atomistica netCDF: no')
endif


atomistica_dependencies = [fftw3, netcdf]

subdir('src')
