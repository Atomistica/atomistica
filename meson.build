project(
    'Atomistica', # Project name
    'c', 'cpp', 'fortran', # Project type. We need a C, C++ and Fortran compiler.
    default_options: ['cpp_std=c++11'],
    version: run_command('python3', '-c', '"import versioneer; print(versioneer.get_version())"', check: true).stdout().strip()
)

pymod = import('python')
python = pymod.find_installation('python3',
    required: true,
)

git = find_program('git', required: false)

atomistica_target = 'standalone'

#add_global_arguments('-ffree-form', language: 'fortran')

# Store global configuration to dictionary
atomistica_date = run_command(python, '-c', '"import versioneer; print(versioneer.get_versions()[\'date\'])"').stdout().strip()
if git.found()
    atomistica_url = run_command(git, 'config', '--get', 'remote.origin.url').stdout().strip()
else
    atomistica_url = 'N/A'
endif

fortran_compiler = meson.get_compiler('fortran')

if fortran_compiler.get_id() == 'gcc'
    add_project_arguments('-fallow-invalid-boz', language: 'fortran')  # Only standalone coe, for r250.f
    add_project_arguments('-ffree-line-length-none', language: 'fortran')
    add_project_arguments('-x', 'f95-cpp-input', language: 'fortran')
    add_project_arguments('-DNO_BIND_C_OPTIONAL', language: 'fortran')  # Only gfortran
    # Disable some warnings
    add_project_arguments('-Wno-unused-dummy-argument', language: 'fortran')
endif

atomistica_configuration = configuration_data({
    'atomistica_veersion': meson.project_version(),
    'atomistica_data': atomistica_date,
    'atomistica_url': atomistica_url,
    'hostname': 'meson',
    'arch': target_machine.cpu_family(),
    'compileroptions': fortran_compiler.get_id(),
    'compilerversion': fortran_compiler.version()
})

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
