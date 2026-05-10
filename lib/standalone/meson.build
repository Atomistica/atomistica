standalone_deps = [atomistica_dep]
standalone_cpp_args = []

if netcdf_dep.found()
  standalone_deps += netcdf_dep
  standalone_cpp_args += '-DATOMISTICA_HAVE_NETCDF'
endif

executable('mdcore',
  ['mdcore.cpp', 'registry.cpp'],
  include_directories: inc,
  dependencies: standalone_deps,
  cpp_args: standalone_cpp_args,
  override_options: ['cpp_std=c++17'],
  install: false,
)
