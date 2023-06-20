incdirs = include_directories('../..', '../../support')

executable(
    'main',
    'main.f90',
    include_directories: incdirs,
    link_with: [atomistica_supplib, atomistica_corelib, potentials_dispersion, potentials_eam],
    dependencies: atomistica_dependencies
)
