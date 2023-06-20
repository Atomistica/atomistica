special_sources = [
    'extrapolation.f90'
]

incdirs = include_directories('..', '../support')

atomistica_special = static_library(
    'special',
    special_sources,
    include_directories: incdirs,
    link_with: [atomistica_supplib, atomistica_corelib],
    dependencies: atomistica_dependencies,
    install: true  # needed, otherwise it is not copied to wheel
)
