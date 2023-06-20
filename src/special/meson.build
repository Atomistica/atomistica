special_sources = [
    'anderson_mixer.f90',
    'extrapolation.f90',
    'table2d.f90',
    'table3d.f90',
    'table4d.f90'
]

incdirs = include_directories('..')

special = static_library(
    'special',
    special_sources,
    include_directories: incdirs,
    link_with: supplib,
    dependencies: atomistica_dependencies,
    install: true  # needed, otherwise it is not copied to wheel
)
