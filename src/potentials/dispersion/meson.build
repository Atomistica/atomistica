dispersion_sources = [
    'dispdftd3.f90'
]

incdirs = include_directories('../..', '../../support')

dispersion = static_library(
    'dispersion',
    dispersion_sources,
    include_directories: incdirs,
    link_with: supplib,
    dependencies: atomistica_dependencies,
    install: true  # needed, otherwise it is not copied to wheel
)
