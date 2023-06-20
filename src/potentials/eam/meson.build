eam_sources = [
    'tabulated_eam.f90',
    'tabulated_alloy_eam.f90'
]

incdirs = include_directories('../..', '../../support')

eam = static_library(
    'eam',
    eam_sources,
    include_directories: incdirs,
    link_with: supplib,
    dependencies: atomistica_dependencies,
    install: true  # needed, otherwise it is not copied to wheel
)
