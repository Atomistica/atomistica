juslin_sources = [
    'juslin.f90',
    'juslin_scr.f90'
]

incdirs = include_directories('../../..', '../../../support')

potentials_juslin = static_library(
    'juslin',
    juslin_sources,
    include_directories: incdirs,
    link_with: [atomistica_supplib, atomistica_corelib],
    dependencies: atomistica_dependencies,
    install: true  # needed, otherwise it is not copied to wheel
)

all_potentials_sources += juslin_sources
