rebo2_sources = [
    'rebo2_default_tables.f90',
    'rebo2.f90',
    'rebo2_scr.f90'
]

incdirs = include_directories('../../..', '../../../support')

potentials_rebo2 = static_library(
    'rebo2',
    rebo2_sources,
    include_directories: incdirs,
    link_with: [atomistica_supplib, atomistica_corelib],
    dependencies: atomistica_dependencies,
    install: true  # needed, otherwise it is not copied to wheel
)

all_potentials_sources += rebo2_sources
