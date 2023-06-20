tersoff_sources = [
    'tersoff.f90',
    'tersoff_scr.f90'
]

incdirs = include_directories('../../..', '../../../support')

potentials_tersoff = static_library(
    'tersoff',
    tersoff_sources,
    include_directories: incdirs,
    link_with: [atomistica_supplib, atomistica_corelib],
    dependencies: atomistica_dependencies,
    install: true  # needed, otherwise it is not copied to wheel
)

all_potentials_sources += tersoff_sources
