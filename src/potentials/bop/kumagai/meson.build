kumagai_sources = [
    'kumagai.f90',
    'kumagai_scr.f90'
]

incdirs = include_directories('../../..', '../../../support')

potentials_kumagai = static_library(
    'kumagai',
    kumagai_sources,
    include_directories: incdirs,
    link_with: [atomistica_supplib, atomistica_corelib],
    dependencies: atomistica_dependencies,
    install: true  # needed, otherwise it is not copied to wheel
)

all_potentials_sources += kumagai_sources
