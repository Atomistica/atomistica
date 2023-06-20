brenner_sources = [
    'brenner.f90',
    'brenner_scr.f90'
]

incdirs = include_directories('../../..', '../../../support')

potentials_brenner = static_library(
    'brenner',
    brenner_sources,
    include_directories: incdirs,
    link_with: [atomistica_supplib, atomistica_corelib],
    dependencies: atomistica_dependencies,
    install: true  # needed, otherwise it is not copied to wheel
)

all_potentials_sources += brenner_sources