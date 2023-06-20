pair_sources = [
    'born_mayer.f90',
    'double_harmonic.f90',
    'harmonic.f90',
    'lj_cut.f90',
    'r6.f90',
]

incdirs = include_directories('../..', '../../support')

potentials_pair = static_library(
    'pair',
    pair_sources,
    include_directories: incdirs,
    link_with: [atomistica_supplib, atomistica_corelib],
    dependencies: atomistica_dependencies,
    install: true  # needed, otherwise it is not copied to wheel
)
