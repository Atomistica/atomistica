coulomb_sources = [
    'coulomb_short_gamma.f90',
    'cutoff_coulomb.f90',
    'damp_short_gamma.f90',
    'direct_coulomb.f90',
    'gaussian_charges.f90',
    'pme_kernel.f90',
    'pme.f90',
    'slater_charges.f90'
]

if not atomistica_with_fftw3
    coulomb_sources += [
        'fft3-public.f',
        'fft_wrap.f',
    ]
endif

incdirs = include_directories('../..', '../../support')

potentials_coulomb = static_library(
    'coulomb',
    pair_sources,
    include_directories: incdirs,
    link_with: [atomistica_supplib, atomistica_corelib],
    dependencies: atomistica_dependencies,
    install: true  # needed, otherwise it is not copied to wheel
)
