old_path = path
path = join_paths(path, 'coulomb')

coulomb_sources = [
    join_paths(path, 'coulomb_short_gamma.f90'),
    join_paths(path, 'cutoff_coulomb.f90'),
    join_paths(path, 'damp_short_gamma.f90'),
    join_paths(path, 'direct_coulomb.f90'),
    join_paths(path, 'gaussian_charges.f90'),
    join_paths(path, 'pme_kernel.f90'),
    join_paths(path, 'pme.f90'),
    join_paths(path, 'slater_charges.f90')
]

if not atomistica_with_fftw3
    coulomb_sources += [
        join_paths(path, 'fft3-public.f'),
        join_paths(path, 'fft_wrap.f')
    ]
endif

all_potentials_sources += coulomb_sources
