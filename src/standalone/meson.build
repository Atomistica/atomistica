standalone_sources = [
	'signal_handler.f90',
	'r250.f',
	'rng.f90',
	'particles.f90',
	'domain_decomposition.f90',
	'neighbors.f90',
	'filter.f90',
	'molecules.f90',
	'dynamics.f90',
	'cyclic.f90',
	'symmetry.f90',
	'versioninfo.f90',
#	'atomistica.f90',
	'native_io.f90',
	'cfg.f90',
	'pdb.f90',
	'vtk.f90',
	'nc.f90',
	'xyz_f90.f90',
	'lammps_data.f90',
	'input_trajectory.f90',
	'lucy.f90',
	'square.f90',
	'interpolation_kernels_dispatch.f90'
]

incdirs = include_directories('..')

standalone = static_library(
    'standalone',
    standalone_sources,
    include_directories: incdirs,
    link_with: supplib,
    dependencies: atomistica_dependencies,
    install: true  # needed, otherwise it is not copied to wheel
)
