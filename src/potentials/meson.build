all_potentials_sources = []

potentials_sources = []
subdir('bop')
subdir('coulomb')
subdir('dispersion')
subdir('eam')
subdir('pair')


message('=== LISTCLASSES ===')

r = run_command(python, '../../tools/listclasses.py', all_potentials_sources, check: truew)
message(r.stdout())
message(r.stderr())

