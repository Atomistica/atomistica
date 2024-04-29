all_potentials_sources = []
potentials_sources = []
path = ''
subdir('bop')
subdir('coulomb')
subdir('dispersion')
subdir('eam')
subdir('pair')


message('=== LISTCLASSES ===')

r = run_command(python, '../../tools/discover_classes.py', all_potentials_sources, check: true)
message(r.stdout())
message(r.stderr())

