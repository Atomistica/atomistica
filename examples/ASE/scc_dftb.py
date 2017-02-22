import os

from ase.structure import molecule

import atomistica.native as native
from atomistica import Atomistica

###

a = molecule('H2O')
a.center(vacuum=10.0)

###

database_folder = os.getenv('MIO')
if database_folder is None:
    raise RuntimeError('Please use environment variable MIO to specify path to mio Slater-Koster tables.')

tb = native.TightBinding(
    database_folder = database_folder,
    SolverLAPACK = dict(electronic_T=0.001),
    SCC = dict(dq_crit = 1e-4,
               mixing = 0.2, # 0.2
               andersen_memory = 3, # 3
               maximum_iterations = 100,
               log = True)
    )
calc = Atomistica(
    [ tb,
      native.DirectCoulomb(),
      native.SlaterCharges(cutoff=10.0) ],
    avgn = 1000
    )

a.set_calculator(calc)
print('potential energy =', a.get_potential_energy())

print('eigenvalues:')
for ev, occ in zip(tb.eigenvalues, tb.occupation):
    print(ev, occ)
