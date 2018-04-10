import os

import numpy as np
np.set_printoptions(precision=2, linewidth=120, suppress=True)

from ase.build import molecule

import atomistica.native as native
from atomistica import Atomistica

###

molstr = 'CO2'
a = molecule(molstr)
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
print('{}, potential energy: {} eV'.format(molstr, a.get_potential_energy()))

### Extract electronic structure and print some properties

potential = calc.get_electrostatic_potential()
charges = calc.get_charges()
charges0 = tb.neutral_charges
norb = tb.number_of_orbitals

print()
print('Atom, i  Number of orbitals  Charge, q  q of neutral atom, q0  Potential, phi')
print('=======  ==================  =========  =====================  ==============')
for i, (n, q, q0, phi) in enumerate(zip(norb, charges, charges0, potential)):
    print('{:7} {:19} {:10.2f} {:22.2f} {:15.2f}'.format(i, n, q, q0, phi))
print('                                sum(q)                           sum(q*phi)/2')
print('                            {:10.2f}                        {:15.2f} eV'.format(charges.sum(), (charges*potential).sum()/2))

# Get electronic structure information from the tight-binding calculator
f = tb.occupation
dm = tb.density_matrix
evecs = tb.eigenvectors
evals = tb.eigenvalues
H0 = tb.Hamiltonian_matrix
S = tb.overlap_matrix

print()
print('Eigenvalue, e  Occupation, f')
print('=============  =============')
for ev, occ in zip(evals, f):
    print('{:13.2f} {:14.2f}'.format(ev, occ))
print('     sum(f*e)         sum(f)')
print('{:13.2f} eV {:11.2f}'.format((f*evals).sum(), f.sum()))

print()
print('Hamiltonian matrix (*without* electrostatic contribution), H0')
print('=============================================================')
print(H0)
print()
print('Overlap matrix, S')
print('=================')
print(S)
print()
print('Density matrix, rho')
print('===================')
print(dm)
print()

# Construct effective Hamiltonian
phiarr = np.ravel(list([_phi]*_n for _phi, _n in zip(potential, norb)))
phimat = phiarr.reshape(1, -1) + phiarr.reshape(-1, 1)
Heff = H0 - S*phimat/2

print('Derived quantities')
print('==================')
print('sum(f)                                    = {:10.2f}'.format(f.sum()))
print('tr(S.rho)                                 = {:10.2f}      [should equal sum(f)]'.format(np.trace(S.dot(dm))))
print('tr(H0.rho)                                = {:10.2f} eV'.format(np.trace(H0.dot(dm))))
print('tr(H0.rho) + sum(q*phi)/2                 = {:10.2f} eV   [this is the potential energy]'.format(np.trace(H0.dot(dm)) + (charges*potential).sum()/2))
print('sum(f*e)                                  = {:10.2f}'.format((f*evals).sum()))
print('tr(Heff.rho)                              = {:10.2f} eV   [should equal sum(f*e)]'.format(np.trace(Heff.dot(dm))))
print('tr(Heff.rho) - sum(q*phi)/2 - sum(q0*phi) = {:10.2f} eV.  [should equal tr(H0.rho) + sum(q*phi)/2]' \
  .format(np.trace(Heff.dot(dm)) - (charges*potential).sum()/2 - (charges0*potential).sum()))
