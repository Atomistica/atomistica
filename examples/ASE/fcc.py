#! /usr/bin/env python

import numpy as np

from ase.constraints import UnitCellFilter
from ase.lattice.cubic import FaceCenteredCubic
from ase.optimize import FIRE
from ase.utils.eos import EquationOfState

from atomistica import TabulatedAlloyEAM

###

n = 2
a = FaceCenteredCubic('Au', size=[n, n, n])
x0 = a.cell[0, 0]/n
c = TabulatedAlloyEAM(fn='Au-Grochola-JCP05.eam.alloy')
a.set_calculator(c)

# Vary volume and fit minimum
def en(a, x):
	a.set_cell([x, x, x], scale_atoms=True)
	return a.get_potential_energy()
x = np.linspace(0.9*x0, 1.1*x0, 101)
e = [ en(a, n*_x)/(n**3) for _x in x ]
eos = EquationOfState(x**3, e)
v0, e0, B = eos.fit()

print 'lattice constant (from equation of state) =', v0**(1./3)

# Keep cell rectangular during optimization
FIRE(UnitCellFilter(a, mask=[1,1,1,0,0,0]), logfile=None).run(fmax=0.0001)

print 'lattice constant (from stress equilibration) =', a.cell[0, 0]/n