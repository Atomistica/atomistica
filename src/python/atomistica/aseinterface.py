"""ASE interface to Atomistica.
"""

import copy
from math import sqrt, log

import _atomistica

import numpy as np

from ase.data import atomic_numbers


class Atomistica:
    """Atomistica ASE calculator.
    """

    CELL_TOL = 1e-16
    POSITIONS_TOL = 1e-16

    def __init__(self, potential, avgn, **kwargs):
        """
        Initialize a potential. *potential* is the a native Atomistica
        potential class, in which case the arguments are given by the keywords
        of this function. Alternatively, it can be a list of tuples
        [ ( pot1, args1 ), ... ] in which case *pot1* is the name of the first
        potential and *args1* a dictionary containing its arguments.
        """
        self.avgn                   = avgn

        self.particles              = None
        self.nl                     = None
        self.potential              = potential(**kwargs)

        self.force_update           = True

        self.kwargs                 = kwargs

        self.lees_edwards_dx        = None
        self.leed_edwards_dv        = None

        self.compute_epot_per_at    = False
        self.compute_epot_per_bond  = False
        self.compute_f_per_bond     = False
        self.compute_wpot_per_at    = False
        self.compute_wpot_per_bond  = False

        self.epot_per_at            = None
        self.epot_per_bond          = None
        self.f_per_bond             = None
        self.wpot_per_at            = None
        self.wpot_per_bond          = None


    def initialize(self, atoms):
        # For now, length and elements are fixed
        # FIXME! Override pbc since code exits when atom moves out of the box!
        pbc = np.array( [ True, True, True ] )
        self.particles  = _atomistica.Particles()
        self.potential.register_data(self.particles)
        self.particles.allocate(len(atoms))
        self.particles.set_cell(atoms.get_cell(), pbc)

        if self.lees_edwards_dx is not None:
            self.particles.set_lees_edwards(self.lees_edwards_dx, self.lees_edwards_dv)

        Z  = self.particles.Z

        for i, at in enumerate(atoms):
            Z[i]   = atomic_numbers[at.symbol]

        self.particles.coordinates[:, :]  = atoms.get_positions()[:, :]
        # Notify the Particles object of a change
        self.particles.I_changed_positions()

        self.particles.update_elements()

        # Initialize and set neighbor list
        self.nl = _atomistica.Neighbors(self.avgn)

        # Tell the potential about the new Particles and Neighbors object
        self.potential.bind_to(self.particles, self.nl)

        # Force re-computation of energies/forces the next time these are requested
        self.force_update = True


    def set_per_at(self, epot=None, f=None, wpot=None):
        if epot is not None:
            self.compute_epot_per_at  = epot
        if wpot is not None:
            self.compute_wpot_per_at  = wpot
        self.force_update = True

    def set_per_bond(self, epot=None, f=None, wpot=None):
        if epot is not None:
            self.compute_epot_per_bond  = epot
        if f is not None:
            self.compute_f_per_bond     = f
        if wpot is not None:
            self.compute_wpot_per_bond  = wpot
        self.force_update = True


    def update(self, atoms, force_update=False):
        # Number of particles changed? -> Reinitialize potential
        if self.particles is None or len(self.particles.Z) != len(atoms):
            self.initialize(atoms)
        # Type of particles changed? -> Reinitialize potential
        elif np.any(self.particles.Z != atoms.get_atomic_numbers()):
            self.initialize(atoms)

        # Cell changed? FIXME! Add PBC,LEBC changed
        cell_chgd  = False
        cell       = self.particles.cell
        pbc        = atoms.get_pbc()
        #if np.any(np.abs(cell - atoms.get_cell()) > self.CELL_TOL):
        if np.any(cell != atoms.get_cell()):
            cell[:, :]  = atoms.get_cell()[:, :]
            self.particles.set_cell(cell, pbc)
            cell_chgd  = True

        pos_chgd   = False
        positions  = self.particles.coordinates
        scaled     = np.linalg.solve(cell.T, positions.T).T
        for i in range(3):
            if pbc[i]:
                # Yes, we need to do it twice.
                # See the scaled_positions.py test
                scaled[:, i] %= 1.0
                scaled[:, i] %= 1.0
        #if np.any(np.abs(scaled - atoms.get_scaled_positions()) > self.POSITIONS_TOL):
        if np.any(scaled != atoms.get_scaled_positions()):
            positions[:, :]  = atoms.get_positions()[:, :]
            pos_chgd         = True
            # Notify the Particles object of a change
            self.particles.I_changed_positions()

        if pos_chgd or cell_chgd or self.force_update or force_update:
            self.calculate()
            self.force_update  = False


    def get_potential_energy(self, atoms, force_consistent=False):
        self.update(atoms)
        return self.epot


    def get_forces(self, atoms):
        self.update(atoms)

        return self.forces.copy()


    def get_stress(self, atoms):
        self.update(atoms)

        st = self.wpot/atoms.get_volume()

        return np.array( [ st[0,0], st[1,1], st[2,2], (st[1,2]+st[2,1])/2, (st[0,2]+st[2,0])/2, (st[0,1]+st[1,0])/2 ] )


    def calculate(self):
        self.epot    = 0.0
        self.forces  = None
        self.wpot    = 0.0
        e, f, w, self.epot_per_at, self.epot_per_bond, self.f_per_bond, \
            self.wpot_per_at, self.wpot_per_bond  = \
            self.potential.energy_and_forces(
            self.particles, self.nl,
            epot_per_at    = self.compute_epot_per_at,
            epot_per_bond  = self.compute_epot_per_bond,
            f_per_bond     = self.compute_f_per_bond,
            wpot_per_at    = self.compute_wpot_per_at,
            wpot_per_bond  = self.compute_wpot_per_bond)

        self.epot += e
        if self.forces is None:
            self.forces  = f
        else:
            self.forces += f
        self.wpot += w


    ### Atomistica features
    def set_lees_edwards(self, dx, dv=None):
        self.lees_edwards_dx  = dx
        if dv is None:
            self.lees_edwards_dv  = [ 0.0, 0.0, 0.0 ]
        else:
            self.lees_edwards_dv  = dv

        if self.particles is not None and self.lees_edwards_dx is not None:
            self.particles.set_lees_edwards(self.lees_edwards_dx, self.lees_edwards_dv)
            self.force_update = True



    def get_atomic_stress(self):
        r = np.zeros( [ 3, 3, len(self.particles) ] )
        for i, a in enumerate(self.particles):
            r[:, :, i] = a.w
        return r


    def get_neighbors(self):
        return self.nl.get_neighbors(self.particles)


    def __str__(self):
        return self.potential.__str__()

    

### Short-cuts

# Simple r^6 pair potential

if hasattr(_atomistica, 'r6'):
    class r6(Atomistica):
        """r^6 potential.
        """

        def __init__(self, **kwargs):
            apply(Atomistica.__init__, (self, _atomistica.r6, 100,), kwargs)


# Tersoff potential

# Parameters
Tersoff_PRB_39_5566_Si_C = {
  "__ref__":  "Tersoff J., Phys. Rev. B 39, 5566 (1989)",
  "el":       [  "C",   "Si"  ],
  "A":        [  1.3936e3,    sqrt(1.3936e3*1.8308e3),  1.8308e3  ],
  "B":        [  3.4674e2,    sqrt(3.4674e2*4.7118e2),  4.7118e2  ],
  "xi":       [  1.0,         0.9776e0,                 1.0       ],
  "lambda":   [  3.4879e0,    (3.4879e0+2.4799e0)/2,    2.4799e0  ],
  "mu":       [  2.2119e0,    (2.2119e0+1.7322e0)/2,    1.7322e0  ],
  "mubo":     [  0.0,         0.0,                      0.0       ],
  "m":        [  1,           1,                        1         ],
  "beta":     [  1.5724e-7,   1.1000e-6   ],
  "n":        [  7.2751e-1,   7.8734e-1   ],
  "c":        [  3.8049e4,    1.0039e5    ],
  "d":        [  4.3484e0,    1.6217e1    ],
  "h":        [  -5.7058e-1,  -5.9825e-1  ],
  "r1":       [  1.80,        sqrt(1.80*2.70),          2.70      ],
  "r2":       [  2.10,        sqrt(2.10*3.00),          3.00      ],
  }

Tersoff_PRB_39_5566_Si_C__Scr = {
  "__ref__":  "Tersoff J., Phys. Rev. B 39, 5566 (1989)",
  "el":       [  "C",   "Si"  ],
  "A":        [  1.3936e3,    sqrt(1.3936e3*1.8308e3),  1.8308e3  ],
  "B":        [  3.4674e2,    sqrt(3.4674e2*4.7118e2),  4.7118e2  ],
  "xi":       [  1.0,         0.9776e0,                 1.0       ],
  "lambda":   [  3.4879e0,    (3.4879e0+2.4799e0)/2,    2.4799e0  ],
  "mu":       [  2.2119e0,    (2.2119e0+1.7322e0)/2,    1.7322e0  ],
  "mubo":     [  0,           0,                        0         ],
  "m":        [  3,           3,                        3         ],
  "beta":     [  1.5724e-7,   1.1000e-6   ],
  "n":        [  7.2751e-1,   7.8734e-1   ],
  "c":        [  3.8049e4,    1.0039e5    ],
  "d":        [  4.3484e0,    1.6217e1    ],
  "h":        [  -5.7058e-1,  -5.9825e-1  ],
  "r1":       [  2.00,        sqrt(2.00*2.50),          2.50      ],
  "r2":       [  2.00*1.2,    sqrt(2.00*2.50)*1.2,      2.50*1.2  ],
  "or1":      [  2.00,        sqrt(2.00*3.00),          3.00      ],
  "or2":      [  2.00*2.0,    sqrt(2.00*3.00)*2.0,      3.00*2.0  ],
  "bor1":     [  2.00,        sqrt(2.00*3.00),          3.00      ],
  "bor2":     [  2.00*2.0,    sqrt(2.00*3.00)*2.0,      3.00*2.0  ],
  "Cmin":     [  1.00,        1.00,                     1.00      ],
  "Cmax":     [  3.00,        3.00,                     3.00      ],
  }
# mubo is 1/dimer length
p = Tersoff_PRB_39_5566_Si_C__Scr
for i in range(3):
    p['mubo'][i] = (p['lambda'][i]-p['mu'][i])/ \
        log((p['lambda'][i]*p['A'][i])/(p['mu'][i]*p['B'][i]))

Goumri_Said_ChemPhys_302_135_Al_N = {
  "__ref__":  "Goumri-Said S., Kanoun M.B., Merad A.E., Merad G., Aourag H., Chem. Phys. 302, 135 (2004)",
  "el":       [  "Al",   "N"  ],
  "r1":       [  3.20,       2.185,     1.60     ],
  "r2":       [  3.60,       2.485,     2.00     ],
  "A":        [  746.698,    3000.214,  636.814  ],
  "B":        [  40.451,     298.81,    511.76   ],
  "xi":       [  1.0,        1.0,       1.0      ],
  "lambda":   [  2.4647,     3.53051,   5.43673  ],
  "mu":       [  0.9683,     1.99995,   2.7      ],
  "beta":     [  1.094932,   5.2938e-3  ],
  "n":        [  6.085605,   1.33041    ],
  "c":        [  0.074836,   2.0312e4   ],
  "d":        [  19.569127,  20.312     ],
  "h":        [  -0.659266,  -0.56239   ]
  }


if hasattr(_atomistica, 'Tersoff'):
    class Tersoff(Atomistica):
        """Tersoff potential.
        Refs: J. Tersoff, Phys. Rev. B 39, 5566 (1989)
        """

        def __init__(self, **kwargs):
            apply(Atomistica.__init__, (self, _atomistica.Tersoff, 100,), kwargs)


if hasattr(_atomistica, 'TersoffScr'):
    class TersoffScr(Atomistica):
        """Tersoff potential including screening functions.
        Refs: J. Tersoff, Phys. Rev. B 39, 5566 (1989)
        """

        def __init__(self, **kwargs):
            apply(Atomistica.__init__, (self, _atomistica.TersoffScr, 1000,), kwargs)


# Karsten Albe's BOP

# Parameters
Erhart_PRB_71_035211_SiC = {
    "__ref__":  "Erhart P., Albe K., Phys. Rev. B 71, 035211 (2005)",
    "el":       [  "C",   "Si" ],
    "D0":       [  6.00,      4.36,       3.24      ],
    "r0":       [  1.4276,    1.79,       2.232     ],
    "S":        [  2.167,     1.847,      1.842     ],
    "beta":     [  2.0099,    1.6991,     1.4761    ],
    "gamma":    [  0.11233,   0.011877,   0.114354  ],
    "c":        [  181.910,   273987.0,   2.00494   ],
    "d":        [  6.28433,   180.314,    0.81472   ],
    "h":        [  0.5556,    0.68,       0.259     ],
    "mu":       [  0.0,       0.0,        0.0       ],
    "n":        [  1.0,       1.0,        1.0       ],
    "m":        [  1,         1,          1         ],
    "r1":       [  1.85,      2.20,       2.68      ],
    "r2":       [  2.15,      2.60,       2.96      ]
    }


Erhart_PRB_71_035211_SiC__Scr = {
    "__ref__":  "Erhart P., Albe K., Phys. Rev. B 71, 035211 (2005)",
    "el":       [  "C",   "Si" ],
    "D0":       [  6.00,        4.36,      3.24       ],
    "r0":       [  1.4276,      1.79,      2.232      ],
    "S":        [  2.167,       1.847,     1.842      ],
    "beta":     [  2.0099,      1.6991,    1.4761     ],
    "gamma":    [  0.11233,     0.011877,  0.114354   ],
    "c":        [  181.910,     273987.0,  2.00494    ],
    "d":        [  6.28433,     180.314,   0.81472    ],
    "h":        [  0.5556,      0.68,      0.259      ],
    "mu":       [  1.0/1.4276,  1.0/1.79,  1.0/1.842  ],
    "n":        [  1.0,         1.0,       1.0        ],
    "m":        [  3,           3,         3          ],
    "r1":       [  2.00,        2.40,      2.50       ],
    "r2":       [  2.00*1.2,    2.40*1.2,  2.50*1.2   ],
    "or1":      [  2.00,        2.40,      3.00       ],
    "or2":      [  2.00*2.0,    2.40*2.0,  3.00*2.0   ],
    "bor1":     [  2.00,        2.40,      3.00       ],
    "bor2":     [  2.00*2.0,    2.40*2.0,  3.00*2.0   ],
    "Cmin":     [  1.00,        1.00,      1.00       ],
    "Cmax":     [  3.00,        3.00,      3.00       ]
    }


Albe_PRB_65_195124_PtC = {
    "__ref__":  "Albe K., Nordlund K., Averback R. S., Phys. Rev. B 65, 195124 (2002)",
    "el":       [  "Pt",   "C" ],
    "D0":       [  3.683,     5.3,        6.0       ],
    "r0":       [  2.384,     1.84,       1.39      ],
    "S":        [  2.24297,   1.1965,     1.22      ],
    "beta":     [  1.64249,   1.836,      2.1       ],
    "gamma":    [  8.542e-4,  9.7e-3,     2.0813e-4 ],
    "c":        [  34.0,      1.23,       330.0     ],
    "d":        [  1.1,       0.36,       3.5       ],
    "h":        [  1.0,       1.0,        1.0       ],
    "mu":       [  1.335,     0.0,        0.0       ],
    "n":        [  1.0,       1.0,        1.0       ],
    "m":        [  1,         1,          1         ],
    "r1":       [  2.9,       2.5,        1.7       ],
    "r2":       [  3.3,       2.8,        2.0       ]
    }


Henriksson_PRB_79_114107_FeC = dict(
    __ref__ = "Henriksson K.O.E., Nordlund K., Phys. Rev. B 79, 144107 (2009)",
    el      = [  "Fe", "C"  ],
    D0      = [  1.5,         4.82645134,   6.0          ],
    r0      = [  2.29,        1.47736510,   1.39         ],
    S       = [  2.0693109,   1.43134755,   1.22         ],
    beta    = [  1.4,         1.63208170,   2.1          ],
    gamma   = [  0.0115751,   0.00205862,   0.00020813   ],
    c       = [  1.2898716,   8.95583221,   330.0        ],
    d       = [  0.3413219,   0.72062047,   3.5          ],
    h       = [ -0.26,        0.87099874,   1.0          ],
    mu      = [  0.0,         0.0,          0.0          ],
    n       = [  1.0,         1.0,          1.0          ],
    m       = [  1,           1,            1            ],
    r1      = [  2.95,        2.3,          1.70         ],
    r2      = [  3.35,        2.7,          2.00         ]
    )


Kioseoglou_PSSb_245_1118_AlN = {
    "__ref__":  "Kioseoglou J., Komninou Ph., Karakostas Th., Phys. Stat. Sol. (b) 245, 1118 (2008)",
    "el":       [  "N",   "Al" ],
    "D0":       [  9.9100,     3.3407,     1.5000   ],
    "r0":       [  1.1100,     1.8616,     2.4660   ],
    "S":        [  1.4922,     1.7269,     2.7876   ],
    "beta":     [  2.05945,    1.7219,     1.0949   ],
    "gamma":    [  0.76612,    1.1e-6,     0.3168   ],
    "c":        [  0.178493,   100390,     0.0748   ],
    "d":        [  0.20172,    16.2170,    19.5691  ],
    "h":        [  0.045238,   0.5980,     0.6593   ],
    "mu":       [  0.0,        0.0,        0.0      ],
    "n":        [  1.0,        0.7200,     6.0865   ],
    "m":        [  1,          1,          1        ],
    "r1":       [  2.00,       2.19,       3.40     ],
    "r2":       [  2.40,       2.49,       3.60     ]
    }


# Brenner's original parameter set without the tables
Brenner_PRB_42_9458_C_I = {
    "__ref__":  "Brenner D., Phys. Rev. B 42, 9458 (1990) [potential I]",
    "el":       [  "C"  ],
    "D0":       [  6.325            ],
    "r0":       [  1.315            ],
    "S":        [  1.29             ],
    "beta":     [  1.5              ],
    "gamma":    [  0.011304         ],
    "c":        [  19.0             ],
    "d":        [  2.5              ],
    "h":        [  1.0              ],
    "mu":       [  0.0              ],
    "n":        [  1.0/(2*0.80469)  ],
    "m":        [  1                ],
    "r1":       [  1.70             ],
    "r2":       [  2.00             ]
    }

Brenner_PRB_42_9458_C_II = {
    "__ref__":  "Brenner D., Phys. Rev. B 42, 9458 (1990) [potential II]",
    "el":       [  "C"  ],
    "D0":       [  6.0          ],
    "r0":       [  1.39         ],
    "S":        [  1.22         ],
    "beta":     [  2.1          ],
    "gamma":    [  0.00020813   ],
    "c":        [  330.0        ],
    "d":        [  3.5          ],
    "h":        [  1.0          ],
    "mu":       [  0.0          ],
    "n":        [  1.0/(2*0.5)  ],
    "m":        [  1            ],
    "r1":       [  1.70         ],
    "r2":       [  2.00         ]
    }

if hasattr(_atomistica, 'Brenner'):
    class Brenner(Atomistica):
        """Karsten Albe-type bond order potential.
           (This is actually the same form used by Brenner, maybe rename.)
        """

        def __init__(self, **kwargs):
            apply(Atomistica.__init__, (self, _atomistica.Brenner, 100,), kwargs)


if hasattr(_atomistica, 'BrennerScr'):
    class BrennerScr(Atomistica):
        """Karsten Albe-type bond order potential.
           (This is actually the same form used by Brenner, maybe rename.)
        """

        def __init__(self, **kwargs):
            apply(Atomistica.__init__, (self, _atomistica.BrennerScr, 1000,), kwargs)


# Kumagai's Si potential
# Parameters
Kumagai_CompMaterSci_39_457_Si = {
  "__ref__":  "Kumagai T., Izumi S., Hara S., Sakai S., "
  "Comp. Mater. Sci. 39, 457 (2007)",
  "el":       [  "Si"        ],
  "A":        [  3281.5905   ],
  "B":        [  121.00047   ],
  "lambda1":  [  3.2300135   ],
  "lambda2":  [  1.3457970   ],
  "eta":      [  1.0000000   ],
  "delta":    [  0.53298909  ],
  "alpha":    [  2.3890327   ],
  "beta":     [  1           ],
  "c1":       [  0.20173476  ],
  "c2":       [  730418.72   ],
  "c3":       [  1000000.0   ],
  "c4":       [  1.0000000   ],
  "c5":       [  26.000000   ],
  "h":        [  -0.36500000 ],
  "r1":       [  2.70        ],
  "r2":       [  3.30        ],
  }

Kumagai_CompMaterSci_39_457_Si__Scr = \
    copy.deepcopy(Kumagai_CompMaterSci_39_457_Si)
Kumagai_CompMaterSci_39_457_Si__Scr.update({
        'r1':      [ 2.50     ],
        'r2':      [ 2.50*1.2 ],
        'or1':     [ 3.00     ],
        'or2':     [ 3.00*2.0 ],
        'bor1':    [ 3.00     ],
        'bor2':    [ 3.00*2.0 ],
        'Cmin':    [ 1.00     ],
        'Cmax':    [ 3.00     ],
        })

if hasattr(_atomistica, 'Kumagai'):
    class Kumagai(Atomistica):
        """Kumagai bond order potential.
        """

        def __init__(self, **kwargs):
            apply(Atomistica.__init__, (self, _atomistica.Kumagai, 100,), kwargs)

if hasattr(_atomistica, 'KumagaiScr'):
    class KumagaiScr(Atomistica):
        """Kumagai bond order potential.
        """

        def __init__(self, **kwargs):
            apply(Atomistica.__init__, (self, _atomistica.KumagaiScr, 1000,), kwargs)


# Tabulated EAM potentials
if hasattr(_atomistica, 'TabulatedEAM'):
    class TabulatedEAM(Atomistica):
        """Tabulated single-element EAM potential.
        """

        def __init__(self, **kwargs):
            apply(Atomistica.__init__, (self, _atomistica.TabulatedEAM, 1000,), kwargs)


if hasattr(_atomistica, 'TabulatedAlloyEAM'):
    class TabulatedAlloyEAM(Atomistica):
        """Tabulated multi-element (alloy) EAM potential.
        """

        def __init__(self, **kwargs):
            apply(Atomistica.__init__, (self, _atomistica.TabulatedAlloyEAM, 1000,),
                  kwargs)

