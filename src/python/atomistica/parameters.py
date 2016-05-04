# ======================================================================
# Atomistica - Interatomic potential library and molecular dynamics code
# https://github.com/Atomistica/atomistica
#
# Copyright (2005-2015) Lars Pastewka <lars.pastewka@kit.edu> and others
# See the AUTHORS file in the top-level Atomistica directory.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
# ======================================================================

"""
Database with parameters for specific parameterizations of the potentials
"""

from __future__ import division

import copy
from math import log, sqrt

# These should correspond to PAIR_INDEX and TRIPLET_INDEX_NS of src/macros.inc.
# Note that Python indices start at 0 while Fortran indices start at 1!
def pair_index(i, j, maxval):
    return min(i+j*maxval, j+i*maxval)-min(i*(i+1)//2, j*(j+1)//2)

def triplet_index(i, j, k, maxval):
    return k+maxval*(j+maxval*i)

# Mixing rules
def mix(p, key, rule):
    nel = len(p['el'])
    for i in range(nel):
        for j in range(i+1,nel):
            ii = pair_index(i,i,nel)
            jj = pair_index(j,j,nel)
            ij = pair_index(i,j,nel)
            p[key][ij] = rule(p[key][ii], p[key][jj])

def mix_arithmetic(p, key):
    mix(p, key, lambda x,y: (x+y)/2)

def mix_geometric(p, key):
    mix(p, key, lambda x,y: sqrt(x*y))



# Tersoff potential --- to be used with the *Tersoff* potential

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

Matsunaga_Fisher_Matsubara_Jpn_J_Appl_Phys_39_48_B_C_N = {
  "__ref__":  "Matsunaga K., Fisher C., Matsubara H., Jpn. J. Appl. Phys. 39, 48 (2000)",
  "el":       [  "C", "N", "B" ],
  "A":        [  1.3936e3,    -1.0,     -1.0,     1.1e4,      -1.0,     2.7702e2  ],
  "B":        [  3.4674e2,    -1.0,     -1.0,     2.1945e2,   -1.0,     1.8349e2  ],
  "xi":       [  1.0,         0.9685,   1.0025,   1.0,        1.1593,   1.0       ],
  "lambda":   [  3.4879,      -1.0,     -1.0,     5.7708,     -1.0,     1.9922    ],
  "mu":       [  2.2119,      -1.0,     -1.0,     2.5115,     -1.0,     1.5856    ],
  "omega":    [  1.0,         0.6381,   1.0,      1.0,        1.0,      1.0       ],
  "mubo":     [  0.0,         0.0,      0.0,      0.0,        0.0,      0.0       ],
  "m":        [  1,           1,        1,        1,          1,        1         ],
  "r1":       [  1.80,        -1.0,     -1.0,     2.0,         -1.0,    1.8       ],
  "r2":       [  2.10,        -1.0,     -1.0,     2.3,         -1.0,    2.1       ],
  "beta":     [  1.5724e-7,   1.0562e-1,   1.6e-6     ],
  "n":        [  7.2751e-1,   12.4498,     3.9929     ],
  "c":        [  3.8049e4,    7.9934e4,    5.2629e-1  ],
  "d":        [  4.3484e0,    1.3432e2,    1.5870e-3  ],
  "h":        [  -5.7058e-1,  -0.9973,     0.5        ],
  }
# Apply mixing rules
mix_geometric(Matsunaga_Fisher_Matsubara_Jpn_J_Appl_Phys_39_48_B_C_N, 'A')
mix_geometric(Matsunaga_Fisher_Matsubara_Jpn_J_Appl_Phys_39_48_B_C_N, 'B')
mix_arithmetic(Matsunaga_Fisher_Matsubara_Jpn_J_Appl_Phys_39_48_B_C_N, 'lambda')
mix_arithmetic(Matsunaga_Fisher_Matsubara_Jpn_J_Appl_Phys_39_48_B_C_N, 'mu')
mix_geometric(Matsunaga_Fisher_Matsubara_Jpn_J_Appl_Phys_39_48_B_C_N, 'r1')
mix_geometric(Matsunaga_Fisher_Matsubara_Jpn_J_Appl_Phys_39_48_B_C_N, 'r2')

Matsunaga_Fisher_Matsubara_Jpn_J_Appl_Phys_39_48_B_C_N__Scr = Matsunaga_Fisher_Matsubara_Jpn_J_Appl_Phys_39_48_B_C_N.copy()
Matsunaga_Fisher_Matsubara_Jpn_J_Appl_Phys_39_48_B_C_N__Scr.update({
    "m":        [  3,          3,      3,      3,          3,      3        ],
    "r1":       [  2.00,       -1.0,   -1.0,   2.00,       -1.0,   1.8      ],
    "r2":       [  2.00*1.2,   -1.0,   -1.0,   2.00*1.2,   -1.0,   1.8*1.2  ],
    "or1":      [  2.00,       -1.0,   -1.0,   3.00,       -1.0,   1.8      ],
    "or2":      [  2.00*2.0,   -1.0,   -1.0,   3.00*2.0,   -1.0,   1.8*2    ],
    "bor1":     [  2.00,       -1.0,   -1.0,   3.00,       -1.0,   1.8      ],
    "bor2":     [  2.00*2.0,   -1.0,   -1.0,   3.00*2.0,   -1.0,   1.8*2    ],
    "Cmin":     [  1.00,       1.00,   1.00,   1.00,       1.00,   1.00     ],
    "Cmax":     [  3.00,       3.00,   3.00,   3.00,       3.00,   3.00     ],
    })
mix_geometric(Matsunaga_Fisher_Matsubara_Jpn_J_Appl_Phys_39_48_B_C_N__Scr, 'r1')
mix_geometric(Matsunaga_Fisher_Matsubara_Jpn_J_Appl_Phys_39_48_B_C_N__Scr, 'r2')
mix_geometric(Matsunaga_Fisher_Matsubara_Jpn_J_Appl_Phys_39_48_B_C_N__Scr, 'or1')
mix_geometric(Matsunaga_Fisher_Matsubara_Jpn_J_Appl_Phys_39_48_B_C_N__Scr, 'or2')
mix_geometric(Matsunaga_Fisher_Matsubara_Jpn_J_Appl_Phys_39_48_B_C_N__Scr, 'bor1')
mix_geometric(Matsunaga_Fisher_Matsubara_Jpn_J_Appl_Phys_39_48_B_C_N__Scr, 'bor2')



# Karsten Albe's BOP --- to be used with the *Brenner* potential

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


# Juslin's W-C-H parameterization
Juslin_JAP_98_123520_WCH = {
    '__ref__': 'Juslin N., Erhart P., Traskelin P., Nord J., Henriksson K.O.E, Nordlund K., Salonen E., Albe K., J. Appl. Phys. 98, 123520 (2005)',
    'el':      [ 'W', 'C', 'H' ],
    'D0':      [ 5.41861,     6.64,      2.748,    0.0,  6.0,         3.6422,      0.0,  3.642,    4.7509  ],
    'r0':      [ 2.34095,     1.90547,   1.727,   -1.0,  1.39,        1.1199,     -1.0,  1.1199,   0.74144 ],
    'S':       [ 1.92708,     2.96149,   1.2489,   0.0,  1.22,        1.69077,     0.0,  1.69077,  2.3432  ],
    'beta':    [ 1.38528,     1.80370,   1.52328,  0.0,  2.1,         1.9583,      0.0,  1.9583,   1.9436  ],
    'gamma':   [ 0.00188227,  0.072855,  0.0054,   0.0,  0.00020813,  0.00020813,  0.0,  12.33,    12.33   ],
    'c':       [ 2.14969,     1.10304,   1.788,    0.0,  330.0,       330.0,       0.0,  0.0,      0.0     ],
    'd':       [ 0.17126,     0.33018,   0.8255,   0.0,  3.5,         3.5,         0.0,  1.0,      1.0     ],
    'h':       [-0.27780,     0.75107,   0.38912,  0.0,  1.0,         1.0,         0.0,  1.0,      1.0     ],
    'n':       [ 1.0,         1.0,       1.0,      0.0,  1.0,         1.0,         0.0,  1.0,      1.0     ],
    'alpha':   [ 0.45876, 0.0, 0.0, 0.45876, 0.0, 0.0, 0.45876, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0, 0.0, 4.0, 4.0, 0.0, 0.0, 0.0, 0.0, 4.0, 4.0, 0.0, 4.0, 4.0 ],
    'omega':   [ 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 2.94586, 4.54415, 1.0, 1.0, 1.0, 1.0, 0.33946, 0.22006, 1.0, 1.0, 1.0 ],
    'r1':      [ 3.20,        2.60,      2.68,     0.0,  1.70,        1.30,        0.0,  1.30,     1.10    ],
    'r2':      [ 3.80,        3.00,      2.96,     0.0,  2.00,        1.80,        0.0,  1.80,     1.70    ],
    }

Juslin_JAP_98_123520_WCH__Scr = copy.deepcopy(Juslin_JAP_98_123520_WCH)
Juslin_JAP_98_123520_WCH__Scr.update({
        'r1':      [ 3.20,        2.60,      2.68,     0.0,  1.70,        1.30,        0.0,  1.30,     1.10    ],
        'r2':      [ 3.80,        3.00,      2.96,     0.0,  2.00,        1.80,        0.0,  1.80,     1.70    ],
        'or1':     [ 3.20,        2.60,      2.68,     0.0,  1.70,        1.30,        0.0,  1.30,     1.10    ],
        'or2':     [ 3.80,        3.00,      2.96,     0.0,  2.00,        1.80,        0.0,  1.80,     1.70    ],
        'bor1':    [ 3.20,        2.60,      2.68,     0.0,  1.70,        1.30,        0.0,  1.30,     1.10    ],
        'bor2':    [ 3.80,        3.00,      2.96,     0.0,  2.00,        1.80,        0.0,  1.80,     1.70    ],
        'Cmin':    [ 3.20,        2.60,      2.68,     0.0,  1.70,        1.30,        0.0,  1.30,     1.10    ],
        'Cmax':    [ 3.80,        3.00,      2.96,     0.0,  2.00,        1.80,        0.0,  1.80,     1.70    ],
        'm':       [ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 ],
        })


# Kuopanportti's Fe-C-H parameterization
Kuopanportti_CMS_111_525_FeCH = {
    '__ref__' : 'Kuopanportti P., Hayward, N., Fu C., Kuronen A., Nordlund K., Comp. Mat. Sci. 111, 525 (2016)',
    'el':     [ 'Fe', 'C', 'H'],
    'D0':     [ 1.5,      4.82645134,  1.630,    0.0,  6.0,         3.6422,      0.0,  3.642,    4.7509  ],
    'r0':     [ 2.29,     1.47736510,  1.589,   -1.0,  1.39,        1.1199,     -1.0,  1.1199,   0.74144 ],
    'S':      [ 2.0693,   1.43134755,  4.000,    0.0,  1.22,        1.69077,     0.0,  1.69077,  2.3432  ],
    'beta':   [ 1.4,      1.63208170,  1.875,    0.0,  2.1,         1.9583,      0.0,  1.9583,   1.9436  ],
    'gamma':  [ 0.01158,  0.00205862,  0.01332,  0.0,  0.00020813,  0.00020813,  0.0,  12.33,    12.33   ],
    'c':      [ 1.2899,   8.95583221,  424.5,    0.0,  330.0,       330.0,       0.0,  0.0,      0.0     ],
    'd':      [ 0.3413,   0.72062047,  7.282,    0.0,  3.5,         3.5,         0.0,  1.0,      1.0     ],
    'h':      [-0.26,     0.87099874, -0.1091,   0.0,  1.0,         1.0,         0.0,  1.0,      1.0     ],
    'n':      [ 1.0,      1.0,         1.0,      0.0,  1.0,         1.0,         0.0,  1.0,      1.0     ],
    'alpha':  [ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0, 0.0, 4.0, 4.0, 0.0, 0.0, 0.0, 0.0, 4.0, 4.0, 0.0, 4.0, 4.0 ],
    'omega':  [ 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 2.94586, 4.54415, 1.0, 1.0, 1.0, 1.0, 0.33946, 0.22006, 1.0, 1.0, 1.0 ],
    'r1':     [ 2.95,     2.30,        2.2974,   0.0,  1.70,        1.30,        0.0,  1.30,     1.10    ],
    'r2':     [ 3.35,     2.70,        2.6966,   0.0,  2.00,        1.80,        0.0,  1.80,     1.70    ],
    }


Kuopanportti_CMS_111_525_FeCH__Scr = copy.deepcopy(Kuopanportti_CMS_111_525_FeCH)
Kuopanportti_CMS_111_525_FeCH__Scr.update({
    'r1':     [ 2.95,   2.30,   2.2974,  0.0,  1.70,  1.30,  0.0,  1.30,  1.10 ],
    'r2':     [ 3.35,   2.70,   2.6966,  0.0,  2.00,  1.80,  0.0,  1.80,  1.70 ],
    'or1':    [ 2.95,   2.30,   2.2974,  0.0,  1.70,  1.30,  0.0,  1.30,  1.10 ],
    'or2':    [ 3.35,   2.70,   2.6966,  0.0,  2.00,  1.80,  0.0,  1.80,  1.70 ],
    'bor1':   [ 2.95,   2.30,   2.2974,  0.0,  1.70,  1.30,  0.0,  1.30,  1.10 ],
    'bor2':   [ 3.35,   2.70,   2.6966,  0.0,  2.00,  1.80,  0.0,  1.80,  1.70 ],
    'Cmin':   [ 1.00,   1.00,   1.00,    0.0,  1.00,  1.00,  0.0,  1.00,  1.00 ],
    'Cmax':   [ 3.00,   3.00,   3.00,    0.0,  3.00,  3.00,  0.0,  3.00,  3.00 ],
    'm':      [ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 ],
    })


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




# Kumagai's Si potential --- to be used with the *Kumagai* potential
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

