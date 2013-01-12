#! /usr/bin/env python

import sys

import numpy as np

from ase.io import write
from ase.optimize import FIRE
from mdcore.io import read

import mdcore

###

# Default parameters
potstr = 'Tersoff'
fmax = 1e-3
outfn = 'fire.traj'

###

import getopt
optlist, args = getopt.getopt(sys.argv[1:], '',
                              [ 'pot=', 'fmax=', 'outfn=' ])

assert len(args) == 1
infn = args[0]
for key, value in optlist:
    if key == '--pot':
        potstr = value
    elif key == '--fmax':
        fmax = float(value)
    elif key == '--outfn':
        outfn = value

###

print '# infn = ', infn
print '# outfn = ', outfn
print '# pot = ', potstr
print '# fmax = ', fmax

###

a = read(infn)

print '{0} atoms.'.format(len(a))

potclass = getattr(mdcore, potstr)
a.set_calculator(potclass())

FIRE(a).run(fmax=fmax)

write(outfn, a)
