#! /usr/bin/env python

import sys

import numpy as np

import ase
from ase.io import read, write
from ase.optimize import FIRE
from ase.md import Langevin
from ase.units import mol, fs, kB

from atomistica.logger import MDLogger

from liquid_tools import *

###

# For coordination counting
#densities  = [ 2.0, 2.3, 2.5, 2.7, 2.9, 3.1, 3.3, 3.5 ]

# Default parameters
T1            = 10000
T2            = 5000
T3            = 300
nat           = 4001
g2_cutoff     = 5.0
coord_cutoff  = 1.85
nbins         = 100
time          = 1e4
fmax          = 10.0

import getopt
optlist, args = getopt.getopt(sys.argv[1:], '',
                              [ 'T1=', 'T2=', 'T3=', 'nat=', 'g2_cutoff=',
                                'coord_cutoff=', 'nbins=', 'time=' ])

assert len(args) == 0
for key, value in optlist:
    if key == '--T1':
        T1 = float(value)
    elif key == '--T2':
        T2 = float(value)
    elif key == '--T3':
        T3 = float(value)
    elif key == '--nat':
        nat = int(value)
    elif key == '--g2_cutoff':
        g2_cutoff = float(value)
    elif key == '--coord_cutoff':
        coord_cutoff = float(value)
    elif key == '--nbins':
        nbins = int(value)
    elif key == '--time':
        time = float(value)

###

print '# T1 = ', T1
print '# T2 = ', T2
print '# T3 = ', T3
print '# nat = ', nat
print '# g2_cutoff = ', g2_cutoff
print '# coord_cutoff = ', coord_cutoff
print '# nbins = ', nbins
print '# time = ', time

###

sys.path += [ '.' ]
from calcs import el, dt, quick_calc, calc, densities

els = [(el, nat)]

###

f = open('hyb.out', 'w')
for density in densities:
    print 'Running for density {0}...'.format(density)
    if isinstance(density, str):
        a = read(density)
    else:
        a = random_solid(els, density)

    # Relax with the unscreened potential
    print 'Relax with quick potential...'
    a.set_calculator(quick_calc)
    FIRE(a).run(fmax=fmax, steps=10000)

    # Relax with the screened potential
    print 'Relax with proper potential...'
    a.set_calculator(calc)
    FIRE(a).run(fmax=fmax, steps=10000)

    # Langevin quench to T1
    print 'Langevin quench to {0}...'.format(T1)
    Langevin(a, dt*fs, T1*kB, 1.0/(500*fs),
             logfile='-', loginterval=int(100/dt)).run(int(time/dt))

    # Langevin quench to T2
    print 'Langevin quench to {0}...'.format(T2)
    Langevin(a, dt*fs, T2*kB, 1.0/(500*fs),
             logfile='-', loginterval=int(100/dt)).run(int(time/dt))

    # Langevin quench to T3
    print 'Langevin quench to {0}...'.format(T3)
    dyn = Langevin(a, dt*fs, T3*kB, 1.0/(500*fs), logfile='-',
                   loginterval=int(100/dt))
    dyn.run(int(time/dt))

    # Collect pair distribution function
    print 'Collect pair distribution function...'
    p = PairAndAngleDistribution(a.get_calculator(), g2_cutoff, coord_cutoff,
                                 npairbins=nbins, nanglebins=nbins)
    dyn.attach(p, interval=int(100/dt))
    dyn.run(int(time/dt))

    print 'Writing files...'

    # Write snapshot
    write('rho_{0}.traj'.format(density), a)

    # Write pair distribution function
    r = (np.arange(nbins)+0.5)*g2_cutoff/nbins
    hist = p.get_pair_hist()
    variance = p.get_pair_variance()
    np.savetxt('g2_{0}.out'.format(density), np.transpose([r, hist, variance]))

    # Write angle distribution function
    r = (np.arange(nbins)+0.5)*pi/nbins
    hist = p.get_angle_hist()
    variance = p.get_angle_variance()
    np.savetxt('angle_{0}.out'.format(density), np.transpose([r, hist, variance]))

    # Count coordination numbers
    print 'Calculation coordination numbers...'
    c = np.zeros(12, dtype=int)
    for i in range(len(a)):
        c[calc.nl.coordination(calc.particles, i, coord_cutoff)] += 1
    assert np.sum(c) == len(a)
    if isinstance(density, str):
        density = np.sum(a.get_masses())/a.get_volume() * 1e24/mol
    np.savetxt(f, [ np.append([density], np.array(c, dtype=float)/len(a)) ])
    f.flush()
f.close()
