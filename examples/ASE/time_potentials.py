#! /usr/bin/env python

import timeit

import ase
import ase.io

from atomistica import *

###

pots = [ Brenner, BrennerScr, Tersoff, TersoffScr, Rebo2, Rebo2Scr ]

for pot in pots:
    a = ase.io.read('rho_2.9.traj')
    a.set_calculator(pot())
    t = timeit.timeit('a.get_potential_energy()',
    	              setup='from __main__ import a',
    			      number=100)

    print '{} {}'.format(pot, t)
