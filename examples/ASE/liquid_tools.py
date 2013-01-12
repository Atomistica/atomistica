from math import pi

import numpy as np

import ase
from ase.units import mol, fs

import mdcore.native as native

###

def random_solid(els, density):
    syms = [ ]
    nat = 0
    for sym, n in els:
        syms += n*[ sym ]
        nat += n
    r = np.random.rand(nat, 3)
    a = ase.Atoms(syms, positions=r, cell=[1,1,1])

    mass = np.sum(a.get_masses())
    a0 = ( 1e24*mass/(density*mol) )**(1./3)
    a.set_cell([a0,a0,a0], scale_atoms=True)

    return a

###

class PairAndAngleDistribution:
    def __init__(self, calc, paircutoff, anglecutoff, npairbins=100,
                 nanglebins=100, avgn=1000):
        self.calc = calc
        # Instantiate an MDCORE neighbor list object, *avgn* is the
        # average number of neighbors per atom.
        self.nl = native.Neighbors(avgn)
        # Ask the neighbor list to compute neighbors at least up to
        # cutoff *cutoff*. Neighbor list may contain farther atoms.
        self.nl.request_interaction_range(max(paircutoff, anglecutoff))

        self.paircutoff = paircutoff
        self.anglecutoff = anglecutoff
        self.npairbins = npairbins
        self.nanglebins = nanglebins

        self.pair_hist = None
        self.pair_hist_var = None
        self.angle_hist = None
        self.angle_hist_var = None
        self.navg = 0

    def __call__(self):
        # Construct neighbor list and return neighbors. Return are a
        # list of integers *i* and *j* that denote the neighbor pair.
        # List *r* contains the distance. If you want to loop over all
        # neighbors and get the distance use
        #    for ii, jj, rr in zip(i,j,r):
        i, j, dr, abs_dr = self.nl.get_neighbors(self.calc.particles, vec=True)
        pavg, pvar = native.pair_distribution(i, abs_dr, self.npairbins,
                                              self.paircutoff)
        aavg, avar = native.angle_distribution(i, j, dr, self.nanglebins,
                                               self.anglecutoff)

        if self.pair_hist is None:
            self.pair_hist = pavg
            self.pair_hist_var = pvar
            self.angle_hist = aavg
            self.angle_hist_var = avar
            self.navg = 1
        else:
            self.pair_hist += pavg
            self.pair_hist_var += pvar
            self.angle_hist += aavg
            self.angle_hist_var += avar
            self.navg += 1

    def get_pair_hist(self):
        return self.pair_hist/self.navg

    def get_angle_hist(self):
        return self.angle_hist/self.navg

    def get_pair_variance(self):
        return self.pair_hist_var/self.navg

    def get_angle_variance(self):
        return self.angle_hist_var/self.navg

###

class CoordinationCount:
    def __init__(self, dyn, calc, cutoff, maxc=10, logfile=None):
        self.dyn = dyn
        self.calc = calc

        self.cutoff = cutoff
        self.maxc = maxc

        self.hist = np.zeros(maxc, dtype=int)
        self.navg = 0

        self.logfile = logfile
        if isinstance(logfile, str):
            self.logfile = open(logfile, 'w')

        self.fmt = '{0:.2f}'
        for i in range(maxc+1):
            self.fmt += ' {{{0}}}'.format(i+1)
        self.fmt += '\n'


    def __call__(self):
        c = self.calc.nl.get_coordination_numbers(self.calc.particles,
                                                  self.cutoff)
        hist = np.bincount(c)
        if len(hist) < self.maxc:
            hist = np.append(hist, np.zeros(self.maxc-len(hist), dtype=int))
        self.hist += hist
        self.navg += 1

        if self.logfile is not None:
            self.hist /= self.navg
            self.hist = np.append(self.hist, np.sum(self.hist))
            self.logfile.write(self.fmt.format(self.dyn.get_time()/(1000*fs),
                                               *self.hist))
            self.hist = np.zeros(self.maxc, dtype=int)
            self.navg = 0
                              

    def get_hist(self):
        return self.hist/self.navg


