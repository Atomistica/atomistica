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
This file contains input and output filters for a deprecated file format used
at Fraunhofer IWM.
"""

from __future__ import print_function

import os

import numpy as np

import ase
from ase.data import atomic_masses
from ase.parallel import paropen

### Input

def read_atoms(fn, cycfn=None, pos_only=False, conv=1.0):
    """
    Read atom information from an atoms.dat file (i.e., tblmd, MDCORE input file)
    """
    f = paropen(fn, "r")

    l = f.readline().lstrip()
    while len(l) > 0 and ( l[0] == '#' or l[0] == '<' ):
        l = f.readline().lstrip()

    n_atoms = int(l)

    l = f.readline().lstrip()
    while len(l) > 0 and ( l[0] == '#' or l[0] == '<' ):
        l = f.readline().lstrip()

    l = f.readline().lstrip()
    while len(l) > 0 and ( l[0] == '#' or l[0] == '<' ):
        l = f.readline().lstrip()

    #
    # Read positions
    #

    forces = np.zeros( [ n_atoms, 3 ] )
    groups = np.zeros( [ n_atoms ] )
    gamma  = np.zeros( [ n_atoms ] )
    T      = np.zeros( [ n_atoms ] )

    ats = [ ]
    for i in range(n_atoms):
        s = l.split()
        #      type   x            y            z
        sym = None
        try:
            Z    = int(s[0])
            sym  = ase.data.chemical_symbols[Z]
        except:
            sym  = s[0]
        a = ase.Atom(sym, ( float(s[2])*conv, float(s[3])*conv, float(s[4])*conv ) )
        groups[i] = int(s[5])
        gamma[i]  = float(s[6])
        T[i]      = float(s[7])

        ats += [ a ]
        l = f.readline()
    this = ase.Atoms(ats, pbc=True)

    if not pos_only:
        while l and l == "":
            l = f.readline().strip()

        while l:
            key = l.strip(" <-#\r\n")

            if key.upper() == "VELOCITIES":
                for i in range(n_atoms):
                    s = f.readline().split()
                    m = this[i].mass
                    if m is None:
                        m = ase.data.atomic_masses[ase.data.chemical_symbols.index(this[i].symbol)]
                    this[i].momentum = ( m*float(s[0]), m*float(s[1]), m*float(s[2]) )

                l = None

            elif key.upper() == "FORCES":
                for i in range(n_atoms):
                    s = f.readline().split()
                    forces[i] =  np.array( [ float(s[0]), float(s[1]), float(s[2]) ] )

                l = None

            elif key.upper() == "CHARGES":
                for i in this:
                    l = f.readline()
                    if l and len(l.split()) == 1:
                        i.charge = float(l)

                l = None

            elif key.upper() == "CELL" or key.upper().split()[0:2] == ("BOX", "VECTORS" ):
                l1 = f.readline()
                l2 = f.readline()
                l3 = f.readline()

                this.set_cell( [ map(float, l1.split()), map(float, l2.split()), map(float, l3.split()) ] )

                l = None

            else:
                aux = [ ]
                l = f.readline().strip()
                while l and l[0] not in [ '<', '#' ]:
                    s = l.split()

                    aux += [ map(float, s) ]

                    l = f.readline().strip()

                if len(aux) == n_atoms:
                    this.set_array(key, np.asarray(aux))
                else:
                    print("Warning: Encountered field '%s' which does not seem to be per-atom data." % key)

            if l is None:
                l = f.readline().strip()
            while l and l == "":
                l = f.readline().strip()


    f.close()

    this.set_array("forces", forces)
    this.set_array("groups", groups)
    this.set_array("gamma", gamma)
    this.set_array("T", T)

    if cycfn:
        read_cyc(this, cycfn, conv=conv)

    return this


def read_cyc(this, fn, conv=1.0):
    """ Read the lattice information from a cyc.dat file (i.e., tblmd input file)
    """
    f = paropen(fn, "r")
    f.readline()
    f.readline()
    f.readline()
    f.readline()
    cell = np.array( [ [ 0.0, 0.0, 0.0 ], [ 0.0, 0.0, 0.0 ], [ 0.0, 0.0, 0.0 ] ] )
    l = f.readline()
    s = map(float, l.split())
    cell[0, 0] = s[0]*conv
    cell[1, 0] = s[1]*conv
    cell[2, 0] = s[2]*conv
    l = f.readline()
    s = map(float, l.split())
    cell[0, 1] = s[0]*conv
    cell[1, 1] = s[1]*conv
    cell[2, 1] = s[2]*conv
    l = f.readline()
    s = map(float, l.split())
    cell[0, 2] = s[0]*conv
    cell[1, 2] = s[1]*conv
    cell[2, 2] = s[2]*conv
    this.set_cell(cell)
    this.set_pbc(True)
    f.close()


### Output

atoms_default_fields = np.array( [ "positions", "momenta", "numbers", "magmoms", "groups", "gamma", "T", "masses" ] )

def write_atoms(fn, this, cycfn=None, conv=1.0, symbols=True):
    """
    Write atom information to an atoms.dat file (i.e., tblmd, MDCORE input file)
    """
    f = paropen(fn, "w")

    f.write("<--- Number of atoms\n")
    f.write("%i\n" % len(this))
    f.write("<--- Number of occupied orbitals\n")
    f.write("%f\n" % 0.0)
    f.write("<--- Atom positions\n")

    groups = None
    if this.has("groups"):
        groups = this.get_array("groups")

    gamma = None
    if this.has("gamma"):
        gamma = this.get_array("gamma")

    T = None
    if this.has("T"):
        T = this.get_array("T")

    for idx, i in enumerate(this):
        group_str = "1"
        if groups is not None:
            group_str = "%i" % groups[idx]
        gamma_str = "0.0"
        if gamma is not None:
            gamma_str = "%20.10e" % gamma[idx]
        T_str = "0.0"
        if T is not None:
            T_str = "%20.10e" % T[idx]

        sym = i.symbol

        r = i.position

        if symbols:
            f.write("%s   %f   %20.10e  %20.10e  %20.10e   %s %s %s\n" % (sym, ase.data.atomic_masses[ase.data.chemical_symbols.index(sym)], r[0]*conv, r[1]*conv, r[2]*conv, group_str, gamma_str, T_str))
        else:
            f.write("%i   %f   %20.10e  %20.10e  %20.10e   %s %s %s\n" % (ase.data.chemical_symbols.index(sym), ase.data.atomic_masses[ase.data.chemical_symbols.index(sym)], r[0]*conv, r[1]*conv, r[2]*conv, group_str, gamma_str, T_str))


    f.write("<--- Velocities\n")
    for i in this:
        m = i.mass
        if m is None:
            m = ase.data.atomic_masses[ase.data.chemical_symbols.index(i.symbol)]
        if i.momentum is not None:
            v = i.momentum/m
        else:
            v = [ 0.0, 0.0, 0.0 ]
        f.write("%20.10e %20.10e %20.10e\n" % ( v[0], v[1], v[2] ))

    f.write("<--- Forces\n")
    for i in this:
        f.write("0 0 0\n")

    f.write("<--- cell\n")
    cell = this.get_cell()
    f.write("%f %f %f\n" % tuple(cell[0, :]))
    f.write("%f %f %f\n" % tuple(cell[1, :]))
    f.write("%f %f %f\n" % tuple(cell[2, :]))

    for name, aux in this.arrays.items():
        if not name in atoms_default_fields:
            f.write("<--- %s\n" % name)
            if aux.dtype == int:
                if len(aux.shape) == 1:
                    for i in this:
                        f.write(" %i\n" % aux[i.index])
                else:
                    for i in this:
                        f.write(( aux.shape[1]*" %i" + "\n" ) % tuple(aux[i.index].tolist()))
            else:
                if len(aux.shape) == 1:
                    for i in this:
                        f.write(" %e\n" % aux[i.index])
                else:
                    for i in this:
                        f.write(( aux.shape[1]*" %e" + "\n" ) % tuple(aux[i.index].tolist()))

    f.close()

    if cycfn:
        write_cyc(cycfn, this, conv=conv)



def write_cyc(fn, this, conv=1.0):
    """ Write the lattice information to a cyc.dat file (i.e., tblmd input file)
    """

    lattice = this.get_cell()

    f = paropen(fn, "w")
    f.write("<------- Simulation box definition\n")
    f.write("<------- Barostat (on = 1, off = 0)\n")
    f.write("  0\n")
    f.write("<------- Box vectors (start)\n")
    f.write("  %20.10f %20.10f %20.10f\n" % (lattice[0][0]*conv, lattice[1][0]*conv, lattice[2][0]*conv))
    f.write("  %20.10f %20.10f %20.10f\n" % (lattice[0][1]*conv, lattice[1][1]*conv, lattice[2][1]*conv))
    f.write("  %20.10f %20.10f %20.10f\n" % (lattice[0][2]*conv, lattice[1][2]*conv, lattice[2][2]*conv))
    f.write("<------- Box vectors (end)\n")
    f.write("  %20.10f %20.10f %20.10f\n" % (lattice[0][0]*conv, lattice[1][0]*conv, lattice[2][0]*conv))
    f.write("  %20.10f %20.10f %20.10f\n" % (lattice[0][1]*conv, lattice[1][1]*conv, lattice[2][1]*conv))
    f.write("  %20.10f %20.10f %20.10f\n" % (lattice[0][2]*conv, lattice[1][2]*conv, lattice[2][2]*conv))
    f.write("<------- Mass and gamma of the box (used in connection with the barostat)\n")
    f.write("  240 0.005\n")
    f.write("<------- Stress tensor (start)\n")
    f.write("  0 0 0\n")
    f.write("  0 0 0\n")
    f.write("  0 0 0\n")
    f.write("<------- Stress tensor (end)\n")
    f.write("  0 0 0\n")
    f.write("  0 0 0\n")
    f.write("  0 0 0\n")
    f.close()
