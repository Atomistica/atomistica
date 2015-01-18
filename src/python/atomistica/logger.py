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

"""Logging for molecular dynamics."""

import weakref
import sys
import ase.units as units
# ase.parallel imported in __init__

class MDLogger:
    """Class for logging molecular dynamics simulations with some
    extended functionality.

    Parameters:
    dyn:           The dynamics.  Only a weak reference is kept.

    atoms:         The atoms.

    logfile:       File name or open file, "-" meaning standart output.

    stress=False:  Include stress in log.

    cell=True:     Include cell in log.

    volume=True:   Include volume in log.

    peratom=False: Write energies per atom.

    mode="a":      How the file is opened if logfile is a filename.
    """
    def __init__(self, dyn, atoms, logfile, header=True, stress=False,
                 cell=False, volume=False, peratom=False, hiprec=False,
                 mode="a"):
        import ase.parallel
        if ase.parallel.rank > 0:
            logfile="/dev/null"  # Only log on master
        if hasattr(dyn, "get_time"):
            self.dyn = weakref.proxy(dyn)
        else:
            self.dyn = None
        self.atoms = atoms
        self.natoms = atoms.get_number_of_atoms()
        if logfile == "-":
            self.logfile = sys.stdout
            self.ownlogfile = False
        elif hasattr(logfile, "write"):
            self.logfile = logfile
            self.ownlogfile = False
        else:
            self.logfile = open(logfile, mode)
            self.ownlogfile = True
        self.stress = stress
        self.cell = cell
        self.volume = volume
        self.peratom = peratom
        if hiprec:
            nf = '%20.12e'
        else:
            nf = '%12.4f '
        i = 1
        if self.dyn is not None:
            self.hdr = "# {0}:Time[ps]".format(i)
            self.fmt = nf
            i += 1
        else:
            self.hdr = "# "
            self.fmt = ""
        if self.peratom:
            self.hdr += "{0}:Etot/N[eV] {1}:Epot/N[eV]" \
                        "{2}:Ekin/N[eV] {3}:T[K]".format(i,i+1,i+2,i+3)
            self.fmt += 4*nf
            i += 4
        else:
            self.hdr += "{0}:Etot[eV] {1}:Epot[eV]" \
                        "{2}:Ekin[eV] {3}:T[K]".format(i,i+1,i+2,i+3)
            self.fmt += 4*nf
            i += 4
        if self.stress:
            self.hdr += "{0}:stress(xx) {1}:stress(yy) {2}:stress(zz)" \
                        "{3}:stress(xy) {4}:stress(yz) {5}:stress(zx)".format(i,i+1,i+2,i+3,i+4,i+5)
            self.fmt += 6*nf
            i += 6
        if self.cell:
            self.hdr += "{0}:cell(xx) {1}:cell(yy) {2}:cell(zz)" \
                        "{3}:cell(xy) {4}:cell(yz) {5}:cell(zx)".format(i,i+1,i+2,i+3,i+4,i+5)
            self.fmt += 6*nf
            i += 6
        if self.volume:
            self.hdr += "{0}:Vol [A^3]".format(i)
            self.fmt += nf
            i += 1
        self.fmt += "\n"
        if header:
            self.logfile.write(self.hdr+"\n")
            
    def __del__(self):
        self.close()

    def close(self):
        if self.ownlogfile:
            self.logfile.close()

    def __call__(self):
        epot = self.atoms.get_potential_energy()
        ekin = self.atoms.get_kinetic_energy()
        temp = ekin / (1.5 * units.kB * self.natoms)
        if self.peratom:
            epot /= self.natoms
            ekin /= self.natoms
        if self.dyn is not None:
            t = self.dyn.get_time() / (1000*units.fs)
            dat = (t,)
        else:
            dat = ()
        dat += (epot+ekin, epot, ekin, temp)
        if self.stress:
            dat += tuple(self.atoms.get_stress() / units.GPa)
        if self.cell:
            cell = self.atoms.get_cell()
            dat += ( cell[0,0], cell[1,1], cell[2,2],
                     (cell[1,2]+cell[2,1])/2,
                     (cell[2,0]+cell[0,2])/2,
                     (cell[0,1]+cell[1,0])/2 )
        if self.volume:
            dat += ( self.atoms.get_volume(), )
        self.logfile.write(self.fmt % dat)
        self.logfile.flush()
        
