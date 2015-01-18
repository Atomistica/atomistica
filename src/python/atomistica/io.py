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
I/O convenience functions.
"""

import os

import ase.io
try:
    from ase.io import NetCDFTrajectory
except:
    pass

try:
    from ase.calculators.lammps import write_lammps_data
except:
    from ase.calculators.lammpsrun import write_lammps_data

from atomistica.mdcore_io import read_atoms, write_atoms

###

def read(fn, **kwargs):
    """
    Convenience function: Detect file extension and read via Atomistica or ASE.
    If reading a NetCDF files, frame numbers can be appended via '@'.
    e.g., a = read('traj.nc@5')
    """
    ext = fn[fn.rfind('.'):].split('@')
    if len(ext) == 1:
        if ext[0] == '.out' or ext[0] == '.dat':
            dirname = os.path.dirname(fn)
            if len(dirname) == 0:
                dirname = '.'
            cycfn = dirname+'/cyc.dat'
            if os.path.exists(cycfn):
                return read_atoms(fn, cycfn=cycfn)
            return read_atoms(fn)
        elif ext[0] == '.nc':
            traj = NetCDFTrajectory(fn, **kwargs)
            return traj[-1]
        else:
            return ase.io.read(fn, **kwargs)
    elif len(ext) == 2:
        if ext[0] == '.nc':
            frame = int(ext[1])
            fn = fn[:fn.rfind('@')]
            traj = NetCDFTrajectory(fn)
            return traj[frame]
        else:
            return ase.io.read(fn, **kwargs)
    else:
        return ase.io.read(fn, **kwargs)


def write(fn, a, **kwargs):
    """
    Convenience function: Detect file extension and write via Atomistica or ASE.
    Has support for writing LAMMPS data files.
    """
    ext = fn[fn.rfind('.'):].split('@')
    if ext[0] == '.out' or ext[0] == '.dat':
        return write_atoms(fn, a)
    elif ext[0] == '.lammps':
        return write_lammps_data(fn, a, velocities=True, **kwargs)
    elif ext[0] == '.nc':
        return NetCDFTrajectory(fn, 'w').write(a)
    else:
        return ase.io.write(fn, a, **kwargs)

