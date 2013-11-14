#! /usr/bin/env python

# ======================================================================
# Atomistica - Interatomic potential library
# https://github.com/pastewka/atomistica
# Lars Pastewka, lars.pastewka@iwm.fraunhofer.de, and others.
# See the AUTHORS file in the top-level Atomistica directory.
#
# Copyright (2005-2013) Fraunhofer IWM
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
Join a number of NetCDF trajectory files into a single one.
"""

import os
import sys

import numpy as np

from netCDF4 import Dataset

###

TIME_TOL = 1e-6
FRAME_DIM = 'frame'
ATOM_DIM = 'atom'

###

if os.path.exists('traj.nc'):
    raise RuntimeError('traj.nc exists already.')

odata = Dataset('traj.nc', 'w', clobber=False, format='NETCDF3_64BIT')
idata_f = map(Dataset, sys.argv[1:])
idata_f.sort(key=lambda x: x.variables['time'][0])

### Global attributes

for attr_str in idata_f[0].ncattrs():
    print "> creating global attribute '%s'..." % attr_str
    odata.setncattr(attr_str, idata_f[0].getncattr(attr_str))

### Copy stuff

cursor = 0
last_time = -1.0
last_data = None
for i, idata in enumerate(idata_f):
    print "Appending '%s' starting at frame %i..." % ( sys.argv[i+1], cursor )
    time = idata.variables['time'][:]
    j = 0
    while time[j] < last_time+TIME_TOL:
        j += 1

    for var_str, var in idata.variables.iteritems():
        if var_str not in odata.variables:
            print "> creating variable '%s'..." % var_str
            for dim_str in var.dimensions:
                if dim_str not in odata.dimensions:
                    print "> creating dimension '%s'......" % dim_str
                    dim = idata.dimensions[dim_str]
                    if dim.isunlimited():
                        odata.createDimension(dim_str)
                    else:
                        odata.createDimension(dim_str, len(dim))
            odata.createVariable(var_str, var.dtype, var.dimensions)
            ovar = odata.variables[var_str]
            for attr_str in var.ncattrs():
                print "> creating attribute '%s' of variable '%s'..." % ( attr_str, var_str )
                ovar.setncattr(attr_str, var.getncattr(attr_str))

        print "Copying variable '%s'..." % var_str
        if var.dimensions[0] == FRAME_DIM:
            odata.variables[var_str][cursor:] = var[j:]
        else:
            if not last_data:
                odata.variables[var_str][:] = var[:]
            else:
                if np.any(last_data.variables[var_str][:] != var[:]):
                    raise RuntimeError("Data for per-file variable '%s' "
                                       "differs in '%s' and '%s'." % 
                                       ( var_str, sys.argv[i], sys.argv[i+1] ))

    cursor += len(idata.dimensions[FRAME_DIM])-j
    if last_data:
        last_data.close()
    last_data = idata
last_data.close()
odata.close()
