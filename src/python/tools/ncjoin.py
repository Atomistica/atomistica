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
Join a individual NetCDF trajectory files into a single one.
"""

import os
import sys
from optparse import OptionParser

import numpy as np

from netCDF4 import Dataset

###

FRAME_DIM = 'frame'
ATOM_DIM = 'atom'

###

def get_nearest_indices(time, every):
    """
    Given a list of times, return indices that are most closely to evenly
    spaced time intervals of length *every*.
    """

    time = np.asarray(time)
    r = [ ]
    m = int(np.max(time)/every)+1
    last_j = -1
    for i in range(m):
        difftime = np.abs(time-i*every)
        j = np.argmin(difftime)
        if difftime[j] < 1.0 and j != last_j:
            r += [ j ]
            last_j = j
    return np.array(r)

###

def open_trajs(trajfns, time_var='time', test_var='coordinates', test_tol=1e-6):
    """
    Open multiple NetCDF trajectory files and check that they are in order.
    Remove overlap if files overlap. Returns list of tuples. Tuple contains
        filename, netCDF4 Dataset object, first frame, last frame
    """

    data_f = zip(trajfns, map(Dataset, trajfns))
    filtered_data_f = [ ]

    fn2, data2 = data_f[0]
    last_time = None
    for i in range(len(data_f)-1):
        fn1, data1 = data_f[i]
        fn2, data2 = data_f[i+1]

        test_first2 = data2.variables[test_var][0]

        last1 = len(data1.variables[time_var])-1
        maxdiff = np.max(np.abs(data1.variables[test_var][last1] - test_first2))
        while last1 >= 0 and maxdiff > test_tol:
            print 'Frame %i of %s does not agree with first frame of %s ' \
                  '(maximum difference %e). Checking frame %i.' % \
                  ( last1, fn1, fn2, maxdiff, last1-1 )
            last1 -= 1
            maxdiff = np.max(np.abs(data1.variables[test_var][last1] -
                                    test_first2))

        if last1 < 0:
            raise RuntimeError('%s and %s are not consecutive. It may help to '
                               'increase *test_tol*.' % ( fn1, fn2 ))

        data_slice = slice(0, last1)
        time = data1.variables[time_var][data_slice]
        if last_time is not None:
            time += last_time - time[0]
        filtered_data_f += [ ( fn1, data1, data_slice, time ) ]
        last_time = data1.variables[time_var][last1]

    time = data2.variables[time_var]
    if last_time is not None:
        time += last_time - time[0]
    filtered_data_f += [ ( fn2, data2, slice(0, len(time)), time ) ]

    return filtered_data_f

###

def filter_trajs(idata_f, every):
    # Create a list of all frames
    idata_oframe = reduce(lambda a,b: a+b,
                          [ [ ( fn, data, i, time[i] )
                              for i in range(len(time))[data_slice] ]
                            for fn, data, data_slice, time in idata_f ],
                          [])

    # Get indices that corresponds to roughly equally spaced time slots
    i = get_nearest_indices([time for fn, data, i, time in idata_oframe], every)
    if len(i) == 0:
        raise RuntimeError('No frames left after filtering.')
    else:
        idata_oframe = [ idata_oframe[_i] for _i in i ]

    # Consolidate into a list that contains per-file information
    filtered_idata_f = [ ]
    cur_fn, cur_data, cur_slice, cur_time = idata_oframe[0]
    cur_slice = [ ]
    cur_time = [ ]
    for fn, data, i, time in idata_oframe:
        if data is not cur_data:
            filtered_idata_f += [ ( cur_fn, cur_data, np.array(cur_slice),
                                    np.array(cur_time) ) ]
            cur_fn = fn
            cur_data = data
            cur_slice = [ ]
            cur_time = [ ]
        cur_slice += [ i ]
        cur_time += [ time ]
    filtered_idata_f += [ ( cur_fn, cur_data, np.array(cur_slice),
                            np.array(cur_time) ) ]
    return filtered_idata_f


### Sanity check

if os.path.exists('traj.nc'):
    raise RuntimeError('traj.nc exists already.')


### Parse command line options

parser = OptionParser()
parser.add_option('-e', '--every', dest='every', type='float',
                  help='copy only frames at times that are multiples of EVERY',
                  metavar='EVERY')
parser.add_option('-v', '--test_var', dest='test_var', default='coordinates',
                  help='use variable VAR to test if two frames are identical',
                  metavar='VAR')
parser.add_option('-t', '--test_tol', dest='test_tol', type='float', default=1e-6,
                  help='use tolerance TOL to test if two frames are identical',
                  metavar='TOL')
options, trajfns = parser.parse_args()
print 'every =', options.every, ', test_var =', options.test_var, \
      ', test_tol =', options.test_tol

if len(trajfns) == 0:
    raise RuntimeError('Please provide one or more files to concatenate.')


### Open input files and filter if requested

print 'Opening files and checking file order...'
idata_f = open_trajs(trajfns, test_var=options.test_var,
                     test_tol=options.test_tol)
if options.every is not None:
    idata_f = filter_trajs(idata_f, options.every)


# Create output file
odata = Dataset('traj.nc', 'w', clobber=False, format='NETCDF3_64BIT')


### Copy global attributes

for attr_str in idata_f[0][1].ncattrs():
    print "> creating global attribute '%s'..." % attr_str
    odata.setncattr(attr_str, idata_f[0][1].getncattr(attr_str))


### Copy everything else

cursor = 0
last_data = None
for trajfn, idata, data_slice, time in idata_f:
    print "Appending '%s' starting at frame %i..." % ( trajfn, cursor )
    print 'File contains %i relevant time slots: ' % len(time), time

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
                print "> creating attribute '%s' of variable '%s'..." % \
                      ( attr_str, var_str )
                ovar.setncattr(attr_str, var.getncattr(attr_str))

        if var.dimensions[0] == FRAME_DIM:
            print "Copying variable '%s'..." % var_str
            if var_str == 'time':
                var_data = time
            else:
                var_data = var[data_slice]
            odata.variables[var_str][cursor:] = var_data
        else:
            if not last_data:
                print "Copying variable '%s'..." % var_str
                odata.variables[var_str][:] = var[:]
            else:
                print "Checking variable '%s' for consistency across files..."%\
                      var_str
                if np.any(last_data.variables[var_str][:] != var[:]):
                    raise RuntimeError("Data for per-file variable '%s' "
                                       "differs in '%s' and '%s'." % 
                                       ( var_str, trajfns[i-1], trajfns[i] ))

    cursor += len(time)
    last_time = time[-1]
    if last_data:
        last_data.close()
    last_data = idata
last_data.close()
odata.close()
