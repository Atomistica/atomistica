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

#! /usr/bin/env python

"""Random analysis tools."""

import os

import numpy as np

from atomistica.io import write

###

VOROPP_PATH = 'voro++'

###

def get_enclosing_orthorhombic_box(cell):
    """
    Return lower and upper bounds of the orthorhombic box that encloses
    the parallelepiped spanned by the three cell vectors of cell.
    """
    
    # Cell vectors
    cx, cy, cz = cell
    
    # The cell has eight corners, one is at the origin, three at cx, cy, cz
    # and the last ones are...
    c1 = cx+cy
    c2 = cx+cz
    c3 = cy+cz
    c4 = cx+cy+cz
    
    corners = np.array([[0,0,0],cx,cy,cz,c1,c2,c3,c4])
    lower = np.min(corners, axis=0)
    upper = np.max(corners, axis=0)
    
    return lower, upper


def voropp_for_non_orthorhombic_cells(_a, q='%v', voropp_path=VOROPP_PATH,
                                      fast=False, dump=None):
    """
    Run voro++ on current configuration and return selected quantities.
    Parameter *q* can be a list of voro++ output quantities.
    Run 'voro++ -hc' to see options. Will take care of Lees-Edwards boundary
    conditions by embedding a sheared cell in its periodic images and then
    throwing away the border (will hence be slower).
    """
    
    # Make a copy because we will modify the Atoms object
    a = _a.copy()
    nat = len(a)
    # Wrap into cell
    a.set_scaled_positions(a.get_scaled_positions()%1.0)

    # shear_dx should go into the cell
    if 'shear_dx' in a.info:
        lx, ly, lz = a.get_cell().diagonal()
        cx, cy, cz = a.get_cell()
        assert abs(cz[0]) < 1e-12
        assert abs(cz[1]) < 1e-12
        shear_dx = a.info['shear_dx']
        cz[0] = shear_dx[0]
        cz[1] = shear_dx[1]
        a.set_cell([cx, cy, cz], scale_atoms=False)
        a.set_scaled_positions(a.get_scaled_positions()%1.0)
        
    cx, cy, cz = a.get_cell()
   
    if fast:
        # Make 2 copies of the box in each direction. Could fail!
        a *= ( 2, 2, 2 )

        # Translate such that the box with the lowest indices sits in the middle
        a.translate(cx/2+cy/2+cz/2)
    else:
        # Make 3 copies of the box in each direction
        a *= ( 3, 3, 3 )
       
        # Translate such that the box with the lowest indices sits in the middle
        a.translate(cx+cy+cz)
        
    # Wrap back to box
    a.set_scaled_positions(a.get_scaled_positions()%1.0)

    # Get enclosing box
    lower, upper = get_enclosing_orthorhombic_box(a.get_cell())
    elx, ely, elz = upper-lower
    
    # Shift and set cell such that the general system is enclosed in the
    # orthorhombic box
    a.translate(-lower)
    a.set_cell([elx,ely,elz], scale_atoms=False)
    
    # Dump snapshot for debugging purposes
    if dump:
        write(dump, a)
    
    # Do Voronoi analysis    
    x, y, z = a.get_positions().T
    f = open('tmp.voronoi', 'w')
    for jn, ( jx, jy, jz ) in enumerate(zip(x, y, z)):
        print >> f, jn, jx, jy, jz
    f.close()
    if isinstance(q, str) or isinstance(q, unicode):
        c = '%s' % format(q)
    else:
        c = reduce(lambda x,y: '{0} {1}'.format(x,y), map(lambda x: '%'+x, q))
    os.system('{0} -o -p -c "%i {1}" 0 {2} 0 {3} 0 {4} tmp.voronoi' \
              .format(voropp_path, c, elx, ely, elz))
    r = np.loadtxt('tmp.voronoi.vol', unpack=True)
    os.remove('tmp.voronoi')
    os.remove('tmp.voronoi.vol')
    
    # Sort particles according to their ids
    r = r[:,np.array(r[0,:], dtype=int)]
    # Use only the lowest indices (i.e. the box in the middle)
    if r.shape[0] == 2:
        return r[1,:nat]
    else:
        return r[1:,:nat]

###

def voropp(a, q='%v', voropp_path=VOROPP_PATH, fast=False, dump=None):
    """
    Run voro++ on current configuration and return selected quantities.
    Parameter *q* can be a list of voro++ output quantities.
    Run 'voro++ -hc' to see options.
    """
    cx, cy, cz = a.get_cell()
    lx, ly, lz = np.linalg.norm(cx), np.linalg.norm(cy), np.linalg.norm(cz)
    if abs(lx*ly*lz - a.get_volume()) > 1e-6 or 'shear_dx' in a.info:
        return voropp_for_non_orthorhombic_cells(a, q, voropp_path, fast, dump)

    x, y, z = a.get_positions().T
    f = open('tmp.voronoi', 'w')
    for jn, ( jx, jy, jz ) in enumerate(zip(x, y, z)):
        print >> f, jn, jx, jy, jz
    f.close()
    if isinstance(q, str) or isinstance(q, unicode):
        c = '%s' % format(q)
    else:
        c = reduce(lambda x,y: '{0} {1}'.format(x,y), map(lambda x: '%'+x, q))
    os.system('{0} -o -p -c "%i {1}" 0 {2} 0 {3} 0 {4} tmp.voronoi' \
              .format(voropp_path, c, lx, ly, lz))
    r = np.loadtxt('tmp.voronoi.vol', unpack=True)
    os.remove('tmp.voronoi')
    os.remove('tmp.voronoi.vol')
    
    # Sort particles according to their ids
    r = r[:,np.array(r[0,:], dtype=int)]
    if r.shape[0] == 2:
        return r[1,:]
    else:
        return r[1:,:]

###

def stress_invariants(s):
    """Receives a list of stress tensors and returns the three invariants.
       Return hydrostatic pressure, octahedral shear stress and J3
    """
    s = np.asarray(s)
    if s.shape == (6,):
        s = s.reshape(1,-1)
    elif s.shape == (3,3):
        s = s.reshape(1,-1,-1)
    if len(s.shape) == 3:
        s = np.transpose([s[:,0,0],s[:,1,1],s[:,2,2],
                          (s[:,0,1]+s[:,1,0])/2,
                          (s[:,1,2]+s[:,2,1])/2,
                          (s[:,2,0]+s[:,0,2])/2])
    I1 = s[:,0]+s[:,1]+s[:,2]
    I2 = s[:,0]*s[:,1]+s[:,1]*s[:,2]+s[:,2]*s[:,0]-s[:,3]**2-s[:,4]**2-s[:,5]**2
    I3 = s[:,0]*s[:,1]*s[:,2]+2*s[:,3]*s[:,4]*s[:,5]-s[:,3]**2*s[:,2]-s[:,4]**2*s[:,0]-s[:,5]**2*s[:,1]

    J2 = I1**2/3-I2
    J3 = 2*I1**3/27-I1*I2/3+I3

    # Return hydrostatic pressure, octahedral shear stress and J3
    return -I1/3, np.sqrt(2*J2/3), J3
