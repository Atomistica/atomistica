#! /usr/bin/env python

"""Random analysis tools."""

import os

import numpy as np

###

VOROPP_PATH = 'voro++'

###

def voropp(a, q='v', voropp_path=VOROPP_PATH):
    """Run voro++ on current configuration and return selected quantities.
       Parameter *q* can be a list of voro++ output quantities.
       Run 'voro++ -hc' to see options.
    """
    lx, ly, lz = a.get_cell().diagonal()
    if abs(lx*ly*lz - a.get_volume()) > 1e-6:
        raise RuntimeError('Atomic volume computation is only supported for orthogonal cells.')
    n = a.get_atomic_numbers()
    x, y, z = a.get_positions().T
    f = open('tmp.voronoi', 'w')
    for jn, jx, jy, jz in zip(n, x, y, z):
        print >> f, jn, jx, jy, jz
    f.close()
    if isinstance(q, str) or isinstance(q, unicode):
        c = '%{0}'.format(q)
    else:
        c = reduce(lambda x,y: '{0} {1}'.format(x,y), map(lambda x: '%'+x, q))
    os.system('{0} -o -p -c "{1}" 0 {2} 0 {3} 0 {4} tmp.voronoi' \
              .format(voropp_path, c, lx, ly, lz))
    r = np.loadtxt('tmp.voronoi.vol', unpack=True)
    os.remove('tmp.voronoi')
    os.remove('tmp.voronoi.vol')
    return r

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
