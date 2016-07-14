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

from __future__ import absolute_import, print_function

import glob
import os
import sys
import re

sys.path += [ 'src/python', 'tools' ]

import numpy as np

from numpy.distutils.core import setup, Extension

import versioneer

from meta import scanallmeta
from listclasses import get_module_list
from gen_factory import write_factory_f90, write_factory_c
from gen_factory import write_coulomb_factory_f90, write_coulomb_factory_c

###

cwd = os.getcwd()
srcdir = '%s/src' % cwd

###

#
# Module configuration
#

inc_dirs = [ ]
lib_dirs = [ ]
# atomisticalib is a dependency, build below by the setup command
libs = [ 'atomisticalib' ]
extra_link_args = [ ]

###

#
# LAPACK configuration from numpy
#

for k, v in np.__config__.__dict__.items():
    if re.match('lapack_.*_info', k):
        if 'library_dirs' in v:
            print("* Using LAPACK information from '%s' dictionary in " \
                "numpy.__config__" % k)
            print("    library_dirs = '%s'" % v['library_dirs'])
            print("    libraries = '%s'" % v['libraries'])
            lib_dirs += v['library_dirs']
            libs += v['libraries']
            # We use whichever lapack_*_info comes first
            break

###

#
# Build list of source files
#

lib_srcs = [ ]
mod_srcs = [ ]

# MDCORE support modules
lib_srcs += [ 'build/versioninfo.f90',
            ]

lib_srcs += [ ('%s/support/' % srcdir)+i for i in
              [ 'c_f.f90',
                'error.f90',
                'System.f90',
                'MPI_context.f90',
                'Units.f90',
                'PeriodicTable.f90',
                'c_linearalgebra.cpp',
                'f_linearalgebra.f90',
                'f_ptrdict.f90',
                'c_ptrdict.c',
                'io.f90',
                'f_logging.f90',
                'c_logging.c',
                'timer.f90',
                'tls.f90',
                'misc.f90',
                'data.f90',
                'simple_spline.f90',
                'nonuniform_spline.f90',
                'cutoff.f90',
                'histogram1d.f90',
                'supplib.f90',
                'atomistica.f90',
                ]
              ]

lib_srcs += [ ('%s/special/' % srcdir)+i for i in
              [ 'table2d.f90',
                'table3d.f90',
                'table4d.f90',
                'anderson_mixer.f90',
                ]
              ]

lib_srcs += [ ('%s/python/f90/' % srcdir)+i for i in
              [ 'python_particles.f90',
                'python_neighbors.f90',
                'particles_wrap.f90',
                'neighbors_wrap.f90',
                'python_helper.f90',
                ]
              ]

lib_srcs += [ ('%s/core/' % srcdir)+i for i in
              [ 'filter.f90',
                ]
              ]

# Generate versioninfo
os.system('sh src/gen_versioninfo.sh src build Python')

# Other stuff
mod_srcs += [ '%s/python/c/py_f.c' % (srcdir),
              '%s/python/c/particles.c' % (srcdir),
              '%s/python/c/neighbors.c' % (srcdir),
              '%s/python/c/coulomb.c' % (srcdir),
              '%s/python/c/coulomb_callback.c' % (srcdir),
              '%s/python/c/potential.c' % (srcdir),
              '%s/python/c/analysis.c' % (srcdir),
              '%s/python/c/atomisticamodule.c' % (srcdir),
              ]

###

#
# C and Fortran compile flags
#

inc_dirs += [ np.get_include(),
              'build',
              'src',
              'src/support',
              'src/potentials',
              'src/notb',
              'src/notb/dense',
              ]

lib_macros = [ ( 'NO_BIND_C_OPTIONAL', None ),
               ( 'MDCORE_PYTHON', None ),
               ]

###

#
# Scan all metadata and get list of potentials
#

#print 'Scanning f90 metadata in directory {0}...'.format(srcdir)
metadata = scanallmeta(['{0}/{1}'.format(srcdir, x) for x in 
                        ['notb', 'potentials', 'potentials_nonfree']])

#print 'Writing factories...'

# Coulomb modules
mods1, fns1 = get_module_list(metadata, 'coulomb',
                              include_list = inc_dirs)
if sys.version_info >= (3,0):
    lib_srcs += fns1
else:
    lib_srcs += [ fn.encode('ascii') for fn in fns1 ]

print('* Found the following Coulomb modules:')
for f90name, f90type, name, features, methods in mods1:
    print('    {0}'.format(name))


# Write coulomb factory
write_coulomb_factory_f90(mods1, 'coulomb', 'build/coulomb_factory_f90.f90')
write_coulomb_factory_c(mods1, 'coulomb',
                        '%s/python/c/coulomb_factory.template.c' \
                             % (srcdir),
                        'build/coulomb_factory_c.c',
                        '%s/python/c/coulomb_factory.template.h' \
                             % (srcdir),
                        'build/coulomb_factory_c.h')
lib_srcs += [ '%s/python/f90/coulomb_dispatch.f90' % (srcdir),
              'build/coulomb_factory_f90.f90',
              'build/coulomb_factory_c.c',
              ]

# Potential modules
mods2, fns2 = get_module_list(metadata, 'potentials',
                              include_list = inc_dirs)
if sys.version_info >= (3,0):
    lib_srcs += fns2
else:
    lib_srcs += [ fn.encode('ascii') for fn in fns2 ]

print('* Found the following potential modules:')
for f90name, f90type, name, features, methods in mods2:
    print('    {0}'.format(name))

write_factory_f90(mods2, 'potential', 'build/potentials_factory_f90.f90')
write_factory_c(mods2, 'potential', 
                '%s/python/c/factory.template.c' % srcdir,
                'build/potentials_factory_c.c',
                '%s/python/c/factory.template.h' % srcdir,
                'build/potentials_factory_c.h')
lib_srcs += [ 'build/potentials_factory_f90.f90',
              'build/potentials_factory_c.c',
              ]

f = open('build/have.inc', 'w')
print('#ifndef __HAVE_INC', file=f)
print('#define __HAVE_INC', file=f)
for classabbrev, classname, classtype, classfeatures, methods in mods1:
    print('#define HAVE_%s' % (classabbrev.upper()), file=f)
for classabbrev, classname, classtype, classfeatures, methods in mods2:
    print('#define HAVE_%s' % (classabbrev.upper()), file=f)
print('#endif', file=f)
f.close()

###

scripts = glob.glob('src/python/tools/*.py')

###

setup(
    name = 'atomistica',
    version = versioneer.get_version(),
    cmdclass = versioneer.get_cmdclass(),
    description = 'Atomistica is a library of interatomic potentials that is compatible with ASE and LAMMPS',
    maintainer = 'Lars Pastewka',
    maintainer_email = 'lars.pastewka@kit.edu',
    url = 'https://github.com/Atomistica/atomistica',
    download_url = 'https://github.com/Atomistica/atomistica/tarball/'+versioneer.get_version(),
    license = 'GPLv3',
    packages = [
        'atomistica'
        ],
    package_dir = {
        'atomistica': 'src/python/atomistica'
        },
    libraries = [
        ( 'atomisticalib', dict(
                sources = lib_srcs,
                include_dirs = inc_dirs,
                macros = lib_macros,
                extra_compiler_args = [ '-fPIC' ],
                )
          )
        ],
    ext_modules = [
        Extension(
            '_atomistica',
            mod_srcs,
            include_dirs = inc_dirs,
            library_dirs = lib_dirs,
            libraries = libs,
            extra_compile_args = [ '-fPIC' ],
            extra_link_args = extra_link_args,
            )
        ],
    scripts = scripts
    )
