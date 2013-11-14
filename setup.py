import glob
import os
import sys

sys.path += [ 'src/python', 'tools' ]

import numpy as np

from numpy.distutils.core import setup, Extension

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
# Get MKL configuration
#

mklroot = os.getenv('MKLROOT')
if mklroot is None:
    mklroot = os.getenv('MKL_ROOT')
if mklroot is not None:
    for lib_dir in [ '%s/lib/em64t' % mklroot, 
                     '%s/lib/intel64' % mklroot ]:
        if os.path.exists(lib_dir):
            lib_dirs += [ lib_dir ]
    #libs += [ 'mkl_intel_lp64', 'mkl_intel_thread', 'mkl_lapack',
    #          'mkl_core', 'mkl_def', 'irc_s', 'iomp5', 'ifcore', 'ifport',
    #          'stdc++' ]
    libs += [ 'iomp5', 'ifcore', 'ifport', 'mkl_rt' ]
    #extra_link_args += [ '-mkl=sequential' ]
else:
    libs += [ 'blas', 'lapack' ]

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
                'python_filter.f90',
                'python_neighbors.f90',
                'particles_wrap.f90',
                'neighbors_wrap.f90',
                'python_helper.f90',
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
# Scan all metadata and get list of potentials
#

#print 'Scanning f90 metadata in directory {0}...'.format(srcdir)
metadata = scanallmeta(srcdir)

#print 'Writing factories...'

# Coulomb modules
mods1, fns1 = get_module_list(metadata, 'coulomb',
                              finterface_list = [ 'register_data',
                                                  'set_Hubbard_U' ])
lib_srcs += fns1


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
                              finterface_list = [ 'register_data',
                                                  'set_Coulomb' ],
                              )
lib_srcs += fns2

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
print >> f, '#ifndef __HAVE_INC'
print >> f, '#define __HAVE_INC'
for classabbrev, classname, classtype, s, s2 in mods1:
    print >> f, '#define HAVE_%s' % (classabbrev.upper())
for classabbrev, classname, classtype, s1, s2 in mods2:
    print >> f, '#define HAVE_%s' % (classabbrev.upper())
print >> f, '#endif'
f.close()

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

scripts = glob.glob('src/python/tools/*.py')

###

setup(
    name = 'atomistica',
    version = '0.1',
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
