#!/usr/bin/env python3
"""
Helper script for Meson build to generate factory code.
Called during build time to scan metadata and generate factory files.
"""

import sys
import os

# Add paths for imports
srcroot = sys.argv[1]
builddir = sys.argv[2]

sys.path.insert(0, os.path.join(srcroot, 'tools'))
sys.path.insert(0, os.path.join(srcroot, 'src/python'))

from meta import scanallmeta
from listclasses import get_module_list
from gen_factory import write_factory_f90, write_factory_c
from gen_factory import write_coulomb_factory_f90, write_coulomb_factory_c

# Include directories for module scanning
inc_dirs = [
    os.path.join(builddir),
    os.path.join(srcroot, 'src'),
    os.path.join(srcroot, 'src/support'),
    os.path.join(srcroot, 'src/potentials'),
    os.path.join(srcroot, 'src/notb'),
    os.path.join(srcroot, 'src/notb/dense'),
]

# Scan all metadata
metadata = scanallmeta([
    os.path.join(srcroot, 'src/notb'),
    os.path.join(srcroot, 'src/potentials'),
    os.path.join(srcroot, 'src/potentials_nonfree')
])

# Coulomb modules
mods1, fns1 = get_module_list(metadata, 'coulomb', include_list=inc_dirs)

print('* Found the following Coulomb modules:')
for f90name, f90type, name, features, methods in mods1:
    print('    {0}'.format(name))

# Write coulomb factory
write_coulomb_factory_f90(
    mods1, 'coulomb',
    os.path.join(builddir, 'coulomb_factory_f90.f90')
)
write_coulomb_factory_c(
    mods1, 'coulomb',
    os.path.join(srcroot, 'src/python/c/coulomb_factory.template.c'),
    os.path.join(builddir, 'coulomb_factory_c.c'),
    os.path.join(srcroot, 'src/python/c/coulomb_factory.template.h'),
    os.path.join(builddir, 'coulomb_factory_c.h')
)

# Potential modules
mods2, fns2 = get_module_list(metadata, 'potentials', include_list=inc_dirs)

print('* Found the following potential modules:')
for f90name, f90type, name, features, methods in mods2:
    print('    {0}'.format(name))

write_factory_f90(
    mods2, 'potential',
    os.path.join(builddir, 'potentials_factory_f90.f90')
)
write_factory_c(
    mods2, 'potential',
    os.path.join(srcroot, 'src/python/c/factory.template.c'),
    os.path.join(builddir, 'potentials_factory_c.c'),
    os.path.join(srcroot, 'src/python/c/factory.template.h'),
    os.path.join(builddir, 'potentials_factory_c.h')
)

# Write have.inc file
with open(os.path.join(builddir, 'have.inc'), 'w') as f:
    print('#ifndef __HAVE_INC', file=f)
    print('#define __HAVE_INC', file=f)
    for classabbrev, classname, classtype, classfeatures, methods in mods1:
        print('#define HAVE_%s' % (classabbrev.upper()), file=f)
    for classabbrev, classname, classtype, classfeatures, methods in mods2:
        print('#define HAVE_%s' % (classabbrev.upper()), file=f)
    print('#endif', file=f)

# Also need to collect the source files for the library
# Write them to a file that Meson can read
with open(os.path.join(builddir, 'potential_sources.txt'), 'w') as f:
    for fn in fns1 + fns2:
        f.write(fn + '\n')

print('Factory generation complete!')
