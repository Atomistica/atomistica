#! /usr/bin/env python

#
# Parse all .f90 files in a certain directory (and subdirectories thereof)
# and scan metadata.
#

from __future__ import print_function

import getopt
import io
import os
import re
import sys

from meta import scanallmeta
from functools import reduce

###

def get_finterfaces(fn, include_list=None, tmpfilename='_cpp.tmp'):
    include_str = ''
    if include_list is not None:
        include_str = reduce(lambda x,y: x+' -I'+y, include_list, '')
    os.system('gfortran -x f95-cpp-input -E {0} {1} > {2}'.format(fn,
                                                                  include_str,
                                                                  tmpfilename))

    iface = re.compile('^\ *interface\ ',re.IGNORECASE)

    finterfaces = []

    f = io.open(tmpfilename, mode='r', encoding='latin-1')
    l = f.readline()
    while l:
        l = f.readline()
        if re.match(iface, l):
            l = iface.sub('', l).strip().lower()
            finterfaces += [l]
    f.close()

    os.remove(tmpfilename)
    # gfortran generates and empty .s file when just preprocessing
    fnroot, fnext = os.path.splitext(fn)
    if os.path.exists(fnroot+'.s'):
        os.remove(fnroot+'.s')

    return [finterface.lower() for finterface in finterfaces]

###

def get_module_list(metadata, interface, finterface_list=[], exclude_list=[],
                    include_list=[]):
    mods = []
    fns = []
    depalready = []

    # Loop over all files and find modules
    for path, metapath in metadata.items():
        for fn, meta in metapath.items():
            if 'interface' in meta:
                if meta['interface'] == interface:
                    classtype = meta['classtype']
                    classname = meta['classname']
                    try:
                        features = meta['features']
                    except:
                        features = ''                   
                    if not classname in exclude_list:
                        s = []
                        finterfaces_present = get_finterfaces(path+'/'+fn,
                                                              include_list)
                        mods += [ ( classtype[:-2], classtype, classname,
                                    features, finterfaces_present ) ]
                        if 'dependencies' in meta:
                            dependencies = meta['dependencies'].split(',')
                            for depfn in dependencies:
                                if not depfn in depalready:
                                    fns += [ path+'/'+depfn ]
                                    depalready += [ depfn ]
                        fns += [ path+'/'+fn ]

    return mods, fns

###

def write_interface_info(metadata, interface, finterface_list, exclude_list,
                         include_list, deffn, mkfn, cfgfn):
    fns = []

    deff = io.open(deffn, mode='a', encoding='latin-1')
    mkf = io.open(mkfn, mode='a', encoding='latin-1')
    cfgf = io.open(cfgfn, mode='a', encoding='latin-1')

    print(u'%s_MODS += \\' % interface.upper(), file=mkf)

    depalready = []
    for path, metapath in metadata.items():
        for fn, meta in metapath.items():
            if 'interface' in meta:
                if meta['interface'] == interface:
                    classtype = meta['classtype']
                    classname = meta['classname']
                    try:
                        features = meta['features']
                    except:
                        features = ''
                    if not classname in exclude_list:
                        s = ''
                        finterfaces_present = get_finterfaces(path+'/'+fn,
                                                              include_list)
                        if len(finterfaces_present) > 0:
                            s = reduce(lambda x,y: x+','+y, finterfaces_present[1:],
                                       finterfaces_present[0])
                        print('%s:%s:%s:%s:%s' % (classtype[:-2],
                                                          classtype, classname,
                                                          features, s), file=deff)
                        if 'dependencies' in meta:
                            dependencies = meta['dependencies'].split(',')
                            for depfn in dependencies:
                                depfn = os.path.basename(depfn)
                                if not depfn in depalready:
                                    print('\t%s \\' % depfn, file=mkf)
                                    depalready += [ depfn ]
                        print(u'\t%s \\' % fn, file=mkf)

                        print('#define HAVE_%s' % \
                            (classtype[:-2].upper()), file=cfgf)

    deff.close()
    print(u'', file=mkf)
    mkf.close()
    cfgf.close()

###

if __name__ == '__main__':
    optlist, args = getopt.getopt(sys.argv[1:], '',
                                  ['exclude=', 'has_finterface='])

    if len(args) < 5:
        raise RuntimeError('Syntax: listclasses.py <path> <interface> '
                           '<definition-file> <makefile> <config-file> '
                           '[-I<include directory>] '
                           '[--exclude=<exclude list>] '
                           '[--has_finterface=<interface list>]')

    path = args[0]
    interface = args[1]
    deffn = args[2]
    mkfn = args[3]
    cfgfn = args[4]

    exclude_list = []
    finterface_list = []
    include_list = []

    for key, value in optlist:
        if key == '--exclude':
            exclude_list = value.split(',')
        elif key == '--has_finterface':
            finterface_list = value.split(',')

    for key in args[5:]:
        if key[:2] == '-I':
            include_list += [key[2:]]
        else:
            raise RuntimeError('Unknown comand line argument: {0}'.format(key))

    print('Scanning metadata of all source files...')
    metadata = scanallmeta(path)
   
    print("Dumping information for classes that implement '{0}' interface..." \
        .format(interface))
    write_interface_info(metadata, interface, finterface_list, exclude_list,
                         include_list, deffn, mkfn, cfgfn)
