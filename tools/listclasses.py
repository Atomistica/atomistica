#! /usr/bin/env python

#
# Parse all .f90 files in a certain directory (and subdirectories thereof)
# and scan metadata.
#

import getopt
import os
import re
import sys

from meta import scanallmeta

###

def get_finterfaces(fn):
    iface = re.compile('^\ *interface\ ',re.IGNORECASE)

    finterfaces = [ ]

    f = open(fn, 'r')
    l = f.readline()
    while l:
        l = f.readline()
        if re.match(iface, l):
            l = iface.sub('', l).strip().lower()
            finterfaces += [ l ]
    f.close()

    return finterfaces

###

def get_module_list(metadata, interface, finterface_list=[], exclude_list=[]):
    mods = [ ]
    fns = [ ]
    depalready = [ ]

    # Loop over all files and find modules
    for path, metapath in metadata.iteritems():
        for fn, meta in metapath.iteritems():
            if 'interface' in meta:
                if meta['interface'] == interface:
                    classtype = meta['classtype']
                    classname = meta['classname']
                    try:
                        features = meta['features']
                    except:
                        features = ''                   
                    if not classname in exclude_list:
                        s = [ ]
                        if len(finterface_list) > 0:
                            finterfaces_present = get_finterfaces(path+'/'+fn)
                            for finterface in finterface_list:
                                if finterface.lower() in finterfaces_present:
                                    s += [ True ]
                                else:
                                    s += [ False ]
                        mods += [ [ classtype[:-2], classtype, classname,
                                    features ] + s ]
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
                         deffn, mkfn, cfgfn):
    fns = [ ]

    deff = open(deffn, 'a')
    mkf = open(mkfn, 'a')
    cfgf = open(cfgfn, 'a')

    print >> mkf, '%s_MODS = \\' % interface.upper()

    depalready = [ ]
    for path, metapath in metadata.iteritems():
        for fn, meta in metapath.iteritems():
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
                        if len(finterface_list) > 0:
                            finterfaces_present = get_finterfaces(path+'/'+fn)
                            for finterface in finterface_list:
                                if finterface in finterfaces_present:
                                    s += ':1'
                                else:
                                    s += ':0'
                        print >> deff, '%s:%s:%s:%s%s' % (classtype[:-2],
                                                          classtype, classname,
                                                          features, s)
                        if 'dependencies' in meta:
                            dependencies = meta['dependencies'].split(',')
                            for depfn in dependencies:
                                depfn = os.path.basename(depfn)
                                if not depfn in depalready:
                                    print >> mkf, '\t%s \\' % depfn
                                    depalready += [ depfn ]
                        print >> mkf, '\t%s \\' % fn

                        print >> cfgf, '#define HAVE_%s' % \
                            (classtype[:-2].upper())

    deff.close()
    print >> mkf
    mkf.close()
    cfgf.close()

###

if __name__ == '__main__':
    optlist, args = getopt.getopt(sys.argv[1:], '',
                                  ['exclude=', 'has_finterface='])

    if len(args) != 5:
        raise RuntimeError('Syntax: listclasses.py <path> <interface> '
                           '<definition-file> <makefile> <config-file> '
                           '[--exclude=<exclude list>] '
                           '[--has_finterface=<interface list>]')

    path = args[0]
    interface = args[1]
    deffn = args[2]
    mkfn = args[3]
    cfgfn = args[4]

    exclude_list = [ ]
    finterface_list = [ ]

    for key, value in optlist:
        if key == '--exclude':
            exclude_list = value.split(',')
        elif key == '--has_finterface':
            finterface_list = value.split(',')

    metadata = scanallmeta(path)
   
    write_interface_info(metadata, interface, finterface_list, exclude_list,
                         deffn, mkfn, cfgfn)
