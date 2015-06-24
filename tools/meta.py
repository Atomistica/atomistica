"""Meta information is stored at the beginning of a Fortran file. This module
   is intended to scan directories for this information.
"""

from __future__ import print_function

import io
import os


srcexts = [ 'f', 'f90', 'f95' ]


def scanmeta(f):
    """Scan file headers for @meta ... @endmeta information and store that into
       a dictionary.
    """
    print(f)
    if isinstance(f, str):
        f = io.open(f, mode='r', encoding='latin-1')

    done = False

    l = f.readline()
    s = None
    while l and s is None:
        i = l.find('!')
        if i >= 0:
            l = l[i+1:]
            i = l.find('@meta')
            if i >= 0:
                l = l[i+5:]
                i = l.find('@endmeta')
                if i >= 0:
                    s = l[:i]
                    done = True
                else:
                    s = l
        l = f.readline()
    
    if not done and not l:
        return { }

    while l and not done:
        i = l.find('!')
        if i >= 0:
            l = l[i+1:]
            i = l.find('@endmeta')
            if i >= 0:
                s += ' '+l[:i]
                done = True
            else:
                s += ' '+l

        l = f.readline()

    s = map(lambda x: x.split(':'), s.split())
    d = { }
    for x in s:
        if len(x) > 2 or len(x) == 0:
            raise RuntimeError('Syntax error in meta information.')
        elif len(x) == 2:
            d[x[0]] = x[1]
        else:
            d[x[0]] = None

    return d


def _scanallmeta(dirname, fns):
    d = {}
    for fn in fns:
        fullfn = dirname+'/'+fn
        if os.path.isfile(fullfn):
            if fn.split('.')[-1] in srcexts:
                d[fn] = scanmeta(fullfn)
    return d


def scanallmeta(path):
    d = {}
    if isinstance(path, str):
        path = [path]
    for p in path:
        for dirpath, dirnames, filenames in os.walk(p, followlinks=True):
            d[dirpath] = _scanallmeta(dirpath, filenames)
    return d
