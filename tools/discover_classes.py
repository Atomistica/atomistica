#! /usr/bin/env python

#
# Parse all .f90 files in a certain directory (and subdirectories thereof)
# and scan metadata.
#

import sys

from meta import scanmeta

if __name__ == '__main__':
    for fn in sys.argv[1:]:
        metadata = scanmeta(fn)
