"""Main MDCore module.
"""

import atexit

import _atomistica

from aseinterface import *

#
# Enabled logging
#

_atomistica.startup()

#
# Close logfile upon exit
#

atexit.register(_atomistica.shutdown)


