"""Main MDCore module.
"""

import atexit

import _mdcore

from aseinterface import *

#
# Enabled logging
#

_mdcore.startup()

#
# Close logfile upon exit
#

atexit.register(_mdcore.shutdown)


