from math import log

from mdcore import *

el = 'C'
dt = 0.5 # fs

densities  = [ 2.0, 2.3, 2.5, 2.7, 2.9, 3.1, 3.3, 3.5 ]

quick_calc = Tersoff(**Tersoff_PRB_39_5566_Si_C)
calc = TersoffScr(**Tersoff_PRB_39_5566_Si_C__Scr)

