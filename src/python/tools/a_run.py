#!/usr/bin/env python
# ======================================================================
# Atomistica - Interatomic potential library and molecular dynamics code
# https://github.com/Atomistica/atomistica
#
# Copyright (2005-2020) Lars Pastewka <lars.pastewka@imtek.uni-freiburg.de>
# and others. See the AUTHORS file in the top-level Atomistica directory.
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
# Emacs: treat this as -*- python -*-

from __future__ import print_function
import os
import stat
import sys
import re
from optparse import OptionParser

from atomistica.hardware import ComputeCluster

defaults = { 'err'   : None,
             'mail'  : None,
             'mem'   : None,
             'name'  : None,
             'depth' : 1,
             'cores' : 2,
             'out'   : None,
             'arch'  : 'pbs',
             'script': 'run.py',
             'time'  : 86400, # one day in seconds
             'wd'    : None,
             'queue' : None,
             'smt'   : False,
             }

set = defaults

#......................................................
# functions

def s_from_dhms(time):
    """return seconds from dhms"""
    dhms_s = { 's' : 1, 'm' : 60, 'h' : 3600, 'd' : 86400 }
    time = time.lower()
    word_list = re.findall('\d*[^\d]*',time)
    seconds=0 
    for word in word_list:
        if word != '':
            sec = 1
            for t in list(dhms_s.keys()):
                nw = word.replace(t,'')
                if nw != word:
                    sec = dhms_s[t]
                    word = nw
                    break
            try:
                seconds += int(word) * sec
            except:
                raise RuntimeError('unknown format in timestring ' + time)
    return seconds

def minutes(secs):
    return int(secs // 60)

def unique_name(name):
    import string
    letters = list(string.letters)
    n = name
    while os.path.exists(n):
        n = name + letters.pop(0)
    return n

# handle command line options

parser = OptionParser(usage='%prog [options] [script ncores]')
parser.add_option("-m", "--mail", dest='mail',
                  help='Where to send an email about starting/ending of the job (def: read from environment variable GPAW_MAIL)')
parser.add_option("-M", "--Memory", action='count', default=None,
                  help='request large memory cores (host specific)')
parser.add_option("-n", "--name", dest='name',
                  help='Name of the job (def: name of parent directory)')
parser.add_option("-d", "--depth", dest='depth',
                  help='depth of directories for naming (def: 1)')
parser.add_option("-o", "--outfile", dest='outfile',
                  help='Name of the output file (def: script.out)')
parser.add_option("-p", "--parameters", dest='parameters',
                  help='Parameters to give to the script (def: empty)')
parser.add_option("-a", "--arch", dest='arch',
                  help='architecture (def: try to guess)')
parser.add_option("-t", "--time", dest='time',
                  help='Time (def: 86400=1140m=24h=1d=one day)')
parser.add_option("-q", "--queue", dest='queue',
                  help='queue to use (host specific)')
parser.add_option("-s", "--smt", dest='smt', action='count',
                  help='Simultaneous Multi-Threading (host specific)')
opt, args = parser.parse_args()
##print "opt=",opt
##print "args=",args

if opt.mail:
    set['mail'] = str(opt.mail)

if opt.Memory is not None:
    set['mem'] = True

if opt.name:
    set['name'] = str(opt.name)

if opt.depth:
    set['depth'] = int(opt.depth)

if opt.outfile:
    set['out'] = str(opt.outfile)

if opt.parameters:
    set['parameters'] = str(opt.parameters)

if opt.time:
    set['time'] = s_from_dhms(opt.time)

if opt.smt:
    set['smt'] = True

if opt.queue:
    set['queue'] = str(opt.queue)

if len(args):
    set['script'] = args[0]
    if len(args) > 1:
        try:
            set['cores'] = int(args[1])
        except ValueError:
            raise ValueError('Number of cores must be integer. ' +
                             'See gpaw-runscript -h')

# ............................................................

try:
    cc = ComputeCluster(opt.arch)
except Exception as ex:
    raise
    print(ex.message, end='')
    sys.exit()

print('using', cc.arch)

# ............................................................

if set['mail'] is None and 'GPAW_MAIL' in os.environ:
    set['mail'] = str(os.environ['GPAW_MAIL'])

parameter_ext = ''
if 'parameters' in set:
     parameter_ext += '_' + set['parameters'].replace(' ','_')

# set output files
if set['out'] is None:
    set['out'] = set['script'] + parameter_ext + ".out"
if set['err'] is None:
    set['err'] = set['script'] + parameter_ext + ".err"

# get the name from current working directory
if set['wd'] is None:
    set['wd'] = os.getcwd()
if set['name'] is None:
    nl = os.getcwd().split('/')[-set['depth']:]
    name = nl[0]
    for string in nl[1:]:
        name += '_' + string
    # avoid beginning with a number
    if name[0].isdigit():
        name = 'j' + name 
    set['name'] = name + parameter_ext 

print(cc.write(**set), 'written')



