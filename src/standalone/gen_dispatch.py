# ======================================================================
# Atomistica - Interatomic potential library and molecular dynamics code
# https://github.com/Atomistica/atomistica
#
# Copyright (2005-2015) Lars Pastewka <lars.pastewka@kit.edu> and others
# See the AUTHORS file in the top-level Atomistica directory.
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
#! /usr/bin/env python

import sys

###

def read_module_list(fn):
    mods = []

    f = open(fn, 'r')
    l = f.readline()
    while l:
        l  = l.strip()
        if len(l) > 0 and l[0] != '!' and l[0] != '#':
            f90name, f90class, name, features, methods = l.split(':')[0:5]
            if f90class.strip() == '':
                f90class = f90name + '_t'
            methods = [s.lower() for s in methods.split(',')]
            mods += [[f90name, f90class, name, features, methods]]
        l = f.readline()
    f.close()

    return mods

###

def write_dispatch_f90(mods, str, template_fn, fn):
    ftmp = open(template_fn, 'r')
    fout = open(fn, 'w')

    method = None

    l = ftmp.readline()
    while l:
        if '#define' in l:
            s = l.split(' ')
            s = s[1].split('(')
            method = s[0].lower()
        elif '#undef' in l:
            method = None
        if '{' in l and '}' in l:
            for f90name, f90class, name, features, methods in mods:
                if method is None or method in methods:
                    fout.write(l.format(classname=f90name, classtype=f90class))
        else:
            fout.write(l)
        l = ftmp.readline()

    ftmp.close()
    fout.close()

###

mods  = read_module_list('integrators.classes')
write_dispatch_f90(mods, 'integrators',
                   '../src/standalone/integrators_dispatch.template.f90',
                   'integrators_dispatch.f90')

mods  = read_module_list('callables.classes')
write_dispatch_f90(mods, 'callables',
                   '../src/standalone/callables_dispatch.template.f90',
                   'callables_dispatch.f90')

mods  = read_module_list('potentials.classes')
write_dispatch_f90(mods, 'potentials',
                   '../src/standalone/potentials_dispatch.template.f90',
                   'potentials_dispatch.f90')

mods  = read_module_list('coulomb.classes')
write_dispatch_f90(mods, 'coulomb',
                   '../src/standalone/coulomb_dispatch.template.f90',
                   'coulomb_dispatch.f90')
