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

def file_from_template(templatefn, fn, keywords):
    ftmp = open(templatefn, 'r')
    fout = open(fn, 'w')

    l = ftmp.readline()
    while l:
        fout.write(l % keywords)
        l = ftmp.readline()

    ftmp.close()
    fout.close()

###

def read_module_list(fn):
    mods = [ ]

    f = open(fn, 'r')
    l = f.readline()
    while l:
        l  = l.strip()
        if len(l) > 0 and l[0] != '!' and l[0] != '#':
            f90name, f90class, name = l.split(':')[0:3]
            if f90class.strip() == '':
                f90class = f90name + '_t'
            mods += [ [ f90name, f90class, name ] ]
        l = f.readline()
    f.close()

    return mods

###

def write_1to1_factory_f90(mods, str, fn):
    f = open(fn, 'w')

    f.write('#include \"macros.inc\"\n\n' +
            'module %s_factory\n' % str +
            '  use supplib\n' +
            '#include "%s.inc"\n' % str)
    #for f90name, f90class, name in mods:
    #    f.write('  use %s\n' % f90name)
    f.write('  use %s, only: %s_t\n' % ( str, str ) +
            '  implicit none\n\n' +
            'contains\n\n')

    for f90name, f90class, name in mods:
        f.write('subroutine %s_new(this_cptr, cfg, m) bind(C)\n' % f90name +
                '  use, intrinsic :: iso_c_binding\n\n' +
                '  implicit none\n\n' +
                '  type(c_ptr), value  :: this_cptr\n' +
                '  type(c_ptr), value  :: cfg\n' +
                '  type(c_ptr), intent(out)  :: m\n' +
                '  type(%s_t), pointer  :: this\n' % str +
                '  call c_f_pointer(this_cptr, this)\n' +
                '  !if (.not. allocated(this%%%s)) then\n' % f90name +
                '     allocate(this%%%s)\n' % f90name +
                '  !else\n' +
                '  !   stop 999\n' +
                '  !endif\n' +
                '  call register(this%%%s, cfg, m)\n' % f90name +
                'endsubroutine %s_new\n\n\n' % f90name)
    
    f.write('endmodule %s_factory\n' % str)
    f.close()

###
    
def write_1ton_factory_f90(mods, str, fn):
    f = open(fn, 'w')

    f.write('#include \"macros.inc\"\n\n' +
            'module %s_factory\n' % str +
            '  use supplib\n' +
            '#include "%s.inc"\n' % str)
    #for f90name, f90class, name in mods:
    #    f.write('  use %s\n' % f90name)
    f.write('  use %s, only: %s_t\n' % ( str, str ) +
            '  implicit none\n\n' +
            'contains\n\n')

    for f90name, f90class, name in mods:
        f.write('subroutine %s_new(this_cptr, cfg, m) bind(C)\n' % f90name +
                '  use, intrinsic :: iso_c_binding\n\n' +
                '  implicit none\n\n' +
                '  type(c_ptr), value  :: this_cptr\n' +
                '  type(c_ptr), value  :: cfg\n' +
                '  type(c_ptr), intent(out)  :: m\n\n' +
                '  type(%s), allocatable  :: tmp(:)\n\n' % f90class +
                '  type(%s_t), pointer  :: this\n' % str +
                '  call c_f_pointer(this_cptr, this)\n' +
                '  if (.not. allocated(this%%%s)) then\n' % f90name +
                '     allocate(this%%%s(1))\n' % f90name +
                '     call register(this%%%s(1), cfg, m)\n' % f90name +
                '  else\n' +
                '     allocate(tmp(size(this%%%s)))\n' % f90name +
                '     tmp  = this%%%s\n' % f90name +
                '     deallocate(this%%%s)\n' % f90name +
                '     allocate(this%%%s(size(tmp)+1))\n' % f90name +
                '     this%%%s(1:size(tmp))  = tmp\n' % f90name +
                '     call register(this%%%s(size(tmp)+1), cfg, m)\n' % f90name +
                '     deallocate(tmp)\n' +
                '  endif\n' +
                'endsubroutine %s_new\n\n\n' % f90name)
    
    f.write('endmodule %s_factory\n' % str)
    f.close()

###

def write_factory_c(mods, str, c_dispatch_template, c_dispatch_file, h_dispatch_template, h_dispatch_file, compiler, system):
    d = { }

    d['disclaimer'] = 'This file has been autogenerated. DO NOT MODIFY.'
    d['name'] = str
    d['n_classes'] = len(mods)

    #
    # Prototypes
    #

    s = ''
    for f90name, f90class, name in mods:
        s += 'void %s_new(void *, section_t *, section_t **);\n' % f90name

    d['prototypes'] = s

    #
    # Classes
    #

    s = '%s_class_t %s_classes[N_CLASSES] = {\n' % ( str, str )
    for f90name, f90class, name in mods:
        s += '  {\n'
        s += '    \"%s\",\n' % name
        s += '    %s_new,\n' % f90name
        s += '  },\n'
            
    s = s[:-2] + '\n};\n'

    d['classes'] = s  

    #
    # Write the dispatch module
    #

    d['dispatch_header'] = h_dispatch_file.split('/')[-1]

    file_from_template(c_dispatch_template, c_dispatch_file, d)
    file_from_template(h_dispatch_template, h_dispatch_file, d)

###

srcdir, compiler, machine, system  = sys.argv[1:5]


mods  = read_module_list('integrators.classes')
write_1to1_factory_f90(mods, 'integrators', 'integrators_factory_f90.f90')
write_factory_c(mods, 'integrators', 
                srcdir + '/factory.template.c', 'integrators_factory_c.c',
                srcdir + '/factory.template.h', 'integrators_factory_c.h',
                compiler, system)


mods  = read_module_list('potentials.classes')
write_1ton_factory_f90(mods, 'potentials', 'potentials_factory_f90.f90')
write_factory_c(mods, 'potentials', 
                srcdir + '/factory.template.c', 'potentials_factory_c.c',
                srcdir + '/factory.template.h', 'potentials_factory_c.h',
                compiler, system)


mods  = read_module_list('coulomb.classes')
write_1to1_factory_f90(mods, 'coulomb', 'coulomb_factory_f90.f90')
write_factory_c(mods, 'coulomb', 
                srcdir + '/factory.template.c', 'coulomb_factory_c.c',
                srcdir + '/factory.template.h', 'coulomb_factory_c.h',
                compiler, system)


mods  = read_module_list('callables.classes')
write_1ton_factory_f90(mods, 'callables', 'callables_factory_f90.f90')
write_factory_c(mods, 'callables', 
                srcdir + '/factory.template.c', 'callables_factory_c.c',
                srcdir + '/factory.template.h', 'callables_factory_c.h',
                compiler, system)

